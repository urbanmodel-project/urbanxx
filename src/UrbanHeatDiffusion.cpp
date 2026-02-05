// Urban Heat Diffusion Implementation
// 1D heat diffusion dynamics for Urban Surfaces
// Based on ELM SoilTemperature.F90

#include "private/UrbanHeatDiffusionImpl.h"
#include "private/UrbanThermalFunctions.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// Crank-Nicolson weighting factor for implicit time-stepping
constexpr Real CNFAC = 0.5;

// Compute ground net energy flux for a surface
KOKKOS_INLINE_FUNCTION
Real ComputeGroundNetEnergyFlux(Real netSw, Real netLw, Real eflxShGrnd,
                                Real qflxEvapSoil = 0.0,
                                Real qflxTranEvap = 0.0) {
  return netSw - netLw -
         (eflxShGrnd + qflxEvapSoil * SHR_CONST_LATVAP +
          qflxTranEvap * SHR_CONST_LATVAP);
}

// Compute derivative of ground net energy flux w.r.t. temperature
KOKKOS_INLINE_FUNCTION
Real ComputeGroundNetEnergyFluxDerivative(Real cgrnds, Real cgrndl, Real emiss,
                                          Real temp) {
  const Real cgrnd = cgrnds + cgrndl * SHR_CONST_LATVAP;
  const Real dlwrdDTemp = 4.0 * emiss * STEBOL * Kokkos::pow(temp, 3.0);
  return -cgrnd - dlwrdDTemp;
}

// Compute effective thermal conductivity for pervious road soil layers
// Based on ELM SoilTemperatureMod.F90 (Johansen 1975 model)
template <typename ViewType>
KOKKOS_INLINE_FUNCTION void ComputeSoilThermalConductivity(
    const int l, const int numLayers, const ViewType &tk_minerals,
    const ViewType &tk_dry, const ViewType &tkLayer, const ViewType &watsat,
    const ViewType &water_liquid, const ViewType &water_ice, const ViewType &dz,
    const ViewType &temp) {

  constexpr Real tfrz = SHR_CONST_TKFRZ;
  constexpr Real denh2o = SHR_CONST_RHOWATER;
  constexpr Real denice = SHR_CONST_RHOICE;
  constexpr Real tkwat = TKWATER;
  constexpr Real tkice = TKICE;

  for (int k = 0; k < numLayers; ++k) {
    // Compute degree of saturation
    // satw = (liquid_vol + ice_vol) / porosity
    const Real satw_raw =
        (water_liquid(l, k) / denh2o + water_ice(l, k) / denice) /
        (dz(l, k) * watsat(l, k));
    const Real satw = Kokkos::fmin(1.0, satw_raw);

    if (satw > 1.0e-6) {
      // Compute Kersten number (thermal conductivity enhancement factor)
      Real dke;
      if (temp(l, k) >= tfrz) {
        // Unfrozen soil: Kersten number based on log of saturation
        dke = Kokkos::fmax(0.0, Kokkos::log10(satw) + 1.0);
      } else {
        // Frozen soil: Kersten number equals saturation
        dke = satw;
      }

      // Compute liquid fraction
      const Real liq_vol = water_liquid(l, k) / (denh2o * dz(l, k));
      const Real ice_vol = water_ice(l, k) / (denice * dz(l, k));
      const Real fl = liq_vol / (liq_vol + ice_vol);

      // Saturated thermal conductivity (Johansen model)
      // Geometric mean of components weighted by porosity and phase
      const Real dksat = tk_minerals(l, k) *
                         Kokkos::pow(tkwat, fl * watsat(l, k)) *
                         Kokkos::pow(tkice, (1.0 - fl) * watsat(l, k));

      // Effective thermal conductivity
      // Linear interpolation between dry and saturated based on Kersten number
      tkLayer(l, k) = dke * dksat + (1.0 - dke) * tk_dry(l, k);
    } else {
      // Dry soil
      tkLayer(l, k) = tk_dry(l, k);
    }
    if (k > NUM_LAYERS_ABV_BEDROCK - 1) {
      tkLayer(l, k) = TK_BEDROCK;
    }
  }
}

// Compute heat capacity times layer thickness for pervious road soil layers
// Based on ELM SoilTemperatureMod.F90
// cv(c,j) = csol(c,j)*(1-watsat(c,j))*dz(c,j) + (h2osoi_ice(c,j)*cpice +
// h2osoi_liq(c,j)*cpliq)
template <typename ViewType>
KOKKOS_INLINE_FUNCTION void ComputeSoilHeatCapacityTimesDz(
    const int l, const int numLayers, const ViewType &cv_solids,
    const ViewType &watsat, const ViewType &water_liquid,
    const ViewType &water_ice, const ViewType &dz, const ViewType &cvTimesDz) {

  constexpr Real cpice = SHR_CONST_CPICE; // Specific heat of ice (J/kg/K)
  constexpr Real cpliq =
      SHR_CONST_CPFW; // Specific heat of liquid water (J/kg/K)

  for (int k = 0; k < numLayers; ++k) {
    // Heat capacity of soil solids component
    const Real cv_solid_component =
        cv_solids(l, k) * (1.0 - watsat(l, k)) * dz(l, k);

    // Heat capacity of water components (ice + liquid)
    const Real cv_water_component =
        water_ice(l, k) * cpice + water_liquid(l, k) * cpliq;

    // Total heat capacity times layer thickness
    if (k < 10) {
      cvTimesDz(l, k) = cv_solid_component + cv_water_component;
    } else {
      cvTimesDz(l, k) = cv_solid_component;
    }
  }
}

// Struct to hold geometry and thermal properties for heat diffusion
template <typename TempView, typename GeomView, typename TkView,
          typename CvView>
struct SurfaceProperties {
  int l;               // Landunit index
  int numLayers;       // Number of layers
  TempView &temp;      // Temperature
  GeomView &zc;        // Layer center depth
  GeomView &zi;        // Layer interface depth
  GeomView &dz;        // Layer thickness
  TkView &tkInterface; // Thermal conductivity at interfaces
  TkView &tkLayer;     // Thermal conductivity at layer centers
  CvView &cv_times_dz; // Heat capacity times layer thickness

  KOKKOS_INLINE_FUNCTION
  SurfaceProperties(int landunit, int layers, TempView &temperature,
                    GeomView &layer_center, GeomView &layer_interface,
                    GeomView &layer_thickness, TkView &tk_interface,
                    TkView &tk_layer, CvView &cv_dz)
      : l(landunit), numLayers(layers), temp(temperature), zc(layer_center),
        zi(layer_interface), dz(layer_thickness), tkInterface(tk_interface),
        tkLayer(tk_layer), cv_times_dz(cv_dz) {}
};

// Struct to hold boundary condition parameters
struct BoundaryConditions {
  Real EflxGnet;              // Ground net energy flux
  Real DEflxGnet_DTemp;       // Derivative of flux w.r.t. temperature
  bool useTopLayerAdjustment; // Use dz_eff adjustment for top layer
  Real capr;               // Turing factor (if useTopLayerAdjustment is true)
  bool hasBottomBoundary;  // Has bottom boundary temperature coupling
  Real bottomBoundaryTemp; // Bottom boundary temperature (if hasBottomBoundary
                           // is true)

  KOKKOS_INLINE_FUNCTION
  BoundaryConditions(Real net_flux, Real flux_derivative, bool top_adjustment,
                     Real turing_factor, bool bottom_bc, Real bottom_temp)
      : EflxGnet(net_flux), DEflxGnet_DTemp(flux_derivative),
        useTopLayerAdjustment(top_adjustment), capr(turing_factor),
        hasBottomBoundary(bottom_bc), bottomBoundaryTemp(bottom_temp) {}
};

// Solve tridiagonal system using Thomas algorithm
// Solves: a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = r[i]
// for i = 0 to n-1
// Note: a[0] and c[n-1] are not used
KOKKOS_INLINE_FUNCTION
void SolveTridiagonal(int n, const Real *a, const Real *b, const Real *c,
                      const Real *r, Real *x) {
  KOKKOS_ASSERT(n <= NUM_SOIL_LAYERS &&
                "SolveTridiagonal: n exceeds NUM_SOIL_LAYERS");

  // Working arrays for modified coefficients
  Real cp[NUM_SOIL_LAYERS]; // Modified upper diagonal
  Real rp[NUM_SOIL_LAYERS]; // Modified right-hand side

  // Forward elimination
  cp[0] = c[0] / b[0];
  rp[0] = r[0] / b[0];

  for (int i = 1; i < n; i++) {
    const Real denom = b[i] - a[i] * cp[i - 1];
    cp[i] = c[i] / denom;
    rp[i] = (r[i] - a[i] * rp[i - 1]) / denom;
  }

  // Back substitution
  x[n - 1] = rp[n - 1];
  for (int i = n - 2; i >= 0; i--) {
    x[i] = rp[i] - cp[i] * x[i + 1];
  }
}

// Solve 1D heat diffusion using Crank-Nicolson implicit time-stepping
// Handles both road surfaces (with optional top layer adjustment and zero flux
// bottom BC) and building surfaces (no top adjustment, heat flux to building
// interior at bottom)
template <typename TempView, typename GeomView, typename TkView,
          typename CvView>
KOKKOS_INLINE_FUNCTION void Solve1DHeatDiffusion(
    Real dtime, SurfaceProperties<TempView, GeomView, TkView, CvView> &surf,
    const BoundaryConditions &bc) {
  int level;

  // Working arrays for heat diffusion computation
  Real fact[NUM_SOIL_LAYERS];
  Real fn[NUM_SOIL_LAYERS];
  Real a[NUM_SOIL_LAYERS];
  Real b[NUM_SOIL_LAYERS];
  Real c[NUM_SOIL_LAYERS];
  Real r[NUM_SOIL_LAYERS];
  Real newTemp[NUM_SOIL_LAYERS];

  // Compute factors and interface fluxes
  // Top layer
  level = 0;
  if (bc.useTopLayerAdjustment) {
    // Road surfaces: use dz_eff adjustment (Turing factor)
    const Real dz1 = surf.zc(surf.l, level) - surf.zi(surf.l, level);
    const Real dz2 = surf.zc(surf.l, level + 1) - surf.zi(surf.l, level);
    const Real dz_eff = 0.5 * (dz1 + bc.capr * dz2);
    fact[level] = dtime / surf.cv_times_dz(surf.l, level) *
                  surf.dz(surf.l, level) / dz_eff;
  } else {
    // Building surfaces: direct calculation
    fact[level] = dtime / surf.cv_times_dz(surf.l, level);
  }
  fn[level] = surf.tkInterface(surf.l, level) *
              (surf.temp(surf.l, level + 1) - surf.temp(surf.l, level)) /
              (surf.zc(surf.l, level + 1) - surf.zc(surf.l, level));

  // Internal layers
  for (level = 1; level < surf.numLayers - 1; ++level) {
    fact[level] = dtime / surf.cv_times_dz(surf.l, level);
    fn[level] = surf.tkInterface(surf.l, level) *
                (surf.temp(surf.l, level + 1) - surf.temp(surf.l, level)) /
                (surf.zc(surf.l, level + 1) - surf.zc(surf.l, level));
  }

  // Bottom layer
  level = surf.numLayers - 1;
  fact[level] = dtime / surf.cv_times_dz(surf.l, level);
  if (bc.hasBottomBoundary) {
    // Building surfaces: heat flux to building interior
    fn[level] = surf.tkInterface(surf.l, level) *
                (bc.bottomBoundaryTemp - surf.temp(surf.l, level)) /
                (surf.zi(surf.l, level + 1) - surf.zc(surf.l, level));
  } else {
    // Road surfaces: zero flux boundary condition
    fn[level] = 0.0;
  }

  // Compute tridiagonal system (RHS vector and coefficient matrix)
  // Top layer
  level = 0;
  Real dzp, dzm;

  r[level] = surf.temp(surf.l, level) +
             fact[level] *
                 (bc.EflxGnet - bc.DEflxGnet_DTemp * surf.temp(surf.l, level) +
                  CNFAC * fn[level]);

  dzp = surf.zc(surf.l, level + 1) - surf.zc(surf.l, level);
  a[level] = 0.0;
  b[level] =
      1.0 +
      (1.0 - CNFAC) * fact[level] * surf.tkInterface(surf.l, level) / dzp -
      fact[level] * bc.DEflxGnet_DTemp;
  c[level] =
      -(1.0 - CNFAC) * fact[level] * surf.tkInterface(surf.l, level) / dzp;

  // Internal layers
  for (level = 1; level < surf.numLayers - 1; ++level) {
    r[level] = surf.temp(surf.l, level) +
               fact[level] * CNFAC * (fn[level] - fn[level - 1]);
    dzm = surf.zc(surf.l, level) - surf.zc(surf.l, level - 1);
    dzp = surf.zc(surf.l, level + 1) - surf.zc(surf.l, level);
    a[level] = -(1.0 - CNFAC) * fact[level] *
               surf.tkInterface(surf.l, level - 1) / dzm;
    b[level] = 1.0 + (1.0 - CNFAC) * fact[level] *
                         (surf.tkInterface(surf.l, level) / dzp +
                          surf.tkInterface(surf.l, level - 1) / dzm);
    c[level] =
        -(1.0 - CNFAC) * fact[level] * surf.tkInterface(surf.l, level) / dzp;
  }

  // Bottom layer
  level = surf.numLayers - 1;
  dzm = surf.zc(surf.l, level) - surf.zc(surf.l, level - 1);
  r[level] = surf.temp(surf.l, level) +
             fact[level] * CNFAC * (fn[level] - fn[level - 1]);
  a[level] =
      -(1.0 - CNFAC) * fact[level] * surf.tkInterface(surf.l, level - 1) / dzm;

  if (bc.hasBottomBoundary) {
    // Building surfaces: include coupling to building interior
    dzp = surf.zi(surf.l, level + 1) - surf.zc(surf.l, level);
    r[level] += (1.0 - CNFAC) * fact[level] * surf.tkInterface(surf.l, level) /
                dzp * bc.bottomBoundaryTemp;
    b[level] = 1.0 + (1.0 - CNFAC) * fact[level] *
                         (surf.tkInterface(surf.l, level - 1) / dzm +
                          surf.tkInterface(surf.l, level) / dzp);
  } else {
    // Road surfaces: zero flux at bottom
    b[level] = 1.0 + (1.0 - CNFAC) * fact[level] *
                         surf.tkInterface(surf.l, level - 1) / dzm;
  }
  c[level] = 0.0;

  // Solve tridiagonal system for new temperatures
  SolveTridiagonal(surf.numLayers, a, b, c, r, newTemp);

  // Update temperature
  for (int j = 0; j < surf.numLayers; ++j) {
    surf.temp(surf.l, j) = newTemp[j];
  }
}

// Compute 1D heat diffusion for all urban surfaces
void ComputeHeatDiffusion(URBANXX::_p_UrbanType &urban) {

  const int numLandunits = urban.numLandunits;
  const int numSoilLayers = urban.numSoilLayers;
  const int numUrbanLayers = urban.numUrbanLayers;

  // Access pervious road soil property views
  auto &perv_tk_minerals = urban.perviousRoad.soil.TkMinerals;
  auto &perv_tk_dry = urban.perviousRoad.soil.TkDry;
  auto &perv_tkLayer = urban.perviousRoad.TkLayer;
  auto &perv_tkInterface = urban.perviousRoad.TkInterface;
  auto &perv_watsat = urban.perviousRoad.soil.WatSat;
  auto &perv_water_liquid = urban.perviousRoad.soil.WaterLiquid;
  auto &perv_water_ice = urban.perviousRoad.soil.WaterIce;
  auto &perv_dz = urban.perviousRoad.Dz;
  auto &perv_zc = urban.perviousRoad.Zc;
  auto &perv_zi = urban.perviousRoad.Zi;
  auto &perv_temp = urban.perviousRoad.Temperature;
  auto &perv_cv_solids = urban.perviousRoad.soil.CvSolids;
  auto &perv_cv_times_dz = urban.perviousRoad.CvTimesDz;
  auto &perv_netLw = urban.perviousRoad.NetLongRad;
  auto &perv_netSw = urban.perviousRoad.NetShortRad;
  auto &perv_EflxShGrnd = urban.perviousRoad.EflxShGrnd;
  auto &perv_QflxEvapSoil = urban.perviousRoad.QflxEvapSoil;
  auto &perv_QflxTranEvap = urban.perviousRoad.QflxTranEvap;
  auto &perv_Cgrndl = urban.perviousRoad.Cgrndl;
  auto &perv_Cgrnds = urban.perviousRoad.Cgrnds;
  auto &perv_Temp = urban.perviousRoad.EffectiveSurfTemp;
  auto &perv_emiss = urban.urbanParams.emissivity.PerviousRoad;

  // Access roof property views
  auto &roof_netLw = urban.roof.NetLongRad;
  auto &roof_netSw = urban.roof.NetShortRad;
  auto &roof_EflxShGrnd = urban.roof.EflxShGrnd;
  auto &roof_QflxEvapSoil = urban.roof.QflxEvapSoil;
  auto &roof_QflxTranEvap = urban.roof.QflxTranEvap;
  auto &roof_Cgrnds = urban.roof.Cgrnds;
  auto &roof_Cgrndl = urban.roof.Cgrndl;
  auto &roof_Temp = urban.roof.EffectiveSurfTemp;
  auto &roof_emiss = urban.urbanParams.emissivity.Roof;
  auto &roof_temp = urban.roof.Temperature;
  auto &roof_tkInterface = urban.roof.TkInterface;
  auto &roof_tkLayer = urban.roof.TkLayer;
  auto &roof_cv_times_dz = urban.roof.CvTimesDz;
  auto &roof_dz = urban.roof.Dz;
  auto &roof_zc = urban.roof.Zc;
  auto &roof_zi = urban.roof.Zi;

  // Access impervious road property views
  auto &imperv_tkLayer = urban.imperviousRoad.TkLayer;
  auto &imperv_tkInterface = urban.imperviousRoad.TkInterface;
  auto &imperv_dz = urban.imperviousRoad.Dz;
  auto &imperv_zc = urban.imperviousRoad.Zc;
  auto &imperv_zi = urban.imperviousRoad.Zi;
  auto &imperv_temp = urban.imperviousRoad.Temperature;
  auto &imperv_cv_times_dz = urban.imperviousRoad.CvTimesDz;
  auto &imperv_netLw = urban.imperviousRoad.NetLongRad;
  auto &imperv_netSw = urban.imperviousRoad.NetShortRad;
  auto &imperv_EflxShGrnd = urban.imperviousRoad.EflxShGrnd;
  auto &imperv_QflxEvapSoil = urban.imperviousRoad.QflxEvapSoil;
  auto &imperv_QflxTranEvap = urban.imperviousRoad.QflxTranEvap;
  auto &imperv_Cgrnds = urban.imperviousRoad.Cgrnds;
  auto &imperv_Cgrndl = urban.imperviousRoad.Cgrndl;
  auto &imperv_Temp = urban.imperviousRoad.EffectiveSurfTemp;
  auto &imperv_emiss = urban.urbanParams.emissivity.ImperviousRoad;

  // Access sunlit wall property views
  auto &sunwall_netLw = urban.sunlitWall.NetLongRad;
  auto &sunwall_netSw = urban.sunlitWall.NetShortRad;
  auto &sunwall_EflxShGrnd = urban.sunlitWall.EflxShGrnd;
  auto &sunwall_Cgrnds = urban.sunlitWall.Cgrnds;
  auto &sunwall_Cgrndl = urban.sunlitWall.Cgrndl;
  auto &sunwall_Temp = urban.sunlitWall.EffectiveSurfTemp;
  auto &sunwall_temp = urban.sunlitWall.Temperature;
  auto &sunwall_tkInterface = urban.sunlitWall.TkInterface;
  auto &sunwall_tkLayer = urban.sunlitWall.TkLayer;
  auto &sunwall_cv_times_dz = urban.sunlitWall.CvTimesDz;
  auto &sunwall_dz = urban.sunlitWall.Dz;
  auto &sunwall_zc = urban.sunlitWall.Zc;
  auto &sunwall_zi = urban.sunlitWall.Zi;
  auto &wall_emiss = urban.urbanParams.emissivity.Wall;

  // Access shaded wall property views
  auto &shadewall_netLw = urban.shadedWall.NetLongRad;
  auto &shadewall_netSw = urban.shadedWall.NetShortRad;
  auto &shadewall_EflxShGrnd = urban.shadedWall.EflxShGrnd;
  auto &shadewall_Cgrnds = urban.shadedWall.Cgrnds;
  auto &shadewall_Cgrndl = urban.shadedWall.Cgrndl;
  auto &shadewall_Temp = urban.shadedWall.EffectiveSurfTemp;
  auto &shadewall_temp = urban.shadedWall.Temperature;
  auto &shadewall_tkInterface = urban.shadedWall.TkInterface;
  auto &shadewall_tkLayer = urban.shadedWall.TkLayer;
  auto &shadewall_cv_times_dz = urban.shadedWall.CvTimesDz;
  auto &shadewall_dz = urban.shadedWall.Dz;
  auto &shadewall_zc = urban.shadedWall.Zc;
  auto &shadewall_zi = urban.shadedWall.Zi;

  // Access building property views
  auto &building_temp = urban.building.Temperature;

  // Single parallel kernel over all landunits
  Kokkos::parallel_for(
      "ComputeHeatDiffusion", numLandunits, KOKKOS_LAMBDA(int l) {
        // Step 1: Compute thermal conductivity for pervious road soil layers
        ComputeSoilThermalConductivity(
            l, numSoilLayers, perv_tk_minerals, perv_tk_dry, perv_tkLayer,
            perv_watsat, perv_water_liquid, perv_water_ice, perv_dz, perv_temp);

        // Step 2: Compute thermal conductivity at layer interfaces
        // (Note: Interface thermal conductivity for impervious road, walls, and
        // roof
        //  are computed during initialization and don't need to be recomputed)

        // Pervious road: all soil layers are active
        ComputeInterfaceThermalConductivity(l, numSoilLayers, numSoilLayers,
                                            perv_tkLayer, perv_tkInterface,
                                            perv_zc, perv_zi, 0.0);

        // Step 3: Compute heat capacity times layer thickness for pervious road
        ComputeSoilHeatCapacityTimesDz(
            l, numSoilLayers, perv_cv_solids, perv_watsat, perv_water_liquid,
            perv_water_ice, perv_dz, perv_cv_times_dz);

        // Step 4: Compute ground net energy flux and its derivative for all
        // surfaces

        // Pervious road (with evaporation)
        const Real perv_EflxGnet = ComputeGroundNetEnergyFlux(
            perv_netSw(l), perv_netLw(l), perv_EflxShGrnd(l),
            perv_QflxEvapSoil(l), perv_QflxTranEvap(l));
        const Real perv_DEflxGnet_DTemp = ComputeGroundNetEnergyFluxDerivative(
            perv_Cgrnds(l), perv_Cgrndl(l), perv_emiss(l), perv_Temp(l));

        // Roof (with evaporation)
        const Real roof_EflxGnet = ComputeGroundNetEnergyFlux(
            roof_netSw(l), roof_netLw(l), roof_EflxShGrnd(l),
            roof_QflxEvapSoil(l), roof_QflxTranEvap(l));
        const Real roof_DEflxGnet_DTemp = ComputeGroundNetEnergyFluxDerivative(
            roof_Cgrnds(l), roof_Cgrndl(l), roof_emiss(l), roof_Temp(l));

        // Impervious road (with evaporation)
        const Real imperv_EflxGnet = ComputeGroundNetEnergyFlux(
            imperv_netSw(l), imperv_netLw(l), imperv_EflxShGrnd(l),
            imperv_QflxEvapSoil(l), imperv_QflxTranEvap(l));
        const Real imperv_DEflxGnet_DTemp =
            ComputeGroundNetEnergyFluxDerivative(
                imperv_Cgrnds(l), imperv_Cgrndl(l), imperv_emiss(l),
                imperv_Temp(l));

        // Sunlit wall (no evaporation - defaults to 0)
        const Real sunwall_EflxGnet = ComputeGroundNetEnergyFlux(
            sunwall_netSw(l), sunwall_netLw(l), sunwall_EflxShGrnd(l));
        const Real sunwall_DEflxGnet_DTemp =
            ComputeGroundNetEnergyFluxDerivative(
                sunwall_Cgrnds(l), sunwall_Cgrndl(l), wall_emiss(l),
                sunwall_Temp(l));

        // Shaded wall (no evaporation - defaults to 0)
        const Real shadewall_EflxGnet = ComputeGroundNetEnergyFlux(
            shadewall_netSw(l), shadewall_netLw(l), shadewall_EflxShGrnd(l));
        const Real shadewall_DEflxGnet_DTemp =
            ComputeGroundNetEnergyFluxDerivative(
                shadewall_Cgrnds(l), shadewall_Cgrndl(l), wall_emiss(l),
                shadewall_Temp(l));

        const Real dtime = 30.0 * 60.0; // seconds
        const Real capr =
            0.34; // Turing factor to turn first layer T into surface T

        // Boundary condition flags for different surface types
        const bool useTopAdjustment_road = true;
        const bool useTopAdjustment_building = false;
        const bool hasBottomBC_road = false;
        const bool hasBottomBC_building = true;
        const Real noBottomTemp = 0.0; // Not used for road surfaces

        // Solve heat diffusion for pervious road
        SurfaceProperties perv_surf(l, numSoilLayers, perv_temp, perv_zc,
                                    perv_zi, perv_dz, perv_tkInterface,
                                    perv_tkLayer, perv_cv_times_dz);
        BoundaryConditions perv_bc(perv_EflxGnet, perv_DEflxGnet_DTemp,
                                   useTopAdjustment_road, capr,
                                   hasBottomBC_road, noBottomTemp);
        Solve1DHeatDiffusion(dtime, perv_surf, perv_bc);

        // Solve heat diffusion for impervious road
        SurfaceProperties imperv_surf(l, numSoilLayers, imperv_temp, imperv_zc,
                                      imperv_zi, imperv_dz, imperv_tkInterface,
                                      imperv_tkLayer, imperv_cv_times_dz);
        BoundaryConditions imperv_bc(imperv_EflxGnet, imperv_DEflxGnet_DTemp,
                                     useTopAdjustment_road, capr,
                                     hasBottomBC_road, noBottomTemp);
        Solve1DHeatDiffusion(dtime, imperv_surf, imperv_bc);

        // Solve heat diffusion for roof
        SurfaceProperties roof_surf(l, numUrbanLayers, roof_temp, roof_zc,
                                    roof_zi, roof_dz, roof_tkInterface,
                                    roof_tkLayer, roof_cv_times_dz);
        BoundaryConditions roof_bc(roof_EflxGnet, roof_DEflxGnet_DTemp,
                                   useTopAdjustment_building, capr,
                                   hasBottomBC_building, building_temp(l));
        Solve1DHeatDiffusion(dtime, roof_surf, roof_bc);

        // Solve heat diffusion for sunlit wall
        SurfaceProperties sunwall_surf(
            l, numUrbanLayers, sunwall_temp, sunwall_zc, sunwall_zi, sunwall_dz,
            sunwall_tkInterface, sunwall_tkLayer, sunwall_cv_times_dz);
        BoundaryConditions sunwall_bc(sunwall_EflxGnet, sunwall_DEflxGnet_DTemp,
                                      useTopAdjustment_building, capr,
                                      hasBottomBC_building, building_temp(l));
        Solve1DHeatDiffusion(dtime, sunwall_surf, sunwall_bc);

        // Solve heat diffusion for shaded wall
        SurfaceProperties shadewall_surf(
            l, numUrbanLayers, shadewall_temp, shadewall_zc, shadewall_zi,
            shadewall_dz, shadewall_tkInterface, shadewall_tkLayer,
            shadewall_cv_times_dz);
        BoundaryConditions shadewall_bc(shadewall_EflxGnet,
                                        shadewall_DEflxGnet_DTemp,
                                        useTopAdjustment_building, capr,
                                        hasBottomBC_building, building_temp(l));
        Solve1DHeatDiffusion(dtime, shadewall_surf, shadewall_bc);
      });
  Kokkos::fence();
}

} // namespace URBANXX
