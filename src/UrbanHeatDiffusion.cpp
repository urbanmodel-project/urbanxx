// Urban Heat Diffusion Implementation
// 1D heat diffusion dynamics for Urban Surfaces
// Based on ELM SoilTemperature.F90

#include "private/UrbanHeatDiffusionImpl.h"
#include "private/UrbanThermalFunctions.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// Crank-Nicolson factor (same value as CRANK_NICONSON_FACTOR = 0.5)
// Used in building surface boundary condition
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

// Solve tridiagonal system using Thomas algorithm
// Solves: a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = r[i]
// for i = 0 to n-1
// Note: a[0] and c[n-1] are not used
KOKKOS_INLINE_FUNCTION
void SolveTridiagonal(int n, const Real *a, const Real *b, const Real *c,
                      const Real *r, Real *x) {
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

// Solve 1D heat diffusion for a road surface (pervious or impervious)
// using Crank-Nicolson implicit time-stepping
template <typename TempView, typename GeomView, typename TkView,
          typename CvView>
KOKKOS_INLINE_FUNCTION void SolveRoadHeatDiffusion(
    int l, int numLayers, Real dtime, Real capr, TempView &temp, GeomView &zc,
    GeomView &zi, GeomView &dz, TkView &tkInterface, TkView &tkLayer,
    CvView &cv_times_dz, Real EflxGnet, Real DEflxGnet_DTemp) {
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
  const Real dz1 = zc(l, level) - zi(l, level);
  const Real dz2 = zc(l, level + 1) - zi(l, level);
  const Real dz_eff = 0.5 * (dz1 + capr * dz2);

  fact[level] = dtime / cv_times_dz(l, level) * dz(l, level) / dz_eff;
  fn[level] = tkInterface(l, level) * (temp(l, level + 1) - temp(l, level)) /
              (zc(l, level + 1) - zc(l, level));

  // Internal layers
  for (level = 1; level < numLayers - 1; ++level) {
    fact[level] = dtime / cv_times_dz(l, level);
    fn[level] = tkInterface(l, level) * (temp(l, level + 1) - temp(l, level)) /
                (zc(l, level + 1) - zc(l, level));
  }

  // Bottom layer
  level = numLayers - 1;
  fact[level] = dtime / cv_times_dz(l, level);
  fn[level] = 0.0;

  // Compute tridiagonal system (RHS vector and coefficient matrix)
  // Top layer
  level = 0;
  Real dzp, dzm;

  r[level] = temp(l, level) +
             fact[level] * (EflxGnet - DEflxGnet_DTemp * temp(l, level) +
                            CRANK_NICONSON_FACTOR * fn[level]);

  dzp = zc(l, level + 1) - zc(l, level);
  a[level] = 0.0;
  b[level] = 1.0 +
             (1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
                 tkInterface(l, level) / dzp -
             fact[level] * DEflxGnet_DTemp;
  c[level] = -(1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
             tkInterface(l, level) / dzp;

  // Internal layers
  for (level = 1; level < numLayers - 1; ++level) {
    r[level] = temp(l, level) + fact[level] * CRANK_NICONSON_FACTOR *
                                    (fn[level] - fn[level - 1]);
    dzm = zc(l, level) - zc(l, level - 1);
    dzp = zc(l, level + 1) - zc(l, level);
    a[level] = -(1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
               tkInterface(l, level - 1) / dzm;
    b[level] = 1.0 + (1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
                         (tkInterface(l, level) / dzp +
                          tkInterface(l, level - 1) / dzm);
    c[level] = -(1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
               tkInterface(l, level) / dzp;
  }

  // Bottom layer
  level = numLayers - 1;
  r[level] = temp(l, level) -
             CRANK_NICONSON_FACTOR * fact[level - 1] * fn[level] +
             fact[level] * fn[level];
  dzm = zc(l, level) - zc(l, level - 1);
  a[level] = -(1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
             tkLayer(l, level - 1) / dzm;
  b[level] = 1.0 + (1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
                       tkLayer(l, level - 1) / dzm;
  c[level] = 0.0;

  // Solve tridiagonal system for new temperatures
  SolveTridiagonal(numLayers, a, b, c, r, newTemp);

  // Update temperature
  for (int j = 0; j < numLayers; ++j) {
    temp(l, j) = newTemp[j];
  }
}

// Solve 1D heat diffusion for a building surface (roof or wall)
// using Crank-Nicolson implicit time-stepping
// Differs from road surfaces: no dz_eff adjustment for top layer,
// and bottom layer has heat flux to building interior
template <typename TempView, typename GeomView, typename TkView,
          typename CvView>
KOKKOS_INLINE_FUNCTION void SolveBuildingSurfaceHeatDiffusion(
    int l, int numLayers, Real dtime, TempView &temp, GeomView &zc,
    GeomView &zi, GeomView &dz, TkView &tkInterface, TkView &tkLayer,
    CvView &cv_times_dz, Real EflxGnet, Real DEflxGnet_DTemp,
    Real buildingTemp) {
  int level;

  // Working arrays for heat diffusion computation
  Real fact[NUM_URBAN_LAYERS];
  Real fn[NUM_URBAN_LAYERS];
  Real a[NUM_URBAN_LAYERS];
  Real b[NUM_URBAN_LAYERS];
  Real c[NUM_URBAN_LAYERS];
  Real r[NUM_URBAN_LAYERS];
  Real newTemp[NUM_URBAN_LAYERS];

  // Compute factors and interface fluxes
  // Top layer (no dz_eff adjustment for building surfaces)
  level = 0;
  fact[level] = dtime / cv_times_dz(l, level);
  fn[level] = tkInterface(l, level) * (temp(l, level + 1) - temp(l, level)) /
              (zc(l, level + 1) - zc(l, level));

  // Internal layers
  for (level = 1; level < numLayers - 1; ++level) {
    fact[level] = dtime / cv_times_dz(l, level);
    fn[level] = tkInterface(l, level) * (temp(l, level + 1) - temp(l, level)) /
                (zc(l, level + 1) - zc(l, level));
  }

  // Bottom layer (heat flux to building interior)
  level = numLayers - 1;
  fact[level] = dtime / cv_times_dz(l, level);
  fn[level] = tkInterface(l, level) * (buildingTemp - CNFAC * temp(l, level)) /
              (zi(l, level) - zc(l, level));

  // Compute tridiagonal system (RHS vector and coefficient matrix)
  // Top layer
  level = 0;
  Real dzp, dzm;

  r[level] = temp(l, level) +
             fact[level] * (EflxGnet - DEflxGnet_DTemp * temp(l, level) +
                            CRANK_NICONSON_FACTOR * fn[level]);

  dzp = zc(l, level + 1) - zc(l, level);
  a[level] = 0.0;
  b[level] = 1.0 +
             (1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
                 tkInterface(l, level) / dzp -
             fact[level] * DEflxGnet_DTemp;
  c[level] = -(1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
             tkInterface(l, level) / dzp;

  // Internal layers
  for (level = 1; level < numLayers - 1; ++level) {
    r[level] = temp(l, level) + fact[level] * CRANK_NICONSON_FACTOR *
                                    (fn[level] - fn[level - 1]);
    dzm = zc(l, level) - zc(l, level - 1);
    dzp = zc(l, level + 1) - zc(l, level);
    a[level] = -(1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
               tkInterface(l, level - 1) / dzm;
    b[level] = 1.0 + (1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
                         (tkInterface(l, level) / dzp +
                          tkInterface(l, level - 1) / dzm);
    c[level] = -(1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
               tkInterface(l, level) / dzp;
  }

  // Bottom layer (different from roads - uses building boundary condition)
  level = numLayers - 1;
  r[level] = temp(l, level) +
             fact[level] * (fn[level] - CRANK_NICONSON_FACTOR * fn[level - 1]);
  dzm = zc(l, level) - zc(l, level - 1);
  dzp = zi(l, level + 1) - zc(l, level);
  a[level] = -(1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
             tkInterface(l, level - 1) / dzm;
  b[level] =
      1.0 + (1.0 - CRANK_NICONSON_FACTOR) * fact[level] *
                (tkInterface(l, level - 1) / dzm + tkInterface(l, level) / dzp);
  c[level] = 0.0;

  // Solve tridiagonal system for new temperatures
  SolveTridiagonal(numLayers, a, b, c, r, newTemp);

  // Update temperature
  for (int j = 0; j < numLayers; ++j) {
    temp(l, j) = newTemp[j];
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

        // Solve heat diffusion for pervious road
        SolveRoadHeatDiffusion(l, numSoilLayers, dtime, capr, perv_temp,
                               perv_zc, perv_zi, perv_dz, perv_tkInterface,
                               perv_tkLayer, perv_cv_times_dz, perv_EflxGnet,
                               perv_DEflxGnet_DTemp);

        // Solve heat diffusion for impervious road
        SolveRoadHeatDiffusion(
            l, numSoilLayers, dtime, capr, imperv_temp, imperv_zc, imperv_zi,
            imperv_dz, imperv_tkInterface, imperv_tkLayer, imperv_cv_times_dz,
            imperv_EflxGnet, imperv_DEflxGnet_DTemp);

        // Solve heat diffusion for roof
        SolveBuildingSurfaceHeatDiffusion(
            l, numUrbanLayers, dtime, roof_temp, roof_zc, roof_zi, roof_dz,
            roof_tkInterface, roof_tkLayer, roof_cv_times_dz, roof_EflxGnet,
            roof_DEflxGnet_DTemp, building_temp(l));

        // Solve heat diffusion for sunlit wall
        SolveBuildingSurfaceHeatDiffusion(
            l, numUrbanLayers, dtime, sunwall_temp, sunwall_zc, sunwall_zi,
            sunwall_dz, sunwall_tkInterface, sunwall_tkLayer,
            sunwall_cv_times_dz, sunwall_EflxGnet, sunwall_DEflxGnet_DTemp,
            building_temp(l));

        // Solve heat diffusion for shaded wall
        SolveBuildingSurfaceHeatDiffusion(
            l, numUrbanLayers, dtime, shadewall_temp, shadewall_zc,
            shadewall_zi, shadewall_dz, shadewall_tkInterface,
            shadewall_tkLayer, shadewall_cv_times_dz, shadewall_EflxGnet,
            shadewall_DEflxGnet_DTemp, building_temp(l));
      });
  Kokkos::fence();
}

} // namespace URBANXX
