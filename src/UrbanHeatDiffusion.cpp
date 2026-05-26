// Urban Heat Diffusion Implementation
// 1D heat diffusion dynamics for Urban Surfaces
// Based on ELM SoilTemperature.F90

#include "Urban.h"
#include "private/UrbanHeatDiffusionImpl.h"
#include "private/UrbanThermalFunctions.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

using namespace URBANXX;

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
    // For bedrock layers (k >= NUM_LAYERS_ABV_BEDROCK), water content is zero
    // so these terms vanish, but we include them for formula consistency with
    // ELM.
    const Real cv_water_component =
        water_ice(l, k) * cpice + water_liquid(l, k) * cpliq;

    // Total heat capacity times layer thickness
    cvTimesDz(l, k) = cv_solid_component + cv_water_component;
  }
}

// Update heat capacity times layer thickness for impervious road deep soil
// layers (layers beyond the active road construction layers).  Those layers
// are natural soil whose water content is maintained as persistent URBANxx
// state (WaterLiquid/WaterIce), initialized once from ELM at startup and
// updated each timestep via internal phase change.
// Matches ELM's cv formula: cv = csol*(1-watsat)*dz + ice*cpice + liq*cpliq
//
// Active road construction layers (k < numActiveLayers) keep their fixed
// cv_road_params*dz value and are NOT touched here.
// Bedrock layers (k >= NUM_LAYERS_ABV_BEDROCK) also keep their fixed value.
template <typename ViewType, typename IntView>
KOKKOS_INLINE_FUNCTION void UpdateImperviousRoadHeatCapacityTimesDz(
    const int l, const int numLayers, const IntView &numActiveLayers,
    const ViewType &cv_solids, const ViewType &watsat,
    const ViewType &water_liquid, const ViewType &water_ice, const ViewType &dz,
    const ViewType &cvTimesDz) {

  constexpr Real cpice = SHR_CONST_CPICE;
  constexpr Real cpliq = SHR_CONST_CPFW;

  const int nActive = numActiveLayers(l);
  for (int k = nActive; k < NUM_LAYERS_ABV_BEDROCK; ++k) {
    const Real cv_solid_component =
        cv_solids(l, k) * (1.0 - watsat(l, k)) * dz(l, k);
    const Real cv_water_component =
        water_ice(l, k) * cpice + water_liquid(l, k) * cpliq;
    cvTimesDz(l, k) = cv_solid_component + cv_water_component;
  }
}

// Update thermal conductivity (layer-center and interface) for impervious road
// deep soil layers (k = nActive..NUM_LAYERS_ABV_BEDROCK-1).
//
// ELM uses moisture-dependent Johansen conductivity for these layers every
// timestep.  URBANxx was only setting dry conductivity at initialization.
// This function applies the same Johansen formula as
// ComputeSoilThermalConductivity and then recomputes the affected interface
// conductivities.
template <typename ViewType, typename IntView>
KOKKOS_INLINE_FUNCTION void UpdateImperviousRoadThermalProperties(
    const int l, const int numLayers, const IntView &numActiveLayers,
    const ViewType &tk_minerals, const ViewType &tk_dry,
    const ViewType &water_liquid, const ViewType &water_ice,
    const ViewType &watsat, const ViewType &dz, const ViewType &temp,
    const ViewType &zc, const ViewType &zi, const ViewType &tkLayer,
    const ViewType &tkInterface) {

  constexpr Real tfrz = SHR_CONST_TKFRZ;
  constexpr Real denh2o = SHR_CONST_RHOWATER;
  constexpr Real denice = SHR_CONST_RHOICE;
  constexpr Real tkwat = TKWATER;
  constexpr Real tkice = TKICE;

  const int nActive = numActiveLayers(l);

  // --- Step A: Update TkLayer for deep natural-soil layers ---
  for (int k = nActive; k < NUM_LAYERS_ABV_BEDROCK; ++k) {
    const Real satw_raw =
        (water_liquid(l, k) / denh2o + water_ice(l, k) / denice) /
        (dz(l, k) * watsat(l, k));
    const Real satw = Kokkos::fmin(1.0, satw_raw);

    if (satw > 1.0e-6) {
      Real dke;
      if (temp(l, k) >= tfrz) {
        dke = Kokkos::fmax(0.0, Kokkos::log10(satw) + 1.0);
      } else {
        dke = satw;
      }
      const Real liq_vol = water_liquid(l, k) / (denh2o * dz(l, k));
      const Real ice_vol = water_ice(l, k) / (denice * dz(l, k));
      const Real fl = liq_vol / (liq_vol + ice_vol);
      const Real dksat = tk_minerals(l, k) *
                         Kokkos::pow(tkwat, fl * watsat(l, k)) *
                         Kokkos::pow(tkice, (1.0 - fl) * watsat(l, k));
      tkLayer(l, k) = dke * dksat + (1.0 - dke) * tk_dry(l, k);
    } else {
      tkLayer(l, k) = tk_dry(l, k);
    }
  }

  // --- Step B: Recompute TkInterface for affected interfaces ---
  // Interface at k connects layer k and k+1.
  // Layers k=nActive..NUM_LAYERS_ABV_BEDROCK-1 changed, so interfaces
  // k=(nActive-1)..NUM_LAYERS_ABV_BEDROCK-1 must be recomputed.
  // tkLayer for k >= NUM_LAYERS_ABV_BEDROCK is TK_BEDROCK (set at init,
  // unchanged), so the bedrock side of the boundary is already correct.
  const int k_start = (nActive > 0) ? nActive - 1 : 0;
  for (int k = k_start; k < NUM_LAYERS_ABV_BEDROCK; ++k) {
    if (k < numLayers - 1) {
      tkInterface(l, k) = tkLayer(l, k) * tkLayer(l, k + 1) *
                          (zc(l, k + 1) - zc(l, k)) /
                          (tkLayer(l, k) * (zc(l, k + 1) - zi(l, k + 1)) +
                           tkLayer(l, k + 1) * (zi(l, k + 1) - zc(l, k)));
    } else {
      tkInterface(l, k) = tkLayer(l, k);
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

  // initialize arrays
  for (int j = 0; j < NUM_SOIL_LAYERS; ++j) {
    a[j] = 0.0;
    b[j] = 1.0;
    c[j] = 0.0;
    r[j] = 0.0;
  }

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

// Apply thin-snow phase change (ELM Phasechange_beta, snl=0, h2osno>0 path)
// for roof, imperviousRoad, and perviousRoad top layer.
// Mirrors ELM: if T[0] > tfrz and H2OSno > 0, clamp T[0] to tfrz,
// compute energy available for melting, reduce h2osno (and snow_depth),
// then if energy remains melt top-layer ice, then apply residual temperature
// correction.
// useTopLayerAdjustment: false for roof (building surface), true for roads.
template <typename TempView, typename ScalarView1D, typename CvView,
          typename GeomView>
KOKKOS_INLINE_FUNCTION void ApplyThinSnowPhaseChange(
    int l, TempView &temp, ScalarView1D &h2osno, ScalarView1D &snow_depth,
    Real &top_wice, Real &top_wliq, const CvView &cv_times_dz,
    const GeomView &dz, const GeomView &zc, const GeomView &zi, Real cgrnds,
    Real cgrndl, Real emiss, Real tgnd0, bool useTopLayerAdjustment, Real capr,
    Real dtime) {

  constexpr Real tfrz = SHR_CONST_TKFRZ;
  constexpr Real hfus = SHR_CONST_LATICE;

  // No thin-snow phase change if no snow water or temperature already at/below
  // freezing
  if (h2osno(l) <= 0.0 || temp(l, 0) <= tfrz)
    return;

  // Clamp temperature to freezing point; tinc is the temperature decrement
  // (negative, since T was above tfrz)
  const Real tinc = tfrz - temp(l, 0); // < 0
  temp(l, 0) = tfrz;

  // Compute fact_top (ELM's fact(c,j) for the top layer)
  Real fact_top;
  if (useTopLayerAdjustment) {
    // Road surfaces: Turing factor (ELM SoilTemperatureMod fact formula)
    const Real dz1 = zc(l, 0) - zi(l, 0);
    const Real dz2 = zc(l, 1) - zi(l, 0);
    const Real dz_eff = 0.5 * (dz1 + capr * dz2);
    fact_top = dtime / cv_times_dz(l, 0) * dz(l, 0) / dz_eff;
  } else {
    // Building surfaces (roof): direct ratio
    fact_top = dtime / cv_times_dz(l, 0);
  }

  // dhsdT = -(cgrnds + cgrndl*hvap + 4*emiss*STEBOL*tgnd0^3)
  // (matches ComputeGroundNetEnergyFluxDerivative)
  const Real cgrnd = cgrnds + cgrndl * SHR_CONST_LATVAP;
  const Real dhsdT = -cgrnd - 4.0 * emiss * STEBOL * Kokkos::pow(tgnd0, 3.0);

  // Energy available for melting (j == snl+1 case, j > 0 in ELM indexing)
  // hm = dhsdT*tinc - tinc/fact_top  (> 0 when T > tfrz)
  const Real hm = dhsdT * tinc - tinc / fact_top;
  if (hm <= 0.0)
    return;

  const Real xm = hm * dtime / hfus;

  // Melt h2osno first
  const Real h2osno_old = h2osno(l);
  h2osno(l) = Kokkos::fmax(0.0, h2osno_old - xm);
  if (h2osno_old > 0.0) {
    const Real propor = h2osno(l) / h2osno_old;
    snow_depth(l) = propor * snow_depth(l);
  }
  // Energy remaining after melting h2osno
  Real heatr = hm - hfus * (h2osno_old - h2osno(l)) / dtime;

  if (heatr <= 0.0)
    return; // All energy consumed; T stays at tfrz

  // Remaining energy: melt top-layer soil ice
  const Real xm_rem = heatr * dtime / hfus;
  const Real wice0 = top_wice;
  const Real wmass0 = top_wice + top_wliq;
  top_wice = Kokkos::fmax(0.0, wice0 - xm_rem);
  top_wliq = Kokkos::fmax(0.0, wmass0 - top_wice);
  const Real heatr2 = heatr - hfus * (wice0 - top_wice) / dtime;

  // Residual temperature correction (j == snl+1, j > 0, frac_h2osfc=0)
  if (Kokkos::fabs(heatr2) > 0.0) {
    const Real denom = 1.0 - fact_top * dhsdT;
    if (Kokkos::fabs(denom) > 1.0e-30) {
      temp(l, 0) = tfrz + fact_top * heatr2 / denom;
    }
  }
}

// Apply soil phase change for pervious road layers (all layers k=0..n-1).
// Mirrors ELM Phasechange_beta for is_soil / icol_road_perv columns.
// Handles both melt (T > tfrz, wice > 0) and freeze (T < tfrz, wliq >
// supercool). For k=0 when thin-snow phase change already clamped T to tfrz,
// imelt is not triggered (T == tfrz → no change).
template <typename TempView, typename WaterView, typename CvView,
          typename SoilView, typename DzView>
KOKKOS_INLINE_FUNCTION void ApplyPerviousRoadSoilPhaseChange(
    int l, int numSoilLayers, TempView &temp, WaterView &water_ice,
    WaterView &water_liq, const CvView &cv_times_dz, const SoilView &watsat,
    const SoilView &bsw, const SoilView &sucsat, const DzView &dz,
    const DzView &zc, const DzView &zi, Real cpar, Real cgrnds, Real cgrndl,
    Real emiss, Real tgnd0, Real dtime) {

  constexpr Real tfrz = SHR_CONST_TKFRZ;
  constexpr Real hfus = SHR_CONST_LATICE;
  constexpr Real grav = GRAVITY;

  for (int k = 0; k < numSoilLayers; ++k) {
    const Real T = temp(l, k);
    const Real wice = water_ice(l, k);
    const Real wliq = water_liq(l, k);
    const Real wmass0 = wice + wliq;
    Real fact_k = dtime / cv_times_dz(l, k);
    if (k == 0) {
      const Real dz1 = zc(l, 0) - zi(l, 0);
      const Real dz2 = zc(l, 1) - zi(l, 0);
      const Real dz_eff = 0.5 * (dz1 + cpar * dz2);
      fact_k = fact_k * dz(l, 0) / dz_eff;
    }

    // Supercool threshold (ELM: for is_soil and icol_road_perv)
    Real supercool = 0.0;
    if (T < tfrz) {
      const Real smp = hfus * (tfrz - T) / (grav * T) * 1000.0; // [mm H2O head]
      if (smp > 0.0 && sucsat(l, k) > 0.0 && bsw(l, k) > 0.0) {
        const Real sc = watsat(l, k) *
                        Kokkos::pow(smp / sucsat(l, k), -1.0 / bsw(l, k)) *
                        dz(l, k) * 1000.0; // [kg/m²]
        supercool = Kokkos::fmax(0.0, sc);
      }
    }

    int imelt = 0;
    Real tinc = 0.0;

    if (wice > 0.0 && T > tfrz) {
      imelt = 1;
      tinc = tfrz - T; // negative
      temp(l, k) = tfrz;
    } else if (wliq > supercool && T < tfrz) {
      imelt = 2;
      tinc = tfrz - T; // positive
      temp(l, k) = tfrz;
    }

    if (imelt == 0)
      continue;

    // Energy: internal layers use hm = -tinc/fact.
    // For top layer (k=0), include dhsdT*tinc, matching ELM top-layer path.
    Real hm;
    if (k == 0) {
      const Real cgrnd = cgrnds + cgrndl * SHR_CONST_LATVAP;
      const Real dhsdT =
          -cgrnd - 4.0 * emiss * STEBOL * Kokkos::pow(tgnd0, 3.0);
      hm = dhsdT * tinc - tinc / fact_k;
    } else {
      hm = -tinc / fact_k;
    }

    // Guard against sign errors (ELM does the same check)
    if (imelt == 1 && hm < 0.0) {
      temp(l, k) = tfrz + tinc;
      continue;
    }
    if (imelt == 2 && hm > 0.0) {
      temp(l, k) = tfrz + tinc;
      continue;
    }

    const Real xm = hm * dtime / hfus;

    Real wice_new;
    if (xm > 0.0) {
      wice_new = Kokkos::fmax(0.0, wice - xm);
    } else if (xm < 0.0) {
      if (wmass0 < supercool) {
        wice_new = 0.0;
      } else {
        wice_new = Kokkos::fmin(wmass0 - supercool, wice - xm);
      }
    } else {
      wice_new = wice;
    }

    const Real heatr = hm - hfus * (wice - wice_new) / dtime;
    water_liq(l, k) = Kokkos::fmax(0.0, wmass0 - wice_new);
    water_ice(l, k) = wice_new;

    if (Kokkos::fabs(heatr) > 0.0) {
      temp(l, k) = tfrz + fact_k * heatr;
    }
  }
}

// Apply soil phase change for deep impervious road layers (k=nActive..
// NUM_LAYERS_ABV_BEDROCK-1). Same algorithm as ApplyPerviousRoadSoilPhaseChange
// but operates on separate WaterLiquid/WaterIce views for impervious road.
// The soil hydraulic properties (watsat, bsw, sucsat) are shared with the
// pervious road (same underlying natural soil).
template <typename TempView, typename WaterView, typename CvView,
          typename SoilView, typename DzView, typename IntView>
KOKKOS_INLINE_FUNCTION void ApplyDeepImperviousRoadPhaseChange(
    int l, int numSoilLayers, const IntView &numActiveLayers, TempView &temp,
    WaterView &water_ice, WaterView &water_liq, const CvView &cv_times_dz,
    const SoilView &watsat, const SoilView &bsw, const SoilView &sucsat,
    const DzView &dz, const DzView &zc, const DzView &zi, Real cpar,
    Real cgrnds, Real cgrndl, Real emiss, Real tgnd0, Real dtime) {

  constexpr Real tfrz = SHR_CONST_TKFRZ;
  constexpr Real hfus = SHR_CONST_LATICE;
  constexpr Real grav = GRAVITY;

  const int nActive = numActiveLayers(l);
  for (int k = nActive; k < NUM_LAYERS_ABV_BEDROCK; ++k) {
    const Real T = temp(l, k);
    const Real wice = water_ice(l, k);
    const Real wliq = water_liq(l, k);
    const Real wmass0 = wice + wliq;
    Real fact_k = dtime / cv_times_dz(l, k);
    if (k == 0) {
      const Real dz1 = zc(l, 0) - zi(l, 0);
      const Real dz2 = zc(l, 1) - zi(l, 0);
      const Real dz_eff = 0.5 * (dz1 + cpar * dz2);
      fact_k = fact_k * dz(l, 0) / dz_eff;
    }

    // Supercool threshold
    Real supercool = 0.0;
    if (T < tfrz) {
      const Real smp = hfus * (tfrz - T) / (grav * T) * 1000.0;
      if (smp > 0.0 && sucsat(l, k) > 0.0 && bsw(l, k) > 0.0) {
        const Real sc = watsat(l, k) *
                        Kokkos::pow(smp / sucsat(l, k), -1.0 / bsw(l, k)) *
                        dz(l, k) * 1000.0;
        supercool = Kokkos::fmax(0.0, sc);
      }
    }

    int imelt = 0;
    Real tinc = 0.0;

    if (wice > 0.0 && T > tfrz) {
      imelt = 1;
      tinc = tfrz - T;
      temp(l, k) = tfrz;
    } else if (wliq > supercool && T < tfrz) {
      imelt = 2;
      tinc = tfrz - T;
      temp(l, k) = tfrz;
    }

    if (imelt == 0)
      continue;

    Real hm = -tinc / fact_k;
    if (k == 0) {
      const Real cgrnd = cgrnds + cgrndl * SHR_CONST_LATVAP;
      const Real dhsdT =
          -cgrnd - 4.0 * emiss * STEBOL * Kokkos::pow(tgnd0, 3.0);
      hm += dhsdT * tinc;
    }

    if (imelt == 1 && hm < 0.0) {
      temp(l, k) = tfrz + tinc;
      continue;
    }
    if (imelt == 2 && hm > 0.0) {
      temp(l, k) = tfrz + tinc;
      continue;
    }

    const Real xm = hm * dtime / hfus;

    Real wice_new;
    if (xm > 0.0) {
      wice_new = Kokkos::fmax(0.0, wice - xm);
    } else if (xm < 0.0) {
      if (wmass0 < supercool) {
        wice_new = 0.0;
      } else {
        wice_new = Kokkos::fmin(wmass0 - supercool, wice - xm);
      }
    } else {
      wice_new = wice;
    }

    const Real heatr = hm - hfus * (wice - wice_new) / dtime;
    water_liq(l, k) = Kokkos::fmax(0.0, wmass0 - wice_new);
    water_ice(l, k) = wice_new;

    if (Kokkos::fabs(heatr) > 0.0) {
      temp(l, k) = tfrz + fact_k * heatr;
    }
  }
}

// Compute 1D heat diffusion for all urban surfaces
void ComputeHeatDiffusion(URBANXX::_p_UrbanType &urban) {

  const int numLandunits = urban.numLandunits;
  const int numSoilLayers = urban.numSoilLayers;
  const int numUrbanLayers = urban.numUrbanLayers;

  // Access pervious road soil property views
  auto perv_tk_minerals = urban.perviousRoad.soil.TkMinerals;
  auto perv_tk_dry = urban.perviousRoad.soil.TkDry;
  auto perv_tkLayer = urban.perviousRoad.TkLayer;
  auto perv_tkInterface = urban.perviousRoad.TkInterface;
  auto perv_watsat = urban.perviousRoad.soil.WatSat;
  auto perv_water_liquid = urban.perviousRoad.soil.WaterLiquid;
  auto perv_water_ice = urban.perviousRoad.soil.WaterIce;
  auto perv_dz = urban.perviousRoad.Dz;
  auto perv_zc = urban.perviousRoad.Zc;
  auto perv_zi = urban.perviousRoad.Zi;
  auto perv_temp = urban.perviousRoad.Temperature;
  auto perv_cv_solids = urban.perviousRoad.soil.CvSolids;
  auto perv_cv_times_dz = urban.perviousRoad.CvTimesDz;
  auto perv_netLw = urban.perviousRoad.NetLongRad;
  auto perv_netSw = urban.perviousRoad.NetShortRad;
  auto perv_EflxShGrnd = urban.perviousRoad.EflxShGrnd;
  auto perv_QflxEvapSoil = urban.perviousRoad.QflxEvapSoil;
  auto perv_QflxTranEvap = urban.perviousRoad.QflxTranEvap;
  auto perv_Cgrndl = urban.perviousRoad.Cgrndl;
  auto perv_Cgrnds = urban.perviousRoad.Cgrnds;
  auto perv_Temp = urban.perviousRoad.EffectiveSurfTemp;
  auto perv_TGrnd0 = urban.perviousRoad.TGrnd0;
  auto perv_emiss = urban.urbanParams.emissivity.PerviousRoad;
  auto perv_H2OSno = urban.perviousRoad.H2OSno;
  auto perv_SnowDepth = urban.perviousRoad.SnowDepth;
  auto perv_bsw = urban.perviousRoad.soil.Bsw;
  auto perv_sucsat = urban.perviousRoad.soil.SucSat;

  // Access roof property views
  auto roof_netLw = urban.roof.NetLongRad;
  auto roof_netSw = urban.roof.NetShortRad;
  auto roof_EflxShGrnd = urban.roof.EflxShGrnd;
  auto roof_QflxEvapSoil = urban.roof.QflxEvapSoil;
  auto roof_QflxTranEvap = urban.roof.QflxTranEvap;
  auto roof_Cgrnds = urban.roof.Cgrnds;
  auto roof_Cgrndl = urban.roof.Cgrndl;
  auto roof_Temp = urban.roof.EffectiveSurfTemp;
  auto roof_TGrnd0 = urban.roof.TGrnd0;
  auto roof_emiss = urban.urbanParams.emissivity.Roof;
  auto roof_temp = urban.roof.Temperature;
  auto roof_tkInterface = urban.roof.TkInterface;
  auto roof_tkLayer = urban.roof.TkLayer;
  auto roof_cv_times_dz = urban.roof.CvTimesDz;
  auto roof_dz = urban.roof.Dz;
  auto roof_zc = urban.roof.Zc;
  auto roof_zi = urban.roof.Zi;
  auto roof_H2OSno = urban.roof.H2OSno;
  auto roof_SnowDepth = urban.roof.SnowDepth;
  auto roof_TopH2OSoiLiq = urban.roof.TopH2OSoiLiq;
  auto roof_TopH2OSoiIce = urban.roof.TopH2OSoiIce;

  // Access impervious road property views
  auto imperv_tkLayer = urban.imperviousRoad.TkLayer;
  auto imperv_tkInterface = urban.imperviousRoad.TkInterface;
  auto imperv_dz = urban.imperviousRoad.Dz;
  auto imperv_zc = urban.imperviousRoad.Zc;
  auto imperv_zi = urban.imperviousRoad.Zi;
  auto imperv_temp = urban.imperviousRoad.Temperature;
  auto imperv_cv_times_dz = urban.imperviousRoad.CvTimesDz;
  auto imperv_netLw = urban.imperviousRoad.NetLongRad;
  auto imperv_netSw = urban.imperviousRoad.NetShortRad;
  auto imperv_EflxShGrnd = urban.imperviousRoad.EflxShGrnd;
  auto imperv_QflxEvapSoil = urban.imperviousRoad.QflxEvapSoil;
  auto imperv_QflxTranEvap = urban.imperviousRoad.QflxTranEvap;
  auto imperv_Cgrnds = urban.imperviousRoad.Cgrnds;
  auto imperv_Cgrndl = urban.imperviousRoad.Cgrndl;
  auto imperv_Temp = urban.imperviousRoad.EffectiveSurfTemp;
  auto imperv_TGrnd0 = urban.imperviousRoad.TGrnd0;
  auto imperv_emiss = urban.urbanParams.emissivity.ImperviousRoad;
  // Per-layer water for dynamic CvTimesDz recomputation
  auto imperv_water_liquid = urban.imperviousRoad.WaterLiquid;
  auto imperv_water_ice = urban.imperviousRoad.WaterIce;
  auto imperv_numActiveLayers = urban.imperviousRoad.NumberOfActiveLayers;
  // Reuse pervious road soil properties for deep impervious road layers
  auto imperv_cv_solids = urban.perviousRoad.soil.CvSolids;
  auto imperv_watsat = urban.perviousRoad.soil.WatSat;
  auto imperv_tk_minerals = urban.perviousRoad.soil.TkMinerals;
  auto imperv_tk_dry = urban.perviousRoad.soil.TkDry;
  auto imperv_H2OSno = urban.imperviousRoad.H2OSno;
  auto imperv_SnowDepth = urban.imperviousRoad.SnowDepth;
  auto imperv_TopH2OSoiLiq = urban.imperviousRoad.TopH2OSoiLiq;
  auto imperv_TopH2OSoiIce = urban.imperviousRoad.TopH2OSoiIce;

  // Access sunlit wall property views
  auto sunwall_netLw = urban.sunlitWall.NetLongRad;
  auto sunwall_netSw = urban.sunlitWall.NetShortRad;
  auto sunwall_EflxShGrnd = urban.sunlitWall.EflxShGrnd;
  auto sunwall_Cgrnds = urban.sunlitWall.Cgrnds;
  auto sunwall_Cgrndl = urban.sunlitWall.Cgrndl;
  auto sunwall_Temp = urban.sunlitWall.EffectiveSurfTemp;
  auto sunwall_TGrnd0 = urban.sunlitWall.TGrnd0;
  auto sunwall_temp = urban.sunlitWall.Temperature;
  auto sunwall_tkInterface = urban.sunlitWall.TkInterface;
  auto sunwall_tkLayer = urban.sunlitWall.TkLayer;
  auto sunwall_cv_times_dz = urban.sunlitWall.CvTimesDz;
  auto sunwall_dz = urban.sunlitWall.Dz;
  auto sunwall_zc = urban.sunlitWall.Zc;
  auto sunwall_zi = urban.sunlitWall.Zi;
  auto wall_emiss = urban.urbanParams.emissivity.Wall;

  // Access shaded wall property views
  auto shadewall_netLw = urban.shadedWall.NetLongRad;
  auto shadewall_netSw = urban.shadedWall.NetShortRad;
  auto shadewall_EflxShGrnd = urban.shadedWall.EflxShGrnd;
  auto shadewall_Cgrnds = urban.shadedWall.Cgrnds;
  auto shadewall_Cgrndl = urban.shadedWall.Cgrndl;
  auto shadewall_Temp = urban.shadedWall.EffectiveSurfTemp;
  auto shadewall_TGrnd0 = urban.shadedWall.TGrnd0;
  auto shadewall_temp = urban.shadedWall.Temperature;
  auto shadewall_tkInterface = urban.shadedWall.TkInterface;
  auto shadewall_tkLayer = urban.shadedWall.TkLayer;
  auto shadewall_cv_times_dz = urban.shadedWall.CvTimesDz;
  auto shadewall_dz = urban.shadedWall.Dz;
  auto shadewall_zc = urban.shadedWall.Zc;
  auto shadewall_zi = urban.shadedWall.Zi;

  // Access building property views
  auto building_temp = urban.building.Temperature;

  // Access geometry and building temperature bounds for internal compute
  auto ht_roof = urban.urbanParams.heights.HtRoof;
  auto canyon_hwr = urban.urbanParams.CanyonHwr;
  auto wt_roof = urban.urbanParams.WtRoof;
  auto bld_min_temp = urban.urbanParams.building.MinTemperature;
  auto bld_max_temp = urban.urbanParams.building.MaxTemperature;

  // Persistent AC heat flux storage (per road area, W/m²), computed from the
  // previous timestep's building heat fluxes and applied to road surfaces
  // in the current timestep.  Mirrors ELM's eflx_heat_from_ac_patch timing.
  auto bldg_eflux_ac = urban.building.EFluxForAC;
  auto bldg_eflux_ac_prev = urban.building.EFluxForAC_Prev;

  // Single parallel kernel over all landunits
  Kokkos::parallel_for(
      "ComputeHeatDiffusion", numLandunits, KOKKOS_LAMBDA(int l) {
        // --- Compute internal building temperature (same formula as ELM
        // UrbanFluxesMod) ---
        // Also detect whether AC or heating is on (building temp out of range).
        bool cool_on = false;
        bool heat_on = false;
        {
          const Real T_roof_inner = roof_temp(l, numUrbanLayers - 1);
          const Real T_sunwall_inner = sunwall_temp(l, numUrbanLayers - 1);
          const Real T_shadwall_inner = shadewall_temp(l, numUrbanLayers - 1);
          const Real lngth_roof =
              (ht_roof(l) / canyon_hwr(l)) * wt_roof(l) / (1.0 - wt_roof(l));
          Real t_bld = (ht_roof(l) * (T_shadwall_inner + T_sunwall_inner) +
                        lngth_roof * T_roof_inner) /
                       (2.0 * ht_roof(l) + lngth_roof);
          // Detect AC/heating BEFORE clamping (mirrors ELM SoilTemperatureMod)
          cool_on = (t_bld > bld_max_temp(l));
          heat_on = (t_bld < bld_min_temp(l));
          // Clamp to [min, max] (same as ELM SoilTemperatureMod L306-321)
          if (t_bld > bld_max_temp(l))
            t_bld = bld_max_temp(l);
          if (t_bld < bld_min_temp(l))
            t_bld = bld_min_temp(l);
          building_temp(l) = t_bld;
        }
        // --- end building temperature computation ---

        // Read AC heat from PREVIOUS timestep to add to road boundary flux.
        // This mirrors ELM's eflx_heat_from_ac_patch timing: UrbanFluxes uses
        // eflx_urban_ac from the previous SolveTemperature, then
        // ComputeGroundHeatFluxAndDeriv adds it to the road eflx_gnet.
        const Real ac_heat_prev = bldg_eflux_ac(l);
        // Save previous value for use in ComputeSoilFluxes (EflxSoilGrnd needs
        // the same AC heat value that ELM's SoilFluxes uses — from UrbanFluxes,
        // which reads the previous timestep's eflx_urban_ac).
        bldg_eflux_ac_prev(l) = ac_heat_prev;

        // Step 1: Compute thermal conductivity for pervious road soil layers
        ComputeSoilThermalConductivity(
            l, numSoilLayers, perv_tk_minerals, perv_tk_dry, perv_tkLayer,
            perv_watsat, perv_water_liquid, perv_water_ice, perv_dz, perv_temp);

        // Step 1b: Update thermal conductivity for impervious road deep soil
        // layers.  ELM uses moisture-dependent (Johansen) conductivity every
        // timestep for layers j > nlev_improad.  URBANxx was using dry
        // conductivity set at init.  This also recomputes the affected
        // interface conductivities.
        UpdateImperviousRoadThermalProperties(
            l, numSoilLayers, imperv_numActiveLayers, imperv_tk_minerals,
            imperv_tk_dry, imperv_water_liquid, imperv_water_ice, imperv_watsat,
            imperv_dz, imperv_temp, imperv_zc, imperv_zi, imperv_tkLayer,
            imperv_tkInterface);

        // Step 2: Compute thermal conductivity at layer interfaces
        // (Note: Interface thermal conductivity for impervious road deep layers
        //  is now handled by UpdateImperviousRoadThermalProperties above.
        //  Walls and roof are computed during initialization only.)

        // Pervious road: all soil layers are active
        ComputeInterfaceThermalConductivity(l, numSoilLayers, numSoilLayers,
                                            perv_tkLayer, perv_tkInterface,
                                            perv_zc, perv_zi, 0.0);

        // Step 3: Compute heat capacity times layer thickness for pervious road
        ComputeSoilHeatCapacityTimesDz(
            l, numSoilLayers, perv_cv_solids, perv_watsat, perv_water_liquid,
            perv_water_ice, perv_dz, perv_cv_times_dz);

        // Step 3b: Update heat capacity times layer thickness for impervious
        // road deep soil layers (beyond the active construction material
        // layers).  The active layers keep fixed cv_road_params*dz from init;
        // only the underlying natural soil layers are moisture-dependent and
        // must be refreshed each timestep to match ELM's cv computation.
        UpdateImperviousRoadHeatCapacityTimesDz(
            l, numSoilLayers, imperv_numActiveLayers, imperv_cv_solids,
            imperv_watsat, imperv_water_liquid, imperv_water_ice, imperv_dz,
            imperv_cv_times_dz);

        // Step 3c: Path A thin-snow heat capacity correction.
        // When snow water is present but no discrete snow layers have formed
        // (thin-snow regime), ELM augments the heat capacity of the top
        // soil/urban layer with cpice * h2osno (SoilTemperatureMod.F90 ~1070).
        // URBANxx must apply the same correction to match ELM's tridiagonal
        // system; otherwise the top-layer temperature responds too strongly to
        // the surface flux, producing a ~0.15 K mismatch on snowy timesteps.
        // cpice = SHR_CONST_CPICE [J/(kg·K)]; H2OSno [kg/m²].
        if (roof_H2OSno(l) > 0.0) {
          roof_cv_times_dz(l, 0) += SHR_CONST_CPICE * roof_H2OSno(l);
        }
        if (imperv_H2OSno(l) > 0.0) {
          imperv_cv_times_dz(l, 0) += SHR_CONST_CPICE * imperv_H2OSno(l);
        }
        if (perv_H2OSno(l) > 0.0) {
          perv_cv_times_dz(l, 0) += SHR_CONST_CPICE * perv_H2OSno(l);
        }

        // Step 4: Compute ground net energy flux and its derivative for all
        // surfaces.
        // For road surfaces, use ac_heat_prev (AC heat from the previous
        // timestep).  This mirrors ELM's 1-timestep lag: UrbanFluxes reads
        // eflx_urban_ac set by the *previous* SoilTemperature call and uses it
        // in the current SoilTemperature call's road boundary condition.

        // Pervious road (with evaporation + previous-timestep AC heat)
        const Real perv_EflxGnet =
            ComputeGroundNetEnergyFlux(perv_netSw(l), perv_netLw(l),
                                       perv_EflxShGrnd(l), perv_QflxEvapSoil(l),
                                       perv_QflxTranEvap(l)) +
            ac_heat_prev;
        const Real perv_DEflxGnet_DTemp = ComputeGroundNetEnergyFluxDerivative(
            perv_Cgrnds(l), perv_Cgrndl(l), perv_emiss(l), perv_Temp(l));

        // Roof (with evaporation)
        const Real roof_EflxGnet = ComputeGroundNetEnergyFlux(
            roof_netSw(l), roof_netLw(l), roof_EflxShGrnd(l),
            roof_QflxEvapSoil(l), roof_QflxTranEvap(l));
        const Real roof_DEflxGnet_DTemp = ComputeGroundNetEnergyFluxDerivative(
            roof_Cgrnds(l), roof_Cgrndl(l), roof_emiss(l), roof_Temp(l));

        // Impervious road (with evaporation + previous-timestep AC heat)
        const Real imperv_EflxGnet =
            ComputeGroundNetEnergyFlux(
                imperv_netSw(l), imperv_netLw(l), imperv_EflxShGrnd(l),
                imperv_QflxEvapSoil(l), imperv_QflxTranEvap(l)) +
            ac_heat_prev;
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

        // Save pre-solve surface temperature for use in UrbanComputeSoilFluxes
        // (tinc calculation).
        roof_TGrnd0(l) = roof_temp(l, 0);
        sunwall_TGrnd0(l) = sunwall_temp(l, 0);
        shadewall_TGrnd0(l) = shadewall_temp(l, 0);
        imperv_TGrnd0(l) = imperv_temp(l, 0);
        perv_TGrnd0(l) = perv_temp(l, 0);

        // Save inner-surface OLD temperatures for building heat flux
        // computation (needed for Crank-Nicolson eflx_building_heat after the
        // solve).
        const int inner_k = numUrbanLayers - 1;
        const Real T_roof_inner_old = roof_temp(l, inner_k);
        const Real T_sunwall_inner_old = sunwall_temp(l, inner_k);
        const Real T_shadewall_inner_old = shadewall_temp(l, inner_k);

        // Solve heat diffusion for roof
        SurfaceProperties roof_surf(l, numUrbanLayers, roof_temp, roof_zc,
                                    roof_zi, roof_dz, roof_tkInterface,
                                    roof_tkLayer, roof_cv_times_dz);
        BoundaryConditions roof_bc(roof_EflxGnet, roof_DEflxGnet_DTemp,
                                   useTopAdjustment_building, capr,
                                   hasBottomBC_building, building_temp(l));
        Solve1DHeatDiffusion(dtime, roof_surf, roof_bc);

        // Phase change (thin-snow): mirrors ELM Phasechange_beta for the
        // snl=0, h2osno>0 case on the roof top layer.
        ApplyThinSnowPhaseChange(
            l, roof_temp, roof_H2OSno, roof_SnowDepth, roof_TopH2OSoiIce(l),
            roof_TopH2OSoiLiq(l), roof_cv_times_dz, roof_dz, roof_zc, roof_zi,
            roof_Cgrnds(l), roof_Cgrndl(l), roof_emiss(l), roof_TGrnd0(l),
            /*useTopLayerAdjustment=*/false, capr, dtime);

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

        // --- Compute building heat flux and update AC heat for next timestep
        // --- This mirrors ELM: after SolveTemperature for walls/roof, compute
        // eflx_building_heat = cnfac*fn_old + (1-cnfac)*fn1_new, then
        // eflx_urban_ac = |eflx_building_heat| when cool_on.
        // The result is stored and used in the road EflxGnet at the NEXT step.
        {
          const Real t_bld = building_temp(l);
          const Real wtr = wt_roof(l);
          const Real hwr = canyon_hwr(l);

          // Helper: compute Crank-Nicolson building heat flux at inner surface
          // fn_old = tk*(t_bld - T_old) / (zi_inner - zc_inner)
          // fn_new = tk*(t_bld - T_new) / (zi_inner - zc_inner)
          // eflx_bldg = CNFAC*fn_old + (1-CNFAC)*fn_new
          auto bldg_heat_flux = [&](Real tk_inner, Real zi_inner, Real zc_inner,
                                    Real T_inner_old,
                                    Real T_inner_new) -> Real {
            const Real inv_dz = 1.0 / (zi_inner - zc_inner);
            const Real fn_old = tk_inner * (t_bld - T_inner_old) * inv_dz;
            const Real fn_new = tk_inner * (t_bld - T_inner_new) * inv_dz;
            return CNFAC * fn_old + (1.0 - CNFAC) * fn_new;
          };

          const Real eflx_bldg_roof = bldg_heat_flux(
              roof_tkInterface(l, inner_k), roof_zi(l, numUrbanLayers),
              roof_zc(l, inner_k), T_roof_inner_old, roof_temp(l, inner_k));

          const Real eflx_bldg_sunwall = bldg_heat_flux(
              sunwall_tkInterface(l, inner_k), sunwall_zi(l, numUrbanLayers),
              sunwall_zc(l, inner_k), T_sunwall_inner_old,
              sunwall_temp(l, inner_k));

          const Real eflx_bldg_shadewall = bldg_heat_flux(
              shadewall_tkInterface(l, inner_k),
              shadewall_zi(l, numUrbanLayers), shadewall_zc(l, inner_k),
              T_shadewall_inner_old, shadewall_temp(l, inner_k));

          // Compute new AC heat per road area (mirrors ELM UrbanFluxesMod):
          //   eflx_heat_from_ac(l) = wt_roof*ac_roof +
          //   (1-wt_roof)*hwr*(sun+shade)
          //   eflx_heat_from_ac_patch = eflx_heat_from_ac(l) / (1 - wt_roof)
          Real ac_heat_new = 0.0;
          if (cool_on) {
            const Real ac_roof = Kokkos::fabs(eflx_bldg_roof);
            const Real ac_sunwall = Kokkos::fabs(eflx_bldg_sunwall);
            const Real ac_shadewall = Kokkos::fabs(eflx_bldg_shadewall);
            const Real eflx_heat_from_ac_l =
                wtr * ac_roof + (1.0 - wtr) * hwr * (ac_sunwall + ac_shadewall);
            ac_heat_new = eflx_heat_from_ac_l / (1.0 - wtr);
          }
          // Store the computed AC heat so it is available as ac_heat_prev
          // in the next timestep (mirrors ELM's 1-timestep lag for road BCs).
          bldg_eflux_ac(l) = ac_heat_new;
        }
        // --- end AC heat update ---

        // Solve heat diffusion for impervious road
        SurfaceProperties imperv_surf(l, numSoilLayers, imperv_temp, imperv_zc,
                                      imperv_zi, imperv_dz, imperv_tkInterface,
                                      imperv_tkLayer, imperv_cv_times_dz);
        BoundaryConditions imperv_bc(imperv_EflxGnet, imperv_DEflxGnet_DTemp,
                                     useTopAdjustment_road, capr,
                                     hasBottomBC_road, noBottomTemp);
        Solve1DHeatDiffusion(dtime, imperv_surf, imperv_bc);

        // Phase change (thin-snow): mirrors ELM Phasechange_beta for the
        // snl=0, h2osno>0 case on the impervious road top layer.
        ApplyThinSnowPhaseChange(
            l, imperv_temp, imperv_H2OSno, imperv_SnowDepth,
            imperv_TopH2OSoiIce(l), imperv_TopH2OSoiLiq(l), imperv_cv_times_dz,
            imperv_dz, imperv_zc, imperv_zi, imperv_Cgrnds(l), imperv_Cgrndl(l),
            imperv_emiss(l), imperv_TGrnd0(l), /*useTopLayerAdjustment=*/true,
            capr, dtime);

        // Phase change for deep impervious road layers (beyond active road
        // construction material).  Replaces the former simplified algorithm
        // with the ELM-matching formulation (supercool threshold + residual
        // temperature correction).
        ApplyDeepImperviousRoadPhaseChange(
            l, numSoilLayers, imperv_numActiveLayers, imperv_temp,
            imperv_water_ice, imperv_water_liquid, imperv_cv_times_dz,
            imperv_watsat, perv_bsw, perv_sucsat, imperv_dz, imperv_zc,
            imperv_zi, capr, imperv_Cgrnds(l), imperv_Cgrndl(l),
            imperv_emiss(l), imperv_TGrnd0(l), dtime);

        // Solve heat diffusion for pervious road
        SurfaceProperties perv_surf(l, numSoilLayers, perv_temp, perv_zc,
                                    perv_zi, perv_dz, perv_tkInterface,
                                    perv_tkLayer, perv_cv_times_dz);
        BoundaryConditions perv_bc(perv_EflxGnet, perv_DEflxGnet_DTemp,
                                   useTopAdjustment_road, capr,
                                   hasBottomBC_road, noBottomTemp);
        Solve1DHeatDiffusion(dtime, perv_surf, perv_bc);

        // Phase change (thin-snow): mirrors ELM Phasechange_beta for the
        // snl=0, h2osno>0 case on the pervious road top layer.
        ApplyThinSnowPhaseChange(l, perv_temp, perv_H2OSno, perv_SnowDepth,
                                 perv_water_ice(l, 0), perv_water_liquid(l, 0),
                                 perv_cv_times_dz, perv_dz, perv_zc, perv_zi,
                                 perv_Cgrnds(l), perv_Cgrndl(l), perv_emiss(l),
                                 perv_TGrnd0(l),
                                 /*useTopLayerAdjustment=*/true, capr, dtime);

        // Soil phase change for all pervious road layers (melt/freeze within
        // natural soil): mirrors ELM Phasechange_beta is_soil/icol_road_perv.
        ApplyPerviousRoadSoilPhaseChange(
            l, numSoilLayers, perv_temp, perv_water_ice, perv_water_liquid,
            perv_cv_times_dz, perv_watsat, perv_bsw, perv_sucsat, perv_dz,
            perv_zc, perv_zi, capr, perv_Cgrnds(l), perv_Cgrndl(l),
            perv_emiss(l), perv_TGrnd0(l), dtime);

        // Update EffectiveSurfTemp from the newly computed first-layer
        // temperature for each surface. This keeps the radiation derivative
        // term (4·ε·σ·T³) in sync without requiring an explicit setter call
        // from the host (ELM/Fortran) before the next timestep.
        roof_Temp(l) = roof_temp(l, 0);
        sunwall_Temp(l) = sunwall_temp(l, 0);
        shadewall_Temp(l) = shadewall_temp(l, 0);
        imperv_Temp(l) = imperv_temp(l, 0);
        perv_Temp(l) = perv_temp(l, 0);
      });
  Kokkos::fence();
}

} // namespace URBANXX

// ============================================================================
// Public C API
// ============================================================================

extern "C" {

void UrbanComputeHeatDiffusion(UrbanType urban, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    URBANXX::ComputeHeatDiffusion(*urban);
    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
