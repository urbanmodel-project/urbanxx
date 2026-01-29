// Urban Thermal Property Functions
// Shared thermal conductivity calculations for urban surfaces
// Based on ELM SoilTemperatureMod.F90

#ifndef URBAN_THERMAL_FUNCTIONS_H
#define URBAN_THERMAL_FUNCTIONS_H

#include "private/UrbanConstants.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

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

// Compute thermal conductivity at layer interfaces using harmonic averaging
// Based on ELM SoilTemperatureMod.F90
template <typename ViewType>
KOKKOS_INLINE_FUNCTION void ComputeInterfaceThermalConductivity(
    const int l, const int numLayers, const int numActiveLayers,
    const ViewType &tkLayer, const ViewType &tkInterface, const ViewType &zc,
    const ViewType &zi, const Real tkFillValue) {

  for (int k = 0; k < numLayers - 1; ++k) {
    if (k < numActiveLayers - 1) {
      // Harmonic average of layer thermal conductivities
      // tk(j) = thk(j)*thk(j+1)*(z(j+1)-z(j)) /
      //         (thk(j)*(z(j+1)-zi(j))+thk(j+1)*(zi(j)-z(j)))
      tkInterface(l, k) = tkLayer(l, k) * tkLayer(l, k + 1) *
                          (zc(l, k + 1) - zc(l, k)) /
                          (tkLayer(l, k) * (zc(l, k + 1) - zi(l, k + 1)) +
                           tkLayer(l, k + 1) * (zi(l, k + 1) - zc(l, k)));
    } else {
      // Inactive layers set to fill value
      tkInterface(l, k) = tkFillValue;
    }
  }

  // Last interface value set to fill value
  tkInterface(l, numLayers - 1) = tkFillValue;
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
    cvTimesDz(l, k) = cv_solid_component + cv_water_component;
  }
}

} // namespace URBANXX

#endif // URBAN_THERMAL_FUNCTIONS_H
