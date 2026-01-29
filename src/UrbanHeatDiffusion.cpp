// Urban Heat Diffusion Implementation
// 1D heat diffusion dynamics for Urban Surfaces
// Based on ELM SoilTemperature.F90

#include "private/UrbanConstants.h"
#include "private/UrbanHeatDiffusionImpl.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// Compute effective thermal conductivity for pervious road soil layers
// Based on ELM SoilTemperatureMod.F90 (Johansen 1975 model)
template <typename ViewType>
KOKKOS_INLINE_FUNCTION void ComputeSoilThermalConductivity(
    const int l, const int numSoilLayers, const ViewType &tk_minerals,
    const ViewType &tk_dry, const ViewType &tk_layer, const ViewType &watsat,
    const ViewType &water_liquid, const ViewType &water_ice, const ViewType &dz,
    const ViewType &temp) {

  constexpr Real tfrz = SHR_CONST_TKFRZ;
  constexpr Real denh2o = SHR_CONST_RHOWATER;
  constexpr Real denice = SHR_CONST_RHOICE;
  constexpr Real tkwat = TKWATER;
  constexpr Real tkice = TKICE;

  for (int k = 0; k < numSoilLayers; ++k) {
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
      tk_layer(l, k) = dke * dksat + (1.0 - dke) * tk_dry(l, k);
    } else {
      // Dry soil
      tk_layer(l, k) = tk_dry(l, k);
    }
    if (k > NUM_LAYERS_ABV_BEDROCK - 1) {
      tk_layer(l, k) = TK_BEDROCK;
    }
  }
}

// Compute 1D heat diffusion for all urban surfaces
void ComputeHeatDiffusion(URBANXX::_p_UrbanType &urban) {

  const int numLandunits = urban.numLandunits;
  const int numSoilLayers = urban.numSoilLayers;

  // Access soil property views
  auto tk_minerals = urban.perviousRoad.soil.TkMinerals;
  auto tk_dry = urban.perviousRoad.soil.TkDry;
  auto tk_layer = urban.perviousRoad.soil.TkLayer;
  auto watsat = urban.perviousRoad.soil.WatSat;
  auto water_liquid = urban.perviousRoad.soil.WaterLiquid;
  auto water_ice = urban.perviousRoad.soil.WaterIce;
  auto dz = urban.perviousRoad.Dz;
  auto temp = urban.perviousRoad.Temperature;

  // Single parallel kernel over all landunits
  Kokkos::parallel_for(
      "ComputeHeatDiffusion", numLandunits, KOKKOS_LAMBDA(int l) {
        // Step 1: Compute thermal conductivity for pervious road soil layers
        ComputeSoilThermalConductivity(l, numSoilLayers, tk_minerals, tk_dry,
                                       tk_layer, watsat, water_liquid,
                                       water_ice, dz, temp);

        // TODO: Add remaining heat diffusion steps:
        // Step 2: Setup tridiagonal system for heat conduction equation
        // Step 3: Apply boundary conditions (surface flux, bottom temperature)
        // Step 4: Solve tridiagonal system for new temperatures
        // Step 5: Compute ground heat flux for roof, wall, and road surfaces
      });
  Kokkos::fence();
}

} // namespace URBANXX
