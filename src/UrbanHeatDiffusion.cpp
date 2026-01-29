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
    if (l == 0)
      printf("tkInterface(%d,%d) = %18.16f\n", l, k, tkInterface(l, k));
  }

  // Last interface value set to fill value
  tkInterface(l, numLayers - 1) = tkFillValue;
}

// Compute 1D heat diffusion for all urban surfaces
void ComputeHeatDiffusion(URBANXX::_p_UrbanType &urban) {

  const int numLandunits = urban.numLandunits;
  const int numSoilLayers = urban.numSoilLayers;

  // Access soil property views
  auto tk_minerals = urban.perviousRoad.soil.TkMinerals;
  auto tk_dry = urban.perviousRoad.soil.TkDry;
  auto tkLayer = urban.perviousRoad.TkLayer;
  auto tkInterface = urban.perviousRoad.TkInterface;
  auto watsat = urban.perviousRoad.soil.WatSat;
  auto water_liquid = urban.perviousRoad.soil.WaterLiquid;
  auto water_ice = urban.perviousRoad.soil.WaterIce;
  auto dz = urban.perviousRoad.Dz;
  auto zc = urban.perviousRoad.Zc;
  auto zi = urban.perviousRoad.Zi;
  auto temp = urban.perviousRoad.Temperature;

  // Single parallel kernel over all landunits
  Kokkos::parallel_for(
      "ComputeHeatDiffusion", numLandunits, KOKKOS_LAMBDA(int l) {
        // Step 1: Compute thermal conductivity for pervious road soil layers
        ComputeSoilThermalConductivity(l, numSoilLayers, tk_minerals, tk_dry,
                                       tkLayer, watsat, water_liquid, water_ice,
                                       dz, temp);

        // Step 2: Compute thermal conductivity at layer interfaces
        // For pervious road, all soil layers are active (numActiveLayers =
        // numSoilLayers) and inactive interfaces are set to 0.0
        ComputeInterfaceThermalConductivity(l, numSoilLayers, numSoilLayers,
                                            tkLayer, tkInterface, zc, zi, 0.0);

        // TODO: Add remaining heat diffusion steps:
        // Step 3: Setup tridiagonal system for heat conduction equation
        // Step 4: Apply boundary conditions (surface flux, bottom temperature)
        // Step 5: Solve tridiagonal system for new temperatures
        // Step 6: Compute ground heat flux for roof, wall, and road surfaces
      });
  Kokkos::fence();
}

} // namespace URBANXX
