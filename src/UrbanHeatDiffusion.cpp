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
  }

  // Last interface value set to fill value
  tkInterface(l, numLayers - 1) = tkFillValue;
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

  // Access impervious road views
  auto imperv_tkLayer = urban.imperviousRoad.TkLayer;
  auto imperv_tkInterface = urban.imperviousRoad.TkInterface;
  auto imperv_zc = urban.imperviousRoad.Zc;
  auto imperv_zi = urban.imperviousRoad.Zi;
  auto imperv_numActiveLayers = urban.imperviousRoad.NumberOfActiveLayers;

  // Access sunlit wall views
  auto sunlit_tkLayer = urban.sunlitWall.TkLayer;
  auto sunlit_tkInterface = urban.sunlitWall.TkInterface;
  auto sunlit_zc = urban.sunlitWall.Zc;
  auto sunlit_zi = urban.sunlitWall.Zi;

  // Access shaded wall views
  auto shaded_tkLayer = urban.shadedWall.TkLayer;
  auto shaded_tkInterface = urban.shadedWall.TkInterface;
  auto shaded_zc = urban.shadedWall.Zc;
  auto shaded_zi = urban.shadedWall.Zi;

  // Access roof views
  auto roof_tkLayer = urban.roof.TkLayer;
  auto roof_tkInterface = urban.roof.TkInterface;
  auto roof_zc = urban.roof.Zc;
  auto roof_zi = urban.roof.Zi;

  // Single parallel kernel over all landunits
  Kokkos::parallel_for(
      "ComputeHeatDiffusion", numLandunits, KOKKOS_LAMBDA(int l) {
        // Step 1: Compute thermal conductivity for pervious road soil layers
        ComputeSoilThermalConductivity(
            l, numSoilLayers, perv_tk_minerals, perv_tk_dry, perv_tkLayer,
            perv_watsat, perv_water_liquid, perv_water_ice, perv_dz, perv_temp);

        // Step 2: Compute thermal conductivity at layer interfaces

        // Pervious road: all soil layers are active
        ComputeInterfaceThermalConductivity(l, numSoilLayers, numSoilLayers,
                                            perv_tkLayer, perv_tkInterface,
                                            perv_zc, perv_zi, 0.0);

        // Impervious road: variable active layers based on density class
        ComputeInterfaceThermalConductivity(
            l, numUrbanLayers, imperv_numActiveLayers(l), imperv_tkLayer,
            imperv_tkInterface, imperv_zc, imperv_zi,
            imperv_tkLayer(l, numUrbanLayers - 1));

        // Sunlit wall: all urban layers are active
        ComputeInterfaceThermalConductivity(
            l, numUrbanLayers, numUrbanLayers, sunlit_tkLayer,
            sunlit_tkInterface, sunlit_zc, sunlit_zi,
            sunlit_tkLayer(l, numUrbanLayers - 1));

        // Shaded wall: all urban layers are active
        ComputeInterfaceThermalConductivity(
            l, numUrbanLayers, numUrbanLayers, shaded_tkLayer,
            shaded_tkInterface, shaded_zc, shaded_zi,
            shaded_tkLayer(l, numUrbanLayers - 1));

        // Roof: all urban layers are active
        ComputeInterfaceThermalConductivity(
            l, numUrbanLayers, numUrbanLayers, roof_tkLayer, roof_tkInterface,
            roof_zc, roof_zi, roof_tkLayer(l, numUrbanLayers - 1));

        // TODO: Add remaining heat diffusion steps:
        // Step 3: Setup tridiagonal system for heat conduction equation
        // Step 4: Apply boundary conditions (surface flux, bottom temperature)
        // Step 5: Solve tridiagonal system for new temperatures
        // Step 6: Compute ground heat flux for roof, wall, and road surfaces
      });
  Kokkos::fence();
}

} // namespace URBANXX
