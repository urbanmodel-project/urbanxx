// Urban Heat Diffusion Implementation
// 1D heat diffusion dynamics for Urban Surfaces
// Based on ELM SoilTemperature.F90

#include "private/UrbanHeatDiffusionImpl.h"
#include "private/UrbanThermalFunctions.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

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

        // TODO: Add remaining heat diffusion steps:
        // Step 4: Setup tridiagonal system for heat conduction equation
        // Step 5: Apply boundary conditions (surface flux, bottom temperature)
        // Step 6: Solve tridiagonal system for new temperatures
        // Step 7: Compute ground heat flux for roof, wall, and road surfaces
      });
  Kokkos::fence();
}

} // namespace URBANXX
