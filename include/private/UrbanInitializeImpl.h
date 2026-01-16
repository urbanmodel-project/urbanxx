#ifndef URBAN_INITIALIZE_IMPL_H
#define URBAN_INITIALIZE_IMPL_H

#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// Class to handle temperature initialization with KOKKOS_CLASS_LAMBDA
class UrbanTemperatureInitializer {
private:
  // References to temperature views
  Array1DR8 roofTemp;
  Array1DR8 imperviousRoadTemp;
  Array1DR8 perviousRoadTemp;
  Array1DR8 sunlitWallTemp;
  Array1DR8 shadedWallTemp;

  // References to canyon air properties
  Array1DR8 taf;
  Array1DR8 qaf;

  int numLandunits;

  // Temperature initialization constants
  static constexpr Real TEMP_ROOF_INIT = 292.0;
  static constexpr Real TEMP_WALL_INIT = 292.0;
  static constexpr Real TEMP_ROAD_INIT = 274.0;
  static constexpr Real TEMP_CANYON_AIR_INIT = 283.0;
  static constexpr Real QAF_INIT = 1.e-4; // kg/kg

public:
  UrbanTemperatureInitializer(_p_UrbanType *urban)
      : roofTemp(urban->roof.Temperature),
        imperviousRoadTemp(urban->imperviousRoad.Temperature),
        perviousRoadTemp(urban->perviousRoad.Temperature),
        sunlitWallTemp(urban->sunlitWall.Temperature),
        shadedWallTemp(urban->shadedWall.Temperature),
        taf(urban->urbanCanyon.Taf), qaf(urban->urbanCanyon.Qaf),
        numLandunits(urban->numLandunits) {}

  // Method to initialize temperatures - called from device code
  KOKKOS_INLINE_FUNCTION
  void initializeAtLandunit(const int l) const {
    // Initialize temperatures
    roofTemp(l) = TEMP_ROOF_INIT;
    imperviousRoadTemp(l) = TEMP_ROAD_INIT;
    perviousRoadTemp(l) = TEMP_ROAD_INIT;
    sunlitWallTemp(l) = TEMP_WALL_INIT;
    shadedWallTemp(l) = TEMP_WALL_INIT;

    // Initialize canyon air properties
    taf(l) = TEMP_CANYON_AIR_INIT;
    qaf(l) = QAF_INIT;
  }

  // Main method to run initialization
  void run() {
    Kokkos::parallel_for(
        "InitializeSurfaceTemperatures", numLandunits,
        KOKKOS_CLASS_LAMBDA(const int l) { initializeAtLandunit(l); });
    Kokkos::fence();
  }
};

} // namespace URBANXX

#endif // URBAN_INITIALIZE_IMPL_H
