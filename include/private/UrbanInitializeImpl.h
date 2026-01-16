#ifndef URBAN_INITIALIZE_IMPL_H
#define URBAN_INITIALIZE_IMPL_H

#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// POD struct to hold temperature views - trivially copyable for CUDA lambdas
struct TemperatureViews {
  Array1DR8 roofTemp;
  Array1DR8 imperviousRoadTemp;
  Array1DR8 perviousRoadTemp;
  Array1DR8 sunlitWallTemp;
  Array1DR8 shadedWallTemp;
  Array1DR8 taf;
  Array1DR8 qaf;
};

// Class to handle temperature initialization with KOKKOS_CLASS_LAMBDA
class UrbanTemperatureInitializer {
private:
  // POD struct containing all views
  TemperatureViews views;
  int numLandunits;

  // Temperature initialization constants
  static constexpr Real TEMP_ROOF_INIT = 292.0;
  static constexpr Real TEMP_WALL_INIT = 292.0;
  static constexpr Real TEMP_ROAD_INIT = 274.0;
  static constexpr Real TEMP_CANYON_AIR_INIT = 283.0;
  static constexpr Real QAF_INIT = 1.e-4; // kg/kg

public:
  UrbanTemperatureInitializer(_p_UrbanType *urban)
      : views{urban->roof.Temperature,
              urban->imperviousRoad.Temperature,
              urban->perviousRoad.Temperature,
              urban->sunlitWall.Temperature,
              urban->shadedWall.Temperature,
              urban->urbanCanyon.Taf,
              urban->urbanCanyon.Qaf},
        numLandunits(urban->numLandunits) {}

  // Method to initialize temperatures - called from device code
  KOKKOS_INLINE_FUNCTION
  void initializeAtLandunit(const int l) const {
    // Initialize temperatures
    views.roofTemp(l) = TEMP_ROOF_INIT;
    views.imperviousRoadTemp(l) = TEMP_ROAD_INIT;
    views.perviousRoadTemp(l) = TEMP_ROAD_INIT;
    views.sunlitWallTemp(l) = TEMP_WALL_INIT;
    views.shadedWallTemp(l) = TEMP_WALL_INIT;

    // Initialize canyon air properties
    views.taf(l) = TEMP_CANYON_AIR_INIT;
    views.qaf(l) = QAF_INIT;
  }

  // Main method to run initialization
  void run() {
    // Use KOKKOS_CLASS_LAMBDA to capture *this (which contains POD struct)
    Kokkos::parallel_for(
        "InitializeSurfaceTemperatures", numLandunits,
        KOKKOS_CLASS_LAMBDA(const int l) { initializeAtLandunit(l); });
    Kokkos::fence();
  }
};

} // namespace URBANXX

#endif // URBAN_INITIALIZE_IMPL_H

