#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

extern "C" {

using namespace URBANXX;

void UrbanInitializeTemperature(UrbanType urban, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    // Get references to temperature views
    auto &roofTemp = urban->roof.Temperature;
    auto &imperviousRoadTemp = urban->imperviousRoad.Temperature;
    auto &perviousRoadTemp = urban->perviousRoad.Temperature;
    auto &sunlitWallTemp = urban->sunlitWall.Temperature;
    auto &shadedWallTemp = urban->shadedWall.Temperature;

    // Get references to canyon air properties
    auto &taf = urban->urbanCanyon.Taf;
    auto &qaf = urban->urbanCanyon.Qaf;

    const int numLandunits = urban->numLandunits;

    // Temperature initialization constants
    constexpr Real TEMP_ROOF_INIT = 292.0;
    constexpr Real TEMP_WALL_INIT = 292.0;
    constexpr Real TEMP_ROAD_INIT = 274.0;
    constexpr Real TEMP_CANYON_AIR_INIT = 283.0;
    constexpr Real QAF_INIT = 1.e-4; // kg/kg

    // Initialize surface temperatures and canyon air properties
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    Kokkos::parallel_for(
        "InitializeSurfaceTemperatures",
        Kokkos::RangePolicy<ExecSpace>(0, numLandunits),
        KOKKOS_LAMBDA(int l) {
          // Initialize temperatures
          roofTemp(l) = TEMP_ROOF_INIT;
          imperviousRoadTemp(l) = TEMP_ROAD_INIT;
          perviousRoadTemp(l) = TEMP_ROAD_INIT;
          sunlitWallTemp(l) = TEMP_WALL_INIT;
          shadedWallTemp(l) = TEMP_WALL_INIT;

          // Initialize canyon air properties
          taf(l) = TEMP_CANYON_AIR_INIT; // Initialize to canyon air temperature
          qaf(l) = QAF_INIT; // Initialize to reasonable humidity value (kg/kg)
        });
    Kokkos::fence();

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
