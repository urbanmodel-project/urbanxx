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

    const int numLandunits = urban->numLandunits;

    // Temperature initialization constants
    constexpr Real TEMP_ROOF_INIT = 274.0;
    constexpr Real TEMP_WALL_INIT = 292.0;
    constexpr Real TEMP_ROAD_INIT = 274.0;

    // Initialize all surface temperatures in a single parallel loop
    Kokkos::parallel_for(
        "InitializeSurfaceTemperatures", numLandunits, KOKKOS_LAMBDA(int l) {
          roofTemp(l) = TEMP_ROOF_INIT;
          imperviousRoadTemp(l) = TEMP_ROAD_INIT;
          perviousRoadTemp(l) = TEMP_ROAD_INIT;
          sunlitWallTemp(l) = TEMP_WALL_INIT;
          shadedWallTemp(l) = TEMP_WALL_INIT;
        });
    Kokkos::fence();

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
