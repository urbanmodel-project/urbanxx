#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

using namespace URBANXX;

// Internal function to initialize temperatures
// Note: Parameters are pre-validated by caller (UrbanSetup)
// Exception handling is done by the caller
static void UrbanInitializeTemperature(UrbanType urban) {
  // Get references to temperature views
  auto &roofTemp = urban->roof.EffectiveSurfTemp;
  auto &imperviousRoadTemp = urban->imperviousRoad.EffectiveSurfTemp;
  auto &perviousRoadTemp = urban->perviousRoad.EffectiveSurfTemp;
  auto &sunlitWallTemp = urban->sunlitWall.EffectiveSurfTemp;
  auto &shadedWallTemp = urban->shadedWall.EffectiveSurfTemp;

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
  Kokkos::parallel_for(
      "InitializeSurfaceTemperatures", numLandunits, KOKKOS_LAMBDA(int l) {
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
}

extern "C" {

// Public setup function that performs all initialization steps
void UrbanSetup(UrbanType urban, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    // Initialize surface temperatures
    UrbanInitializeTemperature(urban);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
