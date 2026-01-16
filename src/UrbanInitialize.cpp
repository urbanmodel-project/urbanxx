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
    printf("DEBUG: Starting UrbanInitializeTemperature\n");
    fflush(stdout);
    
    const int numLandunits = urban->numLandunits;
    printf("DEBUG: numLandunits = %d\n", numLandunits);
    fflush(stdout);
    
    // Copy views for device access - capture by value
    printf("DEBUG: Copying roofTemp...\n");
    fflush(stdout);
    Array1DR8 roofTemp = urban->roof.Temperature;
    
    printf("DEBUG: Copying imperviousRoadTemp...\n");
    fflush(stdout);
    Array1DR8 imperviousRoadTemp = urban->imperviousRoad.Temperature;
    
    printf("DEBUG: Copying perviousRoadTemp...\n");
    fflush(stdout);
    Array1DR8 perviousRoadTemp = urban->perviousRoad.Temperature;
    
    printf("DEBUG: Copying sunlitWallTemp...\n");
    fflush(stdout);
    Array1DR8 sunlitWallTemp = urban->sunlitWall.Temperature;
    
    printf("DEBUG: Copying shadedWallTemp...\n");
    fflush(stdout);
    Array1DR8 shadedWallTemp = urban->shadedWall.Temperature;
    
    printf("DEBUG: Copying taf...\n");
    fflush(stdout);
    Array1DR8 taf = urban->urbanCanyon.Taf;
    
    printf("DEBUG: Copying qaf...\n");
    fflush(stdout);
    Array1DR8 qaf = urban->urbanCanyon.Qaf;
    
    printf("DEBUG: All views copied successfully\n");
    fflush(stdout);

    // Temperature initialization constants
    constexpr Real TEMP_ROOF_INIT = 292.0;
    constexpr Real TEMP_WALL_INIT = 292.0;
    constexpr Real TEMP_ROAD_INIT = 274.0;
    constexpr Real TEMP_CANYON_AIR_INIT = 283.0;
    constexpr Real QAF_INIT = 1.e-4; // kg/kg

    printf("DEBUG: About to launch parallel_for...\n");
    fflush(stdout);
    
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
    
    printf("DEBUG: parallel_for launched, calling fence...\n");
    fflush(stdout);
    Kokkos::fence();
    
    printf("DEBUG: fence completed successfully\n");
    fflush(stdout);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
