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
    
    // Create Kokkos::View objects by copying, not accessing through pointer
    printf("DEBUG: Creating local View copies...\n");
    fflush(stdout);
    
    // Use decltype to get exact type
    using RoofTempView = decltype(urban->roof.Temperature);
    using RoadTempView = decltype(urban->imperviousRoad.Temperature);
    using WallTempView = decltype(urban->sunlitWall.Temperature);
    using CanyonView = decltype(urban->urbanCanyon.Taf);
    
    RoofTempView roofTemp(urban->roof.Temperature);
    RoadTempView imperviousRoadTemp(urban->imperviousRoad.Temperature);
    RoadTempView perviousRoadTemp(urban->perviousRoad.Temperature);
    WallTempView sunlitWallTemp(urban->sunlitWall.Temperature);
    WallTempView shadedWallTemp(urban->shadedWall.Temperature);
    CanyonView taf(urban->urbanCanyon.Taf);
    CanyonView qaf(urban->urbanCanyon.Qaf);
    
    printf("DEBUG: Local View copies created\n");
    printf("DEBUG: View labels and sizes:\n");
    printf("  roofTemp: '%s' size=%zu data=%p\n", roofTemp.label().c_str(), roofTemp.size(), roofTemp.data());
    printf("  imperviousRoadTemp: '%s' size=%zu data=%p\n", imperviousRoadTemp.label().c_str(), imperviousRoadTemp.size(), imperviousRoadTemp.data());
    printf("  perviousRoadTemp: '%s' size=%zu data=%p\n", perviousRoadTemp.label().c_str(), perviousRoadTemp.size(), perviousRoadTemp.data());
    printf("  sunlitWallTemp: '%s' size=%zu data=%p\n", sunlitWallTemp.label().c_str(), sunlitWallTemp.size(), sunlitWallTemp.data());
    printf("  shadedWallTemp: '%s' size=%zu data=%p\n", shadedWallTemp.label().c_str(), shadedWallTemp.size(), shadedWallTemp.data());
    printf("  taf: '%s' size=%zu data=%p\n", taf.label().c_str(), taf.size(), taf.data());
    printf("  qaf: '%s' size=%zu data=%p\n", qaf.label().c_str(), qaf.size(), qaf.data());
    fflush(stdout);

    // Temperature initialization constants
    constexpr Real TEMP_ROOF_INIT = 292.0;
    constexpr Real TEMP_WALL_INIT = 292.0;
    constexpr Real TEMP_ROAD_INIT = 274.0;
    constexpr Real TEMP_CANYON_AIR_INIT = 283.0;
    constexpr Real QAF_INIT = 1.e-4; // kg/kg

    printf("DEBUG: About to launch parallel_for...\n");
    fflush(stdout);
    
    // Get raw pointers from Views - these should work in CUDA lambdas
    Real* roofTempPtr = roofTemp.data();
    Real* imperviousRoadTempPtr = imperviousRoadTemp.data();
    Real* perviousRoadTempPtr = perviousRoadTemp.data();
    Real* sunlitWallTempPtr = sunlitWallTemp.data();
    Real* shadedWallTempPtr = shadedWallTemp.data();
    Real* tafPtr = taf.data();
    Real* qafPtr = qaf.data();
    
    printf("DEBUG: Got raw pointers\n");
    fflush(stdout);
    
    // Initialize surface temperatures and canyon air properties
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    
    // Test 1: Try just roofTemp first using raw pointer
    printf("DEBUG: Testing roofTemp with raw pointer...\n");
    fflush(stdout);
    Kokkos::parallel_for(
        "InitRoofTemp",
        Kokkos::RangePolicy<ExecSpace>(0, numLandunits),
        KOKKOS_LAMBDA(int l) {
          roofTempPtr[l] = TEMP_ROOF_INIT;
        });
    Kokkos::fence();
    printf("DEBUG: roofTemp initialization completed\n");
    fflush(stdout);
    
    // Test 2: Try imperviousRoadTemp
    printf("DEBUG: Testing imperviousRoadTemp with raw pointer...\n");
    fflush(stdout);
    Kokkos::parallel_for(
        "InitImperviousRoadTemp",
        Kokkos::RangePolicy<ExecSpace>(0, numLandunits),
        KOKKOS_LAMBDA(int l) {
          imperviousRoadTempPtr[l] = TEMP_ROAD_INIT;
        });
    Kokkos::fence();
    printf("DEBUG: imperviousRoadTemp initialization completed\n");
    fflush(stdout);
    
    // Continue with the rest...
    Kokkos::parallel_for(
        "InitializeSurfaceTemperatures",
        Kokkos::RangePolicy<ExecSpace>(0, numLandunits),
        KOKKOS_LAMBDA(int l) {
          roofTempPtr[l] = TEMP_ROOF_INIT;
          imperviousRoadTempPtr[l] = TEMP_ROAD_INIT;
          perviousRoadTempPtr[l] = TEMP_ROAD_INIT;
          sunlitWallTempPtr[l] = TEMP_WALL_INIT;
          shadedWallTempPtr[l] = TEMP_WALL_INIT;
          sunlitWallTempPtr[l] = TEMP_WALL_INIT;
          shadedWallTempPtr[l] = TEMP_WALL_INIT;

          // Initialize canyon air properties
          tafPtr[l] = TEMP_CANYON_AIR_INIT; // Initialize to canyon air temperature
          qafPtr[l] = QAF_INIT; // Initialize to reasonable humidity value (kg/kg)
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
