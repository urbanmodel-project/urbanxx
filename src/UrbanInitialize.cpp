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

// Helper function to compute vertical discretization for a surface with uniform
// layers
template <typename ViewType>
KOKKOS_INLINE_FUNCTION void ComputeVertDiscretizationForRoofOrWall(
    const Real thickness, const int numLevels, const int l, const ViewType &zc,
    const ViewType &dz, const ViewType &zi) {
  // Compute cell centers
  for (int k = 0; k < numLevels; ++k) {
    zc(l, k) = (k + 0.5) * (thickness / numLevels);
  }

  // Compute layer thickness
  dz(l, 0) = 0.5 * (zc(l, 0) + zc(l, 1));
  for (int k = 1; k < numLevels - 1; ++k) {
    dz(l, k) = 0.5 * (zc(l, k + 1) - zc(l, k - 1));
  }
  dz(l, numLevels - 1) = zc(l, numLevels - 1) - zc(l, numLevels - 2);

  // Compute interface depths
  zi(l, 0) = 0.0;
  for (int k = 1; k < numLevels; ++k) {
    zi(l, k) = 0.5 * (zc(l, k - 1) + zc(l, k));
  }
  zi(l, numLevels) = zc(l, numLevels - 1) + 0.5 * dz(l, numLevels - 1);
}

// Helper function to compute vertical discretization for roads with
// exponential layers
template <typename ViewType>
KOKKOS_INLINE_FUNCTION void ComputeVertDiscretizationForRoad(
    const int numLevels, const int l, const ViewType &zc, const ViewType &dz,
    const ViewType &zi, const Real scalez, const Real zecoeff) {
  // Compute cell centers (node depths) with exponential spacing
  for (int k = 0; k < numLevels; ++k) {
    zc(l, k) = scalez * (Kokkos::exp(zecoeff * ((double)k + 0.5)) - 1.0);
  }

  // Compute layer thickness
  dz(l, 0) = 0.5 * (zc(l, 0) + zc(l, 1));
  for (int k = 1; k < numLevels - 1; ++k) {
    dz(l, k) = 0.5 * (zc(l, k + 1) - zc(l, k - 1));
  }
  dz(l, numLevels - 1) = zc(l, numLevels - 1) - zc(l, numLevels - 2);

  // Compute interface depths
  zi(l, 0) = 0.0;
  for (int k = 1; k < numLevels; ++k) {
    zi(l, k) = 0.5 * (zc(l, k - 1) + zc(l, k));
  }
  zi(l, numLevels) = zc(l, numLevels - 1) + 0.5 * dz(l, numLevels - 1);
}

static void UrbanInitializeVerticalDiscretization(UrbanType urban) {
  const int numLandunits = urban->numLandunits;
  const int numLevels = urban->numLevels;

  // Access vertical discretization views
  auto &thick_wall = urban->urbanParams.building.WallThickness;
  auto &thick_roof = urban->urbanParams.building.RoofThickness;

  auto &zc_sunlit_wall = urban->sunlitWall.Zc;
  auto &dz_sunlit_wall = urban->sunlitWall.Dz;
  auto &zi_sunlit_wall = urban->sunlitWall.Zi;
  auto &depth_sunlit_wall = urban->sunlitWall.TotalDepth;

  auto &zc_shaded_wall = urban->shadedWall.Zc;
  auto &dz_shaded_wall = urban->shadedWall.Dz;
  auto &zi_shaded_wall = urban->shadedWall.Zi;
  auto &depth_shaded_wall = urban->shadedWall.TotalDepth;

  auto &zc_pervious_road = urban->perviousRoad.Zc;
  auto &dz_pervious_road = urban->perviousRoad.Dz;
  auto &zi_pervious_road = urban->perviousRoad.Zi;
  auto &depth_pervious_road = urban->perviousRoad.TotalDepth;

  auto &zc_impervious_road = urban->imperviousRoad.Zc;
  auto &dz_impervious_road = urban->imperviousRoad.Dz;
  auto &zi_impervious_road = urban->imperviousRoad.Zi;
  auto &depth_impervious_road = urban->imperviousRoad.TotalDepth;

  auto &zc_roof = urban->roof.Zc;
  auto &dz_roof = urban->roof.Dz;
  auto &zi_roof = urban->roof.Zi;
  auto &depth_roof = urban->roof.TotalDepth;

  // Initialize vertical discretization
  Kokkos::parallel_for(
      "UrbanInitializeVerticalDiscretization", numLandunits,
      KOKKOS_LAMBDA(int l) {
        // Road discretization parameters
        constexpr Real scalez = 0.025;
        constexpr Real zecoeff = 0.5;

        // Sunlit wall - uniform discretization
        ComputeVertDiscretizationForRoofOrWall(thick_wall(l), numLevels, l,
                                               zc_sunlit_wall, dz_sunlit_wall,
                                               zi_sunlit_wall);
        depth_sunlit_wall(l) = thick_wall(l);

        // Shaded wall - uniform discretization
        ComputeVertDiscretizationForRoofOrWall(thick_wall(l), numLevels, l,
                                               zc_shaded_wall, dz_shaded_wall,
                                               zi_shaded_wall);
        depth_shaded_wall(l) = thick_wall(l);

        // Roof - uniform discretization
        ComputeVertDiscretizationForRoofOrWall(thick_roof(l), numLevels, l,
                                               zc_roof, dz_roof, zi_roof);
        depth_roof(l) = thick_roof(l);

        // Pervious road - exponential discretization
        ComputeVertDiscretizationForRoad(numLevels, l, zc_pervious_road,
                                         dz_pervious_road, zi_pervious_road,
                                         scalez, zecoeff);
        depth_pervious_road(l) = zc_pervious_road(l, numLevels);

        // Impervious road - exponential discretization
        ComputeVertDiscretizationForRoad(numLevels, l, zc_impervious_road,
                                         dz_impervious_road, zi_impervious_road,
                                         scalez, zecoeff);
        depth_impervious_road(l) = zc_impervious_road(l, numLevels);
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
    UrbanInitializeVerticalDiscretization(urban);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
