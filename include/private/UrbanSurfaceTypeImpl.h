#ifndef URBAN_SURFACE_TYPE_IMPL_H
#define URBAN_SURFACE_TYPE_IMPL_H
#include "private/DataTypesImpl.h"
#include "private/KokkosViewMacros.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {
// Composite data for combined pervious and impervious road surfaces
struct CompositeRoadSurfaceData {
  DECLARE_DEVICE_VIEW(3DR8, DownwellingShortRad) // downwelling shortwave
                                                 // radiation per unit road area
  DECLARE_DEVICE_VIEW(3DR8, SnowAlbedo)          // snow albedo
  DECLARE_DEVICE_VIEW(
      3DR8, AlbedoWithSnowEffects) // albedo of road including snow effects

  CompositeRoadSurfaceData(int numLandunits, int numRadBands, int numRadTypes) {
    ALLOCATE_DEVICE_VIEW(DownwellingShortRad, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
    ALLOCATE_DEVICE_VIEW(SnowAlbedo, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
    ALLOCATE_DEVICE_VIEW(AlbedoWithSnowEffects, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
  }
};

// Base struct for all surface types with common radiative properties
struct SurfaceDataBase {
  DECLARE_DEVICE_VIEW(3DR8,
                      ReflectedShortRad) // reflected solar radiation per unit
                                         // surface area per unit incident flux
  DECLARE_DEVICE_VIEW(3DR8,
                      AbsorbedShortRad)  // absorbed solar radiation per unit
                                         // surface area per unit incident flux
  DECLARE_DEVICE_VIEW(1DR8, Emissivity)  // emissivity
  DECLARE_DEVICE_VIEW(1DR8, Temperature) // temperature
  DECLARE_DEVICE_VIEW(1DR8, NetLongRad)  // net longwave radiation
  DECLARE_DEVICE_VIEW(1DR8, UpwardLongRad) // upward longwave radiation

  SurfaceDataBase(int numLandunits, int numRadBands, int numRadTypes) {
    ALLOCATE_DEVICE_VIEW(ReflectedShortRad, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
    ALLOCATE_DEVICE_VIEW(AbsorbedShortRad, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
    ALLOCATE_DEVICE_VIEW(Emissivity, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Temperature, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(NetLongRad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(UpwardLongRad, Array1DR8, numLandunits)
  }
};

// For surfaces with snow coverage (roads and roofs)
struct SnowCoveredSurfaceData : SurfaceDataBase {
  DECLARE_DEVICE_VIEW(3DR8, SnowAlbedo) // snow albedo
  DECLARE_DEVICE_VIEW(3DR8,
                      AlbedoWithSnowEffects) // albedo including snow effects

  SnowCoveredSurfaceData(int numLandunits, int numRadBands, int numRadTypes)
      : SurfaceDataBase(numLandunits, numRadBands, numRadTypes) {
    ALLOCATE_DEVICE_VIEW(SnowAlbedo, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
    ALLOCATE_DEVICE_VIEW(AlbedoWithSnowEffects, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
  }
};

struct RoadDataType : SnowCoveredSurfaceData {
  // Inherits all fields from SnowCoveredSurfaceData and SurfaceDataBase
  RoadDataType(int numLandunits, int numRadBands, int numRadTypes)
      : SnowCoveredSurfaceData(numLandunits, numRadBands, numRadTypes) {}
};

struct WallDataType : SurfaceDataBase {
  DECLARE_DEVICE_VIEW(
      3DR8,
      DownwellingShortRad) // downwelling shortwave radiation per unit wall area
  // Inherits common radiative fields from SurfaceDataBase

  WallDataType(int numLandunits, int numRadBands, int numRadTypes)
      : SurfaceDataBase(numLandunits, numRadBands, numRadTypes) {
    ALLOCATE_DEVICE_VIEW(DownwellingShortRad, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
  }
};

struct RoofDataType : SnowCoveredSurfaceData {
  // Inherits all fields from SnowCoveredSurfaceData and SurfaceDataBase

  RoofDataType(int numLandunits, int numRadBands, int numRadTypes)
      : SnowCoveredSurfaceData(numLandunits, numRadBands, numRadTypes) {}
};
} // namespace URBANXX

#endif // URBAN_SURFACE_TYPE_IMPL_H