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

struct BuildingDataType {
  DECLARE_DEVICE_VIEW(1DR8, Temperature;)    // building temperature (K)
  DECLARE_DEVICE_VIEW(1DR8, EFlxForHeating;) // building heat flux (W/m^2)
  DECLARE_DEVICE_VIEW(1DR8, EFluxForAC;)     // building cool flux (W/m^2)

  BuildingDataType(int numLandunits) {
    ALLOCATE_DEVICE_VIEW(Temperature, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(EFlxForHeating, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(EFluxForAC, Array1DR8, numLandunits)
  }
};

//
// struct hierarchy:
//
// SurfaceDataBase (base)
// ├── SnowCoveredSurfaceData
// │   ├── RoadDataType
// │   └── RoofDataType
// └── WallDataType
//

// Base struct for all surface types with common radiative properties
struct SurfaceDataBase {
  DECLARE_DEVICE_VIEW(3DR8,
                      ReflectedShortRad) // reflected solar radiation per unit
                                         // surface area per unit incident flux
  DECLARE_DEVICE_VIEW(3DR8,
                      AbsorbedShortRad) // absorbed solar radiation per unit
                                        // surface area per unit incident flux
  DECLARE_DEVICE_VIEW(1DR8, Emissivity) // emissivity
  DECLARE_DEVICE_VIEW(1DR8, EffectiveSurfTemp) // effective temperature
  DECLARE_DEVICE_VIEW(1DR8, NetLongRad)        // net longwave radiation
  DECLARE_DEVICE_VIEW(1DR8, UpwardLongRad)     // upward longwave radiation

  // Saturation humidity variables
  DECLARE_DEVICE_VIEW(1DR8, Es)   // saturation vapor pressure (Pa)
  DECLARE_DEVICE_VIEW(1DR8, EsdT) // d(es)/dT (Pa/K)
  DECLARE_DEVICE_VIEW(1DR8, Qs)   // saturation specific humidity (kg/kg)
  DECLARE_DEVICE_VIEW(1DR8, QsdT) // d(qs)/dT (1/K)

  // Information about layers
  DECLARE_DEVICE_VIEW(2DR8, Zc)         // depth at the center of layer (m)
  DECLARE_DEVICE_VIEW(2DR8, Zi)         // depth at the layer interface (m)
  DECLARE_DEVICE_VIEW(2DR8, Dz)         // layer thickness (m)
  DECLARE_DEVICE_VIEW(1DR8, TotalDepth) // total depth of surface layers (m)

  // Thermal properties
  DECLARE_DEVICE_VIEW(2DR8, ThermalConductivity) // thermal conductivity (W/m/K)
  DECLARE_DEVICE_VIEW(2DR8, HeatCapacity)        // volumetric heat capacity (J/m^3/K)

  SurfaceDataBase(int numLandunits, int numRadBands, int numRadTypes,
                  int numLayers) {
    ALLOCATE_DEVICE_VIEW(ReflectedShortRad, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
    ALLOCATE_DEVICE_VIEW(AbsorbedShortRad, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
    ALLOCATE_DEVICE_VIEW(Emissivity, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(EffectiveSurfTemp, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(NetLongRad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(UpwardLongRad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Es, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(EsdT, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Qs, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(QsdT, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Zc, Array2DR8, numLandunits, numLayers)
    ALLOCATE_DEVICE_VIEW(Zi, Array2DR8, numLandunits, numLayers + 1)
    ALLOCATE_DEVICE_VIEW(Dz, Array2DR8, numLandunits, numLayers)
    ALLOCATE_DEVICE_VIEW(TotalDepth, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ThermalConductivity, Array2DR8, numLandunits, numLayers)
    ALLOCATE_DEVICE_VIEW(HeatCapacity, Array2DR8, numLandunits, numLayers)
  }
};

// For surfaces with snow coverage (roads and roofs)
struct SnowCoveredSurfaceData : SurfaceDataBase {
  DECLARE_DEVICE_VIEW(3DR8, SnowAlbedo) // snow albedo
  DECLARE_DEVICE_VIEW(3DR8,
                      AlbedoWithSnowEffects) // albedo including snow effects

  SnowCoveredSurfaceData(int numLandunits, int numRadBands, int numRadTypes,
                         int numLayers)
      : SurfaceDataBase(numLandunits, numRadBands, numRadTypes, numLayers) {
    ALLOCATE_DEVICE_VIEW(SnowAlbedo, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
    ALLOCATE_DEVICE_VIEW(AlbedoWithSnowEffects, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
  }
};

struct RoadDataType : SnowCoveredSurfaceData {
  // Inherits all fields from SnowCoveredSurfaceData and SurfaceDataBase
  RoadDataType(int numLandunits, int numRadBands, int numRadTypes,
               int numLayers)
      : SnowCoveredSurfaceData(numLandunits, numRadBands, numRadTypes,
                               numLayers) {}
};

struct WallDataType : SurfaceDataBase {
  DECLARE_DEVICE_VIEW(
      3DR8,
      DownwellingShortRad) // downwelling shortwave radiation per unit wall area
  // Inherits common radiative fields from SurfaceDataBase

  WallDataType(int numLandunits, int numRadBands, int numRadTypes,
               int numLayers)
      : SurfaceDataBase(numLandunits, numRadBands, numRadTypes, numLayers) {
    ALLOCATE_DEVICE_VIEW(DownwellingShortRad, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
  }
};

struct RoofDataType : SnowCoveredSurfaceData {
  // Inherits all fields from SnowCoveredSurfaceData and SurfaceDataBase

  RoofDataType(int numLandunits, int numRadBands, int numRadTypes,
               int numLayers)
      : SnowCoveredSurfaceData(numLandunits, numRadBands, numRadTypes,
                               numLayers) {}
};
} // namespace URBANXX

#endif // URBAN_SURFACE_TYPE_IMPL_H