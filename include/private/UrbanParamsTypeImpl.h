#ifndef URBAN_PARAMS_TYPE_IMPL_H
#define URBAN_PARAMS_TYPE_IMPL_H

#include "private/DataTypesImpl.h"
#include "private/KokkosViewMacros.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

struct Albedo {
  DECLARE_DEVICE_VIEW(
      3DR8, PerviousRoad) // albedo for pervious road [landunit, band, type]
  DECLARE_DEVICE_VIEW(
      3DR8, ImperviousRoad) // albedo for impervious road [landunit, band, type]
  DECLARE_DEVICE_VIEW(
      3DR8, SunlitWall) // albedo for sunlit wall [landunit, band, type]
  DECLARE_DEVICE_VIEW(
      3DR8, ShadedWall) // albedo for shaded wall [landunit, band, type]
  DECLARE_DEVICE_VIEW(3DR8, Roof) // albedo for roof [landunit, band, type]

  Albedo(int numLandunits, int numRadBands, int numRadTypes) {
    ALLOCATE_DEVICE_VIEW(PerviousRoad, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
    ALLOCATE_DEVICE_VIEW(ImperviousRoad, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
    ALLOCATE_DEVICE_VIEW(SunlitWall, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
    ALLOCATE_DEVICE_VIEW(ShadedWall, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
    ALLOCATE_DEVICE_VIEW(Roof, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
  }
};

struct Emissivity {
  DECLARE_DEVICE_VIEW(1DR8, PerviousRoad)   // emissivity for road material (-)
  DECLARE_DEVICE_VIEW(1DR8, ImperviousRoad) // emissivity for road material (-)
  DECLARE_DEVICE_VIEW(1DR8, Wall)           // emissivity for wall material (-)
  DECLARE_DEVICE_VIEW(1DR8, Roof)           // emissivity for roof material (-)

  Emissivity(int numLandunits) {
    ALLOCATE_DEVICE_VIEW(PerviousRoad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ImperviousRoad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Wall, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Roof, Array1DR8, numLandunits)
  }
};

struct CommonSurfaceProperties {
  DECLARE_DEVICE_VIEW(
      2DR8, Road) // property value for road material (landunit, level)
  DECLARE_DEVICE_VIEW(
      2DR8, Wall) // property value for wall material (landunit, level)
  DECLARE_DEVICE_VIEW(
      2DR8, Roof) // property value for roof material (landunit, level)

  CommonSurfaceProperties(int numLandunits, int numUrbanLayers) {
    ALLOCATE_DEVICE_VIEW(Road, Array2DR8, numLandunits, numUrbanLayers)
    ALLOCATE_DEVICE_VIEW(Wall, Array2DR8, numLandunits, numUrbanLayers)
    ALLOCATE_DEVICE_VIEW(Roof, Array2DR8, numLandunits, numUrbanLayers)
  }
};

struct ViewFactor {
  DECLARE_DEVICE_VIEW(1DR8, SkyFrmRoad)  // view factor of sky from road (-)
  DECLARE_DEVICE_VIEW(1DR8, SkyFrmWall)  // view factor of sky from wall (-)
  DECLARE_DEVICE_VIEW(1DR8, WallFrmRoad) // view factor of wall from road (-)
  DECLARE_DEVICE_VIEW(1DR8,
                      RoadFrmWall) // view factor of road from wall (-)
  DECLARE_DEVICE_VIEW(
      1DR8, OtherWallFrmWall) // view factor of other wall from wall (-)

  ViewFactor(int numLandunits) {
    ALLOCATE_DEVICE_VIEW(SkyFrmRoad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(SkyFrmWall, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(RoadFrmWall, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(WallFrmRoad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(OtherWallFrmWall, Array1DR8, numLandunits)
  }
};

struct HeightParameters {
  DECLARE_DEVICE_VIEW(1DR8,
                      ForcHgtT) // observational height for temperature (m)
  DECLARE_DEVICE_VIEW(1DR8, ForcHgtU) // observational height for wind (m)
  DECLARE_DEVICE_VIEW(1DR8, ZDTown)   // displacement height (m)
  DECLARE_DEVICE_VIEW(1DR8, Z0Town)   // momentum roughness length (m)
  DECLARE_DEVICE_VIEW(1DR8, HtRoof)   // height of roof (m)
  DECLARE_DEVICE_VIEW(
      1DR8,
      WindHgtCanyon) // height above road at which canyon wind is computed (m)

  HeightParameters(int numLandunits) {
    ALLOCATE_DEVICE_VIEW(ForcHgtT, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ForcHgtU, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ZDTown, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Z0Town, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(HtRoof, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(WindHgtCanyon, Array1DR8, numLandunits)
  }
};

struct BuildingParameters {
  DECLARE_DEVICE_VIEW(
      1DR8, MaxTemperature) // maximum allowable interior temperature (K)
  DECLARE_DEVICE_VIEW(
      1DR8, MinTemperature) // minimum allowable interior temperature (K)
  DECLARE_DEVICE_VIEW(1DR8, WallThickness) // wall thickness (m)
  DECLARE_DEVICE_VIEW(1DR8, RoofThickness) // roof thickness (m)

  BuildingParameters(int numLandunits) {
    ALLOCATE_DEVICE_VIEW(MaxTemperature, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(MinTemperature, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(WallThickness, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(RoofThickness, Array1DR8, numLandunits)
  }
};

struct UrbanParamsType {
  DECLARE_DEVICE_VIEW(1DR8, CanyonHwr) // canyon height-to-width ratio (-)
  DECLARE_DEVICE_VIEW(1DR8,
                      FracPervRoadOfTotalRoad) // fraction of pervious road
                                               // w.r.t. total road (-)
  DECLARE_DEVICE_VIEW(1DR8,
                      WtRoof) // weight of roof w.r.t. total urban area (-)

  ViewFactor viewFactor;
  CommonSurfaceProperties tk; // thermal conductivity (W/m/K)
  CommonSurfaceProperties cv; // heat capacity (J/m^3/K)
  Albedo albedo;              // albedo for various urban surfaces
  Emissivity emissivity;      // emissivity for various urban surfaces
  HeightParameters heights;   // height parameters for surface flux calculations
  BuildingParameters building; // building parameters

  UrbanParamsType(int numLandunits, int numRadBands, int numRadTypes,
                  int numUrbanLayers)
      : viewFactor(numLandunits), tk(numLandunits, numUrbanLayers),
        cv(numLandunits, numUrbanLayers),
        albedo(numLandunits, numRadBands, numRadTypes),
        emissivity(numLandunits), heights(numLandunits),
        building(numLandunits) {
    ALLOCATE_DEVICE_VIEW(CanyonHwr, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(FracPervRoadOfTotalRoad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(WtRoof, Array1DR8, numLandunits)
  }
};
} // namespace URBANXX

#endif // URBAN_PARAMS_TYPE_IMPL_H