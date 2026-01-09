#ifndef URBAN_PARAMS_TYPE_IMPL_H
#define URBAN_PARAMS_TYPE_IMPL_H

#include "private/DataTypesImpl.h"
#include "private/UrbanMacros.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

struct Albedo {
  DECLARE_DEVICE_VIEW(2DR8, PerviousRoadDir) // albedo for pervious road for
                                             // direct beam solar radiation (-)
  DECLARE_DEVICE_VIEW(2DR8, PerviousRoadDif) // albedo for pervious road for
                                             // diffuse solar radiation (-)
  DECLARE_DEVICE_VIEW(2DR8,
                      ImperviousRoadDir) // albedo for impervious road for
                                         // direct beam solar radiation (-)
  DECLARE_DEVICE_VIEW(2DR8, ImperviousRoadDif) // albedo for impervious road for
                                               // diffuse solar radiation (-)
  DECLARE_DEVICE_VIEW(2DR8, SunlitWallDir) // albedo for sunlit wall for direct
                                           // beam solar radiation (-)
  DECLARE_DEVICE_VIEW(
      2DR8,
      SunlitWallDif) // albedo for sunlit wall for diffuse solar radiation (-)
  DECLARE_DEVICE_VIEW(2DR8, ShadedWallDir) // albedo for shaded wall for direct
                                           // beam solar radiation (-)
  DECLARE_DEVICE_VIEW(
      2DR8,
      ShadedWallDif) // albedo for shaded wall for diffuse solar radiation (-)
  DECLARE_DEVICE_VIEW(
      2DR8, RoofDir) // albedo for roof for direct beam solar radiation (-)
  DECLARE_DEVICE_VIEW(
      2DR8, RoofDif) // albedo for roof for diffuse solar radiation (-)

  Albedo(int numLandunits, int numRadBands) {
    ALLOCATE_DEVICE_VIEW(PerviousRoadDir, Array2DR8, numLandunits, numRadBands)
    ALLOCATE_DEVICE_VIEW(PerviousRoadDif, Array2DR8, numLandunits, numRadBands)
    ALLOCATE_DEVICE_VIEW(ImperviousRoadDir, Array2DR8, numLandunits,
                         numRadBands)
    ALLOCATE_DEVICE_VIEW(ImperviousRoadDif, Array2DR8, numLandunits,
                         numRadBands)
    ALLOCATE_DEVICE_VIEW(SunlitWallDir, Array2DR8, numLandunits, numRadBands)
    ALLOCATE_DEVICE_VIEW(SunlitWallDif, Array2DR8, numLandunits, numRadBands)
    ALLOCATE_DEVICE_VIEW(ShadedWallDir, Array2DR8, numLandunits, numRadBands)
    ALLOCATE_DEVICE_VIEW(ShadedWallDif, Array2DR8, numLandunits, numRadBands)
    ALLOCATE_DEVICE_VIEW(RoofDir, Array2DR8, numLandunits, numRadBands)
    ALLOCATE_DEVICE_VIEW(RoofDif, Array2DR8, numLandunits, numRadBands)
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

struct PropertiesForSurfaces {
  DECLARE_DEVICE_VIEW(1DR8, Road) // property value for road material
  DECLARE_DEVICE_VIEW(1DR8, Wall) // property value for wall material
  DECLARE_DEVICE_VIEW(1DR8, Roof) // property value for roof material

  PropertiesForSurfaces(int numLandunits) {
    ALLOCATE_DEVICE_VIEW(Road, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Wall, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Roof, Array1DR8, numLandunits)
  }
};

struct ViewFactor {
  DECLARE_DEVICE_VIEW(1DR8, SkyFrmRoad)  // view factor of sky from road (-)
  DECLARE_DEVICE_VIEW(1DR8, SkyFrmWall)  // view factor of sky from wall (-)
  DECLARE_DEVICE_VIEW(1DR8, WallFrmRoad) // view factor of wall from road (-)
  DECLARE_DEVICE_VIEW(
      1DR8, OthrWallFrmWall) // view factor of other wall from wall (-)

  ViewFactor(int numLandunits) {
    ALLOCATE_DEVICE_VIEW(SkyFrmRoad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(SkyFrmWall, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(WallFrmRoad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(OthrWallFrmWall, Array1DR8, numLandunits)
  }
};

struct UrbanParamsType {
  DECLARE_DUAL_VIEWS(1DR8, CanyonHwr) // canyon height-to-width ratio (-)

  ViewFactor viewFactor;
  PropertiesForSurfaces tk; // thermal conductivity (W/m/K)
  PropertiesForSurfaces cv; // heat capacity (J/m^3/K)
  Albedo albedo;            // albedo for various urban surfaces
  Emissivity emissivity;    // emissivity for various urban surfaces

  UrbanParamsType(int numLandunits, int numRadBands)
      : viewFactor(numLandunits), tk(numLandunits), cv(numLandunits),
        albedo(numLandunits, numRadBands), emissivity(numLandunits) {
    ALLOCATE_DUAL_VIEWS(CanyonHwr, 1DR8, numLandunits)
  }
};
} // namespace URBANXX

#endif // URBAN_PARAMS_TYPE_IMPL_H