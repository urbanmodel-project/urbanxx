#ifndef URBAN_TYPE_IMPL_H
#define URBAN_TYPE_IMPL_H
#include <private/AtmosphereTypeImpl.h>
#include <private/DataTypesImpl.h>
#include <private/UrbanConstants.h>
#include <private/UrbanParamsTypeImpl.h>
#include <private/UrbanSurfaceTypeImpl.h>

namespace URBANXX {

struct _p_UrbanType {
  int numLandunits;
  int numRadBands;
  int numRadTypes;

  AtmosphereType atmosphereData;
  UrbanParamsType urbanParams;
  RoofDataType roof;
  RoadDataType imperviousRoad;
  RoadDataType perviousRoad;
  CompositeRoadSurfaceData compositeRoadSurface;
  WallDataType sunlitWall;
  WallDataType shadedWall;

  _p_UrbanType(int numLandunits_)
      : numLandunits(numLandunits_), numRadBands(NUM_RAD_BANDS),
        numRadTypes(NUM_RAD_TYPES),
        atmosphereData(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES),
        urbanParams(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES),
        roof(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES),
        imperviousRoad(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES),
        perviousRoad(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES),
        compositeRoadSurface(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES),
        sunlitWall(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES),
        shadedWall(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES) {}
};

} // namespace URBANXX

#endif // URBAN_TYPE_IMPL_H