#ifndef URBAN_TYPE_IMPL_H
#define URBAN_TYPE_IMPL_H
#include <private/AtmosphereTypeImpl.h>
#include <private/DataTypesImpl.h>
#include <private/UrbanCanyonTypeImpl.h>
#include <private/UrbanConstants.h>
#include <private/UrbanParamsTypeImpl.h>
#include <private/UrbanSurfaceTypeImpl.h>

namespace URBANXX {

struct _p_UrbanType {
  int numLandunits;
  int numRadBands;
  int numRadTypes;
  int numLevels;

  AtmosphereType atmosphereData;
  UrbanParamsType urbanParams;
  UrbanCanyonType urbanCanyon;
  RoofDataType roof;
  RoadDataType imperviousRoad;
  RoadDataType perviousRoad;
  CompositeRoadSurfaceData compositeRoadSurface;
  WallDataType sunlitWall;
  WallDataType shadedWall;
  BuildingDataType building;

  _p_UrbanType(int numLandunits_)
      : numLandunits(numLandunits_), numRadBands(NUM_RAD_BANDS),
        numRadTypes(NUM_RAD_TYPES), numLevels(NUM_VERTICAL_LEVELS),
        atmosphereData(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES),
        urbanParams(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES,
                    NUM_VERTICAL_LEVELS),
        urbanCanyon(numLandunits_),
        roof(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES, NUM_VERTICAL_LEVELS),
        imperviousRoad(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES,
                       NUM_VERTICAL_LEVELS),
        perviousRoad(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES,
                     NUM_VERTICAL_LEVELS),
        compositeRoadSurface(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES),
        sunlitWall(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES,
                   NUM_VERTICAL_LEVELS),
        shadedWall(numLandunits_, NUM_RAD_BANDS, NUM_RAD_TYPES,
                   NUM_VERTICAL_LEVELS),
        building(numLandunits_) {}
};

} // namespace URBANXX

#endif // URBAN_TYPE_IMPL_H