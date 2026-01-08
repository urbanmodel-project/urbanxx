#ifndef URBAN_TYPE_IMPL_H
#define URBAN_TYPE_IMPL_H
#include <private/AtmosphereTypeImpl.h>
#include <private/DataTypesImpl.h>

namespace URBANXX {

// Constants
constexpr int NUM_RAD_BANDS = 2;

struct _p_UrbanType {
  int numLandunits;
  int numRadBands;

  AtmosphereType atmosphereData;

  _p_UrbanType(int numLandunits_)
      : numLandunits(numLandunits_), numRadBands(NUM_RAD_BANDS),
        atmosphereData(numLandunits_, NUM_RAD_BANDS) {}
};

} // namespace URBANXX

#endif // URBAN_TYPE_IMPL_H