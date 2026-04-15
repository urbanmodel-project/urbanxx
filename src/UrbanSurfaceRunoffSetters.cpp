#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanSetterHelpers.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

extern "C" {

using namespace URBANXX;

void UrbanSetWtfactPerviousRoad(UrbanType urban, const double *values,
                                int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.Wtfact, values, length, status);
}

void UrbanSetFoverPerviousRoad(UrbanType urban, const double *values,
                               int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.Fover, values, length, status);
}

void UrbanSetFrostTablePerviousRoad(UrbanType urban, const double *values,
                                    int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.FrostTable, values, length, status);
}

void UrbanSetZwtPerchedPerviousRoad(UrbanType urban, const double *values,
                                    int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.ZwtPerched, values, length, status);
}

} // extern "C"
