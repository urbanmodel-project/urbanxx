#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanSetterHelpers.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

void UrbanSetNMeltRoof(UrbanType urban, const double *values, int length,
                       UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->roof.NMelt, values, length, status);
}

void UrbanSetNMeltImperviousRoad(UrbanType urban, const double *values,
                                 int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->imperviousRoad.NMelt, values, length, status);
}

void UrbanSetNMeltPerviousRoad(UrbanType urban, const double *values,
                               int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.NMelt, values, length, status);
}
