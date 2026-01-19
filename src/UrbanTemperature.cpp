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

// Temperature setter functions
void UrbanSetTemperatureRoof(UrbanType urban, const double *values, int length,
                             UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->roof.Temperature, values, length, status);
}

void UrbanSetTemperatureImperviousRoad(UrbanType urban, const double *values,
                                       int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->imperviousRoad.Temperature, values, length, status);
}

void UrbanSetTemperaturePerviousRoad(UrbanType urban, const double *values,
                                     int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.Temperature, values, length, status);
}

void UrbanSetTemperatureSunlitWall(UrbanType urban, const double *values,
                                   int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->sunlitWall.Temperature, values, length, status);
}

void UrbanSetTemperatureShadedWall(UrbanType urban, const double *values,
                                   int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->shadedWall.Temperature, values, length, status);
}

} // extern "C"
