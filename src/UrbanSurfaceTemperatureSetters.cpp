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

void UrbanSetEffectiveSurfTempRoof(UrbanType urban, const double *values,
                                   int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->roof.EffectiveSurfTemp, values, length, status);
}

void UrbanSetEffectiveSurfTempImperviousRoad(UrbanType urban,
                                             const double *values, int length,
                                             UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->imperviousRoad.EffectiveSurfTemp, values, length, status);
}

void UrbanSetEffectiveSurfTempPerviousRoad(UrbanType urban,
                                           const double *values, int length,
                                           UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.EffectiveSurfTemp, values, length, status);
}

void UrbanSetEffectiveSurfTempSunlitWall(UrbanType urban, const double *values,
                                         int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->sunlitWall.EffectiveSurfTemp, values, length, status);
}

void UrbanSetEffectiveSurfTempShadedWall(UrbanType urban, const double *values,
                                         int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->shadedWall.EffectiveSurfTemp, values, length, status);
}

} // extern "C"
