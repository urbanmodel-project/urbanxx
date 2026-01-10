#include "Urban.h"
#include "private/AtmosphereTypeImpl.h"
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

void UrbanSetAtmTemp(UrbanType urban, const double *values, int length,
                     UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcTemp, values, length, status);
}

void UrbanSetAtmPotTemp(UrbanType urban, const double *values, int length,
                        UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcPotTemp, values, length, status);
}

void UrbanSetAtmRho(UrbanType urban, const double *values, int length,
                    UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcRho, values, length, status);
}

void UrbanSetAtmSpcHumd(UrbanType urban, const double *values, int length,
                        UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcSpcHumd, values, length, status);
}

void UrbanSetAtmPress(UrbanType urban, const double *values, int length,
                      UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcPress, values, length, status);
}

void UrbanSetAtmWindU(UrbanType urban, const double *values, int length,
                      UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcWindU, values, length, status);
}

void UrbanSetAtmWindV(UrbanType urban, const double *values, int length,
                      UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcWindV, values, length, status);
}

void UrbanSetAtmCoszen(UrbanType urban, const double *values, int length,
                       UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.Coszen, values, length, status);
}

void UrbanSetAtmFracSnow(UrbanType urban, const double *values, int length,
                         UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.FracSnow, values, length, status);
}

void UrbanSetAtmLongwaveDown(UrbanType urban, const double *values, int length,
                             UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcLRad, values, length, status);
}

void UrbanSetAtmShortwaveDown(UrbanType urban, const double *values,
                              const int size[3], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView3D(urban->atmosphereData.ForcSRad, values, size, status);
}

} // extern "C"
