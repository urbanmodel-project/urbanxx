#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanGetterHelpers.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

extern "C" {

using namespace URBANXX;

// Surface runoff (QflxSurf) getter functions — all five surfaces

void UrbanGetQflxSurfRoof(UrbanType urban, double *values, int length,
                          UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->roof.QflxSurf, values, length, status);
}

void UrbanGetQflxSurfImperviousRoad(UrbanType urban, double *values, int length,
                                    UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->imperviousRoad.QflxSurf, values, length, status);
}

void UrbanGetQflxSurfPerviousRoad(UrbanType urban, double *values, int length,
                                  UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->perviousRoad.QflxSurf, values, length, status);
}

void UrbanGetQflxSurfSunlitWall(UrbanType urban, double *values, int length,
                                UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->sunlitWall.QflxSurf, values, length, status);
}

void UrbanGetQflxSurfShadedWall(UrbanType urban, double *values, int length,
                                UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->shadedWall.QflxSurf, values, length, status);
}

// Updated TopH2OSoiLiq after ponding correction (roof and impervious road)

void UrbanGetTopH2OSoiLiqRoof(UrbanType urban, double *values, int length,
                              UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->roof.TopH2OSoiLiq, values, length, status);
}

void UrbanGetTopH2OSoiLiqImperviousRoad(UrbanType urban, double *values,
                                        int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->imperviousRoad.TopH2OSoiLiq, values, length, status);
}

} // extern "C"
