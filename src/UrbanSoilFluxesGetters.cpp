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

// ============================================================================
// EflxSoilGrnd (ground heat flux) getter functions — all five surfaces
// ============================================================================

void UrbanGetEflxSoilGrndRoof(UrbanType urban, double *values, int length,
                              UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->roof.EflxSoilGrnd, values, length, status);
}

void UrbanGetEflxSoilGrndImperviousRoad(UrbanType urban, double *values,
                                        int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->imperviousRoad.EflxSoilGrnd, values, length, status);
}

void UrbanGetEflxSoilGrndPerviousRoad(UrbanType urban, double *values,
                                      int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->perviousRoad.EflxSoilGrnd, values, length, status);
}

void UrbanGetEflxSoilGrndSunlitWall(UrbanType urban, double *values, int length,
                                    UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->sunlitWall.EflxSoilGrnd, values, length, status);
}

void UrbanGetEflxSoilGrndShadedWall(UrbanType urban, double *values, int length,
                                    UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->shadedWall.EflxSoilGrnd, values, length, status);
}

// ============================================================================
// QflxEvapGrnd (liquid ground evaporation) getter functions
// ============================================================================

void UrbanGetQflxEvapGrndRoof(UrbanType urban, double *values, int length,
                              UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->roof.QflxEvapGrnd, values, length, status);
}

void UrbanGetQflxEvapGrndImperviousRoad(UrbanType urban, double *values,
                                        int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->imperviousRoad.QflxEvapGrnd, values, length, status);
}

void UrbanGetQflxEvapGrndPerviousRoad(UrbanType urban, double *values,
                                      int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->perviousRoad.QflxEvapGrnd, values, length, status);
}

// ============================================================================
// QflxSubSnow (sublimation from ice) getter functions
// ============================================================================

void UrbanGetQflxSubSnowRoof(UrbanType urban, double *values, int length,
                             UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->roof.QflxSubSnow, values, length, status);
}

void UrbanGetQflxSubSnowImperviousRoad(UrbanType urban, double *values,
                                       int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->imperviousRoad.QflxSubSnow, values, length, status);
}

void UrbanGetQflxSubSnowPerviousRoad(UrbanType urban, double *values,
                                     int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->perviousRoad.QflxSubSnow, values, length, status);
}

// ============================================================================
// QflxDewSnow (dew deposited to snow pack) getter functions
// ============================================================================

void UrbanGetQflxDewSnowRoof(UrbanType urban, double *values, int length,
                             UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->roof.QflxDewSnow, values, length, status);
}

void UrbanGetQflxDewSnowImperviousRoad(UrbanType urban, double *values,
                                       int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->imperviousRoad.QflxDewSnow, values, length, status);
}

void UrbanGetQflxDewSnowPerviousRoad(UrbanType urban, double *values,
                                     int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->perviousRoad.QflxDewSnow, values, length, status);
}

// ============================================================================
// QflxDewGrnd (dew on bare ground) getter functions
// ============================================================================

void UrbanGetQflxDewGrndRoof(UrbanType urban, double *values, int length,
                             UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->roof.QflxDewGrnd, values, length, status);
}

void UrbanGetQflxDewGrndImperviousRoad(UrbanType urban, double *values,
                                       int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->imperviousRoad.QflxDewGrnd, values, length, status);
}

void UrbanGetQflxDewGrndPerviousRoad(UrbanType urban, double *values,
                                     int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->perviousRoad.QflxDewGrnd, values, length, status);
}

} // extern "C"
