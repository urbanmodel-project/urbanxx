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

// Hydrology getter functions (pervious road outputs of ComputeHydrology)

void UrbanGetSoilLiquidWaterPerviousRoad(UrbanType urban, double *values,
                                         const int size[2],
                                         UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  GetView2D(urban->perviousRoad.H2OSoiLiq, values, size, status);
}

void UrbanGetSoilVolumetricWaterPerviousRoad(UrbanType urban, double *values,
                                             const int size[2],
                                             UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  GetView2D(urban->perviousRoad.H2OSoiVol, values, size, status);
}

void UrbanGetAquiferRechargeRatePerviousRoad(UrbanType urban, double *values,
                                             int length,
                                             UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->perviousRoad.Qcharge, values, length, status);
}

void UrbanGetWaterDeficitFluxPerviousRoad(UrbanType urban, double *values,
                                          int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->perviousRoad.QflxDeficit, values, length, status);
}

} // extern "C"
