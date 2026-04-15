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

// Top-layer soil water setter functions (for soil flux partitioning on
// roof and impervious road surfaces).  Pervious road does not need these
// setters because UrbanComputeSoilFluxes reads H2OSoiLiq(l,0)/H2OSoiIce(l,0)
// directly from the 2D hydrology array set by UrbanSetSoilLiquidWater.

void UrbanSetTopH2OSoiLiqRoof(UrbanType urban, const double *values, int length,
                              UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->roof.TopH2OSoiLiq, values, length, status);
}

void UrbanSetTopH2OSoiIceRoof(UrbanType urban, const double *values, int length,
                              UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->roof.TopH2OSoiIce, values, length, status);
}

void UrbanSetTopH2OSoiLiqImperviousRoad(UrbanType urban, const double *values,
                                        int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->imperviousRoad.TopH2OSoiLiq, values, length, status);
}

void UrbanSetTopH2OSoiIceImperviousRoad(UrbanType urban, const double *values,
                                        int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->imperviousRoad.TopH2OSoiIce, values, length, status);
}

} // extern "C"
