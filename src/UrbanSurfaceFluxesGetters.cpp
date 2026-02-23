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

// Sensible heat flux (EflxShGrnd) getter functions
void UrbanGetSensibleHeatFluxRoof(UrbanType urban, double *values, int length,
                                  UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->roof.EflxShGrnd, values, length, status);
}

void UrbanGetSensibleHeatFluxImperviousRoad(UrbanType urban, double *values,
                                            int length,
                                            UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->imperviousRoad.EflxShGrnd, values, length, status);
}

void UrbanGetSensibleHeatFluxPerviousRoad(UrbanType urban, double *values,
                                          int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->perviousRoad.EflxShGrnd, values, length, status);
}

void UrbanGetSensibleHeatFluxSunlitWall(UrbanType urban, double *values,
                                        int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->sunlitWall.EflxShGrnd, values, length, status);
}

void UrbanGetSensibleHeatFluxShadedWall(UrbanType urban, double *values,
                                        int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->shadedWall.EflxShGrnd, values, length, status);
}

// Soil evaporation flux (QflxEvapSoil) getter functions
void UrbanGetEvapFluxRoof(UrbanType urban, double *values, int length,
                          UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->roof.QflxEvapSoil, values, length, status);
}

void UrbanGetEvapFluxImperviousRoad(UrbanType urban, double *values, int length,
                                    UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->imperviousRoad.QflxEvapSoil, values, length, status);
}

void UrbanGetEvapFluxPerviousRoad(UrbanType urban, double *values, int length,
                                  UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->perviousRoad.QflxEvapSoil, values, length, status);
}

} // extern "C"
