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

// Layer temperature getter functions (output of ComputeHeatDiffusion)
void UrbanGetLayerTempRoof(UrbanType urban, double *values, const int size[2],
                           UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  GetView2D(urban->roof.Temperature, values, size, status);
}

void UrbanGetLayerTempImperviousRoad(UrbanType urban, double *values,
                                     const int size[2],
                                     UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  GetView2D(urban->imperviousRoad.Temperature, values, size, status);
}

void UrbanGetLayerTempPerviousRoad(UrbanType urban, double *values,
                                   const int size[2], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  GetView2D(urban->perviousRoad.Temperature, values, size, status);
}

void UrbanGetLayerTempSunlitWall(UrbanType urban, double *values,
                                 const int size[2], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  GetView2D(urban->sunlitWall.Temperature, values, size, status);
}

void UrbanGetLayerTempShadedWall(UrbanType urban, double *values,
                                 const int size[2], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  GetView2D(urban->shadedWall.Temperature, values, size, status);
}

void UrbanGetBuildingTemperature(UrbanType urban, double *values, int length,
                                 UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->building.Temperature, values, length, status);
}

} // extern "C"
