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

// Absorbed shortwave radiation getter functions
void UrbanGetAbsorbedShortwaveRoof(UrbanType urban, double *values,
                                   const int size[3], UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView3D(urban->roof.AbsorbedShortRad, values, size, status);
}

void UrbanGetAbsorbedShortwaveImperviousRoad(UrbanType urban, double *values,
                                             const int size[3],
                                             UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView3D(urban->imperviousRoad.AbsorbedShortRad, values, size, status);
}

void UrbanGetAbsorbedShortwavePerviousRoad(UrbanType urban, double *values,
                                           const int size[3],
                                           UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView3D(urban->perviousRoad.AbsorbedShortRad, values, size, status);
}

void UrbanGetAbsorbedShortwaveSunlitWall(UrbanType urban, double *values,
                                         const int size[3],
                                         UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView3D(urban->sunlitWall.AbsorbedShortRad, values, size, status);
}

void UrbanGetAbsorbedShortwaveShadedWall(UrbanType urban, double *values,
                                         const int size[3],
                                         UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView3D(urban->shadedWall.AbsorbedShortRad, values, size, status);
}

// Reflected shortwave radiation getter functions
void UrbanGetReflectedShortwaveRoof(UrbanType urban, double *values,
                                    const int size[3], UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView3D(urban->roof.ReflectedShortRad, values, size, status);
}

void UrbanGetReflectedShortwaveImperviousRoad(UrbanType urban, double *values,
                                              const int size[3],
                                              UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView3D(urban->imperviousRoad.ReflectedShortRad, values, size, status);
}

void UrbanGetReflectedShortwavePerviousRoad(UrbanType urban, double *values,
                                            const int size[3],
                                            UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView3D(urban->perviousRoad.ReflectedShortRad, values, size, status);
}

void UrbanGetReflectedShortwaveSunlitWall(UrbanType urban, double *values,
                                          const int size[3],
                                          UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView3D(urban->sunlitWall.ReflectedShortRad, values, size, status);
}

void UrbanGetReflectedShortwaveShadedWall(UrbanType urban, double *values,
                                          const int size[3],
                                          UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView3D(urban->shadedWall.ReflectedShortRad, values, size, status);
}

} // extern "C"
