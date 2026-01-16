#include "Urban.h"
#include "private/AtmosphereTypeImpl.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanLongwaveRadImpl.h"
#include "private/UrbanShortwaveRadImpl.h"
#include "private/UrbanSurfaceFluxesImpl.h"
#include "private/UrbanTypeImpl.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

// Implementation of Urban library functions
// This will be compiled into liburban.a

extern "C" {

void UrbanCreate(int numLandunits, UrbanType *urban, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }
  if (numLandunits <= 0) {
    *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    *urban = new _p_UrbanType(numLandunits);
    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

void UrbanDestroy(UrbanType *urban, UrbanErrorCode *status) {
  if (urban == nullptr || *urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    delete *urban;
    *urban = nullptr;
    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

void UrbanAdvance(UrbanType urban, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    // Compute net longwave radiation
    printf("DEBUG: About to call ComputeNetLongwave...\n");
    fflush(stdout);
    URBANXX::ComputeNetLongwave(*urban);
    printf("DEBUG: ComputeNetLongwave completed\n");
    fflush(stdout);

    // Compute net shortwave radiation
    printf("DEBUG: About to call ComputeNetShortwave...\n");
    fflush(stdout);
    URBANXX::ComputeNetShortwave(*urban);
    printf("DEBUG: ComputeNetShortwave completed\n");
    fflush(stdout);

    // Compute surface fluxes
    printf("DEBUG: About to call ComputeSurfaceFluxes...\n");
    fflush(stdout);
    URBANXX::ComputeSurfaceFluxes(*urban);
    printf("DEBUG: ComputeSurfaceFluxes completed\n");
    fflush(stdout);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
