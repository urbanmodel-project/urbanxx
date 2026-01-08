#include "Urban.h"
#include "private/AtmosphereTypeImpl.h"
#include "private/DataTypesImpl.h"
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

} // extern "C"
