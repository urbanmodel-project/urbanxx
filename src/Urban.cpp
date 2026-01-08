#include "Urban.h"
#include "private/UrbanTypeImpl.h"
#include "private/AtmosphereTypeImpl.h"
#include "private/DataTypesImpl.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

// Implementation of Urban library functions
// This will be compiled into liburban.a

extern "C" {

UrbanErrorCode UrbanCreate(int numLandunits, UrbanType *urban) {
  if (urban == nullptr) {
    return URBAN_ERR_INVALID_ARGUMENT;
  }
  if (numLandunits <= 0) {
    return URBAN_ERR_INVALID_ARGUMENT;
  }

  try {
    *urban = new _p_UrbanType(numLandunits);
    return URBAN_SUCCESS;
  } catch (...) {
    return URBAN_ERR_INTERNAL;
  }
}

UrbanErrorCode UrbanDestroy(UrbanType *urban) {
  if (urban == nullptr || *urban == nullptr) {
    return URBAN_ERR_INVALID_ARGUMENT;
  }

  try {
    delete *urban;
    *urban = nullptr;
    return URBAN_SUCCESS;
  } catch (...) {
    return URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
