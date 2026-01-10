#ifndef URBAN_VALIDATION_H
#define URBAN_VALIDATION_H

#include "Urban.h"

// Inline validation functions to reduce boilerplate
// These are header-only inline functions, so no separate .cpp file needed

static inline void ValidateInputs(UrbanType urban, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
  }
}

static inline bool ValidateInputsWithData(UrbanType urban, const void *data,
                                          UrbanErrorCode *status) {
  if (urban == nullptr || data == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return false;
  }
  return true;
}

static inline bool ValidateInputsWithSize(UrbanType urban, const void *data,
                                          const int *size,
                                          UrbanErrorCode *status) {
  if (urban == nullptr || data == nullptr || size == nullptr ||
      status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return false;
  }
  return true;
}

#endif // URBAN_VALIDATION_H
