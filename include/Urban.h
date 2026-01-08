// Public C API for the Urban library
#ifndef URBAN_H
#define URBAN_H

#include <stddef.h> // size_t

#ifdef __cplusplus
extern "C" {
#endif

// Visibility / linkage macro for C API functions
#ifdef __cplusplus
#define URBAN_EXTERN extern "C"
#else
#define URBAN_EXTERN extern
#endif

#ifndef __cplusplus
#include <stdbool.h>
#endif

// Opaque handle type (pointer to private struct)
typedef struct _p_UrbanType *UrbanType;

// Error codes for API calls
typedef enum {
  URBAN_SUCCESS = 0,
  URBAN_ERR_INVALID_ARGUMENT = 1,
  URBAN_ERR_NOT_INITIALIZED = 2,
  URBAN_ERR_INTERNAL = 3
} UrbanErrorCode;

// API functions
URBAN_EXTERN void UrbanCreate(int numLandunits, UrbanType *urban,
                              UrbanErrorCode *status);
URBAN_EXTERN void UrbanDestroy(UrbanType *urban, UrbanErrorCode *status);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // URBAN_H
