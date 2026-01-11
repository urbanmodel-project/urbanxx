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
  URBAN_ERR_INTERNAL = 3,
  URBAN_ERR_SIZE_MISMATCH = 4
} UrbanErrorCode;

// API functions
URBAN_EXTERN void UrbanCreate(int numLandunits, UrbanType *urban,
                              UrbanErrorCode *status);
URBAN_EXTERN void UrbanDestroy(UrbanType *urban, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetCanyonHwr(UrbanType urban, const double *values,
                                    int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetFracPervRoadOfTotalRoad(UrbanType urban,
                                                  const double *values,
                                                  int length,
                                                  UrbanErrorCode *status);

// Albedo setter functions
URBAN_EXTERN void UrbanSetAlbedoPerviousRoad(UrbanType urban,
                                             const double *values,
                                             const int size[3],
                                             UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAlbedoImperviousRoad(UrbanType urban,
                                               const double *values,
                                               const int size[3],
                                               UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAlbedoSunlitWall(UrbanType urban,
                                           const double *values,
                                           const int size[3],
                                           UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAlbedoShadedWall(UrbanType urban,
                                           const double *values,
                                           const int size[3],
                                           UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAlbedoRoof(UrbanType urban, const double *values,
                                     const int size[3], UrbanErrorCode *status);

// Emissivity setter functions
URBAN_EXTERN void UrbanSetEmissivityPerviousRoad(UrbanType urban,
                                                 const double *values,
                                                 int length,
                                                 UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetEmissivityImperviousRoad(UrbanType urban,
                                                   const double *values,
                                                   int length,
                                                   UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetEmissivityWall(UrbanType urban, const double *values,
                                         int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetEmissivityRoof(UrbanType urban, const double *values,
                                         int length, UrbanErrorCode *status);

// Thermal conductivity setter functions
URBAN_EXTERN void UrbanSetThermalConductivityRoad(UrbanType urban,
                                                  const double *values,
                                                  int length,
                                                  UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetThermalConductivityWall(UrbanType urban,
                                                  const double *values,
                                                  int length,
                                                  UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetThermalConductivityRoof(UrbanType urban,
                                                  const double *values,
                                                  int length,
                                                  UrbanErrorCode *status);

// Heat capacity setter functions
URBAN_EXTERN void UrbanSetHeatCapacityRoad(UrbanType urban,
                                           const double *values, int length,
                                           UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetHeatCapacityWall(UrbanType urban,
                                           const double *values, int length,
                                           UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetHeatCapacityRoof(UrbanType urban,
                                           const double *values, int length,
                                           UrbanErrorCode *status);

// Initialization functions
URBAN_EXTERN void UrbanInitializeTemperature(UrbanType urban,
                                             UrbanErrorCode *status);

// Time-stepping functions
URBAN_EXTERN void UrbanAdvance(UrbanType urban, UrbanErrorCode *status);

// Atmospheric forcing setter functions
URBAN_EXTERN void UrbanSetAtmTemp(UrbanType urban, const double *values,
                                  int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAtmPotTemp(UrbanType urban, const double *values,
                                     int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAtmRho(UrbanType urban, const double *values,
                                 int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAtmSpcHumd(UrbanType urban, const double *values,
                                     int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAtmPress(UrbanType urban, const double *values,
                                   int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAtmWindU(UrbanType urban, const double *values,
                                   int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAtmWindV(UrbanType urban, const double *values,
                                   int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAtmCoszen(UrbanType urban, const double *values,
                                    int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAtmFracSnow(UrbanType urban, const double *values,
                                      int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAtmLongwaveDown(UrbanType urban, const double *values,
                                          int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetAtmShortwaveDown(UrbanType urban,
                                           const double *values,
                                           const int size[3],
                                           UrbanErrorCode *status);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // URBAN_H
