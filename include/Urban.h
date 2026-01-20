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
URBAN_EXTERN void UrbanSetWtRoof(UrbanType urban, const double *values,
                                 int length, UrbanErrorCode *status);

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

// Temperature setter functions
URBAN_EXTERN void UrbanSetTemperatureRoof(UrbanType urban, const double *values,
                                          int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetTemperatureImperviousRoad(UrbanType urban,
                                                    const double *values,
                                                    int length,
                                                    UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetTemperaturePerviousRoad(UrbanType urban,
                                                  const double *values,
                                                  int length,
                                                  UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetTemperatureSunlitWall(UrbanType urban,
                                                const double *values,
                                                int length,
                                                UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetTemperatureShadedWall(UrbanType urban,
                                                const double *values,
                                                int length,
                                                UrbanErrorCode *status);

// Height parameter setter functions
URBAN_EXTERN void UrbanSetForcHgtT(UrbanType urban, const double *values,
                                   int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetForcHgtU(UrbanType urban, const double *values,
                                   int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetZDTown(UrbanType urban, const double *values,
                                 int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetZ0Town(UrbanType urban, const double *values,
                                 int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetHtRoof(UrbanType urban, const double *values,
                                 int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetWindHgtCanyon(UrbanType urban, const double *values,
                                        int length, UrbanErrorCode *status);

// Initialization functions
URBAN_EXTERN void UrbanInitializeTemperature(UrbanType urban,
                                             UrbanErrorCode *status);

// Time-stepping functions
URBAN_EXTERN void UrbanAdvance(UrbanType urban, UrbanErrorCode *status);

// Physics computation functions
URBAN_EXTERN void UrbanComputeNetLongwave(UrbanType urban,
                                          UrbanErrorCode *status);
URBAN_EXTERN void UrbanComputeNetShortwave(UrbanType urban,
                                           UrbanErrorCode *status);
URBAN_EXTERN void UrbanComputeSurfaceFluxes(UrbanType urban,
                                            UrbanErrorCode *status);

// Kokkos utility functions
URBAN_EXTERN bool UrbanKokkosIsLayoutRight(void);
URBAN_EXTERN bool UrbanKokkosIsLayoutLeft(void);

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

// Shortwave radiation getter functions - Absorbed
URBAN_EXTERN void UrbanGetAbsorbedShortwaveRoof(UrbanType urban, double *values,
                                                const int size[3],
                                                UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetAbsorbedShortwaveImperviousRoad(
    UrbanType urban, double *values, const int size[3], UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetAbsorbedShortwavePerviousRoad(UrbanType urban,
                                                        double *values,
                                                        const int size[3],
                                                        UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetAbsorbedShortwaveSunlitWall(UrbanType urban,
                                                      double *values,
                                                      const int size[3],
                                                      UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetAbsorbedShortwaveShadedWall(UrbanType urban,
                                                      double *values,
                                                      const int size[3],
                                                      UrbanErrorCode *status);

// Shortwave radiation getter functions - Reflected
URBAN_EXTERN void UrbanGetReflectedShortwaveRoof(UrbanType urban,
                                                 double *values,
                                                 const int size[3],
                                                 UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetReflectedShortwaveImperviousRoad(
    UrbanType urban, double *values, const int size[3], UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetReflectedShortwavePerviousRoad(
    UrbanType urban, double *values, const int size[3], UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetReflectedShortwaveSunlitWall(UrbanType urban,
                                                       double *values,
                                                       const int size[3],
                                                       UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetReflectedShortwaveShadedWall(UrbanType urban,
                                                       double *values,
                                                       const int size[3],
                                                       UrbanErrorCode *status);

// Longwave radiation getter functions - Net
URBAN_EXTERN void UrbanGetNetLongwaveRoof(UrbanType urban, double *values,
                                          int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetNetLongwaveImperviousRoad(UrbanType urban,
                                                    double *values, int length,
                                                    UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetNetLongwavePerviousRoad(UrbanType urban,
                                                  double *values, int length,
                                                  UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetNetLongwaveSunlitWall(UrbanType urban, double *values,
                                                int length,
                                                UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetNetLongwaveShadedWall(UrbanType urban, double *values,
                                                int length,
                                                UrbanErrorCode *status);

// Longwave radiation getter functions - Upward
URBAN_EXTERN void UrbanGetUpwardLongwaveRoof(UrbanType urban, double *values,
                                             int length,
                                             UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetUpwardLongwaveImperviousRoad(UrbanType urban,
                                                       double *values,
                                                       int length,
                                                       UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetUpwardLongwavePerviousRoad(UrbanType urban,
                                                     double *values, int length,
                                                     UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetUpwardLongwaveSunlitWall(UrbanType urban,
                                                   double *values, int length,
                                                   UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetUpwardLongwaveShadedWall(UrbanType urban,
                                                   double *values, int length,
                                                   UrbanErrorCode *status);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // URBAN_H
