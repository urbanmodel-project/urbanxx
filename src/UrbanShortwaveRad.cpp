#include "private/AtmosphereTypeImpl.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include "private/UrbanDebugUtils.h"
#include "private/UrbanParamsTypeImpl.h"
#include "private/UrbanShortwaveRadImpl.h"
#include "private/UrbanSurfaceTypeImpl.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>
#include <cmath>
#include <iostream>

namespace URBANXX {

// Convergence parameters for shortwave radiation iteration
constexpr int SHORTWAVE_MAX_ITERATIONS = 50;
constexpr Real SHORTWAVE_CONVERGENCE_THRESHOLD = 0.00001;

// Snow albedo constants
constexpr Real SNOW_ALBEDO_VIS = 0.66;
constexpr Real SNOW_ALBEDO_NIR = 0.56;

// Structure to hold shortwave radiation components for a surface
struct SurfaceShortwaveFluxes {
  Real absorbed;          // absorbed shortwave radiation
  Real reflected;         // reflected shortwave radiation
  Real absorbedWeighted;  // absorbed * weight (for fractional surfaces)
  Real reflectedWeighted; // reflected * weight
};

// Structure to hold reflected radiation components from road
struct ReflectedShortwaveFromRoad {
  Real toSky;            // reflected radiation to sky
  Real toSunlitWall;     // reflected radiation to sunlit wall
  Real toShadedWall;     // reflected radiation to shaded wall
  Real toSkyByWt;        // reflected to sky * weight
  Real toSunlitWallByWt; // reflected to sunlit wall * weight
  Real toShadedWallByWt; // reflected to shaded wall * weight
};

// Structure to hold all radiation components for a road surface
struct RoadShortwaveRadiation {
  SurfaceShortwaveFluxes flux;
  ReflectedShortwaveFromRoad ref;
};

// Structure to hold reflected radiation components from wall
struct ReflectedShortwaveFromWall {
  Real toSky;           // reflected radiation to sky
  Real toRoad;          // reflected radiation to road
  Real toOtherWall;     // reflected radiation to other wall
  Real toSkyByWt;       // reflected to sky * weight
  Real toRoadByWt;      // reflected to road * weight
  Real toOtherWallByWt; // reflected to other wall * weight
};

// Structure to hold all radiation components for a wall surface
struct WallShortwaveRadiation {
  SurfaceShortwaveFluxes flux;
  ReflectedShortwaveFromWall ref;
};

// Helper function to compute shortwave radiation components for a surface
KOKKOS_INLINE_FUNCTION
SurfaceShortwaveFluxes ShortwaveFluxes(const Real albedo,
                                       const Real StotForSurface,
                                       const Real weight) {
  SurfaceShortwaveFluxes fluxes;

  fluxes.absorbed = (1.0 - albedo) * StotForSurface;
  fluxes.reflected = albedo * StotForSurface;

  fluxes.absorbedWeighted = fluxes.absorbed * weight;
  fluxes.reflectedWeighted = fluxes.reflected * weight;

  return fluxes;
}

// Helper function to distribute radiation from road to other surfaces
KOKKOS_INLINE_FUNCTION
ReflectedShortwaveFromRoad DistributeShortwaveFromRoad(const Real radiation,
                                                       const Real vf_sky,
                                                       const Real vf_wall,
                                                       const Real weight) {
  ReflectedShortwaveFromRoad ref;

  // Distribute radiation by view factors
  ref.toSky = radiation * vf_sky;
  ref.toSunlitWall = radiation * vf_wall;
  ref.toShadedWall = radiation * vf_wall;

  // Apply weight (fraction of total road)
  ref.toSkyByWt = ref.toSky * weight;
  ref.toSunlitWallByWt = ref.toSunlitWall * weight;
  ref.toShadedWallByWt = ref.toShadedWall * weight;

  return ref;
}

// Helper function to compute reflected radiation components from road
KOKKOS_INLINE_FUNCTION
ReflectedShortwaveFromRoad
ReflectShortwaveRoad(const Real incomingRad, const Real albedo,
                     const Real vf_sky, const Real vf_wall, const Real weight) {
  // Compute reflected radiation
  const Real reflectedRad = albedo * incomingRad;

  // Distribute using common function
  return DistributeShortwaveFromRoad(reflectedRad, vf_sky, vf_wall, weight);
}

// Helper function to initialize a single road surface (impervious or pervious)
KOKKOS_INLINE_FUNCTION
RoadShortwaveRadiation InitializeSingleRoadShortwave(const Real StotForRoad,
                                                     const Real albedo,
                                                     const Real vf_sky,
                                                     const Real vf_wall,
                                                     const Real fraction) {

  RoadShortwaveRadiation rad;
  rad.flux = ShortwaveFluxes(albedo, StotForRoad, fraction);
  rad.ref =
      ReflectShortwaveRoad(StotForRoad, albedo, vf_sky, vf_wall, fraction);
  return rad;
}

// Helper function to distribute radiation from wall to other surfaces
KOKKOS_INLINE_FUNCTION
ReflectedShortwaveFromWall DistributeShortwaveFromWall(const Real radiation,
                                                       const Real vf_sky,
                                                       const Real vf_road,
                                                       const Real vf_wall) {
  ReflectedShortwaveFromWall ref;

  // Distribute radiation by view factors
  ref.toSky = radiation * vf_sky;
  ref.toRoad = radiation * vf_road;
  ref.toOtherWall = radiation * vf_wall;

  // For walls, weight is always 1.0 (no fractional surfaces)
  ref.toSkyByWt = ref.toSky;
  ref.toRoadByWt = ref.toRoad;
  ref.toOtherWallByWt = ref.toOtherWall;

  return ref;
}

// Helper function to compute reflected radiation components from wall
KOKKOS_INLINE_FUNCTION
ReflectedShortwaveFromWall ReflectShortwaveWall(const Real incomingRad,
                                                const Real albedo,
                                                const Real vf_sky,
                                                const Real vf_road,
                                                const Real vf_wall) {
  // Compute reflected radiation
  const Real reflectedRad = albedo * incomingRad;

  // Distribute using common function
  return DistributeShortwaveFromWall(reflectedRad, vf_sky, vf_road, vf_wall);
}

// Helper function to initialize a single wall surface (sunlit or shaded)
KOKKOS_INLINE_FUNCTION
WallShortwaveRadiation InitializeSingleWallShortwave(const Real StotForWall,
                                                     const Real albedo,
                                                     const Real vf_sky,
                                                     const Real vf_road,
                                                     const Real vf_wall) {

  WallShortwaveRadiation rad;
  rad.flux = ShortwaveFluxes(albedo, StotForWall, 1.0);
  rad.ref = ReflectShortwaveWall(StotForWall, albedo, vf_sky, vf_road, vf_wall);
  return rad;
}

// Compute snow albedo for surfaces with snow coverage
KOKKOS_INLINE_FUNCTION
void ComputeSnowAlbedo(const int l, const Real coszen, Real *roof_snowAlbPtr,
                       Real *impRoad_snowAlbPtr, Real *perRoad_snowAlbPtr,
                       const int numLandunits, const int numRadBands) {

  if (coszen > 0) {
    // Set snow albedo for both VIS and NIR bands, both direct and diffuse
    for (int ib = 0; ib < NUM_RAD_BANDS; ++ib) {
      const Real snowAlb = (ib == VIS) ? SNOW_ALBEDO_VIS : SNOW_ALBEDO_NIR;

      // Index calculation: l + numLandunits * (ib + numRadBands * it)
      // Roof snow albedo
      roof_snowAlbPtr[l + numLandunits * (ib + numRadBands * DIRECT)] = snowAlb;
      roof_snowAlbPtr[l + numLandunits * (ib + numRadBands * DIFFUSE)] =
          snowAlb;

      // Impervious road snow albedo
      impRoad_snowAlbPtr[l + numLandunits * (ib + numRadBands * DIRECT)] =
          snowAlb;
      impRoad_snowAlbPtr[l + numLandunits * (ib + numRadBands * DIFFUSE)] =
          snowAlb;

      // Pervious road snow albedo
      perRoad_snowAlbPtr[l + numLandunits * (ib + numRadBands * DIRECT)] =
          snowAlb;
      perRoad_snowAlbPtr[l + numLandunits * (ib + numRadBands * DIFFUSE)] =
          snowAlb;
    }
  }
}

// Compute effective albedo combining base albedo with snow effects
KOKKOS_INLINE_FUNCTION
void ComputeCombinedAlbedo(const int l, const Real frac_sno,
                           Real *roof_snowAlbPtr, Real *roof_baseAlbPtr,
                           Real *roof_albWithSnowPtr, Real *impRoad_snowAlbPtr,
                           Real *impRoad_baseAlbPtr,
                           Real *impRoad_albWithSnowPtr,
                           Real *perRoad_snowAlbPtr, Real *perRoad_baseAlbPtr,
                           Real *perRoad_albWithSnowPtr, const int numLandunits,
                           const int numRadBands) {

  for (int ib = 0; ib < NUM_RAD_BANDS; ++ib) {
    // Index calculation: l + numLandunits * (ib + numRadBands * it)
    // Roof: direct radiation uses direct snow albedo, diffuse uses direct snow
    // albedo
    int idx_direct = l + numLandunits * (ib + numRadBands * DIRECT);
    int idx_diffuse = l + numLandunits * (ib + numRadBands * DIFFUSE);

    roof_albWithSnowPtr[idx_direct] =
        roof_baseAlbPtr[idx_direct] * (1.0 - frac_sno) +
        roof_snowAlbPtr[idx_direct] * frac_sno;
    roof_albWithSnowPtr[idx_diffuse] =
        roof_baseAlbPtr[idx_diffuse] * (1.0 - frac_sno) +
        roof_snowAlbPtr[idx_diffuse] * frac_sno;

    // Impervious road: direct uses direct snow, diffuse uses diffuse snow
    impRoad_albWithSnowPtr[idx_direct] =
        impRoad_baseAlbPtr[idx_direct] * (1.0 - frac_sno) +
        impRoad_snowAlbPtr[idx_direct] * frac_sno;
    impRoad_albWithSnowPtr[idx_diffuse] =
        impRoad_baseAlbPtr[idx_diffuse] * (1.0 - frac_sno) +
        impRoad_snowAlbPtr[idx_diffuse] * frac_sno;

    // Pervious road: direct uses direct snow, diffuse uses diffuse snow
    perRoad_albWithSnowPtr[idx_direct] =
        perRoad_baseAlbPtr[idx_direct] * (1.0 - frac_sno) +
        perRoad_snowAlbPtr[idx_direct] * frac_sno;
    perRoad_albWithSnowPtr[idx_diffuse] =
        perRoad_baseAlbPtr[idx_diffuse] * (1.0 - frac_sno) +
        perRoad_snowAlbPtr[idx_diffuse] * frac_sno;
  }
}

// Compute incident direct and diffuse radiation in VIS and NIR bands
KOKKOS_INLINE_FUNCTION
void ComputeIncidentRadiation(const int l, const Real coszen, const Real hwr,
                              const Real vf_skyFromRoad,
                              const Real vf_skyFromWall,
                              Real *sunlitWall_downRadPtr,
                              Real *shadedWall_downRadPtr,
                              Real *road_downRadPtr, const int numLandunits,
                              const int numRadBands) {

  constexpr Real rpi = M_PI;
  constexpr Real tiny = 1.0e-6;

  // Incident direct and diffuse radiation for VIS and NIR bands
  // (assumed to be unity - per unit incident flux)
  Real sdir[NUM_RAD_BANDS];
  Real sdif[NUM_RAD_BANDS];

  for (int ib = 0; ib < NUM_RAD_BANDS; ++ib) {
    sdir[ib] = 1.0;
    sdif[ib] = 1.0;
  }

  if (coszen > 0) {
    const Real zen = Kokkos::acos(coszen);
    const Real z = Kokkos::fmax(zen, tiny);
    const Real val = Kokkos::fmin(1.0 / (hwr * Kokkos::tan(z)), 1.0);
    const Real theta0 = Kokkos::asin(val);
    const Real tanzen = Kokkos::tan(zen);
    const Real costheta0 = Kokkos::cos(theta0);
    const Real theta0OverPi = theta0 / rpi;

    for (int ib = 0; ib < NUM_RAD_BANDS; ++ib) {
      // Direct radiation
      // Index calculation: l + numLandunits * (ib + numRadBands * it)
      shadedWall_downRadPtr[l + numLandunits * (ib + numRadBands * DIRECT)] =
          0.0; // eqn. 2.15

      road_downRadPtr[l + numLandunits * (ib + numRadBands * DIRECT)] =
          sdir[ib] * (2.0 * theta0OverPi -
                      2.0 / rpi * hwr * tanzen * (1.0 - costheta0)); // eqn 2.17

      sunlitWall_downRadPtr[l + numLandunits * (ib + numRadBands * DIRECT)] =
          2.0 * sdir[ib] *
          ((1.0 / hwr) * (0.5 - theta0OverPi) +
           (1.0 / rpi) * tanzen * (1.0 - costheta0)); // eqn. 2.16

      // Diffuse radiation
      road_downRadPtr[l + numLandunits * (ib + numRadBands * DIFFUSE)] =
          sdif[ib] * vf_skyFromRoad; // eqn 2.30
      sunlitWall_downRadPtr[l + numLandunits * (ib + numRadBands * DIFFUSE)] =
          sdif[ib] * vf_skyFromWall; // eqn 2.32
      shadedWall_downRadPtr[l + numLandunits * (ib + numRadBands * DIFFUSE)] =
          sdif[ib] * vf_skyFromWall; // eqn 2.31
    }
  }
}

// Compute net shortwave radiation for all urban surfaces
void ComputeNetShortwave(URBANXX::_p_UrbanType &urban) {
  const int numLandunits = urban.numLandunits;

  // Get raw pointers to avoid CUDA extended lambda issues (for 1D arrays)
  Real *coszenPtr = urban.atmosphereData.Coszen.data();
  Real *hwrPtr = urban.urbanParams.CanyonHwr.data();
  Real *vf_skyFromRoadPtr = urban.urbanParams.viewFactor.SkyFrmRoad.data();
  Real *vf_skyFromWallPtr = urban.urbanParams.viewFactor.SkyFrmWall.data();
  Real *vf_wallFromRoadPtr = urban.urbanParams.viewFactor.WallFrmRoad.data();
  Real *vf_roadFromWallPtr = urban.urbanParams.viewFactor.RoadFrmWall.data();
  Real *vf_wallFromWallPtr =
      urban.urbanParams.viewFactor.OtherWallFrmWall.data();
  Real *fracPervRoadPtr = urban.urbanParams.FracPervRoadOfTotalRoad.data();

  // For 3D arrays, use raw pointers to avoid CUDA extended lambda issues
  // Layout is LayoutLeft (column-major): index = l + numLandunits * (ib +
  // numRadBands * it)
  Real *sunlitWall_downRadPtr = urban.sunlitWall.DownwellingShortRad.data();
  Real *shadedWall_downRadPtr = urban.shadedWall.DownwellingShortRad.data();
  Real *road_downRadPtr = urban.compositeRoadSurface.DownwellingShortRad.data();

  Real *roof_snowAlbPtr = urban.roof.SnowAlbedo.data();
  Real *impRoad_snowAlbPtr = urban.imperviousRoad.SnowAlbedo.data();
  Real *perRoad_snowAlbPtr = urban.perviousRoad.SnowAlbedo.data();

  Real *roof_baseAlbPtr = urban.urbanParams.albedo.Roof.data();
  Real *impRoad_baseAlbPtr = urban.urbanParams.albedo.ImperviousRoad.data();
  Real *perRoad_baseAlbPtr = urban.urbanParams.albedo.PerviousRoad.data();
  Real *sunlitWall_baseAlbPtr = urban.urbanParams.albedo.SunlitWall.data();
  Real *shadedWall_baseAlbPtr = urban.urbanParams.albedo.ShadedWall.data();

  Real *roof_albWithSnowPtr = urban.roof.AlbedoWithSnowEffects.data();
  Real *impRoad_albWithSnowPtr =
      urban.imperviousRoad.AlbedoWithSnowEffects.data();
  Real *perRoad_albWithSnowPtr =
      urban.perviousRoad.AlbedoWithSnowEffects.data();

  // Access absorbed and reflected shortwave radiation fields (to be updated)
  Real *absImpRoadPtr = urban.imperviousRoad.AbsorbedShortRad.data();
  Real *absPerRoadPtr = urban.perviousRoad.AbsorbedShortRad.data();
  Real *absSunlitWallPtr = urban.sunlitWall.AbsorbedShortRad.data();
  Real *absShadedWallPtr = urban.shadedWall.AbsorbedShortRad.data();

  Real *refImpRoadPtr = urban.imperviousRoad.ReflectedShortRad.data();
  Real *refPerRoadPtr = urban.perviousRoad.ReflectedShortRad.data();
  Real *refSunlitWallPtr = urban.sunlitWall.ReflectedShortRad.data();
  Real *refShadedWallPtr = urban.shadedWall.ReflectedShortRad.data();

  // Get dimensions for 3D array indexing
  const int numRadBands = urban.numRadBands;
  const int numRadTypes = urban.numRadTypes;

  // Snow fraction (placeholder - will be computed from snow model later)
  constexpr Real frac_sno = 0.0;

  // Compute snow albedo, combined albedo, and incident radiation for each
  // landunit
  using ExecSpace = Kokkos::DefaultExecutionSpace;
  Kokkos::parallel_for(
      "ComputeShortwaveRadiation",
      Kokkos::RangePolicy<ExecSpace>(0, numLandunits),
      KOKKOS_LAMBDA(const int l) {
        // Call helper functions with raw pointers
        ComputeIncidentRadiation(l, coszenPtr[l], hwrPtr[l],
                                 vf_skyFromRoadPtr[l], vf_skyFromWallPtr[l],
                                 sunlitWall_downRadPtr, shadedWall_downRadPtr,
                                 road_downRadPtr, numLandunits, numRadBands);
        ComputeSnowAlbedo(l, coszenPtr[l], roof_snowAlbPtr, impRoad_snowAlbPtr,
                          perRoad_snowAlbPtr, numLandunits, numRadBands);
        ComputeCombinedAlbedo(
            l, frac_sno, roof_snowAlbPtr, roof_baseAlbPtr, roof_albWithSnowPtr,
            impRoad_snowAlbPtr, impRoad_baseAlbPtr, impRoad_albWithSnowPtr,
            perRoad_snowAlbPtr, perRoad_baseAlbPtr, perRoad_albWithSnowPtr,
            numLandunits, numRadBands);

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Compute net shortwave radiation with multiple reflections
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        for (int ib = 0; ib < NUM_RAD_BANDS; ++ib) {
          for (int it = 0; it < NUM_RAD_TYPES; ++it) {
            // Index calculation for 3D arrays: l + numLandunits * (ib +
            // numRadBands * it)
            const int idx = l + numLandunits * (ib + numRadBands * it);

            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // Computations for roads
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            // Total shortwave downwelling to road for this band and type
            Real StotForRoad = road_downRadPtr[idx];

            // Impervious road (weight = 1 - fraction of pervious road)
            const Real fracImpRoad = 1.0 - fracPervRoadPtr[l];

            // Initialize impervious and pervious roads
            auto impRoad = InitializeSingleRoadShortwave(
                StotForRoad, impRoad_albWithSnowPtr[idx], vf_skyFromRoadPtr[l],
                vf_wallFromRoadPtr[l], fracImpRoad);
            auto perRoad = InitializeSingleRoadShortwave(
                StotForRoad, perRoad_albWithSnowPtr[idx], vf_skyFromRoadPtr[l],
                vf_wallFromRoadPtr[l], fracPervRoadPtr[l]);

            // Combine both roads
            Real RoadAbs =
                impRoad.flux.absorbedWeighted + perRoad.flux.absorbedWeighted;
            Real RoadRef =
                impRoad.flux.reflectedWeighted + perRoad.flux.reflectedWeighted;

            Real RoadRefToSky = impRoad.ref.toSkyByWt + perRoad.ref.toSkyByWt;
            Real RoadRefToSunlitWall =
                impRoad.ref.toSunlitWallByWt + perRoad.ref.toSunlitWallByWt;
            Real RoadRefToShadedWall =
                impRoad.ref.toShadedWallByWt + perRoad.ref.toShadedWallByWt;

            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // Computations for walls
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            // Total shortwave downwelling to sunlit wall for this band and
            // type
            Real StotForSunlitWall = sunlitWall_downRadPtr[idx];

            // Initialize sunlit wall
            auto sunlitWall = InitializeSingleWallShortwave(
                StotForSunlitWall, sunlitWall_baseAlbPtr[idx],
                vf_skyFromWallPtr[l], vf_roadFromWallPtr[l],
                vf_wallFromWallPtr[l]);

            // Total shortwave downwelling to shaded wall for this band and
            // type
            Real StotForShadedWall = shadedWall_downRadPtr[idx];

            // Initialize shaded wall
            auto shadedWall = InitializeSingleWallShortwave(
                StotForShadedWall, shadedWall_baseAlbPtr[idx],
                vf_skyFromWallPtr[l], vf_roadFromWallPtr[l],
                vf_wallFromWallPtr[l]);

            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // Initialize cumulative absorbed and reflected radiation for all
            // surfaces
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            absImpRoadPtr[idx] = impRoad.flux.absorbed;
            absPerRoadPtr[idx] = perRoad.flux.absorbed;
            absSunlitWallPtr[idx] = sunlitWall.flux.absorbed;
            absShadedWallPtr[idx] = shadedWall.flux.absorbed;

            refImpRoadPtr[idx] = impRoad.ref.toSky;
            refPerRoadPtr[idx] = perRoad.ref.toSky;
            refSunlitWallPtr[idx] = sunlitWall.ref.toSky;
            refShadedWallPtr[idx] = shadedWall.ref.toSky;

            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // Iteration loop for multiple reflections between surfaces
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            for (int iter = 0; iter < SHORTWAVE_MAX_ITERATIONS; ++iter) {

              // step(1): Compute incoming radiation from wall-to-road and
              // wall-to-wall reflections

              // For roads: incoming from walls
              StotForRoad =
                  (sunlitWall.ref.toRoad + shadedWall.ref.toRoad) * hwrPtr[l];

              impRoad.flux = ShortwaveFluxes(impRoad_albWithSnowPtr[idx],
                                             StotForRoad, fracImpRoad);
              perRoad.flux = ShortwaveFluxes(perRoad_albWithSnowPtr[idx],
                                             StotForRoad, fracPervRoadPtr[l]);

              RoadAbs =
                  impRoad.flux.absorbedWeighted + perRoad.flux.absorbedWeighted;
              RoadRef = impRoad.flux.reflectedWeighted +
                        perRoad.flux.reflectedWeighted;

              // For sunlit wall: incoming from roads and shaded wall
              StotForSunlitWall =
                  RoadRefToSunlitWall / hwrPtr[l] + shadedWall.ref.toOtherWall;
              sunlitWall.flux = ShortwaveFluxes(sunlitWall_baseAlbPtr[idx],
                                                StotForSunlitWall, 1.0);

              // For shaded wall: incoming from roads and sunlit wall
              StotForShadedWall =
                  RoadRefToShadedWall / hwrPtr[l] + sunlitWall.ref.toOtherWall;
              shadedWall.flux = ShortwaveFluxes(shadedWall_baseAlbPtr[idx],
                                                StotForShadedWall, 1.0);

              // step(2): Update cumulative absorbed radiation
              absImpRoadPtr[idx] += impRoad.flux.absorbed;
              absPerRoadPtr[idx] += perRoad.flux.absorbed;
              absSunlitWallPtr[idx] += sunlitWall.flux.absorbed;
              absShadedWallPtr[idx] += shadedWall.flux.absorbed;

              // step(3): Compute reflected radiation components for this
              // iteration
              impRoad.ref = ReflectShortwaveRoad(
                  StotForRoad, impRoad_albWithSnowPtr[idx],
                  vf_skyFromRoadPtr[l], vf_wallFromRoadPtr[l], fracImpRoad);
              perRoad.ref = ReflectShortwaveRoad(
                  StotForRoad, perRoad_albWithSnowPtr[idx],
                  vf_skyFromRoadPtr[l], vf_wallFromRoadPtr[l],
                  fracPervRoadPtr[l]);

              RoadRefToSky = impRoad.ref.toSkyByWt + perRoad.ref.toSkyByWt;
              RoadRefToSunlitWall =
                  impRoad.ref.toSunlitWallByWt + perRoad.ref.toSunlitWallByWt;
              RoadRefToShadedWall =
                  impRoad.ref.toShadedWallByWt + perRoad.ref.toShadedWallByWt;

              sunlitWall.ref = ReflectShortwaveWall(
                  StotForSunlitWall, sunlitWall_baseAlbPtr[idx],
                  vf_skyFromWallPtr[l], vf_roadFromWallPtr[l],
                  vf_wallFromWallPtr[l]);

              shadedWall.ref = ReflectShortwaveWall(
                  StotForShadedWall, shadedWall_baseAlbPtr[idx],
                  vf_skyFromWallPtr[l], vf_roadFromWallPtr[l],
                  vf_wallFromWallPtr[l]);

              // step(4): Update cumulative reflected radiation to sky
              refImpRoadPtr[idx] += impRoad.ref.toSky;
              refPerRoadPtr[idx] += perRoad.ref.toSky;
              refSunlitWallPtr[idx] += sunlitWall.ref.toSky;
              refShadedWallPtr[idx] += shadedWall.ref.toSky;

              // Check convergence
              Real convergence_criteria =
                  Kokkos::max(RoadAbs, Kokkos::max(sunlitWall.flux.absorbed,
                                                   shadedWall.flux.absorbed));
              if (convergence_criteria < SHORTWAVE_CONVERGENCE_THRESHOLD) {
                break;
              }
            }
          }
        }
      });

  Kokkos::fence();

  // Debug output disabled - Views have been converted to raw pointers
  // If needed, access views directly from urban structure instead
}

} // namespace URBANXX
