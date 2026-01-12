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
void ComputeSnowAlbedo(const int l, const Real coszen, Array3DR8 roof_snowAlb,
                       Array3DR8 impRoad_snowAlb, Array3DR8 perRoad_snowAlb) {

  if (coszen > 0) {
    // Set snow albedo for both VIS and NIR bands, both direct and diffuse
    for (int ib = 0; ib < NUM_RAD_BANDS; ++ib) {
      const Real snowAlb = (ib == VIS) ? SNOW_ALBEDO_VIS : SNOW_ALBEDO_NIR;

      // Roof snow albedo
      roof_snowAlb(l, ib, DIRECT) = snowAlb;
      roof_snowAlb(l, ib, DIFFUSE) = snowAlb;

      // Impervious road snow albedo
      impRoad_snowAlb(l, ib, DIRECT) = snowAlb;
      impRoad_snowAlb(l, ib, DIFFUSE) = snowAlb;

      // Pervious road snow albedo
      perRoad_snowAlb(l, ib, DIRECT) = snowAlb;
      perRoad_snowAlb(l, ib, DIFFUSE) = snowAlb;
    }
  }
}

// Compute effective albedo combining base albedo with snow effects
KOKKOS_INLINE_FUNCTION
void ComputeCombinedAlbedo(const int l, const Real frac_sno,
                           Array3DR8 roof_snowAlb, Array3DR8 roof_baseAlb,
                           Array3DR8 roof_albWithSnow,
                           Array3DR8 impRoad_snowAlb, Array3DR8 impRoad_baseAlb,
                           Array3DR8 impRoad_albWithSnow,
                           Array3DR8 perRoad_snowAlb, Array3DR8 perRoad_baseAlb,
                           Array3DR8 perRoad_albWithSnow) {

  for (int ib = 0; ib < NUM_RAD_BANDS; ++ib) {
    // Roof: direct radiation uses direct snow albedo, diffuse uses direct snow
    // albedo
    roof_albWithSnow(l, ib, DIRECT) =
        roof_baseAlb(l, ib, DIRECT) * (1.0 - frac_sno) +
        roof_snowAlb(l, ib, DIRECT) * frac_sno;
    roof_albWithSnow(l, ib, DIFFUSE) =
        roof_baseAlb(l, ib, DIFFUSE) * (1.0 - frac_sno) +
        roof_snowAlb(l, ib, DIRECT) * frac_sno;

    // Impervious road: direct uses direct snow, diffuse uses diffuse snow
    impRoad_albWithSnow(l, ib, DIRECT) =
        impRoad_baseAlb(l, ib, DIRECT) * (1.0 - frac_sno) +
        impRoad_snowAlb(l, ib, DIRECT) * frac_sno;
    impRoad_albWithSnow(l, ib, DIFFUSE) =
        impRoad_baseAlb(l, ib, DIFFUSE) * (1.0 - frac_sno) +
        impRoad_snowAlb(l, ib, DIFFUSE) * frac_sno;

    // Pervious road: direct uses direct snow, diffuse uses diffuse snow
    perRoad_albWithSnow(l, ib, DIRECT) =
        perRoad_baseAlb(l, ib, DIRECT) * (1.0 - frac_sno) +
        perRoad_snowAlb(l, ib, DIFFUSE) * frac_sno;
    perRoad_albWithSnow(l, ib, DIFFUSE) =
        perRoad_baseAlb(l, ib, DIRECT) * (1.0 - frac_sno) +
        perRoad_snowAlb(l, ib, DIFFUSE) * frac_sno;
  }
}

// Compute incident direct and diffuse radiation in VIS and NIR bands
KOKKOS_INLINE_FUNCTION
void ComputeIncidentRadiation(const int l, const Real coszen, const Real hwr,
                              const Real vf_skyFromRoad,
                              const Real vf_skyFromWall,
                              Array3DR8 sunlitWall_downRad,
                              Array3DR8 shadedWall_downRad,
                              Array3DR8 road_downRad) {

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
      shadedWall_downRad(l, ib, DIRECT) = 0.0; // eqn. 2.15

      road_downRad(l, ib, DIRECT) =
          sdir[ib] * (2.0 * theta0OverPi -
                      2.0 / rpi * hwr * tanzen * (1.0 - costheta0)); // eqn 2.17

      sunlitWall_downRad(l, ib, DIRECT) =
          2.0 * sdir[ib] *
          ((1.0 / hwr) * (0.5 - theta0OverPi) +
           (1.0 / rpi) * tanzen * (1.0 - costheta0)); // eqn. 2.16

      // Diffuse radiation
      road_downRad(l, ib, DIFFUSE) = sdif[ib] * vf_skyFromRoad; // eqn 2.30
      sunlitWall_downRad(l, ib, DIFFUSE) =
          sdif[ib] * vf_skyFromWall; // eqn 2.32
      shadedWall_downRad(l, ib, DIFFUSE) =
          sdif[ib] * vf_skyFromWall; // eqn 2.31
    }
  }
}

// Compute net shortwave radiation for all urban surfaces
void ComputeNetShortwave(URBANXX::_p_UrbanType &urban) {
  const int numLandunits = urban.numLandunits;

  // Get references to data
  auto coszen = urban.atmosphereData.Coszen;
  auto hwr = urban.urbanParams.CanyonHwr;
  auto vf_skyFromRoad = urban.urbanParams.viewFactor.SkyFrmRoad;
  auto vf_skyFromWall = urban.urbanParams.viewFactor.SkyFrmWall;
  auto vf_wallFromRoad = urban.urbanParams.viewFactor.WallFrmRoad;
  auto vf_roadFromWall = urban.urbanParams.viewFactor.RoadFrmWall;
  auto vf_wallFromWall = urban.urbanParams.viewFactor.OtherWallFrmWall;
  auto fracPervRoad = urban.urbanParams.FracPervRoadOfTotalRoad;

  auto sunlitWall_downRad = urban.sunlitWall.DownwellingShortRad;
  auto shadedWall_downRad = urban.shadedWall.DownwellingShortRad;
  auto road_downRad = urban.compositeRoadSurface.DownwellingShortRad;

  auto roof_snowAlb = urban.roof.SnowAlbedo;
  auto impRoad_snowAlb = urban.imperviousRoad.SnowAlbedo;
  auto perRoad_snowAlb = urban.perviousRoad.SnowAlbedo;

  auto roof_baseAlb = urban.urbanParams.albedo.Roof;
  auto impRoad_baseAlb = urban.urbanParams.albedo.ImperviousRoad;
  auto perRoad_baseAlb = urban.urbanParams.albedo.PerviousRoad;
  auto sunlitWall_baseAlb = urban.urbanParams.albedo.SunlitWall;
  auto shadedWall_baseAlb = urban.urbanParams.albedo.ShadedWall;

  auto roof_albWithSnow = urban.roof.AlbedoWithSnowEffects;
  auto impRoad_albWithSnow = urban.imperviousRoad.AlbedoWithSnowEffects;
  auto perRoad_albWithSnow = urban.perviousRoad.AlbedoWithSnowEffects;

  // Snow fraction (placeholder - will be computed from snow model later)
  constexpr Real frac_sno = 0.0;

  // Compute snow albedo, combined albedo, and incident radiation for each
  // landunit
  Kokkos::parallel_for(
      "ComputeShortwaveRadiation", numLandunits, KOKKOS_LAMBDA(const int l) {
        ComputeIncidentRadiation(l, coszen(l), hwr(l), vf_skyFromRoad(l),
                                 vf_skyFromWall(l), sunlitWall_downRad,
                                 shadedWall_downRad, road_downRad);
        ComputeSnowAlbedo(l, coszen(l), roof_snowAlb, impRoad_snowAlb,
                          perRoad_snowAlb);
        ComputeCombinedAlbedo(
            l, frac_sno, roof_snowAlb, roof_baseAlb, roof_albWithSnow,
            impRoad_snowAlb, impRoad_baseAlb, impRoad_albWithSnow,
            perRoad_snowAlb, perRoad_baseAlb, perRoad_albWithSnow);

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Compute net shortwave radiation with multiple reflections
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        for (int ib = 0; ib < NUM_RAD_BANDS; ++ib) {
          for (int it = 0; it < NUM_RAD_TYPES; ++it) {

            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // Computations for roads
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            // Total shortwave downwelling to road for this band and type
            Real StotForRoad = road_downRad(l, ib, it);

            // Impervious road (weight = 1 - fraction of pervious road)
            const Real fracImpRoad = 1.0 - fracPervRoad(l);

            // Initialize impervious and pervious roads
            auto impRoad = InitializeSingleRoadShortwave(
                StotForRoad, impRoad_albWithSnow(l, ib, it), vf_skyFromRoad(l),
                vf_wallFromRoad(l), fracImpRoad);
            auto perRoad = InitializeSingleRoadShortwave(
                StotForRoad, perRoad_albWithSnow(l, ib, it), vf_skyFromRoad(l),
                vf_wallFromRoad(l), fracPervRoad(l));

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
            Real StotForSunlitWall = sunlitWall_downRad(l, ib, it);

            // Initialize sunlit wall
            auto sunlitWall = InitializeSingleWallShortwave(
                StotForSunlitWall, sunlitWall_baseAlb(l, ib, it),
                vf_skyFromWall(l), vf_roadFromWall(l), vf_wallFromWall(l));

            // Total shortwave downwelling to shaded wall for this band and
            // type
            Real StotForShadedWall = shadedWall_downRad(l, ib, it);

            // Initialize shaded wall
            auto shadedWall = InitializeSingleWallShortwave(
                StotForShadedWall, shadedWall_baseAlb(l, ib, it),
                vf_skyFromWall(l), vf_roadFromWall(l), vf_wallFromWall(l));

            // TODO: Add iteration loop
          }
        }
      });

  Kokkos::fence();

  // Debug output (disabled by default)
  if (0) {
    std::cout << "ShadedWall.DownwellingShortRad:" << std::endl;
    print_view_3d(shadedWall_downRad);
    std::cout << "SunlitWall.DownwellingShortRad:" << std::endl;
    print_view_3d(sunlitWall_downRad);
    std::cout << "CompositeRoad.DownwellingShortRad:" << std::endl;
    print_view_3d(road_downRad);

    std::cout << "Roof.BaseAlbedo:" << std::endl;
    print_view_3d(urban.urbanParams.albedo.Roof, "Roof.BaseAlbedo");
    std::cout << "ImperviousRoad.BaseAlbedo:" << std::endl;
    print_view_3d(urban.urbanParams.albedo.ImperviousRoad,
                  "ImperviousRoad.BaseAlbedo");
    std::cout << "PerviousRoad.BaseAlbedo:" << std::endl;
    print_view_3d(urban.urbanParams.albedo.PerviousRoad,
                  "PerviousRoad.BaseAlbedo");

    std::cout << "Roof.AlbedoWithSnowEffects:" << std::endl;
    print_view_3d(urban.roof.AlbedoWithSnowEffects);
    std::cout << "ImperviousRoad.AlbedoWithSnowEffects:" << std::endl;
    print_view_3d(urban.imperviousRoad.AlbedoWithSnowEffects);
    std::cout << "PerviousRoad.AlbedoWithSnowEffects:" << std::endl;
    print_view_3d(urban.perviousRoad.AlbedoWithSnowEffects);
  }
}

} // namespace URBANXX
