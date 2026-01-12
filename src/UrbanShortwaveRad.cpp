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

  auto sunlitWall_downRad = urban.sunlitWall.DownwellingShortRad;
  auto shadedWall_downRad = urban.shadedWall.DownwellingShortRad;
  auto road_downRad = urban.compositeRoadSurface.DownwellingShortRad;

  auto roof_snowAlb = urban.roof.SnowAlbedo;
  auto impRoad_snowAlb = urban.imperviousRoad.SnowAlbedo;
  auto perRoad_snowAlb = urban.perviousRoad.SnowAlbedo;

  auto roof_baseAlb = urban.urbanParams.albedo.Roof;
  auto impRoad_baseAlb = urban.urbanParams.albedo.ImperviousRoad;
  auto perRoad_baseAlb = urban.urbanParams.albedo.PerviousRoad;

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
