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

  // Compute incident radiation for each landunit
  Kokkos::parallel_for(
      "ComputeIncidentRadiation", numLandunits, KOKKOS_LAMBDA(const int l) {
        ComputeIncidentRadiation(l, coszen(l), hwr(l), vf_skyFromRoad(l),
                                 vf_skyFromWall(l), sunlitWall_downRad,
                                 shadedWall_downRad, road_downRad);
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
  }
}

} // namespace URBANXX
