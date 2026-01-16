#ifndef URBAN_VIEW_FACTOR_IMPL_H
#define URBAN_VIEW_FACTOR_IMPL_H

#include "private/UrbanParamsTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// Class to compute view factors with KOKKOS_CLASS_LAMBDA
class ViewFactorComputer {
private:
  Array1DR8 CanyonHwr;
  Array1DR8 SkyFrmRoad;
  Array1DR8 SkyFrmWall;
  Array1DR8 RoadFrmWall;
  Array1DR8 WallFrmRoad;
  Array1DR8 OtherWallFrmWall;
  int numLandunits;

public:
  ViewFactorComputer(_p_UrbanType *urban)
      : CanyonHwr(urban->urbanParams.CanyonHwr),
        SkyFrmRoad(urban->urbanParams.viewFactor.SkyFrmRoad),
        SkyFrmWall(urban->urbanParams.viewFactor.SkyFrmWall),
        RoadFrmWall(urban->urbanParams.viewFactor.RoadFrmWall),
        WallFrmRoad(urban->urbanParams.viewFactor.WallFrmRoad),
        OtherWallFrmWall(urban->urbanParams.viewFactor.OtherWallFrmWall),
        numLandunits(urban->numLandunits) {}

  // Compute view factors at a single landunit - called from device code
  KOKKOS_INLINE_FUNCTION
  void computeAtLandunit(const int l) const {
    const Real hwr = CanyonHwr(l);
    const Real sqrt_term = Kokkos::sqrt(hwr * hwr + 1.0);

    SkyFrmRoad(l) = sqrt_term - hwr;                            // eqn 2.25
    WallFrmRoad(l) = 0.5 * (1.0 - SkyFrmRoad(l));               // eqn 2.27
    SkyFrmWall(l) = 0.5 * (hwr + 1.0 - sqrt_term) / hwr;        // eqn 2.24
    RoadFrmWall(l) = SkyFrmWall(l);                             // eqn 2.27
    OtherWallFrmWall(l) = 1.0 - SkyFrmWall(l) - RoadFrmWall(l); // eqn 2.28
  }

  // Main method to run computation
  void run() {
    Kokkos::parallel_for(
        "ComputingViewFactor", numLandunits,
        KOKKOS_CLASS_LAMBDA(const int l) { computeAtLandunit(l); });
    Kokkos::fence();
  }
};

} // namespace URBANXX

#endif // URBAN_VIEW_FACTOR_IMPL_H
