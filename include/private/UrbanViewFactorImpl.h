#ifndef URBAN_VIEW_FACTOR_IMPL_H
#define URBAN_VIEW_FACTOR_IMPL_H

#include "private/UrbanParamsTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// POD struct to hold view factor views - trivially copyable for CUDA lambdas
struct ViewFactorViews {
  Array1DR8 CanyonHwr;
  Array1DR8 SkyFrmRoad;
  Array1DR8 SkyFrmWall;
  Array1DR8 RoadFrmWall;
  Array1DR8 WallFrmRoad;
  Array1DR8 OtherWallFrmWall;
};

// Class to compute view factors with KOKKOS_CLASS_LAMBDA
class ViewFactorComputer {
private:
  ViewFactorViews views;
  int numLandunits;

public:
  ViewFactorComputer(_p_UrbanType *urban)
      : views{urban->urbanParams.CanyonHwr,
              urban->urbanParams.viewFactor.SkyFrmRoad,
              urban->urbanParams.viewFactor.SkyFrmWall,
              urban->urbanParams.viewFactor.RoadFrmWall,
              urban->urbanParams.viewFactor.WallFrmRoad,
              urban->urbanParams.viewFactor.OtherWallFrmWall},
        numLandunits(urban->numLandunits) {}

  // Compute view factors at a single landunit - called from device code
  KOKKOS_INLINE_FUNCTION
  static void computeAtLandunit(const int l, const ViewFactorViews &v) {
    const Real hwr = v.CanyonHwr(l);
    const Real sqrt_term = Kokkos::sqrt(hwr * hwr + 1.0);

    v.SkyFrmRoad(l) = sqrt_term - hwr;                                  // eqn 2.25
    v.WallFrmRoad(l) = 0.5 * (1.0 - v.SkyFrmRoad(l));               // eqn 2.27
    v.SkyFrmWall(l) = 0.5 * (hwr + 1.0 - sqrt_term) / hwr;              // eqn 2.24
    v.RoadFrmWall(l) = v.SkyFrmWall(l);                             // eqn 2.27
    v.OtherWallFrmWall(l) = 1.0 - v.SkyFrmWall(l) - v.RoadFrmWall(l); // eqn 2.28
  }

  // Main method to run computation
  void run() {
    // Capture POD struct by value (not *this)
    ViewFactorViews v = views;
    Kokkos::parallel_for(
        "ComputingViewFactor", numLandunits,
        KOKKOS_LAMBDA(const int l) { computeAtLandunit(l, v); });
    Kokkos::fence();
  }
};

} // namespace URBANXX

#endif // URBAN_VIEW_FACTOR_IMPL_H

