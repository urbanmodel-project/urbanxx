#ifndef URBAN_LONGWAVE_RAD_IMPL_H
#define URBAN_LONGWAVE_RAD_IMPL_H

#include "private/UrbanConstants.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// Convergence parameters
constexpr int LONGWAVE_MAX_ITERATIONS = 50;
constexpr Real LONGWAVE_CONVERGENCE_THRESHOLD = 0.0001;

// Forward declarations of helper structures
struct SurfaceLongwaveFluxes;
struct ReflectedRadFromRoad;
struct ReflectedRadFromWall;
struct RoadRadiation;
struct WallRadiation;
struct RoadViewFactors;
struct WallViewFactors;

// Class to compute net longwave radiation with KOKKOS_CLASS_LAMBDA
class LongwaveRadiationComputer {
private:
  // Atmospheric forcing
  Array1DR8 forcLRad;

  // View factors
  Array1DR8 vf_skyFromRoad;
  Array1DR8 vf_skyFromWall;
  Array1DR8 vf_wallFromRoad;
  Array1DR8 vf_roadFromWall;
  Array1DR8 vf_otherWallFromWall;
  Array1DR8 canyonHwr;
  Array1DR8 fracPervRoad;

  // Emissivities
  Array1DR8 emissRoof;
  Array1DR8 emissWall;
  Array1DR8 emissImpRoad;
  Array1DR8 emissPerRoad;

  // Temperatures
  Array1DR8 tempRoof;
  Array1DR8 tempSunlitWall;
  Array1DR8 tempShadedWall;
  Array1DR8 tempImpRoad;
  Array1DR8 tempPerRoad;

  // Output: net longwave
  Array1DR8 netLwSunlitWall;
  Array1DR8 netLwShadedWall;
  Array1DR8 netLwImpRoad;
  Array1DR8 netLwPerRoad;

  // Output: upward longwave
  Array1DR8 upLwSunlitWall;
  Array1DR8 upLwShadedWall;
  Array1DR8 upLwImpRoad;
  Array1DR8 upLwPerRoad;

  int numLandunits;

public:
  LongwaveRadiationComputer(_p_UrbanType *urban)
      : forcLRad(urban->atmosphereData.ForcLRad),
        vf_skyFromRoad(urban->urbanParams.viewFactor.SkyFrmRoad),
        vf_skyFromWall(urban->urbanParams.viewFactor.SkyFrmWall),
        vf_wallFromRoad(urban->urbanParams.viewFactor.WallFrmRoad),
        vf_roadFromWall(urban->urbanParams.viewFactor.RoadFrmWall),
        vf_otherWallFromWall(urban->urbanParams.viewFactor.OtherWallFrmWall),
        canyonHwr(urban->urbanParams.CanyonHwr),
        fracPervRoad(urban->urbanParams.FracPervRoadOfTotalRoad),
        emissRoof(urban->urbanParams.emissivity.Roof),
        emissWall(urban->urbanParams.emissivity.Wall),
        emissImpRoad(urban->urbanParams.emissivity.ImperviousRoad),
        emissPerRoad(urban->urbanParams.emissivity.PerviousRoad),
        tempRoof(urban->roof.Temperature),
        tempSunlitWall(urban->sunlitWall.Temperature),
        tempShadedWall(urban->shadedWall.Temperature),
        tempImpRoad(urban->imperviousRoad.Temperature),
        tempPerRoad(urban->perviousRoad.Temperature),
        netLwSunlitWall(urban->sunlitWall.NetLongRad),
        netLwShadedWall(urban->shadedWall.NetLongRad),
        netLwImpRoad(urban->imperviousRoad.NetLongRad),
        netLwPerRoad(urban->perviousRoad.NetLongRad),
        upLwSunlitWall(urban->sunlitWall.UpwardLongRad),
        upLwShadedWall(urban->shadedWall.UpwardLongRad),
        upLwImpRoad(urban->imperviousRoad.UpwardLongRad),
        upLwPerRoad(urban->perviousRoad.UpwardLongRad),
        numLandunits(urban->numLandunits) {}

  // Helper method - performs computation for one landunit
  KOKKOS_INLINE_FUNCTION
  void computeAtLandunit(const int l) const;

  // Main computation method
  void run();
};

// Original function interface (now uses class internally)
void ComputeNetLongwave(URBANXX::_p_UrbanType &urban);

} // namespace URBANXX

#endif // URBAN_LONGWAVE_RAD_IMPL_H
