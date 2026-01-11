#include "private/AtmosphereTypeImpl.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include "private/UrbanLongwaveRadImpl.h"
#include "private/UrbanParamsTypeImpl.h"
#include "private/UrbanSurfaceTypeImpl.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// Structure to hold longwave radiation components for a surface
struct SurfaceLongwaveFluxes {
  Real absorbed;          // absorbed longwave radiation
  Real reflected;         // reflected longwave radiation
  Real emitted;           // emitted longwave radiation
  Real absorbedWeighted;  // absorbed * weight (for fractional surfaces)
  Real reflectedWeighted; // reflected * weight
  Real emittedWeighted;   // emitted * weight
};

// Helper function to compute longwave radiation components for a surface
KOKKOS_INLINE_FUNCTION
SurfaceLongwaveFluxes ComputeAbsRefEmiRadiation(const Real emissivity,
                                                const Real temperature,
                                                const Real LtotForSurface,
                                                const Real weight) {
  SurfaceLongwaveFluxes fluxes;

  fluxes.absorbed = emissivity * LtotForSurface;
  fluxes.reflected = (1.0 - emissivity) * LtotForSurface;
  fluxes.emitted = emissivity * STEBOL * Kokkos::pow(temperature, 4.0);

  fluxes.absorbedWeighted = fluxes.absorbed * weight;
  fluxes.reflectedWeighted = fluxes.reflected * weight;
  fluxes.emittedWeighted = fluxes.emitted * weight;

  return fluxes;
}

void ComputeNetLongwave(URBANXX::_p_UrbanType &urban) {
  // Get number of landunits for parallel execution
  const int numLandunits = urban.numLandunits;

  // Access atmospheric forcing data
  auto &forcTemp = urban.atmosphereData.ForcTemp;
  auto &forcLRad = urban.atmosphereData.ForcLRad;

  // Access urban parameters - view factors and canyon height-to-width ratio
  auto &vf_sr = urban.urbanParams.viewFactor.SkyFrmRoad;
  auto &vf_sw = urban.urbanParams.viewFactor.SkyFrmWall;
  auto &hwr = urban.urbanParams.CanyonHwr;
  auto &fracPervRoad = urban.urbanParams.FracPervRoadOfTotalRoad;

  // Access urban parameters - emissivities
  auto &emissRoof = urban.urbanParams.emissivity.Roof;
  auto &emissWall = urban.urbanParams.emissivity.Wall;
  auto &emissImpRoad = urban.urbanParams.emissivity.ImperviousRoad;
  auto &emissPerRoad = urban.urbanParams.emissivity.PerviousRoad;

  // Access surface temperatures
  auto &tempRoof = urban.roof.Temperature;
  auto &tempSunlitWall = urban.sunlitWall.Temperature;
  auto &tempShadedWall = urban.shadedWall.Temperature;
  auto &tempImpRoad = urban.imperviousRoad.Temperature;
  auto &tempPerRoad = urban.perviousRoad.Temperature;

  // Access net longwave radiation fields (to be updated)
  auto &netLwRoof = urban.roof.NetLongRad;
  auto &netLwSunlitWall = urban.sunlitWall.NetLongRad;
  auto &netLwShadedWall = urban.shadedWall.NetLongRad;
  auto &netLwImpRoad = urban.imperviousRoad.NetLongRad;
  auto &netLwPerRoad = urban.perviousRoad.NetLongRad;

  Kokkos::parallel_for(
      "ComputeNetLongwave", numLandunits, KOKKOS_LAMBDA(const int l) {
        // Net longwave calculation will go here
        // For now, set to zero as placeholder
        netLwRoof(l) = 0.0;
        netLwSunlitWall(l) = 0.0;
        netLwShadedWall(l) = 0.0;
        netLwImpRoad(l) = 0.0;
        netLwPerRoad(l) = 0.0;

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Computations for roads
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        // Total longwave downwelling to road
        const Real LtotForRoad = forcLRad(l) * vf_sr(l);

        // Impervious road (weight = 1 - fraction of pervious road)
        const Real fracImpRoad = 1.0 - fracPervRoad(l);
        auto fluxImpRoad = ComputeAbsRefEmiRadiation(
            emissImpRoad(l), tempImpRoad(l), LtotForRoad, fracImpRoad);

        // Pervious road
        auto fluxPerRoad = ComputeAbsRefEmiRadiation(
            emissPerRoad(l), tempPerRoad(l), LtotForRoad, fracPervRoad(l));

        // Combinging data from the roads
        Real RoadAbs =
            fluxImpRoad.absorbedWeighted + fluxPerRoad.absorbedWeighted;
        Real RoadRef =
            fluxImpRoad.reflectedWeighted + fluxPerRoad.reflectedWeighted;
        Real RoadEmi =
            fluxImpRoad.emittedWeighted + fluxPerRoad.emittedWeighted;

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Computations for walls
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        // Total longwave downwelling to wall
        const Real LtotForWall = forcLRad(l) * vf_sw(l);

        // Sunlit wall
        auto fluxSunlitWall = ComputeAbsRefEmiRadiation(
            emissWall(l), tempSunlitWall(l), LtotForWall, 1.0);

        // Shaded wall
        auto fluxShadedWall = ComputeAbsRefEmiRadiation(
            emissWall(l), tempShadedWall(l), LtotForWall, 1.0);
      });

  Kokkos::fence();
}

} // namespace URBANXX
