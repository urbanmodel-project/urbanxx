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

// Structure to hold reflected radiation components from road
struct ReflectedRadFromRoad {
  Real toSky;           // reflected radiation to sky
  Real toSunwall;       // reflected radiation to sunlit wall
  Real toShadewall;     // reflected radiation to shaded wall
  Real toSkyByWt;       // reflected to sky * weight
  Real toSunwallByWt;   // reflected to sunlit wall * weight
  Real toShadewallByWt; // reflected to shaded wall * weight
};

// Structure to hold reflected radiation components from wall
struct ReflectedRadFromWall {
  Real toSky;           // reflected radiation to sky
  Real toRoad;          // reflected radiation to road
  Real toOtherwall;     // reflected radiation to other wall
  Real toSkyByWt;       // reflected to sky * weight
  Real toRoadByWt;      // reflected to road * weight
  Real toOtherwallByWt; // reflected to other wall * weight
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

// Helper function to distribute radiation from road to other surfaces
KOKKOS_INLINE_FUNCTION
ReflectedRadFromRoad DistributeRadiationFromRoad(const Real radiation,
                                                 const Real vf_sr,
                                                 const Real vf_wr,
                                                 const Real weight) {
  ReflectedRadFromRoad ref;

  // Distribute radiation by view factors
  ref.toSky = radiation * vf_sr;
  ref.toSunwall = radiation * vf_wr;
  ref.toShadewall = radiation * vf_wr;

  // Apply weight (fraction of total road)
  ref.toSkyByWt = ref.toSky * weight;
  ref.toSunwallByWt = ref.toSunwall * weight;
  ref.toShadewallByWt = ref.toShadewall * weight;

  return ref;
}

// Helper function to compute reflected radiation components from road
KOKKOS_INLINE_FUNCTION
ReflectedRadFromRoad ComputeReflectedRadFromRoad(const Real incomingRad,
                                                 const Real emissivity,
                                                 const Real vf_sr,
                                                 const Real vf_wr,
                                                 const Real weight) {
  // Compute reflected radiation
  const Real reflectedRad = (1.0 - emissivity) * incomingRad;

  // Distribute using common function
  return DistributeRadiationFromRoad(reflectedRad, vf_sr, vf_wr, weight);
}

// Helper function to compute emitted radiation components from road
KOKKOS_INLINE_FUNCTION
ReflectedRadFromRoad ComputeEmittedRadFromRoad(const Real emissivity,
                                               const Real temperature,
                                               const Real vf_sr,
                                               const Real vf_wr,
                                               const Real weight) {
  // Compute emitted radiation
  const Real emittedRad = emissivity * STEBOL * Kokkos::pow(temperature, 4.0);

  // Distribute using common function
  return DistributeRadiationFromRoad(emittedRad, vf_sr, vf_wr, weight);
}

// Helper function to distribute radiation from wall to other surfaces
KOKKOS_INLINE_FUNCTION
ReflectedRadFromWall DistributeRadiationFromWall(const Real radiation,
                                                 const Real vf_sw,
                                                 const Real vf_rw,
                                                 const Real vf_ww) {
  ReflectedRadFromWall ref;

  // Distribute radiation by view factors
  ref.toSky = radiation * vf_sw;
  ref.toRoad = radiation * vf_rw;
  ref.toOtherwall = radiation * vf_ww;

  // For walls, weight is always 1.0 (no fractional surfaces)
  ref.toSkyByWt = ref.toSky;
  ref.toRoadByWt = ref.toRoad;
  ref.toOtherwallByWt = ref.toOtherwall;

  return ref;
}

// Helper function to compute reflected radiation components from wall
KOKKOS_INLINE_FUNCTION
ReflectedRadFromWall ComputeReflectedRadFromWall(const Real incomingRad,
                                                 const Real emissivity,
                                                 const Real vf_sw,
                                                 const Real vf_rw,
                                                 const Real vf_ww) {
  // Compute reflected radiation
  const Real reflectedRad = (1.0 - emissivity) * incomingRad;

  // Distribute using common function
  return DistributeRadiationFromWall(reflectedRad, vf_sw, vf_rw, vf_ww);
}

// Helper function to compute emitted radiation components from wall
KOKKOS_INLINE_FUNCTION
ReflectedRadFromWall ComputeEmittedRadFromWall(const Real emissivity,
                                               const Real temperature,
                                               const Real vf_sw,
                                               const Real vf_rw,
                                               const Real vf_ww) {
  // Compute emitted radiation
  const Real emittedRad = emissivity * STEBOL * Kokkos::pow(temperature, 4.0);

  // Distribute using common function
  return DistributeRadiationFromWall(emittedRad, vf_sw, vf_rw, vf_ww);
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
  auto &vf_wr = urban.urbanParams.viewFactor.WallFrmRoad;
  auto &vf_rw = urban.urbanParams.viewFactor.RoadFrmWall;
  auto &vf_ww = urban.urbanParams.viewFactor.OtherWallFrmWall;
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

  // Access upward longwave radiation fields (to be updated)
  auto &upLwRoof = urban.roof.UpwardLongRad;
  auto &upLwSunlitWall = urban.sunlitWall.UpwardLongRad;
  auto &upLwShadedWall = urban.shadedWall.UpwardLongRad;
  auto &upLwImpRoad = urban.imperviousRoad.UpwardLongRad;
  auto &upLwPerRoad = urban.perviousRoad.UpwardLongRad;

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
        Real LtotForRoad = forcLRad(l) * vf_sr(l);

        // Impervious road (weight = 1 - fraction of pervious road)
        const Real fracImpRoad = 1.0 - fracPervRoad(l);

        auto fluxImpRoad = ComputeAbsRefEmiRadiation(
            emissImpRoad(l), tempImpRoad(l), LtotForRoad, fracImpRoad);

        auto refImpRoad = ComputeReflectedRadFromRoad(
            LtotForRoad, emissImpRoad(l), vf_sr(l), vf_wr(l), fracImpRoad);

        auto emiImpRoad = ComputeEmittedRadFromRoad(
            emissImpRoad(l), tempImpRoad(l), vf_sr(l), vf_wr(l), fracImpRoad);

        // Pervious road
        auto fluxPerRoad = ComputeAbsRefEmiRadiation(
            emissPerRoad(l), tempPerRoad(l), LtotForRoad, fracPervRoad(l));

        auto refPerRoad = ComputeReflectedRadFromRoad(
            LtotForRoad, emissPerRoad(l), vf_sr(l), vf_wr(l), fracPervRoad(l));

        auto emiPerRoad =
            ComputeEmittedRadFromRoad(emissPerRoad(l), tempPerRoad(l), vf_sr(l),
                                      vf_wr(l), fracPervRoad(l));

        // Combinging data from the roads
        Real RoadAbs =
            fluxImpRoad.absorbedWeighted + fluxPerRoad.absorbedWeighted;
        Real RoadRef =
            fluxImpRoad.reflectedWeighted + fluxPerRoad.reflectedWeighted;
        Real RoadEmi =
            fluxImpRoad.emittedWeighted + fluxPerRoad.emittedWeighted;

        // Reflected radiation from roads to sky, sunwall, and shadewall
        Real RoadRefToSky = refImpRoad.toSkyByWt + refPerRoad.toSkyByWt;
        Real RoadRefToSunwall =
            refImpRoad.toSunwallByWt + refPerRoad.toSunwallByWt;
        Real RoadRefToShadewall =
            refImpRoad.toShadewallByWt + refPerRoad.toShadewallByWt;

        // Emitted radiation from roads to sky, sunwall, and shadewall
        Real RoadEmiToSky = emiImpRoad.toSkyByWt + emiPerRoad.toSkyByWt;
        Real RoadEmiToSunwall =
            emiImpRoad.toSunwallByWt + emiPerRoad.toSunwallByWt;
        Real RoadEmiToShadewall =
            emiImpRoad.toShadewallByWt + emiPerRoad.toShadewallByWt;

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Computations for walls
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        // Total longwave downwelling to wall
        Real LtotForWall = forcLRad(l) * vf_sw(l);

        // Sunlit wall
        auto fluxSunlitWall = ComputeAbsRefEmiRadiation(
            emissWall(l), tempSunlitWall(l), LtotForWall, 1.0);

        auto refSunlitWall = ComputeReflectedRadFromWall(
            LtotForWall, emissWall(l), vf_sw(l), vf_rw(l), vf_ww(l));

        auto emiSunlitWall = ComputeEmittedRadFromWall(
            emissWall(l), tempSunlitWall(l), vf_sw(l), vf_rw(l), vf_ww(l));

        // Shaded wall
        auto fluxShadedWall = ComputeAbsRefEmiRadiation(
            emissWall(l), tempShadedWall(l), LtotForWall, 1.0);

        auto refShadedWall = ComputeReflectedRadFromWall(
            LtotForWall, emissWall(l), vf_sw(l), vf_rw(l), vf_ww(l));

        auto emiShadedWall = ComputeEmittedRadFromWall(
            emissWall(l), tempShadedWall(l), vf_sw(l), vf_rw(l), vf_ww(l));

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // initialize the net longwave radiation for each surface
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        netLwImpRoad(l) = fluxImpRoad.emitted - fluxImpRoad.absorbed;
        netLwPerRoad(l) = fluxPerRoad.emitted - fluxPerRoad.absorbed;
        netLwSunlitWall(l) = fluxSunlitWall.emitted - fluxSunlitWall.absorbed;
        netLwShadedWall(l) = fluxShadedWall.emitted - fluxShadedWall.absorbed;

        upLwImpRoad(l) = refImpRoad.toSky + emiImpRoad.toSky;
        upLwPerRoad(l) = refPerRoad.toSky + emiPerRoad.toSky;
        upLwSunlitWall(l) = refSunlitWall.toSky + emiSunlitWall.toSky;
        upLwShadedWall(l) = refShadedWall.toSky + emiShadedWall.toSky;

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Iteration loop for multiple reflections between surfaces
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        const int max_iter = 50;
        for (int iter = 0; iter < max_iter; iter++) {

          // step(1): Compute incoming radiation from wall-to-road and
          // wall-to-wall reflections

          // For roads: incoming from walls
          LtotForRoad = (refSunlitWall.toRoad + emiSunlitWall.toRoad +
                         refShadedWall.toRoad + emiShadedWall.toRoad) *
                        hwr(l);

          fluxImpRoad = ComputeAbsRefEmiRadiation(
              emissImpRoad(l), tempImpRoad(l), LtotForRoad, fracImpRoad);
          fluxPerRoad = ComputeAbsRefEmiRadiation(
              emissPerRoad(l), tempPerRoad(l), LtotForRoad, fracPervRoad(l));

          RoadAbs = fluxImpRoad.absorbedWeighted + fluxPerRoad.absorbedWeighted;
          RoadRef =
              fluxImpRoad.reflectedWeighted + fluxPerRoad.reflectedWeighted;

          // For sunlit wall: incoming from roads and shaded wall
          LtotForWall = (RoadRefToSunwall + RoadEmiToSunwall) / hwr(l) +
                        refShadedWall.toOtherwall + emiShadedWall.toOtherwall;
          fluxSunlitWall = ComputeAbsRefEmiRadiation(
              emissWall(l), tempSunlitWall(l), LtotForWall, 1.0);

          // For shaded wall: incoming from roads and sunlit wall
          LtotForWall = (RoadRefToShadewall + RoadEmiToShadewall) / hwr(l) +
                        refSunlitWall.toOtherwall + emiSunlitWall.toOtherwall;
          fluxShadedWall = ComputeAbsRefEmiRadiation(
              emissWall(l), tempShadedWall(l), LtotForWall, 1.0);

          // Set emitted values to zero so they are not counted multiple times
          emiSunlitWall.toRoad = 0.0;
          emiSunlitWall.toOtherwall = 0.0;
          emiShadedWall.toRoad = 0.0;
          emiShadedWall.toOtherwall = 0.0;
          RoadEmiToSunwall = 0.0;
          RoadEmiToShadewall = 0.0;

          // step(2): Update net longwave by subtracting newly absorbed
          // radiation
          netLwImpRoad(l) -= fluxImpRoad.absorbed;
          netLwPerRoad(l) -= fluxPerRoad.absorbed;
          netLwSunlitWall(l) -= fluxSunlitWall.absorbed;
          netLwShadedWall(l) -= fluxShadedWall.absorbed;

          // step(3): Compute reflected radiation components for this iteration
          refImpRoad = ComputeReflectedRadFromRoad(
              LtotForRoad, emissImpRoad(l), vf_sr(l), vf_wr(l), fracImpRoad);
          refPerRoad =
              ComputeReflectedRadFromRoad(LtotForRoad, emissPerRoad(l),
                                          vf_sr(l), vf_wr(l), fracPervRoad(l));

          RoadRefToSky = refImpRoad.toSkyByWt + refPerRoad.toSkyByWt;
          RoadRefToSunwall =
              refImpRoad.toSunwallByWt + refPerRoad.toSunwallByWt;
          RoadRefToShadewall =
              refImpRoad.toShadewallByWt + refPerRoad.toShadewallByWt;

          refSunlitWall = ComputeReflectedRadFromWall(
              LtotForWall, emissWall(l), vf_sw(l), vf_rw(l), vf_ww(l));
          refShadedWall = ComputeReflectedRadFromWall(
              LtotForWall, emissWall(l), vf_sw(l), vf_rw(l), vf_ww(l));

          // step(4): Update upward longwave radiation
          upLwImpRoad(l) += refImpRoad.toSky;
          upLwPerRoad(l) += refPerRoad.toSky;
          upLwSunlitWall(l) += refSunlitWall.toSky;
          upLwShadedWall(l) += refShadedWall.toSky;

          // Check convergence
          Real crit =
              Kokkos::max(RoadAbs, Kokkos::max(fluxSunlitWall.absorbed,
                                               fluxShadedWall.absorbed));
          const Real errcrit = 0.001;
          if (crit < errcrit)
            break;
        }
      });

  Kokkos::fence();
}

} // namespace URBANXX
