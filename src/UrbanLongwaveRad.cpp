#include "private/AtmosphereTypeImpl.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include "private/UrbanDebugUtils.h"
#include "private/UrbanLongwaveRadImpl.h"
#include "private/UrbanParamsTypeImpl.h"
#include "private/UrbanSurfaceTypeImpl.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>
#include <iostream>

namespace URBANXX {

// Convergence parameters for longwave radiation iteration
constexpr int LONGWAVE_MAX_ITERATIONS = 50;
constexpr Real LONGWAVE_CONVERGENCE_THRESHOLD = 0.001;

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
  Real toSky;            // reflected radiation to sky
  Real toSunlitWall;     // reflected radiation to sunlit wall
  Real toShadedWall;     // reflected radiation to shaded wall
  Real toSkyByWt;        // reflected to sky * weight
  Real toSunlitWallByWt; // reflected to sunlit wall * weight
  Real toShadedWallByWt; // reflected to shaded wall * weight
};

// Structure to hold reflected radiation components from wall
struct ReflectedRadFromWall {
  Real toSky;           // reflected radiation to sky
  Real toRoad;          // reflected radiation to road
  Real toOtherWall;     // reflected radiation to other wall
  Real toSkyByWt;       // reflected to sky * weight
  Real toRoadByWt;      // reflected to road * weight
  Real toOtherWallByWt; // reflected to other wall * weight
};

// Structure to hold all radiation components for a road surface
struct RoadRadiation {
  SurfaceLongwaveFluxes flux;
  ReflectedRadFromRoad ref;
  ReflectedRadFromRoad emi;
};

// Structure to hold all radiation components for a wall surface
struct WallRadiation {
  SurfaceLongwaveFluxes flux;
  ReflectedRadFromWall ref;
  ReflectedRadFromWall emi;
};

// Structure to hold view factors for road surfaces
struct RoadViewFactors {
  Real sky;  // view factor from road to sky
  Real wall; // view factor from road to wall
};

// Structure to hold view factors for wall surfaces
struct WallViewFactors {
  Real sky;  // view factor from wall to sky
  Real road; // view factor from wall to road
  Real wall; // view factor from wall to other wall
};

// Helper function to compute longwave radiation components for a surface
KOKKOS_INLINE_FUNCTION
SurfaceLongwaveFluxes Fluxes(const Real emissivity, const Real temperature,
                             const Real LtotForSurface, const Real weight) {
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
                                                 const RoadViewFactors &vf,
                                                 const Real weight) {
  ReflectedRadFromRoad ref;

  // Distribute radiation by view factors
  ref.toSky = radiation * vf.sky;
  ref.toSunlitWall = radiation * vf.wall;
  ref.toShadedWall = radiation * vf.wall;

  // Apply weight (fraction of total road)
  ref.toSkyByWt = ref.toSky * weight;
  ref.toSunlitWallByWt = ref.toSunlitWall * weight;
  ref.toShadedWallByWt = ref.toShadedWall * weight;

  return ref;
}

// Helper function to compute reflected radiation components from road
KOKKOS_INLINE_FUNCTION
ReflectedRadFromRoad ReflectRoad(const Real incomingRad, const Real emissivity,
                                 const RoadViewFactors &vf, const Real weight) {
  // Compute reflected radiation
  const Real reflectedRad = (1.0 - emissivity) * incomingRad;

  // Distribute using common function
  return DistributeRadiationFromRoad(reflectedRad, vf, weight);
}

// Helper function to compute emitted radiation components from road
KOKKOS_INLINE_FUNCTION
ReflectedRadFromRoad EmitRoad(const Real emissivity, const Real temperature,
                              const RoadViewFactors &vf, const Real weight) {
  // Compute emitted radiation
  const Real emittedRad = emissivity * STEBOL * Kokkos::pow(temperature, 4.0);

  // Distribute using common function
  return DistributeRadiationFromRoad(emittedRad, vf, weight);
}

// Helper function to distribute radiation from wall to other surfaces
KOKKOS_INLINE_FUNCTION
ReflectedRadFromWall DistributeRadiationFromWall(const Real radiation,
                                                 const WallViewFactors &vf) {
  ReflectedRadFromWall ref;

  // Distribute radiation by view factors
  ref.toSky = radiation * vf.sky;
  ref.toRoad = radiation * vf.road;
  ref.toOtherWall = radiation * vf.wall;

  // For walls, weight is always 1.0 (no fractional surfaces)
  ref.toSkyByWt = ref.toSky;
  ref.toRoadByWt = ref.toRoad;
  ref.toOtherWallByWt = ref.toOtherWall;

  return ref;
}

// Helper function to compute reflected radiation components from wall
KOKKOS_INLINE_FUNCTION
ReflectedRadFromWall ReflectWall(const Real incomingRad, const Real emissivity,
                                 const WallViewFactors &vf) {
  // Compute reflected radiation
  const Real reflectedRad = (1.0 - emissivity) * incomingRad;

  // Distribute using common function
  return DistributeRadiationFromWall(reflectedRad, vf);
}

// Helper function to compute emitted radiation components from wall
KOKKOS_INLINE_FUNCTION
ReflectedRadFromWall EmitWall(const Real emissivity, const Real temperature,
                              const WallViewFactors &vf) {
  // Compute emitted radiation
  const Real emittedRad = emissivity * STEBOL * Kokkos::pow(temperature, 4.0);

  // Distribute using common function
  return DistributeRadiationFromWall(emittedRad, vf);
}

// Helper function to initialize a single road surface (impervious or pervious)
KOKKOS_INLINE_FUNCTION
RoadRadiation InitializeSingleRoad(const Real LtotForRoad,
                                   const Real emissivity,
                                   const Real temperature,
                                   const RoadViewFactors &vf,
                                   const Real fraction) {

  RoadRadiation rad;
  rad.flux = Fluxes(emissivity, temperature, LtotForRoad, fraction);
  rad.ref = ReflectRoad(LtotForRoad, emissivity, vf, fraction);
  rad.emi = EmitRoad(emissivity, temperature, vf, fraction);
  return rad;
}

// Helper function to combine radiation from two road surfaces
KOKKOS_INLINE_FUNCTION
void CombineRoadRadiation(const RoadRadiation &impRoad,
                          const RoadRadiation &perRoad, Real &RoadAbs,
                          Real &RoadRef, Real &RoadEmi, Real &RoadRefToSky,
                          Real &RoadRefToSunlitWall, Real &RoadRefToShadedWall,
                          Real &RoadEmiToSky, Real &RoadEmiToSunlitWall,
                          Real &RoadEmiToShadedWall) {

  // Combine absorbed, reflected, and emitted
  RoadAbs = impRoad.flux.absorbedWeighted + perRoad.flux.absorbedWeighted;
  RoadRef = impRoad.flux.reflectedWeighted + perRoad.flux.reflectedWeighted;
  RoadEmi = impRoad.flux.emittedWeighted + perRoad.flux.emittedWeighted;

  // Combine reflected radiation to other surfaces
  RoadRefToSky = impRoad.ref.toSkyByWt + perRoad.ref.toSkyByWt;
  RoadRefToSunlitWall =
      impRoad.ref.toSunlitWallByWt + perRoad.ref.toSunlitWallByWt;
  RoadRefToShadedWall =
      impRoad.ref.toShadedWallByWt + perRoad.ref.toShadedWallByWt;

  // Combine emitted radiation to other surfaces
  RoadEmiToSky = impRoad.emi.toSkyByWt + perRoad.emi.toSkyByWt;
  RoadEmiToSunlitWall =
      impRoad.emi.toSunlitWallByWt + perRoad.emi.toSunlitWallByWt;
  RoadEmiToShadedWall =
      impRoad.emi.toShadedWallByWt + perRoad.emi.toShadedWallByWt;
}

// Helper function to initialize a single wall surface (sunlit or shaded)
KOKKOS_INLINE_FUNCTION
WallRadiation InitializeSingleWall(const Real LtotForWall,
                                   const Real emissivity,
                                   const Real temperature,
                                   const WallViewFactors &vf) {

  WallRadiation rad;
  rad.flux = Fluxes(emissivity, temperature, LtotForWall, 1.0);
  rad.ref = ReflectWall(LtotForWall, emissivity, vf);
  rad.emi = EmitWall(emissivity, temperature, vf);
  return rad;
}

void ComputeNetLongwave(URBANXX::_p_UrbanType &urban) {
  // Get number of landunits for parallel execution
  const int numLandunits = urban.numLandunits;

  // Access atmospheric forcing data
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
  auto &tempRoof = urban.roof.EffectiveSurfTemp;
  auto &tempSunlitWall = urban.sunlitWall.EffectiveSurfTemp;
  auto &tempShadedWall = urban.shadedWall.EffectiveSurfTemp;
  auto &tempImpRoad = urban.imperviousRoad.EffectiveSurfTemp;
  auto &tempPerRoad = urban.perviousRoad.EffectiveSurfTemp;

  // Access net longwave radiation fields (to be updated)
  auto &netLwRoof = urban.roof.NetLongRad;
  auto &netLwSunlitWall = urban.sunlitWall.NetLongRad;
  auto &netLwShadedWall = urban.shadedWall.NetLongRad;
  auto &netLwImpRoad = urban.imperviousRoad.NetLongRad;
  auto &netLwPerRoad = urban.perviousRoad.NetLongRad;
  auto &netLwRoof = urban.roof.NetLongRad;

  // Access upward longwave radiation fields (to be updated)
  auto &upLwRoof = urban.roof.UpwardLongRad;
  auto &upLwSunlitWall = urban.sunlitWall.UpwardLongRad;
  auto &upLwShadedWall = urban.shadedWall.UpwardLongRad;
  auto &upLwImpRoad = urban.imperviousRoad.UpwardLongRad;
  auto &upLwPerRoad = urban.perviousRoad.UpwardLongRad;
  auto &upLwRoof = urban.roof.UpwardLongRad;

  // View to track non-converged landunits
  Kokkos::View<int *> nonConvergedCount("nonConvergedCount", 1);

  Kokkos::parallel_for(
      "ComputeNetLongwave", numLandunits, KOKKOS_LAMBDA(const int l) {
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Computations for roof
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        upLwRoof(l) = emissRoof(l) * STEBOL * Kokkos::pow(tempRoof(l), 4.0) +
                      (1.0 - emissRoof(l)) * forcLRad(l);
        netLwRoof(l) = upLwRoof(l) - forcLRad(l);

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Computations for roads
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        // Total longwave downwelling to road
        Real LtotForRoad = forcLRad(l) * vf_sr(l);

        // Impervious road (weight = 1 - fraction of pervious road)
        const Real fracImpRoad = 1.0 - fracPervRoad(l);

        // View factors for roads
        const RoadViewFactors roadVF = {vf_sr(l), vf_wr(l)};

        // Initialize impervious and pervious roads
        auto impRoad = InitializeSingleRoad(
            LtotForRoad, emissImpRoad(l), tempImpRoad(l), roadVF, fracImpRoad);
        auto perRoad =
            InitializeSingleRoad(LtotForRoad, emissPerRoad(l), tempPerRoad(l),
                                 roadVF, fracPervRoad(l));

        // Combine both roads
        Real RoadAbs, RoadRef, RoadEmi;
        Real RoadRefToSky, RoadRefToSunlitWall, RoadRefToShadedWall;
        Real RoadEmiToSky, RoadEmiToSunlitWall, RoadEmiToShadedWall;
        CombineRoadRadiation(impRoad, perRoad, RoadAbs, RoadRef, RoadEmi,
                             RoadRefToSky, RoadRefToSunlitWall,
                             RoadRefToShadedWall, RoadEmiToSky,
                             RoadEmiToSunlitWall, RoadEmiToShadedWall);

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Computations for walls
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        // Total longwave downwelling to wall
        Real LtotForWall = forcLRad(l) * vf_sw(l);

        // View factors for walls
        const WallViewFactors wallVF = {vf_sw(l), vf_rw(l), vf_ww(l)};

        // Initialize sunlit and shaded walls
        auto sunlitWall = InitializeSingleWall(LtotForWall, emissWall(l),
                                               tempSunlitWall(l), wallVF);
        auto shadedWall = InitializeSingleWall(LtotForWall, emissWall(l),
                                               tempShadedWall(l), wallVF);

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // initialize the net longwave radiation for each surface
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        netLwImpRoad(l) = impRoad.flux.emitted - impRoad.flux.absorbed;
        netLwPerRoad(l) = perRoad.flux.emitted - perRoad.flux.absorbed;
        netLwSunlitWall(l) = sunlitWall.flux.emitted - sunlitWall.flux.absorbed;
        netLwShadedWall(l) = shadedWall.flux.emitted - shadedWall.flux.absorbed;

        upLwImpRoad(l) = impRoad.ref.toSky + impRoad.emi.toSky;
        upLwPerRoad(l) = perRoad.ref.toSky + perRoad.emi.toSky;
        upLwSunlitWall(l) = sunlitWall.ref.toSky + sunlitWall.emi.toSky;
        upLwShadedWall(l) = shadedWall.ref.toSky + shadedWall.emi.toSky;

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Iteration loop for multiple reflections between surfaces
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        bool converged = false;
        for (int iter = 0; iter < LONGWAVE_MAX_ITERATIONS; iter++) {

          // step(1): Compute incoming radiation from wall-to-road and
          // wall-to-wall reflections

          // For roads: incoming from walls
          LtotForRoad = (sunlitWall.ref.toRoad + sunlitWall.emi.toRoad +
                         shadedWall.ref.toRoad + shadedWall.emi.toRoad) *
                        hwr(l);

          impRoad.flux =
              Fluxes(emissImpRoad(l), tempImpRoad(l), LtotForRoad, fracImpRoad);
          perRoad.flux = Fluxes(emissPerRoad(l), tempPerRoad(l), LtotForRoad,
                                fracPervRoad(l));

          RoadAbs =
              impRoad.flux.absorbedWeighted + perRoad.flux.absorbedWeighted;
          RoadRef =
              impRoad.flux.reflectedWeighted + perRoad.flux.reflectedWeighted;

          // For sunlit wall: incoming from roads and shaded wall
          Real LtotForSunWall =
              (RoadRefToSunlitWall + RoadEmiToSunlitWall) / hwr(l) +
              shadedWall.ref.toOtherWall + shadedWall.emi.toOtherWall;
          sunlitWall.flux =
              Fluxes(emissWall(l), tempSunlitWall(l), LtotForSunWall, 1.0);

          // For shaded wall: incoming from roads and sunlit wall
          Real LtotForShadeWall =
              (RoadRefToShadedWall + RoadEmiToShadedWall) / hwr(l) +
              sunlitWall.ref.toOtherWall + sunlitWall.emi.toOtherWall;
          shadedWall.flux =
              Fluxes(emissWall(l), tempShadedWall(l), LtotForShadeWall, 1.0);

          // Set emitted values to zero so they are not counted multiple times
          sunlitWall.emi.toRoad = 0.0;
          sunlitWall.emi.toOtherWall = 0.0;
          shadedWall.emi.toRoad = 0.0;
          shadedWall.emi.toOtherWall = 0.0;
          RoadEmiToSunlitWall = 0.0;
          RoadEmiToShadedWall = 0.0;

          // step(2): Update net longwave by subtracting newly absorbed
          // radiation
          netLwImpRoad(l) -= impRoad.flux.absorbed;
          netLwPerRoad(l) -= perRoad.flux.absorbed;
          netLwSunlitWall(l) -= sunlitWall.flux.absorbed;
          netLwShadedWall(l) -= shadedWall.flux.absorbed;

          // step(3): Compute reflected radiation components for this iteration
          impRoad.ref =
              ReflectRoad(LtotForRoad, emissImpRoad(l), roadVF, fracImpRoad);
          perRoad.ref = ReflectRoad(LtotForRoad, emissPerRoad(l), roadVF,
                                    fracPervRoad(l));

          RoadRefToSky = impRoad.ref.toSkyByWt + perRoad.ref.toSkyByWt;
          RoadRefToSunlitWall =
              impRoad.ref.toSunlitWallByWt + perRoad.ref.toSunlitWallByWt;
          RoadRefToShadedWall =
              impRoad.ref.toShadedWallByWt + perRoad.ref.toShadedWallByWt;

          sunlitWall.ref = ReflectWall(LtotForSunWall, emissWall(l), wallVF);
          shadedWall.ref = ReflectWall(LtotForShadeWall, emissWall(l), wallVF);

          // step(4): Update upward longwave radiation
          upLwImpRoad(l) += impRoad.ref.toSky;
          upLwPerRoad(l) += perRoad.ref.toSky;
          upLwSunlitWall(l) += sunlitWall.ref.toSky;
          upLwShadedWall(l) += shadedWall.ref.toSky;

          // Check convergence
          Real convergence_criteria =
              Kokkos::max(RoadAbs, Kokkos::max(sunlitWall.flux.absorbed,
                                               shadedWall.flux.absorbed));
          if (convergence_criteria < LONGWAVE_CONVERGENCE_THRESHOLD) {
            converged = true;
            break;
          }
        }
        if (!converged) {
          Kokkos::atomic_fetch_add(&nonConvergedCount(0), 1);
        }
      });

  Kokkos::fence();

  // Check for non-convergence on host
  auto nonConvergedCount_h = Kokkos::create_mirror_view(nonConvergedCount);
  Kokkos::deep_copy(nonConvergedCount_h, nonConvergedCount);

  if (nonConvergedCount_h(0) > 0) {
    std::cerr << "WARNING: Longwave radiation iteration did not converge for "
              << nonConvergedCount_h(0) << " landunit(s)" << std::endl;
    // TODO: Consider adding error handling mechanism here
  }

  if (0) {
    print_view_1d(urban.imperviousRoad.NetLongRad, "imperviousRoad.NetLongRad");
    print_view_1d(urban.perviousRoad.NetLongRad, "perviousRoad.NetLongRad");
    print_view_1d(urban.sunlitWall.NetLongRad, "sunlitWall.NetLongRad");
    print_view_1d(urban.shadedWall.NetLongRad, "shadedWall.NetLongRad");
  }
}

} // namespace URBANXX
