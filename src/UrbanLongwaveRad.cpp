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
constexpr Real LONGWAVE_CONVERGENCE_THRESHOLD = 0.0001;

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

  // Get raw pointers from Views to avoid CUDA extended lambda issues
  Real *forcLRadPtr = urban.atmosphereData.ForcLRad.data();

  // View factors and canyon parameters
  Real *vf_srPtr = urban.urbanParams.viewFactor.SkyFrmRoad.data();
  Real *vf_swPtr = urban.urbanParams.viewFactor.SkyFrmWall.data();
  Real *vf_wrPtr = urban.urbanParams.viewFactor.WallFrmRoad.data();
  Real *vf_rwPtr = urban.urbanParams.viewFactor.RoadFrmWall.data();
  Real *vf_wwPtr = urban.urbanParams.viewFactor.OtherWallFrmWall.data();
  Real *hwrPtr = urban.urbanParams.CanyonHwr.data();
  Real *fracPervRoadPtr = urban.urbanParams.FracPervRoadOfTotalRoad.data();

  // Emissivities
  Real *emissRoofPtr = urban.urbanParams.emissivity.Roof.data();
  Real *emissWallPtr = urban.urbanParams.emissivity.Wall.data();
  Real *emissImpRoadPtr = urban.urbanParams.emissivity.ImperviousRoad.data();
  Real *emissPerRoadPtr = urban.urbanParams.emissivity.PerviousRoad.data();

  // Surface temperatures
  Real *tempRoofPtr = urban.roof.Temperature.data();
  Real *tempSunlitWallPtr = urban.sunlitWall.Temperature.data();
  Real *tempShadedWallPtr = urban.shadedWall.Temperature.data();
  Real *tempImpRoadPtr = urban.imperviousRoad.Temperature.data();
  Real *tempPerRoadPtr = urban.perviousRoad.Temperature.data();

  // Net longwave radiation fields
  Real *netLwSunlitWallPtr = urban.sunlitWall.NetLongRad.data();
  Real *netLwShadedWallPtr = urban.shadedWall.NetLongRad.data();
  Real *netLwImpRoadPtr = urban.imperviousRoad.NetLongRad.data();
  Real *netLwPerRoadPtr = urban.perviousRoad.NetLongRad.data();

  // Upward longwave radiation fields
  Real *upLwSunlitWallPtr = urban.sunlitWall.UpwardLongRad.data();
  Real *upLwShadedWallPtr = urban.shadedWall.UpwardLongRad.data();
  Real *upLwImpRoadPtr = urban.imperviousRoad.UpwardLongRad.data();
  Real *upLwPerRoadPtr = urban.perviousRoad.UpwardLongRad.data();

  // View to track non-converged landunits
  Kokkos::View<int *> nonConvergedCount("nonConvergedCount", 1);

  using ExecSpace = Kokkos::DefaultExecutionSpace;
  Kokkos::parallel_for(
      "ComputeNetLongwave", Kokkos::RangePolicy<ExecSpace>(0, numLandunits),
      KOKKOS_LAMBDA(const int l) {
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Computations for roads
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        // Total longwave downwelling to road
        Real LtotForRoad = forcLRadPtr[l] * vf_srPtr[l];

        // Impervious road (weight = 1 - fraction of pervious road)
        const Real fracImpRoad = 1.0 - fracPervRoadPtr[l];

        // View factors for roads
        const RoadViewFactors roadVF = {vf_srPtr[l], vf_wrPtr[l]};

        // Initialize impervious and pervious roads
        auto impRoad =
            InitializeSingleRoad(LtotForRoad, emissImpRoadPtr[l],
                                 tempImpRoadPtr[l], roadVF, fracImpRoad);
        auto perRoad =
            InitializeSingleRoad(LtotForRoad, emissPerRoadPtr[l],
                                 tempPerRoadPtr[l], roadVF, fracPervRoadPtr[l]);

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
        Real LtotForWall = forcLRadPtr[l] * vf_swPtr[l];

        // View factors for walls
        const WallViewFactors wallVF = {vf_swPtr[l], vf_rwPtr[l], vf_wwPtr[l]};

        // Initialize sunlit and shaded walls
        auto sunlitWall = InitializeSingleWall(LtotForWall, emissWallPtr[l],
                                               tempSunlitWallPtr[l], wallVF);
        auto shadedWall = InitializeSingleWall(LtotForWall, emissWallPtr[l],
                                               tempShadedWallPtr[l], wallVF);

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // initialize the net longwave radiation for each surface
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        netLwImpRoadPtr[l] = impRoad.flux.emitted - impRoad.flux.absorbed;
        netLwPerRoadPtr[l] = perRoad.flux.emitted - perRoad.flux.absorbed;
        netLwSunlitWallPtr[l] =
            sunlitWall.flux.emitted - sunlitWall.flux.absorbed;
        netLwShadedWallPtr[l] =
            shadedWall.flux.emitted - shadedWall.flux.absorbed;

        upLwImpRoadPtr[l] = impRoad.ref.toSky + impRoad.emi.toSky;
        upLwPerRoadPtr[l] = perRoad.ref.toSky + perRoad.emi.toSky;
        upLwSunlitWallPtr[l] = sunlitWall.ref.toSky + sunlitWall.emi.toSky;
        upLwShadedWallPtr[l] = shadedWall.ref.toSky + shadedWall.emi.toSky;

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
                        hwrPtr[l];

          impRoad.flux = Fluxes(emissImpRoadPtr[l], tempImpRoadPtr[l],
                                LtotForRoad, fracImpRoad);
          perRoad.flux = Fluxes(emissPerRoadPtr[l], tempPerRoadPtr[l],
                                LtotForRoad, fracPervRoadPtr[l]);

          RoadAbs =
              impRoad.flux.absorbedWeighted + perRoad.flux.absorbedWeighted;
          RoadRef =
              impRoad.flux.reflectedWeighted + perRoad.flux.reflectedWeighted;

          // For sunlit wall: incoming from roads and shaded wall
          Real LtotForSunWall =
              (RoadRefToSunlitWall + RoadEmiToSunlitWall) / hwrPtr[l] +
              shadedWall.ref.toOtherWall + shadedWall.emi.toOtherWall;
          sunlitWall.flux = Fluxes(emissWallPtr[l], tempSunlitWallPtr[l],
                                   LtotForSunWall, 1.0);

          // For shaded wall: incoming from roads and sunlit wall
          Real LtotForShadeWall =
              (RoadRefToShadedWall + RoadEmiToShadedWall) / hwrPtr[l] +
              sunlitWall.ref.toOtherWall + sunlitWall.emi.toOtherWall;
          shadedWall.flux = Fluxes(emissWallPtr[l], tempShadedWallPtr[l],
                                   LtotForShadeWall, 1.0);

          // Set emitted values to zero so they are not counted multiple times
          sunlitWall.emi.toRoad = 0.0;
          sunlitWall.emi.toOtherWall = 0.0;
          shadedWall.emi.toRoad = 0.0;
          shadedWall.emi.toOtherWall = 0.0;
          RoadEmiToSunlitWall = 0.0;
          RoadEmiToShadedWall = 0.0;

          // step(2): Update net longwave by subtracting newly absorbed
          // radiation
          netLwImpRoadPtr[l] -= impRoad.flux.absorbed;
          netLwPerRoadPtr[l] -= perRoad.flux.absorbed;
          netLwSunlitWallPtr[l] -= sunlitWall.flux.absorbed;
          netLwShadedWallPtr[l] -= shadedWall.flux.absorbed;

          // step(3): Compute reflected radiation components for this iteration
          impRoad.ref =
              ReflectRoad(LtotForRoad, emissImpRoadPtr[l], roadVF, fracImpRoad);
          perRoad.ref = ReflectRoad(LtotForRoad, emissPerRoadPtr[l], roadVF,
                                    fracPervRoadPtr[l]);

          RoadRefToSky = impRoad.ref.toSkyByWt + perRoad.ref.toSkyByWt;
          RoadRefToSunlitWall =
              impRoad.ref.toSunlitWallByWt + perRoad.ref.toSunlitWallByWt;
          RoadRefToShadedWall =
              impRoad.ref.toShadedWallByWt + perRoad.ref.toShadedWallByWt;

          sunlitWall.ref = ReflectWall(LtotForSunWall, emissWallPtr[l], wallVF);
          shadedWall.ref =
              ReflectWall(LtotForShadeWall, emissWallPtr[l], wallVF);

          // step(4): Update upward longwave radiation
          upLwImpRoadPtr[l] += impRoad.ref.toSky;
          upLwPerRoadPtr[l] += perRoad.ref.toSky;
          upLwSunlitWallPtr[l] += sunlitWall.ref.toSky;
          upLwShadedWallPtr[l] += shadedWall.ref.toSky;

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