#include "private/AtmosphereTypeImpl.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include "private/UrbanDebugUtils.h"
#include "private/UrbanLongwaveRadImpl.h"
#include "private/UrbanParamsTypeImpl.h"
#include "private/UrbanSurfaceTypeImpl.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>
#include <cmath>
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

// Update internal snow state (H2OSno, IntSno, FracSno) for all snow-covered
// horizontal urban surfaces (roof, impervious road, pervious road) using only
// atmospheric forcing already stored inside the urban struct.  Walls are
// excluded — vertical surfaces do not accumulate snow.
void ComputeSnowCover(URBANXX::_p_UrbanType &urban) {
  const double dtime = 1800.0; // s — hardcoded URBANxx timestep
  // ELM's accum_factor is declared real(r8) but assigned literal 0.1 (no _r8
  // suffix). Without -fdefault-real-8, gfortran widens the nearest float32 to
  // float64: 0.100000001490116119..., NOT the true double
  // 0.1000000000000000056. Using (double)(0.1f) matches ELM's CanopyHydrology
  // exactly.
  const double accum_factor = (double)(0.1f);
  const double n_melt = 1.0;          // flat urban surfaces (roof/road)
  const double large_intsnow = 1.0e8; // ELM cap on int_snow

  auto forc_snow = urban.atmosphereData.ForcSnow; // kg/m²/s
  auto forc_temp = urban.atmosphereData.ForcTemp; // air temperature (K)

  auto roofH2OSno = urban.roof.H2OSno;
  auto roofIntSno = urban.roof.IntSno;
  auto roofFracSno = urban.roof.FracSno;
  auto roofSubSnow = urban.roof.QflxSubSnow;
  auto roofSnowDepth = urban.roof.SnowDepth;

  auto impH2OSno = urban.imperviousRoad.H2OSno;
  auto impIntSno = urban.imperviousRoad.IntSno;
  auto impFracSno = urban.imperviousRoad.FracSno;
  auto impSubSnow = urban.imperviousRoad.QflxSubSnow;
  auto impSnowDepth = urban.imperviousRoad.SnowDepth;

  auto perH2OSno = urban.perviousRoad.H2OSno;
  auto perIntSno = urban.perviousRoad.IntSno;
  auto perFracSno = urban.perviousRoad.FracSno;
  auto perSubSnow = urban.perviousRoad.QflxSubSnow;
  auto perSnowDepth = urban.perviousRoad.SnowDepth;

  Kokkos::parallel_for(
      "ComputeSnowCover", urban.numLandunits, KOKKOS_LAMBDA(const int l) {
        // Compute bulk density of new snow from air temperature
        // (port of ELM CanopyHydrologyMod.F90, urban path, oldfflag=0).
        // tfrz = 273.15 K
        constexpr double tfrz = 273.15;
        const double T = forc_temp(l);
        double bifall;
        if (T > tfrz + 2.0)
          bifall = 50.0 + 1.7 * std::pow(17.0, 1.5);
        else if (T > tfrz - 15.0)
          bifall = 50.0 + 1.7 * std::pow(T - tfrz + 15.0, 1.5);
        else
          bifall = 50.0;

        // Lambda: update one snow-covered surface given its views.
        // forc_snow and bifall are shared across all three surfaces.
        auto update = [&](Kokkos::View<double *> h2osno_v,
                          Kokkos::View<double *> int_sno_v,
                          Kokkos::View<double *> frac_sno_v,
                          Kokkos::View<double *> sub_snow_v,
                          Kokkos::View<double *> snow_depth_v) {
          const double newsnow = forc_snow(l) * dtime; // kg/m²
          const double sub_loss =
              sub_snow_v(l) * dtime; // previous-step sublimation

          double h2osno = h2osno_v(l);
          double int_sno = int_sno_v(l);
          double frac_sno = frac_sno_v(l);

          // Save SWE before sublimation removal for proportional depth scaling
          const double h2osno_old = h2osno;

          // Remove sublimation from existing pack
          h2osno = Kokkos::max(0.0, h2osno - sub_loss);

          if (h2osno <= 0.0) {
            // No existing snowpack
            if (newsnow > 0.0) {
              frac_sno = tanh(accum_factor * newsnow);
              double frac_safe = Kokkos::max(frac_sno, 1.0e-6);
              double denom =
                  0.5 *
                  (cos(M_PI * Kokkos::pow(1.0 - frac_safe, 1.0 / n_melt)) +
                   1.0);
              double temp_intsnow = (denom > 0.0) ? (newsnow / denom) : newsnow;
              int_sno = Kokkos::min(large_intsnow, temp_intsnow) + newsnow;
              h2osno = newsnow;
            } else {
              frac_sno = 0.0;
              int_sno = 0.0;
            }
          } else {
            // Existing snowpack present
            if (newsnow > 0.0) {
              // Accumulation branch
              frac_sno =
                  1.0 - (1.0 - tanh(accum_factor * newsnow)) * (1.0 - frac_sno);
              double frac_safe = Kokkos::max(frac_sno, 1.0e-6);
              double denom =
                  0.5 *
                  (cos(M_PI * Kokkos::pow(1.0 - frac_safe, 1.0 / n_melt)) +
                   1.0);
              double temp_intsnow = (denom > 0.0) ? ((h2osno + newsnow) / denom)
                                                  : (h2osno + newsnow);
              int_sno = Kokkos::min(large_intsnow, temp_intsnow) + newsnow;
              h2osno += newsnow;
            }
            // Depletion branch (Niu-Yang 2007): recompute frac_sno from
            // ratio of current h2osno to the max-season int_sno
            if (int_sno > 0.0) {
              double smr = Kokkos::min(1.0, h2osno / int_sno);
              frac_sno =
                  1.0 -
                  Kokkos::pow(acos(Kokkos::min(1.0, 2.0 * smr - 1.0)) / M_PI,
                              n_melt);
            }
          }

          // --- Update geometric snow depth ---
          // Scale existing depth by the fraction of SWE remaining after
          // sublimation, then add new-snow depth.
          double snow_depth = snow_depth_v(l);
          if (h2osno_old > 0.0)
            snow_depth *=
                (Kokkos::max(0.0, h2osno_old - sub_loss) / h2osno_old);
          if (snow_depth < 0.0 || h2osno_old <= 0.0)
            snow_depth = 0.0;
          snow_depth += newsnow / bifall;
          if (h2osno <= 0.0)
            snow_depth = 0.0;
          snow_depth_v(l) = snow_depth;

          h2osno_v(l) = h2osno;
          int_sno_v(l) = int_sno;
          frac_sno_v(l) = Kokkos::max(0.0, Kokkos::min(1.0, frac_sno));
        };

        update(roofH2OSno, roofIntSno, roofFracSno, roofSubSnow, roofSnowDepth);
        update(impH2OSno, impIntSno, impFracSno, impSubSnow, impSnowDepth);
        update(perH2OSno, perIntSno, perFracSno, perSubSnow, perSnowDepth);
      });
  Kokkos::fence();
}

void ComputeNetLongwave(URBANXX::_p_UrbanType &urban) {
  // Get number of landunits for parallel execution
  const int numLandunits = urban.numLandunits;

  // Access atmospheric forcing data
  auto forcLRad = urban.atmosphereData.ForcLRad;

  // Access urban parameters - view factors and canyon height-to-width ratio
  auto vf_sr = urban.urbanParams.viewFactor.SkyFrmRoad;
  auto vf_sw = urban.urbanParams.viewFactor.SkyFrmWall;
  auto vf_wr = urban.urbanParams.viewFactor.WallFrmRoad;
  auto vf_rw = urban.urbanParams.viewFactor.RoadFrmWall;
  auto vf_ww = urban.urbanParams.viewFactor.OtherWallFrmWall;
  auto hwr = urban.urbanParams.CanyonHwr;
  auto fracPervRoad = urban.urbanParams.FracPervRoadOfTotalRoad;

  // Access urban parameters - emissivities
  auto emissRoof = urban.urbanParams.emissivity.Roof;
  auto emissWall = urban.urbanParams.emissivity.Wall;
  auto emissImpRoad = urban.urbanParams.emissivity.ImperviousRoad;
  auto emissPerRoad = urban.urbanParams.emissivity.PerviousRoad;

  // Snow-covered fraction views for emissivity blending (walls excluded)
  auto fracSnoRoof = urban.roof.FracSno;
  auto fracSnoImpRoad = urban.imperviousRoad.FracSno;
  auto fracSnoPerRoad = urban.perviousRoad.FracSno;
  constexpr double snoem = 0.97; // ELM snow emissivity constant

  // Access surface temperatures
  auto tempRoof = urban.roof.EffectiveSurfTemp;
  auto tempSunlitWall = urban.sunlitWall.EffectiveSurfTemp;
  auto tempShadedWall = urban.shadedWall.EffectiveSurfTemp;
  auto tempImpRoad = urban.imperviousRoad.EffectiveSurfTemp;
  auto tempPerRoad = urban.perviousRoad.EffectiveSurfTemp;

  // Access net longwave radiation fields (to be updated)
  auto netLwSunlitWall = urban.sunlitWall.NetLongRad;
  auto netLwShadedWall = urban.shadedWall.NetLongRad;
  auto netLwImpRoad = urban.imperviousRoad.NetLongRad;
  auto netLwPerRoad = urban.perviousRoad.NetLongRad;
  auto netLwRoof = urban.roof.NetLongRad;

  // Access upward longwave radiation fields (to be updated)
  auto upLwSunlitWall = urban.sunlitWall.UpwardLongRad;
  auto upLwShadedWall = urban.shadedWall.UpwardLongRad;
  auto upLwImpRoad = urban.imperviousRoad.UpwardLongRad;
  auto upLwPerRoad = urban.perviousRoad.UpwardLongRad;
  auto upLwRoof = urban.roof.UpwardLongRad;

  // View to track non-converged landunits
  Kokkos::View<int *> nonConvergedCount("nonConvergedCount", 1);

  Kokkos::parallel_for(
      "ComputeNetLongwave", numLandunits, KOKKOS_LAMBDA(const int l) {
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Computations for roof
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Blend bare-surface emissivity with snow emissivity (walls excluded)
        const Real emissRoofS =
            emissRoof(l) * (1.0 - fracSnoRoof(l)) + snoem * fracSnoRoof(l);
        const Real emissImpRoadS = emissImpRoad(l) * (1.0 - fracSnoImpRoad(l)) +
                                   snoem * fracSnoImpRoad(l);
        const Real emissPerRoadS = emissPerRoad(l) * (1.0 - fracSnoPerRoad(l)) +
                                   snoem * fracSnoPerRoad(l);

        upLwRoof(l) = emissRoofS * STEBOL * Kokkos::pow(tempRoof(l), 4.0) +
                      (1.0 - emissRoofS) * forcLRad(l);
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

        // Initialize impervious and pervious roads (using snow-blended
        // emissivities)
        auto impRoad = InitializeSingleRoad(
            LtotForRoad, emissImpRoadS, tempImpRoad(l), roadVF, fracImpRoad);
        auto perRoad =
            InitializeSingleRoad(LtotForRoad, emissPerRoadS, tempPerRoad(l),
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
              Fluxes(emissImpRoadS, tempImpRoad(l), LtotForRoad, fracImpRoad);
          perRoad.flux = Fluxes(emissPerRoadS, tempPerRoad(l), LtotForRoad,
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
              ReflectRoad(LtotForRoad, emissImpRoadS, roadVF, fracImpRoad);
          perRoad.ref =
              ReflectRoad(LtotForRoad, emissPerRoadS, roadVF, fracPervRoad(l));

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
