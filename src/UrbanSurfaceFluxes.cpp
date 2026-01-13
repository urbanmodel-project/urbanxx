#include "private/AtmosphereTypeImpl.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include "private/UrbanDebugUtils.h"
#include "private/UrbanParamsTypeImpl.h"
#include "private/UrbanSurfaceFluxesImpl.h"
#include "private/UrbanSurfaceTypeImpl.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>
#include <cmath>
#include <iostream>

namespace URBANXX {

#define GRAVITY 9.80616
#define VKC 0.4
#define ZETAM 1.574 // transition point of flux-gradient relation (wind profile)
#define ZETAT                                                                  \
  0.465 // transition point of flux-gradient relation (temperature profile)
#define CPAIR 1004.64 // specific heat of dry air [J/kg/K]
#define SHR_CONST_TKFRZ 273.15

KOKKOS_INLINE_FUNCTION
void QSat(Real T, Real p, Real &es, Real &esdT, Real &qs, Real &qsdT) {
  // For water vapor (temperature range 0C-100C)
  Real a0 = 6.11213476;
  Real a1 = 0.444007856;
  Real a2 = 0.143064234e-01;
  Real a3 = 0.264461437e-03;
  Real a4 = 0.305903558e-05;
  Real a5 = 0.196237241e-07;
  Real a6 = 0.892344772e-10;
  Real a7 = -0.373208410e-12;
  Real a8 = 0.209339997e-15;
  // For derivative:water vapor
  Real b0 = 0.444017302;
  Real b1 = 0.286064092e-01;
  Real b2 = 0.794683137e-03;
  Real b3 = 0.121211669e-04;
  Real b4 = 0.103354611e-06;
  Real b5 = 0.404125005e-09;
  Real b6 = -0.788037859e-12;
  Real b7 = -0.114596802e-13;
  Real b8 = 0.381294516e-16;
  // For ice (temperature range -75C-0C)
  Real c0 = 6.11123516;
  Real c1 = 0.503109514;
  Real c2 = 0.188369801e-01;
  Real c3 = 0.420547422e-03;
  Real c4 = 0.614396778e-05;
  Real c5 = 0.602780717e-07;
  Real c6 = 0.387940929e-09;
  Real c7 = 0.149436277e-11;
  Real c8 = 0.262655803e-14;
  // For derivative:ice
  Real d0 = 0.503277922;
  Real d1 = 0.377289173e-01;
  Real d2 = 0.126801703e-02;
  Real d3 = 0.249468427e-04;
  Real d4 = 0.313703411e-06;
  Real d5 = 0.257180651e-08;
  Real d6 = 0.133268878e-10;
  Real d7 = 0.394116744e-13;
  Real d8 = 0.498070196e-16;

  Real T_limit = T - SHR_CONST_TKFRZ;
  if (T_limit > 100.0)
    T_limit = 100.0;
  if (T_limit < -75.0)
    T_limit = -75.0;

  Real td = T_limit;

  if (td >= 0.0) {
    es =
        a0 +
        td * (a1 +
              td * (a2 +
                    td * (a3 +
                          td * (a4 +
                                td * (a5 + td * (a6 + td * (a7 + td * a8)))))));

    esdT =
        b0 +
        td * (b1 +
              td * (b2 +
                    td * (b3 +
                          td * (b4 +
                                td * (b5 + td * (b6 + td * (b7 + td * b8)))))));

  } else {
    es =
        c0 +
        td * (c1 +
              td * (c2 +
                    td * (c3 +
                          td * (c4 +
                                td * (c5 + td * (c6 + td * (c7 + td * c8)))))));

    esdT =
        d0 +
        td * (d1 +
              td * (d2 +
                    td * (d3 +
                          td * (d4 +
                                td * (d5 + td * (d6 + td * (d7 + td * d8)))))));
  }

  es = es * 100.0;     // pa
  esdT = esdT * 100.0; // pa/K

  Real vp = 1.0 / (p - 0.378 * es);
  Real vp1 = 0.622 * vp;
  Real vp2 = vp1 * vp;

  qs = es * vp1;         // kg/kg
  qsdT = esdT * vp2 * p; // 1 / K
}

KOKKOS_INLINE_FUNCTION
void MoninObukIni(Real ur, Real thv, Real dthv, Real zldis, Real z0m, Real &um,
                  Real &obu) {

  // const Real ustar = 0.06;
  const Real wc = 0.5;

  if (dthv > 0) {
    um = Kokkos::max(ur, 1.0);
  } else {
    um = std::pow(std::pow(ur, 2.0) + std::pow(wc, 2.0), 0.5);
  }

  const Real rib = GRAVITY * zldis * dthv / (thv * std::pow(um, 2.0));

  Real zeta;
  if (rib >= 0) {
    zeta = rib * std::log(zldis / z0m) / (1.0 - 0.5 * Kokkos::min(rib, 0.19));
    zeta = Kokkos::min(2.0, Kokkos::max(zeta, 0.01));
  } else {
    zeta = rib * std::log(zldis / z0m);
    zeta = Kokkos::max(-100.0, Kokkos::min(zeta, -0.01));
  }

  obu = zldis / zeta;
}

KOKKOS_INLINE_FUNCTION
void StabilityFunc1(Real zeta, Real &value) {

  Real chik2 = std::pow(1.0 - 16.0 * zeta, 0.5);
  Real chik = std::pow(chik2, 0.5);

  Real term1 = 2.0 * std::log((1.0 + chik) * 0.5);
  Real term2 = std::log((1.0 + chik2) * 0.5);
  Real term3 = 2.0 * std::atan(chik);
  value = term1 + term2 + term3 + M_PI * 0.5;
}

KOKKOS_INLINE_FUNCTION
void StabilityFunc2(Real zeta, Real &value) {

  Real chik2 = std::pow(1.0 - 16.0 * zeta, 0.5);
  value = 2.0 * std::log((1.0 + chik2) * 0.5);
}

KOKKOS_INLINE_FUNCTION
void ComputeUstar(Real zldis, Real obu, Real z0m, Real um, Real &ustar) {

  const Real zeta = zldis / obu;

  if (zeta < -ZETAM) {
    const Real term1 = std::log(ZETAM * obu / z0m);
    const Real term4 =
        1.14 * (std::pow(-zeta, 0.333) - std::pow(-ZETAM, 0.333));

    Real term2, term3;
    StabilityFunc1(-ZETAM, term2);
    StabilityFunc1(z0m / obu, term3);

    ustar = VKC * um / (term1 - term2 + term3 + term4);

  } else if (zeta < 0.0) {

    const Real term1 = std::log(zldis / z0m);

    Real term2, term3;
    StabilityFunc1(zeta, term2);
    StabilityFunc1(z0m / obu, term3);

    ustar = VKC * um / (term1 - term2 + term3);

  } else if (zeta <= 1.0) {

    const Real denom = std::log(zldis / z0m) + 5.0 * zeta - 5.0 * z0m / obu;
    ustar = VKC * um / denom;

  } else {

    const Real denom = std::log(obu / z0m) + 5.0 - 5.0 * z0m / obu +
                       (5.0 * std::log(zeta) + zeta - 1.0);
    ustar = VKC * um / denom;
  }
}

KOKKOS_INLINE_FUNCTION
void ComputeU10m(Real zldis, Real obu, Real z0m, Real um, Real ustar,
                 Real &u10) {

  const Real zeta = zldis / obu;

  if (zldis - z0m <= 10.0) {

    u10 = um;

  } else {

    if (zeta < -ZETAM) {

      const Real term1 = std::log(-ZETAM * obu / (10.0 + z0m));
      const Real term4 =
          1.14 * (std::pow(-zeta, 0.333) - std::pow(-ZETAM, 0.333));

      Real term2, term3;
      StabilityFunc1(-ZETAM, term2);
      StabilityFunc1((10.0 + z0m) / obu, term3);

      u10 = um - ustar / VKC * (term1 - term2 + term3 + term4);

    } else if (zeta < 0.0) {

      const Real term1 = std::log(zldis / (10.0 + z0m));

      Real term2, term3;
      StabilityFunc1(zeta, term2);
      StabilityFunc1((10.0 + z0m) / obu, term3);

      u10 = um - ustar / VKC * (term1 - term2 + term3);

    } else if (zeta <= 1.0) {

      const Real term1 = std::log(zldis / (10.0 + z0m));
      const Real term2 = -5.0 * zeta;
      const Real term3 = -5.0 * (10.0 + z0m) / obu;

      u10 = um - ustar / VKC * (term1 - term2 + term3);

    } else {
      const Real term1 = std::log(obu / (10.0 + z0m)) + 5.0;
      const Real term2 = 5.0 * (10.0 + z0m) / obu;
      const Real term3 = (5.0 * std::log(zeta) + zeta - 1.0);

      u10 = um - ustar / VKC * (term1 - term2 + term3);
    }
  }
}

KOKKOS_INLINE_FUNCTION
void ComputeRelationForOtherScalarProfiles(Real zldis, Real obu, Real z0,
                                           Real &value) {

  const Real zeta = zldis / obu;

  if (zeta < -ZETAT) {

    const Real term1 = std::log(-zeta * obu / z0);
    const Real term4 = 0.8 * (std::pow(-ZETAT, 0.333) - std::pow(zeta, 0.0333));

    Real term2, term3;
    StabilityFunc2(-ZETAT, term2);
    StabilityFunc2(z0 / obu, term3);

    value = VKC / (term1 - term2 + term3 + term4);

  } else if (zeta < 0.0) {

    const Real term1 = std::log(zldis / z0);

    Real term2, term3;
    StabilityFunc2(zeta, term2);
    StabilityFunc2(z0 / obu, term3);

    value = VKC / (term1 - term2 + term3);

  } else if (zeta < 1.0) {

    const Real term1 = std::log(zldis / z0) + 5.0 * zeta;
    const Real term2 = 5.0 * z0 / obu;

    value = VKC / (term1 - term2);

  } else {

    const Real term1 = std::log(obu / z0) + 5.0;
    const Real term2 = (5.0 * z0 / obu);
    const Real term3 = 5.0 * std::log(zeta) + zeta - 1.0;

    value = VKC / (term1 - term2 + term3);
  }
}

KOKKOS_INLINE_FUNCTION
void FrictionVelocity(int iter, Real forc_hgt_u, Real displa, Real z0, Real obu,
                      Real ur, Real um, Real &ustar, Real &temp1, Real &temp12m,
                      Real &temp2, Real &temp22m, Real &fm) {

  const Real zldis = forc_hgt_u - displa;
  const Real zeta = zldis / obu;

  const Real z0m = z0;
  const Real z0h = z0;
  // const Real z0q = z0;

  ComputeUstar(zldis, obu, z0m, um, ustar);

  /*
  Real vds;
  if (zeta < 0) {
    vds = 2.0e-3 * ustar * (1.0 + std::pow(300.0 / (-obu), 0.666));
  } else {
    vds = 2.0e-3 * ustar;
  }
  */

  Real u10;
  ComputeU10m(zldis, obu, z0m, um, ustar, u10);

  // Calcualte temp1 for the temperature profile
  ComputeRelationForOtherScalarProfiles(zldis, obu, z0h, temp1);

  // Since z0q == z0h, temp2 for the humidity profile is same as
  // that for temperature profile
  temp2 = temp1;

  // Calcualte temp1 at 2m for the temperature profile
  ComputeRelationForOtherScalarProfiles(2.0 + z0h, obu, z0h, temp12m);

  // similarly, set temp2 at 2m for humidity profile to be same that
  // for the temperature profile
  temp22m = temp12m;

  Real fmnew;
  if (Kokkos::min(zeta, 1.0) < 0.0) {
    const Real x = std::pow((1.0 - 16.0 * Kokkos::min(zeta, 1.0)), 0.25);
    const Real tmp2 = std::log((1.0 + x * x) / 2.0);
    const Real tmp3 = std::log((!.0 + x) / 2.0);
    fmnew = 2.0 * tmp3 + tmp2 - std::atan(x) + M_PI / 2.0;
  } else {
    fmnew = -5.0 * Kokkos::min(zeta, 1.0);
  }

  if (iter == 0)
    fm = fmnew;
  else
    fm = 0.5 * (fm + fmnew);
}

KOKKOS_INLINE_FUNCTION
void ComputeCanyonUWind(Real ht_roof, Real z_d_town, Real z_0_town,
                        Real forc_hgt_u, Real wind_hgt_canyon, Real hwr,
                        Real ur, Real &canyon_u_wind) {

  // wind at cannyon top
  const Real term1 = std::log((ht_roof - z_d_town) / z_0_town);
  const Real term2 = std::log((forc_hgt_u - z_d_town) / z_0_town);
  const Real canyontop_wind = ur * term1 / term2;

  const Real factor = std::exp(-0.5 * hwr * (1.0 - wind_hgt_canyon / ht_roof));

  if (hwr < 0.5) {

    canyon_u_wind = canyontop_wind * factor; // eqn 3.61

  } else if (hwr < 1.0) {

    const Real factor2 = 1.0 + 2.0 * (2.0 / M_PI - 1.0) * (hwr - 0.5);
    canyon_u_wind = canyontop_wind * factor2 * factor; // eqn 3.62

  } else {

    canyon_u_wind = canyontop_wind * 2.0 / M_PI * factor; // eqn 3.60
  }
}

// Helper function to compute saturation humidity for a single surface
KOKKOS_INLINE_FUNCTION
void ComputeSurfaceQsat(int l, Real p, const Array1DR8 &temp,
                        const Array1DR8 &es, const Array1DR8 &esdT,
                        const Array1DR8 &qs, const Array1DR8 &qsdT) {
  Real T = temp(l);
  Real es_val, esdT_val, qs_val, qsdT_val;
  QSat(T, p, es_val, esdT_val, qs_val, qsdT_val);
  es(l) = es_val;
  esdT(l) = esdT_val;
  qs(l) = qs_val;
  qsdT(l) = qsdT_val;
}

// Compute saturation humidity for all urban surfaces
KOKKOS_INLINE_FUNCTION
void ComputeQsatForSurfaces(int l, Real forcP,
                            const URBANXX::_p_UrbanType &urban) {
  // Roof
  ComputeSurfaceQsat(l, forcP, urban.roof.Temperature, urban.roof.Es,
                     urban.roof.EsdT, urban.roof.Qs, urban.roof.QsdT);

  // Sunlit wall
  ComputeSurfaceQsat(l, forcP, urban.sunlitWall.Temperature,
                     urban.sunlitWall.Es, urban.sunlitWall.EsdT,
                     urban.sunlitWall.Qs, urban.sunlitWall.QsdT);

  // Shaded wall
  ComputeSurfaceQsat(l, forcP, urban.shadedWall.Temperature,
                     urban.shadedWall.Es, urban.shadedWall.EsdT,
                     urban.shadedWall.Qs, urban.shadedWall.QsdT);

  // Impervious road
  ComputeSurfaceQsat(l, forcP, urban.imperviousRoad.Temperature,
                     urban.imperviousRoad.Es, urban.imperviousRoad.EsdT,
                     urban.imperviousRoad.Qs, urban.imperviousRoad.QsdT);

  // Pervious road
  ComputeSurfaceQsat(l, forcP, urban.perviousRoad.Temperature,
                     urban.perviousRoad.Es, urban.perviousRoad.EsdT,
                     urban.perviousRoad.Qs, urban.perviousRoad.QsdT);
}

// Compute new air temperature and humidity in the canyon
KOKKOS_INLINE_FUNCTION
void ComputeNewTafAndQaf(int l, Real canyonWind, Real thm, Real rahu, Real rawu,
                         const URBANXX::_p_UrbanType &urban, Real qaf,
                         Real &tafNew, Real &qafNew) {

  const Real hwrVal = urban.urbanParams.CanyonHwr(l);

  const Real qRoof = urban.roof.Qs(l);
  const Real qRoadImperv = urban.imperviousRoad.Qs(l);
  const Real qRoadPerv = urban.perviousRoad.Qs(l);
  const Real qSunwall = 0.0;
  const Real qShadewall = 0.0;

  const Real tRoof = urban.roof.Temperature(l);
  const Real tRoadImperv = urban.imperviousRoad.Temperature(l);
  const Real tRoadPerv = urban.perviousRoad.Temperature(l);
  const Real tSunwall = urban.sunlitWall.Temperature(l);
  const Real tShadewall = urban.shadedWall.Temperature(l);

  const Real forcQ = urban.atmosphereData.ForcSpcHumd(l);
  const Real forcRho = urban.atmosphereData.ForcRho(l);
  Real canyonResistance = CPAIR * forcRho / (11.8 + 4.2 * canyonWind);

  const Real wtRoadPerv = urban.urbanParams.FracPervRoadOfTotalRoad(l);
  const Real wtRoof = urban.urbanParams.WtRoof(l);

  const Real fwetRoof = 0.0;
  const Real wtusRoof = wtRoof / canyonResistance;
  const Real wtuqRoof = fwetRoof * wtRoof / canyonResistance;

  const Real wtusRoadPerv = wtRoadPerv * (1.0 - wtRoof) / canyonResistance;
  const Real wtuqRoadPerv = wtRoadPerv * (1.0 - wtRoof) / canyonResistance;

  const Real fwetRoadImperv = (qaf > qRoadImperv) ? 1.0 : 0.0;
  const Real wtusRoadImperv =
      (1.0 - wtRoadPerv) * (1.0 - wtRoof) / canyonResistance;
  const Real wtuqRoadImperv =
      fwetRoadImperv * (1.0 - wtRoadPerv) * (1.0 - wtRoof) / canyonResistance;

  const Real wtusSunwall = hwrVal * (1.0 - wtRoof) / canyonResistance;
  const Real wtuqSunwall = 0.0;

  const Real wtusShadewall = hwrVal * (1.0 - wtRoof) / canyonResistance;
  const Real wtuqShadewall = 0.0;

  const Real tafNumer = thm / rahu + tRoof * wtusRoof +
                        tRoadPerv * wtusRoadPerv +
                        tRoadImperv * wtusRoadImperv + tSunwall * wtusSunwall +
                        tShadewall * wtusShadewall;

  const Real tafDenom = 1.0 / rahu + wtusRoof + wtusRoadPerv + wtusRoadImperv +
                        wtusSunwall + wtusShadewall;

  const Real qafNumer = forcQ / rawu + qRoof * wtuqRoof +
                        qRoadPerv * wtuqRoadPerv +
                        qRoadImperv * wtuqRoadImperv + qSunwall * wtuqSunwall +
                        qShadewall * wtuqShadewall;

  const Real qafDenom = 1.0 / rawu + wtuqRoof + wtuqRoadPerv + wtuqRoadImperv +
                        wtuqSunwall + wtuqShadewall;

  tafNew = tafNumer / tafDenom;
  qafNew = qafNumer / qafDenom;
}

// Compute surface fluxes for all urban surfaces
void ComputeSurfaceFluxes(URBANXX::_p_UrbanType &urban) {
  const int numLandunits = urban.numLandunits;

  // Get references to atmospheric data
  auto forcTemp = urban.atmosphereData.ForcTemp;
  auto forcPotTemp = urban.atmosphereData.ForcPotTemp;
  auto forcSpcHumd = urban.atmosphereData.ForcSpcHumd;
  auto forcPress = urban.atmosphereData.ForcPress;
  auto forcRho = urban.atmosphereData.ForcRho;
  auto forcU = urban.atmosphereData.ForcWindU;
  auto forcV = urban.atmosphereData.ForcWindV;

  // Get references to geometric parameters
  auto hwr = urban.urbanParams.CanyonHwr;

  // Get references to urban canyon air properties
  auto Taf = urban.urbanCanyon.Taf;
  auto Qaf = urban.urbanCanyon.Qaf;

  // Get references to height parameters
  auto forcHgtT = urban.urbanParams.heights.ForcHgtT;
  auto forcHgtU = urban.urbanParams.heights.ForcHgtU;
  auto zDTown = urban.urbanParams.heights.ZDTown;
  auto z0Town = urban.urbanParams.heights.Z0Town;
  auto htRoof = urban.urbanParams.heights.HtRoof;
  auto windHgtCanyon = urban.urbanParams.heights.WindHgtCanyon;

  // Constants
  const Real lapseRate = 0.0098; // dry adiabatic lapse rate (K/m)

  // Compute surface fluxes for each landunit
  Kokkos::parallel_for(
      "ComputeSurfaceFluxes", numLandunits, KOKKOS_LAMBDA(const int l) {
        // Get atmospheric forcing data
        Real taf = Taf(l);
        Real qaf = Qaf(l);

        const Real forcUVal = forcU(l);
        const Real forcVVal = forcV(l);
        const Real u2PlusV2 = std::pow(forcUVal, 2.0) + std::pow(forcVVal, 2.0);
        const Real velocity = std::pow(u2PlusV2, 0.5);
        const Real ur = Kokkos::max(1.0, velocity);

        const Real hwrVal = hwr(l);

        // Get height parameters for this landunit
        const Real forcHgtTVal = forcHgtT(l);
        const Real forcHgtUVal = forcHgtU(l);
        const Real zDTownVal = zDTown(l);
        const Real z0TownVal = z0Town(l);
        const Real htRoofVal = htRoof(l);
        const Real windHgtCanyonVal = windHgtCanyon(l);

        // Initialize Monin-Obukhov variables
        Real um, obu;
        Real thm, thv, zldis;
        const Real forcQ = forcSpcHumd(l);
        const Real forcTh = forcPotTemp(l);

        {
          const Real forcT = forcTemp(l);

          thm = forcT + lapseRate * forcHgtTVal;
          thv = forcTh * (1.0 + 0.61 * forcQ);
          const Real dth = thm - taf;
          const Real dqh = forcQ - qaf;
          const Real dthv = dth * (1.0 + 0.61 * forcQ) + 0.61 * forcTh * dqh;
          zldis = forcHgtUVal - zDTownVal;

          MoninObukIni(ur, thv, dthv, zldis, z0TownVal, um, obu);
        }

        // Get atmospheric pressure and compute saturation humidity for surfaces
        const Real forcP = forcPress(l);
        ComputeQsatForSurfaces(l, forcP, urban);

        // Compute canyon wind speed
        Real canyonUWind;
        ComputeCanyonUWind(htRoofVal, zDTownVal, z0TownVal, forcHgtUVal,
                           windHgtCanyonVal, hwrVal, ur, canyonUWind);

        // Iteration loop to compute friction velocity and surface fluxes
        Real fm = 0.0;
        for (int iter = 0; iter < 3; ++iter) {
          Real ustar;
          Real temp1, temp12m;
          Real temp2, temp22m;
          FrictionVelocity(iter, forcHgtUVal, zDTownVal, z0TownVal, obu, ur, um,
                           ustar, temp1, temp12m, temp2, temp22m, fm);

          // Real ramu = 1.0 / (ustar * ustar / um);
          Real rahu = 1.0 / (temp1 * ustar);
          Real rawu = 1.0 / (temp2 * ustar);

          Real canyonWindPow2 =
              std::pow(canyonUWind, 2.0) + std::pow(ustar, 2.0);
          Real canyonWind = std::pow(canyonWindPow2, 0.5);

          Real tafNew, qafNew;
          ComputeNewTafAndQaf(l, canyonWind, thm, rahu, rawu, urban, qaf,
                              tafNew, qafNew);
          taf = tafNew;
          qaf = qafNew;

          const Real dth = thm - taf;
          const Real dqh = forcQ - qaf;
          const Real tstar = temp1 * dth;
          const Real qstar = temp2 * dqh;
          const Real thvstar =
              tstar * (1.0 + 0.61 * forcQ) + 0.61 * forcTh * qstar;
          Real zeta =
              zldis * VKC * GRAVITY * thvstar / (std::pow(ustar, 2.0) * thv);

          if (zeta > 0.0) {
            zeta = Kokkos::min(2.0, Kokkos::max(zeta, 0.01));
            um = Kokkos::max(ur, 0.1);
          } else {
            const Real beta = 1.0;
            const Real zii = 1000.0;
            const Real wc =
                beta * std::pow(-GRAVITY * ustar * thvstar * zii / thv, 0.333);
            um = std::pow(ur * ur + wc * wc, 0.5);
          }
          obu = zldis / zeta;
        }

        // Store taf and qaf back to arrays
        Taf(l) = taf;
        Qaf(l) = qaf;
      });

  Kokkos::fence();

  // Debug output (disabled by default)
  if (0) {
    std::cout << "ComputeSurfaceFluxes completed for " << numLandunits
              << " landunits" << std::endl;
  }
}

} // namespace URBANXX
