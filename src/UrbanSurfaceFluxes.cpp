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

// Lightweight struct to group saturation humidity views for a surface
struct SurfaceQsatData {
  const Array1DR8 &temp;
  const Array1DR8 &es;
  const Array1DR8 &esdt;
  const Array1DR8 &qs;
  const Array1DR8 &qsdt;
};

// Lightweight struct to group canyon air and urban geometry data
struct CanyonAirData {
  Real hwrVal;
  Real forcQ;
  Real forcRho;
  Real wtRoadPerv;
  Real wtRoof;
};

// Lightweight struct to group surface temperature and humidity scalars
struct SurfaceTempHumidData {
  Real qRoof, qRoadImperv, qRoadPerv;
  Real dQdTRoof, dQdTRoadImperv, dQdTRoadPerv;
  Real tRoof, tRoadImperv, tRoadPerv, tSunwall, tShadewall;
};

// Lightweight struct to store surface flux derivatives for implicit solver
// These are derivatives of heat fluxes w.r.t. surface temperature (Jacobian
// terms)
struct SurfaceFluxDerivatives {
  Real cgrndsRoof, cgrndlRoof;
  Real cgrndsRoadPerv, cgrndlRoadPerv;
  Real cgrndsRoadImperv, cgrndlRoadImperv;
  Real cgrndsSunwall, cgrndlSunwall;
  Real cgrndsShadewall, cgrndlShadewall;
};

// Lightweight struct to store flux conductances for a single surface
struct SurfaceConductances {
  Real wtus;      // Sensible heat conductance (scaled) (m/s)
  Real wtuq;      // Latent heat conductance (scaled) (m/s)
  Real wtusUnscl; // Sensible heat conductance (unscaled) (m/s)
  Real wtuqUnscl; // Latent heat conductance (unscaled) (m/s)
};

// Lightweight struct to store intermediate flux conductances for all surfaces
// These are saved from the iteration and used to compute derivatives
struct SurfaceFluxConductances {
  SurfaceConductances roof;
  SurfaceConductances roadPerv;
  SurfaceConductances roadImperv;
  SurfaceConductances sunwall;
  SurfaceConductances shadewall;
};

constexpr Real VKC = 0.4;
constexpr Real ZETAM =
    1.574; // transition point of flux-gradient relation (wind profile)
constexpr Real ZETAT =
    0.465; // transition point of flux-gradient relation (temperature profile)

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
    const Real term1 = std::log(-ZETAM * obu / z0m);
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
    const Real term4 = 0.8 * (std::pow(-ZETAT, 0.333) - std::pow(-zeta, 0.333));

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

  ComputeUstar(zldis, obu, z0m, um, ustar);

  Real u10;
  ComputeU10m(zldis, obu, z0m, um, ustar, u10);

  // Calculate temp1 for the temperature profile
  ComputeRelationForOtherScalarProfiles(zldis, obu, z0h, temp1);

  // Since z0q == z0h, temp2 for the humidity profile is same as
  // that for temperature profile
  temp2 = temp1;

  // Calculate temp1 at 2m for the temperature profile
  ComputeRelationForOtherScalarProfiles(2.0 + z0h, obu, z0h, temp12m);

  // similarly, set temp2 at 2m for humidity profile to be same that
  // for the temperature profile
  temp22m = temp12m;

  Real fmnew;
  if (Kokkos::min(zeta, 1.0) < 0.0) {
    const Real x = std::pow((1.0 - 16.0 * Kokkos::min(zeta, 1.0)), 0.25);
    const Real tmp2 = std::log((1.0 + x * x) / 2.0);
    const Real tmp3 = std::log((1.0 + x) / 2.0);
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
void ComputeQsatForSurfaces(int l, Real forcP, const SurfaceQsatData &roof,
                            const SurfaceQsatData &sunwall,
                            const SurfaceQsatData &shadewall,
                            const SurfaceQsatData &improad,
                            const SurfaceQsatData &perroad) {
  // Roof
  ComputeSurfaceQsat(l, forcP, roof.temp, roof.es, roof.esdt, roof.qs,
                     roof.qsdt);

  // Sunlit wall
  ComputeSurfaceQsat(l, forcP, sunwall.temp, sunwall.es, sunwall.esdt,
                     sunwall.qs, sunwall.qsdt);

  // Shaded wall
  ComputeSurfaceQsat(l, forcP, shadewall.temp, shadewall.es, shadewall.esdt,
                     shadewall.qs, shadewall.qsdt);

  // Impervious road
  ComputeSurfaceQsat(l, forcP, improad.temp, improad.es, improad.esdt,
                     improad.qs, improad.qsdt);

  // Pervious road
  ComputeSurfaceQsat(l, forcP, perroad.temp, perroad.es, perroad.esdt,
                     perroad.qs, perroad.qsdt);
}

// Compute new air temperature and humidity in the canyon
KOKKOS_INLINE_FUNCTION
void ComputeNewTafAndQaf(Real canyonWind, Real thm, Real rahu, Real rawu,
                         const CanyonAirData &canyon,
                         const SurfaceTempHumidData &surfaces, Real qaf,
                         Real fwetRoof, Real fwetRoadImperv, Real &tafNew,
                         Real &qafNew, SurfaceFluxConductances &condcs) {

  const Real qSunwall = 0.0;
  const Real qShadewall = 0.0;

  Real canyonResistance = CPAIR * canyon.forcRho / (11.8 + 4.2 * canyonWind);

  // fwetRoof is now passed as parameter from FractionWet view
  Real fwetRoofLocal = fwetRoof;
  if (qaf > surfaces.qRoof) {
    fwetRoofLocal = 1.0;
  }
  const Real wtusRoof = canyon.wtRoof / canyonResistance;
  const Real wtusRoofUnscl = 1.0 / canyonResistance;
  const Real wtuqRoof = fwetRoofLocal * canyon.wtRoof / canyonResistance;
  const Real wtuqRoofUnscl = fwetRoofLocal / canyonResistance;
  const Real wtusRoadPerv =
      canyon.wtRoadPerv * (1.0 - canyon.wtRoof) / canyonResistance;
  const Real wtusRoadPervUnscl = 1.0 / canyonResistance;
  const Real wtuqRoadPerv =
      canyon.wtRoadPerv * (1.0 - canyon.wtRoof) / canyonResistance;
  const Real wtuqRoadPervUnscl = 1.0 / canyonResistance;

  // fwetRoadImperv is now passed as parameter from FractionWet view
  Real fwetRoadImpervLocal = fwetRoadImperv;
  if (qaf > surfaces.qRoadImperv) {
    fwetRoadImpervLocal = 1.0;
  }
  const Real wtusRoadImperv =
      (1.0 - canyon.wtRoadPerv) * (1.0 - canyon.wtRoof) / canyonResistance;
  const Real wtusRoadImpervUnscl = 1.0 / canyonResistance;
  const Real wtuqRoadImperv = fwetRoadImpervLocal * (1.0 - canyon.wtRoadPerv) *
                              (1.0 - canyon.wtRoof) / canyonResistance;
  const Real wtuqRoadImpervUnscl = fwetRoadImpervLocal * 1.0 / canyonResistance;

  const Real wtusSunwall =
      canyon.hwrVal * (1.0 - canyon.wtRoof) / canyonResistance;
  const Real wtusSunwallUnscl = 1.0 / canyonResistance;
  const Real wtuqSunwall = 0.0;
  const Real wtuqSunwallUnscl = 0.0;

  const Real wtusShadewall =
      canyon.hwrVal * (1.0 - canyon.wtRoof) / canyonResistance;
  const Real wtusShadewallUnscl = 1.0 / canyonResistance;
  const Real wtuqShadewall = 0.0;
  const Real wtuqShadewallUnscl = 0.0;

  const Real tafNumer = thm / rahu + surfaces.tRoof * wtusRoof +
                        surfaces.tRoadPerv * wtusRoadPerv +
                        surfaces.tRoadImperv * wtusRoadImperv +
                        surfaces.tSunwall * wtusSunwall +
                        surfaces.tShadewall * wtusShadewall;

  const Real tafDenom = 1.0 / rahu + wtusRoof + wtusRoadPerv + wtusRoadImperv +
                        wtusSunwall + wtusShadewall;

  const Real qafNumer = canyon.forcQ / rawu + surfaces.qRoof * wtuqRoof +
                        surfaces.qRoadPerv * wtuqRoadPerv +
                        surfaces.qRoadImperv * wtuqRoadImperv +
                        qSunwall * wtuqSunwall + qShadewall * wtuqShadewall;

  const Real qafDenom = 1.0 / rawu + wtuqRoof + wtuqRoadPerv + wtuqRoadImperv +
                        wtuqSunwall + wtuqShadewall;

  tafNew = tafNumer / tafDenom;
  qafNew = qafNumer / qafDenom;

  // Store conductances for later derivative computation
  condcs.roof.wtus = wtusRoof;
  condcs.roof.wtuq = wtuqRoof;
  condcs.roof.wtusUnscl = wtusRoofUnscl;
  condcs.roof.wtuqUnscl = wtuqRoofUnscl;

  condcs.roadPerv.wtus = wtusRoadPerv;
  condcs.roadPerv.wtuq = wtuqRoadPerv;
  condcs.roadPerv.wtusUnscl = wtusRoadPervUnscl;
  condcs.roadPerv.wtuqUnscl = wtuqRoadPervUnscl;

  condcs.roadImperv.wtus = wtusRoadImperv;
  condcs.roadImperv.wtuq = wtuqRoadImperv;
  condcs.roadImperv.wtusUnscl = wtusRoadImpervUnscl;
  condcs.roadImperv.wtuqUnscl = wtuqRoadImpervUnscl;

  condcs.sunwall.wtus = wtusSunwall;
  condcs.sunwall.wtuq = wtuqSunwall;
  condcs.sunwall.wtusUnscl = wtusSunwallUnscl;
  condcs.sunwall.wtuqUnscl = wtuqSunwallUnscl;

  condcs.shadewall.wtus = wtusShadewall;
  condcs.shadewall.wtuq = wtuqShadewall;
  condcs.shadewall.wtusUnscl = wtusShadewallUnscl;
  condcs.shadewall.wtuqUnscl = wtuqShadewallUnscl;
}

// Compute flux derivatives after iteration has converged
KOKKOS_INLINE_FUNCTION
void ComputeFluxDerivatives(const SurfaceFluxConductances &condcs,
                            const SurfaceTempHumidData &surfaces, Real forcRho,
                            Real rahu, Real rawu,
                            SurfaceFluxDerivatives &derivs) {
  // Compute resistance terms and sums from conductances
  const Real wtas = 1.0 / rahu;
  const Real wtaq = 1.0 / rawu;
  const Real wts_sum = wtas + condcs.roof.wtus + condcs.roadPerv.wtus +
                       condcs.roadImperv.wtus + condcs.sunwall.wtus +
                       condcs.shadewall.wtus;
  const Real wtq_sum = wtaq + condcs.roof.wtuq + condcs.roadPerv.wtuq +
                       condcs.roadImperv.wtuq + condcs.sunwall.wtuq +
                       condcs.shadewall.wtuq;

  const Real dQdTRoof = surfaces.dQdTRoof;
  derivs.cgrndsRoof = forcRho * CPAIR *
                      (wtas + condcs.roadPerv.wtus + condcs.roadImperv.wtus +
                       condcs.sunwall.wtus + condcs.shadewall.wtus) *
                      (condcs.roof.wtusUnscl / wts_sum);
  derivs.cgrndlRoof = forcRho *
                      (wtaq + condcs.roadPerv.wtuq + condcs.roadImperv.wtuq +
                       condcs.sunwall.wtuq + condcs.shadewall.wtuq) *
                      (condcs.roof.wtuqUnscl / wtq_sum) * dQdTRoof;

  const Real dQdTRoadPerv = surfaces.dQdTRoadPerv;
  derivs.cgrndsRoadPerv = forcRho * CPAIR *
                          (wtas + condcs.roof.wtus + condcs.roadImperv.wtus +
                           condcs.sunwall.wtus + condcs.shadewall.wtus) *
                          (condcs.roadPerv.wtusUnscl / wts_sum);
  derivs.cgrndlRoadPerv = forcRho *
                          (wtaq + condcs.roof.wtuq + condcs.roadImperv.wtuq +
                           condcs.sunwall.wtuq + condcs.shadewall.wtuq) *
                          (condcs.roadPerv.wtuqUnscl / wtq_sum) * dQdTRoadPerv;

  const Real dQdTRoadImperv = surfaces.dQdTRoadImperv;
  derivs.cgrndsRoadImperv = forcRho * CPAIR *
                            (wtas + condcs.roof.wtus + condcs.roadPerv.wtus +
                             condcs.sunwall.wtus + condcs.shadewall.wtus) *
                            (condcs.roadImperv.wtusUnscl / wts_sum);
  derivs.cgrndlRoadImperv = forcRho *
                            (wtaq + condcs.roof.wtuq + condcs.roadPerv.wtuq +
                             condcs.sunwall.wtuq + condcs.shadewall.wtuq) *
                            (condcs.roadImperv.wtuqUnscl / wtq_sum) *
                            dQdTRoadImperv;

  derivs.cgrndsSunwall = forcRho * CPAIR *
                         (wtas + condcs.roof.wtus + condcs.roadPerv.wtus +
                          condcs.roadImperv.wtus + condcs.shadewall.wtus) *
                         (condcs.sunwall.wtusUnscl / wts_sum);
  derivs.cgrndlSunwall = forcRho *
                         (wtaq + condcs.roof.wtuq + condcs.roadPerv.wtuq +
                          condcs.roadImperv.wtuq + condcs.shadewall.wtuq) *
                         (condcs.sunwall.wtuqUnscl / wtq_sum);

  derivs.cgrndsShadewall = forcRho * CPAIR *
                           (wtas + condcs.roof.wtus + condcs.roadPerv.wtus +
                            condcs.roadImperv.wtus + condcs.sunwall.wtus) *
                           (condcs.shadewall.wtusUnscl / wts_sum);
  derivs.cgrndlShadewall = forcRho *
                           (wtaq + condcs.roof.wtuq + condcs.roadPerv.wtuq +
                            condcs.roadImperv.wtuq + condcs.sunwall.wtuq) *
                           (condcs.shadewall.wtuqUnscl / wtq_sum);
}

// Compute sensible and latent heat fluxes for a surface
KOKKOS_INLINE_FUNCTION
void ComputeSurfaceHeatFluxes(Real taf, Real qaf, Real tSurf, Real qSurf,
                              Real forcRho, Real wtusUnscl, Real wtuqUnscl,
                              Real &eflxShGrnd, Real &qflxEvapSoil,
                              Real &qflxTranEvap) {
  // Temperature and humidity differences (canyon air - surface)
  const Real dth = taf - tSurf;
  const Real dqh = qaf - qSurf;

  // Sensible heat flux (W/m²) [+ to atmosphere]
  eflxShGrnd = -forcRho * CPAIR * wtusUnscl * dth;

  // Latent heat flux (mm H₂O/s) [+ to atmosphere]
  qflxEvapSoil = -forcRho * wtuqUnscl * dqh;

  // Default: no transpiration (overridden for pervious road)
  qflxTranEvap = 0.0;
}

// Compute sensible and latent heat fluxes for pervious road with
// partitioning between soil evaporation and transpiration
KOKKOS_INLINE_FUNCTION
void ComputePerviousRoadHeatFluxes(Real taf, Real qaf, Real tSurf, Real qSurf,
                                   Real forcRho, Real wtusUnscl, Real wtuqUnscl,
                                   Real &eflxShGrnd, Real &qflxEvapSoil,
                                   Real &qflxTranEvap) {
  // Temperature and humidity differences (canyon air - surface)
  const Real dth = taf - tSurf;
  const Real dqh = qaf - qSurf;

  // Sensible heat flux (W/m²) [+ to atmosphere]
  eflxShGrnd = -forcRho * CPAIR * wtusUnscl * dth;

  // Latent heat flux - partition between soil evaporation and transpiration
  // If dew (dqh > 0), assign to soil evaporation
  // For now, we don't have snow fraction or soil moisture data,
  // so we use simplified logic
  if (dqh > 0.0) {
    // Condensation/dew - assign to soil
    qflxEvapSoil = -forcRho * wtuqUnscl * dqh;
    qflxTranEvap = 0.0;
  } else {
    // Evaporation - for now assign to soil
    // TODO: Add soil moisture availability check to partition to transpiration
    qflxEvapSoil = -forcRho * wtuqUnscl * dqh;
    qflxTranEvap = 0.0;
  }
}

// Compute surface fluxes for all urban surfaces
void ComputeSurfaceFluxes(URBANXX::_p_UrbanType &urban) {
  const int numLandunits = urban.numLandunits;

  // Get references to atmospheric data
  auto &forcTemp = urban.atmosphereData.ForcTemp;
  auto &forcPotTemp = urban.atmosphereData.ForcPotTemp;
  auto &forcSpcHumd = urban.atmosphereData.ForcSpcHumd;
  auto &forcPress = urban.atmosphereData.ForcPress;
  auto &forcRho = urban.atmosphereData.ForcRho;
  auto &forcU = urban.atmosphereData.ForcWindU;
  auto &forcV = urban.atmosphereData.ForcWindV;

  // Get references to geometric parameters
  auto &hwr = urban.urbanParams.CanyonHwr;

  // Get references to urban canyon air properties
  auto &Taf = urban.urbanCanyon.Taf;
  auto &Qaf = urban.urbanCanyon.Qaf;

  // Get references to height parameters
  auto &forcHgtT = urban.urbanParams.heights.ForcHgtT;
  auto &forcHgtU = urban.urbanParams.heights.ForcHgtU;
  auto &zDTown = urban.urbanParams.heights.ZDTown;
  auto &z0Town = urban.urbanParams.heights.Z0Town;
  auto &htRoof = urban.urbanParams.heights.HtRoof;
  auto &windHgtCanyon = urban.urbanParams.heights.WindHgtCanyon;

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

        // Prepare surface data for qsat computation
        SurfaceQsatData roofQsat = {urban.roof.EffectiveSurfTemp, urban.roof.Es,
                                    urban.roof.EsdT, urban.roof.Qs,
                                    urban.roof.QsdT};
        SurfaceQsatData sunwallQsat = {
            urban.sunlitWall.EffectiveSurfTemp, urban.sunlitWall.Es,
            urban.sunlitWall.EsdT, urban.sunlitWall.Qs, urban.sunlitWall.QsdT};
        SurfaceQsatData shadewallQsat = {
            urban.shadedWall.EffectiveSurfTemp, urban.shadedWall.Es,
            urban.shadedWall.EsdT, urban.shadedWall.Qs, urban.shadedWall.QsdT};
        SurfaceQsatData improadQsat = {
            urban.imperviousRoad.EffectiveSurfTemp, urban.imperviousRoad.Es,
            urban.imperviousRoad.EsdT, urban.imperviousRoad.Qs,
            urban.imperviousRoad.QsdT};
        SurfaceQsatData perroadQsat = {
            urban.perviousRoad.EffectiveSurfTemp, urban.perviousRoad.Es,
            urban.perviousRoad.EsdT, urban.perviousRoad.Qs,
            urban.perviousRoad.QsdT};

        ComputeQsatForSurfaces(l, forcP, roofQsat, sunwallQsat, shadewallQsat,
                               improadQsat, perroadQsat);

        // Compute canyon wind speed
        Real canyonUWind;
        ComputeCanyonUWind(htRoofVal, zDTownVal, z0TownVal, forcHgtUVal,
                           windHgtCanyonVal, hwrVal, ur, canyonUWind);

        // Iteration loop to compute friction velocity and surface fluxes
        Real fm = 0.0;
        SurfaceFluxConductances condcs;
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

          // Prepare canyon and surface data
          CanyonAirData canyonData = {
              hwrVal, forcQ, forcRho(l),
              urban.urbanParams.FracPervRoadOfTotalRoad(l),
              urban.urbanParams.WtRoof(l)};
          SurfaceTempHumidData surfaceData = {
              urban.roof.Qs(l),
              urban.imperviousRoad.Qs(l),
              urban.perviousRoad.Qs(l),
              urban.roof.QsdT(l),
              urban.imperviousRoad.QsdT(l),
              urban.perviousRoad.QsdT(l),
              urban.roof.EffectiveSurfTemp(l),
              urban.imperviousRoad.EffectiveSurfTemp(l),
              urban.perviousRoad.EffectiveSurfTemp(l),
              urban.sunlitWall.EffectiveSurfTemp(l),
              urban.shadedWall.EffectiveSurfTemp(l)};

          Real tafNew, qafNew;
          const Real fwetRoof = urban.roof.FractionWet(l);
          const Real fwetRoadImperv = urban.imperviousRoad.FractionWet(l);
          ComputeNewTafAndQaf(canyonWind, thm, rahu, rawu, canyonData,
                              surfaceData, qaf, fwetRoof, fwetRoadImperv,
                              tafNew, qafNew, condcs);
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

        // Compute flux derivatives after iteration has converged
        // Need final rahu and rawu from last iteration
        Real ustar;
        Real temp1, temp12m;
        Real temp2, temp22m;
        FrictionVelocity(2, forcHgtUVal, zDTownVal, z0TownVal, obu, ur, um,
                         ustar, temp1, temp12m, temp2, temp22m, fm);
        Real rahu = 1.0 / (temp1 * ustar);
        Real rawu = 1.0 / (temp2 * ustar);

        SurfaceTempHumidData surfaceData = {
            urban.roof.Qs(l),
            urban.imperviousRoad.Qs(l),
            urban.perviousRoad.Qs(l),
            urban.roof.QsdT(l),
            urban.imperviousRoad.QsdT(l),
            urban.perviousRoad.QsdT(l),
            urban.roof.EffectiveSurfTemp(l),
            urban.imperviousRoad.EffectiveSurfTemp(l),
            urban.perviousRoad.EffectiveSurfTemp(l),
            urban.sunlitWall.EffectiveSurfTemp(l),
            urban.shadedWall.EffectiveSurfTemp(l)};
        SurfaceFluxDerivatives derivs;
        ComputeFluxDerivatives(condcs, surfaceData, forcRho(l), rahu, rawu,
                               derivs);

        // Store flux derivatives to views (for implicit heat diffusion solver)
        urban.roof.Cgrnds(l) = derivs.cgrndsRoof;
        urban.roof.Cgrndl(l) = derivs.cgrndlRoof;
        urban.perviousRoad.Cgrnds(l) = derivs.cgrndsRoadPerv;
        urban.perviousRoad.Cgrndl(l) = derivs.cgrndlRoadPerv;
        urban.imperviousRoad.Cgrnds(l) = derivs.cgrndsRoadImperv;
        urban.imperviousRoad.Cgrndl(l) = derivs.cgrndlRoadImperv;
        urban.sunlitWall.Cgrnds(l) = derivs.cgrndsSunwall;
        urban.sunlitWall.Cgrndl(l) = derivs.cgrndlSunwall;
        urban.shadedWall.Cgrnds(l) = derivs.cgrndsShadewall;
        urban.shadedWall.Cgrndl(l) = derivs.cgrndlShadewall;

        // Compute and store surface heat fluxes
        // Roof
        ComputeSurfaceHeatFluxes(
            taf, qaf, urban.roof.EffectiveSurfTemp(l), urban.roof.Qs(l),
            forcRho(l), condcs.roof.wtusUnscl, condcs.roof.wtuqUnscl,
            urban.roof.EflxShGrnd(l), urban.roof.QflxEvapSoil(l),
            urban.roof.QflxTranEvap(l));

        // Impervious road
        ComputeSurfaceHeatFluxes(
            taf, qaf, urban.imperviousRoad.EffectiveSurfTemp(l),
            urban.imperviousRoad.Qs(l), forcRho(l), condcs.roadImperv.wtusUnscl,
            condcs.roadImperv.wtuqUnscl, urban.imperviousRoad.EflxShGrnd(l),
            urban.imperviousRoad.QflxEvapSoil(l),
            urban.imperviousRoad.QflxTranEvap(l));

        // Pervious road (with partitioning logic)
        ComputePerviousRoadHeatFluxes(
            taf, qaf, urban.perviousRoad.EffectiveSurfTemp(l),
            urban.perviousRoad.Qs(l), forcRho(l), condcs.roadPerv.wtusUnscl,
            condcs.roadPerv.wtuqUnscl, urban.perviousRoad.EflxShGrnd(l),
            urban.perviousRoad.QflxEvapSoil(l),
            urban.perviousRoad.QflxTranEvap(l));

        // Sunlit wall (no evaporation - walls don't have evap views)
        {
          Real qflxEvapSoilDummy, qflxTranEvapDummy;
          ComputeSurfaceHeatFluxes(
              taf, qaf, urban.sunlitWall.EffectiveSurfTemp(l), 0.0, forcRho(l),
              condcs.sunwall.wtusUnscl, condcs.sunwall.wtuqUnscl,
              urban.sunlitWall.EflxShGrnd(l), qflxEvapSoilDummy,
              qflxTranEvapDummy);
        }

        // Shaded wall (no evaporation - walls don't have evap views)
        {
          Real qflxEvapSoilDummy, qflxTranEvapDummy;
          ComputeSurfaceHeatFluxes(
              taf, qaf, urban.shadedWall.EffectiveSurfTemp(l), 0.0, forcRho(l),
              condcs.shadewall.wtusUnscl, condcs.shadewall.wtuqUnscl,
              urban.shadedWall.EflxShGrnd(l), qflxEvapSoilDummy,
              qflxTranEvapDummy);
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
  if (0) {
    print_view_1d(urban.urbanCanyon.Taf, "urbanCanyon.Taf");
    print_view_1d(urban.urbanCanyon.Qaf, "urbanCanyon.Qaf");
  }
}

} // namespace URBANXX
