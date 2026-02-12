#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include "private/UrbanHydrologyFunctions.h"
#include "private/UrbanHydrologyImpl.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"
#include <Kokkos_Core.hpp>

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

using namespace URBANXX;

namespace URBANXX {

// ============================================================================
// Phase 3.1: Preprocessing - Compute Hydraulic Properties
// ============================================================================

void ComputeHydraulicProperties(UrbanType urban) {
  const int nlandunits = urban->numLandunits;
  const int nlevbed = NUM_SOIL_LAYERS;

  auto h2osoi_liq = urban->perviousRoad.H2OSoiLiq;
  auto h2osoi_ice = urban->perviousRoad.H2OSoiIce;
  auto h2osoi_vol = urban->perviousRoad.H2OSoiVol;
  auto eff_porosity = urban->perviousRoad.EffPorosity;
  auto hk = urban->perviousRoad.Hk;
  auto smp = urban->perviousRoad.Smp;
  auto zwt = urban->perviousRoad.Zwt;
  auto jwt = urban->perviousRoad.Jwt;

  auto watsat = urban->perviousRoad.soil.WatSat;
  auto hksat = urban->perviousRoad.soil.HkSat;
  auto bsw = urban->perviousRoad.soil.Bsw;
  auto sucsat = urban->perviousRoad.soil.SucSat;
  auto dz = urban->perviousRoad.Dz;
  auto zi = urban->perviousRoad.Zi;

  Kokkos::parallel_for(
      "ComputeHydraulicProperties", nlandunits, KOKKOS_LAMBDA(const int l) {
        // Convert water table depth from m to mm
        const Real zwtmm = zwt(l) * 1000.0;

        // Find jwt: layer index right above water table
        jwt(l) = nlevbed - 1; // Initialize to bottom layer
        for (int j = 0; j < nlevbed; ++j) {
          if (zwtmm <= zi(l, j) * 1000.0) {
            jwt(l) = (j == 0) ? 0 : j - 1;
            break;
          }
        }

        // Compute properties for each layer
        for (int j = 0; j < nlevbed; ++j) {
          // Compute ice fraction
          const Real vol_ice = Kokkos::fmin(
              watsat(l, j), h2osoi_ice(l, j) / (dz(l, j) * SHR_CONST_RHOICE));
          const Real icefrac = Kokkos::fmin(1.0, vol_ice / watsat(l, j));

          // Compute volumetric liquid water content
          const Real vol_liq = Kokkos::fmax(h2osoi_liq(l, j), 1.0e-6) /
                               (dz(l, j) * SHR_CONST_RHOWATER);

          // Update volumetric water content and effective porosity
          h2osoi_vol(l, j) = vol_liq + vol_ice;
          eff_porosity(l, j) = Kokkos::fmax(0.01, watsat(l, j) - vol_ice);

          // Compute relative saturation for next layer interface
          // (average of current and next layer for interface properties)
          const int j_next = Kokkos::fmin(j + 1, NUM_LAYERS_ABV_BEDROCK - 1);
          const Real s_interface = (h2osoi_vol(l, j) + h2osoi_vol(l, j_next)) /
                                   (watsat(l, j) + watsat(l, j_next));

          // Compute ice impedance for interface
          const Real icefrac_next = Kokkos::fmin(
              1.0, Kokkos::fmin(watsat(l, j_next),
                                h2osoi_ice(l, j_next) /
                                    (dz(l, j_next) * SHR_CONST_RHOICE)) /
                       watsat(l, j_next));
          const Real imped =
              ComputeIceImpedance(0.5 * (icefrac + icefrac_next));

          // Compute hydraulic conductivity at interface
          hk(l, j) = ComputeHydraulicConductivity(s_interface, hksat(l, j),
                                                  bsw(l, j), imped);
          if (j >= NUM_LAYERS_ABV_BEDROCK) {
            hk(l, j) = 0.0; // Hydrologically inactive layers below layer
                            // NUM_LAYERS_ABV_BEDROCK
          }

          // Compute soil matric potential at node
          const Real s_node = Kokkos::fmax(vol_liq / watsat(l, j), 0.01);
          smp(l, j) =
              ComputeMatricPotential(s_node, sucsat(l, j), bsw(l, j), SMP_MIN);
        }
      });
  Kokkos::fence();
}

// ============================================================================
// Phase 3.2: Build Tridiagonal System
// ============================================================================

void SetupHydrologyTridiagonal(UrbanType urban, Real dtime) {
  const int nlandunits = urban->numLandunits;
  const int nlevbed = NUM_LAYERS_ABV_BEDROCK;

  auto h2osoi_liq = urban->perviousRoad.H2OSoiLiq;
  auto h2osoi_ice = urban->perviousRoad.H2OSoiIce;
  auto h2osoi_vol = urban->perviousRoad.H2OSoiVol;
  auto hk = urban->perviousRoad.Hk;
  auto smp = urban->perviousRoad.Smp;
  auto qflx_infl = urban->perviousRoad.QflxInfl;
  auto qflx_tran = urban->perviousRoad.QflxTran;
  auto zwt = urban->perviousRoad.Zwt;
  auto jwt = urban->perviousRoad.Jwt;

  auto watsat = urban->perviousRoad.soil.WatSat;
  auto hksat = urban->perviousRoad.soil.HkSat;
  auto bsw = urban->perviousRoad.soil.Bsw;
  auto sucsat = urban->perviousRoad.soil.SucSat;
  auto dz = urban->perviousRoad.Dz;
  auto zi = urban->perviousRoad.Zi;
  auto Zc = urban->perviousRoad.Zc;

  // Tridiagonal matrix arrays (stored in soil data structure for now)
  auto amx = urban->perviousRoad.soil.Amx;
  auto bmx = urban->perviousRoad.soil.Bmx;
  auto cmx = urban->perviousRoad.soil.Cmx;
  auto rmx = urban->perviousRoad.soil.Rmx;

  Kokkos::parallel_for(
      "SetupHydrologyTridiagonal", nlandunits, KOKKOS_LAMBDA(const int l) {
        const Real zwtmm = zwt(l) * 1000.0;
        const int jwt_l = jwt(l);

        // Compute equilibrium matric potentials for each layer
        Real zq[NUM_SOIL_LAYERS + 1];
        for (int j = 0; j < nlevbed; ++j) {
          const Real z_top = zi(l, j) * 1000.0;
          const Real z_bot = zi(l, j + 1) * 1000.0;

          const Real vol_eq = ComputeEquilibriumWaterContent(
              zwtmm, z_top, z_bot, watsat(l, j), sucsat(l, j), bsw(l, j));

          zq[j] = ComputeEquilibriumMatricPotential(
              vol_eq, watsat(l, j), sucsat(l, j), bsw(l, j), SMP_MIN);
        }

        // Aquifer layer (if water table below soil column)
        if (jwt_l == nlevbed) {
          const Real z_top = zi(l, nlevbed) * 1000.0;
          const Real delta_z_zwt = Kokkos::fmax(zwtmm - z_top, 1.0);

          // Cannot use ComputeEquilibriumWaterContent for aquifer layer since
          // the algorithm here is slightly different than what is used for soil
          // layers
          const Real zwt_loc = zwtmm;
          const Real watsat_loc = watsat(l, nlevbed - 1);
          const Real sucsat_loc = sucsat(l, nlevbed - 1);
          const Real bsw_loc = bsw(l, nlevbed - 1);
          const Real temp_i = 1.0;
          const Real temp_0 = Kokkos::pow(
              (sucsat_loc + zwt_loc - z_top) / sucsat_loc, 1.0 - 1.0 / bsw_loc);
          const Real vol_eq1 = -sucsat_loc * watsat_loc /
                               (1.0 - 1.0 / bsw_loc) / (delta_z_zwt) *
                               (temp_i - temp_0);
          const Real vol_eq =
              Kokkos::fmin(watsat_loc, Kokkos::fmax(vol_eq1, 0.0));

          zq[nlevbed] = ComputeEquilibriumMatricPotential(
              vol_eq, watsat(l, nlevbed - 1), sucsat(l, nlevbed - 1),
              bsw(l, nlevbed - 1), SMP_MIN);
        }

        // Compute derivatives of hydraulic conductivity and matric potential
        Real dhkdw[NUM_SOIL_LAYERS];
        Real dsmpdw[NUM_SOIL_LAYERS];

        for (int j = 0; j < nlevbed; ++j) {
          const Real vol_liq = Kokkos::fmax(h2osoi_liq(l, j), 1.0e-6) /
                               (dz(l, j) * SHR_CONST_RHOWATER);
          const int j_next = Kokkos::fmin(j + 1, nlevbed - 1);
          const Real vol_liq_next =
              Kokkos::fmax(h2osoi_liq(l, j_next), 1.0e-6) /
              (dz(l, j_next) * SHR_CONST_RHOWATER);

          const Real s_interface =
              0.5 * (vol_liq / watsat(l, j) + vol_liq_next / watsat(l, j_next));

          // Compute ice impedance
          const Real vol_ice = Kokkos::fmin(
              watsat(l, j), h2osoi_ice(l, j) / (dz(l, j) * SHR_CONST_RHOICE));
          const Real icefrac = Kokkos::fmin(1.0, vol_ice / watsat(l, j));
          const Real vol_ice_next = Kokkos::fmin(
              watsat(l, j_next),
              h2osoi_ice(l, j_next) / (dz(l, j_next) * SHR_CONST_RHOICE));
          const Real icefrac_next =
              Kokkos::fmin(1.0, vol_ice_next / watsat(l, j_next));
          const Real imped =
              ComputeIceImpedance(0.5 * (icefrac + icefrac_next));

          // NOTE: 0.5 * (watsat(l, j) + watsat(l, j_next)) is not used
          // to be consistent with code in ELM. It is possible that the 0.5
          // neglected after performing algebraic manupilations in the
          // descrtized equations in ELM, but this needs to be verified.
          dhkdw[j] = ComputeHydraulicConductivityDerivative(
              s_interface, hksat(l, j), bsw(l, j), imped,
              (watsat(l, j) + watsat(l, j_next)));

          dsmpdw[j] =
              ComputeMatricPotentialDerivative(smp(l, j), vol_liq, bsw(l, j));
        }

        // Initialize tridiagonal matrix coefficients for layers below soil
        // column
        for (int j = nlevbed; j < NUM_SOIL_LAYERS; ++j) {
          amx(l, j) = 0.0;
          bmx(l, j) = 1.0;
          cmx(l, j) = 0.0;
          rmx(l, j) = 0.0;
        }

        // Layer j=0 (top layer)
        {
          const int j = 0;
          const Real den = (Zc(l, j + 1) - Zc(l, j)) * 1000.0;
          const Real dzq = zq[j + 1] - zq[j];
          const Real num = (smp(l, j + 1) - smp(l, j)) - dzq;

          const Real qin = qflx_infl(l);
          const Real qout = -hk(l, j) * num / den;

          const Real dqodw1 = -(-hk(l, j) * dsmpdw[j] + num * dhkdw[j]) / den;
          const Real dqodw2 =
              -(hk(l, j) * dsmpdw[j + 1] + num * dhkdw[j]) / den;

          rmx(l, j) = qin - qout - qflx_tran(l, j);
          amx(l, j) = 0.0;
          bmx(l, j) = dz(l, j) * 1000.0 / dtime + dqodw1;
          cmx(l, j) = dqodw2;
        }

        // Interior layers j=1 to nlevbed-2
        for (int j = 1; j < nlevbed - 1; ++j) {
          // Inflow from above
          Real den = (Zc(l, j) - Zc(l, j - 1)) * 1000.0;
          Real dzq = zq[j] - zq[j - 1];
          Real num = (smp(l, j) - smp(l, j - 1)) - dzq;

          const Real qin = -hk(l, j - 1) * num / den;
          const Real dqidw0 =
              -(-hk(l, j - 1) * dsmpdw[j - 1] + num * dhkdw[j - 1]) / den;
          const Real dqidw1 =
              -(hk(l, j - 1) * dsmpdw[j] + num * dhkdw[j - 1]) / den;

          // Outflow to below
          den = (Zc(l, j + 1) - Zc(l, j)) * 1000.0;
          dzq = zq[j + 1] - zq[j];
          num = (smp(l, j + 1) - smp(l, j)) - dzq;

          const Real qout = -hk(l, j) * num / den;
          const Real dqodw1 = -(-hk(l, j) * dsmpdw[j] + num * dhkdw[j]) / den;
          const Real dqodw2 =
              -(hk(l, j) * dsmpdw[j + 1] + num * dhkdw[j]) / den;

          rmx(l, j) = qin - qout - qflx_tran(l, j);
          amx(l, j) = -dqidw0;
          bmx(l, j) = dz(l, j) * 1000.0 / dtime - dqidw1 + dqodw1;
          cmx(l, j) = dqodw2;
        }

        // Bottom layer j=nlevbed-1
        {
          const int j = nlevbed - 1;

          if (jwt_l < nlevbed - 1) {
            // Water table is in soil column - zero flow bottom boundary
            Real den = (Zc(l, j) - Zc(l, j - 1)) * 1000.0;
            Real dzq = zq[j] - zq[j - 1];
            Real num = (smp(l, j) - smp(l, j - 1)) - dzq;

            const Real qin = -hk(l, j - 1) * num / den;
            const Real dqidw0 =
                -(-hk(l, j - 1) * dsmpdw[j - 1] + num * dhkdw[j - 1]) / den;
            const Real dqidw1 =
                -(hk(l, j - 1) * dsmpdw[j] + num * dhkdw[j - 1]) / den;

            const Real qout = 0.0;
            const Real dqodw1 = 0.0;

            rmx(l, j) = qin - qout - qflx_tran(l, j);
            amx(l, j) = -dqidw0;
            bmx(l, j) = dz(l, j) * 1000.0 / dtime - dqidw1 + dqodw1;
            cmx(l, j) = 0.0;

            // Aquifer layer is hydrologically inactive
            rmx(l, j + 1) = 0.0;
            amx(l, j + 1) = 0.0;
            bmx(l, j + 1) = dz(l, j) * 1000.0 / dtime;
            cmx(l, j + 1) = 0.0;

          } else {
            // Water table below soil column - include aquifer layer
            // Bottom soil layer
            Real den = (Zc(l, j) - Zc(l, j - 1)) * 1000.0;
            Real dzq = zq[j] - zq[j - 1];
            Real num = (smp(l, j) - smp(l, j - 1)) - dzq;

            const Real qin = -hk(l, j - 1) * num / den;
            const Real dqidw0 =
                -(-hk(l, j - 1) * dsmpdw[j - 1] + num * dhkdw[j - 1]) / den;
            const Real dqidw1 =
                -(hk(l, j - 1) * dsmpdw[j] + num * dhkdw[j - 1]) / den;

            // Compute aquifer layer properties
            const Real vwc_zwt = watsat(l, j);
            const Real vol_liq = Kokkos::fmax(h2osoi_liq(l, j), 1.0e-6) /
                                 (dz(l, j) * SHR_CONST_RHOWATER);
            const Real s_aquifer =
                Kokkos::fmax(0.5 * (vwc_zwt + vol_liq) / watsat(l, j), 0.01);
            const Real smp_aquifer = ComputeMatricPotential(
                s_aquifer, sucsat(l, j), bsw(l, j), SMP_MIN);
            const Real dsmpdw_aquifer = ComputeMatricPotentialDerivative(
                smp_aquifer, s_aquifer * watsat(l, j), bsw(l, j));

            const Real zmm_aquifer = 0.5 * (zwtmm + Zc(l, j) * 1000.0);
            den = (zmm_aquifer - Zc(l, j) * 1000.0);
            dzq = zq[j + 1] - zq[j];
            num = (smp_aquifer - smp(l, j)) - dzq;

            const Real qout = -hk(l, j) * num / den;
            const Real dqodw1 = -(-hk(l, j) * dsmpdw[j] + num * dhkdw[j]) / den;
            const Real dqodw2 =
                -(hk(l, j) * dsmpdw_aquifer + num * dhkdw[j]) / den;

            rmx(l, j) = qin - qout - qflx_tran(l, j);
            amx(l, j) = -dqidw0;
            bmx(l, j) = dz(l, j) * 1000.0 / dtime - dqidw1 + dqodw1;
            cmx(l, j) = dqodw2;

            // Aquifer layer
            const Real qin_aquifer = qout;
            const Real dqidw0_aquifer =
                -(-hk(l, j) * dsmpdw[j] + num * dhkdw[j]) / den;
            const Real dqidw1_aquifer =
                -(hk(l, j) * dsmpdw_aquifer + num * dhkdw[j]) / den;

            const Real qout_aquifer = 0.0; // Zero-flow bottom boundary
            const Real dqodw1_aquifer = 0.0;

            const Real dzmm_aquifer = zwtmm - Zc(l, j) * 1000.0;
            rmx(l, j + 1) = qin_aquifer - qout_aquifer;
            amx(l, j + 1) = -dqidw0_aquifer;
            bmx(l, j + 1) =
                dzmm_aquifer / dtime - dqidw1_aquifer + dqodw1_aquifer;
            cmx(l, j + 1) = 0.0;
          }
        }
      });
  Kokkos::fence();
}

// ============================================================================
// Phase 3.3: Solve Tridiagonal System
// ============================================================================

void SolveHydrologyTridiagonal(UrbanType urban) {
  const int nlandunits = urban->numLandunits;
  const int nlevbed = NUM_SOIL_LAYERS;

  auto jwt = urban->perviousRoad.Jwt;
  auto amx = urban->perviousRoad.soil.Amx;
  auto bmx = urban->perviousRoad.soil.Bmx;
  auto cmx = urban->perviousRoad.soil.Cmx;
  auto rmx = urban->perviousRoad.soil.Rmx;
  auto dwat = urban->perviousRoad.soil.Dwat;

  Kokkos::parallel_for(
      "SolveHydrologyTridiagonal", nlandunits, KOKKOS_LAMBDA(const int l) {
        // Determine number of active layers (including aquifer if needed)
        const int nlayers = NUM_SOIL_LAYERS;

        // Extract local arrays for this landunit
        Real a_local[NUM_SOIL_LAYERS];
        Real b_local[NUM_SOIL_LAYERS];
        Real c_local[NUM_SOIL_LAYERS];
        Real r_local[NUM_SOIL_LAYERS];
        Real x_local[NUM_SOIL_LAYERS];

        for (int j = 0; j < nlayers; ++j) {
          a_local[j] = amx(l, j);
          b_local[j] = bmx(l, j);
          c_local[j] = cmx(l, j);
          r_local[j] = rmx(l, j);
        }

        // Solve tridiagonal system using Thomas algorithm
        // Solves: a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = r[i]
        Real cp[NUM_SOIL_LAYERS]; // Modified upper diagonal
        Real rp[NUM_SOIL_LAYERS]; // Modified right-hand side

        // Forward elimination
        cp[0] = c_local[0] / b_local[0];
        rp[0] = r_local[0] / b_local[0];

        for (int i = 1; i < nlayers; i++) {
          const Real denom = b_local[i] - a_local[i] * cp[i - 1];
          cp[i] = c_local[i] / denom;
          rp[i] = (r_local[i] - a_local[i] * rp[i - 1]) / denom;
        }

        // Back substitution
        x_local[nlayers - 1] = rp[nlayers - 1];
        for (int i = nlayers - 2; i >= 0; i--) {
          x_local[i] = rp[i] - cp[i] * x_local[i + 1];
        }

        // Store solution
        for (int j = 0; j < nlayers; ++j) {
          dwat(l, j) = x_local[j];
        }
      });
  Kokkos::fence();
}

// ============================================================================
// Phase 3.4: Update State
// ============================================================================

void UpdateSoilWater(UrbanType urban, Real dtime) {
  const int nlandunits = urban->numLandunits;
  const int nlevbed = NUM_SOIL_LAYERS;

  auto h2osoi_liq = urban->perviousRoad.H2OSoiLiq;
  auto h2osoi_ice = urban->perviousRoad.H2OSoiIce;
  auto h2osoi_vol = urban->perviousRoad.H2OSoiVol;
  auto qcharge = urban->perviousRoad.Qcharge;
  auto qflx_deficit = urban->perviousRoad.QflxDeficit;
  auto jwt = urban->perviousRoad.Jwt;
  auto zwt = urban->perviousRoad.Zwt;
  auto hk = urban->perviousRoad.Hk;
  auto smp = urban->perviousRoad.Smp;

  auto watsat = urban->perviousRoad.soil.WatSat;
  auto hksat = urban->perviousRoad.soil.HkSat;
  auto bsw = urban->perviousRoad.soil.Bsw;
  auto sucsat = urban->perviousRoad.soil.SucSat;
  auto dz = urban->perviousRoad.Dz;
  auto Zc = urban->perviousRoad.Zc;
  auto zi = urban->perviousRoad.Zi;
  auto dwat = urban->perviousRoad.soil.Dwat;

  Kokkos::parallel_for(
      "UpdateSoilWater", nlandunits, KOKKOS_LAMBDA(const int l) {
        const int jwt_l = jwt(l);

        // Update liquid water content for soil layers
        for (int j = 0; j < nlevbed; ++j) {
          h2osoi_liq(l, j) += dwat(l, j) * dz(l, j) * 1000.0;
        }

        // Compute aquifer recharge rate
        if (jwt_l < nlevbed - 1) {
          // Water table is in soil column
          const Real wh_zwt = 0.0; // At water table: smp = -sucsat and zq =
                                   // -sucsat, so wh_zwt = 0

          // Compute hydraulic conductivity at jwt layer
          const Real s_node = Kokkos::fmax(
              h2osoi_vol(l, jwt_l + 1) / watsat(l, jwt_l + 1), 0.01);
          const Real s1 = Kokkos::fmin(1.0, s_node);

          // Compute ice impedance
          const Real vol_ice = Kokkos::fmin(
              watsat(l, jwt_l + 1),
              h2osoi_ice(l, jwt_l + 1) / (dz(l, jwt_l + 1) * SHR_CONST_RHOICE));
          const Real icefrac =
              Kokkos::fmin(1.0, vol_ice / watsat(l, jwt_l + 1));
          const Real imped = ComputeIceImpedance(icefrac);

          const Real ka = imped * hksat(l, jwt_l + 1) *
                          Kokkos::pow(s1, 2.0 * bsw(l, jwt_l + 1) + 3.0);

          // Compute equilibrium matric potential at jwt layer
          const Real z_top = (jwt_l == 0) ? 0.0 : zi(l, jwt_l - 1) * 1000.0;
          const Real z_bot = zi(l, jwt_l) * 1000.0;
          const Real zwtmm = zwt(l) * 1000.0;

          const Real vol_eq = ComputeEquilibriumWaterContent(
              zwtmm, z_top, z_bot, watsat(l, jwt_l), sucsat(l, jwt_l),
              bsw(l, jwt_l));

          const Real zq = ComputeEquilibriumMatricPotential(
              vol_eq, watsat(l, jwt_l), sucsat(l, jwt_l), bsw(l, jwt_l),
              SMP_MIN);

          const int jwt_max = Kokkos::fmax(0, jwt_l);
          const Real smp1 = Kokkos::fmax(SMP_MIN, smp(l, jwt_max));
          const Real wh = smp1 - zq;

          // Compute recharge rate
          if (jwt_l == 0) {
            qcharge(l) = -ka * (wh_zwt - wh) / ((zwt(l) + 1.e-3) * 1000.0);
          } else {
            qcharge(l) =
                -ka * (wh_zwt - wh) / ((zwt(l) - Zc(l, jwt_l)) * 1000.0 * 2.0);
          }

          // Limit qcharge for numerical stability
          qcharge(l) = Kokkos::fmax(-10.0 / dtime, qcharge(l));
          qcharge(l) = Kokkos::fmin(10.0 / dtime, qcharge(l));

        } else {
          // Water table below soil column - compute from aquifer layer
          qcharge(l) = dwat(l, nlevbed) *
                       (zwt(l) * 1000.0 - Zc(l, nlevbed - 1) * 1000.0) / dtime;
        }

        // Check for negative water content and compute deficit
        qflx_deficit(l) = 0.0;
        for (int j = 0; j < nlevbed; ++j) {
          if (h2osoi_liq(l, j) < 0.0) {
            qflx_deficit(l) -= h2osoi_liq(l, j) / dtime;
            h2osoi_liq(l, j) = 0.0;
          }
        }

        // Update volumetric water content
        for (int j = 0; j < nlevbed; ++j) {
          const Real vol_liq =
              h2osoi_liq(l, j) / (dz(l, j) * SHR_CONST_RHOWATER);
          const Real vol_ice = Kokkos::fmin(
              watsat(l, j), h2osoi_ice(l, j) / (dz(l, j) * SHR_CONST_RHOICE));
          h2osoi_vol(l, j) = vol_liq + vol_ice;
        }
      });
  Kokkos::fence();
}

// ============================================================================
// Internal wrapper for time-stepping integration
// ============================================================================

void ComputeHydrology(URBANXX::_p_UrbanType &urban) {
  // Use 30-minute timestep (1800 seconds) - typical for land surface models
  const Real dtime = 1800.0; // seconds

  ComputeHydraulicProperties(static_cast<UrbanType>(&urban));
  SetupHydrologyTridiagonal(static_cast<UrbanType>(&urban), dtime);
  SolveHydrologyTridiagonal(static_cast<UrbanType>(&urban));
  UpdateSoilWater(static_cast<UrbanType>(&urban), dtime);
}

} // namespace URBANXX

// ============================================================================
// Public C API
// ============================================================================

extern "C" {

void UrbanComputeHydrology(UrbanType urban, Real dtime,
                           UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    ComputeHydraulicProperties(urban);
    SetupHydrologyTridiagonal(urban, dtime);
    SolveHydrologyTridiagonal(urban);
    UpdateSoilWater(urban, dtime);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
