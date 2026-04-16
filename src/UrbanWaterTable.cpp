#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"
#include <Kokkos_Core.hpp>

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

using namespace URBANXX;

// ============================================================================
// UrbanComputeWaterTable
//
// Ports ELM's WaterTable subroutine (SoilHydrologyMod.F90:746-1082)
// for the pervious road column (urbpoi, use_vsfm=false, origflag=0,
// snl=0, qflx_grnd_irrig=0, zengdecker_2009_with_var_soil_thick=false).
//
// Inputs (already set via setters or computed by UrbanComputeHydrology):
//   Qcharge, H2OSoiLiq, H2OSoiIce, H2OSoiVol, Zwt,
//   FracH2osfc, QflxDewGrnd, QflxDewSnow, QflxSubSnow
//   soil.WatSat, soil.SucSat, soil.Bsw, Dz, Zi, Zc, Temperature
//   Wa (in/out: aquifer water storage)
//
// Outputs (written into views):
//   Zwt, Wa, ZwtPerched, FrostTable, QflxSubSnow,
//   QflxDrain (=0, use_vsfm=false), QflxRsubSat (=0),
//   H2OSoiLiq[0], H2OSoiIce[0] (updated by dew), Jwt
// ============================================================================

extern "C" {

void UrbanComputeWaterTable(UrbanType urban, double dtime,
                            UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  const int nlandunits = urban->numLandunits;
  const int nlevbed = NUM_LAYERS_ABV_BEDROCK; // 10 — matches ELM nlevbed

  auto zwt = urban->perviousRoad.Zwt;
  auto wa = urban->perviousRoad.Wa;
  auto jwt_view = urban->perviousRoad.Jwt;
  auto qcharge = urban->perviousRoad.Qcharge;
  auto frost_table = urban->perviousRoad.FrostTable;
  auto zwt_perched = urban->perviousRoad.ZwtPerched;
  auto qflx_drain = urban->perviousRoad.QflxDrain;
  auto qflx_rsub_sat = urban->perviousRoad.QflxRsubSat;
  auto qflx_sub_snow = urban->perviousRoad.QflxSubSnow;
  auto qflx_dew_grnd = urban->perviousRoad.QflxDewGrnd;
  auto qflx_dew_snow = urban->perviousRoad.QflxDewSnow;
  auto frac_h2osfc = urban->perviousRoad.FracH2osfc;
  auto h2osoi_liq = urban->perviousRoad.H2OSoiLiq;
  auto h2osoi_ice = urban->perviousRoad.H2OSoiIce;
  auto h2osoi_vol = urban->perviousRoad.H2OSoiVol;
  auto temperature = urban->perviousRoad.Temperature;

  auto watsat = urban->perviousRoad.soil.WatSat;
  auto sucsat = urban->perviousRoad.soil.SucSat;
  auto bsw = urban->perviousRoad.soil.Bsw;
  auto dz = urban->perviousRoad.Dz;
  auto zi = urban->perviousRoad.Zi;
  auto zc = urban->perviousRoad.Zc;

  Kokkos::parallel_for(
      "UrbanComputeWaterTable", nlandunits, KOKKOS_LAMBDA(const int l) {
        // ----------------------------------------------------------------
        // Initialize outputs (use_vsfm = false)
        // ----------------------------------------------------------------
        qflx_drain(l) = 0.0;
        qflx_rsub_sat(l) = 0.0;

        // ----------------------------------------------------------------
        // Compute jwt: index of the first unsaturated layer above water
        // table.  Mirrors ELM WaterTable lines 847-860.
        // Convention: jwt == nlevbed → water table below all soil layers.
        //             jwt == k       → layers 0..k-1 are above water table
        //             (0-based, so matches ELM jwt=j-1 where j is 1-based).
        // ----------------------------------------------------------------
        int jwt_l = nlevbed; // default: water table below all soil
        for (int k = 0; k < nlevbed; ++k) {
          // zi(l, k+1) = bottom interface of layer k = ELM zi(c, k+1)
          if (zwt(l) <= zi(l, k + 1)) {
            jwt_l = k;
            break;
          }
        }

        // ----------------------------------------------------------------
        // QCHARGE section (ELM lines 862-960).
        // Save qcharge; it is read-only here.
        // ----------------------------------------------------------------
        const Real qcharge_val = qcharge(l);

        // Specific yield at bottom bedrock layer (0-based index nlevbed-1)
        Real rous =
            watsat(l, nlevbed - 1) *
            (1.0 - Kokkos::pow(1.0 + 1.0e3 * zwt(l) / sucsat(l, nlevbed - 1),
                               -1.0 / bsw(l, nlevbed - 1)));
        rous = Kokkos::fmax(rous, 0.02);

        if (jwt_l == nlevbed) {
          // Water table below soil column
          wa(l) += qcharge_val * dtime;
          zwt(l) -= (qcharge_val * dtime) / 1000.0 / rous;
        } else {
          // Water table within soil layers
          Real qcharge_tot = qcharge_val * dtime;

          if (qcharge_tot > 0.0) {
            // Rising water table: loop from jwt_l down to 0
            // ELM: do j = jwt+1, 1, -1  → C++: k = jwt_l down to 0
            for (int k = jwt_l; k >= 0; --k) {
              Real s_y = watsat(l, k) *
                         (1.0 - Kokkos::pow(1.0 + 1.0e3 * zwt(l) / sucsat(l, k),
                                            -1.0 / bsw(l, k)));
              s_y = Kokkos::fmax(s_y, 0.02);

              // Top interface of layer k is zi(l, k)
              Real qcharge_layer =
                  Kokkos::fmin(qcharge_tot, s_y * (zwt(l) - zi(l, k)) * 1.0e3);
              qcharge_layer = Kokkos::fmax(qcharge_layer, 0.0);

              if (s_y > 0.0)
                zwt(l) -= qcharge_layer / s_y / 1000.0;

              qcharge_tot -= qcharge_layer;
              if (qcharge_tot <= 0.0)
                break;
            }
          } else {
            // Deepening water table: loop from jwt_l+1 to nlevbed-1
            // ELM: do j = jwt+1, nlevbed  → C++: k = jwt_l+1 to nlevbed-1
            bool done = false;
            for (int k = jwt_l + 1; k < nlevbed; ++k) {
              Real s_y = watsat(l, k) *
                         (1.0 - Kokkos::pow(1.0 + 1.0e3 * zwt(l) / sucsat(l, k),
                                            -1.0 / bsw(l, k)));
              s_y = Kokkos::fmax(s_y, 0.02);

              // Bottom interface of layer k is zi(l, k+1)
              Real qcharge_layer = Kokkos::fmax(
                  qcharge_tot, -(s_y * (zi(l, k + 1) - zwt(l)) * 1.0e3));
              qcharge_layer = Kokkos::fmin(qcharge_layer, 0.0);
              qcharge_tot -= qcharge_layer;

              if (qcharge_tot >= 0.0) {
                zwt(l) -= qcharge_layer / s_y / 1000.0;
                done = true;
                break;
              } else {
                zwt(l) = zi(l, k + 1);
              }
            }
            if (!done && qcharge_tot > 0.0)
              zwt(l) -= qcharge_tot / 1000.0 / rous;
          }

          // Recompute jwt after moving water table
          jwt_l = nlevbed;
          for (int k = 0; k < nlevbed; ++k) {
            if (zwt(l) <= zi(l, k + 1)) {
              jwt_l = k;
              break;
            }
          }
        }
        jwt_view(l) = jwt_l;

        // ----------------------------------------------------------------
        // Frost table and perched water table (ELM lines 970-1040).
        // origflag == 0, so we may compute perched WT.
        // Temperature(l, k) is 0-based = t_soisno(c, k+1) in ELM.
        // ----------------------------------------------------------------
        const Real tfrz = SHR_CONST_TKFRZ;

        // Find k_frz: first frozen layer with unfrozen layer above (0-based)
        int k_frz = nlevbed - 1; // default: all unfrozen
        if (temperature(l, 0) <= tfrz) {
          k_frz = 0; // top layer is frozen
        } else {
          for (int k = 1; k < nlevbed; ++k) {
            if (temperature(l, k - 1) > tfrz && temperature(l, k) <= tfrz) {
              k_frz = k;
              break;
            }
          }
        }

        frost_table(l) = zc(l, k_frz);
        zwt_perched(l) = frost_table(l); // initialize perched WT to frost table

        // If zwt < frost_table AND frost layer is truly frozen: skip perched WT
        if (zwt(l) < frost_table(l) && temperature(l, k_frz) <= tfrz) {
          // zwt_perched already = frost_table; no further computation
        } else {
          // Locate perched water table from bottom up starting at frost table
          const Real sat_lev = 0.9;
          int k_perch =
              0; // ELM initializes to k_perch=1 (1-based) = 0 (0-based)

          for (int k = k_frz; k >= 0; --k) {
            const Real vol =
                h2osoi_liq(l, k) / (dz(l, k) * SHR_CONST_RHOWATER) +
                h2osoi_ice(l, k) / (dz(l, k) * SHR_CONST_RHOICE);
            if (vol / watsat(l, k) <= sat_lev) {
              k_perch = k;
              break;
            }
          }

          // If frost table layer is actually thawed, no perched WT
          if (temperature(l, k_frz) > tfrz)
            k_perch = k_frz;

          // If perched water table exists between k_perch and k_perch+1
          if (k_frz > k_perch) {
            const Real s1 =
                (h2osoi_liq(l, k_perch) /
                     (dz(l, k_perch) * SHR_CONST_RHOWATER) +
                 h2osoi_ice(l, k_perch) / (dz(l, k_perch) * SHR_CONST_RHOICE)) /
                watsat(l, k_perch);
            const Real s2 = (h2osoi_liq(l, k_perch + 1) /
                                 (dz(l, k_perch + 1) * SHR_CONST_RHOWATER) +
                             h2osoi_ice(l, k_perch + 1) /
                                 (dz(l, k_perch + 1) * SHR_CONST_RHOICE)) /
                            watsat(l, k_perch + 1);

            const Real m = (zc(l, k_perch + 1) - zc(l, k_perch)) / (s2 - s1);
            const Real b = zc(l, k_perch + 1) - m * s2;
            zwt_perched(l) = Kokkos::fmax(0.0, m * sat_lev + b);
          }
        }

        // ----------------------------------------------------------------
        // Dew update (ELM lines 1042-1060).
        // snl(c) == 0 so snl+1 = 1 >= 1: apply dew/sublimation.
        // Updates layer-1 (0-based layer 0) liquid and ice.
        // ----------------------------------------------------------------
        h2osoi_liq(l, 0) += (1.0 - frac_h2osfc(l)) * qflx_dew_grnd(l) * dtime;
        h2osoi_ice(l, 0) += (1.0 - frac_h2osfc(l)) * qflx_dew_snow(l) * dtime;

        if (qflx_sub_snow(l) * dtime > h2osoi_ice(l, 0)) {
          qflx_sub_snow(l) = h2osoi_ice(l, 0) / dtime;
          h2osoi_ice(l, 0) = 0.0;
        } else {
          h2osoi_ice(l, 0) -= (1.0 - frac_h2osfc(l)) * qflx_sub_snow(l) * dtime;
        }

        // ----------------------------------------------------------------
        // Clamp water table depth (ELM lines 1068+)
        // ----------------------------------------------------------------
        zwt(l) = Kokkos::fmax(0.0, Kokkos::fmin(80.0, zwt(l)));
      });
  Kokkos::fence();

  *status = URBAN_SUCCESS;
}

// ============================================================================
// UrbanComputeDewCondensationRoofImperviousRoad
//
// Ports ELM's dew condensation block inside WaterTable
// (SoilHydrologyMod.F90:1062-1078) for the roof and impervious road columns.
//
// Assumption: snl(c) == 0 for both surfaces (no snow layers), so the
// condition snl+1 >= 1 is always satisfied.  Unlike the pervious road path
// there is no (1 - frac_h2osfc) factor.
//
// Inputs (set via setters before calling this function):
//   roof.TopH2OSoiLiq, roof.TopH2OSoiIce
//   roof.QflxDewGrnd, roof.QflxDewSnow, roof.QflxSubSnow
//   imperviousRoad.TopH2OSoiLiq, imperviousRoad.TopH2OSoiIce
//   imperviousRoad.QflxDewGrnd, imperviousRoad.QflxDewSnow,
//   imperviousRoad.QflxSubSnow
//
// Outputs (written into views):
//   roof.TopH2OSoiLiq, roof.TopH2OSoiIce, roof.QflxSubSnow
//   imperviousRoad.TopH2OSoiLiq, imperviousRoad.TopH2OSoiIce,
//   imperviousRoad.QflxSubSnow
// ============================================================================

void UrbanComputeDewCondensationRoofImperviousRoad(UrbanType urban,
                                                   double dtime,
                                                   UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  const int nlandunits = urban->numLandunits;

  // --- Roof ---
  auto roof_liq = urban->roof.TopH2OSoiLiq;
  auto roof_ice = urban->roof.TopH2OSoiIce;
  auto roof_dew_grnd = urban->roof.QflxDewGrnd;
  auto roof_dew_snow = urban->roof.QflxDewSnow;
  auto roof_sub_snow = urban->roof.QflxSubSnow;

  Kokkos::parallel_for(
      "UrbanComputeDewCondensationRoof", nlandunits,
      KOKKOS_LAMBDA(const int l) {
        roof_liq(l) += roof_dew_grnd(l) * dtime;
        roof_ice(l) += roof_dew_snow(l) * dtime;
        if (roof_sub_snow(l) * dtime > roof_ice(l)) {
          roof_sub_snow(l) = roof_ice(l) / dtime;
          roof_ice(l) = 0.0;
        } else {
          roof_ice(l) -= roof_sub_snow(l) * dtime;
        }
      });

  // --- Impervious Road ---
  auto imperv_liq = urban->imperviousRoad.TopH2OSoiLiq;
  auto imperv_ice = urban->imperviousRoad.TopH2OSoiIce;
  auto imperv_dew_grnd = urban->imperviousRoad.QflxDewGrnd;
  auto imperv_dew_snow = urban->imperviousRoad.QflxDewSnow;
  auto imperv_sub_snow = urban->imperviousRoad.QflxSubSnow;

  Kokkos::parallel_for(
      "UrbanComputeDewCondensationImperviousRoad", nlandunits,
      KOKKOS_LAMBDA(const int l) {
        imperv_liq(l) += imperv_dew_grnd(l) * dtime;
        imperv_ice(l) += imperv_dew_snow(l) * dtime;
        if (imperv_sub_snow(l) * dtime > imperv_ice(l)) {
          imperv_sub_snow(l) = imperv_ice(l) / dtime;
          imperv_ice(l) = 0.0;
        } else {
          imperv_ice(l) -= imperv_sub_snow(l) * dtime;
        }
      });

  Kokkos::fence();

  *status = URBAN_SUCCESS;
}

} // extern "C"
