#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include "private/UrbanSetterHelpers.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"
#include <Kokkos_Core.hpp>

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

using namespace URBANXX;

// ============================================================================
// UrbanComputeDrainage
//
// Ports ELM's Drainage subroutine (SoilHydrologyMod.F90:1085-1773)
// for the pervious road column only (lun_pp%urbpoi=1, origflag=0,
// use_vsfm=false, use_vichydro=false, snl=0,
// zengdecker_2009_with_var_soil_thick=false).
//
// Inputs (set via setters before this call):
//   H2OSoiLiq, H2OSoiIce, soil.WatSat, soil.Bsw, soil.SucSat, soil.HkSat
//   EffPorosity, Temperature, Dz, Zi, Zc, Zwt, Wa, HkDepth, TopoSlope
//   Params: RsubTopGlobalMax, Pondmx, Watmin, EIce
//
// Outputs (written into views):
//   QflxDrain, QflxRsubSat
//   H2OSoiLiq (updated in-place), Zwt, Wa (updated in-place, not copied to ELM)
// ============================================================================

extern "C" {

// ============================================================================
// Drainage init-time constant setters
// ============================================================================

void UrbanSetRsubTopGlobalMax(UrbanType urban, double value,
                              UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }
  urban->urbanParams.RsubTopGlobalMax = value;
  *status = URBAN_SUCCESS;
}

void UrbanSetPondmax(UrbanType urban, double value, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }
  urban->urbanParams.Pondmx = value;
  *status = URBAN_SUCCESS;
}

void UrbanSetWatmin(UrbanType urban, double value, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }
  urban->urbanParams.Watmin = value;
  *status = URBAN_SUCCESS;
}

void UrbanSetEice(UrbanType urban, double value, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }
  urban->urbanParams.EIce = value;
  *status = URBAN_SUCCESS;
}

// ============================================================================
// Drainage timestep setter functions
// ============================================================================

void UrbanSetHkDepthForPerviousRoad(UrbanType urban, const double *values,
                                    int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.HkDepth, values, length, status);
}

void UrbanSetTopoSlopeForPerviousRoad(UrbanType urban, const double *values,
                                      int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.TopoSlope, values, length, status);
}

// ============================================================================
// UrbanComputeDrainage
// ============================================================================

void UrbanComputeDrainage(UrbanType urban, double dtime,
                          UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  const int nlandunits = urban->numLandunits;
  const int nlevbed = NUM_LAYERS_ABV_BEDROCK; // 10 — matches ELM nlevbed

  // Drainage constants (set at initialization)
  const Real rsub_top_globalmax =
      static_cast<Real>(urban->urbanParams.RsubTopGlobalMax);
  const Real pondmx = static_cast<Real>(urban->urbanParams.Pondmx);
  const Real watmin = static_cast<Real>(urban->urbanParams.Watmin);
  const Real e_ice = static_cast<Real>(urban->urbanParams.EIce);

  // Physical constants
  constexpr Real tfrz = SHR_CONST_TKFRZ;
  constexpr Real denice = SHR_CONST_RHOICE;   // density of ice [kg/m^3]
  constexpr Real denh2o = SHR_CONST_RHOWATER; // density of water [kg/m^3]
  constexpr Real rpi = 3.14159265358979323846;

  auto zwt = urban->perviousRoad.Zwt;
  auto wa = urban->perviousRoad.Wa;
  auto qflx_drain = urban->perviousRoad.QflxDrain;
  auto qflx_rsub_sat = urban->perviousRoad.QflxRsubSat;
  auto h2osoi_liq = urban->perviousRoad.H2OSoiLiq;
  auto h2osoi_ice = urban->perviousRoad.H2OSoiIce;
  auto eff_porosity = urban->perviousRoad.EffPorosity;
  auto temperature = urban->perviousRoad.Temperature;
  auto hkdepth = urban->perviousRoad.HkDepth;
  auto topo_slope = urban->perviousRoad.TopoSlope;

  auto watsat = urban->perviousRoad.soil.WatSat;
  auto bsw = urban->perviousRoad.soil.Bsw;
  auto sucsat = urban->perviousRoad.soil.SucSat;
  auto hksat = urban->perviousRoad.soil.HkSat;
  auto dz = urban->perviousRoad.Dz;
  auto zi = urban->perviousRoad.Zi;
  auto zc = urban->perviousRoad.Zc;

  Kokkos::parallel_for(
      "UrbanComputeDrainage", nlandunits, KOKKOS_LAMBDA(const int l) {
        // ----------------------------------------------------------------
        // Initialize outputs
        // ----------------------------------------------------------------
        qflx_drain(l) = 0.0;
        qflx_rsub_sat(l) = 0.0;
        Real rsub_top = 0.0;

        // ----------------------------------------------------------------
        // Step 0: Compute dzmm[j] and icefrac[j] for all soil layers
        // ----------------------------------------------------------------
        Real dzmm[NUM_LAYERS_ABV_BEDROCK];
        Real icefrac[NUM_LAYERS_ABV_BEDROCK];
        for (int j = 0; j < nlevbed; ++j) {
          dzmm[j] = dz(l, j) * 1.0e3;
          Real vol_ice = Kokkos::fmin(watsat(l, j),
                                      h2osoi_ice(l, j) / (dz(l, j) * denice));
          icefrac[j] = Kokkos::fmin(1.0, vol_ice / watsat(l, j));
        }

        // ----------------------------------------------------------------
        // Step 1: Compute jwt — index of first unsaturated layer above WT
        // Convention: jwt_l == nlevbed means water table below all soil
        //             jwt_l == k means water table is within layer k
        //             (0-based, matches ELM's jwt = j-1 for 1-based j)
        // ----------------------------------------------------------------
        int jwt_l = nlevbed;
        for (int k = 0; k < nlevbed; ++k) {
          if (zwt(l) <= zi(l, k + 1)) {
            jwt_l = k;
            break;
          }
        }

        // ----------------------------------------------------------------
        // Step 2: q_perch_max
        // ----------------------------------------------------------------
        const Real q_perch_max =
            1.0e-5 * Kokkos::sin(topo_slope(l) * (rpi / 180.0));

        // ----------------------------------------------------------------
        // Step 3: Find k_frz (0-based) — first frozen layer with unfrozen
        //         layer above it. frost_table = z(l, k_frz).
        //         Initialize zwt_perched = frost_table.
        // ----------------------------------------------------------------
        int k_frz = nlevbed - 1; // default: all unfrozen
        if (temperature(l, 0) <= tfrz) {
          k_frz = 0;
        } else {
          for (int k = 1; k < nlevbed; ++k) {
            if (temperature(l, k - 1) > tfrz && temperature(l, k) <= tfrz) {
              k_frz = k;
              break;
            }
          }
        }
        const Real frost_table_val = zc(l, k_frz);
        Real zwt_perched_val = frost_table_val;
        // (zwt_perched and frost_table are updated below but not output)

        // ----------------------------------------------------------------
        // Step 4 (Branch A): water table above frost table AND frost
        //         layer is frozen (origflag == 0 path in ELM).
        //         Compute perched drainage and remove from soil layers
        //         jwt_l..k_frz (0-indexed), then recompute jwt_l.
        // ----------------------------------------------------------------
        if (zwt(l) < frost_table_val && temperature(l, k_frz) <= tfrz) {
          // Compute impedance-weighted hydraulic conductivity for perched
          // saturated layers [ELM lines 1322-1332]
          Real wtsub = 0.0;
          Real q_perch = 0.0;
          // ELM: do k = jwt(c)+1, k_frz (1-indexed) = jwt_l..k_frz (0-indexed)
          for (int k = jwt_l; k <= k_frz; ++k) {
            int k_next = (k < nlevbed - 1) ? k + 1 : k;
            Real imped = Kokkos::pow(
                10.0, -e_ice * (0.5 * (icefrac[k] + icefrac[k_next])));
            q_perch += imped * hksat(l, k) * dzmm[k];
            wtsub += dzmm[k];
          }
          if (wtsub > 0.0)
            q_perch /= wtsub;

          Real qflx_drain_perched =
              q_perch_max * q_perch * (frost_table_val - zwt(l));

          // Remove drainage from perched layers [ELM lines 1335-1356]
          Real rsub_top_tot = -qflx_drain_perched * dtime;
          for (int k = jwt_l; k <= k_frz; ++k) {
            Real rsub_top_layer =
                Kokkos::fmax(rsub_top_tot, -(h2osoi_liq(l, k) - watmin));
            rsub_top_layer = Kokkos::fmin(rsub_top_layer, 0.0);
            rsub_top_tot -= rsub_top_layer;
            h2osoi_liq(l, k) += rsub_top_layer;
            if (rsub_top_tot >= 0.0) {
              zwt(l) -= rsub_top_layer / eff_porosity(l, k) / 1000.0;
              break;
            } else {
              zwt(l) = zi(l, k + 1); // bottom interface of layer k (0-based)
            }
          }
          // Correct qflx_drain_perched for residual
          qflx_drain_perched += rsub_top_tot / dtime;

          // Recompute jwt_l [ELM lines 1359-1368]
          jwt_l = nlevbed;
          for (int k = 0; k < nlevbed; ++k) {
            if (zwt(l) <= zi(l, k + 1)) {
              jwt_l = k;
              break;
            }
          }
          // rsub_top stays 0.0; qflx_drain_perched is internal only

        } else {
          // ----------------------------------------------------------------
          // Step 5 (Branch B): water table at or below frost table.
          //   a. Locate perched water table at sat_lev=0.9 from bottom up.
          //   b. Compute topographic baseflow (rsub_top).
          //      Remove water from soil layers jwt_l..nlevbed-1 (0-indexed).
          // ----------------------------------------------------------------

          // 5a: Find perched WT [ELM lines 1375-1421]
          // (Used internally — zwt_perched not output from this kernel)
          const Real sat_lev = 0.9;
          int k_perch = 0;
          for (int k = k_frz; k >= 0; --k) {
            Real h2osoi_vol_k = h2osoi_liq(l, k) / (dz(l, k) * denh2o) +
                                h2osoi_ice(l, k) / (dz(l, k) * denice);
            if (h2osoi_vol_k / watsat(l, k) <= sat_lev) {
              k_perch = k;
              break;
            }
          }
          // If frost layer is unfrozen, no perched table
          if (temperature(l, k_frz) > tfrz)
            k_perch = k_frz;

          if (k_frz > k_perch) {
            // Interpolate to find perched WT depth
            Real s1 = (h2osoi_liq(l, k_perch) / (dz(l, k_perch) * denh2o) +
                       h2osoi_ice(l, k_perch) / (dz(l, k_perch) * denice)) /
                      watsat(l, k_perch);
            Real s2 =
                (h2osoi_liq(l, k_perch + 1) / (dz(l, k_perch + 1) * denh2o) +
                 h2osoi_ice(l, k_perch + 1) / (dz(l, k_perch + 1) * denice)) /
                watsat(l, k_perch + 1);
            Real m_slope = (zc(l, k_perch + 1) - zc(l, k_perch)) / (s2 - s1);
            Real b_int = zc(l, k_perch + 1) - m_slope * s2;
            zwt_perched_val = Kokkos::fmax(0.0, m_slope * sat_lev + b_int);

            // Compute perched drainage and remove from layers
            Real wtsub = 0.0;
            Real q_perch = 0.0;
            for (int k = k_perch; k <= k_frz; ++k) {
              int k_next = (k < nlevbed - 1) ? k + 1 : k;
              Real imped = Kokkos::pow(
                  10.0, -e_ice * (0.5 * (icefrac[k] + icefrac[k_next])));
              q_perch += imped * hksat(l, k) * dzmm[k];
              wtsub += dzmm[k];
            }
            if (wtsub > 0.0)
              q_perch /= wtsub;

            Real qflx_drain_perched =
                q_perch_max * q_perch * (frost_table_val - zwt_perched_val);

            Real rsub_top_tot = -qflx_drain_perched * dtime;
            for (int k = k_perch + 1; k <= k_frz; ++k) {
              Real rsub_top_layer =
                  Kokkos::fmax(rsub_top_tot, -(h2osoi_liq(l, k) - watmin));
              rsub_top_layer = Kokkos::fmin(rsub_top_layer, 0.0);
              rsub_top_tot -= rsub_top_layer;
              h2osoi_liq(l, k) += rsub_top_layer;
              if (rsub_top_tot >= 0.0) {
                zwt_perched_val -= rsub_top_layer / eff_porosity(l, k) / 1000.0;
                break;
              } else {
                zwt_perched_val = zi(l, k + 1);
              }
            }
            qflx_drain_perched += rsub_top_tot / dtime;
          }

          // 5b: Topographic baseflow [ELM lines 1437-1498]
          const Real fff = 1.0 / hkdepth(l);

          // Compute ice impedance sum over max(jwt_l-1,0)..nlevbed-1
          // [ELM: do j = max(jwt(c),1), nlevbed (1-indexed)]
          Real dzsum = 0.0;
          Real icefracsum = 0.0;
          const int k_start = (jwt_l > 0) ? jwt_l - 1 : 0;
          for (int k = k_start; k < nlevbed; ++k) {
            dzsum += dzmm[k];
            icefracsum += icefrac[k] * dzmm[k];
          }

          // use_vichydro=false, origflag=0 path:
          Real imped = 0.0;
          if (dzsum > 0.0) {
            imped = Kokkos::pow(10.0, -e_ice * (icefracsum / dzsum));
          }
          Real rsub_top_max =
              Kokkos::fmin(10.0 * Kokkos::sin((rpi / 180.0) * topo_slope(l)),
                           rsub_top_globalmax);
          rsub_top = imped * rsub_top_max * Kokkos::exp(-fff * zwt(l));

          // Aquifer yield for potential water table update
          Real rous =
              watsat(l, nlevbed - 1) *
              (1.0 - Kokkos::pow(1.0 + 1.0e3 * zwt(l) / sucsat(l, nlevbed - 1),
                                 -1.0 / bsw(l, nlevbed - 1)));
          rous = Kokkos::fmax(rous, 0.02);

          if (jwt_l == nlevbed) {
            // Water table below all soil layers [ELM lines 1500-1535]
            wa(l) -= rsub_top * dtime;
            zwt(l) += (rsub_top * dtime) / 1000.0 / rous;
            // Cap aquifer storage at 5000 mm; excess goes to bottom layer
            Real excess = Kokkos::fmax(0.0, wa(l) - 5000.0);
            h2osoi_liq(l, nlevbed - 1) += excess;
            wa(l) = Kokkos::fmin(wa(l), 5000.0);
          } else {
            // Water table within soil layers [ELM lines 1536-1601]
            Real rsub_top_tot = -rsub_top * dtime;
            // rsub_top_tot will be <= 0 (deepening water table)
            for (int j = jwt_l; j < nlevbed; ++j) {
              Real s_y = watsat(l, j) *
                         (1.0 - Kokkos::pow(1.0 + 1.0e3 * zwt(l) / sucsat(l, j),
                                            -1.0 / bsw(l, j)));
              s_y = Kokkos::fmax(s_y, 0.02);

              Real rsub_top_layer = Kokkos::fmax(
                  rsub_top_tot, -(s_y * (zi(l, j + 1) - zwt(l)) * 1.0e3));
              rsub_top_layer = Kokkos::fmin(rsub_top_layer, 0.0);
              h2osoi_liq(l, j) += rsub_top_layer;
              rsub_top_tot -= rsub_top_layer;

              if (rsub_top_tot >= 0.0) {
                zwt(l) -= rsub_top_layer / s_y / 1000.0;
                break;
              } else {
                zwt(l) = zi(l, j + 1);
              }
            }
            // Remove residual from aquifer / adjust zwt
            // (zengdecker_2009_with_var_soil_thick = false)
            zwt(l) -= rsub_top_tot / 1000.0 / rous;
            wa(l) += rsub_top_tot;

            // Recompute jwt_l
            jwt_l = nlevbed;
            for (int k = 0; k < nlevbed; ++k) {
              if (zwt(l) <= zi(l, k + 1)) {
                jwt_l = k;
                break;
              }
            }
          }

          // Clamp zwt
          zwt(l) = Kokkos::fmax(0.0, zwt(l));
          zwt(l) = Kokkos::fmin(80.0, zwt(l));

        } // end Branch A / Branch B

        // ----------------------------------------------------------------
        // Step 6: Saturation cascade — propagate excess water upward
        // [ELM lines 1628-1641]
        // ----------------------------------------------------------------
        for (int j = nlevbed - 1; j >= 1; --j) {
          Real xsi = Kokkos::fmax(
              h2osoi_liq(l, j) - eff_porosity(l, j) * dzmm[j], 0.0);
          h2osoi_liq(l, j) =
              Kokkos::fmin(eff_porosity(l, j) * dzmm[j], h2osoi_liq(l, j));
          h2osoi_liq(l, j - 1) += xsi;
        }

        // ----------------------------------------------------------------
        // Step 7: Compute xs1 and qflx_rsub_sat [ELM lines 1643-1679]
        // (urban path: lun_pp%urbpoi == 1)
        // ----------------------------------------------------------------
        Real xs1 =
            Kokkos::fmax(Kokkos::fmax(h2osoi_liq(l, 0) - watmin, 0.0) -
                             Kokkos::fmax(pondmx + watsat(l, 0) * dzmm[0] -
                                              h2osoi_ice(l, 0) - watmin,
                                          0.0),
                         0.0);
        qflx_rsub_sat(l) = xs1 / dtime;
        h2osoi_liq(l, 0) -= xs1;

        // ----------------------------------------------------------------
        // Step 8: watmin redistribution [ELM lines 1703-1760]
        // ----------------------------------------------------------------
        // Layers 0..nlevbed-2: get water from layer below if below watmin
        for (int j = 0; j < nlevbed - 1; ++j) {
          Real xs = 0.0;
          if (h2osoi_liq(l, j) < watmin) {
            xs = watmin - h2osoi_liq(l, j);
            if (j == jwt_l) {
              zwt(l) += xs / eff_porosity(l, j) / 1000.0;
            }
          }
          h2osoi_liq(l, j) += xs;
          h2osoi_liq(l, j + 1) -= xs;
        }

        // Bottom layer: search upward for water if below watmin
        int bot = nlevbed - 1;
        if (h2osoi_liq(l, bot) < watmin) {
          Real xs = watmin - h2osoi_liq(l, bot);
          for (int i = nlevbed - 2; i >= 0; --i) {
            Real avail = Kokkos::fmax(h2osoi_liq(l, i) - watmin - xs, 0.0);
            if (avail >= xs) {
              h2osoi_liq(l, bot) += xs;
              h2osoi_liq(l, i) -= xs;
              xs = 0.0;
              break;
            } else {
              h2osoi_liq(l, bot) += avail;
              h2osoi_liq(l, i) -= avail;
              xs -= avail;
            }
          }
          // Any remaining xs is taken from rsub drainage (water balance)
          rsub_top -= xs / dtime;
          h2osoi_liq(l, bot) += xs;
        }

        // ----------------------------------------------------------------
        // Step 9: Sub-surface runoff
        // ----------------------------------------------------------------
        qflx_drain(l) = qflx_rsub_sat(l) + rsub_top;
      });

  Kokkos::fence();
  *status = URBAN_SUCCESS;
}

} // extern "C"
