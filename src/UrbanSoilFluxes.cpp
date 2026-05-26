#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include "private/UrbanSurfaceTypeImpl.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"
#include <Kokkos_Core.hpp>
#include <cmath>

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

namespace URBANXX {

// ============================================================================
// ComputeSoilFluxes
//
// Ports the ELM SoilFluxes urban path (lun_pp%urbpoi = .true.,
// use_finetop_rad = .false.) to URBANxx.
//
// Assumptions / simplifications (see ELM_URBANxx_SoilFluxes_Plan.md):
//   - frac_h2osfc   = 0 for all surfaces
//   - egirat applied: topsoil evap clamped to (topLiq+topIce)/dtime
//   - dlrad         = 0 for all urban surfaces
//   - do_capsnow    = false  =>  qflx_snwcp update skipped
//   - eflx_wasteheat / heat_from_ac / traffic = 0 for now
//   - htvp = hvap if T > tfrz, else hvap + hfus
// ============================================================================
void ComputeSoilFluxes(URBANXX::_p_UrbanType &urban) {
  const int numLandunits = urban.numLandunits;

  // ---- SurfaceDataBase views for all five surfaces ----
  auto roof_TGrnd0 = urban.roof.TGrnd0;
  auto roof_Temp = urban.roof.EffectiveSurfTemp;
  auto roof_NetSw = urban.roof.NetShortRad;
  auto roof_NetLw = urban.roof.NetLongRad;
  auto roof_Emiss = urban.urbanParams.emissivity.Roof;
  auto roof_Cgrnds = urban.roof.Cgrnds;
  auto roof_Cgrndl = urban.roof.Cgrndl;
  auto roof_EflxShGrnd = urban.roof.EflxShGrnd;
  auto roof_QflxEvapSoil = urban.roof.QflxEvapSoil;
  auto roof_QflxTranEvap = urban.roof.QflxTranEvap;
  auto roof_EflxSoilGrnd = urban.roof.EflxSoilGrnd;
  auto roof_QflxEvapGrnd = urban.roof.QflxEvapGrnd;
  auto roof_QflxSubSnow = urban.roof.QflxSubSnow;
  auto roof_QflxDewSnow = urban.roof.QflxDewSnow;
  auto roof_QflxDewGrnd = urban.roof.QflxDewGrnd;
  auto roof_TopLiq = urban.roof.TopH2OSoiLiq;
  auto roof_TopIce = urban.roof.TopH2OSoiIce;

  auto imperv_TGrnd0 = urban.imperviousRoad.TGrnd0;
  auto imperv_Temp = urban.imperviousRoad.EffectiveSurfTemp;
  auto imperv_NetSw = urban.imperviousRoad.NetShortRad;
  auto imperv_NetLw = urban.imperviousRoad.NetLongRad;
  auto imperv_Emiss = urban.urbanParams.emissivity.ImperviousRoad;
  auto imperv_Cgrnds = urban.imperviousRoad.Cgrnds;
  auto imperv_Cgrndl = urban.imperviousRoad.Cgrndl;
  auto imperv_EflxShGrnd = urban.imperviousRoad.EflxShGrnd;
  auto imperv_QflxEvapSoil = urban.imperviousRoad.QflxEvapSoil;
  auto imperv_QflxTranEvap = urban.imperviousRoad.QflxTranEvap;
  auto imperv_EflxSoilGrnd = urban.imperviousRoad.EflxSoilGrnd;
  auto imperv_QflxEvapGrnd = urban.imperviousRoad.QflxEvapGrnd;
  auto imperv_QflxSubSnow = urban.imperviousRoad.QflxSubSnow;
  auto imperv_QflxDewSnow = urban.imperviousRoad.QflxDewSnow;
  auto imperv_QflxDewGrnd = urban.imperviousRoad.QflxDewGrnd;
  auto imperv_TopLiq = urban.imperviousRoad.TopH2OSoiLiq;
  auto imperv_TopIce = urban.imperviousRoad.TopH2OSoiIce;

  auto perv_TGrnd0 = urban.perviousRoad.TGrnd0;
  auto perv_Temp = urban.perviousRoad.EffectiveSurfTemp;
  auto perv_NetSw = urban.perviousRoad.NetShortRad;
  auto perv_NetLw = urban.perviousRoad.NetLongRad;
  auto perv_Emiss = urban.urbanParams.emissivity.PerviousRoad;
  auto perv_Cgrnds = urban.perviousRoad.Cgrnds;
  auto perv_Cgrndl = urban.perviousRoad.Cgrndl;
  auto perv_EflxShGrnd = urban.perviousRoad.EflxShGrnd;
  auto perv_QflxEvapSoil = urban.perviousRoad.QflxEvapSoil;
  auto perv_QflxTranEvap = urban.perviousRoad.QflxTranEvap;
  auto perv_EflxSoilGrnd = urban.perviousRoad.EflxSoilGrnd;
  auto perv_QflxEvapGrnd = urban.perviousRoad.QflxEvapGrnd;
  auto perv_QflxSubSnow = urban.perviousRoad.QflxSubSnow;
  auto perv_QflxDewSnow = urban.perviousRoad.QflxDewSnow;
  auto perv_QflxDewGrnd = urban.perviousRoad.QflxDewGrnd;
  // Pervious road uses layer-resolved H2OSoiLiq/Ice (index 0 = top layer)
  auto perv_H2OSoiLiq = urban.perviousRoad.H2OSoiLiq;
  auto perv_H2OSoiIce = urban.perviousRoad.H2OSoiIce;

  // Walls: only SurfaceDataBase views (no QflxEvapSoil/TranEvap/dew views)
  auto sunwall_TGrnd0 = urban.sunlitWall.TGrnd0;
  auto sunwall_Temp = urban.sunlitWall.EffectiveSurfTemp;
  auto sunwall_NetSw = urban.sunlitWall.NetShortRad;
  auto sunwall_NetLw = urban.sunlitWall.NetLongRad;
  auto sunwall_Cgrnds = urban.sunlitWall.Cgrnds;
  auto sunwall_Cgrndl = urban.sunlitWall.Cgrndl;
  auto sunwall_EflxShGrnd = urban.sunlitWall.EflxShGrnd;
  auto sunwall_EflxSoilGrnd = urban.sunlitWall.EflxSoilGrnd;

  auto shadewall_TGrnd0 = urban.shadedWall.TGrnd0;
  auto shadewall_Temp = urban.shadedWall.EffectiveSurfTemp;
  auto shadewall_NetSw = urban.shadedWall.NetShortRad;
  auto shadewall_NetLw = urban.shadedWall.NetLongRad;
  auto shadewall_Cgrnds = urban.shadedWall.Cgrnds;
  auto shadewall_Cgrndl = urban.shadedWall.Cgrndl;
  auto shadewall_EflxShGrnd = urban.shadedWall.EflxShGrnd;
  auto shadewall_EflxSoilGrnd = urban.shadedWall.EflxSoilGrnd;

  auto wall_Emiss = urban.urbanParams.emissivity.Wall;

  // AC heat from previous timestep — same value ELM's SoilFluxes uses for
  // eflx_heat_from_ac_patch(p) (set by UrbanFluxes before SoilTemperature,
  // saved by UrbanComputeHeatDiffusion into EFluxForAC_Prev).
  auto bldg_eflux_ac_prev = urban.building.EFluxForAC_Prev;

  Kokkos::parallel_for(
      "ComputeSoilFluxes", numLandunits, KOKKOS_LAMBDA(int l) {
        // ------------------------------------------------------------------
        // Combined lambda for SnowCoveredSurfaceData surfaces.
        // Order matches ELM SoilFluxesMod.F90 for urban columns:
        //   Step 1 — temperature correction to eflxShGrnd and qflxEvapSoil
        //   egirat  — clamp evap demand to topsoil availability; correct
        //             eflxShGrnd using scaled qflxEvapSoil
        //   Step 2 — eflxSoilGrnd (using egirat-corrected fluxes)
        //   Step 3 — four evap/dew fluxes (using egirat-scaled qflxEvapSoil)
        // ------------------------------------------------------------------
        auto computeSteps123_snow = [&](int l, Real tgnd0, Real &effectiveT,
                                        Real netSw, Real netLw, Real emiss,
                                        Real cgrnds, Real cgrndl, Real topLiq,
                                        Real topIce, Real &eflxShGrnd,
                                        Real &qflxEvapSoil, Real qflxTranEvap,
                                        Real &eflxSoilGrnd, Real &qflxEvapGrnd,
                                        Real &qflxSubSnow, Real &qflxDewSnow,
                                        Real &qflxDewGrnd) {
          // --- Step 1: temperature-correction to fluxes ---
          const Real tinc = effectiveT - tgnd0;
          eflxShGrnd += tinc * cgrnds;
          qflxEvapSoil += tinc * cgrndl;

          // htvp computed once; used consistently throughout
          // const Real htvp = (effectiveT > SHR_CONST_TKFRZ)
          //                      ? SHR_CONST_LATVAP
          //                      : (SHR_CONST_LATVAP + SHR_CONST_LATICE);
          const Real htvp = SHR_CONST_LATVAP;

          // --- egirat: clamp total evaporation to available topsoil water ---
          // (ELM SoilFluxesMod: computed and applied before eflx_soil_grnd)
          const Real dtime = 30.0 * 60.0; // seconds
          const Real total = topLiq + topIce;
          const Real egsmax = total / dtime;
          if (qflxEvapSoil > egsmax) {
            const Real egirat = egsmax / qflxEvapSoil;
            const Real save_qflxEvapSoil = qflxEvapSoil;
            qflxEvapSoil *= egirat;
            printf("(l = %d) eflxShGrnd = %18.16f; correction: (%18.16f - "
                   "%18.16f) * %18.16f = %18.16f\n",
                   l, eflxShGrnd, save_qflxEvapSoil, qflxEvapSoil, htvp,
                   (save_qflxEvapSoil - qflxEvapSoil) * htvp);
            eflxShGrnd += (save_qflxEvapSoil - qflxEvapSoil) * htvp;
          }

          // --- Step 2: EflxSoilGrnd (after egirat correction) ---
          // dlrad = 0 for all urban surfaces
          // wasteheat / ac / traffic = 0 for now
          const Real eflx_lwrad_del =
              4.0 * emiss * STEBOL * tgnd0 * tgnd0 * tgnd0 * tinc;

          eflxSoilGrnd = netSw - netLw - eflx_lwrad_del - eflxShGrnd -
                         qflxEvapSoil * htvp - qflxTranEvap * SHR_CONST_LATVAP;

          // --- Step 3: four evap/dew fluxes (after eflxSoilGrnd) ---
          // Uses egirat-scaled qflxEvapSoil; matches ELM order.
          qflxEvapGrnd = 0.0;
          qflxSubSnow = 0.0;
          qflxDewSnow = 0.0;
          qflxDewGrnd = 0.0;

          if (qflxEvapSoil >= 0.0) {
            if (total > 0.0) {
              qflxEvapGrnd = Kokkos::max(qflxEvapSoil * (topLiq / total), 0.0);
            }
            qflxSubSnow = qflxEvapSoil - qflxEvapGrnd;
          } else {
            if (effectiveT < SHR_CONST_TKFRZ) {
              qflxDewSnow = Kokkos::abs(qflxEvapSoil);
            } else {
              qflxDewGrnd = Kokkos::abs(qflxEvapSoil);
            }
          }
        };

        // ------------------------------------------------------------------
        // Helper lambda for WallDataType surfaces (Steps 1 and 2).
        // Walls have no QflxEvapSoil / QflxTranEvap views; both are zero.
        // Cgrndl is also 0 for walls so the qflxEvapSoil update is a no-op.
        // ------------------------------------------------------------------
        auto computeSteps12_wall = [&](Real tgnd0, Real &effectiveT, Real netSw,
                                       Real netLw, Real emiss, Real cgrnds,
                                       Real cgrndl, Real &eflxShGrnd,
                                       Real &eflxSoilGrnd) {
          // --- Step 1: temperature-correction (qflxEvapSoil = 0 for walls) ---
          const Real tinc = effectiveT - tgnd0;
          eflxShGrnd += tinc * cgrnds;
          // cgrndl = 0 for walls; qflxEvapSoil update is a no-op
          (void)cgrndl;

          // --- Step 2: EflxSoilGrnd (qflxEvapSoil = qflxTranEvap = 0) ---
          const Real eflx_lwrad_del =
              4.0 * emiss * STEBOL * tgnd0 * tgnd0 * tgnd0 * tinc;

          eflxSoilGrnd = netSw - netLw - eflx_lwrad_del - eflxShGrnd;
        };

        // ---- Roof ----
        computeSteps123_snow(
            l, roof_TGrnd0(l), roof_Temp(l), roof_NetSw(l), roof_NetLw(l),
            roof_Emiss(l), roof_Cgrnds(l), roof_Cgrndl(l), roof_TopLiq(l),
            roof_TopIce(l), roof_EflxShGrnd(l), roof_QflxEvapSoil(l),
            roof_QflxTranEvap(l), roof_EflxSoilGrnd(l), roof_QflxEvapGrnd(l),
            roof_QflxSubSnow(l), roof_QflxDewSnow(l), roof_QflxDewGrnd(l));

        // ---- Impervious road ----
        computeSteps123_snow(l, imperv_TGrnd0(l), imperv_Temp(l),
                             imperv_NetSw(l), imperv_NetLw(l), imperv_Emiss(l),
                             imperv_Cgrnds(l), imperv_Cgrndl(l),
                             imperv_TopLiq(l), imperv_TopIce(l),
                             imperv_EflxShGrnd(l), imperv_QflxEvapSoil(l),
                             imperv_QflxTranEvap(l), imperv_EflxSoilGrnd(l),
                             imperv_QflxEvapGrnd(l), imperv_QflxSubSnow(l),
                             imperv_QflxDewSnow(l), imperv_QflxDewGrnd(l));
        imperv_EflxSoilGrnd(l) += bldg_eflux_ac_prev(l);

        // ---- Pervious road (top layer = index 0) ----
        computeSteps123_snow(
            l, perv_TGrnd0(l), perv_Temp(l), perv_NetSw(l), perv_NetLw(l),
            perv_Emiss(l), perv_Cgrnds(l), perv_Cgrndl(l), perv_H2OSoiLiq(l, 0),
            perv_H2OSoiIce(l, 0), perv_EflxShGrnd(l), perv_QflxEvapSoil(l),
            perv_QflxTranEvap(l), perv_EflxSoilGrnd(l), perv_QflxEvapGrnd(l),
            perv_QflxSubSnow(l), perv_QflxDewSnow(l), perv_QflxDewGrnd(l));
        perv_EflxSoilGrnd(l) += bldg_eflux_ac_prev(l);

        // ---- Sunlit wall (Steps 1-2 only; walls have no evap/tran views) ----
        computeSteps12_wall(sunwall_TGrnd0(l), sunwall_Temp(l),
                            sunwall_NetSw(l), sunwall_NetLw(l), wall_Emiss(l),
                            sunwall_Cgrnds(l), sunwall_Cgrndl(l),
                            sunwall_EflxShGrnd(l), sunwall_EflxSoilGrnd(l));

        // ---- Shaded wall (Steps 1-2 only; walls have no evap/tran views) ----
        computeSteps12_wall(shadewall_TGrnd0(l), shadewall_Temp(l),
                            shadewall_NetSw(l), shadewall_NetLw(l),
                            wall_Emiss(l), shadewall_Cgrnds(l),
                            shadewall_Cgrndl(l), shadewall_EflxShGrnd(l),
                            shadewall_EflxSoilGrnd(l));
      });
  Kokkos::fence();
}

} // namespace URBANXX

// ============================================================================
// Public C API
// ============================================================================

extern "C" {

void UrbanComputeSoilFluxes(UrbanType urban, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    URBANXX::ComputeSoilFluxes(*urban);
    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
