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

extern "C" {

void UrbanComputeSurfaceRunoff(UrbanType urban, Real dtime,
                               UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    const int nlandunits = urban->numLandunits;

    auto atm_rain = urban->atmosphereData.ForcRain;

    // Roof views
    auto roof_top_liq = urban->roof.TopH2OSoiLiq;
    auto roof_evap = urban->roof.QflxEvapGrnd;
    auto roof_surf = urban->roof.QflxSurf;

    // Impervious road views
    auto imperv_top_liq = urban->imperviousRoad.TopH2OSoiLiq;
    auto imperv_evap = urban->imperviousRoad.QflxEvapGrnd;
    auto imperv_surf = urban->imperviousRoad.QflxSurf;

    // Pervious road views
    auto perv_wtfact = urban->perviousRoad.Wtfact;
    auto perv_fover = urban->perviousRoad.Fover;
    auto perv_zwt = urban->perviousRoad.Zwt;
    auto perv_frost = urban->perviousRoad.FrostTable;
    auto perv_zwtperch = urban->perviousRoad.ZwtPerched;
    auto perv_surf = urban->perviousRoad.QflxSurf;

    // Wall views
    auto sunlit_surf = urban->sunlitWall.QflxSurf;
    auto shaded_surf = urban->shadedWall.QflxSurf;

    Kokkos::parallel_for(
        "UrbanComputeSurfaceRunoff", nlandunits, KOKKOS_LAMBDA(const int l) {
          const Real qflx_top = atm_rain(l); // qflx_top_soil = ForcRain

          // ----------------------------------------------------------------
          // Roof  (snl == 0 assumed — no snow layers)
          // ----------------------------------------------------------------
          {
            const Real top_liq = roof_top_liq(l);
            const Real evap = roof_evap(l);
            const Real xs = Kokkos::fmax(0.0, top_liq / dtime + qflx_top -
                                                  evap - PONDMX_URBAN / dtime);
            if (xs > 0.0) {
              roof_top_liq(l) = PONDMX_URBAN;
            } else {
              roof_top_liq(l) =
                  Kokkos::fmax(0.0, top_liq + (qflx_top - evap) * dtime);
            }
            roof_surf(l) = xs;
          }

          // ----------------------------------------------------------------
          // Impervious road  (snl == 0 assumed — no snow layers)
          // ----------------------------------------------------------------
          {
            const Real top_liq = imperv_top_liq(l);
            const Real evap = imperv_evap(l);
            const Real xs = Kokkos::fmax(0.0, top_liq / dtime + qflx_top -
                                                  evap - PONDMX_URBAN / dtime);
            if (xs > 0.0) {
              imperv_top_liq(l) = PONDMX_URBAN;
            } else {
              imperv_top_liq(l) =
                  Kokkos::fmax(0.0, top_liq + (qflx_top - evap) * dtime);
            }
            imperv_surf(l) = xs;
          }

          // ----------------------------------------------------------------
          // Pervious road
          // ----------------------------------------------------------------
          {
            const Real fff = perv_fover(l);
            Real fsat = perv_wtfact(l) * Kokkos::exp(-0.5 * fff * perv_zwt(l));
            if (perv_frost(l) > perv_zwtperch(l)) {
              fsat =
                  perv_wtfact(l) * Kokkos::exp(-0.5 * fff * perv_zwtperch(l));
            }
            perv_surf(l) = fsat * qflx_top;
          }

          // ----------------------------------------------------------------
          // Walls — always zero
          // ----------------------------------------------------------------
          sunlit_surf(l) = 0.0;
          shaded_surf(l) = 0.0;
        });

    Kokkos::fence();
    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
