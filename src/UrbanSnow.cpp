#include "private/UrbanSnowImpl.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>
#include <cmath>

namespace URBANXX {

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
  const double large_intsnow = 1.0e8; // ELM cap on int_snow

  auto forc_snow = urban.atmosphereData.ForcSnow; // kg/m²/s
  auto forc_temp = urban.atmosphereData.ForcTemp; // air temperature (K)

  auto roofH2OSno = urban.roof.H2OSno;
  auto roofIntSno = urban.roof.IntSno;
  auto roofFracSno = urban.roof.FracSno;
  auto roofSubSnow = urban.roof.QflxSubSnow;
  auto roofSnowMelt = urban.roof.QflxSnowMelt;
  auto roofSnowDepth = urban.roof.SnowDepth;
  auto roofNMelt = urban.roof.NMelt;

  auto impH2OSno = urban.imperviousRoad.H2OSno;
  auto impIntSno = urban.imperviousRoad.IntSno;
  auto impFracSno = urban.imperviousRoad.FracSno;
  auto impSubSnow = urban.imperviousRoad.QflxSubSnow;
  auto impSnowMelt = urban.imperviousRoad.QflxSnowMelt;
  auto impSnowDepth = urban.imperviousRoad.SnowDepth;
  auto impNMelt = urban.imperviousRoad.NMelt;

  auto perH2OSno = urban.perviousRoad.H2OSno;
  auto perIntSno = urban.perviousRoad.IntSno;
  auto perFracSno = urban.perviousRoad.FracSno;
  auto perSubSnow = urban.perviousRoad.QflxSubSnow;
  auto perSnowMelt = urban.perviousRoad.QflxSnowMelt;
  auto perSnowDepth = urban.perviousRoad.SnowDepth;
  auto perNMelt = urban.perviousRoad.NMelt;

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
                          Kokkos::View<double *> snow_melt_v,
                          Kokkos::View<double *> snow_depth_v,
                          Kokkos::View<double *> n_melt_v) {
          const double n_melt = n_melt_v(l);
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

            // ELM CanopyHydrology lines 557-563: frac_sno depletion driven by
            // thin-snow melt from the previous timestep.  Must run before the
            // accumulation branch (matches ELM ordering).
            const double snowmelt = snow_melt_v(l) * dtime;
            if (snowmelt > 0.0 && int_sno > 0.0) {
              double smr = Kokkos::min(1.0, h2osno / int_sno);
              frac_sno =
                  1.0 -
                  Kokkos::pow(acos(Kokkos::min(1.0, 2.0 * smr - 1.0)) / M_PI,
                              n_melt);
            }

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

        update(roofH2OSno, roofIntSno, roofFracSno, roofSubSnow, roofSnowMelt,
               roofSnowDepth, roofNMelt);
        update(impH2OSno, impIntSno, impFracSno, impSubSnow, impSnowMelt,
               impSnowDepth, impNMelt);
        update(perH2OSno, perIntSno, perFracSno, perSubSnow, perSnowMelt,
               perSnowDepth, perNMelt);
      });
  Kokkos::fence();
}

void ComputeUpdateSnowFraction(URBANXX::_p_UrbanType &urban) {
  const int numLandunits = urban.numLandunits;

  auto roofFracSno = urban.roof.FracSno;
  auto impFracSno = urban.imperviousRoad.FracSno;
  auto perFracSno = urban.perviousRoad.FracSno;

  auto roofSnowDepth = urban.roof.SnowDepth;
  auto impSnowDepth = urban.imperviousRoad.SnowDepth;
  auto perSnowDepth = urban.perviousRoad.SnowDepth;

  Kokkos::parallel_for(
      "ComputeUpdateSnowFraction", numLandunits, KOKKOS_LAMBDA(const int l) {
        // Bonan 1996 (LSM Technical Note) — same formula as ELM elm_driver.F90
        // step 7
        roofFracSno(l) = Kokkos::min(roofSnowDepth(l) / 0.05, 1.0);
        impFracSno(l) = Kokkos::min(impSnowDepth(l) / 0.05, 1.0);
        perFracSno(l) = Kokkos::min(perSnowDepth(l) / 0.05, 1.0);
      });
  Kokkos::fence();
}

} // namespace URBANXX
