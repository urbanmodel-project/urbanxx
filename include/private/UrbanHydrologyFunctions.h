#ifndef URBAN_HYDROLOGY_FUNCTIONS_H
#define URBAN_HYDROLOGY_FUNCTIONS_H

#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// Ice impedance constant (Zhao et al. 1997)
constexpr Real E_ICE = 6.0;

// Minimum soil matric potential [mm]
constexpr Real SMP_MIN = -1.e8;

// ============================================================================
// Hydraulic Conductivity Functions
// ============================================================================

// Compute hydraulic conductivity using Clapp-Hornberger formulation
// hk = imped * hksat * s^(2b+3)
// where s = relative saturation (vol_liq / watsat)
KOKKOS_INLINE_FUNCTION Real ComputeHydraulicConductivity(
    const Real s,       // relative saturation [-]
    const Real hksat,   // saturated hydraulic conductivity [mm/s]
    const Real bsw,     // Clapp-Hornberger b parameter [-]
    const Real imped) { // ice impedance factor [-]

  const Real s_clamped = Kokkos::fmin(1.0, Kokkos::fmax(s, 0.01));
  const Real s2 = hksat * Kokkos::pow(s_clamped, 2.0 * bsw + 2.0);
  return imped * s_clamped * s2;
}

// Compute derivative of hydraulic conductivity with respect to volumetric water
// content dhk/dw = imped * (2b+3) * hksat * s^(2b+2) / watsat
KOKKOS_INLINE_FUNCTION Real ComputeHydraulicConductivityDerivative(
    const Real s,        // relative saturation [-]
    const Real hksat,    // saturated hydraulic conductivity [mm/s]
    const Real bsw,      // Clapp-Hornberger b parameter [-]
    const Real imped,    // ice impedance factor [-]
    const Real watsat) { // saturated water content [-]

  const Real s_clamped = Kokkos::fmin(1.0, Kokkos::fmax(s, 0.01));
  const Real s2 = hksat * Kokkos::pow(s_clamped, 2.0 * bsw + 2.0);
  return imped * (2.0 * bsw + 3.0) * s2 / watsat;
}

// Compute ice impedance factor
// imped = 10^(-e_ice * icefrac)
KOKKOS_INLINE_FUNCTION Real ComputeIceImpedance(const Real icefrac) {
  return Kokkos::pow(10.0, -E_ICE * icefrac);
}

// ============================================================================
// Matric Potential Functions
// ============================================================================

// Compute soil matric potential using Clapp-Hornberger formulation
// smp = -sucsat * s^(-b)
// where s = relative saturation (vol_liq / watsat)
KOKKOS_INLINE_FUNCTION Real
ComputeMatricPotential(const Real s,        // relative saturation [-]
                       const Real sucsat,   // saturated suction [mm]
                       const Real bsw,      // Clapp-Hornberger b parameter [-]
                       const Real smpmin) { // minimum matric potential [mm]

  const Real s_clamped = Kokkos::fmax(s, 0.01);
  const Real smp = -sucsat * Kokkos::pow(s_clamped, -bsw);
  return Kokkos::fmax(smpmin, smp);
}

// Compute derivative of matric potential with respect to volumetric water
// content dsmp/dw = -b * smp / (s * watsat)
KOKKOS_INLINE_FUNCTION Real ComputeMatricPotentialDerivative(
    const Real smp,     // soil matric potential [mm]
    const Real vol_liq, // volumetric liquid water content [-]
    const Real bsw) {   // Clapp-Hornberger b parameter [-]

  return -bsw * smp / vol_liq;
}

// ============================================================================
// Equilibrium Water Content Functions
// ============================================================================

// Compute equilibrium volumetric water content based on water table depth
// Three cases:
//   1. Layer fully saturated (depth > zwt): vol_eq = watsat
//   2. Layer partially saturated (zwt within layer): weighted average
//   3. Layer fully unsaturated (depth < zwt): equilibrium profile
KOKKOS_INLINE_FUNCTION Real ComputeEquilibriumWaterContent(
    const Real zwt,    // water table depth [mm]
    const Real z_top,  // top of layer depth [mm]
    const Real z_bot,  // bottom of layer depth [mm]
    const Real watsat, // saturated water content [-]
    const Real sucsat, // saturated suction [mm]
    const Real bsw) {  // Clapp-Hornberger b parameter [-]

  Real vol_eq;

  if (zwt <= z_top) {
    // Layer is fully saturated
    vol_eq = watsat;

  } else if (zwt > z_top && zwt < z_bot) {
    // Layer is partially saturated
    const Real temp_i = 1.0;
    const Real temp_0 =
        Kokkos::pow((sucsat + zwt - z_top) / sucsat, 1.0 - 1.0 / bsw);
    const Real vol_eq1 = -sucsat * watsat / (1.0 - 1.0 / bsw) / (zwt - z_top) *
                         (temp_i - temp_0);

    // Weighted average of saturated and unsaturated parts
    vol_eq =
        (vol_eq1 * (zwt - z_top) + watsat * (z_bot - zwt)) / (z_bot - z_top);
    vol_eq = Kokkos::fmin(watsat, Kokkos::fmax(vol_eq, 0.0));

  } else {
    // Layer is fully unsaturated
    const Real temp_i =
        Kokkos::pow((sucsat + zwt - z_bot) / sucsat, 1.0 - 1.0 / bsw);
    const Real temp_0 =
        Kokkos::pow((sucsat + zwt - z_top) / sucsat, 1.0 - 1.0 / bsw);
    vol_eq = -sucsat * watsat / (1.0 - 1.0 / bsw) / (z_bot - z_top) *
             (temp_i - temp_0);
    vol_eq = Kokkos::fmin(watsat, Kokkos::fmax(vol_eq, 0.0));
  }

  return vol_eq;
}

// Compute equilibrium matric potential for a given equilibrium water content
// zq = -sucsat * (vol_eq/watsat)^(-b)
KOKKOS_INLINE_FUNCTION Real ComputeEquilibriumMatricPotential(
    const Real vol_eq,   // equilibrium volumetric water content [-]
    const Real watsat,   // saturated water content [-]
    const Real sucsat,   // saturated suction [mm]
    const Real bsw,      // Clapp-Hornberger b parameter [-]
    const Real smpmin) { // minimum matric potential [mm]

  const Real s = Kokkos::fmax(vol_eq / watsat, 0.01);
  const Real zq = -sucsat * Kokkos::pow(s, -bsw);
  return Kokkos::fmax(smpmin, zq);
}

} // namespace URBANXX

#endif // URBAN_HYDROLOGY_FUNCTIONS_H
