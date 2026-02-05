// Urban Thermal Property Functions
// Shared thermal conductivity calculations for urban surfaces
// Based on ELM SoilTemperatureMod.F90

#ifndef URBAN_THERMAL_FUNCTIONS_H
#define URBAN_THERMAL_FUNCTIONS_H

#include "private/UrbanConstants.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// Compute thermal conductivity at layer interfaces using harmonic averaging
// Based on ELM SoilTemperatureMod.F90
template <typename ViewType>
KOKKOS_INLINE_FUNCTION void ComputeInterfaceThermalConductivity(
    const int l, const int numLayers, const int numActiveLayers,
    const ViewType &tkLayer, const ViewType &tkInterface, const ViewType &zc,
    const ViewType &zi, const Real tkFillValue) {

  for (int k = 0; k < numLayers; ++k) {
    if (k < numActiveLayers) {
      // Harmonic average of layer thermal conductivities
      // tk(j) = thk(j)*thk(j+1)*(z(j+1)-z(j)) /
      //         (thk(j)*(z(j+1)-zi(j))+thk(j+1)*(zi(j)-z(j)))
      if (k < numLayers - 1) {
        tkInterface(l, k) = tkLayer(l, k) * tkLayer(l, k + 1) *
                            (zc(l, k + 1) - zc(l, k)) /
                            (tkLayer(l, k) * (zc(l, k + 1) - zi(l, k + 1)) +
                             tkLayer(l, k + 1) * (zi(l, k + 1) - zc(l, k)));
      } else {
        tkInterface(l, k) = tkLayer(l, k);
      }
    } else {
      // Inactive layers set to fill value
      tkInterface(l, k) = tkFillValue;
    }
  }
}

} // namespace URBANXX

#endif // URBAN_THERMAL_FUNCTIONS_H
