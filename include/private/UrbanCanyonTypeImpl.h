#ifndef URBAN_CANYON_TYPE_IMPL_H
#define URBAN_CANYON_TYPE_IMPL_H

#include "private/DataTypesImpl.h"
#include "private/KokkosViewMacros.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

struct UrbanCanyonType {
  DECLARE_DEVICE_VIEW(1DR8, Taf) // urban canopy air temperature (K)
  DECLARE_DEVICE_VIEW(1DR8, Qaf) // urban canopy air specific humidity (kg/kg)

  UrbanCanyonType(int numLandunits) {
    ALLOCATE_DEVICE_VIEW(Taf, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Qaf, Array1DR8, numLandunits)
  }
};

} // namespace URBANXX

#endif // URBAN_CANYON_TYPE_IMPL_H
