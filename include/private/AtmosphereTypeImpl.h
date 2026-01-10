#ifndef ATMOSPHERE_TYPE_IMPL_H
#define ATMOSPHERE_TYPE_IMPL_H

#include "private/DataTypesImpl.h"
#include "private/KokkosViewMacros.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

struct AtmosphereType {
  DECLARE_DEVICE_VIEW(1DR8, Coszen) // cosine solar zenith angle (-)
  DECLARE_DEVICE_VIEW(
      3DR8, ForcSRad) // solar radiation on horizontal surface (W/m**2)
                      // [landunit, band, type] where type is direct/diffuse
  DECLARE_DEVICE_VIEW(1DR8, FracSnow) // fraction of ground covered by snow (-)
  DECLARE_DEVICE_VIEW(1DR8, ForcLRad) // downwelling longwave radiation (W/m**2)
  DECLARE_DEVICE_VIEW(1DR8, ForcTemp) // air temperature (K)
  DECLARE_DEVICE_VIEW(1DR8, ForcPotTemp) // potential temperature (Pa)
  DECLARE_DEVICE_VIEW(1DR8, ForcRho)     // air density (kg/m**3)
  DECLARE_DEVICE_VIEW(1DR8, ForcSpcHumd) // specific humidity (kg/kg)
  DECLARE_DEVICE_VIEW(1DR8, ForcPress)   // atmospheric pressure (Pa)
  DECLARE_DEVICE_VIEW(1DR8, ForcWindU)   // wind speed in east direction (m/s)
  DECLARE_DEVICE_VIEW(1DR8, ForcWindV)   // wind speed in north direction (m/s)

  AtmosphereType(int numLandunits, int numRadBands, int numRadTypes) {
    ALLOCATE_DEVICE_VIEW(Coszen, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ForcSRad, Array3DR8, numLandunits, numRadBands, numRadTypes)
    ALLOCATE_DEVICE_VIEW(FracSnow, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ForcLRad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ForcTemp, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ForcPotTemp, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ForcRho, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ForcSpcHumd, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ForcPress, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ForcWindU, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(ForcWindV, Array1DR8, numLandunits)
  }
};

} // namespace URBANXX

#endif // ATMOSPHERE_TYPE_IMPL_H