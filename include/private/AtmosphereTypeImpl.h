#ifndef ATMOSPHERE_TYPE_IMPL_H
#define ATMOSPHERE_TYPE_IMPL_H

#include "private/DataTypesImpl.h"
#include "private/UrbanMacros.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

struct AtmosphereType {
  DECLARE_DUAL_VIEWS(1DR8, Coszen) // cosine solar zenith angle (-)
  DECLARE_DUAL_VIEWS(
      3DR8, ForcSRad) // solar radiation on horizontal surface (W/m**2)
                      // [landunit, band, type] where type is direct/diffuse
  DECLARE_DUAL_VIEWS(1DR8, FracSnow) // fraction of ground covered by snow (-)
  DECLARE_DUAL_VIEWS(1DR8, ForcLRad) // downwelling longwave radiation (W/m**2)
  DECLARE_DUAL_VIEWS(1DR8, ForcTemp) // air temperature (K)
  DECLARE_DUAL_VIEWS(1DR8, ForcPotTemp) // potential temperature (Pa)
  DECLARE_DUAL_VIEWS(1DR8, ForcRho)     // air density (kg/m**3)
  DECLARE_DUAL_VIEWS(1DR8, ForcSpcHumd) // specific humidity (kg/kg)
  DECLARE_DUAL_VIEWS(1DR8, ForcPress)   // atmospheric pressure (Pa)
  DECLARE_DUAL_VIEWS(1DR8, ForcWindU)   // wind speed in east direction (m/s)
  DECLARE_DUAL_VIEWS(1DR8, ForcWindV)   // wind speed in north direction (m/s)

  AtmosphereType(int numLandunits, int numRadBands, int numRadTypes) {
    ALLOCATE_DUAL_VIEWS(Coszen, 1DR8, numLandunits)
    ALLOCATE_DUAL_VIEWS(ForcSRad, 3DR8, numLandunits, numRadBands, numRadTypes)
    ALLOCATE_DUAL_VIEWS(FracSnow, 1DR8, numLandunits)
    ALLOCATE_DUAL_VIEWS(ForcLRad, 1DR8, numLandunits)
    ALLOCATE_DUAL_VIEWS(ForcTemp, 1DR8, numLandunits)
    ALLOCATE_DUAL_VIEWS(ForcPotTemp, 1DR8, numLandunits)
    ALLOCATE_DUAL_VIEWS(ForcRho, 1DR8, numLandunits)
    ALLOCATE_DUAL_VIEWS(ForcSpcHumd, 1DR8, numLandunits)
    ALLOCATE_DUAL_VIEWS(ForcPress, 1DR8, numLandunits)
    ALLOCATE_DUAL_VIEWS(ForcWindU, 1DR8, numLandunits)
    ALLOCATE_DUAL_VIEWS(ForcWindV, 1DR8, numLandunits)
  }
};

} // namespace URBANXX

#endif // ATMOSPHERE_TYPE_IMPL_H