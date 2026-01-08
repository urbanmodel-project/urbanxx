#ifndef ATMOSPHERE_TYPE_IMPL_H
#define ATMOSPHERE_TYPE_IMPL_H

#include "private/DataTypesImpl.h"
#include "private/UrbanMacros.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

struct AtmosphereType {
  // cosine solar zenith angle (-)
  DECLARE_DUAL_ARRAY(1DR8, Coszen);

  // direct beam solar radiation on horizontal surface (W/m**2)
  DECLARE_DUAL_ARRAY(2DR8, SdirHoriz);

  // diffuse solar radiation on horizontal surface (W/m**2)
  DECLARE_DUAL_ARRAY(2DR8, SdifHoriz);

  // fraction of ground covered by snow (-)
  DECLARE_DUAL_ARRAY(1DR8, FracSnow);

  // downwelling longwave radiation (W/m**2)
  DECLARE_DUAL_ARRAY(1DR8, DownwellingLongRad);

  // air temperature (K)
  DECLARE_DUAL_ARRAY(1DR8, ForcTemp);

  // potential temperature (Pa)
  DECLARE_DUAL_ARRAY(1DR8, ForcPotTemp);

  // air density (kg/m**3)
  DECLARE_DUAL_ARRAY(1DR8, ForcRho);

  // specific humidity (kg/kg)
  DECLARE_DUAL_ARRAY(1DR8, ForcSpcHumd);

  // atomspheric pressure (Pa)
  DECLARE_DUAL_ARRAY(1DR8, ForcPress);

  // wind speed in east direction (m/s)
  DECLARE_DUAL_ARRAY(1DR8, ForcWindU);

  // wind speed in north direction (m/s)
  DECLARE_DUAL_ARRAY(1DR8, ForcWindV);

  AtmosphereType(int numLandunits, int numRadBands) {
    ALLOCATE_DUAL_VIEWS(Coszen, 1DR8, numLandunits)
    ALLOCATE_DUAL_VIEWS(SdirHoriz, 2DR8, numLandunits, numRadBands)
    ALLOCATE_DUAL_VIEWS(SdifHoriz, 2DR8, numLandunits, numRadBands)
    ALLOCATE_DUAL_VIEWS(FracSnow, 1DR8, numLandunits)
    ALLOCATE_DUAL_VIEWS(DownwellingLongRad, 1DR8, numLandunits)
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