#include "private/AtmosphereTypeImpl.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include "private/UrbanDebugUtils.h"
#include "private/UrbanParamsTypeImpl.h"
#include "private/UrbanShortwaveRadImpl.h"
#include "private/UrbanSurfaceTypeImpl.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>
#include <iostream>

namespace URBANXX {

// Compute net shortwave radiation for all urban surfaces
void ComputeNetShortwave(URBANXX::_p_UrbanType &urban) {
  // Placeholder implementation
  // TODO: Add shortwave radiation physics

  const int numLandunits = urban.numLandunits;

  std::cout << "ComputeNetShortwave: Processing " << numLandunits
            << " landunits (placeholder)" << std::endl;
}

} // namespace URBANXX
