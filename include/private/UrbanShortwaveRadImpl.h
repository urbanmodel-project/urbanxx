#ifndef URBAN_SHORTWAVE_RAD_IMPL_H
#define URBAN_SHORTWAVE_RAD_IMPL_H

#include "private/UrbanTypeImpl.h"

namespace URBANXX {

// Compute net shortwave radiation for all urban surfaces
void ComputeNetShortwave(URBANXX::_p_UrbanType &urban);

// Compute net shortwave radiation (scalar per landunit) for all five surfaces
// by multiplying AbsorbedShortRad with atmospheric forcing ForcSRad and
// summing over bands (VIS, NIR) and types (DIRECT, DIFFUSE).
void ComputeNetShortwaveRadiation(URBANXX::_p_UrbanType &urban);

} // namespace URBANXX

#endif // URBAN_SHORTWAVE_RAD_IMPL_H
