#ifndef URBAN_LONGWAVE_RAD_IMPL_H
#define URBAN_LONGWAVE_RAD_IMPL_H

#include "private/UrbanTypeImpl.h"

namespace URBANXX {

// Compute and update internal snow state (H2OSno, IntSno, FracSno) for all
// snow-covered horizontal urban surfaces (roof, impervious road, pervious
// road).  Must be called before ComputeNetLongwave each timestep.
void ComputeSnowCover(URBANXX::_p_UrbanType &urban);

// Compute net longwave radiation for all urban surfaces
void ComputeNetLongwave(URBANXX::_p_UrbanType &urban);

} // namespace URBANXX

#endif // URBAN_LONGWAVE_RAD_IMPL_H
