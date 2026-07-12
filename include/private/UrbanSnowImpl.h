#ifndef URBAN_SNOW_IMPL_H
#define URBAN_SNOW_IMPL_H

#include "private/UrbanTypeImpl.h"

namespace URBANXX {

// Update internal snow state (H2OSno, IntSno, FracSno, SnowDepth) for all
// snow-covered horizontal urban surfaces (roof, impervious road, pervious
// road).  Walls are excluded — vertical surfaces do not accumulate snow.
// Must be called before ComputeNetLongwave each timestep.
void ComputeSnowCover(URBANXX::_p_UrbanType &urban);

// Overwrite FracSno with the Bonan 1996 formula (snow_depth/0.05, capped at 1)
// for roof, imperviousRoad, and perviousRoad.  Must be called after hydrology
// updates SnowDepth and before urbanxx_netShortwave (ELM driver step 8).
void ComputeUpdateSnowFraction(URBANXX::_p_UrbanType &urban);

} // namespace URBANXX

#endif // URBAN_SNOW_IMPL_H
