// Urban Heat Diffusion Module
// 1D heat diffusion dynamics for Urban Surfaces
// Based on ELM SoilTemperature.F90

#ifndef URBAN_HEAT_DIFFUSION_IMPL_H
#define URBAN_HEAT_DIFFUSION_IMPL_H

#include "private/UrbanTypeImpl.h"

namespace URBANXX {

// Compute 1D heat diffusion for all urban surfaces
void ComputeHeatDiffusion(URBANXX::_p_UrbanType &urban);

} // namespace URBANXX

#endif // URBAN_HEAT_DIFFUSION_IMPL_H
