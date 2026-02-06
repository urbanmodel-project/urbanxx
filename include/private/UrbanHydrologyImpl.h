#ifndef URBAN_HYDROLOGY_IMPL_H
#define URBAN_HYDROLOGY_IMPL_H

#include "private/DataTypesImpl.h"
#include "private/UrbanHydrologyFunctions.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// Forward declarations of main hydrology functions
void ComputeHydraulicProperties(UrbanType urban);
void SetupHydrologyTridiagonal(UrbanType urban, Real dtime);
void SolveHydrologyTridiagonal(UrbanType urban);
void UpdateSoilWater(UrbanType urban, Real dtime);

} // namespace URBANXX

#endif // URBAN_HYDROLOGY_IMPL_H
