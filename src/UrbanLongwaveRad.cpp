#include "private/AtmosphereTypeImpl.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include "private/UrbanLongwaveRadImpl.h"
#include "private/UrbanParamsTypeImpl.h"
#include "private/UrbanSurfaceTypeImpl.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

void ComputeNetLongwave(URBANXX::_p_UrbanType &urban) {
  // Get number of landunits for parallel execution
  const int numLandunits = urban.numLandunits;

  // Access atmospheric forcing data
  auto& forcTemp = urban.atmosphereData.ForcTemp;
  auto& forcLRad = urban.atmosphereData.ForcLRad;

  // Access urban parameters - emissivities
  auto& emissRoof = urban.urbanParams.emissivity.Roof;
  auto& emissWall = urban.urbanParams.emissivity.Wall;
  auto& emissImpRoad = urban.urbanParams.emissivity.ImperviousRoad;
  auto& emissPerRoad = urban.urbanParams.emissivity.PerviousRoad;

  // Access surface temperatures
  auto& tempRoof = urban.roof.Temperature;
  auto& tempSunlitWall = urban.sunlitWall.Temperature;
  auto& tempShadedWall = urban.shadedWall.Temperature;
  auto& tempImpRoad = urban.imperviousRoad.Temperature;
  auto& tempPerRoad = urban.perviousRoad.Temperature;

  // Access net longwave radiation fields (to be updated)
  auto& netLwRoof = urban.roof.NetLongRad;
  auto& netLwSunlitWall = urban.sunlitWall.NetLongRad;
  auto& netLwShadedWall = urban.shadedWall.NetLongRad;
  auto& netLwImpRoad = urban.imperviousRoad.NetLongRad;
  auto& netLwPerRoad = urban.perviousRoad.NetLongRad;

  Kokkos::parallel_for(
      "ComputeNetLongwave", numLandunits, KOKKOS_LAMBDA(const int l) {
        // Net longwave calculation will go here
        // For now, set to zero as placeholder
        netLwRoof(l) = 0.0;
        netLwSunlitWall(l) = 0.0;
        netLwShadedWall(l) = 0.0;
        netLwImpRoad(l) = 0.0;
        netLwPerRoad(l) = 0.0;
      });

  Kokkos::fence();
}

} // namespace URBANXX
