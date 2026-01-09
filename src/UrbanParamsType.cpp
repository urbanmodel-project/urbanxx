#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanParamsTypeImpl.h"
#include "private/UrbanTypeImpl.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

extern "C" {

using namespace URBANXX;

void UrbanSetCanyonHwr(UrbanType urban, const double *values, int length,
                       UrbanErrorCode *status) {
  if (urban == nullptr || values == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    // Check if the length matches the view extent
    if (length != urban->urbanParams.CanyonHwr.extent(0)) {
      *status = URBAN_ERR_SIZE_MISMATCH;
      return;
    }

    // Create an unmanaged host view from the input array
    auto values_view =
        Kokkos::View<const double *, Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>(values, length);

    // Deep copy from the temporary host view to the device view
    Kokkos::deep_copy(urban->urbanParams.CanyonHwr, values_view);

    auto &CanyonHwr = urban->urbanParams.CanyonHwr;
    auto &sr = urban->urbanParams.viewFactor.SkyFrmRoad;
    auto &sw = urban->urbanParams.viewFactor.SkyFrmWall;
    auto &rw = urban->urbanParams.viewFactor.RoadFromWall;
    auto &wr = urban->urbanParams.viewFactor.WallFrmRoad;
    auto &ww = urban->urbanParams.viewFactor.OtherWallFrmWall;

    Kokkos::parallel_for(
        "ComputingViewFactor", urban->numLandunits, KOKKOS_LAMBDA(int l) {
          const Real hwr = CanyonHwr(l);
          const Real sqrt_term = sqrtf(hwr * hwr + 1.0);

          sr(l) = sqrt_term - hwr;                     // eqn 2.25
          wr(l) = 0.5 * (1.0 - sr(l));                 // eqn 2.27
          sw(l) = 0.5 * (hwr + 1.0 - sqrt_term) / hwr; // eqn 2.24
          rw(l) = sw(l);                               // eqn 2.27
          ww(l) = 1.0 - sw(l) - rw(l);                 // eqn 2.28
        });
    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
