#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanParamsTypeImpl.h"
#include "private/UrbanTypeImpl.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

extern "C" {

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

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
