#ifndef URBAN_GETTER_HELPERS_H
#define URBAN_GETTER_HELPERS_H

#include "Urban.h"
#include "private/DataTypesImpl.h"
#include <Kokkos_Core.hpp>

// Helper template functions for getting Kokkos views to C arrays
// These perform deep copies from device to host memory

// Base template function for getting 3D views
template <typename ViewType>
static void GetView3D(const ViewType &view, double *values, const int size[3],
                      UrbanErrorCode *status) {
  using namespace URBANXX;

  if (values == nullptr || size == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    // Check if each dimension matches the view extent
    if (size[0] != static_cast<int>(view.extent(0)) ||
        size[1] != static_cast<int>(view.extent(1)) ||
        size[2] != static_cast<int>(view.extent(2))) {
      *status = URBAN_ERR_SIZE_MISMATCH;
      return;
    }

    // Create an unmanaged host view wrapping the output array
    // Use the same layout as the source view for compatibility
    using target_layout = typename ViewType::array_layout;
    auto values_view =
        Kokkos::View<double ***, target_layout, Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>(values, size[0],
                                                              size[1], size[2]);

    // Deep copy from the device view to the host view
    Kokkos::deep_copy(values_view, view);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

// Base template function for getting 1D views
template <typename ViewType>
static void GetView1D(const ViewType &view, double *values, int length,
                      UrbanErrorCode *status) {
  using namespace URBANXX;

  if (values == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    // Check if the length matches the view extent
    if (length != static_cast<int>(view.extent(0))) {
      *status = URBAN_ERR_SIZE_MISMATCH;
      return;
    }

    // Create an unmanaged host view wrapping the output array
    // Use the same layout as the source view for compatibility
    using target_layout = typename ViewType::array_layout;
    auto values_view =
        Kokkos::View<double *, target_layout, Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>(values, length);

    // Deep copy from the device view to the host view
    Kokkos::deep_copy(values_view, view);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

#endif // URBAN_GETTER_HELPERS_H
