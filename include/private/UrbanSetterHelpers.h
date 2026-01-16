#ifndef URBAN_SETTER_HELPERS_H
#define URBAN_SETTER_HELPERS_H

#include "Urban.h"
#include "private/DataTypesImpl.h"
#include <Kokkos_Core.hpp>
#include <iostream>

// Helper template functions for setting Kokkos views from C arrays
// These are used by both UrbanParamsType and AtmosphereType setters

// Base template function for setting 3D views
template <typename ViewType>
static void SetView3D(ViewType &view, const double *values, const int size[3],
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

    // Create an unmanaged host view from the input array
    // Use the same layout as the target view for compatibility
    using target_layout = typename ViewType::array_layout;
    auto values_view =
        Kokkos::View<const double ***, target_layout, Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>(values, size[0],
                                                              size[1], size[2]);

    // Deep copy from the temporary host view to the target view
    Kokkos::deep_copy(view, values_view);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

// Base template function for setting 1D views
template <typename ViewType>
static void SetView1D(ViewType &view, const double *values, int length,
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

    // Create an unmanaged host view from the input array
    // Use the same layout as the target view for compatibility
    using target_layout = typename ViewType::array_layout;
    auto values_view =
        Kokkos::View<const double *, target_layout, Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>(values, length);

    // Deep copy from the temporary host view to the target view
    Kokkos::deep_copy(view, values_view);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

#endif // URBAN_SETTER_HELPERS_H
