#ifndef KOKKOS_VIEW_MACROS_H
#define KOKKOS_VIEW_MACROS_H

namespace URBANXX {

#define DECLARE_HOST_VIEW(suffix, name) HostArray##suffix name##H;
#define DECLARE_DEVICE_VIEW(suffix, name) Array##suffix name;

// Macro to declare dual arrays (host and device)
#define DECLARE_DUAL_VIEWS(suffix, name)                                       \
  HostArray##suffix name##H;                                                   \
  Array##suffix name;

// Macro to allocate a single view (no initialization)
#define ALLOCATE_VIEW_NO_INIT(viewname, type, ...)                             \
  viewname = type(#viewname, __VA_ARGS__);

// Macro to allocate a single floating-point view with zero initialization
#define ALLOCATE_VIEW(viewname, type, ...)                                     \
  ALLOCATE_VIEW_NO_INIT(viewname, type, __VA_ARGS__)                           \
  Kokkos::deep_copy(viewname, 0.0);

// Macro to allocate a single integer view with zero initialization
#define ALLOCATE_VIEW_INT(viewname, type, ...)                                 \
  ALLOCATE_VIEW_NO_INIT(viewname, type, __VA_ARGS__)                           \
  Kokkos::deep_copy(viewname, 0);

#define ALLOCATE_DEVICE_VIEW(viewname, type, ...)                              \
  ALLOCATE_VIEW(viewname, type, __VA_ARGS__)
#define ALLOCATE_HOST_VIEW(viewname, type, ...)                                \
  ALLOCATE_VIEW(viewname##H, type, __VA_ARGS__)

#define ALLOCATE_DEVICE_VIEW_INT(viewname, type, ...)                          \
  ALLOCATE_VIEW_INT(viewname, type, __VA_ARGS__)
#define ALLOCATE_HOST_VIEW_INT(viewname, type, ...)                            \
  ALLOCATE_VIEW_INT(viewname##H, type, __VA_ARGS__)

// Macro to allocate dual views (host and device)
#define ALLOCATE_DUAL_VIEWS(viewname, suffix, ...)                             \
  ALLOCATE_VIEW(viewname##H, HostArray##suffix, __VA_ARGS__)                   \
  ALLOCATE_VIEW(viewname, Array##suffix, __VA_ARGS__)

} // namespace URBANXX

#endif // KOKKOS_VIEW_MACROS_H
