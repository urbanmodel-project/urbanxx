#ifndef URBAN_MACROS_H
#define URBAN_MACROS_H

namespace URBANXX {

#define DECLARE_HOST_VIEW(suffix, name) HostArray##suffix name##H;
#define DECLARE_DEVICE_VIEW(suffix, name) Array##suffix name;

// Macro to declare dual arrays (host and device)
#define DECLARE_DUAL_VIEWS(suffix, name)                                       \
  HostArray##suffix name##H;                                                   \
  Array##suffix name;

// Macro to allocate a single view
#define ALLOCATE_VIEW(viewname, type, ...)                                     \
  viewname = type(#viewname, __VA_ARGS__);

#define ALLOCATE_DEVICE_VIEW(viewname, type, ...)                              \
  ALLOCATE_VIEW(viewname, type, __VA_ARGS__);
#define ALLOCATE_HOST_VIEW(viewname, type, ...)                                \
  ALLOCATE_VIEW(viewname##H, type, __VA_ARGS__);

// Macro to allocate dual views (host and device)
#define ALLOCATE_DUAL_VIEWS(viewname, suffix, ...)                             \
  ALLOCATE_VIEW(viewname##H, HostArray##suffix, __VA_ARGS__);                  \
  ALLOCATE_VIEW(viewname, Array##suffix, __VA_ARGS__);

} // namespace URBANXX

#endif // URBAN_MACROS_H
