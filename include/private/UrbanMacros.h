#ifndef URBAN_MACROS_H
#define URBAN_MACROS_H

namespace URBANXX {

// Macro to declare dual arrays (host and device)
#define DECLARE_DUAL_ARRAY(suffix, name)                                       \
  HostArray##suffix name##H;                                                   \
  Array##suffix name

// Macro to allocate a single view
#define ALLOCATE_VIEW(viewname, type, ...)                                     \
  viewname = type(#viewname, __VA_ARGS__);

// Macro to allocate dual views (host and device)
#define ALLOCATE_DUAL_VIEWS(viewname, suffix, ...)                             \
  ALLOCATE_VIEW(viewname##H, HostArray##suffix, __VA_ARGS__);                  \
  ALLOCATE_VIEW(viewname, Array##suffix, __VA_ARGS__);

} // namespace URBANXX

#endif // URBAN_MACROS_H
