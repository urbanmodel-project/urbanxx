# Urban library source files
set(URBAN_SOURCES
    Urban.cpp
    UrbanParamsType.cpp
    UrbanInitialize.cpp
    AtmosphereType.cpp
    UrbanLongwaveRad.cpp
    UrbanShortwaveRad.cpp
    KokkosWrappers.cpp
)

# Fortran module
set(URBAN_FORTRAN_SOURCES
    ${CMAKE_SOURCE_DIR}/include/urban_mod.F90
    f90-mod/Urban.F90
)

# Create static library
add_library(urban STATIC ${URBAN_SOURCES})

# Link against Kokkos
target_link_libraries(urban 
    PUBLIC 
    Kokkos::kokkos
)

# Public include directory (only Urban.h visible to users)
target_include_directories(urban 
    PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE
    ${CMAKE_SOURCE_DIR}/include/private
)

# Set C++ standard
target_compile_features(urban PUBLIC cxx_std_17)

# Create Fortran module library
add_library(urban_fortran STATIC ${URBAN_FORTRAN_SOURCES})
set_target_properties(urban_fortran PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
)

# Installation rules (optional)
install(TARGETS urban urban_fortran
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
install(FILES ${CMAKE_SOURCE_DIR}/include/Urban.h DESTINATION include)
install(DIRECTORY ${CMAKE_BINARY_DIR}/mod/ DESTINATION include/mod)
