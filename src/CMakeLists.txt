# Urban library source files
set(URBAN_SOURCES
    Urban.cpp
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

# Installation rules (optional)
install(TARGETS urban 
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
install(FILES ${CMAKE_SOURCE_DIR}/include/Urban.h DESTINATION include)
