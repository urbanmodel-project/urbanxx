# Driver executable
add_executable(urbanxx_driver driver.c)

# Link against liburban.a, Kokkos, and MPI
target_link_libraries(urbanxx_driver 
    PRIVATE 
    urban
    Kokkos::kokkos
    MPI::MPI_C
    MPI::MPI_CXX
)

# Include directories (only public Urban.h visible)
target_include_directories(urbanxx_driver 
    PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
)

# Set C++ compiler for the driver since it includes C++ headers
set_source_files_properties(driver.c PROPERTIES LANGUAGE CXX)

# Installation rules (optional)
install(TARGETS urbanxx_driver DESTINATION bin)
