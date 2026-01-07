# Driver executable
add_executable(urbanxx_driver driver.c)

# Link against Kokkos
target_link_libraries(urbanxx_driver 
    PRIVATE 
    Kokkos::kokkos
    MPI::MPI_C
    MPI::MPI_CXX
)

# Include directories
target_include_directories(urbanxx_driver 
    PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/private
)

# Set C++ compiler for the driver since it includes C++ headers
set_source_files_properties(driver.c PROPERTIES LANGUAGE CXX)

# Installation rules (optional)
install(TARGETS urbanxx_driver DESTINATION bin)
