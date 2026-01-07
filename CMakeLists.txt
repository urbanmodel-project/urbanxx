cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(UrbanXX LANGUAGES CXX C)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options for Kokkos backends
option(ENABLE_CUDA "Enable CUDA backend" OFF)
option(ENABLE_HIP "Enable HIP backend" OFF)
option(ENABLE_SYCL "Enable SYCL backend" OFF)
option(ENABLE_OPENMP "Enable OpenMP backend" ON)
option(ENABLE_SERIAL "Enable Serial backend" ON)
option(SINGLE_PRECISION "Use single precision" OFF)
option(LAYOUT_RIGHT "Use LayoutRight" OFF)
option(LAYOUT_LEFT "Use LayoutLeft" ON)

# Find MPI
find_package(MPI REQUIRED)

# Add Kokkos as a subdirectory
set(Kokkos_ENABLE_SERIAL ${ENABLE_SERIAL} CACHE BOOL "Enable Serial backend")
set(Kokkos_ENABLE_OPENMP ${ENABLE_OPENMP} CACHE BOOL "Enable OpenMP backend")
set(Kokkos_ENABLE_CUDA ${ENABLE_CUDA} CACHE BOOL "Enable CUDA backend")
set(Kokkos_ENABLE_HIP ${ENABLE_HIP} CACHE BOOL "Enable HIP backend")
set(Kokkos_ENABLE_SYCL ${ENABLE_SYCL} CACHE BOOL "Enable SYCL backend")

add_subdirectory(externals/kokkos)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/externals/kokkos/core/src)

# Add compile definitions based on options
if(ENABLE_CUDA)
    add_compile_definitions(ENABLE_CUDA)
endif()

if(ENABLE_HIP)
    add_compile_definitions(ENABLE_HIP)
endif()

if(ENABLE_SYCL)
    add_compile_definitions(ENABLE_SYCL)
endif()

if(ENABLE_OPENMP)
    add_compile_definitions(ENABLE_OPENMP)
endif()

if(ENABLE_SERIAL)
    add_compile_definitions(ENABLE_SERIAL)
endif()

if(SINGLE_PRECISION)
    add_compile_definitions(SINGLE_PRECISION)
endif()

if(LAYOUT_RIGHT)
    add_compile_definitions(LAYOUT_RIGHT)
elseif(LAYOUT_LEFT)
    add_compile_definitions(LAYOUT_LEFT)
endif()

# Add driver subdirectory
add_subdirectory(driver)
