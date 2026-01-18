# Test configuration for urbanxx
cmake_minimum_required(VERSION 3.16)

# Google Test is included with Kokkos
include_directories(${CMAKE_SOURCE_DIR}/externals/kokkos/tpls/gtest)

# Add unit tests subdirectory
add_subdirectory(unit)
