# Unit tests
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/externals/kokkos/tpls/gtest/gtest)

# Find threading library for gtest
find_package(Threads REQUIRED)

# Build gtest from source (included with Kokkos)
add_library(gtest STATIC
    ${CMAKE_SOURCE_DIR}/externals/kokkos/tpls/gtest/gtest/gtest-all.cc
)
target_include_directories(gtest PUBLIC
    ${CMAKE_SOURCE_DIR}/externals/kokkos/tpls/gtest
)
target_link_libraries(gtest PUBLIC
    Threads::Threads
)

# Test executable for basic functionality
add_executable(test_urban_basic
    test_urban_create_destroy.cpp
)

target_link_libraries(test_urban_basic
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME UrbanBasicTests COMMAND test_urban_basic)

# Test executable for SetView helper functions
add_executable(test_setview_helpers
    test_setview_helpers.cpp
)

target_link_libraries(test_setview_helpers
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME SetViewHelperTests COMMAND test_setview_helpers)

# Test executable for parameter setters
add_executable(test_parameter_setters
    test_parameter_setters.cpp
)

target_link_libraries(test_parameter_setters
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME ParameterSetterTests COMMAND test_parameter_setters)

# Test executable for albedo setters
add_executable(test_albedo_setters
    test_albedo_setters.cpp
)

target_link_libraries(test_albedo_setters
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME AlbedoSetterTests COMMAND test_albedo_setters)

# Test executable for initialization functions
add_executable(test_initialization
    test_initialization.cpp
)

target_link_libraries(test_initialization
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME InitializationTests COMMAND test_initialization)

# Test executable for atmospheric shortwave setter
add_executable(test_atm_shortwave
    test_atm_shortwave.cpp
)

target_link_libraries(test_atm_shortwave
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME AtmShortwaveTests COMMAND test_atm_shortwave)
