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

# Test executable for radiation getters
add_executable(test_radiation_getters
    test_radiation_getters.cpp
)

target_link_libraries(test_radiation_getters
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME RadiationGetterTests COMMAND test_radiation_getters)

# Test executable for canyon air getters
add_executable(test_canyon_air_getters
    test_canyon_air_getters.cpp
)

target_link_libraries(test_canyon_air_getters
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME CanyonAirGetterTests COMMAND test_canyon_air_getters)

# Test executable for temperature setters
add_executable(test_temperature_setters
    test_temperature_setters.cpp
)

target_link_libraries(test_temperature_setters
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME TemperatureSetterTests COMMAND test_temperature_setters)

# Test executable for building temperature setter
add_executable(test_building_temperature_setter
    test_building_temperature_setter.cpp
)

target_link_libraries(test_building_temperature_setter
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME BuildingTemperatureSetterTests COMMAND test_building_temperature_setter)

# Test executable for net shortwave getters
add_executable(test_net_shortwave_getters
    test_net_shortwave_getters.cpp
)

target_link_libraries(test_net_shortwave_getters
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME NetShortwaveGetterTests COMMAND test_net_shortwave_getters)

# Test executable for surface flux getters
add_executable(test_surface_flux_getters
    test_surface_flux_getters.cpp
)

target_link_libraries(test_surface_flux_getters
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME SurfaceFluxGetterTests COMMAND test_surface_flux_getters)

# Test executable for layer temperature getters
add_executable(test_layer_temp_getters
    test_layer_temp_getters.cpp
)

target_link_libraries(test_layer_temp_getters
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME LayerTempGetterTests COMMAND test_layer_temp_getters)

# Test executable for hydrology getters
add_executable(test_hydrology_getters
    test_hydrology_getters.cpp
)

target_link_libraries(test_hydrology_getters
    urban
    kokkoscore
    gtest
)

# Register the test with CTest
add_test(NAME HydrologyGetterTests COMMAND test_hydrology_getters)
