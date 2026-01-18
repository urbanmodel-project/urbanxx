#include "Urban.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

// Test fixture for initialization tests
class InitializationTest : public ::testing::Test {
protected:
  UrbanType urban;
  const int numLandunits = 5;

  void SetUp() override {
    UrbanErrorCode ierr;
    UrbanCreate(numLandunits, &urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS) << "Failed to create Urban object in SetUp";
  }

  void TearDown() override {
    UrbanErrorCode ierr;
    UrbanDestroy(&urban, &ierr);
  }

  // Helper to set minimum required parameters for initialization
  void SetMinimalParameters() {
    UrbanErrorCode status;
    
    // Set basic geometric parameters
    double canyon_hwr[5] = {2.0, 2.0, 2.0, 2.0, 2.0};
    UrbanSetCanyonHwr(urban, canyon_hwr, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);

    double wt_roof[5] = {0.5, 0.5, 0.5, 0.5, 0.5};
    UrbanSetWtRoof(urban, wt_roof, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);

    // Set thermal properties
    double emissivity[5] = {0.90, 0.90, 0.90, 0.90, 0.90};
    UrbanSetEmissivityPerviousRoad(urban, emissivity, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetEmissivityImperviousRoad(urban, emissivity, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetEmissivityWall(urban, emissivity, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetEmissivityRoof(urban, emissivity, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);

    double thermal_cond[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
    UrbanSetThermalConductivityRoad(urban, thermal_cond, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetThermalConductivityWall(urban, thermal_cond, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetThermalConductivityRoof(urban, thermal_cond, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);

    double heat_cap[5] = {1500.0, 1500.0, 1500.0, 1500.0, 1500.0};
    UrbanSetHeatCapacityRoad(urban, heat_cap, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetHeatCapacityWall(urban, heat_cap, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetHeatCapacityRoof(urban, heat_cap, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);

    // Set height parameters
    double heights[5] = {10.0, 10.0, 10.0, 10.0, 10.0};
    UrbanSetForcHgtT(urban, heights, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetForcHgtU(urban, heights, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetZDTown(urban, heights, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetZ0Town(urban, heights, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetHtRoof(urban, heights, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetWindHgtCanyon(urban, heights, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);

    // Set albedos
    const int size[3] = {numLandunits, 2, 2}; // landunits x bands x rad_types
    double albedo[20];
    for (int i = 0; i < 20; ++i) albedo[i] = 0.20;
    
    UrbanSetAlbedoPerviousRoad(urban, albedo, size, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAlbedoImperviousRoad(urban, albedo, size, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAlbedoSunlitWall(urban, albedo, size, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAlbedoShadedWall(urban, albedo, size, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAlbedoRoof(urban, albedo, size, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);

    // Set atmospheric forcing
    double atm_values[5] = {300.0, 300.0, 300.0, 300.0, 300.0};
    UrbanSetAtmTemp(urban, atm_values, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAtmPotTemp(urban, atm_values, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAtmRho(urban, atm_values, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAtmSpcHumd(urban, atm_values, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAtmPress(urban, atm_values, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAtmWindU(urban, atm_values, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAtmWindV(urban, atm_values, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAtmCoszen(urban, atm_values, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAtmFracSnow(urban, atm_values, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAtmLongwaveDown(urban, atm_values, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);

    double shortwave[20];
    for (int i = 0; i < 20; ++i) shortwave[i] = 500.0;
    UrbanSetAtmShortwaveDown(urban, shortwave, size, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);

    double frac_perv[5] = {0.2, 0.2, 0.2, 0.2, 0.2};
    UrbanSetFracPervRoadOfTotalRoad(urban, frac_perv, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
  }
};

// Test: UrbanInitializeTemperature with valid parameters
TEST_F(InitializationTest, InitializeTemperature_WithValidParameters) {
  UrbanErrorCode status;
  
  // Set all required parameters
  SetMinimalParameters();
  
  // Initialize temperature
  UrbanInitializeTemperature(urban, &status);
  
  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanInitializeTemperature should succeed with valid parameters";
}

// Test: UrbanInitializeTemperature with null urban object
TEST_F(InitializationTest, InitializeTemperature_NullUrban) {
  UrbanErrorCode status;
  
  UrbanInitializeTemperature(nullptr, &status);
  
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "UrbanInitializeTemperature should fail with null urban object";
}

// Test: UrbanInitializeTemperature with null status
TEST_F(InitializationTest, InitializeTemperature_NullStatus) {
  // Set all required parameters
  SetMinimalParameters();
  
  // Should not crash with null status
  UrbanInitializeTemperature(urban, nullptr);
}

// Test: Initialize temperature multiple times (should work)
TEST_F(InitializationTest, InitializeTemperature_MultipleTimes) {
  UrbanErrorCode status;
  
  // Set all required parameters
  SetMinimalParameters();
  
  // First initialization
  UrbanInitializeTemperature(urban, &status);
  ASSERT_EQ(status, URBAN_SUCCESS);
  
  // Second initialization (should also succeed)
  UrbanInitializeTemperature(urban, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "Re-initializing temperature should succeed";
}

// Test: Full workflow - Create, Set Parameters, Initialize
TEST_F(InitializationTest, FullWorkflow_CreateSetInitialize) {
  UrbanErrorCode status;
  
  // Set all parameters
  SetMinimalParameters();
  
  // Initialize
  UrbanInitializeTemperature(urban, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
}

// Test: Initialize immediately after create (without setting parameters)
TEST_F(InitializationTest, InitializeTemperature_WithoutParameters) {
  UrbanErrorCode status;
  
  // Try to initialize without setting any parameters
  // This may succeed or fail depending on implementation
  // We're testing that it doesn't crash
  UrbanInitializeTemperature(urban, &status);
  
  // Just check that we get a valid error code (not necessarily success)
  EXPECT_TRUE(status == URBAN_SUCCESS || 
              status == URBAN_ERR_NOT_INITIALIZED || 
              status == URBAN_ERR_INVALID_ARGUMENT)
      << "Should return a valid error code, got: " << status;
}

// Main function to initialize Kokkos
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  
  return result;
}
