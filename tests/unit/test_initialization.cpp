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

    const int numLevels = 15;
    const int totalSize = numLandunits * numLevels;
    double thermal_cond[75];  // 5 * 15
    double heat_cap[75];      // 5 * 15
    for (int i = 0; i < totalSize; ++i) {
      thermal_cond[i] = 1.0;
      heat_cap[i] = 1500.0;
    }
    const int size2D[2] = {numLandunits, numLevels};
    
    UrbanSetThermalConductivityRoad(urban, thermal_cond, size2D, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetThermalConductivityWall(urban, thermal_cond, size2D, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetThermalConductivityRoof(urban, thermal_cond, size2D, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);

    UrbanSetHeatCapacityRoad(urban, heat_cap, size2D, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetHeatCapacityWall(urban, heat_cap, size2D, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetHeatCapacityRoof(urban, heat_cap, size2D, &status);
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
    const int size3D[3] = {numLandunits, 2, 2}; // landunits x bands x rad_types
    double albedo[20];
    for (int i = 0; i < 20; ++i) albedo[i] = 0.20;
    
    UrbanSetAlbedoPerviousRoad(urban, albedo, size3D, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAlbedoImperviousRoad(urban, albedo, size3D, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAlbedoSunlitWall(urban, albedo, size3D, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAlbedoShadedWall(urban, albedo, size3D, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
    UrbanSetAlbedoRoof(urban, albedo, size3D, &status);
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
    UrbanSetAtmShortwaveDown(urban, shortwave, size3D, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);

    double frac_perv[5] = {0.2, 0.2, 0.2, 0.2, 0.2};
    UrbanSetFracPervRoadOfTotalRoad(urban, frac_perv, numLandunits, &status);
    ASSERT_EQ(status, URBAN_SUCCESS);
  }
};

// Test: UrbanSetup with valid parameters
TEST_F(InitializationTest, Setup_WithValidParameters) {
  UrbanErrorCode status;
  
  // Set all required parameters
  SetMinimalParameters();
  
  // Setup urban model
  UrbanSetup(urban, &status);
  
  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetup should succeed with valid parameters";
}

// Test: UrbanSetup with null urban object
TEST_F(InitializationTest, Setup_NullUrban) {
  UrbanErrorCode status;
  
  UrbanSetup(nullptr, &status);
  
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "UrbanSetup should fail with null urban object";
}

// Test: UrbanSetup with null status
TEST_F(InitializationTest, Setup_NullStatus) {
  // Set all required parameters
  SetMinimalParameters();
  
  // Should not crash with null status
  UrbanSetup(urban, nullptr);
}

// Test: Setup multiple times (should work)
TEST_F(InitializationTest, Setup_MultipleTimes) {
  UrbanErrorCode status;
  
  // Set all required parameters
  SetMinimalParameters();
  
  // First setup
  UrbanSetup(urban, &status);
  ASSERT_EQ(status, URBAN_SUCCESS);
  
  // Second setup (should also succeed)
  UrbanSetup(urban, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "Re-running setup should succeed";
}

// Test: Full workflow - Create, Set Parameters, Setup
TEST_F(InitializationTest, FullWorkflow_CreateSetSetup) {
  UrbanErrorCode status;
  
  // Set all parameters
  SetMinimalParameters();
  
  // Setup
  UrbanSetup(urban, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
}

// Test: Setup immediately after create (without setting parameters)
TEST_F(InitializationTest, Setup_WithoutParameters) {
  UrbanErrorCode status;
  
  // Try to setup without setting any parameters
  // This may succeed or fail depending on implementation
  // We're testing that it doesn't crash
  UrbanSetup(urban, &status);
  
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
