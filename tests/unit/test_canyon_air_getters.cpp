#include "Urban.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <cmath>

// Test fixture for canyon air property getter tests
class CanyonAirGetterTest : public ::testing::Test {
protected:
  UrbanType urban;
  const int numLandunits = 5;

  void SetUp() override {
    UrbanErrorCode ierr;
    UrbanCreate(numLandunits, &urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS) << "Failed to create Urban object in SetUp";
    
    // Set required parameters
    SetTestParameters();
    
    // Set atmospheric forcing
    SetTestAtmosphericForcing();
    
    // Run surface fluxes computation to generate canyon air data
    UrbanComputeSurfaceFluxes(urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS) << "Failed to compute surface fluxes";
  }

  void TearDown() override {
    UrbanErrorCode ierr;
    UrbanDestroy(&urban, &ierr);
  }

  void SetTestParameters() {
    UrbanErrorCode ierr;
    
    // Set canyon height-to-width ratio
    double canyon_hwr[5] = {1.0, 1.5, 2.0, 2.5, 3.0};
    UrbanSetCanyonHwr(urban, canyon_hwr, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    
    // Set fraction of pervious road
    double frac_perv[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
    UrbanSetFracPervRoadOfTotalRoad(urban, frac_perv, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    
    // Set roof weight
    double wt_roof[5] = {0.3, 0.35, 0.4, 0.45, 0.5};
    UrbanSetWtRoof(urban, wt_roof, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    
    // Set emissivities
    double emissivity[5] = {0.90, 0.91, 0.92, 0.93, 0.94};
    UrbanSetEmissivityRoof(urban, emissivity, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetEmissivityWall(urban, emissivity, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetEmissivityImperviousRoad(urban, emissivity, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetEmissivityPerviousRoad(urban, emissivity, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    
    // Set temperatures
    double temp[5] = {290.0, 291.0, 292.0, 293.0, 294.0};
    UrbanSetTemperatureRoof(urban, temp, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetTemperatureSunlitWall(urban, temp, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetTemperatureShadedWall(urban, temp, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetTemperatureImperviousRoad(urban, temp, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetTemperaturePerviousRoad(urban, temp, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    
    // Set height parameters
    double heights[5] = {10.0, 11.0, 12.0, 13.0, 14.0};
    UrbanSetForcHgtT(urban, heights, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetForcHgtU(urban, heights, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    
    double z0[5] = {0.1, 0.15, 0.2, 0.25, 0.3};
    UrbanSetZ0Town(urban, z0, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    
    double zd[5] = {1.0, 1.5, 2.0, 2.5, 3.0};
    UrbanSetZDTown(urban, zd, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    
    double ht_roof[5] = {8.0, 9.0, 10.0, 11.0, 12.0};
    UrbanSetHtRoof(urban, ht_roof, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    
    double wind_hgt[5] = {4.0, 4.5, 5.0, 5.5, 6.0};
    UrbanSetWindHgtCanyon(urban, wind_hgt, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
  }

  void SetTestAtmosphericForcing() {
    UrbanErrorCode ierr;
    
    // Set atmospheric forcing
    double atm_temp[5] = {300.0, 301.0, 302.0, 303.0, 304.0};
    double atm_pot_temp[5] = {301.0, 302.0, 303.0, 304.0, 305.0};
    double atm_press[5] = {101325.0, 101325.0, 101325.0, 101325.0, 101325.0};
    double atm_rho[5] = {1.2, 1.2, 1.2, 1.2, 1.2};
    double atm_spc_humd[5] = {0.01, 0.01, 0.01, 0.01, 0.01};
    double atm_wind_u[5] = {3.0, 3.5, 4.0, 4.5, 5.0};
    double atm_wind_v[5] = {1.0, 1.5, 2.0, 2.5, 3.0};
    
    UrbanSetAtmTemp(urban, atm_temp, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAtmPotTemp(urban, atm_pot_temp, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAtmPress(urban, atm_press, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAtmRho(urban, atm_rho, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAtmSpcHumd(urban, atm_spc_humd, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAtmWindU(urban, atm_wind_u, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAtmWindV(urban, atm_wind_v, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
  }
};

// ============================================================================
// Tests for Canyon Air Temperature Getter
// ============================================================================

TEST_F(CanyonAirGetterTest, GetCanyonAirTemperature_ValidData) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetCanyonAirTemperature(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i])) << "Value at index " << i << " is NaN";
    EXPECT_GT(values[i], 200.0) << "Temperature should be > 200K";
    EXPECT_LT(values[i], 400.0) << "Temperature should be < 400K";
  }
}

TEST_F(CanyonAirGetterTest, GetCanyonAirTemperature_LengthMismatch) {
  double values[10];
  UrbanErrorCode status;

  UrbanGetCanyonAirTemperature(urban, values, 10, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(CanyonAirGetterTest, GetCanyonAirTemperature_NullPointer) {
  UrbanErrorCode status;
  double values[5];

  // Null values pointer
  UrbanGetCanyonAirTemperature(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  // Null urban object
  UrbanGetCanyonAirTemperature(nullptr, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  // Null status pointer - should not crash
  UrbanGetCanyonAirTemperature(urban, values, numLandunits, nullptr);
}

// ============================================================================
// Tests for Canyon Air Humidity Getter
// ============================================================================

TEST_F(CanyonAirGetterTest, GetCanyonAirHumidity_ValidData) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetCanyonAirHumidity(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i])) << "Value at index " << i << " is NaN";
    EXPECT_GE(values[i], 0.0) << "Specific humidity should be non-negative";
    EXPECT_LE(values[i], 0.1) << "Specific humidity should be <= 0.1 kg/kg";
  }
}

TEST_F(CanyonAirGetterTest, GetCanyonAirHumidity_LengthMismatch) {
  double values[10];
  UrbanErrorCode status;

  UrbanGetCanyonAirHumidity(urban, values, 10, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(CanyonAirGetterTest, GetCanyonAirHumidity_NullPointer) {
  UrbanErrorCode status;
  double values[5];

  // Null values pointer
  UrbanGetCanyonAirHumidity(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  // Null urban object
  UrbanGetCanyonAirHumidity(nullptr, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  // Null status pointer - should not crash
  UrbanGetCanyonAirHumidity(urban, values, numLandunits, nullptr);
}

// ============================================================================
// Integration Tests
// ============================================================================

TEST_F(CanyonAirGetterTest, BothGetters_Sequential) {
  UrbanErrorCode status;
  
  double temperature[5], humidity[5];
  
  UrbanGetCanyonAirTemperature(urban, temperature, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetCanyonAirHumidity(urban, humidity, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  // Both should have valid data
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(temperature[i]));
    EXPECT_FALSE(std::isnan(humidity[i]));
  }
}

TEST_F(CanyonAirGetterTest, MultipleCallsSameGetter) {
  UrbanErrorCode status;
  
  double temp1[5], temp2[5];
  
  // Call the same getter twice
  UrbanGetCanyonAirTemperature(urban, temp1, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetCanyonAirTemperature(urban, temp2, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  // Values should be identical
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_DOUBLE_EQ(temp1[i], temp2[i]) 
        << "Multiple calls should return same values";
  }
}

TEST_F(CanyonAirGetterTest, TemperatureReasonableRange) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetCanyonAirTemperature(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  // Check values are in a physically reasonable range for urban canyon air
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_GT(values[i], 250.0) << "Canyon air temp should be > 250K";
    EXPECT_LT(values[i], 350.0) << "Canyon air temp should be < 350K";
  }
}

TEST_F(CanyonAirGetterTest, HumidityReasonableRange) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetCanyonAirHumidity(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  // Check values are in a physically reasonable range
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_GE(values[i], 0.0) << "Specific humidity should be >= 0";
    EXPECT_LE(values[i], 0.05) << "Specific humidity should be <= 0.05 for typical conditions";
  }
}

// Main function to initialize Kokkos
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  
  return result;
}
