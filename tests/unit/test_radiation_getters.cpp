#include "Urban.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <cmath>

// Test fixture for radiation getter tests
class RadiationGetterTest : public ::testing::Test {
protected:
  UrbanType urban;
  const int numLandunits = 5;
  const int numBands = 2;
  const int numRadTypes = 2; // Direct and diffuse
  int size[3];

  void SetUp() override {
    UrbanErrorCode ierr;
    UrbanCreate(numLandunits, &urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS) << "Failed to create Urban object in SetUp";
    
    size[0] = numLandunits;
    size[1] = numBands;
    size[2] = numRadTypes;
    
    // Set required parameters
    SetTestParameters();
    
    // Set atmospheric forcing
    SetTestAtmosphericForcing();
    
    // Run physics computations to generate data
    UrbanComputeNetShortwave(urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS) << "Failed to compute net shortwave";
    
    UrbanComputeNetLongwave(urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS) << "Failed to compute net longwave";
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
    
    // Set albedos (3D arrays)
    const int total_size = numLandunits * numBands * numRadTypes;
    double albedo_roof[20], albedo_wall[20], albedo_road[20];
    for (int i = 0; i < total_size; ++i) {
      albedo_roof[i] = 0.15 + (i % 5) * 0.01;
      albedo_wall[i] = 0.25 + (i % 5) * 0.01;
      albedo_road[i] = 0.10 + (i % 5) * 0.01;
    }
    UrbanSetAlbedoRoof(urban, albedo_roof, size, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAlbedoSunlitWall(urban, albedo_wall, size, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAlbedoShadedWall(urban, albedo_wall, size, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAlbedoImperviousRoad(urban, albedo_road, size, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAlbedoPerviousRoad(urban, albedo_road, size, &ierr);
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
  }

  void SetTestAtmosphericForcing() {
    UrbanErrorCode ierr;
    
    // Set atmospheric forcing (1D)
    double atm_temp[5] = {300.0, 301.0, 302.0, 303.0, 304.0};
    double atm_press[5] = {101325.0, 101325.0, 101325.0, 101325.0, 101325.0};
    double atm_longwave[5] = {400.0, 410.0, 420.0, 430.0, 440.0};
    double coszen[5] = {0.5, 0.6, 0.7, 0.8, 0.9};
    
    UrbanSetAtmTemp(urban, atm_temp, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAtmPotTemp(urban, atm_temp, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAtmPress(urban, atm_press, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAtmLongwaveDown(urban, atm_longwave, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAtmCoszen(urban, coszen, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    
    // Set atmospheric shortwave forcing (3D)
    const int total_size = numLandunits * numBands * numRadTypes;
    double shortwave[20];
    for (int i = 0; i < total_size; ++i) {
      shortwave[i] = 500.0 + (i % 10) * 10.0;
    }
    UrbanSetAtmShortwaveDown(urban, shortwave, size, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
  }
};

// ============================================================================
// Tests for 3D Shortwave Radiation Getters - Absorbed
// ============================================================================

TEST_F(RadiationGetterTest, GetAbsorbedShortwaveRoof_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  UrbanErrorCode status;

  UrbanGetAbsorbedShortwaveRoof(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < total_size; ++i) {
    EXPECT_FALSE(std::isnan(values[i])) << "Value at index " << i << " is NaN";
    EXPECT_GE(values[i], 0.0) << "Absorbed radiation should be non-negative";
  }
}

TEST_F(RadiationGetterTest, GetAbsorbedShortwaveRoof_SizeMismatch) {
  double values[40];
  int wrong_size[3] = {10, 2, 2};
  UrbanErrorCode status;

  UrbanGetAbsorbedShortwaveRoof(urban, values, wrong_size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(RadiationGetterTest, GetAbsorbedShortwaveRoof_NullPointer) {
  UrbanErrorCode status;

  UrbanGetAbsorbedShortwaveRoof(urban, nullptr, size, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  double values[20];
  UrbanGetAbsorbedShortwaveRoof(nullptr, values, size, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);
}

TEST_F(RadiationGetterTest, GetAbsorbedShortwaveImperviousRoad_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  UrbanErrorCode status;

  UrbanGetAbsorbedShortwaveImperviousRoad(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < total_size; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

TEST_F(RadiationGetterTest, GetAbsorbedShortwavePerviousRoad_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  UrbanErrorCode status;

  UrbanGetAbsorbedShortwavePerviousRoad(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < total_size; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

TEST_F(RadiationGetterTest, GetAbsorbedShortwaveSunlitWall_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  UrbanErrorCode status;

  UrbanGetAbsorbedShortwaveSunlitWall(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < total_size; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

TEST_F(RadiationGetterTest, GetAbsorbedShortwaveShadedWall_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  UrbanErrorCode status;

  UrbanGetAbsorbedShortwaveShadedWall(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < total_size; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

// ============================================================================
// Tests for 3D Shortwave Radiation Getters - Reflected
// ============================================================================

TEST_F(RadiationGetterTest, GetReflectedShortwaveRoof_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  UrbanErrorCode status;

  UrbanGetReflectedShortwaveRoof(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < total_size; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

TEST_F(RadiationGetterTest, GetReflectedShortwaveRoof_NullPointer) {
  UrbanErrorCode status;

  UrbanGetReflectedShortwaveRoof(urban, nullptr, size, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);
}

TEST_F(RadiationGetterTest, GetReflectedShortwaveImperviousRoad_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  UrbanErrorCode status;

  UrbanGetReflectedShortwaveImperviousRoad(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < total_size; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

TEST_F(RadiationGetterTest, GetReflectedShortwavePerviousRoad_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  UrbanErrorCode status;

  UrbanGetReflectedShortwavePerviousRoad(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < total_size; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

TEST_F(RadiationGetterTest, GetReflectedShortwaveSunlitWall_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  UrbanErrorCode status;

  UrbanGetReflectedShortwaveSunlitWall(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < total_size; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

TEST_F(RadiationGetterTest, GetReflectedShortwaveShadedWall_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  UrbanErrorCode status;

  UrbanGetReflectedShortwaveShadedWall(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < total_size; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

// ============================================================================
// Tests for 1D Longwave Radiation Getters - Net
// ============================================================================

TEST_F(RadiationGetterTest, GetNetLongwaveRoof_ValidData) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetNetLongwaveRoof(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
  }
}

TEST_F(RadiationGetterTest, GetNetLongwaveRoof_LengthMismatch) {
  double values[10];
  UrbanErrorCode status;

  UrbanGetNetLongwaveRoof(urban, values, 10, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(RadiationGetterTest, GetNetLongwaveRoof_NullPointer) {
  UrbanErrorCode status;

  UrbanGetNetLongwaveRoof(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  double values[5];
  UrbanGetNetLongwaveRoof(nullptr, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);
}

TEST_F(RadiationGetterTest, GetNetLongwaveImperviousRoad_ValidData) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetNetLongwaveImperviousRoad(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
  }
}

TEST_F(RadiationGetterTest, GetNetLongwavePerviousRoad_ValidData) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetNetLongwavePerviousRoad(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
  }
}

TEST_F(RadiationGetterTest, GetNetLongwaveSunlitWall_ValidData) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetNetLongwaveSunlitWall(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
  }
}

TEST_F(RadiationGetterTest, GetNetLongwaveShadedWall_ValidData) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetNetLongwaveShadedWall(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
  }
}

// ============================================================================
// Tests for 1D Longwave Radiation Getters - Upward
// ============================================================================

TEST_F(RadiationGetterTest, GetUpwardLongwaveRoof_ValidData) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetUpwardLongwaveRoof(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GT(values[i], 0.0) << "Upward longwave should be positive";
  }
}

TEST_F(RadiationGetterTest, GetUpwardLongwaveRoof_LengthMismatch) {
  double values[10];
  UrbanErrorCode status;

  UrbanGetUpwardLongwaveRoof(urban, values, 10, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(RadiationGetterTest, GetUpwardLongwaveImperviousRoad_ValidData) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetUpwardLongwaveImperviousRoad(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GT(values[i], 0.0);
  }
}

TEST_F(RadiationGetterTest, GetUpwardLongwavePerviousRoad_ValidData) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetUpwardLongwavePerviousRoad(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GT(values[i], 0.0);
  }
}

TEST_F(RadiationGetterTest, GetUpwardLongwaveSunlitWall_ValidData) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetUpwardLongwaveSunlitWall(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GT(values[i], 0.0);
  }
}

TEST_F(RadiationGetterTest, GetUpwardLongwaveShadedWall_ValidData) {
  double values[5];
  UrbanErrorCode status;

  UrbanGetUpwardLongwaveShadedWall(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GT(values[i], 0.0);
  }
}

// ============================================================================
// Integration Tests
// ============================================================================

TEST_F(RadiationGetterTest, MultipleGetters_Sequential) {
  UrbanErrorCode status;
  const int total_size = numLandunits * numBands * numRadTypes;
  
  double absorbed_roof[20], absorbed_impRoad[20];
  double reflected_roof[20], reflected_impRoad[20];
  double net_lw_roof[5], upward_lw_roof[5];
  
  // Call all getters in sequence
  UrbanGetAbsorbedShortwaveRoof(urban, absorbed_roof, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetAbsorbedShortwaveImperviousRoad(urban, absorbed_impRoad, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetReflectedShortwaveRoof(urban, reflected_roof, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetReflectedShortwaveImperviousRoad(urban, reflected_impRoad, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetNetLongwaveRoof(urban, net_lw_roof, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetUpwardLongwaveRoof(urban, upward_lw_roof, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
}

TEST_F(RadiationGetterTest, AllAbsorbedShortwaveGetters_Sequential) {
  UrbanErrorCode status;
  const int total_size = numLandunits * numBands * numRadTypes;
  
  double roof[20], imp_road[20], perv_road[20], sun_wall[20], shade_wall[20];
  
  UrbanGetAbsorbedShortwaveRoof(urban, roof, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetAbsorbedShortwaveImperviousRoad(urban, imp_road, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetAbsorbedShortwavePerviousRoad(urban, perv_road, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetAbsorbedShortwaveSunlitWall(urban, sun_wall, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetAbsorbedShortwaveShadedWall(urban, shade_wall, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
}

TEST_F(RadiationGetterTest, AllReflectedShortwaveGetters_Sequential) {
  UrbanErrorCode status;
  const int total_size = numLandunits * numBands * numRadTypes;
  
  double roof[20], imp_road[20], perv_road[20], sun_wall[20], shade_wall[20];
  
  UrbanGetReflectedShortwaveRoof(urban, roof, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetReflectedShortwaveImperviousRoad(urban, imp_road, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetReflectedShortwavePerviousRoad(urban, perv_road, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetReflectedShortwaveSunlitWall(urban, sun_wall, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetReflectedShortwaveShadedWall(urban, shade_wall, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
}

TEST_F(RadiationGetterTest, AllNetLongwaveGetters_Sequential) {
  UrbanErrorCode status;
  
  double roof[5], imp_road[5], perv_road[5], sun_wall[5], shade_wall[5];
  
  UrbanGetNetLongwaveRoof(urban, roof, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetNetLongwaveImperviousRoad(urban, imp_road, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetNetLongwavePerviousRoad(urban, perv_road, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetNetLongwaveSunlitWall(urban, sun_wall, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetNetLongwaveShadedWall(urban, shade_wall, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
}

TEST_F(RadiationGetterTest, AllUpwardLongwaveGetters_Sequential) {
  UrbanErrorCode status;
  
  double roof[5], imp_road[5], perv_road[5], sun_wall[5], shade_wall[5];
  
  UrbanGetUpwardLongwaveRoof(urban, roof, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetUpwardLongwaveImperviousRoad(urban, imp_road, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetUpwardLongwavePerviousRoad(urban, perv_road, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetUpwardLongwaveSunlitWall(urban, sun_wall, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  
  UrbanGetUpwardLongwaveShadedWall(urban, shade_wall, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
}

// Main function to initialize Kokkos
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  
  return result;
}
