#include "Urban.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <vector>

// Test fixture for UrbanSetBuildingTemperature tests
class BuildingTemperatureSetterTest : public ::testing::Test {
protected:
  UrbanType urban;
  const int numLandunits = 10;

  void SetUp() override {
    UrbanErrorCode ierr;
    UrbanCreate(numLandunits, &urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS) << "Failed to create Urban object in SetUp";
  }

  void TearDown() override {
    UrbanErrorCode ierr;
    UrbanDestroy(&urban, &ierr);
  }
};

// =============================================================================
// Tests for UrbanSetBuildingTemperature
// =============================================================================

// Test: Set building temperature with valid data
TEST_F(BuildingTemperatureSetterTest, SetBuildingTemperature_ValidData) {
  std::vector<double> values(numLandunits, 293.15); // 20°C in Kelvin
  UrbanErrorCode status;

  UrbanSetBuildingTemperature(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetBuildingTemperature should succeed with valid data";
}

// Test: Set building temperature with a range of physically reasonable values
TEST_F(BuildingTemperatureSetterTest, SetBuildingTemperature_ReasonableValues) {
  std::vector<double> values(numLandunits);
  // Temperatures between 273 K (0°C) and 313 K (40°C)
  for (int i = 0; i < numLandunits; ++i) {
    values[i] = 273.0 + i * 4.0;
  }
  UrbanErrorCode status;

  UrbanSetBuildingTemperature(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetBuildingTemperature should succeed with varying temperature values";
}

// Test: Set building temperature with size mismatch
TEST_F(BuildingTemperatureSetterTest, SetBuildingTemperature_SizeMismatch) {
  const int wrongLength = numLandunits / 2;
  std::vector<double> values(wrongLength, 293.15);
  UrbanErrorCode status;

  UrbanSetBuildingTemperature(urban, values.data(), wrongLength, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH)
      << "UrbanSetBuildingTemperature should fail with size mismatch";
}

// Test: Set building temperature with null data pointer
TEST_F(BuildingTemperatureSetterTest, SetBuildingTemperature_NullDataPointer) {
  UrbanErrorCode status;

  UrbanSetBuildingTemperature(urban, nullptr, numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetBuildingTemperature should fail with null data pointer";
}

// Test: Set building temperature with null urban object
TEST_F(BuildingTemperatureSetterTest, SetBuildingTemperature_NullUrban) {
  std::vector<double> values(numLandunits, 293.15);
  UrbanErrorCode status;

  UrbanSetBuildingTemperature(nullptr, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetBuildingTemperature should fail with null urban object";
}

// Test: Set building temperature with null status pointer (should not crash)
TEST_F(BuildingTemperatureSetterTest, SetBuildingTemperature_NullStatus) {
  std::vector<double> values(numLandunits, 293.15);

  // Should not crash
  UrbanSetBuildingTemperature(urban, values.data(), numLandunits, nullptr);
}

// Test: Set building temperature multiple times (last write wins)
TEST_F(BuildingTemperatureSetterTest, SetBuildingTemperature_MultipleSets) {
  std::vector<double> values1(numLandunits, 280.0);
  std::vector<double> values2(numLandunits, 300.0);
  UrbanErrorCode status;

  UrbanSetBuildingTemperature(urban, values1.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanSetBuildingTemperature(urban, values2.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetBuildingTemperature should succeed on repeated calls";
}

// Test: Set building temperature with zero value
TEST_F(BuildingTemperatureSetterTest, SetBuildingTemperature_ZeroValue) {
  std::vector<double> values(numLandunits, 0.0);
  UrbanErrorCode status;

  UrbanSetBuildingTemperature(urban, values.data(), numLandunits, &status);

  // Zero is an unusual value but the API should accept any double
  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetBuildingTemperature should accept zero temperature";
}

// Test: Consistency with UrbanSetBuildingMinTemperature and
// UrbanSetBuildingMaxTemperature; both variants must succeed
TEST_F(BuildingTemperatureSetterTest,
       SetBuildingTemperature_ConsistentWithMinMax) {
  std::vector<double> minTemp(numLandunits, 285.0);
  std::vector<double> maxTemp(numLandunits, 310.0);
  std::vector<double> initTemp(numLandunits, 293.0);
  UrbanErrorCode status;

  UrbanSetBuildingMinTemperature(urban, minTemp.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanSetBuildingMaxTemperature(urban, maxTemp.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanSetBuildingTemperature(urban, initTemp.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetBuildingTemperature should succeed alongside min/max setters";
}

// Main function to initialize Kokkos
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();

  return result;
}
