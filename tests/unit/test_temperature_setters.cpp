#include "Urban.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <vector>

// Test fixture for temperature setter tests
class TemperatureSetterTest : public ::testing::Test {
protected:
  UrbanType urban;
  const int numLandunits = 10;
  const int numUrbanLayers = 5;
  const int numSoilLayers = 15;

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
// Surface Temperature Setter Tests (1D)
// =============================================================================

// Test: UrbanSetEffectiveSurfTempRoof with valid data
TEST_F(TemperatureSetterTest, SetEffectiveSurfTempRoof_ValidData) {
  std::vector<double> values(numLandunits, 292.0);
  UrbanErrorCode status;

  UrbanSetEffectiveSurfTempRoof(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetEffectiveSurfTempRoof should succeed with valid data";
}

// Test: UrbanSetEffectiveSurfTempRoof with size mismatch
TEST_F(TemperatureSetterTest, SetEffectiveSurfTempRoof_SizeMismatch) {
  std::vector<double> values(5, 292.0);
  UrbanErrorCode status;

  UrbanSetEffectiveSurfTempRoof(urban, values.data(), 5, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH)
      << "UrbanSetEffectiveSurfTempRoof should fail with size mismatch";
}

// Test: UrbanSetEffectiveSurfTempRoof with null pointer
TEST_F(TemperatureSetterTest, SetEffectiveSurfTempRoof_NullPointer) {
  UrbanErrorCode status;

  UrbanSetEffectiveSurfTempRoof(urban, nullptr, numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetEffectiveSurfTempRoof should fail with null data pointer";
}

// Test: UrbanSetEffectiveSurfTempRoof with null urban
TEST_F(TemperatureSetterTest, SetEffectiveSurfTempRoof_NullUrban) {
  std::vector<double> values(numLandunits, 292.0);
  UrbanErrorCode status;

  UrbanSetEffectiveSurfTempRoof(nullptr, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetEffectiveSurfTempRoof should fail with null urban object";
}

// Test: UrbanSetEffectiveSurfTempImperviousRoad with valid data
TEST_F(TemperatureSetterTest, SetEffectiveSurfTempImperviousRoad_ValidData) {
  std::vector<double> values(numLandunits, 274.0);
  UrbanErrorCode status;

  UrbanSetEffectiveSurfTempImperviousRoad(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetEffectiveSurfTempImperviousRoad should succeed with valid data";
}

// Test: UrbanSetEffectiveSurfTempImperviousRoad with null pointer
TEST_F(TemperatureSetterTest, SetEffectiveSurfTempImperviousRoad_NullPointer) {
  UrbanErrorCode status;

  UrbanSetEffectiveSurfTempImperviousRoad(urban, nullptr, numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetEffectiveSurfTempImperviousRoad should fail with null data pointer";
}

// Test: UrbanSetEffectiveSurfTempPerviousRoad with valid data
TEST_F(TemperatureSetterTest, SetEffectiveSurfTempPerviousRoad_ValidData) {
  std::vector<double> values(numLandunits, 274.0);
  UrbanErrorCode status;

  UrbanSetEffectiveSurfTempPerviousRoad(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetEffectiveSurfTempPerviousRoad should succeed with valid data";
}

// Test: UrbanSetEffectiveSurfTempPerviousRoad with size mismatch
TEST_F(TemperatureSetterTest, SetEffectiveSurfTempPerviousRoad_SizeMismatch) {
  std::vector<double> values(5, 274.0);
  UrbanErrorCode status;

  UrbanSetEffectiveSurfTempPerviousRoad(urban, values.data(), 5, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH)
      << "UrbanSetEffectiveSurfTempPerviousRoad should fail with size mismatch";
}

// Test: UrbanSetEffectiveSurfTempSunlitWall with valid data
TEST_F(TemperatureSetterTest, SetEffectiveSurfTempSunlitWall_ValidData) {
  std::vector<double> values(numLandunits, 292.0);
  UrbanErrorCode status;

  UrbanSetEffectiveSurfTempSunlitWall(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetEffectiveSurfTempSunlitWall should succeed with valid data";
}

// Test: UrbanSetEffectiveSurfTempSunlitWall with null urban
TEST_F(TemperatureSetterTest, SetEffectiveSurfTempSunlitWall_NullUrban) {
  std::vector<double> values(numLandunits, 292.0);
  UrbanErrorCode status;

  UrbanSetEffectiveSurfTempSunlitWall(nullptr, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetEffectiveSurfTempSunlitWall should fail with null urban object";
}

// Test: UrbanSetEffectiveSurfTempShadedWall with valid data
TEST_F(TemperatureSetterTest, SetEffectiveSurfTempShadedWall_ValidData) {
  std::vector<double> values(numLandunits, 292.0);
  UrbanErrorCode status;

  UrbanSetEffectiveSurfTempShadedWall(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetEffectiveSurfTempShadedWall should succeed with valid data";
}

// Test: UrbanSetEffectiveSurfTempShadedWall with null pointer
TEST_F(TemperatureSetterTest, SetEffectiveSurfTempShadedWall_NullPointer) {
  UrbanErrorCode status;

  UrbanSetEffectiveSurfTempShadedWall(urban, nullptr, numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetEffectiveSurfTempShadedWall should fail with null data pointer";
}

// =============================================================================
// Layer Temperature Setter Tests (2D)
// =============================================================================

// Test: UrbanSetLayerTempRoof with valid data
TEST_F(TemperatureSetterTest, SetLayerTempRoof_ValidData) {
  std::vector<double> values(numLandunits * numUrbanLayers, 292.0);
  int size[2] = {numLandunits, numUrbanLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempRoof(urban, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetLayerTempRoof should succeed with valid data";
}

// Test: UrbanSetLayerTempRoof with size mismatch (wrong number of landunits)
TEST_F(TemperatureSetterTest, SetLayerTempRoof_SizeMismatchLandunits) {
  std::vector<double> values(5 * numUrbanLayers, 292.0);
  int size[2] = {5, numUrbanLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempRoof(urban, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH)
      << "UrbanSetLayerTempRoof should fail with wrong number of landunits";
}

// Test: UrbanSetLayerTempRoof with size mismatch (wrong number of layers)
TEST_F(TemperatureSetterTest, SetLayerTempRoof_SizeMismatchLayers) {
  std::vector<double> values(numLandunits * 3, 292.0);
  int size[2] = {numLandunits, 3};
  UrbanErrorCode status;

  UrbanSetLayerTempRoof(urban, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH)
      << "UrbanSetLayerTempRoof should fail with wrong number of layers";
}

// Test: UrbanSetLayerTempRoof with null pointer
TEST_F(TemperatureSetterTest, SetLayerTempRoof_NullPointer) {
  int size[2] = {numLandunits, numUrbanLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempRoof(urban, nullptr, size, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetLayerTempRoof should fail with null data pointer";
}

// Test: UrbanSetLayerTempRoof with null urban
TEST_F(TemperatureSetterTest, SetLayerTempRoof_NullUrban) {
  std::vector<double> values(numLandunits * numUrbanLayers, 292.0);
  int size[2] = {numLandunits, numUrbanLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempRoof(nullptr, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetLayerTempRoof should fail with null urban object";
}

// Test: UrbanSetLayerTempImperviousRoad with valid data
TEST_F(TemperatureSetterTest, SetLayerTempImperviousRoad_ValidData) {
  std::vector<double> values(numLandunits * numSoilLayers, 274.0);
  int size[2] = {numLandunits, numSoilLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempImperviousRoad(urban, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetLayerTempImperviousRoad should succeed with valid data";
}

// Test: UrbanSetLayerTempImperviousRoad with size mismatch
TEST_F(TemperatureSetterTest, SetLayerTempImperviousRoad_SizeMismatch) {
  std::vector<double> values(numLandunits * 10, 274.0);
  int size[2] = {numLandunits, 10};
  UrbanErrorCode status;

  UrbanSetLayerTempImperviousRoad(urban, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH)
      << "UrbanSetLayerTempImperviousRoad should fail with wrong number of layers";
}

// Test: UrbanSetLayerTempImperviousRoad with null pointer
TEST_F(TemperatureSetterTest, SetLayerTempImperviousRoad_NullPointer) {
  int size[2] = {numLandunits, numSoilLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempImperviousRoad(urban, nullptr, size, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetLayerTempImperviousRoad should fail with null data pointer";
}

// Test: UrbanSetLayerTempPerviousRoad with valid data
TEST_F(TemperatureSetterTest, SetLayerTempPerviousRoad_ValidData) {
  std::vector<double> values(numLandunits * numSoilLayers, 274.0);
  int size[2] = {numLandunits, numSoilLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempPerviousRoad(urban, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetLayerTempPerviousRoad should succeed with valid data";
}

// Test: UrbanSetLayerTempPerviousRoad with null urban
TEST_F(TemperatureSetterTest, SetLayerTempPerviousRoad_NullUrban) {
  std::vector<double> values(numLandunits * numSoilLayers, 274.0);
  int size[2] = {numLandunits, numSoilLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempPerviousRoad(nullptr, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetLayerTempPerviousRoad should fail with null urban object";
}

// Test: UrbanSetLayerTempSunlitWall with valid data
TEST_F(TemperatureSetterTest, SetLayerTempSunlitWall_ValidData) {
  std::vector<double> values(numLandunits * numUrbanLayers, 292.0);
  int size[2] = {numLandunits, numUrbanLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempSunlitWall(urban, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetLayerTempSunlitWall should succeed with valid data";
}

// Test: UrbanSetLayerTempSunlitWall with size mismatch
TEST_F(TemperatureSetterTest, SetLayerTempSunlitWall_SizeMismatch) {
  std::vector<double> values(5 * numUrbanLayers, 292.0);
  int size[2] = {5, numUrbanLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempSunlitWall(urban, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH)
      << "UrbanSetLayerTempSunlitWall should fail with size mismatch";
}

// Test: UrbanSetLayerTempSunlitWall with null pointer
TEST_F(TemperatureSetterTest, SetLayerTempSunlitWall_NullPointer) {
  int size[2] = {numLandunits, numUrbanLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempSunlitWall(urban, nullptr, size, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetLayerTempSunlitWall should fail with null data pointer";
}

// Test: UrbanSetLayerTempShadedWall with valid data
TEST_F(TemperatureSetterTest, SetLayerTempShadedWall_ValidData) {
  std::vector<double> values(numLandunits * numUrbanLayers, 292.0);
  int size[2] = {numLandunits, numUrbanLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempShadedWall(urban, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetLayerTempShadedWall should succeed with valid data";
}

// Test: UrbanSetLayerTempShadedWall with size mismatch
TEST_F(TemperatureSetterTest, SetLayerTempShadedWall_SizeMismatch) {
  std::vector<double> values(numLandunits * 3, 292.0);
  int size[2] = {numLandunits, 3};
  UrbanErrorCode status;

  UrbanSetLayerTempShadedWall(urban, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH)
      << "UrbanSetLayerTempShadedWall should fail with wrong number of layers";
}

// Test: UrbanSetLayerTempShadedWall with null urban
TEST_F(TemperatureSetterTest, SetLayerTempShadedWall_NullUrban) {
  std::vector<double> values(numLandunits * numUrbanLayers, 292.0);
  int size[2] = {numLandunits, numUrbanLayers};
  UrbanErrorCode status;

  UrbanSetLayerTempShadedWall(nullptr, values.data(), size, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetLayerTempShadedWall should fail with null urban object";
}

// =============================================================================
// Canyon Air Property Setter Tests (1D)
// =============================================================================

// Test: UrbanSetCanyonAirTemperature with valid data
TEST_F(TemperatureSetterTest, SetCanyonAirTemperature_ValidData) {
  std::vector<double> values(numLandunits, 283.0);
  UrbanErrorCode status;

  UrbanSetCanyonAirTemperature(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetCanyonAirTemperature should succeed with valid data";
}

// Test: UrbanSetCanyonAirTemperature with size mismatch
TEST_F(TemperatureSetterTest, SetCanyonAirTemperature_SizeMismatch) {
  std::vector<double> values(5, 283.0);
  UrbanErrorCode status;

  UrbanSetCanyonAirTemperature(urban, values.data(), 5, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH)
      << "UrbanSetCanyonAirTemperature should fail with size mismatch";
}

// Test: UrbanSetCanyonAirTemperature with null pointer
TEST_F(TemperatureSetterTest, SetCanyonAirTemperature_NullPointer) {
  UrbanErrorCode status;

  UrbanSetCanyonAirTemperature(urban, nullptr, numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetCanyonAirTemperature should fail with null data pointer";
}

// Test: UrbanSetCanyonAirTemperature with null urban
TEST_F(TemperatureSetterTest, SetCanyonAirTemperature_NullUrban) {
  std::vector<double> values(numLandunits, 283.0);
  UrbanErrorCode status;

  UrbanSetCanyonAirTemperature(nullptr, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetCanyonAirTemperature should fail with null urban object";
}

// Test: UrbanSetCanyonSpecificHumidity with valid data
TEST_F(TemperatureSetterTest, SetCanyonSpecificHumidity_ValidData) {
  std::vector<double> values(numLandunits, 1.e-4);
  UrbanErrorCode status;

  UrbanSetCanyonSpecificHumidity(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS)
      << "UrbanSetCanyonSpecificHumidity should succeed with valid data";
}

// Test: UrbanSetCanyonSpecificHumidity with size mismatch
TEST_F(TemperatureSetterTest, SetCanyonSpecificHumidity_SizeMismatch) {
  std::vector<double> values(5, 1.e-4);
  UrbanErrorCode status;

  UrbanSetCanyonSpecificHumidity(urban, values.data(), 5, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH)
      << "UrbanSetCanyonSpecificHumidity should fail with size mismatch";
}

// Test: UrbanSetCanyonSpecificHumidity with null pointer
TEST_F(TemperatureSetterTest, SetCanyonSpecificHumidity_NullPointer) {
  UrbanErrorCode status;

  UrbanSetCanyonSpecificHumidity(urban, nullptr, numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetCanyonSpecificHumidity should fail with null data pointer";
}

// Test: UrbanSetCanyonSpecificHumidity with null urban
TEST_F(TemperatureSetterTest, SetCanyonSpecificHumidity_NullUrban) {
  std::vector<double> values(numLandunits, 1.e-4);
  UrbanErrorCode status;

  UrbanSetCanyonSpecificHumidity(nullptr, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT)
      << "UrbanSetCanyonSpecificHumidity should fail with null urban object";
}

// =============================================================================
// Main function
// =============================================================================
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
