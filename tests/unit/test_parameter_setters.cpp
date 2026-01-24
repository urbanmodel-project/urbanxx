#include "Urban.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

// Test fixture for parameter setter tests
class ParameterSetterTest : public ::testing::Test {
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

// Test: UrbanSetCanyonHwr with valid data
TEST_F(ParameterSetterTest, SetCanyonHwr_ValidData) {
  double values[10] = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5};
  UrbanErrorCode status;

  UrbanSetCanyonHwr(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetCanyonHwr should succeed with valid data";
}

// Test: UrbanSetCanyonHwr with size mismatch
TEST_F(ParameterSetterTest, SetCanyonHwr_SizeMismatch) {
  double values[5] = {1.0, 1.5, 2.0, 2.5, 3.0};
  UrbanErrorCode status;

  UrbanSetCanyonHwr(urban, values, 5, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH) 
      << "UrbanSetCanyonHwr should fail with size mismatch";
}

// Test: UrbanSetCanyonHwr with null pointer
TEST_F(ParameterSetterTest, SetCanyonHwr_NullPointer) {
  UrbanErrorCode status;

  UrbanSetCanyonHwr(urban, nullptr, numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "UrbanSetCanyonHwr should fail with null data pointer";
}

// Test: UrbanSetCanyonHwr with null urban
TEST_F(ParameterSetterTest, SetCanyonHwr_NullUrban) {
  double values[10] = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5};
  UrbanErrorCode status;

  UrbanSetCanyonHwr(nullptr, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "UrbanSetCanyonHwr should fail with null urban object";
}

// Test: UrbanSetFracPervRoadOfTotalRoad with valid data
TEST_F(ParameterSetterTest, SetFracPervRoadOfTotalRoad_ValidData) {
  double values[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  UrbanErrorCode status;

  UrbanSetFracPervRoadOfTotalRoad(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetFracPervRoadOfTotalRoad should succeed with valid data";
}

// Test: UrbanSetWtRoof with valid data
TEST_F(ParameterSetterTest, SetWtRoof_ValidData) {
  double values[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  UrbanErrorCode status;

  UrbanSetWtRoof(urban, values, numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetWtRoof should succeed with valid data";
}

// Test: Emissivity setters
TEST_F(ParameterSetterTest, SetEmissivity_ValidData) {
  double values[10] = {0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99};
  UrbanErrorCode status;

  UrbanSetEmissivityPerviousRoad(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetEmissivityPerviousRoad should succeed";

  UrbanSetEmissivityImperviousRoad(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetEmissivityImperviousRoad should succeed";

  UrbanSetEmissivityWall(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetEmissivityWall should succeed";

  UrbanSetEmissivityRoof(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetEmissivityRoof should succeed";
}

// Test: Thermal conductivity setters
TEST_F(ParameterSetterTest, SetThermalConductivity_ValidData) {
  const int numLevels = 15;
  const int totalSize = numLandunits * numLevels;
  double values[150];
  for (int i = 0; i < totalSize; ++i) {
    values[i] = 0.5 + (i % numLevels) * 0.1;
  }
  const int size[2] = {numLandunits, numLevels};
  UrbanErrorCode status;

  UrbanSetThermalConductivityRoad(urban, values, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetThermalConductivityRoad should succeed";

  UrbanSetThermalConductivityWall(urban, values, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetThermalConductivityWall should succeed";

  UrbanSetThermalConductivityRoof(urban, values, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetThermalConductivityRoof should succeed";
}

// Test: Heat capacity setters
TEST_F(ParameterSetterTest, SetHeatCapacity_ValidData) {
  const int numLevels = 15;
  const int totalSize = numLandunits * numLevels;
  double values[150];
  for (int i = 0; i < totalSize; ++i) {
    values[i] = 1000.0 + (i % numLevels) * 100.0;
  }
  const int size[2] = {numLandunits, numLevels};
  UrbanErrorCode status;

  UrbanSetHeatCapacityRoad(urban, values, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetHeatCapacityRoad should succeed";

  UrbanSetHeatCapacityWall(urban, values, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetHeatCapacityWall should succeed";

  UrbanSetHeatCapacityRoof(urban, values, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetHeatCapacityRoof should succeed";
}

// Test: Height parameter setters
TEST_F(ParameterSetterTest, SetHeightParameters_ValidData) {
  double values[10] = {10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0};
  UrbanErrorCode status;

  UrbanSetForcHgtT(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetForcHgtT should succeed";

  UrbanSetForcHgtU(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetForcHgtU should succeed";

  UrbanSetZDTown(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetZDTown should succeed";

  UrbanSetZ0Town(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetZ0Town should succeed";

  UrbanSetHtRoof(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetHtRoof should succeed";

  UrbanSetWindHgtCanyon(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetWindHgtCanyon should succeed";
}

// Test: Atmospheric forcing setters
TEST_F(ParameterSetterTest, SetAtmosphericForcing_ValidData) {
  double values[10] = {300.0, 301.0, 302.0, 303.0, 304.0, 305.0, 306.0, 307.0, 308.0, 309.0};
  UrbanErrorCode status;

  UrbanSetAtmTemp(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetAtmTemp should succeed";

  UrbanSetAtmPotTemp(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetAtmPotTemp should succeed";

  UrbanSetAtmRho(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetAtmRho should succeed";

  UrbanSetAtmSpcHumd(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetAtmSpcHumd should succeed";

  UrbanSetAtmPress(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetAtmPress should succeed";

  UrbanSetAtmWindU(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetAtmWindU should succeed";

  UrbanSetAtmWindV(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetAtmWindV should succeed";

  UrbanSetAtmCoszen(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetAtmCoszen should succeed";

  UrbanSetAtmFracSnow(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetAtmFracSnow should succeed";

  UrbanSetAtmLongwaveDown(urban, values, numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) << "UrbanSetAtmLongwaveDown should succeed";
}

// Test: Multiple setters in sequence
TEST_F(ParameterSetterTest, MultipleSetters_Sequential) {
  double canyon_hwr[10] = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5};
  double emissivity[10] = {0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99};
  double heights[10] = {10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0};
  UrbanErrorCode status;

  // Set multiple parameters in sequence
  UrbanSetCanyonHwr(urban, canyon_hwr, numLandunits, &status);
  ASSERT_EQ(status, URBAN_SUCCESS);

  UrbanSetEmissivityRoof(urban, emissivity, numLandunits, &status);
  ASSERT_EQ(status, URBAN_SUCCESS);

  UrbanSetHtRoof(urban, heights, numLandunits, &status);
  ASSERT_EQ(status, URBAN_SUCCESS);
}

// Main function to initialize Kokkos
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  
  return result;
}
