#include "Urban.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

// Test fixture for 3D albedo setter tests
class AlbedoSetterTest : public ::testing::Test {
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
  }

  void TearDown() override {
    UrbanErrorCode ierr;
    UrbanDestroy(&urban, &ierr);
  }
};

// Test: UrbanSetAlbedoPerviousRoad with valid data
TEST_F(AlbedoSetterTest, SetAlbedoPerviousRoad_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20]; // 5*2*2 = 20
  
  // Fill with test albedo values (0.0 to 0.9)
  for (int i = 0; i < total_size; ++i) {
    values[i] = 0.1 + (i % 9) * 0.1;
  }

  UrbanErrorCode status;
  UrbanSetAlbedoPerviousRoad(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetAlbedoPerviousRoad should succeed with valid data";
}

// Test: UrbanSetAlbedoImperviousRoad with valid data
TEST_F(AlbedoSetterTest, SetAlbedoImperviousRoad_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  
  for (int i = 0; i < total_size; ++i) {
    values[i] = 0.15;
  }

  UrbanErrorCode status;
  UrbanSetAlbedoImperviousRoad(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetAlbedoImperviousRoad should succeed with valid data";
}

// Test: UrbanSetAlbedoSunlitWall with valid data
TEST_F(AlbedoSetterTest, SetAlbedoSunlitWall_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  
  for (int i = 0; i < total_size; ++i) {
    values[i] = 0.25;
  }

  UrbanErrorCode status;
  UrbanSetAlbedoSunlitWall(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetAlbedoSunlitWall should succeed with valid data";
}

// Test: UrbanSetAlbedoShadedWall with valid data
TEST_F(AlbedoSetterTest, SetAlbedoShadedWall_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  
  for (int i = 0; i < total_size; ++i) {
    values[i] = 0.20;
  }

  UrbanErrorCode status;
  UrbanSetAlbedoShadedWall(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetAlbedoShadedWall should succeed with valid data";
}

// Test: UrbanSetAlbedoRoof with valid data
TEST_F(AlbedoSetterTest, SetAlbedoRoof_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  
  for (int i = 0; i < total_size; ++i) {
    values[i] = 0.30;
  }

  UrbanErrorCode status;
  UrbanSetAlbedoRoof(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetAlbedoRoof should succeed with valid data";
}

// Test: Size mismatch in first dimension
TEST_F(AlbedoSetterTest, SetAlbedoRoof_SizeMismatch_Dim0) {
  const int wrong_size[3] = {numLandunits + 2, numBands, numRadTypes}; // Wrong first dimension
  const int total_size = wrong_size[0] * wrong_size[1] * wrong_size[2];
  double values[28]; // 7*2*2 = 28
  
  for (int i = 0; i < total_size; ++i) {
    values[i] = 0.30;
  }

  UrbanErrorCode status;
  UrbanSetAlbedoRoof(urban, values, wrong_size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH) 
      << "UrbanSetAlbedoRoof should fail with dimension mismatch";
}

// Test: Size mismatch in second dimension
TEST_F(AlbedoSetterTest, SetAlbedoRoof_SizeMismatch_Dim1) {
  const int wrong_size[3] = {numLandunits, numBands + 1, numRadTypes}; // Wrong second dimension
  const int total_size = wrong_size[0] * wrong_size[1] * wrong_size[2];
  double values[30]; // 5*3*2 = 30
  
  for (int i = 0; i < total_size; ++i) {
    values[i] = 0.30;
  }

  UrbanErrorCode status;
  UrbanSetAlbedoRoof(urban, values, wrong_size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH) 
      << "UrbanSetAlbedoRoof should fail with dimension mismatch";
}

// Test: Null pointer
TEST_F(AlbedoSetterTest, SetAlbedoRoof_NullPointer) {
  UrbanErrorCode status;

  UrbanSetAlbedoRoof(urban, nullptr, size, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "UrbanSetAlbedoRoof should fail with null data pointer";
}

// Test: Null size array
TEST_F(AlbedoSetterTest, SetAlbedoRoof_NullSize) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];

  UrbanErrorCode status;
  UrbanSetAlbedoRoof(urban, values, nullptr, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "UrbanSetAlbedoRoof should fail with null size array";
}

// Test: Null urban object
TEST_F(AlbedoSetterTest, SetAlbedoRoof_NullUrban) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];

  UrbanErrorCode status;
  UrbanSetAlbedoRoof(nullptr, values, size, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "UrbanSetAlbedoRoof should fail with null urban object";
}

// Test: All albedo setters in sequence
TEST_F(AlbedoSetterTest, AllAlbedoSetters_Sequential) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double perv_road[20], imperv_road[20], sunlit_wall[20], shaded_wall[20], roof[20];
  
  // Fill with different values for each surface
  for (int i = 0; i < total_size; ++i) {
    perv_road[i] = 0.15;
    imperv_road[i] = 0.10;
    sunlit_wall[i] = 0.25;
    shaded_wall[i] = 0.20;
    roof[i] = 0.30;
  }

  UrbanErrorCode status;

  UrbanSetAlbedoPerviousRoad(urban, perv_road, size, &status);
  ASSERT_EQ(status, URBAN_SUCCESS);

  UrbanSetAlbedoImperviousRoad(urban, imperv_road, size, &status);
  ASSERT_EQ(status, URBAN_SUCCESS);

  UrbanSetAlbedoSunlitWall(urban, sunlit_wall, size, &status);
  ASSERT_EQ(status, URBAN_SUCCESS);

  UrbanSetAlbedoShadedWall(urban, shaded_wall, size, &status);
  ASSERT_EQ(status, URBAN_SUCCESS);

  UrbanSetAlbedoRoof(urban, roof, size, &status);
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
