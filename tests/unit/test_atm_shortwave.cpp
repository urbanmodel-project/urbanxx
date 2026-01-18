#include "Urban.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

// Test fixture for atmospheric shortwave setter tests
class AtmShortwaveTest : public ::testing::Test {
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

// Test: UrbanSetAtmShortwaveDown with valid data
TEST_F(AtmShortwaveTest, SetAtmShortwaveDown_ValidData) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20]; // 5*2*2 = 20
  
  // Fill with test shortwave radiation values (W/m^2)
  for (int i = 0; i < total_size; ++i) {
    values[i] = 400.0 + i * 10.0; // Range from 400 to 590 W/m^2
  }

  UrbanErrorCode status;
  UrbanSetAtmShortwaveDown(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetAtmShortwaveDown should succeed with valid data";
}

// Test: UrbanSetAtmShortwaveDown with zeros (valid - nighttime)
TEST_F(AtmShortwaveTest, SetAtmShortwaveDown_ZeroValues) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  
  // All zeros (nighttime)
  for (int i = 0; i < total_size; ++i) {
    values[i] = 0.0;
  }

  UrbanErrorCode status;
  UrbanSetAtmShortwaveDown(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetAtmShortwaveDown should succeed with zero values (nighttime)";
}

// Test: UrbanSetAtmShortwaveDown with high values
TEST_F(AtmShortwaveTest, SetAtmShortwaveDown_HighValues) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  
  // High values (intense sunlight)
  for (int i = 0; i < total_size; ++i) {
    values[i] = 1000.0; // ~1000 W/m^2 is typical peak
  }

  UrbanErrorCode status;
  UrbanSetAtmShortwaveDown(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetAtmShortwaveDown should succeed with high values";
}

// Test: UrbanSetAtmShortwaveDown with varying values across bands
TEST_F(AtmShortwaveTest, SetAtmShortwaveDown_VaryingBands) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  
  // Different values for visible vs near-infrared bands
  int idx = 0;
  for (int k = 0; k < numRadTypes; ++k) {
    for (int j = 0; j < numBands; ++j) {
      for (int i = 0; i < numLandunits; ++i) {
        if (j == 0) {
          values[idx] = 300.0; // Visible band
        } else {
          values[idx] = 500.0; // Near-infrared band
        }
        idx++;
      }
    }
  }

  UrbanErrorCode status;
  UrbanSetAtmShortwaveDown(urban, values, size, &status);

  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetAtmShortwaveDown should succeed with varying band values";
}

// Test: Size mismatch in first dimension
TEST_F(AtmShortwaveTest, SetAtmShortwaveDown_SizeMismatch_Dim0) {
  const int wrong_size[3] = {numLandunits + 2, numBands, numRadTypes};
  const int total_size = wrong_size[0] * wrong_size[1] * wrong_size[2];
  double values[28]; // 7*2*2 = 28
  
  for (int i = 0; i < total_size; ++i) {
    values[i] = 400.0;
  }

  UrbanErrorCode status;
  UrbanSetAtmShortwaveDown(urban, values, wrong_size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH) 
      << "UrbanSetAtmShortwaveDown should fail with dimension mismatch";
}

// Test: Size mismatch in second dimension (bands)
TEST_F(AtmShortwaveTest, SetAtmShortwaveDown_SizeMismatch_Dim1) {
  const int wrong_size[3] = {numLandunits, numBands + 1, numRadTypes};
  const int total_size = wrong_size[0] * wrong_size[1] * wrong_size[2];
  double values[30]; // 5*3*2 = 30
  
  for (int i = 0; i < total_size; ++i) {
    values[i] = 400.0;
  }

  UrbanErrorCode status;
  UrbanSetAtmShortwaveDown(urban, values, wrong_size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH) 
      << "UrbanSetAtmShortwaveDown should fail with band dimension mismatch";
}

// Test: Size mismatch in third dimension (radiation types)
TEST_F(AtmShortwaveTest, SetAtmShortwaveDown_SizeMismatch_Dim2) {
  const int wrong_size[3] = {numLandunits, numBands, numRadTypes + 1};
  const int total_size = wrong_size[0] * wrong_size[1] * wrong_size[2];
  double values[30]; // 5*2*3 = 30
  
  for (int i = 0; i < total_size; ++i) {
    values[i] = 400.0;
  }

  UrbanErrorCode status;
  UrbanSetAtmShortwaveDown(urban, values, wrong_size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH) 
      << "UrbanSetAtmShortwaveDown should fail with radiation type dimension mismatch";
}

// Test: Null data pointer
TEST_F(AtmShortwaveTest, SetAtmShortwaveDown_NullPointer) {
  UrbanErrorCode status;

  UrbanSetAtmShortwaveDown(urban, nullptr, size, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "UrbanSetAtmShortwaveDown should fail with null data pointer";
}

// Test: Null size array
TEST_F(AtmShortwaveTest, SetAtmShortwaveDown_NullSize) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  
  for (int i = 0; i < total_size; ++i) {
    values[i] = 400.0;
  }

  UrbanErrorCode status;
  UrbanSetAtmShortwaveDown(urban, values, nullptr, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "UrbanSetAtmShortwaveDown should fail with null size array";
}

// Test: Null urban object
TEST_F(AtmShortwaveTest, SetAtmShortwaveDown_NullUrban) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  
  for (int i = 0; i < total_size; ++i) {
    values[i] = 400.0;
  }

  UrbanErrorCode status;
  UrbanSetAtmShortwaveDown(nullptr, values, size, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "UrbanSetAtmShortwaveDown should fail with null urban object";
}

// Test: Null status pointer
TEST_F(AtmShortwaveTest, SetAtmShortwaveDown_NullStatus) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values[20];
  
  for (int i = 0; i < total_size; ++i) {
    values[i] = 400.0;
  }

  // Should not crash with null status
  UrbanSetAtmShortwaveDown(urban, values, size, nullptr);
}

// Test: Set shortwave multiple times (overwrite behavior)
TEST_F(AtmShortwaveTest, SetAtmShortwaveDown_MultipleTimes) {
  const int total_size = numLandunits * numBands * numRadTypes;
  double values1[20], values2[20];
  
  // First set
  for (int i = 0; i < total_size; ++i) {
    values1[i] = 400.0;
  }
  
  UrbanErrorCode status;
  UrbanSetAtmShortwaveDown(urban, values1, size, &status);
  ASSERT_EQ(status, URBAN_SUCCESS);
  
  // Second set (overwrite)
  for (int i = 0; i < total_size; ++i) {
    values2[i] = 600.0;
  }
  
  UrbanSetAtmShortwaveDown(urban, values2, size, &status);
  EXPECT_EQ(status, URBAN_SUCCESS) 
      << "UrbanSetAtmShortwaveDown should succeed when called multiple times";
}

// Main function to initialize Kokkos
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  
  return result;
}
