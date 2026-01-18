#include "Urban.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

// Test fixture for Urban tests
class UrbanBasicTest : public ::testing::Test {
protected:
  void SetUp() override {
    // Kokkos is initialized in main()
  }

  void TearDown() override {
    // Cleanup if needed
  }
};

// Test: Create and destroy Urban object
TEST_F(UrbanBasicTest, CreateAndDestroy) {
  UrbanType urban = nullptr;
  UrbanErrorCode ierr;
  const int numLandunits = 10;

  // Create Urban object
  UrbanCreate(numLandunits, &urban, &ierr);
  
  ASSERT_EQ(ierr, URBAN_SUCCESS) << "UrbanCreate should succeed";
  ASSERT_NE(urban, nullptr) << "Urban object should not be null after creation";

  // Destroy Urban object
  UrbanDestroy(&urban, &ierr);
  
  ASSERT_EQ(ierr, URBAN_SUCCESS) << "UrbanDestroy should succeed";
  ASSERT_EQ(urban, nullptr) << "Urban object should be null after destruction";
}

// Test: Create with invalid arguments
TEST_F(UrbanBasicTest, CreateWithInvalidArguments) {
  UrbanType urban = nullptr;
  UrbanErrorCode ierr;

  // Test with zero landunits
  UrbanCreate(0, &urban, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT) << "Should fail with zero landunits";

  // Test with negative landunits
  UrbanCreate(-5, &urban, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT) << "Should fail with negative landunits";

  // Test with null status pointer
  UrbanCreate(10, &urban, nullptr);
  // Should not crash, just checking it handles null status gracefully
}

// Test: Destroy with invalid arguments
TEST_F(UrbanBasicTest, DestroyWithInvalidArguments) {
  UrbanErrorCode ierr;

  // Test with null urban pointer
  UrbanDestroy(nullptr, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT) << "Should fail with null urban pointer";

  // Test with null urban object
  UrbanType urban = nullptr;
  UrbanDestroy(&urban, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT) << "Should fail with already null urban";
}

// Main function to initialize Kokkos
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  
  return result;
}
