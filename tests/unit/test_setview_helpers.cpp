#include "Urban.h"
#include "private/UrbanSetterHelpers.h"
#include "private/DataTypesImpl.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

using namespace URBANXX;

// Test fixture for SetView tests
class SetViewTest : public ::testing::Test {
protected:
  void SetUp() override {
    // Kokkos is initialized in main()
  }

  void TearDown() override {
    // Cleanup if needed
  }
};

// Test: SetView1D with valid data
TEST_F(SetViewTest, SetView1D_ValidData) {
  const int length = 5;
  double input_data[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
  UrbanErrorCode status;

  // Create a 1D view
  Array1DR8 view("test_view", length);

  // Set the view using SetView1D
  SetView1D(view, input_data, length, &status);

  ASSERT_EQ(status, URBAN_SUCCESS) << "SetView1D should succeed with valid data";

  // Create a host mirror to verify the data
  auto view_host = Kokkos::create_mirror_view(view);
  Kokkos::deep_copy(view_host, view);

  // Check that all values were copied correctly
  for (int i = 0; i < length; ++i) {
    EXPECT_DOUBLE_EQ(view_host(i), input_data[i]) 
        << "Value mismatch at index " << i;
  }
}

// Test: SetView1D with size mismatch
TEST_F(SetViewTest, SetView1D_SizeMismatch) {
  const int view_length = 5;
  const int data_length = 3;
  double input_data[3] = {1.0, 2.0, 3.0};
  UrbanErrorCode status;

  // Create a 1D view with different size
  Array1DR8 view("test_view", view_length);

  // Try to set with wrong size
  SetView1D(view, input_data, data_length, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH) 
      << "SetView1D should fail with size mismatch";
}

// Test: SetView1D with null pointer
TEST_F(SetViewTest, SetView1D_NullPointer) {
  const int length = 5;
  UrbanErrorCode status;

  Array1DR8 view("test_view", length);

  // Test with null data pointer
  SetView1D(view, nullptr, length, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "SetView1D should fail with null data pointer";
}

// Test: SetView1D with null status
TEST_F(SetViewTest, SetView1D_NullStatus) {
  const int length = 5;
  double input_data[5] = {1.0, 2.0, 3.0, 4.0, 5.0};

  Array1DR8 view("test_view", length);

  // Should not crash with null status
  SetView1D(view, input_data, length, nullptr);
}

// Test: SetView3D with valid data
TEST_F(SetViewTest, SetView3D_ValidData) {
  const int size[3] = {2, 3, 4};
  const int total_size = size[0] * size[1] * size[2];
  double input_data[24]; // 2*3*4 = 24
  
  // Fill with test data
  for (int i = 0; i < total_size; ++i) {
    input_data[i] = static_cast<double>(i);
  }

  UrbanErrorCode status;

  // Create a 3D view
  Array3DR8 view("test_view", size[0], size[1], size[2]);

  // Set the view using SetView3D
  SetView3D(view, input_data, size, &status);

  ASSERT_EQ(status, URBAN_SUCCESS) << "SetView3D should succeed with valid data";

  // Create a host mirror to verify the data
  auto view_host = Kokkos::create_mirror_view(view);
  Kokkos::deep_copy(view_host, view);

  // Check that all values were copied correctly
  // Determine memory layout and iterate accordingly
  using layout_type = typename Array3DR8::array_layout;
  int idx = 0;
  
  if (std::is_same<layout_type, Kokkos::LayoutLeft>::value) {
    // LayoutLeft = column-major: stride-1 on first index
    for (int k = 0; k < size[2]; ++k) {
      for (int j = 0; j < size[1]; ++j) {
        for (int i = 0; i < size[0]; ++i) {
          EXPECT_DOUBLE_EQ(view_host(i, j, k), input_data[idx]) 
              << "Value mismatch at (" << i << "," << j << "," << k << ") with LayoutLeft";
          idx++;
        }
      }
    }
  } else if (std::is_same<layout_type, Kokkos::LayoutRight>::value) {
    // LayoutRight = row-major: stride-1 on last index
    for (int i = 0; i < size[0]; ++i) {
      for (int j = 0; j < size[1]; ++j) {
        for (int k = 0; k < size[2]; ++k) {
          EXPECT_DOUBLE_EQ(view_host(i, j, k), input_data[idx]) 
              << "Value mismatch at (" << i << "," << j << "," << k << ") with LayoutRight";
          idx++;
        }
      }
    }
  } else {
    // For other layouts, check element-wise without assuming order
    FAIL() << "Unsupported layout type for this test";
  }
}

// Test: SetView3D with size mismatch in first dimension
TEST_F(SetViewTest, SetView3D_SizeMismatch_Dim0) {
  const int view_size[3] = {2, 3, 4};
  const int data_size[3] = {3, 3, 4}; // Wrong first dimension
  double input_data[36]; // 3*3*4 = 36

  UrbanErrorCode status;
  Array3DR8 view("test_view", view_size[0], view_size[1], view_size[2]);

  SetView3D(view, input_data, data_size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH) 
      << "SetView3D should fail with dimension mismatch";
}

// Test: SetView3D with size mismatch in second dimension
TEST_F(SetViewTest, SetView3D_SizeMismatch_Dim1) {
  const int view_size[3] = {2, 3, 4};
  const int data_size[3] = {2, 5, 4}; // Wrong second dimension
  double input_data[40]; // 2*5*4 = 40

  UrbanErrorCode status;
  Array3DR8 view("test_view", view_size[0], view_size[1], view_size[2]);

  SetView3D(view, input_data, data_size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH) 
      << "SetView3D should fail with dimension mismatch";
}

// Test: SetView3D with size mismatch in third dimension
TEST_F(SetViewTest, SetView3D_SizeMismatch_Dim2) {
  const int view_size[3] = {2, 3, 4};
  const int data_size[3] = {2, 3, 6}; // Wrong third dimension
  double input_data[36]; // 2*3*6 = 36

  UrbanErrorCode status;
  Array3DR8 view("test_view", view_size[0], view_size[1], view_size[2]);

  SetView3D(view, input_data, data_size, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH) 
      << "SetView3D should fail with dimension mismatch";
}

// Test: SetView3D with null pointer
TEST_F(SetViewTest, SetView3D_NullPointer) {
  const int size[3] = {2, 3, 4};
  UrbanErrorCode status;

  Array3DR8 view("test_view", size[0], size[1], size[2]);

  // Test with null data pointer
  SetView3D(view, nullptr, size, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "SetView3D should fail with null data pointer";
}

// Test: SetView3D with null size array
TEST_F(SetViewTest, SetView3D_NullSize) {
  double input_data[24];
  UrbanErrorCode status;

  Array3DR8 view("test_view", 2, 3, 4);

  // Test with null size array
  SetView3D(view, input_data, nullptr, &status);

  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT) 
      << "SetView3D should fail with null size array";
}

// Test: SetView3D with null status
TEST_F(SetViewTest, SetView3D_NullStatus) {
  const int size[3] = {2, 3, 4};
  double input_data[24];

  Array3DR8 view("test_view", size[0], size[1], size[2]);

  // Should not crash with null status
  SetView3D(view, input_data, size, nullptr);
}

// Main function to initialize Kokkos
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  
  return result;
}
