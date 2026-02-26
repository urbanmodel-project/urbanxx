#include "Urban.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <cmath>
#include <vector>

// Test fixture for net shortwave radiation getter tests.
// These 1D getters read per-landunit net shortwave values computed by
// UrbanComputeNetShortwaveRadiation.
class NetShortwaveGetterTest : public ::testing::Test {
protected:
  UrbanType urban;
  const int numLandunits = 5;
  const int numBands = 2;
  const int numRadTypes = 2;
  int size3D[3];

  void SetUp() override {
    UrbanErrorCode ierr;
    UrbanCreate(numLandunits, &urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS) << "Failed to create Urban object in SetUp";

    size3D[0] = numLandunits;
    size3D[1] = numBands;
    size3D[2] = numRadTypes;

    SetTestParameters();
    SetTestAtmosphericForcing();

    // Compute net shortwave radiation so getter data is available
    UrbanComputeNetShortwaveRadiation(urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS)
        << "Failed to compute net shortwave radiation";
  }

  void TearDown() override {
    UrbanErrorCode ierr;
    UrbanDestroy(&urban, &ierr);
  }

  void SetTestParameters() {
    UrbanErrorCode ierr;

    double canyon_hwr[5] = {1.0, 1.5, 2.0, 2.5, 3.0};
    UrbanSetCanyonHwr(urban, canyon_hwr, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);

    double frac_perv[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
    UrbanSetFracPervRoadOfTotalRoad(urban, frac_perv, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);

    double wt_roof[5] = {0.3, 0.35, 0.4, 0.45, 0.5};
    UrbanSetWtRoof(urban, wt_roof, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);

    // Albedos (3D: landunits x bands x rad_types)
    const int total3D = numLandunits * numBands * numRadTypes;
    std::vector<double> albedo_roof(total3D), albedo_wall(total3D),
        albedo_road(total3D);
    for (int i = 0; i < total3D; ++i) {
      albedo_roof[i] = 0.15 + (i % 5) * 0.01;
      albedo_wall[i] = 0.25 + (i % 5) * 0.01;
      albedo_road[i] = 0.10 + (i % 5) * 0.01;
    }
    UrbanSetAlbedoRoof(urban, albedo_roof.data(), size3D, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAlbedoSunlitWall(urban, albedo_wall.data(), size3D, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAlbedoShadedWall(urban, albedo_wall.data(), size3D, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAlbedoImperviousRoad(urban, albedo_road.data(), size3D, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanSetAlbedoPerviousRoad(urban, albedo_road.data(), size3D, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
  }

  void SetTestAtmosphericForcing() {
    UrbanErrorCode ierr;

    double coszen[5] = {0.5, 0.6, 0.7, 0.8, 0.9};
    UrbanSetAtmCoszen(urban, coszen, numLandunits, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);

    // Shortwave down (3D)
    const int total3D = numLandunits * numBands * numRadTypes;
    std::vector<double> shortwave(total3D);
    for (int i = 0; i < total3D; ++i) {
      shortwave[i] = 300.0 + (i % 10) * 20.0;
    }
    UrbanSetAtmShortwaveDown(urban, shortwave.data(), size3D, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
  }
};

// =============================================================================
// Tests for UrbanGetNetShortwaveRoof
// =============================================================================

TEST_F(NetShortwaveGetterTest, GetNetShortwaveRoof_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetNetShortwaveRoof(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i])) << "Value at index " << i << " is NaN";
    EXPECT_GE(values[i], 0.0) << "Net shortwave should be non-negative";
  }
}

TEST_F(NetShortwaveGetterTest, GetNetShortwaveRoof_SizeMismatch) {
  std::vector<double> values(numLandunits * 2);
  UrbanErrorCode status;

  UrbanGetNetShortwaveRoof(urban, values.data(), numLandunits * 2, &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(NetShortwaveGetterTest, GetNetShortwaveRoof_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetNetShortwaveRoof(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetNetShortwaveRoof(nullptr, values.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  // Null status should not crash
  UrbanGetNetShortwaveRoof(urban, values.data(), numLandunits, nullptr);
}

// =============================================================================
// Tests for UrbanGetNetShortwaveImperviousRoad
// =============================================================================

TEST_F(NetShortwaveGetterTest, GetNetShortwaveImperviousRoad_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetNetShortwaveImperviousRoad(urban, values.data(), numLandunits,
                                     &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

TEST_F(NetShortwaveGetterTest, GetNetShortwaveImperviousRoad_SizeMismatch) {
  std::vector<double> values(numLandunits * 2);
  UrbanErrorCode status;

  UrbanGetNetShortwaveImperviousRoad(urban, values.data(), numLandunits * 2,
                                     &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(NetShortwaveGetterTest, GetNetShortwaveImperviousRoad_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetNetShortwaveImperviousRoad(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetNetShortwaveImperviousRoad(nullptr, values.data(), numLandunits,
                                     &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  // Null status should not crash
  UrbanGetNetShortwaveImperviousRoad(urban, values.data(), numLandunits,
                                     nullptr);
}

// =============================================================================
// Tests for UrbanGetNetShortwavePerviousRoad
// =============================================================================

TEST_F(NetShortwaveGetterTest, GetNetShortwavePerviousRoad_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetNetShortwavePerviousRoad(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

TEST_F(NetShortwaveGetterTest, GetNetShortwavePerviousRoad_SizeMismatch) {
  std::vector<double> values(numLandunits * 2);
  UrbanErrorCode status;

  UrbanGetNetShortwavePerviousRoad(urban, values.data(), numLandunits * 2,
                                   &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(NetShortwaveGetterTest, GetNetShortwavePerviousRoad_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetNetShortwavePerviousRoad(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetNetShortwavePerviousRoad(nullptr, values.data(), numLandunits,
                                   &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  // Null status should not crash
  UrbanGetNetShortwavePerviousRoad(urban, values.data(), numLandunits, nullptr);
}

// =============================================================================
// Tests for UrbanGetNetShortwaveSunlitWall
// =============================================================================

TEST_F(NetShortwaveGetterTest, GetNetShortwaveSunlitWall_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetNetShortwaveSunlitWall(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

TEST_F(NetShortwaveGetterTest, GetNetShortwaveSunlitWall_SizeMismatch) {
  std::vector<double> values(numLandunits * 2);
  UrbanErrorCode status;

  UrbanGetNetShortwaveSunlitWall(urban, values.data(), numLandunits * 2,
                                  &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(NetShortwaveGetterTest, GetNetShortwaveSunlitWall_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetNetShortwaveSunlitWall(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetNetShortwaveSunlitWall(nullptr, values.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  // Null status should not crash
  UrbanGetNetShortwaveSunlitWall(urban, values.data(), numLandunits, nullptr);
}

// =============================================================================
// Tests for UrbanGetNetShortwaveShadedWall
// =============================================================================

TEST_F(NetShortwaveGetterTest, GetNetShortwaveShadedWall_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetNetShortwaveShadedWall(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_FALSE(std::isnan(values[i]));
    EXPECT_GE(values[i], 0.0);
  }
}

TEST_F(NetShortwaveGetterTest, GetNetShortwaveShadedWall_SizeMismatch) {
  std::vector<double> values(numLandunits * 2);
  UrbanErrorCode status;

  UrbanGetNetShortwaveShadedWall(urban, values.data(), numLandunits * 2,
                                  &status);

  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(NetShortwaveGetterTest, GetNetShortwaveShadedWall_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetNetShortwaveShadedWall(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetNetShortwaveShadedWall(nullptr, values.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  // Null status should not crash
  UrbanGetNetShortwaveShadedWall(urban, values.data(), numLandunits, nullptr);
}

// =============================================================================
// Integration Tests
// =============================================================================

TEST_F(NetShortwaveGetterTest, AllNetShortwaveGetters_Sequential) {
  UrbanErrorCode status;
  std::vector<double> roof(numLandunits), imp_road(numLandunits),
      perv_road(numLandunits), sun_wall(numLandunits), sha_wall(numLandunits);

  UrbanGetNetShortwaveRoof(urban, roof.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanGetNetShortwaveImperviousRoad(urban, imp_road.data(), numLandunits,
                                     &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanGetNetShortwavePerviousRoad(urban, perv_road.data(), numLandunits,
                                   &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanGetNetShortwaveSunlitWall(urban, sun_wall.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanGetNetShortwaveShadedWall(urban, sha_wall.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
}

TEST_F(NetShortwaveGetterTest, MultipleCallsSameGetter_ConsistentResults) {
  UrbanErrorCode status;
  std::vector<double> vals1(numLandunits), vals2(numLandunits);

  UrbanGetNetShortwaveRoof(urban, vals1.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanGetNetShortwaveRoof(urban, vals2.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  for (int i = 0; i < numLandunits; ++i) {
    EXPECT_DOUBLE_EQ(vals1[i], vals2[i])
        << "Repeated getter calls should return identical values";
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
