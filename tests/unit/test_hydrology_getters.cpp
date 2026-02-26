#include "Urban.h"
#include "test_full_model_setup.h"
#include <Kokkos_Core.hpp>
#include <cmath>
#include <vector>

class HydrologyGetterTest : public ::testing::Test {
protected:
  UrbanType urban;
  const int numLandunits = 4;
  const int numSoilLayers = 15;
  int size1D;
  int size2D_road[2];

  void SetUp() override {
    UrbanErrorCode ierr;
    UrbanCreate(numLandunits, &urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    size1D          = numLandunits;
    size2D_road[0]  = numLandunits;
    size2D_road[1]  = numSoilLayers;
    SetFullModelParameters(urban, numLandunits);
    UrbanSetup(urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanComputeNetLongwave(urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanComputeSurfaceFluxes(urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanComputeHeatDiffusion(urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    UrbanComputeHydrology(urban, 1800.0, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
  }

  void TearDown() override {
    UrbanErrorCode ierr;
    UrbanDestroy(&urban, &ierr);
  }
};

// ---------------------------------------------------------------------------
// UrbanGetSoilLiquidWaterPerviousRoad  (2-D)
// ---------------------------------------------------------------------------
TEST_F(HydrologyGetterTest, GetSoilLiquidWaterPerviousRoad_ValidData) {
  std::vector<double> buf(numLandunits * numSoilLayers);
  UrbanErrorCode ierr;
  UrbanGetSoilLiquidWaterPerviousRoad(urban, buf.data(), size2D_road, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  for (double v : buf) {
    EXPECT_FALSE(std::isnan(v)) << "NaN in SoilLiquidWaterPerviousRoad";
    EXPECT_GE(v, 0.0) << "Negative liquid water content";
  }
}

TEST_F(HydrologyGetterTest, GetSoilLiquidWaterPerviousRoad_SizeMismatch) {
  std::vector<double> buf(numLandunits * numSoilLayers);
  int badSize[2] = {numLandunits, numSoilLayers + 1};
  UrbanErrorCode ierr;
  UrbanGetSoilLiquidWaterPerviousRoad(urban, buf.data(), badSize, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(HydrologyGetterTest, GetSoilLiquidWaterPerviousRoad_NullPointer) {
  UrbanErrorCode ierr;
  UrbanGetSoilLiquidWaterPerviousRoad(urban, nullptr, size2D_road, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT);
}

// ---------------------------------------------------------------------------
// UrbanGetSoilVolumetricWaterPerviousRoad  (2-D)
// ---------------------------------------------------------------------------
TEST_F(HydrologyGetterTest, GetSoilVolumetricWaterPerviousRoad_ValidData) {
  std::vector<double> buf(numLandunits * numSoilLayers);
  UrbanErrorCode ierr;
  UrbanGetSoilVolumetricWaterPerviousRoad(urban, buf.data(), size2D_road, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  for (double v : buf) {
    EXPECT_FALSE(std::isnan(v)) << "NaN in SoilVolumetricWaterPerviousRoad";
    EXPECT_GE(v, 0.0) << "Negative volumetric water content";
    EXPECT_LE(v, 1.0) << "Volumetric water content > 1";
  }
}

TEST_F(HydrologyGetterTest, GetSoilVolumetricWaterPerviousRoad_SizeMismatch) {
  std::vector<double> buf(numLandunits * numSoilLayers);
  int badSize[2] = {numLandunits + 1, numSoilLayers};
  UrbanErrorCode ierr;
  UrbanGetSoilVolumetricWaterPerviousRoad(urban, buf.data(), badSize, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(HydrologyGetterTest, GetSoilVolumetricWaterPerviousRoad_NullPointer) {
  UrbanErrorCode ierr;
  UrbanGetSoilVolumetricWaterPerviousRoad(urban, nullptr, size2D_road, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT);
}

// ---------------------------------------------------------------------------
// UrbanGetAquiferRechargeRatePerviousRoad  (1-D)
// ---------------------------------------------------------------------------
TEST_F(HydrologyGetterTest, GetAquiferRechargeRatePerviousRoad_ValidData) {
  std::vector<double> buf(numLandunits);
  UrbanErrorCode ierr;
  UrbanGetAquiferRechargeRatePerviousRoad(urban, buf.data(), size1D, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  for (double v : buf) {
    EXPECT_FALSE(std::isnan(v)) << "NaN in AquiferRechargeRatePerviousRoad";
  }
}

TEST_F(HydrologyGetterTest, GetAquiferRechargeRatePerviousRoad_SizeMismatch) {
  std::vector<double> buf(numLandunits);
  UrbanErrorCode ierr;
  UrbanGetAquiferRechargeRatePerviousRoad(urban, buf.data(), size1D + 1, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(HydrologyGetterTest, GetAquiferRechargeRatePerviousRoad_NullPointer) {
  UrbanErrorCode ierr;
  UrbanGetAquiferRechargeRatePerviousRoad(urban, nullptr, size1D, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT);
}

// ---------------------------------------------------------------------------
// UrbanGetWaterDeficitFluxPerviousRoad  (1-D)
// ---------------------------------------------------------------------------
TEST_F(HydrologyGetterTest, GetWaterDeficitFluxPerviousRoad_ValidData) {
  std::vector<double> buf(numLandunits);
  UrbanErrorCode ierr;
  UrbanGetWaterDeficitFluxPerviousRoad(urban, buf.data(), size1D, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  for (double v : buf) {
    EXPECT_FALSE(std::isnan(v)) << "NaN in WaterDeficitFluxPerviousRoad";
  }
}

TEST_F(HydrologyGetterTest, GetWaterDeficitFluxPerviousRoad_SizeMismatch) {
  std::vector<double> buf(numLandunits);
  UrbanErrorCode ierr;
  UrbanGetWaterDeficitFluxPerviousRoad(urban, buf.data(), size1D + 1, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(HydrologyGetterTest, GetWaterDeficitFluxPerviousRoad_NullPointer) {
  UrbanErrorCode ierr;
  UrbanGetWaterDeficitFluxPerviousRoad(urban, nullptr, size1D, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT);
}

// ---------------------------------------------------------------------------
// Integration tests
// ---------------------------------------------------------------------------
TEST_F(HydrologyGetterTest, AllHydrologyGetters_Sequential) {
  std::vector<double> liq(numLandunits * numSoilLayers);
  std::vector<double> vol(numLandunits * numSoilLayers);
  std::vector<double> rech(numLandunits);
  std::vector<double> deficit(numLandunits);
  UrbanErrorCode ierr;

  UrbanGetSoilLiquidWaterPerviousRoad(urban, liq.data(), size2D_road, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  UrbanGetSoilVolumetricWaterPerviousRoad(urban, vol.data(), size2D_road, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  UrbanGetAquiferRechargeRatePerviousRoad(urban, rech.data(), size1D, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  UrbanGetWaterDeficitFluxPerviousRoad(urban, deficit.data(), size1D, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);

  for (double v : liq)     EXPECT_FALSE(std::isnan(v));
  for (double v : vol)     EXPECT_FALSE(std::isnan(v));
  for (double v : rech)    EXPECT_FALSE(std::isnan(v));
  for (double v : deficit) EXPECT_FALSE(std::isnan(v));
}

TEST_F(HydrologyGetterTest, MultipleCallsSameGetter_ConsistentResults) {
  std::vector<double> first(numLandunits * numSoilLayers);
  std::vector<double> second(numLandunits * numSoilLayers);
  UrbanErrorCode ierr;

  UrbanGetSoilLiquidWaterPerviousRoad(urban, first.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanGetSoilLiquidWaterPerviousRoad(urban, second.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  for (int i = 0; i < numLandunits * numSoilLayers; ++i) {
    EXPECT_DOUBLE_EQ(first[i], second[i])
        << "Inconsistent soil liquid water at index " << i;
  }
}

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
