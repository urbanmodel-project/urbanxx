#include "Urban.h"
#include "test_full_model_setup.h"
#include <Kokkos_Core.hpp>
#include <cmath>
#include <vector>

class LayerTempGetterTest : public ::testing::Test {
protected:
  UrbanType urban;
  const int numLandunits = 4;
  const int numUrbanLayers = 5;
  const int numSoilLayers = 15;
  int size2D_urban[2];
  int size2D_road[2];

  void SetUp() override {
    UrbanErrorCode ierr;
    UrbanCreate(numLandunits, &urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS);
    size2D_urban[0] = numLandunits;
    size2D_urban[1] = numUrbanLayers;
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
  }

  void TearDown() override {
    UrbanErrorCode ierr;
    UrbanDestroy(&urban, &ierr);
  }
};

// ---------------------------------------------------------------------------
// UrbanGetLayerTempRoof
// ---------------------------------------------------------------------------
TEST_F(LayerTempGetterTest, GetLayerTempRoof_ValidData) {
  std::vector<double> buf(numLandunits * numUrbanLayers);
  UrbanErrorCode ierr;
  UrbanGetLayerTempRoof(urban, buf.data(), size2D_urban, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  for (double v : buf) {
    EXPECT_FALSE(std::isnan(v)) << "NaN in LayerTempRoof";
    EXPECT_GT(v, 0.0) << "Non-positive temperature in LayerTempRoof";
  }
}

TEST_F(LayerTempGetterTest, GetLayerTempRoof_SizeMismatch) {
  std::vector<double> buf(numLandunits * numUrbanLayers);
  int badSize[2] = {numLandunits + 1, numUrbanLayers};
  UrbanErrorCode ierr;
  UrbanGetLayerTempRoof(urban, buf.data(), badSize, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(LayerTempGetterTest, GetLayerTempRoof_NullPointer) {
  UrbanErrorCode ierr;
  UrbanGetLayerTempRoof(urban, nullptr, size2D_urban, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT);
}

// ---------------------------------------------------------------------------
// UrbanGetLayerTempImperviousRoad
// ---------------------------------------------------------------------------
TEST_F(LayerTempGetterTest, GetLayerTempImperviousRoad_ValidData) {
  std::vector<double> buf(numLandunits * numSoilLayers);
  UrbanErrorCode ierr;
  UrbanGetLayerTempImperviousRoad(urban, buf.data(), size2D_road, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  for (double v : buf) {
    EXPECT_FALSE(std::isnan(v)) << "NaN in LayerTempImperviousRoad";
    EXPECT_GT(v, 0.0) << "Non-positive temperature in LayerTempImperviousRoad";
  }
}

TEST_F(LayerTempGetterTest, GetLayerTempImperviousRoad_SizeMismatch) {
  std::vector<double> buf(numLandunits * numSoilLayers);
  int badSize[2] = {numLandunits, numSoilLayers + 1};
  UrbanErrorCode ierr;
  UrbanGetLayerTempImperviousRoad(urban, buf.data(), badSize, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(LayerTempGetterTest, GetLayerTempImperviousRoad_NullPointer) {
  UrbanErrorCode ierr;
  UrbanGetLayerTempImperviousRoad(urban, nullptr, size2D_road, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT);
}

// ---------------------------------------------------------------------------
// UrbanGetLayerTempPerviousRoad
// ---------------------------------------------------------------------------
TEST_F(LayerTempGetterTest, GetLayerTempPerviousRoad_ValidData) {
  std::vector<double> buf(numLandunits * numSoilLayers);
  UrbanErrorCode ierr;
  UrbanGetLayerTempPerviousRoad(urban, buf.data(), size2D_road, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  for (double v : buf) {
    EXPECT_FALSE(std::isnan(v)) << "NaN in LayerTempPerviousRoad";
    EXPECT_GT(v, 0.0) << "Non-positive temperature in LayerTempPerviousRoad";
  }
}

TEST_F(LayerTempGetterTest, GetLayerTempPerviousRoad_SizeMismatch) {
  std::vector<double> buf(numLandunits * numSoilLayers);
  int badSize[2] = {numLandunits + 1, numSoilLayers};
  UrbanErrorCode ierr;
  UrbanGetLayerTempPerviousRoad(urban, buf.data(), badSize, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(LayerTempGetterTest, GetLayerTempPerviousRoad_NullPointer) {
  UrbanErrorCode ierr;
  UrbanGetLayerTempPerviousRoad(urban, nullptr, size2D_road, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT);
}

// ---------------------------------------------------------------------------
// UrbanGetLayerTempSunlitWall
// ---------------------------------------------------------------------------
TEST_F(LayerTempGetterTest, GetLayerTempSunlitWall_ValidData) {
  std::vector<double> buf(numLandunits * numUrbanLayers);
  UrbanErrorCode ierr;
  UrbanGetLayerTempSunlitWall(urban, buf.data(), size2D_urban, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  for (double v : buf) {
    EXPECT_FALSE(std::isnan(v)) << "NaN in LayerTempSunlitWall";
    EXPECT_GT(v, 0.0) << "Non-positive temperature in LayerTempSunlitWall";
  }
}

TEST_F(LayerTempGetterTest, GetLayerTempSunlitWall_SizeMismatch) {
  std::vector<double> buf(numLandunits * numUrbanLayers);
  int badSize[2] = {numLandunits, numUrbanLayers + 1};
  UrbanErrorCode ierr;
  UrbanGetLayerTempSunlitWall(urban, buf.data(), badSize, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(LayerTempGetterTest, GetLayerTempSunlitWall_NullPointer) {
  UrbanErrorCode ierr;
  UrbanGetLayerTempSunlitWall(urban, nullptr, size2D_urban, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT);
}

// ---------------------------------------------------------------------------
// UrbanGetLayerTempShadedWall
// ---------------------------------------------------------------------------
TEST_F(LayerTempGetterTest, GetLayerTempShadedWall_ValidData) {
  std::vector<double> buf(numLandunits * numUrbanLayers);
  UrbanErrorCode ierr;
  UrbanGetLayerTempShadedWall(urban, buf.data(), size2D_urban, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  for (double v : buf) {
    EXPECT_FALSE(std::isnan(v)) << "NaN in LayerTempShadedWall";
    EXPECT_GT(v, 0.0) << "Non-positive temperature in LayerTempShadedWall";
  }
}

TEST_F(LayerTempGetterTest, GetLayerTempShadedWall_SizeMismatch) {
  std::vector<double> buf(numLandunits * numUrbanLayers);
  int badSize[2] = {numLandunits + 1, numUrbanLayers};
  UrbanErrorCode ierr;
  UrbanGetLayerTempShadedWall(urban, buf.data(), badSize, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(LayerTempGetterTest, GetLayerTempShadedWall_NullPointer) {
  UrbanErrorCode ierr;
  UrbanGetLayerTempShadedWall(urban, nullptr, size2D_urban, &ierr);
  EXPECT_EQ(ierr, URBAN_ERR_INVALID_ARGUMENT);
}

// ---------------------------------------------------------------------------
// Integration tests
// ---------------------------------------------------------------------------
TEST_F(LayerTempGetterTest, AllLayerTempGetters_Sequential) {
  std::vector<double> roof(numLandunits * numUrbanLayers);
  std::vector<double> imp(numLandunits * numSoilLayers);
  std::vector<double> perv(numLandunits * numSoilLayers);
  std::vector<double> sunlit(numLandunits * numUrbanLayers);
  std::vector<double> shaded(numLandunits * numUrbanLayers);
  UrbanErrorCode ierr;

  UrbanGetLayerTempRoof(urban, roof.data(), size2D_urban, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  UrbanGetLayerTempImperviousRoad(urban, imp.data(), size2D_road, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  UrbanGetLayerTempPerviousRoad(urban, perv.data(), size2D_road, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  UrbanGetLayerTempSunlitWall(urban, sunlit.data(), size2D_urban, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);
  UrbanGetLayerTempShadedWall(urban, shaded.data(), size2D_urban, &ierr);
  EXPECT_EQ(ierr, URBAN_SUCCESS);

  for (double v : roof)   EXPECT_FALSE(std::isnan(v));
  for (double v : imp)    EXPECT_FALSE(std::isnan(v));
  for (double v : perv)   EXPECT_FALSE(std::isnan(v));
  for (double v : sunlit) EXPECT_FALSE(std::isnan(v));
  for (double v : shaded) EXPECT_FALSE(std::isnan(v));
}

TEST_F(LayerTempGetterTest, MultipleCallsSameGetter_ConsistentResults) {
  std::vector<double> first(numLandunits * numUrbanLayers);
  std::vector<double> second(numLandunits * numUrbanLayers);
  UrbanErrorCode ierr;

  UrbanGetLayerTempRoof(urban, first.data(), size2D_urban, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanGetLayerTempRoof(urban, second.data(), size2D_urban, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  for (int i = 0; i < numLandunits * numUrbanLayers; ++i) {
    EXPECT_DOUBLE_EQ(first[i], second[i])
        << "Inconsistent roof layer temp at index " << i;
  }
}

TEST_F(LayerTempGetterTest, WallLayerTemps_SunlitAndShaded_SameInit) {
  // Both walls initialised with the same surface temperature so their
  // layer profiles should be equal after heat diffusion.
  std::vector<double> sunlit(numLandunits * numUrbanLayers);
  std::vector<double> shaded(numLandunits * numUrbanLayers);
  UrbanErrorCode ierr;

  UrbanGetLayerTempSunlitWall(urban, sunlit.data(), size2D_urban, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanGetLayerTempShadedWall(urban, shaded.data(), size2D_urban, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  for (int i = 0; i < numLandunits * numUrbanLayers; ++i) {
    EXPECT_DOUBLE_EQ(sunlit[i], shaded[i])
        << "Sunlit and shaded wall layer temps differ at index " << i;
  }
}

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
