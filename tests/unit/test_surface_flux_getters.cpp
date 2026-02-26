#include "Urban.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <cmath>
#include <vector>

// Helper to set up the minimum parameters required before running
// UrbanComputeSurfaceFluxes.  Mirrors the setup in test_canyon_air_getters.cpp.
static void SetSurfaceFluxParameters(UrbanType urban, int numLandunits) {
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

  double emissivity[5] = {0.90, 0.91, 0.92, 0.93, 0.94};
  UrbanSetEmissivityRoof(urban, emissivity, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetEmissivityWall(urban, emissivity, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetEmissivityImperviousRoad(urban, emissivity, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetEmissivityPerviousRoad(urban, emissivity, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  double heights[5] = {10.0, 11.0, 12.0, 13.0, 14.0};
  UrbanSetForcHgtT(urban, heights, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetForcHgtU(urban, heights, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  double z0[5] = {0.1, 0.15, 0.2, 0.25, 0.3};
  UrbanSetZ0Town(urban, z0, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  double zd[5] = {1.0, 1.5, 2.0, 2.5, 3.0};
  UrbanSetZDTown(urban, zd, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  double ht_roof[5] = {8.0, 9.0, 10.0, 11.0, 12.0};
  UrbanSetHtRoof(urban, ht_roof, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  double wind_hgt[5] = {4.0, 4.5, 5.0, 5.5, 6.0};
  UrbanSetWindHgtCanyon(urban, wind_hgt, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // Atmospheric forcing
  double atm_temp[5] = {300.0, 301.0, 302.0, 303.0, 304.0};
  double atm_pot_temp[5] = {301.0, 302.0, 303.0, 304.0, 305.0};
  double atm_press[5] = {101325.0, 101325.0, 101325.0, 101325.0, 101325.0};
  double atm_rho[5] = {1.2, 1.2, 1.2, 1.2, 1.2};
  double atm_spc_humd[5] = {0.01, 0.01, 0.01, 0.01, 0.01};
  double atm_wind_u[5] = {3.0, 3.5, 4.0, 4.5, 5.0};
  double atm_wind_v[5] = {1.0, 1.5, 2.0, 2.5, 3.0};

  UrbanSetAtmTemp(urban, atm_temp, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmPotTemp(urban, atm_pot_temp, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmPress(urban, atm_press, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmRho(urban, atm_rho, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmSpcHumd(urban, atm_spc_humd, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmWindU(urban, atm_wind_u, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmWindV(urban, atm_wind_v, numLandunits, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
}

// Test fixture for surface flux getter tests.
// Covers sensible heat flux, evaporative flux, cgrnds, and cgrndl getters.
class SurfaceFluxGetterTest : public ::testing::Test {
protected:
  UrbanType urban;
  const int numLandunits = 5;

  void SetUp() override {
    UrbanErrorCode ierr;
    UrbanCreate(numLandunits, &urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS) << "Failed to create Urban object in SetUp";

    SetSurfaceFluxParameters(urban, numLandunits);

    UrbanComputeSurfaceFluxes(urban, &ierr);
    ASSERT_EQ(ierr, URBAN_SUCCESS) << "Failed to compute surface fluxes";
  }

  void TearDown() override {
    UrbanErrorCode ierr;
    UrbanDestroy(&urban, &ierr);
  }
};

// =============================================================================
// Sensible Heat Flux Getter Tests
// =============================================================================

TEST_F(SurfaceFluxGetterTest, GetSensibleHeatFluxRoof_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetSensibleHeatFluxRoof(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i])) << "NaN at index " << i;
}

TEST_F(SurfaceFluxGetterTest, GetSensibleHeatFluxRoof_SizeMismatch) {
  std::vector<double> values(numLandunits * 2);
  UrbanErrorCode status;

  UrbanGetSensibleHeatFluxRoof(urban, values.data(), numLandunits * 2, &status);
  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(SurfaceFluxGetterTest, GetSensibleHeatFluxRoof_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetSensibleHeatFluxRoof(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetSensibleHeatFluxRoof(nullptr, values.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  // Null status should not crash
  UrbanGetSensibleHeatFluxRoof(urban, values.data(), numLandunits, nullptr);
}

TEST_F(SurfaceFluxGetterTest, GetSensibleHeatFluxImperviousRoad_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetSensibleHeatFluxImperviousRoad(urban, values.data(), numLandunits,
                                         &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetSensibleHeatFluxImperviousRoad_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetSensibleHeatFluxImperviousRoad(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetSensibleHeatFluxImperviousRoad(nullptr, values.data(), numLandunits,
                                          &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetSensibleHeatFluxImperviousRoad(urban, values.data(), numLandunits,
                                          nullptr);
}

TEST_F(SurfaceFluxGetterTest, GetSensibleHeatFluxPerviousRoad_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetSensibleHeatFluxPerviousRoad(urban, values.data(), numLandunits,
                                       &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetSensibleHeatFluxPerviousRoad_SizeMismatch) {
  std::vector<double> values(numLandunits * 2);
  UrbanErrorCode status;

  UrbanGetSensibleHeatFluxPerviousRoad(urban, values.data(), numLandunits * 2,
                                        &status);
  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(SurfaceFluxGetterTest, GetSensibleHeatFluxSunlitWall_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetSensibleHeatFluxSunlitWall(urban, values.data(), numLandunits,
                                     &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetSensibleHeatFluxShadedWall_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetSensibleHeatFluxShadedWall(urban, values.data(), numLandunits,
                                     &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetSensibleHeatFluxShadedWall_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetSensibleHeatFluxShadedWall(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetSensibleHeatFluxShadedWall(nullptr, values.data(), numLandunits,
                                     &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetSensibleHeatFluxShadedWall(urban, values.data(), numLandunits,
                                     nullptr);
}

// =============================================================================
// Evaporation Flux Getter Tests
// =============================================================================

TEST_F(SurfaceFluxGetterTest, GetEvapFluxRoof_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetEvapFluxRoof(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i])) << "NaN at index " << i;
}

TEST_F(SurfaceFluxGetterTest, GetEvapFluxRoof_SizeMismatch) {
  std::vector<double> values(numLandunits * 2);
  UrbanErrorCode status;

  UrbanGetEvapFluxRoof(urban, values.data(), numLandunits * 2, &status);
  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(SurfaceFluxGetterTest, GetEvapFluxRoof_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetEvapFluxRoof(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetEvapFluxRoof(nullptr, values.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  // Null status should not crash
  UrbanGetEvapFluxRoof(urban, values.data(), numLandunits, nullptr);
}

TEST_F(SurfaceFluxGetterTest, GetEvapFluxImperviousRoad_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetEvapFluxImperviousRoad(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetEvapFluxImperviousRoad_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetEvapFluxImperviousRoad(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetEvapFluxImperviousRoad(nullptr, values.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetEvapFluxImperviousRoad(urban, values.data(), numLandunits, nullptr);
}

TEST_F(SurfaceFluxGetterTest, GetEvapFluxPerviousRoad_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetEvapFluxPerviousRoad(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetEvapFluxPerviousRoad_SizeMismatch) {
  std::vector<double> values(numLandunits * 2);
  UrbanErrorCode status;

  UrbanGetEvapFluxPerviousRoad(urban, values.data(), numLandunits * 2, &status);
  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

// =============================================================================
// Cgrnds (d(SH)/dT) Getter Tests
// =============================================================================

TEST_F(SurfaceFluxGetterTest, GetCgrndsRoof_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetCgrndsRoof(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetCgrndsRoof_SizeMismatch) {
  std::vector<double> values(numLandunits * 2);
  UrbanErrorCode status;

  UrbanGetCgrndsRoof(urban, values.data(), numLandunits * 2, &status);
  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(SurfaceFluxGetterTest, GetCgrndsRoof_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetCgrndsRoof(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetCgrndsRoof(nullptr, values.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetCgrndsRoof(urban, values.data(), numLandunits, nullptr);
}

TEST_F(SurfaceFluxGetterTest, GetCgrndsImperviousRoad_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetCgrndsImperviousRoad(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetCgrndsPerviousRoad_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetCgrndsPerviousRoad(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetCgrndsSunlitWall_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetCgrndsSunlitWall(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetCgrndsSunlitWall_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetCgrndsSunlitWall(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetCgrndsSunlitWall(nullptr, values.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetCgrndsSunlitWall(urban, values.data(), numLandunits, nullptr);
}

TEST_F(SurfaceFluxGetterTest, GetCgrndsShadedWall_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetCgrndsShadedWall(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

// =============================================================================
// Cgrndl (d(LH)/dT) Getter Tests
// =============================================================================

TEST_F(SurfaceFluxGetterTest, GetCgrndlRoof_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetCgrndlRoof(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetCgrndlRoof_SizeMismatch) {
  std::vector<double> values(numLandunits * 2);
  UrbanErrorCode status;

  UrbanGetCgrndlRoof(urban, values.data(), numLandunits * 2, &status);
  EXPECT_EQ(status, URBAN_ERR_SIZE_MISMATCH);
}

TEST_F(SurfaceFluxGetterTest, GetCgrndlRoof_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetCgrndlRoof(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetCgrndlRoof(nullptr, values.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetCgrndlRoof(urban, values.data(), numLandunits, nullptr);
}

TEST_F(SurfaceFluxGetterTest, GetCgrndlImperviousRoad_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetCgrndlImperviousRoad(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetCgrndlPerviousRoad_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetCgrndlPerviousRoad(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetCgrndlSunlitWall_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetCgrndlSunlitWall(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetCgrndlSunlitWall_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetCgrndlSunlitWall(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetCgrndlSunlitWall(nullptr, values.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetCgrndlSunlitWall(urban, values.data(), numLandunits, nullptr);
}

TEST_F(SurfaceFluxGetterTest, GetCgrndlShadedWall_ValidData) {
  std::vector<double> values(numLandunits);
  UrbanErrorCode status;

  UrbanGetCgrndlShadedWall(urban, values.data(), numLandunits, &status);

  EXPECT_EQ(status, URBAN_SUCCESS);
  for (int i = 0; i < numLandunits; ++i)
    EXPECT_FALSE(std::isnan(values[i]));
}

TEST_F(SurfaceFluxGetterTest, GetCgrndlShadedWall_NullPointer) {
  UrbanErrorCode status;
  std::vector<double> values(numLandunits);

  UrbanGetCgrndlShadedWall(urban, nullptr, numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetCgrndlShadedWall(nullptr, values.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_ERR_INVALID_ARGUMENT);

  UrbanGetCgrndlShadedWall(urban, values.data(), numLandunits, nullptr);
}

// =============================================================================
// Integration Tests
// =============================================================================

TEST_F(SurfaceFluxGetterTest, AllSensibleHeatFluxGetters_Sequential) {
  UrbanErrorCode status;
  std::vector<double> roof(numLandunits), imp_road(numLandunits),
      perv_road(numLandunits), sun_wall(numLandunits), sha_wall(numLandunits);

  UrbanGetSensibleHeatFluxRoof(urban, roof.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanGetSensibleHeatFluxImperviousRoad(urban, imp_road.data(), numLandunits,
                                         &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanGetSensibleHeatFluxPerviousRoad(urban, perv_road.data(), numLandunits,
                                       &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanGetSensibleHeatFluxSunlitWall(urban, sun_wall.data(), numLandunits,
                                     &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanGetSensibleHeatFluxShadedWall(urban, sha_wall.data(), numLandunits,
                                     &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
}

TEST_F(SurfaceFluxGetterTest, AllEvapFluxGetters_Sequential) {
  UrbanErrorCode status;
  std::vector<double> roof(numLandunits), imp_road(numLandunits),
      perv_road(numLandunits);

  UrbanGetEvapFluxRoof(urban, roof.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanGetEvapFluxImperviousRoad(urban, imp_road.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanGetEvapFluxPerviousRoad(urban, perv_road.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
}

TEST_F(SurfaceFluxGetterTest, AllCgrndsGetters_Sequential) {
  UrbanErrorCode status;
  std::vector<double> roof(numLandunits), imp_road(numLandunits),
      perv_road(numLandunits), sun_wall(numLandunits), sha_wall(numLandunits);

  UrbanGetCgrndsRoof(urban, roof.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  UrbanGetCgrndsImperviousRoad(urban, imp_road.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  UrbanGetCgrndsPerviousRoad(urban, perv_road.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  UrbanGetCgrndsSunlitWall(urban, sun_wall.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  UrbanGetCgrndsShadedWall(urban, sha_wall.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
}

TEST_F(SurfaceFluxGetterTest, AllCgrndlGetters_Sequential) {
  UrbanErrorCode status;
  std::vector<double> roof(numLandunits), imp_road(numLandunits),
      perv_road(numLandunits), sun_wall(numLandunits), sha_wall(numLandunits);

  UrbanGetCgrndlRoof(urban, roof.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  UrbanGetCgrndlImperviousRoad(urban, imp_road.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  UrbanGetCgrndlPerviousRoad(urban, perv_road.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  UrbanGetCgrndlSunlitWall(urban, sun_wall.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
  UrbanGetCgrndlShadedWall(urban, sha_wall.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);
}

TEST_F(SurfaceFluxGetterTest, MultipleCallsSameGetter_ConsistentResults) {
  UrbanErrorCode status;
  std::vector<double> vals1(numLandunits), vals2(numLandunits);

  UrbanGetSensibleHeatFluxRoof(urban, vals1.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  UrbanGetSensibleHeatFluxRoof(urban, vals2.data(), numLandunits, &status);
  EXPECT_EQ(status, URBAN_SUCCESS);

  for (int i = 0; i < numLandunits; ++i)
    EXPECT_DOUBLE_EQ(vals1[i], vals2[i])
        << "Repeated getter calls should return identical values";
}

// Main function to initialize Kokkos
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();

  return result;
}
