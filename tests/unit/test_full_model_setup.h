// Shared setup helper used by layer temperature and hydrology getter tests.
// Provides a single function that sets all parameters required for a full
// UrbanSetup + UrbanAdvance workflow.
//
// NOTE: This header uses ASSERT_EQ which only works inside a Google Test
// context (test body or SetUp/TearDown).  Include it only from gtest .cpp
// files.
#pragma once

#include "Urban.h"
#include <gtest/gtest.h>
#include <vector>

// Sets the full set of parameters needed before calling UrbanSetup and the
// downstream compute functions.  N is the number of landunits.
static void SetFullModelParameters(UrbanType urban, int N) {
  UrbanErrorCode ierr;
  const int numUrbanLayers = 5;
  const int numSoilLayers = 15;
  const int numBands = 2;
  const int numRadTypes = 2;

  // --- Geometric parameters ---
  std::vector<double> canyon_hwr(N, 2.0);
  UrbanSetCanyonHwr(urban, canyon_hwr.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  std::vector<double> frac_perv(N, 0.2);
  UrbanSetFracPervRoadOfTotalRoad(urban, frac_perv.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  std::vector<double> wt_roof(N, 0.4);
  UrbanSetWtRoof(urban, wt_roof.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Height parameters ---
  // Values must satisfy: forc_hgt_u > ht_roof > z_d_town > z_0_town
  // and wind_hgt_canyon < ht_roof, to avoid log(0) in ComputeCanyonUWind.
  std::vector<double> forc_hgt(N, 30.0); // forcing height above surface (m)
  UrbanSetForcHgtT(urban, forc_hgt.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetForcHgtU(urban, forc_hgt.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  std::vector<double> zd(N, 7.0); // zero-plane displacement (~0.7 * ht_roof)
  UrbanSetZDTown(urban, zd.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  std::vector<double> z0(N, 0.5); // roughness length (m)
  UrbanSetZ0Town(urban, z0.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  std::vector<double> ht_roof(N, 10.0);
  UrbanSetHtRoof(urban, ht_roof.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  std::vector<double> wind_hgt(N, 5.0); // mid-canyon height (< ht_roof)
  UrbanSetWindHgtCanyon(urban, wind_hgt.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Emissivities ---
  std::vector<double> emiss(N, 0.90);
  UrbanSetEmissivityRoof(urban, emiss.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetEmissivityWall(urban, emiss.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetEmissivityImperviousRoad(urban, emiss.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetEmissivityPerviousRoad(urban, emiss.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Albedos (3D: landunits x bands x rad_types) ---
  const int size3D[3] = {N, numBands, numRadTypes};
  const int total3D = N * numBands * numRadTypes;
  std::vector<double> albedo(total3D, 0.20);
  UrbanSetAlbedoRoof(urban, albedo.data(), size3D, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAlbedoSunlitWall(urban, albedo.data(), size3D, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAlbedoShadedWall(urban, albedo.data(), size3D, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAlbedoImperviousRoad(urban, albedo.data(), size3D, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAlbedoPerviousRoad(urban, albedo.data(), size3D, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Thermal conductivity (2D) ---
  const int size2D_urban[2] = {N, numUrbanLayers};
  const int size2D_road[2] = {N, numSoilLayers};
  std::vector<double> tk_urban(N * numUrbanLayers, 1.0);
  std::vector<double> tk_road(N * numSoilLayers, 1.0);
  UrbanSetThermalConductivityWall(urban, tk_urban.data(), size2D_urban, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetThermalConductivityRoof(urban, tk_urban.data(), size2D_urban, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetThermalConductivityRoad(urban, tk_road.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Heat capacity (2D) ---
  std::vector<double> cv_urban(N * numUrbanLayers, 1.5e6);
  std::vector<double> cv_road(N * numSoilLayers, 1.5e6);
  UrbanSetHeatCapacityWall(urban, cv_urban.data(), size2D_urban, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetHeatCapacityRoof(urban, cv_urban.data(), size2D_urban, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetHeatCapacityRoad(urban, cv_road.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Building parameters ---
  std::vector<double> bld_min(N, 285.0), bld_max(N, 310.0);
  UrbanSetBuildingMinTemperature(urban, bld_min.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetBuildingMaxTemperature(urban, bld_max.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  // Internal building temperature (bottom BC for roof and wall heat diffusion)
  std::vector<double> bld_temp(N, 295.0);
  UrbanSetBuildingTemperature(urban, bld_temp.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  std::vector<double> wall_thick(N, 0.286), roof_thick(N, 0.217);
  UrbanSetBuildingWallThickness(urban, wall_thick.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetBuildingRoofThickness(urban, roof_thick.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Soil properties (2D, pervious road) ---
  std::vector<double> sand(N * numSoilLayers, 40.0);
  std::vector<double> clay(N * numSoilLayers, 40.0);
  std::vector<double> organic(N * numSoilLayers, 0.0);
  UrbanSetSandPerviousRoad(urban, sand.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetClayPerviousRoad(urban, clay.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetOrganicPerviousRoad(urban, organic.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Soil water content (2D, pervious road) ---
  // Needed by ComputeSoilThermalConductivity and ComputeSoilHeatCapacityTimesDz
  // inside UrbanComputeHeatDiffusion.
  std::vector<double> soilLiq(N * numSoilLayers, 0.1);
  UrbanSetSoilLiquidWaterForPerviousRoad(urban, soilLiq.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  std::vector<double> soilIce(N * numSoilLayers, 0.0);
  UrbanSetSoilIceContentForPerviousRoad(urban, soilIce.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  std::vector<double> soilVol(N * numSoilLayers, 0.3);
  UrbanSetSoilVolumetricWaterForPerviousRoad(urban, soilVol.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Number of active layers ---
  std::vector<double> nActLayers(N, 2.0);
  UrbanSetNumberOfActiveLayersImperviousRoad(urban, nActLayers.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Layer temperatures (2D) ---
  std::vector<double> ltemp_urban(N * numUrbanLayers, 292.0);
  std::vector<double> ltemp_road(N * numSoilLayers, 274.0);
  UrbanSetLayerTempRoof(urban, ltemp_urban.data(), size2D_urban, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetLayerTempSunlitWall(urban, ltemp_urban.data(), size2D_urban, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetLayerTempShadedWall(urban, ltemp_urban.data(), size2D_urban, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetLayerTempImperviousRoad(urban, ltemp_road.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetLayerTempPerviousRoad(urban, ltemp_road.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Canyon air state ---
  std::vector<double> taf(N, 283.0), qaf(N, 1.e-4);
  UrbanSetCanyonAirTemperature(urban, taf.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetCanyonSpecificHumidity(urban, qaf.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Atmospheric forcing ---
  std::vector<double> atm_temp(N, 297.0), atm_press(N, 98260.0),
      atm_rho(N, 1.14), atm_spc_humd(N, 0.019), atm_wind_u(N, 0.5),
      atm_wind_v(N, 0.0), atm_coszen(N, 0.01), atm_snow(N, 0.0),
      atm_lw(N, 430.0);
  UrbanSetAtmTemp(urban, atm_temp.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmPotTemp(urban, atm_temp.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmPress(urban, atm_press.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmRho(urban, atm_rho.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmSpcHumd(urban, atm_spc_humd.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmWindU(urban, atm_wind_u.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmWindV(urban, atm_wind_v.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmCoszen(urban, atm_coszen.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmFracSnow(urban, atm_snow.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetAtmLongwaveDown(urban, atm_lw.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  std::vector<double> sw(N * numBands * numRadTypes, 0.0);
  UrbanSetAtmShortwaveDown(urban, sw.data(), size3D, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  // --- Hydrology boundary conditions ---
  std::vector<double> zwt(N, 4.8);
  UrbanSetWaterTableDepth(urban, zwt.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  std::vector<double> qinfl(N, 0.0);
  UrbanSetInfiltrationFluxForPerviousRoad(urban, qinfl.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  std::vector<double> qtran(N * numSoilLayers, 0.0);
  UrbanSetTranspirationFluxForPerviousRoad(urban, qtran.data(), size2D_road, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);

  std::vector<double> fwet(N, 0.0);
  UrbanSetFractionWetImperviousRoad(urban, fwet.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
  UrbanSetFractionWetRoof(urban, fwet.data(), N, &ierr);
  ASSERT_EQ(ierr, URBAN_SUCCESS);
}
