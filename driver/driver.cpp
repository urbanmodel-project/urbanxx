#include <Kokkos_Core.hpp>
#include <Urban.h>
#include <UrbanMacros.h>
#include <iostream>
#include <mpi.h>
#include <stdlib.h>

// Constants
const int NUM_LEVELS = 5;
const int NUM_URBAN_DENSITY_CLASSES = 3;

// Helper function to allocate memory with error checking
static double *AllocateArray(int size, const char *name) {
  double *arr = (double *)malloc(size * sizeof(double));
  if (arr == NULL) {
    fprintf(stderr, "Failed to allocate memory for %s\n", name);
    exit(1);
  }
  return arr;
}

// Helper function to safely free multiple arrays
static void FreeArrays(double **arrays, int count) {
  for (int i = 0; i < count; ++i) {
    free(arrays[i]);
    arrays[i] = NULL;
  }
}

void SetCanyonHwr(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  double *canyonHwr = AllocateArray(numLandunits, "canyonHwr");
  for (int i = 0; i < numLandunits; ++i) {
    canyonHwr[i] = 4.80000019073486;
  }
  UrbanCall(UrbanSetCanyonHwr(urban, canyonHwr, numLandunits, &ierr), &ierr);
  free(canyonHwr);

  if (mpi_rank == 0) {
    std::cout << "Set canyon height-to-width ratio" << std::endl;
  }
}

void SetFracPervRoadOfTotalRoad(UrbanType urban, int numLandunits,
                                int mpi_rank) {
  UrbanErrorCode ierr;

  double *fracPervRoadOfTotalRoad =
      AllocateArray(numLandunits, "fracPervRoadOfTotalRoad");
  for (int i = 0; i < numLandunits; ++i) {
    fracPervRoadOfTotalRoad[i] = 0.16666667163372040;
  }
  UrbanCall(UrbanSetFracPervRoadOfTotalRoad(urban, fracPervRoadOfTotalRoad,
                                            numLandunits, &ierr),
            &ierr);
  free(fracPervRoadOfTotalRoad);

  if (mpi_rank == 0) {
    std::cout << "Set fraction of pervious road w.r.t. total road" << std::endl;
  }
}

void SetHeightParameters(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  // Height parameter default values (in meters)
  const double FORC_HGT_T_DEFAULT = 144.44377627618979;
  const double FORC_HGT_U_DEFAULT = 144.44377627618979; // Same as T
  const double Z_D_TOWN_DEFAULT = 113.96331622200367;
  const double Z_0_TOWN_DEFAULT = 0.48046005418613641;
  const double HT_ROOF_DEFAULT = 120.0;
  const double WIND_HGT_CANYON_DEFAULT = 60.0;

  double *forcHgtT = AllocateArray(numLandunits, "forcHgtT");
  double *forcHgtU = AllocateArray(numLandunits, "forcHgtU");
  double *zDTown = AllocateArray(numLandunits, "zDTown");
  double *z0Town = AllocateArray(numLandunits, "z0Town");
  double *htRoof = AllocateArray(numLandunits, "htRoof");
  double *windHgtCanyon = AllocateArray(numLandunits, "windHgtCanyon");

  for (int i = 0; i < numLandunits; ++i) {
    forcHgtT[i] = FORC_HGT_T_DEFAULT;
    forcHgtU[i] = FORC_HGT_U_DEFAULT;
    zDTown[i] = Z_D_TOWN_DEFAULT;
    z0Town[i] = Z_0_TOWN_DEFAULT;
    htRoof[i] = HT_ROOF_DEFAULT;
    windHgtCanyon[i] = WIND_HGT_CANYON_DEFAULT;
  }

  UrbanCall(UrbanSetForcHgtT(urban, forcHgtT, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetForcHgtU(urban, forcHgtU, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetZDTown(urban, zDTown, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetZ0Town(urban, z0Town, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetHtRoof(urban, htRoof, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetWindHgtCanyon(urban, windHgtCanyon, numLandunits, &ierr),
            &ierr);

  double *heightArrays[] = {forcHgtT, forcHgtU, zDTown,
                            z0Town,   htRoof,   windHgtCanyon};
  FreeArrays(heightArrays, 6);

  if (mpi_rank == 0) {
    std::cout << "Set height parameters:" << std::endl;
    std::cout << "  Forcing height (T): " << FORC_HGT_T_DEFAULT << " m"
              << std::endl;
    std::cout << "  Forcing height (U): " << FORC_HGT_U_DEFAULT << " m"
              << std::endl;
    std::cout << "  Zero displacement (town): " << Z_D_TOWN_DEFAULT << " m"
              << std::endl;
    std::cout << "  Roughness length (town): " << Z_0_TOWN_DEFAULT << " m"
              << std::endl;
    std::cout << "  Roof height: " << HT_ROOF_DEFAULT << " m" << std::endl;
    std::cout << "  Wind height (canyon): " << WIND_HGT_CANYON_DEFAULT << " m"
              << std::endl;
  }
}

void SetWtRoof(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  const double WT_ROOF_DEFAULT = 0.69999998807907104;

  double *wtRoof = AllocateArray(numLandunits, "wtRoof");
  for (int i = 0; i < numLandunits; ++i) {
    wtRoof[i] = WT_ROOF_DEFAULT;
  }
  UrbanCall(UrbanSetWtRoof(urban, wtRoof, numLandunits, &ierr), &ierr);
  free(wtRoof);

  if (mpi_rank == 0) {
    std::cout << "Set roof weight: " << WT_ROOF_DEFAULT << std::endl;
  }
}

void SetAlbedo(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  int numBands = 2; // VIS, NIR
  int numTypes = 2; // Direct, Diffuse
  int size3D[3] = {numLandunits, numBands, numTypes};
  int totalSize3D = numLandunits * numBands * numTypes;

  const double ALB_IMPROAD = 0.230000004172325;
  const double ALB_PERROAD = 0.0799999982118607;
  const double ALB_ROOF = 0.254999995231628;
  const double ALB_WALL = 0.200000002980232;

  double *albedoPerviousRoad = AllocateArray(totalSize3D, "albedoPerviousRoad");
  double *albedoImperviousRoad =
      AllocateArray(totalSize3D, "albedoImperviousRoad");
  double *albedoSunlitWall = AllocateArray(totalSize3D, "albedoSunlitWall");
  double *albedoShadedWall = AllocateArray(totalSize3D, "albedoShadedWall");
  double *albedoRoof = AllocateArray(totalSize3D, "albedoRoof");

  for (int ilandunit = 0; ilandunit < numLandunits; ++ilandunit) {
    for (int iband = 0; iband < numBands; ++iband) {
      for (int itype = 0; itype < numTypes; ++itype) {
        int idx = ilandunit * numBands * numTypes + iband * numTypes + itype;
        albedoPerviousRoad[idx] = ALB_PERROAD;
        albedoImperviousRoad[idx] = ALB_IMPROAD;
        albedoSunlitWall[idx] = ALB_WALL;
        albedoShadedWall[idx] = ALB_WALL;
        albedoRoof[idx] = ALB_ROOF;
      }
    }
  }

  UrbanCall(
      UrbanSetAlbedoPerviousRoad(urban, albedoPerviousRoad, size3D, &ierr),
      &ierr);
  UrbanCall(
      UrbanSetAlbedoImperviousRoad(urban, albedoImperviousRoad, size3D, &ierr),
      &ierr);
  UrbanCall(UrbanSetAlbedoSunlitWall(urban, albedoSunlitWall, size3D, &ierr),
            &ierr);
  UrbanCall(UrbanSetAlbedoShadedWall(urban, albedoShadedWall, size3D, &ierr),
            &ierr);
  UrbanCall(UrbanSetAlbedoRoof(urban, albedoRoof, size3D, &ierr), &ierr);

  double *albedoArrays[] = {albedoPerviousRoad, albedoImperviousRoad,
                            albedoSunlitWall, albedoShadedWall, albedoRoof};
  FreeArrays(albedoArrays, 5);

  if (mpi_rank == 0) {
    std::cout << "Set albedo values for all surfaces" << std::endl;
  }
}

void SetEmissivity(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  double *emissivityPerviousRoad =
      AllocateArray(numLandunits, "emissivityPerviousRoad");
  double *emissivityImperviousRoad =
      AllocateArray(numLandunits, "emissivityImperviousRoad");
  double *emissivityWall = AllocateArray(numLandunits, "emissivityWall");
  double *emissivityRoof = AllocateArray(numLandunits, "emissivityRoof");

  const double EMISS_ROOF = 0.90600001811981201;
  const double EMISS_IMPROAD = 0.87999999523162842;
  const double EMISS_PERROAD = 0.94999998807907104;
  const double EMISS_WALL = 0.90200001001358032;

  for (int i = 0; i < numLandunits; ++i) {
    emissivityPerviousRoad[i] = EMISS_PERROAD;
    emissivityImperviousRoad[i] = EMISS_IMPROAD;
    emissivityWall[i] = EMISS_WALL;
    emissivityRoof[i] = EMISS_ROOF;
  }

  UrbanCall(UrbanSetEmissivityPerviousRoad(urban, emissivityPerviousRoad,
                                           numLandunits, &ierr),
            &ierr);
  UrbanCall(UrbanSetEmissivityImperviousRoad(urban, emissivityImperviousRoad,
                                             numLandunits, &ierr),
            &ierr);
  UrbanCall(UrbanSetEmissivityWall(urban, emissivityWall, numLandunits, &ierr),
            &ierr);
  UrbanCall(UrbanSetEmissivityRoof(urban, emissivityRoof, numLandunits, &ierr),
            &ierr);

  double *emissivityArrays[] = {emissivityPerviousRoad,
                                emissivityImperviousRoad, emissivityWall,
                                emissivityRoof};
  FreeArrays(emissivityArrays, 4);

  if (mpi_rank == 0) {
    std::cout << "Set emissivity values for all surfaces" << std::endl;
  }
}

void SetNumberOfActiveLayersImperviousRoad(UrbanType urban, int numLandunits,
                                           int mpi_rank) {
  UrbanErrorCode ierr;

  // Number of active layers for each urban density class
  // 0=Tall Building District, 1=High Density, 2=Medium Density
  int nlevImproadClasses[NUM_URBAN_DENSITY_CLASSES] = {3, 2, 2};

  double *numActiveLayers = AllocateArray(numLandunits, "numActiveLayers");
  for (int i = 0; i < numLandunits; ++i) {
    const int urban_density_class = i % NUM_URBAN_DENSITY_CLASSES;
    numActiveLayers[i] =
        static_cast<double>(nlevImproadClasses[urban_density_class]);
  }
  UrbanCall(UrbanSetNumberOfActiveLayersImperviousRoad(urban, numActiveLayers,
                                                       numLandunits, &ierr),
            &ierr);
  free(numActiveLayers);

  if (mpi_rank == 0) {
    std::cout << "Set number of active layers for impervious road (3, 2, 2 for "
                 "density classes)"
              << std::endl;
  }
}

void SetThermalConductivity(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  const int NUM_ROAD_LEVELS = 15;
  const int NUM_URBAN_LEVELS = 5;
  int size2D_road[2] = {numLandunits, NUM_ROAD_LEVELS};
  int size2D_urban[2] = {numLandunits, NUM_URBAN_LEVELS};
  int totalSizeRoad = numLandunits * NUM_ROAD_LEVELS;
  int totalSizeUrban = numLandunits * NUM_URBAN_LEVELS;

  double *tkRoad = AllocateArray(totalSizeRoad, "tkRoad");
  double *tkWall = AllocateArray(totalSizeUrban, "tkWall");
  double *tkRoof = AllocateArray(totalSizeUrban, "tkRoof");

  // Thermal conductivity values for road (15 levels) across
  // NUM_URBAN_DENSITY_CLASSES urban density classes
  // [layer][urban_density_class]: 0=Tall Building District, 1=High Density,
  // 2=Medium Density Active layers use specified values, inactive layers use
  // bedrock thermal conductivity (3.0 W/m-K)
  const double TK_BEDROCK = 3.0; // thermal conductivity of bedrock [W/m-K]
  double tkRoadLevels[NUM_ROAD_LEVELS][NUM_URBAN_DENSITY_CLASSES] = {
      {1.8999999761581421, 1.6699999570846558, 1.6699999570846558},
      {0.56000000238418579, 0.56000000238418579, 0.56000000238418579},
      {0.36000001430511475, 0.21664454245689402, 0.21664454245689402},
      {0.21572236500511191, 0.21572236500511191, 0.21572236500511191},
      {0.21389214262493184, 0.21389214262493184, 0.21389214262493184},
      {0.21208052418425166, 0.21208052418425166, 0.21208052418425166},
      {0.21028722745802339, 0.21028722745802339, 0.21028722745802339},
      {0.21028722745802339, 0.21028722745802339, 0.21028722745802339},
      {0.21298402568603625, 0.21298402568603625, 0.21298402568603625},
      {0.21664454245689402, 0.21664454245689402, 0.21664454245689402},
      {TK_BEDROCK, TK_BEDROCK, TK_BEDROCK},
      {TK_BEDROCK, TK_BEDROCK, TK_BEDROCK},
      {TK_BEDROCK, TK_BEDROCK, TK_BEDROCK},
      {TK_BEDROCK, TK_BEDROCK, TK_BEDROCK},
      {TK_BEDROCK, TK_BEDROCK, TK_BEDROCK}};
  double tkWallLevels[NUM_URBAN_LEVELS][NUM_URBAN_DENSITY_CLASSES] = {
      {1.44716906547546, 1.06582415103912, 0.970157384872437},
      {1.44716906547546, 1.06582415103912, 0.970157384872437},
      {1.44716906547546, 1.06582415103912, 0.970157384872437},
      {1.44716906547546, 1.06582415103912, 0.970157384872437},
      {1.44716906547546, 1.06582415103912, 0.970157384872437}};
  double tkRoofLevels[NUM_URBAN_LEVELS][NUM_URBAN_DENSITY_CLASSES] = {
      {0.503093481063843, 0.094768725335598, 0.127733826637268},
      {0.503093481063843, 0.094768725335598, 0.127733826637268},
      {0.503093481063843, 0.094768725335598, 0.127733826637268},
      {0.503093481063843, 0.094768725335598, 0.127733826637268},
      {0.503093481063843, 0.094768725335598, 0.127733826637268}};

  // Fill road array
  int idx = 0;
  for (int layer = 0; layer < NUM_ROAD_LEVELS; ++layer) {
    for (int i = 0; i < numLandunits; ++i) {
      const int urban_density_class = i % NUM_URBAN_DENSITY_CLASSES;
      tkRoad[idx] = tkRoadLevels[layer][urban_density_class];
      ++idx;
    }
  }

  // Fill wall and roof arrays
  idx = 0;
  for (int layer = 0; layer < NUM_URBAN_LEVELS; ++layer) {
    for (int i = 0; i < numLandunits; ++i) {
      const int urban_density_class = i % NUM_URBAN_DENSITY_CLASSES;
      tkWall[idx] = tkWallLevels[layer][urban_density_class];
      tkRoof[idx] = tkRoofLevels[layer][urban_density_class];
      ++idx;
    }
  }

  UrbanCall(UrbanSetThermalConductivityRoad(urban, tkRoad, size2D_road, &ierr),
            &ierr);
  UrbanCall(UrbanSetThermalConductivityWall(urban, tkWall, size2D_urban, &ierr),
            &ierr);
  UrbanCall(UrbanSetThermalConductivityRoof(urban, tkRoof, size2D_urban, &ierr),
            &ierr);

  double *tkArrays[] = {tkRoad, tkWall, tkRoof};
  FreeArrays(tkArrays, 3);

  if (mpi_rank == 0) {
    std::cout << "Set thermal conductivity values for all surfaces"
              << std::endl;
  }
}

void SetHeatCapacity(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  const int NUM_ROAD_LEVELS = 15;
  const int NUM_URBAN_LEVELS = 5;
  int size2D_road[2] = {numLandunits, NUM_ROAD_LEVELS};
  int size2D_urban[2] = {numLandunits, NUM_URBAN_LEVELS};
  int totalSizeRoad = numLandunits * NUM_ROAD_LEVELS;
  int totalSizeUrban = numLandunits * NUM_URBAN_LEVELS;

  double *cvRoad = AllocateArray(totalSizeRoad, "cvRoad");
  double *cvWall = AllocateArray(totalSizeUrban, "cvWall");
  double *cvRoof = AllocateArray(totalSizeUrban, "cvRoof");

  // Heat capacity values for road (15 levels) across NUM_URBAN_DENSITY_CLASSES
  // urban density classes [layer][urban_density_class]: 0=Tall Building
  // District, 1=High Density, 2=Medium Density Active layers use specified
  // values, inactive layers use bedrock heat capacity (2.0e6 J/m^3/K)
  const double CV_BEDROCK = 2.0e6; // heat capacity of bedrock [J/m^3/K]
  double cvRoadLevels[NUM_ROAD_LEVELS][NUM_URBAN_DENSITY_CLASSES] = {
      {2100000.0, 2060470.625, 2060470.625},
      {1773000.0, 1712294.75, 1712294.75},
      {1545600.0, CV_BEDROCK, CV_BEDROCK},  // Layer 2: only active for density
                                            // class 0
      {CV_BEDROCK, CV_BEDROCK, CV_BEDROCK}, // Remaining layers use bedrock
      {CV_BEDROCK, CV_BEDROCK, CV_BEDROCK},
      {CV_BEDROCK, CV_BEDROCK, CV_BEDROCK},
      {CV_BEDROCK, CV_BEDROCK, CV_BEDROCK},
      {CV_BEDROCK, CV_BEDROCK, CV_BEDROCK},
      {CV_BEDROCK, CV_BEDROCK, CV_BEDROCK},
      {CV_BEDROCK, CV_BEDROCK, CV_BEDROCK},
      {CV_BEDROCK, CV_BEDROCK, CV_BEDROCK},
      {CV_BEDROCK, CV_BEDROCK, CV_BEDROCK},
      {CV_BEDROCK, CV_BEDROCK, CV_BEDROCK},
      {CV_BEDROCK, CV_BEDROCK, CV_BEDROCK},
      {CV_BEDROCK, CV_BEDROCK, CV_BEDROCK}};
  double cvWallLevels[NUM_URBAN_LEVELS][NUM_URBAN_DENSITY_CLASSES] = {
      {1079394.75, 957632.8125, 899827.1875},
      {1079394.75, 957632.8125, 899827.1875},
      {1079394.75, 957632.8125, 899827.1875},
      {1079394.75, 957632.8125, 899827.1875},
      {1079394.75, 957632.8125, 899827.1875}};
  double cvRoofLevels[NUM_URBAN_LEVELS][NUM_URBAN_DENSITY_CLASSES] = {
      {570998.0, 646213.375, 862451.375},
      {570998.0, 646213.375, 862451.375},
      {570998.0, 646213.375, 862451.375},
      {570998.0, 646213.375, 862451.375},
      {570998.0, 646213.375, 862451.375}};

  // Fill road array
  int idx = 0;
  for (int layer = 0; layer < NUM_ROAD_LEVELS; ++layer) {
    for (int i = 0; i < numLandunits; ++i) {
      const int urban_density_class = i % NUM_URBAN_DENSITY_CLASSES;
      cvRoad[idx] = cvRoadLevels[layer][urban_density_class];
      ++idx;
    }
  }

  // Fill wall and roof arrays
  idx = 0;
  for (int layer = 0; layer < NUM_URBAN_LEVELS; ++layer) {
    for (int i = 0; i < numLandunits; ++i) {
      const int urban_density_class = i % NUM_URBAN_DENSITY_CLASSES;
      cvWall[idx] = cvWallLevels[layer][urban_density_class];
      cvRoof[idx] = cvRoofLevels[layer][urban_density_class];
      ++idx;
    }
  }

  UrbanCall(UrbanSetHeatCapacityRoad(urban, cvRoad, size2D_road, &ierr), &ierr);
  UrbanCall(UrbanSetHeatCapacityWall(urban, cvWall, size2D_urban, &ierr),
            &ierr);
  UrbanCall(UrbanSetHeatCapacityRoof(urban, cvRoof, size2D_urban, &ierr),
            &ierr);

  double *cvArrays[] = {cvRoad, cvWall, cvRoof};
  FreeArrays(cvArrays, 3);

  if (mpi_rank == 0) {
    std::cout << "Set heat capacity values for all surfaces" << std::endl;
  }
}

void SetSoilProperties(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  const int NUM_SOIL_LEVELS = 15;
  int size2D[2] = {numLandunits, NUM_SOIL_LEVELS};
  int totalSize = numLandunits * NUM_SOIL_LEVELS;

  double *sand = AllocateArray(totalSize, "sand");
  double *clay = AllocateArray(totalSize, "clay");
  double *organic = AllocateArray(totalSize, "organic");

  // Soil property values for first 10 layers (layers 11-15 use layer 10 values)
  double sandLevels[10] = {46.0, 46.0, 44.0, 43.0, 41.0,
                           39.0, 37.0, 37.0, 40.0, 44.0};
  double clayLevels[10] = {35.0, 35.0, 37.0, 39.0, 42.0,
                           44.0, 46.0, 45.0, 41.0, 42.0};
  double organicLevels[10] = {25.2229820327902,
                              25.700711396596,
                              22.091324741929,
                              18.1150405358844,
                              14.5211498497041,
                              11.4998502546828,
                              9.04501744160207,
                              7.08594278159189,
                              0.0,
                              0.0};
  int idx = 0;
  for (int layer = 0; layer < NUM_SOIL_LEVELS; ++layer) {
    for (int i = 0; i < numLandunits; ++i) {
      // Use layer 9 (10th layer) values for layers 10-14 (11th-15th)
      int srcLayer = (layer < 10) ? layer : 9;
      sand[idx] = sandLevels[srcLayer];
      clay[idx] = clayLevels[srcLayer];
      organic[idx] = organicLevels[srcLayer];
      ++idx;
    }
  }

  UrbanCall(UrbanSetSandPerviousRoad(urban, sand, size2D, &ierr), &ierr);
  UrbanCall(UrbanSetClayPerviousRoad(urban, clay, size2D, &ierr), &ierr);
  UrbanCall(UrbanSetOrganicPerviousRoad(urban, organic, size2D, &ierr), &ierr);

  double *soilArrays[] = {sand, clay, organic};
  FreeArrays(soilArrays, 3);

  if (mpi_rank == 0) {
    std::cout << "Set soil properties for pervious road (sand, clay, organic)"
              << std::endl;
  }
}

void SetAtmosphericForcing(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  // Atmospheric forcing constants
  const double TEMP_AIR = 297.26422743678319;
  const double TH_AIR = 297.26422743678319;
  const double RHO_AIR = 1.1382761848551157;
  const double Q_AIR = 1.9217052569985755E-002;
  const double PBOT_AIR = 98260.450580263219;
  const double WIND_U = 0.52482489069830152;
  const double WIND_V = 0.0;
  const double COSZEN = 7.9054122593736065E-003;
  const double SNOW = 0.0;
  const double LWDOWN = 432.79580327766450;
  const double SWDOWN = 0.0;

  // Allocate arrays
  double *atmTemp = AllocateArray(numLandunits, "atmTemp");
  double *atmPotTemp = AllocateArray(numLandunits, "atmPotTemp");
  double *atmRho = AllocateArray(numLandunits, "atmRho");
  double *atmSpcHumd = AllocateArray(numLandunits, "atmSpcHumd");
  double *atmPress = AllocateArray(numLandunits, "atmPress");
  double *atmWindU = AllocateArray(numLandunits, "atmWindU");
  double *atmWindV = AllocateArray(numLandunits, "atmWindV");
  double *atmCoszen = AllocateArray(numLandunits, "atmCoszen");
  double *atmFracSnow = AllocateArray(numLandunits, "atmFracSnow");
  double *atmLongwave = AllocateArray(numLandunits, "atmLongwave");

  int numBands = 2; // VIS, NIR
  int numTypes = 2; // Direct, Diffuse
  int size3D[3] = {numLandunits, numBands, numTypes};
  int totalSize3D = numLandunits * numBands * numTypes;
  double *atmShortwave = AllocateArray(totalSize3D, "atmShortwave");

  // Fill arrays with constant values
  for (int i = 0; i < numLandunits; ++i) {
    atmTemp[i] = TEMP_AIR;
    atmPotTemp[i] = TH_AIR;
    atmRho[i] = RHO_AIR;
    atmSpcHumd[i] = Q_AIR;
    atmPress[i] = PBOT_AIR;
    atmWindU[i] = WIND_U;
    atmWindV[i] = WIND_V;
    atmCoszen[i] = COSZEN;
    atmFracSnow[i] = SNOW;
    atmLongwave[i] = LWDOWN;
  }

  for (int i = 0; i < totalSize3D; ++i) {
    atmShortwave[i] = SWDOWN;
  }

  // Set atmospheric forcing
  UrbanCall(UrbanSetAtmTemp(urban, atmTemp, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetAtmPotTemp(urban, atmPotTemp, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetAtmRho(urban, atmRho, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetAtmSpcHumd(urban, atmSpcHumd, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetAtmPress(urban, atmPress, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetAtmWindU(urban, atmWindU, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetAtmWindV(urban, atmWindV, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetAtmCoszen(urban, atmCoszen, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetAtmFracSnow(urban, atmFracSnow, numLandunits, &ierr),
            &ierr);
  UrbanCall(UrbanSetAtmLongwaveDown(urban, atmLongwave, numLandunits, &ierr),
            &ierr);
  UrbanCall(UrbanSetAtmShortwaveDown(urban, atmShortwave, size3D, &ierr),
            &ierr);

  // Free arrays
  double *atmArrays[] = {atmTemp,     atmPotTemp,  atmRho,      atmSpcHumd,
                         atmPress,    atmWindU,    atmWindV,    atmCoszen,
                         atmFracSnow, atmLongwave, atmShortwave};
  FreeArrays(atmArrays, 11);

  if (mpi_rank == 0) {
    std::cout << "Set atmospheric forcing:" << std::endl;
    std::cout << "  Temperature: " << TEMP_AIR << " K" << std::endl;
    std::cout << "  Pressure: " << PBOT_AIR << " Pa" << std::endl;
    std::cout << "  Density: " << RHO_AIR << " kg/m^3" << std::endl;
    std::cout << "  Specific humidity: " << Q_AIR << " kg/kg" << std::endl;
    std::cout << "  Wind U: " << WIND_U << " m/s" << std::endl;
    std::cout << "  Wind V: " << WIND_V << " m/s" << std::endl;
    std::cout << "  Cosine zenith: " << COSZEN << std::endl;
    std::cout << "  Snow fraction: " << SNOW << std::endl;
    std::cout << "  Longwave down: " << LWDOWN << " W/m^2" << std::endl;
    std::cout << "  Shortwave down: " << SWDOWN << " W/m^2" << std::endl;
  }
}

void SetBuildingTemperature(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  const double MIN_TEMP = 285.0; // K
  const double MAX_TEMP = 310.0; // K

  double *minTemp = AllocateArray(numLandunits, "minTemp");
  double *maxTemp = AllocateArray(numLandunits, "maxTemp");

  for (int i = 0; i < numLandunits; ++i) {
    minTemp[i] = MIN_TEMP;
    maxTemp[i] = MAX_TEMP;
  }

  UrbanCall(UrbanSetBuildingMinTemperature(urban, minTemp, numLandunits, &ierr),
            &ierr);
  UrbanCall(UrbanSetBuildingMaxTemperature(urban, maxTemp, numLandunits, &ierr),
            &ierr);

  double *tempArrays[] = {minTemp, maxTemp};
  FreeArrays(tempArrays, 2);

  if (mpi_rank == 0) {
    std::cout << "Set building temperature limits:" << std::endl;
    std::cout << "  Min temperature: " << MIN_TEMP << " K" << std::endl;
    std::cout << "  Max temperature: " << MAX_TEMP << " K" << std::endl;
  }

  // Set building thickness parameters
  const double THICK_WALL = 0.286199986934662; // m
  const double THICK_ROOF = 0.217099994421005; // m

  double *wallThickness = AllocateArray(numLandunits, "wallThickness");
  double *roofThickness = AllocateArray(numLandunits, "roofThickness");

  for (int i = 0; i < numLandunits; ++i) {
    wallThickness[i] = THICK_WALL;
    roofThickness[i] = THICK_ROOF;
  }

  UrbanCall(
      UrbanSetBuildingWallThickness(urban, wallThickness, numLandunits, &ierr),
      &ierr);
  UrbanCall(
      UrbanSetBuildingRoofThickness(urban, roofThickness, numLandunits, &ierr),
      &ierr);

  double *thicknessArrays[] = {wallThickness, roofThickness};
  FreeArrays(thicknessArrays, 2);

  if (mpi_rank == 0) {
    std::cout << "Set building thickness parameters:" << std::endl;
    std::cout << "  Wall thickness: " << THICK_WALL << " m" << std::endl;
    std::cout << "  Roof thickness: " << THICK_ROOF << " m" << std::endl;
  }
}

void SetUrbanParameters(UrbanType urban, int numLandunits, int mpi_rank) {
  SetCanyonHwr(urban, numLandunits, mpi_rank);
  SetFracPervRoadOfTotalRoad(urban, numLandunits, mpi_rank);
  SetWtRoof(urban, numLandunits, mpi_rank);
  SetHeightParameters(urban, numLandunits, mpi_rank);
  SetBuildingTemperature(urban, numLandunits, mpi_rank);
  SetAlbedo(urban, numLandunits, mpi_rank);
  SetEmissivity(urban, numLandunits, mpi_rank);
  SetNumberOfActiveLayersImperviousRoad(urban, numLandunits, mpi_rank);
  SetThermalConductivity(urban, numLandunits, mpi_rank);
  SetHeatCapacity(urban, numLandunits, mpi_rank);
  SetSoilProperties(urban, numLandunits, mpi_rank);
  SetAtmosphericForcing(urban, numLandunits, mpi_rank);
}

int main(int argc, char *argv[]) {

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int mpi_rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  Kokkos::initialize(argc, argv);
  {
    if (mpi_rank == 0) {
      std::cout << "=== Kokkos Configuration (rank 0 only) ===" << std::endl;
      Kokkos::print_configuration(std::cout);
    }

    // Create Urban object
    int numLandunits = 10;
    UrbanType urban = nullptr;
    UrbanErrorCode ierr;
    UrbanCall(UrbanCreate(numLandunits, &urban, &ierr), &ierr);

    if (mpi_rank == 0) {
      std::cout << "Successfully created Urban object with " << numLandunits
                << " landunits" << std::endl;
    }

    // Set all urban parameters
    SetUrbanParameters(urban, numLandunits, mpi_rank);

    // Setup urban model (initialize temperatures and other setup tasks)
    UrbanCall(UrbanSetup(urban, &ierr), &ierr);
    if (mpi_rank == 0) {
      std::cout << "Completed urban model setup" << std::endl;
    }

    // Advance the model one time step
    UrbanCall(UrbanAdvance(urban, &ierr), &ierr);
    if (mpi_rank == 0) {
      std::cout << "Advanced model one time step" << std::endl;
    }

    // Destroy Urban object
    UrbanCall(UrbanDestroy(&urban, &ierr), &ierr);
  }

  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}
