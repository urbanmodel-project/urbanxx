#include <Kokkos_Core.hpp>
#include <iostream>
#include <mpi.h>
#include <stdlib.h>
#include <Urban.h>
#include <UrbanMacros.h>

// Helper function to allocate memory with error checking
static double* AllocateArray(int size, const char* name) {
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

void SetFracPervRoadOfTotalRoad(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  double *fracPervRoadOfTotalRoad = AllocateArray(numLandunits, "fracPervRoadOfTotalRoad");
  for (int i = 0; i < numLandunits; ++i) {
    fracPervRoadOfTotalRoad[i] = 0.16666667163372040;
  }
  UrbanCall(UrbanSetFracPervRoadOfTotalRoad(urban, fracPervRoadOfTotalRoad, numLandunits, &ierr), &ierr);
  free(fracPervRoadOfTotalRoad);
  
  if (mpi_rank == 0) {
    std::cout << "Set fraction of pervious road w.r.t. total road" << std::endl;
  }
}

void SetHeightParameters(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  // Height parameter default values (in meters)
  const double FORC_HGT_T_DEFAULT = 144.44377627618979;
  const double FORC_HGT_U_DEFAULT = 144.44377627618979;  // Same as T
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
  UrbanCall(UrbanSetWindHgtCanyon(urban, windHgtCanyon, numLandunits, &ierr), &ierr);

  double *heightArrays[] = {forcHgtT, forcHgtU, zDTown, z0Town, htRoof, windHgtCanyon};
  FreeArrays(heightArrays, 6);

  if (mpi_rank == 0) {
    std::cout << "Set height parameters:" << std::endl;
    std::cout << "  Forcing height (T): " << FORC_HGT_T_DEFAULT << " m" << std::endl;
    std::cout << "  Forcing height (U): " << FORC_HGT_U_DEFAULT << " m" << std::endl;
    std::cout << "  Zero displacement (town): " << Z_D_TOWN_DEFAULT << " m" << std::endl;
    std::cout << "  Roughness length (town): " << Z_0_TOWN_DEFAULT << " m" << std::endl;
    std::cout << "  Roof height: " << HT_ROOF_DEFAULT << " m" << std::endl;
    std::cout << "  Wind height (canyon): " << WIND_HGT_CANYON_DEFAULT << " m" << std::endl;
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

  int numBands = 2;    // VIS, NIR
  int numTypes = 2;    // Direct, Diffuse
  int size3D[3] = {numLandunits, numBands, numTypes};
  int totalSize3D = numLandunits * numBands * numTypes;

  const double ALB_IMPROAD = 0.230000004172325;
  const double ALB_PERROAD = 0.0799999982118607;
  const double ALB_ROOF = 0.254999995231628;
  const double ALB_WALL = 0.200000002980232;

  double *albedoPerviousRoad = AllocateArray(totalSize3D, "albedoPerviousRoad");
  double *albedoImperviousRoad = AllocateArray(totalSize3D, "albedoImperviousRoad");
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

  UrbanCall(UrbanSetAlbedoPerviousRoad(urban, albedoPerviousRoad, size3D, &ierr), &ierr);
  UrbanCall(UrbanSetAlbedoImperviousRoad(urban, albedoImperviousRoad, size3D, &ierr), &ierr);
  UrbanCall(UrbanSetAlbedoSunlitWall(urban, albedoSunlitWall, size3D, &ierr), &ierr);
  UrbanCall(UrbanSetAlbedoShadedWall(urban, albedoShadedWall, size3D, &ierr), &ierr);
  UrbanCall(UrbanSetAlbedoRoof(urban, albedoRoof, size3D, &ierr), &ierr);

  double *albedoArrays[] = {albedoPerviousRoad, albedoImperviousRoad, albedoSunlitWall, albedoShadedWall, albedoRoof};
  FreeArrays(albedoArrays, 5);

  if (mpi_rank == 0) {
    std::cout << "Set albedo values for all surfaces" << std::endl;
  }
}

void SetEmissivity(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  double *emissivityPerviousRoad = AllocateArray(numLandunits, "emissivityPerviousRoad");
  double *emissivityImperviousRoad = AllocateArray(numLandunits, "emissivityImperviousRoad");
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

  UrbanCall(UrbanSetEmissivityPerviousRoad(urban, emissivityPerviousRoad, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetEmissivityImperviousRoad(urban, emissivityImperviousRoad, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetEmissivityWall(urban, emissivityWall, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetEmissivityRoof(urban, emissivityRoof, numLandunits, &ierr), &ierr);

  double *emissivityArrays[] = {emissivityPerviousRoad, emissivityImperviousRoad, emissivityWall, emissivityRoof};
  FreeArrays(emissivityArrays, 4);

  if (mpi_rank == 0) {
    std::cout << "Set emissivity values for all surfaces" << std::endl;
  }
}

void SetThermalConductivity(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  double *tkRoad = AllocateArray(numLandunits, "tkRoad");
  double *tkWall = AllocateArray(numLandunits, "tkWall");
  double *tkRoof = AllocateArray(numLandunits, "tkRoof");

  for (int i = 0; i < numLandunits; ++i) {
    tkRoad[i] = 1.0;
    tkWall[i] = 0.8;
    tkRoof[i] = 0.9;
  }

  UrbanCall(UrbanSetThermalConductivityRoad(urban, tkRoad, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetThermalConductivityWall(urban, tkWall, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetThermalConductivityRoof(urban, tkRoof, numLandunits, &ierr), &ierr);

  double *tkArrays[] = {tkRoad, tkWall, tkRoof};
  FreeArrays(tkArrays, 3);

  if (mpi_rank == 0) {
    std::cout << "Set thermal conductivity values for all surfaces" << std::endl;
  }
}

void SetHeatCapacity(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  double *cvRoad = AllocateArray(numLandunits, "cvRoad");
  double *cvWall = AllocateArray(numLandunits, "cvWall");
  double *cvRoof = AllocateArray(numLandunits, "cvRoof");

  for (int i = 0; i < numLandunits; ++i) {
    cvRoad[i] = 2.0e6;
    cvWall[i] = 1.8e6;
    cvRoof[i] = 1.9e6;
  }

  UrbanCall(UrbanSetHeatCapacityRoad(urban, cvRoad, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetHeatCapacityWall(urban, cvWall, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetHeatCapacityRoof(urban, cvRoof, numLandunits, &ierr), &ierr);

  double *cvArrays[] = {cvRoad, cvWall, cvRoof};
  FreeArrays(cvArrays, 3);

  if (mpi_rank == 0) {
    std::cout << "Set heat capacity values for all surfaces" << std::endl;
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
  const double SWDOWN = 1.0;

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

  int numBands = 2;    // VIS, NIR
  int numTypes = 2;    // Direct, Diffuse
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
  UrbanCall(UrbanSetAtmFracSnow(urban, atmFracSnow, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetAtmLongwaveDown(urban, atmLongwave, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetAtmShortwaveDown(urban, atmShortwave, size3D, &ierr), &ierr);

  // Free arrays
  double *atmArrays[] = {atmTemp, atmPotTemp, atmRho, atmSpcHumd, atmPress,
                         atmWindU, atmWindV, atmCoszen, atmFracSnow,
                         atmLongwave, atmShortwave};
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

void SetUrbanParameters(UrbanType urban, int numLandunits, int mpi_rank) {
  SetCanyonHwr(urban, numLandunits, mpi_rank);
  SetFracPervRoadOfTotalRoad(urban, numLandunits, mpi_rank);
  SetWtRoof(urban, numLandunits, mpi_rank);
  SetHeightParameters(urban, numLandunits, mpi_rank);
  SetAlbedo(urban, numLandunits, mpi_rank);
  SetEmissivity(urban, numLandunits, mpi_rank);
  SetThermalConductivity(urban, numLandunits, mpi_rank);
  SetHeatCapacity(urban, numLandunits, mpi_rank);
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

    // Initialize temperatures
    UrbanCall(UrbanInitializeTemperature(urban, &ierr), &ierr);
    if (mpi_rank == 0) {
      std::cout << "Initialized surface temperatures" << std::endl;
    }

    // Set all urban parameters
    SetUrbanParameters(urban, numLandunits, mpi_rank);

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
