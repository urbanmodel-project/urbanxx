#include <Kokkos_Core.hpp>
#include <iostream>
#include <mpi.h>
#include <Urban.h>
#include <UrbanMacros.h>

void SetCanyonHwr(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  double canyonHwr[numLandunits];
  for (int i = 0; i < numLandunits; ++i) {
    canyonHwr[i] = 4.80000019073486;
  }
  UrbanCall(UrbanSetCanyonHwr(urban, canyonHwr, numLandunits, &ierr), &ierr);
  
  if (mpi_rank == 0) {
    std::cout << "Set canyon height-to-width ratio" << std::endl;
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

  double albedoPerviousRoad[totalSize3D];
  double albedoImperviousRoad[totalSize3D];
  double albedoSunlitWall[totalSize3D];
  double albedoShadedWall[totalSize3D];
  double albedoRoof[totalSize3D];

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

  if (mpi_rank == 0) {
    std::cout << "Set albedo values for all surfaces" << std::endl;
  }
}

void SetEmissivity(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  double emissivityPerviousRoad[numLandunits];
  double emissivityImperviousRoad[numLandunits];
  double emissivityWall[numLandunits];
  double emissivityRoof[numLandunits];

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

  if (mpi_rank == 0) {
    std::cout << "Set emissivity values for all surfaces" << std::endl;
  }
}

void SetThermalConductivity(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  double tkRoad[numLandunits];
  double tkWall[numLandunits];
  double tkRoof[numLandunits];

  for (int i = 0; i < numLandunits; ++i) {
    tkRoad[i] = 1.0;
    tkWall[i] = 0.8;
    tkRoof[i] = 0.9;
  }

  UrbanCall(UrbanSetThermalConductivityRoad(urban, tkRoad, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetThermalConductivityWall(urban, tkWall, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetThermalConductivityRoof(urban, tkRoof, numLandunits, &ierr), &ierr);

  if (mpi_rank == 0) {
    std::cout << "Set thermal conductivity values for all surfaces" << std::endl;
  }
}

void SetHeatCapacity(UrbanType urban, int numLandunits, int mpi_rank) {
  UrbanErrorCode ierr;

  double cvRoad[numLandunits];
  double cvWall[numLandunits];
  double cvRoof[numLandunits];

  for (int i = 0; i < numLandunits; ++i) {
    cvRoad[i] = 2.0e6;
    cvWall[i] = 1.8e6;
    cvRoof[i] = 1.9e6;
  }

  UrbanCall(UrbanSetHeatCapacityRoad(urban, cvRoad, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetHeatCapacityWall(urban, cvWall, numLandunits, &ierr), &ierr);
  UrbanCall(UrbanSetHeatCapacityRoof(urban, cvRoof, numLandunits, &ierr), &ierr);

  if (mpi_rank == 0) {
    std::cout << "Set heat capacity values for all surfaces" << std::endl;
  }
}

void SetUrbanParameters(UrbanType urban, int numLandunits, int mpi_rank) {
  SetCanyonHwr(urban, numLandunits, mpi_rank);
  SetAlbedo(urban, numLandunits, mpi_rank);
  SetEmissivity(urban, numLandunits, mpi_rank);
  SetThermalConductivity(urban, numLandunits, mpi_rank);
  SetHeatCapacity(urban, numLandunits, mpi_rank);
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

    // Destroy Urban object
    UrbanCall(UrbanDestroy(&urban, &ierr), &ierr);

  }

  Kokkos::finalize();
  MPI_Finalize();
  return 0;

}
