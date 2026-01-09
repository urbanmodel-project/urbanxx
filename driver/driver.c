#include <Kokkos_Core.hpp>
#include <iostream>
#include <mpi.h>
#include <Urban.h>
#include <UrbanMacros.h>

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

    // Set CanyonHwr values
    double canyonHwr[numLandunits];
    for (int i = 0; i < numLandunits; ++i) {
      canyonHwr[i] = 4.80000019073486;
    }
    UrbanCall(UrbanSetCanyonHwr(urban, canyonHwr, numLandunits, &ierr), &ierr);

    // Destroy Urban object
    UrbanCall(UrbanDestroy(&urban, &ierr), &ierr);

  }

  Kokkos::finalize();
  MPI_Finalize();
  return 0;

}
