#include <Kokkos_Core.hpp>
#include <iostream>
#include <mpi.h>
#include <Urban.h>

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
    UrbanErrorCode ierr = UrbanCreate(numLandunits, &urban);
    
    if (ierr != URBAN_SUCCESS) {
      if (mpi_rank == 0) {
        std::cerr << "Error creating Urban object: " << ierr << std::endl;
      }
      Kokkos::finalize();
      MPI_Finalize();
      return 1;
    }

    if (mpi_rank == 0) {
      std::cout << "Successfully created Urban object with " << numLandunits 
                << " landunits" << std::endl;
    }

    // Destroy Urban object
    ierr = UrbanDestroy(&urban);
    if (ierr != URBAN_SUCCESS) {
      if (mpi_rank == 0) {
        std::cerr << "Error destroying Urban object: " << ierr << std::endl;
      }
    }

  }

  Kokkos::finalize();
  MPI_Finalize();
  return 0;

}
