#include <Kokkos_Core.hpp>
#include <iostream>
#include <mpi.h>
#include "private/DataTypesImpl.h"

using namespace URBANXX;

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

  }

  Kokkos::finalize();
  MPI_Finalize();
  return 0;

}
