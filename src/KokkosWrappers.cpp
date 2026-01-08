#include <Kokkos_Core.hpp>
#include <iostream>

// C wrapper functions for Kokkos to be called from Fortran
extern "C" {

void UrbanKokkosInitialize() {
  int argc = 0;
  char** argv = nullptr;
  Kokkos::initialize(argc, argv);
}

void UrbanKokkosFinalize() {
  Kokkos::finalize();
}

void UrbanKokkosPrintConfiguration() {
  Kokkos::print_configuration(std::cout);
}

} // extern "C"
