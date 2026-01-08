program urbanxx_driver_f
#include "finclude/UrbanMacros.h"
  use iso_c_binding
  use mpi
  use urban_mod
  use urban_kokkos_interface
  implicit none

  integer :: ierr, mpi_rank, mpi_size
  type(UrbanType) :: urban
  integer(c_int) :: numLandunits, status

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr)

  ! Initialize Kokkos (via C++ wrapper)
  call UrbanKokkosInitialize()

  if (mpi_rank == 0) then
    write(*,*) '=== Fortran Driver with Kokkos ==='
    write(*,*) 'MPI Size: ', mpi_size
    call UrbanKokkosPrintConfiguration()
  end if

  ! Create Urban object
  numLandunits = 10
  CallA(UrbanCreate(numLandunits, urban%ptr, status))

  if (mpi_rank == 0) then
    write(*,*) 'Successfully created Urban object with ', numLandunits, ' landunits'
  end if

  ! Destroy Urban object
  CallA(UrbanDestroy(urban%ptr, status))

  ! Finalize Kokkos
  call UrbanKokkosFinalize()

  ! Finalize MPI
  call MPI_Finalize(ierr)

end program urbanxx_driver_f
