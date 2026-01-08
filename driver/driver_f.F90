program urbanxx_driver_f
  use iso_c_binding
  use mpi
  use urban_mod
  use urban_kokkos_interface
  implicit none

  integer :: ierr, mpi_rank, mpi_size

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

  ! Finalize Kokkos
  call UrbanKokkosFinalize()

  ! Finalize MPI
  call MPI_Finalize(ierr)

end program urbanxx_driver_f
