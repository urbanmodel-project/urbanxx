module urban_kokkos_interface
  use iso_c_binding
  implicit none

  ! C++ wrapper interfaces for Kokkos functions
  interface
    subroutine UrbanKokkosInitialize() bind(C, name="UrbanKokkosInitialize")
      use iso_c_binding
    end subroutine UrbanKokkosInitialize

    subroutine UrbanKokkosFinalize() bind(C, name="UrbanKokkosFinalize")
      use iso_c_binding
    end subroutine UrbanKokkosFinalize

    subroutine UrbanKokkosPrintConfiguration() bind(C, name="UrbanKokkosPrintConfiguration")
      use iso_c_binding
    end subroutine UrbanKokkosPrintConfiguration
  end interface

end module urban_kokkos_interface
