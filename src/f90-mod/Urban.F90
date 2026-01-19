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

    function UrbanKokkosIsLayoutRight() bind(C, name="UrbanKokkosIsLayoutRight")
      use iso_c_binding
      logical(c_bool) :: UrbanKokkosIsLayoutRight
    end function UrbanKokkosIsLayoutRight

    function UrbanKokkosIsLayoutLeft() bind(C, name="UrbanKokkosIsLayoutLeft")
      use iso_c_binding
      logical(c_bool) :: UrbanKokkosIsLayoutLeft
    end function UrbanKokkosIsLayoutLeft
  end interface

end module urban_kokkos_interface
