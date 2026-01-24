! Fortran module interface for Urban library
module urban_mod
  use iso_c_binding
  implicit none

  ! Opaque handle type
  type :: UrbanType
    type(c_ptr) :: ptr = c_null_ptr
  end type UrbanType

  ! Error codes
  integer(c_int), parameter :: URBAN_SUCCESS = 0
  integer(c_int), parameter :: URBAN_ERR_INVALID_ARGUMENT = 1
  integer(c_int), parameter :: URBAN_ERR_NOT_INITIALIZED = 2
  integer(c_int), parameter :: URBAN_ERR_INTERNAL = 3
  integer(c_int), parameter :: URBAN_ERR_SIZE_MISMATCH = 4

  ! Interface declarations for C API functions
  interface
    subroutine UrbanCreate(numLandunits, urban, status) bind(C, name="UrbanCreate")
      import :: c_ptr, c_int
      integer(c_int), value :: numLandunits
      type(c_ptr) :: urban
      integer(c_int) :: status
    end subroutine UrbanCreate

    subroutine UrbanDestroy(urban, status) bind(C, name="UrbanDestroy")
      import :: c_ptr, c_int
      type(c_ptr) :: urban
      integer(c_int) :: status
    end subroutine UrbanDestroy

    subroutine UrbanSetCanyonHwr(urban, values, length, status) bind(C, name="UrbanSetCanyonHwr")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetCanyonHwr

    subroutine UrbanSetFracPervRoadOfTotalRoad(urban, values, length, status) bind(C, name="UrbanSetFracPervRoadOfTotalRoad")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetFracPervRoadOfTotalRoad

    subroutine UrbanSetWtRoof(urban, values, length, status) bind(C, name="UrbanSetWtRoof")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetWtRoof

    ! Albedo setter functions
    subroutine UrbanSetAlbedoPerviousRoad(urban, values, size, status) bind(C, name="UrbanSetAlbedoPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanSetAlbedoPerviousRoad

    subroutine UrbanSetAlbedoImperviousRoad(urban, values, size, status) bind(C, name="UrbanSetAlbedoImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanSetAlbedoImperviousRoad

    subroutine UrbanSetAlbedoSunlitWall(urban, values, size, status) bind(C, name="UrbanSetAlbedoSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanSetAlbedoSunlitWall

    subroutine UrbanSetAlbedoShadedWall(urban, values, size, status) bind(C, name="UrbanSetAlbedoShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanSetAlbedoShadedWall

    subroutine UrbanSetAlbedoRoof(urban, values, size, status) bind(C, name="UrbanSetAlbedoRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanSetAlbedoRoof

    ! Emissivity setter functions
    subroutine UrbanSetEmissivityPerviousRoad(urban, values, length, status) bind(C, name="UrbanSetEmissivityPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetEmissivityPerviousRoad

    subroutine UrbanSetEmissivityImperviousRoad(urban, values, length, status) bind(C, name="UrbanSetEmissivityImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetEmissivityImperviousRoad

    subroutine UrbanSetEmissivityWall(urban, values, length, status) bind(C, name="UrbanSetEmissivityWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetEmissivityWall

    subroutine UrbanSetEmissivityRoof(urban, values, length, status) bind(C, name="UrbanSetEmissivityRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetEmissivityRoof

    ! Thermal conductivity setter functions
    subroutine UrbanSetThermalConductivityRoad(urban, values, length, status) bind(C, name="UrbanSetThermalConductivityRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetThermalConductivityRoad

    subroutine UrbanSetThermalConductivityWall(urban, values, length, status) bind(C, name="UrbanSetThermalConductivityWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetThermalConductivityWall

    subroutine UrbanSetThermalConductivityRoof(urban, values, length, status) bind(C, name="UrbanSetThermalConductivityRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetThermalConductivityRoof

    ! Heat capacity setter functions
    subroutine UrbanSetHeatCapacityRoad(urban, values, length, status) bind(C, name="UrbanSetHeatCapacityRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetHeatCapacityRoad

    subroutine UrbanSetHeatCapacityWall(urban, values, length, status) bind(C, name="UrbanSetHeatCapacityWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetHeatCapacityWall

    subroutine UrbanSetHeatCapacityRoof(urban, values, length, status) bind(C, name="UrbanSetHeatCapacityRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetHeatCapacityRoof

    ! Height parameter setter functions
    subroutine UrbanSetForcHgtT(urban, values, length, status) bind(C, name="UrbanSetForcHgtT")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetForcHgtT

    subroutine UrbanSetForcHgtU(urban, values, length, status) bind(C, name="UrbanSetForcHgtU")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetForcHgtU

    subroutine UrbanSetZDTown(urban, values, length, status) bind(C, name="UrbanSetZDTown")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetZDTown

    subroutine UrbanSetZ0Town(urban, values, length, status) bind(C, name="UrbanSetZ0Town")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetZ0Town

    subroutine UrbanSetHtRoof(urban, values, length, status) bind(C, name="UrbanSetHtRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetHtRoof

    subroutine UrbanSetWindHgtCanyon(urban, values, length, status) bind(C, name="UrbanSetWindHgtCanyon")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetWindHgtCanyon

    ! Building parameter setter functions
    subroutine UrbanSetBuildingMaxTemperature(urban, values, length, status) bind(C, name="UrbanSetBuildingMaxTemperature")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetBuildingMaxTemperature

    subroutine UrbanSetBuildingMinTemperature(urban, values, length, status) bind(C, name="UrbanSetBuildingMinTemperature")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetBuildingMinTemperature

    ! Initialization functions
    subroutine UrbanInitializeTemperature(urban, status) bind(C, name="UrbanInitializeTemperature")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      integer(c_int) :: status
    end subroutine UrbanInitializeTemperature

    ! Time-stepping functions
    subroutine UrbanAdvance(urban, status) bind(C, name="UrbanAdvance")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      integer(c_int) :: status
    end subroutine UrbanAdvance

    ! Atmospheric forcing setter functions
    subroutine UrbanSetAtmTemp(urban, values, length, status) bind(C, name="UrbanSetAtmTemp")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmTemp

    subroutine UrbanSetAtmPotTemp(urban, values, length, status) bind(C, name="UrbanSetAtmPotTemp")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmPotTemp

    subroutine UrbanSetAtmRho(urban, values, length, status) bind(C, name="UrbanSetAtmRho")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmRho

    subroutine UrbanSetAtmSpcHumd(urban, values, length, status) bind(C, name="UrbanSetAtmSpcHumd")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmSpcHumd

    subroutine UrbanSetAtmPress(urban, values, length, status) bind(C, name="UrbanSetAtmPress")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmPress

    subroutine UrbanSetAtmWindU(urban, values, length, status) bind(C, name="UrbanSetAtmWindU")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmWindU

    subroutine UrbanSetAtmWindV(urban, values, length, status) bind(C, name="UrbanSetAtmWindV")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmWindV

    subroutine UrbanSetAtmCoszen(urban, values, length, status) bind(C, name="UrbanSetAtmCoszen")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmCoszen

    subroutine UrbanSetAtmFracSnow(urban, values, length, status) bind(C, name="UrbanSetAtmFracSnow")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmFracSnow

    subroutine UrbanSetAtmLongwaveDown(urban, values, length, status) bind(C, name="UrbanSetAtmLongwaveDown")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmLongwaveDown

    subroutine UrbanSetAtmShortwaveDown(urban, values, size, status) bind(C, name="UrbanSetAtmShortwaveDown")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanSetAtmShortwaveDown
  end interface

  contains

  ! Error handling subroutine
  subroutine UrbanError(rank, line, status)
    integer, intent(in) :: rank, line, status
    
    if (rank == 0) then
      write(*,*) 'Urban Error:'
      write(*,*) '  Line: ', line
      write(*,*) '  Status: ', status
    end if
    stop 1
  end subroutine UrbanError

end module urban_mod
