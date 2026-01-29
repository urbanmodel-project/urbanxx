! Fortran module interface for Urban library
module urban_mod
  use iso_c_binding
  implicit none

  ! Opaque handle type
  type, bind(C) :: UrbanType
    type(c_ptr) :: ptr = c_null_ptr
  end type UrbanType

  ! Error codes
  integer(c_int), parameter :: URBAN_SUCCESS = 0
  integer(c_int), parameter :: URBAN_ERR_INVALID_ARGUMENT = 1
  integer(c_int), parameter :: URBAN_ERR_NOT_INITIALIZED = 2
  integer(c_int), parameter :: URBAN_ERR_INTERNAL = 3
  integer(c_int), parameter :: URBAN_ERR_SIZE_MISMATCH = 4

  ! Low-level C interface declarations (private)
  interface
    subroutine UrbanCreate_C(numLandunits, urban, status) bind(C, name="UrbanCreate")
      import :: c_ptr, c_int
      integer(c_int), value :: numLandunits
      type(c_ptr) :: urban
      integer(c_int) :: status
    end subroutine UrbanCreate_C

    subroutine UrbanDestroy_C(urban, status) bind(C, name="UrbanDestroy")
      import :: c_ptr, c_int
      type(c_ptr) :: urban
      integer(c_int) :: status
    end subroutine UrbanDestroy_C

    subroutine UrbanSetCanyonHwr_C(urban, values, length, status) bind(C, name="UrbanSetCanyonHwr")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetCanyonHwr_C

    subroutine UrbanSetFracPervRoadOfTotalRoad_C(urban, values, length, status) bind(C, name="UrbanSetFracPervRoadOfTotalRoad")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetFracPervRoadOfTotalRoad_C

    subroutine UrbanSetWtRoof_C(urban, values, length, status) bind(C, name="UrbanSetWtRoof")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetWtRoof_C

    ! Albedo setter functions
    subroutine UrbanSetAlbedoPerviousRoad_C(urban, values, size, status) bind(C, name="UrbanSetAlbedoPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanSetAlbedoPerviousRoad_C

    subroutine UrbanSetAlbedoImperviousRoad_C(urban, values, size, status) bind(C, name="UrbanSetAlbedoImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanSetAlbedoImperviousRoad_C

    subroutine UrbanSetAlbedoSunlitWall_C(urban, values, size, status) bind(C, name="UrbanSetAlbedoSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanSetAlbedoSunlitWall_C

    subroutine UrbanSetAlbedoShadedWall_C(urban, values, size, status) bind(C, name="UrbanSetAlbedoShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanSetAlbedoShadedWall_C

    subroutine UrbanSetAlbedoRoof_C(urban, values, size, status) bind(C, name="UrbanSetAlbedoRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanSetAlbedoRoof_C

    ! Emissivity setter functions
    subroutine UrbanSetEmissivityPerviousRoad_C(urban, values, length, status) bind(C, name="UrbanSetEmissivityPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetEmissivityPerviousRoad_C

    subroutine UrbanSetEmissivityImperviousRoad_C(urban, values, length, status) bind(C, name="UrbanSetEmissivityImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetEmissivityImperviousRoad_C

    subroutine UrbanSetEmissivityWall_C(urban, values, length, status) bind(C, name="UrbanSetEmissivityWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetEmissivityWall_C

    subroutine UrbanSetEmissivityRoof_C(urban, values, length, status) bind(C, name="UrbanSetEmissivityRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetEmissivityRoof_C

    ! Thermal conductivity setter functions
    subroutine UrbanSetThermalConductivityRoad_C(urban, values, size, status) bind(C, name="UrbanSetThermalConductivityRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetThermalConductivityRoad_C

    subroutine UrbanSetThermalConductivityWall_C(urban, values, size, status) bind(C, name="UrbanSetThermalConductivityWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetThermalConductivityWall_C

    subroutine UrbanSetThermalConductivityRoof_C(urban, values, size, status) bind(C, name="UrbanSetThermalConductivityRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetThermalConductivityRoof_C

    ! Heat capacity setter functions
    subroutine UrbanSetHeatCapacityRoad_C(urban, values, size, status) bind(C, name="UrbanSetHeatCapacityRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetHeatCapacityRoad_C

    subroutine UrbanSetHeatCapacityWall_C(urban, values, size, status) bind(C, name="UrbanSetHeatCapacityWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetHeatCapacityWall_C

    subroutine UrbanSetHeatCapacityRoof_C(urban, values, size, status) bind(C, name="UrbanSetHeatCapacityRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetHeatCapacityRoof_C

    ! Temperature setter functions
    subroutine UrbanSetTemperatureRoof_C(urban, values, length, status) bind(C, name="UrbanSetTemperatureRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetTemperatureRoof_C

    subroutine UrbanSetTemperatureImperviousRoad_C(urban, values, length, status) bind(C, name="UrbanSetTemperatureImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetTemperatureImperviousRoad_C

    subroutine UrbanSetTemperaturePerviousRoad_C(urban, values, length, status) bind(C, name="UrbanSetTemperaturePerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetTemperaturePerviousRoad_C

    subroutine UrbanSetTemperatureSunlitWall_C(urban, values, length, status) bind(C, name="UrbanSetTemperatureSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetTemperatureSunlitWall_C

    subroutine UrbanSetTemperatureShadedWall_C(urban, values, length, status) bind(C, name="UrbanSetTemperatureShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetTemperatureShadedWall_C

    ! Soil property setter functions for pervious road
    subroutine UrbanSetSandPerviousRoad(urban, values, size, status) bind(C, name="UrbanSetSandPerviousRoad")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetSandPerviousRoad

    subroutine UrbanSetClayPerviousRoad(urban, values, size, status) bind(C, name="UrbanSetClayPerviousRoad")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetClayPerviousRoad

    subroutine UrbanSetOrganicPerviousRoad(urban, values, size, status) bind(C, name="UrbanSetOrganicPerviousRoad")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetOrganicPerviousRoad

    ! Height parameter setter functions
    subroutine UrbanSetForcHgtT_C(urban, values, length, status) bind(C, name="UrbanSetForcHgtT")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetForcHgtT_C

    subroutine UrbanSetForcHgtU_C(urban, values, length, status) bind(C, name="UrbanSetForcHgtU")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetForcHgtU_C

    subroutine UrbanSetZDTown_C(urban, values, length, status) bind(C, name="UrbanSetZDTown")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetZDTown_C

    subroutine UrbanSetZ0Town_C(urban, values, length, status) bind(C, name="UrbanSetZ0Town")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetZ0Town_C

    subroutine UrbanSetHtRoof_C(urban, values, length, status) bind(C, name="UrbanSetHtRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetHtRoof_C

    subroutine UrbanSetWindHgtCanyon_C(urban, values, length, status) bind(C, name="UrbanSetWindHgtCanyon")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetWindHgtCanyon_C

    ! Building parameter setter functions
    subroutine UrbanSetBuildingMaxTemperature(urban, values, length, status) bind(C, name="UrbanSetBuildingMaxTemperature")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetBuildingMaxTemperature

    subroutine UrbanSetBuildingMinTemperature(urban, values, length, status) bind(C, name="UrbanSetBuildingMinTemperature")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetBuildingMinTemperature

    subroutine UrbanSetBuildingWallThickness(urban, values, length, status) bind(C, name="UrbanSetBuildingWallThickness")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetBuildingWallThickness

    subroutine UrbanSetBuildingRoofThickness(urban, values, length, status) bind(C, name="UrbanSetBuildingRoofThickness")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetBuildingRoofThickness

    ! Setup and initialization functions
    subroutine UrbanSetup(urban, status) bind(C, name="UrbanSetup")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      integer(c_int) :: status
    end subroutine UrbanSetup

    ! Time-stepping functions
    subroutine UrbanAdvance_C(urban, status) bind(C, name="UrbanAdvance")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      integer(c_int) :: status
    end subroutine UrbanAdvance_C

    ! Physics computation functions
    subroutine UrbanComputeNetLongwave_C(urban, status) bind(C, name="UrbanComputeNetLongwave")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      integer(c_int) :: status
    end subroutine UrbanComputeNetLongwave_C

    subroutine UrbanComputeNetShortwave_C(urban, status) bind(C, name="UrbanComputeNetShortwave")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      integer(c_int) :: status
    end subroutine UrbanComputeNetShortwave_C

    subroutine UrbanComputeSurfaceFluxes_C(urban, status) bind(C, name="UrbanComputeSurfaceFluxes")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      integer(c_int) :: status
    end subroutine UrbanComputeSurfaceFluxes_C

    ! Atmospheric forcing setter functions
    subroutine UrbanSetAtmTemp_C(urban, values, length, status) bind(C, name="UrbanSetAtmTemp")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmTemp_C

    subroutine UrbanSetAtmPotTemp_C(urban, values, length, status) bind(C, name="UrbanSetAtmPotTemp")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmPotTemp_C

    subroutine UrbanSetAtmRho_C(urban, values, length, status) bind(C, name="UrbanSetAtmRho")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmRho_C

    subroutine UrbanSetAtmSpcHumd_C(urban, values, length, status) bind(C, name="UrbanSetAtmSpcHumd")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmSpcHumd_C

    subroutine UrbanSetAtmPress_C(urban, values, length, status) bind(C, name="UrbanSetAtmPress")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmPress_C

    subroutine UrbanSetAtmWindU_C(urban, values, length, status) bind(C, name="UrbanSetAtmWindU")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmWindU_C

    subroutine UrbanSetAtmWindV_C(urban, values, length, status) bind(C, name="UrbanSetAtmWindV")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmWindV_C

    subroutine UrbanSetAtmCoszen_C(urban, values, length, status) bind(C, name="UrbanSetAtmCoszen")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmCoszen_C

    subroutine UrbanSetAtmFracSnow_C(urban, values, length, status) bind(C, name="UrbanSetAtmFracSnow")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmFracSnow_C

    subroutine UrbanSetAtmLongwaveDown_C(urban, values, length, status) bind(C, name="UrbanSetAtmLongwaveDown")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAtmLongwaveDown_C

    subroutine UrbanSetAtmShortwaveDown_C(urban, values, size, status) bind(C, name="UrbanSetAtmShortwaveDown")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanSetAtmShortwaveDown_C

    ! Shortwave radiation getter functions - Absorbed
    subroutine UrbanGetAbsorbedShortwaveRoof_C(urban, values, size, status) bind(C, name="UrbanGetAbsorbedShortwaveRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanGetAbsorbedShortwaveRoof_C

    subroutine UrbanGetAbsorbedShortwaveImperviousRoad_C(urban, values, size, status) &
      bind(C, name="UrbanGetAbsorbedShortwaveImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanGetAbsorbedShortwaveImperviousRoad_C

    subroutine UrbanGetAbsorbedShortwavePerviousRoad_C(urban, values, size, status) &
      bind(C, name="UrbanGetAbsorbedShortwavePerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanGetAbsorbedShortwavePerviousRoad_C

    subroutine UrbanGetAbsorbedShortwaveSunlitWall_C(urban, values, size, status) &
      bind(C, name="UrbanGetAbsorbedShortwaveSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanGetAbsorbedShortwaveSunlitWall_C

    subroutine UrbanGetAbsorbedShortwaveShadedWall_C(urban, values, size, status) &
      bind(C, name="UrbanGetAbsorbedShortwaveShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanGetAbsorbedShortwaveShadedWall_C

    ! Shortwave radiation getter functions - Reflected
    subroutine UrbanGetReflectedShortwaveRoof_C(urban, values, size, status) &
      bind(C, name="UrbanGetReflectedShortwaveRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanGetReflectedShortwaveRoof_C

    subroutine UrbanGetReflectedShortwaveImperviousRoad_C(urban, values, size, status) &
      bind(C, name="UrbanGetReflectedShortwaveImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanGetReflectedShortwaveImperviousRoad_C

    subroutine UrbanGetReflectedShortwavePerviousRoad_C(urban, values, size, status) &
      bind(C, name="UrbanGetReflectedShortwavePerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanGetReflectedShortwavePerviousRoad_C

    subroutine UrbanGetReflectedShortwaveSunlitWall_C(urban, values, size, status) &
      bind(C, name="UrbanGetReflectedShortwaveSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanGetReflectedShortwaveSunlitWall_C

    subroutine UrbanGetReflectedShortwaveShadedWall_C(urban, values, size, status) &
      bind(C, name="UrbanGetReflectedShortwaveShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(3) :: size
      integer(c_int) :: status
    end subroutine UrbanGetReflectedShortwaveShadedWall_C

    ! Net longwave radiation getter functions
    subroutine UrbanGetNetLongwaveRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetNetLongwaveRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetNetLongwaveRoof_C

    subroutine UrbanGetNetLongwaveImperviousRoad_C(urban, values, length, &
      status) bind(C, name="UrbanGetNetLongwaveImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetNetLongwaveImperviousRoad_C

    subroutine UrbanGetNetLongwavePerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetNetLongwavePerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetNetLongwavePerviousRoad_C

    subroutine UrbanGetNetLongwaveSunlitWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetNetLongwaveSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetNetLongwaveSunlitWall_C

    subroutine UrbanGetNetLongwaveShadedWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetNetLongwaveShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetNetLongwaveShadedWall_C

    ! Upward longwave radiation getter functions
    subroutine UrbanGetUpwardLongwaveRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetUpwardLongwaveRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetUpwardLongwaveRoof_C

    subroutine UrbanGetUpwardLongwaveImperviousRoad_C(urban, values, length, &
      status) bind(C, name="UrbanGetUpwardLongwaveImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetUpwardLongwaveImperviousRoad_C

    subroutine UrbanGetUpwardLongwavePerviousRoad_C(urban, values, length, &
      status) bind(C, name="UrbanGetUpwardLongwavePerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetUpwardLongwavePerviousRoad_C

    subroutine UrbanGetUpwardLongwaveSunlitWall_C(urban, values, length, &
      status) bind(C, name="UrbanGetUpwardLongwaveSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetUpwardLongwaveSunlitWall_C

    subroutine UrbanGetUpwardLongwaveShadedWall_C(urban, values, length, &
      status) bind(C, name="UrbanGetUpwardLongwaveShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetUpwardLongwaveShadedWall_C

    ! Urban canyon air properties getter functions
    subroutine UrbanGetCanyonAirTemperature_C(urban, values, length, status) &
      bind(C, name="UrbanGetCanyonAirTemperature")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetCanyonAirTemperature_C

    subroutine UrbanGetCanyonAirHumidity_C(urban, values, length, status) &
      bind(C, name="UrbanGetCanyonAirHumidity")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetCanyonAirHumidity_C
  end interface

  contains

  ! High-level wrapper subroutines that accept UrbanType

  subroutine UrbanCreate(numLandunits, urban, status)
    integer(c_int), intent(in) :: numLandunits
    type(UrbanType), intent(inout) :: urban
    integer(c_int), intent(out) :: status
    call UrbanCreate_C(numLandunits, urban%ptr, status)
  end subroutine UrbanCreate

  subroutine UrbanDestroy(urban, status)
    type(UrbanType), intent(inout) :: urban
    integer(c_int), intent(out) :: status
    call UrbanDestroy_C(urban%ptr, status)
  end subroutine UrbanDestroy

  subroutine UrbanSetCanyonHwr(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetCanyonHwr_C(urban%ptr, values, length, status)
  end subroutine UrbanSetCanyonHwr

  subroutine UrbanSetFracPervRoadOfTotalRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetFracPervRoadOfTotalRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetFracPervRoadOfTotalRoad

  subroutine UrbanSetWtRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetWtRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanSetWtRoof

  subroutine UrbanSetAlbedoPerviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetAlbedoPerviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanSetAlbedoPerviousRoad

  subroutine UrbanSetAlbedoImperviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetAlbedoImperviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanSetAlbedoImperviousRoad

  subroutine UrbanSetAlbedoSunlitWall(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetAlbedoSunlitWall_C(urban%ptr, values, size, status)
  end subroutine UrbanSetAlbedoSunlitWall

  subroutine UrbanSetAlbedoShadedWall(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetAlbedoShadedWall_C(urban%ptr, values, size, status)
  end subroutine UrbanSetAlbedoShadedWall

  subroutine UrbanSetAlbedoRoof(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetAlbedoRoof_C(urban%ptr, values, size, status)
  end subroutine UrbanSetAlbedoRoof

  subroutine UrbanSetEmissivityPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetEmissivityPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetEmissivityPerviousRoad

  subroutine UrbanSetEmissivityImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetEmissivityImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetEmissivityImperviousRoad

  subroutine UrbanSetEmissivityWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetEmissivityWall_C(urban%ptr, values, length, status)
  end subroutine UrbanSetEmissivityWall

  subroutine UrbanSetEmissivityRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetEmissivityRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanSetEmissivityRoof

  subroutine UrbanSetThermalConductivityRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetThermalConductivityRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanSetThermalConductivityRoad

  subroutine UrbanSetThermalConductivityWall(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetThermalConductivityWall_C(urban%ptr, values, size, status)
  end subroutine UrbanSetThermalConductivityWall

  subroutine UrbanSetThermalConductivityRoof(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetThermalConductivityRoof_C(urban%ptr, values, size, status)
  end subroutine UrbanSetThermalConductivityRoof

  subroutine UrbanSetHeatCapacityRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetHeatCapacityRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanSetHeatCapacityRoad

  subroutine UrbanSetHeatCapacityWall(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetHeatCapacityWall_C(urban%ptr, values, size, status)
  end subroutine UrbanSetHeatCapacityWall

  subroutine UrbanSetHeatCapacityRoof(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetHeatCapacityRoof_C(urban%ptr, values, size, status)
  end subroutine UrbanSetHeatCapacityRoof

  subroutine UrbanSetTemperatureRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetTemperatureRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanSetTemperatureRoof

  subroutine UrbanSetTemperatureImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetTemperatureImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetTemperatureImperviousRoad

  subroutine UrbanSetTemperaturePerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetTemperaturePerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetTemperaturePerviousRoad

  subroutine UrbanSetTemperatureSunlitWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetTemperatureSunlitWall_C(urban%ptr, values, length, status)
  end subroutine UrbanSetTemperatureSunlitWall

  subroutine UrbanSetTemperatureShadedWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetTemperatureShadedWall_C(urban%ptr, values, length, status)
  end subroutine UrbanSetTemperatureShadedWall

  subroutine UrbanSetForcHgtT(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetForcHgtT_C(urban%ptr, values, length, status)
  end subroutine UrbanSetForcHgtT

  subroutine UrbanSetForcHgtU(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetForcHgtU_C(urban%ptr, values, length, status)
  end subroutine UrbanSetForcHgtU

  subroutine UrbanSetZDTown(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetZDTown_C(urban%ptr, values, length, status)
  end subroutine UrbanSetZDTown

  subroutine UrbanSetZ0Town(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetZ0Town_C(urban%ptr, values, length, status)
  end subroutine UrbanSetZ0Town

  subroutine UrbanSetHtRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetHtRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanSetHtRoof

  subroutine UrbanSetWindHgtCanyon(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetWindHgtCanyon_C(urban%ptr, values, length, status)
  end subroutine UrbanSetWindHgtCanyon

  subroutine UrbanAdvance(urban, status)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(out) :: status
    call UrbanAdvance_C(urban%ptr, status)
  end subroutine UrbanAdvance

  subroutine UrbanComputeNetLongwave(urban, status)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(out) :: status
    call UrbanComputeNetLongwave_C(urban%ptr, status)
  end subroutine UrbanComputeNetLongwave

  subroutine UrbanComputeNetShortwave(urban, status)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(out) :: status
    call UrbanComputeNetShortwave_C(urban%ptr, status)
  end subroutine UrbanComputeNetShortwave

  subroutine UrbanComputeSurfaceFluxes(urban, status)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(out) :: status
    call UrbanComputeSurfaceFluxes_C(urban%ptr, status)
  end subroutine UrbanComputeSurfaceFluxes

  subroutine UrbanSetAtmTemp(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAtmTemp_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAtmTemp

  subroutine UrbanSetAtmPotTemp(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAtmPotTemp_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAtmPotTemp

  subroutine UrbanSetAtmRho(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAtmRho_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAtmRho

  subroutine UrbanSetAtmSpcHumd(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAtmSpcHumd_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAtmSpcHumd

  subroutine UrbanSetAtmPress(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAtmPress_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAtmPress

  subroutine UrbanSetAtmWindU(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAtmWindU_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAtmWindU

  subroutine UrbanSetAtmWindV(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAtmWindV_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAtmWindV

  subroutine UrbanSetAtmCoszen(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAtmCoszen_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAtmCoszen

  subroutine UrbanSetAtmFracSnow(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAtmFracSnow_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAtmFracSnow

  subroutine UrbanSetAtmLongwaveDown(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAtmLongwaveDown_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAtmLongwaveDown

  subroutine UrbanSetAtmShortwaveDown(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetAtmShortwaveDown_C(urban%ptr, values, size, status)
  end subroutine UrbanSetAtmShortwaveDown

  ! Shortwave radiation getter functions - Absorbed
  subroutine UrbanGetAbsorbedShortwaveRoof(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetAbsorbedShortwaveRoof_C(urban%ptr, values, size, status)
  end subroutine UrbanGetAbsorbedShortwaveRoof

  subroutine UrbanGetAbsorbedShortwaveImperviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetAbsorbedShortwaveImperviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanGetAbsorbedShortwaveImperviousRoad

  subroutine UrbanGetAbsorbedShortwavePerviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetAbsorbedShortwavePerviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanGetAbsorbedShortwavePerviousRoad

  subroutine UrbanGetAbsorbedShortwaveSunlitWall(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetAbsorbedShortwaveSunlitWall_C(urban%ptr, values, size, status)
  end subroutine UrbanGetAbsorbedShortwaveSunlitWall

  subroutine UrbanGetAbsorbedShortwaveShadedWall(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetAbsorbedShortwaveShadedWall_C(urban%ptr, values, size, status)
  end subroutine UrbanGetAbsorbedShortwaveShadedWall

  ! Shortwave radiation getter functions - Reflected
  subroutine UrbanGetReflectedShortwaveRoof(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetReflectedShortwaveRoof_C(urban%ptr, values, size, status)
  end subroutine UrbanGetReflectedShortwaveRoof

  subroutine UrbanGetReflectedShortwaveImperviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetReflectedShortwaveImperviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanGetReflectedShortwaveImperviousRoad

  subroutine UrbanGetReflectedShortwavePerviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetReflectedShortwavePerviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanGetReflectedShortwavePerviousRoad

  subroutine UrbanGetReflectedShortwaveSunlitWall(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetReflectedShortwaveSunlitWall_C(urban%ptr, values, size, status)
  end subroutine UrbanGetReflectedShortwaveSunlitWall

  subroutine UrbanGetReflectedShortwaveShadedWall(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(3), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetReflectedShortwaveShadedWall_C(urban%ptr, values, size, status)
  end subroutine UrbanGetReflectedShortwaveShadedWall

  ! Net longwave radiation getter subroutines
  subroutine UrbanGetNetLongwaveRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetNetLongwaveRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetNetLongwaveRoof

  subroutine UrbanGetNetLongwaveImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetNetLongwaveImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetNetLongwaveImperviousRoad

  subroutine UrbanGetNetLongwavePerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetNetLongwavePerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetNetLongwavePerviousRoad

  subroutine UrbanGetNetLongwaveSunlitWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetNetLongwaveSunlitWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetNetLongwaveSunlitWall

  subroutine UrbanGetNetLongwaveShadedWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetNetLongwaveShadedWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetNetLongwaveShadedWall

  ! Upward longwave radiation getter subroutines
  subroutine UrbanGetUpwardLongwaveRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetUpwardLongwaveRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetUpwardLongwaveRoof

  subroutine UrbanGetUpwardLongwaveImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetUpwardLongwaveImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetUpwardLongwaveImperviousRoad

  subroutine UrbanGetUpwardLongwavePerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetUpwardLongwavePerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetUpwardLongwavePerviousRoad

  subroutine UrbanGetUpwardLongwaveSunlitWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetUpwardLongwaveSunlitWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetUpwardLongwaveSunlitWall

  subroutine UrbanGetUpwardLongwaveShadedWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetUpwardLongwaveShadedWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetUpwardLongwaveShadedWall

  ! Urban canyon air properties getter subroutines
  subroutine UrbanGetCanyonAirTemperature(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetCanyonAirTemperature_C(urban%ptr, values, length, status)
  end subroutine UrbanGetCanyonAirTemperature

  subroutine UrbanGetCanyonAirHumidity(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetCanyonAirHumidity_C(urban%ptr, values, length, status)
  end subroutine UrbanGetCanyonAirHumidity

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
