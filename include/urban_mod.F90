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
        subroutine UrbanSetAtmRain_C(urban, values, length, status) bind(C, name="UrbanSetAtmRain")
          import :: c_ptr, c_int
          type(c_ptr), value :: urban
          type(c_ptr), value :: values
          integer(c_int), value :: length
          integer(c_int) :: status
        end subroutine UrbanSetAtmRain_C

        subroutine UrbanSetAtmSnow_C(urban, values, length, status) bind(C, name="UrbanSetAtmSnow")
          import :: c_ptr, c_int
          type(c_ptr), value :: urban
          type(c_ptr), value :: values
          integer(c_int), value :: length
          integer(c_int) :: status
        end subroutine UrbanSetAtmSnow_C
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

    ! Number of active layers setter function
    subroutine UrbanSetNumberOfActiveLayersImperviousRoad(urban, values, &
      length, status) bind(C, name="UrbanSetNumberOfActiveLayersImperviousRoad")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetNumberOfActiveLayersImperviousRoad

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

    subroutine UrbanSetBuildingTemperature(urban, values, length, status) bind(C, name="UrbanSetBuildingTemperature")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetBuildingTemperature

    subroutine UrbanSetBuildingWallThickness(urban, values, length, status) bind(C, name="UrbanSetBuildingWallThickness")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetBuildingWallThickness

    subroutine UrbanSetBuildingRoofThickness(urban, values, length, status) &
      bind(C, name="UrbanSetBuildingRoofThickness")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetBuildingRoofThickness

    ! Layer temperature setter functions (2D)
    subroutine UrbanSetLayerTempRoof(urban, values, size, status) &
      bind(C, name="UrbanSetLayerTempRoof")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetLayerTempRoof

    subroutine UrbanSetLayerTempImperviousRoad(urban, values, size, status) &
      bind(C, name="UrbanSetLayerTempImperviousRoad")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetLayerTempImperviousRoad

    subroutine UrbanSetLayerTempPerviousRoad(urban, values, size, status) &
      bind(C, name="UrbanSetLayerTempPerviousRoad")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetLayerTempPerviousRoad

    subroutine UrbanSetLayerTempSunlitWall(urban, values, size, status) &
      bind(C, name="UrbanSetLayerTempSunlitWall")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetLayerTempSunlitWall

    subroutine UrbanSetLayerTempShadedWall(urban, values, size, status) &
      bind(C, name="UrbanSetLayerTempShadedWall")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetLayerTempShadedWall

    ! Canyon air property setter functions (1D)
    subroutine UrbanSetCanyonAirTemperature(urban, values, length, status) &
      bind(C, name="UrbanSetCanyonAirTemperature")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetCanyonAirTemperature

    subroutine UrbanSetCanyonSpecificHumidity(urban, values, length, status) &
      bind(C, name="UrbanSetCanyonSpecificHumidity")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetCanyonSpecificHumidity

    ! Top-layer soil water setter functions (for soil flux partitioning)
    subroutine UrbanSetTopH2OSoiLiqRoof(urban, values, length, &
      status) bind(C, name="UrbanSetTopH2OSoiLiqRoof")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetTopH2OSoiLiqRoof

    subroutine UrbanSetTopH2OSoiIceRoof(urban, values, length, &
      status) bind(C, name="UrbanSetTopH2OSoiIceRoof")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetTopH2OSoiIceRoof

    subroutine UrbanSetTopH2OSoiLiqImperviousRoad(urban, values, length, &
      status) bind(C, name="UrbanSetTopH2OSoiLiqImperviousRoad")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetTopH2OSoiLiqImperviousRoad

    subroutine UrbanSetTopH2OSoiIceImperviousRoad(urban, values, length, &
      status) bind(C, name="UrbanSetTopH2OSoiIceImperviousRoad")
      import :: c_ptr, c_int, UrbanType
      type(UrbanType), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetTopH2OSoiIceImperviousRoad

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

    subroutine UrbanComputeHydrology_C(urban, dtime, status) bind(C, name="UrbanComputeHydrology")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      real(c_double), value :: dtime
      integer(c_int) :: status
    end subroutine UrbanComputeHydrology_C

    subroutine UrbanComputeHeatDiffusion_C(urban, status) bind(C, name="UrbanComputeHeatDiffusion")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      integer(c_int) :: status
    end subroutine UrbanComputeHeatDiffusion_C

    subroutine UrbanComputeNetShortwaveRadiation_C(urban, status) bind(C, name="UrbanComputeNetShortwaveRadiation")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      integer(c_int) :: status
    end subroutine UrbanComputeNetShortwaveRadiation_C

    subroutine UrbanComputeSoilFluxes_C(urban, status) bind(C, name="UrbanComputeSoilFluxes")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      integer(c_int) :: status
    end subroutine UrbanComputeSoilFluxes_C

    ! Soil ground heat flux (EflxSoilGrnd) getter C interfaces
    subroutine UrbanGetEflxSoilGrndRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetEflxSoilGrndRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetEflxSoilGrndRoof_C

    subroutine UrbanGetEflxSoilGrndImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetEflxSoilGrndImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetEflxSoilGrndImperviousRoad_C

    subroutine UrbanGetEflxSoilGrndPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetEflxSoilGrndPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetEflxSoilGrndPerviousRoad_C

    subroutine UrbanGetEflxSoilGrndSunlitWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetEflxSoilGrndSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetEflxSoilGrndSunlitWall_C

    subroutine UrbanGetEflxSoilGrndShadedWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetEflxSoilGrndShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetEflxSoilGrndShadedWall_C

    ! QflxEvapGrnd getter C interfaces
    subroutine UrbanGetQflxEvapGrndRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxEvapGrndRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxEvapGrndRoof_C

    subroutine UrbanGetQflxEvapGrndImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxEvapGrndImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxEvapGrndImperviousRoad_C

    subroutine UrbanGetQflxEvapGrndPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxEvapGrndPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxEvapGrndPerviousRoad_C

    ! QflxSubSnow getter C interfaces
    subroutine UrbanGetQflxSubSnowRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxSubSnowRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxSubSnowRoof_C

    subroutine UrbanGetQflxSubSnowImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxSubSnowImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxSubSnowImperviousRoad_C

    subroutine UrbanGetQflxSubSnowPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxSubSnowPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxSubSnowPerviousRoad_C

    ! QflxDewSnow getter C interfaces
    subroutine UrbanGetQflxDewSnowRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxDewSnowRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxDewSnowRoof_C

    subroutine UrbanGetQflxDewSnowImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxDewSnowImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxDewSnowImperviousRoad_C

    subroutine UrbanGetQflxDewSnowPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxDewSnowPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxDewSnowPerviousRoad_C

    ! QflxDewGrnd getter C interfaces
    subroutine UrbanGetQflxDewGrndRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxDewGrndRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxDewGrndRoof_C

    subroutine UrbanGetQflxDewGrndImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxDewGrndImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxDewGrndImperviousRoad_C

    subroutine UrbanGetQflxDewGrndPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxDewGrndPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxDewGrndPerviousRoad_C

    ! Surface runoff setter functions (pervious road inputs)
    subroutine UrbanSetWtfactPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetWtfactPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetWtfactPerviousRoad_C

    subroutine UrbanSetFoverPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetFoverPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetFoverPerviousRoad_C

    subroutine UrbanSetFrostTablePerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetFrostTablePerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetFrostTablePerviousRoad_C

    subroutine UrbanSetZwtPerchedPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetZwtPerchedPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetZwtPerchedPerviousRoad_C

    ! Surface runoff compute function
    subroutine UrbanComputeSurfaceRunoff_C(urban, dtime, status) &
      bind(C, name="UrbanComputeSurfaceRunoff")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      real(c_double), value :: dtime
      integer(c_int) :: status
    end subroutine UrbanComputeSurfaceRunoff_C

    ! Surface runoff getter functions (QflxSurf — all five surfaces)
    subroutine UrbanGetQflxSurfRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxSurfRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxSurfRoof_C

    subroutine UrbanGetQflxSurfImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxSurfImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxSurfImperviousRoad_C

    subroutine UrbanGetQflxSurfPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxSurfPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxSurfPerviousRoad_C

    subroutine UrbanGetQflxSurfSunlitWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxSurfSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxSurfSunlitWall_C

    subroutine UrbanGetQflxSurfShadedWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetQflxSurfShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetQflxSurfShadedWall_C

    subroutine UrbanGetTopH2OSoiLiqRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetTopH2OSoiLiqRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetTopH2OSoiLiqRoof_C

    subroutine UrbanGetTopH2OSoiLiqImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetTopH2OSoiLiqImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetTopH2OSoiLiqImperviousRoad_C

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

    ! Hydrology boundary condition setter functions
    subroutine UrbanSetInfiltrationFluxForPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetInfiltrationFluxForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetInfiltrationFluxForPerviousRoad_C

    subroutine UrbanSetTranspirationFluxForPerviousRoad_C(urban, values, size, status) &
      bind(C, name="UrbanSetTranspirationFluxForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetTranspirationFluxForPerviousRoad_C

    subroutine UrbanSetWaterTableDepth_C(urban, values, length, status) bind(C, name="UrbanSetWaterTableDepth")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetWaterTableDepth_C

    ! Soil water content setter functions
    subroutine UrbanSetSoilLiquidWaterForPerviousRoad_C(urban, values, size, status) &
      bind(C, name="UrbanSetSoilLiquidWaterForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetSoilLiquidWaterForPerviousRoad_C

    subroutine UrbanSetSoilIceContentForPerviousRoad_C(urban, values, size, status) &
      bind(C, name="UrbanSetSoilIceContentForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanSetSoilIceContentForPerviousRoad_C

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

    ! Net shortwave radiation getter functions
    subroutine UrbanGetNetShortwaveRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetNetShortwaveRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetNetShortwaveRoof_C

    subroutine UrbanGetNetShortwaveImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetNetShortwaveImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetNetShortwaveImperviousRoad_C

    subroutine UrbanGetNetShortwavePerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetNetShortwavePerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetNetShortwavePerviousRoad_C

    subroutine UrbanGetNetShortwaveSunlitWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetNetShortwaveSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetNetShortwaveSunlitWall_C

    subroutine UrbanGetNetShortwaveShadedWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetNetShortwaveShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetNetShortwaveShadedWall_C

    ! Sensible heat flux getter functions
    subroutine UrbanGetSensibleHeatFluxRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetSensibleHeatFluxRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetSensibleHeatFluxRoof_C

    subroutine UrbanGetSensibleHeatFluxImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetSensibleHeatFluxImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetSensibleHeatFluxImperviousRoad_C

    subroutine UrbanGetSensibleHeatFluxPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetSensibleHeatFluxPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetSensibleHeatFluxPerviousRoad_C

    subroutine UrbanGetSensibleHeatFluxSunlitWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetSensibleHeatFluxSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetSensibleHeatFluxSunlitWall_C

    subroutine UrbanGetSensibleHeatFluxShadedWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetSensibleHeatFluxShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetSensibleHeatFluxShadedWall_C

    ! Soil evaporation flux getter functions
    subroutine UrbanGetEvapFluxRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetEvapFluxRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetEvapFluxRoof_C

    subroutine UrbanGetEvapFluxImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetEvapFluxImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetEvapFluxImperviousRoad_C

    subroutine UrbanGetEvapFluxPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetEvapFluxPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetEvapFluxPerviousRoad_C

    ! Cgrnds (d(sensible heat flux)/dT) getter functions
    subroutine UrbanGetCgrndsRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetCgrndsRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetCgrndsRoof_C

    subroutine UrbanGetCgrndsImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetCgrndsImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetCgrndsImperviousRoad_C

    subroutine UrbanGetCgrndsPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetCgrndsPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetCgrndsPerviousRoad_C

    subroutine UrbanGetCgrndsSunlitWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetCgrndsSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetCgrndsSunlitWall_C

    subroutine UrbanGetCgrndsShadedWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetCgrndsShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetCgrndsShadedWall_C

    ! Cgrndl (d(latent heat flux)/dT) getter functions
    subroutine UrbanGetCgrndlRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetCgrndlRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetCgrndlRoof_C

    subroutine UrbanGetCgrndlImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetCgrndlImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetCgrndlImperviousRoad_C

    subroutine UrbanGetCgrndlPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetCgrndlPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetCgrndlPerviousRoad_C

    subroutine UrbanGetCgrndlSunlitWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetCgrndlSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetCgrndlSunlitWall_C

    subroutine UrbanGetCgrndlShadedWall_C(urban, values, length, status) &
      bind(C, name="UrbanGetCgrndlShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetCgrndlShadedWall_C

    ! Layer temperature getter functions (2D)
    subroutine UrbanGetLayerTempRoof_C(urban, values, size, status) &
      bind(C, name="UrbanGetLayerTempRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanGetLayerTempRoof_C

    subroutine UrbanGetLayerTempImperviousRoad_C(urban, values, size, status) &
      bind(C, name="UrbanGetLayerTempImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanGetLayerTempImperviousRoad_C

    subroutine UrbanGetLayerTempPerviousRoad_C(urban, values, size, status) &
      bind(C, name="UrbanGetLayerTempPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanGetLayerTempPerviousRoad_C

    subroutine UrbanGetLayerTempSunlitWall_C(urban, values, size, status) &
      bind(C, name="UrbanGetLayerTempSunlitWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanGetLayerTempSunlitWall_C

    subroutine UrbanGetLayerTempShadedWall_C(urban, values, size, status) &
      bind(C, name="UrbanGetLayerTempShadedWall")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanGetLayerTempShadedWall_C

    ! Hydrology getter functions (pervious road)
    subroutine UrbanGetSoilLiquidWaterPerviousRoad_C(urban, values, size, status) &
      bind(C, name="UrbanGetSoilLiquidWaterPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanGetSoilLiquidWaterPerviousRoad_C

    subroutine UrbanGetSoilVolumetricWaterPerviousRoad_C(urban, values, size, status) &
      bind(C, name="UrbanGetSoilVolumetricWaterPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanGetSoilVolumetricWaterPerviousRoad_C

    subroutine UrbanGetAquiferRechargeRatePerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetAquiferRechargeRatePerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetAquiferRechargeRatePerviousRoad_C

    subroutine UrbanGetWaterDeficitFluxPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetWaterDeficitFluxPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetWaterDeficitFluxPerviousRoad_C

    ! Infiltration input setter C interfaces
    subroutine UrbanSetSurfaceRunoffForPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetSurfaceRunoffForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetSurfaceRunoffForPerviousRoad_C

    subroutine UrbanSetGroundEvapFluxForPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetGroundEvapFluxForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetGroundEvapFluxForPerviousRoad_C

    ! Infiltration compute C interface
    subroutine UrbanComputeInfiltration_C(urban, status) &
      bind(C, name="UrbanComputeInfiltration")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      integer(c_int) :: status
    end subroutine UrbanComputeInfiltration_C

    ! Infiltration getter C interface
    subroutine UrbanGetInfiltrationFluxPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetInfiltrationFluxPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetInfiltrationFluxPerviousRoad_C

    ! WaterTable setter C interfaces
    subroutine UrbanSetAquiferWaterForPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetAquiferWaterForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetAquiferWaterForPerviousRoad_C

    subroutine UrbanSetFracH2osfcForPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetFracH2osfcForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetFracH2osfcForPerviousRoad_C

    subroutine UrbanSetDewGrndFluxForPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetDewGrndFluxForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetDewGrndFluxForPerviousRoad_C

    subroutine UrbanSetDewSnowFluxForPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetDewSnowFluxForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetDewSnowFluxForPerviousRoad_C

    subroutine UrbanSetSubSnowFluxForPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetSubSnowFluxForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetSubSnowFluxForPerviousRoad_C

    subroutine UrbanSetQchargeForPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetQchargeForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetQchargeForPerviousRoad_C

    ! WaterTable compute C interface
    subroutine UrbanComputeWaterTable_C(urban, dtime, status) &
      bind(C, name="UrbanComputeWaterTable")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      real(c_double), value :: dtime
      integer(c_int) :: status
    end subroutine UrbanComputeWaterTable_C

    ! WaterTable getter C interfaces (1D)
    subroutine UrbanGetWaterTableDepthPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetWaterTableDepthPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetWaterTableDepthPerviousRoad_C

    subroutine UrbanGetAquiferWaterPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetAquiferWaterPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetAquiferWaterPerviousRoad_C

    subroutine UrbanGetZwtPerchedPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetZwtPerchedPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetZwtPerchedPerviousRoad_C

    subroutine UrbanGetSubSnowFluxPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetSubSnowFluxPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetSubSnowFluxPerviousRoad_C

    subroutine UrbanGetDrainFluxPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetDrainFluxPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetDrainFluxPerviousRoad_C

    subroutine UrbanGetRsubSatPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetRsubSatPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetRsubSatPerviousRoad_C

    ! WaterTable getter C interfaces (2D)
    subroutine UrbanGetSoilIceContentPerviousRoad_C(urban, values, size, status) &
      bind(C, name="UrbanGetSoilIceContentPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), dimension(2) :: size
      integer(c_int) :: status
    end subroutine UrbanGetSoilIceContentPerviousRoad_C

    ! DewCondensation setter C interfaces (roof)
    subroutine UrbanSetQflxDewGrndRoof_C(urban, values, length, status) &
      bind(C, name="UrbanSetQflxDewGrndRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetQflxDewGrndRoof_C

    subroutine UrbanSetQflxDewSnowRoof_C(urban, values, length, status) &
      bind(C, name="UrbanSetQflxDewSnowRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetQflxDewSnowRoof_C

    subroutine UrbanSetQflxSubSnowRoof_C(urban, values, length, status) &
      bind(C, name="UrbanSetQflxSubSnowRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetQflxSubSnowRoof_C

    ! DewCondensation setter C interfaces (impervious road)
    subroutine UrbanSetQflxDewGrndImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetQflxDewGrndImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetQflxDewGrndImperviousRoad_C

    subroutine UrbanSetQflxDewSnowImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetQflxDewSnowImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetQflxDewSnowImperviousRoad_C

    subroutine UrbanSetQflxSubSnowImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetQflxSubSnowImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetQflxSubSnowImperviousRoad_C

    ! DewCondensation compute C interface
    subroutine UrbanComputeDewCondensationRoofImperviousRoad_C(urban, dtime, status) &
      bind(C, name="UrbanComputeDewCondensationRoofImperviousRoad")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      real(c_double), value :: dtime
      integer(c_int) :: status
    end subroutine UrbanComputeDewCondensationRoofImperviousRoad_C

    ! DewCondensation getter C interfaces (TopH2OSoiIce for roof and impervious road)
    subroutine UrbanGetTopH2OSoiIceRoof_C(urban, values, length, status) &
      bind(C, name="UrbanGetTopH2OSoiIceRoof")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetTopH2OSoiIceRoof_C

    subroutine UrbanGetTopH2OSoiIceImperviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanGetTopH2OSoiIceImperviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanGetTopH2OSoiIceImperviousRoad_C

    ! Drainage init-time constant setter C interfaces
    subroutine UrbanSetRsubTopGlobalMax_C(urban, value, status) &
      bind(C, name="UrbanSetRsubTopGlobalMax")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      real(c_double), value :: value
      integer(c_int) :: status
    end subroutine UrbanSetRsubTopGlobalMax_C

    subroutine UrbanSetPondmax_C(urban, value, status) &
      bind(C, name="UrbanSetPondmax")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      real(c_double), value :: value
      integer(c_int) :: status
    end subroutine UrbanSetPondmax_C

    subroutine UrbanSetWatmin_C(urban, value, status) &
      bind(C, name="UrbanSetWatmin")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      real(c_double), value :: value
      integer(c_int) :: status
    end subroutine UrbanSetWatmin_C

    subroutine UrbanSetEice_C(urban, value, status) &
      bind(C, name="UrbanSetEice")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      real(c_double), value :: value
      integer(c_int) :: status
    end subroutine UrbanSetEice_C

    ! Drainage timestep setter C interfaces
    subroutine UrbanSetHkDepthForPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetHkDepthForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetHkDepthForPerviousRoad_C

    subroutine UrbanSetTopoSlopeForPerviousRoad_C(urban, values, length, status) &
      bind(C, name="UrbanSetTopoSlopeForPerviousRoad")
      import :: c_ptr, c_int
      type(c_ptr), value :: urban
      type(c_ptr), value :: values
      integer(c_int), value :: length
      integer(c_int) :: status
    end subroutine UrbanSetTopoSlopeForPerviousRoad_C

    ! Drainage compute C interface
    subroutine UrbanComputeDrainage_C(urban, dtime, status) &
      bind(C, name="UrbanComputeDrainage")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: urban
      real(c_double), value :: dtime
      integer(c_int) :: status
    end subroutine UrbanComputeDrainage_C

  end interface

  contains

  subroutine UrbanSetAtmRain(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAtmRain_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAtmRain

  subroutine UrbanSetAtmSnow(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAtmSnow_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAtmSnow

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

  subroutine UrbanComputeHydrology(urban, dtime, status)
    type(UrbanType), intent(in) :: urban
    real(c_double), value :: dtime
    integer(c_int), intent(out) :: status
    call UrbanComputeHydrology_C(urban%ptr, dtime, status)
  end subroutine UrbanComputeHydrology

  subroutine UrbanComputeHeatDiffusion(urban, status)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(out) :: status
    call UrbanComputeHeatDiffusion_C(urban%ptr, status)
  end subroutine UrbanComputeHeatDiffusion

  subroutine UrbanComputeNetShortwaveRadiation(urban, status)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(out) :: status
    call UrbanComputeNetShortwaveRadiation_C(urban%ptr, status)
  end subroutine UrbanComputeNetShortwaveRadiation

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

  ! Hydrology boundary condition setter functions
  subroutine UrbanSetInfiltrationFluxForPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetInfiltrationFluxForPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetInfiltrationFluxForPerviousRoad

  subroutine UrbanSetTranspirationFluxForPerviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetTranspirationFluxForPerviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanSetTranspirationFluxForPerviousRoad

  subroutine UrbanSetWaterTableDepth(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetWaterTableDepth_C(urban%ptr, values, length, status)
  end subroutine UrbanSetWaterTableDepth

  ! Soil water content setter functions
  subroutine UrbanSetSoilLiquidWaterForPerviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetSoilLiquidWaterForPerviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanSetSoilLiquidWaterForPerviousRoad

  subroutine UrbanSetSoilIceContentForPerviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanSetSoilIceContentForPerviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanSetSoilIceContentForPerviousRoad

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

  ! Net shortwave radiation getter subroutines
  subroutine UrbanGetNetShortwaveRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetNetShortwaveRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetNetShortwaveRoof

  subroutine UrbanGetNetShortwaveImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetNetShortwaveImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetNetShortwaveImperviousRoad

  subroutine UrbanGetNetShortwavePerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetNetShortwavePerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetNetShortwavePerviousRoad

  subroutine UrbanGetNetShortwaveSunlitWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetNetShortwaveSunlitWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetNetShortwaveSunlitWall

  subroutine UrbanGetNetShortwaveShadedWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetNetShortwaveShadedWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetNetShortwaveShadedWall

  ! Sensible heat flux getter subroutines
  subroutine UrbanGetSensibleHeatFluxRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetSensibleHeatFluxRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetSensibleHeatFluxRoof

  subroutine UrbanGetSensibleHeatFluxImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetSensibleHeatFluxImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetSensibleHeatFluxImperviousRoad

  subroutine UrbanGetSensibleHeatFluxPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetSensibleHeatFluxPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetSensibleHeatFluxPerviousRoad

  subroutine UrbanGetSensibleHeatFluxSunlitWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetSensibleHeatFluxSunlitWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetSensibleHeatFluxSunlitWall

  subroutine UrbanGetSensibleHeatFluxShadedWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetSensibleHeatFluxShadedWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetSensibleHeatFluxShadedWall

  ! Soil evaporation flux getter subroutines
  subroutine UrbanGetEvapFluxRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetEvapFluxRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetEvapFluxRoof

  subroutine UrbanGetEvapFluxImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetEvapFluxImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetEvapFluxImperviousRoad

  subroutine UrbanGetEvapFluxPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetEvapFluxPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetEvapFluxPerviousRoad

  ! Cgrnds (d(sensible heat flux)/dT) getter subroutines
  subroutine UrbanGetCgrndsRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetCgrndsRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetCgrndsRoof

  subroutine UrbanGetCgrndsImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetCgrndsImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetCgrndsImperviousRoad

  subroutine UrbanGetCgrndsPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetCgrndsPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetCgrndsPerviousRoad

  subroutine UrbanGetCgrndsSunlitWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetCgrndsSunlitWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetCgrndsSunlitWall

  subroutine UrbanGetCgrndsShadedWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetCgrndsShadedWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetCgrndsShadedWall

  ! Cgrndl (d(latent heat flux)/dT) getter subroutines
  subroutine UrbanGetCgrndlRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetCgrndlRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetCgrndlRoof

  subroutine UrbanGetCgrndlImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetCgrndlImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetCgrndlImperviousRoad

  subroutine UrbanGetCgrndlPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetCgrndlPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetCgrndlPerviousRoad

  subroutine UrbanGetCgrndlSunlitWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetCgrndlSunlitWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetCgrndlSunlitWall

  subroutine UrbanGetCgrndlShadedWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetCgrndlShadedWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetCgrndlShadedWall

  ! Layer temperature getter subroutines (2D)
  subroutine UrbanGetLayerTempRoof(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetLayerTempRoof_C(urban%ptr, values, size, status)
  end subroutine UrbanGetLayerTempRoof

  subroutine UrbanGetLayerTempImperviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetLayerTempImperviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanGetLayerTempImperviousRoad

  subroutine UrbanGetLayerTempPerviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetLayerTempPerviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanGetLayerTempPerviousRoad

  subroutine UrbanGetLayerTempSunlitWall(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetLayerTempSunlitWall_C(urban%ptr, values, size, status)
  end subroutine UrbanGetLayerTempSunlitWall

  subroutine UrbanGetLayerTempShadedWall(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetLayerTempShadedWall_C(urban%ptr, values, size, status)
  end subroutine UrbanGetLayerTempShadedWall

  ! Hydrology getter subroutines (pervious road)
  subroutine UrbanGetSoilLiquidWaterPerviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetSoilLiquidWaterPerviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanGetSoilLiquidWaterPerviousRoad

  subroutine UrbanGetSoilVolumetricWaterPerviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetSoilVolumetricWaterPerviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanGetSoilVolumetricWaterPerviousRoad

  subroutine UrbanGetAquiferRechargeRatePerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetAquiferRechargeRatePerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetAquiferRechargeRatePerviousRoad

  subroutine UrbanGetWaterDeficitFluxPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetWaterDeficitFluxPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetWaterDeficitFluxPerviousRoad

  ! Soil fluxes compute subroutine
  subroutine UrbanComputeSoilFluxes(urban, status)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(out) :: status
    call UrbanComputeSoilFluxes_C(urban%ptr, status)
  end subroutine UrbanComputeSoilFluxes

  ! EflxSoilGrnd getter subroutines
  subroutine UrbanGetEflxSoilGrndRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetEflxSoilGrndRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetEflxSoilGrndRoof

  subroutine UrbanGetEflxSoilGrndImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetEflxSoilGrndImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetEflxSoilGrndImperviousRoad

  subroutine UrbanGetEflxSoilGrndPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetEflxSoilGrndPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetEflxSoilGrndPerviousRoad

  subroutine UrbanGetEflxSoilGrndSunlitWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetEflxSoilGrndSunlitWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetEflxSoilGrndSunlitWall

  subroutine UrbanGetEflxSoilGrndShadedWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetEflxSoilGrndShadedWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetEflxSoilGrndShadedWall

  ! QflxEvapGrnd getter subroutines
  subroutine UrbanGetQflxEvapGrndRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxEvapGrndRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxEvapGrndRoof

  subroutine UrbanGetQflxEvapGrndImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxEvapGrndImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxEvapGrndImperviousRoad

  subroutine UrbanGetQflxEvapGrndPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxEvapGrndPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxEvapGrndPerviousRoad

  ! QflxSubSnow getter subroutines
  subroutine UrbanGetQflxSubSnowRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxSubSnowRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxSubSnowRoof

  subroutine UrbanGetQflxSubSnowImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxSubSnowImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxSubSnowImperviousRoad

  subroutine UrbanGetQflxSubSnowPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxSubSnowPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxSubSnowPerviousRoad

  ! QflxDewSnow getter subroutines
  subroutine UrbanGetQflxDewSnowRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxDewSnowRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxDewSnowRoof

  subroutine UrbanGetQflxDewSnowImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxDewSnowImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxDewSnowImperviousRoad

  subroutine UrbanGetQflxDewSnowPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxDewSnowPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxDewSnowPerviousRoad

  ! QflxDewGrnd getter subroutines
  subroutine UrbanGetQflxDewGrndRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxDewGrndRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxDewGrndRoof

  subroutine UrbanGetQflxDewGrndImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxDewGrndImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxDewGrndImperviousRoad

  subroutine UrbanGetQflxDewGrndPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxDewGrndPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxDewGrndPerviousRoad

  ! Surface runoff setter functions (pervious road inputs)
  subroutine UrbanSetWtfactPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetWtfactPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetWtfactPerviousRoad

  subroutine UrbanSetFoverPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetFoverPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetFoverPerviousRoad

  subroutine UrbanSetFrostTablePerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetFrostTablePerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetFrostTablePerviousRoad

  subroutine UrbanSetZwtPerchedPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetZwtPerchedPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetZwtPerchedPerviousRoad

  ! Surface runoff compute function
  subroutine UrbanComputeSurfaceRunoff(urban, dtime, status)
    type(UrbanType), intent(in) :: urban
    real(c_double), value :: dtime
    integer(c_int), intent(out) :: status
    call UrbanComputeSurfaceRunoff_C(urban%ptr, dtime, status)
  end subroutine UrbanComputeSurfaceRunoff

  ! Surface runoff getter functions (QflxSurf — all five surfaces)
  subroutine UrbanGetQflxSurfRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxSurfRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxSurfRoof

  subroutine UrbanGetQflxSurfImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxSurfImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxSurfImperviousRoad

  subroutine UrbanGetQflxSurfPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxSurfPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxSurfPerviousRoad

  subroutine UrbanGetQflxSurfSunlitWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxSurfSunlitWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxSurfSunlitWall

  subroutine UrbanGetQflxSurfShadedWall(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetQflxSurfShadedWall_C(urban%ptr, values, length, status)
  end subroutine UrbanGetQflxSurfShadedWall

  subroutine UrbanGetTopH2OSoiLiqRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetTopH2OSoiLiqRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetTopH2OSoiLiqRoof

  subroutine UrbanGetTopH2OSoiLiqImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetTopH2OSoiLiqImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetTopH2OSoiLiqImperviousRoad

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

  ! Infiltration input setter subroutines
  subroutine UrbanSetSurfaceRunoffForPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetSurfaceRunoffForPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetSurfaceRunoffForPerviousRoad

  subroutine UrbanSetGroundEvapFluxForPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetGroundEvapFluxForPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetGroundEvapFluxForPerviousRoad

  ! Infiltration compute subroutine
  subroutine UrbanComputeInfiltration(urban, status)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(out) :: status
    call UrbanComputeInfiltration_C(urban%ptr, status)
  end subroutine UrbanComputeInfiltration

  ! Infiltration getter subroutine
  subroutine UrbanGetInfiltrationFluxPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetInfiltrationFluxPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetInfiltrationFluxPerviousRoad

  ! WaterTable setter subroutines
  subroutine UrbanSetAquiferWaterForPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetAquiferWaterForPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetAquiferWaterForPerviousRoad

  subroutine UrbanSetFracH2osfcForPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetFracH2osfcForPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetFracH2osfcForPerviousRoad

  subroutine UrbanSetDewGrndFluxForPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetDewGrndFluxForPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetDewGrndFluxForPerviousRoad

  subroutine UrbanSetDewSnowFluxForPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetDewSnowFluxForPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetDewSnowFluxForPerviousRoad

  subroutine UrbanSetSubSnowFluxForPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetSubSnowFluxForPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetSubSnowFluxForPerviousRoad

  subroutine UrbanSetQchargeForPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetQchargeForPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetQchargeForPerviousRoad

  ! WaterTable compute subroutine
  subroutine UrbanComputeWaterTable(urban, dtime, status)
    type(UrbanType), intent(in) :: urban
    real(c_double), value :: dtime
    integer(c_int), intent(out) :: status
    call UrbanComputeWaterTable_C(urban%ptr, dtime, status)
  end subroutine UrbanComputeWaterTable

  ! WaterTable getter subroutines (1D)
  subroutine UrbanGetWaterTableDepthPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetWaterTableDepthPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetWaterTableDepthPerviousRoad

  subroutine UrbanGetAquiferWaterPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetAquiferWaterPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetAquiferWaterPerviousRoad

  subroutine UrbanGetZwtPerchedPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetZwtPerchedPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetZwtPerchedPerviousRoad

  subroutine UrbanGetSubSnowFluxPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetSubSnowFluxPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetSubSnowFluxPerviousRoad

  subroutine UrbanGetDrainFluxPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetDrainFluxPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetDrainFluxPerviousRoad

  subroutine UrbanGetRsubSatPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetRsubSatPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetRsubSatPerviousRoad

  ! WaterTable getter subroutine (2D)
  subroutine UrbanGetSoilIceContentPerviousRoad(urban, values, size, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), dimension(2), intent(in) :: size
    integer(c_int), intent(out) :: status
    call UrbanGetSoilIceContentPerviousRoad_C(urban%ptr, values, size, status)
  end subroutine UrbanGetSoilIceContentPerviousRoad

  ! DewCondensation setter subroutines (roof)
  subroutine UrbanSetQflxDewGrndRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetQflxDewGrndRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanSetQflxDewGrndRoof

  subroutine UrbanSetQflxDewSnowRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetQflxDewSnowRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanSetQflxDewSnowRoof

  subroutine UrbanSetQflxSubSnowRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetQflxSubSnowRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanSetQflxSubSnowRoof

  ! DewCondensation setter subroutines (impervious road)
  subroutine UrbanSetQflxDewGrndImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetQflxDewGrndImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetQflxDewGrndImperviousRoad

  subroutine UrbanSetQflxDewSnowImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetQflxDewSnowImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetQflxDewSnowImperviousRoad

  subroutine UrbanSetQflxSubSnowImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetQflxSubSnowImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetQflxSubSnowImperviousRoad

  ! DewCondensation compute subroutine
  subroutine UrbanComputeDewCondensationRoofImperviousRoad(urban, dtime, status)
    type(UrbanType), intent(in) :: urban
    real(c_double), value :: dtime
    integer(c_int), intent(out) :: status
    call UrbanComputeDewCondensationRoofImperviousRoad_C(urban%ptr, dtime, status)
  end subroutine UrbanComputeDewCondensationRoofImperviousRoad

  ! DewCondensation getter subroutines (TopH2OSoiIce for roof and impervious road)
  subroutine UrbanGetTopH2OSoiIceRoof(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetTopH2OSoiIceRoof_C(urban%ptr, values, length, status)
  end subroutine UrbanGetTopH2OSoiIceRoof

  subroutine UrbanGetTopH2OSoiIceImperviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), intent(in) :: length
    integer(c_int), intent(out) :: status
    call UrbanGetTopH2OSoiIceImperviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanGetTopH2OSoiIceImperviousRoad

  ! Drainage init-time constant setter subroutines
  subroutine UrbanSetRsubTopGlobalMax(urban, value, status)
    type(UrbanType), intent(in) :: urban
    real(c_double), value :: value
    integer(c_int), intent(out) :: status
    call UrbanSetRsubTopGlobalMax_C(urban%ptr, value, status)
  end subroutine UrbanSetRsubTopGlobalMax

  subroutine UrbanSetPondmax(urban, value, status)
    type(UrbanType), intent(in) :: urban
    real(c_double), value :: value
    integer(c_int), intent(out) :: status
    call UrbanSetPondmax_C(urban%ptr, value, status)
  end subroutine UrbanSetPondmax

  subroutine UrbanSetWatmin(urban, value, status)
    type(UrbanType), intent(in) :: urban
    real(c_double), value :: value
    integer(c_int), intent(out) :: status
    call UrbanSetWatmin_C(urban%ptr, value, status)
  end subroutine UrbanSetWatmin

  subroutine UrbanSetEice(urban, value, status)
    type(UrbanType), intent(in) :: urban
    real(c_double), value :: value
    integer(c_int), intent(out) :: status
    call UrbanSetEice_C(urban%ptr, value, status)
  end subroutine UrbanSetEice

  ! Drainage timestep setter subroutines
  subroutine UrbanSetHkDepthForPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetHkDepthForPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetHkDepthForPerviousRoad

  subroutine UrbanSetTopoSlopeForPerviousRoad(urban, values, length, status)
    type(UrbanType), intent(in) :: urban
    type(c_ptr), value :: values
    integer(c_int), value :: length
    integer(c_int), intent(out) :: status
    call UrbanSetTopoSlopeForPerviousRoad_C(urban%ptr, values, length, status)
  end subroutine UrbanSetTopoSlopeForPerviousRoad

  ! Drainage compute subroutine
  subroutine UrbanComputeDrainage(urban, dtime, status)
    type(UrbanType), intent(in) :: urban
    real(c_double), value :: dtime
    integer(c_int), intent(out) :: status
    call UrbanComputeDrainage_C(urban%ptr, dtime, status)
  end subroutine UrbanComputeDrainage

end module urban_mod
