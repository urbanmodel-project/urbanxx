program urbanxx_driver_f
#include "finclude/UrbanFortranMacros.h"
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
  CallA(UrbanCreate(numLandunits, urban, status))

  if (mpi_rank == 0) then
    write(*,*) 'Successfully created Urban object with ', numLandunits, ' landunits'
  end if

  ! Set all urban parameters
  call SetUrbanParameters(urban, numLandunits, mpi_rank)

  ! Setup urban model (initialize temperatures and other setup tasks)
  call UrbanSetup(urban, status)
  if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
  if (mpi_rank == 0) then
    write(*,*) 'Completed urban model setup'
  end if

  ! Set hydrology boundary conditions after setup to override initialization defaults
  call SetHydrologyBoundaryConditions(urban, numLandunits, mpi_rank)

  ! Advance the model one time step
  call UrbanAdvance(urban, status)
  if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
  if (mpi_rank == 0) then
    write(*,*) 'Advanced model one time step'
  end if

  ! Destroy Urban object
  CallA(UrbanDestroy(urban, status))

  ! Finalize Kokkos
  call UrbanKokkosFinalize()

  ! Finalize MPI
  call MPI_Finalize(ierr)

contains

  subroutine SetCanyonHwr(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i
    real(c_double), allocatable, target :: canyonHwr(:)

    allocate(canyonHwr(numLandunits))
    do i = 1, numLandunits
      canyonHwr(i) = 4.80000019073486d0
    end do
    call UrbanSetCanyonHwr(urban, c_loc(canyonHwr), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set canyon height-to-width ratio'
    end if

    deallocate(canyonHwr)
  end subroutine SetCanyonHwr

  subroutine SetFracPervRoadOfTotalRoad(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i
    real(c_double), allocatable, target :: fracPervRoadOfTotalRoad(:)

    allocate(fracPervRoadOfTotalRoad(numLandunits))
    do i = 1, numLandunits
      fracPervRoadOfTotalRoad(i) = 0.16666667163372040d0
    end do
    call UrbanSetFracPervRoadOfTotalRoad(urban, c_loc(fracPervRoadOfTotalRoad), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set fraction of pervious road w.r.t. total road'
    end if

    deallocate(fracPervRoadOfTotalRoad)
  end subroutine SetFracPervRoadOfTotalRoad

  subroutine SetHeightParameters(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i
    real(c_double), allocatable, target :: forcHgtT(:)
    real(c_double), allocatable, target :: forcHgtU(:)
    real(c_double), allocatable, target :: zDTown(:)
    real(c_double), allocatable, target :: z0Town(:)
    real(c_double), allocatable, target :: htRoof(:)
    real(c_double), allocatable, target :: windHgtCanyon(:)

    ! Height parameter default values (in meters)
    real(c_double), parameter :: FORC_HGT_T_DEFAULT = 144.44377627618979d0
    real(c_double), parameter :: FORC_HGT_U_DEFAULT = 144.44377627618979d0  ! Same as T
    real(c_double), parameter :: Z_D_TOWN_DEFAULT = 113.96331622200367d0
    real(c_double), parameter :: Z_0_TOWN_DEFAULT = 0.48046005418613641d0
    real(c_double), parameter :: HT_ROOF_DEFAULT = 120.0d0
    real(c_double), parameter :: WIND_HGT_CANYON_DEFAULT = 60.0d0

    allocate(forcHgtT(numLandunits))
    allocate(forcHgtU(numLandunits))
    allocate(zDTown(numLandunits))
    allocate(z0Town(numLandunits))
    allocate(htRoof(numLandunits))
    allocate(windHgtCanyon(numLandunits))

    do i = 1, numLandunits
      forcHgtT(i) = FORC_HGT_T_DEFAULT
      forcHgtU(i) = FORC_HGT_U_DEFAULT
      zDTown(i) = Z_D_TOWN_DEFAULT
      z0Town(i) = Z_0_TOWN_DEFAULT
      htRoof(i) = HT_ROOF_DEFAULT
      windHgtCanyon(i) = WIND_HGT_CANYON_DEFAULT
    end do

    call UrbanSetForcHgtT(urban, c_loc(forcHgtT), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetForcHgtU(urban, c_loc(forcHgtU), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetZDTown(urban, c_loc(zDTown), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetZ0Town(urban, c_loc(z0Town), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetHtRoof(urban, c_loc(htRoof), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetWindHgtCanyon(urban, c_loc(windHgtCanyon), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set height parameters:'
      write(*,*) '  Forcing height (T):', FORC_HGT_T_DEFAULT, 'm'
      write(*,*) '  Forcing height (U):', FORC_HGT_U_DEFAULT, 'm'
      write(*,*) '  Zero displacement (town):', Z_D_TOWN_DEFAULT, 'm'
      write(*,*) '  Roughness length (town):', Z_0_TOWN_DEFAULT, 'm'
      write(*,*) '  Roof height:', HT_ROOF_DEFAULT, 'm'
      write(*,*) '  Wind height (canyon):', WIND_HGT_CANYON_DEFAULT, 'm'
    end if

    deallocate(forcHgtT)
    deallocate(forcHgtU)
    deallocate(zDTown)
    deallocate(z0Town)
    deallocate(htRoof)
    deallocate(windHgtCanyon)
  end subroutine SetHeightParameters

  subroutine SetBuildingTemperature(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i
    real(c_double), allocatable, target :: minTemp(:)
    real(c_double), allocatable, target :: maxTemp(:)
    real(c_double), allocatable, target :: wallThickness(:)
    real(c_double), allocatable, target :: roofThickness(:)
    real(c_double), parameter :: MIN_TEMP = 285.0d0  ! K
    real(c_double), parameter :: MAX_TEMP = 310.0d0  ! K
    real(c_double), parameter :: THICK_WALL = 0.286199986934662d0  ! m
    real(c_double), parameter :: THICK_ROOF = 0.217099994421005d0  ! m

    allocate(minTemp(numLandunits))
    allocate(maxTemp(numLandunits))

    do i = 1, numLandunits
      minTemp(i) = MIN_TEMP
      maxTemp(i) = MAX_TEMP
    end do

    call UrbanSetBuildingMinTemperature(urban, c_loc(minTemp), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetBuildingMaxTemperature(urban, c_loc(maxTemp), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set building temperature limits:'
      write(*,*) '  Min temperature:', MIN_TEMP, 'K'
      write(*,*) '  Max temperature:', MAX_TEMP, 'K'
    end if

    deallocate(minTemp)
    deallocate(maxTemp)

    ! Set building thickness parameters
    allocate(wallThickness(numLandunits))
    allocate(roofThickness(numLandunits))

    do i = 1, numLandunits
      wallThickness(i) = THICK_WALL
      roofThickness(i) = THICK_ROOF
    end do

    call UrbanSetBuildingWallThickness(urban, c_loc(wallThickness), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetBuildingRoofThickness(urban, c_loc(roofThickness), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set building thickness parameters:'
      write(*,*) '  Wall thickness:', THICK_WALL, 'm'
      write(*,*) '  Roof thickness:', THICK_ROOF, 'm'
    end if

    deallocate(wallThickness)
    deallocate(roofThickness)
  end subroutine SetBuildingTemperature

  subroutine SetSurfaceTemperatures(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i
    real(c_double), allocatable, target :: tempRoof(:)
    real(c_double), allocatable, target :: tempImpRoad(:)
    real(c_double), allocatable, target :: tempPervRoad(:)
    real(c_double), allocatable, target :: tempSunlitWall(:)
    real(c_double), allocatable, target :: tempShadedWall(:)
    real(c_double), parameter :: TEMP_ROOF_INIT = 292.0d0      ! K
    real(c_double), parameter :: TEMP_WALL_INIT = 292.0d0      ! K
    real(c_double), parameter :: TEMP_ROAD_INIT = 274.0d0      ! K

    ! Allocate arrays for surface temperatures
    allocate(tempRoof(numLandunits))
    allocate(tempImpRoad(numLandunits))
    allocate(tempPervRoad(numLandunits))
    allocate(tempSunlitWall(numLandunits))
    allocate(tempShadedWall(numLandunits))

    ! Initialize surface temperatures
    do i = 1, numLandunits
      tempRoof(i) = TEMP_ROOF_INIT
      tempImpRoad(i) = TEMP_ROAD_INIT
      tempPervRoad(i) = TEMP_ROAD_INIT
      tempSunlitWall(i) = TEMP_WALL_INIT
      tempShadedWall(i) = TEMP_WALL_INIT
    end do

    ! Set roof surface temperature
    call UrbanSetEffectiveSurfTempRoof(urban, c_loc(tempRoof), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    ! Set impervious road surface temperature
    call UrbanSetEffectiveSurfTempImperviousRoad(urban, c_loc(tempImpRoad), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    ! Set pervious road surface temperature
    call UrbanSetEffectiveSurfTempPerviousRoad(urban, c_loc(tempPervRoad), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    ! Set sunlit wall surface temperature
    call UrbanSetEffectiveSurfTempSunlitWall(urban, c_loc(tempSunlitWall), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    ! Set shaded wall surface temperature
    call UrbanSetEffectiveSurfTempShadedWall(urban, c_loc(tempShadedWall), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set surface temperatures:'
      write(*,*) '  Roof:', TEMP_ROOF_INIT, 'K'
      write(*,*) '  Impervious road:', TEMP_ROAD_INIT, 'K'
      write(*,*) '  Pervious road:', TEMP_ROAD_INIT, 'K'
      write(*,*) '  Sunlit wall:', TEMP_WALL_INIT, 'K'
      write(*,*) '  Shaded wall:', TEMP_WALL_INIT, 'K'
    end if

    deallocate(tempRoof)
    deallocate(tempImpRoad)
    deallocate(tempPervRoad)
    deallocate(tempSunlitWall)
    deallocate(tempShadedWall)
  end subroutine SetSurfaceTemperatures

  subroutine SetLayerTemperatures(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i, j, idx
    integer(c_int), parameter :: NUM_URBAN_LAYERS = 5
    integer(c_int), parameter :: NUM_SOIL_LAYERS = 15
    integer(c_int), dimension(2) :: size2D_urban, size2D_soil
    logical(c_bool) :: isLayoutLeft
    real(c_double), allocatable, target :: tempRoof(:)
    real(c_double), allocatable, target :: tempImpRoad(:)
    real(c_double), allocatable, target :: tempPervRoad(:)
    real(c_double), allocatable, target :: tempSunlitWall(:)
    real(c_double), allocatable, target :: tempShadedWall(:)
    real(c_double), parameter :: TEMP_ROOF_INIT = 292.0d0      ! K
    real(c_double), parameter :: TEMP_WALL_INIT = 292.0d0      ! K
    real(c_double), parameter :: TEMP_ROAD_INIT = 274.0d0      ! K

    ! Check memory layout
    isLayoutLeft = UrbanKokkosIsLayoutLeft()

    ! Allocate arrays for layer temperatures
    ! Urban surfaces: (numLandunits, NUM_URBAN_LAYERS)
    allocate(tempRoof(numLandunits * NUM_URBAN_LAYERS))
    allocate(tempSunlitWall(numLandunits * NUM_URBAN_LAYERS))
    allocate(tempShadedWall(numLandunits * NUM_URBAN_LAYERS))
    
    ! Road surfaces: (numLandunits, NUM_SOIL_LAYERS)
    allocate(tempImpRoad(numLandunits * NUM_SOIL_LAYERS))
    allocate(tempPervRoad(numLandunits * NUM_SOIL_LAYERS))

    ! Fill arrays based on memory layout
    if (isLayoutLeft) then
      ! LayoutLeft: iterate layers in outer loop, landunits in inner loop
      idx = 0
      do j = 1, NUM_URBAN_LAYERS
        do i = 1, numLandunits
          idx = idx + 1
          tempRoof(idx) = TEMP_ROOF_INIT
          tempSunlitWall(idx) = TEMP_WALL_INIT
          tempShadedWall(idx) = TEMP_WALL_INIT
        end do
      end do
      
      idx = 0
      do j = 1, NUM_SOIL_LAYERS
        do i = 1, numLandunits
          idx = idx + 1
          tempImpRoad(idx) = TEMP_ROAD_INIT
          tempPervRoad(idx) = TEMP_ROAD_INIT
        end do
      end do
    else
      ! LayoutRight: iterate landunits in outer loop, layers in inner loop
      idx = 0
      do i = 1, numLandunits
        do j = 1, NUM_URBAN_LAYERS
          idx = idx + 1
          tempRoof(idx) = TEMP_ROOF_INIT
          tempSunlitWall(idx) = TEMP_WALL_INIT
          tempShadedWall(idx) = TEMP_WALL_INIT
        end do
      end do
      
      idx = 0
      do i = 1, numLandunits
        do j = 1, NUM_SOIL_LAYERS
          idx = idx + 1
          tempImpRoad(idx) = TEMP_ROAD_INIT
          tempPervRoad(idx) = TEMP_ROAD_INIT
        end do
      end do
    end if

    ! Set size arrays
    size2D_urban(1) = numLandunits
    size2D_urban(2) = NUM_URBAN_LAYERS
    size2D_soil(1) = numLandunits
    size2D_soil(2) = NUM_SOIL_LAYERS

    ! Set roof layer temperatures
    call UrbanSetLayerTempRoof(urban, c_loc(tempRoof), size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    ! Set impervious road layer temperatures
    call UrbanSetLayerTempImperviousRoad(urban, c_loc(tempImpRoad), &
      size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    ! Set pervious road layer temperatures
    call UrbanSetLayerTempPerviousRoad(urban, c_loc(tempPervRoad), &
      size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    ! Set sunlit wall layer temperatures
    call UrbanSetLayerTempSunlitWall(urban, c_loc(tempSunlitWall), &
      size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    ! Set shaded wall layer temperatures
    call UrbanSetLayerTempShadedWall(urban, c_loc(tempShadedWall), &
      size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      if (isLayoutLeft) then
        write(*,*) 'Set layer temperatures (LayoutLeft):'
      else
        write(*,*) 'Set layer temperatures (LayoutRight):'
      end if
      write(*,*) '  Roof layers:', TEMP_ROOF_INIT, 'K (', NUM_URBAN_LAYERS, 'layers)'
      write(*,*) '  Impervious road layers:', TEMP_ROAD_INIT, 'K (', NUM_SOIL_LAYERS, 'layers)'
      write(*,*) '  Pervious road layers:', TEMP_ROAD_INIT, 'K (', NUM_SOIL_LAYERS, 'layers)'
      write(*,*) '  Sunlit wall layers:', TEMP_WALL_INIT, 'K (', NUM_URBAN_LAYERS, 'layers)'
      write(*,*) '  Shaded wall layers:', TEMP_WALL_INIT, 'K (', NUM_URBAN_LAYERS, 'layers)'
    end if

    deallocate(tempRoof)
    deallocate(tempImpRoad)
    deallocate(tempPervRoad)
    deallocate(tempSunlitWall)
    deallocate(tempShadedWall)
  end subroutine SetLayerTemperatures

  subroutine SetWtRoof(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i
    real(c_double), allocatable, target :: wtRoof(:)
    real(c_double), parameter :: WT_ROOF_DEFAULT = 0.69999998807907104d0

    allocate(wtRoof(numLandunits))
    do i = 1, numLandunits
      wtRoof(i) = WT_ROOF_DEFAULT
    end do

    call UrbanSetWtRoof(urban, c_loc(wtRoof), &
      numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set roof weight:', WT_ROOF_DEFAULT
    end if

    deallocate(wtRoof)
  end subroutine SetWtRoof

  subroutine SetAlbedo(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i, ilandunit, iband, itype, idx
    integer(c_int) :: numBands, numTypes, totalSize3D
    integer(c_int), dimension(3) :: size3D
    logical(c_bool) :: isLayoutLeft
    real(c_double), allocatable, target :: albedoPerviousRoad(:)
    real(c_double), allocatable, target :: albedoImperviousRoad(:)
    real(c_double), allocatable, target :: albedoSunlitWall(:)
    real(c_double), allocatable, target :: albedoShadedWall(:)
    real(c_double), allocatable, target :: albedoRoof(:)
    real(c_double), parameter :: ALB_IMPROAD = 0.230000004172325d0
    real(c_double), parameter :: ALB_PERROAD = 0.0799999982118607d0
    real(c_double), parameter :: ALB_ROOF = 0.254999995231628d0
    real(c_double), parameter :: ALB_WALL = 0.200000002980232d0

    numBands = 2  ! VIS, NIR
    numTypes = 2  ! Direct, Diffuse
    size3D = [numLandunits, numBands, numTypes]
    totalSize3D = numLandunits * numBands * numTypes

    allocate(albedoPerviousRoad(totalSize3D))
    allocate(albedoImperviousRoad(totalSize3D))
    allocate(albedoSunlitWall(totalSize3D))
    allocate(albedoShadedWall(totalSize3D))
    allocate(albedoRoof(totalSize3D))

    ! Check Kokkos memory layout
    isLayoutLeft = UrbanKokkosIsLayoutLeft()

    if (isLayoutLeft) then
      ! LayoutLeft: First dimension (landunits) varies fastest
      ! Iterate: itype (outer), iband (middle), landunits (inner)
      idx = 1
      do itype = 0, numTypes - 1
        do iband = 0, numBands - 1
          do ilandunit = 0, numLandunits - 1
            albedoPerviousRoad(idx) = ALB_PERROAD
            albedoImperviousRoad(idx) = ALB_IMPROAD
            albedoSunlitWall(idx) = ALB_WALL
            albedoShadedWall(idx) = ALB_WALL
            albedoRoof(idx) = ALB_ROOF
            idx = idx + 1
          end do
        end do
      end do
    else
      ! LayoutRight: Last dimension (types) varies fastest
      ! Iterate: landunits (outer), iband (middle), itype (inner)
      do ilandunit = 0, numLandunits - 1
        do iband = 0, numBands - 1
          do itype = 0, numTypes - 1
            idx = ilandunit * numBands * numTypes + iband * numTypes + itype + 1  ! +1 for Fortran 1-indexing
            albedoPerviousRoad(idx) = ALB_PERROAD
            albedoImperviousRoad(idx) = ALB_IMPROAD
            albedoSunlitWall(idx) = ALB_WALL
            albedoShadedWall(idx) = ALB_WALL
            albedoRoof(idx) = ALB_ROOF
          end do
        end do
      end do
    end if

    call UrbanSetAlbedoPerviousRoad(urban, &
      c_loc(albedoPerviousRoad), size3D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAlbedoImperviousRoad(urban, &
      c_loc(albedoImperviousRoad), size3D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAlbedoSunlitWall(urban, &
      c_loc(albedoSunlitWall), size3D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAlbedoShadedWall(urban, &
      c_loc(albedoShadedWall), size3D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAlbedoRoof(urban, c_loc(albedoRoof), &
      size3D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set albedo values for all surfaces'
    end if

    deallocate(albedoPerviousRoad)
    deallocate(albedoImperviousRoad)
    deallocate(albedoSunlitWall)
    deallocate(albedoShadedWall)
    deallocate(albedoRoof)
  end subroutine SetAlbedo

  subroutine SetEmissivity(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i
    real(c_double), allocatable, target :: emissivityPerviousRoad(:)
    real(c_double), allocatable, target :: emissivityImperviousRoad(:)
    real(c_double), allocatable, target :: emissivityWall(:)
    real(c_double), allocatable, target :: emissivityRoof(:)
    real(c_double), parameter :: EMISS_ROOF = 0.90600001811981201d0
    real(c_double), parameter :: EMISS_IMPROAD = 0.87999999523162842d0
    real(c_double), parameter :: EMISS_PERROAD = 0.94999998807907104d0
    real(c_double), parameter :: EMISS_WALL = 0.90200001001358032d0

    allocate(emissivityPerviousRoad(numLandunits))
    allocate(emissivityImperviousRoad(numLandunits))
    allocate(emissivityWall(numLandunits))
    allocate(emissivityRoof(numLandunits))

    do i = 1, numLandunits
      emissivityPerviousRoad(i) = EMISS_PERROAD
      emissivityImperviousRoad(i) = EMISS_IMPROAD
      emissivityWall(i) = EMISS_WALL
      emissivityRoof(i) = EMISS_ROOF
    end do

    call UrbanSetEmissivityPerviousRoad(urban, &
      c_loc(emissivityPerviousRoad), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetEmissivityImperviousRoad(urban, &
      c_loc(emissivityImperviousRoad), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetEmissivityWall(urban, &
      c_loc(emissivityWall), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetEmissivityRoof(urban, &
      c_loc(emissivityRoof), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set emissivity values for all surfaces'
    end if

    deallocate(emissivityPerviousRoad)
    deallocate(emissivityImperviousRoad)
    deallocate(emissivityWall)
    deallocate(emissivityRoof)
  end subroutine SetEmissivity

  subroutine SetNumberOfActiveLayersImperviousRoad(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i, urban_density_class
    integer(c_int), parameter :: NUM_URBAN_DENSITY_CLASSES = 3
    integer(c_int), dimension(3) :: nlevImproadClasses
    real(c_double), allocatable, target :: numActiveLayers(:)

    ! Number of active layers for each urban density class
    ! 0=Tall Building District, 1=High Density, 2=Medium Density
    nlevImproadClasses = (/ 3, 2, 2 /)

    allocate(numActiveLayers(numLandunits))
    do i = 1, numLandunits
      urban_density_class = mod(i-1, NUM_URBAN_DENSITY_CLASSES) + 1
      numActiveLayers(i) = real(nlevImproadClasses(urban_density_class), c_double)
    end do
    call UrbanSetNumberOfActiveLayersImperviousRoad(urban, &
      c_loc(numActiveLayers), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set number of active layers for impervious road (3, 2, 2 for density classes)'
    end if

    deallocate(numActiveLayers)
  end subroutine SetNumberOfActiveLayersImperviousRoad

  subroutine SetThermalConductivity(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i, k, urban_density_class, idx
    integer(c_int), parameter :: NUM_ROAD_LEVELS = 15
    integer(c_int), parameter :: NUM_URBAN_LEVELS = 5
    integer(c_int) :: totalSizeRoad, totalSizeUrban
    integer(c_int), dimension(2) :: size2D_road, size2D_urban
    real(c_double), allocatable, target :: tkRoad(:)
    real(c_double), allocatable, target :: tkWall(:)
    real(c_double), allocatable, target :: tkRoof(:)
    logical(c_bool) :: isLayoutLeft
    ! Thermal conductivity values for road (15 levels) across 3 urban density classes
    ! [layer][urban_density_class]: 0=Tall Building District, 1=High Density, 2=Medium Density
    real(c_double), parameter :: TK_BEDROCK = 3.0d0  ! thermal conductivity of bedrock [W/m-K]
    real(c_double), dimension(15,3) :: tkRoadLevels
    real(c_double), dimension(5,3) :: tkWallLevels
    real(c_double), dimension(5,3) :: tkRoofLevels

    tkRoadLevels = reshape((/ &
      1.8999999761581421d0, 1.6699999570846558d0, 1.6699999570846558d0, &
      0.56000000238418579d0, 0.56000000238418579d0, 0.56000000238418579d0, &
      0.36000001430511475d0, 0.21664454245689402d0, 0.21664454245689402d0, &
      0.21572236500511191d0, 0.21572236500511191d0, 0.21572236500511191d0, &
      0.21389214262493184d0, 0.21389214262493184d0, 0.21389214262493184d0, &
      0.21208052418425166d0, 0.21208052418425166d0, 0.21208052418425166d0, &
      0.21028722745802339d0, 0.21028722745802339d0, 0.21028722745802339d0, &
      0.21028722745802339d0, 0.21028722745802339d0, 0.21028722745802339d0, &
      0.21298402568603625d0, 0.21298402568603625d0, 0.21298402568603625d0, &
      0.21664454245689402d0, 0.21664454245689402d0, 0.21664454245689402d0, &
      TK_BEDROCK, TK_BEDROCK, TK_BEDROCK, &
      TK_BEDROCK, TK_BEDROCK, TK_BEDROCK, &
      TK_BEDROCK, TK_BEDROCK, TK_BEDROCK, &
      TK_BEDROCK, TK_BEDROCK, TK_BEDROCK, &
      TK_BEDROCK, TK_BEDROCK, TK_BEDROCK /), shape(tkRoadLevels))
    
    tkWallLevels = reshape((/ &
      1.44716906547546d0, 1.06582415103912d0, 0.970157384872437d0, &
      1.44716906547546d0, 1.06582415103912d0, 0.970157384872437d0, &
      1.44716906547546d0, 1.06582415103912d0, 0.970157384872437d0, &
      1.44716906547546d0, 1.06582415103912d0, 0.970157384872437d0, &
      1.44716906547546d0, 1.06582415103912d0, 0.970157384872437d0 /), shape(tkWallLevels))
    
    tkRoofLevels = reshape((/ &
      0.503093481063843d0, 0.094768725335598d0, 0.127733826637268d0, &
      0.503093481063843d0, 0.094768725335598d0, 0.127733826637268d0, &
      0.503093481063843d0, 0.094768725335598d0, 0.127733826637268d0, &
      0.503093481063843d0, 0.094768725335598d0, 0.127733826637268d0, &
      0.503093481063843d0, 0.094768725335598d0, 0.127733826637268d0 /), shape(tkRoofLevels))

    totalSizeRoad = numLandunits * NUM_ROAD_LEVELS
    totalSizeUrban = numLandunits * NUM_URBAN_LEVELS
    size2D_road(1) = numLandunits
    size2D_road(2) = NUM_ROAD_LEVELS
    size2D_urban(1) = numLandunits
    size2D_urban(2) = NUM_URBAN_LEVELS

    allocate(tkRoad(totalSizeRoad))
    allocate(tkWall(totalSizeUrban))
    allocate(tkRoof(totalSizeUrban))

    ! Check Kokkos memory layout
    isLayoutLeft = UrbanKokkosIsLayoutLeft()

    if (isLayoutLeft) then
      ! LayoutLeft: First dimension (landunits) varies fastest
      ! Iterate: layer (outer), landunits (inner)
      idx = 1
      do k = 1, NUM_ROAD_LEVELS
        do i = 1, numLandunits
          urban_density_class = mod(i-1, 3) + 1
          tkRoad(idx) = tkRoadLevels(k, urban_density_class)
          idx = idx + 1
        end do
      end do

      idx = 1
      do k = 1, NUM_URBAN_LEVELS
        do i = 1, numLandunits
          urban_density_class = mod(i-1, 3) + 1
          tkWall(idx) = tkWallLevels(k, urban_density_class)
          tkRoof(idx) = tkRoofLevels(k, urban_density_class)
          idx = idx + 1
        end do
      end do
    else
      ! LayoutRight: Last dimension (layers) varies fastest
      ! Iterate: landunits (outer), layer (inner)
      idx = 1
      do i = 1, numLandunits
        urban_density_class = mod(i-1, 3) + 1
        do k = 1, NUM_ROAD_LEVELS
          tkRoad(idx) = tkRoadLevels(k, urban_density_class)
          idx = idx + 1
        end do
      end do

      idx = 1
      do i = 1, numLandunits
        urban_density_class = mod(i-1, 3) + 1
        do k = 1, NUM_URBAN_LEVELS
          tkWall(idx) = tkWallLevels(k, urban_density_class)
          tkRoof(idx) = tkRoofLevels(k, urban_density_class)
          idx = idx + 1
        end do
      end do
    end if

    call UrbanSetThermalConductivityRoad(urban, &
      c_loc(tkRoad), size2D_road, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetThermalConductivityWall(urban, &
      c_loc(tkWall), size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetThermalConductivityRoof(urban, &
      c_loc(tkRoof), size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set thermal conductivity values for all surfaces'
    end if

    deallocate(tkRoad)
    deallocate(tkWall)
    deallocate(tkRoof)
  end subroutine SetThermalConductivity

  subroutine SetHeatCapacity(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i, k, urban_density_class, idx
    integer(c_int), parameter :: NUM_ROAD_LEVELS = 15
    integer(c_int), parameter :: NUM_URBAN_LEVELS = 5
    integer(c_int) :: totalSizeRoad, totalSizeUrban
    integer(c_int), dimension(2) :: size2D_road, size2D_urban
    real(c_double), allocatable, target :: cvRoad(:)
    real(c_double), allocatable, target :: cvWall(:)
    real(c_double), allocatable, target :: cvRoof(:)
    logical(c_bool) :: isLayoutLeft
    ! Heat capacity values for road (15 levels) across 3 urban density classes
    ! [layer][urban_density_class]: 0=Tall Building District, 1=High Density, 2=Medium Density
    real(c_double), parameter :: CV_BEDROCK = 2.0d6  ! heat capacity of bedrock [J/m^3/K]
    real(c_double), dimension(15,3) :: cvRoadLevels
    real(c_double), dimension(5,3) :: cvWallLevels
    real(c_double), dimension(5,3) :: cvRoofLevels

    cvRoadLevels = reshape((/ &
      2100000.0d0, 2060470.625d0, 2060470.625d0, &
      1773000.0d0, 1712294.75d0, 1712294.75d0, &
      1545600.0d0, CV_BEDROCK, CV_BEDROCK, &
      CV_BEDROCK, CV_BEDROCK, CV_BEDROCK, &
      CV_BEDROCK, CV_BEDROCK, CV_BEDROCK, &
      CV_BEDROCK, CV_BEDROCK, CV_BEDROCK, &
      CV_BEDROCK, CV_BEDROCK, CV_BEDROCK, &
      CV_BEDROCK, CV_BEDROCK, CV_BEDROCK, &
      CV_BEDROCK, CV_BEDROCK, CV_BEDROCK, &
      CV_BEDROCK, CV_BEDROCK, CV_BEDROCK, &
      CV_BEDROCK, CV_BEDROCK, CV_BEDROCK, &
      CV_BEDROCK, CV_BEDROCK, CV_BEDROCK, &
      CV_BEDROCK, CV_BEDROCK, CV_BEDROCK, &
      CV_BEDROCK, CV_BEDROCK, CV_BEDROCK, &
      CV_BEDROCK, CV_BEDROCK, CV_BEDROCK /), shape(cvRoadLevels))
    
    cvWallLevels = reshape((/ &
      1079394.75d0, 957632.8125d0, 899827.1875d0, &
      1079394.75d0, 957632.8125d0, 899827.1875d0, &
      1079394.75d0, 957632.8125d0, 899827.1875d0, &
      1079394.75d0, 957632.8125d0, 899827.1875d0, &
      1079394.75d0, 957632.8125d0, 899827.1875d0 /), shape(cvWallLevels))
    
    cvRoofLevels = reshape((/ &
      570998.0d0, 646213.375d0, 862451.375d0, &
      570998.0d0, 646213.375d0, 862451.375d0, &
      570998.0d0, 646213.375d0, 862451.375d0, &
      570998.0d0, 646213.375d0, 862451.375d0, &
      570998.0d0, 646213.375d0, 862451.375d0 /), shape(cvRoofLevels))

    totalSizeRoad = numLandunits * NUM_ROAD_LEVELS
    totalSizeUrban = numLandunits * NUM_URBAN_LEVELS
    size2D_road(1) = numLandunits
    size2D_road(2) = NUM_ROAD_LEVELS
    size2D_urban(1) = numLandunits
    size2D_urban(2) = NUM_URBAN_LEVELS

    allocate(cvRoad(totalSizeRoad))
    allocate(cvWall(totalSizeUrban))
    allocate(cvRoof(totalSizeUrban))

    ! Check Kokkos memory layout
    isLayoutLeft = UrbanKokkosIsLayoutLeft()

    if (isLayoutLeft) then
      ! LayoutLeft: First dimension (landunits) varies fastest
      ! Iterate: layer (outer), landunits (inner)
      idx = 1
      do k = 1, NUM_ROAD_LEVELS
        do i = 1, numLandunits
          urban_density_class = mod(i-1, 3) + 1
          cvRoad(idx) = cvRoadLevels(k, urban_density_class)
          idx = idx + 1
        end do
      end do

      idx = 1
      do k = 1, NUM_URBAN_LEVELS
        do i = 1, numLandunits
          urban_density_class = mod(i-1, 3) + 1
          cvWall(idx) = cvWallLevels(k, urban_density_class)
          cvRoof(idx) = cvRoofLevels(k, urban_density_class)
          idx = idx + 1
        end do
      end do
    else
      ! LayoutRight: Last dimension (layers) varies fastest
      ! Iterate: landunits (outer), layer (inner)
      idx = 1
      do i = 1, numLandunits
        urban_density_class = mod(i-1, 3) + 1
        do k = 1, NUM_ROAD_LEVELS
          cvRoad(idx) = cvRoadLevels(k, urban_density_class)
          idx = idx + 1
        end do
      end do

      idx = 1
      do i = 1, numLandunits
        urban_density_class = mod(i-1, 3) + 1
        do k = 1, NUM_URBAN_LEVELS
          cvWall(idx) = cvWallLevels(k, urban_density_class)
          cvRoof(idx) = cvRoofLevels(k, urban_density_class)
          idx = idx + 1
        end do
      end do
    end if

    call UrbanSetHeatCapacityRoad(urban, &
      c_loc(cvRoad), size2D_road, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetHeatCapacityWall(urban, &
      c_loc(cvWall), size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetHeatCapacityRoof(urban, &
      c_loc(cvRoof), size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set heat capacity values for all surfaces'
    end if

    deallocate(cvRoad)
    deallocate(cvWall)
    deallocate(cvRoof)
  end subroutine SetHeatCapacity

  subroutine SetSoilProperties(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i, k, srcLayer, idx
    integer(c_int), parameter :: NUM_SOIL_LEVELS = 15
    integer(c_int) :: totalSize
    integer(c_int), dimension(2) :: size2D
    logical(c_bool) :: isLayoutLeft
    real(c_double), allocatable, target :: sand(:)
    real(c_double), allocatable, target :: clay(:)
    real(c_double), allocatable, target :: organic(:)
    ! Soil property values for first 10 layers (layers 11-15 use layer 10 values)
    real(c_double), dimension(10) :: sandLevels
    real(c_double), dimension(10) :: clayLevels
    real(c_double), dimension(10) :: organicLevels

    sandLevels = (/ 46.0d0, 46.0d0, 44.0d0, 43.0d0, 41.0d0, 39.0d0, 37.0d0, 37.0d0, 40.0d0, 44.0d0 /)
    clayLevels = (/ 35.0d0, 35.0d0, 37.0d0, 39.0d0, 42.0d0, 44.0d0, 46.0d0, 45.0d0, 41.0d0, 42.0d0 /)
    organicLevels = (/ 25.2229820327902d0, 25.700711396596d0, 22.091324741929d0, 18.1150405358844d0, &
                       14.5211498497041d0, 11.4998502546828d0, 9.04501744160207d0, 7.08594278159189d0, &
                       0.0d0, 0.0d0 /)

    totalSize = numLandunits * NUM_SOIL_LEVELS
    size2D(1) = numLandunits
    size2D(2) = NUM_SOIL_LEVELS

    allocate(sand(totalSize))
    allocate(clay(totalSize))
    allocate(organic(totalSize))

    ! Check Kokkos memory layout
    isLayoutLeft = UrbanKokkosIsLayoutLeft()

    if (isLayoutLeft) then
      ! LayoutLeft: First dimension (landunits) varies fastest
      ! Iterate: layer (outer), landunits (inner)
      idx = 1
      do k = 1, NUM_SOIL_LEVELS
        do i = 1, numLandunits
          ! Use layer 10 values for layers 11-15
          if (k <= 10) then
            srcLayer = k
          else
            srcLayer = 10
          end if
          sand(idx) = sandLevels(srcLayer)
          clay(idx) = clayLevels(srcLayer)
          organic(idx) = organicLevels(srcLayer)
          idx = idx + 1
        end do
      end do
    else
      ! LayoutRight: Last dimension (layers) varies fastest
      ! Iterate: landunits (outer), layer (inner)
      do i = 1, numLandunits
        do k = 1, NUM_SOIL_LEVELS
          ! Use layer 10 values for layers 11-15
          if (k <= 10) then
            srcLayer = k
          else
            srcLayer = 10
          end if
          sand((i-1) * NUM_SOIL_LEVELS + k) = sandLevels(srcLayer)
          clay((i-1) * NUM_SOIL_LEVELS + k) = clayLevels(srcLayer)
          organic((i-1) * NUM_SOIL_LEVELS + k) = organicLevels(srcLayer)
        end do
      end do
    end if

    call UrbanSetSandPerviousRoad(urban, c_loc(sand), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetClayPerviousRoad(urban, c_loc(clay), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetOrganicPerviousRoad(urban, c_loc(organic), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set soil properties for pervious road (sand, clay, organic)'
    end if

    deallocate(sand)
    deallocate(clay)
    deallocate(organic)
  end subroutine SetSoilProperties

  subroutine SetAtmosphericForcing(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i, iband, itype, idx
    integer(c_int) :: numBands, numTypes, totalSize3D
    integer(c_int), dimension(3) :: size3D
    logical(c_bool) :: isLayoutLeft
    real(c_double), allocatable, target :: atmTemp(:)
    real(c_double), allocatable, target :: atmPotTemp(:)
    real(c_double), allocatable, target :: atmRho(:)
    real(c_double), allocatable, target :: atmSpcHumd(:)
    real(c_double), allocatable, target :: atmPress(:)
    real(c_double), allocatable, target :: atmWindU(:)
    real(c_double), allocatable, target :: atmWindV(:)
    real(c_double), allocatable, target :: atmCoszen(:)
    real(c_double), allocatable, target :: atmFracSnow(:)
    real(c_double), allocatable, target :: atmLongwave(:)
    real(c_double), allocatable, target :: atmShortwave(:)

    ! Atmospheric forcing constants
    real(c_double), parameter :: TEMP_AIR = 297.26422743678319d0
    real(c_double), parameter :: TH_AIR = 297.26422743678319d0
    real(c_double), parameter :: RHO_AIR = 1.1382761848551157d0
    real(c_double), parameter :: Q_AIR = 1.9217052569985755d-2
    real(c_double), parameter :: PBOT_AIR = 98260.450580263219d0
    real(c_double), parameter :: WIND_U = 0.52482489069830152d0
    real(c_double), parameter :: WIND_V = 0.0d0
    real(c_double), parameter :: COSZEN = 7.9054122593736065d-3
    real(c_double), parameter :: SNOW = 0.0d0
    real(c_double), parameter :: LWDOWN = 432.79580327766450d0
    real(c_double), parameter :: SWDOWN = 0.0d0

    ! Allocate arrays
    allocate(atmTemp(numLandunits))
    allocate(atmPotTemp(numLandunits))
    allocate(atmRho(numLandunits))
    allocate(atmSpcHumd(numLandunits))
    allocate(atmPress(numLandunits))
    allocate(atmWindU(numLandunits))
    allocate(atmWindV(numLandunits))
    allocate(atmCoszen(numLandunits))
    allocate(atmFracSnow(numLandunits))
    allocate(atmLongwave(numLandunits))

    numBands = 2  ! VIS, NIR
    numTypes = 2  ! Direct, Diffuse
    size3D = [numLandunits, numBands, numTypes]
    totalSize3D = numLandunits * numBands * numTypes
    allocate(atmShortwave(totalSize3D))

    ! Fill arrays with constant values
    do i = 1, numLandunits
      atmTemp(i) = TEMP_AIR
      atmPotTemp(i) = TH_AIR
      atmRho(i) = RHO_AIR
      atmSpcHumd(i) = Q_AIR
      atmPress(i) = PBOT_AIR
      atmWindU(i) = WIND_U
      atmWindV(i) = WIND_V
      atmCoszen(i) = COSZEN
      atmFracSnow(i) = SNOW
      atmLongwave(i) = LWDOWN
    end do

    ! Check Kokkos memory layout
    isLayoutLeft = UrbanKokkosIsLayoutLeft()

    if (isLayoutLeft) then
      ! LayoutLeft: First dimension (landunits) varies fastest
      ! Iterate: itype (outer), iband (middle), landunits (inner)
      idx = 1
      do itype = 1, numTypes
        do iband = 1, numBands
          do i = 1, numLandunits
            atmShortwave(idx) = SWDOWN
            idx = idx + 1
          end do
        end do
      end do
    else
      ! LayoutRight: Last dimension (types) varies fastest
      ! Iterate: landunits (outer), iband (middle), itype (inner)
      do i = 1, totalSize3D
        atmShortwave(i) = SWDOWN
      end do
    end if

    ! Set atmospheric forcing
    call UrbanSetAtmTemp(urban, c_loc(atmTemp), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAtmPotTemp(urban, c_loc(atmPotTemp), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAtmRho(urban, c_loc(atmRho), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAtmSpcHumd(urban, c_loc(atmSpcHumd), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAtmPress(urban, c_loc(atmPress), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAtmWindU(urban, c_loc(atmWindU), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAtmWindV(urban, c_loc(atmWindV), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAtmCoszen(urban, c_loc(atmCoszen), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAtmFracSnow(urban, c_loc(atmFracSnow), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAtmLongwaveDown(urban, c_loc(atmLongwave), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    call UrbanSetAtmShortwaveDown(urban, c_loc(atmShortwave), size3D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    if (mpi_rank == 0) then
      write(*,*) 'Set atmospheric forcing:'
      write(*,*) '  Temperature:', TEMP_AIR, 'K'
      write(*,*) '  Pressure:', PBOT_AIR, 'Pa'
      write(*,*) '  Density:', RHO_AIR, 'kg/m^3'
      write(*,*) '  Specific humidity:', Q_AIR, 'kg/kg'
      write(*,*) '  Wind U:', WIND_U, 'm/s'
      write(*,*) '  Wind V:', WIND_V, 'm/s'
      write(*,*) '  Cosine zenith:', COSZEN
      write(*,*) '  Snow fraction:', SNOW
      write(*,*) '  Longwave down:', LWDOWN, 'W/m^2'
      write(*,*) '  Shortwave down:', SWDOWN, 'W/m^2'
    end if

    ! Free arrays
    deallocate(atmTemp)
    deallocate(atmPotTemp)
    deallocate(atmRho)
    deallocate(atmSpcHumd)
    deallocate(atmPress)
    deallocate(atmWindU)
    deallocate(atmWindV)
    deallocate(atmCoszen)
    deallocate(atmFracSnow)
    deallocate(atmLongwave)
    deallocate(atmShortwave)
  end subroutine SetAtmosphericForcing

  subroutine SetHydrologyBoundaryConditions(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status, i, j, idx
    integer(c_int), parameter :: NUM_SOIL_LAYERS = 15
    integer(c_int) :: size2D(2)
    integer(c_int) :: totalSize
    real(c_double), allocatable, target :: zwt(:)
    real(c_double), allocatable, target :: qflxInfl(:)
    real(c_double), allocatable, target :: qflxTran(:)
    real(c_double), allocatable, target :: fwet(:)
    logical(c_bool) :: isLayoutLeft

    ! Constants
    real(c_double), parameter :: ZWT_INITIAL = 4.8018819123227204d0  ! m (initial water table depth)
    real(c_double), parameter :: QFLX_INFL = 0.0d0                   ! mm/s (infiltration flux)
    real(c_double), parameter :: FWET_INITIAL = 0.0d0                ! fraction of surface that is wet [-]

    ! Set water table depth (1D: per landunit)
    allocate(zwt(numLandunits))
    do i = 1, numLandunits
      zwt(i) = ZWT_INITIAL
    end do
    call UrbanSetWaterTableDepth(urban, c_loc(zwt), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    deallocate(zwt)

    ! Set infiltration flux (1D: per landunit)
    allocate(qflxInfl(numLandunits))
    do i = 1, numLandunits
      qflxInfl(i) = QFLX_INFL
    end do
    call UrbanSetInfiltrationFlux(urban, c_loc(qflxInfl), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    deallocate(qflxInfl)

    ! Set transpiration flux (2D: per landunit and soil layer)
    size2D(1) = numLandunits
    size2D(2) = NUM_SOIL_LAYERS
    totalSize = numLandunits * NUM_SOIL_LAYERS
    allocate(qflxTran(totalSize))

    ! Check memory layout
    isLayoutLeft = UrbanKokkosIsLayoutLeft()

    ! Initialize to zero for all landunits and layers
    idx = 1
    if (isLayoutLeft) then
      ! LayoutLeft: landunits vary fastest
      do j = 1, NUM_SOIL_LAYERS
        do i = 1, numLandunits
          qflxTran(idx) = 0.0d0
          idx = idx + 1
        end do
      end do
    else
      ! LayoutRight: layers vary fastest
      do i = 1, numLandunits
        do j = 1, NUM_SOIL_LAYERS
          qflxTran(idx) = 0.0d0
          idx = idx + 1
        end do
      end do
    end if

    call UrbanSetTranspirationFlux(urban, c_loc(qflxTran), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    deallocate(qflxTran)

    ! Set fraction wet for impervious road (1D: per landunit)
    allocate(fwet(numLandunits))
    fwet(:) = FWET_INITIAL
    call UrbanSetFractionWetImperviousRoad(urban, c_loc(fwet), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
    deallocate(fwet)

    if (mpi_rank == 0) then
      write(*,*) 'Set hydrology boundary conditions:'
      write(*,*) '  Initial water table depth:', ZWT_INITIAL, 'm'
      write(*,*) '  Infiltration flux:', QFLX_INFL, 'mm/s'
      write(*,*) '  Transpiration flux: 0.0 mm/s (all layers)'
      write(*,*) '  Fraction wet (impervious road):', FWET_INITIAL
    end if
  end subroutine SetHydrologyBoundaryConditions

  subroutine SetCanyonAirStates(urban, numLandunits, mpi_rank)
    use iso_c_binding, only: c_int, c_double, c_loc
    implicit none
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank
    integer(c_int) :: status
    real(c_double), allocatable, target :: tempCanyonAir(:)
    real(c_double), allocatable, target :: qafCanyonAir(:)
    real(c_double), parameter :: TEMP_CANYON_AIR_INIT = 283.0_c_double
    real(c_double), parameter :: QAF_INIT = 1.e-4_c_double
    integer :: i

    allocate(tempCanyonAir(numLandunits))
    allocate(qafCanyonAir(numLandunits))

    do i = 1, numLandunits
      tempCanyonAir(i) = TEMP_CANYON_AIR_INIT
      qafCanyonAir(i) = QAF_INIT
    end do

    call UrbanSetCanyonAirTemperature(urban, c_loc(tempCanyonAir), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    call UrbanSetCanyonSpecificHumidity(urban, c_loc(qafCanyonAir), numLandunits, status)
    if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)

    deallocate(tempCanyonAir)
    deallocate(qafCanyonAir)

    if (mpi_rank == 0) then
      print *, 'Set canyon air properties:'
      print *, '  Canyon air temperature: ', TEMP_CANYON_AIR_INIT, ' K'
      print *, '  Canyon specific humidity: ', QAF_INIT, ' kg/kg'
    end if
  end subroutine SetCanyonAirStates

  subroutine SetUrbanParameters(urban, numLandunits, mpi_rank)
    type(UrbanType), intent(in) :: urban
    integer(c_int), intent(in) :: numLandunits
    integer, intent(in) :: mpi_rank

    call SetCanyonHwr(urban, numLandunits, mpi_rank)
    call SetFracPervRoadOfTotalRoad(urban, numLandunits, mpi_rank)
    call SetWtRoof(urban, numLandunits, mpi_rank)
    call SetHeightParameters(urban, numLandunits, mpi_rank)
    call SetBuildingTemperature(urban, numLandunits, mpi_rank)
    call SetSurfaceTemperatures(urban, numLandunits, mpi_rank)
    call SetLayerTemperatures(urban, numLandunits, mpi_rank)
    call SetCanyonAirStates(urban, numLandunits, mpi_rank)
    call SetAlbedo(urban, numLandunits, mpi_rank)
    call SetEmissivity(urban, numLandunits, mpi_rank)
    call SetNumberOfActiveLayersImperviousRoad(urban, numLandunits, mpi_rank)
    call SetThermalConductivity(urban, numLandunits, mpi_rank)
    call SetHeatCapacity(urban, numLandunits, mpi_rank)
    call SetSoilProperties(urban, numLandunits, mpi_rank)
    call SetAtmosphericForcing(urban, numLandunits, mpi_rank)
  end subroutine SetUrbanParameters

end program urbanxx_driver_f
