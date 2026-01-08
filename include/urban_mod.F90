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
