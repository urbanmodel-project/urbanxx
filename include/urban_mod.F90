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

  ! Interface declarations for C functions
  ! Add your C API function interfaces here as needed
  ! Example:
  ! interface
  !   function UrbanCreate(urban) bind(C, name="UrbanCreate")
  !     import :: c_ptr, c_int
  !     type(c_ptr) :: urban
  !     integer(c_int) :: UrbanCreate
  !   end function UrbanCreate
  ! end interface

end module urban_mod
