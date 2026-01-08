! Fortran preprocessor macros for Urban library
! This file should be included with #include "finclude/UrbanMacros.h"

#define CallA(func) \
  call func ;\
  if (status /= URBAN_SUCCESS) call UrbanError(mpi_rank, __LINE__, status)
