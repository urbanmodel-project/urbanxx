// C/C++ macros for Urban library
#ifndef URBAN_MACROS_H
#define URBAN_MACROS_H

#include <stdio.h>
#include <stdlib.h>

// Macro for error checking similar to PETSc's PetscCall
// For functions that return status via output parameter
#define UrbanCall(func, ierr_ptr)                                              \
  do {                                                                         \
    func;                                                                      \
    if (*(ierr_ptr) != URBAN_SUCCESS) {                                        \
      fprintf(stderr, "UrbanCall failed at %s:%d\n", __FILE__, __LINE__);      \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#endif // URBAN_MACROS_H
