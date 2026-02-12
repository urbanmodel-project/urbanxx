# URBANXX

[![Build Status](https://github.com/urbanmodel-project/urbanxx/actions/workflows/auto_test.yml/badge.svg)](https://github.com/urbanmodel-project/urbanxx/actions)
[![Code Coverage](https://codecov.io/github/urbanmodel-project/urbanxx/branch/main/graph/badge.svg)](https://codecov.io/github/urbanmodel-project/urbanxx)

## Prerequisites

- CMake
- gfortran
- MPI library (e.g., OpenMPI)
- C++ compiler with OpenMP support
- Git (for submodules)

## Building

Clone the repository with submodules (required for Kokkos):

```bash
git clone --recursive https://github.com/urbanmodel-project/urbanxx.git
cd urbanxx
```

Configure the build:

```bash
cmake -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-DURBANXX_ENABLE_OPENMP" \
  -DKokkos_ENABLE_OPENMP=ON \
  -DKokkos_ARCH_NATIVE=ON \
  -DKokkos_ENABLE_DEPRECATED_CODE_4=OFF
```

Build the project:

```bash
cmake --build build
```

## Running the Model

URBANXX provides two driver executables:

**C++ driver:**
```bash
./build/driver/urbanxx_driver
```

**Fortran driver:**
```bash
./build/driver/urbanxx_driver_f
```

## Testing

Run the test suite:

```bash
cd build
ctest --output-on-failure
```

## Build Configuration Options

- `CMAKE_BUILD_TYPE`: Build type (`Release`, `Debug`, or `RelWithDebInfo`)
- `URBANXX_ENABLE_OPENMP`: Enable OpenMP parallelization
- `Kokkos_ENABLE_OPENMP`: Enable Kokkos OpenMP backend
- `Kokkos_ARCH_NATIVE`: Optimize for native architecture
- `Kokkos_ENABLE_DEPRECATED_CODE_4`: Disable deprecated Kokkos 4.x code

## Project Structure

- `src/`: Core C++ source files
- `driver/`: Driver applications (C++ and Fortran)
- `tests/`: Unit tests
- `include/`: Public headers
- `externals/`: External dependencies (Kokkos)
- `doc/` : Documentation (work in progress)
