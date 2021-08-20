# stella

![Github Actions badge](https://github.com/stellaGK/stella/actions/workflows/tests.yml/badge.svg)

stella solves the gyrokinetic-Poisson system of equations in the local limit
using an operator-split, implicit-explicit numerical scheme. It is capable of
evolving electrostatic fluctuations with fully kinetic electrons and an
arbitrary number of ion species in general magnetic geometry, including
stellarators.

## Dependencies

stella has several optional dependencies:

- MPI
- netCDF Fortran
- FFTW3

## Installation

1. Set `GK_SYSTEM='system'`, with `system` replaced by the appropriate system on
   which you are running. See the `Makefiles` directory for a list of supported
   systems.
2. Optionally, set the following environment variables to override the locations
   in the `GK_SYSTEM` Makefile:
   - `FFTW_LIB_DIR`: directory containing libfftw3
   - `FFTW_INC_DIR`: directory including fftw3.f
   - `NETCDF_LIB_DIR`: directory containing libnetcdff
   - `NETCDF_INC_DIR`: directory including netcdf.inc
4. Set the environment variable `MAKEFLAGS=-IMakefiles`, or set `-IMakefiles`
   when you run `make`
5. Run `make`

For example, to compile on Ubuntu:

```bash
# using bash:
export GK_SYSTEM=gnu_ubuntu
export MAKEFLAGS=-IMakefiles
make

# or in one line:
make -IMakefiles GK_SYSTEM=gnu_ubuntu
```
