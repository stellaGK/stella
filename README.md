# stella

[![Check stella](https://github.com/stellaGK/stella/actions/workflows/check_stella.yml/badge.svg)](https://github.com/stellaGK/stella/actions/workflows/check_stella.yml)

$\texttt{stella}$ solves the gyrokinetic-Poisson system of equations in the local limit
using an operator-split, implicit-explicit numerical scheme. It is capable of
evolving electrostatic fluctuations with fully kinetic electrons and an
arbitrary number of ion species in general magnetic geometry, including
stellarators.

Automatic documentation is created based on the `stella` code, however this is a work in progress.

https://stellagk.github.io/stella/

<br>

## Table of contents
  * [Dependencies](#dependencies)
  * [Installation and Compilation](#installation-and-compilation)
    + [CMake](#cmake)
    + [Make](#make)
  * [Verification of stella output](#verification-of-stella-output)
    + [Set-up](#set-up)
    + [Numerical tests](#numerical-tests)


<br>

## Dependencies

stella requires MPI, and has several optional dependencies:

- netCDF Fortran
- FFTW3
- LAPACK


<br>

## Installation and Compilation

There are two ways to build stella: with CMake (experimental) or with
plain `make`.


### CMake

**Note**: If you have previously built stella with plain `make` you
_must_ run `make clean` before attempting to build with CMake, or the
existing built objects will interfere with the CMake build.

Building stella with CMake requires CMake >= 3.16. You can download
the latest version from the [CMake
website](https://cmake.org/download/), but it is often easier to
install with `pip`:

```
pip install cmake
```

Building stella is then a matter of first configuring the build:

```
cmake . -B COMPILATION/build_cmake
```

and then building proper:

```
cmake --build COMPILATION/build_cmake
```

You may need to pass a few flags to the first `cmake` command to tell
it where to find some dependencies:

```
cmake . -B build \
  -DnetCDFFortran_ROOT=/path/to/netcdf/fortran
  -DFFTW_ROOT=/path/to/fftw
```

There are a few build options:

- `STELLA_ENABLE_LAPACK`: Enable LAPACK (default: on)
- `STELLA_ENABLE_FFT`: Enable FFTs (default: on)
- `STELLA_ENABLE_NETCDF`: Enable NetCDF (default: on)
- `STELLA_ENABLE_DOUBLE`: Promotes precisions of real and complex to double
  (default: on)
- `STELLA_ENABLE_LOCAL_SPFUNC`: Enable local special functions" (default: off)
- `STELLA_ENABLE_NAGLIB`: Use the NAG library (default: off)
- `STELLA_ENABLE_POSIX`: Enable POSIX functions for command line functionality
  (default: off)
- `STELLA_ENABLE_F200X`: Enable use of F2003/F2008 functionality (default: on)

You can turn these on or off with `-D<option name>=ON/OFF`. You can
get a complete list of options by running the following in a build
directory:

```
cmake -LH
```


### Make

The other build system uses plain `make`:

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
If the exports of `GK_SYSTEM` and `MAKEFLAGS` are set, compiling stella is achieved through:
```
make
```
To `clean` the directory, the following commands exist:
```
make clean              # removes compiled stella files, utils files and mini_libstell files
make clean-quick        # only removes the compiled stella files, not the utils and mini_libstell files
make clean-submodules   # clean + remove git_version, neasyf and pFUnit folders
make distclean          # clean + remove stelle executable + invoke clean on pFUnit
```

<br>

## Verification of stella output

Automated python tests are available to test the output of stella. 

### Set-up
The first time you want to run these tests, you need to install the python virtual environment:
```
make create-test-virtualenv
```
Next, activate the virtual environment:
```
source AUTOMATIC_TESTS/venv/bin/activate
```

### Numerical tests
Run the automated python tests, which tests the numerical output of stella:
```
make numerical-tests
```


