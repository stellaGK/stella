# stella

[![Check stella](https://github.com/stellaGK/stella/actions/workflows/check_stella.yml/badge.svg)](https://github.com/stellaGK/stella/actions/workflows/check_stella.yml)

$\texttt{stella}$ solves the gyrokinetic-Poisson system of equations in the local limit
using an operator-split, implicit-explicit numerical scheme. It is capable of
evolving electrostatic fluctuations with fully kinetic electrons and an
arbitrary number of ion species in general magnetic geometry, including
stellarators.

A detailed $\texttt{stella}$ manual will be released in the coming months, which will be available on the GitHub repository, as well as on arXiv.

<br><br>

## Table of contents
- [Dependencies](#dependencies)
- [Installation and Compilation](#installation-and-compilation)
   * [CMake](#cmake)
   * [Make](#make)
- [New input variables](#new-input-variables)
- [Verification of stella output](#verification-of-stella-output)
   * [Set-up the automatic tests](#set-up-the-automatic-tests)
   * [Run numerical tests](#run-numerical-tests)
- [Acknowledgments](#acknowledgments)
   * [Gyrokinetic model & physics extensions](#gyrokinetic-model-and-physics-extensions)
   * [Magnetic geometries](#magnetic-geometries)
   * [Code infrastructure](#code-infrastructure)
   * [Diagnostics](#diagnostics)
   * [Small features](#small-features)
   * [Bug fixes](#bug-fixes)


<br><br>

## Dependencies

In order to compile and run $\texttt{stella}$, install the following libraries:
- A Fortran compiler, typically $\texttt{gfortan}$ on Ubuntu, or $\texttt{gcc}$ on macOS.
- NetCDF libraries in order to write the diagnostics.
- Libraries for linear algebra routines, such as $\texttt{blas}$ and $\texttt{lapack}$.
- An implementation of $\texttt{MPI}$ for parallel runs.
- FFTW for computing discrete Fourier transformations.
- The $\texttt{make}$ or $\texttt{cmake}$ library, for building the code.
- Python utilities for the automatic tests.

Note that the following libraries are optional:
- NetCDF Fortran
- FFTW3
- LAPACK


<br><br>

## Installation and Compilation

There are two ways to build stella: using $\texttt{CMake}$ or plain $\texttt{make}$.


<br>

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


<br>

### Make

The other build system uses plain `make`:

1. Set `STELLA_SYSTEM='system'`, with `system` replaced by the appropriate system on
   which you are running. See the `COMPILATION/Makefiles/` directory for a list of supported
   systems.
2. Optionally, set the following environment variables to override the locations
   in the `STELLA_SYSTEM` Makefile:
   - `FFTW_LIB_DIR`: directory containing libfftw3
   - `FFTW_INC_DIR`: directory including fftw3.f
   - `NETCDF_LIB_DIR`: directory containing libnetcdff
   - `NETCDF_INC_DIR`: directory including netcdf.inc
4. Set the environment variable `MAKEFLAGS=-IMakefiles`, or set `-IMakefiles`
   when you run `make`
5. Run `make`

For example, to compile on Ubuntu:

```bash
# Using bash:
export STELLA_SYSTEM=gnu_ubuntu
export MAKEFLAGS=-IMakefiles
make

# Or in one line:
make -IMakefiles STELLA_SYSTEM=gnu_ubuntu
```
In summary the exports of `STELLA_SYSTEM` and `MAKEFLAGS` are set, compiling stella is achieved through:
```
make clean
make depend
make
```
To `clean` the directory, the following commands exist:
```
make clean              # Removes compiled stella files, utils files and mini_libstell files
make clean-quick        # Only removes the compiled stella files, not the utils and mini_libstell files
make clean-submodules   # Clean + Remove git_version, neasyf and pFUnit folders
make distclean          # Clean + Remove stelle executable + Invoke clean on pFUnit
```




<br><br>

## New input variables

Many of the namelists and variable names have been changed throughout the years to make them more intuitive. An up-to-date default $\texttt{stella}$ input file can be found at:

```
STELLA_CODE/read_namelists_from_input_file/default_input_file.in
```

To convert old $\texttt{stella}$ input files to the new format, run the following command in the folder containing the input file:

```
python3 $STELLA/AUTOMATIC_TESTS/convert_input_files/convert_inputFile.py
```




<br><br>

## Verification of stella output

To ensure that the $\texttt{stella}$ code is functioning correctly, a suite of numerical tests has been implemented and is automatically run on every push to GitHub. These tests cover a wide range of routines and options within $\texttt{stella}$. The numerical tests are organized into eight categories. Tests 1-5 use the electrostatic flux-tube version of $\texttt{stella}$, while tests 6 and 7 target the full-flux-surface and electromagnetic versions, respectively. Finally, test 8 checks whether simulations can be succesfully restarted.
- **Test 1**: Confirms that the $\texttt{stella}$ executable exists and that $\texttt{stella}$ runs correctly by checking that the executable **produces output files**.
- **Test 2**: Verifies that the **magnetic geometries** are implemented correctly. This includes geometries based on Miller parameters, $\texttt{VMEC}$ equilibria, and slab geometry.
- **Test 3**: Checks each term of the **gyrokinetic equation** independently. It confirms that the electrostatic potential remains constant when no terms (nor collisions or dissipation) are included; validates the initialization options for the distribution function, and verifies the correct evolution of the potential for each term. The $(k_x,k_y)$ grid options ("box" and "range") are also tested.
- **Test 4**: Tests the **parallel boundary conditions**: (1) standard twist-and-shift, (2) stellarator-symmetric, (3) periodic, and (4) zero boundary conditions. Additionally, multiple input flags are tested simultaneously. In the future, each input parameter would be tested individually.
- **Test 5**: Validates the **diagnostics**, including growth rates, fluxes, density and temperature, distribution function, and electrostatic potential.
- **Test 6**: Checks the electrostatic **full-flux-surface** version of $\texttt{stella}$. Geometric quantities are verified first, followed by each term of the gyrokinetic equation.
- **Test 7**: Tests the **electromagnetic** flux-tube version, verifying each term of the gyrokinetic equation and simulating a KBM and TAE instability.
- **Test 8**: Checks whether a simulation can succesfully be **restarted**.


<br>

### Set-up the automatic tests
The first time you want to run these tests, you need to install the python virtual environment:
```
 make create-test-virtualenv
```
Next, activate the virtual environment:
```
source AUTOMATIC_TESTS/venv/bin/activate
```


<br>

### Run numerical tests
Run the automated python tests, which tests the numerical output of stella:
```
make numerical-tests
```
When debugging the code, the tests can be run in chunks instead:
```
make numerical-tests-1
make numerical-tests-2
make numerical-tests-3
make numerical-tests-4
make numerical-tests-5
make numerical-tests-6
make numerical-tests-7
make numerical-tests-8
```




<br><br>

## Acknowledgments

Disclaimer: This section is incomplete.
It is the responsibility of the respective authors to list their contributions to the $\texttt{stella}$ code.

### Gyrokinetic model and physics extensions

- The electrostatic flux-tube version of the $\texttt{stella}$ gyrokinetic code has been developed by M. Barnes in 2018. The details of the code can be found in [[2019 - Barnes - Journal of Computational Physics](https://doi.org/10.1016/j.jcp.2019.01.025)].
- The "parallel nonlinearity" term has been added to the gyrokinetic equation by M. Barnes in 2018, as detailed in [[2018 - Barnes - Plasma Physics and Controlled Fusion](https://doi.org/10.1088/1361-6587/aaeb69)].
- A Linearised Fokker-Planck collision model has been implemented by A. von Boetticher in August 2021, as detailed in [[2024 - A. von Boetticher - Plasma Physics and Controlled Fusion](https://doi.org/10.1088/1361-6587/ad6c7c)].
- Stellarator-symmetric parallel boundary conditions have been implemented by A. González-Jerez in November 2021, based on [[2018 - Martin - Plasma Physics and Controlled Fusion](https://doi.org/10.1088/1361-6587/aad38a)].
- The radially global version of the $\texttt{stella}$ code has been developed by D. A. St-Onge in 2022, as detailed in [[2022 - D. A. St-Onge - Journal of Computational Physics](https://doi.org/10.1016/j.jcp.2022.111498)]. This paper used [stella release 0.5.1](https://github.com/stellaGK/stella/releases/tag/v0.5.1).
- Hyper dissipation was added by S. Stroteich in March 2024.
- The electromagnetic extension of the $\texttt{stella}$ code has been implemented by M. Hardman in July 2024, with the help of G. Acton, M. Barnes and R. Davies.
- The full-flux-surface version of $\texttt{stella}$ will be implemented by G. Acton in 2026.

### Magnetic geometries

- The $\texttt{VMEC}$ interface has been written by M. Landreman and has been incorporated into $\texttt{stella}$ by M. Barnes in 2018.
- The z-pinch magnetic geometry has been added by L. Podavini and M. Barnes in September 2024.

### Code infrastructure

- The CMake compilation has been set-up by P. Hill in September 2021.
- The construction of the response matrix now uses an LU decomposition, speeding up this routine drastically, this has been implemented by D. A. St-Onge in April 2022.
- The automatic testing infrastructure has been implemented by H. Thienpondt in July 2024.
- The $\texttt{stella}$ code has been reorganized and cleaned up by H. Thienpondt and G. Acton in October 2025.
- The shared-memory domains for the response matrix can now be parallelized over NUMA domains rather than over nodes. This is particularly useful on supercomputers with multiple sockets per node, where inter-socket communication can be relatively slow. To enable this feature, set `SPLIT_BY_NUMA = on` in the makefile and ensure that `--ntasks-per-socket` is specified in the sbatch script. This has been implemented by H. Thienpondt in December 2025.

### Diagnostics

- Expanded the `gvmu` and `gvpas` diagnostics to print various diagnostics related to the distribution function, i.e. `|g|^2(t, mu, vpa, s)`; `|g|^2(t, z, vpa, s)`; `|g|^2(t, z, mu, s)`; `|g|^2(t, kx, ky, z, s)` and `|g|^2(t, z, vpa, mu, s)`, as well as their non-zonal components. Moreover, these diagnostics can be written for the perturbed distribution function $f$, the gyro-averaged distribution function $g$ and the non-adiabatic part of the distribution function $h$ (H. Thienpondt, July 2024).

### Small features

- Print GitHub commit number, brach and date to the NetCDF file and header (H. Thienpondt and P. Hill, January 2022).
- Add more detailed code timers to the code (H. Thienpondt, January 2022).
- Improved the writing of the diagnostics to the NetCDF file using `neasyf` (P. Hill, May 2022).
- The CFL condition is now capable of increasing the time step as well, which is controlled by `cfl_cushion_upper`, `cfl_cushion_middle` and `cfl_cushion_lower` (H. Thienpondt, March 2023).
- Automatically stop linear simulations if `(gamma, omega)` has saturated. This is controlled by the `autostop` and `navg` input variables (H. Thienpondt, September 2024).
- Extended the `CMake` setup to make sure stella builds on newer compilers, since `ifort` has been deprecated (V. Seitz, October 2025).

### Bug fixes

- September 2024: Fixed masking of the `omega` data in the NetCDF output (H. Thienpondt).
- December 2025: Fixed an integer-overflow bug in the MPI shared-memory implementation for the response matrix. The shared-window size and associated pointers were previously defined using 32-bit integers, limiting the usable shared memory to approximately 8–9 GB and leading to segmentation faults for larger problems. These quantities are now defined using 64-bit integers, removing this artificial limitation and enabling reliable handling of response matrices of at least ~40 GB (and likely larger). The maximum supported size is now determined by the available node RAM rather than by internal stella integer limits (H. Thienpondt).

