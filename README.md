# stella

[![Check stella](https://github.com/stellaGK/stella/actions/workflows/check_stella.yml/badge.svg)](https://github.com/stellaGK/stella/actions/workflows/check_stella.yml)

$\texttt{stella}$ solves the gyrokinetic-Poisson system of equations in the local limit using an operator-split, implicit-explicit numerical scheme. It is capable of evolving electrostatic and electromagnetic fluctuations with fully kinetic electrons and an arbitrary number of ion species in general magnetic geometry, including stellarators.

A detailed $\texttt{stella}$ manual will be released in the coming months, which will be available on the GitHub repository, as well as on arXiv.

<br><br>

## Table of contents
- [Quick start](#quick-start)
- [Dependencies](#dependencies)
- [Installation and Compilation](#installation-and-compilation)
   * [Make](#make)
   * [CMake](#cmake)
- [Running stella](#running-stella)
- [Input files](#input-files)
- [Verification of stella output](#verification-of-stella-output)
   * [Set-up the automatic tests](#set-up-the-automatic-tests)
   * [Run numerical tests](#run-numerical-tests)
- [Acknowledgements](#acknowledgements)
   * [Gyrokinetic model & physics extensions](#gyrokinetic-model-and-physics-extensions)
   * [Magnetic geometries](#magnetic-geometries)
   * [Code infrastructure](#code-infrastructure)
   * [Diagnostics](#diagnostics)
   * [Small features](#small-features)
   * [Bug fixes](#bug-fixes)
      * [Major bug fixes](#major-bug-fixes)
      * [Minor bug fixes](#minor-bug-fixes)


<br><br>

## Quick start

On Ubuntu, the following commands download, compile and run $\texttt{stella}$:

```bash
git clone https://github.com/stellaGK/stella.git
cd stella
sudo apt install gfortran make libfftw3-dev libnetcdf-dev libnetcdff-dev netcdf-bin openmpi-bin libopenmpi-dev libblas-dev liblapack-dev
export STELLA_SYSTEM=gnu_ubuntu
make -j 4
mkdir -p ../stella_runs/first_run
cp DOCUMENTATION/example_input_files/miller_linear_CBC.in ../stella_runs/first_run/
cd ../stella_runs/first_run
mpirun -n 4 ../../stella/stella miller_linear_CBC.in
```

This runs a linear simulation of the Cyclone Base Case ITG mode (Miller geometry, adiabatic electrons) at four values of $k_y$, which takes about a minute on 4 cores. The growth rates and frequencies are written to `miller_linear_CBC.omega`, and the expected values, which agree with GENE, GS2, GX and CGYRO, are listed at the top of the input file; see [Running stella](#running-stella) for the other output files. On macOS, see [Make](#make) for the Homebrew packages and paths.


<br><br>

## Dependencies

In order to compile and run $\texttt{stella}$, install the following:
- A Fortran compiler, typically $\texttt{gfortran}$ (on macOS it comes with the $\texttt{gcc}$ package of Homebrew).
- An implementation of $\texttt{MPI}$, such as OpenMPI.
- The NetCDF Fortran library, to write the diagnostics.
- FFTW3, to compute the discrete Fourier transforms.
- BLAS and LAPACK, for the linear algebra routines.
- $\texttt{make}$ or $\texttt{cmake}$, to build the code.
- Python 3, for the automatic tests and the post-processing.

NetCDF, FFTW3 and LAPACK can be switched off at compile time (see [CMake](#cmake)), but nonlinear simulations need FFTW3 and the diagnostics need NetCDF, so in practice all of them are required.

On Ubuntu:

```bash
sudo apt install gfortran make libfftw3-dev libnetcdf-dev libnetcdff-dev netcdf-bin openmpi-bin libopenmpi-dev libblas-dev liblapack-dev
```

On macOS (with [Homebrew](https://brew.sh)):

```bash
brew install gcc make openmpi fftw netcdf-fortran
```

The submodules of $\texttt{stella}$ (in `EXTERNALS/`) are downloaded automatically during the first build. They can also be downloaded by hand with `git submodule update --init --recursive`.


<br><br>

## Installation and Compilation

There are two ways to build $\texttt{stella}$: using plain $\texttt{make}$ or $\texttt{CMake}$. Both create the executable `stella` in the main folder of the repository.


<br>

### Make

1. Set `STELLA_SYSTEM` to the system on which you are compiling, e.g. `gnu_ubuntu` or `macos`, which are tested automatically on GitHub. The supported systems are the files `COMPILATION/Makefiles/Makefile.<system>`.
2. Optionally, set the following environment variables to override the locations of the libraries in the `Makefile.<system>` file:
   - `FFTW_LIB_DIR`: directory containing libfftw3
   - `FFTW_INC_DIR`: directory containing fftw3.f
   - `NETCDF_LIB_DIR`: directory containing libnetcdff
   - `NETCDF_INC_DIR`: directory containing netcdf.inc
3. Run `make`, or `make -j 4` to compile on 4 cores.

On Ubuntu:

```bash
export STELLA_SYSTEM=gnu_ubuntu
make -j 4
```

On macOS with Homebrew:

```bash
export STELLA_SYSTEM=macos
export FFTW_LIB_DIR=$(brew --prefix fftw)/lib/
export FFTW_INC_DIR=$(brew --prefix)/include/
export NETCDF_LIB_DIR="$(brew --prefix)/lib -L$(brew --prefix netcdf-fortran)/lib/"
export NETCDF_INC_DIR=$(brew --prefix netcdf-fortran)/include/
make -j 4
```

Add the `export` lines to your `~/.bashrc` (or `~/.zshrc`) so that they are set in every new terminal. After changing the source code, `make` recompiles only what is needed; run `make depend` first when you add a new Fortran file or `use` statement. To clean the directory, the following commands exist:

```
make clean              # Removes compiled stella files, utils files and mini_libstell files
make clean-quick        # Only removes the compiled stella files, not the utils and mini_libstell files
make clean-submodules   # Clean + Remove git_version, neasyf and pFUnit folders
make distclean          # Clean + Remove stella executable + Invoke clean on pFUnit
```


<br>

### CMake

**Note**: If you have previously built $\texttt{stella}$ with plain `make` you _must_ run `make clean` before attempting to build with CMake, or the existing built objects will interfere with the CMake build.

Building $\texttt{stella}$ with CMake requires CMake >= 3.16. You can download the latest version from the [CMake website](https://cmake.org/download/), but it is often easier to install with `pip`:

```
pip install cmake
```

Building $\texttt{stella}$ is then a matter of first configuring the build:

```
cmake . -B COMPILATION/build_cmake
```

and then building proper:

```
cmake --build COMPILATION/build_cmake -j 4
```

The executable is written to `COMPILATION/build_cmake/stella`, with a symbolic link `stella` in the main folder. You may need to pass a few flags to the first `cmake` command to tell it where to find some dependencies:

```
cmake . -B COMPILATION/build_cmake \
  -DnetCDFFortran_ROOT=/path/to/netcdf/fortran \
  -DFFTW_ROOT=/path/to/fftw
```

There are a few build options:

- `STELLA_ENABLE_LAPACK`: Enable LAPACK (default: on)
- `STELLA_ENABLE_FFT`: Enable FFTs (default: on)
- `STELLA_ENABLE_NETCDF`: Enable NetCDF (default: on)
- `STELLA_ENABLE_DOUBLE`: Promotes the precision of reals and complex numbers to double (default: on)
- `STELLA_ENABLE_LOCAL_SPFUNC`: Enable local special functions (default: off)
- `STELLA_ENABLE_NAGLIB`: Use the NAG library (default: off)
- `STELLA_ENABLE_POSIX`: Enable POSIX functions for command line functionality (default: off)
- `STELLA_ENABLE_F200X`: Enable use of F2003/F2008 functionality (default: on)

You can turn these on or off with `-D<option name>=ON/OFF`. You can get a complete list of options by running the following in a build directory:

```
cmake -LH
```


<br><br>

## Running stella

$\texttt{stella}$ reads one input file, and is run with MPI from the folder that contains it:

```bash
mpirun -n 4 /path/to/stella/stella input.in
```

All output files start with the name of the input file:
- `input.out.nc`: all diagnostics (fields, fluxes, moments, geometry, grids) in NetCDF format.
- `input.out`: the time traces of $|\phi|^2$, $|A_\parallel|^2$ and $|B_\parallel|^2$ that are also written to the screen.
- `input.fluxes`: the particle, momentum and heat fluxes of each species versus time.
- `input.omega`: the frequency $\omega$ and growth rate $\gamma$ of every $(k_x, k_y)$ mode versus time (linear simulations).
- `input.geometry`: the geometric coefficients along the field line.
- `input.final_fields`: the electrostatic potential at the end of the simulation.
- `input.error`: warnings and errors; check this file when a simulation stops early.
- `input.input`: the input file with all variables, including the default values that were not set.

Good starting points for new simulations are the linear Cyclone Base Case example `DOCUMENTATION/example_input_files/miller_linear_CBC.in` and the default input file, which lists every input variable (see [Input files](#input-files)). The input files of the automatic tests in `AUTOMATIC_TESTS/numerical_tests/` show how to switch on other options, but they use very coarse grids and only a few time steps. Stellarator simulations need a $\texttt{VMEC}$ equilibrium (`vmec_filename`); the W7-X standard configuration and the NCSX equilibrium `wout_161s1.nc` are included in `DOCUMENTATION/VMEC_equilibria/`. The output can be plotted with the python package in `POST_PROCESSING/stellapy` (see its `README`) or read directly with `netCDF4` or `xarray` in python.


<br><br>

## Input files

Many of the namelists and variable names have been changed throughout the years to make them more intuitive. An up-to-date default $\texttt{stella}$ input file, which lists every input variable with its default value, can be found at:

```
STELLA_CODE/read_namelists_from_input_file/default_input_file.in
```

To convert input files of older $\texttt{stella}$ versions (v0.5 to v0.8) to the new format, run the following command in the folder containing the input files, where `/path/to/stella` is the folder of this repository:

```
python3 /path/to/stella/AUTOMATIC_TESTS/convert_input_files/convert_inputFile.py
```

Changes to the input variables since stella v1.0:
- `boundary_option = 'non-twisting'` (namelist `z_boundary_condition`) selects the non-twisting flux tube of Ball and Brunner (2021).
- `flip_flop = .false.` is the default (namelist `numerical_algorithms`); see the major bug fixes below.
- `zed_upwind = 0` is the default when `include_apar = .true.` (0.02 otherwise), and a warning is printed when `zed_upwind > 0` is set with $A_\parallel$.
- `no_j1l1 = .false.` is the default, so the $j = 1$, $l = 1$ term of the Fokker–Planck field-particle part is included.
- `cfac2` has been removed; `cfac` now scales all finite-Larmor-radius terms of the Fokker–Planck operator.
- `include_radial_variation` and `include_full_flux_annulus` (namelist `gyrokinetic_terms`) have been removed, since their values were always overwritten by `simulation_domain`: use `simulation_domain = 'multibox'` for radial variation and `simulation_domain = 'full_flux_annulus'` for a full flux surface.




<br><br>

## Verification of stella output

To ensure that the $\texttt{stella}$ code is functioning correctly, a suite of numerical tests has been implemented and is automatically run on every push to GitHub. These tests cover a wide range of routines and options within $\texttt{stella}$. The numerical tests are organised into nine categories. Tests 1, 2, 4, 5 and 6 use the electrostatic flux-tube version of $\texttt{stella}$, test 3 guards the numerical reproducibility of $\texttt{stella}$, while tests 7 and 8 target the full-flux-surface and electromagnetic versions, respectively. Finally, test 9 checks whether simulations can be successfully restarted.
- **Test 1**: Confirms that the $\texttt{stella}$ executable exists and that $\texttt{stella}$ runs correctly by checking that the executable **produces output files**.
- **Test 2**: Verifies that the **magnetic geometries** are implemented correctly. This includes geometries based on Miller parameters, $\texttt{VMEC}$ equilibria, and slab geometry.
- **Test 3**: Checks the **numerical reproducibility** of short nonlinear simulations, which switch on different numerical schemes (collisions, electromagnetic fields, adiabatic electrons, explicit algorithms, ...). The full state of the simulation (fields, moments, fluxes, frequency and distribution functions) is compared with a tight tolerance, it is verified that the results do not depend on the number of MPI processes or the data layouts, and that writing diagnostics does not change the time evolution. Run these tests before and after optimising the memory usage or speed of $\texttt{stella}$.
- **Test 4**: Checks each term of the **gyrokinetic equation** independently. It confirms that the electrostatic potential remains constant when no terms (nor collisions or dissipation) are included; validates the initialisation options for the distribution function, and verifies the correct evolution of the potential for each term. The $(k_x,k_y)$ grid options ("box" and "range") are also tested.
- **Test 5**: Tests the **parallel boundary conditions**: (1) standard twist-and-shift, (2) stellarator-symmetric, (3) periodic, and (4) zero boundary conditions. Additionally, multiple input flags are tested simultaneously. In the future, each input parameter would be tested individually.
- **Test 6**: Validates the **diagnostics**, including growth rates, fluxes, density and temperature, distribution function, and electrostatic potential.
- **Test 7**: Checks the electrostatic **full-flux-surface** version of $\texttt{stella}$. Geometric quantities are verified first, followed by each term of the gyrokinetic equation.
- **Test 8**: Tests the **electromagnetic** flux-tube version, verifying each term of the gyrokinetic equation. Short runs of a KBM and TAE instability check the corresponding code paths, while the full KBM and TAE benchmarks (against GS2) are slow, hence they are only run with `make numerical-tests-slow`.
- **Test 9**: Checks whether a simulation can successfully be **restarted**.


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
make numerical-tests-9
```
The slow physics benchmarks (KBM and TAE) are skipped by default, they can be run with:
```
make numerical-tests-slow
```




<br><br>

## Acknowledgements

Disclaimer: This section is incomplete.
It is the responsibility of the respective authors to list their contributions to the $\texttt{stella}$ code.

### Gyrokinetic model and physics extensions

- The electrostatic flux-tube version of the $\texttt{stella}$ gyrokinetic code has been developed by M. Barnes in 2018. The details of the code can be found in [[2019 - Barnes - Journal of Computational Physics](https://doi.org/10.1016/j.jcp.2019.01.025)].
- The "parallel nonlinearity" term has been added to the gyrokinetic equation by M. Barnes in 2018, as detailed in [[2018 - Barnes - Plasma Physics and Controlled Fusion](https://doi.org/10.1088/1361-6587/aaeb69)].
- A Linearised Fokker-Planck collision model has been implemented by A. von Boetticher in August 2021, as detailed in [[2024 - A. von Boetticher - Plasma Physics and Controlled Fusion](https://doi.org/10.1088/1361-6587/ad6c7c)].
- Stellarator-symmetric parallel boundary conditions have been implemented by A. González-Jerez in November 2021, based on [[2018 - Martin - Plasma Physics and Controlled Fusion](https://doi.org/10.1088/1361-6587/aad38a)].
- The radially global version of the $\texttt{stella}$ code has been developed by D. A. St-Onge in 2022, as detailed in [[2022 - D. A. St-Onge - Journal of Computational Physics](https://doi.org/10.1016/j.jcp.2022.111498)]. This paper used [stella release 0.5.1](https://github.com/stellaGK/stella/releases/tag/v0.5.1).
- Hyper dissipation in (vpa, z) was added by S. Stroteich in March 2024. These options can be controlled by the `D_zed` and `D_vpa` input variables in the `hyper_dissipation` namelist.
- The electromagnetic extension of the $\texttt{stella}$ code has been implemented by M. Hardman in July 2024, with the help of G. Acton, M. Barnes and R. Davies.
- The non-twisting flux tube has been implemented by H. Thienpondt in October 2026, based on [[2021 - Ball - Plasma Physics and Controlled Fusion](https://doi.org/10.1088/1361-6587/abf8f4)].
- The full-flux-surface version of $\texttt{stella}$ is being implemented by G. Acton in 2026 (see test 7).

### Magnetic geometries

- The $\texttt{VMEC}$ interface has been written by M. Landreman and has been incorporated into $\texttt{stella}$ by M. Barnes in 2018.
- The z-pinch magnetic geometry has been added by L. Podavini and M. Barnes in September 2024.

### Code infrastructure

- The CMake compilation has been set-up by P. Hill in September 2021.
- The construction of the response matrix now uses an LU decomposition, speeding up this routine drastically, this has been implemented by D. A. St-Onge in April 2022.
- The automatic testing infrastructure has been implemented by H. Thienpondt in July 2024.
- The $\texttt{stella}$ code has been reorganised and cleaned up by H. Thienpondt and G. Acton in October 2025.
- The shared-memory domains for the response matrix can now be parallelised over NUMA domains rather than over nodes. This is particularly useful on supercomputers with multiple sockets per node, where inter-socket communication can be relatively slow. To enable this feature, set `SPLIT_BY_NUMA = on` in the makefile and ensure that `--ntasks-per-socket` is specified in the sbatch script. This has been implemented by H. Thienpondt in December 2025.
- More automatic tests have been added by H. Thienpondt in October 2026 to test the full numerical reproducibility of the code. The electromagnetic tests (test 8) now compare with the new $\texttt{stella}$ data.

### Diagnostics

- Expanded the `gvmu` and `gvpas` diagnostics to print various diagnostics related to the distribution function, i.e. `|g|^2(t, mu, vpa, s)`; `|g|^2(t, z, vpa, s)`; `|g|^2(t, z, mu, s)`; `|g|^2(t, kx, ky, z, s)` and `|g|^2(t, z, vpa, mu, s)`, as well as their non-zonal components. Moreover, these diagnostics can be written for the perturbed distribution function $f$, the gyro-averaged distribution function $g$ and the non-adiabatic part of the distribution function $h$ (H. Thienpondt, July 2024).
- Added turbulent stress diagnostics for studying zonal flow physics (Y. Zhang, August 2026).
- The distribution diagnostics now sum over the processors once, instead of with 4 MPI reductions for every $(z, \mathrm{tube})$ point, which makes writing them much faster; the results are unchanged (H. Thienpondt, October 2026).

### Small features

- Print GitHub commit number, branch and date to the NetCDF file and header (H. Thienpondt and P. Hill, January 2022).
- Add more detailed code timers to the code (H. Thienpondt, January 2022).
- Improved the writing of the diagnostics to the NetCDF file using `neasyf` (P. Hill, May 2022).
- The CFL condition is now capable of increasing the time step as well, which is controlled by `cfl_cushion_upper`, `cfl_cushion_middle` and `cfl_cushion_lower` (H. Thienpondt, March 2023).
- Automatically stop linear simulations if `(gamma, omega)` has saturated. This is controlled by the `autostop` and `navg` input variables (H. Thienpondt, September 2024).
- Extended the `CMake` setup to make sure stella builds on newer compilers, since `ifort` has been deprecated (V. Seitz, October 2025).

### Bug fixes

#### Major bug fixes

These fixes change the results of $\texttt{stella}$ simulations, so results obtained before the fix may be affected.

- May 2026: Fixed the dependence of the implicit collision operators on the number of MPI processes. The response matrix is only solved on the processor that owns `is = 1` for each `(kx, ky, z, tube)` point, but the other processors kept using the fields from before the solve, and `phi` was summed over all processors. The solution is now collected from all processors, so the results no longer depend on the number of processes (G. Acton).
- October 2026: Fixed a bug in electromagnetic $\texttt{stella}$: the $B_\parallel$ contribution to the explicit $x$-component of the magnetic drift was gyro-averaged with $J_0$ instead of $J_1/a$. This affected all simulations with `include_bpar = .true.`, explicit magnetic drifts and nonzero $k_x$ (H. Thienpondt).
- October 2026: The explicit `rk3` time advance now uses the Shu–Osher SSPRK3 scheme. The previous scheme was only third order for linear terms; for the nonlinear term it was second order. Linear results do not change (H. Thienpondt).
- October 2026: Fixed the parallel boundary conditions: no twist-and-shift phase shift at the joins of unconnected chains; the zero-incoming boundary condition now acts in the same way with explicit as with implicit parallel streaming; the two copies of every twist-and-shift join are re-imposed during the explicit time advance; the response matrix applied the periodic end point to every segment of a periodic chain; `boundary_option = 'stellarator'` with Miller geometry now falls back to periodic boundary conditions with a warning (G. Acton).
- October 2026: Periodic modes now have the same $k_\perp^2$ at both copies of their end point $z = \pm z_\mathrm{max}$, which differed when the geometry is not up-down symmetric (e.g. stellarators). Before this fix, periodic boundary conditions did not work in stellarators: nonlinear W7-X simulations with kinetic electrons and `boundary_option = 'periodic'` always blew up, with the time step collapsing from 0.05 to $5\times10^{-5}$, while with this fix the same simulation runs normally with a constant time step (G. Acton).
- October 2026: The $v_\parallel$ hyper-dissipation is now applied; it was calculated but never added to the gyrokinetic equation (H. Thienpondt).
- October 2026: The explicit Dougherty energy-conservation term now gyro-averages the whole factor $(v_\parallel^2 + v_\perp^2 - 3/2)$ with $J_0$ instead of only $v_\parallel^2$ (H. Thienpondt).
- October 2026: The matrices of the implicit collision operators are now rebuilt when the CFL condition changes the time step; before, the collision frequency was effectively multiplied by `dt_old/dt_new` (R. Nies).
- October 2026: Fixed the Fokker–Planck field-particle part: the `psi` moments are zeroed before the sum over processors (NaN with implicit collisions), the field-particle part with several species now works on more than one processor, it stops with a clear message with explicit collisions (not implemented), and `cfac` now also scales its $J_1$ terms, so that `cfac = 0` gives the drift-kinetic operator (`cfac2` is removed). The perpendicular momentum moment of the implicit Dougherty operator now uses $v_\perp^2 J_1/a$ (H. Thienpondt).
- October 2026: Fixed implicit collisions with the implicit mirror term and `flip_flop = .true.`: on every other time step the mirror term acted on an outdated velocity-space distribution function, which made the results wrong even for vanishing collision frequency (H. Thienpondt).
- November 2026: `flip_flop = .false.` remains the default, and we have now investigated in detail why `flip_flop = .true.` should not be used. Alternating the order of the explicit and implicit parts of the operator splitting every time step should make the time advance second order accurate instead of first order, but with kinetic electrons `flip_flop = .true.` is unstable at time steps that `flip_flop = .false.` handles: the nearly time-centred implicit scheme flips the sign of the stiff electron components every time step, and reversing the order of the operator splitting every step makes the splitting errors add up. In nonlinear stellarator runs with kinetic electrons a grid-scale mode at the smallest $k_y$, at a maximum of $|B|$, burst the heat flux to up to 370 $Q_{gB}$ from the first time steps. The comments in `stella.f90` explain what is required to make the time advance second order (H. Thienpondt).
- November 2026: `zed_upwind = 0` is now the default when `include_apar = .true.` (0.02 otherwise), and stella prints a warning when `zed_upwind > 0` is set with $A_\parallel$. The upwinding of the parallel streaming acts as a numerical diffusion of the electrons along the field line; in a kinetic-ballooning mode the charge balance is the small remainder, of order $k_\perp^2\rho^2$, of much larger electron and ion responses, so this diffusion lowers the growth rate by an amount that depends only on the product of `zed_upwind` and $\Delta z/k_y^2$. With the old default the growth rate of the kinetic-ballooning mode fell to zero for $k_y \to 0$ in W7-X. Electromagnetic simulations that relied on the default `zed_upwind` should be checked at low $k_y$. The problem was reported in [issue #206](https://github.com/stellaGK/stella/issues/206) by [Toby Adkins](https://github.com/physycola), who found that the eigenfunctions and the electron/ion flux partition of an electromagnetic CBC disagreed with GS2, GX and CGYRO at low $k_y$ when $A_\parallel$ is included with implicit parallel streaming; [Yujia Zhang](https://github.com/Scottyujia) identified the upwinding in $z$ as the cause, showing that `zed_upwind = 0` alone makes the eigenfunctions agree with GS2 (Toby Adkins, Yujia Zhang and H. Thienpondt).
- November 2026: Removed the density error of the implicit Fokker–Planck test-particle solve on the default (non-uniform) $\mu$ grid. The solve changed the density by about $10^{-5}$ per time step, which quasineutrality amplifies by $1/(k_\perp\rho_i)^2$ with kinetic electrons, so that collisions alone were unstable ($|\phi|^2$ grew by a factor 90 by $t = 2$, and $\gamma = 0.30$ instead of 0.12 at $\nu_{ii} = 0.3$ in the CBC benchmark with kinetic electrons). The error is now removed after the solve with a Maxwellian without momentum and energy moments, keeping the density change of the classical (gyro-)diffusion; this is applied with kinetic electrons unless `density_conservation = .true.` (H. Thienpondt).
- November 2026: Removed the density moment of the energy-exchange ($l = 0$, $j \ge 1$) terms of the Fokker–Planck field-particle part. For electrons colliding with ions, these terms depend on the ion speed, which lies inside the first cells of the electron velocity grid, so on the grid their density moment was as large as their energy moment: the energy exchange created electron density, and with collisions only the potential grew by a factor 6.5 by $t = 0.5$. The energy exchange itself is unchanged (H. Thienpondt).
- November 2026: The $j = 1$, $l = 1$ term of the Fokker–Planck field-particle part is now included by default (`no_j1l1 = .false.`). It was switched off by default without a documented reason; without the finite-Larmor-radius terms (`cfac = 0`) it hardly matters, but with them (`cfac = 1`) the growth rates were too low at higher collisionality compared with GS2 and CGYRO, and with it $\texttt{stella}$ agrees with both (H. Thienpondt).
- November 2026: The $v_\parallel$ advection speed of the parallel nonlinearity (`include_parallel_nonlinearity = .true.`) now keeps its dependence on $v_\parallel$. Before, the speed at the first grid point ($v_\parallel = -v_{\parallel,\max}$) was used for every $v_\parallel$, although its magnetic-drift and curvature parts are proportional to $v_\parallel$; each point is now upwinded with the sign of its own speed, and the CFL condition uses the largest speed along $v_\parallel$ (H. Thienpondt).
- November 2026: The parallel nonlinearity now calculates $\partial\langle\phi\rangle/\partial z$ with the same centred stencil for both signs of $v_\parallel$. It used the routine of the parallel streaming term, which at the ends of every chain of connected $k_x$ upwinds with the sign of $v_\parallel$ (since October 2026 with the row $\langle\phi\rangle/\Delta z$ at the incoming end), so the parallel acceleration differed between $v_\parallel < 0$ and $v_\parallel > 0$ there. In the nonlinear CBC at $a/L_{T_i} = 3.6$ and $\rho_* = 0.05$, this part alone drove a toroidal momentum flux of $+0.22\,\Pi_{gB}$, raised the heat flux by 14% and halved the time step, while in GENE it gives no momentum flux (H. Thienpondt).

#### Minor bug fixes

These fixes do not change the results of earlier simulations (e.g. crashes, debug builds, timers and diagnostics).

- September 2024: Fixed masking of the `omega` data in the NetCDF output (H. Thienpondt).
- December 2025: Fixed an integer-overflow bug in the MPI shared-memory implementation for the response matrix. The shared-window size and associated pointers were previously defined using 32-bit integers, limiting the usable shared memory to approximately 8–9 GB and leading to segmentation faults for larger problems. These quantities are now defined using 64-bit integers, removing this artificial limitation and enabling reliable handling of response matrices of at least ~40 GB (and likely larger). The maximum supported size is now determined by the available node RAM rather than by internal stella integer limits (H. Thienpondt).
- October 2026: Fixed the code timers `calculate_phi_and_bpar` and the `int_dv_g` part of the field solve in the kxkyz layout, which were started but never stopped, so that they also measured the time between two calls and reported more time than the whole field solve (H. Thienpondt).
- October 2026: Writing all diagnostics (`write_all = .true.`) no longer crashes simulations without radial variation, since the radial moments and fluxes are only written with radial variation; radial variation with VMEC geometry, which used zeros for the radial derivatives of the geometry, now stops with a clear error message (H. Thienpondt).
- October 2026: Fixed the errors found by running the automatic tests with a debug build: uninitialised `initialised_*` flags of the grids and arrays, the electrostatic implicit mirror term reading `apar` beyond its size, `jtwist`, `jtwistfac` and `ikx_twist_shift` uninitialised for the range grid (`ikx_twist_shift` is now also declared as an integer instead of a real), geometry quantities that were not broadcast to all processors, and a partly initialised collision frequency (H. Thienpondt).
- October 2026: The input file converter no longer drops `include_apar`, `include_bpar`, `explicit_option` and the diagnostics flags of old input files, and the default input file no longer makes $\texttt{stella}$ abort (H. Thienpondt).
- October 2026: Fixed the gyro-averages of the flow shear and stress diagnostics, which indexed `apar` and `bpar` in electrostatic simulations (H. Thienpondt).
- October 2026: With `flip_flop = .true.`, $(\omega, \gamma)$ are now calculated over two time steps; before, they alternated between two wrong values every time step, which also prevented `autostop` (H. Thienpondt).
- October 2026: Fixed a division by zero for the $k_x = k_y = 0$ mode of the modified Boltzmann response, which made debug builds stop in nonlinear simulations with adiabatic electrons (H. Thienpondt).
- October 2026: Fixed a first-order error in time at the parallel boundaries: the end points of periodic (e.g. zonal) modes are now made consistent on every Runge–Kutta stage, and the parallel boundary conditions are imposed on the initial distribution function. Before, the explicit time advance was first order in $\Delta t$ for zero and linked boundary conditions, and in linear simulations with periodic boundary conditions; the results only change at the level of the time-stepping error (H. Thienpondt).
- October 2026: The CFL condition of explicit Dougherty collisions now bounds the combined stiffness of the $v_\parallel$ operator, the $\mu$ operator and the classical (gyro-)diffusion; before, it took the minimum of separate limits for the $v_\parallel$ and $\mu$ operators, and with long field lines or many $\mu$ points (e.g. `nperiod = 4`, `nmu = 24`) the time step it chose made the explicit collisions blow up (H. Thienpondt).
- November 2026: Linear simulations with implicit mirror and parallel streaming terms and with both magnetic drifts switched off (`include_xdrift = include_ydrift = .false.`) were flagged as fully implicit, so the explicit time advance was skipped. Since the drive term is advanced in the explicit step (unless `drifts_implicit = .true.`), switching off the drifts also silently removed the drive. The explicit step is now kept whenever the drive is included (H. Thienpondt).
- November 2026: Fixed the toroidal angular momentum flux (`vflux`) of Miller geometry: its parallel part used $R_0 B_\mathrm{tor}/B$ with the constant major radius $R_0$ (`rmaj`) instead of $R(z) B_\mathrm{tor}/B = I/B$, which made it too small on the outboard and too large on the inboard side ($R_0/R(z)$, about 15% in the CBC). The $A_\parallel$ and $B_\parallel$ contributions to the momentum flux now also use $R^2\,\hat{b}\cdot\nabla\zeta$; before, they multiplied by `rmaj*btor/bmag`, which was wrong in the same way for Miller geometry and meaningless for VMEC geometry, where `btor` is not defined (set to $-1000$). The radial-variation fluxes use $I = R(z) B_\mathrm{tor}$ in the same way. Only the momentum flux diagnostic changes (H. Thienpondt).
- November 2026: Fixed the double Maxwellian factor of the drive term with neoclassical terms: with `include_neoclassical_terms = .true.`, `wstar` was multiplied by $F_0$ a second time (introduced in v1.0), so the gradient drive was proportional to $F_0^2$ and `dfneo_drho` picked up a spurious factor $F_0$. Neoclassical terms are currently disabled at input, so no v1.0 results are affected (H. Thienpondt).
