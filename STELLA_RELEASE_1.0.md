
The $\texttt{stella}$ code has been extensively cleaned and reorganized by **H. Thienpondt** and **G. Acton**. Numerous files have been moved, renamed, split, and streamlined. Additional comments and headers have been introduced throughout the code to enhance readability and clarity. Each folder now includes a dedicated README file providing an overview of its contents. Furthermore, new automated tests have been added to strengthen reliability.

## Restructuring of the code in July 2024

The cleanup of the $\texttt{stella}$ code began in July 2024. To preserve the pre-cleanup version, it was released as `stella release v0.6`. Subsequently, the electromagnetic extension was implemented and released as `stella release v0.7`, which introduced the following restructurings and improvements:

- The code has been organized into subfolders:
   - AUTOMATIC_TESTS
   - COMPILATION
   - DOCUMENTATION
   - EXTERNALS
   - POST_PROCESSING
   - STELLA_CODE
      - calculations
      - diagnostics
      - dissipation
      - fields
      - geometry
      - grids
      - gyrokinetic_terms
      - neoclassical
      - parameters
      - radial_variation
      - utils
      - stella.f90

- Most importantly, the main $\texttt{stella}$ code has been centralized in the "STELLA_CODE" folder and further organized into subfolders. This structure makes it easier for users to search for keywords, to locate key sections, and to focus on the core scripts.

- The compilation process using `make` has been fully revamped and moved to the "COMPILATION" folder, preventing compiled object files from cluttering the main code.

- An automatic testing infrastructure has been implemented to ensure code reliability.

This initial cleaned version of $\texttt{stella}$ has been saved and released as `stella release v0.8`.


## Restructuring of the code in September 2025

The cleanup of the $\texttt{stella}$ code was resumed in the summer of 2025. Key updates include:

- The input variables have been fully reorganized, with most variables and namelists renamed. 

- The main stella code within "STELLA_CODE" has been split up further into the following subfolders:
   - arrays
   - calculations
   - diagnostics
   - dissipation
   - field_equations
   - geometry
   - grids
   - gyrokinetic_equation
   - neoclassical
   - parallelisation
   - parameters
   - radial_variation
   - read_namelists_from_input_file
   - utils
   - init_stella.f90
   - stella.f90
   
- Many variable and routine names have been updated to improve clarity and readability.

- Additional `abort` statements have been added throughout the code to prevent incorrect usage.

- Added an `Acknowledgments` section to the main `README.md` file to keep track of contributions to the $\texttt{stella}$ code. This list is incomplete and it is the responsibility of the various authors to add their contributions.

- Various bugs have been fixed.


### New input variables

The default $\texttt{stella}$ input file can be found at:

```
STELLA_CODE/read_namelists_from_input_file/default_input_file.in
```

To convert old $\texttt{stella}$ input files to the new format, run the following command in the folder containing the input file:

```
python3 $STELLA/AUTOMATIC_TESTS/convert_input_files/convert_inputFile.py
```


### Restructuring of the $\texttt{stella}$ scripts

The main $\texttt{stella}$ script has been cleaned up and divided into two separate files:: 
- stella.f90
- init_stella.f90

Added a new directory called “read_namelists_from_input_file” which reads in the input file. All namelists are read here now.

The script that evolves the distribution function has been split up and renamed:
- The time_advance.f90 script has been replaced by:
   - gyrokinetic_equation/gyrokinetic_equation_initialisation.f90
   - gyrokinetic_equation/gyrokinetic_equation_explicit.f90
   - gyrokinetic_equation/gyrokinetic_equation_implicit.f90
   - gyrokinetic_equation/gk_magnetic_drift.f90
   - gyrokinetic_equation/gk_drive.f90
   - gyrokinetic_equation/gk_nonlinearity.f90
   
- Several routines from time_advance.f90 have been separated into:
   - calculations/calculations_checksim.f90
   - calculations/calculations_add_explicit_terms.f90
   - calculations/calculations_timestep.f90
   
The following scripts have been renamed:
- ffs_solve.f90 → gyrokinetic_equation/gk_ffs_solve.f90
- flow_shear.f90 → gyrokinetic_equation/gk_flow_shear.f90
- parallel_streaming.f90 → gyrokinetic_equation/gk_parallel_streaming.f90
- implicit_solve.f90 → gyrokinetic_equation/gk_implicit_terms.f90

The gyro_averages.f90 script has been split up: 
- calculations_gyro_averages.f90
- arrays_gyro_averages.f90

The fields scripts have been renamed to field_equations, and the following files have been moved out of "fields":
- fields/dist_fn.f90 → arrays/initialise_arrays.f90
- fields/init_g.f90 → arrays/initialise_distribution_function.f90
 
The vpamu_grids.f90 script has been renamed and split into:
- calculations/calculations_velocity_integrals.f90
- grids/grids_velocity.f90

Other moves and name changes:
- grids/stella_time.f90 → grids/grids_time.f90
- grids/arrays_dist_fn.f90 → arrays/arrays_distribution_function.f90
- grids/arrays_fields.f90 → arrays/arrays_fields.f90
- grids/common_types.f90 → parallelisation/common_types.f90
- grids/stella_layouts.f90 → parallelisation/parallelisation_layouts.f90
- utils/redistribute.f90 → parallelisation/redistribute.f90
- calculations/dist_redistribute.f90 → parallelisation/initialise_redistribute.f90

New files:
- parallelisation/timers.f90

The scripts related to the radially global version of $\texttt{stella}$ have been isolated:
- radial_variation/gk_radial_variation.f90
- radial_variation/gk_sources.fpp
- radial_variation/multibox.f90

In addition, within the other scripts, code specific to the radially global version has been moved into subroutines, ensuring a clearer separation from the main $\texttt{stella}$ code. Specifically,
- dist_fn.f90/init_gxyz() → initialise_distribution_function/add_corrections_to_g_for_radial_variation()


### Important nane changes

- The geometry variables have been redefined:
   - <bmag>(ia,iz) = $B / B_ref$
   - <gradx_dot_gradx>(ia,iz) = $|\nabla x|^2$
   - <grady_dot_grady>(ia,iz) = $|\nabla y|^2$
   - <gradx_dot_grady>(ia,iz) = $\nabla x \cdot \nabla y$
   - <B_times_gradB_dot_gradx>(ia,iz) = $B × \nabla B · \nabla x (a*B_ref/B^3)$
   - <B_times_gradB_dot_grady>(ia,iz) = $B × \nabla B · \nabla y (a*B_ref/B^3)$
   - <B_times_kappa_dot_gradx>(ia,iz) = $B × \kappa · \nabla x (a*B_ref/B^2)$
   - <B_times_kappa_dot_grady>(ia,iz) = $B × \kappa · \nabla y (a*B_ref/B^2)$
   - <b_dot_gradz>(ia,iz) = $b \cdot \nabla z$
   - <b_dot_gradz_avg>(iz) = $\sum_\alpha b · \nabla z \mathcal{J} d \alpha / \sum_\alpha \mathcal{J} d \alpha$
   
- Therefore, the following name changes have been implemented:
   - gds22 → gradx_dot_gradx * shat * shat
   - gds2 → grady_dot_grady
   - gds21 → gradx_dot_grady * shat
   - gbdrift0  → B_times_gradB_dot_gradx * 2. * shat
   - gbdrift → B_times_gradB_dot_grady * 2.
   - cvdrift0 → B_times_kappa_dot_gradx * 2. * shat
   - cvdrift →  B_times_kappa_dot_grady * 2.
   - gradpar → b_dot_gradz_avg
   - dgradpardrho  → d_bdotgradz_drho
   - gradpar_c → b_dot_gradz_centredinz (in gk_parallel_streaming.f90)
   - d_gradydotgrady_drho → d_gradydotgrady_drho
   - d_gradxdotgradx_drho → d_gradxdotgradx_drho
   - d_gradxdotgrady_drho → d_gradxdotgrady_drho

## Automatic tests

To ensure that the stella code is functioning correctly, a suite of numerical tests has been implemented and is automatically run on every push to GitHub. These tests cover a wide range of routines and options within stella. The automatic testing infrastructure was implemented by **H. Thienpondt**, who also designed most of the tests, while the full-flux-surface tests were added by **G. Acton** and the electromagnetic tests were contributed by **M. Hardman**. The numerical tests are organized into seven categories. Tests 1–5 use the electrostatic flux-tube version of stella, while tests 6 and 7 target the full-flux-surface and electromagnetic versions, respectively.

- **Test 1**: Confirms that the stella executable exists and that stella runs correctly by
checking that the executable produces output files.
- **Test 2**: Verifies that the magnetic geometries are implemented correctly. This includes
geometries based on Miller parameters, VMEC equilibria, and slab geometry.
- **Test 3**: Checks each term of the gyrokinetic equation independently. It confirms that
the electrostatic potential remains constant when no terms (nor collisions or dissipation)
are included; validates the initialization options for the distribution function, and verifies
the correct evolution of the potential for each term. The (kx , ky ) grid options (“box” and
“range”) are also tested.
- **Test 4**: Tests the parallel boundary conditions: (1) standard twist-and-shift, (2) stellarator-
symmetric, (3) periodic, and (4) zero boundary conditions. Additionally, multiple input flags
are tested simultaneously. In the future, each input parameter would be tested individually.
- **Test 5**: Validates the diagnostics, including growth rates, fluxes, density and temperature,
distribution function, and electrostatic potential.
- **Test 6**: Checks the electrostatic full-flux-surface version of stella. Geometric quantities
are verified first, followed by each term of the gyrokinetic equation.
- **Test 7**: Tests the electromagnetic flux-tube version, verifying each term of the gyrokinetic
equation and simulating a KBM and TAE instability.
- **Test 8**: Checks whether a simulation can be correctly restarted.

The numerical tests can be run locally by executing,

```
cd $STELLA
make create - test - virtualenv
source $STELLA / AUTOMATIC_TESTS / venv / bin / activate
make numerical - tests
```

When debugging, the tests can be run in chunks through,

```
make numerical - tests -1
make numerical - tests -2
make numerical - tests -3
make numerical - tests -4
make numerical - tests -5
make numerical - tests -6
make numerical - tests -7
make numerical - tests -8
```

The input files for the numerical tests are located in,

```
cd AUTOMATIC_TESTS/numerical_tests/
```

## Fixed Bugs

- Broadcast the phase shift (it was only defined on one processor).
- The `theta0` grid was badly defined when using the "range" mode when selecting a range in both kx and ky. This lead to issues in `kperp2` and the Bessel functions.
- Restarted simulations were off by 1 time step.

## Contributors

- **H. Thienpondt** – Code cleanup, restructuring, input reorganization, automatic testing framework, and extensive commenting within the code and documentation, including README files.
- **G. Acton** – Code cleanup, restructuring, input reorganization, and extensive commenting within the code and documentation, including README files.



