Numerical reproducibility tests (test 3)
========================================

These tests guard the numerics of stella while its memory usage or speed is being
optimised. Such changes should not alter the results, so these tests compare the
**full state** of short nonlinear simulations with a tight tolerance, instead of only
the time trace of |phi|^2 as most other numerical tests do.

Run them from the main stella directory with:

    make numerical-tests-3


What is tested
--------------

- `test_3a_regression_of_full_state.py`: runs every input file on 1 processor and
  compares the fields (phi, apar, bpar) on the full (t, tube, z, kx, ky) grid, the
  moments, the fluxes, the frequency, the distribution functions g and h on the
  (z, vpa, mu) grid, and the geometry (bmag, kperp2) with `EXPECTED_OUTPUT.*.out.nc`.
- `test_3b_mpi_invariance.py`: runs every input file on 2, 3 and 4 processors, and with
  different `xyzs_layout` and `vms_layout`, and compares the full state with the
  single-processor run. No expected output is needed.
- `test_3c_diagnostics_do_not_affect_evolution.py`: checks that writing the diagnostics
  every time step (nwrite = 1) gives the same evolution as nwrite = 10.

The input files all derive from `es_nonlinear.in` (nonlinear, shaped Miller geometry,
kinetic ions and electrons, all gyrokinetic terms, hyper dissipation). Each other input
file switches on one numerical scheme, listed at the top of the file: electromagnetic
fields, adiabatic electrons, implicit/explicit Dougherty and Fokker-Planck collisions,
the rk4 and euler explicit algorithms, flip_flop, the mirror term without semi-Lagrange,
and explicit parallel streaming and mirror terms.


Tolerance
---------

Every quantity must satisfy `max|local - expected| <= rtol * max|expected| + atol` with
`rtol = 1e-8` (see `reproducibility_settings.py`). For reference, changing the number of
processors or the layouts changes the results by less than 1e-10, and the results are
bit-for-bit identical when a simulation is repeated. Perturbing a single input parameter
(delt, vnew_ref, tprim or kappa) by a relative amount of 1e-6 makes these tests fail.


Known bugs
----------

If the results of an input file depend on the number of processors because of a known
bug, add the input file to `known_mpi_bugs` in `reproducibility_settings.py`. Its MPI
invariance tests are then marked as strict expected failures (xfail). Once the bug is
fixed, the test reports an XPASS which makes it fail, as a reminder to remove the input
file from `known_mpi_bugs` (and to recreate its expected output if the single-processor
result has changed). There are currently no known bugs; the implicit Dougherty and
Fokker-Planck collision operators used to depend on the number of processors.


Intended changes of the numerics
--------------------------------

If a change of the numerics is intended, document it in the commit message and recreate
the expected output (on a single processor) with:

    cd AUTOMATIC_TESTS/numerical_tests/test_3_numerical_reproducibility
    python3 create_expected_output.py                  # all input files
    python3 create_expected_output.py flip_flop.in     # or only some of them
