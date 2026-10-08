Automated stella tests
======================

An extensive package of automated stella tests have been developed to ensure
that the majority of bugs introduced by developers can be caught. There are four
main sets of test:
    - Quick Python numerical tests of stella, testing the exact time evolution
    - Quick Python bench-mark tests, testing whether stella matches roughly with other codes
    - Slow Python bench-mark tests, testing whether stella matches accurately with other codes
    - Fortran tests, testing specific Fortran routines, this is underdeveloped
    
The numerical tests are performed in a logical order. First it is tested whether stella runs
and produces output files, since if this test fails all subsequent tests will fail as well. 
Next, the different geometry options (Miller, VMEC) are tested, since errors in the geometry
will make all subsequent tests fails. Then the numerical reproducibility is tested: the full 
state of short nonlinear simulations (fields, moments, fluxes and distribution functions) is
compared with a tight tolerance, and it is checked that the results do not depend on the number
of MPI processes. Finally it is tested whether the potential remains 
constant if none of the gyrokinetic terms nor dissipation are included, and whether the 
potential (and distribution function) are initialized the same in the current stella version
compared to a previous run, seeing that the time evolution will differ if the initialization has
been changed. After these initial checks, each term in the gyrokinetic equation is tested
separately, as well as their implicit/explicit implementation, to catch any bugs introduced
into the gyrokinetic terms. Note that if a valid change was made to a term or algorithm, this
change should be clearly documented in the tests and the `EXPECTED_OUTPUT_*_out.nc` should 
be updated.


Running automated stella tests on Github
----------------------------------------
All these test are run automatically on Github at each push and pull request. This allows
us to catch bugs (or valid changes) automatically on Github. It is very important that
each user ensures that all tests are passed successfully on Github. These tests are performed
automatically on Github by workflows defined in ['.github/workflows/tests.yml'](test_1_whether_stella_runs/test_whether_stella_runs.py).


Running automated stella tests with python
------------------------------------------

The automated stella tests here are written using the Python [pytest][pytest] package,
and make use of [xarray][xarray] for reading in the data. To install these
packages in a standalone environment, you can run:

    export STELLA_HEAD_DIR=your_path_to_stella/stella
    export STELLA_EXE_PATH=your_path_to_stella/stella/stella
    make create-test-virtualenv
    source $STELLA_HEAD_DIR/AUTOMATIC_TESTS/venv/bin/activate

This will create a Python [virtual environment][venv] with the packages needed
for running the tests. You can then run all the tests in the main stella directory using:
    
    make numerical-tests
    
If you would like to see more information while running the tests run: 
    
    make numerical-tests-verbose
    
Slow tests, such as the KBM and TAE benchmarks of electromagnetic stella, are marked with
`@pytest.mark.slow` and are skipped by default. Run them with:
    
    make numerical-tests-slow
    
Before and after optimising the memory usage or speed of stella, run the numerical 
reproducibility tests, which compare the full state of the simulations (fields, moments, 
fluxes and distribution functions) with a tight tolerance, and check that the results do 
not depend on the number of MPI processes, see 
[`numerical_tests/test_3_numerical_reproducibility`](numerical_tests/test_3_numerical_reproducibility/README.md):

    make numerical-tests-3
    
(TODO-HT) Besides the numerical tests create a package for quick
and slow physics tests, used as benchmarks.


Number of threads
-----------------
The simulations can be run on as many threads as your local computer supports. The number of threads can be set in config.ini.
    

Writing new Python tests
------------------------

The testing framework is setup to automatically find files called `test_*.py`
and run any functions it finds in them called `test_*`. Writing a new automated
test for `stella` is as "simple" as writing a new function that starts with `test_`.
Each test folder contains the input files, and the `EXPECTED_OUTPUT.*` files obtained 
by running a bench-marked stella version. If the results change due to a valid change 
of stella, the `EXPECTED_OUTPUT.*` files should be updated, and the change documented.

Every test module first loads [`run_local_stella_simulation.py`](run_local_stella_simulation.py), 
which contains the functions to run stella inside a temporary folder (copying the input 
file and VMEC file, and running `mpirun -np <nproc>` with `nproc` from `config.ini`), 
and the functions to compare the output with the expected output:

```python
# Package to run stella 
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")
```

A single test runs stella inside the temporary folder `tmp_path` provided by `pytest`, 
and compares the output with the expected output, e.g.:

```python
def test_whether_miller_linear_evolves_correctly(tmp_path, stella_version):
    run_data = run_local_stella_simulation('miller_geometry_linear.in', tmp_path, stella_version)
    compare_local_potential_with_expected_potential(run_data=run_data)
```

If several tests check the same stella simulation, run it in a fixture with 
`scope="module"`, rather than in the first test. The simulation is then run once and 
shared by all the tests that request it, and each test can still be run on its own 
(e.g. with `pytest -k <name>`), see 
[`test_1_whether_stella_runs/test_whether_stella_runs.py`](numerical_tests/test_1_whether_stella_runs/test_whether_stella_runs.py):

```python
@pytest.fixture(scope="module")
def local_stella_run_directory(tmp_path_factory, stella_version):
    tmp_path = tmp_path_factory.mktemp('stella_run')
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    return tmp_path

def test_whether_all_output_files_are_gerenated_when_running_stella(local_stella_run_directory):
    local_files = os.listdir(local_stella_run_directory)
    ...
```

Some things to keep in mind when comparing data:

- `np.allclose()` broadcasts arrays, so a simulation which wrote a single time step 
  would match any time trace. Check the shapes first, e.g. with `assert_same_shape()`.
- Compare full arrays with `compare_netcdf_quantities_normwise()`, which requires 
  `max|local - expected| <= rtol * max|expected| + atol`, rather than element by element, 
  since elements which are zero up to round-off errors would otherwise fail on noise.
- Text files are written with a limited number of digits, so compare them with a 
  relative tolerance which matches the printed precision, rather than rounding to a 
  fixed number of decimals (which hides errors in quantities smaller than one).
- Make sure that a check can actually fail, e.g. by perturbing an input parameter 
  slightly and checking that the test catches it.
- Slow tests should be marked with `@pytest.mark.slow`, they are skipped by default.

[pytest]: https://pytest.org
[xarray]: http://xarray.pydata.org
[venv]: https://docs.python.org/3/library/venv.html
