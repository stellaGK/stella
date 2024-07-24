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
will make all subsequent tests fails. Finally it is tested whether the potential remains 
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

    export GK_HEAD_DIR=your_path_to_stella/stella
    export STELLA_EXE_PATH=your_path_to_stella/stella/stella
    cd $GK_HEAD_DIR/
    make create-test-virtualenv
    source $GK_HEAD_DIR/tests/venv/bin/activate

This will create a Python [virtual environment][venv] with the packages needed
for running the tests. You can then run all the tests:
    
    make run-automated-numerical-tests-for-stella
    
If you would like to see more information while running the tests run: 
    
    make run-automated-numerical-tests-for-stella-verbose
    
(TODO-HT) Besides the numerical tests create a package for quick
and slow physics tests, used as benchmarks.
    

Writing new Python tests
------------------------

The testing framework is setup to automatically find files called `test_*.py`
and run any functions it finds in them called `test_*`. Writing a new automated
test for `stella` is as "simple" as writing a new function that starts with `test_`.

Let's look at a basic test, `test_whether_stella_runs` in 
[`test_1_whether_stella_runs/test_whether_stella_runs.py`](test_1_whether_stella_runs/test_whether_stella_runs.py) to see how to write 
a test that runs `stella` and checks the output netCDF file.
The test module is split up into three tests, first we simply
run a local stella simulation inside the folder 'tmp_path':

```python

# Package to run stella locally
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Global variables  
input_filename = 'miller_nonlinear_CBC.in'  
local_stella_run_directory = 'Not/Run/Yet'

def test_whether_we_can_run_a_local_stella_simulation(tmp_path):
    '''Run a local stella simulation in a temporary folder <tmp_path>.'''  
    
    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global local_stella_run_directory
    local_stella_run_directory = tmp_path
    
    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path)
    print('\n  -->  Successfully ran a local stella simulation.')
    return 
    
```

First, we need to import the package [`run_local_stella_simulation.py`](run_local_stella_simulation.py) which
contains the functions that are needed to set-up and run stella inside `tmp_path`. 
Specifically, it will copy the input file (and VMEC file) to `tmp_path` and it
will run the stella simulation using `mpirun` through `subprocess.run`.

```python
    # Package to run stella locally
    module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
    with open(module_path, 'r') as file: exec(file.read())
```

By default stella will run on 2 processors, to ensure that the parallelization 
is implemented correctly. It is very common to forget a `broadcast(variable)` 
statement, hence it is important to run tests on more than 1 processor to catch 
these errors. If a user wishes to perform the test on more processors they can 
change the '2' in the following line of code inside [`run_local_stella_simulation.py`](run_local_stella_simulation.py):

```python
    subprocess.run(['mpirun', '-np', '2', stella_path, input_filename], check=True)  
```

We define two global variables so that they are accessible to all three tests.

```python
    # Global variables  
    input_filename = 'miller_nonlinear_CBC.in'  
    local_stella_run_directory = 'Not/Run/Yet
```

Each function starting with `test_*.py` will be executed the `pytest` framework.
The function takes one argument called `tmp_path`: this is a path to a temporary
directory done magically by the test framework `pytest`. This is used so we can run
the test in an isolated environment, and not worry about other tests potentially
interfering with output files, for example. We can help out future developers by 
providing a brief description of the aim of the test, how it works, and so on, in the docstring.

```python
def test_whether_we_can_run_a_local_stella_simulation(tmp_path):
    '''Run a local stella simulation in a temporary folder <tmp_path>.'''      
```

If we want to perform multiple tests on the same stella simulation, we will
run the stella simulation as a first test, and we will let the other tests
know about the location of the stella simulation by setting the global variable:

```python
    global local_stella_run_directory
    local_stella_run_directory = tmp_path
```

Next, we run the local stella simulation, by using the [`run_local_stella_simulation.py`](run_local_stella_simulation.py) package. 
This function will also check that `stella` completed successfully, hence this 
is a useful first test to see whether stella runs successfully.

```python
    run_local_stella_simulation(input_filename, tmp_path)
```

Next we run a second test, to check whether all the expected output files 
have been generated within the `tmp_path` folder. If an output file is missing
an error will be printed to the command prompt, to help developers track
down the specific error. If one (or more) of the output files are missing, 
the `assert` statement will get triggered and an `AssertionError: <error message>` 
which will label the test as `Failed` in the command prompt. 

```python
def test_whether_all_output_files_are_gerenated_when_running_stella(error=False):  
    '''To ensure that stella ran correctly, check that all output files are generated.'''
    
    # Gather the output files generated during the local stella run inside <tmp_path>
    local_files = os.listdir(local_stella_run_directory)
    
    # Create a list of the output files we expect when stella has been run
    stem = input_filename.replace('.in','')
    expected_files = ['in', 'out.nc', 'geometry', 'fluxes', 'omega', 'out', 'final_fields', 'species.input', 'error']
    expected_files = [stem+'.'+extension for extension in expected_files]
    
    # Check whether all these output files are present
    for expected_file in expected_files:
        if not (expected_file in local_files):
            print(f'ERROR: The "{expected_file}" output file was not generated when running stella.'); error = True
            
    # The <pytest> module will check whether all <assert> statements are true,
    # if it runs into a statement which is false, the test will be labeled as
    # "Failed" and the string in the second argument of the <assert> statement 
    # will be printed to the command prompt as "AssertionError: <error message>"
    assert (not error), f'Some output files were not generated when running stella.'
    
    # The test will stop as soon as an <assert> error was triggered, hence we
    # will only reach this final line of code if the test ran successfully
    print(f'  -->  All the expected files (.in, .geometry, .out.nc, .fluxes, .omega, ...) are generated.')
    return 
```

The third test can be found within [`test_1_whether_stella_runs/test_whether_stella_runs.py`](test_1_whether_stella_runs/test_whether_stella_runs.py), which will not be discussed in detail here.

If you wish to set-up a single test the following lines of code are sufficient. Where
you need to add the `EXPECTED_OUTPUT_*_out.nc` file (or any other output file) to the test folder,
obtained by running a bench-marked stella version, so that the output files of the local stella 
version can be compared against the expected output. Note that if the results are expected to change, 
due to a valid change to stella, the `EXPECTED_OUTPUT_*_out.nc` file should be updated on Github.

```python
    # Python modules
    import os, sys
    import pathlib 
    import numpy as np
    import xarray as xr  

    # Package to run stella 
    module_path = str(pathlib.Path(__file__).parent.parent / 'run_local_stella_simulation.py')
    with open(module_path, 'r') as file: exec(file.read())
         
    # Check a specific kind of simulation, use a sensible name for the 
    # test, so that it is clear in the command prompt what is being tested.
    def test_whether_miller_linear_evolves_correctly(tmp_path):

        # Input file name  
        input_filename = 'miller_geometry_linear.in'  
        
        # Run a local stella simulation
        run_local_stella_simulation(input_filename, tmp_path)
         
        # File names  
        local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc') 
        expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   
        
        # Check whether the potential data matches in the netcdf file
        with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        
            # Check whether all the keys are present
            assert set(local_netcdf.keys()) == set(expected_netcdf.keys()), f'The netcdf file contains different quantities.'
            
            # Read the time axis
            local_time = local_netcdf['t']
            expected_time = expected_netcdf['t']
            
            # Read the potential axis
            local_phi2 = local_netcdf['phi2']
            expected_phi2 = expected_netcdf['phi2'] 
            
            # Check whether we have the same time and potential data
            if not (np.allclose(local_time, expected_time, equal_nan=True)):
                print('\nERROR: The time axis does not match in the netCDF files.'); error = True
                print('\nCompare the time arrays in the local and expected netCDF files:')
                compare_local_array_with_expected_array(local_time, expected_time)  
            if not (np.allclose(local_phi2, expected_phi2, equal_nan=True)):
                print('\nERROR: The potential data does not match in the netCDF files.'); error = True 
                print('\nCompare the potential arrays in the local and expected netCDF files:')
                compare_local_array_with_expected_array(local_phi2, expected_phi2) 
            assert (not error), f'The potential data does not match in the netCDF files.'  
                    
        print(f'\n  -->  The potential is evolving correctly in a linear flux-tube simulation using Miller geometry ({int(local_netcdf["nproc"])} CPUs).')
        return
```

[pytest]: https://pytest.org
[xarray]: http://xarray.pydata.org
[venv]: https://docs.python.org/3/library/venv.html
