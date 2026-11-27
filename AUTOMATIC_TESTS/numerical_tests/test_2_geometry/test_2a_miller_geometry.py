################################################################################
#                           Check the Miller geometry                          #
################################################################################
# Test all the geometry arrays when using a Miller equilibrium. 
################################################################################

# Python modules
import pytest
import os, sys
import pathlib 
import numpy as np
import xarray as xr

# Package to run stella 
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Global variables  
input_filename = 'miller_geometry.in'

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                         Run local stella simulation                          #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def stella_run(tmp_path_factory, stella_version):
    '''Run a local stella simulation, which is shared by all the tests in this 
    module, so that each test can also be run on its own.'''
    filename = input_filename if stella_version=='master' else input_filename.replace('.in', f'_v{stella_version}.in')
    run_data = run_local_stella_simulation(filename, tmp_path_factory.mktemp('miller_geometry'), stella_version)
    run_data['input_file_stem'] = filename.replace('.in','')
    # Older stella versions named the Miller files <millerlocal.*> instead of <geometry_miller.*>
    new_name_exists = (run_data['tmp_path'] / f'geometry_miller.{run_data["input_file_stem"]}.input').exists()
    run_data['miller_file_name'] = 'geometry_miller' if new_name_exists else 'millerlocal'
    return run_data

#-------------------------------------------------------------------------------
#                    Check whether output files are present                    #
#-------------------------------------------------------------------------------
def test_whether_miller_output_files_are_present(stella_run, error=False):  
    
    # Gather the output files generated during the local stella run
    local_files = os.listdir(stella_run['tmp_path'])
    input_file, miller_file_name = stella_run['input_file_stem'], stella_run['miller_file_name']
    
    # Check whether all the output files we expect are present
    expected_files = [f'{miller_file_name}.{input_file}.input', f'{miller_file_name}.{input_file}.output', f'{input_file}.geometry']
    for expected_file in expected_files:
        if not (expected_file in local_files):
            print(f'ERROR: The "{expected_file}" output file was not generated when running stella.'); error = True
    assert (not error), f'Some output files were not generated when running stella.'
    print(f'  -->  All the expected files ({miller_file_name}.input, {miller_file_name}.output, .geometry) are generated.')
    return 

#-------------------------------------------------------------------------------
#                    Check whether Miller output files match                   #
#-------------------------------------------------------------------------------
def test_whether_miller_output_files_are_correct(stella_run):
    '''Check that the results are identical to a previous run.'''
    
    # File names
    local_directory = stella_run['tmp_path']
    input_file, miller_file_name = stella_run['input_file_stem'], stella_run['miller_file_name']
    local_geometry_file = local_directory / f'{input_file}.geometry' 
    expected_geometry_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.miller_geometry.geometry' 
    local_miller_input_file = local_directory / f'{miller_file_name}.{input_file}.input' 
    expected_miller_input_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.miller_geometry.millerlocal.input' 
    local_miller_output_file = local_directory / f'{miller_file_name}.{input_file}.output' 
    expected_miller_output_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.miller_geometry.millerlocal.output'
    
    # Compare text files (first check the input file to save <shat>)
    compare_geometry_files(local_geometry_file, expected_geometry_file, error=False)
    shat = compare_miller_input_files(local_miller_input_file, expected_miller_input_file, error=False)
    compare_miller_output_files(local_miller_output_file, expected_miller_output_file, shat=shat, error=False)
    print(f'  -->  Geometry output file matches.')
    return

#-------------------------------------------------------------------------------
#              Check whether the data in the netcdf file matches               #
#-------------------------------------------------------------------------------
def test_whether_miller_geometry_data_in_netcdf_file_is_correct(stella_run, error=False): 
    compare_geometry_in_netcdf_files(stella_run, error=False)  
    print('  -->  All Miller geometry data in the netcdf file matches the expected output.')
    return
    


#-------------------------------------------------------------------------------
#          Check whether the .geometry file follows from the netcdf file       #
#-------------------------------------------------------------------------------
# Every variable of the <.geometry> file is also written to the netcdf file, so
# that the geometry can be rebuilt from the netcdf output alone (#109).
#-------------------------------------------------------------------------------
def test_whether_geometry_file_follows_from_netcdf_file(stella_run, error=False):

    # Read the <.geometry> file: two lines with the scalars, then the arrays per field line
    local_directory = stella_run['tmp_path']
    input_file = stella_run['input_file_stem']
    with open(local_directory / f'{input_file}.geometry', 'r') as file: lines = file.readlines()
    scalar_names = lines[0].replace('#', '').split()
    scalar_values = [float(value) for value in lines[1].replace('#', '').split()]
    array_names = lines[3].replace('#', '').split()
    array_values = np.array([[float(value) for value in line.split()] for line in lines[4:] if line.strip()])

    # The names of the variables in the netcdf file
    netcdf_names = {'qinp': 'q', '1/Grho': 'one_over_nablarho', 'b.Gz': 'b_dot_gradz', 'Gy.Gy': 'grady_dot_grady',
        'Gx.Gy': 'gradx_dot_grady', 'Gx.Gx': 'gradx_dot_gradx', 'BxGB.Gy': 'B_times_gradB_dot_grady',
        'Bxkappa.Gy': 'B_times_kappa_dot_grady', 'BxGB.Gx': 'B_times_gradB_dot_gradx'}

    # Compare them with the netcdf file, to the 4 digits of the <.geometry> file
    with xr.open_dataset(local_directory / f'{input_file}.out.nc') as netcdf:
        nalpha = netcdf.sizes['alpha']
        for name, value in zip(scalar_names, scalar_values):
            assert np.isclose(float(netcdf[netcdf_names.get(name, name)]), value, rtol=1e-3, atol=1e-4), f'{name} differs from the netcdf file.'
        for i, name in enumerate(array_names):
            if name == 'alpha': netcdf_array = np.repeat(netcdf['alpha'].values, netcdf.sizes['zed'])
            elif name == 'zed': netcdf_array = np.tile(netcdf['zed'].values, nalpha)
            else: netcdf_array = netcdf[netcdf_names.get(name, name)].values.flatten()
            assert np.allclose(netcdf_array, array_values[:, i], rtol=1e-3, atol=1e-4), f'{name} differs from the netcdf file.'
    print('  -->  Every variable of the .geometry file is also in the netcdf file.')
    return
