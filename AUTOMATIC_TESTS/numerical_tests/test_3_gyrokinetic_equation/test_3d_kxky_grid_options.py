
################################################################################
#               Check all the terms in the gyrokinetic equation                #
################################################################################
# Test the options for the (kx,ky) grid which are 'box' and 'range'
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

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                                BOX NONLINEAR                                 #
#-------------------------------------------------------------------------------
def test_kxky_grid_box_nonlinear(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('kxky_grid_box_nonlinear.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files   
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)    
    print('  -->  The grid_option="box" is evolving correctly during a nonlinear simulation.')
    return
    
#-------------------------------------------------------------------------------
#                                 BOX LINEAR                                   #
#-------------------------------------------------------------------------------
def test_kxky_grid_box_linear(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('kxky_grid_box_linear.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files   
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)    
    print('  -->  The grid_option="box" is evolving correctly during a linear simulation.')
    return
    
#-------------------------------------------------------------------------------
#                                RANGE LINEAR                                  #
#-------------------------------------------------------------------------------
def test_kxky_grid_range_linear_ky1_rangekx(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('kxky_grid_ky1_rangekx_linear.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The grid_option="range" is evolving correctly during a linear simulation.')
    return

#-------------------------------------------------------------------------------
def test_kxky_grid_range_linear_rangeky_kx0(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('kxky_grid_rangeky_kx0_linear.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The grid_option="range" is evolving correctly during a linear simulation.')
    return

#-------------------------------------------------------------------------------
def test_kxky_grid_range_linear_rangeky_kx0_withoutky0(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('kxky_grid_rangeky_kx0_withoutky0_linear.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The grid_option="range" is evolving correctly during a linear simulation.')
    return

#-------------------------------------------------------------------------------
def test_kxky_grid_range_linear_rangeky_rangekx(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('kxky_grid_rangeky_rangekx_linear.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The grid_option="range" is evolving correctly during a linear simulation.')
    return

#-------------------------------------------------------------------------------
#                                RANGE LINEAR                                  #
#-------------------------------------------------------------------------------
def test_kxky_grid_singlemode_linear(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('kxky_grid_singlemode_linear.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files   
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)    
    print('  -->  The grid_option="range" choosing a single (kx,ky) mode is evolving correctly during a linear simulation.')
    return
    

#-------------------------------------------------------------------------------
#            RANGE LINEAR, PARALLEL STREAMING WITH LINKED BOUNDARY             #
#-------------------------------------------------------------------------------
# The modes of a (kx,ky) range are unconnected. Check that the implicit parallel 
# streaming is still applied to them when boundary_option = 'linked', by comparing 
# the singlemode run with a run where the parallel streaming is switched off.
#-------------------------------------------------------------------------------
def test_kxky_grid_singlemode_linear_applies_implicit_parallel_streaming(tmp_path, stella_version, error=False):

    # Run the input file with and without parallel streaming
    input_file = 'kxky_grid_singlemode_linear.in'
    phi2 = {}
    for name, include_streaming in [('with_streaming', True), ('without_streaming', False)]:
        run_path = tmp_path / name
        run_path.mkdir()
        path_input_file = copy_input_file(input_file, run_path)
        if not include_streaming:
            text = path_input_file.read_text().replace('include_nonlinear = .false.', 'include_nonlinear = .false.\n  include_parallel_streaming = .false.')
            path_input_file.write_text(text)
        os.chdir(run_path)
        run_stella(get_stella_path(stella_version), input_file)
        with xr.open_dataset(run_path / input_file.replace('.in', '.out.nc')) as netcdf:
            phi2[name] = netcdf['phi2'].values.copy()

    # The parallel streaming has to change the potential
    assert not np.allclose(phi2['with_streaming'][1:], phi2['without_streaming'][1:], rtol=1e-3), \
        'The implicit parallel streaming has no effect on a (kx,ky) range with linked boundary conditions.'
    print('  -->  The implicit parallel streaming is applied to the unconnected modes of a (kx,ky) range.')
    return
