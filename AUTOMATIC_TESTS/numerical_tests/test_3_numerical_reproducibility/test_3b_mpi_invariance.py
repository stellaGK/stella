################################################################################
#   Check that the results do not depend on the parallelisation of stella      #
################################################################################
# Optimising the memory usage or speed of stella often means changing how the
# data is distributed over the MPI processes. The physical result should not
# depend on this. Here each input file is run on 1 processor (the reference)
# and on several other decompositions (2, 3 and 4 processors, where 3 gives an
# uneven split of the grids, and different xyzs_layout and vms_layout), and the
# full state of the simulation is compared between them.
#
# These tests do not need any expected output, so they keep working when a
# change of the numerics is intended. The relevant parameters are:
#       &parallelisation
#         xyzs_layout = 'yxzs'
#         vms_layout = 'vms'
#       /
################################################################################

# Python modules
import pytest
import os, sys
import shutil
import pathlib
import numpy as np
import xarray as xr

# Package to run stella
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Settings of the reproducibility tests
module_path = str(pathlib.Path(__file__).parent / 'reproducibility_settings.py')
with open(module_path, 'r') as file: exec(file.read())

# Decompositions which are compared to the single-processor reference
decompositions = {
    'nproc2': (2, {}),
    'nproc3': (3, {}),
    'nproc4_xzys': (4, {'xyzs_layout': "'xzys'"}),
    'nproc3_zyxs_mvs': (3, {'xyzs_layout': "'zyxs'", 'vms_layout': "'mvs'"}),
}

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
def run_stella_with_decomposition(input_file, folder, stella_version, nproc, layouts):
    '''Run <input_file> in <folder> on <nproc> processors with the given layouts.'''
    folder.mkdir(parents=True, exist_ok=True)
    path_input_file = copy_input_file(input_file, folder)
    if layouts:
        with open(path_input_file, 'a') as file:
            file.write('\n&parallelisation\n')
            for key, value in layouts.items(): file.write(f'  {key} = {value}\n')
            file.write('/\n')
    os.chdir(folder)
    run_stella(get_stella_path(stella_version), input_file, nproc=nproc)
    return folder / input_file.replace('.in', '.out.nc')

#-------------------------------------------------------------------------------
# The single-processor reference is only run once per input file
reference_netcdf_files = {}
def get_reference_netcdf_file(input_file, tmp_path_factory, stella_version):
    if input_file not in reference_netcdf_files:
        folder = tmp_path_factory.mktemp(input_file.replace('.in', '_nproc1'))
        reference_netcdf_files[input_file] = run_stella_with_decomposition(input_file, folder, stella_version, 1, {})
    return reference_netcdf_files[input_file]

#-------------------------------------------------------------------------------
def mark_known_bugs(input_file):
    if input_file in known_mpi_bugs:
        return pytest.param(input_file, marks=pytest.mark.xfail(strict=True, reason=known_mpi_bugs[input_file]))
    return input_file

#-------------------------------------------------------------------------------
#           Compare the full state between different decompositions           #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('decomposition', decompositions.keys())
@pytest.mark.parametrize('input_file', [mark_known_bugs(f) for f in input_files])
def test_whether_results_are_independent_of_mpi_decomposition(input_file, decomposition, tmp_path, tmp_path_factory, stella_version):

    # These input files only exist for the current stella version
    if stella_version != 'master': pytest.skip('Only implemented for the master branch of stella.')

    # Run stella on a single processor, and with the given decomposition
    nproc, layouts = decompositions[decomposition]
    reference_netcdf_file = get_reference_netcdf_file(input_file, tmp_path_factory, stella_version)
    local_netcdf_file = run_stella_with_decomposition(input_file, tmp_path, stella_version, nproc, layouts)

    # Only compare the keys which are present in the reference
    with xr.open_dataset(reference_netcdf_file) as reference_netcdf:
        keys = [key for key in state_keys if key in reference_netcdf.variables]

    # Compare the full arrays
    compare_netcdf_quantities_normwise(local_netcdf_file, reference_netcdf_file, keys, rtol=rtol, atol=atol, label='nproc1')
    print(f'  -->  The results of {input_file} do not depend on the MPI decomposition ({decomposition}).')
    return
