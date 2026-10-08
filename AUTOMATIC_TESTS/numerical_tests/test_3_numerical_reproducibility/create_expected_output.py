#!/usr/bin/python3
################################################################################
#          Create the expected output for the reproducibility tests            #
################################################################################
# Only run this script when a change of the numerics of stella is intended, and
# document the reason in the commit message. It runs every input file of test 3
# on a single processor and stores the compared quantities, compressed, in
# EXPECTED_OUTPUT.<input_file>.out.nc. Usage, from this folder:
#     python3 create_expected_output.py [input_file.in ...]
# The stella executable can be set with the STELLA_EXE_PATH environment variable.
################################################################################

import os, sys
import shutil
import pathlib
import tempfile
import numpy as np
import xarray as xr

# Package to run stella
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Settings of the reproducibility tests
module_path = str(pathlib.Path(__file__).parent / 'reproducibility_settings.py')
with open(module_path, 'r') as file: exec(file.read())

test_directory = pathlib.Path(__file__).parent.absolute()
for input_file in (sys.argv[1:] or input_files):
    with tempfile.TemporaryDirectory() as tmp_path:
        tmp_path = pathlib.Path(tmp_path)
        run_local_stella_simulation(input_file, tmp_path, 'master', nproc=nproc_expected_output)
        with xr.open_dataset(tmp_path / input_file.replace('.in', '.out.nc')) as netcdf:
            keys = [key for key in regression_keys if key in netcdf.variables and key != 't']
            dataset = netcdf[keys].load()
        expected_netcdf_file = test_directory / f'EXPECTED_OUTPUT.{input_file.replace(".in","")}.out.nc'
        dataset.to_netcdf(expected_netcdf_file, encoding={key: {'zlib': True, 'complevel': 9} for key in keys})
        print(f'Created {expected_netcdf_file.name} ({os.path.getsize(expected_netcdf_file)/1e3:.0f} kB)')
    os.chdir(test_directory)
