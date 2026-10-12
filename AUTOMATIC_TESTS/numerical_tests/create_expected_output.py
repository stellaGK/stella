#!/usr/bin/python3
################################################################################
#                 Create the expected output of a numerical test               #
################################################################################
# Only run this script when a change of the numerics of stella is intended, and
# document the reason in the commit message. It runs the given input files and
# stores the time traces of the fields together with the full state of the
# simulation (see <full_state_keys> in run_local_stella_simulation.py), compressed,
# in EXPECTED_OUTPUT.<input_file>.out.nc next to the input file. Usage:
#     python3 create_expected_output.py test_4_gyrokinetic_equation/mirror_implicit.in [...]
#     python3 create_expected_output.py test_5_fluxtube/test_useless_flags1.in --name test_useless_flags
# The number of processors can be set with --nproc (default: 4), and the stella
# executable with the STELLA_EXE_PATH environment variable. Note that the numerical
# reproducibility tests (test 3) have their own create_expected_output.py.
################################################################################

import os, re, sys
import shutil
import pathlib
import argparse
import tempfile
import xarray as xr

# Package to run stella. Rename this script first, since convert_inputFile.py (which is executed
# by run_local_stella_simulation.py) converts all input files in the current folder when it is
# executed as __main__.
# It also locates files relative to a test module, which lives one folder deeper than this script.
__name__ = 'create_expected_output'
__file__ = str(pathlib.Path(__file__).absolute().parent / 'test_1_whether_stella_runs' / 'create_expected_output.py')
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Time traces which are read by the tests, besides the full state
time_trace_keys = ['phi2', 'apar2', 'bpar2']

parser = argparse.ArgumentParser(description='Create the expected output of numerical tests.')
parser.add_argument('input_files', nargs='+', type=pathlib.Path, help='input files of the tests')
parser.add_argument('--name', help='stem of the expected output file, if it differs from the input file')
parser.add_argument('--nproc', type=int, default=4, help='number of processors (default: 4)')
args = parser.parse_args()
if args.name and len(args.input_files) > 1: sys.exit('ERROR: --name can only be used with a single input file.')

stella_path = get_stella_path('master')
input_files = [input_file.absolute() for input_file in args.input_files]  # before changing directories
for input_file in input_files:
    test_directory = input_file.parent
    with tempfile.TemporaryDirectory() as tmp_path:
        tmp_path = pathlib.Path(tmp_path)
        shutil.copyfile(input_file, tmp_path / input_file.name)

        # Copy the VMEC equilibrium if the input file uses one
        vmec_file = re.search(r"^\s*vmec_filename\s*=\s*['\"]([^'\"]+)['\"]", input_file.read_text(), re.M)
        if vmec_file: shutil.copyfile(test_directory / vmec_file.group(1), tmp_path / vmec_file.group(1))

        # Run stella and keep the time traces and the full state
        os.chdir(tmp_path)
        run_stella(stella_path, input_file.name, nproc=args.nproc)
        with xr.open_dataset(tmp_path / input_file.name.replace('.in', '.out.nc')) as netcdf:
            keys = [key for key in time_trace_keys + full_state_keys if key in netcdf.variables]
            dataset = netcdf[keys].load()
        dataset.attrs['full_state_reference'] = 1  # compare the full state in the tests
        os.chdir(test_directory)

    # Write the expected output
    stem = args.name if args.name else input_file.stem
    expected_netcdf_file = test_directory / f'EXPECTED_OUTPUT.{stem}.out.nc'
    dataset.to_netcdf(expected_netcdf_file, encoding={key: {'zlib': True, 'complevel': 9} for key in keys})
    print(f'Created {expected_netcdf_file.relative_to(test_directory.parent)} ({os.path.getsize(expected_netcdf_file)/1e3:.0f} kB)')
