#!/usr/bin/env python3
################################################################################
#          Make the equilibrium files of test_2g_gvec_vs_vmec_geometry         #
################################################################################
# Makes the three equilibrium files that test_2g_gvec_vs_vmec_geometry.py
# compares, in AUTOMATIC_TESTS/numerical_tests/test_2_geometry:
#
#   wout_tok_circ.nc   VMEC++ equilibrium of the circular tokamak, from
#                      DOCUMENTATION/GVEC_equilibria/circular_tok/vmec/input.tok_circ
#   gvec_tok_circ.nc   stella geometry file (gvec_to_stella.py, s = 0.25) of the
#                      GVEC equilibrium of the same tokamak, solved by GVEC from
#                      DOCUMENTATION/GVEC_equilibria/circular_tok/gvec/parameter.ini
#   gvec_w7x_kjm.nc    stella geometry file (gvec_to_stella.py, s = 0.49, one
#                      field period) of the GVEC state converted from the VMEC
#                      equilibrium wout_w7x_kjm.nc with pygvec convert-wout
#   gvec_*_Bminus.nc   the same geometry files with the magnetic field reversed
#                      (reverse_magnetic_field.py)
#   GVEC_*_State.dat, GVEC_*_parameter.ini
#                      the GVEC states and parameter files, for the direct path
#
# The tokamak is solved independently by both codes, from the same boundary and
# profiles; the W7-X state is the VMEC equilibrium represented in GVEC, so it
# tests the 3-D terms of the geometry with the same equilibrium in both codes.
# Both geometry files use the minor radius Aminor_p of the VMEC file as the
# reference length, as stella's VMEC geometry does.
#
# Requirements (made with these versions on 25/11/2026):
#   pip install gvec==1.5.1 vmecpp==0.9.4 netCDF4
#
# Usage, from any folder (the GVEC and VMEC runs are made in <work_folder>):
#   python make_gvec_vs_vmec_files.py <work_folder>
################################################################################

import sys
import shutil
import pathlib
import subprocess
import netCDF4 as nc

here = pathlib.Path(__file__).resolve().parent
test_folder = here.parent
stella_folder = test_folder.parents[2]
documentation_folder = stella_folder / 'DOCUMENTATION' / 'GVEC_equilibria'
sys.path.insert(0, str(documentation_folder))
sys.path.insert(0, str(here))
from gvec_to_stella import gvec_to_stella
from reverse_magnetic_field import reverse_magnetic_field_of_parameter_file, reverse_magnetic_field_of_state_file

# The z grid of the test inputs (nzed = 48) and the surfaces of the test inputs
nzgrid = 24
s_tokamak = 0.25
s_w7x = 0.49


def aminor(wout_file):
    with nc.Dataset(wout_file) as wout:
        return float(wout['Aminor_p'][...])


def write_geometry_files(name, state_file, parameter_file, **options):
    '''Write the geometry file of <state_file> and of its reversed magnetic field,
    and copy the state and parameter files to the test folder.'''
    shutil.copyfile(state_file, test_folder / f'GVEC_{name}_State.dat')
    shutil.copyfile(parameter_file, test_folder / f'GVEC_{name}_parameter.ini')
    state_file_Bminus = state_file.parent / f'GVEC_{name}_State_Bminus.dat'
    parameter_file_Bminus = parameter_file.parent / f'GVEC_{name}_parameter_Bminus.ini'
    reverse_magnetic_field_of_state_file(state_file, state_file_Bminus)
    reverse_magnetic_field_of_parameter_file(parameter_file, parameter_file_Bminus)
    for suffix, state, parameters in [('', state_file, parameter_file), ('_Bminus', state_file_Bminus, parameter_file_Bminus)]:
        gvec_to_stella(state, parameter_file=parameters, output_file=test_folder / f'gvec_{name}{suffix}.nc',
                       verbose=False, **options)


def make_tokamak_files(work_folder):
    import vmecpp

    # VMEC++ equilibrium
    vmec_folder = work_folder / 'tok_circ_vmec'
    vmec_folder.mkdir(parents=True, exist_ok=True)
    vmec_input = vmecpp.VmecInput.from_file(documentation_folder / 'circular_tok' / 'vmec' / 'input.tok_circ')
    vmec_output = vmecpp.run(vmec_input, verbose=False)
    vmec_output.wout.save(str(vmec_folder / 'wout_tok_circ.nc'))
    shutil.copyfile(vmec_folder / 'wout_tok_circ.nc', test_folder / 'wout_tok_circ.nc')

    # GVEC equilibrium, solved from the boundary
    gvec_folder = work_folder / 'tok_circ_gvec'
    gvec_folder.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(documentation_folder / 'circular_tok' / 'gvec' / 'parameter.ini', gvec_folder / 'parameter.ini')
    subprocess.run(['pygvec', 'run', 'parameter.ini', '-q'], cwd=gvec_folder, check=True)
    write_geometry_files('tok_circ', gvec_folder / 'GVEC_TOK_CIRC_State_final.dat',
                         gvec_folder / 'parameter_GVEC_TOK_CIRC_final.ini',
                         a_minor_override=aminor(test_folder / 'wout_tok_circ.nc'),
                         s_input=s_tokamak, nalpha=1, nzgrid=nzgrid)


def make_w7x_file(work_folder):

    # GVEC state of the VMEC equilibrium, with a finer radial grid than the
    # default of convert-wout, so that the radial derivatives are accurate
    gvec_folder = work_folder / 'w7x_kjm_gvec'
    gvec_folder.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(test_folder / 'wout_w7x_kjm.nc', gvec_folder / 'wout_w7x_kjm.nc')
    subprocess.run(['pygvec', 'convert-wout', 'wout_w7x_kjm.nc', '.', '--param', 'sgrid_nElems=40', '-q'],
                   cwd=gvec_folder, check=True)
    write_geometry_files('w7x_kjm', gvec_folder / 'GVEC_State_0000_00000000.dat',
                         gvec_folder / 'parameter.ini',
                         a_minor_override=aminor(test_folder / 'wout_w7x_kjm.nc'),
                         s_input=s_w7x, nalpha=1, nzgrid=nzgrid, nfield_periods=1)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('Usage: python make_gvec_vs_vmec_files.py <work_folder>')
    work_folder = pathlib.Path(sys.argv[1]).resolve()
    make_tokamak_files(work_folder)
    make_w7x_file(work_folder)
    print(f'Written the VMEC and GVEC files of test_2g in {test_folder}')
