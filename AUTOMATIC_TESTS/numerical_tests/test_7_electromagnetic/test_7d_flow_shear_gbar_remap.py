################################################################################
#          Check the Hammett remap of gbar in electromagnetic flow shear       #
################################################################################
# With perpendicular flow shear (g_exb /= 0) the Hammett remap moves each k_y
# by one k_x cell at a time. In electromagnetic runs the option
#       &flow_shear
#         remap_em = 0 / 1
#       /
# chooses whether the remap moves g (0, the original) or gbar = g + c J0 A_par
# (1), with A_par solved again from gbar afterwards. These tests check that
#   1. remap_em = 1 reproduces the stored time traces of phi, A_par and B_par;
#   2. remap_em = 1 changes A_par relative to remap_em = 0, so the option acts;
#   3. without A_par (gbar = g) both options give the same result.
# See DOCUMENTATION/flow_shear_gbar_remap.md for why gbar is moved.
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
input_filename = 'EM_flow_shear_remap_gbar.in'

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#        Run a copy of the test input file with some lines replaced           #
#-------------------------------------------------------------------------------
def run_modified_input_file(tmp_path, stella_version, name, replacements):
    text = (get_stella_expected_run_directory() / input_filename).read_text()
    for old, new in replacements:
        assert text.count(old) == 1, f'"{old}" should appear once in {input_filename}'
        text = text.replace(old, new)
    (tmp_path / name).write_text(text)
    os.chdir(tmp_path)
    run_stella(get_stella_path(stella_version), name)
    return xr.open_dataset(tmp_path / name.replace('.in', '.out.nc'))

#-------------------------------------------------------------------------------
#      Check that the gbar remap reproduces the stored phi, A_par, B_par       #
#-------------------------------------------------------------------------------
def test_flow_shear_gbar_remap(tmp_path, stella_version):

    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path, stella_version)

    # Compare the time traces of |phi|^2, |A_par|^2 and |B_par|^2 with the expected output
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc')
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        for name in ['t', 'phi2', 'apar2', 'bpar2']:
            assert np.allclose(local_netcdf[name], expected_netcdf[name], equal_nan=True, atol=1e-20), \
                f'<{name}> does not match the expected output for remap_em = 1.'

    print('  -->  The Hammett remap of gbar reproduces the expected phi, A_par and B_par.')
    return

#-------------------------------------------------------------------------------
#         Check that remap_em = 1 changes A_par relative to remap_em = 0       #
#-------------------------------------------------------------------------------
def test_flow_shear_gbar_remap_differs_from_g_remap(tmp_path, stella_version):

    gbar = run_modified_input_file(tmp_path, stella_version, 'remap_gbar.in', [])
    g = run_modified_input_file(tmp_path, stella_version, 'remap_g.in', [('  remap_em = 1\n', '  remap_em = 0\n')])

    # The two remaps solve A_par differently at every remap, so A_par must differ
    apar2_gbar, apar2_g = gbar['apar2'].values, g['apar2'].values
    relative_difference = np.max(np.abs(apar2_gbar - apar2_g)) / np.max(np.abs(apar2_g))
    assert relative_difference > 1e-2, \
        f'remap_em = 1 should change A_par relative to remap_em = 0 (relative difference {relative_difference:.1e}).'

    print(f'  -->  remap_em = 1 changes |A_par|^2 by up to {relative_difference:.0%} relative to remap_em = 0.')
    return

#-------------------------------------------------------------------------------
#         Check that both remaps agree when there is no A_par (gbar = g)       #
#-------------------------------------------------------------------------------
def test_flow_shear_gbar_remap_electrostatic(tmp_path, stella_version):

    electrostatic = [('include_apar = .true.', 'include_apar = .false.'),
                     ('include_bpar = .true.', 'include_bpar = .false.')]
    gbar = run_modified_input_file(tmp_path, stella_version, 'ES_remap_gbar.in', electrostatic)
    g = run_modified_input_file(tmp_path, stella_version, 'ES_remap_g.in', electrostatic + [('  remap_em = 1\n', '  remap_em = 0\n')])

    # Without A_par, gbar = g and both options take the same code path
    assert np.array_equal(gbar['phi2'].values, g['phi2'].values), \
        'Without A_par, remap_em = 1 and remap_em = 0 should give identical results.'

    print('  -->  Without A_par, remap_em = 1 and remap_em = 0 give identical results.')
    return
