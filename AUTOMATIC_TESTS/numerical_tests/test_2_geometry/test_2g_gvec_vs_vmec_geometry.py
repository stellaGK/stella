################################################################################
#                    Compare the GVEC geometry with the VMEC geometry          #
################################################################################
# Run stella with the VMEC geometry and with the GVEC geometry of the same
# equilibrium, and check that the geometry arrays, the flux-tube constants and
# a short linear ITG simulation agree. Two equilibria are compared:
#
#   tokamak  The circular tokamak of DOCUMENTATION/GVEC_equilibria/circular_tok,
#            solved independently by VMEC++ (wout_tok_circ.nc) and by GVEC
#            (gvec_tok_circ.nc), from the same boundary and profiles.
#   w7x      The W7-X equilibrium wout_w7x_kjm.nc, read by the VMEC geometry,
#            and the same equilibrium converted to GVEC with pygvec convert-wout
#            (gvec_w7x_kjm.nc), which tests the 3-D terms of the geometry.
#
# The GVEC geometry is read from the files of gvec_to_stella.py and, when stella
# is linked to the GVEC library, directly from the GVEC states (these runs are
# skipped otherwise). Both are also run with the magnetic field reversed (B -> -B
# with the same flux surfaces and iota, see gvec_vs_vmec/reverse_magnetic_field.py),
# which must keep |B|, the metric, the Jacobian, |phi|^2(t), the growth rate and
# the heat flux, and reverse b.grad z, the drifts, the Clebsch factor (in
# exb_nonlin_fac and flux_fac) and the real frequency, as for VMEC.
#
# The files are made by gvec_vs_vmec/make_gvec_vs_vmec_files.py, so these tests
# do not need GVEC or VMEC. The tolerances are relative to the largest value of
# each array along the field line. They are set by the differences between the
# two representations of the equilibrium: for the tokamak, both codes agree to
# better than 0.6%; for W7-X, |B| and the metric agree to 0.05%, the binormal
# drifts to 0.5% and the radial drifts to 5.5%, which is a difference of the
# poloidal harmonics of |B| between the VMEC file and its GVEC representation (it
# does not decrease with the radial resolution of GVEC, and the tokamak solved by
# both codes agrees to 0.2%).
#
# Quantities that are odd under a reversal of the handedness of the frame (the
# Clebsch factor, the Jacobian, the radial drifts and grad x . grad y) are
# compared with their signs, since a wrong handedness reverses the drive relative
# to the drifts while |B| and the binormal drifts stay the same; the linear
# simulation of W7-X then differs by ~90%.
################################################################################

# Python modules
import pytest
import os, sys
import shutil
import pathlib
import subprocess
import numpy as np
import netCDF4 as nc

# Package to run stella
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Routines to reverse the magnetic field of a GVEC state
sys.path.insert(0, str(pathlib.Path(__file__).parent / 'gvec_vs_vmec'))
from reverse_magnetic_field import reverse_magnetic_field_of_parameter_file, reverse_magnetic_field_of_state_file

# Global variables
cases = {'tokamak': 'tok_circ', 'w7x': 'w7x_kjm'}
gvec_paths = ['file', 'direct']

# Tolerances of the geometry arrays, relative to the largest value of each array
array_tolerances = {
    'tokamak': {
        'bmag': 1e-2, 'b_dot_gradz': 1e-2, 'gradx_dot_gradx': 1e-2, 'gradx_dot_grady': 1e-2,
        'grady_dot_grady': 1e-2, 'grho': 1e-2, 'jacob': 1e-2,
        'B_times_gradB_dot_grady': 1e-2, 'B_times_kappa_dot_grady': 1e-2,
        'B_times_gradB_dot_gradx': 1e-2, 'B_times_kappa_dot_gradx': 1e-2,
    },
    'w7x': {
        'bmag': 2e-3, 'b_dot_gradz': 2e-3, 'gradx_dot_gradx': 2e-3, 'gradx_dot_grady': 2e-3,
        'grady_dot_grady': 2e-3, 'grho': 2e-3, 'jacob': 2e-3,
        'B_times_gradB_dot_grady': 1e-2, 'B_times_kappa_dot_grady': 1e-2,
        'B_times_gradB_dot_gradx': 8e-2, 'B_times_kappa_dot_gradx': 8e-2,
    },
}

# Arrays that keep or reverse their sign when the magnetic field is reversed
arrays_even_in_B = ['bmag', 'gradx_dot_gradx', 'gradx_dot_grady', 'grady_dot_grady', 'grho', 'jacob']
arrays_odd_in_B = ['b_dot_gradz', 'B_times_gradB_dot_grady', 'B_times_kappa_dot_grady',
                   'B_times_gradB_dot_gradx', 'B_times_kappa_dot_gradx']

# Constants of the flux tube in the header of the .geometry file
constants = ['qinp', 'shat', 'rhotor', 'aref', 'bref', 'dxdpsi', 'dydalpha', 'exb_nonlin', 'flux_fac']
constants_odd_in_B = ['exb_nonlin', 'flux_fac']

# Tolerance of the constants of the flux tube (q, shat, exb_nonlin_fac, ...)
constant_tolerance = 1e-3

# Tolerance of |phi|^2(t) of the linear simulation, relative to |phi|^2
phi2_tolerance = 1e-2

# Tolerance of the comparison of B and -B, which only differ by round-off
reversal_tolerance = 1e-10

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
def run_gvec_direct(case, tmp_path, stella_version, reverse_field):
    '''Run the direct path of the GVEC geometry, which reads the GVEC state, and
    return None if stella has not been linked to the GVEC library.'''
    name = cases[case]
    input_filename = f'gvec_vs_vmec_{case}_gvec_direct.in'
    copy_input_file(input_filename, tmp_path)
    state_file, parameter_file = f'GVEC_{name}_State.dat', f'GVEC_{name}_parameter.ini'
    if reverse_field:
        reverse_magnetic_field_of_state_file(get_stella_expected_run_directory() / state_file, tmp_path / state_file)
        reverse_magnetic_field_of_parameter_file(get_stella_expected_run_directory() / parameter_file, tmp_path / parameter_file)
    else:
        shutil.copyfile(get_stella_expected_run_directory() / state_file, tmp_path / state_file)
        shutil.copyfile(get_stella_expected_run_directory() / parameter_file, tmp_path / parameter_file)
    os.chdir(tmp_path)
    run = subprocess.run(['mpirun', '--oversubscribe', '-np', f'{read_nproc()}', get_stella_path(stella_version), input_filename],
                         capture_output=True, text=True)
    if run.returncode != 0:
        if 'compiled without GVEC support' in run.stdout + run.stderr: return None
        raise RuntimeError(f'stella stopped with the direct path of the GVEC geometry:\n{run.stdout}\n{run.stderr}')
    return tmp_path / input_filename.replace('.in', '')

#-------------------------------------------------------------------------------
#                         Run local stella simulations                         #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def stella_runs(tmp_path_factory, stella_version):
    '''Run stella with the VMEC geometry and with both paths of the GVEC geometry,
    with B and -B, of both equilibria; the runs are shared by all the tests in this
    module. The key is (case, geometry, magnetic field), the value the path of the
    output files without extension, or None for a direct path that is not available.'''
    if stella_version!='master': pytest.skip('The GVEC geometry has not been implemented for older stella versions.')
    runs = {}
    for case, name in cases.items():
        input_filename = f'gvec_vs_vmec_{case}_vmec.in'
        tmp_path = tmp_path_factory.mktemp(f'gvec_vs_vmec_{case}_vmec')
        run_local_stella_simulation(input_filename, tmp_path, stella_version, f'wout_{name}.nc')
        runs[(case, 'vmec', '+B')] = tmp_path / input_filename.replace('.in', '')
        for field, suffix in [('+B', ''), ('-B', '_Bminus')]:
            input_filename = f'gvec_vs_vmec_{case}_gvec.in'
            tmp_path = tmp_path_factory.mktemp(f'gvec_vs_vmec_{case}_gvec{suffix}')
            geometry_file = f'gvec_{name}{suffix}.nc'
            copy_input_file(input_filename, tmp_path)
            with open(tmp_path / input_filename) as file: text = file.read()
            with open(tmp_path / input_filename, 'w') as file: file.write(text.replace(f"'gvec_{name}.nc'", f"'{geometry_file}'"))
            copy_vmec_file(geometry_file, tmp_path)
            os.chdir(tmp_path)
            run_stella(get_stella_path(stella_version), input_filename)
            runs[(case, 'file', field)] = tmp_path / input_filename.replace('.in', '')
            tmp_path = tmp_path_factory.mktemp(f'gvec_vs_vmec_{case}_gvec_direct{suffix}')
            runs[(case, 'direct', field)] = run_gvec_direct(case, tmp_path, stella_version, field=='-B')
    return runs

#-------------------------------------------------------------------------------
def get_run(stella_runs, case, geometry, field='+B'):
    run = stella_runs[(case, geometry, field)]
    if run is None: pytest.skip('stella has not been linked to the GVEC library (USE_GVEC), so the direct path is not tested.')
    return run

#-------------------------------------------------------------------------------
def read_flux_tube_constants(path):
    '''Read the constants of the flux tube from the header of the .geometry file.'''
    with open(str(path) + '.geometry') as geometry_file:
        names = geometry_file.readline().lstrip('#').split()
        values = [float(value) for value in geometry_file.readline().lstrip('#').split()]
    return dict(zip(names, values))

#-------------------------------------------------------------------------------
def read_netcdf(path, name):
    with nc.Dataset(str(path) + '.out.nc') as data:
        return np.asarray(data[name][...], dtype=float)

#-------------------------------------------------------------------------------
def relative_difference(array, reference):
    return np.max(np.abs(array - reference)) / np.max(np.abs(reference))

#-------------------------------------------------------------------------------
#                      Check the geometry arrays along z                       #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('gvec_path', gvec_paths)
@pytest.mark.parametrize('case', list(cases))
def test_whether_gvec_and_vmec_geometry_arrays_agree(stella_runs, case, gvec_path, error=False):
    vmec, gvec = get_run(stella_runs, case, 'vmec'), get_run(stella_runs, case, gvec_path)
    for name, tolerance in array_tolerances[case].items():
        difference = relative_difference(read_netcdf(gvec, name), read_netcdf(vmec, name))
        if difference > tolerance:
            print(f'ERROR: {name} of the GVEC geometry ({gvec_path}) differs from the VMEC geometry by {difference:.2e} (tolerance {tolerance:.0e}).'); error = True
        else:
            print(f'  {name}: relative difference {difference:.2e} (tolerance {tolerance:.0e})')
    assert (not error), f'The GVEC ({gvec_path}) and VMEC geometry arrays of the {case} do not agree.'
    print(f'  -->  The GVEC ({gvec_path}) and VMEC geometry arrays of the {case} agree.')
    return

#-------------------------------------------------------------------------------
#                    Check the constants of the flux tube                      #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('gvec_path', gvec_paths)
@pytest.mark.parametrize('case', list(cases))
def test_whether_gvec_and_vmec_flux_tube_constants_agree(stella_runs, case, gvec_path, error=False):
    vmec = read_flux_tube_constants(get_run(stella_runs, case, 'vmec'))
    gvec = read_flux_tube_constants(get_run(stella_runs, case, gvec_path))
    for name in constants:
        if abs(gvec[name] - vmec[name]) > constant_tolerance * abs(vmec[name]):
            print(f'ERROR: {name} = {gvec[name]} for the GVEC geometry ({gvec_path}), but {vmec[name]} for the VMEC geometry.'); error = True
    assert (not error), f'The GVEC ({gvec_path}) and VMEC constants of the flux tube of the {case} do not agree.'
    print(f'  -->  The GVEC ({gvec_path}) and VMEC constants of the flux tube of the {case} agree.')
    return

#-------------------------------------------------------------------------------
#                      Check the linear ITG simulation                         #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('gvec_path', gvec_paths)
@pytest.mark.parametrize('case', list(cases))
def test_whether_gvec_and_vmec_linear_simulations_agree(stella_runs, case, gvec_path):
    vmec, gvec = get_run(stella_runs, case, 'vmec'), get_run(stella_runs, case, gvec_path)
    vmec_phi2, gvec_phi2 = read_netcdf(vmec, 'phi2'), read_netcdf(gvec, 'phi2')
    difference = np.max(np.abs(gvec_phi2 - vmec_phi2) / np.abs(vmec_phi2))
    assert difference < phi2_tolerance, f'|phi|^2(t) of the {case} differs by {difference:.2e} between the GVEC ({gvec_path}) and VMEC geometries (tolerance {phi2_tolerance:.0e}).'
    print(f'  -->  |phi|^2(t) of the {case} agrees to {difference:.2e} between the GVEC ({gvec_path}) and VMEC geometries.')
    return

#-------------------------------------------------------------------------------
#                    Check the reversal of the magnetic field                  #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('gvec_path', gvec_paths)
@pytest.mark.parametrize('case', list(cases))
def test_whether_gvec_geometry_handles_the_sign_of_the_magnetic_field(stella_runs, case, gvec_path, error=False):
    plus, minus = get_run(stella_runs, case, gvec_path, '+B'), get_run(stella_runs, case, gvec_path, '-B')

    # Geometry arrays
    for name in arrays_even_in_B + arrays_odd_in_B:
        sign = -1 if name in arrays_odd_in_B else 1
        difference = relative_difference(sign * read_netcdf(minus, name), read_netcdf(plus, name))
        if difference > reversal_tolerance:
            print(f'ERROR: {name} with -B is not {"-" if sign<0 else ""}{name} with +B (relative difference {difference:.2e}).'); error = True

    # Constants of the flux tube
    constants_plus, constants_minus = read_flux_tube_constants(plus), read_flux_tube_constants(minus)
    for name in constants:
        sign = -1 if name in constants_odd_in_B else 1
        if abs(sign * constants_minus[name] - constants_plus[name]) > constant_tolerance * abs(constants_plus[name]):
            print(f'ERROR: {name} = {constants_minus[name]} with -B, but {constants_plus[name]} with +B.'); error = True

    # Linear simulation: same |phi|^2(t), growth rate and heat flux, opposite real frequency
    for name, sign in [('phi2', 1), ('qflux_vs_s', 1)]:
        difference = relative_difference(read_netcdf(minus, name), read_netcdf(plus, name))
        if difference > reversal_tolerance:
            print(f'ERROR: {name} with -B differs from {name} with +B by {difference:.2e}.'); error = True
    omega_plus, omega_minus = read_netcdf(plus, 'omega')[-1].ravel(), read_netcdf(minus, 'omega')[-1].ravel()
    if np.max(np.abs(omega_minus * np.array([-1, 1]) - omega_plus)) > reversal_tolerance * np.max(np.abs(omega_plus)):
        print(f'ERROR: omega = {omega_minus} with -B, but {omega_plus} with +B (expected the opposite real frequency).'); error = True

    assert (not error), f'The GVEC geometry ({gvec_path}) of the {case} does not handle the reversal of the magnetic field.'
    print(f'  -->  The GVEC geometry ({gvec_path}) of the {case} handles the reversal of the magnetic field.')
    return
