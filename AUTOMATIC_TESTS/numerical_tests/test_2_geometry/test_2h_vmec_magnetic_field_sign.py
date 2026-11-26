################################################################################
#              Check the sign of the magnetic field in the VMEC geometry       #
################################################################################
# Run stella with two VMEC equilibria that only differ by the sign of the
# magnetic field (B+ and B-, the same flux surfaces, iota and pressure), for a
# circular tokamak and for W7-X. Reversing B must keep |B|, the metric, the
# Jacobian, |phi|^2(t), the growth rate, the heat flux and the toroidal angular
# momentum flux of a linear ITG simulation, and reverse b.grad z, the magnetic
# drifts, the Clebsch factor (in exb_nonlin_fac and flux_fac) and the real
# frequency, since the mode then propagates in the opposite direction of y.
#
# The equilibria are made by vmec_magnetic_field_sign/make_vmec_files.py with
# VMEC++, so the B- file is the B+ file with the field reversed, and the two
# simulations agree to round-off.
#
# The B- solution is the complex conjugate of the B+ solution at the opposite
# vpa, so both the parallel part (vpa b.grad(zeta) R^2) and the perpendicular
# part ((b x k).grad(zeta) R^2/B) of the toroidal angular momentum flux keep
# their sign, since both contain the unit vector b. This agrees with the
# symmetry of the momentum flux in tokamaks, which changes sign when either the
# toroidal field or the plasma current is reversed, but not when both are
# (Camenen et al., Phys. Plasmas 16, 062501 (2009)). Before November 2026 the
# perpendicular part missed the sign of b, so that the momentum flux changed
# with the sign of B (6.2e-8 with B+ and 4.5e-9 with B- for the tokamak).
################################################################################

# Python modules
import pytest
import os, sys
import pathlib
import numpy as np
import netCDF4 as nc

# Package to run stella
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Global variables
cases = {'tokamak': 'tok_circ', 'w7x': 'w7x_kjm'}

# Arrays that keep or reverse their sign when the magnetic field is reversed
arrays_even_in_B = ['bmag', 'gradx_dot_gradx', 'gradx_dot_grady', 'grady_dot_grady', 'grho', 'jacob']
arrays_odd_in_B = ['b_dot_gradz', 'B_times_gradB_dot_grady', 'B_times_kappa_dot_grady',
                   'B_times_gradB_dot_gradx', 'B_times_kappa_dot_gradx']

# Constants of the flux tube in the header of the .geometry file
constants = ['qinp', 'shat', 'rhotor', 'aref', 'bref', 'dxdpsi', 'dydalpha', 'exb_nonlin', 'flux_fac']
constants_odd_in_B = ['exb_nonlin', 'flux_fac']

# The B+ and B- simulations only differ by round-off
tolerance = 1e-10

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                         Run local stella simulations                         #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def stella_runs(tmp_path_factory, stella_version):
    '''Run stella with the B+ and the B- equilibrium of both cases; the runs are
    shared by all the tests in this module.'''
    if stella_version!='master': pytest.skip('This test has not been implemented for older stella versions.')
    runs = {}
    for case, name in cases.items():
        input_filename = f'vmec_magnetic_field_sign_{case}.in'
        for field in ['Bplus', 'Bminus']:
            tmp_path = tmp_path_factory.mktemp(f'vmec_magnetic_field_sign_{case}_{field}')
            vmec_filename = f'wout_{name}_{field}.nc'
            copy_input_file(input_filename, tmp_path)
            with open(tmp_path / input_filename) as file: text = file.read()
            with open(tmp_path / input_filename, 'w') as file: file.write(text.replace(f"'wout_{name}_Bplus.nc'", f"'{vmec_filename}'"))
            copy_vmec_file(vmec_filename, tmp_path)
            os.chdir(tmp_path)
            run_stella(get_stella_path(stella_version), input_filename)
            runs[(case, field)] = tmp_path / input_filename.replace('.in', '')
    return runs

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
#                 Check the geometry arrays and flux-tube constants            #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('case', list(cases))
def test_whether_vmec_geometry_handles_the_sign_of_the_magnetic_field(stella_runs, case, error=False):
    plus, minus = stella_runs[(case, 'Bplus')], stella_runs[(case, 'Bminus')]
    for name in arrays_even_in_B + arrays_odd_in_B:
        sign = -1 if name in arrays_odd_in_B else 1
        difference = relative_difference(sign * read_netcdf(minus, name), read_netcdf(plus, name))
        if difference > tolerance:
            print(f'ERROR: {name} with -B is not {"-" if sign<0 else ""}{name} with +B (relative difference {difference:.2e}).'); error = True
    constants_plus, constants_minus = read_flux_tube_constants(plus), read_flux_tube_constants(minus)
    for name in constants:
        sign = -1 if name in constants_odd_in_B else 1
        if abs(sign * constants_minus[name] - constants_plus[name]) > 1e-3 * abs(constants_plus[name]):
            print(f'ERROR: {name} = {constants_minus[name]} with -B, but {constants_plus[name]} with +B.'); error = True
    assert (not error), f'The VMEC geometry of the {case} does not handle the reversal of the magnetic field.'
    print(f'  -->  The VMEC geometry of the {case} handles the reversal of the magnetic field.')
    return

#-------------------------------------------------------------------------------
#                      Check the linear ITG simulation                         #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('case', list(cases))
def test_whether_vmec_linear_simulation_handles_the_sign_of_the_magnetic_field(stella_runs, case, error=False):
    plus, minus = stella_runs[(case, 'Bplus')], stella_runs[(case, 'Bminus')]
    for name in ['phi2', 'qflux_vs_s']:
        difference = relative_difference(read_netcdf(minus, name), read_netcdf(plus, name))
        if difference > tolerance:
            print(f'ERROR: {name} with -B differs from {name} with +B by {difference:.2e}.'); error = True
    omega_plus, omega_minus = read_netcdf(plus, 'omega')[-1].ravel(), read_netcdf(minus, 'omega')[-1].ravel()
    if np.max(np.abs(omega_minus * np.array([-1, 1]) - omega_plus)) > tolerance * np.max(np.abs(omega_plus)):
        print(f'ERROR: omega = {omega_minus} with -B, but {omega_plus} with +B (expected the opposite real frequency).'); error = True
    assert (not error), f'The linear simulation of the {case} does not handle the reversal of the magnetic field.'
    print(f'  -->  The linear simulation of the {case} handles the reversal of the magnetic field.')
    return

#-------------------------------------------------------------------------------
#                            Check the momentum flux                           #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('case', list(cases))
def test_whether_vmec_momentum_flux_handles_the_sign_of_the_magnetic_field(stella_runs, case):
    plus, minus = stella_runs[(case, 'Bplus')], stella_runs[(case, 'Bminus')]
    vflux_plus, vflux_minus = read_netcdf(plus, 'vflux_vs_s'), read_netcdf(minus, 'vflux_vs_s')
    qflux_plus = read_netcdf(plus, 'qflux_vs_s')
    assert np.max(np.abs(vflux_plus)) > 1e-3 * np.max(np.abs(qflux_plus)), f'The momentum flux of the {case} is too small to test its sign.'
    difference = relative_difference(vflux_minus, vflux_plus)
    assert difference < tolerance, f'vflux of the {case} changes by {difference:.2e} when the magnetic field is reversed.'
    print(f'  -->  The momentum flux of the {case} does not change when the magnetic field is reversed.')
    return
