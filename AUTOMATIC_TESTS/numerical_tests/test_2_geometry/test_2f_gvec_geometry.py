################################################################################
#                            Check the GVEC geometry                           #
################################################################################
# Tests for <geometry_option> = 'gvec', which reads a pre-computed geometry
# netCDF file written by gvec_to_stella.py.
#
# These tests deliberately do NOT require GVEC (or pygvec) to be installed.
# Instead they synthesise a geometry file whose contents are known analytically,
# so that every geometry array stella reports can be predicted in closed form.
# That makes the tests fast, hermetic, and - most importantly - sensitive to the
# index/ordering bugs that netCDF readers are prone to.
#
# Why the synthetic file matters
# ------------------------------
# gvec_to_stella.py declares its 2-D variables with dimensions
# ("nalpha", "nz_file").  The Fortran netCDF API presents dimensions in the
# REVERSE order, so the receiving buffer must be (nz_file, nalpha).  Getting
# this wrong makes every 2-D read fail - including for nalpha = 1 - which is
# exactly the failure these tests are designed to catch.
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
input_filename = 'gvec_geometry.in'
input_file = input_filename.replace('.in', '')
geometry_filename = 'gvec_geometry.nc'
stella_local_run_directory = 'Not/Run/Yet'
run_data = {}

# Parameters of the synthetic equilibrium. These must stay consistent with
# <gvec_geometry.in>: nzed = 2*NZGRID and nalpha = 1 (flux tube).
S_TORFLUX     = 0.25          # normalised toroidal flux -> rho = 0.5
RHO           = np.sqrt(S_TORFLUX)
NZGRID        = 8             # nzed = 16
NALPHA        = 1
NFP           = 5
SIGN_TORFLUX  = 1
SAFETY_Q      = 2.0
SHAT          = 0.8
L_REFERENCE   = 1.0
B_REFERENCE   = 1.0

# The .geometry text file is written with 4 significant figures, so a relative
# tolerance of 2e-3 is the tightest meaningful comparison. Any index-ordering or
# rho-scaling error produces an O(1) discrepancy, far above this.
RTOL = 2e-3


#-------------------------------------------------------------------------------
#                       Build a synthetic GVEC geometry file                   #
#-------------------------------------------------------------------------------
def analytic_arrays(nzgrid=NZGRID, nalpha=NALPHA, constant_b_dot_gradzeta=None):
    '''Analytic geometry arrays on the (alpha, zeta) field-line grid.

    The zeta-dependence is deliberately non-symmetric and different for each
    quantity, so that a transposed, reversed or mis-shifted read is visible.
    '''
    nz_file = 2 * nzgrid + 1
    zeta = np.linspace(-np.pi, np.pi, nz_file)
    alpha = np.zeros(nalpha)

    # Broadcast to (nalpha, nz_file)
    z = zeta[np.newaxis, :] * np.ones((nalpha, 1))

    arrays = {
        # |B| / B_ref : strictly positive, single well
        'bmag':                  1.0 + 0.20 * np.cos(z),
        # metric coefficients, chosen positive-definite and Cauchy-Schwarz safe
        'grad_alpha_grad_alpha': 2.0 + 0.50 * np.sin(z),
        'grad_psit_grad_psit':   3.0 + 0.70 * np.cos(2.0 * z),
        'grad_alpha_grad_psit':  0.30 * np.sin(z),
        # remaining quantities: smooth, small, physically harmless
        'theta':                 z.copy(),
        'grho':                  1.0 + 0.05 * np.sin(z),
        'b_dot_gradzeta':        (0.50 + 0.05 * np.cos(z)
                                  if constant_b_dot_gradzeta is None
                                  else constant_b_dot_gradzeta * np.ones_like(z)),
        'B_sub_theta':           0.10 * np.ones_like(z),
        'B_sub_zeta':            1.00 * np.ones_like(z),
        'gds23':                 0.02 * np.sin(z),
        'gds24':                 0.03 * np.cos(z),
        'gds25':                 0.04 * np.sin(2.0 * z),
        'gds26':                 0.05 * np.cos(2.0 * z),
        'gbdrift_alpha':              0.10 * np.sin(z),
        'B_times_gradB_dot_gradpsit': 0.11 * np.cos(z),
        'gbdrift_kappa_alpha':        0.12 * np.sin(z),
        'B_times_kappa_dot_gradpsit': 0.13 * np.cos(z),
        'psit_displacement_fac':      0.14 * np.ones_like(z),
        'gradzeta_gradalpha_R2overB2': 0.15 * np.sin(z),
        'gradzeta_gradpsit_R2overB2':  0.16 * np.cos(z),
        'b_dot_gradzeta_RR':           0.17 * np.ones_like(z),
    }
    return alpha, zeta, arrays


def write_synthetic_gvec_file(path, nzgrid=NZGRID, nalpha=NALPHA,
                              dimension_order=('nalpha', 'nz_file'),
                              constant_b_dot_gradzeta=None):
    '''Write a GVEC geometry netCDF file in the format gvec_to_stella.py produces.

    <dimension_order> is exposed so a test can deliberately write the file with
    the wrong dimension order and confirm stella rejects it rather than
    silently reading scrambled data.
    '''
    alpha, zeta, arrays = analytic_arrays(nzgrid, nalpha, constant_b_dot_gradzeta)
    nz_file = 2 * nzgrid + 1

    ds = nc.Dataset(str(path), 'w')
    try:
        ds.createDimension('nalpha', nalpha)
        ds.createDimension('nz_file', nz_file)

        for name, value in [('nfp', NFP), ('sign_torflux', SIGN_TORFLUX),
                            ('nfp_sim', NFP)]:
            ds.createVariable(name, 'i4')[...] = value

        for name, value in [('s', S_TORFLUX), ('safety_factor_q', SAFETY_Q),
                            ('shat', SHAT), ('L_reference', L_REFERENCE),
                            ('B_reference', B_REFERENCE), ('dzdzeta', 1.0),
                            ('alpha0', 0.0), ('zeta_center', 0.0)]:
            ds.createVariable(name, 'f8')[...] = value

        ds.createVariable('alpha', 'f8', ('nalpha',))[:] = alpha
        ds.createVariable('zeta', 'f8', ('nz_file',))[:] = zeta

        for name, arr in arrays.items():
            var = ds.createVariable(name, 'f8', dimension_order)
            var[:] = arr if dimension_order == ('nalpha', 'nz_file') else arr.T
    finally:
        ds.close()
    return path


#-------------------------------------------------------------------------------
#                          Read stella's .geometry file                        #
#-------------------------------------------------------------------------------
GEOMETRY_COLUMNS = ['alpha', 'zed', 'zeta', 'bmag', 'b_dot_gradz',
                    'grady_dot_grady', 'gradx_dot_grady', 'gradx_dot_gradx',
                    'BxgradB_dot_grady', 'Bxkappa_dot_grady',
                    'BxgradB_dot_gradx', 'bmag_psi0']


def read_geometry_file(path):
    '''Parse the <input>.geometry text file into a dict of numpy arrays.'''
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            rows.append([float(value) for value in line.split()])
    data = np.array(rows)
    return {name: data[:, i] for i, name in enumerate(GEOMETRY_COLUMNS)}


#-------------------------------------------------------------------------------
#                    Check whether output files are present                    #
#-------------------------------------------------------------------------------
def test_whether_gvec_output_files_are_present(tmp_path, error=False):
    '''Run stella on a synthetic GVEC geometry file and check it produced output.'''

    global stella_local_run_directory, run_data
    stella_local_run_directory = tmp_path

    # The geometry file must exist in the run directory before stella starts
    write_synthetic_gvec_file(tmp_path / geometry_filename)

    run_data = run_local_stella_simulation(input_filename, tmp_path, 'master')

    local_files = os.listdir(stella_local_run_directory)
    expected_files = [f'{input_file}.geometry', f'{input_file}.out.nc']
    for expected_file in expected_files:
        if expected_file not in local_files:
            print(f'ERROR: The "{expected_file}" output file was not generated when running stella.')
            error = True

    assert not error, f'Some output files were not generated when running stella with GVEC geometry.'
    print('\n  -->  stella ran successfully with geometry_option = "gvec".')


#-------------------------------------------------------------------------------
#            The geometry arrays must match the analytic prediction            #
#-------------------------------------------------------------------------------
def test_whether_gvec_geometry_arrays_match_analytic_values(error=False):
    '''Compare stella's geometry arrays against closed-form expectations.

    For radial_coordinate = 'sgnpsitpsit' (the default) and sign_torflux = +1:

        dxdpsi   = 1/rho          dydalpha = rho
        gradx_dot_gradx = grad_psit_grad_psit / rho^2
        grady_dot_grady = grad_alpha_grad_alpha * rho^2
        gradx_dot_grady = sign_torflux * grad_alpha_grad_psit
        bmag            = bmag (passed through unchanged)

    A transposed or mis-scaled read breaks these at the O(1) level.
    '''
    geometry = read_geometry_file(stella_local_run_directory / f'{input_file}.geometry')
    _, _, arrays = analytic_arrays()

    expected = {
        'bmag':            arrays['bmag'][0],
        'gradx_dot_gradx': arrays['grad_psit_grad_psit'][0] / RHO**2,
        'grady_dot_grady': arrays['grad_alpha_grad_alpha'][0] * RHO**2,
        'gradx_dot_grady': SIGN_TORFLUX * arrays['grad_alpha_grad_psit'][0],
    }

    for name, want in expected.items():
        got = geometry[name]
        if got.shape != want.shape:
            print(f'ERROR: "{name}" has shape {got.shape}, expected {want.shape}.')
            error = True
            continue
        if not np.allclose(got, want, rtol=RTOL, atol=1e-6):
            worst = int(np.argmax(np.abs(got - want)))
            print(f'ERROR: "{name}" does not match the analytic value.')
            print(f'       largest deviation at z index {worst}: '
                  f'stella = {got[worst]:.6e}, expected = {want[worst]:.6e}')
            error = True

    assert not error, 'The GVEC geometry arrays do not match their analytic values.'
    print('\n  -->  The GVEC geometry arrays match the analytic values.')


#-------------------------------------------------------------------------------
#                 The z-grid must be read in the correct order                 #
#-------------------------------------------------------------------------------
def test_whether_gvec_zeta_grid_is_read_in_order(error=False):
    '''The zeta grid must come back monotonically increasing over [-pi, pi].

    This catches a reversed or rolled z-index, which a symmetric test array
    would hide.
    '''
    geometry = read_geometry_file(stella_local_run_directory / f'{input_file}.geometry')
    _, zeta, _ = analytic_arrays()

    if not np.allclose(geometry['zeta'], zeta, rtol=RTOL, atol=1e-6):
        print(f'ERROR: the zeta grid does not match the file.')
        print(f'       stella   = {geometry["zeta"][:5]} ...')
        print(f'       expected = {zeta[:5]} ...')
        error = True

    assert not error, 'The GVEC zeta grid was not read in the correct order.'
    print('\n  -->  The GVEC zeta grid is read in the correct order.')


#-------------------------------------------------------------------------------
#                   Physical invariants of the geometry arrays                 #
#-------------------------------------------------------------------------------
def test_whether_gvec_geometry_satisfies_physical_invariants(error=False):
    '''|B| > 0, the metric is positive definite, and Cauchy-Schwarz holds.

    These hold for any valid equilibrium, so this test also applies to geometry
    files generated from real GVEC states.
    '''
    geometry = read_geometry_file(stella_local_run_directory / f'{input_file}.geometry')

    if not np.all(geometry['bmag'] > 0):
        print('ERROR: bmag is not strictly positive.'); error = True
    if not np.all(geometry['gradx_dot_gradx'] > 0):
        print('ERROR: gradx_dot_gradx is not strictly positive.'); error = True
    if not np.all(geometry['grady_dot_grady'] > 0):
        print('ERROR: grady_dot_grady is not strictly positive.'); error = True

    # Cauchy-Schwarz: (grad x . grad y)^2 <= |grad x|^2 |grad y|^2
    lhs = geometry['gradx_dot_grady']**2
    rhs = geometry['gradx_dot_gradx'] * geometry['grady_dot_grady']
    if not np.all(lhs <= rhs * (1.0 + RTOL)):
        worst = int(np.argmax(lhs - rhs))
        print(f'ERROR: Cauchy-Schwarz violated at z index {worst}: '
              f'{lhs[worst]:.6e} > {rhs[worst]:.6e}')
        error = True

    assert not error, 'The GVEC geometry arrays violate a physical invariant.'
    print('\n  -->  The GVEC geometry satisfies the physical invariants.')


#-------------------------------------------------------------------------------
#              A mismatched grid in the file must be rejected                  #
#-------------------------------------------------------------------------------
def test_whether_gvec_rejects_a_mismatched_grid(tmp_path, error=False):
    '''stella must abort - not silently continue - if the file grid is wrong.

    The reader compares the file's <nalpha> and <nz_file> dimensions against
    what stella asked for. Writing a file with the wrong nzgrid must therefore
    stop the run.
    '''
    write_synthetic_gvec_file(tmp_path / geometry_filename, nzgrid=NZGRID + 2)

    failed_as_expected = False
    try:
        run_local_stella_simulation(input_filename, tmp_path, 'master')
    except Exception:
        failed_as_expected = True

    if not failed_as_expected:
        print('ERROR: stella accepted a GVEC geometry file with the wrong z-grid size.')
        error = True

    assert not error, 'stella did not reject a mismatched GVEC geometry grid.'
    print('\n  -->  stella correctly rejects a mismatched GVEC geometry grid.')


#-------------------------------------------------------------------------------
#                 A missing geometry file must abort cleanly                   #
#-------------------------------------------------------------------------------
def test_whether_gvec_rejects_a_missing_file(tmp_path, error=False):
    '''A missing geometry file must abort with an error, not run with zeros.'''

    # Deliberately do not create <gvec_geometry.nc> in tmp_path
    failed_as_expected = False
    try:
        run_local_stella_simulation(input_filename, tmp_path, 'master')
    except Exception:
        failed_as_expected = True

    if not failed_as_expected:
        print('ERROR: stella did not abort when the GVEC geometry file was missing.')
        error = True

    assert not error, 'stella did not abort on a missing GVEC geometry file.'
    print('\n  -->  stella correctly aborts when the GVEC geometry file is missing.')


#-------------------------------------------------------------------------------
#              The arc-length z-grid must produce finite geometry              #
#-------------------------------------------------------------------------------
ARC_REFINEMENT = 4
arc_input_filename = 'gvec_geometry_arc.in'
arc_input_file = arc_input_filename.replace('.in', '')


def test_whether_gvec_arc_length_grid_is_finite(tmp_path, error=False):
    '''Run with zed_equal_arc = .true. and check the geometry is correct.

    With an arc-length z-grid stella re-interpolates every array from the fine
    (nzgrid_file) grid onto the arc-length grid - a different code path from the
    direct copy the other tests exercise, whose failure mode is NaN rather than
    an error: stella still exits 0, so nothing catches it unless values are
    checked.

    To make the test sharp we give the file a CONSTANT b.grad zeta. The arc
    length is then exactly proportional to zeta, so after stella compresses it
    onto [-pi, pi] the arc-length grid coincides with the zeta grid, and every
    interpolated array must reproduce the analytic function pointwise. A
    mis-scaled arc-length grid samples the data at the wrong positions and shows
    up here even though the values stay finite and in range.
    '''
    b_dot_gradzeta_const = 0.5      # != 1 so a missing rescaling is visible
    write_synthetic_gvec_file(tmp_path / geometry_filename,
                              nzgrid=NZGRID * ARC_REFINEMENT,
                              constant_b_dot_gradzeta=b_dot_gradzeta_const)

    run_local_stella_simulation(arc_input_filename, tmp_path, 'master')
    geometry = read_geometry_file(tmp_path / f'{arc_input_file}.geometry')

    for name in ['bmag', 'gradx_dot_gradx', 'grady_dot_grady',
                 'gradx_dot_grady', 'b_dot_gradz']:
        if not np.all(np.isfinite(geometry[name])):
            n_bad = int(np.sum(~np.isfinite(geometry[name])))
            print(f'ERROR: "{name}" has {n_bad} non-finite entries on the '
                  f'arc-length grid.')
            error = True

    if error:
        assert not error, 'The GVEC arc-length z-grid produced non-finite geometry.'

    # With constant b.grad zeta the arc-length grid equals the zeta grid, so the
    # analytic values evaluated on stella's zed grid are the expected answer.
    zed = geometry['zed']
    expected_bmag = 1.0 + 0.20 * np.cos(zed)
    expected_gxx = (3.0 + 0.70 * np.cos(2.0 * zed)) / RHO**2
    expected_gyy = (2.0 + 0.50 * np.sin(zed)) * RHO**2

    for name, want in [('bmag', expected_bmag),
                       ('gradx_dot_gradx', expected_gxx),
                       ('grady_dot_grady', expected_gyy)]:
        got = geometry[name]
        # Spline interpolation onto a coincident grid: allow a little more than
        # the 4-significant-figure output tolerance, but far below an O(1) error.
        if not np.allclose(got, want, rtol=5e-3, atol=1e-5):
            worst = int(np.argmax(np.abs(got - want)))
            print(f'ERROR: "{name}" on the arc-length grid does not reproduce '
                  f'the analytic value.')
            print(f'       worst at zed = {zed[worst]:.4f}: '
                  f'stella = {got[worst]:.6e}, expected = {want[worst]:.6e}')
            error = True

    if not np.all(geometry['bmag'] > 0):
        print('ERROR: bmag is not strictly positive on the arc-length grid.')
        error = True

    assert not error, 'The GVEC arc-length z-grid produced incorrect geometry.'
    print('\n  -->  The GVEC arc-length z-grid reproduces the analytic geometry.')
