################################################################################
#                       Check the sheared slab geometry                        #
################################################################################
# Test all the geometry arrays when using a sheared slab equilibrium.
#
# Analytic reference values (shat=0.5, qinp=1.4, rmaj=3.0, betaprim=0):
#   b.gradz  = 1/(q*R) = 1/4.2 ≈ 0.2381
#   |∇y|²    = 1 + (shat*z)²   (z-dependent)
#   ∇x.∇y    = -shat*z          (z-dependent)
#   |∇x|²    = 1.0              (constant)
#   drift terms = 0             (flat slab, no curvature)
#
# NOTE on precision: stella writes geometry files in Fortran e13.4 format
# (4-digit mantissa, e.g. 0.3142E+01 for π). This gives ~4 significant
# figures. Analytic checks use atol=5e-3 to accommodate this precision.
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
input_filename = 'slab_geometry.in'
stella_local_run_directory = 'Not/Run/Yet'
run_data = {}

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                    Check whether output files are present                    #
#-------------------------------------------------------------------------------
def test_whether_slab_output_files_are_present(tmp_path, stella_version, error=False):

    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global stella_local_run_directory, run_data
    stella_local_run_directory = tmp_path

    # Run stella inside of <tmp_path> based on <input_filename>
    run_data = run_local_stella_simulation(input_filename, tmp_path, stella_version)

    # Gather the output files generated during the local stella run inside <tmp_path>
    local_files = os.listdir(stella_local_run_directory)

    # Create a list of the output files we expect when stella has been run
    expected_files = ['slab_geometry.geometry']

    # Check whether all these output files are present
    for expected_file in expected_files:
        if not (expected_file in local_files):
            print(f'ERROR: The "{expected_file}" output file was not generated when running stella.'); error = True

    assert (not error), f'Some output files were not generated when running stella.'
    print(f'  -->  All the expected files are generated.')
    return

#-------------------------------------------------------------------------------
#                    Check whether slab output files match                     #
#-------------------------------------------------------------------------------
def test_whether_slab_geometry_output_matches(stella_version):
    '''Check that the results are identical to the analytic slab values.'''

    # File names
    local_geometry_file    = stella_local_run_directory / 'slab_geometry.geometry'
    expected_geometry_file = get_stella_expected_run_directory() / 'EXPECTED_OUTPUT.slab_geometry.geometry'

    # Compare geometry file against expected output
    compare_geometry_files(local_geometry_file, expected_geometry_file, error=False)

    print(f'  -->  Slab geometry output file matches expected values.')
    return

#-------------------------------------------------------------------------------
#           Verify slab metric values against analytic expressions             #
#-------------------------------------------------------------------------------
def test_slab_metric_analytic_values():
    '''Verify that metric coefficients match the analytic slab formulae.

    The geometry file is written in Fortran e13.4 format (4-digit mantissa),
    so values like zed=π are stored as 3.142 and gds2 as 3.467, giving a
    precision of ~0.001. The analytic checks use atol=5e-3 accordingly.
    '''

    import numpy as np

    # Slab parameters used in slab_geometry.in
    shat = 0.5
    qinp = 1.4
    rmaj = 3.0
    betaprim = 0.0

    # Read the geometry file written by stella
    local_geometry_file = stella_local_run_directory / 'slab_geometry.geometry'

    with open(local_geometry_file, 'r') as f:
        lines = f.readlines()

    # Read scalar parameters from line 2
    scalars = [v for v in lines[1].replace('#','').split() if v]
    rhoc_out      = float(scalars[0])
    qinp_out      = float(scalars[1])
    shat_out      = float(scalars[2])
    aref_out      = float(scalars[4])
    bref_out      = float(scalars[5])
    dxdpsi_out    = float(scalars[6])
    dydalpha_out  = float(scalars[7])
    exb_nonlin_out = float(scalars[8])

    # Check scalar parameters (these are exact in e13.4 format)
    assert np.isclose(rhoc_out, 1.0),       f"rhoc should be 1.0, got {rhoc_out}"
    assert np.isclose(qinp_out, qinp),      f"qinp should be {qinp}, got {qinp_out}"
    assert np.isclose(shat_out, shat),      f"shat should be {shat}, got {shat_out}"
    assert np.isclose(aref_out, 1.0),       f"aref should be 1.0, got {aref_out}"
    assert np.isclose(bref_out, 1.0),       f"bref should be 1.0, got {bref_out}"
    assert np.isclose(dxdpsi_out, 1.0),     f"dxdpsi should be 1.0, got {dxdpsi_out}"
    assert np.isclose(dydalpha_out, 1.0),   f"dydalpha should be 1.0, got {dydalpha_out}"
    assert np.isclose(exb_nonlin_out, 0.5), f"exb_nonlin should be 0.5, got {exb_nonlin_out}"

    # Read per-z data; stella writes 12 columns in the new geometry format
    data = np.loadtxt(local_geometry_file, skiprows=4, dtype='float')
    if data.ndim == 1:
        data = data.reshape(1, -1)

    zed_col  = data[:, 1]
    bmag_col = data[:, 3]
    bdot_col = data[:, 4]

    ncols = data.shape[1]
    if ncols == 12:
        # New format columns: alpha, zed, zeta, bmag, b.Gz, Gy.Gy, Gx.Gy, Gx.Gx,
        #                     BxGB.Gy, Bxkappa.Gy, BxGB.Gx, bmag_psi0
        grady_dot_grady         = data[:, 5]
        gradx_dot_grady         = data[:, 6]
        gradx_dot_gradx         = data[:, 7]
        B_times_gradB_dot_grady = data[:, 8]
        B_times_kappa_dot_grady = data[:, 9]
        B_times_gradB_dot_gradx = data[:, 10]
        gds2    = grady_dot_grady
        gds21   = gradx_dot_grady * shat_out
        gbdrift  = B_times_gradB_dot_grady * 2.0
        cvdrift  = B_times_kappa_dot_grady * 2.0
        gbdrift0 = B_times_gradB_dot_gradx * 2.0 * shat_out
    else:
        # Old 15-column format: ..., gds2, gds21, gds22, gds23, gds24, gbdrift, cvdrift, gbdrift0, ...
        gds2     = data[:, 5]
        gds21    = data[:, 6]
        gbdrift  = data[:, 10]
        cvdrift  = data[:, 11]
        gbdrift0 = data[:, 12]

    # --- Analytic checks ---
    # The z-grid has nzed=4 → 5 points at -π, -π/2, 0, π/2, π.
    # We use the EXACT analytic values for comparison (not the rounded file values),
    # with atol=5e-3 to account for the e13.4 output precision (~0.001 uncertainty).

    import math
    nzgrid = 2  # nzed/2
    zed_exact = np.array([-math.pi, -math.pi/2, 0.0, math.pi/2, math.pi])

    # 1. bmag = 1 everywhere in slab
    assert np.allclose(bmag_col, 1.0, atol=1e-6), \
        f"bmag should be 1.0 everywhere, got {bmag_col}"

    # 2. b.grad z = 1/(q*R), exact in the output to 4 sig figs
    expected_bdot = 1.0 / (qinp * rmaj)
    assert np.allclose(bdot_col, expected_bdot, rtol=1e-3), \
        f"b_dot_gradz should be {expected_bdot:.6f}, got {bdot_col}"

    # 3. |∇y|² = 1 + (shat*z)²  (atol=5e-3 for e13.4 precision)
    expected_gds2 = 1.0 + (shat * zed_exact)**2
    assert np.allclose(gds2, expected_gds2, atol=5e-3), \
        f"gds2 (|∇y|²) mismatch (atol=5e-3):\n  expected={expected_gds2}\n  got={gds2}"

    # 4. ∇x·∇y * shat = -shat²*z  (atol=5e-3)
    expected_gds21 = (-shat * zed_exact) * shat
    assert np.allclose(gds21, expected_gds21, atol=5e-3), \
        f"gds21 mismatch (atol=5e-3):\n  expected={expected_gds21}\n  got={gds21}"

    # 5. Drift terms all zero (flat slab, betaprim=0)
    assert np.allclose(gbdrift, 0.0, atol=1e-6), \
        f"gbdrift should be 0 for flat slab, got {gbdrift}"
    assert np.allclose(cvdrift, 0.0, atol=1e-6), \
        f"cvdrift should be 0 for flat slab with betaprim=0, got {cvdrift}"
    assert np.allclose(gbdrift0, 0.0, atol=1e-6), \
        f"gbdrift0 should be 0 for flat slab, got {gbdrift0}"

    print(f'  -->  All analytic metric checks passed.')
    return

#-------------------------------------------------------------------------------
#              Check whether the data in the netcdf file matches               #
#-------------------------------------------------------------------------------
def test_whether_geometry_data_in_netcdf_file_is_correct(error=False):
    '''Compare geometry arrays in the netcdf file against a reference.

    This test is skipped if no reference netcdf file exists yet (the reference
    must be generated by running stella and committing the output file as
    EXPECTED_OUTPUT.slab_geometry.out.nc).
    '''
    expected_nc = get_stella_expected_run_directory() / 'EXPECTED_OUTPUT.slab_geometry.out.nc'
    if not expected_nc.exists():
        pytest.skip('No EXPECTED_OUTPUT.slab_geometry.out.nc found; '
                    'run stella and commit the output to enable this test.')
    compare_geometry_in_netcdf_files(run_data, error=False)
    print('  -->  All geometry data in the netcdf file matches the expected output.')
    return
