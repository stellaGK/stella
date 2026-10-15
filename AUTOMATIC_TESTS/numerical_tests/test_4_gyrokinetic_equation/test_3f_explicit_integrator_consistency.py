
################################################################################
#     Cross-scheme consistency of the new explicit time integrators            #
################################################################################
# This is a CORRECTNESS (not regression) test for the six new explicit schemes:
#     SSP(3,2), SSP(4,2), SSP(5,2), SSP(4,3), SSP(9,3) and RK6-4.
#
# The same minimal linear problem (only the explicit parallel-streaming term
# active, single-ky-mode Gaussian init) is run with EACH of the six schemes at a
# fixed, small time step. The six resulting potential traces <phi2>(t) are
# required to agree within a loose relative tolerance.
#
# Rationale: every scheme integrates the *same* ODE, so at a small time step
# their trajectories must be close. This directly validates the ported
# Runge-Kutta coefficients -- a mistyped coefficient (the main risk in a port)
# would shift a scheme's trajectory and fail the comparison. The tolerance is
# deliberately loose (a consistency bound, not an absolute-accuracy check); the
# spread is dominated by the 2nd-order SSP(3,2) scheme deviating from the
# 3rd/6th-order cluster, measured at ~1.8e-3 relative at delt = 5e-4.
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

# The six new explicit schemes under test
SCHEMES = ['ssp32', 'ssp42', 'ssp52', 'ssp43', 'ssp93', 'rk64']

# Loose relative tolerance for the cross-scheme agreement (see module docstring).
RTOL = 1e-2

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#   Run the consistency problem with every scheme and collect <phi2>(t)       #
#-------------------------------------------------------------------------------
def test_whether_new_explicit_integrators_agree_with_each_other(tmp_path, stella_version, error=False):

    # For each scheme, make a copy of the base input with <explicit_algorithm>
    # set to that scheme, and run stella.
    base_input = 'integrator_consistency.in'
    base_text = (get_stella_expected_run_directory() / base_input).read_text()

    phi2_by_scheme = {}
    for scheme in SCHEMES:
        # Build the per-scheme input file in <tmp_path> by rewriting the
        # <explicit_algorithm> line. (For the base scheme this is a no-op,
        # which is fine: the base input already uses that scheme.)
        scheme_input = f'consistency_{scheme}.in'
        scheme_text = base_text.replace("explicit_algorithm = 'ssp93'",
                                        f"explicit_algorithm = '{scheme}'")
        (tmp_path / scheme_input).write_text(scheme_text)

        # Run stella on the per-scheme input (which already lives in <tmp_path>).
        # We call run_stella directly rather than run_local_stella_simulation,
        # because the latter re-copies the input from the source test directory,
        # where our generated <consistency_<scheme>.in> files do not exist.
        os.chdir(tmp_path)
        run_stella(get_stella_path(stella_version), scheme_input)

        # Read the potential trace <phi2>(t) from the output netCDF file
        local_netcdf_file = tmp_path / scheme_input.replace('.in', '.out.nc')
        with xr.open_dataset(local_netcdf_file) as local_netcdf:
            phi2 = np.asarray(local_netcdf['phi2'].values, dtype=float)
        phi2_by_scheme[scheme] = phi2

    # All schemes must have the same number of time points
    npts = len(phi2_by_scheme[SCHEMES[0]])
    assert all(len(phi2_by_scheme[s]) == npts for s in SCHEMES), \
        f'The number of time points differs between schemes: {[len(phi2_by_scheme[s]) for s in SCHEMES]}'

    # The initial value (t = 0) is the same initialization for every scheme
    t0 = phi2_by_scheme[SCHEMES[0]][0]
    for scheme in SCHEMES[1:]:
        if not (np.isclose(phi2_by_scheme[scheme][0], t0, rtol=1e-10, atol=0.0)):
            print(f'ERROR: <phi2> at t=0 differs between {SCHEMES[0]} and {scheme}.')
            error = True

    # Compare every pair of schemes over the whole trace (excluding t=0, which
    # is the identical initialization).
    ref = np.max(np.abs(phi2_by_scheme[SCHEMES[0]]))
    for i in range(len(SCHEMES)):
        for j in range(i + 1, len(SCHEMES)):
            a, b = SCHEMES[i], SCHEMES[j]
            pa, pb = phi2_by_scheme[a][1:], phi2_by_scheme[b][1:]
            rel = np.max(np.abs(pa - pb)) / ref
            if rel > RTOL:
                print(f'ERROR: {a} and {b} disagree beyond tolerance: max rel diff = {rel:.3e} > {RTOL:.1e}')
                error = True

    if error:
        # Print the traces to aid debugging
        print('\n<phi2>(t) traces for each scheme:')
        for scheme in SCHEMES:
            print(f'  {scheme:6s}: {np.array2string(phi2_by_scheme[scheme], precision=5, suppress_small=True)}')
        assert (not error), f'New explicit integrators do not agree within {RTOL:.1e} relative.'

    print(f'  -->  The six new explicit integrators agree with each other within {RTOL:.1e} relative '
          f'({npts} time points).')
    return
