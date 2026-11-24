################################################################################
#   Check that implicit collisions act when no other term is treated implicitly #
################################################################################
# The implicit part of the time step is skipped when fully_explicit = .true., which
# read_parameters_numerical sets when the mirror term, the parallel streaming and the
# drifts are all explicit (or switched off). Before November 2026, only the Fokker-Planck
# operator and the hyper-dissipation cleared this flag, so implicit Dougherty collisions
# were silently skipped in such runs: the collision timer stayed at 0 s and the parallel
# flow was never damped.
# Here we run a zonal mode with only Dougherty collisions (without momentum conservation)
# and the y drift, which does nothing for the zonal mode, starting from a parallel flow,
# with implicit and with explicit collisions, and check that:
#     - the collisions damp the parallel flow, sum_z |upar|, by about 18% at t = 2
#       (vnew_ref = 0.1), whereas it stayed constant when the collisions were skipped;
#     - implicit and explicit collisions agree (to 8e-5; the tolerance is 1e-3).
# The relevant parameters are:
#       &gyrokinetic_terms
#         include_parallel_streaming = .false.
#         include_mirror = .false.
#       /
#       &dissipation_and_collisions_options
#         include_collisions = .true.
#         collision_model = 'dougherty'
#         collisions_implicit = .true. or .false.
#       /
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

# Input file with implicit Dougherty collisions as the only implicit term
input_file = 'collisions_dougherty_only.in'

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
def get_parallel_flow(netcdf_file):
    """ Return sum_z |upar| at every time, relative to its value at t = 0. """
    with xr.open_dataset(netcdf_file) as netcdf:
        upar = netcdf['upar'].values
        upar = upar[..., 0] + 1j * upar[..., 1]
        flow = np.abs(upar.reshape(upar.shape[0], -1)).sum(axis=1)
    return flow / flow[0]

#-------------------------------------------------------------------------------
def test_whether_implicit_dougherty_collisions_act_without_other_implicit_terms(tmp_path, stella_version):

    # This input file only exists for the current stella version
    if stella_version != 'master': pytest.skip('Only implemented for the master branch of stella.')

    # Run the simulation with implicit and with explicit collisions
    flows = {}
    for implicit in [True, False]:
        folder = tmp_path / ('implicit' if implicit else 'explicit'); folder.mkdir()
        path_input_file = copy_input_file(input_file, folder)
        if not implicit:
            text = open(path_input_file, 'r').read()
            assert 'collisions_implicit = .true.' in text, f'Expected "collisions_implicit = .true." in {input_file}.'
            open(path_input_file, 'w').write(text.replace('collisions_implicit = .true.', 'collisions_implicit = .false.'))
        os.chdir(folder)
        run_stella(get_stella_path(stella_version), input_file, nproc=2)
        flows[implicit] = get_parallel_flow(folder / input_file.replace('.in', '.out.nc'))

    # The collisions must damp the parallel flow
    for implicit in [True, False]:
        label = 'implicit' if implicit else 'explicit'
        print(f'    {label} collisions: sum_z |upar| relative to t = 0: {flows[implicit][-1]:.6f}')
        assert np.all(np.isfinite(flows[implicit])), f'The parallel flow is not finite with {label} collisions.'
        assert flows[implicit][-1] < 0.9, f'The {label} Dougherty collisions do not damp the parallel flow.'

    # Implicit and explicit collisions must agree
    difference = np.max(np.abs(flows[True] - flows[False]))
    assert difference < 1e-3, f'Implicit and explicit Dougherty collisions differ by {difference:.2e} in the parallel flow.'
    print(f'  -->  Implicit Dougherty collisions act when no other term is implicit.')
    return
