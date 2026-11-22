################################################################################
#        Check that implicit and explicit collisions agree with A_parallel     #
################################################################################
# With A_par, the time derivative of the gyrokinetic equation acts on
#     gbar = g + 2 (Z_s/T_s) J_0 vpa A_par F_s,
# so a collisional change of the parallel current changes A_par, and the inductive
# electric field -dA_par/dt acts back on the distribution. The explicit time advance
# includes this (the collisions are added to the right-hand side for gbar). Before
# November 2026, the implicit collision step advanced g with A_par fixed and A_par was
# then recomputed from g, so this inductive response was missing: the electron current
# was damped at nu_ei instead of resistively, which suppressed the microtearing mode.
# Here we run the same linear electromagnetic simulation with implicit and with explicit
# collisions, and check that phi and A_par agree at the end of the run:
#     - Dougherty collisions without momentum conservation (the electron current is
#       damped), nvgrid = 8, nmu = 4: the two agree to 0.7% (first order in delt),
#       while without the inductive response they differed by 20% in phi.
#     - Fokker-Planck test-particle collisions between like species only (the explicit
#       Fokker-Planck operator has no inter-species collisions and no field-particle
#       part), nvgrid = 16, nmu = 8: the two discretisations agree to 3.6% (0.7% with
#       nvgrid = 32, nmu = 16), while without the inductive response they differed by 21% in phi.
# The relevant parameters are:
#       &electromagnetic
#         include_apar = .true.
#         include_bpar = .true.
#       /
#       &dissipation_and_collisions_options
#         include_collisions = .true.
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

# Input file with implicit Dougherty collisions
input_file = 'EM_implicit_collisions.in'

# Collision operator, replacements in the input file, and the tolerance on max|implicit - explicit|/max
collision_operators = {
    'dougherty': ([], 0.02),
    'fokker-planck': ([("collision_model = 'dougherty'", "collision_model = 'fokker-planck'"),
                       ("&collisions_dougherty\n  momentum_conservation = .false.\n/",
                        "&collisions_fokker_planck\n  interspec = .false.\n  fieldpart = .false.\n/"),
                       ('nmu = 4', 'nmu = 8'), ('nvgrid = 8', 'nvgrid = 16')], 0.1),
}

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
def write_input_file(folder, replacements, implicit):
    """ Copy <input_file> to <folder> with the collision operator given by <replacements>
    and with implicit or explicit collisions. """
    path_input_file = copy_input_file(input_file, folder)
    text = open(path_input_file, 'r').read()
    if not implicit: replacements = replacements + [('collisions_implicit = .true.', 'collisions_implicit = .false.')]
    for old, new in replacements:
        assert old in text, f'Expected "{old}" in {input_file}.'
        text = text.replace(old, new)
    open(path_input_file, 'w').write(text)
    return path_input_file

#-------------------------------------------------------------------------------
@pytest.mark.parametrize('collision_model', list(collision_operators))
def test_whether_implicit_collisions_with_apar_agree_with_explicit_collisions(collision_model, tmp_path, stella_version):

    # These input files only exist for the current stella version
    if stella_version != 'master': pytest.skip('Only implemented for the master branch of stella.')
    replacements, rtol = collision_operators[collision_model]

    # Run the simulation with implicit and with explicit collisions
    netcdf_files = {}
    for implicit in [True, False]:
        folder = tmp_path / ('implicit' if implicit else 'explicit'); folder.mkdir()
        write_input_file(folder, replacements, implicit)
        os.chdir(folder)
        run_stella(get_stella_path(stella_version), input_file, nproc=4)
        netcdf_files[implicit] = folder / input_file.replace('.in', '.out.nc')

    # Compare phi and A_par at the last time
    with xr.open_dataset(netcdf_files[True]) as implicit, xr.open_dataset(netcdf_files[False]) as explicit:
        assert np.allclose(implicit['t'], explicit['t']), 'The time axes do not match.'
        failed = []
        for key in ['phi_vs_t', 'apar_vs_t']:
            field_implicit = implicit[key].values[-1]; field_explicit = explicit[key].values[-1]
            scale = np.max(np.abs(field_explicit))
            difference = np.max(np.abs(field_implicit - field_explicit))
            assert np.isfinite(scale) and np.isfinite(difference), f'{key} is not finite: the simulation is unstable.'
            print(f'    {key:<10s} max|implicit - explicit|/max = {difference/scale:.3e}')
            if not difference <= rtol * scale: failed.append(key)
    assert not failed, f'Implicit and explicit {collision_model} collisions with A_par differ in {failed}.'
    print(f'  -->  Implicit and explicit {collision_model} collisions with A_par agree.')
    return
