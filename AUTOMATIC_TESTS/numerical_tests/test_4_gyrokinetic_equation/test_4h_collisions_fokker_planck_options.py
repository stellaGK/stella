################################################################################
#              Check the input options of the Fokker-Planck operator           #
################################################################################
# The options of the Fokker-Planck operator are read from <collisions_fokker_planck>
# (user options) and <collisions_fokker_planck_developers> (developer options).
# The names used before Nov 2026 (fieldpart, cfac, nuxfac, iiknob, no_j1l2, ...)
# are still accepted in <collisions_fokker_planck>; test that they give exactly
# the same simulation as the new names. Combinations that are not implemented,
# or that would silently do something else than asked, must stop stella.
################################################################################

# Python modules
import pytest
import pathlib
import subprocess
import numpy as np
import xarray as xr

# Package to run stella
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Input files: collisions only, with the field-particle part (test 4e), and with exact_conservation (test 4g)
input_conservation = 'collisions_fokker_planck_conservation.in'
input_exact = 'collisions_fokker_planck_exact_conservation.in'

# The <collisions_fokker_planck> namelist of <input_conservation>
namelist_conservation = '&collisions_fokker_planck\n  include_field_particle_part = .true.\n  include_flr_terms = .false.\n/'

# The same options with the new names and with the old names
options_new_names = ('&collisions_fokker_planck\n  include_field_particle_part = .true.\n  include_flr_terms = .false.\n/\n'
    '&collisions_fokker_planck_developers\n  scale_nu_ii = 0.5\n  include_j1_l2_term = .false.\n'
    '  scale_energy_diffusion_ei = 0.8\n  field_particle_jmax = 1\n/')
options_old_names = ('&collisions_fokker_planck\n  fieldpart = .true.\n  cfac = 0.0\n  iiknob = 0.5\n  no_j1l2 = .true.\n'
    '  eiediffknob = 0.8\n  jmax = 1\n  nvel_local = 64\n/')

# Inputs that must stop stella: the input file, replacements in it, the number of processors, and part of the message
developers = '&collisions_fokker_planck_developers\n'
before_dissipation = '&dissipation_and_collisions_options\n'
stops = {
    'explicit_algorithm_is_a_developer_option': (input_conservation,
        [('collisions_implicit = .true.', 'collisions_implicit = .false.'),
         (namelist_conservation, '&collisions_fokker_planck\n  include_field_particle_part = .false.\n/')],
        None, 'is a developer option'),
    'explicit_algorithm_without_unlike_species_collisions': (input_conservation,
        [(namelist_conservation, '&collisions_fokker_planck\n  include_field_particle_part = .false.\n/\n'
          + developers + '  use_explicit_algorithm = .true.\n/')],
        None, 'only has collisions between particles of the same species'),
    'explicit_algorithm_with_field_particle_part': (input_conservation,
        [(namelist_conservation, namelist_conservation + '\n' + developers
          + '  use_explicit_algorithm = .true.\n  include_unlike_species_collisions = .false.\n/')],
        None, 'is only implemented in the implicit algorithm'),
    'old_and_new_name': (input_conservation,
        [(namelist_conservation, namelist_conservation.replace('/', '  cfac = 0.0\n/'))],
        None, 'Both the deprecated option cfac'),
    'removed_advfield_coll': (input_conservation,
        [(namelist_conservation, namelist_conservation.replace('/', '  advfield_coll = .false.\n/'))],
        None, 'has been removed'),
    'exact_conservation_with_mixed_derivative_terms': (input_exact,
        [('  include_mixed_derivative_terms = .false.\n', '')], None, 'without the mixed-derivative terms'),
    'exact_conservation_without_conservative_weights': (input_exact,
        [('  conservative_wgts_vpa = .true.\n', '')], None, 'conservative vpa integration weights'),
    'exact_conservation_without_field_particle_part': (input_exact,
        [('include_field_particle_part = .true.', 'include_field_particle_part = .false.')],
        None, 'needs the field-particle part'),
    'exact_conservation_with_lmax_2': (input_exact,
        [(before_dissipation, developers + '  field_particle_lmax = 2\n/\n' + before_dissipation)],
        None, 'field_particle_lmax = field_particle_jmax = 1'),
    'exact_conservation_gauss_laguerre_on_several_processors': (input_exact,
        [('  equally_spaced_mu_grid = .true.\n', '')], 2, 'not MPI parallel'),
    'unknown_exact_conservation_algorithm': (input_exact,
        [(before_dissipation, developers + "  exact_conservation_algorithm = 'blocks'\n/\n" + before_dissipation)],
        None, 'exact_conservation_algorithm in &collisions_fokker_planck_developers must be'),
}

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
def write_input_file(input_filename, replacements, folder):
    '''Copy <input_filename> to <folder>, with the replacements (old, new) in the text.'''
    text = open(get_stella_expected_run_directory() / input_filename).read()
    for old, new in replacements:
        assert old in text, f'"{old}" is not in {input_filename}.'
        text = text.replace(old, new)
    with open(folder / input_filename, 'w') as file: file.write(text)

#-------------------------------------------------------------------------------
#               The old names give the same simulation as the new names        #
#-------------------------------------------------------------------------------
def test_whether_the_old_option_names_give_the_same_simulation(tmp_path, stella_version):
    if stella_version!='master': pytest.skip('Older stella versions do not have the new option names.')

    # Run the same simulation with the new and with the old names of the options
    netcdf_files = {}
    for names, options in [('new', options_new_names), ('old', options_old_names)]:
        folder = tmp_path / names; folder.mkdir()
        write_input_file(input_conservation, [(namelist_conservation, options)], folder)
        os.chdir(folder)
        run_stella(get_stella_path(stella_version), input_conservation)
        netcdf_files[names] = folder / input_conservation.replace('.in', '.out.nc')

    # The potential and the moments must be the same. They are identical on one thread; with several OpenMP threads,
    # two runs of the same input differ by about 1e-11 (the order of threaded sums, Nov 2026), hence the tolerance
    with xr.open_dataset(netcdf_files['new']) as new, xr.open_dataset(netcdf_files['old']) as old:
        for key in ['phi_vs_t', 'density', 'upar', 'temperature']:
            difference = np.max(np.abs(new[key].values - old[key].values)) / np.max(np.abs(new[key].values))
            assert difference < 1e-10, f'The old option names give another {key} (relative difference {difference:.1e}).'
    print('  -->  The old option names give the same simulation as the new names.')
    return

#-------------------------------------------------------------------------------
#                  Combinations that must stop stella at the start             #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('case', stops.keys())
def test_whether_stella_stops_on_wrong_fokker_planck_options(case, tmp_path, stella_version):
    if stella_version!='master': pytest.skip('Older stella versions do not check these options.')

    input_filename, replacements, nproc, message = stops[case]
    write_input_file(input_filename, replacements, tmp_path)
    os.chdir(tmp_path)
    nproc = nproc if nproc else read_nproc()
    result = subprocess.run(['mpirun', '--oversubscribe', '-np', f'{nproc}', get_stella_path(stella_version),
        input_filename], capture_output=True, text=True)
    output = result.stdout + result.stderr
    assert result.returncode != 0, f'stella did not stop ({case}).'
    assert message in output, f'stella did not stop with the message "{message}" ({case}):\n{output[-2000:]}'
    print(f'  -->  stella stops ({case}).')
    return
