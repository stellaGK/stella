################################################################################
#   Check that exact_conservation conserves momentum and energy exactly       #
################################################################################
# With exact_conservation, the Fokker-Planck collision operator must conserve the
# total parallel momentum sum_s m_s n_s upar_s and the total energy
# sum_s 3/2 (delta T_s + T_s delta n_s) of all species to machine precision, and
# the density of every species. stella chooses the algorithm from the mu grid:
# the block products on the equally spaced grid, the test-particle matrix on the
# Gauss-Laguerre grid. The developer option exact_conservation_algorithm can also
# use the test-particle matrix on the equally spaced grid, to compare the two.
################################################################################

# Python modules
import pytest
import pathlib
import numpy as np
import xarray as xr

# Package to run stella
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Global variables
input_filename = 'collisions_fokker_planck_exact_conservation.in'

# Replacements in the input file for each case, the number of processors (the test-particle-matrix algorithm with
# several species only runs on one processor), and the tolerance on the density: on the Gauss-Laguerre grid the
# density is only conserved to about 1e-8 (Nov 2026)
cases = {
    'equally_spaced_mu_grid': ([], None, 1e-10),
    'gauss_laguerre_mu_grid': ([('  equally_spaced_mu_grid = .true.\n', '')], 1, 1e-7),
    'equally_spaced_mu_grid_test_particle_matrix': ([('&dissipation_and_collisions_options\n',
        "&collisions_fokker_planck_developers\n  exact_conservation_algorithm = 'test_particle_matrix'\n/\n"
        "&dissipation_and_collisions_options\n")], 1, 1e-10),
}

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#             Check the conservation of density, momentum and energy           #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('case', cases.keys())
def test_whether_exact_conservation_conserves_momentum_and_energy(case, tmp_path, stella_version):

    # Write the input file of this case, and run stella inside of <tmp_path>
    replacements, nproc, tolerance_density = cases[case]
    text = open(get_stella_expected_run_directory() / input_filename).read()
    for old, new in replacements:
        assert old in text, f'"{old}" is not in {input_filename}.'
        text = text.replace(old, new)
    with open(tmp_path / input_filename, 'w') as file: file.write(text)
    os.chdir(tmp_path)
    run_stella(get_stella_path(stella_version), input_filename, nproc=nproc)
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc')

    with xr.open_dataset(local_netcdf_file) as data:

        # Moments on the (t, species, tube, z, kx, ky) grid
        complex_moment = lambda name: data[name].values[..., 0] + 1j * data[name].values[..., 1]
        dens, upar, temp = complex_moment('density'), complex_moment('upar'), complex_moment('temperature')
        shape = (1, -1, 1, 1, 1, 1)
        mass = data['mass'].values.reshape(shape)
        n0 = data['dens'].values.reshape(shape)
        T0 = data['temp'].values.reshape(shape)

        # Momentum and energy of every species, and their sum over the species
        momentum_s = mass * n0 * upar
        energy_s = 1.5 * (temp + T0 * dens)
        momentum, energy = momentum_s.sum(axis=1), energy_s.sum(axis=1)

        # The largest change in time, relative to the largest value of a single species
        change_momentum = np.max(np.abs(momentum - momentum[0])) / np.max(np.abs(momentum_s))
        change_energy = np.max(np.abs(energy - energy[0])) / np.max(np.abs(energy_s))
        change_density = max(np.max(np.abs(dens[:, s] - dens[0, s])) / np.max(np.abs(dens[:, s]))
            for s in range(dens.shape[1]))
        # The collisions must exchange momentum and energy between the species, otherwise the test is empty
        exchange = np.max(np.abs(momentum_s[:, 0] - momentum_s[0, 0])) / np.max(np.abs(momentum_s))
        print(f'  Changes: momentum {change_momentum:.1e}, energy {change_energy:.1e}, density {change_density:.1e}; '
              f'momentum exchanged between the species {exchange:.1e}')

        assert exchange > 1e-4, f'The collisions hardly exchange momentum ({exchange:.1e}): the test is empty.'
        assert change_momentum < 1e-11, f'exact_conservation changes the momentum by {change_momentum:.2e}.'
        assert change_energy < 1e-11, f'exact_conservation changes the energy by {change_energy:.2e}.'
        assert change_density < tolerance_density, f'exact_conservation changes the density by {change_density:.2e}.'

    print(f'  -->  exact_conservation conserves momentum, energy and density ({case}).')
    return
