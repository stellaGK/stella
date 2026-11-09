################################################################################
#        Check that the Fokker-Planck collision operator conserves density     #
################################################################################
# With only collisions switched on, the collision operator must conserve the
# density of every species, and it must not drive the potential. This is tested
# for ions and electrons with the real mass ratio at a large collision frequency,
# where the energy-exchange terms of the field-particle part of electrons
# colliding with ions used to create electron density, because part of these
# terms lies inside the first cells of the electron velocity grid. The potential
# then grew by a factor 6.5 by t = 0.5 with collisions only (Nov 2026).
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
input_filename = 'collisions_fokker_planck_conservation.in'

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#      Check that collisions conserve density and do not drive the potential   #
#-------------------------------------------------------------------------------
def test_whether_fokker_planck_collisions_conserve_density(tmp_path, stella_version):

    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    local_netcdf_file = tmp_path / input_filename.replace('.in','.out.nc')

    with xr.open_dataset(local_netcdf_file) as local_netcdf:

        # Density of every species on the (t, species, tube, z, kx, ky) grid
        density = local_netcdf['density'].data
        density = density[..., 0] + 1j * density[..., 1]
        for ispec in range(density.shape[1]):
            change = np.max(np.abs(density[:, ispec] - density[0, ispec])) / np.max(np.abs(density[0, ispec]))
            assert change < 1e-6, f'Collisions change the density of species {ispec+1} by {change:.2e} (relative).'

        # Collisions alone must not make the potential grow
        phi2 = local_netcdf['phi2'].data
        growth = np.max(phi2) / phi2[0]
        assert growth < 1.01, f'Collisions alone make |phi|^2 grow by a factor {growth:.3f}.'

    print(f'  -->  The Fokker-Planck collisions conserve density and do not drive the potential.')
    return
