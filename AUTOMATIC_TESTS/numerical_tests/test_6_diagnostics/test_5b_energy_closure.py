################################################################################
#             Free-energy closure regression test                               #
################################################################################
# Exercises the free-energy closure under several electron models and
# configurations, and checks the closure arithmetic and the NetCDF outputs.
#
# The closure writes an 18-column '.energy' ascii file:
#   time, particle_energy, field_energy, total_energy, dedt_particle,
#   dedt_field, dedt_total, dedt_projected, drive, heat_sources, collisions,
#   diss_vpa, diss_zed, diss_perp, drifts, streaming, mirror, nonlinearity
# where dedt_projected = <c_f|RHS at t_n> (projection identity) and dedt_total =
# (H(t_n)-H(t_{n-1}))/(t_n-t_{n-1}) is the exact two-level change, reported
# SEPARATELY and not equated to the projection.
#
# Cases (each run once, results cached):
#   kinetic         : linear, kinetic electrons (2 species)
#   adiabatic_nombt : linear, adiabatic electrons (periodic / no MBR)
#   adiabatic_mbt   : linear, adiabatic electrons (field-line average / MBR)
#   nonlinear       : nonlinear, adiabatic electrons (small box)
#   per_species     : kinetic + write_energy_species on
#   default_off     : write_energy off (diagnostic must be fully gated)
#   exact_zero      : drive/drifts/streaming/mirror/hyper off (columns = 0)
#   ascii_only      : write_energy on, kxkyz/vmu NetCDF off
#
# Checks per physics case: ascii produced (18 cols), particle/field energy split,
# dedt_total finite-difference consistency, NetCDF closure vars, per-mode identity,
# and kxkyz / vmu cross-representation consistency (both match the ASCII column).
# Extra checks: per-species files, default-off gating, exact-zero guards,
# ASCII-only path.
################################################################################

# Python modules
import pytest
import os, sys
import pathlib
import subprocess
import numpy as np
import xarray as xr
import tempfile

# Package to run stella
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# (case name, input file, nstep, nspec, has_kxkyz_vmu_nc)
CASES = [
    ('kinetic',         'energy_closure.in',          10, 2, True),
    ('adiabatic_nombt', 'energy_adiabatic_nombt.in',  10, 1, True),
    ('adiabatic_mbt',   'energy_adiabatic_mbt.in',    10, 1, True),
    ('nonlinear',       'energy_nonlinear.in',         8, 1, True),
]

# Special-configuration cases (single-file inputs)
SPECIAL_FILES = {
    'per_species': 'energy_per_species.in',
    'default_off': 'energy_default_off.in',
    'exact_zero':  'energy_exact_zero.in',
    'ascii_only':  'energy_ascii_only.in',
}

# Cache of completed runs, shared across all tests
_run_cache = {}

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                           Run a case once (cached)
#-------------------------------------------------------------------------------
def run_case(input_file, key, stella_version):
    if key in _run_cache:
        return _run_cache[key]
    stem = input_file.replace('.in', '')
    tmp_dir = pathlib.Path(tempfile.mkdtemp(prefix=f'energy_{key}_'))
    # Run the (short) energy cases on a single rank for determinism; the small
    # velocity grids can produce NaN under large-rank decompositions.
    run_local_stella_simulation(input_file, tmp_dir, stella_version, nproc=1)
    res = {
        'dir': tmp_dir,
        'energy': tmp_dir / f'{stem}.energy',
        'nc': tmp_dir / f'{stem}.out.nc',
    }
    _run_cache[key] = res
    return res

#-------------------------------------------------------------------------------
#   Base checks applied to each physics case (parametrized)
#-------------------------------------------------------------------------------
@pytest.fixture(scope="module", params=CASES, ids=[c[0] for c in CASES])
def result(request, stella_version):
    name, f, nstep, nspec, has_nc = request.param
    res = run_case(f, name, stella_version)
    res['nstep'] = nstep
    res['nspec'] = nspec
    res['has_nc'] = has_nc
    return res

#-------------------------------------------------------------------------------
def test_energy_ascii_file_produced(result):
    """'.energy' exists, is finite, and has 18 columns (and nstep+1 rows)."""
    energy = result['energy']
    assert energy.exists(), f'.energy file {energy} not produced.'
    data = np.loadtxt(energy, dtype='float', comments='#')
    assert data.ndim == 2, f'expected 2D data, got shape {data.shape}'
    nrows, ncols = data.shape
    assert ncols == 18, f'expected 18 columns, got {ncols}'
    assert nrows >= 2, f'expected >= 2 rows, got {nrows}'
    assert np.all(np.isfinite(data)), 'non-finite values in .energy file'

#-------------------------------------------------------------------------------
def test_energy_split(result):
    """particle_energy (col 2) and field_energy (col 3) are finite and >= 0, and
    total_energy (col 4) = particle_energy + field_energy."""
    data = np.loadtxt(result['energy'], dtype='float', comments='#')
    particle = data[:, 1]
    field = data[:, 2]
    total = data[:, 3]
    assert np.all(np.isfinite(particle)) and np.all(np.isfinite(field)) and np.all(np.isfinite(total))
    assert np.all(particle >= -1e-12), 'particle energy should be >= 0'
    assert np.all(field >= -1e-12), 'field energy should be >= 0'
    assert np.allclose(total, particle + field, atol=1e-12), 'total_energy != particle + field'

#-------------------------------------------------------------------------------
def test_dedt_total_is_separate(result):
    """col 7 (dedt_total) is the finite difference of (particle+field) between
    saved steps: dedt_total[t] = (H(t)-H(t-1))/(t(t)-t(t-1)). This is the actual
    two-level change, reported separately from the projection (col 8)."""
    data = np.loadtxt(result['energy'], dtype='float', comments='#')
    t = data[:, 0]
    H = data[:, 1] + data[:, 2]
    dedt_total = data[:, 6]
    assert np.all(np.isfinite(dedt_total))
    # skip the first (degenerate) startup row
    with np.errstate(divide='ignore', invalid='ignore'):
        fd = np.diff(H) / np.diff(t)
    assert np.allclose(fd, dedt_total[1:], rtol=1e-2, atol=1e-9), \
        f'dedt_total not a finite difference of H: max |fd-dd| = {np.max(np.abs(fd - dedt_total[1:]))}'

#-------------------------------------------------------------------------------
def test_netcdf_energy_variables(result):
    """The closure NetCDF variables are present, shaped and finite."""
    if not result['has_nc']:
        pytest.skip('NetCDF energy output not enabled for this case')
    ds = xr.open_dataset(result['nc'])
    kxkyz_vars = ['free_energy_kxkyz', 'field_energy_kxkyz', 'dedt_projected_kxkyz',
                  'drive_kxkyz', 'collisions_kxkyz', 'diss_vpa_kxkyz', 'diss_zed_kxkyz',
                  'diss_perp_kxkyz', 'nonlinearity_kxkyz', 'streaming_kxkyz',
                  'mirror_kxkyz', 'drifts_kxkyz']
    for v in kxkyz_vars:
        assert v in ds.variables, f'NetCDF var {v} missing.'
        assert np.all(np.isfinite(ds[v].values)), f'NetCDF var {v} not finite.'

    vmu_vars = ['free_energy_vmu', 'field_energy_vmu', 'dedt_projected_vmu',
                'drive_vmu', 'collisions_vmu', 'diss_vpa_vmu', 'diss_zed_vmu',
                'diss_perp_vmu', 'nonlinearity_vmu', 'streaming_vmu',
                'mirror_vmu', 'drifts_vmu']
    for v in vmu_vars:
        assert v in ds.variables, f'NetCDF var {v} missing.'
        assert np.all(np.isfinite(ds[v].values)), f'NetCDF var {v} not finite.'

    # Per-mode projection identity in the NetCDF (last saved step).
    dd = ds['dedt_projected_kxkyz'].values[-1]
    sm = (ds['drive_kxkyz'].values[-1] + ds['collisions_kxkyz'].values[-1] +
          ds['diss_vpa_kxkyz'].values[-1] + ds['diss_zed_kxkyz'].values[-1] +
          ds['diss_perp_kxkyz'].values[-1] + ds['nonlinearity_kxkyz'].values[-1] +
          ds['streaming_kxkyz'].values[-1] + ds['mirror_kxkyz'].values[-1] +
          ds['drifts_kxkyz'].values[-1])
    assert np.allclose(dd, sm, atol=1e-12), \
        f'per-mode projection identity fails: max diff = {np.max(np.abs(dd - sm))}'
    ds.close()

#-------------------------------------------------------------------------------
#                      Cross-representation consistency
#-------------------------------------------------------------------------------
def test_kxkyz_consistency(result):
    """Sum of per-mode free_energy_kxkyz weighted by dVolume/volume_total equals
    the ASCII particle_energy (tight tolerance). In a flux tube there is one
    alpha (tube) but many ky harmonics, so all ky modes share the same volume
    element dVolume(alpha=0, kx, z)."""
    if not result['has_nc']:
        pytest.skip('NetCDF energy output not enabled for this case')
    ds = xr.open_dataset(result['nc'])
    # dVolume dims: (zed, kx, alpha); free_energy_kxkyz dims: (t,spec,tube,zed,kx,ky)
    fe = ds['free_energy_kxkyz'].values
    dv = ds['dVolume'].values
    nt, nspec, ntube, nz, nkx, nky = fe.shape
    vol = dv[:, 0, 0].sum()
    integ = np.zeros(nt)
    for t in range(nt):
        for s in range(nspec):
            for tu in range(ntube):
                for z in range(nz):
                    for kx in range(nkx):
                        for ky in range(nky):
                            integ[t] += fe[t, s, tu, z, kx, ky] * dv[z, kx, 0]
    integ /= vol
    data = np.loadtxt(result['energy'], dtype='float', comments='#')
    particle = data[:, 1]
    assert np.allclose(integ[1:], particle[1:], rtol=1e-4, atol=1e-12), \
        f'kxkyz-consistency fails: max rel = {np.max(np.abs(integ[1:] - particle[1:]) / np.abs(particle[1:]))}'
    ds.close()

#-------------------------------------------------------------------------------
def test_vmu_consistency(result):
    """Sum of per-mode free_energy_vmu weighted by wgts_mu_bare * wgts_vpa equals
    the ASCII particle_energy. The vpa weight used by the free-energy quadrature
    is the CONSERVATIVE one (dvpa/sqrt(pi)) set by set_vpa_weights(.true.); the
    out.nc wgts_vpa is the default (conservative flag) so we reconstruct dvpa here
    from the vpa grid. Matches to ~1e-8, confirming kxkyz and vmu give the same H_p."""
    if not result['has_nc']:
        pytest.skip('NetCDF energy output not enabled for this case')
    ds = xr.open_dataset(result['nc'])
    fv = ds['free_energy_vmu'].values   # (nt, nspec, nmu, nvpa)
    wm = ds['wgts_mu_bare'].values      # (nmu,)
    vpa = ds['vpa'].values
    wv = np.gradient(vpa) / np.sqrt(np.pi)   # conservative wgts_vpa = dvpa/sqrt(pi)
    nt = fv.shape[0]
    integ = np.zeros(nt)
    for t in range(nt):
        for s in range(fv.shape[1]):
            integ[t] += np.sum(fv[t, s, :, :] * wm[:, None] * wv[None, :])
    data = np.loadtxt(result['energy'], dtype='float', comments='#')
    particle = data[:, 1]
    assert np.allclose(integ[1:], particle[1:], rtol=1e-4, atol=1e-12), \
        f'vmu-consistency fails: max rel = {np.max(np.abs(integ[1:] - particle[1:]) / np.abs(particle[1:]))}'
    ds.close()

#-------------------------------------------------------------------------------
#                      Special-configuration checks
#-------------------------------------------------------------------------------
def test_per_species_output(stella_version):
    """Per-species files exist; the per-species particle energy sums to the ASCII
    combined particle energy; the field energy is NOT in the per-species files."""
    res = run_case(SPECIAL_FILES['per_species'], 'per_species', stella_version)
    # For the kinetic base (nspec=2) the labels are 'ion_1' and 'electron_1'.
    ion = pathlib.Path(str(res['dir'])) / 'energy_per_species.energy_species_ion_1'
    elec = pathlib.Path(str(res['dir'])) / 'energy_per_species.energy_species_electron_1'
    assert ion.exists() and elec.exists(), 'per-species files not produced.'
    di = np.loadtxt(ion, dtype='float', comments='#')
    de = np.loadtxt(elec, dtype='float', comments='#')
    # 12 columns: time, energy, dedt_proj, drive, coll, diss_vpa, diss_zed,
    #             diss_perp, nonlinearity, streaming, mirror, drifts
    assert di.shape[1] == 12 and de.shape[1] == 12
    comb = np.loadtxt(res['energy'], dtype='float', comments='#')
    # per-species particle energy (col 1) summed = combined particle_energy (col 1)
    assert np.allclose(di[:, 1] + de[:, 1], comb[:, 1], rtol=1e-3, atol=1e-12), \
        'per-species particle energy does not sum to combined'

#-------------------------------------------------------------------------------
def test_default_off_gating(stella_version):
    """write_energy=.false. -> no '.energy' file and no free-energy NetCDF vars."""
    res = run_case(SPECIAL_FILES['default_off'], 'default_off', stella_version)
    assert not res['energy'].exists(), '.energy file should not be produced when write_energy is off.'
    ds = xr.open_dataset(res['nc'])
    energy_vars = [v for v in ds.variables if 'energy' in v or v.startswith(('free_', 'dedt_proj', 'drive_', 'diss_', 'streaming_', 'mirror_', 'drifts_', 'nonlinearity_', 'collisions_', 'field_'))]
    assert not any(('kxkyz' in v or 'vmu' in v) for v in energy_vars), \
        f'free-energy NetCDF vars leaked when off: {energy_vars}'
    ds.close()

#-------------------------------------------------------------------------------
def test_exact_zero_guards(stella_version):
    """When drive/drifts/streaming/mirror/hyper are off, the corresponding term
    columns (cols 9..18) are exactly 0 and the projection identity still closes."""
    res = run_case(SPECIAL_FILES['exact_zero'], 'exact_zero', stella_version)
    data = np.loadtxt(res['energy'], dtype='float', comments='#')
    # cols 9..18 are the 10 term columns (drive...nonlinearity). The exact-change
    # columns (4..7) are not terms and need not be zero.
    assert np.all(data[:, 8:18] == 0.0), 'off term columns should be exactly zero.'
    # identity still closes
    assert np.allclose(data[:, 7], np.sum(data[:, 8:18], axis=1), atol=1e-12)

#-------------------------------------------------------------------------------
def test_ascii_only_path(stella_version):
    """write_energy on but kxkyz/vmu off -> '.energy' written, no energy NetCDF vars."""
    res = run_case(SPECIAL_FILES['ascii_only'], 'ascii_only', stella_version)
    assert res['energy'].exists(), '.energy file not produced.'
    data = np.loadtxt(res['energy'], dtype='float', comments='#')
    assert data.shape[1] == 18
    ds = xr.open_dataset(res['nc'])
    energy_vars = [v for v in ds.variables if 'kxkyz' in v or 'vmu' in v or 'free_energy' in v or 'field_energy' in v]
    assert not energy_vars, f'free-energy NetCDF vars should be absent, got {energy_vars}'
    ds.close()
