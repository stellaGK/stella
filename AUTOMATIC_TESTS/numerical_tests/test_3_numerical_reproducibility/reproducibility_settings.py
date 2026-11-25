################################################################################
#                Settings shared by the numerical reproducibility tests        #
################################################################################
# These tests are designed to guard the numerics of stella while its memory
# usage or speed is being optimised. Such changes should not alter the results,
# hence the tolerances below are tight: the results of a nonlinear simulation
# change by less than 1e-10 (relative to the maximum of each quantity) when the
# number of MPI processes or the data layouts are changed, while a real change
# of the numerics typically gives differences of 1e-7 or larger.
################################################################################

# Input files in this folder, each switching on a different numerical scheme.
input_files = [
    'es_nonlinear.in',
    'em_nonlinear.in',
    'adiabatic_electrons.in',
    'collisions_dougherty_implicit.in',
    'collisions_dougherty_explicit.in',
    'collisions_fokker_planck_implicit.in',
    'collisions_fokker_planck_explicit.in',
    'explicit_algorithm_rk4.in',
    'explicit_algorithm_euler.in',
    'flip_flop.in',
    'mirror_without_semi_lagrange.in',
    'stream_and_mirror_explicit.in',
]

# Quantities that fully describe the state of the simulation: the fields on the
# full (t, tube, z, kx, ky) grid, the moments, the fluxes, the frequency, and the
# distribution functions on the (z, vpa, mu) grid. Keys that do not exist in the
# expected output (e.g. apar_vs_t for electrostatic runs) are skipped.
geometry_keys = ['bmag', 'kperp2']
state_keys = ['phi_vs_t', 'apar_vs_t', 'bpar_vs_t', 'density', 'upar', 'temperature',
              'pflux_vs_kxkys', 'vflux_vs_kxkys', 'qflux_vs_kxkys', 'omega',
              'g2_vs_zvpamus', 'h2_vs_zvpamus']
regression_keys = ['t'] + geometry_keys + state_keys

# Tolerance on max|local - expected| <= rtol * max|expected| + atol, where <atol> only
# matters for quantities which vanish analytically (e.g. the particle flux with
# adiabatic electrons, which is of the order of 1e-28)
rtol = 1e-8
atol = 1e-24

# The expected output is created on a single processor, see create_expected_output.py
nproc_expected_output = 1

# Known bugs where the results depend on the number of MPI processes. The MPI
# invariance tests are marked as strict expected failures for these input files,
# so once the bug is fixed, the test will report an XPASS (and fail) as a reminder
# to remove the input file from this dictionary.
known_mpi_bugs = {}
