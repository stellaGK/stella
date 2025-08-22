# Input parameters

This module will set the default input input parameters for each name list,
and it will read the stella input file per namelist.

For each namelists two (or three) routines will exist:
   - set_default_parameters_<namelist>
   - read_namelist_<namelist>
   - check_inputs_<namelist>

First we will set the default input parameters, and then we will overwrite
any default options with those specified in the input file. Optionally
we can check if any input variables are clashing with each other.

### Namelists

Overview of stella namelists:

GEOMETRY
  geometry_option
  overwrite_geometry
  geometry_vmec (renamed from vmec_parameters)
  geometry_miller (renamed from millergeo_parameters)

PHYSICS
  gyrokinetic_terms
  scale_gyrokinetic_terms
  adiabatic_electron_response
  adiabatic_ion_response
  electromagnetic
  full_flux_surface
  extra_physics

GRIDS
  vpamu_grid
  z_grid (renamed from z_grid_parameters)
  z_boundary_condition
  species_knobs
  species_parameters_1
  species_parameters_2
  kxky_grid_option
  kxky_grid_range
  kxky_grid_box

DIAGNOSTICS
  diagnostics
  diagnostics_potential
  diagnostics_omega
  diagnostics_distribution
  diagnostics_fluxes
  diagnostics_moments

INITIALIZE FIELDS
  initialize_distribution (renamed from init_g_knobs)
  initialize_distribution_maxwellian
  initialize_distribution_noise
  initialize_distribution_kpar
  initialize_distribution_rh
  restart_options

DISSIPATION AND COLLISIONS
  dissipation_and_collisions_options (renamed from &dissipation)
  collisions_dougherty
  collisions_fokker_planck (renamed from &collisions_fp)
  hyper_dissipation

TIME TRACE
  time_trace_options
  time_step

NUMERICS
  numerical_algorithms
  numerical_upwinding_for_derivatives

NEOCLASSICS
  neoclassical_input
  euterpe_parameters
  sources

RADIAL VARIATION
  radial_variation

PARALLELISATION
  parallelisation (renamed from &layout_knobs)

VERBOSE
  debug_flags
