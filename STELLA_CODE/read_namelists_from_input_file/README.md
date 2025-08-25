# Read namelists from input file

This module will set the default input parameters for each name list,
and it will read the stella input file per namelist.

For each namelists two (or three) routines will exist:
   - set_default_parameters_<namelist>
   - read_namelist_<namelist>
   - check_inputs_<namelist>

First we will set the default input parameters, and then we will overwrite
any default options with those specified in the input file. Optionally
we can check if any input variables are clashing with each other.

### Namelists

An overview of all default input parameters of stella can be found in:
   STELLA_CODE/read_namelists_from_input_file/default_input_file.in

Overview of stella namelists:

GEOMETRY
  geometry_options
  geometry_vmec
  geometry_miller
  geometry_zpinch
  geometry_from_txt

PHYSICS
  gyrokinetic_terms
  scale_gyrokinetic_terms
  adiabatic_electron_response
  adiabatic_ion_response
  electromagnetic
  flow_shear
  physics_inputs

KINETIC SPECIES
  species_options
  species_parameters_1
  species_parameters_2
  euterpe_parameters
  
GRIDS
  kxky_grid_option
  kxky_grid_range
  kxky_grid_box
  z_grid
  z_boundary_condition
  velocity_grids

DIAGNOSTICS
  diagnostics
  diagnostics_potential
  diagnostics_omega
  diagnostics_distribution
  diagnostics_fluxes
  diagnostics_moments

INITIALIZE FIELDS
  initialize_distribution
  initialize_distribution_maxwellian
  initialize_distribution_noise
  initialize_distribution_kpar
  initialize_distribution_rh
  restart_options

DISSIPATION AND COLLISIONS
  dissipation_and_collisions_options
  collisions_dougherty
  collisions_fokker_planck
  hyper_dissipation

TIME TRACE
  time_trace_options
  time_step

NUMERICS
  numerical_algorithms
  numerical_upwinding_for_derivatives
  flux_annulus

NEOCLASSICS
  neoclassical_input

RADIAL VARIATION
  multibox_parameters
  sources

PARALLELISATION
  parallelisation

VERBOSE
  debug_flags
