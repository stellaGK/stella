# stella

[![Check stella](https://github.com/stellaGK/stella/actions/workflows/check_stella.yml/badge.svg)](https://github.com/stellaGK/stella/actions/workflows/check_stella.yml)

`stella` solves the gyrokinetic-Poisson system of equations in the local limit
using an operator-split, implicit-explicit numerical scheme. It is capable of
evolving electrostatic fluctuations with fully kinetic electrons and an
arbitrary number of ion species in general magnetic geometry, including
stellarators.

<br>

## Table of contents 
  * [Code](#code)
    + [stella.f90](#stellaf90)
    + [arrays](#arrays)
    + [Calculations](#calculations)
    + [Diagnostics](#diagnostics)
    + [Dissipation](#dissipation)
    + [Field_equations](#field_equations)
    + [Geometry](#geometry)
    + [Grids](#grids)
    + [Gyrokinetic_terms](#gyrokinetic-terms)
    + [Neoclassical](#neoclassical)
    + [Parameters](#parameters)
    + [Radial_variation](#radial-variation)
    + [Read_namelists_from_input_file](#read-namelists-from-input-file)

<br>

## Code

The code is organized in the following folders:
- arrays
- calculations 
- diagnostics  
- dissipation  
- field_equations
- geometry  
- grids  
- gyrokinetic_terms  
- neoclassical  
- parameters  
- radial_variation  
- stella.f90

<br>

### stella.f90

This is the main script.

<br>

### Arrays

The `arrays` folder contains the following scripts:

- arrays_distribution_function.f90
- arrays_fields.f90 
- arrays_gyro_averages.f90 
- arrays.f90 
- initialise_arrays.f90 
- initialise_distribution_function.f90 

<br>

### Calculations

The `calculations` folder contains the following scripts:

- calculations_add_explicit_terms.f90
- calculations_checksum.f90
- calculations_finite_differences.f90 
- calculations_gyro_averages.f90 
- calculations_kxky_derivatives.f90
- calculations_kxky.f90 
- calculations_redistribute.f90 
- calculations_timestep.f90 
- calculations_tofrom_ghf.f90
- calculations_transforms.f90 
- calculations_velocity_integrals.f90
- calculations_volume_averages.f90 

<br>

### Diagnostics

The `diagnostics` folder contains the following scripts:

- diagnostics.f90
- diagnostics_distribution.f90
- diagnostics_fluxes.f90
- diagnostics_fluxes_fluxtube.f90
- diagnostics_fluxes_fullfluxsurface.f90
- diagnostics_fluxes_radialvariation.f90
- diagnostics_moments.f90
- diagnostics_omega.f90
- diagnostics_potential.f90 
- stella_io.fpp
- stella_save.fpp

<br>


### Dissipation  

The `dissipation` folder contains the following scripts:

- collisions_dougherty.f90
- collisions_fokkerplanck.f90
- dissipation_and_collisions.f90
- dissipation_hyper.f90

<br>


### Field_equations

The `field_equations` folder contains the following scripts:

- field_equations_collisions.fpp
- field_equations_electromagnetic.fpp
- field_equations_fluxtube.fpp
- field_equations_fullfluxsurface.fpp
- field_equations_quasineutrality.fpp
- field_equations_radialvariation.fpp

<br>


### Geometry  

The `geometry` folder contains the following scripts:

- geometry.f90
- geometry_inputprofiles_interface.f90
- geometry_miller.f90
- geometry_vmec.f90
- geometry_vmec_read_netCDF_file.f90
- geometry_zpinch.f90

<br>



### Grids  

The `grids` folder contains the following scripts:

- grids_extended_zgrid.f90
- grids_kxky.f90
- grids_species_from_euterpe.f90
- grids_species.f90
- grids_time.f90
- grids_velocity.f90 
- grids_z.f90 

<br>


### Gyrokinetic terms  

The `gyrokinetic_terms` folder contains the following scripts:

- gk_drive.f90
- gk_ffs_solve.f90
- gk_flow_shear.f90
- gk_implicit_terms.f90
- gk_magnetic_drift.f90
- gk_mirror.f90
- gk_nonlinearity.f90
- gk_parallel_streaming.f90
- gyrokinetic_equation_explicit.f90
- gyrokinetic_equation_implicit.f90
- gyrokinetic_equation_initialisation.f90
- response_matrix.fpp

<br>


### Neoclassical  

The `neoclassical` folder contains the following scripts:

- neoclassical_terms.f90
- sfincs_interface.fpp

<br>


### Parameters  

The `parameters` folder contains the following scripts:

- common_type.f90
- debug_flags.f90
- file_units.f90 
- interface_random_number_generator.f90 
- parallelisation_layouts.f90 
- parameters_diagnostics.f90 
- parameters_multibox.f90 
- parameters_numerical.f90
- parameters_physics.f90 

<br>

### Radial variation 

The `radial_variation` folder contains the following scripts:

- gk_radial_variation.f90 
- gk_sources.fpp
- multibox.f90

<br>

### Read Namelists

- namelist_debug.f90
- namelist_diagnostics.f90
- namelist_dissipation.f90
- namelist_flow_shear.f90 
- namelist_geometry.f90 
- namelist_initialise_distribution_function.f90
- namelist_kxky_grid.f90
- namelist_neoclassical_input.f90
- namelist_parameters_numerical.f90
- namelist_parallelisation.f90
- namelist_parameters_physics.f90
- namelist_radial_variation.f90
- namelist_species.f90
- namelist_velocity_grids.f90
- namelist_z_grid.f90
