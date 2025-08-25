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
    + [Calculations](#calculations)
    + [Diagnostics](#diagnostics)
    + [Dissipation](#dissipation)
    + [Fields](#fields)
    + [Geometry](#geometry)
    + [Grids](#grids)
    + [Gyrokinetic_terms](#gyrokinetic-terms)
    + [Neoclassical](#neoclassical)
    + [Parameters](#parameters)
    + [Radial_variation](#radial-variation)

<br>

## Code

The code is organized in the following folders:
- calculations 
- diagnostics  
- dissipation  
- fields  
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

### Calculations

The `calculations` folder contains the following scripts:

- calculations_kxky.f90  
- calculations_redistribute.f90  
- calculations_finite_differences.f90  
- calculations_tofrom_ghf.f90  
- gyro_averages.f90  
- stella_transforms.f90  
- volume_averages.f90

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

- coll_dougherty.f90
- coll_fokkerplanck.f90
- dissipation.f90
- hyper.f90

<br>


### Fields  

The `fields` folder contains the following scripts:

- dist_fn.f90
- fields.fpp
- fields_collisions
- fields_electromagnetic.fpp
- fields_fluxtube.fpp
- fields_fullfluxsurface.fpp
- fields_radialvariation.fpp
- init_g.f90

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

- store_arrays_distribution_fn.f90
- arrays_fields.f90
- common_types.f90
- extended_zgrid.f90
- grids_kxky.f90
- species.f90
- stella_layouts.f90
- stella_time.f90
- velocity_grids.f90
- write_radial_grid.f90
- z_grid.f90

<br>


### Gyrokinetic terms  

The `gyrokinetic_terms` folder contains the following scripts:

- ffs_solve.f90
- flow_shear.f90
- implicit_solve.f90
- mirror_terms.f90
- parallel_streaming.f90
- response_matrix.fpp
- sources.fpp
- time_advance.f90

<br>


### Neoclassical  

The `neoclassical` folder contains the following scripts:

- euterpe_interface.f90
- neoclassical_terms.f90
- sfincs_interface.fpp

<br>


### Parameters  

The `parameters` folder contains the following scripts:

- debug_flags.f90
- parameters_diagnostics.f90  
- parameters_kxky_grid_box.f90  
- parameters_kxky_grid.f90  
- parameters_kxky_grid_range.f90  
- parameters_numerical.f90  
- parameters_physics.f90

<br>

### Radial variation 

The `radial_variation` folder contains the following scripts:

- multibox.f90

<br>


