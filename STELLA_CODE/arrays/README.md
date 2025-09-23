# Arrays

This directory contains modules that store and initialise arrays that are used throughout the code. 

The following modules act as a storage for frequently used arrays:

- arrays_distribution_function.f90 -> Stores distribution sized arrays
- arrays_fields.f90 -> Stores field sized arrays
- arrays_gyro_averages.f90 -> Stores Bessel function arrays
- arrays.f90 -> Stores other useful arrays used throughout stella, such as kperp, vperp etc. 


The following modules initilase arrays:

- arrays_gyro_averages.f90 -> Initialises Bessel function arrays
- initialise_arrays.f90 -> Initialises many of the arrays in arrays.f90 
- initialise_distribution_function.f90 -> Initialises the distribution function arrays