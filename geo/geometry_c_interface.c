/* This file exists to get round the worst of the 
 * name mangling issues when using geo.a with CUDA 
 * Maintainer: edmundhighcock@sourceforge.net */

#include "geometry_c_interface.h"

extern void geometry_set_inputs_(int * equilibrium_type,
		 											 char * eqfile,
													 int * irho,
		 											 double * rhoc, 
													 int * bishop,
													 int * nperiod,
													 int * ntheta_out);
 

/*void allocate_coefficients(int ntgrid){*/
/**/
/*}*/

extern void geometry_get_default_advanced_parameters_(
		struct advanced_parameters_struct * advanced_parameters_out);

void geometry_get_default_advanced_parameters_c(
		struct advanced_parameters_struct * advanced_parameters_out){

				geometry_get_default_advanced_parameters_(
						advanced_parameters_out);
}


extern void geometry_set_advanced_parameters_(
		struct advanced_parameters_struct * advanced_parameters_out);

void geometry_set_advanced_parameters_c(
	struct advanced_parameters_struct * advanced_parameters_out){

	geometry_set_advanced_parameters_(advanced_parameters_out);
}


extern void geometry_get_miller_parameters_(
		struct miller_parameters_struct * miller_parameters_out);

void geometry_get_miller_parameters_c(
		struct miller_parameters_struct * miller_parameters_out){

			geometry_get_miller_parameters_(miller_parameters_out);
}

extern void geometry_set_miller_parameters_(
		struct miller_parameters_struct * miller_parameters_in);

void geometry_set_miller_parameters_c(
		struct miller_parameters_struct * miller_parameters_in){

			geometry_set_miller_parameters_(miller_parameters_in);
}

extern void geometry_get_constant_coefficients_(
		struct constant_coefficients_struct * constant_coefficients_out);

void geometry_get_constant_coefficients_c(
		struct constant_coefficients_struct * constant_coefficients_out){

			geometry_get_constant_coefficients_(constant_coefficients_out);
}

extern void geometry_vary_s_alpha_(
		double * s_hat_input_in, double *beta_prime_input_in);

void geometry_vary_s_alpha_c(double * s_hat_input_in, double * beta_prime_input_in){
	geometry_vary_s_alpha_(s_hat_input_in, beta_prime_input_in);
}

extern void geometry_calculate_coefficients_(
		int * grid_size_out);

void geometry_calculate_coefficients_c(int * grid_size_out){
	geometry_calculate_coefficients_(grid_size_out);
}

extern void geometry_get_coefficients_(int * grid_size_in, struct coefficients_struct * coefficients_out);
void geometry_get_coefficients_c(int * grid_size_in, struct coefficients_struct * coefficients_out){
		geometry_get_coefficients_(grid_size_in, coefficients_out);
}


void geometry_set_inputs_c(int * equilibrium_type,
		 											 char * eqfile,
													 int * irho,
		 											 double * rhoc, 
													 int * bishop,
													 int * nperiod,
													 int * ntheta_out){
	geometry_set_inputs_( equilibrium_type,
		 											  eqfile,
													  irho,
		 											  rhoc, 
													  bishop,
													  nperiod,
													  ntheta_out);
}


