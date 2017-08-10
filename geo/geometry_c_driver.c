/* This file exists to test and exemplify the interface
 * to the geometry module for a C/cuda program.
 * Maintainer: edmundhighcock@sourceforge.net */

#include <stdlib.h>
#include <stdio.h>

extern void geometry_set_inputs_(int * equilibrium_type,
		 											 char * eqfile,
													 int * irho,
		 											 double * rhoc, 
													 int * bishop,
													 int * nperiod,
													 int * ntheta_out);
 
/* These MUST be in the same order as the type declaration in
 * geometry.f90*/
struct advanced_parameters_struct {
	 	int equal_arc;
    double dp_mult;
    double delrho;                                                                
    double rmin;
    double rmax;
    int isym;
    int in_nt;
    int write_lots;
    int itor;
} ;

/* These MUST be in the same order as the type declaration in
 * geometry.f90*/
struct miller_parameters_struct { 
     double rmaj;
     double R_geo;
     double akappa;
     double akappri;
     double tri;
     double tripri;
     double shift;
     double qinp;
     double shat;
     double asym;
     double asympri;
};

struct coefficients_struct {
		 double  grho;   
		 double  bmag;       
		 double  gradpar;    
		 double  cvdrift;    
		 double  cvdrift0;   
		 double  gbdrift;    
		 double  gbdrift0;   
		 double  cdrift;    
		 double  cdrift0;    
		 double  gbdrift_th; 
		 double  cvdrift_th; 
		 double  gds2;       
		 double  gds21;      
		 double  gds22;      
		 double  gds23;      
		 double  gds24;      
		 double  gds24_noq;  
		 double  jacob;      
		 double  Rplot;      
		 double  Zplot;      
		 double  aplot;      
		 double  Rprime;     
		 double  Zprime;     
		 double  aprime;     
		 double  Uk1;        
		 double  Uk2;        
		 double  Bpol;       
};

/*void allocate_coefficients(int ntgrid){*/
/**/
/*}*/

extern void geometry_get_default_advanced_paramters_(
		struct advanced_parameters_struct * advanced_parameters_out);

extern void geometry_set_advanced_paramters_(
		struct advanced_parameters_struct * advanced_parameters_out);

extern void geometry_get_miller_parameters_(
		struct miller_parameters_struct * miller_parameters_out);
		
extern void geometry_set_miller_parameters_(
		struct miller_parameters_struct * miller_parameters_in);

extern void geometry_vary_s_alpha_(
		double * s_hat_input_in, double *beta_prime_input_in);

extern void geometry_get_coefficients_(int * ntheta, struct coefficients_struct * coefficients_out);

void set_inputs (int tprim, int fprim){

}

void get_outputs (double * pflux, double * qflux) {

}

int main (int argv, char * argc){

	double rhoc, delrho;
	int equilibrium_type, ntheta_out, nperiod;
	int bishop, irho;
	char * eqfile;

	eqfile = (char *)malloc(sizeof(char)*800);
	equilibrium_type = 3;
	rhoc = 0.6;
	delrho = 0.01; /* Leave at this value*/
	nperiod = 1;
	eqfile = "ogyropsi.dat"; 
	ntheta_out = 16;
	bishop = 1;
	irho = 2;

	geometry_set_inputs_(&equilibrium_type, 
											 eqfile, 
											 &irho, 
											 &rhoc, 
											 &bishop, 
											 &nperiod, 
											 &ntheta_out);
	printf("ntheta_out was %d\n", ntheta_out);

	struct advanced_parameters_struct advanced_parameters;
	geometry_get_default_advanced_parameters_(&advanced_parameters);

	printf("default delrho was %f and default write_lots was %d\n", 
					advanced_parameters.delrho,
					advanced_parameters.write_lots);  


	/*advanced_parameters.delrho = 2.0;*/
	geometry_set_advanced_parameters_(&advanced_parameters);

	struct miller_parameters_struct miller_parameters;

	geometry_get_miller_parameters_(&miller_parameters);
	printf("The default s_hat was %f\n\n\n", miller_parameters.shat);

	/*miller_parameters.shat = 5.6;*/
	geometry_set_miller_parameters_(&miller_parameters);

	double s_hat_input, beta_prime_input;
	geometry_vary_s_alpha_(&s_hat_input, &beta_prime_input);

	struct coefficients_struct * coefficients_array;

	printf("Got here!");
	int ntheta_out_2;
	geometry_calculate_coefficients_(&ntheta_out_2);

	coefficients_array = (struct coefficients_struct *)malloc(
			sizeof(struct coefficients_struct)*ntheta_out_2);
	geometry_get_coefficients_(&ntheta_out, coefficients_array);


	printf("Got here!");
	printf("\nbmag was %f %f %f %f %f etc\n ", 
			coefficients_array[0].bmag,
			coefficients_array[1].bmag,
			coefficients_array[2].bmag,
			coefficients_array[3].bmag,
			coefficients_array[4].bmag
			);



	return 0;

}
