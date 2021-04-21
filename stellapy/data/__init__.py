#=====================================================================
# Read variables from .txt, .h5 or .nc files
#=====================================================================
''' READ PACKAGE

Read a file from stella and return a dictionary containing the corresponding information.
All read functions require the <folder> and <file_name> as arguments.


Modules for stella files
-------------------------
input: reads the ".in" file
    read.input_species1     -->     dict[z; mass; temp; tprim; frim]
    read.input_species2     -->     {z; mass; temp; tprim; frim}
    read.input_parameters   -->     {torflux; rho; ns; tite; nine; vnew_ref; beta; zeff}
    read.input_grid         -->     {nzed; nperiod; zeta_center; nvgrid; vpa_max; vperp_max; delt; nstep; cfl_cushion; vperp_max; \
                                    nfield_periods; geo_option; boundary_option; adiabatic_option; vmec_filename}
                                    
fluxes: read the ".fluxes" file
    read.fluxes             -->     dict[/;shot][/;rho][time; q_flux; p_flux; v_flux]


vmecgeo: reads the ".vmecgeo" file
    read.vmecgeo            -->     {rhotor; qinp; shat; aref; Bref; z_scalefac}    

netcdf: reads the ".out" or ".h5" file       
    read.netcdf             -->     {dim_species; vec_time; dim_time; vec_kx; dim_kx; vec_ky; dim_ky; vec_density; phi2_vs_kxky; vec_z; dim_z; iz0}
    read.netcdf_geom..      -->     {bmag; gradpar; gbdrift; gbdrift0; cvdrift; cvdrift0; gds2; gds21; gds22; time}

wout: read the ".wout" or ".h5" file
    read.wout               -->     {iota, ns, b0}




Modules for self-made output files: 
------------------------------------
profile:
    read.profile            -->     {s1_dens; s1_temp}

lineardata:         
    read.lineardata[0]      -->     linear_data[shot][rho][kx]['<quantity>'] if nesteddict = True;  else    linear_data['<quantity>']
                                    <quantity> = {'ky', 'omega', 'gamma', 'omega_full', 'gamma_full', 'phi2', 'p_flux', 'q_flux', 'v_flux}

    read.lineardata[1]      -->     input_values[shot][rho]['<quantity>'] if nesteddict = True;     else    input_values['<quantity>']     
                                    <quantity> = {'ns', 'rho', 'tite', 'ref_a', 'ref_B', 's1_z', 's1_mass', 
                                                  's1_dens', 's1_temp', 's1_tprim', 's1_fprim', 'shot_number'}

Functions
---------
check_inputs
    Prints some information about the input parameters to the command prompt.


'''




#=====================================================================
# Make function avaliable as package.func instead of package.mod.func
#=====================================================================

# Import the current modules and sub-package list in the package folder
import os, glob

divider = '\\' if (os.name == 'nt') else '/'
mod_list = [file_name.split(divider)[-1].split('.')[0] for file_name in glob.glob(__path__[0]+'/[!_]*.py')]
sub_pack_list = [folder_name.split(divider)[-2] for folder_name in glob.glob(__path__[0]+'/[!_]*/')]

# Import all functions from the modules
for mod in mod_list:
    exec('from . import ' + mod)
    exec('from .' + mod + ' import *')

# Import all subpackages
for pack in sub_pack_list:
    exec('from . import ' + pack)

# Clean up
del glob
try:
    del mod
except:
    pass
try:
    del pack
except:
    pass
