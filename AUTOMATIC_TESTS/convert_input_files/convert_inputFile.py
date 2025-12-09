''' 

#===============================================================================
#            Upgrade an input file from stella version v0.5 to v1.0            #
#===============================================================================

Use the following command line in a folder with input files:
   >> python3 $STELLA/AUTOMATIC_TESTS/convert_input_files/convert_inputFile.py

Hanne Thienpondt
09/07/2025

'''

import copy
import os, sys
import pathlib
import numpy as np
import f90nml, json

#--------------------------------------------------------------
def update_inputFile(path_input_file='', add_default_variables=False, downgrade=False):
    ''' Update an input file from stella version 5, 6, 7 or 8 to the master version. '''
    
    # Read the input file
    input_parameters = f90nml.read(path_input_file)
    input_parameters = json.loads(json.dumps(input_parameters))  
    
    # We can upgrade the input file from stella version v0.5 to v1.0
    # Or we can downgrade the input file from stella version v1.0 to v0.5
    upgrade = not downgrade
    
    # Make sure everything is in lower case
    for namelist in input_parameters.keys():
        if namelist.lower()!=namelist:
            input_parameters[namelist.lower()] = input_parameters[namelist]
            del input_parameters[namelist]
        namelist = namelist.lower()
        for variable in input_parameters[namelist].keys():
            if variable.lower()!=variable:
                input_parameters[namelist][variable.lower()] = input_parameters[namelist][variable] 
                del input_parameters[namelist][variable] 
            
    #===============================================================================
    #                                   Geometry                                  
    #===============================================================================
    
    renamed_variables = { 
        #------------------- geo_knobs --> geometry_options ------------------
        'geo_knobs:geo_option:local'                : 'geometry_options:geometry_option:local',
        'geo_knobs:q_as_x:False'                    : 'geometry_options:q_as_x:False',
        #------------------- geo_knobs --> geometry_from_txt ------------------
        'geo_knobs:overwrite_bmag:False'            : 'geometry_from_txt:overwrite_bmag:False',
        'geo_knobs:overwrite_gradpar:False'         : 'geometry_from_txt:overwrite_b_dot_grad_zeta:False',
        'geo_knobs:overwrite_gds2:False'            : 'geometry_from_txt:overwrite_grady_dot_grady:False',
        'geo_knobs:overwrite_gds21:False'           : 'geometry_from_txt:overwrite_gradx_dot_grady:False',
        'geo_knobs:overwrite_gds22:False'           : 'geometry_from_txt:overwrite_gradx_dot_gradx:False',
        'geo_knobs:overwrite_gbdrift:False'         : 'geometry_from_txt:overwrite_b_times_gradb_dot_grady:False',
        'geo_knobs:overwrite_cvdrift:False'         : 'geometry_from_txt:overwrite_b_times_kappa_dot_grady:False',
        'geo_knobs:overwrite_gbdrift0:False'        : 'geometry_from_txt:overwrite_b_times_gradb_dot_gradx:False',
        'geo_knobs:overwrite_gds23:False'           : 'geometry_from_txt:overwrite_gds23:False',
        'geo_knobs:overwrite_gds24:False'           : 'geometry_from_txt:overwrite_gds24:False',
        'geo_knobs:set_bmag_const:False'            : 'geometry_from_txt:set_bmag_const:DEPRECATED',
        'geo_knobs:geo_file:input.geometry'         : 'geometry_from_txt:geometry_file:input.geometry',
        #------------------- vmec_parameters --> geometry_vmec ------------------
        'vmec_parameters:vmec_filename:wout*.nc'    : 'geometry_vmec:vmec_filename:wout*.nc',
        'vmec_parameters:alpha0:0.0'                : 'geometry_vmec:alpha0:0.0',
        'vmec_parameters:zeta_center:0.0'           : 'geometry_vmec:zeta_center:0.0',
        'vmec_parameters:nfield_periods:-1.0'       : 'geometry_vmec:nfield_periods:-1.0',
        'vmec_parameters:torflux:0.6354167'         : 'geometry_vmec:torflux:0.6354167',
        'vmec_parameters:surface_option:0'          : 'geometry_vmec:surface_option:0',
        'vmec_parameters:verbose:True'              : 'geometry_vmec:verbose:True',
        'vmec_parameters:gradpar_zeta_prefac:1.0'   : 'vmec_parameters:gradpar_zeta_prefac:DEPRECATED',
        'vmec_parameters:zgrid_scalefac:1.0'        : 'vmec_parameters:zgrid_scalefac:DEPRECATED',
        'geometry_vmec:z_grid_refinement_factor:DOESNT EXIST YET' : 'geometry_vmec:z_grid_refinement_factor:1.0',
        'geometry_vmec:rectangular_cross_section:DOESNT EXIST YET' : 'geometry_vmec:rectangular_cross_section:False',
        'geometry_vmec:radial_coordinate:DOESNT EXIST YET' : 'geometry_vmec:radial_coordinate:sgn(psi_t)psi_t',
        'geometry_vmec:n_tolerated_test_arrays_inconsistencies:DOESNT EXIST YET' : 'geometry_vmec:n_tolerated_test_arrays_inconsistencies:0',
        #------------------- millergeo_parameters --> geometry_miller ------------------
        'millergeo_parameters:nzed_local:128'       : 'geometry_miller:nzed_local:128',
        'millergeo_parameters:rhoc:0.5'             : 'geometry_miller:rhoc:0.5',
        'millergeo_parameters:rmaj:3.0'             : 'geometry_miller:rmaj:3.0',
        'millergeo_parameters:rgeo:3.0'             : 'geometry_miller:rgeo:3.0',
        'millergeo_parameters:qinp:1.4'             : 'geometry_miller:qinp:1.4',
        'millergeo_parameters:shat:0.8'             : 'geometry_miller:shat:0.8',
        'millergeo_parameters:shift:0.0'            : 'geometry_miller:shift:0.0',
        'millergeo_parameters:kappa:0.0'            : 'geometry_miller:kappa:0.0',
        'millergeo_parameters:kapprim:0.0'          : 'geometry_miller:kapprim:0.0',
        'millergeo_parameters:tri:0.0'              : 'geometry_miller:tri:0.0',
        'millergeo_parameters:triprim:0.0'          : 'geometry_miller:triprim:0.0',
        'millergeo_parameters:betaprim:0.0'         : 'geometry_miller:betaprim:0.0',
        'millergeo_parameters:betadbprim:0.0'       : 'geometry_miller:betadbprim:0.0',
        'millergeo_parameters:d2qdr2:0.0'           : 'geometry_miller:d2qdr2:0.0',
        'millergeo_parameters:d2psidr2:0.0'         : 'geometry_miller:d2psidr2:0.0',
        'millergeo_parameters:read_profile_variation:False' : 'geometry_miller:read_profile_variation:False',
        'millergeo_parameters:write_profile_variation:False' : 'geometry_miller:write_profile_variation:False',
        'millergeo_parameters:load_psi0_variables:True' : 'millergeo_parameters:load_psi0_variables:DEPRECATED',
        'millergeo_parameters:rhotor:rhoc'          : 'millergeo_parameters:rhotor:DEPRECATED',
        'millergeo_parameters:psitor_lcfs:1.0'      : 'millergeo_parameters:psitor_lcfs:DEPRECATED',
        'millergeo_parameters:drhotordrho:1.0'      : 'millergeo_parameters:drhotordrho:DEPRECATED',
        'millergeo_parameters:rhoc0:0.5'            : 'millergeo_parameters:rhoc0:DEPRECATED',
        #------------------- millergeo_parameters --> geometry_miller ------------------
        'geometry_zpinch:betaprim:DOESNT EXIST YET' : 'geometry_zpinch:betaprim:0.0',
        
        }
    
    # Replace variables
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    # Deal with stella v0.8
    renamed_variables = { 
        'geo_knobs:overwrite_b_dot_grad_zeta:False'  : 'geometry_from_txt:overwrite_b_dot_grad_zeta:False',
        }
        
    # Replace variables
    if upgrade:
        input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    #===============================================================================
    #                                     Physics                                  
    #===============================================================================
    
    renamed_variables = { 
        #------------------- time_advance_knobs --> numerical_algorithms ------------------
        'time_advance_knobs:explicit_option:rk3'    : 'numerical_algorithms:explicit_algorithm:rk3',
        'time_advance_knobs:flip_flop:False'        : 'numerical_algorithms:flip_flop:False',
        'numerical_algorithms:stream_iterative_implicit:DOESNT EXIST YET' : 'numerical_algorithms:stream_iterative_implicit:False',
        'numerical_algorithms:fully_implicit:DOESNT EXIST YET' : 'numerical_algorithms:fully_implicit:False',
        'numerical_algorithms:fully_explicit:DOESNT EXIST YET' : 'numerical_algorithms:fully_explicit:False',
        'numerical_algorithms:split_parallel_dynamics:DOESNT EXIST YET' : 'numerical_algorithms:split_parallel_dynamics:False',
        'numerical_algorithms:use_deltaphi_for_response_matrix:DOESNT EXIST YET' : 'numerical_algorithms:use_deltaphi_for_response_matrix:False',
        'numerical_algorithms:maxwellian_normalization:DOESNT EXIST YET' : 'numerical_algorithms:maxwellian_normalization:False',
        #------------------- time_advance_knobs --> scale_gyrokinetic_terms ------------------
        'time_advance_knobs:xdriftknob:1.0'         : 'scale_gyrokinetic_terms:xdriftknob:1.0',
        'time_advance_knobs:ydriftknob:1.0'         : 'scale_gyrokinetic_terms:ydriftknob:1.0',
        'time_advance_knobs:wstarknob:1.0'          : 'scale_gyrokinetic_terms:wstarknob:1.0',
        'scale_gyrokinetic_terms:suppress_zonal_interaction:DOESNT EXIST YET' : 'scale_gyrokinetic_terms:suppress_zonal_interaction:False',
        #------------------- physics_flags --> gyrokinetic_terms ------------------
        'physics_flags:radial_variation:False'      : 'gyrokinetic_terms:include_radial_variation:False',
        'physics_flags:include_parallel_nonlinearity:False' : 'gyrokinetic_terms:include_parallel_nonlinearity:False',
        'physics_flags:include_parallel_streaming:True' : 'gyrokinetic_terms:include_parallel_streaming:True',
        'physics_flags:include_mirror:True'         : 'gyrokinetic_terms:include_mirror:True',
        'physics_flags:nonlinear:False'             : 'gyrokinetic_terms:include_nonlinear:False',
        'gyrokinetic_terms:include_electromagnetic:DOESNT EXIST YET' : 'gyrokinetic_terms:include_electromagnetic:False',
        #------------------- physics_flags --> debug_flags ------------------
        'physics_flags:const_alpha_geo:False'       : 'debug_flags:const_alpha_geo:False',
        #------------------- physics_flags --> multibox_parameters ------------------
        'physics_flags:include_pressure_variation:True' : 'multibox_parameters:include_pressure_variation:True',
        'physics_flags:include_geometric_variation:True' : 'multibox_parameters:include_geometric_variation:True',
        #------------------- physics_flags --> adiabatic_electron_response ------------------
        'physics_flags:adiabatic_option:no-field-line-average-term' : 'adiabatic_electron_response:adiabatic_option:no-field-line-average-term',
        #------------------- parameters --> adiabatic_electron_response ------------------
        'parameters:tite:1.0'                       : 'adiabatic_electron_response:tite:1.0',
        'parameters:nine:1.0'                       : 'adiabatic_electron_response:nine:1.0',
        #------------------- parameters --> electromagnetic ------------------
        'parameters:beta:0.0'                       : 'electromagnetic:beta:0.0',
        #------------------- parameters --> physics_inputs ------------------
        'parameters:vnew_ref:-1.0'                  : 'dissipation_and_collisions_options:vnew_ref:-1.0',
        'parameters:zeff:1.0'                       : 'dissipation_and_collisions_options:zeff:1.0',
        'parameters:rhostar:-1.0'                   : 'physics_inputs:rhostar:-1.0',
        #------------------- parameters --> flow_shear ------------------
        'parameters:g_exb:0.0'                      : 'flow_shear:g_exb:0.0',
        'parameters:g_exbfac:1.0'                   : 'flow_shear:g_exbfac:1.0',
        'parameters:omprimfac:1.0'                  : 'flow_shear:omprimfac:1.0',
        'flow_shear:prp_shear_enabled:DOESNT EXIST YET' : 'flow_shear:prp_shear_enabled:False',
        'flow_shear:hammett_flow_shear:DOESNT EXIST YET' : 'flow_shear:hammett_flow_shear:True',
        #------------------- parameters --> flux_annulus ------------------
        'flux_annulus:nitt:DOESNT EXIST YET'        : 'flux_annulus:nitt:1',
        #------------------- parameters --> deprecated ------------------
        'parameters:irhostar:-1.0'                  : 'parameters:irhostar:DEPRECATED',
        }
    
    # Replace variables
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    # Add the default variables
    if upgrade:
        if add_default_variables:
            if 'gyrokinetic_terms' not in input_parameters.keys(): input_parameters['gyrokinetic_terms'] = {}
            input_parameters['gyrokinetic_terms']['include_xdrift'] = True
            input_parameters['gyrokinetic_terms']['include_ydrift'] = True
            input_parameters['gyrokinetic_terms']['include_drive'] = True
    
    # Manually change some variables
    if upgrade:
        if 'time_advance_knobs' in input_parameters.keys():
            if 'xdriftknob' in input_parameters['time_advance_knobs'].keys():
                value_old = input_parameters['time_advance_knobs']['xdriftknob']
                del input_parameters['time_advance_knobs']['xdriftknob']
                if 'gyrokinetic_terms' not in input_parameters.keys(): input_parameters['gyrokinetic_terms'] = {}
                if 'gyrokinetic_terms' not in input_parameters.keys(): input_parameters['scale_gyrokinetic_terms'] = {}
                input_parameters['scale_gyrokinetic_terms']['xdriftknob'] = value_old
                value_old = True if value_old > 0 else False
                input_parameters['gyrokinetic_terms']['include_xdrift'] = value_old
            if 'ydriftknob' in input_parameters['time_advance_knobs'].keys():
                value_old = input_parameters['time_advance_knobs']['ydriftknob']
                del input_parameters['time_advance_knobs']['ydriftknob']
                if 'gyrokinetic_terms' not in input_parameters.keys(): input_parameters['gyrokinetic_terms'] = {}
                if 'gyrokinetic_terms' not in input_parameters.keys(): input_parameters['scale_gyrokinetic_terms'] = {}
                input_parameters['scale_gyrokinetic_terms']['ydriftknob'] = value_old
                value_old = True if value_old > 0 else False
                input_parameters['gyrokinetic_terms']['include_ydrift'] = value_old
            if 'wstarknob' in input_parameters['time_advance_knobs'].keys():
                value_old = input_parameters['time_advance_knobs']['wstarknob']
                del input_parameters['time_advance_knobs']['wstarknob']
                if 'gyrokinetic_terms' not in input_parameters.keys(): input_parameters['gyrokinetic_terms'] = {}
                if 'gyrokinetic_terms' not in input_parameters.keys(): input_parameters['scale_gyrokinetic_terms'] = {}
                input_parameters['scale_gyrokinetic_terms']['wstarknob'] = value_old
                value_old = True if value_old > 0 else False
                input_parameters['gyrokinetic_terms']['include_drive'] = value_old
    if downgrade:
        if 'gyrokinetic_terms' in input_parameters.keys():
            if 'include_xdrift' in input_parameters['gyrokinetic_terms'].keys():
                value_old = input_parameters['gyrokinetic_terms']['include_xdrift']
                del input_parameters['gyrokinetic_terms']['include_xdrift']
                if 'time_advance_knobs' not in input_parameters.keys(): input_parameters['time_advance_knobs'] = {}
                value_old = 0 if value_old==False else 1
                input_parameters['time_advance_knobs']['xdriftknob'] = value_old
            if 'include_ydrift' in input_parameters['gyrokinetic_terms'].keys():
                value_old = input_parameters['gyrokinetic_terms']['include_ydrift']
                del input_parameters['gyrokinetic_terms']['include_ydrift']
                if 'time_advance_knobs' not in input_parameters.keys(): input_parameters['time_advance_knobs'] = {}
                value_old = 0 if value_old==False else 1
                input_parameters['time_advance_knobs']['ydriftknob'] = value_old
                input_parameters['gyrokinetic_terms']['include_ydrift'] = value_old
            if 'include_drive' in input_parameters['gyrokinetic_terms'].keys():
                value_old = input_parameters['gyrokinetic_terms']['include_drive']
                del input_parameters['gyrokinetic_terms']['include_drive']
                if 'time_advance_knobs' not in input_parameters.keys(): input_parameters['time_advance_knobs'] = {}
                value_old = True if value_old > 0 else False
                input_parameters['time_advance_knobs']['wstarknob'] = value_old
    
    # Manually change some variables
    if upgrade:
        if 'physics_flags' in input_parameters.keys():
            if 'full_flux_surface' in input_parameters['physics_flags'].keys():
                value_old = input_parameters['physics_flags']['full_flux_surface']
                if value_old==True: value_old = 'full_flux_annulus'
                if value_old==False: value_old = 'fluxtube'
                del input_parameters['physics_flags']['full_flux_surface']
                if 'gyrokinetic_terms' not in input_parameters.keys(): input_parameters['gyrokinetic_terms'] = {}
                input_parameters['gyrokinetic_terms']['simulation_domain'] = value_old
        if 'parameters_physics' in input_parameters.keys():
            if 'full_flux_surface' in input_parameters['parameters_physics'].keys():
                value_old = input_parameters['parameters_physics']['full_flux_surface']
                if value_old==True: value_old = 'full_flux_annulus'
                if value_old==False: value_old = 'fluxtube'
                del input_parameters['parameters_physics']['full_flux_surface']
                if 'gyrokinetic_terms' not in input_parameters.keys(): input_parameters['gyrokinetic_terms'] = {}
                input_parameters['gyrokinetic_terms']['simulation_domain'] = value_old
    if downgrade:
        if 'gyrokinetic_terms' in input_parameters.keys():
            if 'simulation_domain' in input_parameters['gyrokinetic_terms'].keys():
                value_old = input_parameters['gyrokinetic_terms']['simulation_domain']
                if value_old=='full_flux_annulus': value_old = True
                if value_old!='full_flux_annulus': value_old = False
                del input_parameters['gyrokinetic_terms']['simulation_domain']
                if 'physics_flags' not in input_parameters.keys(): input_parameters['physics_flags'] = {}
                input_parameters['physics_flags']['full_flux_surface'] = value_old
                
                
    # Deal with stella v0.8
    renamed_variables = { 
        #-------------------------- parameters_physics -------------------------
        'parameters_physics:explicit_option:rk3'    : 'numerical_algorithms:explicit_algorithm:rk3',
        'parameters_physics:flip_flop:False'        : 'numerical_algorithms:flip_flop:False',
        'parameters_physics:stream_iterative_implicit:DOESNT EXIST YET' : 'numerical_algorithms:stream_iterative_implicit:False',
        'parameters_physics:fully_implicit:DOESNT EXIST YET' : 'numerical_algorithms:fully_implicit:False',
        'parameters_physics:fully_explicit:DOESNT EXIST YET' : 'numerical_algorithms:fully_explicit:False',
        'parameters_physics:split_parallel_dynamics:DOESNT EXIST YET' : 'numerical_algorithms:split_parallel_dynamics:False',
        'parameters_physics:use_deltaphi_for_response_matrix:DOESNT EXIST YET' : 'numerical_algorithms:use_deltaphi_for_response_matrix:False',
        'parameters_physics:maxwellian_normalization:DOESNT EXIST YET' : 'numerical_algorithms:maxwellian_normalization:False',
        'parameters_physics:xdriftknob:1.0'         : 'scale_gyrokinetic_terms:xdriftknob:1.0',
        'parameters_physics:ydriftknob:1.0'         : 'scale_gyrokinetic_terms:ydriftknob:1.0',
        'parameters_physics:wstarknob:1.0'          : 'scale_gyrokinetic_terms:wstarknob:1.0',
        'parameters_physics:suppress_zonal_interaction:DOESNT EXIST YET' : 'scale_gyrokinetic_terms:suppress_zonal_interaction:False',
        'parameters_physics:radial_variation:False'      : 'gyrokinetic_terms:include_radial_variation:False',
        'parameters_physics:include_parallel_nonlinearity:False' : 'gyrokinetic_terms:include_parallel_nonlinearity:False',
        'parameters_physics:include_parallel_streaming:True' : 'gyrokinetic_terms:include_parallel_streaming:True',
        'parameters_physics:include_mirror:True'         : 'gyrokinetic_terms:include_mirror:True',
        'parameters_physics:nonlinear:False'             : 'gyrokinetic_terms:include_nonlinear:False',
        'parameters_physics:include_electromagnetic:DOESNT EXIST YET' : 'gyrokinetic_terms:include_electromagnetic:False',
        'parameters_physics:const_alpha_geo:False'       : 'debug_flags:const_alpha_geo:False',
        'parameters_physics:include_pressure_variation:True' : 'multibox_parameters:include_pressure_variation:True',
        'parameters_physics:include_geometric_variation:True' : 'multibox_parameters:include_geometric_variation:True',
        'parameters_physics:adiabatic_option:no-field-line-average-term' : 'adiabatic_electron_response:adiabatic_option:no-field-line-average-term',
        'parameters_physics:tite:1.0'                       : 'adiabatic_electron_response:tite:1.0',
        'parameters_physics:nine:1.0'                       : 'adiabatic_electron_response:nine:1.0',
        'parameters_physics:beta:0.0'                       : 'electromagnetic:beta:0.0',
        'parameters_physics:vnew_ref:-1.0'                  : 'dissipation_and_collisions_options:vnew_ref:-1.0',
        'parameters_physics:zeff:1.0'                       : 'dissipation_and_collisions_options:zeff:1.0',
        'parameters_physics:rhostar:-1.0'                   : 'physics_inputs:rhostar:-1.0',
        'parameters_physics:g_exb:0.0'                      : 'flow_shear:g_exb:0.0',
        'parameters_physics:g_exbfac:1.0'                   : 'flow_shear:g_exbfac:1.0',
        'parameters_physics:omprimfac:1.0'                  : 'flow_shear:omprimfac:1.0',
        'parameters_physics:prp_shear_enabled:DOESNT EXIST YET' : 'flow_shear:prp_shear_enabled:False',
        'parameters_physics:hammett_flow_shear:DOESNT EXIST YET' : 'flow_shear:hammett_flow_shear:True',
        'parameters_physics:nitt:DOESNT EXIST YET'        : 'flux_annulus:nitt:1',
        'parameters_physics:irhostar:-1.0'                  : 'parameters:irhostar:DEPRECATED',
        }
        
    # Replace variables
    if upgrade:
        input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    #===============================================================================
    #                                 Kinetic species                                  
    #===============================================================================
    
    renamed_variables = { 
        #------------------- species_knobs --> species_options ------------------
        'species_knobs:nspec:2' : 'species_options:nspec:2',
        'species_knobs:read_profile_variation:False' : 'species_options:read_profile_variation:False',
        'species_knobs:write_profile_variation:False' : 'species_options:write_profile_variation:False',
        'species_knobs:species_option:stella' : 'species_options:species_option:stella',
        'species_knobs:ecoll_zeff:False' : 'dissipation_and_collisions_options:ecoll_zeff:False'}
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    # Do this for nspec species
    nspec = 2
    if 'species_knobs' in input_parameters.keys():
        if 'nspec' in input_parameters['species_knobs'].keys():
            nspec = int(input_parameters['species_knobs']['nspec'])
    for i in range(nspec):
        ispec = i+1
        renamed_variables = { 
            #------------------- species_parameters_1 ------------------
            f'species_parameters_{ispec}:z:1' : f'species_parameters_{ispec}:z:1',
            f'species_parameters_{ispec}:mass:1.0' : f'species_parameters_{ispec}:mass:1.0',
            f'species_parameters_{ispec}:dens:1.0' : f'species_parameters_{ispec}:dens:1.0',
            f'species_parameters_{ispec}:temp:1.0' : f'species_parameters_{ispec}:temp:1.0',
            f'species_parameters_{ispec}:tprim:-999.9' : f'species_parameters_{ispec}:tprim:-999.9',
            f'species_parameters_{ispec}:fprim:-999.9' : f'species_parameters_{ispec}:fprim:-999.9',
            f'species_parameters_{ispec}:d2ndr2:0.0' : f'species_parameters_{ispec}:d2ndr2:0.0',
            f'species_parameters_{ispec}:d2tdr2:0.0' : f'species_parameters_{ispec}:d2tdr2:0.0',
            f'species_parameters_{ispec}:bess_fac:1.0' : f'species_parameters_{ispec}:bess_fac:1.0',
            f'species_parameters_{ispec}:type:ion' : f'species_parameters_{ispec}:type:ion',
            }

        # Replace variables
        input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    #===============================================================================
    #                       Discretized (kx,ky,z,mu,vpa) grid                       
    #===============================================================================
    
    renamed_variables = { 
        #------------------- kt_grids_knobs --> kxky_grid_option ------------------
        'kt_grids_knobs:grid_option:range'          : 'kxky_grid_option:grid_option:range',
        #------------------- kt_grids_range_parameters --> kxky_grid_range------------------
        'kt_grids_range_parameters:nalpha:1'        : 'kt_grids_range_parameters:nalpha:DEPRECATED',
        'kt_grids_range_parameters:naky:1'          : 'kxky_grid_range:naky:1',
        'kt_grids_range_parameters:nakx:1'          : 'kxky_grid_range:nakx:1',
        'kt_grids_range_parameters:aky_min:0.0'     : 'kxky_grid_range:aky_min:0.0',
        'kt_grids_range_parameters:aky_max:0.0'     : 'kxky_grid_range:aky_max:0.0',
        'kt_grids_range_parameters:akx_min:0.0'     : 'kxky_grid_range:akx_min:0.0',
        'kt_grids_range_parameters:akx_max:-1.0'    : 'kxky_grid_range:akx_max:-1.0',
        'kt_grids_range_parameters:theta0_min:0.0'  : 'kxky_grid_range:theta0_min:0.0',
        'kt_grids_range_parameters:theta0_max:-1.0' : 'kxky_grid_range:theta0_max:-1.0',
        'kt_grids_range_parameters:phase_shift_angle:0.0' : 'kt_grids_range_parameters:phase_shift_angle:DEPRECATED',
        'kxky_grid_range:kyspacing_option:DOESNT EXIST YET' : 'kxky_grid_range:kyspacing_option:default',
        #------------------- kt_grids_box_parameters ------------------
        'kt_grids_box_parameters:nx:1'              : 'kxky_grid_box:nx:1',
        'kt_grids_box_parameters:ny:1'              : 'kxky_grid_box:ny:1',
        'kt_grids_box_parameters:jtwist:-1'         : 'kxky_grid_box:jtwist:-1',
        'kt_grids_box_parameters:jtwistfac:1.0'     : 'kxky_grid_box:jtwistfac:1.0',
        'kt_grids_box_parameters:x0:-1.0'           : 'kxky_grid_box:x0:-1.0',
        'kt_grids_box_parameters:y0:-1.0'           : 'kxky_grid_box:y0:-1.0',
        'kt_grids_box_parameters:nalpha:1'          : 'kt_grids_box_parameters:nalpha:DEPRECATED',
        'kt_grids_box_parameters:centered_in_rho:True' : 'kxky_grid_box:centered_in_rho:True',
        'kt_grids_box_parameters:randomize_phase_shift:False' : 'kxky_grid_box:randomize_phase_shift:False',
        'kt_grids_box_parameters:periodic_variation:False' : 'kxky_grid_box:periodic_variation:False',
        'kt_grids_box_parameters:phase_shift_angle:0.0' : 'kxky_grid_box:phase_shift_angle:0.0',
        #------------------- zgrid_parameters --> z_grid------------------
        'zgrid_parameters:nzed:24'                  : 'z_grid:nzed:24',
        'zgrid_parameters:nperiod:1'                : 'z_grid:nperiod:1',
        'zgrid_parameters:ntubes:1'                 : 'z_grid:ntubes:1',
        'zgrid_parameters:zed_equal_arc:False'      : 'z_grid:zed_equal_arc:False',
        #------------------- zgrid_parameters --> z_boundary_condition------------------
        'zgrid_parameters:boundary_option:zero'     : 'z_boundary_condition:boundary_option:zero',
        'zgrid_parameters:shat_zero:1e-05'          : 'z_boundary_condition:shat_zero:1e-05',
        'zgrid_parameters:grad_x_grad_y_zero:1e-05' : 'z_boundary_condition:grad_x_grad_y_zero:1e-05',
        'z_boundary_condition:dkx_over_dky:DOESNT EXIST YET' : 'z_boundary_condition:dkx_over_dky:-1.0',
        #------------------- vpamu_grids_parameters --> velocity_grids ------------------
        'vpamu_grids_parameters:nvgrid:24'          : 'velocity_grids:nvgrid:24',
        'vpamu_grids_parameters:vpa_max:3.0'        : 'velocity_grids:vpa_max:3.0',
        'vpamu_grids_parameters:nmu:12'             : 'velocity_grids:nmu:12',
        'vpamu_grids_parameters:vperp_max:3.0'      : 'velocity_grids:vperp_max:3.0',
        'vpamu_grids_parameters:equally_spaced_mu_grid:False' : 'velocity_grids:equally_spaced_mu_grid:False',
        'vpamu_grids_parameters:conservative_wgts_vpa:False' : 'velocity_grids:conservative_wgts_vpa:False',
        }
    
    # Replace variables
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    # Deal with stella version 0.8
    renamed_variables = { 
        #------------------- kt_grids_box_parameters ------------------
        'parameters_kxky_grids_box:nx:1'              : 'kxky_grid_box:nx:1',
        'parameters_kxky_grids_box:ny:1'              : 'kxky_grid_box:ny:1',
        'parameters_kxky_grids_box:jtwist:-1'         : 'kxky_grid_box:jtwist:-1',
        'parameters_kxky_grids_box:jtwistfac:1.0'     : 'kxky_grid_box:jtwistfac:1.0',
        'parameters_kxky_grids_box:x0:-1.0'           : 'kxky_grid_box:x0:-1.0',
        'parameters_kxky_grids_box:y0:-1.0'           : 'kxky_grid_box:y0:-1.0',
        'parameters_kxky_grids_box:nalpha:1'          : 'kt_grids_box_parameters:nalpha:DEPRECATED',
        'parameters_kxky_grids_box:centered_in_rho:True' : 'kxky_grid_box:centered_in_rho:True',
        'parameters_kxky_grids_box:randomize_phase_shift:False' : 'kxky_grid_box:randomize_phase_shift:False',
        'parameters_kxky_grids_box:periodic_variation:False' : 'kxky_grid_box:periodic_variation:False',
        'parameters_kxky_grids_box:phase_shift_angle:0.0' : 'kxky_grid_box:phase_shift_angle:0.0',
        }
        
    # Replace variables
    if upgrade:
        input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    #===============================================================================
    #                                   Diagnostics                                  
    #===============================================================================
    
    renamed_variables = { 
        #------------------- stella_diagnostics_knobs --> diagnostics ------------------
        'stella_diagnostics_knobs:nwrite:50'        : 'diagnostics:nwrite:50',
        'stella_diagnostics_knobs:navg:50'          : 'diagnostics:navg:50',
        'stella_diagnostics_knobs:nsave:-1'         : 'diagnostics:nsave:-1',
        'stella_diagnostics_knobs:nc_mult:1'        : 'diagnostics:nc_mult:1',
        'stella_diagnostics_knobs:save_for_restart:False' : 'diagnostics:save_for_restart:False',
        #------------------- stella_diagnostics_knobs --> diagnostics_potential ------------------
        'stella_diagnostics_knobs:write_phi_vs_time:False' : 'diagnostics_potential:write_phi_vs_kxkyz:False',
        #------------------- stella_diagnostics_knobs --> diagnostics_fluxes ------------------
        'stella_diagnostics_knobs:write_radial_fluxes:False' : 'diagnostics_fluxes:write_radial_fluxes:False',
        'stella_diagnostics_knobs:write_fluxes_kxkyz:False' : 'diagnostics_fluxes:write_fluxes_kxkyz:False',
        'stella_diagnostics_knobs:flux_norm:True'   : 'diagnostics_fluxes:flux_norm:True',
        #------------------- stella_diagnostics_knobs --> diagnostics_distribution ------------------
        'stella_diagnostics_knobs:write_gvmus:False' : 'diagnostics_distribution:write_g2_vs_vpamus:False',
        'stella_diagnostics_knobs:write_gzvs:False' : 'diagnostics_distribution:write_g2_vs_zvpas:False',
        #------------------- stella_diagnostics_knobs --> diagnostics_moments ------------------
        'stella_diagnostics_knobs:write_moments:False' : 'diagnostics_moments:write_moments:False',
        'stella_diagnostics_knobs:write_radial_moments:False' : 'diagnostics_moments:write_radial_moments:False',
        }
    
    # Replace variables
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    # Add the default variables
    if add_default_variables:
        if upgrade:
            if 'diagnostics' not in input_parameters.keys(): input_parameters['diagnostics'] = {}
            if 'diagnostics_potential' not in input_parameters.keys(): input_parameters['diagnostics_potential'] = {}
            if 'diagnostics_omega' not in input_parameters.keys(): input_parameters['diagnostics_omega'] = {}
            if 'diagnostics_distribution' not in input_parameters.keys(): input_parameters['diagnostics_distribution'] = {}
            if 'diagnostics_fluxes' not in input_parameters.keys(): input_parameters['diagnostics_fluxes'] = {}
            if 'diagnostics_moments' not in input_parameters.keys(): input_parameters['diagnostics_moments'] = {}
            input_parameters['diagnostics']['nwrite'] = 50.0
            input_parameters['diagnostics']['navg'] = 50.0
            input_parameters['diagnostics']['nsave'] = -1.0
            input_parameters['diagnostics']['nc_mult'] = 1.0
            input_parameters['diagnostics']['save_for_restart'] = False
            input_parameters['diagnostics']['write_all'] = False
            input_parameters['diagnostics']['write_all_time_traces'] = True
            input_parameters['diagnostics']['write_all_spectra_kxkyz'] = False
            input_parameters['diagnostics']['write_all_spectra_kxky'] = False
            input_parameters['diagnostics']['write_all_velocity_space'] = False
            input_parameters['diagnostics']['write_all_potential'] = False
            input_parameters['diagnostics']['write_all_omega'] = False
            input_parameters['diagnostics']['write_all_distribution'] = False
            input_parameters['diagnostics']['write_all_fluxes'] = False
            input_parameters['diagnostics']['write_all_moments'] = False
            input_parameters['diagnostics_potential']['write_all_potential_time_traces'] = False
            input_parameters['diagnostics_potential']['write_all_potential_spectra'] = False
            input_parameters['diagnostics_potential']['write_phi2_vs_time'] = True
            input_parameters['diagnostics_potential']['write_apar2_vs_time'] = True
            input_parameters['diagnostics_potential']['write_bpar2_vs_time'] = True
            input_parameters['diagnostics_potential']['write_phi_vs_kxkyz'] = False
            input_parameters['diagnostics_potential']['write_apar_vs_kxkyz'] = False
            input_parameters['diagnostics_potential']['write_bpar_vs_kxkyz'] = False
            input_parameters['diagnostics_potential']['write_phi2_vs_kxky'] = False
            input_parameters['diagnostics_potential']['write_apar2_vs_kxky'] = False
            input_parameters['diagnostics_potential']['write_bpar2_vs_kxky'] = False
            input_parameters['diagnostics_omega']['write_omega_vs_kxky'] = False
            input_parameters['diagnostics_omega']['write_omega_avg_vs_kxky'] = False
            input_parameters['diagnostics_distribution']['write_g2_vs_vpamus'] = False
            input_parameters['diagnostics_distribution']['write_g2_vs_zvpas'] = False
            input_parameters['diagnostics_distribution']['write_g2_vs_zmus'] = False
            input_parameters['diagnostics_distribution']['write_g2_vs_kxkyzs'] = False
            input_parameters['diagnostics_distribution']['write_g2_vs_zvpamus'] = False
            input_parameters['diagnostics_distribution']['write_distribution_g'] = True
            input_parameters['diagnostics_distribution']['write_distribution_h'] = False
            input_parameters['diagnostics_distribution']['write_distribution_f'] = False
            input_parameters['diagnostics_fluxes']['flux_norm'] = True
            input_parameters['diagnostics_fluxes']['write_fluxes_vs_time'] = True
            input_parameters['diagnostics_fluxes']['write_radial_fluxes'] = False
            input_parameters['diagnostics_fluxes']['write_fluxes_kxkyz'] = False
            input_parameters['diagnostics_fluxes']['write_fluxes_kxky'] = False
            input_parameters['diagnostics_moments']['write_moments'] = False
            input_parameters['diagnostics_moments']['write_radial_moments'] = False
        if downgrade:
            input_parameters['stella_diagnostics_knobs']['write_omega'] = False
            input_parameters['stella_diagnostics_knobs']['write_phi_vs_time'] = False
            input_parameters['stella_diagnostics_knobs']['write_gvmus'] = False
            input_parameters['stella_diagnostics_knobs']['write_gzvs'] = False
            input_parameters['stella_diagnostics_knobs']['write_kspectra'] = False
            input_parameters['stella_diagnostics_knobs']['write_moments'] = False
            input_parameters['stella_diagnostics_knobs']['write_radial_moments'] = False
            input_parameters['stella_diagnostics_knobs']['write_fluxes_kxkyz'] = False
    
    # Do most of the diagnostics manually
    if upgrade:
        if 'stella_diagnostics_knobs' in input_parameters.keys():
            if 'write_omega' in input_parameters['stella_diagnostics_knobs'].keys():
                value_old = input_parameters['stella_diagnostics_knobs']['write_omega']
                del input_parameters['stella_diagnostics_knobs']['write_omega']
                if 'diagnostics_omega' not in input_parameters.keys(): input_parameters['diagnostics_omega'] = {}
                input_parameters['diagnostics_omega']['write_omega_vs_kxky'] = value_old
                input_parameters['diagnostics_omega']['write_omega_avg_vs_kxky'] = value_old
            if 'write_kspectra' in input_parameters['stella_diagnostics_knobs'].keys():
                value_old = input_parameters['stella_diagnostics_knobs']['write_kspectra']
                del input_parameters['stella_diagnostics_knobs']['write_kspectra']
                if 'diagnostics_potential' not in input_parameters.keys(): input_parameters['diagnostics_potential'] = {}
                input_parameters['diagnostics_potential']['write_phi2_vs_kxky'] = value_old
                input_parameters['diagnostics_potential']['write_apar2_vs_kxky'] = value_old
                input_parameters['diagnostics_potential']['write_bpar2_vs_kxky'] = value_old
                if 'diagnostics_fluxes' not in input_parameters.keys(): input_parameters['diagnostics_fluxes'] = {}
                input_parameters['diagnostics_fluxes']['write_fluxes_kxky'] = value_old
                if 'diagnostics_omega' not in input_parameters.keys(): input_parameters['diagnostics_omega'] = {}
                input_parameters['diagnostics_omega']['write_omega_vs_kxky'] = value_old
                input_parameters['diagnostics_omega']['write_omega_avg_vs_kxky'] = value_old
    if downgrade:
        if 'diagnostics' in input_parameters.keys():
            if 'write_all' in input_parameters['diagnostics'].keys():
                write_all = input_parameters['diagnostics']['write_all']
                input_parameters['stella_diagnostics_knobs']['write_omega'] = write_all
                input_parameters['stella_diagnostics_knobs']['write_phi_vs_time'] = write_all
                input_parameters['stella_diagnostics_knobs']['write_gvmus'] = write_all
                input_parameters['stella_diagnostics_knobs']['write_gzvs'] = write_all
                input_parameters['stella_diagnostics_knobs']['write_kspectra'] = write_all
                input_parameters['stella_diagnostics_knobs']['write_moments'] = write_all
                input_parameters['stella_diagnostics_knobs']['write_radial_moments'] = write_all
                input_parameters['stella_diagnostics_knobs']['write_fluxes_kxkyz'] = write_all
        
    
    #===============================================================================
    #                                Initialise potential                                  
    #===============================================================================
    
    renamed_variables = { 
        #------------------- init_g_knobs --> initialise_distribution------------------
        'init_g_knobs:ginit_option:default'         : 'initialise_distribution:initialise_distribution_option:default',
        'init_g_knobs:phiinit:1.0'                  : 'initialise_distribution:phiinit:1.0',
        'init_g_knobs:scale_to_phiinit:False'       : 'initialise_distribution:scale_to_phiinit:False',  
        #------------------- init_g_knobs --> initialise_distribution_noise------------------
        'init_g_knobs:zf_init:1.0'                  : 'initialise_distribution_noise:zf_init:1.0',
        #------------------- init_g_knobs --> initialise_distribution_kpar------------------
        'init_g_knobs:tpar0:0.0'                    : 'initialise_distribution_kpar:tpar0:0.0',
        'init_g_knobs:tperp0:0.0'                   : 'initialise_distribution_kpar:tperp0:0.0',
        'init_g_knobs:den1:0.0'                     : 'initialise_distribution_kpar:den1:0.0',
        'init_g_knobs:upar1:0.0'                    : 'initialise_distribution_kpar:upar1:0.0',
        'init_g_knobs:tpar1:0.0'                    : 'initialise_distribution_kpar:tpar1:0.0',
        'init_g_knobs:tperp1:0.0'                   : 'initialise_distribution_kpar:tperp1:0.0',
        'init_g_knobs:den2:0.0'                     : 'initialise_distribution_kpar:den2:0.0',
        'init_g_knobs:upar2:0.0'                    : 'initialise_distribution_kpar:upar2:0.0',
        'init_g_knobs:tpar2:0.0'                    : 'initialise_distribution_kpar:tpar2:0.0',
        'init_g_knobs:tperp2:0.0'                   : 'initialise_distribution_kpar:tperp2:0.0',
        #------------------- init_g_knobs --> initialise_distribution_rh------------------
        'init_g_knobs:kxmax:1e+100'                 : 'initialise_distribution_rh:kxmax:1e+100',
        'init_g_knobs:kxmin:0.0'                    : 'initialise_distribution_rh:kxmin:0.0',
        #------------------- init_g_knobs --> restart_options------------------
        'init_g_knobs:restart_file:dummy.nc'        : 'restart_options:restart_file:dummy.nc',
        'init_g_knobs:restart_dir:./'               : 'restart_options:restart_dir:./',
        'init_g_knobs:tstart:0.0'                   : 'restart_options:tstart:0.0',
        'init_g_knobs:scale:1.0'                    : 'restart_options:scale:1.0',
        'restart_options:read_many:DOESNT EXIST YET': 'restart_options:save_many:True',
        'restart_options:save_many:DOESNT EXIST YET': 'restart_options:save_many:True',
        }
    
    # Replace variables
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    
    renamed_variables = { 
        'init_g_knobs:width0:-3.5'                  : 'initialise_distribution_maxwellian:width0:-3.5',
        'init_g_knobs:den0:1.0'                     : 'initialise_distribution_maxwellian:den0:1.0',
        'init_g_knobs:upar0:0.0'                    : 'initialise_distribution_maxwellian:upar0:0.0',
        }
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade, delete_variable=False)
    renamed_variables = { 
        'init_g_knobs:width0:-3.5'                  : 'initialise_distribution_kpar:width0:-3.5',
        'init_g_knobs:den0:1.0'                     : 'initialise_distribution_kpar:den0:1.0',
        'init_g_knobs:upar0:0.0'                    : 'initialise_distribution_kpar:upar0:0.0',
        }
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade, delete_variable=True)
    renamed_variables = { 
        'init_g_knobs:chop_side:False'              : 'initialise_distribution_maxwellian:chop_side:False',
        'init_g_knobs:left:True'                    : 'initialise_distribution_maxwellian:left:True',
        }
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade, delete_variable=False)
    renamed_variables = { 
        'init_g_knobs:chop_side:False'              : 'initialise_distribution_kpar:chop_side:False',
        'init_g_knobs:left:True'                    : 'initialise_distribution_kpar:left:True',
        }
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade, delete_variable=False)
    renamed_variables = { 
        'init_g_knobs:chop_side:False'              : 'initialise_distribution_noise:chop_side:False',
        'init_g_knobs:left:True'                    : 'initialise_distribution_noise:left:True',
        }
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade, delete_variable=True)
    renamed_variables = { 
        'init_g_knobs:refac:1.0'                    : 'initialise_distribution_kpar:refac:1.0',
        'init_g_knobs:imfac:0.0'                    : 'initialise_distribution_kpar:imfac:0.0',
        }
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade, delete_variable=False)
    renamed_variables = { 
        'init_g_knobs:refac:1.0'                    : 'initialise_distribution_rh:refac:1.0',
        'init_g_knobs:imfac:0.0'                    : 'initialise_distribution_rh:imfac:0.0',
        }
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade, delete_variable=True)
    
    # Manually change init_g_knobs:even --> initialise_distribution_maxwellian:oddparity
    if upgrade:
        if 'init_g_knobs' in input_parameters.keys():
            if 'even' in input_parameters['init_g_knobs'].keys():
                value_old = input_parameters['init_g_knobs']['even']
                del input_parameters['init_g_knobs']['even']
                if type(value_old)!=bool: print('This should be a boolean. Abort.'); sys.exit()
                if 'initialise_distribution_maxwellian' not in input_parameters.keys(): input_parameters['initialise_distribution_maxwellian'] = {}
                input_parameters['initialise_distribution_maxwellian']['oddparity'] = not value_old
    if downgrade:
        if 'initialise_distribution_maxwellian' in input_parameters.keys():
            if 'oddparity' in input_parameters['initialise_distribution_maxwellian'].keys():
                value_old = input_parameters['initialise_distribution_maxwellian']['oddparity']
                del input_parameters['initialise_distribution_maxwellian']['oddparity']
                if type(value_old)!=bool: print('This should be a boolean. Abort.'); sys.exit()
                if 'init_g_knobs' not in input_parameters.keys(): input_parameters['init_g_knobs'] = {}
                input_parameters['init_g_knobs']['even'] = not value_old
    
    #===============================================================================
    #                                    Numerics                                  
    #===============================================================================
    
    renamed_variables = {  
        #------------------- knobs --> time_trace_options------------------
        'knobs:avail_cpu_time:10000000000.0'        : 'time_trace_options:avail_cpu_time:10000000000.0',
        'knobs:tend:-1.0'                           : 'time_trace_options:tend:-1.0',
        'knobs:nstep:-1'                            : 'time_trace_options:nstep:-1',
        'time_trace_options:autostop:DOESNT EXIST YET' : 'time_trace_options:autostop:True',
        #------------------- knobs --> time_step------------------
        'knobs:delt:-1'                             : 'time_step:delt:0.03',
        'knobs:delt_option:check_restart'           : 'time_step:delt_option:check_restart',
        'knobs:delt_adjust:2.0'                     : 'knobs:delt_adjust:DEPRECATED',
        'knobs:delt_max:-1'                         : 'time_step:delt_max:-1',
        'time_step:delt_min:DOESNT EXIST YET'       : 'time_step:delt_min:1e-10',
        'knobs:cfl_cushion:0.5'                     : 'time_step:cfl_cushion_upper:0.5', 
        'time_step:cfl_cushion_middle:DOESNT EXIST YET' : 'time_step:cfl_cushion_middle:0.25',
        'time_step:cfl_cushion_lower:DOESNT EXIST YET'  : 'time_step:cfl_cushion_lower:1e-05',
        #------------------- knobs --> scale_gyrokinetic_terms------------------
        'knobs:fphi:1.0'                            : 'scale_gyrokinetic_terms:fphi:1.0',
        #------------------- knobs --> parallelisation------------------
        'knobs:fields_kxkyz:False'                  : 'parallelisation:fields_kxkyz:False',
        'knobs:lu_option:none'                      : 'parallelisation:lu_option:default',
        'knobs:mat_gen:True'                        : 'parallelisation:mat_gen:True',
        'knobs:mat_read:False'                      : 'parallelisation:mat_read:False',
        #------------------- knobs --> numerical_algorithms------------------
        'knobs:stream_implicit:True'                : 'numerical_algorithms:stream_implicit:True',
        'knobs:mirror_implicit:True'                : 'numerical_algorithms:mirror_implicit:True',
        'knobs:drifts_implicit:False'               : 'numerical_algorithms:drifts_implicit:False',
        'knobs:driftkinetic_implicit:False'         : 'numerical_algorithms:driftkinetic_implicit:False',
        'knobs:maxwellian_inside_zed_derivative:False' : 'numerical_algorithms:maxwellian_inside_zed_derivative:False',
        'knobs:mirror_semi_lagrange:True'           : 'numerical_algorithms:mirror_semi_lagrange:True',
        'knobs:mirror_linear_interp:False'          : 'numerical_algorithms:mirror_linear_interp:False',
        'knobs:stream_matrix_inversion:False'       : 'numerical_algorithms:stream_matrix_inversion:False',
        #------------------- knobs --> numerical_upwinding_for_derivatives------------------
        'knobs:zed_upwind:0.02'                     : 'numerical_upwinding_for_derivatives:zed_upwind:0.02',
        'knobs:vpa_upwind:0.02'                     : 'numerical_upwinding_for_derivatives:vpa_upwind:0.02',
        'knobs:time_upwind:0.02'                    : 'numerical_upwinding_for_derivatives:time_upwind:0.02',
        #------------------- knobs --> initialise_distribution_noise------------------
        'knobs:rng_seed:-1'                         : 'initialise_distribution_noise:rng_seed:-1',
        #------------------- knobs --> multibox_parameters------------------
        'knobs:ky_solve_radial:0'                   : 'multibox_parameters:ky_solve_radial:0',
        'knobs:ky_solve_real:False'                 : 'multibox_parameters:ky_solve_real:False',
        }
    
    # Replace variables
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    # Add the default variables
    if add_default_variables:
        if 'electromagnetic' not in input_parameters.keys(): input_parameters['electromagnetic'] = {}
        input_parameters['electromagnetic']['include_apar'] = False
        input_parameters['electromagnetic']['include_bpar'] = False
        
    # Manually change some values
    if upgrade:
        if 'knobs' in input_parameters.keys():
            if 'fapar' in input_parameters['knobs'].keys():
                value_old = input_parameters['knobs']['fapar']
                del input_parameters['knobs']['fapar']
                value_old = True if value_old > 0 else False
                if 'electromagnetic' not in input_parameters.keys(): input_parameters['electromagnetic'] = {}
                input_parameters['electromagnetic']['include_apar'] = value_old
            if 'fbpar' in input_parameters['knobs'].keys():
                value_old = input_parameters['knobs']['fbpar']
                del input_parameters['knobs']['fbpar']
                value_old = True if value_old > 0 else False
                if 'electromagnetic' not in input_parameters.keys(): input_parameters['electromagnetic'] = {}
                input_parameters['electromagnetic']['include_bpar'] = value_old
    if downgrade:
        if 'electromagnetic' in input_parameters.keys():
            if 'include_apar' in input_parameters['electromagnetic'].keys():
                value_old = input_parameters['electromagnetic']['include_apar']
                del input_parameters['electromagnetic']['include_apar']
                value_old = 1 if value_old==True else 0
                if 'knobs' not in input_parameters.keys(): input_parameters['knobs'] = {}
                input_parameters['knobs']['fapar'] = value_old
            if 'include_bpar' in input_parameters['electromagnetic'].keys():
                value_old = input_parameters['electromagnetic']['include_bpar']
                del input_parameters['electromagnetic']['include_bpar']
                value_old = 1 if value_old==True else 0
                if 'knobs' not in input_parameters.keys(): input_parameters['knobs'] = {}
                input_parameters['knobs']['fbpar'] = value_old
                
    # Deal with stella v0.8
    renamed_variables = { 
        #-------------------------- parameters_numerical------------------------
        'parameters_numerical:avail_cpu_time:10000000000.0'        : 'time_trace_options:avail_cpu_time:10000000000.0',
        'parameters_numerical:tend:-1.0'                           : 'time_trace_options:tend:-1.0',
        'parameters_numerical:nstep:-1'                            : 'time_trace_options:nstep:-1',
        'parameters_numerical:autostop:True'                       : 'time_trace_options:autostop:True',
        'parameters_numerical:delt:-1'                             : 'time_step:delt:0.03',
        'parameters_numerical:delt_option:check_restart'           : 'time_step:delt_option:check_restart',
        'parameters_numerical:delt_adjust:2.0'                     : 'knobs:delt_adjust:DEPRECATED',
        'parameters_numerical:delt_max:-1'                         : 'time_step:delt_max:-1',
        'parameters_numerical:delt_min:1.e-10'                     : 'time_step:delt_min:1e-10',
        'parameters_numerical:cfl_cushion:0.5'                     : 'time_step:cfl_cushion_upper:0.5', 
        'parameters_numerical:cfl_cushion_upper:0.5'               : 'time_step:cfl_cushion_upper:0.5', 
        'parameters_numerical:cfl_cushion_middle:0.25'             : 'time_step:cfl_cushion_middle:0.25',
        'parameters_numerical:cfl_cushion_lower:0.00001'           : 'time_step:cfl_cushion_lower:1e-05',
        'parameters_numerical:fphi:1.0'                            : 'scale_gyrokinetic_terms:fphi:1.0',
        'parameters_numerical:fields_kxkyz:False'                  : 'parallelisation:fields_kxkyz:False',
        'parameters_numerical:lu_option:none'                      : 'parallelisation:lu_option:default',
        'parameters_numerical:mat_gen:True'                        : 'parallelisation:mat_gen:True',
        'parameters_numerical:mat_read:False'                      : 'parallelisation:mat_read:False',
        'parameters_numerical:stream_implicit:True'                : 'numerical_algorithms:stream_implicit:True',
        'parameters_numerical:mirror_implicit:True'                : 'numerical_algorithms:mirror_implicit:True',
        'parameters_numerical:drifts_implicit:False'               : 'numerical_algorithms:drifts_implicit:False',
        'parameters_numerical:driftkinetic_implicit:False'         : 'numerical_algorithms:driftkinetic_implicit:False',
        'parameters_numerical:maxwellian_inside_zed_derivative:False' : 'numerical_algorithms:maxwellian_inside_zed_derivative:False',
        'parameters_numerical:mirror_semi_lagrange:True'           : 'numerical_algorithms:mirror_semi_lagrange:True',
        'parameters_numerical:mirror_linear_interp:False'          : 'numerical_algorithms:mirror_linear_interp:False',
        'parameters_numerical:stream_matrix_inversion:False'       : 'numerical_algorithms:stream_matrix_inversion:False',
        'parameters_numerical:zed_upwind:0.02'                     : 'numerical_upwinding_for_derivatives:zed_upwind:0.02',
        'parameters_numerical:vpa_upwind:0.02'                     : 'numerical_upwinding_for_derivatives:vpa_upwind:0.02',
        'parameters_numerical:time_upwind:0.02'                    : 'numerical_upwinding_for_derivatives:time_upwind:0.02',
        'parameters_numerical:rng_seed:-1'                         : 'initialise_distribution_noise:rng_seed:-1',
        'parameters_numerical:ky_solve_radial:0'                   : 'multibox_parameters:ky_solve_radial:0',
        'parameters_numerical:ky_solve_real:False'                 : 'multibox_parameters:ky_solve_real:False',
        }
        
    # Replace variables
    if upgrade:
        input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    #===============================================================================
    #                                   Dissipation                                  
    #===============================================================================
    
    renamed_variables = { 
        #------------------- dissipation --> dissipation ------------------
        'dissipation:include_collisions:False'      : 'dissipation_and_collisions_options:include_collisions:False',
        'dissipation:hyper_dissipation:False'       : 'dissipation_and_collisions_options:hyper_dissipation:False',
        'dissipation:collision_model:dougherty'     : 'dissipation_and_collisions_options:collision_model:dougherty',
        'dissipation:collisions_implicit:True'      : 'dissipation_and_collisions_options:collisions_implicit:True',
        #------------------- dissipation --> hyper_dissipation ------------------
        'hyper:d_hyper:0.05'                        : 'hyper_dissipation:d_hyper:0.05',
        'dissipation:d_hyper:0.05'                  : 'hyper_dissipation:d_hyper:0.05',
        'hyper_dissipation:d_zed:DOESNT EXIST YET'  : 'hyper_dissipation:d_zed:0.05',
        'hyper_dissipation:d_vpa:DOESNT EXIST YET'  : 'hyper_dissipation:d_vpa:0.05',
        'hyper_dissipation:hyp_zed:DOESNT EXIST YET': 'hyper_dissipation:hyp_zed:False',
        'hyper_dissipation:hyp_vpa:DOESNT EXIST YET': 'hyper_dissipation:hyp_vpa:False',
        'hyper_dissipation:use_physical_ksqr:DOESNT EXIST YET' : 'hyper_dissipation:use_physical_ksqr:True',
        'hyper_dissipation:scale_to_outboard:DOESNT EXIST YET' : 'hyper_dissipation:scale_to_outboard:False', 
        #------------------- dissipation --> collisions_dougherty ------------------
        'dissipation:momentum_conservation:True'    : 'collisions_dougherty:momentum_conservation:True',
        'dissipation:energy_conservation:True'      : 'collisions_dougherty:energy_conservation:True',
        #------------------- dissipation --> collisions_fokker_planck ------------------
        'dissipation:testpart:True'                 : 'collisions_fokker_planck:testpart:True',
        'dissipation:fieldpart:False'               : 'collisions_fokker_planck:fieldpart:False',
        'dissipation:intraspec:True'                : 'collisions_fokker_planck:intraspec:True',
        'dissipation:interspec:True'                : 'collisions_fokker_planck:interspec:True',
        'dissipation:iiknob:1.0'                    : 'collisions_fokker_planck:iiknob:1.0',
        'dissipation:ieknob:1.0'                    : 'collisions_fokker_planck:ieknob:1.0',
        'dissipation:eeknob:1.0'                    : 'collisions_fokker_planck:eeknob:1.0',
        'dissipation:eiknob:1.0'                    : 'collisions_fokker_planck:eiknob:1.0',
        'dissipation:eiediffknob:1.0'               : 'collisions_fokker_planck:eiediffknob:1.0',
        'dissipation:deflknob:1.0'                  : 'collisions_fokker_planck:deflknob:1.0',
        'dissipation:eimassr_approx:False'          : 'collisions_fokker_planck:eimassr_approx:False',
        'dissipation:advfield_coll:True'            : 'collisions_fokker_planck:advfield_coll:True',
        'dissipation:density_conservation:False'    : 'collisions_fokker_planck:density_conservation:False',
        'dissipation:density_conservation_tp:False' : 'collisions_fokker_planck:density_conservation_tp:False',
        'dissipation:exact_conservation:False'      : 'collisions_fokker_planck:exact_conservation:False',
        'dissipation:spitzer_problem:False'         : 'collisions_fokker_planck:spitzer_problem:False',
        'dissipation:cfac:1'                        : 'collisions_fokker_planck:cfac:1',
        'dissipation:cfac2:1'                       : 'collisions_fokker_planck:cfac2:1',
        'dissipation:nuxfac:1'                      : 'collisions_fokker_planck:nuxfac:1',
        'dissipation:jmax:1'                        : 'collisions_fokker_planck:jmax:1',
        'dissipation:lmax:1'                        : 'collisions_fokker_planck:lmax:1',
        'dissipation:i1fac:1'                       : 'collisions_fokker_planck:i1fac:1',
        'dissipation:i2fac:0'                       : 'collisions_fokker_planck:i2fac:0',
        'dissipation:no_j1l1:True'                  : 'collisions_fokker_planck:no_j1l1:True',
        'dissipation:no_j1l2:False'                 : 'collisions_fokker_planck:no_j1l2:False',
        'dissipation:no_j0l2:False'                 : 'collisions_fokker_planck:no_j0l2:False',
        'dissipation:nvel_local:512'                : 'collisions_fokker_planck:nvel_local:512',
        'dissipation:density_conservation_field:False' : 'collisions_fokker_planck:density_conservation_field:False',
        'collisions_fokker_planck:eideflknob:DOESNT EXIST YET' : 'collisions_fokker_planck:eideflknob:1.0',
        'collisions_fokker_planck:exact_conservation_tp:DOESNT EXIST YET' : 'collisions_fokker_planck:exact_conservation_tp:False',
        }
    
    # Replace variables
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    # The next variable is used in two separate namelists, so don't delete it on the first run
    renamed_variables = { 
        'dissipation:vpa_operator:True' : 'collisions_fokker_planck:vpa_operator:True',
        'dissipation:mu_operator:True'  : 'collisions_fokker_planck:mu_operator:True',}
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade, delete_variable=False)
    renamed_variables = { 
        'dissipation:vpa_operator:True' : 'collisions_dougherty:vpa_operator:True',
        'dissipation:mu_operator:True'  : 'collisions_dougherty:mu_operator:True',}
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    # Deal with stella v0.8
    renamed_variables = { 
        #------------------- dissipation --> dissipation ------------------
        'dissipation:include_collisions:False'      : 'dissipation_and_collisions_options:include_collisions:False',
        'dissipation:hyper_dissipation:False'       : 'dissipation_and_collisions_options:hyper_dissipation:False',
        'dissipation:collision_model:dougherty'     : 'dissipation_and_collisions_options:collision_model:dougherty',
        'dissipation:collisions_implicit:True'      : 'dissipation_and_collisions_options:collisions_implicit:True',
        #------------------- dissipation --> hyper_dissipation ------------------
        'dissipation:d_hyper:0.05'                  : 'hyper_dissipation:d_hyper:0.05',
        'hyper_dissipation:d_zed:DOESNT EXIST YET'  : 'hyper_dissipation:d_zed:0.05',
        'hyper_dissipation:d_vpa:DOESNT EXIST YET'  : 'hyper_dissipation:d_vpa:0.05',
        'hyper_dissipation:hyp_zed:DOESNT EXIST YET': 'hyper_dissipation:hyp_zed:False',
        'hyper_dissipation:hyp_vpa:DOESNT EXIST YET': 'hyper_dissipation:hyp_vpa:False',
        'hyper_dissipation:use_physical_ksqr:DOESNT EXIST YET' : 'hyper_dissipation:use_physical_ksqr:True',
        'hyper_dissipation:scale_to_outboard:DOESNT EXIST YET' : 'hyper_dissipation:scale_to_outboard:False', 
        #------------------- dissipation --> collisions_dougherty ------------------
        'collisions_dougherty:momentum_conservation:True'  : 'collisions_dougherty:momentum_conservation:True',
        'collisions_dougherty:energy_conservation:True'    : 'collisions_dougherty:energy_conservation:True',
        #------------------- dissipation --> collisions_fokker_planck ------------------
        'collisions_fp:testpart:True'                 : 'collisions_fokker_planck:testpart:True',
        'collisions_fp:fieldpart:False'               : 'collisions_fokker_planck:fieldpart:False',
        'collisions_fp:intraspec:True'                : 'collisions_fokker_planck:intraspec:True',
        'collisions_fp:interspec:True'                : 'collisions_fokker_planck:interspec:True',
        'collisions_fp:iiknob:1.0'                    : 'collisions_fokker_planck:iiknob:1.0',
        'collisions_fp:ieknob:1.0'                    : 'collisions_fokker_planck:ieknob:1.0',
        'collisions_fp:eeknob:1.0'                    : 'collisions_fokker_planck:eeknob:1.0',
        'collisions_fp:eiknob:1.0'                    : 'collisions_fokker_planck:eiknob:1.0',
        'collisions_fp:eiediffknob:1.0'               : 'collisions_fokker_planck:eiediffknob:1.0',
        'collisions_fp:deflknob:1.0'                  : 'collisions_fokker_planck:deflknob:1.0',
        'collisions_fp:eimassr_approx:False'          : 'collisions_fokker_planck:eimassr_approx:False',
        'collisions_fp:advfield_coll:True'            : 'collisions_fokker_planck:advfield_coll:True',
        'collisions_fp:density_conservation:False'    : 'collisions_fokker_planck:density_conservation:False',
        'collisions_fp:density_conservation_tp:False' : 'collisions_fokker_planck:density_conservation_tp:False',
        'collisions_fp:exact_conservation:False'      : 'collisions_fokker_planck:exact_conservation:False',
        'collisions_fp:spitzer_problem:False'         : 'collisions_fokker_planck:spitzer_problem:False',
        'collisions_fp:cfac:1'                        : 'collisions_fokker_planck:cfac:1',
        'collisions_fp:cfac2:1'                       : 'collisions_fokker_planck:cfac2:1',
        'collisions_fp:nuxfac:1'                      : 'collisions_fokker_planck:nuxfac:1',
        'collisions_fp:jmax:1'                        : 'collisions_fokker_planck:jmax:1',
        'collisions_fp:lmax:1'                        : 'collisions_fokker_planck:lmax:1',
        'collisions_fp:i1fac:1'                       : 'collisions_fokker_planck:i1fac:1',
        'collisions_fp:i2fac:0'                       : 'collisions_fokker_planck:i2fac:0',
        'collisions_fp:no_j1l1:True'                  : 'collisions_fokker_planck:no_j1l1:True',
        'collisions_fp:no_j1l2:False'                 : 'collisions_fokker_planck:no_j1l2:False',
        'collisions_fp:no_j0l2:False'                 : 'collisions_fokker_planck:no_j0l2:False',
        'collisions_fp:nvel_local:512'                : 'collisions_fokker_planck:nvel_local:512',
        'collisions_fp:vpa_operator:True'             : 'collisions_fokker_planck:vpa_operator:True',
        'collisions_fp:mu_operator:True'              : 'collisions_fokker_planck:mu_operator:True',
        'collisions_fp:density_conservation_field:False' : 'collisions_fokker_planck:density_conservation_field:False',
        'collisions_fp:eideflknob:1.0' : 'collisions_fokker_planck:eideflknob:1.0',
        'collisions_fp:exact_conservation_tp:False' : 'collisions_fokker_planck:exact_conservation_tp:False',
        }
        
    # Replace variables
    if upgrade:
        input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    #===============================================================================
    #                                     Others                                  
    #===============================================================================
    
    renamed_variables = { 
        #------------------- sfincs_input ------------------
        'sfincs_input:read_sfincs_output_from_file:False' : 'sfincs_input:read_sfincs_output_from_file:False',
        'sfincs_input:nproc_sfincs:1'               : 'sfincs_input:nproc_sfincs:1',
        'sfincs_input:calculate_radial_electric_field:True' : 'sfincs_input:calculate_radial_electric_field:True',
        'sfincs_input:includexdotterm:True'         : 'sfincs_input:includexdotterm:True',
        'sfincs_input:includeelectricfieldterminxidot:True' : 'sfincs_input:includeelectricfieldterminxidot:True',
        'sfincs_input:magneticdriftscheme:0'        : 'sfincs_input:magneticdriftscheme:0',
        'sfincs_input:includephi1:True'             : 'sfincs_input:includephi1:True',
        'sfincs_input:includephi1inkineticequation:False' : 'sfincs_input:includephi1inkineticequation:False',
        'sfincs_input:geometryscheme:1'             : 'sfincs_input:geometryscheme:1',
        'sfincs_input:vmecradialoption:0'           : 'sfincs_input:vmecradialoption:0',
        'sfincs_input:equilibriumfile:wout_161s1.nc' : 'sfincs_input:equilibriumfile:wout_161s1.nc',
        'sfincs_input:coordinatesystem:3'           : 'sfincs_input:coordinatesystem:3',
        'sfincs_input:inputradialcoordinate:3'      : 'sfincs_input:inputradialcoordinate:3',
        'sfincs_input:inputradialcoordinateforgradients:3' : 'sfincs_input:inputradialcoordinateforgradients:3',
        'sfincs_input:ahat:1.0'                     : 'sfincs_input:ahat:1.0',
        'sfincs_input:delta:-1.0'                   : 'sfincs_input:delta:-1.0',
        'sfincs_input:nu_n:-1.0'                    : 'sfincs_input:nu_n:-1.0',
        'sfincs_input:dphihatdrn:-9999.9'           : 'sfincs_input:dphihatdrn:-9999.9',
        'sfincs_input:er_window:0.3'                : 'sfincs_input:er_window:0.3',
        'sfincs_input:nxi:48'                       : 'sfincs_input:nxi:48',
        'sfincs_input:nx:12'                        : 'sfincs_input:nx:12',
        'sfincs_input:ntheta:65'                    : 'sfincs_input:ntheta:65',
        'sfincs_input:nzeta:1'                      : 'sfincs_input:nzeta:1',
        #------------------- layouts_knobs --> parallelisation ------------------
        'layouts_knobs:xyzs_layout:xyzs'            : 'parallelisation:xyzs_layout:xyzs',
        'layouts_knobs:vms_layout:vms'              : 'parallelisation:vms_layout:vms',
        #------------------- sources ------------------
        'sources:include_krook_operator:False'      : 'sources:include_krook_operator:False',
        'sources:exclude_boundary_regions:False'    : 'sources:exclude_boundary_regions:False',
        'sources:exclude_boundary_regions_qn:False' : 'sources:exclude_boundary_regions_qn:False',
        'sources:remove_zero_projection:False'      : 'sources:remove_zero_projection:False',
        'sources:nu_krook:0.05'                     : 'sources:nu_krook:0.05',
        'sources:tcorr_source:0.02'                 : 'sources:tcorr_source:0.02',
        'sources:tcorr_source_qn:0.0'               : 'sources:tcorr_source_qn:0.0',
        'sources:ikxmax_source:1'                   : 'sources:ikxmax_source:1',
        'sources:krook_odd:True'                    : 'sources:krook_odd:True',
        'sources:from_zero:True'                    : 'sources:from_zero:True',
        'sources:conserve_momentum:False'           : 'sources:conserve_momentum:False',
        'sources:conserve_density:False'            : 'sources:conserve_density:False',
        'sources:source_option:DOESNT EXIST YET'    : 'sources:source_option:none',
        #------------------- neoclassical_input ------------------
        'neoclassical_input:include_neoclassical_terms:False' : 'neoclassical_input:include_neoclassical_terms:False',
        'neoclassical_input:nradii:5'               : 'neoclassical_input:nradii:5',
        'neoclassical_input:drho:0.01'              : 'neoclassical_input:drho:0.01',
        'neoclassical_input:neo_option:sfincs'      : 'neoclassical_input:neo_option:sfincs',
        #------------------- multibox_parameters ------------------
        'multibox_parameters:boundary_size:4'       : 'multibox_parameters:boundary_size:4',
        'multibox_parameters:krook_size:0'          : 'multibox_parameters:krook_size:0',
        'multibox_parameters:phi_bound:0'           : 'multibox_parameters:phi_bound:0',
        'multibox_parameters:phi_pow:0'             : 'multibox_parameters:phi_pow:0',
        'multibox_parameters:krook_exponent:0.0'    : 'multibox_parameters:krook_exponent:0.0',
        'multibox_parameters:krook_efold:3.0'       : 'multibox_parameters:krook_efold:3.0',
        'multibox_parameters:nu_krook_mb:0.0'       : 'multibox_parameters:nu_krook_mb:0.0',
        'multibox_parameters:mb_debug_step:1000'    : 'multibox_parameters:mb_debug_step:1000',
        'multibox_parameters:smooth_zfs:False'      : 'multibox_parameters:smooth_zf:False',
        'multibox_parameters:comm_at_init:False'    : 'multibox_parameters:comm_at_init:False',
        'multibox_parameters:rk_step:False'         : 'multibox_parameters:rk_step:False',
        'multibox_parameters:zf_option:default'     : 'multibox_parameters:zf_option:default',
        'multibox_parameters:krook_option:linear'   : 'multibox_parameters:krook_option:linear',
        'multibox_parameters:lr_debug_option:default': 'multibox_parameters:lr_debug_option:default',
        'multibox_parameters:use_dirichlet_bc:False' : 'multibox_parameters:use_dirichlet_bc:False',
        #------------------- euterpe_parameters ------------------
        'euterpe_parameters:nradii:1000'            : 'euterpe_parameters:nradii:1000',
        'euterpe_parameters:data_file:euterpe.dat'  : 'euterpe_parameters:data_file:euterpe.dat',
        #------------------- debug_flags ------------------
        'debug_flags:print_extra_info_to_terminal:DOESNT EXIST YET' : 'debug_flags:print_extra_info_to_terminal:True',
        'debug_flags:debug_all:DOESNT EXIST YET'    : 'debug_flags:debug_all:False',
        'debug_flags:stella_debug:DOESNT EXIST YET' : 'debug_flags:stella_debug:False',
        'debug_flags:ffs_solve_debug:DOESNT EXIST YET' : 'debug_flags:ffs_solve_debug:False',
        'debug_flags:fields_all_debug:DOESNT EXIST YET' : 'debug_flags:fields_all_debug:False',
        'debug_flags:fields_debug:DOESNT EXIST YET' : 'debug_flags:fields_debug:False',
        'debug_flags:fields_fluxtube_debug:DOESNT EXIST YET' : 'debug_flags:fields_fluxtube_debug:False',
        'debug_flags:fields_electromagnetic_debug:DOESNT EXIST YET' : 'debug_flags:fields_electromagnetic_debug:False',
        'debug_flags:fields_ffs_debug:DOESNT EXIST YET' : 'debug_flags:fields_ffs_debug:False',
        'debug_flags:implicit_solve_debug:DOESNT EXIST YET' : 'debug_flags:implicit_solve_debug:False',
        'debug_flags:parallel_streaming_debug:DOESNT EXIST YET' : 'debug_flags:parallel_streaming_debug:False',
        'debug_flags:mirror_terms_debug:DOESNT EXIST YET' : 'debug_flags:mirror_terms_debug:False',
        'debug_flags:neoclassical_terms_debug:DOESNT EXIST YET' : 'debug_flags:neoclassical_terms_debug:False',
        'debug_flags:response_matrix_debug:DOESNT EXIST YET' : 'debug_flags:response_matrix_debug:False',
        'debug_flags:time_advance_debug:DOESNT EXIST YET' : 'debug_flags:time_advance_debug:False',
        'debug_flags:extended_grid_debug:DOESNT EXIST YET' : 'debug_flags:extended_grid_debug:False',
        'debug_flags:diagnostics_all_debug:DOESNT EXIST YET' : 'debug_flags:diagnostics_all_debug:False',
        'debug_flags:diagnostics_parameters:DOESNT EXIST YET' : 'debug_flags:diagnostics_parameters:False',
        'debug_flags:diagnostics_fluxes_fluxtube_debug:DOESNT EXIST YET' : 'debug_flags:diagnostics_fluxes_fluxtube_debug:False',
        'debug_flags:diagnostics_omega_debug:DOESNT EXIST YET' : 'debug_flags:diagnostics_omega_debug:False',
        'debug_flags:diagnostics_debug:DOESNT EXIST YET' : 'debug_flags:diagnostics_debug:False',
        'debug_flags:dist_fn_debug:DOESNT EXIST YET' : 'debug_flags:dist_fn_debug:False',
        'debug_flags:gyro_averages_debug:DOESNT EXIST YET' : 'debug_flags:gyro_averages_debug:False',
        'debug_flags:fluxes_debug:DOESNT EXIST YET' : 'debug_flags:fluxes_debug:False',
        'debug_flags:geometry_debug:DOESNT EXIST YET' : 'debug_flags:geometry_debug:False',
        'debug_flags:const_alpha_geo:DOESNT EXIST YET' : 'debug_flags:const_alpha_geo:False', 
        }
    
    # Replace variables
    input_parameters = replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade)
    
    #===============================================================================
    #                                   Extra rules                                  
    #===============================================================================
    
    # For old stella versions, always turn off apar and radial_moments
    if downgrade:
        if 'knobs' not in input_parameters.keys(): input_parameters['knobs'] = {}
        input_parameters['knobs']['fapar'] = 0
        input_parameters['knobs']['fbpar'] = 0
        if 'stella_diagnostics_knobs' not in input_parameters.keys(): input_parameters['stella_diagnostics_knobs'] = {}
        input_parameters['stella_diagnostics_knobs']['write_radial_moments'] = False
    
    #===============================================================================
    #                                   Write file                                  
    #===============================================================================
    
    # Write the namelist
    path_new_input_file = str(path_input_file).replace('.in', '_updated.in')
    if downgrade: path_new_input_file = str(path_input_file).replace('.in', '_downgraded.in')
    write_dictionaryToNamelist(path_new_input_file, input_parameters, downgrade) 
    return input_parameters

#--------------------------------------------------------------
def replace_variables(input_parameters, renamed_variables, add_default_variables, downgrade=False, delete_variable=True):
    
    # If we want to downgrade an input file from version 1.0 to version 0.5, invert the dictionary
    if downgrade: renamed_variables = {v: k for k, v in renamed_variables.items()}
    
    # Find the name change namelist_old:variable_old --> namelist_new:variable_new
    for namelist_variable_value_old in renamed_variables.keys():
        namelist_variable_value_new = renamed_variables[namelist_variable_value_old]
        namelist_new = namelist_variable_value_new.split(':')[0]
        variable_new = namelist_variable_value_new.split(':')[1]
        value_new = namelist_variable_value_new.split(':')[2]
        namelist_old = namelist_variable_value_old.split(':')[0]
        variable_old = namelist_variable_value_old.split(':')[1]
        value_old = namelist_variable_value_old.split(':')[2]
        
        # If the variable name hasn't changed, continue
        if namelist_new==namelist_old and variable_new==variable_old: 
            input_parameters = add_default_variable(input_parameters, namelist_new, variable_new, value_new, namelist_old, variable_old, add_default_variables)
            continue
        
        # If the variable is present, and the name has changed, update it
        if namelist_old in input_parameters.keys():
            if variable_old in input_parameters[namelist_old].keys(): 
                if namelist_new not in input_parameters.keys(): input_parameters[namelist_new] = {}
                input_parameters[namelist_new][variable_new] = input_parameters[namelist_old][variable_old]
                if delete_variable: del input_parameters[namelist_old][variable_old]
                continue
                
        # If the variable is not present, add the default value
        input_parameters = add_default_variable(input_parameters, namelist_new, variable_new, value_new, namelist_old, variable_old, add_default_variables)
        
    return input_parameters

#--------------------------------------------------------------
def add_default_variable(input_parameters, namelist_new, variable_new, value_new, namelist_old, variable_old, add_default_variables):
    
    # Only add default variables if we want to
    if add_default_variables==False: 
        return input_parameters
    
    # Parse the value to the correct type
    value_numeric = value_new.replace('.','').replace('-','').replace('e','').replace(' ','')
    if value_new=='False': value_new = False
    elif value_new=='True': value_new = True
    elif value_numeric.isnumeric():
        if '.' in value_new: value_new = float(value_new)
        elif 'e' in value_new: value_new = float(value_new)
        else: value_new = int(value_new)
        
    # Add the default value
    if namelist_new not in input_parameters.keys(): 
        input_parameters[namelist_new] = {}
    if namelist_old not in input_parameters.keys(): 
        input_parameters[namelist_new][variable_new] = value_new
        return input_parameters
    if variable_old not in input_parameters[namelist_old].keys(): 
        input_parameters[namelist_old][variable_new] = value_new
    return input_parameters

#--------------------------------------------------------------
def write_dictionaryToNamelist(path, dictionary, downgrade, indent="  ", sort_knobs=False):
    
    # Order of namelists in version 1.0
    ordered_namelists = ['&Geometry', 'geometry_options', 'geometry_miller', 'geometry_vmec', 'geometry_zpinch', 'geometry_from_txt',
        '&Physics', 'gyrokinetic_terms', 'scale_gyrokinetic_terms', 'physics_inputs', 'flux_annulus', 'electromagnetic',  
        '&Diagnostics', 'diagnostics', 'diagnostics_moments', 'diagnostics_omega', 'diagnostics_distribution', 'diagnostics_fluxes', 'diagnostics_potential',
        '&Initialise Distribution', 'initialise_distribution', 'initialise_distribution_maxwellian', 'initialise_distribution_noise', 'initialise_distribution_kpar', 'initialise_distribution_rh', 'restart_options',
        '&Ion, electron and impurity species', 'species_options', 'adiabatic_electron_response', 'adiabatic_ion_response', 'species_parameters_1', 'species_parameters_2', 'species_parameters_3', 'species_parameters_4', 'species_parameters_5', 'euterpe_parameters', 
        '&Discretized (kx,ky,z,mu,vpa) grid', 'kxky_grid_option', 'kxky_grid_range', 'kxky_grid_box', 'z_grid', 'z_boundary_condition', 'velocity_grids',
        '&Numerics', 'time_trace_options', 'time_step', 'numerical_algorithms', 'numerical_upwinding_for_derivatives', 
        '&Dissipation', 'dissipation_and_collisions_options', 'hyper_dissipation', 'collisions_dougherty', 'collisions_fokker_planck',
        '&Neoclassics', 'neoclassical_input', 'sfincs_input', 
        '&Radial variation', 'multibox_parameters', 'sources', 'flow_shear', 
        '&Parallelisation', 'parallelisation', 
        '&Debug', 'debug_flags']
    if downgrade:
        ordered_namelists = ['&Geometry', 'geo_knobs', 'vmec_parameters', 'millergeo_parameters',
            '&Physics', 'parameters', 'physics_flags',
            '&Diagnostics', 'stella_diagnostics_knobs',
            '&Initialise Distribution', 'init_g_knobs',
            '&Ion, electron and impurity species', 'species_knobs', 'species_parameters_1', 'species_parameters_2', 
            'species_parameters_3', 'species_parameters_4', 'species_parameters_5', 'euterpe_parameters', 
            '&Discretized (kx,ky,z,mu,vpa) grid', 'kt_grids_knobs', 'kt_grids_range_parameters', 'kt_grids_box_parameters', 'zgrid_parameters', 'vpamu_grids_parameters',
            '&Numerics', 'knobs', 'time_advance_knobs', 
            '&Dissipation', 'dissipation',
            '&Neoclassics', 'neoclassical_input', 'sfincs_input', 
            '&Radial variation', 'multibox_parameters', 'sources', 
            '&Parallelisation', 'layouts_knobs']
        
    
    # Write dictionary as a namelist
    with open(path, 'w') as f:
        for namelist in ordered_namelists:
            
            # Print headers
            if '&' in namelist:
                f.write('\n!'+''.center(80,'=')+'\n')
                f.write('!'+f' {namelist[1:]} '.center(80,'=')+'\n')
                f.write('!'+''.center(80,'=')+'\n\n')
                continue 
            
            # Check
            if namelist not in dictionary.keys():
                continue
                
            # Print namelists
            if namelist in dictionary.keys():
                
                # Sort the variable names
                variables = list(dictionary[namelist].keys())
                variables.sort()
                
                # Print namelists name
                if len(variables)==0: 
                    print(f'There are no variables in {namelist}')
                    continue
                f.write(f'&{namelist}\n')
                
                # Print variables
                for variable in variables:
                        value = dictionary[namelist][variable]
                        if type(value)==bool:
                            if value==False: value = '.false.'
                            if value==True: value = '.true.'
                        elif type(value)==str:
                            value = "'"+value+"'"
                        f.write(f'  {variable} = {value}\n')
                f.write('/\n\n')

    return 

#--------------------------------------------------------------
def get_differenceInDictionaries(dict1, dict2):
    ''' Returns the difference in keys between two dictionaries. '''

    # Initiate the different dictionary
    dict_difference = {}
    
    # Get the difference of a two layered dictionary
    for namelist_name, namelist in dict1.items():
            
        # Check whether the namelist is present
        if namelist_name not in dict2.keys():
            print(f'WARNING: the {namelist_name} namelist is missing in the new dictionary.')
            dict_difference[namelist_name] = {'ALL OF THEM' : 'ARE MISSING'}
            continue
            
        # Check the variables within the namelist
        for variable, value in namelist.items(): 
            
            # Check whether the variable is present
            if variable not in dict2[namelist_name].keys():
                print(f'WARNING: the {namelist_name}:{variable} variable is missing in the new dictionary.')
                if namelist_name not in dict_difference:
                    dict_difference[namelist_name] = {variable : 'VARIABLE IS MISSING'}
                else:
                    dict_difference[namelist_name][variable] = 'VARIABLE IS MISSING' 
                continue
            
            # Compare two numbers
            if type(value)==float or type(value)==int or isinstance(value, np.floating):
                if value != np.nan and not np.isnan(value) and value != None and dict2[namelist_name][variable] != None:  
                    if dict2[namelist_name][variable] != value:  
                        if namelist_name not in dict_difference:
                            dict_difference[namelist_name] = {variable : f'{value}, {dict2[namelist_name][variable]}'}
                        else:
                            dict_difference[namelist_name][variable] = f'{value}, {dict2[namelist_name][variable]}'
                            
            # Compare two strings
            else:
                if dict2[namelist_name][variable] != value:  
                    if namelist_name not in dict_difference:
                        dict_difference[namelist_name] = {variable : f'{value}, {dict2[namelist_name][variable]}'}
                    else:
                        dict_difference[namelist_name][variable] = f'{value}, {dict2[namelist_name][variable]}'
                
    # Count the number of differences
    numberOfDifferences = 0
    for namelist in dict_difference.keys():
        numberOfDifferences += len(dict_difference[namelist])  
    return dict_difference, numberOfDifferences



#p------------------------------------------------------
def tool_to_write_this_script():

    # Read the default parameters
    path = pathlib.Path(os.path.realpath(__file__))
    path = pathlib.Path(str(path).split('AUTOMATIC_TESTS')[0]) 
    input_file_v10 = path / 'AUTOMATIC_TESTS/convert_input_files/input_stella_v1.0.in' 
    input_file_v5 = path / 'AUTOMATIC_TESTS/convert_input_files/input_stella_v0.5.in' 
    
    # Read default input parameters  
    input_parameters_v10 = f90nml.read(input_file_v10)
    input_parameters_v10 = json.loads(json.dumps(input_parameters_v10))  
    input_parameters_v5 = f90nml.read(input_file_v5)
    input_parameters_v5 = json.loads(json.dumps(input_parameters_v5))  
    
    # Read the text files
    with open(input_file_v10) as f:
        text_input_file_v10 = f.read()
    with open(input_file_v5) as f:
        text_input_file_v5 = f.read()
        
    # Get the namelists
    namelists_stella_v05 = []
    namelists_stella_v10 = []
    for line in text_input_file_v5.split('\n'):
        if '&' in line: namelists_stella_v05.append(line.replace(' ','').replace('&',''))
    for line in text_input_file_v10.split('\n'):
        if '&' in line: namelists_stella_v10.append(line.replace(' ','').replace('&',''))
    namelists_stella_v05 = list(set(namelists_stella_v05))
    namelists_stella_v10 = list(set(namelists_stella_v10))
    if False:
        for i in range(len(namelists_stella_v05)):
            print(namelists_stella_v05[i])
        print('---------')
        for i in range(len(namelists_stella_v10)):
            print(namelists_stella_v10[i])
                
    # Namelists with their old (v0.5) and new (v1.0) names
    geometry_namelists_v05 = ['geo_knobs', 'vmec_parameters', 'millergeo_parameters']
    dissipation_namelists_v05 = ['dissipation']
    physics_namelists_v05 = ['physics_flags', 'time_advance_knobs', 'parameters',]
    numerical_namelists_v05 = ['knobs']
    fields_namelists_v05 = ['init_g_knobs']
    diagnostics_namelists_v05 = ['stella_diagnostics_knobs']
    species_namelists_v05 = [ 'species_knobs',  'species_parameters_1']
    grids_namelists_v05 = ['kt_grids_knobs', 'kt_grids_range_parameters', 'kt_grids_box_parameters', 'zgrid_parameters', 'vpamu_grids_parameters', ]
    remaining_namelists_v05 = ['neoclassical_input',  'sfincs_input', 'layouts_knobs', 'euterpe_parameters', 'multibox_parameters', 'sources']
        
    # Write dicts
    if True:
        print('\n\n')
        input_parameters = input_parameters_v10
        for namelist in input_parameters:
            if 'diagnostics' in namelist:
                print(f'------------------- {namelist} ------------------') 
                for variable in input_parameters[namelist].keys():
                    if False: print(f"'{namelist}:{variable}:{input_parameters[namelist][variable]}' : '{namelist}:{variable}:{input_parameters[namelist][variable]}',")
                    print(f"input_parameters['{namelist}']['{variable}'] = {input_parameters[namelist][variable]}")
    

#===============================================================================
#                                DEBUG SCRIPT                                  #
#===============================================================================
 
if __name__ == '__main__' and False:
    
    # Use the following tool to write this script
    #tool_to_write_this_script()
    
    # Upgrade input file from stella version 0.5 to stella version 1.0
    path = pathlib.Path(os.path.realpath(__file__))
    path_input_file = path.parent / 'input_stella_v0.5.in'
    input_parameters = update_inputFile(path_input_file, add_default_variables=True)
    
    # Upgrade input file from stella version 0.8 to stella version 1.0
    path = pathlib.Path(os.path.realpath(__file__))
    path_input_file = path.parent / 'input_stella_v0.8.in'
    input_parameters = update_inputFile(path_input_file, add_default_variables=True)
    
    # Downgrade input file from stella version 1.0 to stella version 0.5
    path = pathlib.Path(os.path.realpath(__file__))
    path_input_file = path.parent / 'input_stella_v1.0.in'
    input_parameters = update_inputFile(path_input_file, add_default_variables=True, downgrade=True)
    
    # Check whether we can correctly upgrade a file from 0.5 to 1.0
    if False:
        print('\n Test updating from stella v0.5 to stella v1.0\n')
        input_parameters_correct = f90nml.read(path.parent/'input_stella_v1.0.in')
        input_parameters_correct = json.loads(json.dumps(input_parameters_correct))  
        input_parameters_testscript = f90nml.read(path.parent/'input_stella_v0.5_updated.in' )
        input_parameters_testscript = json.loads(json.dumps(input_parameters_testscript))  
        dict_difference, numberOfDifferences = get_differenceInDictionaries(input_parameters_correct, input_parameters_testscript)
        print('\n\n------------ Print differences ------------------', numberOfDifferences)
        for namelist in dict_difference.keys():
            print('\n==================================')
            print(f'{namelist}'.center(34))
            print('==================================')
            for variable, value in dict_difference[namelist].items():
                print(f'{variable}: {value}') 
    
    # Check whether we can correctly upgrade a file from 0.8 to 1.0
    if True:
        print('\n Test updating from stella v0.8 to stella v1.0\n')
        input_parameters_correct = f90nml.read(path.parent/'input_stella_v1.0.in')
        input_parameters_correct = json.loads(json.dumps(input_parameters_correct))  
        input_parameters_testscript = f90nml.read(path.parent/'input_stella_v0.8_updated.in' )
        input_parameters_testscript = json.loads(json.dumps(input_parameters_testscript))  
        dict_difference, numberOfDifferences = get_differenceInDictionaries(input_parameters_correct, input_parameters_testscript)
        print('\n\n------------ Print differences ------------------', numberOfDifferences)
        for namelist in dict_difference.keys():
            print('\n==================================')
            print(f'{namelist}'.center(34))
            print('==================================')
            for variable, value in dict_difference[namelist].items():
                print(f'{variable}: {value}') 
    
    # Check whether we can correctly downgrade a file from 1.0 to 0.5
    if False:
        print('\n Test downgrading from stella v1.0 to stella v0.5\n')
        input_parameters_correct = f90nml.read(path.parent/'input_stella_v0.5.in')
        input_parameters_correct = json.loads(json.dumps(input_parameters_correct))  
        input_parameters_testscript = f90nml.read(path.parent/'input_stella_v1.0_downgraded.in' )
        input_parameters_testscript = json.loads(json.dumps(input_parameters_testscript))  
        dict_difference, numberOfDifferences = get_differenceInDictionaries(input_parameters_correct, input_parameters_testscript)
        print('\n\n------------ Print differences ------------------', numberOfDifferences)
        for namelist in dict_difference.keys():
            print('\n==================================')
            print(f'{namelist}'.center(34))
            print('==================================')
            for variable, value in dict_difference[namelist].items():
                print(f'{variable}: {value}') 

    # Some examples
    if False:
        path_input_file = path.parent / 'miller_nonlinear_CBC.in'
        input_parameters = update_inputFile(path_input_file)
        path_input_file = path.parent / 'W7X_old_input.in'
        input_parameters = update_inputFile(path_input_file)
    
    
#===============================================================================
#                             RUN AS MAIN SCRIPT                               #
#===============================================================================

if __name__ == '__main__' and True:

    # Import modules
    import glob
    
    # Get the input files in the current folder
    folder = pathlib.Path(os.getcwd())
    input_files = [pathlib.Path(i) for i in glob.glob(str(folder)+'/*.in')]
    input_files = [i for i in input_files if '_updated.in' not in i.name]
    input_files = [i for i in input_files if '_downgraded.in' not in i.name]
    
    # Upgrade or downgrade the input file
    for input_file in input_files:
        version_05 = False; version_10 = False; version_08 = False
        with open(input_file) as f:
            text = f.read()
            if '&geo_knobs' in text: version_05 = True
            if '&millergeo_parameters' in text: version_05 = True
            if '&vmec_parameters' in text: version_05 = True
            if '&geometry_options' in text: version_10 = True
            if '&geometry_miller' in text: version_10 = True
            if '&geometry_vmec' in text: version_10 = True
            if version_05:
                if '&parameters_numerical' in text: version_08 = True
                if '&collisions_fp' in text: version_08 = True 
                if '&parameters_physics' in text: version_08 = True
        if version_05:
            print(f'\nUpdating {input_file.name} from stella version 0.5 to stella version 1.0')
            update_inputFile(input_file)
        if version_10:
            print(f'\nDowngrading {input_file.name} from stella version 1.0 to stella version 0.5')
            update_inputFile(input_file, downgrade=True)
        print(f'   --> Converted {input_file.name}')
    
