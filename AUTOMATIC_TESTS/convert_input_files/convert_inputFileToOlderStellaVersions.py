#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Hanne Thienpondt

import copy
import os, sys
import pathlib
import f90nml, json

#--------------------------------------------------------------
def convert_inputFileToOlderStellaVersions(path_input_file='', stella_version='0.5'):
    
    # Read the default parameters of the master branch of stella
    path = pathlib.Path(os.path.realpath(__file__))
    path = pathlib.Path(str(path).split('AUTOMATIC_TESTS/numerical_tests')[0]) 
    default_input_file = path / 'AUTOMATIC_TESTS/convert_input_files/input_stella_master.in' 
    
    # Read default input parameters  
    input_parameters = f90nml.read(default_input_file)
    input_parameters = json.loads(json.dumps(input_parameters)) 
    input_parameters_default =  copy.deepcopy(input_parameters)
    
    # Read the real input file
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))  
    
    # First update the number of species
    try: input_parameters["species_knobs"].update(input_parameters_read["species_knobs"])
    except: pass
     
    # Add more default species if nspec>2
    for i in range(2, int(input_parameters["species_knobs"]["nspec"])+1): 
        input_parameters["species_parameters_"+str(i)] = input_parameters["species_parameters_1"].copy() 

    # Now overwrite with the real input file
    for knob in input_parameters_read.keys():  
        if knob in input_parameters.keys():   
            input_parameters[knob].update(input_parameters_read[knob])
            
    # Remember some values which didn't exist in v0.5 but did exist in v0.6 or v0.7
    saved_values = {
      'cfl_cushion_upper'  : input_parameters['parameters_numerical']['cfl_cushion_upper'],
      'cfl_cushion_middle' : input_parameters['parameters_numerical']['cfl_cushion_middle'],
      'cfl_cushion_lower'  : input_parameters['parameters_numerical']['cfl_cushion_lower'],
      'delt_min' : input_parameters['parameters_numerical']['delt_min'],}

    # Convert to older stella version
    input_parameters = convert_stellaMasterToStellav05(input_parameters, input_parameters_default, stella_version)
    
    # Upgrade again
    if stella_version in ['0.6', '0.7']:
        input_parameters = convert_stellaMasterv05ToStellav06(input_parameters, saved_values)
    return input_parameters


#--------------------------------------------------
def convert_stellaMasterv05ToStellav06(input_parameters, saved_values): 

    # Once the input file has been converted to stella v0.5, it's easy to convert it to v0.6
    input_parameters['knobs']['cfl_cushion_upper'] = saved_values['cfl_cushion_upper']
    input_parameters['knobs']['cfl_cushion_middle'] = saved_values['cfl_cushion_middle']
    input_parameters['knobs']['cfl_cushion_lower'] = saved_values['cfl_cushion_lower']
    input_parameters['knobs']['delt_min'] = saved_values['delt_min']
    del input_parameters['knobs']['cfl_cushion']
    del input_parameters['knobs']['delt_adjust']
    return input_parameters

#--------------------------------------------------
def convert_stellaMasterToStellav05(input_parameters, input_parameters_default, stella_version): 

    # Remember some values
    write_all = input_parameters['stella_diagnostics_knobs']['write_all']
        
    # Obtain the changes between the input files
    renamed_variables, new_variables, deprecated_variables, deprecated_namelist, new_namelists = update_stella05ToStellaMaster()
        
    # If new variables are not set to their default, given a warning
    for namelist_variable in new_variables.keys():
        namelist_new = namelist_variable.split(':')[0]
        variable_new = namelist_variable.split(':')[-1]
        if input_parameters[namelist_new][variable_new]!=input_parameters_default[namelist_new][variable_new]:
            if variable_new not in ['print_extra_info_to_terminal', 'write_all']:
                print('WARNING: a variable which only exists in the newest stella version has been set to a non-default value.')
                print('Older stella versions will not be able to take this into account:')
                print(f'       [{namelist_new}][{variable_new}] = {input_parameters[namelist_new][variable_new]} != {input_parameters_default[namelist_new][variable_new]}\n') 
        
        # Remove new variables
        del input_parameters[namelist_new][variable_new]
        
    # Convert the new variables back to the old variables
    for namelist_variable in renamed_variables.keys():
        namelist_old = namelist_variable.split(':')[0]
        variable_old = namelist_variable.split(':')[-1]
        namelist_new = renamed_variables[namelist_variable].split(':')[0]
        variable_new = renamed_variables[namelist_variable].split(':')[-1] 
        if namelist_old not in input_parameters.keys(): input_parameters[namelist_old] = {}
        
        # Exceptions: in v0.6 the dissipation module had already been split up
        if stella_version in ['0.6', '0.7'] and namelist_old=='dissipation':
            pass
        else:
            input_parameters[namelist_old][variable_old] = input_parameters[namelist_new][variable_new]
        
        # Exceptions when it comes to values
        if namelist_old=='knobs':
            if variable_old=='fapar' or variable_old=='fbpar':
                if input_parameters[namelist_old][variable_old]==False: input_parameters[namelist_old][variable_old] = -1
                if input_parameters[namelist_old][variable_old]==True: input_parameters[namelist_old][variable_old] = 1
        
        # Exceptions when old variable has been parsed to 2 different namelists: 
        if namelist_old=='dissipation' and variable_old=='vpa_operator' and stella_version=='0.5':
            del input_parameters['collisions_fp'][variable_new]
            del input_parameters['collisions_dougherty'][variable_new]
        elif namelist_old=='dissipation' and variable_old=='mu_operator' and stella_version=='0.5':
            del input_parameters['collisions_fp'][variable_new]
            del input_parameters['collisions_dougherty'][variable_new]
        else:
            del input_parameters[namelist_new][variable_new]

    # Remove new namelists
    for namelist_new in new_namelists:
        if stella_version in ['0.6', '0.7'] and namelist_new in ['collisions_dougherty',  'collisions_fp']: continue 
        if len(input_parameters[namelist_new].keys())!=0: 
            print('ERROR: after correct conversion, the new namelist should be empty!')
            print(f'Namelist "{namelist_new}":', input_parameters[namelist_new]); sys.exit()
        del input_parameters[namelist_new]

    # Take care of synonyms
    input_parameters = convert_synonyms(input_parameters, write_all)
    return input_parameters

#--------------------------------------------------
def convert_synonyms(input_parameters, write_all):
    if input_parameters['physics_flags']['adiabatic_option']=='default': 
        input_parameters['physics_flags']['adiabatic_option'] = 'no-field-line-average-term'
    if input_parameters['physics_flags']['adiabatic_option']=='iphi00=2': 
        input_parameters['physics_flags']['adiabatic_option'] = 'field-line-average-term'
        
    # Turn on all diagnostics is <write_all> = True
    # Except for radial diagnostics since they'll give memory errors in stella
    if write_all: 
        input_parameters['stella_diagnostics_knobs']['write_omega'] = True
        input_parameters['stella_diagnostics_knobs']['write_phi_vs_time'] = True
        input_parameters['stella_diagnostics_knobs']['write_gvmus'] = True
        input_parameters['stella_diagnostics_knobs']['write_gzvs'] = True
        input_parameters['stella_diagnostics_knobs']['write_kspectra'] = True
        input_parameters['stella_diagnostics_knobs']['write_moments'] = True
        input_parameters['stella_diagnostics_knobs']['write_fluxes_kxkyz'] = True
        input_parameters['stella_diagnostics_knobs']['write_radial_fluxes'] = False
        input_parameters['stella_diagnostics_knobs']['write_radial_moments'] = False
    return input_parameters

#--------------------------------------------------
def update_stella05ToStellaMaster():
    renamed_variables =  {
        'kt_grids_box_parameters:centered_in_rho' : 'parameters_kxky_grids_box:centered_in_rho',
        'kt_grids_box_parameters:randomize_phase_shift' : 'parameters_kxky_grids_box:randomize_phase_shift',
        'kt_grids_box_parameters:nx' : 'parameters_kxky_grids_box:nx',
        'kt_grids_box_parameters:ny' : 'parameters_kxky_grids_box:ny',
        'kt_grids_box_parameters:jtwist' : 'parameters_kxky_grids_box:jtwist',
        'kt_grids_box_parameters:jtwistfac' : 'parameters_kxky_grids_box:jtwistfac',
        'kt_grids_box_parameters:phase_shift_angle' : 'parameters_kxky_grids_box:phase_shift_angle',
        'kt_grids_box_parameters:x0' : 'parameters_kxky_grids_box:x0',
        'kt_grids_box_parameters:y0' : 'parameters_kxky_grids_box:y0',
        'kt_grids_box_parameters:periodic_variation' : 'parameters_kxky_grids_box:periodic_variation', 
        'dissipation:d_hyper' : 'hyper:d_hyper',
        'dissipation:cfac2' : 'collisions_fp:cfac2',
        'dissipation:eiediffknob' : 'collisions_fp:eiediffknob',
        'dissipation:lmax' : 'collisions_fp:lmax',
        'dissipation:density_conservation_tp' : 'collisions_fp:density_conservation_tp',
        'dissipation:intraspec' : 'collisions_fp:intraspec',
        'dissipation:ieknob' : 'collisions_fp:ieknob',
        'dissipation:spitzer_problem' : 'collisions_fp:spitzer_problem',
        'dissipation:mu_operator' : 'collisions_fp:mu_operator',
        'dissipation:nvel_local' : 'collisions_fp:nvel_local',
        'dissipation:interspec' : 'collisions_fp:interspec',
        'dissipation:vpa_operator' : 'collisions_fp:vpa_operator',
        'dissipation:jmax' : 'collisions_fp:jmax',
        'dissipation:i2fac' : 'collisions_fp:i2fac',
        'dissipation:nuxfac' : 'collisions_fp:nuxfac',
        'dissipation:iiknob' : 'collisions_fp:iiknob',
        'dissipation:advfield_coll' : 'collisions_fp:advfield_coll',
        'dissipation:exact_conservation' : 'collisions_fp:exact_conservation',
        'dissipation:density_conservation_field' : 'collisions_fp:density_conservation_field',
        'dissipation:no_j0l2' : 'collisions_fp:no_j0l2',
        'dissipation:testpart' : 'collisions_fp:testpart',
        'dissipation:testpart' : 'collisions_fp:testpart',
        'dissipation:no_j1l1' : 'collisions_fp:no_j1l1',
        'dissipation:eimassr_approx' : 'collisions_fp:eimassr_approx',
        'dissipation:eiknob' : 'collisions_fp:eiknob',
        'dissipation:deflknob' : 'collisions_fp:deflknob',
        'dissipation:no_j1l2' : 'collisions_fp:no_j1l2',
        'dissipation:cfac' : 'collisions_fp:cfac',
        'dissipation:i1fac' : 'collisions_fp:i1fac',
        'dissipation:iiknob' : 'collisions_fp:iiknob',
        'dissipation:eeknob' : 'collisions_fp:eeknob',
        'dissipation:fieldpart' : 'collisions_fp:fieldpart',
        'dissipation:density_conservation' : 'collisions_fp:density_conservation',
        'dissipation:momentum_conservation' : 'collisions_dougherty:momentum_conservation',
        'dissipation:mu_operator' : 'collisions_dougherty:mu_operator',
        'dissipation:vpa_operator' : 'collisions_dougherty:vpa_operator',
        'dissipation:energy_conservation' : 'collisions_dougherty:energy_conservation',
        'knobs:fapar' : 'parameters_physics:include_apar',
        'knobs:fbpar' : 'parameters_physics:include_bpar',
        'knobs:cfl_cushion' : 'parameters_numerical:cfl_cushion_middle',
        'knobs:delt_adjust' : 'parameters_numerical:cfl_cushion_lower',
        'knobs:lu_option' : 'parameters_numerical:lu_option',
        'knobs:delt' : 'parameters_numerical:delt',
        'knobs:mat_gen' : 'parameters_numerical:mat_gen',
        'knobs:delt_option' : 'parameters_numerical:delt_option',
        'knobs:drifts_implicit' : 'parameters_numerical:drifts_implicit',
        'knobs:stream_matrix_inversion' : 'parameters_numerical:stream_matrix_inversion',
        'knobs:zed_upwind' : 'parameters_numerical:zed_upwind',
        'knobs:ky_solve_radial' : 'parameters_numerical:ky_solve_radial',
        'knobs:nstep' : 'parameters_numerical:nstep',
        'knobs:tend' : 'parameters_numerical:tend',
        'knobs:mirror_semi_lagrange' : 'parameters_numerical:mirror_semi_lagrange',
        'knobs:mirror_implicit' : 'parameters_numerical:mirror_implicit',
        'knobs:maxwellian_inside_zed_derivative' : 'parameters_numerical:maxwellian_inside_zed_derivative',
        'knobs:ky_solve_real' : 'parameters_numerical:ky_solve_real',
        'knobs:fphi' : 'parameters_numerical:fphi',
        'knobs:fields_kxkyz' : 'parameters_numerical:fields_kxkyz',
        'knobs:delt_max' : 'parameters_numerical:delt_max',
        'knobs:mirror_linear_interp' : 'parameters_numerical:mirror_linear_interp',
        'knobs:rng_seed' : 'parameters_numerical:rng_seed',
        'knobs:avail_cpu_time' : 'parameters_numerical:avail_cpu_time',
        'knobs:time_upwind' : 'parameters_numerical:time_upwind',
        'knobs:vpa_upwind' : 'parameters_numerical:vpa_upwind',
        'knobs:stream_implicit' : 'parameters_numerical:stream_implicit',
        'knobs:mat_read' : 'parameters_numerical:mat_read',
        'parameters:g_exb' : 'parameters_physics:g_exb',
        'parameters:omprimfac' : 'parameters_physics:omprimfac',
        'parameters:tite' : 'parameters_physics:tite',
        'parameters:zeff' : 'parameters_physics:zeff',
        'parameters:nine' : 'parameters_physics:nine',
        'parameters:rhostar' : 'parameters_physics:rhostar',
        'parameters:irhostar' : 'parameters_physics:irhostar',
        'parameters:beta' : 'parameters_physics:beta',
        'parameters:vnew_ref' : 'parameters_physics:vnew_ref',
        'parameters:g_exbfac' : 'parameters_physics:g_exbfac',
        'physics_flags:const_alpha_geo' : 'debug_flags:const_alpha_geo',
        'physics_flags:adiabatic_option' : 'parameters_physics:adiabatic_option',
        'physics_flags:include_pressure_variation' : 'parameters_physics:include_pressure_variation',
        'physics_flags:include_parallel_streaming' : 'parameters_physics:include_parallel_streaming',
        'physics_flags:radial_variation' : 'parameters_physics:radial_variation',
        'physics_flags:include_mirror' : 'parameters_physics:include_mirror',
        'physics_flags:nonlinear' : 'parameters_physics:nonlinear',
        'physics_flags:include_geometric_variation' : 'parameters_physics:include_geometric_variation',
        'physics_flags:include_parallel_nonlinearity' : 'parameters_physics:include_parallel_nonlinearity',
        'physics_flags:full_flux_surface' : 'parameters_physics:full_flux_surface',
        'time_advance_knobs:flip_flop' : 'parameters_numerical:flip_flop',
        'time_advance_knobs:explicit_option' : 'parameters_numerical:explicit_option',
        'time_advance_knobs:xdriftknob' : 'parameters_physics:xdriftknob',
        'time_advance_knobs:ydriftknob' : 'parameters_physics:ydriftknob',
        'time_advance_knobs:wstarknob' : 'parameters_physics:wstarknob',
        'stella_diagnostics_knobs:write_phi_vs_time' : 'stella_diagnostics_knobs:write_phi_vs_kxkyz',
        'stella_diagnostics_knobs:write_kspectra' : 'stella_diagnostics_knobs:write_phi2_vs_kxky',
        'stella_diagnostics_knobs:write_gvmus' : 'stella_diagnostics_knobs:write_g2_vs_zvpamus',
        'stella_diagnostics_knobs:write_gzvs' : 'stella_diagnostics_knobs:write_g2_vs_zvpas',
        'geo_knobs:overwrite_gradpar' : 'geo_knobs:overwrite_b_dot_grad_zeta',
        'parameters:tite' : 'parameters_physics:tite',
        'physics_flags:nonlinear' : 'parameters_physics:nonlinear',
        'physics_flags:full_flux_surface' : 'parameters_physics:full_flux_surface',
        'kt_grids_box_parameters:ny' : 'parameters_kxky_grids_box:ny',
        }
    deprecated_variables = {
        'knobs' : 'driftkinetic_implicit',
        'sources' : 'remove_zero_projection',
        'sources' : 'include_krook_operator',
        'init_g_knobs' : 'even',
        'vmec_parameters' : 'gradpar_zeta_prefac',
        }
    deprecated_namelist = [
        'time_advance_knobs', 
        'knobs', 
        'physics_flags', 
        'kt_grids_box_parameters', 
        'parameters']
    new_namelists = [
        'parameters_physics', 
        'parameters_numerical', 
        'collisions_dougherty', 
        'collisions_fp', 
        'hyper', 
        'parameters_kxky_grids_box',
        'debug_flags']
    new_variables = {
        'debug_flags:debug_all': -1,
        'debug_flags:fields_all_debug': -1,
        'debug_flags:diagnostics_all_debug': -1,
        'debug_flags:stella_debug': -1,
        'debug_flags:ffs_solve_debug': -1,
        'debug_flags:fields_debug': -1,
        'debug_flags:fields_fluxtube_debug': -1,
        'debug_flags:fields_electromagnetic_debug': -1,
        'debug_flags:fields_ffs_debug': -1,
        'debug_flags:implicit_solve_debug': -1,
        'debug_flags:mirror_terms_debug': -1,
        'debug_flags:neoclassical_terms_debug': -1,
        'debug_flags:parallel_streaming_debug': -1,
        'debug_flags:response_matrix_debug': -1,
        'debug_flags:time_advance_debug': -1,
        'debug_flags:extended_grid_debug': -1,
        'debug_flags:geometry_debug': -1,
        'debug_flags:dist_fn_debug': -1,
        'debug_flags:gyro_averages_debug': -1,
        'debug_flags:diagnostics_debug': -1,
        'debug_flags:diagnostics_parameters': -1,
        'debug_flags:diagnostics_omega_debug': -1,
        'debug_flags:diagnostics_fluxes_fluxtube_debug': -1,
        'debug_flags:fluxes_debug': -1,
        'hyper:hyp_vpa': -1,
        'hyper:use_physical_ksqr': -1,
        'hyper:scale_to_outboard': -1,
        'hyper:d_zed': -1,
        'hyper:d_vpa': -1,
        'hyper:hyp_zed': -1,
        'collisions_fp:exact_conservation_tp': -1,
        'collisions_fp:eideflknob': -1, 
        'parameters_physics:hammett_flow_shear': -1,
        'parameters_physics:suppress_zonal_interaction': -1,
        'parameters_physics:prp_shear_enabled': -1,
        'parameters_numerical:maxwellian_normalization': -1,
        'parameters_numerical:cfl_cushion_upper': -1,
        'parameters_numerical:print_extra_info_to_terminal': -1,
        'parameters_numerical:split_parallel_dynamics': -1,
        'parameters_numerical:stream_iterative_implicit': -1,
        'parameters_numerical:delt_min': -1,
        'parameters_numerical:autostop': -1,
        'parameters_numerical:nitt': -1,
        'parameters_numerical:use_deltaphi_for_response_matrix': -1,  
        'layouts_knobs:kymus_layout': -1,
        'kt_grids_range_parameters:kyspacing_option': -1,
        'sources:source_option': -1,
        'init_g_knobs:oddparity': -1,
        'init_g_knobs:read_many': -1,
        'zgrid_parameters:dkx_over_dky': -1,
        'vmec_parameters:n_tolerated_test_arrays_inconsistencies': -1,
        'vmec_parameters:rectangular_cross_section': -1,
        'vmec_parameters:radial_coordinate': -1,
        'vmec_parameters:zgrid_refinement_factor': -1,
        'stella_diagnostics_knobs:write_distribution_g': -1,
        'stella_diagnostics_knobs:write_distribution_h': -1,
        'stella_diagnostics_knobs:write_distribution_f': -1,
        'stella_diagnostics_knobs:write_phi2_vs_time': -1,
        'stella_diagnostics_knobs:write_apar2_vs_time': -1,
        'stella_diagnostics_knobs:write_bpar2_vs_time': -1,
        'stella_diagnostics_knobs:write_apar2_vs_kxky': -1,
        'stella_diagnostics_knobs:write_bpar2_vs_kxky': -1,
        'stella_diagnostics_knobs:write_apar_vs_kxkyz': -1,
        'stella_diagnostics_knobs:write_bpar_vs_kxkyz': -1,
        'stella_diagnostics_knobs:write_omega_vs_kxky': -1,
        'stella_diagnostics_knobs:write_omega_avg_vs_kxky': -1,
        'stella_diagnostics_knobs:write_all': -1,
        'stella_diagnostics_knobs:write_fluxes_vs_time': -1,
        'stella_diagnostics_knobs:write_fluxes_kxky': -1,
        'stella_diagnostics_knobs:write_g2_vs_vpamus': -1,
        'stella_diagnostics_knobs:write_g2_vs_zmus': -1,
        'stella_diagnostics_knobs:write_g2_vs_kxkyzs': -1,
        'stella_diagnostics_knobs:autostop': -1,
        }
    return renamed_variables, new_variables, deprecated_variables, deprecated_namelist, new_namelists




#===============================================================================
#                             RUN AS MAIN SCRIPT                               #
#===============================================================================
 
if __name__ == "__main__":
    path = pathlib.Path(os.path.realpath(__file__))
    path_input_file = path.parent / 'miller_nonlinear_CBC.in'
    input_parameters = convert_inputFileToOlderStellaVersions(path_input_file)
    for namelist in input_parameters.keys():
        print('\n==================================')
        print(f'{namelist}'.center(34))
        print('==================================')
        for variable, value in input_parameters[namelist].items():
            print(f'{variable}: {value}')

    
