#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Hanne Thienpondt

import copy
import os, sys
import pathlib
import f90nml, json

#--------------------------------------------------------------
def convert_inputFile(path_input_file='', stella_version_original='v0.5', stella_version_objective='master'):

    # Valid versions
    stella_versions = ['v0.5', 'v0.6', 'v0.7', 'v0.8', 'master']
    if (stella_version_original not in stella_versions) or (stella_version_objective not in stella_versions):
        print('\nABORT: The supplied stella versions are not recognized.')
        print(f'Choose from: {stella_versions}\n')
        sys.exit()
    
    # Read the default parameters of the two mentioned branches of stella
    folder = pathlib.Path(os.path.realpath(__file__)) 
    folder = pathlib.Path(str(folder).split('AUTOMATIC_TESTS')[0]) 
    folder = folder / 'AUTOMATIC_TESTS/convert_input_files'
    input_file_original_stella = folder / f'input_stella_{stella_version_original}.in'
    input_file_objective_stella = folder / f'input_stella_{stella_version_objective}.in'
    
    # Read the text files
    with open(input_file_original_stella) as f:
        text_input_file_original_stella = f.read()
    with open(input_file_objective_stella) as f:
        text_input_file_objective_stella = f.read()
        
    # Get the namelists
    namelists_original_stella = []
    namelists_objective_stella = []
    for line in text_input_file_original_stella.split('\n'):
        if '&' in line: namelists_original_stella.append(line.replace(' ','').replace('&',''))
    for line in text_input_file_objective_stella.split('\n'):
        if '&' in line: namelists_objective_stella.append(line.replace(' ','').replace('&',''))
    namelists_original_stella = list(set(namelists_original_stella))
    namelists_objective_stella = list(set(namelists_objective_stella))
    for i in range(len(namelists_original_stella)):
        print(namelists_original_stella[i])
    print('---------')
    for i in range(len(namelists_objective_stella)):
        print(namelists_objective_stella[i])
    return 
    
    # Read default input parameters
    input_parameters_original = f90nml.read(input_file_original_stella)
    input_parameters_original = json.loads(json.dumps(input_parameters_original))
    input_parameters_objective = f90nml.read(input_file_objective_stella)
    input_parameters_objective = json.loads(json.dumps(input_parameters_objective))
    
    # Read the real input file
    input_parameters = f90nml.read(path_input_file)
    input_parameters = json.loads(json.dumps(input_parameters))
    
    # Iterate over the values in the input file and update them
    for knob in input_parameters.keys():
        for value in input_parameters[knob]:
            pass
        
    # First update the number of species
    try: input_parameters['species_knobs'].update(input_parameters_read['species_knobs'])
    except: pass
     
    # Add more default species if nspec>2
    for i in range(2, int(input_parameters['species_knobs']['nspec'])+1): 
        input_parameters['species_parameters_'+str(i)] = input_parameters['species_parameters_1'].copy() 

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

#===============================================================================
#                             RUN AS MAIN SCRIPT                               #
#===============================================================================
 
if __name__ == '__main__':
    path = pathlib.Path(os.path.realpath(__file__))
    path_input_file = path.parent / 'miller_nonlinear_CBC.in'
    path_input_file = path.parent / 'W7X_old_input.in'
    input_parameters = convert_inputFile(path_input_file)
    for namelist in input_parameters.keys():
        print('\n==================================')
        print(f'{namelist}'.center(34))
        print('==================================')
        for variable, value in input_parameters[namelist].items():
            print(f'{variable}: {value}')

    
