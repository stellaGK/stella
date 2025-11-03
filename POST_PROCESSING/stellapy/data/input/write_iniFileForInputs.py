  
import os, pathlib, configparser 
from datetime import datetime, timedelta
from stellapy.data.input.read_inputFile import read_inputFile 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder  
from stellapy.data.input.read_inputFile import read_extraInputParameters
from stellapy.data.stella.load_defaultInputParameters import load_defaultInputParameters  
from stellapy.data.input.write_listOfMatchingInputFiles import read_inputFilesInDummyInputFiles

#===============================================================================
#                         WRITE THE INI INPUT FILES
#===============================================================================
  
def write_iniFileForInputs(folder=None):  
    
    # Get the input_files inside <folder> that have a *.out.nc file
    input_files = get_filesInFolder(folder, end='.in') 
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    
    # For linear flux tube simulations, read the dummy input file instead 
    dummy_input_files = get_filesInFolder(folder, end='_dummy.in') 
    
    # Remove the input files in the dummy inputs from the list of input files
    input_files_in_dummy_files = read_inputFilesInDummyInputFiles(dummy_input_files)
    input_files = [i for i in input_files if i not in input_files_in_dummy_files]
    input_files += dummy_input_files
    
    # Collect the parent directories (Ignore OLD directories)
    parent_directories = list(set([i.parent if 'OLD' not in i.parent.name else i.parent.parent for i in input_files]))
    
    # Only look at the input files in each parent folder
    for directory in parent_directories:
        
        # Get the files in <directory>
        input_files_directory = [i for i in input_files if (i.parent==directory) or (i.parent.parent==directory and 'OLD' in i.parent.name)]
        
        # Iterate over the input files
        for input_file in input_files_directory: 
            
            # Only write the *.ini file if it doesn't exist or it's older than the *.in file
            if os.path.isfile(input_file.with_suffix('.ini')):
                if datetime.fromtimestamp(input_file.stat().st_mtime) < (datetime.fromtimestamp(input_file.with_suffix('.ini').stat().st_mtime)+timedelta(minutes=5)):
                    continue
            
            # Read the input parameters
            input_parameters = read_inputFile(input_file) 
            input_parameters = read_extraInputParameters(input_parameters, input_file) 
            
            # Write the *.ini file
            write_iniFile(input_parameters, input_file)

    return
    
#===============================================================================
#                             SAVE THE INPUT FILE                              #
#===============================================================================

def write_iniFile(input_parameters, input_file):
     
    # Define the configuration file
    ini_path = input_file.with_suffix('.ini')
    ini_data = configparser.ConfigParser()
    
    # Get the default input parameters 
    default_parameters = load_defaultInputParameters() 
    default_parameters['adiabatic_electron_response']['teti'] = -1

    # Look at the knobs in the following order
    knobs = list(input_parameters.keys()) 
    knobs = ['='*15+' SPECIES '+'='*15] + sorted([s for s in knobs if 'species_' in s]) + ['adiabatic_electron_response', 'species_options']
    knobs += ['='*15+' GEOMETRY '+'='*15] + ['geometry_options', 'geometry_miller', 'geometry_vmec']
    knobs += ['='*15+' GRIDS '+'='*15] + ['kxky_grid_range', 'kxky_grid_box', 'z_grid', 'z_boundary_condition', 'velocity_grids']
    knobs += ['='*15+' SIMULATION '+'='*15] + ['gyrokinetic_terms', 'scale_gyrokinetic_terms', 'electromagnetic', 'physics_inputs', 'diagnostics']
    knobs += ['='*15+' OTHER '+'='*15] + ['time_trace_options', 'time_step', 'kxky_grid_option', 'initialise_distribution', 'numerical_algorithms', 'dissipation_and_collisions_options', 'hyper_dissipation', 'parallelisation']
    
    # Only save the variables that differ from the default parameters 
    for knob in knobs:
        
        # Add some headers to sort the stella knobs
        if '====' in knob: ini_data.add_section(knob)
        
        # If it's not a header but a knob, then add the values
        if '====' not in knob:
            
            # For all knobs except the specie knob, sort the keys
            keys = list(input_parameters[knob].keys())
            if 'species_' not in knob: keys = sorted(keys)
            
            # Iterate over the keys
            for key in keys: 
                if 'species_' not in knob:
                    if default_parameters[knob][key]!=input_parameters[knob][key]:
                        if knob not in ini_data: ini_data.add_section(knob)
                        ini_data[knob][key] = str(input_parameters[knob][key])
                elif 'species_' in knob:
                    if default_parameters[knob]!=input_parameters[knob]:
                        if key in ['d2ndr2', 'd2Tdr2']:
                            if default_parameters[knob][key]!=input_parameters[knob][key]:
                                if knob not in ini_data: ini_data.add_section(knob)
                                ini_data[knob][key] = str(input_parameters[knob][key])
                        else:
                            if knob not in ini_data: ini_data.add_section(knob)
                            ini_data[knob][key] = str(input_parameters[knob][key])
                            
    # Save the configuration file
    with open(ini_path, 'w') as ini_file:
        ini_data.write(ini_file)
        print('    ----> Saved the input file as ' + ini_path.name)
    
    return ini_path

