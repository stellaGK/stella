#!/usr/bin/python3
''' Run check_resolution() as a bash command on the command prompt.'''

# Tell python where to find the personal modules
import sys, os, pathlib

# Tell python where to find the personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.utils import get_filesInFolder
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.data.read_inputParameters import read_floatInput, read_stringInput
from stellapy.config import turnOnVerboseWrapper_configurationFile

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        turnOnVerboseWrapper_configurationFile()
        check_resolution(folder=pathlib.Path(os.getcwd()))

def check_resolution(folder):
    ''' Print the resolution of the simulations in <folder> to the command prompt: 
       {nzed; nfield_periods; nvgrid; nmu; naky; nakx; y0; delt; nstep}
       
    In case the resolution differs between the input files, a list of the
    used resolutions will be shown.
    '''
    # Get the indentation for the text printed to the command prompt
    from stellapy.utils.decorators.verbose_wrapper import indent
    
    # List the parameters of interest: keys must match with those of read_grid()
    resolution = {'nzed' : [], 'nfield_periods' : [], 'nvgrid' : [], 'nmu' : [], 'naky' : [], 'nakx' : [], 'y0' : [], 'delt' : [], 'nstep' : []}
    
    # Go through all the input files present in the <folder>
    input_files = get_filesInFolder(folder, end=".in")
    for input_file in input_files:
        
        # Read the grid parameters
        input_parameters = read_grid(input_file)
        
        # Add them to the dictionary <resolution> if they are unique
        for key in resolution.keys():
            if input_parameters[key] not in resolution[key]:
                resolution[key].append(input_parameters[key])
    
    # Print the resolution parameters
    print()
    print(indent, "Resolution:")
    for key in resolution.keys():
        print(indent, '    ','{0:15}'.format(key),resolution[key])

def read_grid(input_file):
    ''' Read *.in and return input_parameters = {nzed; nperiod; zeta_center; nvgrid; vpa_max; vperp_max; delt; nstep; \
        cfl_cushion; vperp_max; nfield_periods; geo_option; boundary_option; adiabatic_option; vmec_filename} '''    
        
    # Read the *in file and get the parameters for species_1 
    input_file   = open(input_file, 'r')
    input_text   = input_file.read().replace(' ', '') 
    input_parameters                = {} 
    input_parameters['nmu']         = read_floatInput(input_text, 'nmu')
    input_parameters['nvgrid']      = read_floatInput(input_text, 'nvgrid')
    input_parameters['nzed']        = read_floatInput(input_text, 'nzed')
    input_parameters['naky']        = read_floatInput(input_text, 'naky')
    input_parameters['nakx']        = read_floatInput(input_text, 'nakx')
    input_parameters['y0']          = read_floatInput(input_text, 'y0')
    input_parameters['nfield_periods'] = read_floatInput(input_text, 'nfield_periods')
    input_parameters['nperiod']     = read_floatInput(input_text, 'nperiod') 
    input_parameters['zeta_center'] = read_floatInput(input_text, 'zeta_center')
    input_parameters['vpa_max']     = read_floatInput(input_text, 'vpa_max')
    input_parameters['delt']        = read_floatInput(input_text, 'delt')
    input_parameters['nstep']       = read_floatInput(input_text, 'nstep')
    input_parameters['vperp_max']   = read_floatInput(input_text, 'vperp_max')
    input_parameters['cfl_cushion'] = read_floatInput(input_text, 'cfl_cushion')
    input_parameters['geo_option']      = read_stringInput(input_text, 'geo_option')
    input_parameters['boundary_option'] = read_stringInput(input_text, 'boundary_option')
    input_parameters['vmec_filename']   = read_stringInput(input_text, 'vmec_filename')
    input_parameters['adiabatic_option'] = read_stringInput(input_text, 'adiabatic_option')

    # Close the file
    input_file.close()
    return input_parameters

# Execute the python script
print_decoratorOnCommandPrompt("begin")
main_program()
print_decoratorOnCommandPrompt("end") 


    







    




