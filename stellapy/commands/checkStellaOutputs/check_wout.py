#!/usr/bin/python3
''' Run read_wout() as a bash command on the command prompt.'''

# Load the modules
import sys, os, pathlib
import numpy as np
import netCDF4 as nc4 

# Tell python where to find the personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0]) 
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.utils import get_filesInFolder
from stellapy.data import read_wout
from stellapy.config import turnOnVerboseWrapper_configurationFile, turnOffVerboseWrapper_configurationFile


# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
         
        # Get the variables
        turnOffVerboseWrapper_configurationFile()
        wout_variables = read_wout(pathlib.Path(os.getcwd()))
    
        # Get the file path
        wout_file_paths = None
        if get_filesInFolder(pathlib.Path(os.getcwd()), start="wout", end=".nc"): 
            wout_file_paths = get_filesInFolder(pathlib.Path(os.getcwd()), start="wout", end=".nc")
        elif get_filesInFolder(pathlib.Path(os.getcwd()), start="wout", end=".h5"): 
            wout_file_paths = get_filesInFolder(pathlib.Path(os.getcwd()), start="wout", end=".h5")

        # Check all the parameters
        check_netcdfFile(wout_file_paths[0])
            
        # Check the parameters read by us
        turnOnVerboseWrapper_configurationFile()
        check_parameters_wout(wout_variables, wout_file_paths[0])
        
def check_netcdfFile(wout_path):
    
    # Load the file
    wout_file = nc4.Dataset(wout_path)
    
    # Print the dimensions
    print('\n', 'MAGNETIC FIELD DIMENSIONS:')
    print('\n', ' ', '{0:20}'.format("Dimensions"), '{0:20}'.format("Shape"), '{0:35}'.format("Value"), '\n')
    for name, value in wout_file.dimensions.items():
        variable = name
        shape = str(value).split('size =')[1].split('\n')[0]
        try:    
            value = list(wout_file.variables[variable])
        except: 
            value = [np.nan]
    
        if np.array(value).size > 1:     
            value = np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[0:60] + ' ...'
        else:  
            try:                          
                value = str(value[0])
            except: 
                value = 'nan'
        if value == 'nan':               
            value = " /  "
    
        print('    ', '{0:20}'.format(variable+':'), '{0:20}'.format(shape), '{0:35}'.format(value))
    
    # Print the wout variables
    print('\n', 'MAGNETIC FIELD VARIABLES:')
    print('\n', '    ', '{0:20}'.format("Variable"), '{0:20}'.format("Shape"), '{0:35}'.format("Function"), '{0:60}'.format("Long name"), '\n')
    for name, value in wout_file.variables.items():
        variable            = name
        function            = variable + str(value).split(" "+variable)[1].split('\n')[0].replace(" ", "")     
        shape               = str(value).split("shape = ")[1].split('\n')[0].replace(" ", "")
        try:    long_name   = str(value).split("long_name: ")[1].split('\n')[0]
        except: long_name   = ""
        
        print('    ', '{0:20}'.format(variable+':'), '{0:20}'.format(shape), '{0:35}'.format(function), '{0:60}'.format(long_name))
  
        
def check_parameters_wout(wout_variables, wout_file_path):
    ''' Check the parameters related to the magnetic field. '''

    from stellapy.utils.decorators import restrict_sizeList
    from stellapy.utils.decorators.verbose_wrapper import indent

    print()
    print(indent, 'Variables in wout*.nc file:', wout_file_path)
    print(indent, '    '+'div along s=(r/a)^2=rho^2' + ": ", wout_variables['ns'])
    print(indent, '    '+'Aminor_p' + ":                  ", wout_variables['Aminor_p'])
    print(indent, '    '+'iotas' + ":                     ", restrict_sizeList(wout_variables['iotas'],5))
    print(indent, '    '+'iota' + ":                      ", wout_variables['iota']) 
    #print(indent, '    '+'phi' + ":                       ", wout_variables['phi']) 
    print()
    print(indent, 'If we want 1 poloidal turn in w7x:')
    print(indent, '    nfieldperiods(s) = #poloidalturns * #fieldperiods / iota(s)')
    print(indent, '    nfieldperiods(s) = 1 * 5 /', wout_variables['iota'], '=', 5/wout_variables['iota'])
    print()
    
# Execute the python script
turnOffVerboseWrapper_configurationFile()
print_decoratorOnCommandPrompt("begin")
main_program()
print_decoratorOnCommandPrompt("end")

    

