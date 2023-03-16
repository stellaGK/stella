#!/usr/bin/python3
''' Run read_wout() as a bash command on the command prompt.'''

# Load the modules
import sys, os, pathlib
import numpy as np
import netCDF4 as nc4 

# Tell python where to find the personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.utils import get_filesInFolder
from stellapy.data import read_wout
from stellapy.utils.config import turnOnVerboseWrapper_configurationFile, turnOffVerboseWrapper_configurationFile
from stellapy.calculations.calculate_gridDivisionsAndSize import get_dkx, get_shat

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
         
        # Position
        rho = 0.75
        svalue = rho*rho
         
        # Get the variables
        turnOffVerboseWrapper_configurationFile()
        wout_variables = read_wout(pathlib.Path(os.getcwd()), svalue=svalue)
    
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
        check_parameters_wout(wout_variables, wout_file_paths[0], svalue)
        
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
            value = list(wout_file.variables[variable][:])
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
    print('\n', '    ', '{0:20}'.format("Variable"), '{0:10}'.format("Shape"), '{0:30}'.format("     Value"), '{0:35}'.format("Function"), '{0:40}'.format("Long name"), '\n')
    for name, value in wout_file.variables.items():
        variable            = name
        function            = variable + str(value).split(" "+variable)[1].split('\n')[0].replace(" ", "")     
        shape               = str(value).split("shape = ")[1].split('\n')[0].replace(" ", "")
        try:    long_name   = str(value).split("long_name: ")[1].split('\n')[0]
        except: long_name   = "" 
        try:    
            value = list(wout_file.variables[variable][:]) 
        except: 
            value = [wout_file.variables[variable][:]] 
        if np.array(value).size > 1:     
            value = np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[:10].split(" ")[0]\
                     + ' ... ' + \
                    np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[-10:].split(" ")[-1]  
            value = value.replace("\n", "") 
        else:  
            try:                          
                value = str(value[0])
            except: 
                value = 'nan'
        if value == 'nan':               
            value = " /  "
        
        print('    ', '{0:20}'.format(variable+':'), '{0:10}'.format(shape), '{0:30}'.format(value), '{0:35}'.format(function), '{0:40}'.format(long_name))
  
        
def check_parameters_wout(wout_variables, wout_file_path, svalue):
    ''' Check the parameters related to the magnetic field. '''

    from stellapy.utils.decorators import restrict_sizeList
    from stellapy.utils.decorators.verbose import indent

    if "NCSX" in str(wout_file_path) or "ncsx" in str(wout_file_path): y0 = 15; device="NCSX"
    if "LHD"  in str(wout_file_path) or "lhd"  in str(wout_file_path): y0 = 15; device="LHD"
    if "TJ"   in str(wout_file_path) or "tj"   in str(wout_file_path): y0 = 10; device="TJII"
    if "W7"   in str(wout_file_path) or "w7x"  in str(wout_file_path): y0 = 20;   device="W7X"
    print()
    print(indent, 'Variables in wout*.nc file:', wout_file_path)
    print(indent, '    '+'div along s=(r/a)^2=rho^2' + ": ", wout_variables['ns'])
    print(indent, '    '+'Aminor_p' + ":                  ", wout_variables['Aminor_p'])
    print(indent, '    '+'iotas' + ":                     ", restrict_sizeList(wout_variables['iotas'],5))
    print(indent, '    '+'iota at rho = ' + str(round(svalue,2)) + ":        ", wout_variables['iota'])  
    print(indent, '    '+'shat at rho = ' + str(round(svalue,2)) + ":        ", get_shat(wout_file_path, svalue, iota=wout_variables['iota'], diotaOverds=None, ns=wout_variables['ns'], iotaf=wout_variables['iotaf']))  
    print(indent, '    '+'B0 at rho =   ' + str(round(svalue,2)) + ":        ", wout_variables['b0'])  
    print()
    print(indent, 'If we want 1 poloidal turn in this stellarator:')
    print(indent, '    nfieldperiods(s) = #poloidalturns * #fieldperiods / iota(s)')
    print(indent, '    nfieldperiods(s) = 1 * ', wout_variables['nfp'] ,'/', wout_variables['iota'], '=', wout_variables['nfp']/wout_variables['iota']) 
    print()
    print(indent, 'Step size in (kx,ky) for 1 poloidal turn in '+device+':')
    nfield_periods = 1*wout_variables['nfp']/wout_variables['iota']
    print(indent, '    nfp = ', nfield_periods)
    print(indent, '    dky = ', 1/y0)
    print(indent, '    dkx = ', get_dkx(wout_file_path, svalue, nfield_periods, y0, iota=wout_variables['iota'], phi=wout_variables['phi'],\
                                        nfp=wout_variables['nfp'], ns=wout_variables['ns'], iotaf=wout_variables['iotaf'])) 
    print()
    print(indent, 'Step size in (kx,ky) for 3 poloidal turns in '+device+':')
    nfield_periods = 3*wout_variables['nfp']/wout_variables['iota']
    print(indent, '    nfp = ', nfield_periods)
    print(indent, '    dky = ', 1/y0)
    print(indent, '    dkx = ', get_dkx(wout_file_path, svalue, nfield_periods, y0, iota=wout_variables['iota'], phi=wout_variables['phi'],\
                                        nfp=wout_variables['nfp'], ns=wout_variables['ns'], iotaf=wout_variables['iotaf'])) 
    print()
    print(indent, 'Step size in (kx,ky) for x poloidal turn in '+device+':')
    if "NCSX" in str(wout_file_path) or "ncsx" in str(wout_file_path): nfield_periods=5.259
    if "LHD"  in str(wout_file_path) or "lhd"  in str(wout_file_path): nfield_periods=13.1
    if "TJ"   in str(wout_file_path) or "tj"   in str(wout_file_path): nfield_periods=2.5012
    if "W7X"  in str(wout_file_path) or "w7x"  in str(wout_file_path): nfield_periods=5.6
    print(indent, '    nfp = ', nfield_periods)
    print(indent, '    dky = ', 1/y0)
    print(indent, '    dkx = ', get_dkx(wout_file_path, svalue, nfield_periods, y0, iota=wout_variables['iota'], phi=wout_variables['phi'],\
                                        nfp=wout_variables['nfp'], ns=wout_variables['ns'], iotaf=wout_variables['iotaf'])) 

# Execute the python script
turnOffVerboseWrapper_configurationFile()
print_decoratorOnCommandPrompt("begin")
main_program()
print_decoratorOnCommandPrompt("end")

    

