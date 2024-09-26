''' Run read_wout() as a bash command on the command prompt.'''

# Load the modules
import numpy as np
import netCDF4 as nc4  
import sys, os, pathlib  

# Tell python where to find the personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_gridDivisionsAndSize  
from stellapy.utils.files.get_filesInFolder import get_filesInFolder
from stellapy.data.geometry.read_wout import read_woutFile 
from stellapy.utils.commandprompt.bash import Bash 

#===============================================================================
#                          Print geometric coefficients                        #
#===============================================================================

def print_geometry(folder, vmec_filename=None, rho=0.25): 
        
    # Folder 
    if vmec_filename!=None: vmec = get_filesInFolder(folder, start=vmec_filename, end=".nc")[0]
    if vmec_filename==None: vmec = get_filesInFolder(folder, start=vmec_filename, end=".nc")[0]
    
    # Readial location
    svalue = rho*rho  

    # Read the VMEC file   
    path = type('Path', (object,), {'vmec' : vmec})
    woutParameters = read_woutFile(path)     

    # Check all the parameters
    check_netcdfFile(path)
        
    # Check the parameters read by us 
    check_parameters_wout(woutParameters, vmec, svalue)
    return

#------------------------------
def check_netcdfFile(path, print_keys=False):
    
    # Load the file
    wout_file = nc4.Dataset(path.vmec)
    
    # Keys
    if print_keys:
        keys = list(wout_file.dimensions.keys()) + list(wout_file.variables.keys())
        keys.sort(); print("KEYS:")
        for i in range(len(keys)):
            print("   "+keys[i])
    
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
            array1 = np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[:20].split(" ") 
            array2 = np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[-20:].split(" ") 
            value = array1[0] + ', '+ array1[1] + ' ... ' + array2[-2] + ', ' + array2[-1]  
            value = array1[0] + ' ... ' + array2[-1]  
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

    rho = str(round(np.sqrt(svalue),2))  
    args = calculate_gridDivisionsAndSize(np.nan, np.nan, wout_variables, svalue=svalue) 
    args = calculate_gridDivisionsAndSize(15, args['nfp']/args['iota'], wout_variables, svalue=svalue) 
    
    print()
    print(indent, 'Variables in wout*.nc file:',str(wout_file_path))
    print(indent, '    '+'div along s=(r/a)^2=rho^2' + ": ", args['ns'])
    print(indent, '    '+'Aminor_p' + ":                  ", args['Aminor_p'])
    print(indent, '    '+'iotas' + ":                     ", restrict_sizeList(args['iotas'],5))
    print()
    print(indent, 'Variables for rho = '+str(rho)+', y0 = 15 and nfield_periods = '+str(args['nfp']/args['iota']))
    print(indent, '    '+'iota at rho = ' + rho + ":        ", args['iota'])  
    print(indent, '    '+'shat at rho = ' + rho + ":        ", args['shat'])  
    print(indent, '    '+'B0 at rho =   ' + rho + ":        ", args['b0'])   
    print(indent, '    '+'dkx at rho =  ' + rho + ":        ", args['dkx'])  
    print(indent, '    '+'dky at rho =  ' + rho + ":        ", args['dky']) 
    print(indent, '    '+'dkx/dky atrho=' + rho + ":        ", args['dkx']/args['dky'])  
    print(indent, '    '+'field_period_ratio at rho =   ' + rho + ":     ", args['field_period_ratio'])   
    print(indent, '    '+'twist_and_shift_geo_fac at rho = ' + rho + ":  ", args['twist_and_shift_geo_fac'])  
    print(indent, '    '+'jtwist at rho = ' + rho + ":                   ", args['jtwist'])  
    print()
    print(indent, 'If we want 1 poloidal turn in this stellarator:')
    print(indent, '    nfieldperiods(s) = #poloidalturns * #fieldperiods / iota(s)')
    print(indent, '    nfieldperiods(s) = 1 * ', args['nfp'] ,'/', args['iota'], '=', args['nfp']/args['iota'])
    print()     
    
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface
    bash = Bash(print_geometry, __doc__)  
    
    # Define the wout 
    bash.add_option('vmec_filename', 'str', 'v', '', 'Give the vmec file name.')
    
    # Give the radial location
    bash.add_option('rho', 'float', 'r', '', 'Define the radial location rho.')
    
    # Get the bash arguments and execute the script 
    print_geometry(**bash.get_arguments() )   
