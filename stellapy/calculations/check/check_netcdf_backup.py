#!/usr/bin/python3
 
#################################################################
#  STELLAPY.DATA.MAGNETICFIELD AND STELLAPY.PLOT.MAGNETICFIELD
#################################################################

# Load the modules
import sys, os, pathlib, time, copy
import netCDF4 as nc4 
import numpy as np 

# Load the personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)
from stellapy.utils import get_filesInFolder
from stellapy.utils.decorators.verbose import indent

#======================
# Read the NETCDF file
#======================

# Print the size of the vectors
print_size = False

# Get the folder
netcdf_folder = pathlib.Path(os.getcwd())

# Load the file
if get_filesInFolder(netcdf_folder, end="out.nc") or get_filesInFolder(netcdf_folder, start="restart.nc"):
    if get_filesInFolder(netcdf_folder, end="out.nc"):
        netcdf_paths = get_filesInFolder(netcdf_folder, end="out.nc")
    elif get_filesInFolder(netcdf_folder, start="restart.nc")[0]:
        netcdf_paths = get_filesInFolder(netcdf_folder, start="restart.nc") 
        netcdf_paths = netcdf_paths[:1]
        
    # Read the file
    for netcdf_path in netcdf_paths:
        netcdf_file = nc4.Dataset(netcdf_path)
        #print(netcdf_file)
        #print(list(netcdf_file.variables.keys()))
        header = str(netcdf_file.variables["code_info"]).split("stella")[-1]
        date = str(int(header.split("c1: Date: ")[-1].split("\n")[0]))
        date = date[0:4] + "/" + date[4:6]  + "/" + date[6:]
        hour = float(header.split("c2: Time:")[-1].split("+")[0])
        hour = time.strftime('%H:%M:%S', time.gmtime(hour))
        
        # Print the dimensions
        print()
        print(indent, '######################################################')
        print(indent, '                STELLA SIMULATION DATA ')
        print(indent, '######################################################')
        print("    Date:           ", date, hour)
        print("    Netcdf version: ", header.split("c3: ")[-1].split("\n")[0])
        print("    Netcdf file:    ", netcdf_path)
        print()
        print("   ", header.split("c4: ")[-1].split("\n")[0])
        print("   ", header.split("c5: ")[-1].split("\n")[0])
        print("   ", header.split("c6: ")[-1].split("\n")[0])
        print("   ", header.split("c7: ")[-1].split("\n")[0])
        print("   ", header.split("c8: ")[-1].split("\n")[0])
        print("   ", header.split("c9: ")[-1].split("\n")[0])
        print("   ", header.split("c10: ")[-1].split("\n")[0])
        print("   ", header.split("c11: ")[-1].split("\n")[0])
        print()
        print(indent, '======================')
        print(indent, '  NETCDF DIMENSIONS')
        print(indent, '======================', '\n')
        print(indent, ' ', '{0:16}'.format("Dimensions"), '{0:16}'.format("  Shape"), '{0:35}'.format("     Value"))
        print(indent, ' ', '{0:16}'.format("----------"), '{0:16}'.format("  -----"), '{0:35}'.format("     -----"), '\n')
        for name, value in netcdf_file.dimensions.items():
            variable = name
            shape = str(value).split('size =')[1].split('\n')[0] 
            try:   
                value = netcdf_file.variables[variable][:]
                try:
                    value = list(value)
                except:
                    value = [value]
            except: 
                value = [np.nan]
        
            if np.array(value).size > 1: 
                value = np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[:30] + ' ... ' + \
                        np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[-30:]  
                value = value.replace("\n", "") 
            else:  
                try:                          
                    value = str(value[0])
                except: 
                    value = 'nan'
            if value == 'nan':               
                value = " /  "
            if variable=="kx":
                try:
                    vec_kx = copy.deepcopy(np.copy(netcdf_file.variables['kx'][:]))
                    vec_kx = np.abs(vec_kx[np.nonzero(vec_kx)])
                    Lx = np.max(vec_kx)
                    dkx = np.min(vec_kx) 
                    value = value + "   dkx ="+str("{:9.5f}".format(dkx)) + "   Lx ="+str("{:9.5f}".format(Lx))
                except: pass
            if variable=="ky":
                try:
                    vec_ky = copy.deepcopy(np.copy(netcdf_file.variables['ky'][:]))
                    vec_ky = np.abs(vec_ky[np.nonzero(vec_ky)])
                    Ly = np.max(vec_ky)
                    dky = np.min(vec_ky) 
                    value = value + "   dky ="+str("{:9.5f}".format(dky)) + "   Ly ="+str("{:9.5f}".format(Ly))
                except: pass
                
            # Print the dimensions
            print( indent, '    ', '{0:16}'.format(variable+':'), '{0:16}'.format(shape), '{0:50}'.format(value))
            if variable=="vpa" and False:
                vector = list(netcdf_file.variables[variable])
                vector = np.array2string(np.array(vector), formatter={'float_kind':lambda x: "%.2f" % x})
                print(value)
                print(vector)
                
        # Print the netcdf variables
        print()
        print(indent, '======================')
        print(indent, '   NETCDF VARIABLES')
        print(indent, '======================', '\n')
        if print_size:
            print(indent, '    ', '{0:16}'.format("Variable"), '{0:20}'.format("Shape"), '{0:20}'.format("Value"), '{0:10}'.format("      Size"), '{0:35}'.format("     Function"), '{0:15}'.format("Units"), '{0:25}'.format("Long name"))
            print(indent, '    ', '{0:16}'.format("--------"), '{0:20}'.format("-----"), '{0:20}'.format("-----"), '{0:10}'.format("      ----"), '{0:35}'.format("     --------"), '{0:15}'.format("-----"), '{0:25}'.format("---------"), '\n')
        else:
            print(indent, '    ', '{0:16}'.format("Variable"), '{0:20}'.format("Shape"), '{0:20}'.format("Value"), '{0:40}'.format("     Function"), '{0:15}'.format("Units"), '{0:25}'.format("Long name"))
            print(indent, '    ', '{0:16}'.format("--------"), '{0:20}'.format("-----"), '{0:20}'.format("-----"), '{0:40}'.format("     --------"), '{0:15}'.format("-----"), '{0:25}'.format("---------"), '\n')
        for name, value in netcdf_file.variables.items():
            variable            = name
            function            = variable + str(value).split(" "+variable)[1].split('\n')[0].replace(" ", "")     
            shape               = str(value).split("shape = ")[1].split('\n')[0].replace(" ", "") 
            try:    long_name   = str(value).split("long_name: ")[1].split('\n')[0]
            except: long_name   = ""
            try:    units       = str(value).split("units: ")[1].split('\n')[0]
            except: units       = ""
            shapes = shape.split("(")[-1].split(")")[0].split(",")
            shapes = [ int(s) for s in shapes if s!=""]; size=1
            for i in range(len(shapes)):
                size = size*shapes[i]
                
            # Get the value
            try:   
                value = netcdf_file.variables[variable][:]
                try:
                    value = list(value)
                except:
                    value = [value]
            except: 
                value = [np.nan]
        
            if not np.ma.is_masked(value[0]):
                if np.array(value).size > 1:  
                    value = np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[:10].split(" ")[0]\
                             + ' ... ' + \
                            np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[-10:].split(" ")[-1]  
                    value = value.replace("\n", "").replace("[", "").replace("]", "")  
                    if shape != "()": value = "[" + value + "]"
                else:  
                    try:      
                        value = float(value[0])   
                        if variable in ["q", "shat"]: value = round(value, 8)
                        value = str(value)
                    except: 
                        value = 'nan'
            else:
                value = "Masked"
            if value == 'nan':               
                value = " /  "
                
            # Print a line between each section
            if variable[0].isupper():
                print()
            
            # Print the netcdf variables
            if print_size:
                print( indent, '    ', '{0:16}'.format(variable+':'), '{0:20}'.format(shape), '{0:20}'.format(value), '{0:10}'.format(size), '{0:35}'.format("   "+function), '{0:15}'.format(units), '{0:25}'.format(long_name))
            else:
                print( indent, '    ', '{0:16}'.format(variable+':'), '{0:20}'.format(shape), '{0:20}'.format(value), '{0:40}'.format("   "+function), '{0:15}'.format(units), '{0:25}'.format(long_name[:25]))
        

else:
    print("No 'out.nc' file was found.")