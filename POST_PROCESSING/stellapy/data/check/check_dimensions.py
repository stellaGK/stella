#!/usr/bin/python3
import numpy as np 
import netCDF4 as nc4 
import sys, os, pathlib, copy 

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)
from stellapy.data.dimensions.read_dimensionsAndVectors import read_dimensionsFromNcFile
from stellapy.utils.decorators.verbose import indent
from stellapy.utils.files import get_filesInFolder
if __name__ == "__main__":
    
#===============================================================================
#                            READ THE NETCDF FILE                              #
#===============================================================================

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
            
            # Print the header
            print()  
            print(indent,"".center(110,"=")) 
            print(indent," NETCDF DIMENSIONS ".center(110,"="))
            print(indent,"".center(110,"="), "\n") 
            #print(indent, ' ', '{0:16}'.format("Dimensions"), '{0:16}'.format("  Shape"), '{0:16}'.format("  dstart"), '{0:16}'.format("  dstop"), '{0:35}'.format("     Value"))
            #print(indent, ' ', '{0:16}'.format("----------"), '{0:16}'.format("  -----"), '{0:16}'.format("  ------"), '{0:16}'.format("  -----"), '{0:35}'.format("     -----"), '\n')
            print(indent, ' ', '{0:16}'.format("Dimensions"), '{0:10}'.format("  Shape"), '{0:60}'.format("     Value"), '{0:16}'.format("  Max"), '{0:16}'.format("  Step"))
            print(indent, ' ', '{0:16}'.format("----------"), '{0:10}'.format("  -----"), '{0:60}'.format("     -----"), '{0:16}'.format("  ---"), '{0:16}'.format("  ----"), '\n')
                
            # Quick reading
            try: 
                dimensions = read_dimensionsFromNcFile(netcdf_path)
                for variable in ['vec_z', 'vec_kx', 'vec_ky', 'vec_mu', 'vec_vpa', 'vec_time']:
                    value = dimensions[variable]
                    dstart = str(round(value[1]-value[0],2))
                    dstop = str(round(value[-1]-value[-2],2))
                    length = str(len(value))
                    if np.array(value).size > 1: 
                            value = np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[:20] + ' ... ' + \
                                    np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[-20:]  
                            value = value.replace("\n", "") 
                    if variable=="vec_kx": 
                        vec_kx = copy.deepcopy(np.copy(dimensions[variable][:]))
                        vec_kx = np.abs(vec_kx[np.nonzero(vec_kx)])
                        Lx = np.max(vec_kx)
                        dkx = np.min(vec_kx) 
                        step = "dkx ="+str("{:9.5f}".format(dkx))
                        max = "Lkx ="+str("{:6.2f}".format(Lx)) 
                    if variable=="vec_ky": 
                        vec_ky = copy.deepcopy(np.copy(dimensions[variable]))
                        vec_ky = np.abs(vec_ky[np.nonzero(vec_ky)])
                        Ly = np.max(vec_ky)
                        dky = np.min(vec_ky) 
                        step = "dky ="+str("{:9.5f}".format(dky))
                        max = "Lky ="+str("{:6.2f}".format(Ly))  
                    print(indent, '  ', '{0:17}'.format(variable+':'), '{0:10}'.format(length), '{0:60}'.format(value), '{0:16}'.format(max), '{0:16}'.format(step))
                print()
                
            # Slow reading but works on corrupted files
            except: 
                netcdf_file = nc4.Dataset(netcdf_path) 
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
                        value = np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[:20] + ' ... ' + \
                                np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[-20:]  
                        value = value.replace("\n", "") 
                    else:  
                        try:                          
                            value = str(value[0])
                        except: 
                            value = 'nan'
                    if value == 'nan':               
                        value = " /  "
                    if variable=="kx" or variable=="vec_kx":
                        try:
                            vec_kx = copy.deepcopy(np.copy(netcdf_file.variables['kx'][:]))
                            vec_kx = np.abs(vec_kx[np.nonzero(vec_kx)])
                            Lx = np.max(vec_kx)
                            dkx = np.min(vec_kx) 
                            value = value + "   dkx ="+str("{:9.5f}".format(dkx)) + "   Lkx ="+str("{:9.5f}".format(Lx))
                        except: pass
                    if variable=="ky" or variable=="vec_ky":
                        try:
                            vec_ky = copy.deepcopy(np.copy(netcdf_file.variables['ky'][:]))
                            vec_ky = np.abs(vec_ky[np.nonzero(vec_ky)])
                            Ly = np.max(vec_ky)
                            dky = np.min(vec_ky) 
                            value = value + "   dky ="+str("{:9.5f}".format(dky)) + "   Lky ="+str("{:9.5f}".format(Ly))
                        except: pass
                         
                    # Print the dimensions
                    print( indent, '    ', '{0:16}'.format(variable+':'), '{0:16}'.format(shape), '{0:50}'.format(value))
    
    elif get_filesInFolder(netcdf_folder, end="dimensions"):
        import h5py
        dimensions = {}
        dimensions_paths = get_filesInFolder(netcdf_folder, end="dimensions")
        for path in dimensions_paths:
            with h5py.File(path, 'r') as f: 
                for dim in ["z", "kx", "ky", "mu", "vpa", "time", "species"]:
                    for extra in ["dim_", "vec_"]:
                        key = extra + dim 
                        if key in f.keys():  
                            dimensions[key] = f[key][()] 
            print()  
            print(indent,"".center(110,"=")) 
            print(indent," NETCDF DIMENSIONS ".center(110,"="))
            print(indent,"".center(110,"="), "\n") 
            print(indent, ' ', '{0:16}'.format("Dimensions"), '{0:10}'.format("  Shape"), '{0:50}'.format("     Value"), '{0:16}'.format("  Max"), '{0:16}'.format("  Step"))
            print(indent, ' ', '{0:16}'.format("----------"), '{0:10}'.format("  -----"), '{0:50}'.format("     -----"), '{0:16}'.format("  ---"), '{0:16}'.format("  ----"), '\n')
            for key in dimensions.keys():
                if "dim_" not in key:
                    variable = key
                    shape = str(np.shape(dimensions[key]))
                    value = dimensions[key] 
                    value = list(value)
                    step = ""
                    max = ""
                    if np.array(value).size > 1: 
                        value = np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[:20] + ' ... ' + \
                                np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[-20:]  
                        value = value.replace("\n", "") 
                    if variable=="vec_kx": 
                        vec_kx = copy.deepcopy(np.copy(dimensions[key][:]))
                        vec_kx = np.abs(vec_kx[np.nonzero(vec_kx)])
                        Lx = np.max(vec_kx)
                        dkx = np.min(vec_kx) 
                        step = "dkx ="+str("{:9.5f}".format(dkx))
                        max = "Lx ="+str("{:6.2f}".format(Lx)) 
                    if variable=="vec_ky": 
                        vec_ky = copy.deepcopy(np.copy(dimensions[key]))
                        vec_ky = np.abs(vec_ky[np.nonzero(vec_ky)])
                        Ly = np.max(vec_ky)
                        dky = np.min(vec_ky) 
                        step = "dky ="+str("{:9.5f}".format(dky))
                        max = "Ly ="+str("{:6.2f}".format(Ly))  
                    print(indent, '  ', '{0:17}'.format(key+':'), '{0:10}'.format(shape), '{0:50}'.format(value), '{0:16}'.format(max), '{0:16}'.format(step))
            print()
    else:
        print("No 'out.nc' file was found.")
