
import copy
import os, sys
import h5py
import pathlib
import numpy as np  
import netCDF4 as nc4 
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.data.output.read_outputFile import read_outputFile
from stellapy.data.output.read_outputFile import read_netcdfVariables 

#===============================================================================
#                          READ THE GEOMETRY FILE
#===============================================================================
# Structure of the new geometry file (14 columns):
#    alpha; zed; zeta; bmag; gradpar; gds2; gds21; gds22; gds23; 
#    gds24; gbdrift; cvdrift; gbdrift0; bmag_psi0 
#
# Structure of the new geometry file (13 columns):
#    alpha; zed; zeta; bmag; gradpar; gds2; gds21; gds22; gds23; 
#    gds24; gbdrift; cvdrift; gbdrift0
# 
# B_ref*(pi*a^2)=toroidal_flux_for_the_last_closed_flux_surface_over_two_p
#
# TODO: 
#    Calculate the curvature
#===============================================================================

def read_geometry(path):
    ''' Read the "*.geometry" file. ''' 
    
    # Read the *.geo.h5 or *unique.geometryX.h5 file  
    if os.path.isfile(path.geometry):
        return read_geometryFromH5File(path.geometry) 
    
    # Read the *.geometry file 
    elif os.path.isfile(path.geometry_stella):
        return read_geometryFromTxtFile(path.geometry_stella)
    
    # Read the *.out.nc file  
    elif os.path.isfile(path.geometry):
        return read_geometryFromNetcdfFile(path.output_stella) 
    
    # Critical error if we didn't find any data
    exit_reason = "The geometry data can not be found. Make sure that\n" 
    exit_reason += "either the *.geometry or *.geo.h5 file is present.\n"
    exit_reason += str(path.geometry)
    exit_reason += str(path.geometry_stella) 
    exit_program(exit_reason, read_geometry, sys._getframe().f_lineno)
    return 

#-------------------------------------
def read_geometryFromNetcdfFile(path):
    
    # Recall what file we are reading from
    geo_data = {"source" : path.name}

    # Read the wout data from the "*.out.nc" file
    with nc4.Dataset(path) as f:
        for key in ["zed", "alpha", "bmag", "b_dot_grad_z", "gradpar", "gbdrift", "gbdrift0", "cvdrift", "cvdrift0", "gds2", "gds21", "gds22", "grho", "jacob", "djacdrho"]:
            if key in f.variables.keys():
                geo_data[key] = copy.deepcopy(f.variables[key][:])
    return geo_data

#-------------------------------------
def read_geometryFromTxtFile(path): 

    # Store the data for the plots in a nested dictionary
    keys = ["alpha", "zed", "zeta", "bmag", "gradpar", "gds2", "gds21", "gds22", "gds23", "gds24", "gbdrift", "cvdrift", "gbdrift0", "bmag_psi0", "btor"]
    geo_data = {"source" : path.name}

    # Read the scalars at the top of the file
    with open(path, 'r') as vmec_file:
        header_line = vmec_file.readline() # @UnusedVariable
        first_line  = vmec_file.readline().split()
        for i, key in enumerate(["rhoc", "qinp", "shat", "rhotor", "aref", "bref", "dxdXcoord", "dydalpha", "exb_nonlin", "exb_nonlin_p"]):
            geo_data[key] = float(first_line[i+1])
    
    # Read all the relevant arrays from the new geometry file:  
    try: 
        geometry_data = np.loadtxt(path, dtype='float', skiprows=4).reshape(-1, 15) 
        for key in keys[:-1]: geo_data[key] = geometry_data[:,keys.index(key)] 
    except:  
        try: 
            geometry_data = np.loadtxt(path, dtype='float', skiprows=4).reshape(-1, 14)
            for key in keys: geo_data[key] = geometry_data[:,keys.index(key)]
        except:
            try:
                geometry_data = np.loadtxt(path, dtype='float', skiprows=4).reshape(-1, 13)
                for key in keys[:-2]: geo_data[key] = geometry_data[:,keys.index(key)]
            except:
                # Perhaps there's an +100 instead of +E100:
                input_data = open(path, 'r')
                input_text = input_data.read() 
                new_file = open(str(path)+"_corrected", "w" )  
                new_file.write(input_text.replace("+1", "E+1").replace("-3", "E-3").replace("EE", "E"))  
                new_file.close(); input_data.close()
                try: 
                    geometry_data = np.loadtxt(str(path)+"_corrected", dtype='float', skiprows=4).reshape(-1, 14)
                    for key in keys[:-1]: geo_data[key] = geometry_data[:,keys.index(key)] 
                except:  
                    try: 
                        geometry_data = np.loadtxt(str(path)+"_corrected", dtype='float', skiprows=4).reshape(-1, 15)
                        for key in keys: geo_data[key] = geometry_data[:,keys.index(key)]
                    except:
                        try:
                            geometry_data = np.loadtxt(str(path)+"_corrected", dtype='float', skiprows=4).reshape(-1, 13)
                            for key in keys[:-2]: geo_data[key] = geometry_data[:,keys.index(key)]
                        except:
                            exit_reason = "Something went wrong when reading the geometry in:\n"
                            exit_reason += str(path)
                            exit_program(exit_reason, read_geometryFromTxtFile, sys._getframe().f_lineno)   
    
    # Make sure the data is given as a function of (nzed, nalpha)
    alphas = list(set(geo_data["alpha"].reshape(-1))); nalpha = len(alphas)
    zed = list(set(geo_data["zed"].reshape(-1))); nzed = len(zed)
    for key in keys:
        if key in geo_data.keys():
            geo_data[key] = geo_data[key].reshape(nalpha, nzed).T
    
    # Return the geometry data
    return geo_data
        
#-------------------------------------
def read_geometryFromH5File(path):
    
    # Recall what file we are reading from
    geo_data = {"source" : path.name}

    # Read the wout data from the "*.geo.h5" file
    with h5py.File(path, 'r') as f: 
        for key in ["alpha","zed", "zeta", "bmag", "gradpar", "gds2", "gds21", "gds22", "gds23", "gds24", "gbdrift", "cvdrift", "gbdrift0", "bmag_psi0"]:
            if key in f.keys():
                geo_data[key] = f[key][()]  
        for key in ["qinp", "shat", "rhotor", "aref", "bref"]:
            if key in f.keys():
                geo_data[key] = f[key][()]
                
    # Return the geometry data
    return geo_data

#===============================================================================
#                ATTACH THE GEOMETRY DATA TO THE <GEOMETRY> OBJECT
#===============================================================================
def get_geometryData(self):
             
    # Read the geometry file
    geo_data = read_geometry(self.path)
    
    # Attach the geometry scalars
    self.aref   = geo_data['aref']
    self.bref   = geo_data['bref']
    self.rhotor = geo_data['rhotor'] 
    self.qinp   = geo_data['qinp'] 
    return

#-------------------------
def get_geometryDataExtra(self):
             
    # Read the geometry file
    geo_data = read_geometry(self.path)
    
    # Attach the geometry scalars
    self.aref   = geo_data['aref']
    self.bref   = geo_data['bref']
    self.rhotor = geo_data['rhotor'] 
    self.qinp   = geo_data['qinp']
    
    # Attach the geometry data
    self.zeta = geo_data['zeta'] 
    self.gbdrift = geo_data['gbdrift'] 
    self.alpha = geo_data['alpha']
    self.zed = geo_data['zed']  
    self.gradpar = geo_data['gradpar']
    self.gds2 = geo_data['gds2']
    self.gds21 = geo_data['gds21']
    self.gds22 = geo_data['gds22']
    self.gds23 = geo_data['gds23']
    self.gds24 = geo_data['gds24'] 
    self.cvdrift = geo_data['cvdrift']
    self.gbdrift0 = geo_data['gbdrift0']
    self.bmag_psi0 = geo_data['bmag_psi0'] if "bmag_psi0" in geo_data.keys() else np.nan
    self.btor = geo_data['btor'] if "btor" in geo_data.keys() else np.nan
    return

################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    import timeit; start = timeit.timeit()
    input_file = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_GEOMETRY_OBJECT/input_ky3.3125.in") 
    geo_data = read_geometry(input_file)
    print("\n"+"  "+"-"*30+"\n"+" "*10+"GEOMETRY DATA"+"\n"+"  "+"-"*30+"\n")
    for key, value in geo_data.items():
        if isinstance(value, np.ndarray):
            if value.ndim==1: 
                value = "["+str(value[0])+" ,"+str(value[1])+", ..., "+str(value[-2])+" ,"+str(value[-1])+"]"
        print('{:15}'.format("  "+key), value)
    print("\n"+"  Elapsed time: "+str((start-timeit.timeit())*1000)+" ms")
    
    