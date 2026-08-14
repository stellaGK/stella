
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
# Structure of the new geometry file (12 columns):
#    alpha; zed; zeta; bmag; b_dot_gradz; grady_dot_grady; gradx_dot_grady;
#    gradx_dot_gradx; B_times_gradB_dot_grady; B_times_kappa_dot_grady;
#    B_times_gradB_dot_gradx; bmag_psi0
#
# The new netcdf additionally provides b_dot_gradz_avg and B_times_kappa_dot_gradx.
#
# Old variable names like gradpar, gds2, gds21, gds22, gbdrift, cvdrift and
# their derivatives are no longer written by stella; see geometry.f90 for the
# (old = new * factor) mapping. gds23 and gds24 are also dropped from the output.
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

    # Read the geometric quantities written by the new stella netcdf
    keys = ["zed", "alpha", "zeta", "bmag", "bmag_psi0", "grho", "jacob", "djacdrho", "btor", "shat",
            "b_dot_gradz", "b_dot_gradz_avg",
            "grady_dot_grady", "gradx_dot_grady", "gradx_dot_gradx",
            "B_times_gradB_dot_grady", "B_times_gradB_dot_gradx",
            "B_times_kappa_dot_grady", "B_times_kappa_dot_gradx"]
    with nc4.Dataset(path) as f:
        for key in keys:
            if key in f.variables.keys():
                geo_data[key] = copy.deepcopy(f.variables[key][:])
    return geo_data

#-------------------------------------
def read_geometryFromTxtFile(path):

    # Store the data for the plots in a nested dictionary
    keys = ["alpha", "zed", "zeta", "bmag", "b_dot_gradz",
            "grady_dot_grady", "gradx_dot_grady", "gradx_dot_gradx",
            "B_times_gradB_dot_grady", "B_times_kappa_dot_grady", "B_times_gradB_dot_gradx",
            "bmag_psi0"]
    geo_data = {"source" : path.name}

    # Read the scalars at the top of the file
    with open(path, 'r') as vmec_file:
        header_line = vmec_file.readline() # @UnusedVariable
        first_line  = vmec_file.readline().split()
        for i, key in enumerate(["rhoc", "qinp", "shat", "rhotor", "aref", "bref", "dxdpsi", "dydalpha", "exb_nonlin", "flux_fac", "one_over_nablarho"]):
            if i+1 < len(first_line):
                try: geo_data[key] = float(first_line[i+1])
                except ValueError: pass

    # Read all the relevant arrays from the new geometry file (12 columns)
    try:
        geometry_data = np.loadtxt(path, dtype='float', skiprows=4).reshape(-1, len(keys))
    except:
        # Perhaps there's an +100 instead of +E100:
        input_data = open(path, 'r')
        input_text = input_data.read()
        new_file = open(str(path)+"_corrected", "w" )
        new_file.write(input_text.replace("+1", "E+1").replace("-3", "E-3").replace("EE", "E"))
        new_file.close(); input_data.close()
        try:
            geometry_data = np.loadtxt(str(path)+"_corrected", dtype='float', skiprows=4).reshape(-1, len(keys))
        except:
            exit_reason = "Something went wrong when reading the geometry in:\n"
            exit_reason += str(path)
            exit_program(exit_reason, read_geometryFromTxtFile, sys._getframe().f_lineno)

    # Map columns to keys
    for i, key in enumerate(keys):
        geo_data[key] = geometry_data[:, i]

    # Make sure the data is given as a function of (nzed, nalpha)
    alphas = list(set(geo_data["alpha"].reshape(-1))); nalpha = len(alphas)
    zed = list(set(geo_data["zed"].reshape(-1))); nzed = len(zed)
    for key in keys:
        geo_data[key] = geo_data[key].reshape(nalpha, nzed).T

    # Return the geometry data
    return geo_data

#-------------------------------------
def read_geometryFromH5File(path):

    # Recall what file we are reading from
    geo_data = {"source" : path.name}

    # Read the wout data from the "*.geo.h5" file
    with h5py.File(path, 'r') as f:
        for key in ["alpha","zed", "zeta", "bmag", "b_dot_gradz", "b_dot_gradz_avg",
                    "grady_dot_grady", "gradx_dot_grady", "gradx_dot_gradx",
                    "B_times_gradB_dot_grady", "B_times_gradB_dot_gradx",
                    "B_times_kappa_dot_grady", "B_times_kappa_dot_gradx",
                    "bmag_psi0"]:
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

    # Attach the geometry data (new names)
    self.alpha = geo_data['alpha']
    self.zed = geo_data['zed']
    self.zeta = geo_data['zeta']
    self.b_dot_gradz = geo_data.get('b_dot_gradz')
    self.b_dot_gradz_avg = geo_data.get('b_dot_gradz_avg')
    self.grady_dot_grady = geo_data.get('grady_dot_grady')
    self.gradx_dot_grady = geo_data.get('gradx_dot_grady')
    self.gradx_dot_gradx = geo_data.get('gradx_dot_gradx')
    self.B_times_gradB_dot_grady = geo_data.get('B_times_gradB_dot_grady')
    self.B_times_gradB_dot_gradx = geo_data.get('B_times_gradB_dot_gradx')
    self.B_times_kappa_dot_grady = geo_data.get('B_times_kappa_dot_grady')
    self.B_times_kappa_dot_gradx = geo_data.get('B_times_kappa_dot_gradx')
    self.bmag_psi0 = geo_data.get('bmag_psi0')
    self.btor = geo_data.get('btor')

    # gds23 / gds24 are no longer written by stella but kept for old runs
    self.gds23 = geo_data.get('gds23')
    self.gds24 = geo_data.get('gds24')
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


