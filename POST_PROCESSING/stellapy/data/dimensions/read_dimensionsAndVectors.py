""" 

#===============================================================================
#                           Read the dimensions data                           #
#===============================================================================
 
Read the dimensions of the simulation from the *.dimensions or *.out.nc files.

Returns
-------
    {dim_species, dim_z, dim_kx, dim_ky, dim_mu, dim_vpa, dim_time, 
    vec_z, vec_kx, vec_ky, vec_mu, vec_vpa, vec_time}

Hanne Thienpondt
20/01/2023

"""

#!/usr/bin/python3
import os, sys
import pathlib
import copy, h5py 
import numpy as np

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.dimensions.write_h5FileForDimensions import write_h5FileForDimensions 
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables
from stellapy.data.utils.show_progressWhenReadingFiles import show_progressWhenReadingFiles
from stellapy.data.input.read_inputFile import read_nspecFromInputFile
from stellapy.utils.decorators.exit_program import exit_program

#===============================================================================
#                           Read the dimensions data                           #
#===============================================================================

def read_dimensions(path): 
     
    # Try to write the *.dimensions file   
    if not os.path.isfile(path.dimensions): 
        write_h5FileForDimensions(path.folder)
        
    # Read the dimensions files from multiple simulations 
    if path.multiple_input_files: 
        return read_dimensionsFromMultipleFiles(path.paths, path) 
                                     
    # Read from the *.dimensions file 
    if os.path.isfile(path.dimensions):
        return read_dimensionsFromDimensionsFile(path.dimensions) 
     
    # Read from the *.out.h5 file  
    if os.path.isfile(path.output):
        return read_dimensionsFromH5File(path.output)  
    
    # Read from the *.out.nc file 
    if os.path.isfile(path.output_stella):
        return read_dimensionsFromNcFile(path.output_stella)   

    # Critical error if we didn't find any data
    exit_reason = "The dimension data could not be found:\n" 
    exit_reason += "    "+str(path.dimensions)+"\n" 
    exit_program(exit_reason, read_dimensions, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_dimensionsFromDimensionsFile(path):  
    
    # Initiate
    dimensions = {}

    # Read the wout data from the *.dimensions file
    with h5py.File(path, 'r') as f: 
        for dim in ["z", "kx", "ky", "mu", "vpa", "time", "species"]:
            for extra in ["dim_", "vec_"]:
                key = extra + dim 
                if key in f.keys():  
                    dimensions[key] = f[key][()] 

    # Return the dimensions
    return dimensions

#-------------------------------------
def read_dimensionsFromMultipleFiles(paths, path):  
    
    # Initiate
    dimensions = {} 
    dim_time = 0 

    # Read the dimensions from the *.dimensions file 
    with h5py.File(path.dimensions, 'r') as f: 
        for dim in ["z", "kx", "ky", "mu", "vpa", "time", "species"]:
            for extra in ["dim_", "vec_"]:
                key = extra + dim 
                if key in f.keys():  
                    dimensions[key] = f[key][()] 
                    
    # Get vec_kx and vec_ky by combining the simulations
    for path in paths:  
        if os.path.isfile(path.dimensions):
            with h5py.File(path.dimensions, 'r') as f: 
                dimensions["vec_kx"] = list(dimensions["vec_kx"]) + list(f["vec_kx"][()])
                dimensions["vec_ky"] = list(dimensions["vec_ky"]) + list(f["vec_ky"][()])  
                dim_time = np.max([dim_time, dimensions["dim_time"]]) 
    dimensions["vec_kx"] = np.array(sorted(list(set(dimensions["vec_kx"]))))
    dimensions["vec_ky"] = np.array(sorted(list(set(dimensions["vec_ky"]))))
    dimensions["dim_kx"] = len(dimensions["vec_kx"])
    dimensions["dim_ky"] = len(dimensions["vec_ky"])
    dimensions["dim_time"] = dim_time  

    # Return the dimensions   
    return dimensions

#-------------------------------------
def read_dimensionsFromH5File(path): 
    
    # Initiate
    dimensions = {}
    
    # Read the wout data from the *.dimensions file
    with h5py.File(path, 'r') as f:  
        for key in ["vec_z", "vec_kx", "vec_ky", "vec_mu", "vec_vpa", "vec_time", "vec_species"]: 
            if key in f.keys(): 
                dimensions[key] = f[key][()] 
                
    # Get the dimensions
    if 'dim_z' not in dimensions:
        dimensions['dim_z']   = len(dimensions['vec_z'])
        dimensions['dim_kx']  = len(dimensions['vec_kx'])
        dimensions['dim_ky']  = len(dimensions['vec_ky'])
        dimensions['dim_mu']  = len(dimensions['vec_mu'])
        dimensions['dim_vpa'] = len(dimensions['vec_vpa'])
        dimensions['dim_time'] = len(dimensions['vec_time']) 
                    
    # Get the dimension of the species
    if 'dim_species' not in dimensions:  
        dimensions['dim_species'] = read_nspecFromInputFile(str(path).replace(".out.h5",".in"))
        dimensions['vec_species'] = range(dimensions['dim_species'])
        
    # Return the dimensions
    return dimensions

#-------------------------------------
def read_dimensionsFromNcFile(path): 
    
    # Initiate
    dimensions = {}

    # Get the dimensions from the output file
    netcdf_file = read_outputFile(path)   
    dimensions['vec_z']   = read_netcdfVariables('vec_z',   netcdf_file)
    dimensions['vec_kx']  = read_netcdfVariables('vec_kx',  netcdf_file)
    dimensions['vec_ky']  = read_netcdfVariables('vec_ky',  netcdf_file)
    dimensions['vec_mu']  = read_netcdfVariables('vec_mu',  netcdf_file)
    dimensions['vec_vpa'] = read_netcdfVariables('vec_vpa', netcdf_file)
    dimensions['vec_time'] = read_netcdfVariables('vec_time', netcdf_file)
    dimensions['vec_species'] = read_netcdfVariables('vec_species',  netcdf_file) 
    netcdf_file.close()
    
    # Get the size of the dimensions
    dimensions['dim_z']   = len(dimensions['vec_z'])
    dimensions['dim_kx']  = len(dimensions['vec_kx'])
    dimensions['dim_ky']  = len(dimensions['vec_ky'])
    dimensions['dim_mu']  = len(dimensions['vec_mu'])
    dimensions['dim_vpa'] = len(dimensions['vec_vpa'])
    dimensions['dim_time'] = len(dimensions['vec_time'])
    dimensions['dim_species'] = len(dimensions['vec_species'])
    return dimensions
            
#===============================================================================
#                    ATTACH THE DIMENSIONS AND VECTORS
#===============================================================================

def get_dimensionsAndVectors(self):
    """ Add the kx, ky, time and z dimensions as attributes to the self object."""
    
    # Get the dimension and vector objects
    dim = self.dim
    vec = self.vec
    
    # Initiate the attributes  
    vec.z, vec.kx, vec.ky, vec.kx_stella, vec.mu, vec.vpa, vec.time = [], [], [], [], [], [], []
    dim.species, dim.z, dim.kx, dim.ky, dim.mu, dim.vpa, dim.time= 0, 0, 0, 0, 0, 0, 0
        
    # Show the reading progress in the GUI 
    show_progressWhenReadingFiles(self, "Reading dimensions")
    
    # Read the netcdf data to determine the dimensions of the input files
    dimensions = read_dimensions(self.path)   
    
    # Save the dimensions 
    dim.species  = dimensions['dim_species']  
    dim.z        = dimensions['dim_z']  
    dim.kx       = dimensions['dim_kx']  
    dim.ky       = dimensions['dim_ky']  
    dim.mu       = dimensions['dim_mu'] 
    dim.vpa      = dimensions['dim_vpa']   
    dim.time     = dimensions['dim_time']   
    vec.z        = dimensions['vec_z'] 
    vec.kx       = dimensions['vec_kx']   
    vec.ky       = dimensions['vec_ky']    
    vec.mu       = dimensions['vec_mu']     
    vec.vpa      = dimensions['vec_vpa']  
    vec.time     = dimensions['vec_time']  
    
    # Save the vector kx like it occurs in stella
    vec.kx_stella = list(copy.deepcopy(dimensions['vec_kx']))        
           
    # Sort the modes along kx
    vec.kx = sorted(vec.kx)
    
    # Create a species vector
    vec.species = range(dim.species)
    return

#-------------------------------------
def get_extraZVectors(self):
    """ Add the z axis of a nonlinear simulation as an attribute of the simulation object. 
    Calculate three extra z-axis: zeta-axis; poloidal turns axis; toroidal turn aixs. """
    
    # Get the dimension and vector objects
    dim = self.dim
    vec = self.vec
    
    # Get the geometry quantities 
    nfp = self.simulation.geometry.nfp
    iota = self.simulation.geometry.iota 
    nperiod = self.simulation.input.nperiod
    nfield_periods = self.simulation.geometry.nfield_periods
    
    # Initiate the attributes: the time, z and phi2 can differ for each input file
    vec.zeta = np.ones((dim.z)) * np.NaN
    vec.pol  = np.ones((dim.z)) * np.NaN
    vec.tor  = np.ones((dim.z)) * np.NaN

    # Read the length of one poloidal turn around the stellarator
    # nfieldperiods(s) = #poloidalturns * #fieldperiods / iota(s)
    nfield_periods_oneturn = nfp/iota 
    poloidal_turns = nfield_periods/nfield_periods_oneturn
    if np.isnan(poloidal_turns): poloidal_turns = nperiod/2

    # z in the netcdf file ranges from [-pi, pi], normalize this to [-0.5,0.5]
    z_normalized = (vec.z)/(2*np.pi)                          

    # Save the extra z vectors 
    vec.pol = z_normalized*poloidal_turns    
    vec.tor = vec.pol/iota    
    
    # For miller coordinates
    if nfield_periods<0:
        
        # z in the netcdf file ranges from [-(nperiod+1)*pi, (nperiod+1)*pi], normalize this to [-(nperiod+1), (nperiod+1)]
        z_normalized = (vec.z)/(2*np.pi)           
        vec.pol = z_normalized 
        vec.tor = vec.pol/iota    
    return

