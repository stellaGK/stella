 
import os, sys
import copy, h5py 
import numpy as np
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.data.input.read_inputFile import read_nspecFromInputFile
from stellapy.utils.decorators.printWhichFileWeRead import printWhichFileWeRead
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables
from stellapy.data.dimensions.write_h5FileForDimensions import write_h5FileForDimensions 

#===============================================================================
#                        ATTACH THE DIMENSIONS DATA
#===============================================================================

def read_dimensions(path):
    ''' Read the dimensions from *.dimensions or *.out.nc '''
     
    # Try to write the *.dimensions file  
    if not os.path.isfile(path.dimensions):
        write_h5FileForDimensions(path.folder)
                                     
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
    exit_reason = "The dimension data could not be found."
    exit_program(exit_reason, read_dimensions, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_dimensionsFromDimensionsFile(path):  
    printWhichFileWeRead("READ DIMENSIONS FROM THE DIMENSION FILE")
    
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
def read_dimensionsFromH5File(path): 
    printWhichFileWeRead("READ DIMENSIONS FROM THE OUT.H5 FILE")
    
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
    printWhichFileWeRead("READ DIMENSIONS FROM THE OUT.NC FILE")
    
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
    
    # For a linear simulation, gather (kx,ky) of all the modes
    if self.object=="Simulation" and self.simulation.linear:
        vec.kx = sorted(list(set([mode.kx for mode in self.simulation.modes])))
        vec.ky = sorted(list(set([mode.ky for mode in self.simulation.modes])))
        dim.kx = len(vec.kx)
        dim.ky = len(vec.ky)
    
    # Save the vector kx like it accors in stella
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
    if self.object=="Simulation": self_object = self.simulation
    if self.object=="Mode": self_object = self.mode
    nfp = self_object.geometry.nfp
    iota = self_object.geometry.iota
    zeta = self_object.geometry.zeta
    nperiod = self_object.input.nperiod
    nfield_periods = self_object.geometry.nfield_periods
    
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
    vec.zeta = zeta
    vec.pol = z_normalized*poloidal_turns    
    vec.tor = vec.pol/iota    
    
    # For miller coordinates
    if nfield_periods<0:
        
        # z in the netcdf file ranges from [-(nperiod+1)*pi, (nperiod+1)*pi], normalize this to [-(nperiod+1), (nperiod+1)]
        z_normalized = (vec.z)/(2*np.pi)           
        vec.pol = z_normalized 
        vec.tor = vec.pol/iota    
    return

