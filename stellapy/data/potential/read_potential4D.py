
#!/usr/bin/python3 
import pathlib
import numpy as np
import os, sys, h5py  

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.data.potential.write_h5FileForPotential4D import write_h5FileForPotential4D 
from stellapy.data.paths.load_pathObject import get_potential4DPath 
from stellapy.utils.decorators.exit_program import exit_program 
from stellapy.data.utils import Data 

#===============================================================================
#                        READ THE 4D POTENTIAL DATA
#===============================================================================

def read_potential4D(path, dim, vec, linear, nonlinear, nakxnaky): 
    ''' Read the potential from *.potential4D or *.out.nc ''' 
               
    # Make sure that the *.dt10.potential4D file exists
    if not os.path.isfile(path.potential4D) and os.path.isfile(path.output_stella): 
        write_h5FileForPotential4D(path.folder)
        get_potential4DPath(path)  
                              
    # Read from the *.dt10.potential4D file 
    if nonlinear:
        if os.path.isfile(path.potential4D):
            return read_potentialFromH5File(path.potential4D)  
    
    # Read the *.dt1.phi_vs_t h5 file or combine the *.dt1.phi_vs_t txt files
    if linear:
        
        # Read the potential data from multiple simulations 
        if path.multiple_input_files: 
            return read_potentialFromMultipleFiles(dim, vec, path.paths) 
        
        # Read the *.phi_vs_t file  
        if os.path.isfile(path.phi_vs_t):  
            if nakxnaky>1: return read_potentialFromH5File(path.phi_vs_t)  
            if nakxnaky==1: return read_potentialFromTxtFile(path.phi_vs_t)   

    # Critical error if we didn't find any data
    exit_reason = "The potential data could not be found for:\n"
    exit_reason += "     "+str(path.input_file)
    exit_program(exit_reason, read_potential4D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_potentialFromH5File(path):
    potential = {}    
    with h5py.File(path, 'r') as f: 
        for key in ["vec_time", "phi_vs_tkxky", "phi2_vs_tkxky", "phi_vs_tkxky_zeta0"]:  
            if key in f.keys():  
                potential[key] = f[key][()]  
    return potential  

#-------------------------------------
def read_potentialFromTxtFile(path): 
    with open(path) as f:
        first_line = f.readline() 
    if "nozonal" in first_line:
        data = np.loadtxt(path,skiprows=1,dtype='float').reshape(-1, 6) 
        return {"vec_time" : data[:,0], "phi2_vs_tkxky" : (data[:,1])[:,np.newaxis,np.newaxis], "phi_vs_tkxky" : (data[:,4]+1j*data[:,5])[:,np.newaxis,np.newaxis]} 
    if "nozonal" not in first_line: 
        data = np.loadtxt(path, skiprows=1, dtype='float').reshape(-1, 4)[:,:]    
        return {"vec_time" : data[:,0], "phi2_vs_tkxky" : (data[:,1])[:,np.newaxis,np.newaxis], "phi_vs_tkxky" : (data[:,2]+1j*data[:,3])[:,np.newaxis,np.newaxis]}  
    return 

#-----------------------------
def read_potentialFromMultipleFiles(dim, vec, paths): 
        
    # Create matrices (kx,ky)
    time_vs_tkxky = np.ones((0, dim.kx, dim.ky))*np.nan
    phi2_vs_tkxky = np.ones((0, dim.kx, dim.ky))*np.nan
    phi_vs_tkxky = np.ones((0, dim.kx, dim.ky), dtype=complex)*np.nan 

    # Iterate over the simulations (1 mode per simulation)
    for path in paths:
        
        # Don't read the dummy simulation, we need to read each input file inside of it
        if "_dummy.in" not in str(path.input_file):
        
            # Read the potential data from a txt or an h5 file
            if os.path.isfile(path.phi_vs_t):  
                if path.nakxnaky==1: data = read_potentialFromTxtFile(path.phi_vs_t)
                if path.nakxnaky>1: data = read_potentialFromH5File(path.phi_vs_t)   
                        
                # Make sure we have enough time points 
                dim_time = len(data["vec_time"]); dim_time_old = np.shape(time_vs_tkxky)[0]
                if dim_time_old < dim_time:
                    time_vs_tkxky = np.append(time_vs_tkxky, np.ones((dim_time-dim_time_old, dim.kx, dim.ky))*np.nan, axis=0)
                    phi_vs_tkxky  = np.append(phi_vs_tkxky, np.ones((dim_time-dim_time_old, dim.kx, dim.ky))*np.nan, axis=0)
                    phi2_vs_tkxky = np.append(phi2_vs_tkxky, np.ones((dim_time-dim_time_old, dim.kx, dim.ky))*np.nan, axis=0)
    
                # Iterate over the modes (kx,ky) 
                for iikx, kx in enumerate(path.vec_kx):
                    for iiky, ky in enumerate(path.vec_ky): 
                        
                        # Get the mode (kx,ky)   
                        ikx = list(vec.kx).index(kx) 
                        iky = list(vec.ky).index(ky)
                            
                        # Put the omega data in the matrices
                        time_vs_tkxky[:dim_time,ikx,iky] = data["vec_time"] 
                        phi_vs_tkxky[:dim_time,ikx,iky] = data["phi_vs_tkxky"][:,iikx,iiky]
                        phi2_vs_tkxky[:dim_time,ikx,iky] = data["phi2_vs_tkxky"][:,iikx,iiky] 
 
    return {"vec_time" : time_vs_tkxky, "phi_vs_tkxky" : phi_vs_tkxky, "phi2_vs_tkxky" : phi2_vs_tkxky}
            
#===============================================================================
#                  ATTACH THE POTENTIAL TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_potential4D(self): 
    
    # Read the potential
    potential = read_potential4D(self.path, self.dim, self.vec, self.linear, self.nonlinear, self.nakxnaky)  
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
    sorted_indexes = np.argsort(self.vec.kx_stella)  
    potential["phi_vs_tkxky"] = potential["phi_vs_tkxky"][:,sorted_indexes,:]
    potential["phi2_vs_tkxky"] = potential["phi2_vs_tkxky"][:,sorted_indexes,:] 
    if "phi_vs_tkxky_zeta0" in potential:
        potential["phi_vs_tkxky_zeta0"] = potential["phi_vs_tkxky_zeta0"][:,sorted_indexes,:] 
    
    # Save the potential
    self.phi_vs_tkxky = Data(["phi","t","kx","ky"], potential["phi_vs_tkxky"], potential["vec_time"], self.vec.kx, self.vec.ky) 
    self.phi2_vs_tkxky = Data(["phi2","t","kx","ky"], potential["phi2_vs_tkxky"], self.phi_vs_tkxky.t, self.vec.kx, self.vec.ky)    
    if "phi_vs_tkxky_zeta0" in potential:
        self.phi_vs_tkxky_zeta0 = Data(["phi","t","kx","ky"], potential["phi_vs_tkxky_zeta0"], self.phi_vs_tkxky.t, self.vec.kx, self.vec.ky)    
    else:
        self.phi_vs_tkxky_zeta0 = None
    return 
