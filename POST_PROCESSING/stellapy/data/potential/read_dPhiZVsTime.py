
import numpy as np   
import os, sys, h5py
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program 
from stellapy.data.potential.write_txtFileForDPhiZVsTime import write_txtFileForDPhiZVsTime

#===============================================================================
#                    READ THE DPHI(Z) DATA FROM THE TXT FILE                   #
#===============================================================================    

def read_dPhiZVsTime(path, dim, vec, nakxnaky):
    
    # Make sure the *.dphiz_vs_t file exists
    if not os.path.isfile(path.dphiz_vs_t) and os.path.isfile(path.output_stella): 
        write_txtFileForDPhiZVsTime(path.folder)
        
    # Read dphiz(t) from multiple simulations 
    if path.multiple_input_files: 
        return read_dPhiZVsTimeFromMultipleFiles(dim, vec, path.paths) 
    
    # Read the *.dphiz_vs_t file  
    if os.path.isfile(path.dphiz_vs_t):  
        if nakxnaky>1: return read_dPhiZVsTimeFromH5File(path.dphiz_vs_t)  
        if nakxnaky==1: return read_dPhiZVsTimeFromTxtFile(path.dphiz_vs_t)  
    
    # Critical error if we didn't find any data
    exit_reason = "The dphi(z) data versus time can not be found for:\n"
    exit_reason += "     "+str(path.input_file)
    exit_program(exit_reason, read_dPhiZVsTime, sys._getframe().f_lineno)   
    return 

#-----------------------------
def read_dPhiZVsTimeFromH5File(path): 
    with h5py.File(path, 'r') as f:  
        vec_time = f["vec_time"][()]
        dphiz_vs_tkxky = f["dphiz_vs_tkxky"][()]  
    return {"time_vs_tkxky" : vec_time, "dphiz_vs_tkxky" : dphiz_vs_tkxky} 

#-----------------------------
def read_dPhiZVsTimeFromTxtFile(path):  
    data = np.loadtxt(path, skiprows=1, dtype='float').reshape(-1, 2)
    vec_time = data[:,0]; dphiz_vs_tkxky = (data[:,1])[:,np.newaxis,np.newaxis] 
    return {"time" : vec_time, "time_vs_tkxky" : vec_time, "dphiz_vs_tkxky" : dphiz_vs_tkxky}

#-----------------------------
def read_dPhiZVsTimeFromMultipleFiles(dim, vec, paths): 
        
    # Create matrices (kx,ky)
    time_vs_tkxky = np.ones((0, dim.kx, dim.ky))*np.nan
    dphiz_vs_tkxky = np.ones((0, dim.kx, dim.ky))*np.nan 
    
    # Iterate over the simulations (1 mode per simulation)
    for path in paths:
        
        # Don't read the dummy simulation, we need to read each input file inside of it
        if "_dummy.in" not in str(path.input_file):
        
            # Read dphiz(t) from a txt or an h5 file
            if os.path.isfile(path.dphiz_vs_t):  
                if path.nakxnaky==1: data = read_dPhiZVsTimeFromTxtFile(path.dphiz_vs_t)
                if path.nakxnaky>1: data = read_dPhiZVsTimeFromH5File(path.dphiz_vs_t)   
                    
                # Make sure we have enough time points 
                dim_time = len(data["time_vs_tkxky"]); dim_time_old = np.shape(time_vs_tkxky)[0]
                if dim_time_old < dim_time:
                    time_vs_tkxky = np.append(time_vs_tkxky, np.ones((dim_time-dim_time_old, dim.kx, dim.ky))*np.nan, axis=0)
                    dphiz_vs_tkxky = np.append(dphiz_vs_tkxky, np.ones((dim_time-dim_time_old, dim.kx, dim.ky))*np.nan, axis=0) 
    
                # Iterate over the modes (kx,ky) 
                for iikx, kx in enumerate(path.vec_kx):
                    for iiky, ky in enumerate(path.vec_ky): 
                        
                        # Get the mode (kx,ky)   
                        ikx = list(vec.kx).index(kx) 
                        iky = list(vec.ky).index(ky)
                            
                        # Put the omega data in the matrices
                        time_vs_tkxky[:dim_time,ikx,iky] = data["time_vs_tkxky"] 
                        dphiz_vs_tkxky[:dim_time,ikx,iky] = data["dphiz_vs_tkxky"][:,iikx,iiky]  
     
    return {"time_vs_tkxky" : time_vs_tkxky, "dphiz_vs_tkxky" : dphiz_vs_tkxky}
            
#===============================================================================
#                ATTACH THE DPHI(Z) DATA TO THE SIMULATION OBJECT              #
#===============================================================================

def get_dPhiZVsTime(self): 
    
    # Read the data in the fluxes file
    potential = read_dPhiZVsTime(self.path, self.dim, self.vec, self.nakxnaky)  
    
    # Save the distribution data  
    self.dphiz_vs_tkxky = Data(["dphiz", "t", "kx", "ky"], potential["dphiz_vs_tkxky"], potential["time_vs_tkxky"], self.vec.kx, self.vec.ky)
    return 