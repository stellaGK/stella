
import os, sys
import numpy as np  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.GUI.widgets.Progress import show_progressWhenReadingFiles
from stellapy.data.input.read_inputFile import read_numberOfModesFromInputFile
 
#===============================================================================
#                            READ THE OMEGA FILE
#===============================================================================
        
def read_omegaFile(path, ikx=0, iky=0):
    ''' Read the *.omega file: [time,  ky,  kx,  Re(om),  Im(om),  Re(omavg), Im(omavg)]'''
    
    # Read the reduced omega file   
    if os.path.isfile(path.omega):    
        dim_kx, dim_ky = read_numberOfModesFromInputFile(path.input_file)
        multiple_modes_per_file = True if (dim_kx*dim_ky>1) else False
        if not multiple_modes_per_file:
            omega_data = np.loadtxt(path.omega, dtype='float').reshape(-1, 3)[:,:]    
            return {"time" : omega_data[:,0], "omega" : omega_data[:,1], "gamma" : omega_data[:,2]} 
        if multiple_modes_per_file:
            omega_data = np.loadtxt(path.omega, dtype='float').reshape(dim_kx, dim_ky, -1, 5)[ikx,iky,:,:] 
            return {"time" : omega_data[:,2], "omega" : omega_data[:,3], "gamma" : omega_data[:,4]} 
        
    # Read the full omega file
    elif os.path.isfile(path.omega_stella):  
        dim_kx, dim_ky = read_numberOfModesFromInputFile(path.input_file)
        omega_data = np.loadtxt(path.omega_stella, dtype='float').reshape(-1, dim_kx, dim_ky, 7)[:,:,:,:]   
        return {"time" : omega_data[:,ikx,iky,0], "omega" : omega_data[:,ikx,iky,3], "gamma" : omega_data[:,ikx,iky,4]} 
    
    # Critical error if we didn't find any data 
    exit_reason = "The omega data can not be found.\n" 
    exit_program(exit_reason, read_omegaFile, sys._getframe().f_lineno)   
    return
    
#===============================================================================
#                ATTACH THE OMEGA DATA TO THE SIMULATION OBJECT                #
#===============================================================================

def get_omegaAndGamma(self):
    ''' Attach omega and gamma with axis [time, kx, ky] to the simulation object.'''
    
    # Show the reading progress
    show_progressWhenReadingFiles(self, "Reading omega and gamma") 
    
    # Read the data in the omega file  
    omega_data = read_omegaFile(self.path, self.ikx, self.iky)  
     
    # Multiply omega with sign(b0) to make sure it has the correct sign
    omega_data['omega'][:] = omega_data['omega'][:]*self.sign_B  
    # Save the data
    self.omega_vs_t = Data(["omega", "t"], omega_data['omega'], omega_data['time'])
    self.gamma_vs_t = Data(["gamma", "t"], omega_data['gamma'], self.omega_vs_t.t)
    return



