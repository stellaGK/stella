""" 

#===============================================================================
#                            Read the *.omega file                             #
#===============================================================================

Read the *.omega or the *.dt1.omega file, which is a txt file if there is only
one (kx,ky) mode in the simulation, or an h5 file if there are multiple modes. If
multiple linear simulations have the same input parameters but a different 
wavevector (kx,ky) then they are combined in a singular dummy simulation. In this 
case we read the omega data from different data files and combine them in one
output object.

Returns
-------
    {time_vs_tkxky, omega_vs_tkxky, gamma_vs_tkxky}

Hanne Thienpondt
20/01/2023

"""

#!/usr/bin/python3
import h5py
import os, sys
import pathlib
import numpy as np   

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.input.read_inputFile import read_numberOfModesFromInputFile 
from stellapy.utils.commandprompt.print_progressbar import print_progressbar
from stellapy.utils.decorators.exit_program import exit_program  
from stellapy.data.utils.Data import Data  

#===============================================================================
#                            Read the *.omega file                             #
#===============================================================================
        
def read_omegaFile(path, nakxnaky, dim, vec):
    ''' Read the *.omega file: [time,  ky,  kx,  Re(om),  Im(om),  Re(omavg), Im(omavg)]'''

    # Read the omega files from multiple simulations 
    if path.multiple_input_files: 
        return read_omegaFromMultipleFiles(dim, vec, path.paths, path) 
    
    # Read the reduced omega h5 or txt file      
    if os.path.isfile(path.omega):  
        if nakxnaky>1: return read_omegaFromH5File(path)  
        if nakxnaky==1: return read_omegaFromTxtFile(path) 
    
    # Critical error if we didn't find any data 
    exit_reason = "The omega file does not exist:\n" 
    exit_reason += path.omega+"\n" 
    exit_program(exit_reason, read_omegaFile, sys._getframe().f_lineno)   
    return

#-----------------------------
def read_omegaFromH5File(path):
    """ If multiple modes (kx,ky) have been run together, the <omega_vs_tkxky> 
    quantity located in the stella output file is saved to an h5 file, therefore
    the real part is the frequency, and the imaginary part is the growth rate. """
    
    # Read the omega data from the stellapy h5 file
    if os.path.isfile(path.omega):
        with h5py.File(path.omega, 'r') as f: 
            vec_time = f["vec_time"][()]
            omega_vs_tkxky = f["omega_vs_tkxky"][()] 
        return {"time" : vec_time, "omega_vs_tkxky" : omega_vs_tkxky[:,:,:,0], "gamma_vs_tkxky" : omega_vs_tkxky[:,:,:,1]}
        
    # Read the omega data from the stella txt file
    elif os.path.isfile(path.omega_stella): 
        return read_omegaFromStellaTxtFile(path)
        
    # Exit if neither could be found
    else: 
        exit_reason = "The omega data can not be found.\n" 
        exit_program(exit_reason, read_omegaFromH5File, sys._getframe().f_lineno)   
    return 
    
#-----------------------------
def read_omegaFromTxtFile(path): 
    """ If one mode (kx,ky) is run per simulation, than the omega data is 
    written to a txt file containing {time; omega; gamma}. """
        
    # Read the omega data from the stellapy txt file
    if os.path.isfile(path.omega):
        data = np.loadtxt(path.omega, dtype='float').reshape(-1, 3)[:,:] 
        return {"time" : data[:,0], "omega_vs_tkxky" : data[:,np.newaxis,np.newaxis,1], "gamma_vs_tkxky" : data[:,np.newaxis,np.newaxis,2]} 
        
    # Read the omega data from the stella txt file
    elif os.path.isfile(path.omega_stella):
        return read_omegaFromStellaTxtFile(path)
        
    # Exit if neither could be found
    else: 
        exit_reason = "The omega data can not be found.\n" 
        exit_program(exit_reason, read_omegaFromTxtFile, sys._getframe().f_lineno)   
    return 
    
#-----------------------------
def read_omegaFromStellaTxtFile(path): 
    dim_kx, dim_ky = read_numberOfModesFromInputFile(path.input_file)
    data = np.loadtxt(path.omega_stella, dtype='float').reshape(-1, dim_kx, dim_ky, 7)[:,:,:,:]  
    return {"time" : data[:,0,0,0], "time_vs_tkxky" : data[:,:,:,0], "omega_vs_tkxky" : data[:,:,:,3], "gamma_vs_tkxky" : data[:,:,:,4]}  

#-----------------------------
def read_omegaFromMultipleFiles(dim, vec, paths, original_path):
        
    # Create matrices (kx,ky)
    omega_vs_tkxky = np.ones((0, dim.kx, dim.ky))*np.nan
    gamma_vs_tkxky = np.ones((0, dim.kx, dim.ky))*np.nan
    time_vs_tkxky = np.ones((0, dim.kx, dim.ky))*np.nan
        
    # Show reading progress 
    count = print_progressbar(0, len(paths), prefix = '   Reading gamma(t,kx,ky):', suffix = 'Complete', length = 50)+1
    
    # Iterate over the simulations
    for path in paths: 
        
        # Don't read the dummy simulation, we need to read each input file inside of it
        if "_dummy.in" not in str(path.input_file):
        
            # Read the omega data from a txt or an h5 file
            if path.nakxnaky==1: data = read_omegaFromTxtFile(path)
            if path.nakxnaky>1: data = read_omegaFromH5File(path)  
                    
            # Make sure we have enough time points 
            dim_time = len(data["time"]); dim_time_old = np.shape(time_vs_tkxky)[0]
            if dim_time_old < dim_time:
                time_vs_tkxky = np.append(time_vs_tkxky, np.ones((dim_time-dim_time_old, dim.kx, dim.ky))*np.nan, axis=0)
                omega_vs_tkxky = np.append(omega_vs_tkxky, np.ones((dim_time-dim_time_old, dim.kx, dim.ky))*np.nan, axis=0)
                gamma_vs_tkxky = np.append(gamma_vs_tkxky, np.ones((dim_time-dim_time_old, dim.kx, dim.ky))*np.nan, axis=0)

            # Iterate over the modes (kx,ky) 
            for iikx, kx in enumerate(path.vec_kx):
                for iiky, ky in enumerate(path.vec_ky): 
                    
                    # Get the mode (kx,ky)    
                    try:
                        ikx = list(vec.kx).index(kx) 
                        iky = list(vec.ky).index(ky)
                    except:
                        exit_reason = f"The mode with (kx,ky) = ({kx},{ky}) could not be found in the dimensions file.\n" 
                        exit_reason += f"Remove the dimensions file (and perhaps all dummy files) and rewrite the stellapy files.\n" 
                        exit_reason += f"      {original_path.dimensions}.\n"
                        exit_program(exit_reason, read_omegaFromH5File, sys._getframe().f_lineno)   
                        
                        
                    # Put the omega data in the matrices
                    time_vs_tkxky[:dim_time,ikx,iky] = data["time"] 
                    omega_vs_tkxky[:dim_time,ikx,iky] = data["omega_vs_tkxky"][:,iikx,iiky] 
                    gamma_vs_tkxky[:dim_time,ikx,iky] = data["gamma_vs_tkxky"][:,iikx,iiky] 
                          
        # Show the progress
        print_progressbar(count, len(paths), prefix = '   Reading gamma(t,kx,ky):', suffix = 'Complete', length = 50); count += 1      
     
    return {"time" : time_vs_tkxky, "time_vs_tkxky" : time_vs_tkxky, "omega_vs_tkxky" : omega_vs_tkxky, "gamma_vs_tkxky" : gamma_vs_tkxky}
    
#===============================================================================
#             Attach the omega and gamma data to the <omega> object            #
#===============================================================================

def get_omegaAndGamma(self):
    ''' Attach omega and gamma with axis [time, kx, ky] to the simulation object.'''
    
    # Read the data in the omega file  
    omega_data = read_omegaFile(self.path, self.nakxnaky, self.dim, self.vec)  
     
    # Multiply omega with sign(b0) to make sure it has the correct sign
    omega_data['omega_vs_tkxky'][:] = omega_data['omega_vs_tkxky'][:]*self.sign_B  
    
    # Save the data   
    self.omega_vs_tkxky = Data(["omega", "t"], omega_data['omega_vs_tkxky'], omega_data['time'], self.vec.kx, self.vec.ky)
    self.gamma_vs_tkxky = Data(["gamma", "t"], omega_data['gamma_vs_tkxky'], self.omega_vs_tkxky.t, self.vec.kx, self.vec.ky)
    return



