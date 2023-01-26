
import h5py
import os, sys
import numpy as np  
from stellapy.data.utils import Data 
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.data.potential.write_txtFileForPotentialVsTime import write_txtFileForPotentialVsTime
from stellapy.data.output.read_outFile import read_outFile
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile,\
    read_numberOfModesFromInputFile

#===============================================================================
#                    READ THE DPHI(Z) DATA FROM THE TXT FILE                   #
#===============================================================================    

def read_potentialVsTime(path, ikx=0, iky=0): 
    
    # Check whether we have a linear or nonlinear simulation
    nonlinear = read_linearNonlinearFromInputFile(path.input_file)[1]
      
    # Make sure the *.phi_vs_t file exists
    if not os.path.isfile(path.phi_vs_t) and os.path.isfile(path.output_stella):  
        write_txtFileForPotentialVsTime(path.folder) 

    # Read the *.phi_vs_t file
    if os.path.isfile(path.phi_vs_t): 
         
        if nonlinear:
            data = np.loadtxt(path.phi_vs_t,skiprows=1,dtype='float').reshape(-1, 6)
            potential_data = {} 
            potential_data["vec_time"]  = data[:,0]
            potential_data["phi2_vs_t"] = data[:,1]
            potential_data["phi2_vs_t_zonal"] = data[:,2]
            potential_data["phi2_vs_t_nozonal"] = data[:,3]
            potential_data["phi_vs_t"]  = data[:,4] + 1j*data[:,5]
            return potential_data
        
        if not nonlinear:
            dim_kx, dim_ky = read_numberOfModesFromInputFile(path.input_file)
            multiple_modes_per_file = True if (dim_kx*dim_ky>1) else False
            with open(path.phi_vs_t) as f:
                first_line = f.readline() 
            if "nozonal" in first_line:
                data = np.loadtxt(path.phi_vs_t,skiprows=1,dtype='float').reshape(-1, 6) 
                return {"vec_time" : data[:,0], "phi2_vs_t" : data[:,1], "phi_vs_t" : data[:,4]+1j*data[:,5]} 
            if "nozonal" not in first_line:
                if not multiple_modes_per_file:
                    data = np.loadtxt(path.phi_vs_t, skiprows=1, dtype='float').reshape(-1, 4)[:,:]    
                    return {"vec_time" : data[:,0], "phi2_vs_t" : data[:,1], "phi_vs_t" : data[:,2]+1j*data[:,3]} 
                if multiple_modes_per_file:  
                    data = np.loadtxt(path.phi_vs_t, skiprows=1, dtype='float').reshape(dim_kx, dim_ky, -1, 6)[ikx,iky,:,:] 
                    return {"vec_time" : data[:,2], "phi2_vs_t" : data[:,3], "phi_vs_t" : data[:,4]+1j*data[:,5]} 
                                
    # This file doen't exist for old simulations
    elif os.path.isfile(path.output):
        potential_data = read_fromH5File(path)
        return potential_data 

    # Critical error if we didn't find any data
    exit_program("The potential data versus time can not be found.", read_potentialVsTime, sys._getframe().f_lineno) 
    return 

#===============================================================================
#                            READ FROM OLD H5 FILE                             #
#===============================================================================

def read_fromH5File(path, potential_data={}): 
    with h5py.File(path.output, 'r') as f:    
        potential_data["vec_time"] = f["vec_time"][()]
        potential_data["phi2_vs_t"] = f["phi2"][()] if ("phi2" in f.keys()) else f["phi2_vs_t"][()]
        phi2_vs_tkxky = f["phi2_vs_kxky"][()] if ("phi2_vs_kxky" in f.keys()) else f["phi2_vs_tkxky"][()]  
        phi_vs_tri = f["phi_vs_tri"][()]   
        potential_data["phi_vs_t"] = phi_vs_tri[:,0]+1j*phi_vs_tri[:,1]
        potential_data["phi2_vs_t_zonal"] = np.sum(phi2_vs_tkxky[:,:,0], axis=1)
        potential_data["phi2_vs_t_nozonal"] = 2*np.sum(phi2_vs_tkxky[:,:,1:], axis=(1,2)) 
    return potential_data

#===============================================================================
#                ATTACH THE DPHI(Z) DATA TO THE SIMULATION OBJECT              #
#===============================================================================

def get_potentialVsTime(self): 
    
    # Read the data in the fluxes file 
    potential = read_potentialVsTime(self.path, self.ikx, self.iky)  
    
    # Save the distribution data 
    self.phi_vs_t = Data(["phi", "t"], potential["phi_vs_t"], potential["vec_time"]) 
    self.phi2_vs_t = Data(["phi2", "t"], potential["phi2_vs_t"], self.phi_vs_t.t)  
    if self.nonlinear:
        self.phi2_vs_t_zonal = Data(["phi2_zonal", "t"], potential["phi2_vs_t_zonal"], self.phi_vs_t.t)  
        self.phi2_vs_t_nozonal = Data(["phi2_nozonal", "t"], potential["phi2_vs_t_nozonal"], self.phi_vs_t.t) 
    if self.linear:
        self.phi2_vs_t_zonal = np.nan
        self.phi2_vs_t_nozonal = np.nan
    return 