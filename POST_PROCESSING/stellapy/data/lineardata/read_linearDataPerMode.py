
import os, sys
import numpy as np
import configparser  
from stellapy.data.utils.show_progressWhenReadingFiles import show_progressWhenReadingFiles
from stellapy.data.lineardata.write_linearDataPerMode import write_linearDataPerMode
from stellapy.data.lineardata.get_modeIdentifier import get_modeIdentifier
from stellapy.utils.decorators.exit_program import exit_program

#===============================================================================
#                              GET THE LINEAR DATA                             #
#=============================================================================== 

def read_linearDataPerMode(self):
    ''' Read the linear data from the simulation file or calculate it. Grab 
    omega(t_last); gamma(t_last) or np.mean(omega[85%:]) and np.mean(gamma[85%:]). 
    Here t_last is the last time value of the simulation, and [85%:] is the last 15%
    of the time trace, the percentage is defined by <self.percentage> = 0.85. '''  
    
    # The "*.lineardata" file needs to be written for each <range> simulation, or for
    # each <dummy_input_file> which contains multiple <one-mode-per-simulation> input files,
    # <path.input_files> can consist of a combination of both kinds of input files.
    paths = [self.path] if not self.path.multiple_input_files else self.path.paths 
            
    # Check whether the lineardata file exists, or rewrite it if the omega file is more recent
    for path in paths: 
        
        # <one-mode-per-simulation> input files are included in a dummy input file 
        if path.dummy_input_file!=None: input_file = path.dummy_input_file 
        if path.dummy_input_file==None: input_file = path.input_file 
        
        # If the "*.lineardata" file doesn't exist, write it
        if not os.path.isfile(input_file.with_suffix(".lineardata.ini")): 
            write_linearDataPerModeForInputFile(self, input_file)
            
        # If the "*.lineardata" file exists, make sure it's older than the omega file 
        elif not "_dummy.in" in str(path.input_file):
            time_linearDataFile = input_file.with_suffix(".lineardata.ini").stat().st_mtime
            try: time_omegaFile = path.omega.stat().st_mtime
            except: time_omegaFile = path.omega_stella.stat().st_mtime
            if time_omegaFile>time_linearDataFile: 
                write_linearDataPerModeForInputFile(self, input_file)
        
    # Read the "*.lineardata" file
    show_progressWhenReadingFiles(self, "Reading linear data")  
    if self.path.multiple_input_files: 
        return read_linearDataFromMultipleLinearDataFiles(self.dim, self.vec, self.path.paths)  
    if not self.path.multiple_input_files: 
        return read_linearDataFromLinearDataFile(self=self)   
    
#----------------------------------
def write_linearDataPerModeForInputFile(self, input_file):
    show_progressWhenReadingFiles(self, "Calculating linear data") 
    from stellapy.simulations.Simulation import create_simulations 
    dummy_simulation = create_simulations(input_files=[input_file])[0]
    write_linearDataPerMode(self=dummy_simulation.lineardata); del dummy_simulation 
    
#===============================================================================
#                         READ THE LINEAR DATA PER MODE                        #
#===============================================================================  

def read_linearDataFromMultipleLinearDataFiles(dim, vec, paths):
    
    # Store (gamma, omega) in a matrices versus (kx,ky)
    gamma_last_vs_kxky = np.ones((dim.kx, dim.ky)) * np.nan
    omega_last_vs_kxky = np.ones((dim.kx, dim.ky)) * np.nan
    gamma_avg_vs_kxky = np.ones((dim.kx, dim.ky)) * np.nan
    omega_avg_vs_kxky = np.ones((dim.kx, dim.ky)) * np.nan
    gamma_min_vs_kxky = np.ones((dim.kx, dim.ky)) * np.nan
    gamma_max_vs_kxky = np.ones((dim.kx, dim.ky)) * np.nan
    omega_min_vs_kxky = np.ones((dim.kx, dim.ky)) * np.nan
    omega_max_vs_kxky = np.ones((dim.kx, dim.ky)) * np.nan
    converged_vs_kxky = np.ones((dim.kx, dim.ky)) * np.nan
    unstable_vs_kxky = np.ones((dim.kx, dim.ky)) * np.nan
    stable_vs_kxky = np.ones((dim.kx, dim.ky)) * np.nan
    
    # Iterate over the simulations
    for path in paths: 
        
        # Don't read input files which are included in a dummy input file
        if path.dummy_input_file==None:
    
            # Read (gamma, omega) of one simulation
            gamma_last_vs_kxky_temp, omega_last_vs_kxky_temp, gamma_avg_vs_kxky_temp, omega_avg_vs_kxky_temp, \
            gamma_min_vs_kxky_temp, gamma_max_vs_kxky_temp, omega_min_vs_kxky_temp, omega_max_vs_kxky_temp, \
            stable_vs_kxky_temp, unstable_vs_kxky_temp, converged_vs_kxky_temp = \
            read_linearDataFromLinearDataFile(path.input_file, path.dim_kx, path.dim_ky, path.vec_kx, path.vec_ky) 
            
            # Put the (kx,ky) data in the bigger matrices 
            for iikx, kx in enumerate(path.vec_kx):
                for iiky, ky in enumerate(path.vec_ky): 

                    # Get the mode (kx,ky)
                    # WARNING: If an error pops up here, make sure the list of input files
                    # contains the dummy input files and not only the normal input files!    
                    ikx = list(vec.kx).index(kx) 
                    iky = list(vec.ky).index(ky) 
                        
                    # Put (gamma, omega) in the matrices
                    gamma_last_vs_kxky[ikx,iky] = gamma_last_vs_kxky_temp[iikx,iiky]
                    omega_last_vs_kxky[ikx,iky] = omega_last_vs_kxky_temp[iikx,iiky]
                    gamma_avg_vs_kxky[ikx,iky] = gamma_avg_vs_kxky_temp[iikx,iiky]
                    omega_avg_vs_kxky[ikx,iky] = omega_avg_vs_kxky_temp[iikx,iiky]
                    gamma_min_vs_kxky[ikx,iky] = gamma_min_vs_kxky_temp[iikx,iiky]
                    gamma_max_vs_kxky[ikx,iky] = gamma_max_vs_kxky_temp[iikx,iiky]
                    omega_min_vs_kxky[ikx,iky] = omega_min_vs_kxky_temp[iikx,iiky]
                    omega_max_vs_kxky[ikx,iky] = omega_max_vs_kxky_temp[iikx,iiky]
                    converged_vs_kxky[ikx,iky] = converged_vs_kxky_temp[iikx,iiky]
                    unstable_vs_kxky[ikx,iky] = unstable_vs_kxky_temp[iikx,iiky]
                    stable_vs_kxky[ikx,iky] = stable_vs_kxky_temp[iikx,iiky] 
        
    # Return the linear data versus (kx,ky)
    return gamma_last_vs_kxky, omega_last_vs_kxky, gamma_avg_vs_kxky, omega_avg_vs_kxky, \
    gamma_min_vs_kxky, gamma_max_vs_kxky, omega_min_vs_kxky, omega_max_vs_kxky, \
    stable_vs_kxky, unstable_vs_kxky, converged_vs_kxky
            
#-------------------------------------
def read_linearDataFromLinearDataFile(input_file=None, dim_kx=None, dim_ky=None, vec_kx=None, vec_ky=None, self=None):
    
    # Extract the data
    if self!=None:
        input_file = self.input_file 
        dim_kx = self.dim.kx
        dim_ky = self.dim.ky
        vec_kx = self.vec.kx 
        vec_ky = self.vec.ky

    # Read the linear data from "*.lineardata"  
    file = configparser.ConfigParser() 
    path = input_file.with_suffix(".lineardata.ini")
    file.read(path)       
    
    # Store (gamma, omega) in a matrices versus (kx,ky)
    gamma_last_vs_kxky = np.ones((dim_kx, dim_ky)) * np.nan
    omega_last_vs_kxky = np.ones((dim_kx, dim_ky)) * np.nan
    gamma_avg_vs_kxky = np.ones((dim_kx, dim_ky)) * np.nan
    omega_avg_vs_kxky = np.ones((dim_kx, dim_ky)) * np.nan
    gamma_min_vs_kxky = np.ones((dim_kx, dim_ky)) * np.nan
    gamma_max_vs_kxky = np.ones((dim_kx, dim_ky)) * np.nan
    omega_min_vs_kxky = np.ones((dim_kx, dim_ky)) * np.nan
    omega_max_vs_kxky = np.ones((dim_kx, dim_ky)) * np.nan
    converged_vs_kxky = np.ones((dim_kx, dim_ky)) * np.nan
    unstable_vs_kxky = np.ones((dim_kx, dim_ky)) * np.nan
    stable_vs_kxky = np.ones((dim_kx, dim_ky)) * np.nan
    
    # Read the linear data per mode
    for ikx in range(dim_kx):
        for iky in range(dim_ky):
            
            # Get the mode identifier 
            mode_identifier, _, _ = get_modeIdentifier(vec_kx[ikx], vec_ky[iky])
            
            # Make sure it it written 
            if mode_identifier not in file:   
                from stellapy.simulations.Simulation import create_simulations
                dummy_simulation = create_simulations(input_files=[input_file])[0]
                write_linearDataPerMode(self=dummy_simulation.lineardata); del dummy_simulation
                
            # Read the data
            try:
                gamma_last_vs_kxky[ikx, iky] = float(file[mode_identifier]["gamma_last"])
                omega_last_vs_kxky[ikx, iky] = float(file[mode_identifier]["omega_last"])
                gamma_avg_vs_kxky[ikx, iky] = float(file[mode_identifier]["gamma_avg"])
                omega_avg_vs_kxky[ikx, iky] = float(file[mode_identifier]["omega_avg"])
                gamma_min_vs_kxky[ikx, iky] = float(file[mode_identifier]["gamma_min"])
                gamma_max_vs_kxky[ikx, iky] = float(file[mode_identifier]["gamma_max"])
                omega_min_vs_kxky[ikx, iky] = float(file[mode_identifier]["omega_min"])
                omega_max_vs_kxky[ikx, iky] = float(file[mode_identifier]["omega_max"])
                converged_vs_kxky[ikx, iky] = True if (file[mode_identifier]["converged"]=="True") else False
                unstable_vs_kxky[ikx, iky] = True if (file[mode_identifier]["unstable"]=="True") else False
                stable_vs_kxky[ikx, iky] = True if (file[mode_identifier]["stable"]=="True") else False
            except: 
                exit_reason = f"Could not find the linear data for mode (kx,ky)={mode_identifier}.\n"
                exit_reason += "     "+str(path)
                exit_program(exit_reason, read_linearDataFromLinearDataFile, sys._getframe().f_lineno)   
            
    # Return the linear data versus (kx,ky)
    return gamma_last_vs_kxky, omega_last_vs_kxky, gamma_avg_vs_kxky, omega_avg_vs_kxky, \
    gamma_min_vs_kxky, gamma_max_vs_kxky, omega_min_vs_kxky, omega_max_vs_kxky, \
    stable_vs_kxky, unstable_vs_kxky, converged_vs_kxky
            
#===============================================================================
#                              GET THE LINEAR DATA                             #
#===============================================================================

def get_linearDataPerMode(self):
    
    # Read the linear data 
    gamma_last_vs_kxky, omega_last_vs_kxky, gamma_avg_vs_kxky, omega_avg_vs_kxky, gamma_min_vs_kxky, gamma_max_vs_kxky, \
    omega_min_vs_kxky, omega_max_vs_kxky, stable_vs_kxky, unstable_vs_kxky, converged_vs_kxky = read_linearDataPerMode(self) 

    # Save this data  
    self.stable_vs_kxky     = stable_vs_kxky
    self.unstable_vs_kxky   = unstable_vs_kxky
    self.converged_vs_kxky  = converged_vs_kxky
    self.omega_last_vs_kxky = omega_last_vs_kxky
    self.omega_avg_vs_kxky  = omega_avg_vs_kxky
    self.omega_min_vs_kxky  = omega_min_vs_kxky
    self.omega_max_vs_kxky  = omega_max_vs_kxky
    self.gamma_last_vs_kxky = gamma_last_vs_kxky
    self.gamma_avg_vs_kxky  = gamma_avg_vs_kxky
    self.gamma_min_vs_kxky  = gamma_min_vs_kxky
    self.gamma_max_vs_kxky  = gamma_max_vs_kxky
    
    # Remove certain modes 
    for kx in self.removedKxModes:
        ikx = list(self.vec.kx).index(kx) 
        self.stable_vs_kxky[ikx,:] = np.nan
        self.unstable_vs_kxky[ikx,:] = np.nan
        self.converged_vs_kxky[ikx,:] = np.nan
        self.omega_last_vs_kxky[ikx,:] = np.nan
        self.omega_avg_vs_kxky[ikx,:] = np.nan
        self.omega_min_vs_kxky[ikx,:] = np.nan
        self.omega_max_vs_kxky[ikx,:] = np.nan
        self.gamma_last_vs_kxky[ikx,:] = np.nan
        self.gamma_avg_vs_kxky[ikx,:] = np.nan
        self.gamma_min_vs_kxky[ikx,:] = np.nan
        self.gamma_max_vs_kxky[ikx,:] = np.nan
    for ky in self.removedKyModes:
        iky = list(self.vec.ky).index(ky) 
        self.stable_vs_kxky[:,iky] = np.nan
        self.unstable_vs_kxky[:,iky] = np.nan
        self.converged_vs_kxky[:,iky] = np.nan
        self.omega_last_vs_kxky[:, iky] = np.nan
        self.omega_avg_vs_kxky[:,iky] = np.nan
        self.omega_min_vs_kxky[:,iky] = np.nan
        self.omega_max_vs_kxky[:,iky] = np.nan
        self.gamma_last_vs_kxky[:,iky] = np.nan
        self.gamma_avg_vs_kxky[:,iky] = np.nan
        self.gamma_min_vs_kxky[:,iky] = np.nan
        self.gamma_max_vs_kxky[:,iky] = np.nan
    return 
    
    
    