
import os 
import configparser 
from stellapy.data.utils.show_progressWhenReadingFiles import show_progressWhenReadingFiles
from stellapy.data.lineardata.write_linearDataPerMode import write_linearDataPerMode
from stellapy.data.lineardata.write_linearDataForFFS import write_linearDataForFFS

#===============================================================================
#                              GET THE LINEAR DATA                             #
#=============================================================================== 

def get_linearDataForFFS(self):
    ''' Read the linear data from the simulation file or calculate it.  
    Grab omega(t_last); gamma(t_last) and fluxes(t_last). Here t_last is the last
    time value of the simulation, assuming the simulation has converged at this point. '''

    # Check whether the file exists 
    if not os.path.isfile(self.input_file.with_suffix(".lineardata.ini")):
        show_progressWhenReadingFiles(self, "Calculating linear data") 
        write_linearDataPerMode(self=self)  
        write_linearDataForFFS(self=self) 
        return 
    
    # If the linear data file is older than omega(t), then rewrite it
    try: time_omegaFile = self.path.omega.stat().st_mtime
    except: time_omegaFile = self.path.omega_stella.stat().st_mtime
    time_linearDataFile = self.input_file.with_suffix(".lineardata.ini").stat().st_mtime
    if time_omegaFile>time_linearDataFile:
        show_progressWhenReadingFiles(self, "Calculating linear data")  
        write_linearDataPerMode(self=self)  
        write_linearDataForFFS(self=self) 
        return 
    
    # Read the "*.lineardata" file
    return read_linearDataForFFSFromLinearDataFile(self=self)  
    
#===============================================================================
#                         READ THE LINEAR DATA PER MODE                        #
#===============================================================================  
 
def read_linearDataForFFSFromLinearDataFile(omega=None, dim=None, vec=None, input_file=None, 
        full_flux_surface=None, percentage=None, self=None):
    
    # Read the linear data from "*.lineardata"  
    file = configparser.ConfigParser() 
    path = input_file.with_suffix(".lineardata.ini")
    file.read(path)        
    
    # Make sure the data is written
    if "full flux surface" not in file:  
        if self!=None:
            full_flux_surface = self.input.full_flux_surface
            input_file = self.input_file
            percentage = self.percentage
            dim = self.dim; vec = self.vec
            omega = self.omega
        write_linearDataPerMode(omega, dim, vec, input_file, full_flux_surface, file, percentage)
    
    # Read the full flux surface data
    gamma = float(file["full flux surface"]["gamma"]) 
    omega = float(file["full flux surface"]["omega"]) 
    gamma_min = float(file["full flux surface"]["gamma_min"]) 
    gamma_max = float(file["full flux surface"]["gamma_max"]) 
    omega_min = float(file["full flux surface"]["omega_min"]) 
    omega_max = float(file["full flux surface"]["omega_max"]) 
    kymin = float(file["full flux surface"]["kymin"]) 
            
    # Return the linear data per mode
    return gamma, omega, gamma_min, gamma_max, omega_min, omega_max, kymin
 
#===============================================================================
#                              GET THE LINEAR DATA                             #
#===============================================================================

def get_linearDataPerMode(self):
    
    # Read the linear data
    gamma, omega, gamma_min, gamma_max, omega_min, omega_max, kymin = get_linearDataForFFS(self)
    
    # Save this data 
    self.ffs_kymin = kymin 
    self.ffs_gamma = gamma
    self.ffs_omega = omega
    self.ffs_gamma_min = gamma_min
    self.ffs_gamma_max = gamma_max
    self.ffs_omega_min = omega_min
    self.ffs_omega_max = omega_max
    return 
    
            