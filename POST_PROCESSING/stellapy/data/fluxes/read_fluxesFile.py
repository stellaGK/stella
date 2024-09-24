    
import os, sys
import numpy as np 
from datetime import datetime
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.data.utils.show_progressWhenReadingFiles import show_progressWhenReadingFiles
 
#===============================================================================
#                             READ THE FLUXES FILE
#===============================================================================
 
def read_fluxesFile(path, dim_species, sign_B, fluxes={}):
    ''' Read the heat, momentum and particle fluxes calculated by stella.

    Return fluxes[time, qflux_0, pflux_0, vflux_0, ...] 
    The data is read from the ".fluxes" files in <folder>.
    The "wout*" file is required to give the fluxes the correct sign based on sign(B).
    '''
        
    # Get the fluxes file
    if os.path.isfile(path.fluxes): 
        date = datetime.fromtimestamp(path.fluxes.stat().st_mtime) 
        flux_data = np.loadtxt(path.fluxes, dtype='float').reshape(-1, 1+3*dim_species) 
    elif os.path.isfile(path.fluxes_stella): 
        date = datetime.fromtimestamp(path.fluxes_stella.stat().st_mtime)
        flux_data = np.loadtxt(path.fluxes_stella, dtype='float').reshape(-1, 1+3*dim_species) 
    else:  
        exit_reason = "The fluxes data can not be found.\n"
        exit_reason += "    "+str(path.fluxes_stella)
        exit_program(exit_reason, read_fluxesFile, sys._getframe().f_lineno)  

    # Store the data {time, qflux, pflux, vflux} in a dictionary  
    fluxes['time'] = flux_data[:,0]   
    for specie in range(dim_species): 
        fluxes['pflux_'+str(specie)] = flux_data[:,0*dim_species+specie+1]*sign_B
        fluxes['vflux_'+str(specie)] = flux_data[:,1*dim_species+specie+1]*sign_B
        fluxes['qflux_'+str(specie)] = flux_data[:,2*dim_species+specie+1]*sign_B  
         
    # Return fluxes[time, qflux_0, pflux_0, vflux_0, qflux_1, pflux_1, vflux_1, ...]
    return fluxes, date

#===============================================================================
#                ATTACH THE FLUXES DATA TO THE SIMULATION OBJECT               #
#===============================================================================

def get_fluxes(self):
    ''' For nonlinear simulation we need to concatentate the multiple "*.fluxes" files. '''
    
    # Read the data in the fluxes file  
    fluxes, date = read_fluxesFile(self.path, self.dim.species, self.sign_B) 
    self.dim_time = len(fluxes['time'])
    self.date = date 
    
    # Initiate the attributes
    vec_time = np.ones((self.dim_time)) * np.NaN
    pflux_vs_ts = np.ones((self.dim_time, self.dim.species)) * np.NaN
    qflux_vs_ts = np.ones((self.dim_time, self.dim.species)) * np.NaN
    vflux_vs_ts = np.ones((self.dim_time, self.dim.species)) * np.NaN 
        
    # Show the reading progress
    show_progressWhenReadingFiles(self, "Reading fluxes")
        
    # Save the fluxes data
    vec_time[:] = fluxes['time'][:]
    for i in range(self.dim.species):
        pflux_vs_ts[:,i] = fluxes['pflux_'+str(i)][:]
        qflux_vs_ts[:,i] = fluxes['qflux_'+str(i)][:]
        vflux_vs_ts[:,i] = fluxes['vflux_'+str(i)][:] 

    # Save the data 
    self.pflux_vs_ts = Data(["pflux", "t", "s"], pflux_vs_ts, vec_time, range(self.dim.species))
    self.qflux_vs_ts = Data(["qflux", "t", "s"], qflux_vs_ts, self.pflux_vs_ts.t, range(self.dim.species))
    self.vflux_vs_ts = Data(["vflux", "t", "s"], vflux_vs_ts, self.pflux_vs_ts.t, range(self.dim.species)) 
    return 
