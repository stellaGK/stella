#!/usr/bin/python3
import numpy as np

# Stellapy package     
from stellapy.utils.commandprompt.print_progressbar import print_progressbar

#===============================================================================
#                                READ THE DATA                                 #
#=============================================================================== 

def get_gammaVsParameters(research, parameters1, parameters2, knob1, key1, knob2, key2, start=0 , stop=100, count=1):
    
    # Initiate the data
    gamma = np.empty((len(parameters1), len(parameters2)))*[np.NaN]
    omega = np.empty((len(parameters1), len(parameters2)))*[np.NaN]
    ky    = np.empty((len(parameters1), len(parameters2)))*[np.NaN] 

    # Add the (0,0) point since we can't simulate it
    if 0 in parameters1 and 0 in parameters2: 
        gamma[0,0] = 0
        omega[0,0] = 0
        ky[0,0] = 0 
            
    # Iterate over the simulations
    for experiment in research.experiments: 
        for simulation in experiment.simulations: 
            
            # Progress
            print_progressbar(start+(stop-start)*count/research.numberOfSimulations, 100, prefix = '   Progress:', suffix = 'Reading gamma(parameter1,parameter2).', length=50); count += 1
            
            # Get the value of the parameters for this <simulation>
            parameter1 = simulation.inputParameters[knob1][key1]
            parameter2 = simulation.inputParameters[knob2][key2]
            
            # Get the index of the parameters  
            a = parameters1.index(parameter1)
            b = parameters2.index(parameter2)
            
            # Get the linear data
            iky = [mode.ky if mode.lineardata.unstable else 0 for mode in simulation.modes]
            iomega = [mode.lineardata.omega_avg if mode.lineardata.unstable else 0 for mode in simulation.modes]
            igamma = [mode.lineardata.gamma_avg if mode.lineardata.unstable else 0 for mode in simulation.modes]
              
            # Find the most unstable mode
            index = np.nanargmax(igamma) 
            gamma[a,b] = igamma[index]
            omega[a,b] = iomega[index]
            ky[a,b] = iky[index]
            
    return gamma, omega, ky  

