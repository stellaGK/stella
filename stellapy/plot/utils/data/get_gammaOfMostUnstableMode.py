"""

#===============================================================================
#                      GET DATA FOR THE MOST UNSTABLE MODE                     #
#===============================================================================

For <reseach>, the {ky, gamma, omega} data is gathered of the most unstable
for each gamma(ky) spectrum of each <simulation>. The data is sorted per
<experiment>. 

Arguments
---------
    x_quantity : {fprim, tiprim, teprim, rho, ...} 

Hanne Thienpondt
01/09/2022

"""

#!/usr/bin/python3
import sys, os
import numpy as np

# Personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])    
from stellapy.plot.utils.labels import standardParameters

#===============================================================================
#                      GET DATA FOR THE MOST UNSTABLE MODE                     #
#===============================================================================

def get_gammaOfMostUnstableMode(research, x_quantity="tiprim", kx_range=[-999,999], ky_range=[-999,999]):
    
    # Identify the stella knob and key of the chosen parameter <x_quantity>
    knob = standardParameters[x_quantity]["knob"]
    key = standardParameters[x_quantity]["key"]
    
    # Initiate the data
    gamma_max = {}; omega_max = {}; ky_max = {}; parameters = {}
        
    # Gather gamma_max(parameter) or omega_max(parameter) with gamma the most unstable mode of gamma(ky) 
    for experiment in research.experiments:
        
        # For each experiments gather gamma(parameter)
        ky_max[experiment.id] = []
        gamma_max[experiment.id] = []
        omega_max[experiment.id] = []
        parameters[experiment.id] = []
        
        # Iterate over the simulations
        for simulation in experiment.simulations:
            
            # Get the unstable modes that lie within kx_range and ky_range
            modes = [mode for mode in simulation.modes if mode.lineardata.unstable]
            modes = [mode for mode in modes if (round(mode.ky,2) >= ky_range[0] and round(mode.ky,2) <= ky_range[1])]
            modes = [mode for mode in modes if (round(mode.kx,2) >= kx_range[0] and round(mode.kx,2) <= kx_range[1])]
            
            # For each mode get (ky, gamma, omega)
            ky = [mode.ky for mode in modes]  
            gamma = [mode.lineardata.gamma_avg for mode in modes]
            omega = [mode.lineardata.omega_avg for mode in modes]
            
            # Save the values for the most unstable mode
            parameters[experiment.id].append(simulation.modes[0].input.inputParameters[knob][key])
            index_most_unstable_mode = np.argmax(gamma)
            ky_max[experiment.id].append(ky[index_most_unstable_mode])
            gamma_max[experiment.id].append(gamma[index_most_unstable_mode])
            omega_max[experiment.id].append(omega[index_most_unstable_mode]) 
            
        # For each experiment sort the data on <parameter>
        sorted_indexes = list(np.array(parameters[experiment.id]).argsort())
        parameters[experiment.id] = [parameters[experiment.id][i] for i in sorted_indexes]  
        ky_max[experiment.id] = [ky_max[experiment.id][i] for i in sorted_indexes] 
        gamma_max[experiment.id] = [gamma_max[experiment.id][i] for i in sorted_indexes] 
        omega_max[experiment.id] = [omega_max[experiment.id][i] for i in sorted_indexes] 
            
    return gamma_max, omega_max, ky_max, parameters



