# """
# 
# #===============================================================================
# #                      GET DATA FOR THE MOST UNSTABLE MODE                     #
# #===============================================================================
# 
# For <reseach>, the {ky, gamma, omega} data is gathered of the most unstable
# for each gamma(ky) spectrum of each <simulation>. The data is sorted per
# <experiment>. 
# 
# Arguments
# ---------
#     x_quantity : {fprim, tiprim, teprim, rho, ...} 
# 
# Hanne Thienpondt
# 19/10/2022
# 
# """
# 
# #!/usr/bin/python3
# import sys, os
# import numpy as np
# 
# # Personal modules
# sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)    
# from stellapy.data.lineardata.get_filterForSelectedModes import get_filterForSelectedModes
# 
# #===============================================================================
# #                      GET DATA FOR THE MOST UNSTABLE MODE                     #
# #===============================================================================
# 
# def get_gammaOfMostUnstableMode(research, modes_id="unstable", 
#         kx_range=[-9999,9999], ky_range=[-9999,9999]):
#         
#     # Gather gamma_max(parameter) or omega_max(parameter) with gamma the most unstable mode of gamma(ky) 
#     for experiment in research.experiments:
#         for simulation in experiment.simulations: 
#             
#             # Get the modes
#             vec_kx = np.array(simulation.vec.kx); dim_kx = simulation.dim.kx
#             vec_ky = np.array(simulation.vec.ky); dim_ky = simulation.dim.ky
#             
#             # Get the modes which are selected (usually the unstable modes only)
#             selected_modes = get_filterForSelectedModes(simulation, modes_id, kx_range, ky_range)
#             
#             # Get the unstable modes that lie within kx_range and ky_range
#             gamma = simulation.lineardata.gamma_avg_vs_kxky[:,selected_modes] 
#             omega = simulation.lineardata.omega_avg_vs_kxky[:,selected_modes] 
# 
#             
#             # For each mode get (ky, gamma, omega)
#             ky = [mode.ky for mode in modes]  
#             gamma = [mode.lineardata.gamma_avg for mode in modes]
#             omega = [mode.lineardata.omega_avg for mode in modes]
#             gamma_error = [[mode.lineardata.gamma_min, mode.lineardata.gamma_max] for mode in modes]
#             omega_error = [[mode.lineardata.omega_min, mode.lineardata.omega_max] for mode in modes]
#             
#             # Save the values for the most unstable mode
#             parameters[experiment.id].append(simulation.modes[0].input.inputParameters[knob][key])
#             index_most_unstable_mode = np.argmax(gamma)
#             ky_max[experiment.id].append(ky[index_most_unstable_mode])
#             gamma_max[experiment.id].append(gamma[index_most_unstable_mode])
#             omega_max[experiment.id].append(omega[index_most_unstable_mode]) 
#             gamma_error_max[experiment.id].append(gamma_error[index_most_unstable_mode]) 
#             omega_error_max[experiment.id].append(omega_error[index_most_unstable_mode]) 
#             
#         # For each experiment sort the data on <parameter>
#         sorted_indexes = list(np.array(parameters[experiment.id]).argsort())
#         parameters[experiment.id] = [parameters[experiment.id][i] for i in sorted_indexes]  
#         ky_max[experiment.id] = [ky_max[experiment.id][i] for i in sorted_indexes] 
#         gamma_max[experiment.id] = [gamma_max[experiment.id][i] for i in sorted_indexes] 
#         omega_max[experiment.id] = [omega_max[experiment.id][i] for i in sorted_indexes] 
#         gamma_error_max[experiment.id] = [gamma_error_max[experiment.id][i] for i in sorted_indexes] 
#         omega_error_max[experiment.id] = [omega_error_max[experiment.id][i] for i in sorted_indexes] 
#             
#     return gamma_max, omega_max, gamma_error_max, omega_error_max, ky_max, parameters
# 


