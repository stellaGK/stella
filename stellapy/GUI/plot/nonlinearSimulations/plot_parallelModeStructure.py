# 
# # Load the modules 
# import sys
# import numpy as np
# import matplotlib.pyplot as plt
# 
# # Personal modules
# from stellapy.utils.decorators.verbose import noverbose  
# 
# # Plotting modules 
# from stellapy.GUI.plot.utils.Plot import Plot
# from stellapy.GUI.plot.utils import Axis, Legend
# from stellapy.GUI.plot.utils import create_figure
# from stellapy.GUI.utils import DisplayProgressGUI   
# from stellapy.utils.config.species.display_multipleSpecies import get_speciesLineStyle 
# from stellapy.utils.decorators.exit_program import exit_program
# 
# #===============================================================================
# #                    PLOT QUANTITY(Z) FOR NONLINEAR SIMULATIONS
# #===============================================================================
# 
# @noverbose 
# def plot_parallelModeStructure(\
#             # Simulations
#                 research=None,\
#                 experiment_id="All experiments",\
#                 simulation_id="All simulations",\
#                 species=[0],\
#             # Data 
#                 x_quantity="z",\
#                 y_quantity="phi2",\
#                 units="normalized",\
#                 x_range=None,\
#                 y_range=None,\
#             # Labels
#                 x_label=None,\
#                 y_label=None,\
#                 title=None,\
#             # Figure 
#                 ax=None,\
#                 Progress=None,\
#                 show_figure=True,\
#             # Toggles
#                 log=False): 
#      
#     # Only implemented for specific data
#     if y_quantity!="phi2" and x_quantity!="z":
#         exit_program("Not implemented", plot_parallelModeStructure, sys._getframe().f_lineno)  
#     
#     # Save the plotting details to a <plot> object
#     plot = Plot()  
#     plot.update_xyData(x_quantity, y_quantity, x_range, y_range, units)  
#     plot.update_simulations(experiment_id, simulation_id, species=[0]) 
#     plot.update_figure(Progress, ax, show_figure)
#     plot.update_labels(title, x_label, y_label) 
#     plot.update_legend(fontsize=14, loc='upper left') 
#     plot.process_plottingVariables(research)  
#     
#     # Update the progress bar of the GUI 
#     displayProgressGUI = DisplayProgressGUI(research, "Plot "+plot.yname+" versus "+plot.xname)   
#     
#     # Create the figure and a <legend> and <axis> object
#     ax = create_figure(ax, plot, research) 
#     legend = Legend(ax, plot)
#     axis = Axis(ax, plot, overshoot_y=1.4)      
# 
#     # Define a color per mode
#     colors = plt.cm.get_cmap('jet')(np.linspace(0,1,research.numberOfPlottedModes))
#           
#     # Iterate over the experiments, simulations and the species 
#     for experiment in research.plotted_experiments:  
#         for simulation in experiment.plotted_simulations: 
#             for specie in species:
#  
#                 # Update the progress bar of the GUI  
#                 displayProgressGUI.move_progressBar() 
#                 
#                 # The line style depends on the species and the color on yquantity/experiments
#                 style = get_speciesLineStyle(species, specie) 
#                  
#                 # Get the (time,z,quanitity) data for the plots  
#                 vec_phi2 = simulation.potential.phi2_vs_tz.phi2
#                 vec_time = simulation.potential.phi2_vs_tz.data.t
#                 vec_z = simulation.potential.phi2_vs_tz.data.z 
#                            
#                # Keep track of the axis limits 
#                axis.update_axisLimits(vec_z, phi_vs_z)         
# 
#     # Plot the magnetic field
#     B_bmag = simulation.B_bmag-np.min(simulation.B_bmag)
#     B_bmag = B_bmag/np.max(B_bmag)*axis.ylims[1]
#     ax.plot(vec_z, B_bmag, color="red", linestyle=":")
#                      
#     # Add the legend and rescale the axis
#     legend.add_legend()
#     axis.rescale_axis()
#     if averagetimeframe: ax.set_title("Parallel mode structure of the 10 biggest modes at t = ["+"{:.1f}".format(t1)+', '+"{:.1f}".format(t2)+"]")
#     if not averagetimeframe: ax.set_title("Parallel mode structure of the 10 biggest modes at t = "+"{:.1f}".format(t1))
#  
#     # Show the figure 
#     if show_figure: plt.show()
#      
#     # Update the progress bar of the GUI  
#     displayProgressGUI.finilize_progressBar()
#     if True: return
#  
# #################################################################
# #                   GET THE (X,Y) DATA
# #################################################################
#  
# def get_zAxis(simulation, x_quantity):
#      
#     # Get the z-axis
#     if x_quantity == "z":     vec_z = simulation.vec_z 
#     if x_quantity == "zeta":  vec_z = simulation.vec_zeta
#     if x_quantity == "pol":   vec_z = simulation.vec_pol
#     if x_quantity == "tor":   vec_z = simulation.vec_tor
#     return vec_z
#  
# #----------------------------- 
# def get_potentialVersusZ(simulation, y_quantity, units, specie, t_specific, it, averagetimeframe):
#  
#     # Get the potential(t,z,kx,ky) from the netcdf file  
#     vec_time = simulation.parallelModeStructure.vec_time
#     phi_vs_tzkxkyri = simulation.parallelModeStructure.phi_vs_tzkxkyri
#     phi_vs_tzkxkyri_tavg = simulation.parallelModeStructure.phi_vs_tzkxkyri_tavg 
#  
#     # Extract the quantity we need
#     if it==None:
#         if y_quantity=="phiReal": phi_vs_tzkxky = phi_vs_tzkxkyri_tavg[:,:,:,:,0]
#         if y_quantity=="phiImag": phi_vs_tzkxky = phi_vs_tzkxkyri_tavg[:,:,:,:,1]
#         if y_quantity=="phi2":    phi_vs_tzkxky = (np.abs(phi_vs_tzkxkyri_tavg[:,:,:,:,0] + 1j*phi_vs_tzkxkyri_tavg[:,:,:,:,1]))**2
#         phi_vs_zkxky = calculate_timeAverage(simulation, phi_vs_tzkxky, specie, vec_time=vec_time, t_axis=0, t_specific=t_specific) 
#      
#     if it!=None and not averagetimeframe:
#         if y_quantity=="phiReal": phi_vs_zkxky = phi_vs_tzkxkyri[it,:,:,:,0]
#         if y_quantity=="phiImag": phi_vs_zkxky = phi_vs_tzkxkyri[it,:,:,:,1]
#         if y_quantity=="phi2":    phi_vs_zkxky = (np.abs(phi_vs_tzkxkyri[it,:,:,:,0] + 1j*phi_vs_tzkxkyri[it,:,:,:,1]))**2
#          
#     if it!=None and averagetimeframe:
#         if y_quantity=="phiReal": phi_vs_zkxky = phi_vs_tzkxkyri_tavg[it,:,:,:,0]
#         if y_quantity=="phiImag": phi_vs_zkxky = phi_vs_tzkxkyri_tavg[it,:,:,:,1]
#         if y_quantity=="phi2":    phi_vs_zkxky = (np.abs(phi_vs_tzkxkyri_tavg[it,:,:,:,0] + 1j*phi_vs_tzkxkyri_tavg[it,:,:,:,1]))**2
#          
#     # Denormalize the potential
#     if units=="SI":  phi_vs_zkxky  = phi_vs_zkxky*simulation.referenceUnits['phi'] 
#       
#     # Return the y data
#     if it!=None: t1, t2 = vec_time[it], vec_time[it+1]
#     if it==None: t1, t2 = get_timeFrame(simulation, specie)    
#     return phi_vs_zkxky, t1, t2
#    

