# """
# 
# #===============================================================================
# #       Plot Gamma(ky); Omega(ky); Gamma(parameter) and Omega(parameter)       #
# #===============================================================================
# 
# Based on <folder> create a <research> and create a figure with 4 subplots:
#     - ax1: gamma(parameter) based on <subplot_gamma_vs_parameter>
#     - ax2: omega(parameter) based on <subplot_gamma_vs_parameter --omega>
#     - ax3: gamma(ky) based on <subplot_gamma_vs_wavenumber>
#     - ax4: omega(ky) based on <subplot_gamma_vs_wavenumber --omega>
# 
# Hanne Thienpondt
# 10/09/2022
# 
# """
# 
# #!/usr/bin/python3  
# import pathlib
# import sys, os 
# import matplotlib as mpl 
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# 
# # Personal modules
# sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
# from stellapy.data.input.read_inputFile import read_nonlinearFullFluxSurfaceFromInputFile
# from stellapy.plot.linearFFS.gamma_vs_wavenumber import subplot_gamma_vs_wavenumber
# from stellapy.plot.linear.gamma_vs_parameter import subplot_gamma_vs_parameter  
# from stellapy.utils.commandprompt.get_bashArguments import get_bashArguments
# from stellapy.plot.utils.style.create_figure import update_figure_style 
# from stellapy.utils.files.get_firstInputFile import get_firstInputFile  
# from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
# from stellapy.simulations.Research import create_research 
# from stellapy.plot.utils.labels import standardParameters 
# from stellapy.utils.commandprompt.bash import Bash
# 
# #===============================================================================
# #                            PLOT GAMMA(PARAMETER)                             #
# #===============================================================================
# 
# def plot_gamma_overview_specta(folder, parameter="ky max", 
#         experiment_parameter1="-", experiment_parameter2="-",
#         kx_range=[-999,999], ky_range=[-999,999], modes_id="all"): 
# 
#     # Find the stella knob and key of the given parameter
#     knob1 = standardParameters[experiment_parameter1]["knob"]
#     key1 = standardParameters[experiment_parameter1]["key"]
#     knob2 = standardParameters[experiment_parameter2]["knob"]
#     key2 = standardParameters[experiment_parameter2]["key"] 
#     
#     # Create a <research> based on the given <folder> 
#     research = create_research(folders=folder, knob1=knob1, key1=key1, knob2=knob2, key2=key2)  
#     
#     # Create a figure
#     fig = plt.figure(figsize=(18, 9))
#     grid_specifications = gridspec.GridSpec(2, 2)
#     grid_specifications.update(top=0.95, left=0.06, right=0.95, bottom=0.1, hspace=0.25)
#     ax1 = plt.subplot(grid_specifications[0])     
#     ax2 = plt.subplot(grid_specifications[1])     
#     ax3 = plt.subplot(grid_specifications[2])     
#     ax4 = plt.subplot(grid_specifications[3])   
#     update_figure_style(fig, [ax1, ax2, ax3, ax4])
#     
#     # Plot gamma(parameter) and omega(parameter)
#     subplot_gamma_vs_parameter(ax1, research, parameter, "gamma", modes_id, kx_range, ky_range) 
#     subplot_gamma_vs_parameter(ax2, research, parameter, "omega", modes_id, kx_range, ky_range)  
#      
#     # Plot gamma(ky) and omega(ky)
#     subplot_gamma_vs_wavenumber(ax3, research, "ky", "gamma", modes_id, kx_range, ky_range, highlight_dominant_gamma=False) 
#     subplot_gamma_vs_wavenumber(ax4, research, "ky", "omega", modes_id, kx_range, ky_range, highlight_dominant_gamma=False)  
#     
#     # Show the figure  
#     mpl.rcParams["savefig.directory"] = folder   
#     plt.show()  
#     return
# 
# #===============================================================================
# #                             RUN AS BASH COMMAND                              #
# #===============================================================================
#  
# if __name__ == "__main__":    
#     
#     # Similar scripts
#     script_ft = "python3 $STELLAPY/plot/linear/gamma_overview_spectra.py"
#     
#     # Check which script to launch      
#     input_file = get_firstInputFile(pathlib.Path(os.getcwd())) 
#     nonlinear, full_flux_surface = read_nonlinearFullFluxSurfaceFromInputFile(input_file) 
#     if not full_flux_surface: os.system(script_ft+get_bashArguments(sys.argv)); sys.exit() 
#     if nonlinear: print("Plot gamma(ky) and gamma(t) makes no sense for a nonlinear simulation."); sys.exit()
#     
#     # Create a bash-like interface
#     bash = Bash(plot_gamma_overview_specta, __doc__)   
#     
#     # Toggle the quantities to be plotted through --kx --omega  
#     bash.add_toggle('y_quantity', 'ky', 'Plot ky(parameter).') 
#     bash.add_toggle('y_quantity', 'omega', 'Plot omega(parameter).')    
#      
#     # Select the x-quantity and the experiment parameter
#     bash.add_option('parameter', 'str', 'p', '{\
# rho, tiprim, teprim, fprim, tite, delt, nmu, nvgrid, rk, \n\
# dvpa, dmu, nz, nzed, nzgrid, nx, ny, y0, kx max, ky max\n\
# dkx, dky, Lx, Ly, nfield, pol.turns, nperiod, D_hyper, alpha0}') 
#     bash.add_option('experiment_parameter', 'str', 'e', '{\
# rho, tiprim, teprim, fprim, tite, delt, nmu, nvgrid, rk, \n\
# dvpa, dmu, nz, nzed, nzgrid, nx, ny, y0, kx max, ky max\n\
# dkx, dky, Lx, Ly, nfield, pol.turns, nperiod, D_hyper, alpha0}')     
#     
#     # Adjust the range of wavenumbers
#     bash.add_option('kxmin', 'float', '-', 'Minimum kx.')   
#     bash.add_option('kxmax', 'float', '-', 'Maximum kx.')  
#     bash.add_option('kymin', 'float', '-', 'Minimum ky.')   
#     bash.add_option('kymax', 'float', 'k', 'Maximum ky.')   
#     
#     # Choose whether we plot the stable, unstable or all modes
#     bash.add_option('modes_id', 'str', 'm', 'Choose {unstable, stable, all}.')  
#     
#     # Get the bash arguments and execute the script
#     args = bash.get_arguments()
#     args = get_rangesKxKy(args)
#     plot_gamma_overview_specta(**args)   
# 
#     
#     
