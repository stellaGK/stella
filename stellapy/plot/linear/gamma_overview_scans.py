"""
 
#===============================================================================
#       Plot Gamma(ky); Omega(ky); Gamma(x_quantity) and Omega(x_quantity)       #
#===============================================================================
 
Based on <folder> create a <research> and create a figure with 4 subplots:
    - ax1: gamma(x_quantity) based on <subplot_gamma_vs_x_quantity>
    - ax2: omega(x_quantity) based on <subplot_gamma_vs_x_quantity --omega>
    - ax3: gamma(ky) based on <subplot_gamma_vs_wavenumber>
    - ax4: omega(ky) based on <subplot_gamma_vs_wavenumber --omega>
 
Hanne Thienpondt
16/12/2022
 
"""
 
#!/usr/bin/python3  
import pathlib
import sys, os 
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
 
# Personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.plot.utils.data.recognize_varied_parameter import recognize_varied_parameter
from stellapy.data.input.read_inputFile import read_nonlinearFullFluxSurfaceFromInputFile
from stellapy.plot.linear.gamma_vs_wavenumber import subplot_gamma_vs_wavenumber
from stellapy.plot.linear.gamma_vs_parameter import subplot_gamma_vs_parameter
from stellapy.utils.commandprompt.get_bashArguments import get_bashArguments 
from stellapy.plot.utils.style.create_figure import update_figure_style
from stellapy.utils.files.get_firstInputFile import get_firstInputFile   
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.simulations.Research import create_research 
from stellapy.plot.utils.labels import standardParameters 
from stellapy.utils.commandprompt.bash import Bash
 
#===============================================================================
#                            PLOT GAMMA(PARAMETER)                             #
#===============================================================================
 
def plot_gamma_overview_scans(
        folder, 
        # Quantities to be plotted
        x_quantity=None,  
        # Research options
        experiment_parameter1="rho", 
        experiment_parameter2="-", 
        ignore_resolution=True,
        folderIsExperiment=False, 
        # Selected modes
        modes_id="unstable", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999]): 
    
    # Try to recognize the scanned parameter
    x_quantity = recognize_varied_parameter(folder, x_quantity) 
    
    # Find the stella knob and key of the given x_quantity
    knob1 = standardParameters[experiment_parameter1]["knob"]
    key1 = standardParameters[experiment_parameter1]["key"]
    knob2 = standardParameters[experiment_parameter2]["knob"]
    key2 = standardParameters[experiment_parameter2]["key"] 
     
    # Create a <research> based on the given <folder> 
    research = create_research(folders=folder, knob1=knob1, key1=key1, knob2=knob2, key2=key2, ignore_resolution=ignore_resolution, folderIsExperiment=folderIsExperiment)  
     
    # Create a figure
    fig = plt.figure(figsize=(18, 9))
    grid_specifications = gridspec.GridSpec(2, 2)
    grid_specifications.update(top=0.95, left=0.06, right=0.95, bottom=0.1, hspace=0.25)
    ax1 = plt.subplot(grid_specifications[0])     
    ax2 = plt.subplot(grid_specifications[1])     
    ax3 = plt.subplot(grid_specifications[2])     
    ax4 = plt.subplot(grid_specifications[3])   
    update_figure_style(fig, [ax1, ax2, ax3, ax4])
     
    # Plot gamma(x_quantity) and omega(x_quantity)
    subplot_gamma_vs_parameter(ax1, research, x_quantity, "gamma", modes_id, kx_range, ky_range) 
    subplot_gamma_vs_parameter(ax2, research, x_quantity, "omega", modes_id, kx_range, ky_range)  
      
    # Plot gamma(ky) and omega(ky)
    subplot_gamma_vs_wavenumber(ax3, research, "ky", "gamma", modes_id, kx_range, ky_range) 
    subplot_gamma_vs_wavenumber(ax4, research, "ky", "omega", modes_id, kx_range, ky_range)  
     
    # Show the figure  
    mpl.rcParams["savefig.directory"] = folder   
    if len(ax1.get_lines())!=0: plt.show() 
    return  
 
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
  
if __name__ == "__main__":    
     
    # Similar scripts
    script_ffs = "python3 $STELLAPY/plot/linearFFS/gamma_overview_spectra.py"
     
    # Check which script to launch      
    input_file = get_firstInputFile(pathlib.Path(os.getcwd()))
    if input_file!=None: 
        nonlinear, full_flux_surface = read_nonlinearFullFluxSurfaceFromInputFile(input_file) 
        if full_flux_surface: os.system(script_ffs+get_bashArguments(sys.argv)); sys.exit() 
        if nonlinear: print("Plot gamma(ky) and gamma(t) makes no sense for a nonlinear simulation."); sys.exit()
     
    # Create a bash-like interface
    bash = Bash(plot_gamma_overview_scans, __doc__)    
      
    # Adjust the range of wavenumbers
    bash.add_optionheader("wavenumbers")
    bash.add_option('kxmin', 'float', '', '', 'Minimum kx.')   
    bash.add_option('kxmax', 'float', '', '', 'Maximum kx.')  
    bash.add_option('kymin', 'float', '', '', 'Minimum ky.')   
    bash.add_option('kymax', 'float', 'k', '', 'Maximum ky.') 
    bash.add_optionspace()
    
    # Choose the y-quantity
    bash.add_toggleheader("x-quantity")
    bash.add_toggle('x_quantity', 'fprim', '', '', 'Plot the frequency or growthrate vs the density gradient.')  
    bash.add_toggle('x_quantity', 'tprim', '', 'tprim', 'Plot the frequency or growthrate vs the ion temperature gradient.')  
    bash.add_toggle('x_quantity', 'tiprim', '', '', 'Plot the frequency or growthrate vs the ion temperature gradient.')  
    bash.add_toggle('x_quantity', 'teprim', '', '', 'Plot the frequency or growthrate vs the electron temperature gradient.')  
    bash.add_toggle('x_quantity', 'rho', '', '', 'Plot the frequency or growthrate vs the effective minor radius.')   
    bash.add_togglespace()   
    
    # Other options 
    bash.add_toggleheader("other")
    bash.add_optionheader("other")   
      
    # Select the x-quantity and the experiment parameter
    bash.add_option('x_quantity', 'str', 'x', '', '{\
rho, tiprim, teprim, fprim, tite, delt, nmu, nvgrid, rk, \n\
dvpa, dmu, nz, nzed, nzgrid, nx, ny, y0, kx max, ky max\n\
dkx, dky, Lx, Ly, nfield, pol.turns, nperiod, D_hyper}') 
    bash.add_option('experiment_parameter1', 'str', 'e', '', 'See list above.')    
    bash.add_option('experiment_parameter2', 'str', '', '', 'See list above.')       
     
    # Choose whether we plot the stable, unstable or all modes
    bash.add_option('modes_id', 'str', 'm', '', 'Choose {unstable, stable, all}.')
     
    # Research options 
    bash.add_toggle('ignore_resolution', False, 'i', 'include_resolution', 'Include the resolution (delta t, nzed, ...) between simulations.')   
    bash.add_toggle('folderIsExperiment', False, 'f', 'folderIsExperiment', 'Create an experiment for each folder')
     
    # Get the bash arguments and execute the script
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_gamma_overview_scans(**args)   
 
     
     
