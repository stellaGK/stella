"""

#===============================================================================
#                              Plot dphiz(time)                                #
#=============================================================================== 

Plot dphiz(z) versus time for each (kx,ky) mode in the <research> for <folder>. 
In order to keep the plot clean, this can only be done for a single <simulation>. 

The difference in the shape of phi(z) in time is the best measure to judge the 
convergence of linear simulations. Only when this quantity reaches numerical
convergence in time (around 1e-19) can we be sure that the mode has completely
converged, if it does not reach this, there is always a possibility that the mode
looks converged, but at some point it jumps to a different (gamma, omega). These
modes are called "jumpers" due to their evolution of gamma(t) and occur in all 
devices. It is the result of initializing the mode along the field line as a 
Maxwellian, centered at zeta=0, while the most unstable mode has peaks elsewhere, 
thus first the mode converges with a peak at zeta=0, but later it jumps to a mode 
structure with peaks at zeta!=0. 

Arguments
--------- 
    folder : pathlib.Path directory where the command has been executed 
    kx_range : [-9999, 9999]
    ky_range : [-9999, 9999] 

Hanne Thienpondt
19/12/2022

"""

#!/usr/bin/python3
import sys, os
import pathlib
import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt  

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.data.lineardata.load_mostUnstableModeObject import load_mostUnstableModeObject
from stellapy.data.lineardata.get_filterForSelectedModes import get_filterForSelectedModes
from stellapy.plot.utils.labels.standardLabels import standardLabels
from stellapy.utils.files.get_filesInFolder import get_filesInFolder    
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.plot.utils.style.create_figure import create_figure  
from stellapy.simulations.Research import create_research  
from stellapy.utils.commandprompt.bash import Bash
from stellapy.plot.utils.style import Axis, Plot

#===============================================================================
#                              Plot dphiz(time)                                #
#=============================================================================== 

def plot_dphiz_vs_time(
        folder, 
        # Selected modes
        modes_id="all", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999], 
        # Research options
        folderIsExperiment=False, 
        ignore_resolution=True,  
        # Plotting options
        number_of_decimal_places=2, 
        verbose=False): 
    
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder, ignore_resolution=ignore_resolution, folderIsExperiment=folderIsExperiment)   
    simulations = [s for experiment in research.experiments for s in experiment.simulations] 

    # Create a figure for each simulation
    for simulation in simulations: 
        
        # Create a figure
        title = "Check the convergence of the linear simulations"
        if len(simulations)>1: title += ": "+simulation.line_label
        if "k_y" in title: title = title.split("$k_y")[0]+";".join(title.split("$k_y")[-1].split(";")[1:])
        ax = create_figure(title=title)
        
        # Plot dphiz(t)
        subplot_dphiz_vs_time(ax, simulation, modes_id=modes_id, kx_range=kx_range, ky_range=ky_range, 
            verbose=verbose, number_of_decimal_places=number_of_decimal_places)
        
        # Change the appearance
        ax.yaxis.labelpad = 15
        
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0: plt.show()
    return
 
#--------------------------
def subplot_dphiz_vs_time(ax, simulation, 
        modes_id="all", kx_range=[-9999,9999], ky_range=[-9999,9999], 
        fontsize=20, verbose=True, print_all_labels=True, number_of_decimal_places=2):
    
    # Automate the axis limits  
    axis = Axis(ax, Plot(),  xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.1, logy=True) 
    
    # Print information of the convergence
    if verbose:
        print("\n  |"+"{0:^15}".format("ky") + "{0:^15}".format("Tend") + "{0:^15}".format("Tconverged") +"  |")
        print("  |"+"{0:^15}".format("--") + "{0:^15}".format("----")+ "{0:^15}".format("----------") +"  |") 
    
    # Color map for modes    
    selected_modes = get_filterForSelectedModes(simulation, modes_id, kx_range, ky_range) 
    number_of_modes = np.sum(selected_modes.astype(int), axis=(0,1)) 
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,number_of_modes)); plot_i=0
    
    # Calculate the most unstable mode
    load_mostUnstableModeObject(simulation, modes_id, kx_range, ky_range)
            
    # Iterate over the modes
    for ikx, kx in enumerate(simulation.vec.kx):
        for iky, ky in enumerate(simulation.vec.ky):
            
            # Check whether it needs to be plotted
            if selected_modes[ikx, iky]==False: continue
            
            # String of ky 
            ky = simulation.vec.ky[iky] 
            string_ky = "{:.2f}".format(ky) if round(ky,2)==float(ky) else ("{:."+str(number_of_decimal_places)+"f}").format(ky)
            
            # Plot the most unstable mode in black
            if ky_range[0]==ky_range[1]:
                if kx!=simulation.mostunstablemode.kx: color = colors[plot_i]; lw=2; label=""
                if kx==simulation.mostunstablemode.kx: color = "black"; lw=4; label="Most unstable mode: $k_y\\rho_i = "+string_ky+"$"
            if ky_range[0]!=ky_range[1]:
                if ky!=simulation.mostunstablemode.ky: color = colors[plot_i]; lw=2; label=""
                if ky==simulation.mostunstablemode.ky: color = "black"; lw=4; label="Most unstable mode: $k_y\\rho_i = "+string_ky+"$"
            
            # Line label
            if print_all_labels and label !="":
                label = label.replace("Most unstable mode: ","")
            if print_all_labels and label=="":
                label = "$k_y\\rho_i = "+ string_ky + "$"
                
            # Get the dphiz(t) data
            dphiz = simulation.potential.dphiz_vs_tkxky.dphiz[:,ikx,iky]
            
            # The time axis is the same for all modes for a <range> simulation 
            # and it is unique for each (kx,ky) if we launched one mode per simulation
            try: time = simulation.potential.dphiz_vs_tkxky.t[:,ikx,iky]
            except: time = simulation.potential.dphiz_vs_tkxky.t
            
            # Remove nans and zeros (since we have a log y-scale)
            time[dphiz==0] = np.nan
            dphiz = dphiz[~np.isnan(time)]
            time = time[~np.isnan(time)] 
            
            # Plot dphiz(t) 
            ax.plot(time, dphiz, lw=lw, color=color, alpha=0.5)
            axis.update_axisLimits(time, dphiz); plot_i+=1    
            
            # Also plot the moving average to reduce the noise
            if len(time)>1:
                if len(time)<10: continue
                index = [i for i in range(len(time)) if time[i]>10.0][0]
                time = moving_average(time, index)
                dphiz = moving_average(dphiz, index)
                ax.plot(time, dphiz, lw=lw, color=color, label=label)
                 
                # Find the time at which the simulation was saturated 
                small_dphiz = [ False if dphiz[i]>1.E-16 else True for i in range(len(dphiz))] 
                if np.sum(small_dphiz)==0: small_dphiz = [ False if dphiz[i]>2.E-16 else True for i in range(len(dphiz))] 
                if np.sum(small_dphiz)==0: small_dphiz = [ False if dphiz[i]>3.E-16 else True for i in range(len(dphiz))] 
                if np.sum(small_dphiz)==0: small_dphiz = [ False if dphiz[i]>4.E-16 else True for i in range(len(dphiz))] 
                if np.sum(small_dphiz)==0: small_dphiz = [ False if dphiz[i]>5.E-16 else True for i in range(len(dphiz))] 
                if np.sum(small_dphiz)==0: small_dphiz = [ False if dphiz[i]>1.E-15 else True for i in range(len(dphiz))] 
                if np.sum(small_dphiz)>0:   
                    index_converged = [ i for i in range(len(small_dphiz)) if small_dphiz[i]==True][0] 
                    if verbose: print("  |"+"{0:^15}".format(ky)+"{0:^15}".format(round(time[-1]))+"{0:^15}".format(round(time[index_converged]))+"  |") 
                    if verbose: plt.axvline(x=time[index_converged], lw=lw, color=color) 
                else:
                    if verbose: plt.axvline(x=0, lw=lw, color=color)   
     
    # Only continue if we plotted data
    if len(ax.get_lines())!=0:  
        
        # Labels
        ax.set_xlabel(standardLabels["normalized"]["t"])
        ax.set_ylabel("max$(\\varphi_\\mathbf{k}(z,t) - \\varphi_\\mathbf{k}(z,t_{last}))$")
        
        # Limits 
        axis.rescale_axis()
        
        # Legend 
        ax.legend(labelspacing=0.1, handlelength=1, prop={'size':fontsize})
        
    return


#--------------------------
def find_input_files(folder, kymax=9999):
    
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    
    # Make sure there is a ".dtX.dphiz_vs_t" file
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix(".dt1.dphiz_vs_t")) or os.path.isfile(i.with_suffix(".dt0.1.dphiz_vs_t"))]
    
    # Find the (kx,ky) values
    kx_modes = []; ky_modes = []
    for input_file in input_files:
        kx_modes.append(float(input_file.name.split("_kx")[-1].split(".in")[0]) if "_kx" in input_file.name else 0)
        ky_modes.append(float(input_file.name.split("_ky")[-1].split("_kx")[0]) if "_kx" in input_file.name else float(input_file.name.split("_ky")[-1].split(".in")[0]))

    # Apply the kymax limit
    input_files = [input_files[i] for i in range(len(input_files)) if ky_modes[i]<=kymax] 
    return input_files
        
#--------------------------
def moving_average(x, N):
    return np.convolve(x, np.ones(N) / float(N), 'valid') 

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface
    bash = Bash(plot_dphiz_vs_time, __doc__)   
      
    # Adjust the range of wavenumbers
    bash.add_optionheader("wavenumbers")
    bash.add_option('kxmin', 'float', '', '', 'Minimum kx.')   
    bash.add_option('kxmax', 'float', '', '', 'Maximum kx.')  
    bash.add_option('kymin', 'float', '', '', 'Minimum ky.')   
    bash.add_option('kymax', 'float', 'k', '', 'Maximum ky.') 
    bash.add_optionspace()
    
    # Choose the modes
    bash.add_toggleheader("modes")
    bash.add_toggle('modes_id', 'all', '', '', 'Plot all the modes.')  
    bash.add_toggle('modes_id', 'stable', '', '', 'Plot the stable modes or unconverged modes.')  
    bash.add_toggle('modes_id', 'unstable', '', '', 'Plot the unstable converged modes. (DEFAULT)')   
    bash.add_togglespace()  
    
    # Other options  
    bash.add_optionheader("other")
    bash.add_toggleheader("other")
    
    # Choose whether we plot the stable, unstable or all modes
    bash.add_option('modes_id', 'str', 'm', '', 'Choose {unstable, stable, all}.')
    
    # Research options 
    bash.add_toggle('ignore_resolution', False, 'i', 'include_resolution', 'Include the resolution (delta t, nzed, ...) between simulations.')   
    bash.add_toggle('folderIsExperiment', False, 'f', 'folderIsExperiment', 'Create an experiment for each folder')
    
    # Plotting options 
    bash.add_option('number_of_decimal_places', 'int', 'd', 'decimal_places', 'Choose the number of decimal places in the ky label.')  
    
    # Get the bash arguments and execute the script
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_dphiz_vs_time(**args)   
