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

Hanne Thienpondt
01/09/2022

"""

#!/usr/bin/python3
import sys, os
import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt  

# Stellapy package
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])  
from stellapy.plot.utils.labels.standardLabels import standardLabels
from stellapy.utils.files.get_filesInFolder import get_filesInFolder    
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.utils.get_modes import get_modes 
from stellapy.simulations.Research import create_research  
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                              Plot dphiz(time)                                #
#=============================================================================== 

def plot_dphiz_vs_time(folder, kx_range=[-999,999], ky_range=[-999,999], verbose=True):  
    
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder)
    if research.numberOfSimulations>1:
        exit_reason = "Stellapy can only plot dphiz(t) of one simulation at a time.\n" 
        exit_reason += "Please go in to a subfolder which contains a single simulation." 
        exit_program(exit_reason, plot_dphiz_vs_time, sys._getframe().f_lineno)  

    # Create a figure
    title = "Check the convergence of the linear simulations"
    ax = create_figure(title=title)
    
    # Plot dphiz(t)
    subplot_dphiz_vs_time(ax, research, kx_range, ky_range, verbose)
    
    # Change the appearance
    ax.yaxis.labelpad = 15
    
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    plt.show()
    return
 
#--------------------------
def subplot_dphiz_vs_time(ax, research, kx_range=[-999, 999], ky_range=[0,999], fontsize=20, verbose=True):
    
    # Only plot this for the first simulations
    simulation = research.experiments[0].simulations[0] 
    
    # Print information of the convergence
    if verbose:
        print("\n  |"+"{0:^15}".format("ky") + "{0:^15}".format("Tend") + "{0:^15}".format("Tconverged") +"  |")
        print("  |"+"{0:^15}".format("--") + "{0:^15}".format("----")+ "{0:^15}".format("----------") +"  |")
    
    # Get the modes that need to be plotted
    simulation.plotted_modes = get_modes(simulation, kx_range, ky_range, "all")  
    
    # Define a color for each mode that is plotted
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,len(simulation.plotted_modes)))
    
    # Keep track of the axis limits
    plot_i = 0; xright = 1; ybottom = 1.E-17; ytop = 1
    
    # Find the most unstable mode    
    vec_ky = [mode.ky if mode.lineardata.unstable else 0 for mode in simulation.plotted_modes] 
    vec_kx = [mode.kx if mode.lineardata.unstable else 0 for mode in simulation.plotted_modes] 
    vec_gamma = [mode.lineardata.gamma_avg if mode.lineardata.unstable else 0 for mode in simulation.plotted_modes] 
    ky = vec_ky[np.nanargmax(vec_gamma)] 
    kx = vec_kx[np.nanargmax(vec_gamma)] 
    
    # Plot the time evolution
    for mode in simulation.plotted_modes:
        
        # Plot the most unstable mode in black
        if ky_range[0]==ky_range[1]:
            if mode.kx!=kx: color = colors[plot_i]; lw=2; label=""
            if mode.kx==kx: color = "black"; lw=4; label="Most unstable mode: $k_y\\rho_i = "+str(mode.ky)+"$"
        if ky_range[0]!=ky_range[1]:
            if mode.ky!=ky: color = colors[plot_i]; lw=2; label=""
            if mode.ky==ky: color = "black"; lw=4; label="Most unstable mode: $k_y\\rho_i = "+str(mode.ky)+"$"
        
        # Get the data and plot it
        time = mode.potential.dphiz_vs_t.t 
        dphiz = mode.potential.dphiz_vs_t.dphiz 
        ax.plot(time, dphiz, lw=lw, color=color, alpha=0.5)
        
        # Also plot the moving average to reduce the noise
        if len(time)>1:
            if len(time)<10: continue
            index = [i for i in range(len(time)) if time[i]>10.0][0]
            time = moving_average(mode.potential.dphiz_vs_t.t, index)
            dphiz = moving_average(mode.potential.dphiz_vs_t.dphiz, index)
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
                if verbose: print("  |"+"{0:^15}".format(mode.ky)+"{0:^15}".format(round(time[-1]))+"{0:^15}".format(round(time[index_converged]))+"  |") 
                plt.axvline(x=time[index_converged], lw=lw, color=color) 
            else:
                plt.axvline(x=0, lw=lw, color=color) 
            
        # Update the axis limits
        dphiz = list(mode.potential.dphiz_vs_t.dphiz[np.isfinite(mode.potential.dphiz_vs_t.dphiz)]) + [np.nan]
        xright = np.nanmax([xright, np.max(mode.potential.dphiz_vs_t.t)])
        ytop = np.nanmax([ytop, np.max(dphiz)]) 
        plot_i += 1 
        
    # Edit the axis limits 
    ax.set_xlim(left=0, right=xright)
    ax.set_ylim(bottom=ybottom, top=ytop)   
    
    # Change appearance
    ax.legend(labelspacing=0.1, handlelength=1, prop={'size':fontsize})
    ax.set_xlabel(standardLabels["normalized"]["t"])
    ax.set_ylabel("max$(\\varphi_\\mathbf{k}(z,t) - \\varphi_\\mathbf{k}(z,t_{last}))$")
    ax.set_yscale('log')
    print()
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
    bash.add_option('kxmin', 'float', '-', 'Minimum kx.')   
    bash.add_option('kxmax', 'float', '-', 'Maximum kx.')  
    bash.add_option('kymin', 'float', '-', 'Minimum ky.')   
    bash.add_option('kymax', 'float', 'k', 'Maximum ky.')
    
    # Get the bash arguments and execute the script
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_dphiz_vs_time(**args)   
