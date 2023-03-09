"""

#===============================================================================
#                               Plot phi2(time)                                #
#===============================================================================

For <research>, created for <folder>, the time evolution of the average potential 
squared in real space is plotted. If there is only one <simulation>, the zonal 
contributions are shown as well. If the simulations are linear, then 
plot/linear/potential_vs_time.py is called instead.


Average potential squared in real space
---------------------------------------
In stella the Fourier components hat{phi}_{kxky} are related to the convential
Fourier components through hat{phi}_{kxky} = hat{PHI}_{kxky}/NxNy. In other words,
the inverse Fourier transformation is defined as
    
    phi_{xy} = 1/NxNy sum_{kxky} hat{PHI}_{kxky} exp(i*kx*x+i*ky*y)
             = sum_{kxky} hat{phi}_{kxky} exp(i*kx*x+i*ky*y)
             
As a result, Parseval's theorem is given by

    sum_{xy} |phi_{xy}|^2 = 1/NxNy sum_{kxky} |hat{PHI}_{kxky}|^2
                          = NxNy sum_{kxky} |hat{phi}_{kxky}|^2  
                          
Therefore, the average potential squared in real space, is simply given by the sum
over the (kx,ky) Fourier components as calculated by stella. 

    |phi|^2 = 1/NxNy sum_{xy} |phi(x,y)|^2 
            = 1/(NxNy)^2 sum_{kxky} |hat{PHI}_{kxky}|^2
            = sum_{kxky} |hat{phi}_{kxky}|^2
            
Note that phi refers to the perturbed potential delta phi_1 in the stella bible.

t_range
-------
Change the time range where the simulation is considered saturated. 
    >> plot_saturatedflux_vs_parameter -t "[500, 1000]" (or -t "500, 1000" or --trange "[500, 1000]") 

zonal
-----
If none then <phi2> is plotted, if True then <phi2_zonal> is plotted, if False
then <phi2_nozonal> is plotted.

Hanne Thienpondt 
15/12/2022

"""


#!/usr/bin/python3 
import pathlib 
import sys, os
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.data.input.read_inputFile import read_nonlinearFullFluxSurfaceFromInputFile
from stellapy.utils.files.get_firstInputFile import get_firstInputFile
from stellapy.plot.utils.labels.standardLabels import standardLabels 
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.Research import create_research
from stellapy.utils.commandprompt.bash import Bash
from stellapy.plot.utils.style import Axis, Plot

#===============================================================================
#                               Plot phi2(time)                                #
#===============================================================================

def plot_potential_vs_time(
        folder, 
        # Quanities to be plotted
        y_quantity="phi2",
        # Time frame to calculate saturated quantities 
        trange=None, 
        # Research
        folderIsExperiment=False, 
        folderIsSimulation=False,
        # Plot
        xrange=None,
        yrange=None,
        label_units="_norm",
        zonal=None, 
        log=False):  
    
    # Create <simulations> based on the given <folder>
    research = create_research(folders=folder, resolutionScan=False, folderIsExperiment=folderIsExperiment, folderIsSimulation=folderIsSimulation, ignore_resolution=False) 
    simulations = [s for experiment in research.experiments for s in experiment.simulations] 
    research.update_timeFrame(trange); experiments = research.experiments
    
    # Create a figure for each experiment
    for experiment in experiments: 
     
        # Create a figure
        title = "Time evolution of the average potential squared in real space"
        if research.numberOfExperiments>1: title += ": "+experiment.line_label
        ax = create_figure(title=title)    
        
        # Plot phi2(t)
        research.experiments = [experiment]; research.numberOfSimulations = len(experiment.simulations)
        if len(simulations)==1: subplot_potential_vs_time_zonal_contributions(ax, simulations[0], y_quantity=y_quantity, log=log)
        if len(simulations)!=1: subplot_potential_vs_time(ax, research, y_quantity=y_quantity, zonal=zonal, log=log, label_units=label_units)   
        
        # Appearance
        ax.yaxis.labelpad = 15
        ax.xaxis.labelpad = 10
        if log: ax.set_yscale("log"); ax.set_ylim([10e-9, 1000])         
        if xrange: ax.set_xlim(xrange)    
        if yrange: ax.set_ylim(yrange)  
            
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0: plt.show()
    return

#----------------------------------------
def subplot_potential_vs_time_zonal_contributions(ax, simulation, y_quantity="phi2", log=False, fontsize=20):

    # Get the data
    if y_quantity!="phi2": exit_program("Only y_quantity==phi2 has been implemented.", subplot_potential_vs_time_zonal_contributions, sys._getframe().f_lineno)
    phi2 = simulation.potential.phi2_vs_t.phi2 
    phi2_zonal = simulation.potential.phi2_vs_t_zonal.phi2_zonal
    phi2_nozonal = simulation.potential.phi2_vs_t_nozonal.phi2_nozonal
    time = simulation.potential.phi2_vs_t.t
    
    # Get the label
    tstart = simulation.time.tstart; tend = simulation.time.tend
    saturated_potential = np.mean(phi2[(time>=tstart) & (phi2<=tend)]) 
    saturated_potential_zonal = np.mean(phi2_zonal[(time>=tstart) & (phi2_zonal<=tend)]) 
    saturated_potential_nozonal = np.mean(phi2_nozonal[(time>=tstart) & (phi2_nozonal<=tend)]) 
    label = "$\\langle |\\sum_{k_x,k_y} \\hat{\\varphi}_{\\text{sat}}|^2 ({\\rho_i T_i}/{ae})^2 \\rangle_z = "+str("{:.2f}".format(saturated_potential))+"$"
    label_zonal = "$\\langle |\\sum_{k_x,k_y=0} \\hat{\\varphi}_{\\text{Z,sat}}|^2 ({\\rho_i T_i}/{ae})^2 \\rangle_z = "+str("{:.2f}".format(saturated_potential_zonal))+"$"
    label_nozonal = "$\\langle |\\sum_{k_x,k_y\\neq 0} \\hat{\\varphi}_{\\text{NZ,sat}}|^2 ({\\rho_i T_i}/{ae})^2 \\rangle_z = "+str("{:.2f}".format(saturated_potential_nozonal))+"$"
 
    # Plot the data
    ax.plot(time, phi2[:], color="black")  
    ax.plot(time, phi2_zonal[:], color="navy")  
    ax.plot(time, phi2_nozonal[:], color="crimson")  
    ax.plot([tstart, tend], [saturated_potential, saturated_potential], color="black", label=label) 
    ax.plot([tstart, tend], [saturated_potential_zonal, saturated_potential_zonal], color="navy", label=label_zonal) 
    ax.plot([tstart, tend], [saturated_potential_nozonal, saturated_potential_nozonal], color="crimson", label=label_nozonal) 
    
    # Axis limits
    ax.set_xlim(left=0, right=np.max(time))
    phi2 = phi2[int(len(phi2)/2):]
    if log==True:  ax.set_ylim(bottom=np.min(phi2), top=np.max(phi2))
    if log==False: ax.set_ylim(bottom=0, top=np.max(phi2)*1.1)
    
    # legend and labels
    ax.legend(labelspacing=0.5, handlelength=1, shadow=True, prop={'size':fontsize})
    ax.set_xlabel(standardLabels["normalized"]["t"])
    ax.set_ylabel(standardLabels["normalized"]["phi2"])
    return

#----------------------------------------
def subplot_potential_vs_time(ax, research, y_quantity="phi2", zonal=None, log=False, label_units=""):
    
    # Automate the axis limits  
    axis = Axis(ax, Plot(),  xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.1, logy=log, percentage=0.5)
    
    # Colors based on the number of simulations
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,research.numberOfSimulations)); count=0
    
    # Plot phi2(t)
    for experiment in research.experiments:
        for simulation in experiment.simulations:  
        
            # Get the data
            if y_quantity!="phi2": exit_program("Only y_quantity==phi2 has been implemented.", subplot_potential_vs_time_zonal_contributions, sys._getframe().f_lineno)
            if zonal==None: phi2 = simulation.potential.phi2_vs_t.phi2 
            if zonal==True: phi2 = simulation.potential.phi2_vs_t_zonal.phi2_zonal 
            if zonal==False: phi2 = simulation.potential.phi2_vs_t_nozonal.phi2_nozonal 
            time = simulation.potential.phi2_vs_t.t 
            
            # Get the label 
            tstart = simulation.time.tstart; tend = simulation.time.tend
            saturated_potential = np.mean(phi2[(time>=tstart) & (time<=tend)]) 
            label = "$\\langle \\tilde{\\varphi}_{\\text{sat}}^2 \\rangle_z = "+str("{:.2f}".format(saturated_potential))+"$"
            label = simulation.line_label+"; "+label if research.numberOfSimulations!=research.numberOfExperiments else experiment.line_label+"; "+label 
            
            # Plot the time evolution of the potential squared  
            ax.plot(time, phi2[:], color=colors[count])  
            ax.plot([tstart, tend], [saturated_potential, saturated_potential], color=colors[count], label=label) 
            axis.update_axisLimits(time, phi2) ; count +=1  
    
    # Appearance
    if zonal==None: ax.set_ylabel(standardLabels["normalized"]["phi2_fieldlineaverage"+label_units])
    if zonal==True: ax.set_ylabel(standardLabels["normalized"]["phi2_zonal_fieldlineaverage"+label_units])
    if zonal==False: ax.set_ylabel(standardLabels["normalized"]["phi2_nozonal_fieldlineaverage"+label_units])
    ax.legend(labelspacing=0.2, handlelength=1, shadow=True)
    
    # Axis limits
    axis.rescale_axis()
    return
            
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#=============================================================================== 
 
if __name__ == "__main__":
    
    # Check which script to launch      
    input_file = get_firstInputFile(pathlib.Path(os.getcwd()))  
    nonlinear, full_flux_surface = read_nonlinearFullFluxSurfaceFromInputFile(input_file)
    if not nonlinear: os.system("python3 $STELLAPY/plot/linear/potential_vs_time.py "+" ".join(sys.argv[1:])); sys.exit()
    
    # Create a bash-like interface
    bash = Bash(plot_potential_vs_time, __doc__)  
    
    # Choose the y-quantity
    bash.add_toggleheader("y_quantity")
    bash.add_toggle('y_quantity', 'phi2', '', '', 'Plot the time evolution of the potential squared.')  
    bash.add_toggle('zonal', True, 'z', 'zonal', 'Plot only the zonal modes.') 
    bash.add_toggle('zonal', False, 'n', 'nozonal', 'Plot only the non-zonal modes.') 
    bash.add_togglespace()  
    
    # Manipulate the data
    bash.add_optionheader("simulation") 
    bash.add_toggleheader("simulation")
    bash.add_toggle('folderIsExperiment', True, 'f', 'folderIsExp', 'Each folder is an experiment.') 
    bash.add_toggle('folderIsSimulation', True, '', 'all', 'Each folder is a simulation.')  
    bash.add_option('trange', 'range', 't', '', "Select the time frame for the saturated state of the simulations, e.g. '-t [500, 1000]'.") 
    bash.add_option('xrange', 'range', '', 'xrange', "Select the xrange for the figure, e.g. '[0,1000]'.") 
    bash.add_option('yrange', 'range', '', 'yrange', "Select the yrange for the figure, e.g. '[0,10]'.") 
    bash.add_optionspace()
    bash.add_togglespace()  
    
    # Change the appearance of the figure
    bash.add_optionheader("Figure") 
    bash.add_toggleheader("Figure") 
    bash.add_toggle('log', True, 'l', 'log', 'Use a logarithmic y-axis.') 
    bash.add_toggle('label_units', '_norm', '', 'normlabels', 'Put the y-axis in normalized symbols \\tilde{\\varphi}.') 
    bash.add_toggle('label_units', '', '', 'exactlabels', 'Put the y-axis in exact symbols \\varphi^2 ({\\rho_i T_i}/{ae})^2.') 
    bash.add_optionspace()
    bash.add_togglespace()
    
    # Other options 
    bash.add_toggleheader("other")
    bash.add_optionheader("other")    
    
    # Get the arguments and execute the script
    plot_potential_vs_time(**bash.get_arguments())   


