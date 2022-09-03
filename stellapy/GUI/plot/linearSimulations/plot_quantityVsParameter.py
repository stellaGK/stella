

# Load modules 
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

# Personal modules
from stellapy.utils.decorators.verbose import noverbose   
from stellapy.GUI.plot.utils.get_dataForPlots import get_dataForPlots  

# Plotting modules
from stellapy.GUI.plot.utils.Plot import Plot
from stellapy.GUI.plot.utils import Axis, Legend
from stellapy.GUI.plot.utils import create_figure
from stellapy.GUI.utils import DisplayProgressGUI 
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import get_styleForLinesAndMarkers  
from stellapy.plot.utils.devices.get_colorsPerDevice import get_markerPerDevice, get_colorsPerDevice

#===============================================================================
#                PLOT GAMMA(PARAMETER) OF THE MOST UNSTABLE MODE               #
#===============================================================================

@noverbose
def plot_quantityVsParameter(\
            # Simulations
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
            # Parameter
                knob="vmec_parameters",\
                key="rho",\
            # Data on the (x,y) axis  
                x_quantity='rho',\
                y_quantity='omega_avg',\
                x_range=None,\
                y_range=None,\
                units="normalized",\
            # Modes
                modes_id="all",\
                kx_range=[0,0],\
                ky_range=[0,100], 
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
                fontsize_title=None,\
            # Legend
                ls=None,\
                lw=None,\
                ms=None,\
                color=None,\
                label=None,\
                marker=None,\
                removelegend=False,\
                fontsize=20,\
            # Figure
                ax=None, \
                Progress=None,\
                show_figure = False,\
            # Toggles
                interpolate=False,\
                show_error = False,\
                add_dottedLine = False,\
                splitInKineticAdiabatic = False): 
    
    # We can look at a speficic_k
    specific_k = ky_range[0] if (kx_range[0]==kx_range[1] and ky_range[0]==ky_range[1]) else None
 
    # Save the plotting details to a <plot> object
    plot = Plot()  
    plot.update_legend(fontsize, splitInKineticAdiabatic=splitInKineticAdiabatic, removelegend=removelegend)
    plot.update_labels(title, x_label, y_label, fontsize_title=fontsize_title) 
    plot.update_toggles(show_error=show_error, add_dottedLine=add_dottedLine, interpolate=interpolate)
    plot.update_figure(Progress, ax, show_figure, logy=(y_quantity=="phi2"))
    plot.update_simulations(experiment_id, simulation_id, species=[0]) 
    plot.update_xyData(x_quantity, y_quantity, x_range, y_range, units) 
    plot.update_modes(modes_id, kx_range, ky_range, specific_k)
    plot.update_parameter(knob, key)
    plot.process_plottingVariables(research)
        
    # Update the progress bar of the GUI 
    displayProgressGUI1 = DisplayProgressGUI(research, "Read "+plot.yname+" versus "+plot.xname) 
    displayProgressGUI2 = DisplayProgressGUI(research, "Plot "+plot.yname+" versus "+plot.xname) 
         
    # Create the figure and a <legend> and <axis> object
    ax = create_figure(ax, plot, research) 
    legend = Legend(ax, plot); axis = Axis(ax, plot) 
    axis.set_limits(xbot_pos=0, ybot_pos=0) 
 
    # Set some fixed axis limits depending on the plot
    if x_quantity=="rk":  axis.set_limits(xbot=1.5, xtop=3.5, ybot_pos=0)
    if x_quantity=="rho": axis.set_limits(xbot=0, xtop=1, ybot_pos=0)
    
    # For the labels we needed the real xdim but for the data we need "ky"
    plot.x_quantity="ky"; plot.process_plottingVariables(research)
     
    # Iterate over the experiments
    for experiment in research.plotted_experiments: 
         
        # First do the calculations for this experiment
        parameters, simulations = get_parameters(experiment, key, knob)
        quantity, error = get_linearData(experiment, simulations, plot, displayProgressGUI1)
      
        # Now plot the data 
        plot_data(ax, parameters, quantity, error, research, experiment, plot, legend, axis, displayProgressGUI2, color, marker, ls, label)
        
    # Add the legend and rescale the axis 
    if not removelegend: legend.add_legend()
    axis.rescale_axis()
 
    # Show the figure 
    if show_figure: plt.show()
     
    # Update the progress bar of the GUI  
    displayProgressGUI2.finilize_progressBar()
    if True: return
    
#===============================================================================
#                              CALCULATE THE DATA                              #
#=============================================================================== 

def get_parameters(experiment, key, knob):
    
    # Intitiate the parameters
    parameters = [] 
     
    # Get the requested parameter
    for simulation in experiment.plotted_simulations: 
        if key!="explicit_option" and key!="nfield_periods":
            value = float(simulation.input.inputParameters[knob][key])
        if key=="nfield_periods": 
            if simulation.input.vmec:     value = simulation.input.nfield_periods
            if simulation.input.miller:   value = simulation.input.nperiod 
        if key=="explicit_option": 
            if (simulation.input.inputParameters[knob][key]=="rk2"): value = 2 
            if (simulation.input.inputParameters[knob][key]=="rk3"): value = 3
        parameters.append(value) 
     
    # Sort the simulations by this parameter
    sorted_indexes = list(np.array(parameters).argsort(kind="heapsort"))
    simulations = [experiment.plotted_simulations[i] for i in sorted_indexes]
    parameters  = [parameters[i] for i in sorted_indexes]
    return parameters, simulations

#------------------------------------
def get_linearData(experiment, simulations, plot, displayProgressGUI1):
    
    # Initiate the data
    quantity = [np.NaN]*len(experiment.plotted_simulations)
    error = np.ones((2, len(experiment.plotted_simulations)))*np.NaN

    # Iterate over the simulations
    for i, simulation in enumerate(simulations): 
         
        # Update the progress bar of the GUI  
        if displayProgressGUI1: displayProgressGUI1.move_progressBar()
        
        # Get ky**2
        if plot.quantity=="ky**2_vs_ky":
            plot.ydim = "gamma_avg"
            plot.quantity = "gamma_avg_vs_ky"
            simulation.lineardata.get_linearDataPerSimulation(simulation.plotted_modes, plot)  
            _, vec_gamma    = get_dataForPlots(plot, simulation, "gamma_avg_vs_ky")  
            index = np.nanargmax(vec_gamma) 
            quantity[i] = [mode.ky for mode in simulation.plotted_modes][index]   
            plot.ydim = "ky**2"
            plot.quantity = "ky**2_vs_ky"
        
        else:
            
            # Calculate the linear data for these plotted modes and save them as e.g. simulation.lineardata.gamma_avg_vs_ky
            simulation.lineardata.get_linearDataPerSimulation(simulation.plotted_modes, plot)  
                 
            # Get the data
            _, vec_quantity = get_dataForPlots(plot, simulation, plot.quantity) 
            _, vec_error    = get_dataForPlots(plot, simulation, plot.quantity+"_error")  
            
            # Get the gamma vector
            vec_gamma = [mode.lineardata.gamma_avg if mode.lineardata.unstable else 0 for mode in simulation.plotted_modes]
            
            # Find the most unstable mode 
            index = np.nanargmax(vec_gamma) 
            quantity[i] = vec_quantity[index]
            error[:,i] = vec_error[:,index] 
        
    return quantity, error

#===============================================================================
#                                PLOT THE DATA                                 #
#=============================================================================== 

def plot_data(ax, parameters, quantity, error, research, experiment, plot, legend, axis, displayProgressGUI2, color, marker, ls, label):
    
    # Update the progress bar of the GUI  
    displayProgressGUI2.move_progressBar() 
    
    # Get the lines and marker styles for the experiments
    style = get_styleForLinesAndMarkers(plot, legend, research, experiment) 
    style["marker"] = experiment.marker_style if not (get_markerPerDevice(style["label"])) else get_markerPerDevice(style["label"])
    
    # Represent kinetic simulations with full lines, and adiabatic electrons 
    # with dotted lines and adiabatic ions with striped lines  
    if plot.splitInKineticAdiabatic and "nspec" in style["label"]: 
        charge = experiment.simulations[0].input.inputParameters["species_parameters_1"]["z"]
        adiabatic = "\\texttt{nspec}$\\,=\\,$1" in style["label"]
        lss = "--" if (adiabatic and charge>0) else (":" if (adiabatic and charge<0) else "-")
        label_exp_new = "" if "\\texttt{nspec}$\\,=\\,$1" in style["label"] else style["label"]
        label_exp_new = label_exp_new.replace("\\texttt{nspec}$\\,=\\,$1", "").replace("\\texttt{nspec}$\\,=\\,$2", "").replace(";","")
        if label_exp_new!="" and style["label"] in legend.labels_points: 
            legend.labels_points.remove(style["label"])
            if label_exp_new not in legend.labels_points: legend.labels_points.append(label_exp_new)
        if label_exp_new!="" and style["label"] in legend.labels_lines:  
            legend.labels_lines.remove(style["label"])  
            if label_exp_new not in legend.labels_lines: legend.labels_lines.append(label_exp_new)
        style["label"] = label_exp_new 
        style["linestyle"] = lss
    style["marker"] = get_markerPerDevice(experiment.simulations[0].modes[0].input_file)
    style["color"] = get_colorsPerDevice(style["label"])
    
    # Overwrite styles
    style["linestyle"] = ls if (ls!=None) else style["linestyle"]
    style["color"] = color if (color!=None) else style["color"]
    style["label"] = label if (label!=None) else style["label"]
    style["marker"] = marker if (marker!=None) else style["marker"]

    # Keep track of the axis limits 
    axis.update_axisLimits(parameters, quantity) 
    
    # Interpolate
    if plot.interpolate:
        ax.plot(parameters[0], quantity[0], **style)
        lw = style['linewidth']; style['linewidth'] = 0; max = np.max(parameters)
        parameters = np.array(parameters); quantity = np.array(quantity)
        ax.plot(parameters[quantity!=0], quantity[quantity!=0], **style); style['linewidth'] = lw; style['ms'] = 0
        parameters = np.array(parameters)[np.array(quantity)>0]
        quantity = np.array(quantity)[np.array(quantity)>0] 
        f = interpolate.interp1d(parameters, quantity, fill_value="extrapolate")
        parameters = np.linspace(np.min(parameters), max, 100)
        quantity = f(parameters)
        
    # Plot the (parameter, quantity) data
    if not plot.show_error:  
        ax.plot(parameters, quantity, **style)   
     
    # Plot the (parameter, quantity, error) data
    if plot.show_error:  
        del style["marker"]
        ax.errorbar(parameters, quantity, yerr=error, fmt="o", capsize=2, **style) 
     
    # Add a dotted line to show the saturated value
    if plot.add_dottedLine:
        if plot.key not in ["dmu", "dvpa", "delt", "dkx", "dky"]:  
            vec_x = np.linspace(0, parameters[-1], 10) 
            vec_y = np.ones(len(vec_x))*quantity[-1]
        if plot.key in ["dmu", "dvpa", "delt", "dkx", "dky"]:  
            vec_x = np.linspace(parameters[0]*0.0001, parameters[-1]*2, 10) 
            vec_y = np.ones(len(vec_x))*quantity[0]
        ax.plot(vec_x, vec_y, lw=2, linestyle="dotted", color=style['color']) 
    return


################################################################################
#                  USE THIS PLOTTING FUNCTION AS A MAIN SCRIPT                 #
################################################################################ 
    
if __name__ == "__main__":   
    
    import pathlib, time
    from stellapy.simulations.Research import create_research  
    start = time.time()
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI_PREVIOUS/fprim2tprim6/")  
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI_PREVIOUS")    
    research = create_research(folders=folder); research.print_research()  
    plot_quantityVsParameter(research, knob="species_parameters_1", key="tprim", x_quantity="tiprim")
    print("    ---> It took "+str(time.time() - start)+" seconds to plot a quantity versus parameter.")
    plt.show()
    