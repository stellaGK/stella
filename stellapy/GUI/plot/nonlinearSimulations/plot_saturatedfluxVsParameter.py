
# Load the modules
import numpy as np 
import matplotlib.pyplot as plt

# Personal modules
from stellapy.utils.decorators.verbose import noverbose  

# Plotting modules 
from stellapy.GUI.utils import DisplayProgressGUI  
from stellapy.GUI.plot.utils.Plot import Plot
from stellapy.GUI.plot.utils import Axis, Legend
from stellapy.GUI.plot.utils import create_figure
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import get_styleForLinesAndMarkers 
from stellapy.plot.utils.devices.recognize_device import recognize_device
from stellapy.plot.utils.devices.get_colorsPerDevice import get_colorsPerDevice,\
    get_markerPerDevice 
 
#===============================================================================
#                      PLOT SATURATED HEAT FLUX(PARAMETER)                     #
#===============================================================================

@noverbose 
def plot_saturatedfluxVsParameter(\
            # Specify which simulations to plot
                research=None,\
                experiment_id='All experiments',\
                simulation_id='All simulations',\
                species=[0],\
            # Parameter
                knob='vmec_parameters',\
                key='rho',\
            # Data
                x_quantity='rho',\
                y_quantity='qflux',\
                x_range=None,\
                y_range=None,\
                units='normalized',\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # Legend 
                shadow=True,\
                fontsize=20,\
                removelegend = False,\
                splitInKineticAdiabatic = False,\
                splitInKineticAdiabaticTEM = False,\
            # Toggles
                lw=None,\
                ls=None,\
                ms=None,\
                color=None,\
                marker=None,\
                verbose=False,\
                normalize=False,\
                show_errorStd=False,\
                show_errorMinMax=False,\
            # For the GUI the figure object already exists 
                ax=None,\
                Progress=None,\
                show_figure=True): 
    """ Dotted lines at 10, 20 and 50%. """
    
    # Save the plotting details to a <plot> object
    plot = Plot()  
    plot.update_legend(fontsize, shadow=shadow, removelegend=removelegend, splitInKineticAdiabatic=splitInKineticAdiabatic, splitInKineticAdiabaticTEM=splitInKineticAdiabaticTEM)
    plot.update_toggles(show_errorStd=show_errorStd, show_errorMinMax=show_errorMinMax, normalize=normalize)
    plot.update_xyData(x_quantity, y_quantity, x_range, y_range, units)  
    plot.update_simulations(experiment_id, simulation_id, species) 
    plot.update_figure(Progress, ax, show_figure)
    plot.update_labels(title, x_label, y_label)  
    plot.update_parameter(knob, key) 
    plot.process_plottingVariables(research)  

    # Update the progress bar of the GUI 
    displayProgressGUI1 = DisplayProgressGUI(research, "Read "+plot.yname+" versus "+plot.xname) 
    displayProgressGUI2 = DisplayProgressGUI(research, "Plot "+plot.yname+" versus "+plot.xname) 
 
    # Create the figure and a <legend> and <axis> object
    ax = create_figure(ax, plot, research) 
    legend = Legend(ax, plot); axis = Axis(ax, plot) 
    axis.set_limits(ytop_neg=0,  ybot_pos=0, )   
    
    # Save the data to a dictionary
    savedata = {}

    # Iterate over the experiments 
    for experiment in research.plotted_experiments: 
        for specie in species:
         
            # First do the calculations for this experiment
            parameters, simulations = get_parameters(experiment, key, knob)
            quantity, error_std, error_minmax = get_saturatedFlux(experiment, simulations, plot, specie, displayProgressGUI1)
          
            # Now plot the data 
            plot_data(ax, parameters, quantity, error_std, error_minmax, research, experiment, plot, legend, axis, displayProgressGUI2, savedata, color, ls, lw, ms, marker)
        
            # Print some information
            if verbose: print_relative_error_saturated_heat_flux(parameters, quantity, plot)   
            
    # Add the legend and rescale the axis 
    if not removelegend: legend.add_legend()
    axis.rescale_axis()
 
    # Show the figure 
    if show_figure: plt.show()
     
    # Update the progress bar of the GUI  
    displayProgressGUI2.finilize_progressBar()
    return savedata

#===============================================================================
#                              CALCULATE THE DATA                              #
#=============================================================================== 

def get_parameters(experiment, key, knob):
    
    # Initiate the parameters
    parameters = [] 
     
    # Get the requested parameter
    for simulation in experiment.plotted_simulations:
        simulation.input.inputParameters["parameters"]["teti"] = 1/simulation.input.inputParameters["parameters"]["tite"]
        if key!="explicit_option" and key!="nfield_periods":
            value = float(simulation.input.inputParameters[knob][key])
        if key=="nfield_periods": 
            if simulation.input.vmec:     value = simulation.input.nfield_periods
            if simulation.input.miller:   value = simulation.input.nperiod  
        if key=="poloidal_turns": 
            if simulation.input.vmec:     value = float(simulation.input.inputParameters[knob][key])
            if simulation.input.miller:   value = 2*(simulation.input.nperiod-1)+1 
        if key=="explicit_option": 
            if (simulation.input.inputParameters[knob][key]=="rk2"): value = 2 
            if (simulation.input.inputParameters[knob][key]=="rk3"): value = 3
        parameters.append(value)
     
    # Sort the simulations by this parameter
    sorted_indexes = list(np.array(parameters).argsort(kind="heapsort"))
    simulations = [experiment.plotted_simulations[i] for i in sorted_indexes]
    parameters  = [parameters[i] for i in sorted_indexes]
    return parameters, simulations

#------------------
def get_saturatedFlux(experiment, simulations, plot, specie, displayProgressGUI1):
    
    # Initiate the data
    quantity = [np.NaN]*len(experiment.plotted_simulations)
    error_std = [np.NaN]*len(experiment.plotted_simulations)
    error_minmax = np.ones((2, len(experiment.plotted_simulations)))*np.NaN

    # Iterate over the simulations 
    for i, simulation in enumerate(simulations): 
         
        # Update the progress bar of the GUI  
        displayProgressGUI1.move_progressBar() 
             
        # Get the fluxes or potential data
        if "flux" in plot.y_quantity: 
            minimum = simulation.time.satFluxMinimum[plot.y_quantity][specie]
            maximum = simulation.time.satFluxMaximum[plot.y_quantity][specie]
            satflux = simulation.time.saturatedFluxes[plot.y_quantity][specie]
            stderror = simulation.time.satFluxStdErrors[plot.y_quantity][specie]
        if "phi" in plot.y_quantity:
            minimum = simulation.time.satPotMinimum[plot.y_quantity]
            maximum = simulation.time.satPotMaximum[plot.y_quantity]
            satflux = simulation.time.saturatedPotential[plot.y_quantity]
            stderror = simulation.time.satPotStdErrors[plot.y_quantity]
        if "g" in plot.y_quantity:
            minimum = simulation.time.satDistMinimum[plot.y_quantity]
            maximum = simulation.time.satDistMaximum[plot.y_quantity]
            satflux = simulation.time.saturatedDistribution[plot.y_quantity]
            stderror = simulation.time.satDistStdErrors[plot.y_quantity]
        
        # Correct the units
        if plot.units=="SI" and "flux" in plot.y_quantity: 
            minimum = minimum*simulation.referenceunits.ref[plot.y_quantity] 
            maximum = maximum*simulation.referenceunits.ref[plot.y_quantity] 
            satflux = satflux*simulation.referenceunits.ref[plot.y_quantity] 
            stderror = stderror*simulation.referenceunits.ref[plot.y_quantity]   
        if plot.units=="SI" and "phi" in plot.y_quantity: 
            print("ERROR: impement the SI units for potential.")
            import sys; sys.exit(1) 
        
        # Save the data for simulation <i> 
        quantity[i] = satflux
        error_std[i] = stderror
        error_minmax[:,i] = [minimum, maximum]
        
    # Normalize to see only the shape
    if plot.normalize: 
        error_minmax = np.array(error_minmax)/quantity[-1] 
        error_std = np.array(error_std)/quantity[-1] 
        quantity = np.array(quantity)/quantity[-1] 
        
    return quantity, error_std, error_minmax

#===============================================================================
#                                PLOT THE DATA                                 #
#=============================================================================== 

def plot_data(ax, parameters, quantity, error_std, error_minmax, research, experiment, plot, legend, axis, displayProgressGUI2, savedata, color, ls, lw, ms, marker):
    
    # Update the progress bar of the GUI  
    displayProgressGUI2.move_progressBar() 
    
    # Get the lines and marker styles for the experiments
    style = get_styleForLinesAndMarkers(plot, legend, research, experiment) 
    if "marker" in style: del style["marker"]
    if "{Te}" in style["label"] and not marker: marker = "s" if "3.0" in style["label"] else "o" 
    marker = experiment.marker_style if not marker else marker; capsize = 2
    style["color"] = color if (color!=None and "Low resolution" in style["label"]) else style["color"] 
    style["color"] = color if (color!=None) else style["color"] 
    style["linestyle"] = ls if (ls!=None) else style["linestyle"] 
    style["linewidth"] = lw if (lw!=None) else style["linewidth"] 
    style["ms"] = ms if (ms!=None) else style["ms"] 
    
    style["marker"] = get_markerPerDevice(style["label"])
    style["color"] = get_colorsPerDevice(style["label"]) 
    
    # Represent kinetic simulations with full lines, and adiabatic electrons 
    # with dotted lines and adiabatic ions with striped lines
    if (plot.splitInKineticAdiabatic or plot.splitInKineticAdiabaticTEM) and "nspec" in style["label"]: 
        style["color"], style["label"], style["linestyle"] = get_adiabaticKineticData(ax, experiment, style["label"], legend, plot) 
        
    # If we do a resolution scan e.g. nvgrid: 24 -> 48 then don't connect the points
    if "rightarrow" in style["label"]: style["linewidth"] = 0;  capsize = 0
    
    # Put the exponent in the label 
    if plot.units=="SI":
        if plot.y_quantity=="qflux": 
            ax.set_ylabel("$Q^\\text{turb}_{i}$ [$10^5$ W/m$^2$]")
            quantity = np.array(quantity)*1e-5
            error_std = np.array(error_std)*1e-5
            error_minmax = np.array(error_minmax)*1e-5
        if plot.y_quantity=="pflux": 
            ax.set_ylabel("$\Gamma^{\\text{turb}}_s\,$[$10^{19}\,$m$^{-2}$s$^{-1}$]")
            quantity = np.array(quantity)*1e-19
            error_std = np.array(error_std)*1e-19
            error_minmax = np.array(error_minmax)*1e-19 
    
    # Get the error
    error = error_std if plot.show_errorStd else (error_minmax if plot.show_errorMinMax else None) 
     
    # Plot the (parameter, quantity, error) data  
    ax.errorbar(parameters, quantity, yerr=error, fmt=marker, capsize=capsize, **style)  #mew=0.1, 
     
    # Add a dotted line to show the saturated value
    if plot.add_dottedLine:
        if plot.key not in ["dmu", "dvpa", "delt", "dkx", "dky"]:  
            vec_x = np.linspace(0, parameters[-1], 10) 
            vec_y = np.ones(len(vec_x))*quantity[-1]
        if plot.key in ["dmu", "dvpa", "delt", "dkx", "dky"]:  
            vec_x = np.linspace(parameters[0]*0.0001, parameters[-1]*2, 10) 
            vec_y = np.ones(len(vec_x))*quantity[0]
        ax.plot(vec_x, vec_y, lw=2, linestyle="dotted", color=style['color']) 
        
    # If we have a resolution scan of a single (fprim, tprim) then show the errors   
    if len(parameters)==1 and style["label"]=="Low resolution ": 
        ax.hlines(quantity[0]*1.1, -100, 100, ls=":", color="grey")
        ax.hlines(quantity[0]*0.9, -100, 100, ls=":", color="grey")
        ax.hlines(quantity[0]*1.2, -100, 100, ls="--", color="grey")
        ax.hlines(quantity[0]*0.8, -100, 100, ls="--", color="grey")
        ax.hlines(quantity[0]*1.5, -100, 100, ls="-", color="grey")
        ax.hlines(quantity[0]*0.5, -100, 100, ls="-", color="grey")
        axis.set_limits(xbot=parameters[0]-0.5, xtop=parameters[0]+0.1, ybot_pos=0)
             
    # Keep track of the axis limits 
    axis.update_axisLimits(parameters, quantity)
    
    # Save the data
    savedata[style["label"]] = [parameters, quantity]
    return 

#===============================================================================
#                                   METHODS                                    #
#=============================================================================== 
  
def print_relative_error_saturated_heat_flux(parameters, vec_satflux, plot):  
    average = np.mean(vec_satflux)
    print()
    print("==============================================")
    print("         ANALYSIS: RELATIVE ERROR             ")
    print("==============================================")
    print()
    print("   The average saturated heat flux is", str(round(average,2)),"\n")
    for i in range(len(parameters)):
        sentence = plot.key+" = "+str(parameters[i])
        relative_error = str(round(np.abs(average-vec_satflux[i])/average*100,2))+"%"
        satflux = "Qsat = "+str(round(vec_satflux[i],2))
        arrow = "   --->   "
        print("      "+"{:<15}".format(sentence)+arrow+"{:>6}".format(relative_error)+arrow+"{:<10}".format(satflux))
    print()
    print("==============================================")
    print()
    return 
 
#----------------------------
def get_adiabaticKineticData(ax, experiment, label_exp, legend, plot):
     
    # Determine whether the simulation has run with adiabtic electrons/ions
    charge = experiment.simulations[0].input.inputParameters["species_parameters_1"]["z"] 
    adiabatic = "\\texttt{nspec}$\\,=\\,$1" in label_exp
    if (adiabatic and charge>0): legend.adiabaticElectrons = True
    if (adiabatic and charge<0): legend.adiabaticIons = True
     
    # Check whether we are tweaking the ion temperature gradient
    TEM = ("a/L_{T_i}" in label_exp) and ("0.0" in label_exp)
     
    # Represent kinetic both with '-', adiabatic electrons with '--' and adiabatic ions with ':'
    if plot.splitInKineticAdiabatic: 
        ls = "--" if (adiabatic and charge>0) else (":" if (adiabatic and charge<0) else "-")
     
    # Represent kinetic both with '-', adiabatic electrons with '--' and TEM ':'
    if plot.splitInKineticAdiabaticTEM: 
        ls = "--" if ((not TEM) and adiabatic) else (":" if (TEM) else "-")
     
    # If we have only one device, make sure the colors are correct
    device = recognize_device(label_exp)  
    color  = get_colorsPerDevice(device) if (device!=None) else None
    device = recognize_device(str(experiment.simulations[0].input_file)) if (color==None) else device
    color  = get_colorsPerDevice(device) if (device!=None and color==None) else color
    experiment.line_color = color if (color!=None) else experiment.line_color 
     
    # Update the label since we'll have a seperate legend for the species: 
    #     "DEVICE; \\texttt{nspec}$\\,=\\,$1"                 -->  "DEVICE"
    if plot.splitInKineticAdiabatic: 
         
        # Create a new label
        label_exp_new = "" if "\\texttt{nspec}$\\,=\\,$1" in label_exp else label_exp
        label_exp_new = label_exp_new.replace("\\texttt{nspec}$\\,=\\,$1", "").replace("\\texttt{nspec}$\\,=\\,$2", "").replace(";","")
     
        # Update the legend lists
        if label_exp_new!="" and label_exp in legend.labels_points: 
            legend.labels_points.remove(label_exp)
            if label_exp_new not in legend.labels_points: legend.labels_points.append(label_exp_new)
        if label_exp_new!="" and label_exp in legend.labels_lines:  
            legend.labels_lines.remove(label_exp)  
            if label_exp_new not in legend.labels_lines: legend.labels_lines.append(label_exp_new)
        label_exp = label_exp_new
         
    # Update the labels: 
    #     "$a/L_{T_i}$$\,=\,$3.0; \\texttt{nspec}$\\,=\\,$2"   -->  "a/LTi=3 with kinetic electrons"    "-"
    #     "$a/L_{T_i}$$\,=\,$3.0; \\texttt{nspec}$\\,=\\,$1"   -->  "a/LTi=3 with adiabatic electrons"  "--"
    #     "$a/L_{T_i}$$\,=\,$0.0; \\texttt{nspec}$\\,=\\,$2"   -->  "a/LTi=0 with kinetic electrons"    ":"
    if plot.splitInKineticAdiabaticTEM:
         
        # Create a new label
        label_exp_new = label_exp
        label_exp_new = label_exp_new.replace("; \\texttt{nspec}$\\,=\\,$2", " with kin.\ elec.")
        label_exp_new = label_exp_new.replace("; \\texttt{nspec}$\\,=\\,$1", " with adi.\ elec.")
        color = "red" if "0.0" in label_exp_new else ("green" if "adi" in label_exp_new else color) 
        ax.plot(-100, -100, ls=ls, label=label_exp_new, color=color)
        legend.labels_lines.append(label_exp_new)
        if label_exp in legend.labels_points: legend.labels_points.remove(label_exp)
        if label_exp in legend.labels_lines:  legend.labels_lines.remove(label_exp)  
        legend.font_size = 18
        label_exp = ""
     
    # Return
    return color, label_exp, ls 

################################################################################
#                  USE THIS PLOTTING FUNCTION AS A MAIN SCRIPT                 #
################################################################################ 
     
if __name__ == "__main__":   
     
    import pathlib, time
    from stellapy.simulations.Research import create_research  
    start = time.time()
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/fprim2tprim3/")    
    research = create_research(folders=folder) 
    research.print_research()  
    plot_saturatedfluxVsParameter(research, show_figure=False)
    print("    ---> It took "+str(time.time() - start)+" seconds to plot a quantity versus parameter.")
    plt.show()    
     
    

