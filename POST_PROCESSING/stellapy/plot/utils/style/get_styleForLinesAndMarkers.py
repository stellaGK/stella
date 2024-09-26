
import copy
import numpy as np
import matplotlib.pyplot as plt 
from itertools import takewhile
from stellapy.utils.files import initiate_nesteddict
from stellapy.plot.utils.devices.get_colorsPerDevice import replaceWithCustomLineColors 
from stellapy.data.geometry.recognize_device import recognize_device

#----------------------------------
def get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation=None, mode=None): 
    """ Depending on the (experiments; simulations) we will represent the data
    with different styles of lines and markers. """

    # Standard style for lines/markers
    style = {"linewidth" : 3, "linestyle" : "-", "ms" : 6} 

    # Get the labels for the experiment and the simulation
    label_exp, label_sim = get_labels(plot, legend, research, experiment, simulation, mode)

    # Plot lines for each experiment  
    if research.plotting_option["LinesPerExperiment"] or simulation==None: 
        style["color"] = experiment.line_color
        style["label"] = label_exp 
        markers = {}
        
    if research.plotting_option["LinesPerSimulation"] and simulation!=None:  
        style["color"] = simulation.marker_color
        style["label"] = label_sim  
        markers = {}
                
    if research.plotting_option["MarkersPerSimulation"] and simulation!=None:
        style["linewidth"] = 0
        style["color"] = simulation.marker_color
        style["label"] = label_sim   
        markers["marker"] = simulation.marker_style
        markers["mec"] = simulation.marker_color
        markers["mfc"] = simulation.marker_color 
                
    # Return what to plot  
    return style

#----------------------------------
def get_labels(plot, legend, research, experiment, simulation=None, mode=None): 
    
    # Get the standard label of the experiment
    label_exp = experiment.line_label 
    
    # For each simulation, add a label to show what is varied between the simulations
    label_sim = simulation.marker_label if simulation!=None else "" 
    label_sim = label_sim.replace("_"," ") if "$" not in label_sim else label_sim

    # Linearly we add the (kx,ky) mode in the label 
    if mode!=None: 
        if research.numberOfPlottedModes!=research.numberOfPlottedSimulations or research.numberOfPlottedModes!=research.numberOfPlottedExperiments:        
            
            # Check whether we scanned at fixed kx or fixed ky
            scanAtFixedKx = (plot.kx_range[0]==plot.kx_range[1])
            scanAtFixedKy = (plot.ky_range[0]==plot.ky_range[1])
            scanDiffKxKy = (not scanAtFixedKx) and (not scanAtFixedKy)
            
            # The label for the mode (kx,ky) depends on whether kx or ky is fixed
            if scanAtFixedKx:   label_sim_temp = "$k_y\\rho_i = " + "{0:.2f}".format(mode.ky) + "$"  
            elif scanAtFixedKy: label_sim_temp = "$k_x\\rho_i = " + "{0:.2f}".format(mode.kx) + "$" 
            elif scanDiffKxKy:  label_sim_temp = "$(k_x\\rho_i, k_y\\rho_i) = ("+"{0:.2f}".format(mode.kx)+", "+"{0:.2f}".format(mode.ky)+ ")$" 
            
            # If we have multiple experiments/simulations then add that label as well
            if   research.numberOfPlottedSimulations>1: label_sim += ": " + label_sim_temp      
            elif research.numberOfPlottedExperiments>1: label_sim += ": " + label_exp 
            else:                                       label_sim = label_sim_temp
    
    # Custom lines
    if not (("{ne}" in label_exp) and ("{ni}" in label_exp)): label_exp = label_exp.replace("{ne}","{n}")
    if not (("{ne}" in label_sim) and ("{ni}" in label_sim)): label_sim = label_sim.replace("{ne}","{n}")  
    
    # Make sure the labels work with latex
    label_exp = label_exp.replace("wout_", "wout ")
    label_sim = label_sim.replace("wout_", "wout ")

    # Don't use the same label twice
    label_exp = label_exp if (label_exp not in legend.labels_lines) else ""
    label_sim = label_sim if (label_sim not in legend.labels_points) else ""

    # Remember which labels we already plotted for the simulations and experiments
    legend.labels_points.append(label_sim) 
    legend.labels_lines.append(label_exp)  
        
    # Return the experiment and simulation label
    return label_exp, label_sim

#===============================================================================
#                           FIXED FOR THE RESEARCH
#===============================================================================

def get_plottingOption(experiments):    
    """ Depending on the (experiments; simulations) we will represent the data
    with different styles of lines and markers. """
    
    # Get the number of varied values
    from stellapy.simulations.utils.get_numberOfVariedValues import get_numberOfVariedValues
    number_ofVariedValues = get_numberOfVariedValues(experiments)  
    
    # Plot the lines with the same color for each simulation if there is only one varied value  
    if (number_ofVariedValues==1) and (len(experiments)==1):  
        LinesPerExperiment   = False
        LinesPerSimulation   = True 
        MarkersPerSimulation = False
        
    # Plot the lines with the same color for each experiment if there are multiple experiments
    elif (number_ofVariedValues==1) and (len(experiments)>1):  
        LinesPerExperiment   = True
        LinesPerSimulation   = False 
        MarkersPerSimulation = False
    
    # If there is only one experiment and multiple simulations, plot lines per simulation
    elif (number_ofVariedValues!=1) and (len(experiments)==1): 
        LinesPerExperiment   = False
        LinesPerSimulation   = True 
        MarkersPerSimulation = False
    
    # Plot the lines with the same color for each experiment if there are multiple experiments
    # Then add markers to differentiate between different simulations
    elif (number_ofVariedValues!=1) and (len(experiments)>1): 
        LinesPerExperiment   = True
        LinesPerSimulation   = False 
        MarkersPerSimulation = True
    
    return {"LinesPerExperiment" : LinesPerExperiment, "LinesPerSimulation" : LinesPerSimulation, "MarkersPerSimulation" : MarkersPerSimulation}


#----------------------------------
def load_labelsLinesMarkers(experiments, varied_values, return_dictionaries=True):
    """ Load labels, lines, color and markers for the (simulation, groups) combinations.     

    Returns
    -------
    line_style, line_color, marker_color, marker_style : dictionaries with keys [simlutation] and [groups]
    """

    # Initiate the dictionaries
    line_label   = initiate_nesteddict()
    marker_label = initiate_nesteddict()
    line_style   = initiate_nesteddict()
    line_color   = initiate_nesteddict()
    marker_color = initiate_nesteddict()
    marker_style = initiate_nesteddict() 
    
    # Automatically decide on an appropriate color map    
    if len(experiments)<=3:     l_colors = ["navy", "crimson", "black"]
    if len(varied_values)<=3:   m_colors = ["navy", "crimson", "black"]
    if len(experiments)>3:      l_colors = "jet"
    if len(varied_values)>3:    m_colors = "jet"
    #if l_colors==m_colors:      l_colors = ["lime", "magenta", "green"] if isinstance(l_colors, list) else "rainbow"
        
    # Load the color map 
    if isinstance(m_colors, str): m_colors = copy.copy(plt.get_cmap(m_colors))( np.linspace(0,1,len(varied_values)) )
    if isinstance(l_colors, str): l_colors = copy.copy(plt.get_cmap(l_colors))( np.linspace(0,1,len(experiments)) ) 
        
    # Initiate line and marker styles to choose from
    l_styles = ["-", "--", ":","-.","densely dotted","densely dashed","densely"]*100
    m_styles = ["o","X","D","P","s","x","d","p"]*100
    
    # Find the common prefix of the experiments
    common_prefix = "".join(c[0] for c in takewhile(lambda x:  all(x[0] == y for y in x), zip(*experiments))) 
        
    # Add lines/marker styles and colors for each experiment
    for experiment in experiments:
        line_color[experiment]   = l_colors[experiments.index(experiment)]  
        line_style[experiment]   = l_styles[experiments.index(experiment)]
        marker_color[experiment] = l_colors[experiments.index(experiment)]    
        marker_style[experiment] = m_styles[experiments.index(experiment)]  

    # Add lines/marker styles and colors for each simulation
    varied_values_with_common_prefixes = copy.deepcopy(varied_values)
    common_prefix_simulation = "".join(c[0] for c in takewhile(lambda x:  all(x[0] == y for y in x), zip(*varied_values))) 
    if ";" in common_prefix_simulation: 
        common_prefix_simulation = "; ".join(common_prefix_simulation.split('; ')[:-1])+"; "
        for i in range(len(varied_values)): varied_values_with_common_prefixes[i] = varied_values[i].replace(common_prefix_simulation, "");  
    for i, variation in enumerate(varied_values):  
        line_label[variation]    = str(varied_values_with_common_prefixes[i])
        marker_label[variation]  = str(varied_values_with_common_prefixes[i]) 
        line_color[variation]    = m_colors[i]  
        line_style[variation]    = "-"
        marker_color[variation]  = m_colors[i]   
        marker_style[variation]  = m_styles[i] 
        marker_label[variation] = marker_label[variation].replace("_"," ") if "$" not in marker_label[variation] else marker_label[variation]
        line_label[variation] = line_label[variation].replace("_"," ") if "$" not in line_label[variation] else line_label[variation]


    # Edit the line/marker labels for the experiments/simulations
    for experiment in experiments:
        if common_prefix==experiments[0] or common_prefix=="":
            line_label[experiment] = str(experiment).replace("nl_","").replace("th_","").replace("_", " ")
        if "=" in common_prefix:
            line_label[experiment] = str(experiment)
        else:
            line_label[experiment] = str(experiment).replace("_", " ")
        marker_label[experiment] = str(experiment).replace("__", ": ").replace("_", " ") 

    # Replace some custom line colors  
    label = line_label[list(line_label.keys())[0]]
    if ";" not in label and recognize_device(label):
        line_color = replaceWithCustomLineColors(experiments, line_label, line_color) 

    # Return a single choice
    if return_dictionaries==False:
        return line_label[varied_values[0]], marker_label[varied_values[0]], line_style[varied_values[0]], line_color[varied_values[0]], marker_style[varied_values[0]], marker_color[varied_values[0]]
    
    # Return the lines/markers labels/styles/colors 
    return line_label, marker_label, line_style, line_color, marker_style, marker_color
 
