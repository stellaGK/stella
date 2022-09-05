 
#================================================================
# Plot flux(t) for non-linear runs
#================================================================
# TODO: Fluxes aren't normalized.
 
# Load the modules
import numpy as np
import matplotlib.pyplot as plt

# Personal modules
from stellapy.utils.decorators import verbose
from stellapy.simulations.utils.get_simulations import get_simulations
from stellapy.simulations.utils.get_experiments import get_experiments  
from stellapy.utils.files.sort_listByNumbers import sort_listByNumbers
from stellapy.data.geometry import calculate_geometricQuantitiesSTELLA
# from stellapy.calculations.calculate_geometricQuantities import calculate_geometricQuantities
 
@verbose 
def plot_geometryVsZ(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
            # Specify data range 
                x_quantity='zeta',\
                y_quantity='bmag',\
                x_range=None,\
                y_range=None,\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # For the GUI the figure object already exists 
                show_figure = True,\
                ax=None,\
                Progress=None,\
            # Toggles
                tick_style='sci',\
                normalize=False,\
                showOnlyMarkers=False,\
                log=False,\
                units="normalized"):
    
    ''' Note that 
            -gbdrift is probably "$\\mathcal{K}_2$ [a.u.]" 
            gbdrift0 is $\\mathbf{B} \\times \\nabla\\mathbf{B} \\cdot \\nabla \\psi$
            bmag is "$B$ [a.u.]" '''

    # Update the progress bar of the GUI 
    if Progress: Progress.move(0,"Plot geometry versus z.")
     
    # Set the labels, and change the standard labels depending on the units
    label = determine_labels(x_label, y_label, x_quantity, y_quantity, title, units, normalize)
    load_plotbox2d(x_label=label['x'], y_label=label['y'], title=label["title"], ax=ax) 

    # Keep track of the labels that are already used and save the axis limits
    labels = []; xlims=[0,0]; ylims=[0,0]; 
    
    # Get the total number of varied values, to determine whether we need to use markers
    # to diferentiate between different simulations or whether the line colors are sufficient
    number_ofVariedValues = 0
    experiments = get_experiments(research, experiment_id)
    for experiment in experiments:
        number_ofVariedValues = np.max([number_ofVariedValues,len(experiment.variedValues)])
                
    # First plot the lines with the same color for each experiment
    # Then plot the data points with the same color for each varied_value
    # Keep track of the labels to fix the order in the legend
    labels_lines = []; labels_points = []
     
    # Iterate over the experiments
    for experiment in experiments:
            
        # Iterate over the simulations
        simulations = get_simulations(experiment, simulation_id)
        for simulation in simulations:
                     
            # Get the data
            style    = '-'
            vec_z    = getattr(simulation, "B_"+x_quantity)
            vec_geo  = getattr(simulation, "B_"+y_quantity)
            
            # Get the real accurate data 
            if "LHD"  in str(simulation.input_files[0]): wout_path = "/home/hanne/CIEMAT/MAGNETICFIELDS/LHD/wout_lhd0.nc"
            if "TJII" in str(simulation.input_files[0]): wout_path = "/home/hanne/CIEMAT/MAGNETICFIELDS/TJII/wout_tjii_edi.nc"
            if "W7X"  in str(simulation.input_files[0]): wout_path = "/home/hanne/CIEMAT/MAGNETICFIELDS/W7X/wout_w7xr003.nc"
            if "NCSX" in str(simulation.input_files[0]): wout_path = "/home/hanne/CIEMAT/MAGNETICFIELDS/NCSX/wout_ncsx.nc"
            if "CBC"  in str(simulation.input_files[0]): wout_path = "/home/hanne/CIEMAT/MAGNETICFIELDS/CBC/wout_tok_cbc_2T.nc"  
            geo = calculate_geometricQuantitiesSTELLA(wout_path, simulation.nfield_periods, 48, rho=simulation.rho, nperiod=1)
            vec_zReal = geo[x_quantity]
            vec_geoReal = geo[y_quantity][0,:]
            
            # Normalize the y-quantity to min and max is between 0 and 1
            if normalize:
                minimum = np.nanmin(np.abs(vec_geo))
                maximum = np.nanmax(np.abs(vec_geo-minimum))
                minimum = minimum if maximum!=0.0 else 0.0
                maximum = maximum if maximum!=0.0 else 1.0
                vec_geo = (vec_geo-minimum)/maximum
                vec_z = vec_z/np.max(vec_z)*np.pi 
                
                minimum = np.nanmin(np.abs(vec_geoReal))
                maximum = np.nanmax(np.abs(vec_geoReal-minimum))
                minimum = minimum if maximum!=0.0 else 0.0
                maximum = maximum if maximum!=0.0 else 1.0
                vec_geoReal = (vec_geoReal-minimum)/maximum
                vec_zReal = vec_zReal/np.max(vec_zReal)*np.pi 
             
            # Add a label to the legend for the first simulation of this experiment
            label_exp = experiment.line_label if simulations.index(simulation) == 0 else ""
            
            # Add a label to the legend for the first simulation of this experiment
            label_var = simulation.marker_label; lw=2
            label_var = label_var.replace("_"," ") if "$" not in label_var else label_var
            label_var = "" if (label_var in labels) else label_var; labels.append(label_var)  
                    
            # First plot the lines with the same color for each experiment
            # Only if there are multiple experiments or if there is only one varied value
            if (len(experiments) > 1) or (number_ofVariedValues==1):
                lines = {'linewidth' : 3, 'color' : experiment.line_color}; lw = 0 #'linestyle' : style, 
                if showOnlyMarkers: lines['linewidth'] = 0  
                if label_exp!="": labels_lines.append(label_exp)
                ax.plot(vec_z, vec_geo, label=label_exp, **lines, linestyle="-") 
                ax.plot(vec_zReal, vec_geoReal, label=label_exp, **lines, linestyle=":") 
            
            # Then plot the data points with the same color for each varied_value
            # If there is only one experiment, add lines in the same color as the marker 
            if not number_ofVariedValues==1 or showOnlyMarkers:
                if showOnlyMarkers: lines['linewidth'] = 0; lw=0
                if lw!=0: lines = {'linewidth' : lw, 'color' : simulation.marker_color}
                if lw!=0: markers = {}
                if lw==0: lines = {'linewidth' : lw}
                if lw==0: markers = {'ms' : 6, 'marker' : simulation.marker_style, 'mec' : simulation.marker_color, 'mfc' : simulation.marker_color}
                if label_var!="": labels_points.append(label_var)
                ax.plot(vec_z, vec_geo, label=label_var, **markers, **lines, linestyle="-") 
                ax.plot(vec_zReal, vec_geoReal, label=label_var, **markers, **lines, linestyle=":") 
            
            # Keep track of the axis limits 
            xlims = [np.nanmin([np.nanmin(vec_z), xlims[0]]), np.nanmax([np.nanmax(vec_z), xlims[1]])]
            ylims = [np.nanmin([np.nanmin(vec_geo), ylims[0]]), np.nanmax([np.nanmax(vec_geo), ylims[1]])]
                    
    # Change appearance plot 
    ax.autoscale() 
    ax.set_xlim(xmin=0)
    if tick_style=='sci': ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    if tick_style==None:  ax.ticklabel_format(useOffset=False)
    if log == True: ax.set_yscale('log')
    
    # Add the legend 
    handles, labels = ax.get_legend_handles_labels()  
    labels_points = sort_listByNumbers(labels_points) 
    labels_order = labels_lines + labels_points 
    labels_order = [label for label in labels_order if label in labels]
    labels_index = [labels.index(label) for label in labels_order] 
    labels  = [labels[i] for i in labels_index]
    handles = [handles[i] for i in labels_index] 
    if len(labels_order)>1:
        ax.legend(handles, labels, labelspacing=0.0, prop={'size':16}, handlelength=1.5)

    # Rescale the axis
    if x_range==None: x_range = xlims
    if y_range==None: y_range = ylims 
    ax.set_xlim(x_range) 
    ax.set_ylim(y_range)  
    
    # Show the figure 
    if show_figure: plt.show()
    if True: return
    
#################################################################
#                        METHODS
#################################################################

def determine_labels(x_label, y_label, x_quantity, y_quantity, title, units, normalize):
    ''' Set the labels and titels for the (omega, ki) or (gamma, ki) plot '''
    
    # Save the labels in one dictionary
    label = {'x' : x_label, 'y' : y_label, 'title' : title}

    # Determine the label of the x-axis
    if label["x"] is None:
        if x_quantity=="z": 
            if units=="normalized": label["x"]  = '$z/a$' 
            if units=="SI units":   label["x"]  = '$z$ [m]'
        if x_quantity=="zeta":  
            if units=="normalized": label["x"]  = '$\\zeta$' 
            if units=="SI units":   label["x"]  = '$\\zeta$' 
        if x_quantity=="pol":   
            if units=="normalized": label["x"] = 'Poloidal turns' 
            if units=="SI units":   label["x"] = 'Poloidal turns' 
        if x_quantity=="tor":    
            if units=="normalized": label["x"] = 'Toroidal turns' 
            if units=="SI units":   label["x"] = 'Toroidal turns' 
            
    # Determine the label of the y-axis        
    if label["y"] is None:   
        label["y"] = y_quantity.replace("_", "")
        if y_quantity=="gbdrift": 
            label["y"]  = "$-\\mathcal{K}_2$ [a.u.]"   
        if y_quantity=="gbdrift0": 
            label["y"]  = "$\\mathbf{B} \\times \\nabla\\mathbf{B} \\cdot \\nabla \\psi$"  
        if y_quantity=="bmag": 
            label["y"]  = "$B$ [a.u.]"

    return label
    