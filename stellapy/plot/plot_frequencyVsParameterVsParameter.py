
#================================================================
# Plot omega(k) or gamma(k)
#================================================================

# Load the modules
import numpy as np
import configparser, pathlib, os
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Personal modules
from stellapy.plot.utils import load_plotbox2d  
from stellapy.utils.decorators import verbose_wrapper
from stellapy.simulations.utils.get_simulation import get_simulation
from stellapy.simulations.utils.get_experiment import get_experiment
from stellapy.plot.utils import get_axisOfScan  

#===========================
# plot omega(k) or gamma(k)
#===========================

@verbose_wrapper
def plot_frequencyVsParameterVsParameter(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
                parameter_knob1="vmec_parameters",\
                parameter_key1="rho",\
                parameter_knob2="vmec_parameters",\
                parameter_key2="rho",\
            # Specify data range
                z_quantity='gamma',\
                x_range=None,\
                y_range=None,\
            # Details of the modes
                kx_range=-0.0, \
                ky_range=[0,100], \
                k_value=9.0, \
                lineardata="average",\
            # Save time by reading the data elsewhere
                parameters1=None,\
                parameters2=None,\
                gamma_mostUnstableMode=None,\
                omega_mostUnstableMode=None,\
                ky_mostUnstableMode=None,\
                modes=[0,1,2],\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # For the GUI the figure object already exists 
                show_figure = False,\
                ax=None, \
                Progress=None,\
                root=None,\
            # Extra options
                interpolate=False,\
                step=4):
    ''' Plot the frequency/omega or growthrate/gamma versus the wave number. 
    
    Parameters
    ----------
    research: object
        Contains the selected experiments and simulations.
    experiment_id, simulation_id: str
        Is used to select which simulations/experiments to plot.
    kvalue: {"max", integer}
        Determines which modes to plot.
    units: {"normalized", "SI units"}
        Determines the units of the data.
    lineardata: {"average", "last"}
        Either plot the last time value of omega/gamma or the average
        over the last 10% of time.
    '''


    # Update the progress bar of the GUI
    if Progress: Progress.move(0,"Plot growthrate versus parameter.")
         
    # Decide whether to scan modes along the x-axis or y-axis
    scan, k_fixed, k_range = get_axisOfScan(kx_range, ky_range, root) 
    if scan == False: return 
    
    # Create the figure
    if ax is None:
        fig = plt.figure(figsize=(18, 9)); show_figure = True
        grid_specifications = gridspec.GridSpec(1, 1, figure=fig)
        grid_specifications.update(top=0.95, left=0.09, right=0.82, bottom=0.12)
        ax = plt.subplot(grid_specifications)
        
    # Set the labels, and change the standard labels depending on the units
    label = determine_labels(x_label, y_label, kx_range, ky_range, title,  "normalized", k_value, parameter_key1, parameter_key2, z_quantity)
    load_plotbox2d(x_label=label['x'], y_label=label['y'], title=label["title"], ax=ax) 

    # Keep track of the labels that are already used and save the axis limits
    xlims=[0,0]; ylims=[0,0]; 
    
    #==================
    # GET THE DATA
    #==================
    if (parameters1 is None\
    or parameters2 is None\
    or ky_mostUnstableMode is None\
    or gamma_mostUnstableMode is None\
    or omega_mostUnstableMode is None):
    
        # Initiate the parameters
        parameters1 = []
        parameters2 = []
        
        # Iterate over the experiments
        experiments = get_experiment(research, experiment_id)
        for experiment in experiments:
            
            # Get the requested parameters
            simulations = get_simulation(experiment, simulation_id)
            for simulation in simulations:
                parameters1.append(simulation.inputParameters[parameter_knob1][parameter_key1])
                parameters2.append(simulation.inputParameters[parameter_knob2][parameter_key2])
            
        # Sort the parameters
        parameters1 = sorted(list(set(parameters1)))
        parameters2 = sorted(list(set(parameters2)))
            
        # Now initiate the z data: we can plot growth rate, frequency and the ky
        if k_value == "max":
            gamma_mostUnstableMode = np.empty((len(parameters1), len(parameters2))); gamma_mostUnstableMode[:, :] = np.NaN 
            omega_mostUnstableMode = np.empty((len(parameters1), len(parameters2))); omega_mostUnstableMode[:, :] = np.NaN 
            ky_mostUnstableMode    = np.empty((len(parameters1), len(parameters2))); ky_mostUnstableMode[:, :] = np.NaN 
            ky_secondUnstableMode  = np.empty((len(parameters1), len(parameters2))); ky_secondUnstableMode[:, :] = np.NaN 
        if k_value != "max":
            gamma_atSpecificK      = np.empty((len(parameters1), len(parameters2))); gamma_atSpecificK[:, :] = np.NaN 
            omega_atSpecificK      = np.empty((len(parameters1), len(parameters2))); omega_atSpecificK[:, :] = np.NaN 
            ky_atSpecificK         = np.empty((len(parameters1), len(parameters2))); ky_atSpecificK[:, :] = k_value
            
        # Add the (0,0) point since we can't simulate it
        if 0 in parameters1 and 0 in parameters2:
            if k_value == "max":
                gamma_mostUnstableMode[0,0] = 0
                omega_mostUnstableMode[0,0] = 0
                ky_mostUnstableMode[0,0] = 0
                ky_secondUnstableMode[0,0] = 0
            if k_value != "max":
                gamma_atSpecificK[0,0] = 0
                omega_atSpecificK[0,0] = 0
                ky_atSpecificK[0,0] = k_value
                
        # Iterate over the experiments
        experiments = get_experiment(research, experiment_id)
        for experiment in experiments:
                
            # Iterate over the simulations
            simulations = get_simulation(experiment, simulation_id)
            for simulation in simulations:
                
                # Get the values of the parameters
                parameter1 = simulation.inputParameters[parameter_knob1][parameter_key1]
                parameter2 = simulation.inputParameters[parameter_knob2][parameter_key2]
                
                # Get the index of the simulation and the parameters 
                i = simulations.index(simulation)
                a = parameters1.index(parameter1)
                b = parameters2.index(parameter2)
                
                # Update the progress bar of the GUI 
                if Progress: Progress.move(i/len(simulations)*100,"Doing analysis ("+str(i)+"/"+str(len(simulations))+")")
                    
                # Get the modes of this simulation that need to be plotted and their indices in the (kx,ky) matrixes
                modes  = simulation.get_modesForAOneDimensionalScan(scan, k_fixed, k_range, plotted_modes="unstable") 
                i_kxky = simulation.get_indicesOfModes(scan, k_fixed, modes)     
                
                # Get the omega or gamma at the last time value
                if lineardata=="average": 
                    omega_data = simulation.omega_avg[i_kxky]
                    gamma_data = simulation.gamma_avg[i_kxky]
                if lineardata=="last":    
                    omega_data = simulation.omega_last[i_kxky]
                    gamma_data = simulation.gamma_last[i_kxky]
        
                # Get the maximum of omega/gamma or omega/gamma at a specific ky
                if k_value == "max":
                    #try: 
                    index = np.nanargmax(gamma_data)
                    gamma_mostUnstableMode[a,b] = np.nanmax(gamma_data) 
                    omega_mostUnstableMode[a,b] = omega_data[index] 
                    ky_mostUnstableMode[a,b] = modes[index] 
                    
                    # Calculate the second maxima directly
                    x_data = list(modes); y_data = list(gamma_data)
                    maxima = [y for y in y_data[1:-1] if (y>y_data[y_data.index(y)-1] and y>y_data[y_data.index(y)+1]) ]
                    modes =  [x_data[y_data.index(m)] for m in maxima]
                    sorted_indexes = list(np.array(maxima).argsort())
                    sorted_indexes.reverse()
                    maxima = [maxima[i] for i in sorted_indexes]
                    modes  = [modes[i]  for i in sorted_indexes]
                    print()
                    print(modes)
                    print(maxima)
                    if len(modes) > 1: 
                        ky_secondUnstableMode[a,b] = modes[1] 
    
                    #except: pass
                    
                if k_value != "max":
                    mask = np.isin(modes, k_value)
                    if k_value in modes: gamma_atSpecificK[a,b] = gamma_data[mask]
                    if k_value in modes: omega_atSpecificK[a,b] = omega_data[mask]
                    else: pass


    #==================
    # PLOT THE DATA
    #==================
    
    # Update the progress bar of the GUI 
    if Progress: Progress.move(50,"Plotting linear map.")
        
    # Plot the data points with the same color for each varied_value 
    if len(modes)>0:
        
        # Get the z-data
        if k_value == "max":
            if z_quantity=="gamma": y_data = gamma_mostUnstableMode
            if z_quantity=="omega": y_data = omega_mostUnstableMode
            if z_quantity=="ky":    y_data = ky_mostUnstableMode
        if k_value != "max":
            if z_quantity=="gamma": y_data = gamma_atSpecificK
            if z_quantity=="omega": y_data = omega_atSpecificK
            if z_quantity=="ky":    y_data = ky_atSpecificK
            
        # Interpolate the data
        if interpolate:     
            # If we ran without the interpolate we manually added an extra parameter
            if len(y_data)==len(parameters1)-1:
                parameters1 = parameters1[0:-1]
                parameters2 = parameters2[0:-1]            
                
            # Replace nans by the value next to them 
            for index in np.argwhere(np.isnan(y_data)): # At ky,kx=0,0 the values are nan's
                if not index[0]==0: 
                    y_data[index[0], index[1]] = y_data[index[0]-1, index[1]]  
                elif not index[0]==len(y_data[:,1])-1:           
                    y_data[index[0], index[1]] = y_data[index[0]+1, index[1]]           
            function = interp2d(parameters1, parameters2, y_data, kind='linear')
            xnew = np.linspace(0, parameters1[-1], int(len(parameters1)+1)*step)
            ynew = np.linspace(0, parameters2[-1], int(len(parameters2)+1)*step)
            y_data = function(xnew,ynew)
            parameters1, parameters2 = np.meshgrid(xnew, ynew)
            parameters1 = parameters1[0,:]
            parameters2 = parameters2[:,0]
        
        # If we don't interpolate, add an extra parameter because the color is shown from x to x+1
        if not interpolate and len(y_data)==len(parameters1):
            parameters1 = parameters1 + [parameters1[-1] + (parameters1[1]-parameters1[0])]
            parameters2 = parameters2 + [parameters2[-1] + (parameters2[1]-parameters2[0])]
            
        # Add an offset to the ticks if we are not interpolating
        offset = -0.5 if not interpolate else 0
        parameters1_temp = [p + offset for p in parameters1]
        parameters2_temp = [p + offset for p in parameters2]
            
        # Plot y_avrg(kx,ky)
        cmap = plt.get_cmap('jet')
        cmap.set_bad(color='black')
        vmax = np.nanmax(np.nanmax(y_data, axis=0))
        vmin = np.nanmin(np.nanmin(y_data, axis=0))
        vmin = 0 if vmin>0 else vmin
        img = ax.pcolormesh(parameters1_temp, parameters2_temp, y_data.T, cmap=cmap, vmin=vmin, vmax=vmax)  
        cbar = plt.colorbar(img, ax=ax)
        
        # Set a label for the colorbar
        if z_quantity=="gamma": cbar.set_label('$\\gamma a/v_{\\mathrm{th},i}$')
        if z_quantity=="omega": cbar.set_label('$\\omega a/v_{\\mathrm{th},i}$')
        if z_quantity=="ky": cbar.set_label('$k_{y}\\rho_i$')
        
        # Keep track of the axis limits
        xlims = [np.nanmin([np.nanmin(parameters1_temp), xlims[0]]), np.nanmax([np.nanmax(parameters1_temp), xlims[1]])]
        ylims = [np.nanmin([np.nanmin(parameters2_temp), ylims[0]]), np.nanmax([np.nanmax(parameters2_temp), ylims[1]])]
            
    # Change appearance plot 
    ax.autoscale()
    
    # Rescale the axis
    if x_range==None: x_range = xlims
    if y_range==None: y_range = ylims
    ax.set_xlim(x_range) 
    ax.set_ylim(y_range) 
    
    # Set the ticks
    if not interpolate:
        plt.xticks(parameters1[:-1])
        plt.yticks(parameters2[:-1])
    
    # Show the figure
    if show_figure: plt.show()
    
    # Return the data so we can save it to speed up processing
    if interpolate: data = None
    if not interpolate:
        
        # Save the data so we can save it to speed up processing
        try: 
            data = {parameter_key1 : parameters1,\
                parameter_key2 : parameters2,\
                "key1"  : parameter_key1,\
                "key2"  : parameter_key2,\
                "gamma" : gamma_mostUnstableMode,\
                "omega" : omega_mostUnstableMode,\
                "ky"    : ky_mostUnstableMode, \
                "ky2"   : ky_secondUnstableMode}
        
            # Write the data to a configuration file
            write_configurationFile(research, data, k_value)
            
        except: 
            data = {parameter_key1 : parameters1,\
                parameter_key2 : parameters2,\
                "key1"  : parameter_key1,\
                "key2"  : parameter_key2,\
                "gamma" : gamma_mostUnstableMode,\
                "omega" : omega_mostUnstableMode,\
                "ky"    : ky_mostUnstableMode}
    
    return cbar, data

#################################################################
#                        METHODS
#################################################################

def write_configurationFile(research, data, k_value):
    ''' Write a configuration file and save it as "linearmap_default.ini". '''
    
    # Create a configuration object
    linearmap_file = configparser.ConfigParser() 
    linearmap_path = pathlib.Path(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])  
    linearmap_path = linearmap_path / "stellapy/config/linearmap_default.ini"
    
    # Rewrite the input_files and modes to make them look orderly
    experiments = [e.id for e in research.experiments]
    experiments = "\t\n" + ("\t\n").join(experiments)
    
    # Get the parameters
    parameters1 = data[data['key1']]
    parameters2 = data[data['key2']]
    
    # Get which values are extracted
    extraction = "most unstable mode" if k_value == "max" else "mode at ky="+str(k_value)
        
    # Make the sections in the configuration file 
    linearmap_file['GENERAL'] = {
        'research' : "default",\
        'experiments' : experiments,\
        'extraction' : extraction,\
        'parameter 1' :  data['key1'],\
        'parameter 2' :  data['key2'],\
        data['key1'] : parameters1[:-1],\
        data['key2'] : parameters2[:-1],}
    linearmap_file["Growth rate"] = {}
    linearmap_file["Frequency"] = {}
    linearmap_file["Ky of the most unstable mode"] = {}
    linearmap_file["Ky of the second maxima"] = {}
    
    # Save the linear map:  
    for i in range(len(parameters1)-1):
        for j in range(len(parameters2)-1):
            parameter1 = parameters1[i];  parameter2 = parameters2[j]; 
            parameters = '(' + str(parameter1) + ", " + str(parameter2) + ')'
            linearmap_file["Growth rate"][parameters] = str(data['gamma'][i,j])
            linearmap_file["Frequency"][parameters]   = str(data['omega'][i,j])
            linearmap_file["Ky of the most unstable mode"][parameters] = str(data['ky'][i,j])
            linearmap_file["Ky of the second maxima"][parameters] = str(data['ky2'][i,j])
    
    # Write the configuration file
    linearmap_file.write(open(linearmap_path, 'w'))
    return

def get_modes(simulation, kx_range, ky_range, plotted_modes):
    ''' Only plot the stable/unstable/converged/unconverged modes along kx/ky. '''
    
    # Get the modes in the simulation
    vec_kx  = simulation.vec_kx
    vec_ky  = simulation.vec_ky
    
    # Get the stable/unstable or converged/unconverged modes
    if isinstance(kx_range, float):  simulation.get_modesForAOneDimensionalScan(scan="ky", fixed_value=kx_range, k_range=ky_range)
    if isinstance(ky_range, float):  simulation.get_modesForAOneDimensionalScan(scan="kx", fixed_value=ky_range, k_range=kx_range)
    if plotted_modes=="stable":      vec_k = simulation.vec_kStable
    if plotted_modes=="unstable":    vec_k = simulation.vec_kUnstable
    if plotted_modes=="converged":   vec_k = simulation.vec_kConverge
    if plotted_modes=="unconverged": vec_k = simulation.vec_kNotConverge
    
    # Get the modes to plot
    if isinstance(kx_range, float): 
        selected_kx = [kx_range]
        selected_ky = [ ky for ky in vec_ky if (ky >= ky_range[0] and ky <= ky_range[1] and ky in vec_k)]
    if isinstance(ky_range, float): 
        selected_kx = [ kx for kx in vec_kx if (kx >= kx_range[0] and kx <= kx_range[1] and kx in vec_k)]
        selected_ky = [ky_range]
        
    # Return the total amount of modes and those to plot
    return vec_kx, vec_ky, selected_kx, selected_ky

#---------------------------------------------
def determine_labels(x_label, y_label, kx_range, ky_range, title,  units, k_value, parameter_key1, parameter_key2, z_quantity):
    ''' Set the labels and titels for the (omega, ki) or (gamma, ki) plot '''
    
    # Save the labels in one dictionary
    label = {'x' : x_label, 'y' : y_label, 'title' : title}
    
    # Set a title
    if k_value=="max" and z_quantity=="gamma": label['title'] = 'Growth rate of the most unstable mode'
    if k_value=="max" and z_quantity=="omega": label['title'] = 'Frequency of the most unstable mode'
    if k_value=="max" and z_quantity=="ky":    label['title'] = 'Wavenumber of the most unstable mode'
    
    # Find out the label of the x-axis
    for i in range(2):
        parameter_key = [parameter_key1, parameter_key2][i]
        label_key = ["x", "y"][i]
        if   parameter_key=="rho":      label[label_key] = "$\\rho$"
        elif parameter_key=="tprim":    label[label_key] = "$a/L_{T_i}$"
        elif parameter_key=="tiprim":   label[label_key] = "$a/L_{T_i}$"
        elif parameter_key=="teprim":   label[label_key] = "$a/L_{T_e}$"
        elif parameter_key=="fprim":    label[label_key] = "$a/L_{n}$"
        elif parameter_key=="delta t":  label[label_key] = "$\\Delta t$"
        elif parameter_key=="delt":     label[label_key] = "$\\Delta t$"
        elif parameter_key=="nmu":      label[label_key] = "$n_{\\mu}$"
        elif parameter_key=="nvgrid":   label[label_key] = "$n_{v}$"
        else:                           label[label_key] = parameter_key

  
    return label

