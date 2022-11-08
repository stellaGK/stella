
#================================================================
# Plot omega(k) or gamma(k)
#================================================================

# Load the modules
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from scipy.interpolate import interp2d

# Personal modules
from stellapy.GUI.plot.utils import create_figure 
from stellapy.utils.decorators import verbose
from stellapy.simulations.utils.get_simulations import get_simulations
from stellapy.simulations.utils.get_experiments import get_experiments 

#===========================
# plot omega(k) or gamma(k)
#===========================

@verbose
def plot_saturatedFluxVsParameterVsParameter(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
                parameter_knob1="vmec_parameters",\
                parameter_key1="rho",\
                parameter_knob2="vmec_parameters",\
                parameter_key2="rho",\
            # Specify data range
                species=[0],\
                x_quantity='fprim',\
                y_quantity='tprim',\
                z_quantity='qflux',\
                x_range=None,\
                y_range=None,\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # For the GUI the figure object already exists 
                show_figure = False,\
                ax=None, \
                Progress=None,\
            # Extra options
                spectrum_along_axis = "kx",\
                removeZonal = True,\
                log=False,\
                units="normalized",\
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
    
    # Create the figure
    if ax is None:
        fig = plt.figure(figsize=(18, 9)); show_figure = True
        grid_specifications = gridspec.GridSpec(1, 1, figure=fig)
        grid_specifications.update(top=0.95, left=0.09, right=0.82, bottom=0.12)
        ax = plt.subplot(grid_specifications)
        
    # Set the labels, and change the standard labels depending on the units 
    z2_quantity = "para" if z_quantity=="qflux" else spectrum_along_axis
    z2_quantity = x_quantity+"NZ" if (removeZonal==True and z_quantity!="qflux") else z2_quantity
    ax, z_quantity = create_figure(ax, units, title, x_quantity, x_label, y_quantity, y_label, z_quantity, None, z2_quantity=z2_quantity) 

    # Keep track of the labels that are already used and save the axis limits
    xlims=[0,0]; ylims=[0,0]; 
    
    #=======================
    # GET THE XY PARAMETERS
    #=======================
    
    # Initiate the parameters
    parameters1 = []
    parameters2 = []
            
    # Iterate over the experiments
    experiments = get_experiments(research, experiment_id)
    for experiment in experiments:
        
        # Iterate over the simulations
        simulations = get_simulations(experiment, simulation_id)
        for simulation in simulations:
            
            # Get the requested parameters
            parameters1.append(simulation.inputParameters[parameter_knob1][parameter_key1])
            parameters2.append(simulation.inputParameters[parameter_knob2][parameter_key2])
        
    # Sort the parameters
    parameters1 = sorted(list(set(parameters1)))
    parameters2 = sorted(list(set(parameters2)))
    
    #====================
    # GET THE Z-QUANTITY
    #====================
    
    # Initiate the saturated fluxes
    z_data = np.empty((len(parameters1), len(parameters2))); z_data[:, :] = np.NaN 
    
    # Iterate over the experiments
    experiments = get_experiments(research, experiment_id)
    for experiment in experiments:
        
        # Iterate over the simulations
        simulations = get_simulations(experiment, simulation_id)
        for simulation in simulations:
            
            # Get the values of the parameters
            parameter1 = simulation.inputParameters[parameter_knob1][parameter_key1]
            parameter2 = simulation.inputParameters[parameter_knob2][parameter_key2]
                
            # Get the index of the simulation and the parameters 
            a = parameters1.index(parameter1)
            b = parameters2.index(parameter2)
                
            # Get the saturated flux
            if z_quantity=="qflux":
                z_data[a,b] = simulation.calculate_saturatedFlux()[z_quantity][species[0]]
                
            # Get the height of the heat flux spectra along x or along y
            if z_quantity!="qflux":
                
                # Get the spectrum versus time
                if z_quantity == "phi2_range": data = simulation.phi2_kxky[:,:,:] 
                if z_quantity == "q_range":    data = simulation.qflx_kxky[:,:,:,species[0]] 
            
                # Get the average over a time frame   
                t_last  = simulation.vec_t[(~np.isnan(simulation.vec_t)).sum() - 1]
                t_start = simulation.t_range[species[0]][0]
                t_stop  = np.nanmin([t_last, simulation.t_range[species[0]][1]])
                data = np.nanmean(data[(simulation.vec_t > t_start) & (simulation.vec_t < t_stop),:,:], axis=0) 
        
                # Sum the ky components of phi to get the spectrum along kx 
                if spectrum_along_axis=="kx" and removeZonal==False: 
                    data = data[:,0]+2*np.sum(data[:,1:],axis=1) 
                if spectrum_along_axis=="kx" and removeZonal==True: 
                    data = 2*np.sum(data[:,1:],axis=1) 
                if spectrum_along_axis=="ky" and removeZonal==False: 
                    data = np.sum(data[:,:],axis=0)
                if spectrum_along_axis=="ky" and removeZonal==True: 
                    data = np.sum(data[:,:],axis=0)
                    data[0] = np.NaN
                    
                # Now get the range of values
                z_data[a,b] = np.log10(np.nanmax(data)) - np.log10(np.nanmin(data))
    
    #==================
    # PLOT THE DATA
    #==================
    
    # Update the progress bar of the GUI 
    if Progress: Progress.move(50,"Plotting linear map.")
        
    # Interpolate the data
    if interpolate:     
        
        # If we ran without the interpolate we manually added an extra parameter
        if len(z_data)==len(parameters1)-1:
            parameters1 = parameters1[0:-1]
            parameters2 = parameters2[0:-1]            
            
        # Replace nans by the value next to them 
        for index in np.argwhere(np.isnan(z_data)): # At ky,kx=0,0 the values are nan's
            if not index[0]==0: 
                z_data[index[0], index[1]] = z_data[index[0]-1, index[1]]  
            elif not index[0]==len(z_data[:,1])-1:           
                z_data[index[0], index[1]] = z_data[index[0]+1, index[1]]           
        function = interp2d(parameters1, parameters2, z_data, kind='linear')
        xnew = np.linspace(0, parameters1[-1], int(len(parameters1)+1)*step)
        ynew = np.linspace(0, parameters2[-1], int(len(parameters2)+1)*step)
        z_data = function(xnew,ynew)
        parameters1, parameters2 = np.meshgrid(xnew, ynew)
        parameters1 = parameters1[0,:]
        parameters2 = parameters2[:,0]
    
    # If we don't interpolate, add an extra parameter because the color is shown from x to x+1
    if not interpolate and len(z_data)==len(parameters1):
        parameters1 = parameters1 + [parameters1[-1] + (parameters1[1]-parameters1[0])]
        parameters2 = parameters2 + [parameters2[-1] + (parameters2[1]-parameters2[0])]
        
    # Add an offset to the ticks if we are not interpolating
    offset = -0.5 if not interpolate else 0
    parameters1_temp = [p + offset for p in parameters1]
    parameters2_temp = [p + offset for p in parameters2]

    # If z-axis is logaritmic, remove the zero values
    if log: z_data[z_data < 1.E-25] = np.NaN 
        
    # Define the color map and range
    cmap = plt.get_cmap('jet')
    cmap.set_bad(color='black')
    vmax = np.nanmax(np.nanmax(z_data, axis=0))
    vmin = np.nanmin(np.nanmin(z_data, axis=0))
    vmin = 0 if (vmin>0 and not log) else vmin
    norm = LogNorm(vmin=vmin, vmax=vmax) if log else None
    
    # Plot sat_flux(par1, par2)
    img = ax.pcolormesh(parameters1_temp, parameters2_temp, z_data.T, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)  
    cbar = plt.colorbar(img, ax=ax)
    
    # Set a label for the colorbar 
    if len(species)==1:
        if int(species[0]) == 0: s = "_i" 
        if int(species[0]) == 1: s = "_e" 
    else: s = "_s" 
    if z_quantity=="qflux": z_label = "$Q"+s+"/Q_{gB}$" if units=="normalized" else "$Q"+s+"$ [W/m$^2$]"
    if z_quantity=="vflux": z_label = "$\\Gamma"+s+"/\\Gamma_{gB}$" if units=="normalized" else "$\\Gamma"+s+"$ [W/m$^2$]"
    if z_quantity=="pflux": z_label = "$\\Pi"+s+"/\\Pi_{gB}$" if units=="normalized" else "$\\Pi"+s+"$ [W/m$^2$]"
    if z_quantity=="q_range" and spectrum_along_axis=="ky" and removeZonal==False: z_label = "max( log$_{10}(\\sum_{k_x} Q)-$min(  log$_{10}(\\sum_{k_x} Q)$"
    if z_quantity=="q_range" and spectrum_along_axis=="ky" and removeZonal==True:  z_label = "max( log$_{10}(\\sum_{k_x} Q)-$min(  log$_{10}(\\sum_{k_x} Q)$"
    if z_quantity=="q_range" and spectrum_along_axis=="kx" and removeZonal==False: z_label = "max( log$_{10}(\\sum_{k_y} Q)-$min(  log$_{10}(\\sum_{k_y} Q)$"
    if z_quantity=="q_range" and spectrum_along_axis=="kx" and removeZonal==True:  z_label = "max( log$_{10}(\\sum_{k_y} Q)-$min(  log$_{10}(\\sum_{k_y} Q)$"
    if z_quantity=="phi2_range" and spectrum_along_axis=="ky" and removeZonal==False: z_label = "max( log$_{10}(\\sum_{k_x} |\\langle\\hat\\phi\\rangle|^2)-$min(  log$_{10}(\\sum_{k_x} |\\langle\\hat\\phi\\rangle|^2)$"
    if z_quantity=="phi2_range" and spectrum_along_axis=="ky" and removeZonal==True:  z_label = "max( log$_{10}(\\sum_{k_x} |\\langle\\hat\\phi\\rangle|^2)_{k_y \\neq 0}-$min(  log$_{10}(\\sum_{k_x} |\\langle\\hat\\phi\\rangle|^2)_{k_y \\neq 0}$"
    if z_quantity=="phi2_range" and spectrum_along_axis=="kx" and removeZonal==False: z_label = "max( log$_{10}(\\sum_{k_y} |\\langle\\hat\\phi\\rangle|^2)-$min(  log$_{10}(\\sum_{k_y} |\\langle\\hat\\phi\\rangle|^2)$"
    if z_quantity=="phi2_range" and spectrum_along_axis=="kx" and removeZonal==True:  z_label = "max( log$_{10}(\\sum_{k_y > 0} |\\langle\\hat\\phi\\rangle|^2)-$min(  log$_{10}(\\sum_{k_y > 0} |\\langle\\hat\\phi\\rangle|^2)$" 
    cbar.set_label(z_label)
    
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
        
    # Add a title
    if z_quantity=="qflux": title = "Saturated heat flux"
    if z_quantity=="vflux": title = "Saturated momentum flux"
    if z_quantity=="pflux": title = "Saturated particle flux"
    if z_quantity=="q_range" and spectrum_along_axis=="ky" and removeZonal==False: title = "Height of the spectrum of\n the heat flux along $k_y$"
    if z_quantity=="q_range" and spectrum_along_axis=="ky" and removeZonal==True:  title = "Height of the spectrum of\n the heat flux along $k_y$"
    if z_quantity=="q_range" and spectrum_along_axis=="kx" and removeZonal==False: title = "Height of the spectrum of\n the heat flux along $k_x$"
    if z_quantity=="q_range" and spectrum_along_axis=="kx" and removeZonal==True:  title = "Height of the spectrum of\n the heat flux along $k_x$"
    if z_quantity=="phi2_range" and spectrum_along_axis=="ky" and removeZonal==False: title = "Height of the spectrum of\n the potential squared along $k_y$"
    if z_quantity=="phi2_range" and spectrum_along_axis=="ky" and removeZonal==True:  title = "Height of the spectrum of\n the potential squared along $k_y$ for $k_y \\neq 0$"
    if z_quantity=="phi2_range" and spectrum_along_axis=="kx" and removeZonal==False: title = "Height of the spectrum of\n the potential squared along $k_x$"
    if z_quantity=="phi2_range" and spectrum_along_axis=="kx" and removeZonal==True:  title = "Height of the spectrum of\n the potential squared along $k_x$ for $k_y \\neq 0$"
    ax.set_title(title, fontsize=18)
    
    # Show the figure
    if show_figure: plt.show()
    return cbar

#################################################################
#                        METHODS
#################################################################

def determine_labels(x_label, y_label, title, parameter_key1, parameter_key2, z_quantity):
    ''' Set the labels and titels for the (omega, ki) or (gamma, ki) plot '''
    
    # Save the labels in one dictionary
    label = {'x' : x_label, 'y' : y_label, 'title' : title}
    
    # Set a title
    if z_quantity=="qflux": label['title'] = 'Saturated heat flux' 
    if z_quantity=="vflux": label['title'] = 'Saturated momentum flux' 
    if z_quantity=="pflux": label['title'] = 'Saturated particle flux' 
    
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

