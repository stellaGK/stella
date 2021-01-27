
#================================================================
# Plot |phi^2(t)|
#================================================================

# Load the modules 
import numpy as np
import matplotlib.pyplot as plt
# from .plotbox import plotbox_2d
# import matplotlib.gridspec as gridspec
from stellapy.utils.decorators import verbose_wrapper
# from stellapy.lineardata import calculate_lineardata_for_one_kx_ky
# from stellapy.utils import get_NetcdfFilesOrReducedFilesInside, get_filesinfolder, initiate_nesteddict, get_fullpathoffolder

@verbose_wrapper
# stellapy.plt.phi2_vs_t
def plot_potentialVsTime(folder, kmin=0, kmax=1000, thesis=True, plot_heatflux=False, input_file=None):
    ''' Plot phi^2(t) with phi=sum(phi(kx,ky)) '''

    linear_data = initiate_nesteddict()
    plt.rc('font', size=22)
    # Prepare the figure
    if thesis==False:
        fig = plt.figure(figsize=[16, 12])
        gs1 = gridspec.GridSpec(1, 1)
        gs1.update(wspace=0.025, hspace=0.00, top=0.92, right=0.7, bottom = 0.1)
        ax = plt.subplot(gs1[0]) 
        x_label = "$t\, v_{\mathrm{th}}/a$"
        y_label = "$|\\varphi(t)|^2$"    if plot_heatflux==False else "$Q_s/Q_{gB}$"
        plotbox_2d(x_label=x_label, y_label=y_label, title="", fig_size=(8.5, 7.5), ax=ax);
    if thesis==True:
        fig = plt.figure(figsize=[16, 6])
        gs1 = gridspec.GridSpec(1, 1)
        ax = plt.subplot(gs1[0]) 
        x_label = "$t\, v_{\mathrm{th}}/a$"
        y_label = "$|\\varphi|^2$"     if plot_heatflux==False else "$Q_s/Q_{gB}$"
        plotbox_2d(x_label=x_label, y_label=y_label, fig_size=(8.5, 7.5), ax=ax);

    # Two scanarios based on nonlinear = {True or False}
    input_parameters = read.input_parameters(folder, get_filesinfolder(folder, end=".in")[0])
    nonlinear_run = input_parameters['nonlinear']
    ky_scan = not input_parameters['nonlinear']

    # Read the netcdf files to get potential['time', 'phi2_vs_kxky', 'phi2'] or potential[ky]['time', 'phi2_vs_kxky', 'phi2']
    potential = read.phi(folder, kmin, kmax, input_file)
    keys_shot, keys_rho = read.sorted_shot_and_rho_keys(potential, return_type='keys')
    
    # CASE 1: Nonlinear runs with extra runs in run01, run02, ... folders    
    if nonlinear_run == True:
        for shot in keys_shot:
            for rho in keys_rho[shot]:
                for ky in potential[shot][rho].keys():

                    # Read the data
                    heat_flux   = np.abs(potential[shot][rho][ky]['qflx'])
                    time_q      = potential[shot][rho][ky]['time_fluxes']
                    phi         = potential[shot][rho][ky]['phi2']
                    time        = potential[shot][rho][ky]['time'] 
                    phi_log     = np.log10(np.abs(phi))

                    # Plot phi^2(t)
                    if '348' in folder: color='navy'
                    if '169' in folder: color='crimson'
                    ax.plot(time, phi, lw=3., color=color, label="$|\\varphi|^2$") # plot phi
                    if plot_heatflux==True:     ax.plot(time_q, heat_flux, lw=3., color='black', label="$Q_s/Q_{gB}$") # plot heat_flux
                    print(time_q)
                    # Show phi(kx,ky) in the linear growth
                    try: 
                        index_start = list(phi_log).index( list(filter(lambda i: i > -5.0, phi_log))[0] )
                        index_stop  = list(phi_log).index( list(filter(lambda i: i > -1.0, phi_log))[0] )
                        ax.axvspan(time[index_start], time[index_stop], facecolor='crimson', alpha=0.4)
                    except:
                        print()

                    # Show phi(kx,ky) on the stagnation
                    try:
                        index_start = list(phi_log).index( list(filter(lambda i: i > 1.0, phi_log))[0] )
                        index_stop  = list(time).index( list(filter(lambda i: i > time[index_start]+20, time))[0] )                      
                        ax.axvspan(time[index_start], time[index_stop], facecolor='navy', alpha=0.4)
                    except:
                        print()
                    # Show phi(kx,ky) in steady-state
                    ax.axvspan(400, time[-1], facecolor='green', alpha=0.4)
                    if "004" in folder: ax.axvspan(200, time[-1], facecolor='green', alpha=0.4)

    # CASE 2: kyscan with multiple files for different ky
    if ky_scan == True:
        for shot in keys_shot:
            for rho in keys_rho[shot]:

                # The number of lines defines the colors
                ky_vec = list(potential[shot][rho].keys())
                color  = plt.cm.jet(np.linspace(0,1,len(ky_vec)))
                
                # Plot phi(t)
                for ky in ky_vec:
                    ax.plot(potential[shot][rho][ky]['time'], potential[shot][rho][ky]['phi2'], lw=2, color=color[ky_vec.index(ky)],label=round(ky,1))

                # Check the calculated growth rates
                if input_file:
                    if "348" in input_file: wout_file_path = get_fullpathoffolder(folder) + "/wout_w7xr348+252.h5"; B_field = "348"
                    if "169" in input_file: wout_file_path = get_fullpathoffolder(folder) + "/wout_w7xr169+252.h5"; B_field = "169"
                    wout_parameters = read.wout(wout_file_path)
                    b0 = wout_parameters['b0']
                    linear_data_one_kx_ky = calculate_lineardata_for_one_kx_ky(folder, input_file, b0=b0) 
                    for index in range(len(linear_data_one_kx_ky["ky"][:])):
                        ky    = linear_data_one_kx_ky["ky"][index]
                        last_20percent = int(np.size(potential[shot][rho][ky]['time'])*8/10)
                        gamma = linear_data_one_kx_ky["gamma_unstable"][index] 
                        time  = potential[shot][rho][ky]['time'][last_20percent:-1]
                        phi   = potential[shot][rho][ky]['phi2'][last_20percent:-1]
                        if plot_heatflux==False: ax.plot([time[0], time[-1]], [phi[0], phi[-1]*10**gamma], lw=2, ls=':', color=color[ky_vec.index(ky)])
 
    # Finish the figure
    ax.autoscale()
    if plot_heatflux==True:     ax.set_xlim(xmin=0, xmax=time_q[-1]) 
    ax.set_yscale("log")
    #ax.set_ylim(ymax=100000000000000, ymin=1.E-1)
    ax.set_xlim([0,600])
    if "004" in folder: ax.set_xlim([0,400])
    #ax.legend(loc='best', labelspacing=0.0, ncol=1, shadow=True, prop={'size':24})
    if ky_scan == True:
        ax.legend(loc='upper center', bbox_to_anchor=(1.25, 1), labelspacing=0.0, ncol=2, shadow=True, prop={'size':24},title='$k_y$')
    plt.show()


















