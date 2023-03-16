#!/usr/bin/python3 
import h5py 
import pathlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.interpolate as interpolate 


#===============================================================================
#                            READ A PROFILES FILE                              #
#===============================================================================

def calculate_plasmaProfiles(file, Ti_diag='CXRS', rho0=0.7, plot=True, verbose=False):
    """ Calculate the plama profiles, and the plasma gradients at a specific <rho>
    based on a given Profiles.h5 file. A dictionary is returned containing the
    profile data at <rho> with <ni_norm> = 1; <Ti_norm> = 1; <ne_norm> = ne/ni_ref
    and <Te_norm> = Te/Ti_ref which can be written directly to a stella input file. 
    
    The ion temperature can come from 'XICS' or 'CXRS'. """

    # Open the h5 file containing the profile data
    file = pathlib.Path(file)
    h5_file = h5py.File(file, 'r') 
    print_overviewH5File(file, h5_file, verbose) 
        
    # Radial location
    rho = h5_file['Te']['fit'][0,:] 
       
    # Collect the temperature data 
    Te = h5_file['Te']['fit'][1,:]
    Ti = h5_file['Ti']['CXRS_fit'][1,:] if (Ti_diag == 'CXRS') else h5_file['Ti']['XICS_fit'][1,:]
    dTidrho = h5_file['Ti']['CXRS_fit'][2,:] if (Ti_diag == 'CXRS') else h5_file['Ti']['XICS_fit'][2,:]
    dTedrho = h5_file['Te']['fit'][2,:]
    
    # Make the density profiles of the electrons match that of the ions
    ne = ni = h5_file['ne20']['fit'][1,:] 
    dnedrho = dnidrho = h5_file['ne20']['fit'][2,:] 
    
    # Density and temperature ratios
    nine = ni/ne  
    tite = Ti/Te
        
    # Profile gradients: a/L_T= -a/T * dT/dr = -1/T * dT/drho
    tiprim = -dTidrho/Ti; teprim = -dTedrho/Te
    fiprim = -dnidrho/ni; feprim = -dnedrho/ne
    
    # Calulate the profile data at the given radial location
    nine, tite, ni, Ti, ne, Te, tiprim, fiprim, teprim, feprim, ni_norm, Ti_norm, ne_norm, Te_norm = calculate_parametersAtChosenRadialLocation(rho, rho0, nine, tite, ni, ne, Ti, Te, fiprim, feprim, tiprim, teprim)

    # Print the data at the given radial location
    print_overviewDataAtChosenRadialLocation(rho0, verbose, nine, tite, ni, Ti, ne, Te, tiprim, fiprim, teprim, feprim, ni_norm, Ti_norm, ne_norm, Te_norm)

    # Plot the profiles
    plot_profiles(plot, h5_file)
    
    # Return the profile data
    return {"nine" : nine, "tite" : tite, "ni" : ni, "Ti" : Ti, "ne" : ne, "Te" : Te, 
            "tiprim" : tiprim, "fiprim" : fiprim, "teprim" : teprim, "feprim" : feprim, 
            "ni_norm" : ni_norm, "Ti_norm" : Ti_norm, "ne_norm" : ne_norm, "Te_norm" : Te_norm}

#------------------------
def print_overviewH5File(file, h5_file, verbose, length=60):
    if verbose:
        print("") 
        print("     ","".center(length,"=")) 
        print("     ",(" "+file.name+" ").center(length,"="))
        print("     ","".center(length,"=")) 
        for key in list(h5_file.keys()):
            print("\n                  ", key)
            for sub_key in list(h5_file[key].keys()):
                print("                     - ", sub_key) 
        print("")  
    return 

#------------------------
def calculate_parametersAtChosenRadialLocation(rho, rho0, nine, tite, ni, ne, Ti, Te, fiprim, feprim, tiprim, teprim):
    nine = interpolate_toChosenRadialLocation(rho,nine,rho0)
    tite = interpolate_toChosenRadialLocation(rho,tite,rho0)
    ni = interpolate_toChosenRadialLocation(rho,ni,rho0)*10
    Ti = interpolate_toChosenRadialLocation(rho,Ti,rho0)
    tiprim = interpolate_toChosenRadialLocation(rho,tiprim,rho0)
    fiprim = interpolate_toChosenRadialLocation(rho,fiprim,rho0)
    ne = interpolate_toChosenRadialLocation(rho,ne,rho0)*10
    Te = interpolate_toChosenRadialLocation(rho,Te,rho0)
    teprim = interpolate_toChosenRadialLocation(rho,teprim,rho0)
    feprim = interpolate_toChosenRadialLocation(rho,feprim,rho0)
    ni_ref = ni
    Ti_ref = Ti
    ni_norm = 1
    Ti_norm = 1
    ne_norm = ne/ni_ref
    Te_norm = Te/Ti_ref
    return nine, tite, ni, Ti, ne, Te, tiprim, fiprim, teprim, feprim, ni_norm, Ti_norm, ne_norm, Te_norm
    
#------------------------
def print_overviewDataAtChosenRadialLocation(rho0, verbose, nine, tite, ni, Ti, ne, Te, tiprim, fiprim, teprim, feprim, ni_norm, Ti_norm, ne_norm, Te_norm, length=60):
    if rho0!=None and verbose:
        print(""); space = "         "
        print("     ","".center(length,"=")) 
        print("     ",(" Stella input data at rho = "+str(rho0)+" ").center(length,"="))
        print("     ","".center(length,"="),"\n") 
        print(space, "&vmec_parameters")
        print(space, "torflux = "+str("{0:<15}".format(rho0*rho0))+'\n')
        print(space, "&parameters")
        print(space, "nine = "  +str("{0:<15}".format(nine)))
        print(space, "tite = "  +str("{0:<15}".format(tite))+'\n')
        print(space, "&species_parameters")
        print(space, "dens = "  +str("{0:<15}".format(ni)))
        print(space, "temp = "  +str("{0:<15}".format(Ti)))
        print(space, "tprim = " +str("{0:<15}".format(tiprim)))
        print(space, "fprim = " +str("{0:<15}".format(fiprim))+'\n')
        print(space, "&species_parameters_2")        
        print(space, "dens = "  +str("{0:<15}".format(ne))) 
        print(space, "temp = "  +str("{0:<15}".format(Te)))
        print(space, "tprim = " +str("{0:<15}".format(teprim)))
        print(space, "fprim = " +str("{0:<15}".format(feprim))+'\n')
        print(space, "&species_parameters (with densi = tempi = 1)")
        print(space, "dens = "  +str("{0:<15}".format(ni_norm)))
        print(space, "temp = "  +str("{0:<15}".format(Ti_norm)))
        print(space, "tprim = " +str("{0:<15}".format(tiprim)))
        print(space, "fprim = " +str("{0:<15}".format(fiprim))+'\n')
        print(space, "&species_parameters_2 (with densi = tempi = 1)")        
        print(space, "dens = "  +str("{0:<15}".format(ne_norm))) 
        print(space, "temp = "  +str("{0:<15}".format(Te_norm)))
        print(space, "tprim = " +str("{0:<15}".format(teprim)))
        print(space, "fprim = " +str("{0:<15}".format(tiprim))+'\n')
    return

#------------------------   
def interpolate_toChosenRadialLocation(x, y, x0, der=0): 
    tck = interpolate.splrep(x, y, s=0)
    value_y = interpolate.splev(x0,tck,der=der)
    return value_y

#------------------------   
def plot_profiles(plot, h5_file):
    if plot == True: 
        
        # Figure 
        from stellapy.plot.utils.style.load_styleFigures import load_styleFigures
        load_styleFigures() 
        fig = plt.figure(figsize=(18, 5))
        grid_specifications = gridspec.GridSpec(1, 3)
        grid_specifications.update(top=0.95, left=0.05, right=0.98, bottom=0.15, wspace=0.25)
        ax1 = plt.subplot(grid_specifications[0])   
        ax2 = plt.subplot(grid_specifications[1])   
        ax3 = plt.subplot(grid_specifications[2])   
        
        # Temperature
        ax1.set_xlabel('$r/a$'); ax1.set_ylabel('$T_{e}, T_{i}$ [keV]')  
        ax1.set_xlim([0,1]); ax1.set_ylim([0,3])  
        args = {"linewidth" : 2.0,  "markersize" : 10, "markeredgewidth" : 2.0, "mfc" : 'white'}
        args_XICS = {"label" : '$T_{i}^{\\mathrm{XICS}}$', "fmt" : 's', "color" : 'navy', "mec" : 'navy', **args}
        args_CXRS = {"label" : '$T_{i}^{\\mathrm{CXRS}}$', "fmt" : 'o', "color" : 'cornflowerblue', "mec" : 'cornflowerblue', **args}
        args_Te = {"label" : '$T_{e}^{\\mathrm{TS}}$', "fmt" : '^', "color" : 'crimson', "mec" : 'crimson', **args}
        if 'XICS_data' in h5_file['Ti'].keys():
            Ti_data = h5_file['Ti']['XICS_data']; Ti_fit = h5_file['Ti']['XICS_fit']
            ax1.errorbar(Ti_data[0,:], Ti_data[1,:], Ti_data[2,:], **args_XICS)
            ax1.plot(Ti_fit[0,:], Ti_fit[1,:], '-', linewidth=3.0, color='navy')
        if 'CXRS_data' in h5_file['Ti'].keys():
            Ti_data = h5_file['Ti']['CXRS_data']; Ti_fit = h5_file['Ti']['CXRS_fit']
            ax1.errorbar(Ti_data[0,:], Ti_data[1,:], Ti_data[2,:], **args_CXRS)
            ax1.plot(Ti_fit[0,:], Ti_fit[1,:], '-', linewidth=3.0, color='cornflowerblue')
        Te_data = h5_file['Te']['data']; Te_fit = h5_file['Te']['fit']
        ax1.errorbar(Te_data[0,:], Te_data[1,:], Te_data[2,:], **args_Te)
        ax1.plot(Te_fit[0,:], Te_fit[1,:], '-', linewidth=3.0, color='crimson')
        ax1.legend(loc='best', labelspacing=0.1, prop={'size': 20})
         
        # Density
        args = {"label" : '$n_{e}^{\\mathrm{TS}}$', "color" : 'limegreen', "mec" : 'limegreen', "linewidth" : 2.0,  "markersize" : 10, "markeredgewidth" : 2.0, "mfc" : 'white', "fmt" : 's'}
        ax2.set_xlabel('$r/a$'); ax2.set_ylabel('$n_e$ [$10^{20}$ m$^{-3}$]')  
        ax2.set_xlim([0,1]); ax2.set_ylim([0,1])   
        ne_data = h5_file['ne20']['data']; ne_fit = h5_file['ne20']['fit']
        ax2.errorbar(ne_data[0,:], ne_data[1,:], ne_data[2,:], **args)
        ax2.plot(ne_fit[0,:], ne_fit[1,:], '-', linewidth=3.0, color='limegreen')
        ax2.legend(loc='best', labelspacing=0.1, prop={'size':20})
         
        # Gradient length-scales
        args = {"linewidth" : 3.0,  "markersize" : 10, "mfc" : 'white', "markeredgewidth" : 2.0}
        ax3.set_xlabel('$r/a$'); ax3.set_ylabel('$a/L_{X}$')  
        ax3.set_xlim([0,1]); ax3.set_ylim([0,10])    
        ne = h5_file['ne20']['fit']; Te = h5_file['Te']['fit']
        ax3.plot(ne[0,:], -ne[2,:]/ne[1,:], "-s", label='$X=n_{e}^{\\mathrm{TS}}$', color='limegreen', **args)
        ax3.plot(Te[0,:], -Te[2,:]/Te[1,:], '-^', color='crimson', label='$X=T_{e}^{\\mathrm{TS}}$', **args)
        if 'XICS_data' in h5_file['Ti'].keys():
            Ti = h5_file['Ti']['XICS_fit']
            ax3.plot(Ti[0,:], -Ti[2,:]/Ti[1,:], '-o', color='navy', label='$X=T_{i}^{\\mathrm{XICS}}$', **args)
        if 'CXRS_data' in h5_file['Ti'].keys():
            Ti = h5_file['Ti']['CXRS_fit']
            ax3.plot(Ti[0,:], -Ti[2,:]/Ti[1,:], '-o', color='cornflowerblue', label='$X=T_{i}^{\\mathrm{CXRS}}$', **args)
        ax3.legend(loc='best', labelspacing=0.1, prop={'size':20})
        fig.set_tight_layout(False); plt.show()
        figure_folder = "/home/hanne/Pictures/"
        name = "Profiles.png"
        fig.savefig(figure_folder+name, format='png', dpi=300)
        return
    
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__": 
    file = "/home/hanne/CIEMAT/RESEARCH/PARTICLE_FLUXES_PPT/Profiles/Profiles_20180920.017_from_2.00_to_6.00s.h5"
    calculate_plasmaProfiles(file, rho0=0.7, verbose=True, plot=False)
    