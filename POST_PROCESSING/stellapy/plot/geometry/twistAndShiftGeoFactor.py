"""

#===============================================================================
#              Plot the twist and shift geometry factor along zeta             #
#===============================================================================

Plot the twist and shift geometry factor along zeta. If the boundary_option is
chosen to be 'stellarator', then the local shear is used to determine dkx. Ideally,
one chooses an <nfield_periods> such that <twist_and_shift_geo_factor>=1. In order
to find this <nfield_periods>, start with e.g. 1.2 poloidal turns and find the
nfield_periods where <twist_and_shift_geo_factor>=1 close to 1 poloidal turn. 
Then rerun the script with the chosen <nfield_periods> to confirm that the 
<twist_and_shift_geo_factor> has the desired value, the interpolation errors 
are quite big since the <twist_and_shift_geo_factor> can vary very quickly along 
zeta, therefore, it's important to iterate this process.  

Hanne Thienpondt
26/10/2022

"""

#!/usr/bin/python3 
import pathlib
import sys, os
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.calculations.calculate_twistAndShiftGeoFactor import calculate_twistAndShiftGeoFactor 
from stellapy.plot.utils.labels.standardLabels import standardLabels
from stellapy.plot.utils.style.create_figure import create_figure 

#===============================================================================
#              Plot the twist and shift geometry factor along zeta             #
#===============================================================================
 
def plot_twistAndShiftGeoFactor(
        vmec_path = None,\
        poloidal_turns = None,\
        nfield_periods = 9.47408308356441,\
        nzed = 128,\
        rho = 0.25,\
        verbose = True): 
    
    # Calculate the twist and shift geo factor
    zeta, twist_and_shift_geo_factor, poloidal_turns = calculate_twistAndShiftGeoFactor(
        vmec_path = vmec_path,\
        poloidal_turns = poloidal_turns,\
        nfield_periods = nfield_periods,\
        nzed=nzed,\
        rho=rho,\
        verbose=verbose) 
    
    # Interpolate
    zeta_new = np.linspace(zeta[0], zeta[-1], num=nzed*1000, endpoint=True)
    f = interp1d(zeta, twist_and_shift_geo_factor, kind='cubic')
    twist_and_shift_geo_factor_new = f(zeta_new)
    
    # Calculate where we reach one poloidal turn, and where twist_and_shift_geo_factor is one or zero
    zeta_right_oneturn = zeta[-1]/poloidal_turns
    idx_zero = np.argwhere(np.diff(np.sign(twist_and_shift_geo_factor_new - [0]*len(zeta_new)))).flatten()
    idx_one = np.argwhere(np.diff(np.sign(np.abs(twist_and_shift_geo_factor_new) - [1]*len(zeta_new)))).flatten() 
     
    # Plot the twist and shift geo factor   
    ax = create_figure(title="Twist and shift geometry factor")  
    ax.plot(zeta, twist_and_shift_geo_factor, color='black', lw=0, marker="o")  
    ax.plot(zeta_new, twist_and_shift_geo_factor_new, color='black', lw=2)    
    
    # Crossing with zero and one
    ax.axhline(y=1, color='navy', linestyle='-')
    ax.axhline(y=-1, color='navy', linestyle='-')
    ax.axhline(y=0, color='maroon', linestyle='-')
    
    # Highlight the crossings
    ax.axvspan(-zeta_right_oneturn, zeta_right_oneturn, alpha=0.5, color='gray', label="One poloidal turn")
    ax.plot(zeta_new[idx_zero], twist_and_shift_geo_factor_new[idx_zero], marker="x", color="maroon", lw=0)
    ax.plot(zeta_new[idx_one], twist_and_shift_geo_factor_new[idx_one], marker="x", color="blue", lw=0, ms=10, mew=2)
    
    # Appearance
    ax.set_xlabel(standardLabels["normalized"]["zeta"])
    ax.set_ylabel("twist and shift geometry factor")
    ax.set_xlim([np.min(zeta), np.max(zeta)]) 
    ax.set_ylim([np.min(twist_and_shift_geo_factor)*1.1, np.max(twist_and_shift_geo_factor)*1.1])  
    ax.legend(handlelength=1, handletextpad=1, prop={'size':20})
    mpl.rcParams["savefig.directory"] = pathlib.Path(vmec_path).parent 
    plt.show() 
    return
 
#===============================================================================
#                              Main script
#===============================================================================
if __name__ == "__main__": 
    plot_twistAndShiftGeoFactor()
