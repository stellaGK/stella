"""

#===============================================================================
#           Calculate the twist and shift geometry factor along zeta           #
#===============================================================================

In order to use the stellarator symmetric boundary conditions in stella, it is
important to chose the flux tube length such that one obtains a favourable value
of the local shear at the ends of the flux tube. The script plots the twist_and_
shift_geo_factor as a function of zeta, which can be used to choose nfield_periods 
(or poloidal_turns) such that e.g. dkx = dky by picking a length at which twist_
and_shift_geo_factor = 1 since dkx = twist_and_shift_geo_factor*dky. The nzed 
chosen here is arbitrary, a higher nzed takes longer to compute but gives
a better interpolation of the twist_and_shift_geo_factor.


Hanne Thienpondt
30/05/2022

"""

#!/usr/bin/python3 
import os, sys 
import pathlib
import numpy as np 
from scipy.interpolate import interp1d

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.data.geometry.calculate_geometricQuantities import calculate_geometricQuantities  
 
#===============================================================================
#           Calculate the twist and shift geometry factor along zeta           #
#===============================================================================

def calculate_twistAndShiftGeoFactor(
        vmec_path=None,\
        poloidal_turns=None,\
        nfield_periods=None,\
        nzed=128,\
        rho=0.25,\
        verbose=True): 

    # Calculate the geometric quantities
    _, STELLA, VMEC, _ = calculate_geometricQuantities(
        vmec_path, 
        # Radial location
        rho=rho, 
        # Flux tube
        nfield_periods=nfield_periods, 
        poloidal_turns=poloidal_turns,
        nzed=nzed, 
        nperiod=1, 
        # Choose between linked, stellarator, periodic or default (unlinked)
        boundary_option="stellarator",  
        zed_equal_arc=True,
        # If <grad_x_grad_y_end> is smaller than <grad_x_grad_y_zero> then period BC are used 
        grad_x_grad_y_zero=None,
        # Other
        save_as_pickle=True, 
        verbose=False) 

    # Get the data we need
    zeta = STELLA.zeta[0,:]
    twist_and_shift_geo_factor = STELLA.twist_and_shift_geo_fac_full[0,:] 
    
    # Print where the factor crosses through 0 or 1
    if verbose: print_placesWhereTwistAndShiftEqualsZeroAndOne(zeta, twist_and_shift_geo_factor, nzed, VMEC.poloidal_turns, VMEC.nfield_periods)
    return zeta, twist_and_shift_geo_factor, VMEC.poloidal_turns

#-----------------------------------------------------
def print_placesWhereTwistAndShiftEqualsZeroAndOne(zeta, twist_and_shift_geo_factor, nzed, poloidal_turns, nfield_periods):
    
    # Interpolate 
    zeta_new = np.linspace(zeta[0], zeta[-1], num=nzed*1000, endpoint=True)
    f = interp1d(zeta, twist_and_shift_geo_factor, kind='cubic')
    factor_new = f(zeta_new)
    
    # Calculate where we reach one poloidal turn, and where twist_and_shift_geo_factor is one or zero
    zeta_right_oneturn = zeta[-1]/poloidal_turns
    idx_zero = np.argwhere(np.diff(np.sign(factor_new - [0]*len(zeta_new)))).flatten()
    idx_one = np.argwhere(np.diff(np.sign(np.abs(factor_new) - [1]*len(zeta_new)))).flatten() 
    
    # Print this information
    print()
    print('                     ====================== ')
    print('                         Twist and shift     ')
    print('                     ======================\n ') 
    print("{0:>30}".format("nfield_periods"), " ", "{0:<10}".format(nfield_periods))
    print("{0:>30}".format("poloidal_turns"), " ", "{0:<10}".format(poloidal_turns))
    print("{0:>30}".format("twist and shift at edge"), " ", "{0:<10}".format(twist_and_shift_geo_factor[-1]))
    print()
    print('                     ====================== ')
    print('                           Crossing Zero     ')
    print('                     ======================\n ') 
    for idx in idx_zero:
        if zeta_new[idx]>=0:
            print("{0:>30}".format("zeta: "), " ", "{0:<10}".format(zeta_new[idx])) 
    print()
    print('                     ====================== ')
    print('                           Crossing One     ')
    print('                     ======================\n ') 
    for idx in idx_one:
        if zeta_new[idx]>=0:
            print("{0:>30}".format("zeta: "), " ", "{0:<10}".format(zeta_new[idx]))   
    print()
    print('                     ====================== ')
    print('                           Crossing Zero     ')
    print('                     ======================\n ') 
    for idx in idx_zero:
        if zeta_new[idx]>=0:
            print("{0:>30}".format("nfield_periods: "), " ", "{0:<10}".format(zeta_new[idx]/zeta_new[-1]*nfield_periods)) 
    print()
    print('                     ====================== ')
    print('                           Crossing One     ')
    print('                     ======================\n ') 
    for idx in idx_one:
        if zeta_new[idx]>=0:
            print("{0:>30}".format("nfield_periods: "), " ", "{0:<10}".format(zeta_new[idx]/zeta_new[-1]*nfield_periods))  
    try:
        print()
        print('                     ====================== ')
        print('                     Maximum outside 1 turn     ')
        print('                     ======================\n ') 
        factor_outside = factor_new[zeta_new>zeta_right_oneturn]
        zeta_outside = zeta_new[zeta_new>zeta_right_oneturn]
        print("{0:>30}".format("zeta: "), " ", "{0:<10}".format(zeta_outside[np.argmax(np.abs(factor_outside))]))  
        print("{0:>30}".format("nfield_periods: "), " ", "{0:<10}".format(zeta_outside[np.argmax(np.abs(factor_outside))]/zeta_new[-1]*nfield_periods))  
    except: pass
    return
 
#===============================================================================
#                              Main script
#===============================================================================
if __name__ == "__main__": 
    calculate_twistAndShiftGeoFactor()
