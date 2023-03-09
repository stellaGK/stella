""" 

#===============================================================================
#                      Calculate (nx, ny, nfield_periods)                      #
#===============================================================================

For a given <vmec> or <miller> (defined by <shat> and <q>) at <rho> the values of
<nx> and <ny> will be calculated in order to obtain <kxmax> and <kymax> with <y0>. 
These values depend on <poloidal_turns> or <nfield_periods>, the script is also
useful to convert these two between each other. The <twist_and_shift_geo_factor>
depends on the <boundary_option> = {'linked', 'stellarator', 'unlinked', 'periodic'},
and it determines dkx since dkx = dky * twist_and_shift_geo_factor.

If boundary_option=='stellarator', then <nfield_periods> should be chosen very 
carefully, since <twist_and_shift_geo_factor> will depend on the local shear at 
the end of the flux tube. The script <calculate_twistAndShiftGeoFactor()> can be 
used to find <nfield_periods> where <twist_and_shift_geo_factor>=1 so that dkx=dky 
and kymax=kymax for nx=ny. 

One can also choose <jtwist> manually by setting <jtwist_overwrite>. 

Hanne Thienpondt
30/05/2022

"""

#!/usr/bin/python3 
import sys, os
import pathlib
import numpy as np

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)
from stellapy.calculations.calculate_twistAndShiftGeoFactor import calculate_twistAndShiftGeoFactor
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_gridDivisionsAndSize 
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_ny, calculate_nx 
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_iota  
from stellapy.data.geometry.recognize_device import recognize_device
from stellapy.data.geometry.read_wout import read_woutFile
        
#===============================================================================
#                      Calculate (nx, ny, nfield_periods)                      #
#===============================================================================
         
def calculate_NxNyNfieldperiods(
        # Option 1: Give a vmec path
        vmec=None,
        # Option 2: Define a miller equilibrium
        miller=False,
        shat=0.796, 
        q=1.4,
        # Radial location
        rho=0.8,
        svalue=None,
        # Define the desired box 
        kxmax=2, 
        kymax=2, 
        y0=15, 
        # Flux tube length
        poloidal_turns=1, 
        nfield_periods=None, 
        nperiod=1, 
        # Boundary conditions
        boundary_option="linked", 
        # One can manually overwrite jtwist in the input file
        jtwist_overwrite=None):  
        
        # Radial location
        if svalue==None: svalue = rho*rho
        if rho==None: rho = np.sqrt(svalue)
        
        # Get the geometric quantities 
        if vmec: nfield_periods, poloidal_turns, wout_variables = read_vmecData(vmec, svalue, poloidal_turns, nfield_periods)
        if miller: nfield_periods, poloidal_turns, wout_variables = read_millerData(shat, q, nperiod)
            
        # Calculate the box size
        args = calculate_gridDivisionsAndSize(y0, nfield_periods, wout_variables, svalue=svalue, nperiod=nperiod, jtwist=jtwist_overwrite)
        twist_and_shift_geo_factor = args["twist_and_shift_geo_fac"]
        dky = args["dky"]; dkx = args["dkx"]; jtwist = args["jtwist"]
        
        # For the stellarator symmetric boundary conditions, dkx depends on nfield_periods
        if boundary_option=="stellarator":
            _, twist_and_shift_geo_factor, _ = calculate_twistAndShiftGeoFactor(vmec, nfield_periods=nfield_periods, rho=rho, verbose=False)
            twist_and_shift_geo_factor = np.abs(twist_and_shift_geo_factor[-1])
            dkx = twist_and_shift_geo_factor*dky/jtwist
            
        # Calculate ny
        nky = round(kymax/dky,0)+1
        ny = calculate_ny(nky) 
        
        # Calculate nx
        nkx = 2*round(kxmax/dkx,0)+1 
        nx = calculate_nx(nkx)  
        
        # Print the results aminor
        jtwist_overwrite = "None" if jtwist_overwrite==None else jtwist_overwrite 
        device = recognize_device(vmec) if recognize_device(vmec) else "OVERVIEW"
        print()
        print('                     ====================== ')
        print('                              '+device+'       ')
        print('                     ======================\n ')
        print('                     ------- INPUT -------\n')
        print("{0:>30}".format("rho"), " ", "{0:<10}".format(rho))
        print("{0:>30}".format("y0"), " ", "{0:<10}".format(y0))
        print("{0:>30}".format("kxmax"), " ", "{0:<10}".format(kxmax))
        print("{0:>30}".format("kymax"), " ", "{0:<10}".format(kymax))
        print("{0:>30}".format("poloidal_turns"), " ", "{0:<10}".format(poloidal_turns))
        print("{0:>30}".format("boundary option"), " ", "{0:<10}".format(boundary_option))
        print("{0:>30}".format("twist_and_shift"), " ", "{0:<10}".format(twist_and_shift_geo_factor))
        print("{0:>30}".format("jtwist_overwrite"), " ", "{0:<10}".format(jtwist_overwrite))
        print("{0:>30}".format("jtwist"), " ", "{0:<10}".format(args["jtwist"]))
        print('\n                     -------- BOX --------\n')
        print("{0:>30}".format("nx"), " ", "{0:<10}".format(nx))
        print("{0:>30}".format("ny"), " ", "{0:<10}".format(ny))
        print("{0:>30}".format("dkx"), " ", "{0:<10}".format(dkx))
        print("{0:>30}".format("dky"), " ", "{0:<10}".format(dky))
        print("{0:>30}".format("dkx/dky"), " ", "{0:<10}".format(dkx/dky))  
        print("{0:>30}".format("kx max"), " ", "{0:<10}".format((nkx-1)/2*dkx))
        print("{0:>30}".format("ky max"), " ", "{0:<10}".format((nky-1)*dky))
        print('\n                     ------ GEOMETRY ------\n')
        print("{0:>30}".format("nfield_periods"), " ", "{0:<10}".format(nfield_periods))
        print("{0:>30}".format("shat"), " ", "{0:<10}".format(args["shat"]))
        print("{0:>30}".format("iota"), " ", "{0:<10}".format(args["iota"]))
        print("{0:>30}".format("1/iota"), " ", "{0:<10}".format(1/args["iota"]))
        print("{0:>30}".format("jtwist"), " ", "{0:<10}".format(args["jtwist"]))
        print("{0:>30}".format("2*pi*shat"), " ", "{0:<10}".format(2*np.pi*args["shat"]))
        print("{0:>30}".format("twist_and_shift_geo_fac"), " ", "{0:<10}".format(twist_and_shift_geo_factor))
        print("{0:>30}".format("twist_and_shift_geo_fac_linked"), " ", "{0:<10}".format(args["twist_and_shift_geo_fac"]))
        if vmec: 
            print('\n                     -------- VMEC --------\n')
            print("{0:>30}".format("aminor"), " ", "{0:<10}".format(wout_variables["aminor"]))
            print("{0:>30}".format("rmajor"), " ", "{0:<10}".format(wout_variables["rmajor"]))
            print("{0:>30}".format("field_period_ratio"), " ", "{0:<10}".format(args["field_period_ratio"]))
            print("{0:>30}".format("nfp"), " ", "{0:<10}".format(round(wout_variables["nfp"],2)))
            print("{0:>30}".format("nfield/nfp"), " ", "{0:<10}".format(round(nfield_periods/wout_variables["nfp"],2)))
            print()
            
#-----------------------------------------
def read_vmecData(vmec_path, svalue, poloidal_turns, nfield_periods):
    wout_variables = read_woutFile(type('Paths', (object,), {"vmec" : vmec_path, "geometry" : "NotFound"})) 
    iota = calculate_iota(svalue=svalue, ns=wout_variables['ns'], iotas=wout_variables['iotas'])    
    if not nfield_periods: nfield_periods = poloidal_turns*wout_variables['nfp']/iota
    if not poloidal_turns: poloidal_turns = nfield_periods*iota/wout_variables['nfp']
    return nfield_periods, poloidal_turns, wout_variables

#-----------------------------------------
def read_millerData(shat, q, nperiod):
    wout_variables = {"source" : "Miller coordinates"} 
    wout_variables["nfp"] = 1
    wout_variables["shat"] = shat
    wout_variables["iota"] = 1/q  
    nfield_periods = (2*nperiod-1) 
    poloidal_turns = nfield_periods*wout_variables["iota"]
    return nfield_periods, poloidal_turns, wout_variables
 
#-----------------------------------------    
if __name__ == "__main__": 
    calculate_NxNyNfieldperiods()  