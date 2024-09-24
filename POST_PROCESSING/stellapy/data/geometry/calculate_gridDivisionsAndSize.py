'''
################################################################################
                          MODULE GET GRID DIVISIONS 
################################################################################
Obtain the divisions in kx and ky for non-linear runs.

For non-linear runs we simulate a box of size (x0, 2*pi*y0) around the magnetic field line.
The input parameters are (y0, nzed, ny and nx); then (x0, nakx, naky) are calculated by stella.

The resolution of the simulation is:
    - nzed along the magnetic field line
    - ny along the y-axis with length 2*pi*y0  
    - nx along the x-axis with length x0 ~ 1/hat(s)  with  hat(s) = -2s/iota  with  s = (r/a)^2
   

Parameters
----------
from stellapy.data.read_wout import read_wout # Circular dependency 
wout_data = read_wout(woutfile)
    
Theory
------
The size of the grid and it's divisions are based on the standard parallel boundary condition.
Their values are calculated according to the paper: "The parallel boundary condition 
for turbulence simulations in low magnetic shear devices." by M. F. Martin [1]
Or according to "Definitions for GS2 full-flux-surface stellarator geometry" by
M. Landreman which can be found in "stella/geo/vmec_interface/doc" [2]


CASE GEO_OPTION_VMEC
---------------------
The input parameters are read from a vmec file.
Moreover, nalpha may be specified via input file.


Default vmec_parameters
-----------------------
vmec_filename = 'equilibria/wout_w7x_standardConfig.nc'
alpha0 = 0.0
zeta_center = 0.0
nfield_periods = -1.0
torflux = 0.6354167d+0
surface_option = 0
verbose = .true.
zgrid_scalefac = 1.0
zgrid_refinement_factor = 1


Geometry coefficients from vmec
-------------------------------
nzgrid, nalpha, geo_surf, grho, bmag, gradpar, grad_alpha_grad_alpha, &
grad_alpha_grad_psi, grad_psi_grad_psi, field_period_ratio, &
gds23, gds24, gds25, gds26, gbdrift_alpha, gbdrift0_psi, &
cvdrift_alpha, cvdrift0_psi, calculate_signTorFlux, &
theta_vmec, zed_scalefac, aref, bref, alpha, zeta, &
 

Signs
-----
-sign(psi_tor) = sign(psi) = dxdpsi_sign*calculate_signTorFlux
dxdpsi_sign = -1
dydalpha_sign = 1
dxdXcoord_sign = -1


Coordinates
-----------
B = nabla psi x nabla alpha
x = dx/dpsi*(psi-psi_0)
y = dy/dalpha*(alpha-alpha_0)
kx = kpsi*dpsi/dx
ky = kalpha*dalpha/dy


Other coordinates
-------------------
x = a*sqrt(s)
s = psi_toroidal/psi_edge
psi = -psi_tor
rho = rhotor = sqrt(s) = r/a = sqrt(psitor/psitor_LCFS) = sqrt(torflux)  (FOR VMEC)

        
Definitions
-----------
q = 1/iota 
Bref = 2*abs(psi_tor_LCFS)/a^2
field_period_ratio = nfield_periods/nfp  with nfp = 5 for Wendelstein 7-X


Derivatives
-----------
dx/dpsi = -dx/dpsi_tor = 1/(a*Bref)*sqrt(psitor_LCFS/psitor) = -sign(psi_tor)/(rhotor*a*Bref)
dy/dalpha = a*sqrt(psitor/psitor_LCFS) = a*sign(dydalpha)*rhotor


Normalized quantities
---------------------
psi_N = -psitor/(a**2*Bref)         (FOR VMEC)


Normalized derivatives
----------------------
dxdpsi_N   = a*Bref*dx/dpsi = -a*Bref*dx/dpsi_tor 
           = -sign(psi_tor)/rhotor = dxdpsi_sign*calculate_signTorFlux/rhotor
dydalpha_N = (dy/dalpha) / a = sign(dydalpha)*(dpsi/dr)/(a*Bref) = dydalpha_sign*dpsidrho 
           = sign(dydalpha)*rhotor = dydalpha_sign*geo_surf%rhotor
drhodpsi_N = -drho/d(rho**2)*(aref**2*Bref/psitor_lcfs) = -1.0/rho = dxdpsi_sign*calculate_signTorFlux/rhotor    

    
Derivations
-----------
shat = x/q*(dq/dx) 
     = (a*sqrt(s))/(1/iota) * (d(1/iota)/d(a*sqrt(s))
     = (a*sqrt(s)*iota) * (-1/iota^2) * (2*sqrt(s)/a) * diota/ds
     = -2 * s/iota * diota/ds
     
Lx   = 2*pi/dkx = 2*pi*jtwist/(factor*dky) = -jtwist/(shat*iota*field_period_ratio*dky)
     = L/(P*shat*dky) since Lx/Ly = jtwist/(2*np.pi*P*shat)
         
calculate_twistAndShiftGeoFac = dkx/dky*jtwist
                        = -2.*pi*shat*drhodpsi*dydalpha/(qinp*dxdpsi*rhotor)*field_period_ratio
                        = -2.*np.pi*shat*iota*drhodpsi_N*dydalpha_N/(dxdpsi_N*np.sqrt(svalue))*field_period_ratio
                        = -2.*np.pi*shat*iota*field_period_ratio
                        
exb_nonlin_fac = the factor in front of equation (23) (M.Barnes)
               = Bref/2*dy/dalpha*dx/dpsi 
               = 0.5 * a*Bref*dx/dpsi * (dy/dalpha)/a 
               = 0.5 * dxdXcoord * dydalpha
                
gds2 = |grad y|^2 = |grad alpha|^2 * (dy/dalpha)^2
     = grad_alpha_grad_alpha * dydalpha**2

gds21 = shat * grad x . grad y = shat * dx/dpsi_t * dy/dalpha * grad alpha . grad psi_t 
      = -grad_alpha_grad_psi * geo_surf%shat * dxdXcoord * dydalpha

gds22 = shat^2 * |grad x|^2 = shat^2 * |grad psi_t|^2 * (dx/dpsi_t)^2
      = geo_surf%shat**2 * grad_psi_grad_psi * dxdXcoord**2

gbdrift_alpha and cvdrift_alpha contain the grad-B and curvature drifts projected onto
the grad alpha direction and need the projections on grad y
   gbdrift = gbdrift_alpha * dydalpha
   cvdrift = cvdrift_alpha * dydalpha

gbdrift0_psi and cvdrift0_psi contain the grad-B and curvature drifts projected onto
the grad psi direction need the projections on grad x
   gbdrift0 = gbdrift0_psi * dxdXcoord
   cvdrift0 = cvdrift0_psi * dxdXcoord

################################################################################
'''
 
import numpy as np 
from stellapy.data.output.read_outputFile import read_netcdfVariables
from stellapy.data.input.read_inputFile import calculate_stellaVariables

#################################################################
#                        ARGUMENTS
#################################################################
def calculate_gridDivisionsAndSize(y0, nfield_periods, wout_data, svalue=None, rho=None, jtwist=None, nperiod=1):
    ''' The functions rely on {svalue, y0, nfield_periods, wout_data}. ''' 
    args  =  {
        # Variables related to the simulation
        "rho"    : rho,\
        "svalue" : svalue,\
        "y0"     : y0,\
        "nperiod" : nperiod,\
        # Variables that can be given or calculated
        "jtwist": jtwist,\
        # Variables that are calculated
        "P"     : np.nan,\
        "dkx"   : np.nan,\
        "shat"  : np.nan,\
        "iota"  : np.nan,\
        "diotaOverds"  : np.nan,\
        "torflux_sign" : np.nan,\
        "twist_and_shift_geo_fac" : np.nan,\
        # Print information
        "verbose" : False,\
        }
    
    # Determine the position
    if rho==None: args["rho"] = np.sqrt(args["svalue"] )
    if svalue==None: args["svalue"] = args["rho"]*args["rho"] 
    
    # Calculate the variables at the flux tube for VMEC
    if wout_data["source"]!="Miller coordinates": 
        
        # Nfieldperiods
        args["nfield_periods"] = nfield_periods
        args["field_period_ratio"] = nfield_periods/wout_data["nfp"]
        
        # Variables read from the equilibrium
        args["b0"] = wout_data["b0"]
        args["ns"] = wout_data["ns"]
        args["nfp"] = wout_data["nfp"]
        args["phi"] = wout_data["phi"]
        args["iotas"] = wout_data["iotas"]
        args["iotaf"] = wout_data["iotaf"]
        args["Aminor_p"] = wout_data["aminor"]
             
        # Calculations at rho=rho0
        args["torflux_sign"] = calculate_signTorFlux(**args) 
        args["iota"]         = calculate_iota(**args) 
        args["diotaOverds"]  = calculate_diotaOverds(**args) 
        args["shat"]         = calculate_shat(**args)  
        
    # For Miller, the variables are already given
    if wout_data["source"]=="Miller coordinates":
        args["torflux_sign"] = "Miller"
        args["iota"] = wout_data["iota"]
        args["shat"] = wout_data["shat"]
        
    # Only (iota, shat) are needed for the following calculations
    args["dky"]          = calculate_dky(args["y0"])
    args["Ly"]           = calculate_Ly(args["y0"]) 
    args["twist_and_shift_geo_fac"] = calculate_twistAndShiftGeoFac(**args) 
    args["jtwist"]       = calculate_jtwist(**args)  
    args["P"]            = calculate_zetaOverTwoPi(**args)
    args["dkx"]          = calculate_dkx(**args)
    args["Lx"]           = calculate_Lx(**args) 
    return args

#################################################################
#                       GRID DIVISIONS
#################################################################

def calculate_naky(ny):
    ''' Returns the number of modes <ky> along the y-direction, when <ny> is given. 
    Note that only positive ky values are simulated due to the mirror condition. 
    The real amount of ky modes is (2*naky-1).'''
    naky = int((ny-1)/3) + 1
    return int(naky)

#---------------- 
def calculate_nakx(nx):
    ''' Returns the number of modes <kx> along the x-direction, when <nx> is given. '''
    nakx = int(2*(int((nx-1)/3))) + 1
    return int(nakx)

#---------------- 
def calculate_ny(naky): 
    ''' Returns the number of grid divisions along y in real space. '''
    ny = 3*(naky-1)+1   
    return int(ny)

#---------------- 
def calculate_nx(nakx): 
    ''' Returns the number of grid divisions along x in real space. '''
    nx = (nakx-1)*3/2+1  
    return int(nx)

#---------------- 
def calculate_dky(y0):
    ''' Returns the step size between the <ky> modes. '''
    return 1.0/y0

#---------------- 
def calculate_dkx(shat=None, dky=None, jtwist=None, twist_and_shift_geo_fac=None, nperiod=None, **_):
    ''' Returns the step size between the <kx> modes. '''

    # Non-quantized boundary condition assumed to be periodic instead of 
    # linked boundary conditions if zero magnetic shear
    if abs(shat) <= shat_zero(): dkx = dky/jtwist 
        
    # The division along x is based on the twist-and-shift boundary condition
    else: dkx = (2*nperiod-1)*dky*twist_and_shift_geo_fac/jtwist  
    return abs(dkx)

#---------------- 
def calculate_aky(ny, y0):
    ''' Returns a list of the <ky> values simulated by stella, when <ny> and <y0> are given. '''

    # Get the step size between the <ky> values and the number of <ky> values
    dky  = calculate_dky(y0)
    naky = calculate_naky(ny)

    # Create a list with <naky> values of <ky> 
    aky = [ round(iky*dky,2) for iky in np.arange(0, naky) ] 
    return aky

#---------------- 
def calculate_akx(nx, dkx, **_):
    ''' Returns a list of the <kx> values simulated by stella, when <nx>, <y0> and <s> are given. '''
    
    # Get the step size between the <kx> values and the number of <kx> values
    nakx = calculate_nakx(nx)
    
    # get the ikx index corresponding to kx_max
    ikx_max = int(calculate_nakx(nx)/2+1)
    
    # kx goes from zero to kx_max down to zero and then from -kx_max to -|kx_min|    
    akx = [ round(ikx*dkx,2) for ikx in np.arange(0, ikx_max) ]       
    akx = akx + [ round((ikx-nakx)*dkx,2) for ikx in np.arange(ikx_max, nakx) ]
    if True: return akx

#################################################################
#                      GRID SIZE
#################################################################

def calculate_zetaOverTwoPi(iota=None, shat=None, twist_and_shift_geo_fac=None, field_period_ratio=None, verbose=False, **_):
    ''' The end of the zeta domain is defined as zeta = 2*pi*P. '''    

    # Calculate P = zeta/2*pi        
    P_method1 = twist_and_shift_geo_fac/(2*np.pi*shat)

    # Compare with a different calculation method:
    if verbose:   
        P_method2 = -iota*field_period_ratio
        print("Calculation 1: zeta/2*pi = ", P_method1)
        print("Calculation 2: zeta/2*pi = ", P_method2)
        
    # Return the factor P
    return P_method1

#------------------------------------
def calculate_Ly(y0):
    ''' Returns the length of the box along y. '''
    return 2*np.pi*y0

#------------------------------------
def calculate_Lx(dkx=None, P=None, dky=None, jtwist=None, shat=None, verbose=False, **_):
    ''' Returns the length of the box along x.
    
    Definitions
    -----------
    Lx = 2*pi/dkx = 2*pi*jtwist/(factor*dky) = -jtwist/(shat*iota*field_period_ratio*dky)
    Lx = L/(P*shat*dky) since Lx/Ly = jtwist/(2*np.pi*P*shat)
    '''

    # Calculation through the division along kx 
    Lx_method1 = 2*np.pi/dkx
    
    # Calculation through Lx/Ly = jtwist/(2*np.pi*P*shat) so Lx = L/(P*shat*dky)
    if verbose:
        Lx_method2 = jtwist/(shat*P*dky)
        print("Calculation 1: Lx = ", Lx_method1)
        print("Calculation 2: Lx = ", Lx_method2)
        
    # Return the box size along x
    return Lx_method1

#------------------------------------
def calculate_Lx_over_Ly(P=None, shat=None, jtwist=None, Lx=None, Ly=None, verbose=False, **_):
    ''' Returns the perpendicular aspect ratio of the domain (see eq.11 in [1]).  '''
    
    # Official equation  
    LxLy_method1 = jtwist/(2*np.pi*P*shat)
    
    # Direct calculation
    if verbose:
        LxLy_method2 = Lx/Ly
        print("Calculation 1: LxLy = ", LxLy_method1)
        print("Calculation 2: LxLy = ", LxLy_method2)
        
    # Return LxLy
    if True: return LxLy_method1
    
#################################################################
#                   MATHEMATIC FACTORS
#################################################################

def calculate_iota(svalue=None, ns=None, iotas=None, **_):
    ''' Calculate the rotational transform iota at s=svalue if <iotas> and <ns> are given, 
    or read it through read_wout(woutfile, s=svalue).
    
    Definition
    ----------
    iota = 1/q = psi_poloidal/psi_toroidal = R*B_theta/(r*B_phi)
    '''    
    
    # Get the weighing factors
    vmec_radial_i_half, vmec_radial_w_half = calculate_vmecRadialWeightHalf(ns, svalue)
    
    # Calculate iota
    iota = iotas[vmec_radial_i_half[0]] * vmec_radial_w_half[0] + iotas[vmec_radial_i_half[1]] * vmec_radial_w_half[1]    
    return iota

#------------------------------------         
def calculate_diotaOverds(svalue=None, ns=None, iotaf=None, **_):
    ''' Return diota/ds at s=svalue. '''
    
    # Get the weighing factors
    vmec_radial_i_half, vmec_radial_w_half = calculate_vmecRadialWeightHalf(ns, svalue)
    
    # Step size of the divisions along s=(r/a)^2=rho^2 during the VMEC run
    ds = 1.0/float(ns-1)
    
    # Calculate diota/ds
    diotaOverds_on_half_grid = np.empty(ns)*0.0 
    diotaOverds_on_half_grid[1:ns-1] = (iotaf[1:ns-1] - iotaf[0:ns-2]) / ds
    diotaOverds = diotaOverds_on_half_grid[vmec_radial_i_half[0]] * vmec_radial_w_half[0]+\
                  diotaOverds_on_half_grid[vmec_radial_i_half[0]] * vmec_radial_w_half[1]
    return diotaOverds

#------------------------------------
def shat_zero():
    ''' Shat_zero is minimum shat value below which the periodic boundary condition is enforced. '''
    return 1.E-2

#------------------------------------
def calculate_shat(svalue=None, iota=None, diotaOverds=None, **_):
    ''' Returns hat(s) = x/q*(dq/dx) =  -2*s/iota*diota/ds. 
    
    Definitions
    -----------
    q = 1/iota 
    x = a*sqrt(s)
    s = psi_toroidal/psi_edge
    
    Derivation
    ----------
    shat = x/q*(dq/dx) 
         = (a*sqrt(s))/(1/iota) * (d(1/iota)/d(a*sqrt(s))
         = (a*sqrt(s)*iota) * (-1/iota^2) * (2*sqrt(s)/a) * diota/ds
         = -2 * s/iota * diota/ds
    
    '''    
    # Calculate shat
    shat = (-2 * svalue / iota) * diotaOverds
    return shat


#------------------------------------
def calculate_jtwist(twist_and_shift_geo_fac=None, jtwist=None, **_):
    ''' J is a non-zero integer that can be set in the code to 
    potentially achieve a more desirable aspect ratio. ''' 
    if np.isnan(twist_and_shift_geo_fac): return np.nan
    if jtwist==None or jtwist<0: jtwist = max(1,int(abs(twist_and_shift_geo_fac)+0.5))
    return jtwist

#------------------------------------
def calculate_twistAndShiftGeoFac(svalue=None, iota=None, shat=None, torflux_sign=None, field_period_ratio=None, verbose=False, **_):
    ''' A factor used directly in the stella code: 
    
    Signs
    -----
    -sign(psi_tor) = sign(psi) = dxdpsi_sign*calculate_signTorFlux
    dxdpsi_sign = -1
    dydalpha_sign = 1
    
    Definitions
    -----------
    psi = -psi_tor
    rho = rhotor = sqrt(s)  = r/a = sqrt(psitor/psitor_LCFS) = sqrt(torflux)
    Bref = 2*abs(psi_tor_LCFS)/a^2
    a*Bref*dx/dpsi_tor = sign(psi_tor)/rhotor
    dxdpsi_N = a*Bref*dx/dpsi = -a*Bref*dx/dpsi_tor = -sign(psi_tor)/rhotor = dxdpsi_sign*calculate_signTorFlux/rhotor; 
    dydalpha_N = (dy/dalpha) / a = sign(dydalpha)*rhotor
    psi_N = -psitor/(a**2*Bref)
    drhodpsi_N = -drho/d(rho**2)*(aref**2*Bref/psitor_lcfs) = -1.0/rho = dxdpsi_sign*calculate_signTorFlux/rhotor
    field_period_ratio = nfield_periods/nfp  with nfp = 5 for Wendelstein 7-X
    
    Derivation
    -----------
    abs(calculate_twistAndShiftGeoFac) = dkx/dky*jtwist
    calculate_twistAndShiftGeoFac = -2.*pi*shat*drhodpsi*dydalpha/(qinp*dxdpsi*rhotor)*field_period_ratio
    
    '''   
    
    # For miller
    if str(torflux_sign)=="Miller":
        twist_and_shift_geo_fac = 2.0*np.pi*shat
        return twist_and_shift_geo_fac
    
    # Define the signs of the derivatives
    dxdpsi_sign   = -1.
    dydalpha_sign =  1. 
    
    # Get the factors
    dydalpha_N = dydalpha_sign*np.sqrt(svalue)                      
    drhodpsi_N = dxdpsi_sign*torflux_sign/np.sqrt(svalue)  # Opposite sign than in marconi 
    dxdpsi_N   = dxdpsi_sign*torflux_sign/np.sqrt(svalue)  # Opposite sign than in marconi  

    # The factor calculated by marconi is: 
    twist_and_shift_geo_fac_method1 = -2.*np.pi*shat*iota*drhodpsi_N*dydalpha_N/(dxdpsi_N*np.sqrt(svalue))*field_period_ratio
    
    # However this can be simplied too: 
    twist_and_shift_geo_fac_method2 = -2.*np.pi*shat*iota*field_period_ratio
    
    # Compare calculations
    if verbose: 
        print("Calculation 1: calculate_twistAndShiftGeoFac = ", twist_and_shift_geo_fac_method1)
        print("Calculation 2: calculate_twistAndShiftGeoFac = ", twist_and_shift_geo_fac_method2)
    if True: return twist_and_shift_geo_fac_method2

#################################################################
#                         SIGNS
#################################################################

def calculate_signTorFlux(phi=None, **_):
    ''' Returns the sign of the toroidal flux. '''
    return np.sign(edge_toroidal_flux_over_2pi(phi))   

#------------------------------------
def edge_toroidal_flux_over_2pi(phi=None, **_):
    ''' Returns psi_edge/(2*pi). '''
    toroidalFlux = phi
    dimension = np.size(toroidalFlux)
    if True: return toroidalFlux[dimension-1] 

#################################################################
#                   OTHERS
#################################################################

def calculate_vmecRadialWeightHalf(ns, torflux):
    ''' Interpolation for an <s> value between two flux surfaces, 
    since this value is not in the output file. '''
    
    # Initiate the weights
    radial_index_full  = np.empty(2, dtype='int')
    radial_weight_full = np.empty(2, dtype='float')
    radial_index_half  = np.empty(2, dtype='int')
    radial_weight_half = np.empty(2, dtype='float')
    
    # Calculate the normalized toroidal flux on the full and half grids
    psiN_full = [ (interpolation-1)/[ns-1] for interpolation in np.arange(1,ns+1,1) ] 
    psiN_half = [ (psiN_full[i] + psiN_full[i+1])*0.5 for i in range(len(psiN_full)-1)]
    normalized_toroidal_flux_full_grid = psiN_full #@UnusedVariable
    normalized_toroidal_flux_half_grid = psiN_half
        
    # Handle quantities for the full grid
    if (torflux>1):
        print("Error# torflux cannot be >1"); return
    elif (torflux<0):
        print("Error# torflux cannot be <0"); return 
    elif (torflux==1):
        radial_index_full[0] = ns-1
        radial_index_full[1] = ns
        radial_weight_full[0] = 0
          
    # torflux is >= 0 and <1, this is the most common case.
    else:
        radial_index_full[0] = np.floor(torflux*(ns-1))+1
        radial_index_full[1] = radial_index_full[0] + 1
        radial_weight_full[0] = radial_index_full[0] - torflux*(ns-1)
  
    # Derive the second component from the first one
    radial_weight_full[1] = 1 - radial_weight_full[0]
  
    # Handle quantities for the half grid
    if (torflux < normalized_toroidal_flux_half_grid[0]):
          
        # We start at element 2 since element 1 is always 0 for quantities on the half grid.
        print("Warning: extrapolating beyond the end of VMEC's half grid.")
        print("Extrapolating towards the magnetic axis.) Results are likely to be inaccurate.")
        radial_index_half[0] = 2
        radial_index_half[1] = 3
        radial_weight_half[0] = (normalized_toroidal_flux_half_grid[1] - torflux) / (normalized_toroidal_flux_half_grid[1] - normalized_toroidal_flux_half_grid[0])
  
    # We are past the VMECs grid
    elif (torflux > normalized_toroidal_flux_half_grid[ns-2]):
        print("Warning: extrapolating beyond the end of VMEC's half grid.")
        print("(Extrapolating towards the last closed flux surface.) Results may be inaccurate.")
        radial_index_half[0] = ns-1
        radial_index_half[1] = ns
        radial_weight_half[0] = (normalized_toroidal_flux_half_grid[ns-1] - torflux) \
             / (normalized_toroidal_flux_half_grid[ns-1] - normalized_toroidal_flux_half_grid(ns-2))
  
    # We are exactly at the last point of the half grid
    elif (torflux == normalized_toroidal_flux_half_grid[ns-2]):
        radial_index_half[0] = ns-1
        radial_index_half[1] = ns
        radial_weight_half[0] = 0
          
    # torflux is inside the half grid, this is the most common case.
    else:
        radial_index_half[0] = np.floor(torflux*(ns-1) + 0.5) + 1
        if (radial_index_half[0] < 2): radial_index_half[0] = 2 
        radial_index_half[1] = radial_index_half[0] + 1
        radial_weight_half[0] = radial_index_half[0] - torflux*(ns-1) - 0.5
  
    # Derive the second component from the first one
    radial_weight_half[1] = 1-radial_weight_half[0]
      
    # Correct the indexing from fortan to python
    radial_index_full = [ i-1 for i in radial_index_full]
    radial_index_half = [ i-1 for i in radial_index_half] 
    return radial_index_half, radial_weight_half

#################################################################
#                      VELOCITY WEIGTHS
#################################################################

def calculate_vpaWeights(input_file, nvgrid=None, vpa_max=None):
    ''' integer :: iv, idx, iseg, nvpa_seg 
    nvgrid=36, vpa_max=3 gives 72 values:
    [-3.00 -2.92 -2.83 -2.75 -2.66 ... 2.58 2.66 2.75 2.83 2.92 3.00]
    [-3.00 -2.92 -2.83 -2.75 -2.66 -2.58 -2.49 -2.41 -2.32 -2.24 -2.15 -2.07
 -1.99 -1.90 -1.82 -1.73 -1.65 -1.56 -1.48 -1.39 -1.31 -1.23 -1.14 -1.06
 -0.97 -0.89 -0.80 -0.72 -0.63 -0.55 -0.46 -0.38 -0.30 -0.21 -0.13 -0.04
 0.04 0.13 0.21 0.30 0.38 0.46 0.55 0.63 0.72 0.80 0.89 0.97 1.06 1.14
 1.23 1.31 1.39 1.48 1.56 1.65 1.73 1.82 1.90 1.99 2.07 2.15 2.24 2.32
 2.41 2.49 2.58 2.66 2.75 2.83 2.92 3.00]

    '''

    # Read the input file 
    if input_file:
        from stellapy.data.input.read_inputFile import read_inFile
        inputs = read_inFile(input_file) 
        inputs = calculate_stellaVariables(inputs)
        nvgrid = inputs["vpamu_grids_parameters"]["nvgrid"]
        vpa_max = inputs["vpamu_grids_parameters"]["vpa_max"]  
    
    # Define the public variable nvpa
    nvpa = 2*nvgrid 
    
    # Initialize the vectors
    vpa = np.empty((nvpa)); vpa[:] = 0
    wgts_vpa = np.empty((nvpa)); wgts_vpa[:] = 0
    
    # Apply equal grid spacing in vpa
    dvpa = 2*vpa_max/(nvpa-1) 
    
    # Obtain the vpa grid for vpa > 0 and then fill in vpa grid for vpa < 0 by mirroring
    for iv in range(nvgrid+1, nvpa+1): 
        vpa[iv-1] = (iv-nvgrid-0.5)*dvpa
    vpa[:nvgrid] = -vpa[nvgrid:][::-1] 

    # Use the Simpson's 3/8 rule at the lower boundary
    delta=0.375*dvpa
    wgts_vpa[0]= delta
    wgts_vpa[1:3] = 3.*delta
    wgts_vpa[3] = delta 

    # Use the composite simpson for all other points
    nvpa_seg = int((nvpa-4)/2)
    delta = dvpa/3.
    for iseg in range(1, nvpa_seg+1): 
        idx = int(2*(iseg-1))+4  
        wgts_vpa[idx-1] = wgts_vpa[idx-1] + delta 
        wgts_vpa[idx+0] = wgts_vpa[idx+0] + 4.*delta
        wgts_vpa[idx+1] = wgts_vpa[idx+1] + delta 
 
    # For the sake of symmetry, use the Simpson's 3/8 rule at the upper boundary
    delta = 0.375*dvpa
    wgts_vpa[nvpa-4] = wgts_vpa[nvpa-4] + delta
    wgts_vpa[nvpa-3:nvpa-1] = wgts_vpa[nvpa-3:nvpa-1] + 3.*delta 
    wgts_vpa[nvpa-1] = wgts_vpa[nvpa-1] + delta 
 
    # Use the composite simpson for all other points
    nvpa_seg = int((nvpa-4)/2)
    delta = dvpa/3.
    for iseg in range(1, nvpa_seg+1):
        idx = int(2*(iseg-1))+1
        wgts_vpa[idx-1] = wgts_vpa[idx-1] + delta
        wgts_vpa[idx+0] = wgts_vpa[idx+0] + 4.*delta
        wgts_vpa[idx+1] = wgts_vpa[idx+1] + delta 
 
    # Divide by 2 to account for the double counting
    wgts_vpa = 0.5*wgts_vpa
    return wgts_vpa

#------------------------------------
def calculate_muWeights(input_file, nalpha=1, bmag_psi0=[[None]], equally_spaced_mu_grid=None):
    """ Bmag_psi has dimension (nz, nalpha=naky). """
        
    # Read the input file 
    from stellapy.data.input.read_inputFile import read_inFile
    inputs = read_inFile(input_file) 
    inputs = calculate_stellaVariables(inputs)
    equally_spaced_mu_grid = inputs["vpamu_grids_parameters"]["equally_spaced_mu_grid"] if equally_spaced_mu_grid==None else equally_spaced_mu_grid
    nspec = inputs["species_knobs"]["nspec"]
    nmu = inputs["vpamu_grids_parameters"]["nmu"] 
    vperp_max = inputs["vpamu_grids_parameters"]["vperp_max"]
    nzgrid = int(inputs["zgrid_parameters"]["nzgrid"]) 

    # For the vmec option we have bmag_psi0 = bmag(zed,alpha)  
    if bmag_psi0[0,0]==None: 
        from stellapy.data.output.read_outputFile import read_outputFile
        netcdf_file = read_outputFile(input_file.with_suffix(".out.nc")) 
        bmag_psi0 = read_netcdfVariables('bmag', netcdf_file)   
        netcdf_file.close()  
        
    # Initialize the vectors
    mu = np.empty((nmu)); mu[:] = 0
    dmu = np.empty((nmu-1)); dmu[:] = 0
    wgts_mu_tmp = np.empty((nmu)); wgts_mu_tmp[:] = 0
    wgts_mu = np.empty((nalpha,2*nzgrid+1,nmu)); wgts_mu[:] = 0
    maxwell_mu= np.empty((nalpha,2*nzgrid+1,nmu,nspec)) ; maxwell_mu[:,:,:,:] = 0.0
    
    #======================================
    #===== Equally spaced grid in mu ======
    #======================================
    # Get an equally space grid in mu with equal grid divisions dmu
    if equally_spaced_mu_grid:
        
        #  Get the maximum value of mu: mu_max = vperp_max**2/(2*max(bmag))
        #  with vperp_max = 3.0 by default and for VMEC: bmag = bmag_psi0(nalpha,z)
        mu_max = vperp_max**2/(2.*np.max(bmag_psi0))
        
        # The first grid point needs to be at dmu/2 to avoid the special point mu=0
        # We have, dmu/2 + (nmu-1)*dmu = mu_max, therefore, dmu = mu_max/(nmu-1/2)
        # Assign this mu interval to all elements of the vector dmu(nmu-1)
        dmu[:] = mu_max/(nmu-0.5)
        
        # The first grid point of mu is dmu/2
        mu[0] = 0.5*dmu[0]
        # The rest of the grid points is given by the previous point + dmu
        for imu in range(1, nmu): 
            mu[imu] = mu[0]+imu*dmu[0]
                   
        # Do the simplest thing to start
        wgts_mu_tmp[:] = dmu[0]  
        
    #============================================
    #===== Gauss-Laguerre quadrature in mu ======
    #============================================
    # Use Gauss-Laguerre quadrature to divide the space grid in mu
    else: 
 
        # Get the Laguerre grids for mu and its weights
        mu, wgts_mu_tmp = calculate_laguerre_grids(mu, wgts_mu_tmp) 
  
        # Rescale the weights with exp(mu_hat)*T*B/(B_0*m)   Eq. (10)
        #    with B_0 = B_min*mu_hatmax/v_perpmaxN**2 and v_th**2 = 2T/m
        #    and vperp_max = (v_perp/v_ths)**2 = 3.0 is an input parameter
        #       rescale = exp(mu_hat)*T*B/m*(1/B_0)
        #               = exp(mu_hat)*T*B/m*(v_perpmaxN**2/(B_min*mu_hatmax))
        #               = exp(mu_hat)*T/m*(v_perpmaxN**2/(B_minN*mu_hatmax))
        # Equation 10 has units (m/s)**2 so divide by v_th**2 = 2T/m
        #       rescale = exp(mu_hat)*(v_perpmaxN**2/(2*B_minN*mu_hatmax))
        wgts_mu_tmp[:] = wgts_mu_tmp*np.exp(mu[:])/(2.*np.min(bmag_psi0)*mu[nmu-1]/vperp_max**2)
        
        # The Laguerre method returns mu_hat = mu*B_0/T so transform this to mu
        # Recall that the normalized mu still has units of mass so ignore the "m"
        #       mu = mu_hat*T/B_0 = mu_hat*T*(v_perpmaxN**2/(B_min*mu_hatmax))
        #     mu_N = mu*B_r/v_thr**2 = mu*B_r*m/(2T)
        #          = mu_hat*T*v_perpmaxN**2/(B_min*mu_hatmax) * B*m/(2T)
        #          = mu_hat*v_perpmaxN**2/(2*B_minN*mu_hatmax) * m
        mu = mu/(2.*np.min(bmag_psi0)*mu[nmu-1]/vperp_max**2)
        
        # Define the grid spacings as the difference between the mu elements
        # Leave dmu(nmu) uninitialized since it should never be used, this way
        # valgrind or similar will return an error if it is
        dmu[:nmu-2] = mu[1:nmu-1]-mu[:nmu-2] 
  
    # This is the mu part of the v-space Maxwellian: (T_psi0/T)*exp(-2*mu*B)
    #     exp(-2*mu*B) = exp(-2*m_s*v_perp**2*B/(2B)) = exp(-m_s*v_perp**2)  
    temperature = np.empty((nspec)); temperature[:] = 1 # = spec%temp_psi0/spec%temp 
    for alpha in range(nalpha):
        for z in range(2*nzgrid+1):
            for s in range(nspec):   
                maxwell_mu[alpha,z,:,s] = np.exp(-2*mu*bmag_psi0[z,alpha]*temperature[s])
   
    # The factor of 2./sqrt(pi) is necessary to account for the 2pi from the integration
    # over the gyro-angle and the 1/pi^(3/2) normalization of the velocity space Jacobian
    for alpha in range(nalpha):
        for z in range(2*nzgrid+1): 
            wgts_mu[alpha,z,:] = 2./np.sqrt(np.pi)*wgts_mu_tmp*bmag_psi0[z,alpha]
            
    # Return the weights
    if True: return wgts_mu
    
#################################################################
#                   LAGUERRE GRIDS
#################################################################
 
def calculate_laguerre_grids(zero, wgt):

    # Variables
    error = False
    weight_roundoff_correction = False

    # Define the length
    n = len(zero)
    
    # Initialize the vectors
    zz = np.empty((n)); zz[:] = 0 
 
    # Error if n is too high
    if (n > 180):
        print('ERROR: can''t get so many laguerre grid points')
        print('ERROR: size(zero)= ', n)
        error = True 

    # search zero from 0 to xmax using evenly spaced grid points
    if (n==1):
        zero[0] = 1.0
        wgt[0] = 1.0
        return
    else:
        nzero = 0
        pold = laguerre_l(n,0.0)
        delx = 0.001
        for i in range(1, 1000+1):     
            x = delx*i
            pnew = laguerre_l(n,x)
            if (pold*pnew < 0.0):
                nzero = nzero+1
                zz[nzero-1] = find_zero(n, x-delx, x, pold, pnew, zz[nzero-1]) 
            pold=pnew 

        for _ in range(0, 3+1):
            delx = delx * 10.
            for i in range(1, 900+1):
                x = delx*(i+100)
                pnew = laguerre_l(n,x)
                if (pold*pnew < 0.0):
                    nzero = nzero+1
                    zz[nzero-1] = find_zero(n, x-delx, x, pold, pnew, zz[nzero-1]) 
                if (nzero == n): break
                pold=pnew 
            if (nzero == n): break  

    zero = zz
    wgt = zz / (n+1)**2 / laguerre_l(n+1,zz)**2

    # Roundoff correction
    if (weight_roundoff_correction):
        i = np.argmax(wgt)
        wgt[i] = 1.0 - sum(wgt[0:i-1]) - sum(wgt[i+1:]) 

    # Check number of found zeros
    if (nzero < n):
        print('ERROR in laguerre: didn''t find all zeros')

    # Check zeros and weights 
    check_laguerre_zeros(zero)
    check_laguerre_weights(wgt,eps=1.0e-7)

    if (error): print('STOP in calculate_laguerre_grids')
    return zero, wgt

#------------------------------------
def laguerre_l(n, x):

    p1 = 1.0 - x
    p2 = 1.0

    if (n==0):
        laguerre_l = p2
        return laguerre_l
    elif (n==1):
        laguerre_l = p1
        return laguerre_l 

    for k in range(2, n+1):
        p = ((2.0*k-1.0-x)*p1 - (k-1.0)*p2) / k
        p2 = p1
        p1 = p 

    laguerre_l = p
    return laguerre_l

#------------------------------------
def laguerre_lp (n, x):
 
    laguerre_lp = n*(laguerre_l(n,x) - laguerre_l(n-1,x)) / x
    return laguerre_lp
  
#------------------------------------
def find_zero (n, xold, xnew, pold, pnew, zz):
 
    eps = 8.8817841970012523E-016 
    x1 = xold
    x2 = xnew
    p1 = pold
    p2 = pnew
    maxit=100
 
    # bisection
    for i in range(1, maxit+1):
        zz = (x1+x2) * 0.5
        pz = laguerre_l(n,zz)
        if (abs(pz) <= eps): return zz
        if (pz*p1 < 1E-17): p2=pz ; x2=zz
        else: p1=pz ; x1=zz  
 
    if (i==maxit+1):
        # newton-raphson
        if (zz==x1): x1 = x2
        for i in range(1, maxit+1):
            x1 = zz
            p1 = laguerre_l(n,x1)
            zz = x1 - p1 / laguerre_lp(n,x1)
            pz = laguerre_l(n,zz)
            if (min(abs(zz/x1-1.0), abs(pz)) < eps): break 
         
        if (i==maxit+1):
            print('WARNING: too many iterations in calculate_laguerre_grids')
            print('Why do we never use p2?', p2)
            return  
    
    return zz
         
#------------------------------------   
def check_laguerre_zeros (zero):
 
    # Variables
    error = False  
    n = len(zero)
 
    # Check positiveness 
    if any(x<0 for x in zero):
        print('ERROR in laguerre: grid not positive')
        error = True 
 
    # check alignment
    if any(zero[1:n]-zero[0:n-1] < 0.0):
        print('ERROR in laguerre: wrong alignment')
        error = True 
 
    # Check distances are increasing
    for i in range(0, n-2):
        if (zero[i+1]-zero[i] > zero[i+2]-zero[i+1]):
            print('ERROR in laguerre: distances are decreasing at i= ', i)
            error = True 
 
    if (error): print('STOP in check_laguerre_zeros')
    return

#------------------------------------
def check_laguerre_weights (wgt, eps):

    # Variables
    error = False  
    n = len(wgt)
 
    # check if weights are all positive
    if any(wgt <= 0.0):
        print('ERROR in laguerre: weights are not positive at n =', n)
        error = True 
 
    # check if there is a single maximum
    imax = np.argmax(wgt)  
    if (any(wgt[0:imax] > wgt[1:imax+1])):
        print('ERROR in laguerre: weights decreasing before maximum')
        error = True 
        
    if (any(wgt[imax:n-1] < wgt[imax+1:n])):
        print('ERROR in laguerre: weights increasing after maximum')
        error = True 
 
    if (error): print('STOP error in check_laguerre_weights')
    if True: return 
        
#==================================================
#  ATTACH THE NETCDF DATA TO THE SIMULATION OBJECT
#==================================================
 
def calculate_velocityWeights(self):
    """ Add the potential squared averaged over (kx,ky) as an attribute of the simulation object. """
    
    # Initiate the attributes 
    self.wgts_vpa = np.empty((self.nvgrid)); self.wgts_vpa[:] = np.NaN
    self.wgts_mu = np.empty((2*self.nzgrid+1, self.nmu)); self.wgts_mu[:,:] = np.NaN 
    
    # Allocate the attributes
    self.wgts_vpa = calculate_vpaWeights(self.input_files[0])   
    self.wgts_mu = calculate_muWeights(self.input_files[0])[0,:,:] # Only get the first tube
    return
    
  
    
    