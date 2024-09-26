
import numpy as np
import netCDF4 as nc4   
from stellapy.plot.geometry.utils.get_interpolationWeights import get_interpolationWeightsOnHalfGrid
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_iota
    
def get_CosineAndSineTransformOfRZQ(wout_path, svalue, quantity, ntheta=50, nzeta=50):

    # Read the netcdf file
    dataset = nc4.Dataset(wout_path, mode='r')
  
    # Poloidal and toroidal mode numbers
    xm = dataset.variables['xm'][:]             #@UnusedVariable # Poloidal mode numbers  xm(mn_mode)   
    xn = dataset.variables['xn'][:]             #@UnusedVariable # Toroidal mode numbers  xn(mn_mode)  
    xm_nyq = dataset.variables['xm_nyq'][:]     # Poloidal mode numbers (Nyquist) xm_nyq(mn_mode_nyq)
    xn_nyq = dataset.variables['xn_nyq'][:]     # Toroidal mode numbers (Nyquist) xn_nyq(mn_mode_nyq)
    
    # Dimensions
    ns  = int(dataset.variables['ns'][:])               
    nfp = int(dataset.variables['nfp'][:]) 
    mnmax = int(dataset.variables['mnmax'][:]) 
     
    # Read rmnc, zmns and bmnc versus (radius,mn_mode)
    rmnc = dataset.variables['rmnc'][:]     # cosmn component of cylindrical R, full mesh  
    zmns = dataset.variables['zmns'][:]     # sinmn component of cylindrical Z, full mesh 
    bmnc = dataset.variables['bmnc'][:]     # cosmn component of mod-B, half mesh   
    lmns = dataset.variables['lmns'][:]     # cosmn component of lambda, half mesh   
    gmnc = dataset.variables['gmnc'][:]     # cosmn component of jacobian, half mesh   
    
    # Select the chosen quantity
    if quantity=="B": qmnc = bmnc;      qmns = [None]
    if quantity=="J": qmnc = gmnc;      qmns = [None]
    if quantity=="L": qmnc = [None];    qmns = lmns

    # Axis in field-aligned coordinates, in stellarators usually denoted by (u,v)
    theta = np.linspace(0, 2*np.pi, num=ntheta)       # u 
    zeta  = np.linspace(0, 2*np.pi/nfp, num=nzeta)    # v
    
    # Copy the toroidal and poloidal mode numbers to ease the notation
    m = np.array([ int(m) for m in xm_nyq ])
    n = np.array([ int(n) for n in xn_nyq ])
       
    # Get the interpolation weights, for when svalue lies between two flux surfaces   
    surfaces_indexes, surfaces_weights = get_interpolationWeightsOnHalfGrid(ns, svalue)
    
    # Initiate the (R,Z) coordinates and the quantity Q
    R_mn = np.zeros((ntheta, nzeta, mnmax))
    Z_mn = np.zeros((ntheta, nzeta, mnmax))
    Q_mn = np.zeros((ntheta, nzeta, mnmax))
    
    # Consine and Sine Fourier transforms 
    for interpolation in [0,1]:
        for itheta in range(len(theta)):
            for izeta in range(len(zeta)):
                
                # When svalue denotes a flux surface between two flux surfaces
                # that are given by vmec, we need to interpolate the values
                weight = surfaces_weights[interpolation]
                isurf  = surfaces_indexes[interpolation]
                
                # Get the angle 
                angle = m * theta[itheta] - n * zeta[izeta]
                
                # Get the cosinus and sinus
                cos_angle = np.cos(angle)
                sin_angle = np.sin(angle)
                
                # Calculate the (R,Z) coordinates  
                R_mn[itheta, izeta, :] += rmnc[isurf,:]*cos_angle*weight
                Z_mn[itheta, izeta, :] += zmns[isurf,:]*sin_angle*weight 
                
                # Calculate the quantity Q
                if np.all(qmnc[0]!=None): Q_mn[itheta, izeta, :] += qmnc[isurf,:]*cos_angle*weight
                if np.all(qmns[0]!=None): Q_mn[itheta, izeta, :] += qmns[isurf,:]*sin_angle*weight
                
    # Sum away the (m,n) modes
    R = np.sum(R_mn[:, :, :], axis=2)
    Z = np.sum(Z_mn[:, :, :], axis=2)
    Q = np.sum(Q_mn[:, :, :], axis=2) 
    
    # Return R(zeta,theta); Z(zeta, theta) and Q(zeta, theta)
    return R, Z, Q, zeta, theta, nfp

#-------------------------------------
def get_CosineAndSineTransformOfRZQ_forFieldLine(wout_path, svalue, quantity, nzeta=50, toroidal_turns=1, alpha0=0):

    # Read the netcdf file
    dataset = nc4.Dataset(wout_path, mode='r')
  
    # Poloidal and toroidal mode numbers
    xm = dataset.variables['xm'][:]             #@UnusedVariable # Poloidal mode numbers  xm(mn_mode)   
    xn = dataset.variables['xn'][:]             #@UnusedVariable # Toroidal mode numbers  xn(mn_mode)  
    xm_nyq = dataset.variables['xm_nyq'][:]     # Poloidal mode numbers (Nyquist) xm_nyq(mn_mode_nyq)
    xn_nyq = dataset.variables['xn_nyq'][:]     # Toroidal mode numbers (Nyquist) xn_nyq(mn_mode_nyq)
    
    # Dimensions
    ns  = int(dataset.variables['ns'][:])               
    nfp = int(dataset.variables['nfp'][:]) 
    mnmax = int(dataset.variables['mnmax'][:]) 
     
    # Read rmnc, zmns and bmnc versus (radius,mn_mode)
    rmnc = dataset.variables['rmnc'][:]     # cosmn component of cylindrical R, full mesh  
    zmns = dataset.variables['zmns'][:]     # sinmn component of cylindrical Z, full mesh 
    bmnc = dataset.variables['bmnc'][:]     # cosmn component of mod-B, half mesh   
    lmns = dataset.variables['lmns'][:]     # cosmn component of lambda, half mesh   
    gmnc = dataset.variables['gmnc'][:]     # cosmn component of jacobian, half mesh   
    
    # Select the chosen quantity
    if quantity=="B": qmnc = bmnc;      qmns = [None]
    if quantity=="J": qmnc = gmnc;      qmns = [None]
    if quantity=="L": qmnc = [None];    qmns = lmns
    
    # Needed quantities
    iota = calculate_iota(svalue, ns, dataset.variables['iotas'][:]) 
    
    # Get the pest theta and initiate the vmec theta    
    poloidal_turns = toroidal_turns*iota
    theta_pest = np.linspace(-poloidal_turns*np.pi, poloidal_turns*np.pi, nzeta)
    theta_vmec = np.empty((nzeta))
                                 
    # The zeta array is calculted from alpha0 = theta_pest - iota * zeta
    zeta = np.array([(x-alpha0)/iota for x in theta_pest])

    # Copy the toroidal and poloidal mode numbers to ease the notation
    m = np.array([ int(m) for m in xm_nyq ])
    n = np.array([ int(n) for n in xn_nyq ])
       
    # Get the interpolation weights, for when svalue lies between two flux surfaces   
    surfaces_indexes, surfaces_weights = get_interpolationWeightsOnHalfGrid(ns, svalue)
    
    # Initiate the (R,Z) coordinates and the quantity Q
    R_mn = np.zeros((nzeta, mnmax))
    Z_mn = np.zeros((nzeta, mnmax))
    Q_mn = np.zeros((nzeta, mnmax))
        
    # Collect the vmec data needed for <pest_to_vmec>
    vmec_data = {'xm' : xm, 'xn' : xn, 'mnmax' : mnmax, 'lmns' : lmns}
    vmec_data['surfaces_indexes'] = surfaces_indexes; vmec_data['surfaces_weights'] = surfaces_weights
    
    # Consine and Sine Fourier transforms 
    for interpolation in [0,1]:
        for izeta in range(len(zeta)):
            
            # Calculate the theta angle for vmec
            theta_vmec[izeta] = pest_to_vmec(theta_pest[izeta], zeta[izeta], vmec_data)
            
            # When svalue denotes a flux surface between two flux surfaces
            # that are given by vmec, we need to interpolate the values
            weight = surfaces_weights[interpolation]
            isurf  = surfaces_indexes[interpolation]
            
            # Get the angle 
            angle = m * theta_vmec[izeta] - n * zeta[izeta]
            
            # Get the cosinus and sinus
            cos_angle = np.cos(angle)
            sin_angle = np.sin(angle)
            
            # Calculate the (R,Z) coordinates  
            R_mn[izeta, :] += rmnc[isurf,:]*cos_angle*weight
            Z_mn[izeta, :] += zmns[isurf,:]*sin_angle*weight 
            
            # Calculate the quantity Q
            if np.all(qmnc[0]!=None): Q_mn[izeta, :] += qmnc[isurf,:]*cos_angle*weight
            if np.all(qmns[0]!=None): Q_mn[izeta, :] += qmns[isurf,:]*sin_angle*weight
                
    # Sum away the (m,n) modes
    R = np.sum(R_mn[:, :], axis=1)
    Z = np.sum(Z_mn[:, :], axis=1)
    Q = np.sum(Q_mn[:, :], axis=1) 
    
    # Return R(zeta,theta); Z(zeta, theta) and Q(zeta, theta)
    return R, Z, Q, zeta, theta_pest, nfp


#--------------------------
def pest_to_vmec(theta_pest=0.8, zeta0=0.0, vmec_data=None):
    from scipy.optimize import fsolve 
    x_roots = fsolve(lambda x: fzero_residual(x, theta_pest, zeta0, vmec_data), x0=theta_pest)
    return x_roots
 
#--------------------------  
def fzero_residual(x, theta_pest_target=0.0, zeta0=0.0, vmec_data=None):
    ''' Note that lmns and lmnc use the non-Nyquist xm, xn, and mnmax.
    Also note that lmns and lmnc are on the HALF grid.
    residual = (theta_pest based on theta_vmec_try) - theta_pest_target = theta_vmec_try + Lambda - theta_pest_target
    The variable x plays the role of theta_vmec_try. '''
   
    # The data needed from the vmec file
    xm = vmec_data['xm']  
    xn = vmec_data['xn']      
    mnmax = vmec_data['mnmax']   
    lmns = vmec_data['lmns']     
    surfaces_indexes = vmec_data['surfaces_indexes']
    surfaces_weights = vmec_data['surfaces_weights']
    
    # x = theta_vmec_try in Matt's VMEC_2_GS2 interface
    fzero_residual = x - theta_pest_target
    
    # Iterate over the modes
    for imn in range(mnmax):
        
        # Get the angle 
        sin_angle = np.sin(xm[imn]*x - xn[imn]*zeta0)
        
        # Interpolate the calculation between the flux surfaces
        for interpolation in [0,1]:
            
            # When svalue denotes a flux surface between two flux surfaces
            # that are given by vmec, we need to interpolate the values
            weight = surfaces_weights[interpolation]
            isurf  = surfaces_indexes[interpolation]
            
            # Calculate fzero_residual
            fzero_residual = fzero_residual + lmns[isurf, imn]*sin_angle*weight

    return fzero_residual
 
 