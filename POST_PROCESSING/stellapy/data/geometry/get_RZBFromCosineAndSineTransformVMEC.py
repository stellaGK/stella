
#!/usr/bin/python3   
import sys, os
import pathlib
import numpy as np
import pickle
import netCDF4 as nc4   
from stellapy.utils.decorators.exit_program import exit_program

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.data.geometry.get_interpolationWeights import get_interpolationWeightsOnHalfGrid
from stellapy.data.geometry.get_interpolationWeights import get_interpolationWeightsOnFullGrid
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_iota 
from stellapy.utils.files.ensure_dir import ensure_dir

#===============================================================================
#                             Get (R,Z,B) from VMEC                            #
#===============================================================================
    
def get_RZBFromCosineAndSineTransformVMEC(wout_path, svalue, quantity, ntheta=50, nzeta=50):
    
    # Read from pickle 
    pickle_path = pathlib.Path(wout_path).parent/'Pickles'; ensure_dir(pickle_path)
    pickle_path = pickle_path / (pathlib.Path(wout_path).stem+f"_{svalue:.2f}_{quantity}_{nzeta}_{ntheta}.pickle")
    if (os.path.isfile(pickle_path)):  
        file = open(pickle_path, 'rb');  
        DATA = pickle.load(file); file.close()  
        return DATA.R, DATA.Z, DATA.Q, DATA.zeta, DATA.theta, DATA.nfp
    
    # Reac the VMEC file
    xm, xn, xm_nyq, xn_nyq, _, _, ns, nfp, rmnc, zmns, lmns, bmnc, gmnc, _ = read_vmecdata(wout_path)

    # Create a (zeta,theta) grid in cylindrical coordinates with (<ntheta>,<nzeta>) points   
    # In stellarators these coordinates are usually denoted by (u,v)
    theta = np.linspace(0, 2*np.pi, num=ntheta)       # u 
    zeta  = np.linspace(0, 2*np.pi/nfp, num=nzeta)    # v
    
    # Copy the toroidal and poloidal mode numbers to ease the notation
    # Note that in stella we put int(n)/nfp below, however, this is already included above in <zeta> 
    m_nyq = np.array([ int(m) for m in xm_nyq ])
    n_nyq = np.array([ int(n) for n in xn_nyq ])
    m = np.array([ int(m) for m in xm ])
    n = np.array([ int(n) for n in xn ])
       
    # Get the interpolation weights, for when svalue lies between two flux surfaces   
    surfaces_indexes, surfaces_weights = get_interpolationWeightsOnHalfGrid(ns, svalue)

    # Initiate the (R,Z) coordinates and the quantity Q
    R = np.zeros((ntheta, nzeta))
    Z = np.zeros((ntheta, nzeta))
    Q = np.zeros((ntheta, nzeta))
    
    # Select the chosen quantity
    if quantity=="B": qmnc = bmnc;      qmns = [None]
    if quantity=="J": qmnc = gmnc;      qmns = [None]
    if quantity=="L": qmnc = [None];    qmns = lmns
    
    # Consine and Sine Fourier transforms 
    for interpolation in [0,1]:
        for itheta in range(len(theta)):
            for izeta in range(len(zeta)):
    
                # When svalue denotes a flux surface between two flux surfaces
                # that are given by vmec, we need to interpolate the values
                weight = surfaces_weights[interpolation]
                isurf  = surfaces_indexes[interpolation]
    
                # Get the <imn> angles versus (m,n)
                angle = m * theta[itheta] - n * zeta[izeta]
                angle_nyq = m_nyq * theta[itheta] - n_nyq * zeta[izeta]
    
                # Get the cosinus and sinus of these angles
                cos_angle = np.cos(angle); cos_angle_nyq = np.cos(angle_nyq)
                sin_angle = np.sin(angle); sin_angle_nyq = np.sin(angle_nyq)
    
                # Calculate the (R,Z) coordinates by adding each (m,n) contribution 
                R[itheta, izeta] += np.sum(rmnc[isurf,:]*cos_angle[:]*weight)
                Z[itheta, izeta] += np.sum(zmns[isurf,:]*sin_angle[:]*weight)
    
                # Calculate the quantities Q (qmnc), L (lmns) or J (bmnc) from the (m_nyq, n_nyq) contributions
                if np.all(qmnc[0]!=None): Q[itheta, izeta] += np.sum(qmnc[isurf,:]*cos_angle_nyq[:]*weight)
                if np.all(qmns[0]!=None): Q[itheta, izeta] += np.sum(qmns[isurf,:]*sin_angle_nyq[:]*weight)
                
    # Save the data as a pickle
    DATA = DATAOBJECT_RZQ(R, Z, Q, zeta, theta, nfp) 
    file = open(pickle_path,  "wb"); 
    pickle.dump(DATA , file); 
    file.close()      
    return R, Z, Q, zeta, theta, nfp

#--------------------------
class DATAOBJECT_RZQ:
    def __init__(self, R, Z, Q, zeta, theta, nfp):
        self.R = R 
        self.Z = Z 
        self.Q = Q 
        self.zeta = zeta 
        self.theta = theta 
        self.nfp = nfp 
        return
    
#-----------------------------
def get_RZBKappaFromCosineAndSineTransformVMEC(wout_path, svalue, ntheta=50, nzeta=50):
    
    # Read from pickle  
    pickle_path = pathlib.Path(wout_path).parent/'Pickles'; ensure_dir(pickle_path)
    pickle_path = pickle_path / (pathlib.Path(wout_path).stem+f"_{svalue:.2f}_kappa_{nzeta}_{ntheta}.pickle")
    if (os.path.isfile(pickle_path)):  
        file = open(pickle_path, 'rb');  
        DATA = pickle.load(file); file.close()   
        vmec_data = {'xm' : DATA.xm, 'xn' : DATA.xn, 'mnmax' : DATA.mnmax, 'lmns' : DATA.lmns}
        vmec_data['surfaces_indexes'] = DATA.surfaces_indexes; vmec_data['surfaces_weights'] = DATA.surfaces_weights
        return DATA, vmec_data
    
    # Reac the VMEC file
    xm, xn, xm_nyq, xn_nyq, mnmax, mnmax_nyq, ns, nfp, rmnc, zmns, lmns, bmnc, gmnc, iotas = read_vmecdata(wout_path)
    dataset = nc4.Dataset(wout_path, mode='r') 
    bsupumnc = dataset.variables['bsupumnc'][:].data 
    bsupvmnc = dataset.variables['bsupvmnc'][:].data 
    bsubumnc = dataset.variables['bsubumnc'][:].data
    bsubvmnc = dataset.variables['bsubvmnc'][:].data
    bsubsmns = dataset.variables['bsubsmns'][:].data
    Aminor = dataset.variables['Aminor_p'][:].data 
    isigng = dataset.variables['signgs'][:].data
    phi = dataset.variables['phi'][:].data   
    iotaf = dataset.variables['iotaf'][:].data
    iota = calculate_iota(svalue, ns, iotas)  
    edge_toroidal_flux_over_2pi = phi[ns-1] / (2*np.pi)*isigng 
    
    # Create a (zeta,theta) grid in cylindrical coordinates with (<ntheta>,<nzeta>) points   
    # In stellarators these coordinates are usually denoted by (u,v) 
    theta = np.linspace(0, 2*np.pi, num=ntheta)     # u 
    zeta  = np.linspace(0, 2*np.pi/nfp, num=nzeta)  # v
    
    # Copy the toroidal and poloidal mode numbers to ease the notation
    # Note that in stella we put int(n)/nfp below, however, this is already included above in <zeta> 
    m_nyq = np.array([ int(m) for m in xm_nyq ])
    n_nyq = np.array([ int(n) for n in xn_nyq ])
    m = np.array([ int(m) for m in xm ])
    n = np.array([ int(n) for n in xn ])
       
    # Get the interpolation weights, for when svalue lies between two flux surfaces   
    surfaces_indexes, surfaces_weights = get_interpolationWeightsOnHalfGrid(ns, svalue)
    surfaces_indexes_fullgrid, surfaces_weights_fullgrid = get_interpolationWeightsOnFullGrid(ns, svalue)
    
    # Initiate quentities versus r 
    d_R_d_s_mnc = np.zeros((ns, mnmax)) 
    d_Z_d_s_mns = np.zeros((ns, mnmax)) 
    d_Lambda_d_s_mns = np.zeros((ns, mnmax)) 
    d_B_d_s_mnc = np.zeros((ns, mnmax_nyq)) 
    
    # Initiate the coordinates
    Z = np.zeros((ntheta, nzeta)) 
    R = np.zeros((ntheta, nzeta)) 
    d_R_d_s = np.zeros((ntheta, nzeta)) 
    d_R_d_zeta = np.zeros((ntheta, nzeta)) 
    d_R_d_theta = np.zeros((ntheta, nzeta)) 
    d_Z_d_s = np.zeros((ntheta, nzeta)) 
    d_Z_d_zeta = np.zeros((ntheta, nzeta)) 
    d_Z_d_theta = np.zeros((ntheta, nzeta)) 
    d_Lambda_d_s = np.zeros((ntheta, nzeta)) 
    d_Lambda_d_zeta = np.zeros((ntheta, nzeta)) 
    d_Lambda_d_theta = np.zeros((ntheta, nzeta)) 
    
    # Initiate the quantities
    B = np.zeros((ntheta, nzeta)) 
    sqrt_g = np.zeros((ntheta, nzeta)) 
    d_B_d_s = np.zeros((ntheta, nzeta)) 
    d_B_d_zeta = np.zeros((ntheta, nzeta)) 
    d_B_d_theta = np.zeros((ntheta, nzeta)) 
    B_sup_theta = np.zeros((ntheta, nzeta)) 
    B_sup_zeta = np.zeros((ntheta, nzeta)) 
    B_sub_s = np.zeros((ntheta, nzeta)) 
    B_sub_theta = np.zeros((ntheta, nzeta)) 
    B_sub_zeta = np.zeros((ntheta, nzeta)) 
    B_cross_grad_B_dot_grad_alpha = np.zeros((ntheta, nzeta)) 
    B_cross_grad_B_dot_grad_alpha_noshear = np.zeros((ntheta, nzeta)) 
    
    # Step along s=rho**2
    normalized_toroidal_flux_full_grid = [ (i-1)/[ns-1] for i in np.arange(1,ns+1,1) ] 
    ds = normalized_toroidal_flux_full_grid[1] - normalized_toroidal_flux_full_grid[0]
    
    # Consine and Sine Fourier transforms 
    for interpolation in [0,1]:
        for itheta in range(len(theta)):
            for izeta in range(len(zeta)):  
                
                # When svalue denotes a flux surface between two flux surfaces
                # that are given by vmec, we need to interpolate the values
                weight = surfaces_weights[interpolation]
                isurf  = surfaces_indexes[interpolation]
                weight_fullgrid = surfaces_weights_fullgrid[interpolation] 
                isurff_fullgrid  = surfaces_indexes_fullgrid[interpolation]
                        
                # Get the <imn> angles versus (m,n)
                angle = m * theta[itheta] - n * zeta[izeta]
                angle_nyq = m_nyq * theta[itheta] - n_nyq * zeta[izeta]
                
                # Get the cosinus and sinus of these angles
                cos_angle = np.cos(angle); cos_angle_nyq = np.cos(angle_nyq)
                sin_angle = np.sin(angle); sin_angle_nyq = np.sin(angle_nyq)  
            
                # Get the derivative of rmns; zmns and lmns with respect to s 
                # B and Lambda are on the half mesh, so their radial derivatives are on the full mesh.
                # R and Z are on the full mesh, so their radial derivatives are on the half mesh.
                d_R_d_s_mnc[1:ns,:] = (rmnc[1:ns,:] - rmnc[0:ns-1,:]) / ds
                d_Z_d_s_mns[1:ns,:] = (zmns[1:ns,:] - zmns[0:ns-1,:]) / ds
                d_Lambda_d_s_mns[1:ns-1,:] = (lmns[2:ns,:] - lmns[1:ns-1,:]) / ds 
                d_R_d_s_mnc[0] = 0; d_Z_d_s_mns[0] = 0
                d_Lambda_d_s_mns[0] = d_Lambda_d_s_mns[1]
                d_Lambda_d_s_mns[ns-1] = d_Lambda_d_s_mns[ns-2]  
                
                # Get the coordinate grids
                Z[itheta, izeta] += np.sum(zmns[isurff_fullgrid,:]*sin_angle[:])*weight_fullgrid
                R[itheta, izeta] += np.sum(rmnc[isurff_fullgrid,:]*cos_angle[:])*weight_fullgrid
                d_R_d_theta[itheta, izeta] += -np.sum(m[:]*rmnc[isurff_fullgrid,:]*sin_angle[:])*weight_fullgrid 
                d_R_d_zeta[itheta, izeta] += np.sum(n[:]*rmnc[isurff_fullgrid,:]*sin_angle[:])*weight_fullgrid 
                d_Z_d_theta[itheta, izeta] += np.sum(m[:]*zmns[isurff_fullgrid,:]*cos_angle[:])*weight_fullgrid
                d_Z_d_zeta[itheta, izeta] += -np.sum(n[:]*zmns[isurff_fullgrid,:]*cos_angle[:])*weight_fullgrid
                d_Lambda_d_theta[itheta, izeta] += np.sum(m[:]*lmns[isurf,:]*cos_angle[:])*weight
                d_Lambda_d_zeta[itheta, izeta] += -np.sum(n[:]*lmns[isurf,:]*cos_angle[:])*weight       
                  
                # Radial derivatives, since R is on the full mesh, its radial derivative is on the half mesh.
                d_R_d_s[itheta, izeta] += np.sum(d_R_d_s_mnc[isurf]*cos_angle[:])*weight
                d_Z_d_s[itheta, izeta] += np.sum(d_Z_d_s_mns[isurf]*sin_angle[:])*weight
                d_Lambda_d_s[itheta, izeta] += np.sum(d_Lambda_d_s_mns[isurff_fullgrid]*sin_angle[:])*weight_fullgrid
            
                # Get the derivative of bmnc with respect to s and do a simplistic extrapolation at the endpoints (mn_nyq)
                d_B_d_s_mnc[1:ns-1,:] = (bmnc[2:ns,:] - bmnc[1:ns-1,:]) / ds
                d_B_d_s_mnc[0,:] = d_B_d_s_mnc[1,:]; d_B_d_s_mnc[ns-1,:] = d_B_d_s_mnc[ns-2,:]
                
                # Calculate magnetic quantities (derivatives to <s> are clacluated on the full grid
                B[itheta, izeta] += np.sum(bmnc[isurf,:]*cos_angle_nyq[:]*weight)
                sqrt_g[itheta, izeta] += np.sum(gmnc[isurf,:]*cos_angle_nyq[:]*weight)
                
                d_B_d_s[itheta, izeta] += np.sum(d_B_d_s_mnc[isurff_fullgrid,:]*cos_angle_nyq[:]*weight_fullgrid) # d_B_d_s_mnc??
                d_B_d_theta[itheta, izeta] += -np.sum(m_nyq[:]*bmnc[isurf,:]*sin_angle_nyq[:]*weight)
                d_B_d_zeta[itheta, izeta] += np.sum(n_nyq[:]*bmnc[isurf,:]*sin_angle_nyq[:]*weight) # n -> n*nfp???
                
                B_sup_theta[itheta, izeta] += np.sum(bsupumnc[isurf,:]*cos_angle_nyq[:]*weight)
                B_sup_zeta[itheta, izeta] += np.sum(bsupvmnc[isurf,:]*cos_angle_nyq[:]*weight)
                
                B_sub_s[itheta, izeta] += np.sum(bsubsmns[isurff_fullgrid,:]*sin_angle_nyq[:]*weight_fullgrid)
                B_sub_theta[itheta, izeta] += np.sum(bsubumnc[isurf,:]*cos_angle_nyq[:]*weight)
                B_sub_zeta[itheta, izeta] += np.sum(bsubvmnc[isurf,:]*cos_angle_nyq[:]*weight)
                
    # Normalize B
    L_reference = a = Aminor
    B_reference = Bref = 2*abs(edge_toroidal_flux_over_2pi) / (L_reference*L_reference)
    sign_torflux = np.sign(edge_toroidal_flux_over_2pi)
    bmag = B/B_reference 
    
    # Calculate the shear
    d_iota_d_s_on_half_grid = np.zeros(ns)
    d_iota_d_s_on_half_grid[1:ns-1] = (iotaf[1:ns-1] - iotaf[0:ns-2]) / ds
    d_iota_d_s = d_iota_d_s_on_half_grid[surfaces_indexes[0]]*surfaces_weights[0] 
    d_iota_d_s += d_iota_d_s_on_half_grid[surfaces_indexes[1]]*surfaces_weights[1]
    shat = (-2*svalue / iota)*d_iota_d_s
    
    # We need cos(tzeta) and sin(zeta) to transform (R) to (X,Y) as a function of (theta, zeta)
    cos_angle = (np.cos(zeta))[np.newaxis,:]
    sin_angle = (np.sin(zeta))[np.newaxis,:]
    
    # X = R*cos(ζ); dX/dζ = (dR/dζ)*cos(ζ) - R*sin(ζ); dX/dθ = dR/dθ*cos(ζ)
    d_X_d_s = d_R_d_s*cos_angle
    d_X_d_zeta = d_R_d_zeta*cos_angle - R*sin_angle
    d_X_d_theta = d_R_d_theta*cos_angle
    
    # Y = R*sin(ζ); dY/dζ = (dR/dζ)*sin(ζ) + R*cos(ζ); dY/dθ = dR/dθ*sin(ζ)
    d_Y_d_s = d_R_d_s*sin_angle
    d_Y_d_zeta = d_R_d_zeta*sin_angle + R*cos_angle
    d_Y_d_theta = d_R_d_theta*sin_angle
    
    # Use the dual relations to get the Cartesian components of ∇s, ∇θ, and ∇ζ
    #     ∇s = 1/sqrt(g) (e_θ x e_ζ) = 1/sqrt(g) ((dX/dθ)e_X + (dY/dθ)e_Y + (dZ/dθ)e_Z) x ((dX/dζ)e_X + (dY/dζ)e_Y + (dZ/dζ)e_Z)
    #     ∇s . e_X = 1/sqrt(g) ( dY/dθ * dZ/dζ - dZ/dθ * dY/dζ ) 
    #     ∇s . e_Y = 1/sqrt(g) ( dZ/dθ * dX/dζ - dX/dθ * dZ/dζ )
    #     ∇s . e_Z = 1/sqrt(g) ( dX/dθ * dY/dζ - dY/dθ * dX/dζ ) 
    grad_s_X = (d_Y_d_theta * d_Z_d_zeta - d_Z_d_theta * d_Y_d_zeta) / sqrt_g
    grad_s_Y = (d_Z_d_theta * d_X_d_zeta - d_X_d_theta * d_Z_d_zeta) / sqrt_g
    grad_s_Z = (d_X_d_theta * d_Y_d_zeta - d_Y_d_theta * d_X_d_zeta) / sqrt_g
    grad_zeta_X = (d_Y_d_s * d_Z_d_theta - d_Z_d_s * d_Y_d_theta) / sqrt_g
    grad_zeta_Y = (d_Z_d_s * d_X_d_theta - d_X_d_s * d_Z_d_theta) / sqrt_g
    grad_zeta_Z = (d_X_d_s * d_Y_d_theta - d_Y_d_s * d_X_d_theta) / sqrt_g
    grad_theta_X = (d_Y_d_zeta * d_Z_d_s - d_Z_d_zeta * d_Y_d_s) / sqrt_g
    grad_theta_Y = (d_Z_d_zeta * d_X_d_s - d_X_d_zeta * d_Z_d_s) / sqrt_g
    grad_theta_Z = (d_X_d_zeta * d_Y_d_s - d_Y_d_zeta * d_X_d_s) / sqrt_g
    
    # In VMEC s = psi_t/psi_{t,LCFS} so psi_t = s*psi_{t,LCFS}
    # ∇psi_t = (dpsi_t/ds) * ∇s + 0*∇θ + 0*∇ζ = d(psi_{t,LCFS}*s)/ds * ∇s = psi_{t,LCFS} * ∇s 
    grad_psit_X = grad_s_X * edge_toroidal_flux_over_2pi
    grad_psit_Y = grad_s_Y * edge_toroidal_flux_over_2pi
    grad_psit_Z = grad_s_Z * edge_toroidal_flux_over_2pi
    
    # ∇alpha = ∇(theta_pest - iota*zeta) = ∇(theta + Lambda - iota*zeta)
    #        = (dλ/ds - ζ (diota/ds) ∇s + (1 + dλ/dθ) ∇θ + (dλ/dζ - iota) ∇ζ
    # First add the part the final two terms in ∇alpha = (dλ/ds - ζ (diota/ds) ∇s + (1 + dλ/dθ) ∇θ + (dλ/dζ - iota) ∇ζ
    grad_alpha_X = (1 + d_Lambda_d_theta) * grad_theta_X + (-iota + d_Lambda_d_zeta) * grad_zeta_X
    grad_alpha_Y = (1 + d_Lambda_d_theta) * grad_theta_Y + (-iota + d_Lambda_d_zeta) * grad_zeta_Y
    grad_alpha_Z = (1 + d_Lambda_d_theta) * grad_theta_Z + (-iota + d_Lambda_d_zeta) * grad_zeta_Z
    
    # Then add the part which depends on zeta, assuming alpha = theta_pest - iota*zeta
    grad_alpha_X += (d_Lambda_d_s - zeta[np.newaxis,:]*d_iota_d_s)*grad_s_X
    grad_alpha_Y += (d_Lambda_d_s - zeta[np.newaxis,:]*d_iota_d_s)*grad_s_Y
    grad_alpha_Z += (d_Lambda_d_s - zeta[np.newaxis,:]*d_iota_d_s)*grad_s_Z
    
    # Only the secular term, and without the secular term 
    grad_alpha_X_shearingpart = -zeta[np.newaxis,:]*d_iota_d_s*grad_s_X
    grad_alpha_Y_shearingpart = -zeta[np.newaxis,:]*d_iota_d_s*grad_s_Y
    grad_alpha_Z_shearingpart = -zeta[np.newaxis,:]*d_iota_d_s*grad_s_Z
    grad_alpha_X_noshear = grad_alpha_X - grad_alpha_X_shearingpart
    grad_alpha_Y_noshear = grad_alpha_Y - grad_alpha_Y_shearingpart
    grad_alpha_Z_noshear = grad_alpha_Z - grad_alpha_Z_shearingpart
    
    # <grad_alpha_grad_alpha> = (a^2) ∇α . ∇α
    # <grad_alpha_grad_psit> = (1/Bref) ∇α . ∇ψt 
    # <grad_psit_grad_psit> = (1/(a^2*Bref^2)) ∇ψt . ∇ψt 
    grad_alpha_grad_alpha = (a*a) * (grad_alpha_X * grad_alpha_X + grad_alpha_Y * grad_alpha_Y + grad_alpha_Z * grad_alpha_Z) 
    grad_alpha_grad_psit = (1./Bref) * (grad_alpha_X * grad_psit_X + grad_alpha_Y * grad_psit_Y + grad_alpha_Z * grad_psit_Z)   
    grad_psit_grad_psit = 1/(a*a*Bref*Bref) * (grad_psit_X * grad_psit_X + grad_psit_Y * grad_psit_Y + grad_psit_Z * grad_psit_Z) 
    
    # Without shear 
    grad_alpha_grad_alpha_noshear = (a*a) * (grad_alpha_X_noshear * grad_alpha_X_noshear + grad_alpha_Y_noshear * grad_alpha_Y_noshear + grad_alpha_Z_noshear * grad_alpha_Z_noshear) 
    grad_alpha_grad_psit_noshear = (1./Bref) * (grad_alpha_X_noshear * grad_psit_X + grad_alpha_Y_noshear * grad_psit_Y + grad_alpha_Z_noshear * grad_psit_Z)   

    # <B_cross_grad_B_dot_grad_psit> = B x ∇B . ∇ψt = psi_{t,LCFS}/sqrt(g) (- dB/dθ B_zeta + dB/dζ B_theta)  
    B_cross_grad_B_dot_grad_psit = edge_toroidal_flux_over_2pi / sqrt_g * (B_sub_theta * d_B_d_zeta - B_sub_zeta * d_B_d_theta)
    
    # We neglect shear in the calculation of <B_cross_grad_B_dot_grad_alpha>
    # By setting <d_iota_d_s> = 0 or making sure (zeta-zeta_center)=0 at each point
    # So that we have a rectangular cross-section at each point 
    B_cross_grad_B_dot_grad_alpha[:,:] = 0 \
        + (B_sub_s[:,:]*d_B_d_theta[:,:]*(d_Lambda_d_zeta[:,:] - iota) \
        + B_sub_theta[:,:]*d_B_d_zeta[:,:]*(d_Lambda_d_s[:,:] - zeta[:]*d_iota_d_s) \
        + B_sub_zeta[:,:]*d_B_d_s[:,:]*(1 + d_Lambda_d_theta[:,:]) \
        - B_sub_zeta[:,:]*d_B_d_theta[:,:]*(d_Lambda_d_s[:,:] - zeta[:]*d_iota_d_s) \
        - B_sub_theta[:,:]*d_B_d_s[:,:]*(d_Lambda_d_zeta[:,:] - iota) \
        - B_sub_s[:,:]*d_B_d_zeta[:,:]*(1 + d_Lambda_d_theta[:,:])) / sqrt_g[:,:] 
    B_cross_grad_B_dot_grad_alpha_shearingpart = \
        + (B_sub_theta[:,:]*d_B_d_zeta[:,:]*(-d_iota_d_s) \
        - B_sub_zeta[:,:]*d_B_d_theta[:,:]*(-d_iota_d_s)) / sqrt_g[:,:] 
    B_cross_grad_B_dot_grad_alpha_noshear[:,:] = 0 \
        + (B_sub_s[:,:]*d_B_d_theta[:,:]*(d_Lambda_d_zeta[:,:] - iota) \
        + B_sub_theta[:,:]*d_B_d_zeta[:,:]*(d_Lambda_d_s[:,:]) \
        + B_sub_zeta[:,:]*d_B_d_s[:,:]*(1 + d_Lambda_d_theta[:,:]) \
        - B_sub_zeta[:,:]*d_B_d_theta[:,:]*(d_Lambda_d_s[:,:]) \
        - B_sub_theta[:,:]*d_B_d_s[:,:]*(d_Lambda_d_zeta[:,:] - iota) \
        - B_sub_s[:,:]*d_B_d_zeta[:,:]*(1 + d_Lambda_d_theta[:,:])) / sqrt_g[:,:] 
        
    # Calculate grad B (could be wrong!)
    # ∇B = (dB/ds) ∇s + (dB/dtheta) ∇theta + (dB/dzeta) ∇zeta
    gradB_X = d_B_d_s*grad_s_X + d_B_d_theta*grad_theta_X + d_B_d_zeta*grad_zeta_X 
    gradB_Y = d_B_d_s*grad_s_Y + d_B_d_theta*grad_theta_Y + d_B_d_zeta*grad_zeta_Y 
    gradB_Z = d_B_d_s*grad_s_Y + d_B_d_theta*grad_theta_Z + d_B_d_zeta*grad_zeta_Z 
        
    # Collect the vmec data needed for <pest_to_vmec>
    vmec_data = {'xm' : xm, 'xn' : xn, 'mnmax' : mnmax, 'lmns' : lmns}
    vmec_data['surfaces_indexes'] = surfaces_indexes; vmec_data['surfaces_weights'] = surfaces_weights
    
    # Save the data as a pickle
    DATA = DATAOBJECT_KAPPA(zeta, theta, nfp, iota, xm, xn, mnmax, lmns, surfaces_indexes, surfaces_weights, 
        d_R_d_s_mnc, d_Z_d_s_mns, d_Lambda_d_s_mns, d_B_d_s_mnc, Z, R, B, d_R_d_s, d_R_d_zeta, d_R_d_theta, 
        d_Z_d_s, d_Z_d_zeta, d_Z_d_theta, d_Lambda_d_s, d_Lambda_d_zeta, d_Lambda_d_theta, sqrt_g, 
        d_B_d_s, d_B_d_zeta, d_B_d_theta, B_sup_theta, B_sup_zeta, B_sub_s, B_sub_theta, 
        B_sub_zeta, B_cross_grad_B_dot_grad_alpha, ns, normalized_toroidal_flux_full_grid, ds, 
        B_cross_grad_B_dot_grad_alpha_noshear, B_cross_grad_B_dot_grad_alpha_shearingpart, 
        grad_alpha_grad_alpha, grad_alpha_grad_psit, grad_psit_grad_psit, B_cross_grad_B_dot_grad_psit,
        grad_alpha_grad_alpha_noshear, grad_alpha_grad_psit_noshear, shat, sign_torflux, 
        L_reference, B_reference, bmag)
    file = open(pickle_path,  "wb"); 
    pickle.dump(DATA , file); 
    file.close()      
    
    # Return the curvature and the coordinates
    return DATA, vmec_data

#--------------------------
class DATAOBJECT_KAPPA:
    def __init__(self, zeta, theta, nfp, iota, xm, xn, mnmax, lmns, surfaces_indexes, surfaces_weights, 
        d_R_d_s_mnc, d_Z_d_s_mns, d_Lambda_d_s_mns, d_B_d_s_mnc, Z, R, B, d_R_d_s, d_R_d_zeta, d_R_d_theta, 
        d_Z_d_s, d_Z_d_zeta, d_Z_d_theta, d_Lambda_d_s, d_Lambda_d_zeta, d_Lambda_d_theta, sqrt_g, 
        d_B_d_s, d_B_d_zeta, d_B_d_theta, B_sup_theta, B_sup_zeta, B_sub_s, B_sub_theta, 
        B_sub_zeta, B_cross_grad_B_dot_grad_alpha, ns, normalized_toroidal_flux_full_grid, ds, 
        B_cross_grad_B_dot_grad_alpha_noshear, B_cross_grad_B_dot_grad_alpha_shearingpart, 
        grad_alpha_grad_alpha, grad_alpha_grad_psit, grad_psit_grad_psit, B_cross_grad_B_dot_grad_psit,
        grad_alpha_grad_alpha_noshear, grad_alpha_grad_psit_noshear, shat, sign_torflux,
        L_reference, B_reference, bmag):
        self.R = R; self.Z = Z; self.B = B; self.sqrt_g = sqrt_g 
        self.zeta = zeta; self.theta = theta; self.ds = ds; self.ns = ns 
        self.nfp = nfp; self.iota = iota; self.shat = shat; self.sign_torflux = sign_torflux
        self.xm = xm; self.xn = xn; self.mnmax = mnmax 
        self.lmns = lmns; self.surfaces_indexes = surfaces_indexes; self.surfaces_weights = surfaces_weights 
        self.d_R_d_s_mnc = d_R_d_s_mnc; self.d_Z_d_s_mns = d_Z_d_s_mns; self.d_B_d_s_mnc = d_B_d_s_mnc;  self.d_Lambda_d_s_mns = d_Lambda_d_s_mns 
        self.d_R_d_s = d_R_d_s; self.d_Z_d_s = d_Z_d_s; self.d_Lambda_d_s = d_Lambda_d_s 
        self.d_R_d_zeta = d_R_d_zeta; self.d_Z_d_zeta = d_Z_d_zeta; self.d_Lambda_d_zeta = d_Lambda_d_zeta 
        self.d_R_d_theta = d_R_d_theta; self.d_Z_d_theta = d_Z_d_theta; self.d_Lambda_d_theta = d_Lambda_d_theta 
        self.d_B_d_s = d_B_d_s; self.B_sub_s = B_sub_s; self.B_sub_zeta = B_sub_zeta  
        self.d_B_d_zeta = d_B_d_zeta; self.B_sup_zeta = B_sup_zeta 
        self.d_B_d_theta = d_B_d_theta; self.B_sup_theta = B_sup_theta; self.B_sub_theta = B_sub_theta 
        self.B_cross_grad_B_dot_grad_alpha = B_cross_grad_B_dot_grad_alpha 
        self.normalized_toroidal_flux_full_grid = normalized_toroidal_flux_full_grid  
        self.B_cross_grad_B_dot_grad_alpha_noshear = B_cross_grad_B_dot_grad_alpha_noshear
        self.B_cross_grad_B_dot_grad_alpha_shearingpart = B_cross_grad_B_dot_grad_alpha_shearingpart
        self.grad_alpha_grad_alpha = grad_alpha_grad_alpha; self.grad_alpha_grad_psit = grad_alpha_grad_psit
        self.grad_psit_grad_psit = grad_psit_grad_psit; self.B_cross_grad_B_dot_grad_psit = B_cross_grad_B_dot_grad_psit
        self.grad_alpha_grad_alpha_noshear = grad_alpha_grad_alpha_noshear; self.grad_alpha_grad_psit_noshear = grad_alpha_grad_psit_noshear
        self.L_reference = L_reference; self.B_reference = B_reference; self.bmag = bmag
        return

#-------------------------------------
def get_CosineAndSineTransformOfRZQ_forFieldLine(wout_path, svalue, quantity, nzeta=50, toroidal_turns=1, alpha0=0):
    
    # Reac the VMEC file
    xm, xn, xm_nyq, xn_nyq, mnmax, mnmax_nyq, ns, nfp, rmnc, zmns, lmns, bmnc, gmnc, iotas = read_vmecdata(wout_path) 
    
    # Select the quantity
    if quantity=='B': qmnc = bmnc; qmns = [None]
    elif quantity=='lambda': qmnc = gmnc; qmns = [None]
    else: exit_program('No other options are implemented.', get_CosineAndSineTransformOfRZQ_forFieldLine, sys._getframe().f_lineno)
        
    # Needed quantities
    iota = calculate_iota(svalue, ns, iotas) 
    
    # Get the pest theta and initiate the vmec theta    
    poloidal_turns = toroidal_turns*iota
    theta_pest = np.linspace(-poloidal_turns*np.pi, poloidal_turns*np.pi, nzeta)
    theta_vmec = np.empty((nzeta))
                                 
    # The zeta array is calculted from alpha0 = theta_pest - iota * zeta
    zeta = np.array([(x-alpha0)/iota for x in theta_pest]) 
    
    # Copy the toroidal and poloidal mode numbers to ease the notation
    m_nyq = np.array([ int(m) for m in xm_nyq ])
    n_nyq = np.array([ int(n) for n in xn_nyq ])
    m = np.array([ int(m) for m in xm ])
    n = np.array([ int(n) for n in xn ])
       
    # Get the interpolation weights, for when svalue lies between two flux surfaces   
    surfaces_indexes, surfaces_weights = get_interpolationWeightsOnHalfGrid(ns, svalue)
    
    # Initiate the (R,Z) coordinates and the quantity Q
    R_mn = np.zeros((nzeta, mnmax))
    Z_mn = np.zeros((nzeta, mnmax))
    Q_mn = np.zeros((nzeta, mnmax_nyq))
        
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
            angle_nyq = m_nyq * theta_vmec[izeta] - n_nyq * zeta[izeta]
            
            # Get the cosinus and sinus
            cos_angle = np.cos(angle); cos_angle_nyq = np.cos(angle_nyq)
            sin_angle = np.sin(angle); sin_angle_nyq = np.sin(angle_nyq)  
            
            # Calculate the (R,Z) coordinates  
            R_mn[izeta, :] += rmnc[isurf,:]*cos_angle*weight
            Z_mn[izeta, :] += zmns[isurf,:]*sin_angle*weight 
            
            # Calculate the quantity Q
            if np.all(qmnc[0]!=None): Q_mn[izeta, :] += qmnc[isurf,:]*cos_angle_nyq*weight
            if np.all(qmns[0]!=None): Q_mn[izeta, :] += qmns[isurf,:]*sin_angle_nyq*weight 
                
    # Sum away the (m,n) modes
    R = np.sum(R_mn[:, :], axis=1)
    Z = np.sum(Z_mn[:, :], axis=1)
    Q = np.sum(Q_mn[:, :], axis=1) 
    
    # Return R(zeta,theta); Z(zeta, theta) and Q(zeta, theta)
    return R, Z, Q, zeta, theta_pest, nfp

#--------------------------  
def read_vmecdata(wout_path):

    # Read the netcdf file
    dataset = nc4.Dataset(wout_path, mode='r')
  
    # Poloidal and toroidal mode numbers
    # For most VMECs
    xm = dataset.variables['xm'][:]             # Poloidal mode numbers  xm(mn_mode)   
    xn = dataset.variables['xn'][:]             # Toroidal mode numbers  xn(mn_mode)  
    xm_nyq = dataset.variables['xm_nyq'][:]     # Poloidal mode numbers (Nyquist) xm_nyq(mn_mode_nyq)
    xn_nyq = dataset.variables['xn_nyq'][:]     # Toroidal mode numbers (Nyquist) xn_nyq(mn_mode_nyq) 
    
    # Dimensions
    ns  = int(dataset.variables['ns'][:])               
    nfp = int(dataset.variables['nfp'][:]) 
    mnmax = int(dataset.variables['mnmax'][:])  
    mnmax_nyq = int(dataset.variables['mnmax_nyq'][:]) 
     
    # Read rmnc, zmns and bmnc versus (radius,mn_mode)
    rmnc = dataset.variables['rmnc'][:]     # cosmn component of cylindrical R, full mesh  
    zmns = dataset.variables['zmns'][:]     # sinmn component of cylindrical Z, full mesh 
    bmnc = dataset.variables['bmnc'][:]     # cosmn component of mod-B, half mesh   
    lmns = dataset.variables['lmns'][:]     # cosmn component of lambda, half mesh   
    gmnc = dataset.variables['gmnc'][:]     # cosmn component of jacobian, half mesh   
    
    # Rotational transform 
    iotas = dataset.variables['iotas'][:]
    
    # Stellarator-asymmetric terms 
    lasym = dataset.variables['lasym__logical__'][:]
    if lasym: print('ABORT: <lasym> is not implemented yet.'); sys.exit()
    
    # Check grid sizes
    print(f'  xm: {np.shape(xm)}') 
    print(f'  xn: {np.shape(xn)}') 
    print(f'  xm_nyq: {np.shape(xm_nyq)}') 
    print(f'  xn_nyq: {np.shape(xn_nyq)}') 
    print(f'  rmnc: {np.shape(rmnc)}') 
    print(f'  zmns: {np.shape(zmns)}') 
    print(f'  lmns: {np.shape(lmns)}')
    print(f'  bmnc: {np.shape(bmnc)}') 
    print(f'  gmnc: {np.shape(gmnc)}'); print()
    return xm, xn, xm_nyq, xn_nyq, mnmax, mnmax_nyq, ns, nfp, rmnc, zmns, lmns, bmnc, gmnc, iotas

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
 
#===============================================================================
#                                   TEST CODE                                  #
#===============================================================================

if __name__ == "__main__":   
    wout_path = pathlib.Path('/home/hanne/CIEMAT/RESEARCH/MagneticEquilibria/VmecFiles/wout_asdex.nc')
    wout_path = pathlib.Path('/home/hanne/CIEMAT/RESEARCH/MagneticEquilibria/VmecFiles/wout_tok_cbc_2T.nc')
    wout_path = pathlib.Path('/home/hanne/CIEMAT/RESEARCH/MagneticEquilibria/VmecFiles/wout_w7xr003.nc')
    wout_path = pathlib.Path('/home/hanne/CIEMAT/RESEARCH/MagneticEquilibria/VmecFiles/wout_aug_vacuum.nc')
    _, _, Q, zeta, theta, nfp = get_RZBFromCosineAndSineTransformVMEC(wout_path, svalue=0.7, quantity='B', ntheta=20, nzeta=20)
    