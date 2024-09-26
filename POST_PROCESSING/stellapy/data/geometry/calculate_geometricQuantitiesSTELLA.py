 
import numpy as np

class stella: 
    ''' All quantities here have dimensions nzgrid_stella along z. '''
    
    def __init__(self, verbose=True): 
        
        # Toroidal flux 
        self.sign_torflux = None
        self.sign_toroidal_flux = self.sign_torflux
        
        # Other
        self.twist_and_shift_geo_fac = None
        self.field_period_ratio = None
        self.verbose = verbose
        return 
    
    #*********************************************************************
    #      Initiate the geometry arrays with nzeta = nzgrid*2+1
    #*********************************************************************
    def initiate_arraysForStella(self, nalpha, nzgrid, nzeta):  
        
        # Dimensions of the arrays
        self.nalpha = nalpha
        self.nzgrid = nzgrid
        self.nzeta = nzeta
        
        # Dimensions
        self.zeta = np.zeros((nalpha, nzeta))
        self.theta = np.zeros((nalpha, nzeta))
        
        # Geometry arrays used by stella
        self.bmag = np.zeros((nalpha, nzeta))
        self.bpol = np.zeros((nalpha, nzeta))
        self.btor = np.zeros((nalpha, nzeta))
        self.gradpar_zeta = np.zeros((nalpha, nzeta))
        self.grad_alpha_grad_alpha = np.zeros((nalpha, nzeta))
        self.grad_alpha_grad_psi = np.zeros((nalpha, nzeta))
        self.grad_psi_grad_psi = np.zeros((nalpha, nzeta))
        
        self.gds23 = np.zeros((nalpha, nzeta))
        self.gds24 = np.zeros((nalpha, nzeta))
        self.gds25 = np.zeros((nalpha, nzeta))
        self.gds26 = np.zeros((nalpha, nzeta))
        self.gbdrift_alpha = np.zeros((nalpha, nzeta))
        self.gbdrift0_psi  = np.zeros((nalpha, nzeta))
        self.cvdrift_alpha = np.zeros((nalpha, nzeta))
        self.cvdrift0_psi  = np.zeros((nalpha, nzeta))
        self.theta_vmec = np.zeros((nalpha, nzeta))
        
        self.gbdrift = np.zeros((nalpha, nzeta))
        self.cvdrift = np.zeros((nalpha, nzeta))
        self.gbdrift0 = np.zeros((nalpha, nzeta))
        self.cvdrift0 = np.zeros((nalpha, nzeta))
        
        # These can't be calculated for VMEC!
        self.dgds21dr = np.ones((nalpha, nzeta))*np.nan
        self.dgds22dr = np.ones((nalpha, nzeta))*np.nan
        self.dgds2dr = np.ones((nalpha, nzeta))*np.nan
        
    #********************************************************************* 
    #                      Set the field period ratio
    #*********************************************************************
    def add_dataFromVmecAndZgrid(self, VMEC, ZGRID):
        self.nfp = VMEC.nfp
        self.shat = VMEC.shat
        self.rhoc = VMEC.normalized_toroidal_flux_used
        self.qinp = VMEC.safety_factor_q
        self.alpha = VMEC.alpha
        self.sign_torflux = VMEC.sign_toroidal_flux
        self.torflux_sign = VMEC.sign_toroidal_flux
        self.L_reference = VMEC.L_reference
        self.B_reference = VMEC.B_reference
        self.zed = ZGRID.zed

    #********************************************************************* 
    #                      Set the field period ratio
    #*********************************************************************
    def set_fieldPeriodRatio(self, field_period_ratio):
        self.field_period_ratio = field_period_ratio
        
    #*********************************************************************
    #                   Interpolate the geometric quantities
    #*********************************************************************
    def interpolate_geometric_quantities(self, KGRID, STELLA, VMEC, ZGRID, full_flux_surface=False):
        
        # Interpolate if we read more zpoints from VMEC then stella needs
        if (ZGRID.nzgrid_vmec != ZGRID.nzgrid):
         
            # We must interpolate geometric quantities from (zeta,alpha) grid to
            # (zed,alpha) grid, with zed the arc-length
            if self.verbose: print('stella::vmec_geo::get_total_arc_length')
        
            # Unpack the grid sizes since they are used a lot 
            nalpha = KGRID.nalpha
            zed = ZGRID.zed
            nzgrid = ZGRID.nzgrid
            nzgrid_vmec = ZGRID.nzgrid_vmec
            zgrid_refinement_factor = VMEC.zgrid_refinement_factor
            
            # Initiate the arrays
            zed_domain_size = np.ones((nalpha))*np.nan 
            arc_length = np.ones((nalpha,int(2*nzgrid_vmec+1)))*np.nan
         
            # Note that nzgrid*zgrid_refinement_factor gives the index
            # for the max zeta of the nominal zeta grid
            zetamax_idx = nzgrid*zgrid_refinement_factor
             
            # First we need to get zed(zeta,alpha) defined via 
            # 1 = b . grad z = b . grad zeta * dz/dzeta
            for ia in range(nalpha):
                 
                # This is z(zeta_max) - z(zeta_min) for nominal zeta domain 
                zed_domain_size[ia] = get_total_arc_length(VMEC.gradpar_zeta[ia,int(-zetamax_idx+nzgrid_vmec):int(zetamax_idx+nzgrid_vmec+1)], ZGRID.dzeta_vmec)
                
                # Now get z(zeta)
                zmin = -zed_domain_size[ia]*0.5
                arc_length[ia,:] = get_arc_length_grid(zetamax_idx, nzgrid_vmec, zmin, VMEC.gradpar_zeta[ia,:], ZGRID.dzeta_vmec, nzgrid_vmec)
    
            # Now that we know the min/max values of z corresponding to min/max values
            # of the nominal zeta at each of the alphas, construct a regular z grid
            # Make the max z value on this regular grid to the the maximum over all alpha of z(zeta_max,alpha)
            zmax = np.max(zed_domain_size)*0.5
         
            # Scale zed so that it is arc-length compressed (or expanded) to the range [-pi:pi]
            self.zed_scalefac = VMEC.zed_scalefac = np.pi/zmax
            arc_length = arc_length*VMEC.zed_scalefac
             
            # Interpolate the geometric quantities
            if self.verbose: print('stella::vmec_geo::geo_spline')
            for ia in range(nalpha):
                 
                # Now that we have z(alpha,zeta), interpolate from regular zeta grid 
                # (which is irregular in z) to regular zed grid (irregular in zeta)  
                STELLA.zeta[ia,:]                   = geo_spline (arc_length[ia,:], VMEC.zeta, zed)
                STELLA.bmag[ia,:]                   = geo_spline (arc_length[ia,:], VMEC.bmag[ia,:], zed)
                STELLA.gradpar_zeta[ia,:]           = geo_spline (arc_length[ia,:], VMEC.gradpar[ia,:], zed)
                STELLA.grad_psi_grad_psi[ia,:]      = geo_spline (arc_length[ia,:], VMEC.grad_psi_grad_psi[ia,:], zed)
                STELLA.grad_alpha_grad_alpha[ia,:]  = geo_spline (arc_length[ia,:], VMEC.grad_alpha_grad_alpha[ia,:], zed)
                STELLA.grad_alpha_grad_psi[ia,:]    = geo_spline (arc_length[ia,:], VMEC.grad_alpha_grad_psi[ia,:], zed)
                STELLA.grad_psi_grad_psi[ia,:]      = geo_spline (arc_length[ia,:], VMEC.grad_psi_grad_psi[ia,:], zed)
                STELLA.gds23[ia,:]                  = geo_spline (arc_length[ia,:], VMEC.gds23[ia,:], zed)
                STELLA.gds24[ia,:]                  = geo_spline (arc_length[ia,:], VMEC.gds24[ia,:], zed)
                STELLA.gds25[ia,:]                  = geo_spline (arc_length[ia,:], VMEC.gds25[ia,:], zed)
                STELLA.gds26[ia,:]                  = geo_spline (arc_length[ia,:], VMEC.gds26[ia,:], zed)
                STELLA.gbdrift_alpha[ia,:]          = geo_spline (arc_length[ia,:], VMEC.gbdrift_alpha[ia,:], zed)
                STELLA.gbdrift0_psi[ia,:]           = geo_spline (arc_length[ia,:], VMEC.gbdrift0_psi[ia,:], zed)
                STELLA.cvdrift_alpha[ia,:]          = geo_spline (arc_length[ia,:], VMEC.cvdrift_alpha[ia,:], zed)
                STELLA.cvdrift0_psi[ia,:]           = geo_spline (arc_length[ia,:], VMEC.cvdrift0_psi[ia,:], zed)
                STELLA.theta_vmec[ia,:]             = geo_spline (arc_length[ia,:], VMEC.theta[ia,:], zed) 
                    
                STELLA.bpol[ia,:]                   = geo_spline (arc_length[ia,:], VMEC.bpol[ia,:], zed)
                STELLA.btor[ia,:]                   = geo_spline (arc_length[ia,:], VMEC.btor[ia,:], zed)
         
                # gradpar at this point is b . grad zeta
                # but want it to be b . grad z = b . grad zeta * dz/dzeta
                # we have constructed z so that b . grad z = 1
                # so dz/dzeta = 1 / b . grad zeta
                # gds23 and gds24 involve grad z factors
                # but currently calculated in terms of grad zeta
                # so convert via multiplication with dz/dzeta
                STELLA.gds23[ia,:] = STELLA.gds23[ia,:]/STELLA.gradpar_zeta[ia,:]
                STELLA.gds24[ia,:] = STELLA.gds24[ia,:]/STELLA.gradpar_zeta[ia,:]
         
            # The parallel gradient is one
            STELLA.gradpar = np.ones((int(2*ZGRID.nzgrid+1)))
             
            # We now have geometric coefficients on the alpha-grid as we will be multiplying 
            # this with functions of g and phi we must take care to avoid aliasing, this
            # is accomplished by filtering out the highest third of the wavenumber spectra
            if (full_flux_surface):
                print('vmec_geo::get_vmec_geo::filter_geo_coef')
                for iz in range(2*nzgrid+1):
                    filter_geo_coef (nalpha,STELLA.bmag[:,iz])
                    filter_geo_coef (nalpha,STELLA.grad_alpha_grad_alpha[:,iz])
                    filter_geo_coef (nalpha,STELLA.grad_alpha_grad_psi[:,iz])
                    filter_geo_coef (nalpha,STELLA.grad_psi_grad_psi[:,iz])
                    filter_geo_coef (nalpha,STELLA.gds23[:,iz])
                    filter_geo_coef (nalpha,STELLA.gds24[:,iz])
                    filter_geo_coef (nalpha,STELLA.gds25[:,iz])
                    filter_geo_coef (nalpha,STELLA.gds26[:,iz])
                    filter_geo_coef (nalpha,STELLA.gbdrift_alpha[:,iz])
                    filter_geo_coef (nalpha,STELLA.gbdrift0_psi[:,iz])
                    filter_geo_coef (nalpha,STELLA.cvdrift_alpha[:,iz])
                    filter_geo_coef (nalpha,STELLA.cvdrift0_psi[:,iz]) 
        
        else: 
            STELLA.zeta                  = VMEC.zeta*np.ones((1,KGRID.nalpha))
            STELLA.bmag                  = VMEC.bmag
            STELLA.gradpar               = VMEC.gradpar[0,:]
            STELLA.grad_alpha_grad_alpha = VMEC.grad_alpha_grad_alpha
            STELLA.grad_alpha_grad_psi   = VMEC.grad_alpha_grad_psi
            STELLA.grad_psi_grad_psi     = VMEC.grad_psi_grad_psi
            STELLA.gds23                 = VMEC.gds23
            STELLA.gds24                 = VMEC.gds24
            STELLA.gds25                 = VMEC.gds25
            STELLA.gds26                 = VMEC.gds26
            STELLA.gbdrift_alpha         = VMEC.gbdrift_alpha
            STELLA.gbdrift0_psi          = VMEC.gbdrift0_psi
            STELLA.cvdrift_alpha         = VMEC.cvdrift_alpha
            STELLA.cvdrift0_psi          = VMEC.cvdrift0_psi
            STELLA.thetamod_vmec         = VMEC.theta
     
            # Scale zed so that it is zeta compressed (or expanded) to the range [-pi,pi]
            # this is 1/p from stella JCP paper
            self.zed_scalefac = VMEC.zed_scalefac = VMEC.nfp/VMEC.nfield_periods
        return 

    #*********************************************************************
    #                  Finish the <get_vmec_geo> routine
    #*********************************************************************
    def finish_vmecGeoRoutine(self, VMEC):
         
        # Calculate (b . grad zed) with zed = zeta or arc-length scaled to run from -pi to pi
        self.gradpar = self.gradpar*self.zed_scalefac
        self.gds23 = self.gds23*self.zed_scalefac
        self.gds24 = self.gds24*self.zed_scalefac
     
        # vmec_to_stella_geometry_interface returns psitor/psitor_lcfs as rhoc
        # stella uses rhoc = sqrt(psitor/psitor_lcfs) = rhotor
        self.rhoc = np.sqrt(VMEC.rhoc)
        self.rhotor = self.rhoc
     
        # rho = sqrt(psi_t / psi_{t,LCFS})
        # Bref = 2|psi_LCFS|/a^2
     
        # grho = a * |grad rho| = a * |drho/dpsi_t| * |grad psi_t|
        # = |drho/dpsi_t|*(a^2*Bref) * |grad psi_t|/(a*Bref)
        # = a^2*Bref/(2*rho)/|psi_LCFS| * sqrt(grad_psi_grad_psi)
        # = 1/rho * sqrt(grad_psi_grad_psi)
        self.grho = np.sqrt(self.grad_psi_grad_psi)/self.rhotor
     
        # grho = |grad rho| = |drho/dx| * |grad x|
        # |drho/dx| = L_reference
        # gds22 = shat^2 * |grad x|^2
        # grho = sqrt(gds22/surf%shat**2)/L_reference
        self.drhotordrho = 1.0
        self.psitor_lcfs = 0.5*self.sign_torflux
     
        # This is zeta0 that appears everywhere in tandem with zeta-zeta0
        # converted to the zed coordinate (which is possibly arc-length
        # and is compressed to fit on the range -pi,pi)
        self.zed0_fac = -self.zed[self.nzgrid-1]/self.zeta[0,self.nzgrid-1]*self.qinp
     
        # scale the vmec output
        # alpha = theta_pest - iota*zeta
        # theta_pest = theta_vmec + Lambda(psi,alpha,theta_vmec)
        # with theta_pest a straight-field-line angle
        # but not theta_vmec
        # so theta is theta_pest up to constant (alpha)
     
        # theta = zeta/nfp/surf%qinp
        self.theta = self.alpha*np.ones((2,2*self.nzgrid+1))+self.zeta/self.qinp 
         
        # This is the vmec theta (not straight-field-line coordinate) scaled to run between -pi and pi
        self.theta_vmec = self.theta_vmec/self.nfp
        if True: return 
        
    #*********************************************************************
    #                  Finish the <stella_geometry> routine
    #*********************************************************************
    def finish_stellaGeometryRoutine(self, VMEC, ZGRID):  
        
        # Stella routine
        if self.verbose: print('stella_geometry::init_geometry::init_geometry')
        
        # Signs
        self.dxdXcoord_sign = -1
        self.dydalpha_sign = 1
        
        # dxdXcoord = a*Bref*dx/dpsi = -a*Bref*dx/dpsi_tor = -sign(psi_tor)/rhotor
        self.dxdXcoord = self.dxdXcoord_sign*self.torflux_sign/self.rhotor
         
        # dydalpha = (dy/dalpha) / a = sign(dydalpha)*rhotor 
        self.dydalpha = self.dydalpha_sign*self.rhotor
         
        # drho/dpsiN = -drho/d(rho**2)*(aref**2*Bref/psitor_lcfs) = -1.0/rho
        self.drhodpsi = self.dxdXcoord_sign*self.torflux_sign/self.rhotor
        
        # Grad x grad y at the end of the flux tube
        self.grad_x_grad_y_end = self.grad_alpha_grad_psi[0,-1]*(VMEC.L_reference*VMEC.L_reference*VMEC.B_reference)
         
        # abs(twist_and_shift_geo_fac) is dkx/dky*jtwist
        # minus its sign gives the direction of the shift in kx
        # to be used for the twist-and-shift boundary condition
        # The equation actually reduces to f=-2.*np.pi*shat*iota*field_period_ratio
        if ZGRID.boundary_option_switch=="boundary_option_linked":
            if self.verbose: print('Standard twist and shift BC selected')
            self.twist_and_shift_geo_fac = -2.*np.pi*self.shat*self.drhodpsi*self.dydalpha/(self.qinp*self.dxdXcoord*self.rhotor)*self.field_period_ratio
        if ZGRID.boundary_option_switch=="boundary_option_linked_stellarator":
            if self.verbose: print('Stellarator symmetric twist and shift BC selected')
            self.twist_and_shift_geo_fac_full = 2*(self.rhotor*self.rhotor)*(self.grad_alpha_grad_psi)/(self.grad_psi_grad_psi) 
            if (np.abs(self.grad_x_grad_y_end) <= ZGRID.grad_x_grad_y_zero):
                if self.verbose: print('Using periodic boundary conditions as grad_x_grad_y_end < grad_x_grad_y_zero')
            self.twist_and_shift_geo_fac = self.twist_and_shift_geo_fac_full[0,-1]
            
        # gds2 = |grad y|^2 = |grad alpha|^2*(dy/dalpha)^2
        self.gds2 = self.grad_alpha_grad_alpha*self.dydalpha**2
         
        # gds21 = shat*grad x . grad y = shat*dx/dpsi_t*dy/dalpha*grad alpha . grad psi_t
        # Note: psi = -psi_t and so dx/dpsi = - dx/dpsi_t, which is why there is a minus sign here
        self.gds21 = -self.grad_alpha_grad_psi*self.shat*self.dxdXcoord*self.dydalpha
         
        # gds22 = shat^2*|grad x|^2 = shat^2*|grad psi_t|^2*(dx/dpsi_t)^2
        self.gds22 = self.shat**2*self.grad_psi_grad_psi*self.dxdXcoord**2
         
        # gbdrift_alpha and cvdrift_alpha contain the grad-B and curvature drifts
        # projected onto the grad alpha direction, thus we need the projections on grad y
        self.gbdrift = self.gbdrift_alpha*self.dydalpha
        self.cvdrift = self.cvdrift_alpha*self.dydalpha
         
        # gbdrift0_psi and cvdrift0_psi contain the grad-B and curvature drifts
        # projected onto the grad psi direction, thus we need the projections on grad x
        self.gbdrift0 = self.gbdrift0_psi*self.dxdXcoord
        self.cvdrift0 = self.cvdrift0_psi*self.dxdXcoord
        
        # Synonyms
        self.drhodpsi_psi0 = self.drhodpsi 
        self.bmag_psi0 = self.bmag
        
        # In the Clebsch representation (vec_B = grad alpha x grad psi) we have:
        #   J = 1/(delta ell * vec_B) = 1/B
        # Calculate the Jacobian(ia,iz):
        #   J = 1/( nabla s cdot nabla theta xtimes zeta )
        #     = 1/( drho/dpsi * a*nabla_parallel*z * B )
        self.jacob = 1.0/np.abs(self.drhodpsi*self.gradpar*self.bmag)
        
        # write_geometric_coefficients
        if False:
            print("{0:>30}".format("rhoc"), " ", "{0:<10}".format(self.rhoc))
            print("{0:>30}".format("qinp"), " ", "{0:<10}".format(self.qinp))
            print("{0:>30}".format("shat"), " ", "{0:<10}".format(self.shat))
            print("{0:>30}".format("rhotor"), " ", "{0:<10}".format(self.rhotor))
            print("{0:>30}".format("aref"), " ", "{0:<10}".format(VMEC.L_reference))
            print("{0:>30}".format("bref"), " ", "{0:<10}".format(VMEC.B_reference))
            print("{0:>30}".format("dxdXcoord"), " ", "{0:<10}".format(self.dxdXcoord))
            print("{0:>30}".format("dydalpha"), " ", "{0:<10}".format(self.dydalpha)) 
        
        if True: return
        
#===============================================================================
#                                   zgrid.f90
#===============================================================================
 
def get_total_arc_length(gp, dz): 
    return integrate_zed (dz, 1./gp) 

#--------------------------   
def get_arc_length_grid (nz_max, nzext_max, zboundary, gp, dz, nzgrid_vmec):
    ''' Fortran arrays [-nzext_max, nzext_max] correspond to [0, 2*nzext_max+1]. '''

    # Initiate the array
    zarc = np.ones((int(2*nzgrid_vmec+1)))*np.nan
    zarc[int(-nz_max+nzext_max)] = zboundary 
    
    # Fill the array
    if (nz_max != nzext_max):
        for iz in range(int(-nz_max+nzext_max)):
            zarc[iz] = integrate_zed(dz, 1./gp[iz:int(-nz_max+nzext_max+1)])
            zarc[iz] = zarc[int(-nz_max+nzext_max)] - zarc[iz] 
    for iz in np.arange(int(-nz_max+nzext_max+1), int(2*nzext_max+1)):
        zarc[iz] = integrate_zed(dz, 1./gp[int(-nz_max+nzext_max):iz+1])
        zarc[iz] = zarc[int(-nz_max+nzext_max)] + zarc[iz] 
    return zarc
  
#--------------------------   
def integrate_zed(dz, f, intf=0):
    ''' Trapezoidal rule to integrate in zed''' 
    for iz in range(len(f)-1): 
        intf += dz*(f[iz]+f[iz+1]) 
    return 0.5*intf

#--------------------------
def filter_geo_coef(_):
    print("WARNING: <filter_geo_coef> hasn't been copied yet from stella.")
    return 
  
#===============================================================================
#                                   spl.f90
#===============================================================================

def geo_spline(x,y,xint):  
    ''' Use the f2py wrapper for interpolating tensioned splines. This wraps the 
    curv1() and curv2() functions in Alan K Cline's fitpack. 
    
   curv1
   -----
    Required arguments:
        x : input array of strictly increasing x values in strictly increasing order.
        y : input array of y values.
        sigma : tension value is a float where 0.0 generates a standard cubic spline and large values (e.g. 20.0)
                generate increasingly highly tensioned splines. A typical value might be 1.0. A practical upper value
                of about 40.0 seems to be about the maximum that is handled before the derivates can blow up.
    Optional arguments:
        slp1 : start slope (default=0.0)
        slpn : end slope (default=0.0)
        islpsw : switch indicating whether to use slp1 and slpn slope data or whether
                 to estimate the end slopes:
                 = 0 if slp1 and slpn are to be used,
                 = 1 if slp1 is to be used but not slpn,
                 = 2 if slpn is to be used but not slp1,
                 = 3 (default) if both slp1 and slpn are to be estimated internally.
    Return objects:
        yp : array of derivatives at each x,y.
        ierr : error flag,
                 = 0 for normal return,
                 = 1 if n is less than 2,
                 = 2 if x-values are not strictly increasing.
              
   curv2
   -----
    Required arguments:
        xt : evaluate spline at this x-value.
        yp : result of calling curv1 above.
        x, y, sigma : must be the same as used when curv1 was called.
    Return objects:
       yt : y-value evaluated at xt
    '''
    
    # Load the package 
    try: from stellapy.calculations.utils import fitpack as fitpack; loaded_fitpack = True
    except: loaded_fitpack = False
    
    # Use fitpack to exactly mirror stella:
    if loaded_fitpack:
    
        # Initiate yint
        yint = np.ones((len(xint)))*np.nan
        
        # Sigma is one in stella
        sigma = 1.0
        
        # Combines the curv1() and curv2() calls for a real number
        if len(x)==1:
            yp,_ = fitpack.curv1(x,y,sigma,islpsw=3)                #@UndefinedVariable
            yint = fitpack.curv2(xint,x,y,yp,sigma)                 #@UndefinedVariable
        else: 
            yp,_ = fitpack.curv1(x,y,sigma,islpsw=3)                #@UndefinedVariable
            for ix in range(len(xint)):
                yint[ix] = fitpack.curv2(xint[ix],x,y,yp,sigma)     #@UndefinedVariable 
        return yint
    
    # Use python interpolations if fitpack is missing:
    if not loaded_fitpack:  
        from scipy.interpolate import interp1d 
        f = interp1d(x, y, kind='cubic')
        yint = f(xint)
        return yint
        



