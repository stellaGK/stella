
import sys
import numpy as np 
import netCDF4 as nc4  

class vmec:
 
    #*********************************************************************
    #                   Read default vmec parameters
    #*********************************************************************
    def __init__(self, verbose=True):
        
        # Stella routine
        self.verbose = verbose
        if self.verbose: print("stella::vmec_geo::init_vmec_defaults")    
            
        # Constants
        self.mu_0 = 4*np.pi*10**(-7)

        # Default input parameters
        self.vmec_filename = 'equilibria/wout_w7x_standardConfig.nc'
        self.alpha0 = 0.0
        self.zeta_center = 0.0
        self.nfield_periods = -1.0
        self.torflux = 0.6354167 
        self.surface_option = 0 
        return 

    #*********************************************************************
    #               Replace defaults by input parameters
    #*********************************************************************
    def update_inputParameters(self, ZGRID, vmec_filename=None, rho=None, nfield_periods=None, alpha0=None):
    
        # Stella routine
        if self.verbose: print("stella::vmec_geo::read_vmec_parameters")
        
        # Overwrite default parameters
        if vmec_filename!=None:  self.vmec_filename = vmec_filename
        if nfield_periods!=None: self.nfield_periods = nfield_periods
        if alpha0!=None:         self.alpha0 = alpha0
        if rho!=None:            self.torflux = rho*rho

        # Add some synonyms
        self.rho = np.sqrt(self.torflux)
        self.svalue = self.torflux
        self.rhotor = self.torflux
        self.desired_normalized_toroidal_flux = self.torflux 
        
        # When simulating an entire flux surface, we must obtain vmec geo quantities
        # on a zeta grid that is longer than the zeta grid that will ultimately be used
        # in simulation. This is related to need for gradpar to be independent of alpha.
        if (ZGRID.zed_equal_arc):
            self.zgrid_scalefac = 2.0
            self.zgrid_refinement_factor = 4
        else:
            self.zgrid_scalefac = 1.0
            self.zgrid_refinement_factor = 1 
            
        # Check the zgrid scaling factor
        if (self.zgrid_scalefac < 1.0):
            if self.verbose: print('zgrid_scalefac = ', self.zgrid_scalefac)
            if self.verbose: print('zgrid_scalefac should always be >= 1.0.  aborting'); exit
        elif (not ZGRID.zed_equal_arc):
            if (self.zgrid_scalefac > 1.0):
                if self.verbose: print('There is no reason to use zgrid_scalefac different from 1.0 unless zed_equal_arc=T')
                if self.verbose: print('Setting zgrid_scalefac = 1.0')
                self.zgrid_scalefac = 1.0
            elif (self.zgrid_refinement_factor > 1):
                if self.verbose: print('There is no reason to use zgrid_refinement_factor > 1 unless zed_equal_arc=T')
                if self.verbose: print('Setting zgrid_refinement_factor = 1')
                self.zgrid_refinement_factor = 1
        return 

    #*********************************************************************
    #                       Basic calculations
    #*********************************************************************   
    def calculate_toroidalFluxAtEdge(self, ns):
        ''' The toroidal flux at the edge is usually denoted by 2*pi*psi_LCFS.'''
        self.edge_toroidal_flux_over_2pi = self.phi[ns-1] / (2*np.pi)*self.isigng 
        self.sign_toroidal_flux = np.sign(self.edge_toroidal_flux_over_2pi)
        self.sign_torflux = self.sign_toroidal_flux
        return
    
    #*********************************************************************
    #      Reference units
    #*********************************************************************
    def calculate_referenceUnits(self): 
        ''' The reference length in stella is the minor radius.
        The reference magnetic field is 2*psi_LCFS/(2*pi*L**2). 
        Using the choices made by Pavlos Xanthopoulos in GIST. '''
        self.L_reference = self.Aminor
        self.B_reference = 2*abs(self.edge_toroidal_flux_over_2pi) / (self.L_reference*self.L_reference)
        return 
    
    #*********************************************************************
    #      Initiate the geometry arrays with nzeta = nzgrid_vmec*2+1
    #*********************************************************************
    def initiate_arraysForVmec(self, nalpha, nzgrid, nzeta):  
        
        # Dimensions of the arrays
        self.nzgrid = nzgrid
        self.nzeta  = nzeta
        self.nalpha = nalpha         
        
        # Main dimensions 
        self.alpha = np.empty((self.nalpha))
        self.zeta = np.empty((self.nzeta))
        self.theta_vmec = np.empty((self.nalpha, self.nzeta))
        self.theta_pest = np.empty((self.nalpha, self.nzeta))
        
        # Magnetic field along (alpha, zeta)
        self.B = np.zeros((nalpha, nzeta))
        self.B_sub_s = np.zeros((nalpha, nzeta))
        self.B_sub_zeta = np.zeros((nalpha, nzeta))
        self.B_sub_theta_vmec = np.zeros((nalpha, nzeta))
        self.B_sup_theta_vmec = np.zeros((nalpha, nzeta))
        self.B_sup_zeta = np.zeros((nalpha, nzeta))
        
        # Derivatives of the magnetic field along (alpha, zeta)
        self.d_B_d_theta_vmec = np.zeros((nalpha, nzeta))
        self.d_B_d_zeta = np.zeros((nalpha, nzeta))
        self.d_B_d_s = np.zeros((nalpha, nzeta))
        
        # Magnetic field in (X,Y,Z) coordinates along (alpha, zeta)
        self.B_X = np.zeros((nalpha, nzeta))  
        self.B_Y = np.zeros((nalpha, nzeta))  
        self.B_Z = np.zeros((nalpha, nzeta))  
        self.grad_B_X = np.zeros((nalpha, nzeta))  
        self.grad_B_Y = np.zeros((nalpha, nzeta))  
        self.grad_B_Z = np.zeros((nalpha, nzeta))  
        
        # Magnetic field cross products along (alpha, zeta)
        self.B_cross_grad_B_dot_grad_alpha = np.zeros((nalpha, nzeta))  
        self.B_cross_grad_B_dot_grad_alpha_alternate = np.zeros((nalpha, nzeta))  
        self.B_cross_grad_s_dot_grad_alpha = np.zeros((nalpha, nzeta))  
        self.B_cross_grad_s_dot_grad_alpha_alternate = np.zeros((nalpha, nzeta))
    
        # Jacobian 
        self.sqrt_g = np.zeros((nalpha, nzeta))   

        # Gradient products
        self.grad_alpha_grad_alpha = None
        self.grad_alpha_grad_psi = None
        self.grad_psi_grad_psi = None
        
        # Gradients
        self.grad_s_X = np.zeros((nalpha, nzeta)) 
        self.grad_s_Y = np.zeros((nalpha, nzeta)) 
        self.grad_s_Z = np.zeros((nalpha, nzeta)) 
        self.grad_theta_vmec_X = np.zeros((nalpha, nzeta)) 
        self.grad_theta_vmec_Y = np.zeros((nalpha, nzeta)) 
        self.grad_theta_vmec_Z = np.zeros((nalpha, nzeta)) 
        self.grad_theta_pest_X = np.zeros((nalpha, nzeta)) 
        self.grad_theta_pest_Y = np.zeros((nalpha, nzeta)) 
        self.grad_theta_pest_Z = np.zeros((nalpha, nzeta)) 
        self.grad_zeta_X = np.zeros((nalpha, nzeta)) 
        self.grad_zeta_Y = np.zeros((nalpha, nzeta)) 
        self.grad_zeta_Z = np.zeros((nalpha, nzeta)) 
        self.grad_psi_X = np.zeros((nalpha, nzeta)) 
        self.grad_psi_Y = np.zeros((nalpha, nzeta)) 
        self.grad_psi_Z = np.zeros((nalpha, nzeta)) 
        self.grad_alpha_X = np.zeros((nalpha, nzeta)) 
        self.grad_alpha_Y = np.zeros((nalpha, nzeta)) 
        self.grad_alpha_Z = np.zeros((nalpha, nzeta))  

        # Gradient products
        self.gradzeta_grady = np.zeros((nalpha, nzeta))  
        self.gradzeta_gradx = np.zeros((nalpha, nzeta))  
        self.gradtheta_grady = np.zeros((nalpha, nzeta))  
        self.gradtheta_gradx = np.zeros((nalpha, nzeta))   

        # Signs
        self.dxdXcoord_sign = -1
        self.dydalpha_sign = 1
        
        # Coordinate derivatives
        self.dxdXcoord = None
        self.dydalpha = None
        self.drhodpsi = None
        self.drhodpsi_psi0 = self.drhodpsi
        
        # (R,Z,Lambda) coordinates
        self.R = np.zeros((nalpha, nzeta))
        self.d_R_d_theta_vmec = np.zeros((nalpha, nzeta))
        self.d_R_d_zeta = np.zeros((nalpha, nzeta))
        self.d_R_d_s = np.zeros((nalpha, nzeta))
        self.d_Z_d_theta_vmec = np.zeros((nalpha, nzeta))
        self.d_Z_d_zeta = np.zeros((nalpha, nzeta))
        self.d_Z_d_s = np.zeros((nalpha, nzeta))
        self.d_Lambda_d_theta_vmec = np.zeros((nalpha, nzeta))
        self.d_Lambda_d_zeta = np.zeros((nalpha, nzeta))
        self.d_Lambda_d_s = np.zeros((nalpha, nzeta)) 
        
        # (X,Y) coordinates
        self.d_X_d_s = np.zeros((nalpha, nzeta))
        self.d_X_d_theta_vmec = np.zeros((nalpha, nzeta))
        self.d_X_d_zeta = np.zeros((nalpha, nzeta))
        self.d_Y_d_s = np.zeros((nalpha, nzeta))
        self.d_Y_d_theta_vmec = np.zeros((nalpha, nzeta))
        self.d_Y_d_zeta = np.zeros((nalpha, nzeta))  
    
    #*********************************************************************
    #             Set the used normalized toroidal flux
    #*********************************************************************
    def set_normalizedToroidalFluxUsed(self, torflux):
        self.normalized_toroidal_flux_used = torflux
        return
    
    #*********************************************************************
    #          Determine the number of field periods
    #*********************************************************************
    def get_numberOfFieldPeriods(self):
        if self.verbose: print("stella::vmec_geo::get_modified_vmec_zeta_grid-->call get_nominal_vmec_zeta_grid(...)")
        self.number_of_field_periods_device = self.nfp_device = self.nfp
        self.number_of_field_periods_stella = self.nfp_stella = self.nfield_periods
        if self.verbose: print("stella::vmec_geo::get_vmec_geo-->call vmec_to_stella_geometry_interface(...)")
        self.number_of_field_periods_to_include = self.nfp_include = self.nfield_periods*self.zgrid_scalefac
        if self.verbose: print("stella::vmec_to_stella_geometry_interface-->Set up the coordinate grids.")
        self.number_of_field_periods_to_include_final = self.nfp_final = self.number_of_field_periods_to_include
        if (self.number_of_field_periods_to_include <= 0): self.number_of_field_periods_to_include_final = self.nfp_final = self.nfp 
        return

    #*********************************************************************
    #                Set up the coordinate grids for VMEC
    #*********************************************************************
    def get_alphaGridForVmec(self):
        ''' Alpha coordinates. '''
        if self.verbose: print("stella::vmec_to_stella_geometry_interface-->Set up the coordinate grids.")
        self.alpha = [ self.alpha0 + ((i-1)*2*np.pi)/self.nalpha for i in np.arange(1, self.nalpha+1, 1)]
        return
    
    #----------------------
    def get_zetaGridForVmec(self):
        ''' Zeta coordinates. '''   
        if self.verbose: print("stella::vmec_to_stella_geometry_interface-->Set up the coordinate grids.")
        self.zeta = [ self.zeta_center + (np.pi*i*self.nfp_final)/(self.nfp*self.nzgrid) for i in np.arange(-self.nzgrid, self.nzgrid+1, 1)]
        return
    
    #----------------------
    def get_thetaGridForVmec(self):
        ''' We need to determine theta_vmec = theta_pest - Lambda. 
        Knowing that theta_pest = alpha + iota*zeta. '''
        if self.verbose: print("stella::vmec_to_stella_geometry_interface-->Beginning root solves to determine theta_self.")
        for izeta in range(self.nzeta): 
            for ialpha in range(self.nalpha):
                # Current zeta
                zeta0 = self.zeta[izeta]
                # The theta pest angle should equal alpha + iota*zeta
                theta_pest_target = self.alpha[ialpha] + self.iota*zeta0; 
                # Guess that theta_vmec will be within 0.3 radians of theta_pest:
                theta_vmec_min = theta_pest_target - 0.3
                theta_vmec_max = theta_pest_target + 0.3
                # Calculate theta in the vmec file
                self.theta_pest[ialpha, izeta] = theta_pest_target
                self.theta_vmec[ialpha, izeta], theta_converged = self.get_root(theta_vmec_min, theta_vmec_max, theta_pest_target, zeta0)
                if (not theta_converged): print("ERROR: could not find root needed to compute theta_vmec. aborting"); exit
        self.theta = self.theta_vmec


    def get_root(self, a0, b0, theta_pest_target, zeta0): 
        """ Apply the bisection method, which is a root-finding method that applies 
        to any continuous functions for which one knows two values with opposite signs. 
        The method consists of repeatedly bisecting the interval defined by these values 
        and then selecting the subinterval in which the function changes sign, and 
        therefore must contain a root.
        
         We know that (theta_vmec + Lambda) = theta_pest = theta_pest_vmec
         We want to find for which <theta_vmec> we can most accuratly get <theta_pest_target>.
         Therefore calculate the difference between <theta_pest_vmec> and <theta_pest_target>:
            f = residual = (theta_vmec_try + Lambda) - theta_pest_target 
        And apply the bisection method to minimalize this difference.
         """
         
        # We can increase or interval [theta_pest_target-0.3, theta_pest_target+0.3] up to ten times
        # in order to find the first occurance of a root (f=0 within the interval).
        itmax_bracket = 10
        
        # We can narrow down our interval up to 10 times to find the point at which f=0
        # In other words we try to find <theta_pest_vmec> as close as possible to <theta_pest_target> 
        itmax_root = 10
        
        # If <theta_pest_vmec> and <theta_pest_target> differ less than <tol> we are content
        tol = 1.0e-10
        
        # The first two attempts are done with (a0=theta_vmec_min) and (b0=theta_vmec_max)
        a = a0
        b = b0
        
        # Calculate the difference between theta_pest_vmec and theta_pest_target for these two first attemps
        fa = self.fzero_residual(a, theta_pest_target, zeta0)
        fb = self.fzero_residual(b, theta_pest_target, zeta0)
        
        # When f(a), f(b) have opposite signs, it it said that they "bracket the root" 
        # since in a drawing the function f will pass through zero between points (a,b)
        # If however this is not the case, we need to put our original points (a,b) further apart.
        for _ in range(itmax_bracket):
            eps = sys.float_info.epsilon
            if ((fa > 0.0 and fb > 0.0) or (fa < 0.0 and fb < 0.0)):
                if self.verbose: print('\nIn vmec_to_stella_geometry_interface, theta_min='+str(a)+' and theta_max='+str(b)+' do not bracket root.')
                if self.verbose: print('f(theta_min)='+str(fa)+'and f(theta_max)='+str(fb)+'.')
                a = a-0.3
                b = b+0.3
                if self.verbose: print('Trying again with values '+str(a)+' and '+str(b)+' .')
                fa = self.fzero_residual(a, theta_pest_target, zeta0)
                fb = self.fzero_residual(b, theta_pest_target, zeta0)
            else: 
                break
        
        # Assume that now the points (a,b) are such that f(a) and f(b) have opposite
        # sign, and thus the root (or f=0 point) lies in the interval [a,b]. Now
        # narrow down the interval to approach the f=0 point.
        c = b
        fc = fb
        for _ in range(itmax_bracket):
            if ((fb > 0.0 and fc > 0.0) or (fb < 0.0 and fc < 0.0)):
                c = a
                fc = fa
                d = b-a
                e = d 
            if (abs(fc) < abs(fb)):
                a = b
                b = c
                c = a
                fa = fb
                fb = fc
                fc = fa 
            tol1 = 2.0*eps*abs(b)+0.5*tol
            xm = 0.5*(c-b)
            
            # If we found exactly f=0 in point b then we return this theta_vmec 
            # Or if <theta_pest_vmec> and <theta_pest_target> lie close enough to each other
            if (abs(xm) <= tol1 or fb == 0.0):
                root = b
                converged = True
                return root, converged
            
            # Otherwise, narrow down the interval to get closer to f=0
            if (abs(e) >= tol1 and abs(fa) > abs(fb)):
                s = fb/fa
                if (a==c):
                    p = 2.0*xm*s
                    q = 1.0-s
                else:
                    q = fa/fc
                    r = fb/fc
                    p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
                    q = (q-1.0)*(r-1.0)*(s-1.0) 
                if (p > 0.0): q = -q
                p=abs(p)
                if (2.0*p < min(3.0*xm*q-abs(tol1*q),abs(e*q))):
                    e = d
                    d = p/q
                else:
                    d = xm
                    e = d 
            else:
                d=xm
                e=d 
            a = b
            fa = fb
            b = b + d if np.abs(d) > tol1 else b + tol1*np.sign(xm) 
            fb = self.fzero_residual(b, theta_pest_target, zeta0)
        return root, converged

    #--------------------------  
    def fzero_residual(self, theta_vmec_try, theta_pest_target, zeta0):
        ''' Note that lmns and lmnc use the non-Nyquist xm, xn, and mnmax.
        Also note that lmns and lmnc are on the HALF grid.
        residual = (theta_pest based on theta_vmec_try) - theta_pest_target = theta_vmec_try + Lambda - theta_pest_target
        The variable x plays the role of theta_vmec_try. '''

        # Not implemented for assymetric fields yet
        if self.lasym and self.verbose: print("lasym == True has not been implemented yet."); return 
         
        # Residual = (theta_pest based on theta_vmec_try) - theta_pest_target = theta_vmec_try + Lambda - theta_pest_target
        fzero_residual = theta_vmec_try - theta_pest_target
         
        # Iterate over the modes
        for imn in range(self.mnmax):
             
            # Get the angle 
            sin_angle = np.sin(self.xm[imn]*theta_vmec_try - self.xn[imn]*zeta0)
             
            # Interpolate the calculation between the flux surfaces
            for interpolation in [0,1]:
                 
                # When svalue denotes a flux surface between two flux surfaces
                # that are given by vmec, we need to interpolate the values
                weight = self.surfaces_weights[interpolation]
                isurf  = self.surfaces_indexes[interpolation] 
                 
                # Calculate fzero_residual
                fzero_residual = fzero_residual + sin_angle*weight*self.lmns[isurf,imn]
        return fzero_residual

    #*********************************************************************
    #  Calculate the normalized toroidal flux on the full and half grids
    #*********************************************************************
    def calculate_normalizedToroidalFluxOnFullAndHalfGrids(self, ns):
        psiN_full = [ (interpolation-1)/[ns-1] for interpolation in np.arange(1,ns+1,1) ] 
        psiN_half = [ (psiN_full[i] + psiN_full[i+1])*0.5 for i in range(len(psiN_full)-1)]
        self.normalized_toroidal_flux_full_grid = psiN_full
        self.normalized_toroidal_flux_half_grid = psiN_half
        return 
     
    #*********************************************************************
    #                    Read the vmec equilibrium
    #*********************************************************************
    def read_vmec_equilibrium(self):
        ''' <phi> is vmec's array of the toroidal flux (not divided by 2pi#) 
        on vmec's radial grid. '''
        
        # Stella routine
        if self.verbose: print('stella::vmec_geo::read_vmec_equilibrium')
        
        # Read the vmec netcdf file
        dataset = nc4.Dataset(self.vmec_filename, mode='r')
    
        # Basic quantities
        self.Aminor = dataset.variables['Aminor_p'][:] 
             
        # Poloidal and toroidal mode numbers
        self.xm = dataset.variables['xm'][:]             # Poloidal mode numbers  xm(mn_mode)   
        self.xn = dataset.variables['xn'][:]             # Toroidal mode numbers  xn(mn_mode)  
        self.xm_nyq = dataset.variables['xm_nyq'][:]     # Poloidal mode numbers (Nyquist) xm_nyq(mn_mode_nyq)
        self.xn_nyq = dataset.variables['xn_nyq'][:]     # Toroidal mode numbers (Nyquist) xn_nyq(mn_mode_nyq)
           
        # Dimensions
        self.ns  = int(dataset.variables['ns'][:])               
        self.nfp = int(dataset.variables['nfp'][:]) 
        self.mpol = int(dataset.variables['mpol'][:])     
        self.ntor = int(dataset.variables['ntor'][:])     
        self.mnmax = int(dataset.variables['mnmax'][:]) 
        self.mnmax_nyq = int(dataset.variables['mnmax_nyq'][:]) 
          
        # Read the vectors on the <ns> flux surfaces
        self.phi = dataset.variables['phi'][:]           # Toroidal flux on full mesh 
        self.phip = dataset.variables['phips'][:]        # d(phi)/ds: Toroidal flux deriv on full mesh
        self.iotas = dataset.variables['iotas'][:]       # iota on half mesh 
        self.iotaf = dataset.variables['iotaf'][:]       # iota on full mesh 
        self.presf = dataset.variables['presf'][:]
            
        # Read rmnc, zmns and bmnc versus (radius,mn_mode)
        self.rmnc = dataset.variables['rmnc'][:]     # cosmn component of cylindrical R, full mesh  
        self.zmns = dataset.variables['zmns'][:]     # sinmn component of cylindrical Z, full mesh 
        self.bmnc = dataset.variables['bmnc'][:]     # cosmn component of mod-B, half mesh   
        self.lmns = dataset.variables['lmns'][:]     # sinmn component of lambda, half mesh   
        self.gmnc = dataset.variables['gmnc'][:]     # cosmn component of jacobian, half mesh 
        self.bsupumnc = dataset.variables['bsupumnc'][:] 
        self.bsupvmnc = dataset.variables['bsupvmnc'][:] 
        self.bsubumnc = dataset.variables['bsubumnc'][:]
        self.bsubvmnc = dataset.variables['bsubvmnc'][:]
        self.bsubsmns = dataset.variables['bsubsmns'][:]
            
        # Others
        self.lasym = dataset.variables['lasym__logical__'][:]     
        self.isigng = dataset.variables['signgs'][:]
        
        # Initiate the arrays for the derivatives of the sinus/consinus 
        # transformations, which will be filled for each <imn_nyq> and <imn>
        # Sinus/cosinus transformation of derivatives of the (m,n) modes
        self.d_B_d_s_mnc = np.zeros((self.ns))
        self.d_B_d_s_mns = np.zeros((self.ns))
        self.d_R_d_s_mnc = np.zeros((self.ns))
        self.d_R_d_s_mns = np.zeros((self.ns))
        self.d_Z_d_s_mnc = np.zeros((self.ns))
        self.d_Z_d_s_mns = np.zeros((self.ns))
        self.d_Lambda_d_s_mnc = np.zeros((self.ns))
        self.d_Lambda_d_s_mns = np.zeros((self.ns))
        if True: return 
    
#===============================================================================
#                               CALCULATIONS
#===============================================================================

    #*********************************************************************
    #   Evaluate several radial-profile functions at the flux surface
    #*********************************************************************  
    def calculate_radialProfileFunctionsAtFluxSurface(self, ns):
    
        # Stella routine
        if self.verbose: print("stella::vmec_to_stella_geometry_interface-->Evaluate several radial-profile functions")
    
        # Calculate iota and the safety factor for the chosen flux surface
        self.iota  = self.iotas[self.radial_index_half[0]]*self.radial_weight_half[0] 
        self.iota += self.iotas[self.radial_index_half[1]]*self.radial_weight_half[1] 
        self.safety_factor_q = 1/self.iota
      
        # Step along s=rho**2
        self.ds = self.normalized_toroidal_flux_full_grid[1] - self.normalized_toroidal_flux_full_grid[0]
          
        # Derive iota along s=rho**2
        d_iota_d_s_on_half_grid = np.zeros(ns)
        d_iota_d_s_on_half_grid[1:ns-1] = (self.iotaf[1:ns-1] - self.iotaf[0:ns-2]) / self.ds
        self.d_iota_d_s  = d_iota_d_s_on_half_grid[self.radial_index_half[0]]*self.radial_weight_half[0] 
        self.d_iota_d_s += d_iota_d_s_on_half_grid[self.radial_index_half[1]]*self.radial_weight_half[1]
        del d_iota_d_s_on_half_grid
          
        # shat = (r/q)(dq/dr) where r = a sqrt(s) = - (r/iota) (d iota / d r) = -2 (s/iota) (d iota / d s)
        self.shat = (-2*self.normalized_toroidal_flux_used / self.iota)*self.d_iota_d_s
      
        # Derive pressure along s
        d_pressure_d_s_on_half_grid = np.zeros(ns) 
        d_pressure_d_s_on_half_grid[1:ns-1] = (self.presf[1:ns-1] - self.presf[0:ns-2]) / self.ds
        self.d_pressure_d_s  = d_pressure_d_s_on_half_grid[self.radial_index_half[0]]*self.radial_weight_half[0]
        self.d_pressure_d_s += d_pressure_d_s_on_half_grid[self.radial_index_half[1]]*self.radial_weight_half[1]
        del d_pressure_d_s_on_half_grid
        
        # Synonyms
        self.rhoc = self.normalized_toroidal_flux_used
        self.qinp = self.safety_factor_q
        return
        
    #*********************************************************************
    #     Get the radial derivatives of the sine/cosine transformation 
    #*********************************************************************
    def get_radialDerivativesOfMnModes(self, ds, ns, imn, imn_nyq, non_Nyquist_mode_available):

        # Stella routine
        if imn_nyq==0 and imn==0: 
            if self.verbose: print("stella::vmec_to_stella_geometry_interface-->Evaluate the radial derivatives we will need")
        
        # Get the derivative of bmnc with respect to s and do a simplistic extrapolation at the endpoints
        self.d_B_d_s_mnc[1:ns-1] = (self.bmnc[2:ns,imn_nyq] - self.bmnc[1:ns-1,imn_nyq]) / ds
        self.d_B_d_s_mnc[0] = self.d_B_d_s_mnc[1]
        self.d_B_d_s_mnc[ns-1] = self.d_B_d_s_mnc[ns-2]
    
        # Get the derivative of rmns; zmns and lmns with respect to s
        if (non_Nyquist_mode_available):
          
            # R is on the full mesh: 
            self.d_R_d_s_mnc[1:ns] = (self.rmnc[1:ns,imn] - self.rmnc[0:ns-1,imn]) / ds
            self.d_R_d_s_mnc[0] = 0
          
            # Z is on the full mesh:
            self.d_Z_d_s_mns[1:ns] = (self.zmns[1:ns,imn] - self.zmns[0:ns-1,imn]) / ds
            self.d_Z_d_s_mns[0] = 0
          
            # Lambda is on the half mesh:
            self.d_Lambda_d_s_mns[1:ns-1] = (self.lmns[2:ns,imn] - self.lmns[1:ns-1,imn]) / ds 
            self.d_Lambda_d_s_mns[0] = self.d_Lambda_d_s_mns[1]
            self.d_Lambda_d_s_mns[ns-1] = self.d_Lambda_d_s_mns[ns-2]  
          
        else: 
            self.d_R_d_s_mnc[:] = 0
            self.d_Z_d_s_mns[:] = 0
            self.d_Lambda_d_s_mns[:] = 0 
        
        # Now consider the stellarator-asymmetric terms.
        if (self.lasym):
 
            # B and Lambda are on the half mesh, so their radial derivatives are on the full mesh.
            # R and Z are on the full mesh, so their radial derivatives are on the half mesh.
            self.d_B_d_s_mns[1:ns-2] = (self.bmns[2:ns-1,imn_nyq] - self.bmns[1:ns-2,imn_nyq]) / ds 
            self.d_B_d_s_mns[0] = self.d_B_d_s_mns[1]
            self.d_B_d_s_mns[ns-1] = self.d_B_d_s_mns[ns-2]
            
            if (non_Nyquist_mode_available):
                
                # R is on the full mesh:
                self.d_R_d_s_mns[1:ns-1] = (self.rmns[1:ns-1,imn] - self.rmns[0:ns-2,imn]) / ds
                self.d_R_d_s_mns[0] = 0
                
                # Z is on the full mesh:
                self.d_Z_d_s_mnc[1:ns-1] = (self.zmnc[1:ns-1,imn] - self.zmnc[0:ns-2,imn]) / ds
                self.d_Z_d_s_mnc[0] = 0
                
                # Lambda is on the half mesh:
                self.d_Lambda_d_s_mnc[1:ns-2] = (self.lmnc[2:ns-1,imn] - self.lmnc[1:ns-2,imn]) / ds 
                self.d_Lambda_d_s_mnc[0] = self.d_Lambda_d_s_mnc[1]
                self.d_Lambda_d_s_mnc[ns-1] = self.d_Lambda_d_s_mnc[ns-2]
                
            else:
                self.d_R_d_s_mns = 0
                self.d_Z_d_s_mnc = 0
                self.d_Lambda_d_s_mnc = 0

        return  

    #*********************************************************************
    #           Add contributions to the magnetic quantities
    #*********************************************************************
    def add_contributionsToMagneticQuantities(self, ialpha, izeta, m, n, imn_nyq, scale_factor):

        # Stella routine
        if imn_nyq==0 and ialpha==0 and izeta==0: 
            if self.verbose: print("stella::vmec_to_stella_geometry_interface-->! Handle |B|:    (BFIELD)")

        # Get the angle and their cosinus and sinus  
        angle = m*self.theta[ialpha,izeta] - n*self.nfp*self.zeta[izeta]
        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)
                  
        # Iterate over the interpolation weights
        for interpolation in [0,1]:
            
            # When svalue denotes a flux surface between two flux surfaces
            # that are given by vmec, we need to interpolate the values
            weight  = self.radial_weight_half[interpolation]*scale_factor
            isurf   = self.radial_index_half[interpolation]
            weightf = self.radial_weight_full[interpolation]*scale_factor
            isurff  = self.radial_index_full[interpolation]
                
            # Handle |B|:
            temp = self.bmnc[isurf,imn_nyq]*weight 
            self.B[ialpha,izeta] += temp*cos_angle
            self.d_B_d_theta_vmec[ialpha,izeta] -= m*temp*sin_angle
            self.d_B_d_zeta[ialpha,izeta] += n*self.nfp*temp*sin_angle       
              
            # Handle Jacobian:
            self.sqrt_g[ialpha,izeta] += self.gmnc[isurf,imn_nyq]*weight*cos_angle
              
            # Handle B sup theta: 
            self.B_sup_theta_vmec[ialpha,izeta] += self.bsupumnc[isurf,imn_nyq]*weight*cos_angle
              
            # Handle B sup zeta: 
            self.B_sup_zeta[ialpha,izeta] += self.bsupvmnc[isurf,imn_nyq]*weight*cos_angle
              
            # Handle B sub theta: 
            self.B_sub_theta_vmec[ialpha,izeta] += self.bsubumnc[isurf,imn_nyq]*weight*cos_angle
              
            # Handle B sub zeta: 
            self.B_sub_zeta[ialpha,izeta] += self.bsubvmnc[isurf,imn_nyq]*weight*cos_angle
              
            # Handle B sub psi. Unlike the other components of B, this one is on the full mesh.
            self.B_sub_s[ialpha,izeta] += self.bsubsmns[isurff,imn_nyq]*weightf*sin_angle
              
            # Handle d B / d s. Since bmnc is on the half mesh, its radial derivative is on the full mesh. 
            self.d_B_d_s[ialpha,izeta] += self.d_B_d_s_mnc[isurff]*weightf*cos_angle
        return 

    #*********************************************************************
    #           Add contributions to the coordinate quantities
    #*********************************************************************
    def add_contributionsToCoordinateQuantities(self, ialpha, izeta, m, n, imn, scale_factor):

        # Stella routine
        if imn==0 and ialpha==0 and izeta==0: 
            if self.verbose: print("stella::vmec_to_stella_geometry_interface-->! Handle |B|:    (COORD)")

        # Get the angle and their cosinus and sinus  
        angle = m*self.theta[ialpha,izeta] - n*self.nfp*self.zeta[izeta]
        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)
                  
        # Iterate over the interpolation weights
        for interpolation in [0,1]: 
            
            # When svalue denotes a flux surface between two flux surfaces
            # that are given by vmec, we need to interpolate the values
            weight  = self.radial_weight_half[interpolation]*scale_factor
            isurf   = self.radial_index_half[interpolation]
            weightf = self.radial_weight_full[interpolation]*scale_factor
            isurff  = self.radial_index_full[interpolation]
            
            # Handle R, which is on the full mesh
            temp = self.rmnc[isurff,imn]*weightf 
            self.R[ialpha,izeta] += temp*cos_angle
            self.d_R_d_theta_vmec[ialpha,izeta] -= temp*m*sin_angle
            self.d_R_d_zeta[ialpha,izeta] += temp*n*self.nfp*sin_angle
              
            # Handle Z, which is on the full mesh (Z[ialpha,izeta] += temp*sin_angle)
            temp = self.zmns[isurff,imn]*weightf 
            self.d_Z_d_theta_vmec[ialpha,izeta] += temp*m*cos_angle
            self.d_Z_d_zeta[ialpha,izeta] -= temp*n*self.nfp*cos_angle
              
            # Handle Lambda:
            temp = self.lmns[isurf,imn]*weight 
            self.d_Lambda_d_theta_vmec[ialpha,izeta] += m*temp*cos_angle
            self.d_Lambda_d_zeta[ialpha,izeta] -= n*self.nfp*temp*cos_angle       
              
            # Handle d R / d s. Since R is on the full mesh, its radial derivative is on the half mesh.
            self.d_R_d_s[ialpha,izeta] += self.d_R_d_s_mnc[isurf]*weight*cos_angle
              
            # Handle d Z / d s.  Since Z is on the full mesh, its radial derivative is on the half mesh.
            self.d_Z_d_s[ialpha,izeta] += self.d_Z_d_s_mns[isurf]*weight*sin_angle
              
            # Handle d Lambda / d s. Since Lambda is on the half mesh, its radial derivative is on the full mesh.
            self.d_Lambda_d_s[ialpha,izeta] += self.d_Lambda_d_s_mns[isurff]*weightf*sin_angle
        
        # Now consider the stellarator-asymmetric terms
        #add_stellaratorAsymmetricTerms()
        return

    #*********************************************************************
    # Using R(theta,zeta) and Z(theta,zeta), compute the Cartesian
    # components of the gradient basis vectors using the dual relations:
    #*********************************************************************
    def calculate_cartesianComponents(self):
        
        # Stella routine
        if self.verbose: print("stella::vmec_to_stella_geometry_interface-->Use the dual relations to get")
        
        # Construct the coordinate derivatives
        for izeta in range(self.nzeta):
         
            # Get the cosinus and sinus of the angles
            cos_angle = np.cos(self.zeta[izeta])
            sin_angle = np.sin(self.zeta[izeta])
             
            # X = R*cos(zeta)
            self.d_X_d_theta_vmec[:,izeta] = self.d_R_d_theta_vmec[:,izeta]*cos_angle
            self.d_X_d_zeta[:,izeta] = self.d_R_d_zeta[:,izeta]*cos_angle - self.R[:,izeta]*sin_angle
            self.d_X_d_s[:,izeta] = self.d_R_d_s[:,izeta]*cos_angle
             
            # Y = R*sin(zeta)
            self.d_Y_d_theta_vmec[:,izeta] = self.d_R_d_theta_vmec[:,izeta]*sin_angle
            self.d_Y_d_zeta[:,izeta] = self.d_R_d_zeta[:,izeta]*sin_angle + self.R[:,izeta]*cos_angle
            self.d_Y_d_s[:,izeta] = self.d_R_d_s[:,izeta]*sin_angle
        
        # Use the dual relations to get the Cartesian components of grad s 
        self.grad_s_X = (self.d_Y_d_theta_vmec*self.d_Z_d_zeta - self.d_Z_d_theta_vmec*self.d_Y_d_zeta) / self.sqrt_g
        self.grad_s_Y = (self.d_Z_d_theta_vmec*self.d_X_d_zeta - self.d_X_d_theta_vmec*self.d_Z_d_zeta) / self.sqrt_g
        self.grad_s_Z = (self.d_X_d_theta_vmec*self.d_Y_d_zeta - self.d_Y_d_theta_vmec*self.d_X_d_zeta) / self.sqrt_g
      
        # Use the dual relations to get the Cartesian components of grad theta_vmec 
        self.grad_theta_vmec_X = (self.d_Y_d_zeta*self.d_Z_d_s - self.d_Z_d_zeta*self.d_Y_d_s) / self.sqrt_g
        self.grad_theta_vmec_Y = (self.d_Z_d_zeta*self.d_X_d_s - self.d_X_d_zeta*self.d_Z_d_s) / self.sqrt_g
        self.grad_theta_vmec_Z = (self.d_X_d_zeta*self.d_Y_d_s - self.d_Y_d_zeta*self.d_X_d_s) / self.sqrt_g
      
        # Use the dual relations to get the Cartesian components of grad zeta:
        self.grad_zeta_X = (self.d_Y_d_s*self.d_Z_d_theta_vmec - self.d_Z_d_s*self.d_Y_d_theta_vmec) / self.sqrt_g
        self.grad_zeta_Y = (self.d_Z_d_s*self.d_X_d_theta_vmec - self.d_X_d_s*self.d_Z_d_theta_vmec) / self.sqrt_g
        self.grad_zeta_Z = (self.d_X_d_s*self.d_Y_d_theta_vmec - self.d_Y_d_s*self.d_X_d_theta_vmec) / self.sqrt_g 
      
        # Get grad_theta_pest = grad (theta_vmec + Lambda)
        self.grad_theta_pest_X = (1.0 + self.d_Lambda_d_theta_vmec)*self.grad_theta_vmec_X + self.d_Lambda_d_zeta*self.grad_zeta_X + self.d_Lambda_d_s*self.grad_s_X
        self.grad_theta_pest_Y = (1.0 + self.d_Lambda_d_theta_vmec)*self.grad_theta_vmec_Y + self.d_Lambda_d_zeta*self.grad_zeta_Y + self.d_Lambda_d_s*self.grad_s_Y
        self.grad_theta_pest_Z = (1.0 + self.d_Lambda_d_theta_vmec)*self.grad_theta_vmec_Z + self.d_Lambda_d_zeta*self.grad_zeta_Z + self.d_Lambda_d_s*self.grad_s_Z
        
        # The following gradient must be zero 
        self.grad_zeta_Z[:,:] = 0
        return
    
    #*********************************************************************
    # Compute the Cartesian components of other quantities we need:
    #*********************************************************************
    def calculate_moreCartestionComponents(self):

        # Stella routine
        if self.verbose: print("stella::vmec_to_stella_geometry_interface-->Compute the Cartesian components")
        
        # Calculate gradients of the magnetic field
        self.grad_B_X = self.d_B_d_s*self.grad_s_X + self.d_B_d_theta_vmec*self.grad_theta_vmec_X + self.d_B_d_zeta*self.grad_zeta_X
        self.grad_B_Y = self.d_B_d_s*self.grad_s_Y + self.d_B_d_theta_vmec*self.grad_theta_vmec_Y + self.d_B_d_zeta*self.grad_zeta_Y
        self.grad_B_Z = self.d_B_d_s*self.grad_s_Z + self.d_B_d_theta_vmec*self.grad_theta_vmec_Z + self.d_B_d_zeta*self.grad_zeta_Z
      
        # Calculate the components of the magnetic field
        self.B_X = self.edge_toroidal_flux_over_2pi*((1 + self.d_Lambda_d_theta_vmec)*self.d_X_d_zeta + (self.iota - self.d_Lambda_d_zeta)*self.d_X_d_theta_vmec) / self.sqrt_g
        self.B_Y = self.edge_toroidal_flux_over_2pi*((1 + self.d_Lambda_d_theta_vmec)*self.d_Y_d_zeta + (self.iota - self.d_Lambda_d_zeta)*self.d_Y_d_theta_vmec) / self.sqrt_g
        self.B_Z = self.edge_toroidal_flux_over_2pi*((1 + self.d_Lambda_d_theta_vmec)*self.d_Z_d_zeta + (self.iota - self.d_Lambda_d_zeta)*self.d_Z_d_theta_vmec) / self.sqrt_g
      
        # Radial coordinate
        self.sqrt_s = np.sqrt(self.normalized_toroidal_flux_used)
        
        # Calculate gradients
        self.grad_psi_X = self.grad_s_X*self.edge_toroidal_flux_over_2pi
        self.grad_psi_Y = self.grad_s_Y*self.edge_toroidal_flux_over_2pi
        self.grad_psi_Z = self.grad_s_Z*self.edge_toroidal_flux_over_2pi
      
        # Form grad alpha = grad (theta_vmec + Lambda - iota*zeta)
        for izeta in range(self.nzeta):
            self.grad_alpha_X[:,izeta] = (self.d_Lambda_d_s[:,izeta] - self.zeta[izeta]*self.d_iota_d_s)*self.grad_s_X[:,izeta]
            self.grad_alpha_Y[:,izeta] = (self.d_Lambda_d_s[:,izeta] - self.zeta[izeta]*self.d_iota_d_s)*self.grad_s_Y[:,izeta]
            self.grad_alpha_Z[:,izeta] = (self.d_Lambda_d_s[:,izeta] - self.zeta[izeta]*self.d_iota_d_s)*self.grad_s_Z[:,izeta]
     
        # Calculate gradients
        self.grad_alpha_X = self.grad_alpha_X + (1 + self.d_Lambda_d_theta_vmec)*self.grad_theta_vmec_X + (-self.iota + self.d_Lambda_d_zeta)*self.grad_zeta_X
        self.grad_alpha_Y = self.grad_alpha_Y + (1 + self.d_Lambda_d_theta_vmec)*self.grad_theta_vmec_Y + (-self.iota + self.d_Lambda_d_zeta)*self.grad_zeta_Y
        self.grad_alpha_Z = self.grad_alpha_Z + (1 + self.d_Lambda_d_theta_vmec)*self.grad_theta_vmec_Z + (-self.iota + self.d_Lambda_d_zeta)*self.grad_zeta_Z
            
    #*********************************************************************
    # For gbdrift, we need \vect{B} cross grad |B| dot grad alpha.
    # For cvdrift, we also need \vect{B} cross grad s dot grad alpha.
    # Let us compute both of these quantities 2 ways, and make sure the two
    # approaches give the same answer (within some tolerance).
    #*********************************************************************
    def calculate_crossProducts(self):
        
        # Stella routine
        if self.verbose: print("stella::vmec_to_stella_geometry_interface-->For gbdrift, we need \vect{B} cross grad |B| dot grad alpha.")
        
        # Calculate \vect{B} cross grad |B| dot grad alpha
        self.B_cross_grad_s_dot_grad_alpha = (self.B_sub_zeta*(1 + self.d_Lambda_d_theta_vmec) \
             - self.B_sub_theta_vmec*(self.d_Lambda_d_zeta - self.iota) ) / self.sqrt_g
      
        # Calculate \vect{B} cross grad s dot grad alpha
        self.bpol = self.B_sub_theta_vmec
        self.btor = self.B_sub_zeta
        self.brho = self.B_sub_s
        for izeta in range(self.nzeta):
            self.B_cross_grad_B_dot_grad_alpha[:,izeta] = 0 \
                + (self.B_sub_s[:,izeta]*self.d_B_d_theta_vmec[:,izeta]*(self.d_Lambda_d_zeta[:,izeta] - self.iota) \
                + self.B_sub_theta_vmec[:,izeta]*self.d_B_d_zeta[:,izeta]*(self.d_Lambda_d_s[:,izeta] - self.zeta[izeta]*self.d_iota_d_s) \
                + self.B_sub_zeta[:,izeta]*self.d_B_d_s[:,izeta]*(1 +self.d_Lambda_d_theta_vmec[:,izeta]) \
                - self.B_sub_zeta[:,izeta]*self.d_B_d_theta_vmec[:,izeta]*(self.d_Lambda_d_s[:,izeta] - self.zeta[izeta]*self.d_iota_d_s) \
                - self.B_sub_theta_vmec[:,izeta]*self.d_B_d_s[:,izeta]*(self.d_Lambda_d_zeta[:,izeta] - self.iota) \
                - self.B_sub_s[:,izeta]*self.d_B_d_zeta[:,izeta]*(1 + self.d_Lambda_d_theta_vmec[:,izeta])) / self.sqrt_g[:,izeta] 
        return

    #********************************************************************* 
    #        Finally, assemble the quantities needed for stella. 
    #*********************************************************************
    def calculate_quantitiesForStella(self): 
        
        # Stella routine
        if self.verbose: print("stella::vmec_to_stella_geometry_interface-->Finally, assemble the quantities")
        
        # Extract needed quantities
        B_reference = self.B_reference
        L_reference = self.L_reference
        sign_torflux = self.sign_torflux
        torflux = self.normalized_toroidal_flux_used 
        B = self.B  
        
        # Normalized magnetic field
        self.bmag = bmag = self.B/self.B_reference 
        
        # Parallel gradient
        self.gradpar_zeta = self.gradpar = L_reference*self.B_sup_zeta/B
  
        # grad alpha . grad alpha in units of 1/L_ref^2, with alpha = theta_pest - iota*zeta
        self.grad_alpha_grad_alpha = L_reference*L_reference*(self.grad_alpha_X*self.grad_alpha_X \
            + self.grad_alpha_Y*self.grad_alpha_Y + self.grad_alpha_Z*self.grad_alpha_Z)
        
        # grad alpha . grad psi_t in units of B_reference
        self.grad_alpha_grad_psi = (self.grad_alpha_X*self.grad_psi_X \
            + self.grad_alpha_Y*self.grad_psi_Y + self.grad_alpha_Z*self.grad_psi_Z) / B_reference
        
        # grad psi_t . grad psi_t in units of (B_reference*L_reference)^2
        self.grad_psi_grad_psi = (self.grad_psi_X*self.grad_psi_X + self.grad_psi_Y*self.grad_psi_Y \
            + self.grad_psi_Z*self.grad_psi_Z) / (L_reference*L_reference*B_reference*B_reference)
        
        # (grad zeta . grad x_stella) / bmag^2 = (grad zeta . grad psitor)*dx/dpsitor / bmag^2 = (grad zeta . grad psitor)*sign_torflux/rhotor/Lref/Bref/ bmag^2 
        self.gradzeta_gradx = sign_torflux*(self.grad_zeta_X*self.grad_psi_X + \
            self.grad_zeta_Y*self.grad_psi_Y + self.grad_zeta_Z*self.grad_psi_Z) \
            /(L_reference*B_reference*np.sqrt(torflux)*bmag**2) 
        
        # (grad zeta . grad y_stella) / bmag^2 = (grad zeta . grad alpha)*dy/dalpha / bmag^2 = (grad zeta . grad alpha)*Lref*rhotor / bmag^2
        self.gradzeta_grady = (self.grad_zeta_X*self.grad_alpha_X \
            + self.grad_zeta_Y*self.grad_alpha_Y + self.grad_zeta_Z*self.grad_alpha_Z) \
            *L_reference*np.sqrt(torflux) / bmag**2
        
        # (grad theta_pest . grad x_stella) / bmag^2
        self.gradtheta_gradx = (self.grad_theta_pest_X*self.grad_psi_X \
            + self.grad_theta_pest_Y*self.grad_psi_Y + self.grad_theta_pest_Z*self.grad_psi_Z) \
            /(L_reference*B_reference*np.sqrt(torflux)*bmag**2)
        
        # (grad theta_pest . grad y_stella) / bmag^2
        self.gradtheta_grady = (self.grad_theta_pest_X*self.grad_alpha_X \
            + self.grad_theta_pest_Y*self.grad_alpha_Y + self.grad_theta_pest_Z*self.grad_alpha_Z) \
            *L_reference*np.sqrt(torflux) / bmag**2
        
        # psitor/|psitor|*((grad y_stella . grad zeta)*(grad x_stella . grad y_stella) - (grad x_stella . grad zeta)*|grad y_stella|^2) / (B/Bref)^2
        self.gds23 = sign_torflux*(self.gradzeta_grady*sign_torflux*self.grad_alpha_grad_psi \
             - self.gradzeta_gradx*self.grad_alpha_grad_alpha*torflux) 
      
        # psitor/|psitor|*((grad y_stella . grad zeta)*|grad x_stella|^2 - (grad x_stella . grad zeta)*(grad x_stella . grad y_stella)) / (B/Bref)^2
        self.gds24 = (self.gradzeta_grady*self.grad_psi_grad_psi/torflux \
             - self.gradzeta_gradx*sign_torflux*self.grad_alpha_grad_psi)*sign_torflux 
      
        # ((grad y_stella . grad theta_pest)*(grad x_stella . grad y_stella) - (grad x_stella . grad theta_pest)*|grad y_stella|^2) / (B/Bref)^2
        self.gds25 = self.gradtheta_grady*sign_torflux*self.grad_alpha_grad_psi \
             - self.gradtheta_gradx*self.grad_alpha_grad_alpha*torflux 
      
        # ((grad y_stella . grad theta_pest)*|grad x_stella|^2  - (grad x_stella . grad theta_pest)*(grad x_stella . grad y_stella)) / 2 / (B/Bref)^2
        self.gds26 = 0.5*(self.gradtheta_grady*self.grad_psi_grad_psi/torflux \
             - self.gradtheta_gradx*sign_torflux*self.grad_alpha_grad_psi) 
      
        # Other quantities
        self.gbdrift_alpha = 2*B_reference*L_reference*L_reference*self.B_cross_grad_B_dot_grad_alpha / (B*B*B) 
        
        self.gbdrift0_psi = - self.edge_toroidal_flux_over_2pi \
            *(self.B_sub_theta_vmec*self.d_B_d_zeta - self.B_sub_zeta*self.d_B_d_theta_vmec) / self.sqrt_g \
            *2*self.shat / (B*B*B)
            
        self.cvdrift_alpha = self.gbdrift_alpha + 2*B_reference*L_reference*L_reference*self.mu_0*self.d_pressure_d_s \
            *self.B_cross_grad_s_dot_grad_alpha / (B*B*B*B)
            
        self.cvdrift0_psi = self.gbdrift0_psi
        if True: return 

#===============================================================================
#                               SANITY CHECKS
#===============================================================================

    def perform_sanityChecks(self, ns):
        ''' Sanity checks of the VMEC equilibrium. '''
        
        # Stella routine
        if self.verbose: print("stella::vmec_to_stella_geometry_interface-->! Do some sanity checking")
        
        # Do some validation
        if (self.nalpha<1): print("Error! nalpha must be >= 1. Instead it is",self.nalpha); exit
        if (self.nzgrid<1): print("Error! nzgrid must be >= 1. Instead it is",self.nzgrid); exit
        if (self.desired_normalized_toroidal_flux<=0): print("Error! desired_normalized_toroidal_flux must be > 0. Instead it is", self.desired_normalized_toroidal_flux); exit
        if (self.desired_normalized_toroidal_flux>1): print("Error! desired_normalized_toroidal_flux must be <= 1. Instead it is",self.desired_normalized_toroidal_flux); exit
            
        # Check that the vector <phi> starts with a zero
        if (abs(self.phi[0]) > 1e-14):
            print("Error# VMEC phi array does not begin with 0.")
            print("phi:",self.phi); exit 
        
        # Check whether <phi> is uniformly spaced
        dphi = self.phi[1] - self.phi[0]
        for i in range(2,ns):
            if (abs(self.phi[i]-self.phi[i-1]-dphi) > 1e-11):
                print("Error# VMEC phi array is not uniformly spaced.")
                print("phi:",self.phi); exit  
        
        # Check that <phip> containts the derived toroidal flux d(phi)/ds
        # <phip> is given on the half-mesh, so skip first point
        for i in range(1,ns):
            if (abs(self.phip[i]+self.phi[ns-1]/(2*np.pi)) > 1e-11):
                print("Error# VMEC phips array is not constant and equal to -phi(ns)/(2*pi).")
                print("phip(s):",self.phip); exit 
         
        # The first mode in the m and n arrays should be m=n=0:
        if (self.xm[0] != 0): print("First element of xm in the wout file should be 0."); exit
        if (self.xn[0] != 0): print("First element of xn in the wout file should be 0."); exit
        if (self.xm_nyq[0] != 0): print("First element of xm_nyq in the wout file should be 0."); exit
        if (self.xn_nyq[0] != 0): print("First element of xn_nyq in the wout file should be 0."); exit
         
        # Lambda should be on the half mesh, so its value at radial index 1 should be 0 for all (m,n)
        if (np.max(np.abs(self.lmns[:,0])) > 0):
            print("Error# Expected lmns to be on the half mesh, but its value at radial index 1 is nonzero.")
            print("Here comes lmns[:,0]:", self.lmns[:,0]); exit 
            
        # Check the asymmetric lambda as well
        if (self.lasym):
            if (np.max(abs(self.lmnc[:,0])) > 0):
                print("Error# Expected lmnc to be on the half mesh, but its value at radial index 1 is nonzero.")
                print("Here comes lmnc(:,1):", self.lmnc[:,0]); exit
        if True: return 
    
#===============================================================================
#                            INTERPOLATION WEIGHTS
#===============================================================================
    def calculate_interpolationWeights(self, ns):
        ''' In general, we get quantities for stella by linear interpolation, taking a weighted average of the quantity from
        2 surfaces in the VMEC file. Sometimes the weights are 0 and 1, i.e. no interpolation is needed.

        For any VMEC quantity Q on the full grid, the value used in stella will be
        Q_stella = Q(vmec_radial_index_full[0])*vmec_radial_weight_full[0] + Q(vmec_radial_index_full[1])*vmec_radial_weight_full[1]
        
        For any VMEC quantity Q on the half grid, the value used in stella will be
        Q_stella = Q[vmec_radial_index_half[0]]*vmec_radial_weight_half[0] + Q[vmec_radial_index_half[1]]*vmec_radial_weight_half[1]
        '''
        
        # Stella routine
        if self.verbose: print("stella::vmec_to_stella_geometry_interface-->! For any VMEC quantity Q on the full grid")

        # Initiate the weights
        self.radial_index_full  = np.empty(2, dtype='int')
        self.radial_weight_full = np.empty(2, dtype='float')
        self.radial_index_half  = np.empty(2, dtype='int')
        self.radial_weight_half = np.empty(2, dtype='float')
      
        # Handle quantities for the full grid
        if (self.normalized_toroidal_flux_used>1):
            if self.verbose: print("Error# normalized_toroidal_flux_used cannot be >1"); return
        elif (self.normalized_toroidal_flux_used<0):
            if self.verbose: print("Error# normalized_toroidal_flux_used cannot be <0"); return 
        elif (self.normalized_toroidal_flux_used==1):
            self.radial_index_full[0] = ns-1
            self.radial_index_full[1] = ns
            self.radial_weight_full[0] = 0
              
        # Normalized_toroidal_flux_used is >= 0 and <1, this is the most common case.
        else:
            self.radial_index_full[0] = np.floor(self.normalized_toroidal_flux_used*(ns-1))+1
            self.radial_index_full[1] = self.radial_index_full[0] + 1
            self.radial_weight_full[0] = self.radial_index_full[0] - self.normalized_toroidal_flux_used*(ns-1)
      
        # Derive the second component from the first one
        self.radial_weight_full[1] = 1 - self.radial_weight_full[0]
      
        # Handle quantities for the half grid
        if (self.normalized_toroidal_flux_used < self.normalized_toroidal_flux_half_grid[0]):
              
            # We start at element 2 since element 1 is always 0 for quantities on the half grid.
            if self.verbose: print("Warning: extrapolating beyond the end of VMEC's half grid.")
            if self.verbose: print("Extrapolating towards the magnetic axis.) Results are likely to be inaccurate.")
            self.radial_index_half[0] = 2
            self.radial_index_half[1] = 3
            self.radial_weight_half[0] = (self.normalized_toroidal_flux_half_grid[1] - self.normalized_toroidal_flux_used) / (self.normalized_toroidal_flux_half_grid[1] - self.normalized_toroidal_flux_half_grid[0])
      
        # We are past the VMECs grid
        elif (self.normalized_toroidal_flux_used > self.normalized_toroidal_flux_half_grid[ns-2]):
            if self.verbose: print("Warning: extrapolating beyond the end of VMEC's half grid.")
            if self.verbose: print("(Extrapolating towards the last closed flux surface.) Results may be inaccurate.")
            self.radial_index_half[0] = ns-1
            self.radial_index_half[1] = ns
            self.radial_weight_half[0] = (self.normalized_toroidal_flux_half_grid[ns-1] - self.normalized_toroidal_flux_used) \
                 / (self.normalized_toroidal_flux_half_grid[ns-1] - self.normalized_toroidal_flux_half_grid(ns-2))
      
        # We are exactly at the last point of the half grid
        elif (self.normalized_toroidal_flux_used == self.normalized_toroidal_flux_half_grid[ns-2]):
            self.radial_index_half[0] = ns-1
            self.radial_index_half[1] = ns
            self.radial_weight_half[0] = 0
              
        # Normalized_toroidal_flux_used is inside the half grid, this is the most common case.
        else:
            self.radial_index_half[0] = np.floor(self.normalized_toroidal_flux_used*(ns-1) + 0.5) + 1
            if (self.radial_index_half[0] < 2): self.radial_index_half[0] = 2 
            self.radial_index_half[1] = self.radial_index_half[0] + 1
            self.radial_weight_half[0] = self.radial_index_half[0] - self.normalized_toroidal_flux_used*(ns-1) - 0.5
      
        # Derive the second component from the first one
        self.radial_weight_half[1] = 1-self.radial_weight_half[0]
          
        # Correct the indexing from fortan to python
        self.radial_index_full = [ i-1 for i in self.radial_index_full]
        self.radial_index_half = [ i-1 for i in self.radial_index_half] 
        
        # Get the interpolation weights, for when svalue lies between two flux surfaces s
        self.surfaces_indexes = self.radial_index_half
        self.surfaces_weights = self.radial_weight_half
        if True: return 

#===============================================================================
#                               Extra functions
#===============================================================================

def get_stellaThetaAndZetaForVmec(path, alpha0=0, rho=0.7, poloidal_turns=1, nzed=96, nperiod=1, zeta_center=0):
    ''' Use the <vmec> class to get <zeta> and <theta>. '''
    
    # Packages
    import pathlib, os
    from stellapy.data.geometry.read_wout import read_allVariablesFromWoutFile
    from stellapy.data.paths.load_pathObject import create_dummyPathObject 
    from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_iota   
    
    # Read the vmec file to determine <nfield_periods> 
    stella = os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0]
    vmec_filename = pathlib.Path(path).name; svalue = rho*rho
    input_file = pathlib.Path(stella+"/stellapy.utils.config/input.in")
    path = create_dummyPathObject(input_file, vmec_filename, True)
    woutParameters = read_allVariablesFromWoutFile(path.vmec)     
    woutParameters['iota'] = calculate_iota(svalue=svalue, ns=woutParameters['ns'], iotas=woutParameters['iotas'])     
    nfield_periods = poloidal_turns*woutParameters['nfp']/woutParameters['iota'] 
    
    # Create the <vmec> class
    VMEC = vmec() 
    
    # Change the input parameters
    VMEC.alpha0 = alpha0 
    VMEC.nalpha = 1
    VMEC.nzgrid = int(nzed/2) + int((nperiod-1)*nzed)  
    VMEC.nzeta  = VMEC.nzgrid*2+1
    VMEC.woutParameters = woutParameters
    VMEC.nfp = nfp = woutParameters['nfp']
    VMEC.iota = woutParameters['iota']
    
    # Read the netcdf file
    dataset  = nc4.Dataset(path.vmec, mode='r')
    VMEC.normalized_toroidal_flux_used = rho*rho
    VMEC.iotas = dataset.variables['iotas'][:]
    VMEC.iotaf = dataset.variables['iotaf'][:]
    VMEC.presf = dataset.variables['presf'][:]
    VMEC.lasym = dataset.variables['lasym__logical__'][:]  
    VMEC.lmns = dataset.variables['lmns'][:]
    VMEC.xn = dataset.variables['xn'][:]
    VMEC.xm = dataset.variables['xm'][:] 
    VMEC.mnmax = int(dataset.variables['mnmax'][:]) 
    ns = int(dataset.variables['ns'][:])     
    
    # Initiate the arrays
    VMEC.initiate_arraysForVmec(VMEC.nalpha, VMEC.nzgrid, VMEC.nzeta)
    VMEC.alpha = [ VMEC.alpha0 + ((i-1)*2*np.pi)/VMEC.nalpha for i in np.arange(1, VMEC.nalpha+1, 1)]
    VMEC.zeta = [ zeta_center + (np.pi*i*nfield_periods)/(nfp*VMEC.nzgrid) for i in np.arange(-VMEC.nzgrid, VMEC.nzgrid+1, 1)]
    VMEC.calculate_normalizedToroidalFluxOnFullAndHalfGrids(ns)
    VMEC.calculate_interpolationWeights(ns)
    VMEC.calculate_radialProfileFunctionsAtFluxSurface(ns)
    VMEC.get_thetaGridForVmec()
    
    # Make sure we have arrays
    VMEC.theta_vmec = np.array(VMEC.theta_vmec[0,:])
    VMEC.theta_pest = np.array(VMEC.theta_pest[0,:])
    VMEC.zeta = np.array(VMEC.zeta) 
    return VMEC 

if __name__ == "__main__":
    get_stellaThetaAndZetaForVmec("/home/hanne/CIEMAT/RESEARCH/MagneticEquilibria/VmecFiles/wout_w7xr003.nc") 
