 
import sys
import numpy as np


class kgrid:

    #*********************************************************************
    #                   Read default vmec parameters
    #*********************************************************************
    def __init__(self, kt_grids_knobs="range", verbose=True):

        # From the physics flags (not relevant yet for this python module)
        # The first one turns on the full flux surface, the second one the radial variation
        self.full_flux_surface = False
        self.q_as_x = False 
        self.rhostar = -1
        self.runtype_option_switch = "runtype_standalone"
        self.job = 1
        self.kt_grids_knobs = kt_grids_knobs
        self.verbose = verbose

        # Default input parameters for nonlinear simulations
        if self.kt_grids_knobs=="box":
            
            # Stella routine
            if self.verbose: print("stella::kt_grids::read_kt_grids_box")
            
            # Default input parameters
            self.nx = 1
            self.ny = 1
            self.jtwist = -1
            self.jtwistfac = 1.0
            self.y0 = -1.0
            self.nalpha = 1
            self.centered_in_rho = True
            self.periodic_variation = False
            self.reality = True
            
            # Get the number of de-aliased modes in x and y. The number of  
            # de-aliased modes (nak) is 2/3 of the total number of modes (n)
            # This makes the nakx/naky modes free of aliasing when we put them in a
            # vector of length nx/ny where we add zeros to the left and right
            # Along naky we only take half the modes since we consider only ky>0
            self.naky = int((self.ny-1)/3) + 1
            self.nakx = int(2*(int((self.nx-1)/3))) + 1
            
            # Full flux surface
            if (self.full_flux_surface): self.nalpha = self.ny
            
        # Default input parameters for linear simulations
        if self.kt_grids_knobs=="range":
            
            # Stella routine
            if self.verbose: print("stella::kt_grids::read_kt_grids_range")
            
            # Default input parameters
            self.nalpha = 1
            self.naky = 1
            self.nakx = 1
            self.aky_min = 0.0
            self.aky_max = 0.0 
            self.akx_min = 0.0
            self.akx_max = -1.0
            self.theta0_min = 0.0
            self.theta0_max = -1.0    
            
    #*********************************************************************
    #               Replace defaults by input parameters
    #*********************************************************************
    def update_inputParameters(self, nalpha, kx, ky, nx, ny, y0):
        
        # Overwrite default parameters
        if nalpha!=None:   self.nalpha = nalpha
        if kx!=None:       self.akx_min = kx; self.akx_max = kx
        if ky!=None:       self.aky_min = ky; self.aky_max = ky
        if nx!=None:       self.nx = nx 
        if ny!=None:       self.ny = ny
        if y0!=None:       self.y0 = y0
        if nx!=None:       self.naky = int((self.ny-1)/3) + 1
        if ny!=None:       self.nakx = int(2*(int((self.nx-1)/3))) + 1

    #*********************************************************************
    #                    Initiate the kt grids range
    #*********************************************************************
    def init_kt_grids_range(self, STELLA, ZGRID): 
    
        # set default akx and theta0 to 0
        self.akx = np.zeros((self.nakx))
        self.aky = np.zeros((self.nakx))
        self.theta0 = np.zeros((self.naky,self.nakx))
        self.zed0 = np.zeros((self.naky,self.nakx))
        self.ikx_twist_shift = 0
        
        # NB: we are assuming here that all ky are positive when running in range mode
        self.dky = 0.0
        if (self.naky > 1): self.dky = (self.aky_max - self.aky_min)/(self.naky - 1)
        self.aky = [ self.aky_min + self.dky*i for i in range(self.naky) ]
        
        if (self.q_as_x):
            self.tfac = 1.0
        else:
            self.tfac = STELLA.shat 
        
        zero = 100.*sys.float_info.epsilon 
        
        # if theta0_min and theta0_max have been specified,
        # use them to determine akx_min and akx_max
        if (self.theta0_max > self.theta0_min-zero):
            if (STELLA.shat > sys.float_info.epsilon):
                self.akx_min = self.theta0_min * self.tfac * self.aky[0]
                self.akx_max = self.theta0_max * self.tfac * self.aky[0]
            else:
                self.akx_min = self.theta0_max * self.tfac * self.aky[0]
                self.akx_max = self.theta0_min * self.tfac * self.aky[0] 
        
        # shat_zero is minimum shat value below which periodic BC is enforced
        if (abs(STELLA.shat) > ZGRID.shat_zero):  # ie assumes boundary_option .eq. 'linked'
            # if akx_min and akx_max specified in input
            # instead of theta0_min and theta0_max,
            # use them to get theta0_min and theta0_max 
            if (self.theta0_min > self.theta0_max+zero and abs(self.aky[0]) > zero):
                self.theta0_min = self.akx_min/(self.tfac*self.aky[0])
                self.theta0_max = self.akx_max/(self.tfac*self.aky[0])
                self.dtheta0 = 0.0
                 
                if (self.nakx > 1): self.dtheta0 = (self.theta0_max - self.theta0_min)/(self.nakx - 1)
                   
                for j in range(self.naky):
                    self.theta0[j,:] = [ self.theta0_min + self.dtheta0*i for i in range(self.nakx) ] 
                self.akx = self.theta0[0,:] * self.tfac * self.aky[0]
                
            elif (self.akx_max > self.akx_min-zero):
                self.dkx = 0.0
                if (self.nakx > 1): self.dkx = (self.akx_max - self.akx_min)/(self.nakx - 1)
                self.akx = [ self.akx_min + self.dkx*i for i in range(self.nakx) ]
            
                self.dtheta0 = 0.0
                if (self.nakx > 1): self.dtheta0 = (self.theta0_max - self.theta0_min)/(self.nakx - 1)
            
                if (STELLA.shat > sys.float_info.epsilon):
                    for j in range(self.naky):
                        self.theta0[j,:] = [ self.theta0_min + self.dtheta0*i for i in range(self.nakx) ] 
                else:
                    for j in range(self.naky):   
                        self.theta0[j,:] = [ self.theta0_min + self.dtheta0*i for i in np.arange(self.nakx-1,0-1,-1) ] 
            else:
                if self.verbose: print('ky=0 is inconsistent with kx_min different from kx_max. aborting.'); exit 
              
        else:
            # here assume boundary_option .eq. 'periodic'
            # used for periodic finite kx ballooning space runs with shat=0
            self.dkx = 0.0
            if (self.nakx > 1): self.dkx = (self.akx_max - self.akx_min)/(self.nakx - 1)
            self.akx = [ self.akx_min + self.dkx*i for i in range(self.nakx) ] 
        
        self.ikx_max = self.nakx
        self.naky_all = self.naky
        
        # Determine if iky corresponds to zonal mode
        self.zonal_mode = np.zeros((self.naky))  
        if (abs(self.aky[0]) < sys.float_info.epsilon): self.zonal_mode[0] = True
        return
 

    #*********************************************************************
    #                    Initiate the kt grids box
    #*********************************************************************
    def init_kt_grids_box(self, STELLA, ZGRID):
        
        # Stella routine
        if self.verbose: print("stella::kt_grids::init_kt_grids_box")

        # set jtwist and y0 for cases where they have not been specified
        # and for which it makes sense to set them automatically
        if (self.jtwist < 1): 
            self.jtwist = int(np.max([1,int(np.abs(STELLA.twist_and_shift_geo_fac)+0.5)])) 
            self.jtwist = int(np.max([1,int(self.jtwistfac*self.jtwist+0.5)]))

        # signed version of jtwist, with sign determined by, e.g., magnetic shear
        self.ikx_twist_shift = int(-self.jtwist*int(np.sign(STELLA.twist_and_shift_geo_fac)))
        
        # y0 defines the box size along y
        if (self.y0 < 0.):
            if (self.full_flux_surface):
                # if simulating a flux annulus, then y0 determined by the physical extent of the device
                if (self.rhostar > 0.):
                    self.y0 = 1./(self.rhostar*STELLA.rhotor)
                elif self.verbose: print('must set rhostar if simulating a full flux surface. aborting.'); exit 
            else:
                # if simulating a flux tube, then it makes no sense to have y0 < 0.0 so abort
                if self.verbose: print('y0 negative only makes sense when simulating a flux annulus.  aborting.'); exit
 
        # Get the grid spacing in ky and then in kx using twist-and-shift BC
        self.dky = 1./self.y0
        
        # non-quantized b/c assumed to be periodic instead 
        # of linked boundary conditions if zero magnetic shear
        if (abs(STELLA.shat) <= ZGRID.shat_zero):
            self.dkx = self.dky / self.jtwist
        else:
            self.dkx = self.dky * np.abs(STELLA.twist_and_shift_geo_fac) / self.jtwist
      
        # Define the real space box size
        self.x0 = 1./self.dkx
        
        # ky goes from zero to ky_max
        self.aky = [ round(iky*self.dky,2) for iky in np.arange(0, self.naky) ]  
        
        # get the ikx index corresponding to kx_max
        self.ikx_max = int(self.nakx/2)+1
        
        # get the total number of ky values, including negative ky
        self.naky_all = 2*self.naky-1
        
        # kx goes from zero to kx_max down to zero and then from -kx_max to -|kx_min|     
        self.akx = [ round(ikx*self.dkx,2) for ikx in np.arange(0, self.ikx_max) ]        
        self.akx = self.akx + [ round((ikx-self.nakx)*self.dkx,2) for ikx in np.arange(self.ikx_max, self.nakx) ] 
        
        # turn lists into arrays 
        self.akx = np.array(self.akx)
        self.aky = np.array(self.aky)
        
        # set theta0=0 for ky=0
        self.theta0 = np.zeros((self.naky, self.nakx)) 
         
        # theta0 = kx/ky
        if (self.q_as_x):
            for ikx in range(self.nakx):
                self.theta0[1:,ikx] = self.akx[ikx]/self.aky[1:]
                 
        # theta0 = kx/ky/shat
        elif (np.abs(STELLA.shat) > ZGRID.shat_zero):
            for ikx in range(self.nakx): 
                self.theta0[1:,ikx] = self.akx[ikx]/(self.aky[1:]*STELLA.shat) 
              
        # if shat=0, theta0 is meaningless, so be careful
        else: 
            for ikx in range(self.nakx): 
                self.theta0[1:,ikx] = - self.akx[ikx]/self.aky[1:] 
        
        # Allocate arrays for radial variation
        self.x = np.zeros((self.nx))
        self.x_d = np.zeros((self.nakx))
        self.rho = np.zeros((self.nx))
        self.rho_d = np.zeros((self.nakx)) 
         
        # Step size in real space
        self.dx = (2*np.pi*self.x0)/self.nx
        self.dy = (2*np.pi*self.y0)/self.ny
         
#         # Shift
#         self.x_shift = np.pi*self.x0
#         if (self.centered_in_rho):
#             if (self.q_as_x):
#                 self.dqdrho = STELLA.shat*STELLA.qinp/STELLA.rhoc
#                 self.x_shift = np.pi*self.x0*(1.0 - 0.5*self.rhostar*np.pi*self.x0*STELLA.d2qdr2/(self.dqdrho**2*STELLA.dxdXcoord))
#             else:
#                 self.x_shift = np.pi*self.x0*(1.0 - 0.5*self.rhostar*np.pi*self.x0*STELLA.d2psidr2*STELLA.drhodpsi**2/STELLA.dxdXcoord)
#          
#         # x vector
#         for ikx in range(self.nx):
#             if (self.runtype_option_switch=="runtype_multibox" and self.job==1):
#                 self.x[ikx] = (ikx+1-0.5)*self.dx - self.x_shift
#             elif (self.runtype_option_switch=="runtype_multibox" and self.radial_variation):
#                 if (self.periodic_variation):
#                     if (ikx < self.nx/2):
#                         self.x[ikx] = (ikx+1-0.5)*self.dx - 0.5*self.x_shift
#                     else:
#                         self.x[ikx] = self.x[self.nx-ikx+1] 
#                 else:
#                     self.x[ikx] = (ikx-0.5)*self.dx - self.x_shift 
#             else:
#                 self.x[ikx] = [ikx+1-1]*self.dx 
#          
#         # Something
#         self.dx_d = (2*np.pi*self.x0)/self.nakx
#         for ikx in range(self.nakx):
#             if (self.runtype_option_switch=="runtype_multibox" and self.job==1):
#                 self.x_d[ikx] = (ikx+1-0.5)*self.dx_d - self.x_shift
#             elif (self.runtype_option_switch=="runtype_multibox") and self.radial_variation:
#                 if (self.periodic_variation):
#                     if(ikx <= (self.nakx+1)/2):
#                         self.x_d[ikx] = (ikx+1-0.5)*self.dx_d - 0.5*self.x_shift
#                     else:
#                         self.x_d[ikx] = self.x_d[self.nakx-ikx+1] 
#                 else:
#                     self.x_d[ikx] = (ikx+1-0.5)*self.dx_d - self.x_shift 
#             else:
#                 self.x_d[ikx] = (ikx+1-1)*self.dx_d 
#         
#         # x to rho
#         self.rho = self.get_x_to_rho(1,self.x,self.rho,STELLA)
#         self.rho_d = self.get_x_to_rho(1,self.x_d,self.rho_d,STELLA)
#           
#         # Rho clamped
#         self.rho_clamped = self.rho
#         self.rho_d_clamped = self.rho_d
#           
#         # Zed0
#         self.zed0 = self.theta0*STELLA.zed0_fac
#           
#         # Radial varaiation  
#         if (self.radial_variation): self.dump_radial_grid()
#         if (self.radial_variation and (np.any((self.rho+STELLA.rhoc) < 0.0) or np.any((self.rho+STELLA.rhoc) > 1.0))):
#             if self.verbose: print('rho(x) is beyond range [0,1]. Try changing rhostar or q/psi profiles'); exit 
            
        # Determine if iky corresponds to zonal mode
        self.zonal_mode = np.zeros((self.naky)) 
        if (abs(self.aky[0]) < sys.float_info.epsilon): self.zonal_mode[0] = True
        return 

  
    def get_x_to_rho (self, llim, x_in, rho_out, STELLA):
        ''' Rho is the root of the quadratic equation in x. '''

        ulim = np.len(x_in)+llim-1 
    
        # The default option when radial_variation==true
        if (self.q_as_x):
            a = 0.5*STELLA.d2qdr2/self.dqdrho
            b = 1.0
            c = -self.rhostar/(self.dqdrho*STELLA.dxdXcoord)
    
        # The default option when radial_variation==false
        else:
            a = 0.5*STELLA.d2psidr2*STELLA.drhodpsi
            b = 1.0
            c = -STELLA.rhostar*STELLA.drhodpsi/STELLA.dxdXcoord 
    
        for ix in np.arange(llim-1, ulim):
            if (abs(4.0*a*c*x_in(ix)) < 1.e-6):
                rho_out[ix] = -(c*x_in[ix])/b - a*(c*x_in[ix])**2/b**3
            else:
                rho_out[ix] = (-b + np.sqrt(b**2 - 4.*a*c*x_in[ix]))/(2.*a) 
        return rho_out


    def dump_radial_grid(self):
        pass
        return