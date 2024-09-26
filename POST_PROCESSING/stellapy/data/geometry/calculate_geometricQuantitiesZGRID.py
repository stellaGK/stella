'''               ZGRID MODULE OF STELLA

zed_equal_arc
-------------
If zed_equal_arc = T, then zed is chosen to be arc length
If zed_equal_arc = F, then zed is poloidal (axisymmetric) or zeta (toroidal) angle

shat_zero
----------
Set minimum shat value below which we assume periodic BC

'''

import sys
import numpy as np

class zgrid:

    #*********************************************************************
    #                   Read default vmec parameters
    #*********************************************************************
    def __init__(self, verbose=True):
        
        # Stella routine 
        self.verbose = verbose
        if self.verbose: print("stella::zgrid::read_parameters")

        # Default input parameters
        self.nzed = 24
        self.nperiod = 1
        self.ntubes = 1
        self.boundary_option = "linked"
        self.boundary_option_switch = "boundary_option_linked"
        self.zed_equal_arc = False
        self.shat_zero = 1.0e-5
        self.grad_x_grad_y_zero = 1.0e-5 

    #*********************************************************************
    #               Replace defaults by input parameters
    #*********************************************************************
    def update_inputParameters(self, nzed=None, nperiod=None, boundary_option=None, zed_equal_arc=None, grad_x_grad_y_zero=None):

        # Stella routine
        if self.verbose: print("stella::zgrid::read_parameters")
        
        # Overwrite default parameters
        if grad_x_grad_y_zero!=None: self.grad_x_grad_y_zero = grad_x_grad_y_zero
        if boundary_option!=None: self.boundary_option = boundary_option
        if zed_equal_arc!=None:  self.zed_equal_arc = zed_equal_arc
        if nperiod!=None:        self.nperiod = nperiod
        if nzed!=None:           self.nzed = nzed
        
        # Other parameters
        if self.boundary_option=="zero":            self.boundary_option_switch = "boundary_option_zero"
        if self.boundary_option=="default":         self.boundary_option_switch = "boundary_option_zero"
        if self.boundary_option=="unconnected":     self.boundary_option_switch = "boundary_option_zero"
        if self.boundary_option=="periodic":        self.boundary_option_switch = "boundary_option_self_periodic"
        if self.boundary_option=="self-periodic":   self.boundary_option_switch = "boundary_option_self_periodic"
        if self.boundary_option=="linked":          self.boundary_option_switch = "boundary_option_linked"
        if self.boundary_option=="stellarator":     self.boundary_option_switch = "boundary_option_linked_stellarator"
        
        # Calculations 
        self.nzgrid = int(self.nzed/2) + int((self.nperiod-1)*self.nzed)  
        self.nzeta  = self.nzgrid*2+1
        self.zed =[ i*np.pi/(nzed/2) for i in np.arange(-self.nzgrid, self.nzgrid+1,1) ]

        # Get the step size <delzed> at every point of the z-axis 
        self.delzed = np.zeros((int(2*self.nzgrid+1))) 
        for i in range(2*self.nzgrid):
            self.delzed[i] = self.zed[i+1]-self.zed[i] 
        self.delzed[-1] = self.delzed[0]
        
        # Synonyms
        self.nzgrid_stella = self.nzgrid

    #*********************************************************************
    #                      Get <nzgrid_vmec>
    #*********************************************************************
    def get_zgrid_vmec(self, VMEC):
        if self.zed_equal_arc==True:  self.nzgrid_vmec = self.get_modified_vmec_zeta_grid(VMEC)
        if self.zed_equal_arc==False: self.nzgrid_vmec = self.nzgrid 
        
    #*********************************************************************
    #                      Get modified vmec zeta grid
    #*********************************************************************
    def get_modified_vmec_zeta_grid(self, VMEC):
        ''' We need to extend the maximum and minimum zeta values by <zgrid_scalefac>
        to ensure that we have information about geometric coefficients everywhere  
        on a fixed equal-arc grid in zed. First figure out how many extra zeta grid  
        points are required at the nominal grid spacing to get out to the ends of 
        the extended zeta domain. '''
        
        # Stella routine
        if self.verbose: print('stella::vmec_geo::get_modified_vmec_zeta_grid')
     
        # First calculate the nominal zeta grid used for vmec. Note that <nfield_periods>
        # is the number of field periods sampled in stella, while <nfield_periods_device>
        # is the number of field periods in the device. <nfield_periods> may be reasonably 
        # bigger than <nfield_periods_device> as the former is sampled while keeping 
        # alpha fixed (rather than theta)
        zeta = self.get_nominal_vmec_zeta_grid(VMEC)
     
        # Maximum zeta value for nominal zeta grid
        zeta_max = zeta[-1]
         
        # Excess_zeta is difference between expanded zeta_max and nominal zeta_max
        excess_zeta = zeta_max*(VMEC.zgrid_scalefac-1.0)
         
        # Assume equal grid spacing in zeta
        dzeta = zeta[1]-zeta[0]
        tmp = excess_zeta/dzeta
         
        # nzgrid_excess is the number of additional zeta grid points needed to cover at least excess_zeta
        if (abs(tmp - np.round(tmp)) < 0.1):
            nzgrid_excess = np.round(tmp)
        else:
            nzgrid_excess = np.round(tmp) + 1 
     
        # Now refine the zeta grid by desired amount in preparation for interpolation
        nzgrid_modified = (self.nzgrid + nzgrid_excess)*VMEC.zgrid_refinement_factor
        dzeta_modified = dzeta/VMEC.zgrid_refinement_factor
        self.dzeta_vmec = dzeta_modified
        return nzgrid_modified 


    #*********************************************************************
    #                    Get norminal vmec zeta grid
    #*********************************************************************
    def get_nominal_vmec_zeta_grid (self, VMEC):
        ''' 2*nzgrid+1 is the number of zeta grid points for the nominal zeta grid.
        The zeta domain is centered at zeta_center. Setting zeta_center = 2*pi*N/nfp for any integer N should
        yield identical results to setting zeta_center = 0, where nfp is the number of field periods (as in VMEC).
        <number_of_field_periods_stella> is the number of field periods sampled by stella.
        <number_of_field_periods_device> is the number of field periods for the device. 
        On exit, zeta holds the nominal grid points in the toroidal angle zeta. ''' 
        
        # Get the variables
        nzgrid = self.nzgrid
        nfp_device = VMEC.number_of_field_periods_device
        nfp_stella = VMEC.number_of_field_periods_stella
        if (VMEC.nfp_stella < 0.0): nfp_stella = nfp_device
        
        # Calculate the zeta grid used by stella
        zeta = [ VMEC.zeta_center + (np.pi*i*nfp_stella)/(nfp_device*nzgrid) for i in np.arange(-nzgrid, nzgrid+1, 1)]
        return zeta

    #********************************************************************* 
    #    Change variables from <vmec_to_stella_geometry> to <vmec_geo>
    #*********************************************************************
    def updateVariableNamesForVmecGeo(self, VMEC, nalpha, nzeta):
        # Vmec reads <zgrid_vmec> points, so remember that this is a vmec quantity
        self.zeta_vmec = VMEC.zeta 
        # Stella only needs <zgrid> points so prepare for the interpolation
        self.zeta = np.zeros((nalpha, nzeta)) 
        
        
    #********************************************************************* 
    #                      Initiate extended zgrid
    #********************************************************************* 
    def init_extended_zgrid(self, KGRID): 

        if self.verbose: print("stella::extended_zgrid::init_extended_zgrid")
       
        self.ntg = self.nzed/2
        self.neigen = np.zeros((KGRID.naky), dtype=int)
        self.periodic = np.zeros((KGRID.naky))     
         
        if (self.boundary_option_switch=="boundary_option_self_periodic"):
            self.periodic = True
        else:
            self.periodic = [True if np.abs(KGRID.aky[i]) < sys.float_info.epsilon else False for i in range(len(KGRID.aky))]  
          
        if (self.boundary_option_switch=="boundary_option_linked" or self.boundary_option_switch=="boundary_option_linked_stellarator"):
         
            # if linked BC, then iky=1 corresponds to ky=0 which has no connections
            self.neigen[0] = KGRID.nakx
            if (KGRID.naky > 1):
                for iky in np.arange(1,KGRID.naky):
                    # must link different kx values at theta = +/- pi
                    # neigen is the number of independent eigenfunctions along the field line
                    self.neigen[iky] = np.min([int((iky+1-1)*KGRID.jtwist),KGRID.nakx])
            
            self.neigen_max = np.max(self.neigen)
            
            self.ikx_shift_end = np.zeros((self.neigen_max)) 
            self.ikx_shift = np.zeros((KGRID.nakx, KGRID.naky))  
            
            # phi(kx-kx_shift,-nzgrid) = phi(kx,nzgrid) from twist-and-shift BC
            # for positive (negative) magnetic shear, kx_shift is positive (negative),
            # so start at most positive (negative) kx and
            # progress to smaller (larger) kx values as connections occur
            
            # figure out how much to shift ikx by to get to the end of the kx chain
            # for positive (negative) magnetic shear, this is the left-most (right-most) theta-theta0 
            # in each set of connected 2pi segments
            # note that theta0 goes from 0 to theta0_max and then from theta0_min back
            # to -dtheta0
            
            for ikx in range(self.neigen_max):
                # first ikx_max=nakx/2+1 theta0s are 0 and all positive theta0 values
                # remainder are negative theta0s
                # theta_0 = kx / ky / shat
                # if ky > 0, then most positive theta_0 corresponds to most positive kx
                 
                # first consider case where shift in kx is negative (corresponds to positive magnetic shear)
                if (KGRID.ikx_twist_shift < 0):
                    if (ikx <= KGRID.ikx_max):
                        self.ikx_shift_end[ikx] = KGRID.ikx_max-2*(ikx+1)+1
                    else:
                        self.ikx_shift_end[ikx] = 3*KGRID.ikx_max-2*(ikx+1)
                        
                # then consider case where shift in kx is positive
                elif (KGRID.ikx_twist_shift > 0):
                    if (ikx < KGRID.ikx_max):
                        if (ikx + KGRID.ikx_max <= KGRID.nakx):
                            self.ikx_shift_end[ikx] = KGRID.ikx_max
                        else:
                            self.ikx_shift_end[ikx] = ikx - KGRID.nakx 
                    else:
                        self.ikx_shift_end[ikx] = 1 - KGRID.ikx_max 
 
            
            for iky in range(KGRID.naky):
                # ikx_shift is how much to shift each ikx by to connect
                # to the next theta0 (from most positive to most negative for positive magnetic shear
                # and vice versa for negative magnetic shear)
                
                # first consider shift in index for case where shift is negative
                # (corresponds to positive magnetic shear)
                if (KGRID.ikx_twist_shift < 0):
                    # if ky > 0, then going to more negative theta0
                    # corresponds to going to more negative kx
                    for ikx in range(KGRID.ikx_max):
                        # if theta0 is sufficiently positive, shifting to more
                        # negative theta0 corresponds to decreasing ikx
                        if (ikx-self.neigen[iky] > 0):
                            self.ikx_shift[ikx,iky] = -self.neigen[iky]
                            # if a positive theta0 connects to a negative theta0
                            # must do more complicated mapping of ikx
                        elif (ikx-self.neigen[iky]+KGRID.nakx >= KGRID.ikx_max+1):
                            self.ikx_shift[ikx,iky] = KGRID.nakx - self.neigen[iky]
                            
                    # if theta0 is negative, then shifting to more negative
                    # theta0 corresponds to decreasing ikx
                    for ikx in np.arange(KGRID.ikx_max, KGRID.nakx): 
                        # if theta0 is sufficiently negative, it has no
                        # more negative theta0 with which it can connect
                        if (ikx-self.neigen[iky] > KGRID.ikx_max):
                            self.ikx_shift[ikx,iky] = -self.neigen[iky]
                            
                elif (KGRID.ikx_twist_shift > 0):
                    # if ky > 0, then going to more positive theta0
                    # corresponds to going to more positive kx
                    for ikx in range(KGRID.ikx_max):
                        # if shift in kx, kx_shift, is less than kx-kx_max,
                        # then shift by the appropriate amount
                        if (ikx+self.neigen[iky] <= KGRID.ikx_max):
                            self.ikx_shift[ikx,iky] = self.neigen[iky] 
                        
                    for ikx in np.arange(KGRID.ikx_max, KGRID.nakx): 
                        # if kx+kx_shift < 0, then simple shift by neigen
                        if (ikx+self.neigen[iky] <= KGRID.nakx):
                            self.ikx_shift[ikx,iky] = self.neigen[iky]
                        # if 0 < kx+kx_shift <= kx_max, then more complicated shift
                        # to positive set of kx values
                        elif (ikx-KGRID.ikx_max+self.neigen[iky] <= KGRID.nakx):
                            self.ikx_shift[ikx,iky] = self.neigen[iky] - KGRID.nakx 
            
            self.nsegments = np.zeros((self.neigen_max,KGRID.naky), dtype=int)
            
            for iky in range(KGRID.naky):
                if (self.neigen[iky] == 0):
                    self.nsegments[:,iky] = 1
                else:
                    self.nsegments[:,iky] = int((KGRID.nakx-1)/self.neigen[iky])
            
                    for ie in range(int((KGRID.nakx-1)%(self.neigen[iky]+1))):
                        self.nsegments[ie,iky] = int(self.nsegments[ie,iky] + 1)
             
            self.nseg_max = int(np.max(self.nsegments))
            self.iz_low = np.ones((self.nseg_max))*(-self.nzgrid)
            self.iz_mid = np.ones((self.nseg_max))*(0.0)
            self.iz_up = np.ones((self.nseg_max))*(self.nzgrid)
            
        if (self.boundary_option_switch=="default"):       
            
            self.neigen[:] = KGRID.nakx 
            self.neigen_max = KGRID.nakx
            self.ikx_shift_end = np.zeros((self.neigen_max))
            self.ikx_shift = np.zeros((KGRID.nakx,KGRID.naky)) 
            self.nsegments = np.zeros((self.neigen_max,KGRID.naky), dtype=int)
            
            # this is the number of 2pi poloidal segments in the extended theta domain,
            # which is needed in initializing the reponse matrix and doing the implicit sweep
            self.nsegments[:,:] = int(2*(self.nperiod-1) + 1)
            self.nseg_max = np.max(self.nsegments) 
            self.iz_low = np.ones((self.nseg_max))
            self.iz_mid = np.ones((self.nseg_max))
            self.iz_up = np.ones((self.nseg_max)) 
            
            # iz_low(j) is the ig index corresponding to the inboard midplane from below (theta=-pi) within the jth segment
            # iz_mid(j) is the ig index corresponding to the outboard midplane (theta=0) within the jth segment
            for iseg in range(self.nseg_max):
                self.iz_low[iseg] = -self.nzgrid + (iseg+1-1)*self.nzed
                self.iz_mid[iseg] = self.iz_low[iseg] + int(self.nzed/2)
                self.iz_up[iseg] = self.iz_low[iseg] + self.nzed 
         
        # initialize ikxmod to nakx
        # should not be necessary but just in case one tries to access
        # a value beyond nsegments(ie,iky)
        self.ikxmod = np.ones((self.nseg_max,self.neigen_max,KGRID.naky), dtype=int)*(KGRID.nakx) 
         
        for iky in range(KGRID.naky):
            # only do the following once for each independent set of theta0s
            # the assumption here is that all kx are on processor and sequential
            for ie in range(self.neigen[iky]):
                # remap to start at theta0 = theta0_max (theta0_min) for negative (positive) kx shift
                # for this set of connected theta0s   
                self.ikxmod[0,ie,iky] = (ie+1) + self.ikx_shift_end[ie] 
                if (self.nsegments[ie,iky] > 1): 
                    for iseg in range(1, self.nsegments[ie,iky], 1):   
                        self.ikxmod[iseg,ie,iky] = self.ikxmod[iseg-1,ie,iky] + self.ikx_shift[self.ikxmod[iseg-1,ie,iky],iky]
         
        del self.ikx_shift_end
        del self.ikx_shift
         
        self.it_left = np.zeros((self.ntubes))
        self.it_right = np.zeros((self.ntubes)) 
         
        self.it_right[self.ntubes-1] = 1 
        if (self.ntubes > 1):
            for it in range(self.ntubes-1):
                self.it_right[it] = it+1+1 
         
        self.it_left[0] = self.ntubes
        if (self.ntubes > 1):
            for it in np.arange(1, self.ntubes):
                self.it_left[it] = it+1-1 
         
        # this is the number of unique zed values in all segments but the first
        # the first has one extra unique zed value (all others have one grid common
        # with the previous segment due to periodicity)
        self.nzed_segment = self.iz_up[0]-self.iz_low[0]   
        return
 
 
 
 
