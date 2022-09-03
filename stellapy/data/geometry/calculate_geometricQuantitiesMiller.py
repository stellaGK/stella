
import numpy, sys
from stellapy.plot.utils.style.load_styleFigures import load_styleFigures
from stellapy.data.geometry.calculate_geometricQuantitiesSTELLA import geo_spline

#===============================================================================
#              CALCULATE THE MILLER COORDINATES LIKE IN STELLA
#===============================================================================

def calculate_geometricQuantitiesMiller(miller_variables, nzed, nperiod, nalpha):
    ''' Note that <geo_surf>, <local> and <local_out> are a class in stella 
    containing all the geometric quantities from the Miller equilibrium. 

    kperp2 = (n0/a)^2*(drho/dpsiN)^2*(gds2 + 2*theta0*gds21 + theta0^2*gds22)
    theta0 = kx/(ky*shat)
    
    Note that the definitions of gbdrift, gbdrift0, dgbdriftdr and dgbdrift0dr
    are such that it gets multiplied by vperp2, not mu.  This is in contrast to
    Michael's GS3 notes'''
    
    # Create the miller and stella class, and load the default parameters
    LOCAL = local()
    MILLER = miller(LOCAL)
    STELLA = stella(nzed, nperiod)
    
    # Read in the Miller local parameters
    MILLER.read_local_parameters(nzed, nperiod, miller_variables)

    # Allocate the geometry arrays for stella
    STELLA.allocate_arrays(nalpha, 2*MILLER.nzgrid+1)
    
    # Pass the stella variables onto the miller class (arguments of <get_local_geo>)
    MILLER.get_local_geo_arguments(STELLA)
    
    # Use Miller local parameters to get geometric coefficients needed by stella
    MILLER.get_local_geo()
    return MILLER

#===============================================================================
#                              BASIC FUNCTIONS
#===============================================================================
def calculate_Rpos(r, theta, j, self, LOCAL, nz):
    """ Here r = {
            self.LOCAL.rhoc - 1.e-3*(self.LOCAL.rhoc/self.LOCAL.rmaj),
            self.LOCAL.rhoc,
            self.LOCAL.rhoc + 1.e-3*(self.LOCAL.rhoc/self.LOCAL.rmaj)}
    """ 
    

    
    # Get dr = (r-r0) = {0, +- 1.e-3*(LOCAL.rhoc/LOCAL.rmaj)}
    dr = r - LOCAL.rhoc 

    # For Y Xiao: 
    #    g = LOCAL.delp/LOCAL.rhoc + LOCAL.d * sin(theta)**2
    #    Rpos = LOCAL.rmaj*(1.+r*(cos(theta)-g)-g*dr)

    # For others?
    g = numpy.cos(theta + LOCAL.tri * numpy.sin(theta))
    gp = -numpy.sin(theta + LOCAL.tri * numpy.sin(theta)) * LOCAL.triprim*numpy.sin(theta)

    # Allow for strange specification of R_psi
    if j==(2*nz+1): i = 0   # Fortran:  if (j==nz+1): i = -nz
    if j!=(2*nz+1): i = j 
    
    # Calculate R(theta) for dr[0], dr[1] or dr[2]
    Rpos = LOCAL.rmaj + LOCAL.shift*dr + g*LOCAL.rhoc + (g+LOCAL.rhoc*gp)*dr
    
    # If read_profile_variation, then add (1/2)*(r-r0)**2*d2R/dr|_r0, otherwise d2R=0
    Rpos += 0.5*(r-self.rhoc0)**2*self.d2R[i]
    return Rpos

#---------------------------
def calculate_Zpos(r, theta, j, self, LOCAL, nz):

    # Get dr = (r-r0) = {0, +- 1.e-3*(LOCAL.rhoc/LOCAL.rmaj)}
    dr = r - LOCAL.rhoc
    
    # Allow for strange specification of Z_psi
    if j==(2*nz+1): i = 0   # Fortran:  if (j==nz+1): i = -nz
    if j!=(2*nz+1): i = j 

    # Calculate Z(theta) for dr[0], dr[1] or dr[2]
    Zpos = LOCAL.kappa*numpy.sin(theta)*LOCAL.rhoc + (LOCAL.rhoc*LOCAL.kapprim + LOCAL.kappa)*numpy.sin(theta)*dr 
    
    # If read_profile_variation, then add (1/2)*(r-r0)**2*d2R/dr|_r0, otherwise d2R=0
    Zpos += 0.5*(r-self.rhoc0)**2*self.d2Z[i]
    return Zpos

#--------------------------- 
def get_drho(f, LOCAL):
    """ Takes in f(r), with r given at three radial locations and 
    returns df = df/dr at the middle radius. """ 
    return 0.5*(f[2,:]-f[0,:])/LOCAL.dr

#--------------------------- 
def get_dthet(f, df, delthet):
    """ Given function f(theta:-pi->pi), calculate theta derivative second order
    accurate, with equal grid spacing assumed. Assumes periodic in theta -- 
    may need to change this in future. """

    # Assuming equal grid spacing in theta here 
    df[1:-1] = (f[2:]-f[:-2])/(delthet[:-1]+delthet[1:])

    # Use periodicity at boundary
    df[0] = (f[1]-f[-2])/(delthet[0]+delthet[-1])
    df[-1] = df[0]
    return df

#--------------------------- 
def get_d2dthet2(f, d2f, delthet):
    """ Given function f(theta), calculate second derivative of f with respect to 
    theta second order accurate, with equal grid spacing assumed. """ 

    # Assuming equal grid spacing in theta here 
    d2f[1:-1] = (f[:-2]-2*f[1:-1]+f[2:])/delthet[1:]**2
    
    # Use periodicity at boundary
    d2f[0] = (f[-2]-2*f[0]+f[1])/delthet[1]**2
    d2f[-1] = d2f[0] 
    return d2f

#--------------------------- 
def theta_integrate(integrand, delthet, nz2pi, nz): 
    """ Use trapezoidal rule to integrate in theta. """
    imin = int(-nz2pi+nz); imax = int(nz2pi+nz+1) 
    integral = 0.5*numpy.sum(delthet[imin:imax-1]*(integrand[imin:imax-1] + integrand[imin+1:imax]))
    return integral

# Get indefinite integral of integrand
def theta_integrate_indef(integrand, integral, delthet, nz):
    """ Use trapezoidal rule to integrate in theta. """ 
    integral[0] = 0.0
    for i in range(nz+1,2*nz+1,1): integral[i] = integral[i-1]+0.5*delthet[i-1]*(integrand[i-1]+integrand[i])
    for i in range(nz-1, -1, -1):  integral[i] = integral[i+1]-0.5*delthet[i]*(integrand[i+1]+integrand[i]) 
    return integral 

#--------------------------- 
def get_dIdrho (dpsidrho, grho, nz, nz2pi, jacrho, bi, gpsi, Rr, dRdrho, dqdr, \
                dRdth, dZdth, d2Zdrdth, drzdth, d2Rdrdth, delthet, LOCAL):

    # Create a dummy vector
    dum = numpy.zeros((2*nz+1))
    imin = int(-nz2pi+nz); imax = int(nz2pi+nz+1) 

    # Calculate the dummy vector
    dum = jacrho*( 1.0 + (bi/gpsi)**2 ) / Rr[1,:]**2
    denom = theta_integrate(dum[imin:imax], delthet, nz2pi, nz)
    dum = jacrho*( 2.*dRdrho/Rr[1,:] + dqdr/LOCAL.qinp ) / Rr[1,:]**2
    num1 = theta_integrate(dum[imin:imax], delthet, nz2pi, nz)

    # betaprim below is (4*pi*ptot/B0^2)*(-d ln ptot / drho)
    dum = ( -2.*(dRdth*d2Rdrdth + dZdth*d2Zdrdth)/jacrho + drzdth + LOCAL.betaprim*jacrho/dpsidrho**2 ) / grho**2
    num2 = theta_integrate(dum[imin:imax], delthet, nz2pi, nz)

    # Calculate the current
    dIdrho = bi*(num1 + num2)/denom
    return dIdrho
    
#===============================================================================
#               SAVE SOME MILLER COORDINATES TO A <LOCAL> CLASS
#===============================================================================

class local: 
    def __init__(self):
        return

#===============================================================================
#               SAVE THE MILLER COORDINATES TO A <MILLER> CLASS
#===============================================================================

class miller:
    
    #*********************************************************************
    #                   Read default miller parameters
    #*********************************************************************
    def __init__(self, LOCAL):
        self.LOCAL = LOCAL
        self.init_local_defaults()
        return
        
    def init_local_defaults(self):
        
        # Default parameters
        self.nzed_local = 128                # Miller coordinates gets calculated with this nzed
        self.rhoc = 0.5
        self.rmaj = 3.0
        self.rgeo = 3.0
        self.qinp = 1.4                      # Safery factor (inverse of the rotational transform i)
        self.shat = 0.8
        self.shift = 0.0
        self.kappa = 0.0
        self.kapprim = 0.0
        self.tri = 0.0
        self.triprim = 0.0
        self.betaprim = 0.0                  # betaprim = -(4pi/Bref^2)*d(ptot)/drho
        self.betadbprim = 0.0                # betadbprim = -(4pi/Bref^2)*d^2ptot/drho^2
        self.d2qdr2 = 0.0
        self.d2psidr2 = 0.0
        
        # This variable can be changed in the profile variation file
        self.rhoc0 = 0.5                     # Radial location of the data given in the profile variation file 
        
        # Toggles
        self.read_profile_variation = False
        self.write_profile_variation = False
        self.load_psi0_variables = True
        
        # Only needed for sfincs when not using geo info from file
        self.rhotor = self.rhoc
        self.psitor_lcfs = 1.0
        self.drhotordrho = 1.0
        return

    #*********************************************************************
    #                Read miller parameters from input file
    #*********************************************************************
    def read_local_parameters(self, nzed, nperiod, miller_variables):
        
        # Remember nzed and nzgrid
        self.nzed = nzed
        self.nzgrid = int(nzed/2) + int((nperiod-1)*nzed)  
        self.nzeta  = self.nzgrid*2+1
        self.zed = [ i*numpy.pi/(nzed/2) for i in numpy.arange(-self.nzgrid, self.nzgrid+1,1) ]
        
        # The following variables can be modified
        keys = ["rhoc", "rmaj", "shift", "qinp", "shat", "kappa", "kapprim", "tri", \
                "triprim", "rgeo", "betaprim", "betadbprim", "d2qdr2", "d2psidr2", \
                "nzed_local", "read_profile_variation", "write_profile_variation"]
        
        # Overwrite the variables that were given
        for key in keys:
            if key in miller_variables:
                setattr(self, key, miller_variables[key])
                
        # This routine also sets <zed0_fac> in stella
        self.zed0_fac = 1.0
        
        # Derived input variables
        self.dr = 1.e-3*(self.rhoc/self.rmaj)

        # Set more variables for sfincs?
        self.dpsitordrho = 0.0
        self.d2psitordrho2 = 0.0
        
        # The next three variables are for multibox simulations with radial variation
        self.rhoc_psi0 = self.rhoc
        self.qinp_psi0 = self.qinp
        self.shat_psi0 = self.shat
        
        # Save a bunch of these variables to <local>
        keys += ["zed0_fac", "dr", "dpsitordrho", "d2psitordrho2", "rhoc_psi0", "qinp_psi0", "shat_psi0"]
        for key in keys: setattr(self.LOCAL, key, getattr(self, key))
        
        # First get nperiod corresponding to input number of grid points
        self.nz2pi = self.nzed/2
        self.np = (self.nzgrid-self.nz2pi)/self.nzed + 1
        
        # Now switch to using (possible higher resolution) local grid
        self.nz2pi = self.nzed_local/2
        
        # This is the equivalent of nzgrid on the local grid
        self.nz = int(self.nz2pi + self.nzed_local*(self.np-1))
        
        # Initialize arrays to zero, they will be overwritten if reading in from file
        # Only relevant for profile variation tests. These needs to be deallocated somewhere
        self.d2R = numpy.zeros((2*self.nz+1))          # d2R(-nz:nz)
        self.d2Z = numpy.zeros((2*self.nz+1))          # d2R(-nz:nz)
        self.bmag_psi0 = numpy.zeros((2*self.nz+1))    # bmag_psi0(-nz:nz)
        self.d2Rgrho_psi0= numpy.zeros((2*self.nz+1))  # grho_psi0(-nz:nz)
        self.dI = 0.
        return 
    
    
    #*********************************************************************
    #           Pass the stella arguments onto the miller routine
    #*********************************************************************
    def get_local_geo_arguments(self, STELLA):
        self.STELLA = STELLA
        self.zeta_out = STELLA.zeta
        self.nzed = STELLA.nzed
        self.nzgrid = STELLA.nzgrid
        self.zed_in = STELLA.zed
        self.nzed_in = STELLA.nzed
        self.dpsidrho_out = STELLA.dpsidrho
        self.dpsidrho_psi0_out = STELLA.dpsidrho_psi0
        self.dIdrho_out = STELLA.dIdrho
        self.grho_out = STELLA.grho[0,:]
        self.bmag_out = STELLA.bmag[0,:]
        self.bmag_psi0_out = STELLA.bmag_psi0[0,:]
        self.gds22 = STELLA.gds2[0,:]
        self.gds21_out = STELLA.gds21[0,:]
        self.gds22_out = STELLA.gds22[0,:]
        self.gds23_out = STELLA.gds23[0,:]
        self.gds24_out = STELLA.gds24[0,:]
        self.gradpar_out = STELLA.gradpar
        self.gbdrift0_out = STELLA.gbdrift0[0,:]
        self.gbdrift_out = STELLA.gbdrift[0,:]
        self.cvdrift0_out = STELLA.cvdrift0[0,:]
        self.cvdrift_out = STELLA.cvdrift[0,:]
        self.dBdrho_out = STELLA.dBdrho 
        self.d2Bdrdth_out = STELLA.d2Bdrdth 
        self.dgradpardrho_out = STELLA.dgradpardrho 
        self.btor_out = STELLA.btor 
        self.rmajor_out = STELLA.rmajor 
        self.dcvdrift0drho_out = STELLA.dcvdrift0drho[0,:]
        self.dcvdriftdrho_out = STELLA.dcvdriftdrho[0,:]
        self.dgbdrift0drho_out = STELLA.dgbdrift0drho[0,:]
        self.dgbdriftdrho_out = STELLA.dgbdriftdrho[0,:]
        self.dgds2dr_out = STELLA.dgds2dr[0,:]
        self.dgds21dr_out = STELLA.dgds21dr[0,:]
        self.dgds22dr_out = STELLA.dgds22dr[0,:]
        self.djacdrho_out = STELLA.djacdrho[0,:]
        return

    #*********************************************************************
    #                      Allocate arrays for miller
    #*********************************************************************
    def allocate_arrays(self, nr, nzeta):
        """ Periodic quantities can be computed on 2*pi grid and replicated. """
        
        
        # Arrays with dimensions (nr,-nz:nz)
        self.Rr = numpy.zeros((nr,nzeta)) 
        self.Zr = numpy.zeros((nr,nzeta))  
        
        # Arrays with dimensions (-nz:nz)
        self.zeta = numpy.zeros((1,nzeta))
        self.grho = numpy.zeros((nzeta))
        self.bmag = numpy.zeros((nzeta))
        self.gradpar = numpy.zeros((nzeta)) 
        self.gds2 = numpy.zeros((nzeta)) 
        self.gds21 = numpy.zeros((nzeta)) 
        self.gds22 = numpy.zeros((nzeta)) 
        self.gds23 = numpy.zeros((nzeta)) 
        self.gds24 = numpy.zeros((nzeta)) 
        self.gbdrift0 = numpy.zeros((nzeta)) 
        self.gbdrift = numpy.zeros((nzeta)) 
        self.cvdrift0 = numpy.zeros((nzeta)) 
        self.cvdrift = numpy.zeros((nzeta)) 
        self.jacrho = numpy.zeros((nzeta))  
        self.djacdrho = numpy.zeros((nzeta))  
        self.djacrdrho = numpy.zeros((nzeta))  
        self.d2jacdr2 = numpy.zeros((nzeta))  
        self.d2Rdrdth = numpy.zeros((nzeta))  
        self.d2Zdrdth = numpy.zeros((nzeta))  
        self.gpsi = numpy.zeros((nzeta))  
        self.dBdrho = numpy.zeros((nzeta))  
        self.dgradpardrho = numpy.zeros((nzeta))  
        self.d2Bdrdth = numpy.zeros((nzeta))  
        self.dgradparBdrho = numpy.zeros((nzeta))  
        self.dBdth = numpy.zeros((nzeta))  
        self.gradparb = numpy.zeros((nzeta))  
        self.dgradparBdrho = numpy.zeros((nzeta))  
        self.dcvdrift0drho = numpy.zeros((nzeta))  
        self.dgbdrift0drho = numpy.zeros((nzeta))  
        self.theta = numpy.zeros((nzeta))  
        self.dRdrho = numpy.zeros((nzeta))  
        self.dZdrho = numpy.zeros((nzeta))  
        self.dRdth = numpy.zeros((nzeta))  
        self.dZdth = numpy.zeros((nzeta))   
        self.gradrho_gradthet = numpy.zeros((nzeta))   
        self.gradthet2 = numpy.zeros((nzeta))   
        self.dgr2dr = numpy.zeros((nzeta))   
        self.dgpsi2dr = numpy.zeros((nzeta))   
        self.dgrgt = numpy.zeros((nzeta))   
        self.dgt2 = numpy.zeros((nzeta))   
        self.dgagr = numpy.zeros((nzeta))   
        self.dgagt = numpy.zeros((nzeta))   
        self.dga2 = numpy.zeros((nzeta))   
        self.d2Rdr2 = numpy.zeros((nzeta))   
        self.d2Zdr2 = numpy.zeros((nzeta))   
        self.d2Bdr2 = numpy.zeros((nzeta))   
        self.drz = numpy.zeros((nzeta))   
        self.drzdth = numpy.zeros((nzeta))   
        self.d2Rdr2dth = numpy.zeros((nzeta))   
        self.d2Rdth2 = numpy.zeros((nzeta))    
        self.d2Zdth2 = numpy.zeros((nzeta))    
        self.d2gpsidr2 = numpy.zeros((nzeta))    
        self.gradalph_gradthet = numpy.zeros((nzeta))    
        self.gradalph2 = numpy.zeros((nzeta))    
        self.gradrho_gradalph = numpy.zeros((nzeta))    
        self.dgds2dr = numpy.zeros((nzeta))     
        self.dgds21dr = numpy.zeros((nzeta))  
        self.dgds22dr = numpy.zeros((nzeta))   
        self.dcvdriftdrho = numpy.zeros((nzeta))   
        self.dgbdriftdrho = numpy.zeros((nzeta))   
        self.varthet = numpy.zeros((nzeta))   
        self.dvarthdr = numpy.zeros((nzeta))   
        self.d2varthdr2 = numpy.zeros((nzeta))     
        self.cross = numpy.zeros((nzeta))       
        self.dcrossdr = numpy.zeros((nzeta))   
        return
     
    #*********************************************************************
    #                   Calculate Miller coordinates
    #*********************************************************************
    def get_local_geo(self):

        # Number of grid points used for radial derivatives
        nr = 3
    
        # First get nperiod corresponding to input number of grid points
        nz2pi = self.nzed/2
        np = (self.nzgrid-nz2pi)/self.nzed + 1
    
        # Now switch to using (possible higher resolution) local grid
        nz2pi = self.nzed_local/2
        
        # This is the equivalent of nzgrid on the local grid
        nz = int(nz2pi + self.nzed_local*(np-1))
        
        # Save the z-variables on the local grid
        self.nzed = self.nzed_local
        self.nzgrid = int(self.nzed/2) + int((np-1)*self.nzed)  
        self.nzeta  = self.nzgrid*2+1
        self.zed = [ i*numpy.pi/(self.nzed/2) for i in numpy.arange(-self.nzgrid, self.nzgrid+1,1) ]
    
        # Allocate the arrays for the Miller coordinates
        self.allocate_arrays(nr, 2*nz+1)

        # Read the geometry from a text file
        if (self.read_profile_variation):
            
            # Open and read the text file
            f = open('RZ.in','r')
            text = f.readlines()
            text = text.split(" "); i = 0
            text = [t for t in text if (t!=" " and t!="\n" and t!="\\n")]
            
            # Read the variables
            rhoc0           = text[i]; i += 1
            dI              = text[i]; i += 1
            qinp            = text[i]; i += 1
            shat            = text[i]; i += 1
            d2qdr2          = text[i]; i += 1
            kappa           = text[i]; i += 1
            kapprim         = text[i]; i += 1
            tri             = text[i]; i += 1
            triprim         = text[i]; i += 1
            betaprim        = text[i]; i += 1
            betadbprim      = text[i]; i += 1
            dpsidrho_psi0   = text[i]; i += 1
                
            # Read the vectors
            d2R = numpy.zeros((2*nz+1))
            d2Z = numpy.zeros((2*nz+1))
            theta = numpy.zeros((2*nz+1))
            bmag_psi0 = numpy.zeros((2*nz+1))
            grho_psi0 = numpy.zeros((2*nz+1))
            for j in range(2*nz+1):
                theta[j]     = text[i]; i += 1
                d2R[j]       = text[i]; i += 1
                d2Z[j]       = text[i]; i += 1
                bmag_psi0[j] = text[i]; i += 1
                grho_psi0[j] = text[i]; i += 1
                
            # Close the text file
            f.close()
            
            # Interpolate the data at rho=rhoc0 (given in RZ.in) to rho=rhoc (given in input.in)
            self.LOCAL.qinp = qinp + shat*qinp/rhoc0*(self.LOCAL.rhoc-rhoc0) + 0.5*(self.LOCAL.rhoc-rhoc0)**2*d2qdr2
            self.LOCAL.shat = (self.LOCAL.rhoc/self.LOCAL.qinp) * (shat*qinp/rhoc0 + (self.LOCAL.rhoc-rhoc0)*d2qdr2)
            self.LOCAL.kappa = kappa + kapprim*(self.LOCAL.rhoc-rhoc0)
            self.LOCAL.tri = tri + triprim*(self.LOCAL.rhoc-rhoc0)
            self.LOCAL.betaprim = betaprim + betadbprim*(self.LOCAL.rhoc-rhoc0)
            self.LOCAL.rhoc_psi0 = rhoc0
            self.LOCAL.qinp_psi0 = qinp
            self.LOCAL.shat_psi0 = shat
            
            # Save these variables
            self.rhoc0 = rhoc0
            self.dI = dI
            self.qinp = qinp         
            self.shat = shat           
            self.d2qdr2 = d2qdr2         
            self.kappa = kappa          
            self.kapprim = kapprim      
            self.tri = tri            
            self.triprim = triprim        
            self.betaprim = betaprim        
            self.betadbprim = betadbprim      
            self.dpsidrho_psi0 = dpsidrho_psi0 
    
        # Calculate dq/dr
        dqdr = self.LOCAL.shat*self.LOCAL.qinp/self.LOCAL.rhoc
        
        # Something? With self.LOCAL.dr = 1.e-3*(self.LOCAL.rhoc/self.LOCAL.rmaj)
        dr = numpy.zeros((3))
        dr[0] = -self.LOCAL.dr
        dr[1] = 0.
        dr[2] = self.LOCAL.dr
        
        # Calculate R(nr,theta) and Z(nr,theta)
        for j in range(2*nz+1):
            self.theta[j] = (j-nz)*(2*np-1)*numpy.pi/nz
            for i in [0,1,2]:
                rmin = self.LOCAL.rhoc + dr[i]
                self.Rr[i,j] = calculate_Rpos(rmin,self.theta[j],j,self,self.LOCAL,nz)
                self.Zr[i,j] = calculate_Zpos(rmin,self.theta[j],j,self,self.LOCAL,nz)
        Rr = self.Rr
        Zr = self.Zr 

        # Get delta theta as a function of theta
        delthet = numpy.zeros(2*nz) #(-nz:nz-1))
        delthet = self.theta[1:] - self.theta[0:-1] 
     
        # Get dR/drho and dZ/drho
        dRdrho = self.dRdrho = get_drho(Rr, self.LOCAL)
        dZdrho = self.dZdrho = get_drho(Zr, self.LOCAL) 
     
        # Get dR/dtheta and dZ/dtheta
        dRdth = self.dRdth = get_dthet(Rr[1,:], self.dRdth, delthet)
        dZdth = self.dZdth = get_dthet(Zr[1,:], self.dZdth, delthet)
     
        # Get second derivatives of R and Z with respect to theta
        d2Rdth2 = self.d2Rdth2 = get_d2dthet2(Rr[1,:], self.d2Rdth2, delthet)
        d2Zdth2 = self.d2Zdth2 = get_d2dthet2(Zr[1,:], self.d2Zdth2, delthet)
        
        # Get mixed theta and rho derivatives of R and Z
        d2Rdrdth = self.d2Rdrdth = get_dthet(dRdrho, self.d2Rdrdth, delthet)
        d2Zdrdth = self.d2Zdrdth = get_dthet(dZdrho, self.d2Zdrdth, delthet)
     
        # Get the Jacobian of the transformation from (rho,theta,zeta) to (R,Z,zeta)
        # this is what I call jacr or jacrho in following comments, as opposed to jacobian, 
        # which is for tranformation from (psi,theta,zeta) to (R,Z,zeta)
        jacrho = self.jacrho = Rr[1,:]*(dRdrho*dZdth - dRdth*dZdrho)
     
        # <theta_integrate> returns integral from 0 -> 2*pi, note that dpsidrho 
        # here is an intermediary that requires manipulation to get final dpsidrho 
        imin = int(-nz2pi+nz); imax = int(nz2pi+nz+1)  
        dpsidrho = theta_integrate(jacrho[imin:imax]/Rr[1,imin:imax]**2, delthet, nz2pi, nz)
        dpsidrho = dpsidrho/(2*numpy.pi)
     
        # Get dpsinorm/drho = (I/2*pi*q)*int_0^{2*pi} dthet jacrho/R**2
        # If using input.profiles, we are given dpsitordrho and must use it to compute rgeo
        if (abs(self.LOCAL.dpsitordrho) > sys.float_info.epsilon):
            self.LOCAL.rgeo = self.LOCAL.dpsitordrho/dpsidrho
            dpsidrho = self.LOCAL.dpsitordrho/self.LOCAL.qinp
            self.LOCAL.d2psidr2 = (self.LOCAL.d2psitordrho2-self.LOCAL.dpsitordrho*self.LOCAL.shat/self.LOCAL.hoc) / self.LOCAL.qinp
            # I=Btor*R is a flux function so bi = I/(Btor(psi,theta of Rgeo)*a) = Rgeo/a
            bi = self.LOCAL.rgeo
        else:
            # otherwise, we are given rgeo and must use it to compute dpsidrho
            # I=Btor*R is a flux function so bi = I/(Btor(psi,theta of Rgeo)*a) = Rgeo/a
            bi = self.LOCAL.rgeo + self.dI*(self.rhoc-self.rhoc0)
            dpsidrho = dpsidrho*bi/self.LOCAL.qinp 
     
        # Get |grad rho| and |grad psi|
        grho = self.grho = Rr[1,:]*numpy.sqrt(dRdth**2 + dZdth**2)/jacrho 
        gpsi = grho*dpsidrho
     
        # Quantity needed in calculation of dI/drho and djacrho/drho
        drz = (dRdrho*dRdth + dZdrho*dZdth)/jacrho
        drzdth = self.drzdth = get_dthet(drz, self.drzdth, delthet)
     
        # Get dI/drho
        dIdrho = get_dIdrho(dpsidrho, grho, nz, nz2pi, jacrho, bi, gpsi, Rr, dRdrho, dqdr, dRdth, dZdth, d2Zdrdth, drzdth, d2Rdrdth, delthet, self.LOCAL)
        self.dIdrho_out = dIdrho
      
        # Get djacobian/drho*dpsi/drho and djacr/drho
        # this is dpsi/dr * d/dr (jacobian) with betaprim below is (4*pi*ptot/B0^2)*(-d ln ptot / drho)
        djacdrho = self.djacdrho = (Rr[1,:]/grho)**2*(2*(dRdth*d2Rdrdth+dZdth*d2Zdrdth)/jacrho \
            - drzdth + jacrho*(bi*dIdrho/Rr[1,:]**2 - self.LOCAL.betaprim)/dpsidrho**2)
     
        # This is d/dr (jacobian_r)
        djacrdrho = self.djacrdrho = djacdrho + jacrho*self.LOCAL.d2psidr2/dpsidrho
     
        # Get d2R/drho2 and d2Z/drho2: get factor common to both d2R/drho2 and d2Z/drho2
        d2Rdr2 = ((djacrdrho-jacrho*dRdrho/Rr[1,:])/Rr[1,:] - dRdrho*d2Zdrdth + dZdrho*d2Rdrdth)/(dRdth**2+dZdth**2)
        d2Zdr2 = -d2Rdr2*dRdth
        d2Rdr2 = d2Rdr2*dZdth
        d2R = self.d2Rdr2 = d2Rdr2
        d2Z = self.d2Zdr2 = d2Zdr2
         
        # Get theta derivative of d2R/drho2 and d2Z/drho2
        d2Rdr2dth = self.d2Rdr2dth = get_dthet(d2Rdr2, self.d2Rdr2dth, delthet)
        d2Rdr2dth = self.d2Zdr2dth = get_dthet(d2Zdr2, self.d2Rdr2dth, delthet)
      
        # Calculate the magnitude of B (normalized by B(psi,theta corresponding to Rgeo))
        # B/B0 = sqrt(I**2 + |grad psi|**2)/R
        bmag = self.bmag = numpy.sqrt(bi**2 + gpsi**2)/Rr[1,:] 
 
        # the next line is for multibox runs
        if (self.load_psi0_variables):
            dpsidrho_psi0 = dpsidrho
            bmag_psi0 = bmag
            grho_psi0 = grho 
        
        if (self.write_profile_variation):
            file = open('RZ.out', 'w')
            text = self.LOCAL.rhoc, dIdrho, self.LOCAL.qinp, self.LOCAL.shat, self.LOCAL.d2qdr2, \
                self.LOCAL.kappa, self.LOCAL.kapprim, self.LOCAL.tri, self.LOCAL.triprim, \
                self.LOCAL.betaprim, self.LOCAL.betadbprim, dpsidrho, "\n"
            for j in range(2*nz+1):
                text += str(theta[j]) + str(d2Rdr2[j]) + str(d2Zdr2[j]) + str(bmag[j]) + str(grho[j]) + "\n"
            file.write(text)  
            file.close() 
 
        # Calculate dB/dtheta
        dBdth = self.dBdth = get_dthet(bmag, self.dBdth, delthet)
         
        # Calculate b . grad theta
        gradpar = self.gradpar = dpsidrho/(bmag*jacrho)
        
        # Calculate b . grad B
        gradparb = self.gradparb = gradpar*dBdth
     
        # Calculate d|grad rho|^2/drho and d|grad psi|^2/drho
        dgr2dr = self.dgr2dr = 2*(grho**2*(dRdrho/Rr[1,:]-djacrdrho/jacrho) + (Rr[1,:]/jacrho)**2*(dRdth*d2Rdrdth + d2Zdrdth*dZdth))
        dgpsi2dr = self.dgpsi2dr = 2*(gpsi**2*(dRdrho/Rr[1,:]-djacdrho/jacrho) + (Rr[1,:]/jacrho)**2*(dRdth*d2Rdrdth + d2Zdrdth*dZdth)*dpsidrho**2)
    
        # Calculate dB/drho and d2B/drho2 
        dBdrho = self.dBdrho = (bi*dIdrho + 0.5*dgpsi2dr) / (bmag*Rr[1,:]**2) - bmag*dRdrho/Rr[1,:]
     
        # Calculate d (b . grad theta) / drho
        dgradpardrho = self.dgradpardrho = -gradpar*(dBdrho/bmag + djacdrho/jacrho)
     
        # Calculate d/dtheta (dB/drho)
        d2Bdrdth = self.d2Bdrdth = get_dthet(dBdrho, self.d2Bdrdth, delthet)
         
        # Calculate d(b . grad B)/drho
        dgradparBdrho = self.dgradparBdrho = dgradpardrho*dBdth + gradpar*d2Bdrdth
     
        # Calculate varthet = (I/(q*(dpsi/dr)) * int_0^theta dtheta' jacrho/R^2 
        varthet = theta_integrate_indef(jacrho/Rr[1,:]**2, self.varthet, delthet, nz)
        varthet = self.varthet = bi*varthet/(dpsidrho*self.LOCAL.qinp)
     
        # Calculate dvarthet/drho 
        dum = bi*jacrho*( dIdrho/bi - dqdr/self.LOCAL.qinp + djacdrho/jacrho - 2*dRdrho/Rr[1,:] )/Rr[1,:]**2
        dvarthdr = theta_integrate_indef(dum, self.dvarthdr, delthet, nz)
        dvarthdr = self.dvarthdr = dvarthdr/(dpsidrho*self.LOCAL.qinp)
     
        # Calculate |grad theta|^2, grad r . grad theta, grad alpha . grad theta, etc.
        gradthet2 = self.gradthet2 = (Rr[1,:]/jacrho)**2*(dRdrho**2 + dZdrho**2)
        gradrho_gradthet = self.gradrho_gradthet = -(Rr[1,:]/jacrho)**2*(dRdrho*dRdth+dZdrho*dZdth)
        gradalph_gradthet = self.gradalph_gradthet = -(varthet*dqdr + self.LOCAL.qinp*dvarthdr)*gradrho_gradthet - bi*jacrho/(dpsidrho*Rr[1,:]**2)*gradthet2
        gradrho_gradalph = self.gradrho_gradalph = -(varthet*dqdr + self.LOCAL.qinp*dvarthdr)*grho**2 - bi*jacrho/(dpsidrho*Rr[1,:]**2)*gradrho_gradthet
        gradalph2 = self.gradalph2 = (1./Rr[1,:]**2) + ((varthet*dqdr+self.LOCAL.qinp*dvarthdr)*grho)**2 \
                    + 2.*bi*jacrho*(varthet*dqdr+self.LOCAL.qinp*dvarthdr)*gradrho_gradthet/(dpsidrho*Rr[1,:]**2) \
                    + (bi*jacrho/(dpsidrho*Rr[1,:]**2))**2*gradthet2
     
        # Calculate gds2 = |grad alpha|^2 * (dpsiN/drho)^2 (dpsiN/drho factor accounts for ky normalization)
        gds2 = self.gds2 = gradalph2*dpsidrho_psi0**2
        
        # Calculate gds21 = (grad q . grad alpha) * (dpsiN/drho)^2
        gds21 = self.gds21 = gradrho_gradalph*dqdr*dpsidrho_psi0**2
        
        # Calculate gds22 =|grad q|^2 * (dpsiN/drho)^2
        gds22 = self.gds22 = (grho*dpsidrho_psi0*dqdr)**2
        
        # Calculate gds23 = (grad rho . grad theta * |grad alpha|^2 - grad alpha . grad theta * grad rho . grad alpha) * (dpsiN/drho)^2 / B^2
        gds23 = self.gds23 = (gradrho_gradthet*gradalph2-gradalph_gradthet*gradrho_gradalph)*(dpsidrho_psi0/bmag)**2
        
        # Calculate gds24 = (grad rho . grad theta * grad rho . grad alpha - grad alpha . grad theta * |grad rho|^2) * (dpsiN/drho)^2 / B^2 * q/rho
        gds24 = self.gds24 = (gradrho_gradthet*gradrho_gradalph-gradalph_gradthet*grho**2) \
                 *(dpsidrho_psi0/bmag)**2*(self.LOCAL.qinp_psi0/self.LOCAL.rhoc_psi0)
    
        # Calculate cross = (grad alpha x B) . grad theta
        cross = self.cross = dpsidrho*(gradrho_gradalph*gradalph_gradthet - gradalph2*gradrho_gradthet)
        
        # Calculate gbdrift = bhat/B x (grad B/B) . grad alpha * 2 * dpsiN/drho
        gbdrift = self.gbdrift = 2.0*(-dBdrho + cross*dBdth*dpsidrho/bmag**2)/bmag
        
        # Calculate cvdrift = bhat/B x (bhat . grad bhat) . grad alpha * 2 * dpsiN/drho
        # Assuming betaprim = 4*pi*ptot/B0^2 * (-d ln ptot / drho)
        cvdrift = self.cvdrift = (gbdrift + 2.0*self.LOCAL.betaprim/bmag**2)
        
        # Calculate cvdrift0 = 2 *(bhat/B x grad B / B) . (grad q) * dpsiN/drho / (bhat . grad B)
        # same as usual GS2 definition once bhat . grad B is added in below
        cvdrift0 = self.cvdrift0 = -2.*bi*dqdr/bmag**2
        
        # Calculate dcvdrift0drho = 2*dpsiN/drho times the rho derivative (bhat/B x grad B / B) . (grad q)
        dcvdrift0drho = self.dcvdrift0drho = cvdrift0*(dgradparBdrho + gradparb*(dIdrho/bi - 2.*dBdrho/bmag \
            - self.LOCAL.d2psidr2/dpsidrho))  - 2.*bi*gradparb*self.LOCAL.d2qdr2/bmag**2
        
        # Calculate dgbdrift0drho = 2*dpsiN/drho/B times the rho derivative of (bhat x gradB/B) . (grad q)
        # note that there's an extra factor of 1/B that's not expanded due to v_perp -> mu
        dgbdrift0drho = self.dgbdrift0drho = cvdrift0*(dgradparBdrho + gradparb*(dIdrho/bi - dBdrho/bmag \
            - self.LOCAL.d2psidr2/dpsidrho)) - 2.*bi*gradparb*self.LOCAL.d2qdr2/bmag**2
        
        cvdrift0 = self.cvdrift0 = cvdrift0*gradparb
        
        # Calculate gbdrift0 = 2 * dpsiN/drho * (bhat/B x gradB/B) . (grad q)
        gbdrift0 = cvdrift0
        
#         # Calculate d^2I/drho^2 and d^2 Jac / dr^2
#         call get_d2Idr2_d2jacdr2 (grho, dIdrho)
#         
#         # Calculate d^2varhteta/drho^2
#         call get_d2varthdr2 (dpsidrho, dIdrho)
#         
#         # Calculate d2B/drho^2
#         call get_d2Bdr2 (bmag, dIdrho)
#         
#         # Calculate d/dr [(grad alpha x B) . grad theta]
#         call get_dcrossdr (dpsidrho, dIdrho, grho)
#         
#         ! dgbdriftdrho is d/drho [(bhat/B x (grad B) . grad alpha) * 2 * dpsiN/drho] / B
#         ! note that there's an extra factor of 1/B that's not expanded due to v_perp -> mu
#         dgbdriftdrho = 2.0*(local%d2psidr2*dBdrho/dpsidrho - d2Bdr2 &
#              + dpsidrho*(dcrossdr*dBdth+cross*(d2Bdrdth-2.*dBdth*dBdrho/bmag))/bmag**2)/bmag
#         ! dcvdriftdrho is d/drho (bhat/B x [bhat . grad bhat] . grad alpha) * 2 * dpsiN/drho
#         dcvdriftdrho = dgbdriftdrho - gbdrift*dBdrho/bmag &
#              + 2.0*local%betadbprim/bmag**2 - 4.0*local%betaprim*dBdrho/bmag**3 &
#              - 2.0*local%betaprim*local%d2psidr2/dpsidrho
 
        # The next two sets of lines are corrections needed for the side boxes in a multibox simulation
        #  gbdrift  = gbdrift *(dpsidrho_psi0/dpsidrho)*(bmag/bmag_psi0)
        #  gbdrift0 = gbdrift0*(dpsidrho_psi0/dpsidrho)*(bmag/bmag_psi0)
        gbdrift  = gbdrift *(dpsidrho_psi0/dpsidrho)
        gbdrift0 = gbdrift0*(dpsidrho_psi0/dpsidrho)
        cvdrift  = cvdrift *(dpsidrho_psi0/dpsidrho)
        cvdrift0 = cvdrift0*(dpsidrho_psi0/dpsidrho)
 
#         #  dgbdriftdrho  = dgbdriftdrho *(dpsidrho_psi0/dpsidrho)*(bmag/bmag_psi0)
#         #  dgbdrift0drho = dgbdrift0drho*(dpsidrho_psi0/dpsidrho)*(bmag/bmag_psi0)
#         dgbdriftdrho  = dgbdriftdrho *(dpsidrho_psi0/dpsidrho)
#         dgbdrift0drho = dgbdrift0drho*(dpsidrho_psi0/dpsidrho)
#         dcvdriftdrho  = dcvdriftdrho *(dpsidrho_psi0/dpsidrho)
#         dcvdrift0drho = dcvdrift0drho*(dpsidrho_psi0/dpsidrho)

        # Interpolate here 
        self.bmag_out = geo_spline(self.theta, bmag, self.zed_in)
        self.rmajor_out = geo_spline(self.theta, Rr[1,:], self.zed_in) 
        self.grho_out = geo_spline(self.theta, grho_psi0, self.zed_in) 
        self.bmag_psi0_out = geo_spline(self.theta, bmag_psi0, self.zed_in)
        self.gds2_out = geo_spline(self.theta, gds2, self.zed_in)
        self.gds21_out = geo_spline(self.theta, gds21, self.zed_in)
        self.gds22_out = geo_spline(self.theta, gds22, self.zed_in)
        self.gds23_out = geo_spline(self.theta, gds21, self.zed_in)
        self.gds24_out = geo_spline(self.theta, gds21, self.zed_in)
        self.gradpar_out = geo_spline(self.theta, gradpar, self.zed_in)
        self.gbdrift_out = geo_spline(self.theta, gbdrift, self.zed_in)
        self.gbdrift0_out = geo_spline(self.theta, gbdrift0, self.zed_in)
        self.cvdrift_out = geo_spline(self.theta, cvdrift, self.zed_in)
        self.cvdrift0_out = geo_spline(self.theta, cvdrift0, self.zed_in)
        self.dBdrho_out = geo_spline(self.theta, dBdrho, self.zed_in)
        self.d2Bdrdth_out = geo_spline(self.theta, d2Bdrdth, self.zed_in)
        self.dgradpardrho_out = geo_spline(self.theta, dgradpardrho, self.zed_in)
        self.dcvdrift0drho_out = geo_spline(self.theta, dcvdrift0drho, self.zed_in)
        self.dgbdrift0drho_out = geo_spline(self.theta, dgbdrift0drho, self.zed_in)
        self.djacdrho_out = geo_spline(self.theta, djacdrho/dpsidrho, self.zed_in)
#         self.dcvdriftdrho_out = geo_spline(self.theta, dcvdriftdrho, self.zed_in)
#         self.dgbdriftdrho_out = geo_spline(self.theta, dgbdriftdrho, self.zed_in)
#         self.dgds2dr_out = geo_spline(self.theta, dgds2dr, self.zed_in)
#         self.dgds21dr_out = geo_spline(self.theta, dgds21dr, self.zed_in)
#         self.dgds22dr_out = geo_spline(self.theta, dgds22dr, self.zed_in)

        # Get the toroidal component of the magnetic field
        # btor = B_toroidal/Bref = I/R Bref = rgeo * a/R
        self.btor = bi/Rr[1,:]
        self.btor_out = self.Btor = bi/self.rmajor_out 
        self.bpol = gpsi/Rr[1,:] 

        # Create zeta for Miller 
        self.zeta[0,:] = numpy.array(self.zed)*self.LOCAL.qinp
        
        # Create zeta for stella
        self.nzed_out = self.nzed_in
        self.nzgrid_out = int(self.nzed_in/2) + int((np-1)*self.nzed_in)  
        self.nzeta_out  = self.nzgrid_out*2+1
        self.zed_out = [ i*numpy.pi/(self.nzed_in/2) for i in numpy.arange(-self.nzgrid_out, self.nzgrid_out+1,1) ]
        self.zeta_out[0,:] = numpy.array(self.zed_in)*self.LOCAL.qinp
        
        # Check the data with millerlocal.input
        if False:
            self.LOCAL.rhotor = numpy.nan
            self.LOCAL.drhotordrho = numpy.nan
            self.LOCAL.psitor_lcfs = numpy.nan
            print("-------------------------------")
            print('#1.rhoc', '2.rmaj', '3.rgeo', '4.shift', '5.qinp')
            print(self.LOCAL.rhoc, self.LOCAL.rmaj, self.LOCAL.rgeo, self.LOCAL.shift, self.LOCAL.qinp)
            print('#6.shat', '7.kappa', '8.kapprim', '9.tri', '10.triprim')
            print(self.LOCAL.shat, self.LOCAL.kappa, self.LOCAL.kapprim, self.LOCAL.tri, self.LOCAL.triprim)
            print('11.betaprim', '12.dpsitordrho', '13.rhotor', '14.drhotordrho', '15.d2qdr2')
            print(self.LOCAL.betaprim, self.LOCAL.dpsitordrho, self.LOCAL.rhotor, self.LOCAL.drhotordrho, self.LOCAL.d2qdr2)
            print('16.d2psidr2', '17.betadbprim', '18.psitor_lcfs')
            print(self.LOCAL.d2psidr2, self.LOCAL.betadbprim, self.LOCAL.psitor_lcfs)
        
        # CBC stella output
        stella = numpy.loadtxt("/home/hanne/CIEMAT/RUNS/NONLINEAR_CBC/TIPRIM3TEPRIM3/test/test.output", dtype='float').reshape(-1, 64) 
        
        # Check the data with millerlocal.output 
        if False:
            print("-------------------------------")
            print("shape:     ", numpy.shape(stella))
            print("theta:     ", numpy.all(numpy.round(stella[:,0],6)==numpy.round(self.theta,6)))  
            print("R:         ", numpy.all(numpy.round(stella[:,1],6)==numpy.round(self.Rr[1,:],6)))  
            print("Z:         ", numpy.all(numpy.round(stella[:,57],6)==numpy.round(self.Zr[1,:],6))) 
            print("dR/dr:     ", numpy.all(numpy.round(stella[:,2],6)==numpy.round(self.dRdrho,6)))  
            print("dZ/dr:     ", numpy.all(numpy.round(stella[:,6],6)==numpy.round(self.dZdrho,6))) 
            print("dR/dth:    ", numpy.all(numpy.round(stella[:,4],6)==numpy.round(self.dRdth,6)))  
            print("dZ/dth:    ", numpy.all(numpy.round(stella[:,8],6)==numpy.round(self.dZdth,6))) 
            print("d2Rdrdth:  ", numpy.all(numpy.round(stella[:,5],6)==numpy.round(self.d2Rdrdth,6))) 
            print("d2Zdrdth:  ", numpy.all(numpy.round(stella[:,9],6)==numpy.round(self.d2Zdrdth,6)))
            
            print("jacr:      ", numpy.all(numpy.round(stella[:,18],6)==numpy.round(self.jacrho,6))) 
            print("grho2:     ", numpy.all(numpy.round(stella[:,22],6)==numpy.round(self.grho**2,6)))  
    
            print("dpsidrho:  ", numpy.all(numpy.round(stella[0,58],6)==numpy.round(dpsidrho,6)))
            print("bi:        ", numpy.all(numpy.round(stella[0,59],6)==numpy.round(bi,6)))
            print("gpsi:      ", numpy.all(numpy.round(stella[:,60],6)==numpy.round(gpsi,6)))
            print("drz:       ", numpy.all(numpy.round(stella[:,61],6)==numpy.round(drz,6)))
            print("drzdth:    ", numpy.all(numpy.round(stella[:,62],6)==numpy.round(drzdth,6)))
            print("dIdrho:    ", numpy.all(numpy.round(stella[0,63],6)==numpy.round(dIdrho,6)))
            
            print("djacrdr:   ", numpy.all(numpy.round(stella[:,19],5)==numpy.round(self.djacrdrho,5)))  
            print("djacdrho:  ", numpy.all(numpy.round(stella[:,20],5)==numpy.round(self.djacdrho,5)))     
            
            print("d2Rdr2:    ", numpy.all(numpy.round(stella[:,3],5)==numpy.round(self.d2Rdr2,5)))   
            print("d2Zdr2:    ", numpy.all(numpy.round(stella[:,7],5)==numpy.round(self.d2Zdr2,5))) 
             
            print("bmag:      ", numpy.all(numpy.round(stella[:,10],5)==numpy.round(self.bmag,5)))  
            
            print("dBdr:      ", numpy.all(numpy.round(stella[:,11],5)==numpy.round(self.dBdrho,5)))  
            print("d2Bdr2:    ", numpy.all(numpy.round(stella[:,12],5)==numpy.round(self.d2Bdr2,5))) 
            print("dB/dth':   ", numpy.all(numpy.round(stella[:,13],5)==numpy.round(self.dBdth,5))) 
            print("d2Bdrdth:  ", numpy.all(numpy.round(stella[:,14],5)==numpy.round(self.d2Bdrdth,5))) 
            print("varthet:   ", numpy.all(numpy.round(stella[:,15],5)==numpy.round(self.varthet,5))) 
            print("dvarthdr:  ", numpy.all(numpy.round(stella[:,16],5)==numpy.round(self.dvarthdr,5))) 
            print("d2varthdr2:", numpy.all(numpy.round(stella[:,17],5)==numpy.round(self.d2varthdr2,5))) 
            print("d2jacdr2:  ", numpy.all(numpy.round(stella[:,21],5)==numpy.round(self.d2jacdr2,5)))
            print("--------------- bmag ----------------") 
            print(self.bmag[:5])
            print(stella[:5,10])
        return self
#===============================================================================
#               SAVE THE STELLA COORDINATES TO A <STELLA> CLASS
#===============================================================================

class stella:
    
    def __init__(self, nzed, nperiod):
    
        # Remember nzed and nzgrid
        self.nzed = nzed
        self.nzgrid = int(nzed/2) + int((nperiod-1)*nzed)  
        self.nzeta  = self.nzgrid*2+1
        self.zed = [ i*numpy.pi/(nzed/2) for i in numpy.arange(-self.nzgrid, self.nzgrid+1,1) ]
        
        
    #*********************************************************************
    #                       Allocate arrays
    #*********************************************************************
    def allocate_arrays(self, nalpha, nzeta): 
        
        # Arrays with dimensions (nalpha,-nzgrid:nzgrid)
        self.bmag = numpy.zeros((nalpha, nzeta))
        self.bmag_psi0 = numpy.zeros((nalpha, nzeta))
        self.gds2 = numpy.zeros((nalpha, nzeta))
        self.gds21 = numpy.zeros((nalpha, nzeta))
        self.gds22 = numpy.zeros((nalpha, nzeta))
        self.gds22 = numpy.zeros((nalpha, nzeta))
        self.gds23 = numpy.zeros((nalpha, nzeta))
        self.gds24 = numpy.zeros((nalpha, nzeta))
        self.gds25 = numpy.zeros((nalpha, nzeta))
        self.gds26 = numpy.zeros((nalpha, nzeta))
        self.dgds2dr = numpy.zeros((nalpha, nzeta))
        self.dgds21dr = numpy.zeros((nalpha, nzeta))
        self.dgds22dr = numpy.zeros((nalpha, nzeta))
        self.gbdrift = numpy.zeros((nalpha, nzeta))
        self.gbdrift0 = numpy.zeros((nalpha, nzeta))
        self.cvdrift = numpy.zeros((nalpha, nzeta))
        self.cvdrift0 = numpy.zeros((nalpha, nzeta))
        self.dgbdriftdrho = numpy.zeros((nalpha, nzeta))
        self.dcvdriftdrho = numpy.zeros((nalpha, nzeta))
        self.dgbdrift0drho = numpy.zeros((nalpha, nzeta))
        self.dcvdrift0drho = numpy.zeros((nalpha, nzeta))
        self.dbdzed = numpy.zeros((nalpha, nzeta)) 
        self.theta_vmec = numpy.zeros((nalpha, nzeta)) 
        self.nalpha = numpy.zeros((nalpha, nzeta)) 
        self.djacdrho = numpy.zeros((nalpha, nzeta)) 
        self.grho = numpy.zeros((nalpha, nzeta)) 
        self.dl_over_b = numpy.zeros((nalpha, nzeta)) 
        self.d_dl_over_b_drho = numpy.zeros((nalpha, nzeta))  
        
        # Arrays with dimensions (-nzgrid:nzgrid)
        self.gradpar = numpy.zeros((nzeta))  
        self.zed_eqarc = numpy.zeros((nzeta))  
        self.btor = numpy.zeros((nzeta))  
        self.rmajor = numpy.zeros((nzeta))  
        self.dBdrho = numpy.zeros((nzeta))  
        self.d2Bdrdth = numpy.zeros((nzeta))  
        self.dgradpardrho = numpy.zeros((nzeta))   
        
        # Other arrays
        self.alpha = numpy.zeros((nalpha)) 
        self.zeta = numpy.zeros((nalpha, nzeta))  
        
        # Variables 
        self.dpsidrho_psi0 = numpy.nan
        self.dpsidrho = numpy.nan
        self.dIdrho = numpy.nan
        return 

#===============================================================================
#                                 RUN THIS SCRIPT
#===============================================================================
    
def test_calculateGoemetricQuantitiesMiller(device="CBC", rho=0.7): 
    
    # Details per device
    if device=="CBC" and rho==0.7:
            
        # Stella output
        path = "/home/hanne/CIEMAT/RUNS/NONLINEAR_CBC/TIPRIM3/fprim2tprim3/input.geometry"
        keys_geo = ["alpha","zed", "zeta", "bmag", "gradpar", "gds2", "gds21", "gds22", "gds23", "gds24", "gbdrift", "cvdrift", "gbdrift0", "bmag_psi0"]
    
        # Stela input
        nzed = 24
        nalpha = 1
        nperiod = 1
        miller_variables = {
            "nzed_local" : 128,
            "rhoc" : 0.7,
            "shat" : 0.796,
            "qinp" : 1.4,
            "rmaj" : 2.77778,
            "rgeo" : 2.77778,
            "shift" : 0.0,
            "kappa" : 1.0,
            "kapprim" : 0.0,
            "tri" : 0.0,
            "triprim" : 0.0,
            "betaprim" : 0.0,
            "d2qdr2" : 0.0,
            "d2psidr2" : 0.0,
            "betadbprim" : 0.0}

    # Create the figure
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec 
    fig = plt.figure(figsize=(18,9)); load_styleFigures(); axes = []
    grid_specifications = gridspec.GridSpec(2, 6, figure=fig)
    grid_specifications.update(top=0.93, left=0.05, right=0.98, bottom=0.07, wspace=0.4)
    axes.append(plt.subplot(grid_specifications[0:2])); j = 1  
    for i in range(2,12): axes.append(plt.subplot(grid_specifications[i]))   
    for i in range(1,11): axes[i].set_xlabel('$\\zeta$')
    keys_plot = ["bmag", "bmag_psi0", "gds2", "gds21", "gds22", "gds23", "gds24", "gbdrift", "cvdrift", "gbdrift0"]
    axes[0].set_xlabel('$R$ [m]')
    axes[0].set_ylabel('$Z$ [m]')
    axes[1].set_ylabel('$B$ [a.u.]')
    axes[2].set_ylabel('bmag psi0 [a.u.]') 
    axes[3].set_ylabel('gds2 [a.u.]') 
    axes[4].set_ylabel('gds21 [a.u.]')
    axes[5].set_ylabel('gds22 [a.u.]')
    axes[6].set_ylabel('gds23 [a.u.]')
    axes[7].set_ylabel('gds24 [a.u.]')
    axes[8].set_ylabel('gbdrift0 [a.u.]')
    axes[9].set_ylabel('gbdrift [a.u.]') 
    axes[10].set_ylabel('cvdrift [a.u.]') 
    
    # Calculate the Miller coefficients based on this python script
    MILLER = calculate_geometricQuantitiesMiller(miller_variables, nzed, nperiod, nalpha)
    
    # Read the Miller coefficients from the stella output
    data = numpy.loadtxt(path, skiprows=3); STELLA = type('Stella', (object,), {})
    for i in range(len(keys_geo)): setattr(STELLA, keys_geo[i], data[:,i])


    # Plot the cross-section
    Z = MILLER.Zr; R = MILLER.Rr
    axes[0].plot(R[1,:],Z[1,:],color="black",label='Python')
    update_limits(axes[0], R*1.1, Z*1.1, R*0.9, Z*1.1) 
    
    # Plot the geometric coefficients from Miller and Stella
    for key in keys_plot:
        x1 = getattr(MILLER, "zeta_out")[0,:]; y1 = getattr(MILLER, key+"_out")
        x2 = getattr(STELLA, "zeta"); y2 = getattr(STELLA, key)
        print(key+": "+str(y1[10:13])+" vs "+str(y2[10:13]))
        axes[j].plot(x1,y1,color="black",label='Python')
        axes[j].plot(x2,y2,color="red",label='Stella')
        update_limits(axes[j], x1, y1, x2, y2*1.1) 
        j += 1
    
    # Show figure 
    plt.show()
    return

def update_limits(ax, x1, y1, x2=[numpy.nan], y2=[numpy.nan]):
    xmin = 0; xmax = 0; ymin = 0; ymax = 0
    xmax = numpy.nanmax([xmax, numpy.max(x1)]); ymax = numpy.nanmax([ymax, numpy.max(y1)])
    xmin = numpy.nanmin([xmin, numpy.min(x1)]); ymin = numpy.nanmin([ymin, numpy.min(y1)])  
    xmax = numpy.nanmax([xmax, numpy.max(x2)]); ymax = numpy.nanmax([ymax, numpy.max(y2)])
    xmin = numpy.nanmin([xmin, numpy.min(x2)]); ymin = numpy.nanmin([ymin, numpy.min(y2)]) 
    ax.set_xlim(xmin,xmax); ax.set_ylim(ymin,ymax)  
    return 
 
if __name__ == "__main__":
    test_calculateGoemetricQuantitiesMiller()
