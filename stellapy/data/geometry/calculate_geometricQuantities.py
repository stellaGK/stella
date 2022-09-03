''' 

#===============================================================================
#                Calculate geometric quantities used by stella                 #
#===============================================================================

The goal is to reproduce all geometric quantities calculated in stella. 
Moreover, we want to reproduce the geometry file that is printed by stella, 
so that we can calculate it directly from vmec and get a higher resolution. 
This is the main script that will fill the following classes:
    {BIELD; COORD; GEO; GRAD; SURF; VMEC; ZGRID} 
It follows the geometry module of the stellac code. 

Geometry file
-------------
'# alpha', 'zed', 'zeta', 'bmag', 'gradpar', 'gds2', 'gds21', 'gds22', 
'gds23', 'gds24', 'gbdrift', 'cvdrift', 'gbdrift0', 'bmag_psi0


Input parameters
-----------------
nalpha is the number of grid points in the alpha coordinate
alpha0 is the first alpha value to include in the alpha grid
The zeta grid has nzgrid*2+1 points, including the "repeated" point at index -nzgrid and +nzgrid.
The zeta domain is centered at zeta_center. Setting zeta_center = 2*pi*N/nfp for any integer N should
yield identical results to setting zeta_center = 0, where nfp is the number of field periods (as in VMEC).
If number_of_field_periods_to_include is > 0, then this parameter does what you think:
the extent of the toroidal in zeta will be 2*pi*number_of_field_periods_to_include/nfp.
If number_of_field_periods_to_include is <= 0, the entire 2*pi toroidal domain will be included.
The parameter desired_normalized_toroidal_flux determines which flux surface from the VMEC file will be used
for the computation. This parameter should lie in the interval [0,1]. 


Output quantities
-----------------
On exit, normalized_toroidal_flux_used holds the flux surface that was actually 
used for the geometry, as measured by psi_toroidal / psi_{toroidal,edge} 
Safety factor q = 1/iota.
Magnetic shear shat = (x/q)*(d q / d x) where x = Aminor_p*np.sqrt(psi_toroidal / psi_{toroidal,edge})
and Aminor_p is the minor radius calculated by VMEC. 
L_reference is the reference length used for stella normalization, in meters.
B_reference is the reference magnetic field strength used for stella normalization, in Tesla.
nfp is the number of field periods given by VMEC
On exit, alpha holds the grid points in alpha = theta_p - iota*zeta, where theta_p is the PEST toroidal angle
On exit, zeta holds the grid points in the toroidal angle zeta
Zeta is the physical toroidal angle taken to increase in the counter-clockwise direction
<gradpar_zeta>: gradpar_zeta = b . grad zeta
<grad_alpha_grad_alpha>:  grad alpha . grad alpha in units of 1/L_ref^2, with alpha = theta_pest - iota*zeta
<grad_alpha_grad_psi> grad alpha . grad psi_t in units of B_ref, with alpha = theta_pest - iota*zeta
<grad_psi_grad_psi> grad psi_t . grad psi_t in units of (a*B_ref)^2, with alpha = theta_pest - iota*zeta
<gbdrift_alpha> 2*bhat/B x (grad B / B) . grad alpha*B_ref*L_ref^2
<cvdrift_alpha># 2*bhat/B x (bhat . grad bhat) . grad alpha*B_ref*L_ref^2


VMEC variables of interest
--------------------------
All VMEC quantities (B, pressure, etc) are in SI units. The quantities on the
half grid have the same number of array elements (ns) as quantities on the full grid.
but the first array element is 0.
    ns = number of flux surfaces used by VMEC
    nfp = number of field periods, e.g. 5 for W7-X, 4 for HSX
    iotas = rotational transform (1/q) on the half grid.
    iotaf = rotational transform on the full grid.
    presf = pressure on the full grid.
 
Integers
--------
nalpha, nzgrid, vmec_surface_option, sign_toroidal_flux
j, index, izeta, ialpha, interpolation, m, n, imn, imn_nyq
 
'''

#!/usr/bin/python3 
import numpy as np  
import sys, os, pathlib, pickle

# Stellapy package
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.data.geometry.calculate_geometricQuantitiesSTELLA import stella 
from stellapy.data.geometry.calculate_geometricQuantitiesKGRID import kgrid 
from stellapy.data.geometry.calculate_geometricQuantitiesVMEC import vmec
from stellapy.data.geometry.calculate_geometricQuantitiesZGRID import zgrid
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_iota
from stellapy.data.geometry.read_wout import read_woutFile
from stellapy.utils.files.ensure_dir import ensure_dir
from stellapy.utils.config import CONFIG 

#===============================================================================
#                Calculate geometric quantities used by stella                 #
#===============================================================================

def calculate_geometricQuantities(
        vmec_path, 
        # Radial location
        rho=0.7, 
        # Flux tube
        nfield_periods=None, 
        poloidal_turns=None,
        nzed=48, 
        nperiod=1, 
        # Choose between linked, stellarator, periodic or default (unlinked)
        boundary_option="linked",  
        zed_equal_arc=True,
        # If <grad_x_grad_y_end> is smaller than <grad_x_grad_y_zero> then period BC are used 
        grad_x_grad_y_zero=None,
        # For nonlinear simulations select (nx,ny,y0), for linear simulations select (kx,ky)
        kt_grids_knobs="range",  # Box is broken for now
        nx=91, ny=91, y0=15,
        kx=0, ky=2,
        nalpha=1,
        # Other
        save_as_pickle=True, 
        verbose=True): 

    # Read the vmec file to determine <nfield_periods> or <poloidal_turns>  
    boundary_option = boundary_option.replace("'","").replace('"','')
    if vmec_path==None: print("Please specify the path of the VMEC."); sys.exit()
    PATH = type('Path', (object,), {'vmec' : pathlib.Path(vmec_path), 'vmec_filename' : 'test', 'geometry' : "DoesntExist"})
    woutParameters = read_woutFile(PATH)    
    woutParameters['iota'] = calculate_iota(svalue=rho*rho, ns=woutParameters['ns'], iotas=woutParameters['iotas']) 
    if nfield_periods==None: nfield_periods = np.abs(poloidal_turns*woutParameters['nfp']/woutParameters['iota'])
    if poloidal_turns==None: poloidal_turns = np.abs(nfield_periods*woutParameters['iota']/woutParameters['nfp'])
    
    # If the quantities have been calculated before, read the pickle file
    pickly_path = get_picklePath(vmec_path, nzed, rho, poloidal_turns, nperiod, boundary_option, grad_x_grad_y_zero, zed_equal_arc)
    KGRID, STELLA, VMEC, ZGRID = read_pickleFile(vmec_path, pickly_path, verbose)
    if KGRID!=None: 
        if verbose: print("    READING PICKLE FILE:")
        if verbose: print("         "+str(pickly_path))
        return KGRID, STELLA, VMEC, ZGRID
    
    # Warning that this script takes a long time
    print("Calculating the geometric quantities used in stella. ")
    print("This might take a while ... ")
        
    # Initiate the classes that hold the geometric quantities
    KGRID = kgrid(kt_grids_knobs, verbose=verbose)
    STELLA = stella(verbose=verbose)
    VMEC = vmec(verbose=verbose)
    ZGRID = zgrid(verbose=verbose) 
    
    # Fill in the input parameters 
    KGRID.update_inputParameters(nalpha, kx, ky, nx, ny, y0)
    ZGRID.update_inputParameters(nzed, nperiod, boundary_option, zed_equal_arc=zed_equal_arc, grad_x_grad_y_zero=grad_x_grad_y_zero)
    VMEC.update_inputParameters(ZGRID, vmec_path, rho, nfield_periods)
    VMEC.poloidal_turns = poloidal_turns
    
    # Execute the statements in the <stella_geometry.f90> file
    init_geometry(KGRID, STELLA, VMEC, ZGRID) 
    
    # Save the data as a pickly file and return the data
    if save_as_pickle: save_pickleFile(vmec_path, pickly_path, KGRID, STELLA, VMEC, ZGRID)  
    return KGRID, STELLA, VMEC, ZGRID

#===============================================================================
#                                 Pickly files
#===============================================================================
def get_picklePath(vmec_path, nzed, rho, poloidal_turns, nperiod, boundary_option, grad_x_grad_y_zero, zed_equal_arc): 
    if "/" not in str(vmec_path): vmec_path = pathlib.Path(CONFIG['PATHS']['VMEC']+"/"+str(vmec_path)) 
    zero = "_zero"+str(round(grad_x_grad_y_zero,8)) if grad_x_grad_y_zero else ""
    pickly_details = "_nzed"+str(nzed)+"_rho"+str(rho)+"_pol"+str(round(poloidal_turns,8))+"_nperiod"+str(nperiod) 
    pickly_details += "_"+boundary_option+zero+"_arc"+str(zed_equal_arc)+"_KGRID"
    pickly_path = pathlib.Path(vmec_path).with_suffix('.pickle') 
    pickle_directory = pathlib.Path(str(pickly_path.parent)+"/Pickles"); ensure_dir(pickle_directory) 
    pickly_path = pickle_directory / (pickly_path.stem + pickly_details + pickly_path.suffix)
    return pickly_path

def read_pickleFile(vmec_path, pickly_path, verbose):  
    if (os.path.isfile(pickly_path)): 
        pickle_KGRID  = pickly_path
        pickle_STELLA = str(pickly_path).replace('KGRID', 'STELLA')
        pickle_VMEC   = str(pickly_path).replace('KGRID', 'VMEC')
        pickle_ZGRID  = str(pickly_path).replace('KGRID', 'ZGRID') 
        file = open(pickle_KGRID, 'rb');  KGRID  = pickle.load(file); file.close()
        file = open(pickle_STELLA,'rb');  STELLA = pickle.load(file); file.close()
        file = open(pickle_VMEC,  'rb');  VMEC   = pickle.load(file); file.close()
        file = open(pickle_ZGRID, 'rb');  ZGRID  = pickle.load(file); file.close()   
        if verbose: print("Read pickle files for", pathlib.Path(vmec_path).name) 
        if verbose: print("     ", pickly_path)
        return KGRID, STELLA, VMEC, ZGRID
    return None, None, None, None

def save_pickleFile(vmec_path, pickly_path, KGRID, STELLA, VMEC, ZGRID):
    pickle_KGRID  = pickly_path
    pickle_STELLA = str(pickly_path).replace('KGRID', 'STELLA')
    pickle_VMEC   = str(pickly_path).replace('KGRID', 'VMEC')
    pickle_ZGRID  = str(pickly_path).replace('KGRID', 'ZGRID')
    file = open(pickle_KGRID,  "wb"); pickle.dump(KGRID , file); file.close()
    file = open(pickle_STELLA, "wb"); pickle.dump(STELLA, file); file.close()
    file = open(pickle_VMEC,   "wb"); pickle.dump(VMEC,   file); file.close()
    file = open(pickle_ZGRID,  "wb"); pickle.dump(ZGRID,  file); file.close()
    print("\nSaved pickle files for", pathlib.Path(vmec_path).name)
    print("     ", pickly_path)
    return

#===============================================================================
#                             stella_geometry.f90
#===============================================================================

def init_geometry(KGRID, STELLA, VMEC, ZGRID, geo_option_switch='vmec'):
    
    # Stella routine
    if STELLA.verbose: print("stella::stella_geometry::init_geometry")
    
    # VMEC geometry
    if geo_option_switch=='vmec':
        
        # Execute the statements in the <vmec_geo.f90> file 
        get_vmec_geo(KGRID, STELLA, VMEC, ZGRID)
    
        # Finally, assemble the quantities needed for stella.
        STELLA.finish_stellaGeometryRoutine(VMEC, ZGRID)
        
    # Also calculate <theta0> and <kperp2> 
    #init_kperp2(KGRID, STELLA, ZGRID)    
    return
    
#===============================================================================
#                                 vmec_geo.f90
#===============================================================================

    #*********************************************************************
    #                Get the geometry coefficients from vmec
    #*********************************************************************
def get_vmec_geo (KGRID, STELLA, VMEC, ZGRID, verbose=False):

    # First read in equilibrium information from vmec file 
    VMEC.read_vmec_equilibrium()
    
    # We read more field periods from VMEC if (zed_equal_arc==True)
    VMEC.get_numberOfFieldPeriods()
        
    # Calculate <nzgrid_vmec>, which is the number of positive/negative zeta  
    # locations at which to get geometry data from vmec. 
    ZGRID.get_zgrid_vmec(VMEC)
    
    # Vmec to stella geometry interface
    vmec_to_stella_geometry_interface(KGRID, VMEC, ZGRID)
    
    # Copy some VMEC and ZGRID data to STELLA
    STELLA.add_dataFromVmecAndZgrid(VMEC, ZGRID)
     
    # Get ratio of number of simulated field periods to the number of field periods of the device
    STELLA.set_fieldPeriodRatio(VMEC.nfield_periods/VMEC.nfp)
    
    # Initiate the arrays for the stella quantities with nzgrid <= nzgrid_vmec
    STELLA.initiate_arraysForStella(KGRID.nalpha, ZGRID.nzgrid, int(2*ZGRID.nzgrid+1))
    
    # Interpolate the geometric quantities from (zeta,alpha) grid to
    # (zed,alpha) grid, with zed the arc-length
    STELLA.interpolate_geometric_quantities(KGRID, STELLA, VMEC, ZGRID)
     
    # Finish the remaining calculations
    STELLA.finish_vmecGeoRoutine(VMEC)
     
    # Write the 'vmec.geo' file   
    if verbose: 
        with open("/home/hanne/CIEMAT/MAGNETICFIELDS/W7X/geometry_test1.txt", 'w') as f:
            f.write("\n  ------------------------------------------------------------------------------------------------")
            f.write('\n        rhotor           qinp            shat           aref            Bref         z_scalefac')
            f.write("\n  ------------------------------------------------------------------------------------------------")
            f.write('\n' + '{:16.4e}'.format(STELLA.rhoc) +
                    '{:16.4e}'.format(STELLA.qinp) +
                    '{:16.4e}'.format(STELLA.shat) +
                    '{:16.4e}'.format(STELLA.L_reference) +
                    '{:16.4e}'.format(STELLA.B_reference) +
                    '{:16.4e}'.format(VMEC.zed_scalefac))
            f.write("\n\n\n  ------------------------------------------------------------------------------------------------")
            f.write('\n        alpha            zed             zeta           bmag          theta_vmec       gradpar')
            f.write("\n  ------------------------------------------------------------------------------------------------")
            for j in range(2*STELLA.nzgrid+1):
                for i in range(STELLA.nalpha):
                    f.write('\n' + '{:16.4e}'.format(STELLA.alpha[i]) +
                          '{:16.4e}'.format(STELLA.zed[j]) +
                          '{:16.4e}'.format(STELLA.zeta[i,j]) +
                          '{:16.4e}'.format(STELLA.bmag[i,j]) +
                          '{:16.4e}'.format(STELLA.theta_vmec[i,j]) +
                          '{:16.4e}'.format(STELLA.gradpar[j]))
        with open("/home/hanne/CIEMAT/MAGNETICFIELDS/W7X/geometry_test2.txt", 'w') as f:
            f.write("\n  -----------------------------------------------------------------------------------------------------------------------------------------------")
            f.write('\n        rhotor           qinp            shat           aref            Bref         z_scalefac')
            f.write("\n  -----------------------------------------------------------------------------------------------------------------------------------------------")
            f.write('\n' + '{:16.4e}'.format(STELLA.rhoc) +
                    '{:16.4e}'.format(STELLA.qinp) +
                    '{:16.4e}'.format(STELLA.shat) +
                    '{:16.4e}'.format(STELLA.L_reference) +
                    '{:16.4e}'.format(STELLA.B_reference) +
                    '{:16.4e}'.format(VMEC.zed_scalefac))
            f.write("\n\n\n  -----------------------------------------------------------------------------------------------------------------------------------------------")
            f.write('\n         gaga            gagp            gpgp           gds23          gds24          gbdrift         gbdrift0         cvdrift        cvdrift0')
            f.write("\n  -----------------------------------------------------------------------------------------------------------------------------------------------")
            for j in range(2*STELLA.nzgrid+1):
                for i in range(STELLA.nalpha):
                    f.write('\n' + '{:16.4e}'.format(STELLA.grad_alpha_grad_alpha[i,j]) +
                          '{:16.4e}'.format(STELLA.grad_alpha_grad_psi[i,j]) +
                          '{:16.4e}'.format(STELLA.grad_psi_grad_psi[i,j]) +
                          '{:16.4e}'.format(STELLA.gds23[i,j]) +
                          '{:16.4e}'.format(STELLA.gds24[i,j]) +
                          '{:16.4e}'.format(STELLA.gbdrift_alpha[i,j]) +
                          '{:16.4e}'.format(STELLA.gbdrift0_psi[i,j]) +
                          '{:16.4e}'.format(STELLA.cvdrift_alpha[i,j]) +
                          '{:16.4e}'.format(STELLA.cvdrift0_psi[i,j]))
    return
 
#===============================================================================
#                      vmec_to_stella_geometry_interface.f90
#===============================================================================

def vmec_to_stella_geometry_interface(KGRID, VMEC, ZGRID):
    
#*********************************************************************
#                         Get data
#*********************************************************************

    # Stella routine
    if ZGRID.verbose: print('stella::vmec_geo::vmec_to_stella_geometry_interface')
         
    # Unpack the grid sizes since they are used a lot
    ns = VMEC.ns
    nzeta = int(ZGRID.nzgrid_vmec*2+1)
    nzgrid = ZGRID.nzgrid_vmec
    nalpha = KGRID.nalpha 
    
    # We can read more nzgrid points from the VMEC than actually used in stella
    VMEC.desired_normalized_toroidal_flux = VMEC.torflux
    VMEC.initiate_arraysForVmec(KGRID.nalpha, nzgrid, nzeta) 
    
    # Do some sanity checking to ensure the VMEC arrays have some expected properties.
    VMEC.perform_sanityChecks(ns)
    
#*********************************************************************
#                         Basic Calculations
#*********************************************************************

    # Calculate basic quantities for the geometry (edge_toroidal_flux_over_2pi, 
    # sign_toroidal_flux, normalized_toroidal_flux_full_grid, normalized_toroidal_flux_half_grid)
    VMEC.calculate_toroidalFluxAtEdge(ns)
      
    # Calculate the reference units, which is independent of the chosen flux surface
    VMEC.calculate_referenceUnits()
      
    # Calculate the toroidal flux on the full grid and half grid
    VMEC.calculate_normalizedToroidalFluxOnFullAndHalfGrids(ns)

    #  Determine which flux surface to use, based on desired_normalized_toroidal_flux and vmec_surface_option.
    if VMEC.surface_option==0: VMEC.set_normalizedToroidalFluxUsed(VMEC.desired_normalized_toroidal_flux)
    if VMEC.surface_option!=0: print("WARNING: didn't implement the other surface options."); exit
  
    # Get the interpolation weights for the flux surface
    VMEC.calculate_interpolationWeights(ns)
  
#*********************************************************************
# Evaluate several radial-profile functions at the flux surface
#*********************************************************************

    # Evaluate several radial-profile functions at the flux surface
    VMEC.calculate_radialProfileFunctionsAtFluxSurface(ns)

#*********************************************************************
# Set up the coordinate grids for VMEC
#*********************************************************************
  
    # Get alpha and zeta coordinates with nfp_to_include=nfield_periods*zgrid_scalefac
    VMEC.get_alphaGridForVmec()
    VMEC.get_zetaGridForVmec()
   
    # We know theta_pest = alpha + iota*zeta, but we need to determine
    # theta_vmec = theta_pest - Lambda.  
    VMEC.get_thetaGridForVmec()
    
#*********************************************************************
# Now that we know the grid points in theta_vmec, we can evaluate
# all the geometric quantities on the grid points.
#*********************************************************************
      
    # All the quantities we need except R, Z, and Lambda use the "_nyq" mode numbers. 
    for imn_nyq in range(VMEC.mnmax_nyq):
  
        # Get the modes
        m = int(VMEC.xm_nyq[imn_nyq])
        n = int(VMEC.xn_nyq[imn_nyq]/VMEC.nfp)
          
        # Check whether the nyquist modes are available
        non_Nyquist_mode_available, scale_factor, imn = check_presenceNyquistModes(m,n,VMEC)  
        
        # Evaluate the radial derivatives we will need:
        # B and Lambda are on the half mesh, so their radial derivatives are on the full mesh.
        # R and Z are on the full mesh, so their radial derivatives are on the half mesh.
        VMEC.get_radialDerivativesOfMnModes(VMEC.ds, ns, imn, imn_nyq, non_Nyquist_mode_available)
  
        # For this <imn> add the contribution to the geometric quantities
        for ialpha in range(nalpha):
            for izeta in range(nzeta):
                VMEC.add_contributionsToMagneticQuantities(ialpha, izeta, m, n, imn_nyq, scale_factor)
                  
                # Handle arrays that use xm and xn instead of xm_nyq and xn_nyq.
                if (non_Nyquist_mode_available):
                    VMEC.add_contributionsToCoordinateQuantities(ialpha, izeta, m, n, imn, scale_factor)

#*********************************************************************
# Using R(theta,zeta) and Z(theta,zeta), compute the Cartesian
# components of the gradient basis vectors using the dual relations:
#********************************************************************* 
    VMEC.calculate_cartesianComponents() 
    VMEC.calculate_moreCartestionComponents()  
    VMEC.calculate_crossProducts()
 
#********************************************************************* 
#        Finally, assemble the quantities needed for stella. 
#********************************************************************* 
    VMEC.calculate_quantitiesForStella() 
    return


#-----------------------------
def check_presenceNyquistModes(m,n,VMEC):
     
    # Check whether the nyquist modes are available 
    if (np.abs(m) >= VMEC.mpol or np.abs(n) > VMEC.ntor):
        non_Nyquist_mode_available = False
        imn = None
         
    # Find the imn in the non-Nyquist arrays that corresponds to the same m and n.
    else:
        non_Nyquist_mode_available = True
        found_imn = False
        for imn in range(VMEC.mnmax):
            if (VMEC.xm[imn]==m and VMEC.xn[imn]==n*VMEC.nfp):
                found_imn = True
                break 
        if ((VMEC.xm[imn] != m) or (VMEC.xn[imn] != n*VMEC.nfp)):
            print("Something went wrong#"); return
        if (not found_imn):
            print("Error# imn could not be found matching the given imn_nyq."); return

    # All quantities are multiplied by a variable scale_factor which can in principle depend on m and n.
    # For now we just set scale_factor = 1. In the future, scale_factor could be used to lower the
    # symmetry-breaking Fourier components, or filter out certain Fourier components in some way.
    scale_factor = 1
        
    # Return whether the non nyquist modes were available 
    return non_Nyquist_mode_available, scale_factor, imn

#===============================================================================
#                                 dist_fn.f90
#=============================================================================== 
# We have the following geometric quantities in stella:
#     gds2 = |nabla y|**2
#     gds21 = s_hat * (nabla x cdot nabla y)
#     gds22 = s_hat**2 * |nabla x|**2
#     theta0 = 0            for ky=0
#     theta0 = kx/ky/shat   for ky =#= 0
#     theta0 = (shat*kx)/ky for ky =#= 0
#
# The perpendicular wavenumber squared is defined as:
#     kperp2 = k_x**2 * |nabla x|**2 + k_y**2 * |nabla y|**2 + 2*k_x*k_y*(nabla x cdot nabla y)
#            = k_y**2 * theta0**2 * gds22 + k_y**2 * gds2 + 2*k_y**2*theta0*gds21 )
#            = k_y**2 ( theta0**2 *gds22 + gds2 + 2*theta0*gds21 )
#            = aky**2 *(gds2 + 2.0*theta0*gds21 + theta0*theta0*gds22)

def init_kperp2(KGRID, STELLA, ZGRID):

    # Stella routine
    if STELLA.verbose: print("stella::dist_fn::init_kperp2")
    
    # Make sure the KGRID is initiated
    if KGRID.kt_grids_knobs=="box":     KGRID.init_kt_grids_box(STELLA, ZGRID)
    if KGRID.kt_grids_knobs=="range":   KGRID.init_kt_grids_range(STELLA, ZGRID)

    # Allocate the arrays <kperp2> and <dkperp2dr>
    STELLA.kperp2 = np.zeros((KGRID.naky,KGRID.nakx,KGRID.nalpha,2*ZGRID.nzgrid+1))
    STELLA.dkperp2dr = np.zeros((KGRID.naky,KGRID.nakx,KGRID.nalpha,2*ZGRID.nzgrid+1))

    for iky in range(KGRID.naky): 

        # Calculate kperp2 for a zonal mode: if aky(iky)==0
        if (KGRID.zonal_mode[iky]):
            for ikx in range(KGRID.nakx): 

                # If radial_variation then q_as_x=true otherwise q_as_x==false by default
                if (KGRID.q_as_x):
                    STELLA.kperp2[iky,ikx,:,:] = KGRID.akx[ikx]*KGRID.akx[ikx]*STELLA.gds22
                    STELLA.dkperp2dr[iky,ikx,:,:] = KGRID.akx[ikx]*KGRID.akx[ikx]*STELLA.dgds22dr/STELLA.kperp2[iky,ikx,:,:]
                else:
                    # kperp2 = k_x**2 * |nabla x|**2 
                    STELLA.kperp2[iky,ikx,:,:] = KGRID.akx[ikx]*KGRID.akx[ikx]*STELLA.gds22/(STELLA.shat**2)
                    # dkperp2/dr = (n0/a)^2*(drho/dpsiN)^2*(theta0^2*dgds22dr)
                    STELLA.dkperp2dr[iky,ikx,:,:] = KGRID.akx[ikx]*KGRID.akx[ikx]*STELLA.dgds22dr/(STELLA.shat**2*STELLA.kperp2[iky,ikx,:,:]) 

                # Make sure dkperp2dr=0 when kperp2=0
                if (np.any(STELLA.kperp2[iky,ikx,:,:] < sys.float_info.epsilon)): STELLA.dkperp2dr[iky,ikx,:,:] = 0. 

        # Calculate kperp2 for all modes that are not zonal modes
        else:
            for ikx in range(KGRID.nakx): 
                
                # kperp2 = k_x**2 * |nabla x|**2 + k_y**2 * |nabla y|**2 + 2*k_x*k_y*(nabla x cdot nabla y)
                STELLA.kperp2[iky,ikx,:,:] = KGRID.aky[iky]*KGRID.aky[iky] \
                     *(STELLA.gds2 + 2.0*KGRID.theta0[iky,ikx]*STELLA.gds21 \
                     + KGRID.theta0[iky,ikx]*KGRID.theta0[iky,ikx]*STELLA.gds22)
                
                # dkperp2/dr = (n0/a)^2*(drho/dpsiN)^2*(dgds2dr + 2*theta0*dgds21dr + theta0^2*dgds22dr)
                # note that dkperp2dr can not be calculated when reading from a VMEW file
                STELLA.dkperp2dr[iky,ikx,:,:] = KGRID.aky[iky]*KGRID.aky[iky] \
                     *(STELLA.dgds2dr + 2.0*KGRID.theta0[iky,ikx]*STELLA.dgds21dr \
                     + KGRID.theta0[iky,ikx]*KGRID.theta0[iky,ikx]*STELLA.dgds22dr)
                STELLA.dkperp2dr[iky,ikx,:,:] = STELLA.dkperp2dr[iky,ikx,:,:]/STELLA.kperp2[iky,ikx,:,:] 
                
                # Make sure dkperp2dr=0 when kperp2=0
                if (np.any(STELLA.kperp2[iky,ikx,:,:] < sys.float_info.epsilon)): STELLA.dkperp2dr[iky,ikx,:,:] = 0.
    
    # Enforce periodicity: kperp2(-z) = kperp2(z)
    enforce_single_valued_kperp2(KGRID, STELLA, ZGRID)



#------------------------------  
def enforce_single_valued_kperp2(KGRID, STELLA, ZGRID):
    ''' Enforce periodicity: kperp2(-z) = kperp2(z) '''
    
    # Stella routine
    if STELLA.verbose: print("stella::dist_fn::enforce_single_valued_kperp2")
    
    # Make sure <nsegments> and <neigen> are defined
    ZGRID.init_extended_zgrid(KGRID) 

    tmp = np.zeros((KGRID.nalpha)) 
    for iky in range(KGRID.naky):
        for ie in range(int(ZGRID.neigen[iky])):
            if (ZGRID.nsegments[ie,iky] > 1):
                for iseg in np.arange(2, ZGRID.nsegments[ie,iky]):
                    tmp[:] = 0.5*(STELLA.kperp2[iky,ZGRID.ikxmod[iseg-1,ie,iky],:,2*ZGRID.nzgrid+1] + STELLA.kperp2[iky,ZGRID.ikxmod[iseg,ie,iky],:,0])
                    STELLA.kperp2[iky,ZGRID.ikxmod[iseg,ie,iky],:,0] = tmp
                    STELLA.kperp2[iky,ZGRID.ikxmod[iseg-1,ie,iky],:,2*ZGRID.nzgrid+1] = tmp 
    return 
