# -*- coding: utf-8 -*-

## Description: import variables from stella netcdf file

from scipy.io import netcdf
import numpy as np

####### Import variables from netcdf file #########
#infile = input("Path to netcdf file: ")
infile = '/Users/michaelbarnes/codes/stella/runs/jcp_paper_sims/NCSX_lin_benchmark/ky1p0.out.nc'
#infile = '../stella.out.nc'
print()
#outdir = input("Path for output: ")
#outdir = '/Users/michaelbarnes/codes/gs2/runs/flow_shear_tests/nofs01_figs/'
outdir = '/Users/michaelbarnes/codes/stella/runs/jcp_paper_sims/NCSX_lin_benchmark/figures/'
ncfile = netcdf.netcdf_file(infile,'r')   

print()
print('reading data from netcdf file...')
print()

# get kx and ky grids
kx_stella = np.copy(ncfile.variables['kx'][:])
ky_stella = np.copy(ncfile.variables['ky'][:])

# get zed grid
zed = np.copy(ncfile.variables['zed'][:])

# get time grid
time_stella = np.copy(ncfile.variables['t'][:])

# number of kinetic species
nspec = ncfile.dimensions['species']

# get geometric quantities
bmag = np.copy(ncfile.variables['bmag'][:])
gradpar = np.copy(ncfile.variables['gradpar'][:])
gbdrift = np.copy(ncfile.variables['gbdrift'][:])
gbdrift0 = np.copy(ncfile.variables['gbdrift0'][:])
cvdrift = np.copy(ncfile.variables['cvdrift'][:])
cvdrift0 = np.copy(ncfile.variables['cvdrift0'][:])
gds2 = np.copy(ncfile.variables['gds2'][:])
gds21 = np.copy(ncfile.variables['gds21'][:])
gds22 = np.copy(ncfile.variables['gds22'][:])

def read_stella_float(var):

    import numpy as np

    try:
        arr = np.copy(ncfile.variables[var][:])
        flag = True
    except KeyError:
        print('INFO: '+var+' not found in netcdf file')
        arr = np.arange(1,dtype=float)
        flag = False

    return arr, flag

# electrostatic potential as a function of (z,kx,ky,t)
phi_vs_t, phi_vs_t_present \
    = read_stella_float('phi_vs_t')
# |g|^2 averaged over kx, ky, and z
gvmus, gvmus_present \
    = read_stella_float('gvmus')
# |g|^2 averaged over kx, ky, and mu
gzvs, gzvs_present \
    = read_stella_float('gzvs')
# jacobian for transformation to (rho,alpha,z) coordinates
jacob, jacob_present \
    = read_stella_float('jacob')
# gradient of normalized radial coordinate rho
grho, grho_present \
    = read_stella_float('grho')
# modulus squared of electrostatic potential (averaged over space)
phi2_stella, phi2_stella_present \
    = read_stella_float('phi2')
# time-dependent electrostatic particle flux for each species
es_part_flux, es_part_flux_present \
    = read_stella_float('es_part_flux')
# electrostatic heat flux
es_heat_flux, es_heat_flux_present \
    = read_stella_float('es_heat_flux')
# electrostatic momentum flux
es_mom_flux, es_mom_flux_present \
    = read_stella_float('es_mom_flux')
# turbulent energy exchange
es_energy_exchange, es_energy_exchange_present = \
    read_stella_float('es_energy_exchange')
# time-dependent particle flux for each species as a function of (kx,ky)
es_part_by_k, es_part_by_k_present = \
    read_stella_float('es_part_by_k')
if es_part_by_k_present is not True:
    es_part_by_k, es_part_by_k_present = \
        read_stella_float('es_part_flux_by_mode')
# time-dependent heat flux for each species as a function of (kx,ky)
es_heat_by_k, es_heat_by_k_present = \
    read_stella_float('es_heat_by_k')
if es_heat_by_k_present is not True:
    es_heat_by_k, es_heat_by_k_present = \
        read_stella_float('es_heat_flux_by_mode')
# time-dependent momentum flux for each species as a function of (kx,ky)
es_mom_by_k, es_mom_by_k_present = \
    read_stella_float('es_mom_by_k')
if es_mom_by_k_present is not True:
    es_mom_by_k, es_mom_by_k_present = \
        read_stella_float('es_mom_flux_by_mode')
es_energy_exchange_by_k, es_energy_exchange_by_k_present = \
    read_stella_float('es_energy_exchange_by_k')
if es_energy_exchange_by_k_present is not True:
    es_energy_exchange_by_k, es_energy_exchange_by_k_present = \
        read_stella_float('es_energy_exchange_by_mode')
    
es_energy_exchange_by_ky, es_energy_exchange_by_ky_present = \
    read_stella_float('es_energy_exchange_by_ky')
# modulus squared of electrostatic potential as function of (kx,ky,t)
phi2_by_mode, phi2_by_mode_present = \
    read_stella_float('phi2_by_mode')
# modulus squared of electrostatic potential as function of (ky,t)
phi2_by_ky, phi2_by_ky_present = \
    read_stella_float('phi2_by_ky')
# parallel velocity grid
vpa, vpa_present = \
    read_stella_float('vpa')
# mu grid
mu, mu_present = \
    read_stella_float('mu')
# electrostatic particle flux as function of (vpa,z)
es_part_sym, es_part_sym_present = \
    read_stella_float('es_part_sym')
# electrostatic heat flux as function of (vpa,z)
es_heat_sym, es_heat_sym_present = \
    read_stella_float('es_heat_sym')
# electrostatic momentum flux as function of (vpa,z)
es_mom_sym, es_mom_sym_present = \
    read_stella_float('es_mom_sym')
# if vpa not present then do not try to plot es_mom_sym
if vpa_present == False:
    es_mom_sym_present = False
# phi(kx,ky,t) for given z (usually outboard midplane)
# note that name of phi_igomega_by_mode used to be phi0
# so must allows for this case
phi_igomega_by_mode, phi_igomega_by_mode_present  = \
    read_stella_float('phi_igomega_by_mode')
if phi_igomega_by_mode_present == False:
    phi_igomega_by_mode, phi_igomega_by_mode_present  = \
        read_stella_float('phi0')
phip0, phip0_present = \
    read_stella_float('phip0')
xgrid, xgrid_present = \
    read_stella_float('xgrid')
xgrid = np.concatenate((xgrid[kx_stella.shape[0]//2+1:],xgrid[:kx_stella.shape[0]//2+1]))
# density fluctuation (kx,ky,t) for given z (usually outboard midplane)
ntot_igomega_by_mode, ntot_igomega_by_mode_present = \
    read_stella_float('ntot_igomega_by_mode')
# parallel flow fluctuation (kx,ky,t) for given z (usually outboard midplane)
upar_igomega_by_mode, upar_igomega_by_mode_present = \
    read_stella_float('upar_igomega_by_mode')
# parallel temperature fluctuation (kx,ky,t) for given z (usually outboard midplane)
tpar_igomega_by_mode, tpar_igomega_by_mode_present = \
    read_stella_float('tpar_igomega_by_mode')
# perp termperature fluctuation (kx,ky,t) for given z (usually outboard midplane)
tperp_igomega_by_mode, tperp_igomega_by_mode_present = \
    read_stella_float('tperp_igomega_by_mode')

ncfile.close()

print()
print('...finished reading data from netcdf file')
