# -*- coding: utf-8 -*-

## Description: import variables from stella netcdf file

from scipy.io import netcdf
import numpy as np
import netCDF4 as nc4

####### Import variables from netcdf file #########
#infile = input("Path to netcdf file: ")
#input_directory = '/Users/barnesm/codes/stella/runs/cbc_example/SL_test/'
input_directory = '/Users/barnesm/codes/stella/runs/zpinch_om/ExB/'
file_prefix = 'LBLT200_2d_kymax10'
infile = input_directory + file_prefix + '.out.nc'
#infile = '../stella.out.nc'
#print()
#outdir = input("Path for output: ")
outdir = input_directory
ncfile = nc4.Dataset(infile)

print('reading data from netcdf file...')
print()

# get kx and ky grids
kx_stella = np.copy(ncfile.variables['kx'][:])
nakx = len(ncfile.dimensions['kx'])
ky = np.copy(ncfile.variables['ky'][:])
naky = len(ncfile.dimensions['ky'])

# this is the index of the first negative value of kx
# note stella orders kx as (0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx)
nakx_mid = nakx//2+1
kx = np.concatenate((kx_stella[nakx_mid:],kx_stella[:nakx_mid]))

# get zed grid
zed = np.copy(ncfile.variables['zed'][:])
nzed = zed.size
delzed = zed[1]-zed[0]
iz0 = nzed//2+1

# get time grid
time = np.copy(ncfile.variables['t'][:])
ntime = time.size

# number of kinetic species
nspec = len(ncfile.dimensions['species'])

# get geometric quantities
bmag = np.copy(ncfile.variables['bmag'][:])
gradpar = np.copy(ncfile.variables['gradpar'][:])
gbdrift = np.copy(ncfile.variables['gbdrift'][:])
gbdrift0 = np.copy(ncfile.variables['gbdrift0'][:])
cvdrift = np.copy(ncfile.variables['cvdrift'][:])
cvdrift0 = np.copy(ncfile.variables['cvdrift0'][:])
#gds2 = np.copy(ncfile.variables['gds2'][:])
#gds21 = np.copy(ncfile.variables['gds21'][:])
#gds22 = np.copy(ncfile.variables['gds22'][:])

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

# electrostatic potential averaged over z as function of (ky,kx,t)
phi2_vs_kxky_stella, phi2_vs_kxky_present \
    = read_stella_float('phi2_vs_kxky')
if (phi2_vs_kxky_present):
    phi2_vs_kxky_stella[:,0,0] = 0.0
    phi2_vs_kxky = np.concatenate((phi2_vs_kxky_stella[:,nakx_mid:,:],phi2_vs_kxky_stella[:,:nakx_mid,:]),axis=1)
# electrostatic potential as a function of (z,kx,ky,t)
phi_vs_t, phi_vs_t_present \
    = read_stella_float('phi_vs_t')
# |g|^2 averaged over kx, ky, and z
gvmus, gvmus_present \
    = read_stella_float('gvmus')
# |g|^2 averaged over kx, ky, and mu
gzvs, gzvs_present \
    = read_stella_float('g2_vs_zvpas')
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
    es_heat_by_k_stella, es_heat_by_k_present = \
        read_stella_float('qflux_vs_kxkys')
    es_heat_by_k_stella[:,:,0,0] = 0.0
    es_heat_by_k = np.concatenate((es_heat_by_k_stella[:,:,nakx_mid:,:],es_heat_by_k_stella[:,:,:nakx_mid,:]),axis=2)

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
xgrid, xgrid_present = \
    read_stella_float('xgrid')
if (xgrid_present):
    xgrid = np.concatenate((xgrid[kx_stella.shape[0]//2+1:],xgrid[:kx_stella.shape[0]//2+1]))
else:
    xgrid = np.arange(nakx,dtype=float)
    if (kx_stella.shape[0] > 1):
        dx = 2*np.pi/(kx_stella[1]*nakx)
    else:
        if (kx_stella[0] > 0.0):
            dx = 2*np.pi/kx_stella[0]
        else:
            dx = 0.0

    for ikx in range(nakx):
        xgrid[ikx] = ikx*dx
# density fluctuation (kx,ky,z,t)
density, density_present = \
    read_stella_float('density')
# parallel flow fluctuation (kx,ky,z,t)
upar, upar_present = \
    read_stella_float('upar')
# temperature fluctuation (kx,ky,z,t)
temperature, temperature_present = \
    read_stella_float('temperature')

ncfile.close()

print()
print('...finished reading data from netcdf file')
