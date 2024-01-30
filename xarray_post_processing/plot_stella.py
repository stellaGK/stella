# -*- coding: utf-8 -*-
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg') # this line allows plots to be made without using a display environment variable
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import xarray as xr
import argparse
import os

from utils import dimension_sizes
from utils import flux_data
from utils import field_data
from utils import species_data
from utils import plot_fluxes_nc
from utils import plot_fields_nc
from utils import plot_fields_fluxes_spectra_nc
from utils import plot_fields_kyspectra_with_time_nc

# define the command line inputs
parser = argparse.ArgumentParser()
ncfilehelpstr = "the file_name of the stella netCDF4 file_name.out.nc file to be analysed"
parser.add_argument("ncfile", help=ncfilehelpstr)
# read the command line inputs
args = parser.parse_args()
# get the current working directory to construct an absolute path
workdir = os.path.abspath(os.getcwd())

filename = workdir + "/" + args.ncfile

stelladata = xr.open_dataset(filename+".out.nc")

nky, nkx, nzed, nvpa, nmu, nspec, ntubes, nalpha = dimension_sizes(stelladata)
pflx, vflx, qflx, time = flux_data(stelladata)
phi2, apar2, apar2_present, bpar2, bpar2_present = field_data(stelladata)
charge, mass, dens, temp, tprim, fprim, vnew, typeint, typestring = species_data(stelladata,filename)
stelladata.close()

# plot the heat, particle and momentum fluxes
plot_fluxes_nc(filename, pflx, vflx, qflx, time, typestring)
# plot the fields
plot_fields_nc(filename, phi2, apar2, apar2_present, bpar2, bpar2_present, time)
# plot the spectra
plot_fields_fluxes_spectra_nc(filename)
plot_fields_kyspectra_with_time_nc(filename)
