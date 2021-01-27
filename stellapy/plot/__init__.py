#=====================================================================
# Plot data from stella or the self made files
#=====================================================================
''' Plot stella data

BASH GUIDE
------

The netcdf file is needed for some plots, instead use this function on marconi and download the smaller out.h5 file
    -->   reduce_size_netcdf       
Download the small and important files locally for faster plotting
    -->   sync_runs_marconi   
Plot the desired quantities
    -->   plot_omega_vs_t
    -->   plot_phi_vs_z



Stella output files
-------------------

plot.omega_vs_t         -->     (omega,t)          {*.omega; *.out.nc or *out.h5;}
plot.phi_vs_z           -->     (phi,z)            {*.out.nc or *out.h5}




Self-made output files: 
----------------------

plot.



Other functions
---------------


plotbox                 -->     Holds layouts for specific plots

    plotbox.plotbox_2d  -->     Figure kayour for a (x,y) plot


'''


#=====================================================================
# Make function avaliable as package.func instead of package.mod.func
#=====================================================================

# Import the current modules and sub-package list in the package folder
import os, glob

divider = '\\' if (os.name == 'nt') else '/'
mod_list = [file_name.split(divider)[-1].split('.')[0] for file_name in glob.glob(__path__[0]+'/[!_]*.py')]
sub_pack_list = [folder_name.split(divider)[-2] for folder_name in glob.glob(__path__[0]+'/[!_]*/')]

# Import all functions from the modules
for mod in mod_list:
    exec('from . import ' + mod)
    exec('from .' + mod + ' import *')

# Import all subpackages
for pack in sub_pack_list:
    exec('from . import ' + pack)

# Clean up
del glob
try:
    del mod
except:
    pass
try:
    del pack
except:
    pass
