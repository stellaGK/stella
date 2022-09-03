################################################################################
#                       DISTRIBUTION OF THE GUIDING CENTERS
################################################################################
# 
# In order to use this data perform this command on the supercomputer:
#     write_h5FileForDistribution 
#
# The size can be increased to include more data
#     write_h5FileForDistribution -s 0     --> g_vs_smu, g_vs_svpa, g_vs_ts
#     write_h5FileForDistribution -s 1     --> g_vs_tsmu, g_vs_tsvpa
#     write_h5FileForDistribution -s 2     --> g_vs_tsmuvpa, g_vs_tsvpaz
# The default size if <0> for linear and <1> for nonlinear simulations
# 
# The default time step is dt=10, and it can be changed through:
#     write_h5FileForDistribution -t 1
#
################################################################################

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
