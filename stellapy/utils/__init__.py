###################################################################
# Functions to read data from stella output files
###################################################################
''' 
Functions to read data from stella output files. 

Add the following line to the .alias file of your PC:
	export RUNS='/home/hanne/Documents/CIEMAT/RUNS
			
The hierarchy is such that $RUNS/$CODE/$EQUIL_NAME'_'$NUM 
contains the input and output of the run.
Example : $RUNS/stella/w7xr003_0125

Methods
-------
runsdir() : str
	Returns the path name where the cases or runs are located
outdir(case=None) : str
	Returns the path name of a specific [case] folder.
outfile(case=None, quant=None) : str
	Returns the full path of an output file, with extension "quant"
geotxtfile(case=None) : str
	Returns the full path of the *vmec_geo file
infile(case=None) : str
	Returns the full path of the netcdf file *out.nc

Parameters
----------
case : string
    The name of the case folder [$EQUIL_NAME'_'$NUM] plus input file.
	Example: 'w7xr169+252_0001/kyscan_s17_B169_rho0.5_AE_ky_1.4.in'

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
