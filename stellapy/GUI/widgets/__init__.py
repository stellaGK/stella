''' 
#===============================================================================      
###########################  STELLA GUI: WIDGETS ###############################
#===============================================================================
Each tab reuses certain complex widgets which are constructed in this module.


'''

# Global variables
PAD_TITLE = {}
PAD_LABEL = {}
PAD_ENTRY = {}
PAD_TITLE2 = {}
PAD_LABEL2 = {}
PAD_ENTRY2 = {}
PAD_DIVID2 = {}
PAD_ENTRYR = {}

#===============================================================================
# Make the functions avaliable as package.func instead of package.mod.func
#===============================================================================

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
