###################################################################
# Functions to read data from stella output files
###################################################################
''' 

Create and manage the widgets for the first tab of the GUI. Create the three 
main classes and frames for the simulations (see TabSelectedFiles) 
    1. Select simulations and display them through class_simulations in 
       frame_simulations. (see TabSelectedFiles_Simulations)
    2. Show the progress of the computations after simulations have been selected 
       through class_progress in frame_progress. (see TabSelectedFiles_InputParemeters)
    3. Display the input parameters related to the selected simulations through 
       class_inputParameters in 6 frames. (see TabSelectedFiles_Progress)


'''


#=====================================================================
# Make function avaliable as package.func instead of package.mod.func
#=====================================================================

# Import the current modules and sub-package list in the package folder
import os, glob

_divider = '\\' if (os.name == 'nt') else '/'
_mod_list = [file_name.split(_divider)[-1].split('.')[0] for file_name in glob.glob(__path__[0]+'/[!_]*.py')]
_sub_pack_list = [folder_name.split(_divider)[-2] for folder_name in glob.glob(__path__[0]+'/[!_]*/')]

# Import all functions from the modules
for mod in _mod_list:
    exec('from . import ' + mod)
    exec('from .' + mod + ' import *')

# Import all subpackages
for pack in _sub_pack_list:
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
