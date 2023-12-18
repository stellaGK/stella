""" 

#===============================================================================
#                   READ THE CONFIGURATION FILE STELLAPY.INI                   #
#=============================================================================== 

For some functionalities of stellapy, it is required to set some paths/values, 
which have been defined in the '$STELLAPY/stellapy.ini' configuration file.
The default configuration file is defined in 'create_defaultConfigurationFile.py'.
    
Hanne Thienpondt
23/01/2023

"""

#!/usr/bin/python3
import os, sys
import pathlib
import configparser 

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.utils.decorators.exit_program import exit_program 

# Grab the location of the configuration file 
divider = '\\' if (os.name == 'nt') else '/'  
PATH_CONFIGURATIONFILE = pathlib.Path(os.environ.get('STELLAPY')+"/stellapy.ini")  
PATH_SOURCEFILE = pathlib.Path(os.environ.get('STELLAPY')+"/source.sh")

# Create a configuration object 
CONFIG = configparser.ConfigParser(allow_no_value=True)                                                   
CONFIG.optionxform = lambda option: option  # preserve case for letters  
CONFIG.read(PATH_CONFIGURATIONFILE)   

#===============================================================================
#                   READ THE CONFIGURATION FILE STELLAPY.INI                   #
#===============================================================================

def read_configurationFile():
    CONFIG.read(PATH_CONFIGURATIONFILE)
    return

#------------------------------    
def write_configurationFile(): 
    with open(PATH_CONFIGURATIONFILE, 'w') as configfile:
        CONFIG.write(configfile)
    return

#------------------------------      
def check_pathsToCodeAndSimulations(): 
    
    # Read the configuration file
    read_configurationFile()
        
    # Check whether the VMEC path is set correctly
    if 'MANDATORY' in CONFIG['PATHS']['VMEC']:
        exit_reason = "A default configuration file has been created, please fill in [PATHS][VMEC]:\n"
        exit_reason += "    "+os.environ.get('STELLAPY')+"stellapy.ini \n"
        exit_program(exit_reason, read_configurationFile, sys._getframe().f_lineno)  
    return

