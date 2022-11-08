
import os

#===============================================================================
#                 CREATE A DEFAULT CONFIGURATION FILE STELLAPY.INI
#===============================================================================

def create_defaultConfigurationFile(CONFIG, NAME_CONFIGURATIONFILE):
    ''' Write a default configuration file and save it as "stellapy.utils.config/config.ini". '''
    
    # Get the path to the current file
    path_currentFile = os.path.dirname(os.path.abspath(__file__))
    
    # Default section
    CONFIG['GENERAL'] = { 
        '# General information about the stellapy package.' : None,\
        'Package'   : 'Stellapy',\
        'Author'    : 'Hanne Thienpondt',\
        'Email'     : 'Hanne.Thienpondt@outlook.com',\
        'Date'      : '01/09/2022',\
        'Version'   : '5.0'}  
    
    CONFIG['CODE'] = { 
        '# Save the location of the stella and stellapy packages. ' : None,\
        'Stella'    : path_currentFile.split("stellapy")[0],\
        'Stellapy'  : path_currentFile.split("stellapy")[0]+"stellapy/"} 
    
    CONFIG['PATHS'] = {
        '# VMEC: Location of VMEC files ("wout*.nc").' : None,\
        '# Simulations: Simulations that are downloaded from the supercomputer. ' : None,\
        '# GUI: The pickly files and figures will be saved here.' : None,\
        '# User: Name of the user with a corresponding stellapy/<user> folder.' : None,\
        'VMEC' : '/home/user/STELLA/VmecFiles/',\
        'Runs': '/home/user/STELLA/RUNS/',\
        'GUI' : '/home/user/STELLA/GUI/',\
        'User' : 'User'} 

    CONFIG['GUI SETTINGS'] = {
        '# Configure the GUI.' : None,\
        'Theme' : 'awlight',\
        'TextEditor' : 'emacs'}
    
    CONFIG['SUPERCOMPUTER'] = {
        '# Define the email for the run updates.' : None,\
        'Email' : 'dummy_email@gmail.com'} 
    
    CONFIG['MARCONI'] = {
        '# Define the default account, partition and number of nodes used for the simulations.' : None,\
        'VMEC' : '/marconi_scratch/userexternal/user/SCRATCH/EQUIL/',\
        'Account'   : 'FUA36_TRASIMEX',\
        'Partition' : 'skl_fua_prod',\
        'Nonlinear Nodes' : 4,\
        'Linear Nodes' : 1,\
        'UserName' : 'user'} 

    # Write the configuration file
    CONFIG.write(open(NAME_CONFIGURATIONFILE, 'w'))    
    return
