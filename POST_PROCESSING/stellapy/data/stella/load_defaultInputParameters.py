'''

Load default input parameters of the stella code, by reading the 
'input_stella_v1.0.in' file in the 'stellapy/data/stella' folder.

Returns
-------
input_parameters = dict[knobs][variable]
    Returns the namelists from the stella code with their default values.

Hanne Thienpondt
03/10/2025
 
'''

#!/usr/bin/python3
import sys, os
import pathlib
import f90nml
import json

#===============================================================================
#                   Load default input parameters for stella                    
#===============================================================================
def load_defaultInputParameters():
        
    # Read the default input file
    stellapy_folder = pathlib.Path(os.path.abspath(__file__).split('/stellapy/')[0]) / 'stellapy'
    default_input_file = stellapy_folder / 'data/stella/input_stella_v1.0.in'
    
    # Read default input parameters
    input_parameters = f90nml.read(default_input_file)
    input_parameters = json.loads(json.dumps(input_parameters))
    
    # Extra variables that will be set throughout stellapy
    input_parameters['z_grid']['nz'] = -1
    input_parameters['geometry_vmec']['rho'] = -1.0
    input_parameters['geometry_vmec']['poloidal_turns'] = -1.0
    input_parameters['kxky_grid_box']['lx'] = -1.0
    input_parameters['kxky_grid_box']['ly'] = -1.0
    input_parameters['kxky_grid_box']['Lx'] = -1.0
    input_parameters['kxky_grid_box']['Ly'] = -1.0
    input_parameters['kxky_grid_box']['dkx'] = -1.0
    input_parameters['kxky_grid_box']['dky'] = -1.0
    input_parameters['kxky_grid_box']['kx max'] = -1.0
    input_parameters['kxky_grid_box']['ky max'] = -1.0
    input_parameters['kxky_grid_box']['shat'] = -1.0
    input_parameters['velocity_grids']['dmu'] = -1.0
    input_parameters['velocity_grids']['dvpa'] = -1.0
    input_parameters['velocity_grids']['vpa max'] = -1.0
    input_parameters['velocity_grids']['mu max'] = -1.0
    input_parameters['adiabatic_electron_response']['teti'] = 1/input_parameters['adiabatic_electron_response']['tite']
    return input_parameters

