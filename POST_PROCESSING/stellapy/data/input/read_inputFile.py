''' 

#===============================================================================
#                              Read the input file                             #
#===============================================================================

Read the input file located at <path_input_file>. Additionally, calculate some 
of the interal stella variables and some extra parameters such as {Lx, Ly}.
    
Hanne Thienpondt
03/10/2025

'''

#!/usr/bin/python3 
import h5py
import numpy as np
import configparser
import json, f90nml
import os, sys, pathlib

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)
from stellapy.data.stella.load_defaultInputParameters import load_defaultInputParameters
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables
from stellapy.utils.decorators.exit_program import exit_program
    
#===============================================================================
#                             Read the input file                              #
#===============================================================================

def read_inputFile(path_input_file):
    ''' Read the input file and overwrite the default stella parameters.'''
            
    # Read the input file and overwrite the default stella parameters.
    if os.path.isfile(path_input_file):
        if path_input_file.suffix=='.in':  input_parameters = read_inFile(path_input_file)
        if path_input_file.suffix=='.ini': input_parameters = read_iniFile(path_input_file)
        input_parameters['adiabatic_electron_response']['teti'] = 1/input_parameters['adiabatic_electron_response']['tite']
        return input_parameters
    if os.path.isfile(path_input_file.with_suffix('.in')):
        input_parameters = read_inFile(path_input_file.with_suffix('.in'))
        input_parameters['adiabatic_electron_response']['teti'] = 1/input_parameters['adiabatic_electron_response']['tite']
        return input_parameters

    # Critical error if we didn't find any data
    exit_reason = 'The input data can not be found. The following file does not exist:\n'
    exit_reason += '    '+str(path_input_file)+'\n\n'
    exit_reason += 'You can try to write it through the command:\n'
    exit_reason += '    write_dataFiles -s ini \n\n'
    exit_reason += 'If the data was written before October 2022 try:\n'
    exit_reason += '    write_listOfMatchingInputFilesForOldStellapy \n'
    exit_program(exit_reason, read_inputFile, sys._getframe().f_lineno)
    return

#------------------------- 
def read_inFile(path_input_file):
    
    # Read the input namelist and parse it to a nested dict
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))

    # Initiate the dictionary: load the default stella parameters
    input_parameters = load_defaultInputParameters()
    input_parameters['species_options'].update(input_parameters_read['species_options'])

    # Add more default species if nspec>2
    for i in range(2, input_parameters['species_options']['nspec']+1):
        input_parameters['species_parameters_'+str(i)] = input_parameters['species_parameters_1'].copy()

    # Overwrite the default values if they have been changed in the <path_input_file>
    for knob in input_parameters_read.keys():
        input_parameters[knob].update(input_parameters_read[knob])

    # Close the input file
    return input_parameters

#-------------------------
def read_iniFile(path_input_file):
    
    # Read the '*.ini' file
    input_parameters_read = configparser.ConfigParser()
    input_parameters_read.read(path_input_file)
    input_parameters_read = {s:dict(input_parameters_read.items(s)) for s in input_parameters_read.sections()}

    # Initiate the dictionary: load the default stella parameters
    input_parameters = load_defaultInputParameters()
    
    # Parse the ini data strings to int/float/bool
    for knob in input_parameters_read.keys():
        for parameter in input_parameters_read[knob].keys():
            try: input_parameters_read[knob][parameter] = read_parameterFromIniFile(input_parameters_read[knob], parameter, input_parameters[knob][parameter])
            except:
                input_parameters_read[knob][parameter] = read_parameterFromIniFile(input_parameters_read[knob], parameter, input_parameters[knob][parameter])
                    
    # Add more default species if nspec>2
    if 'species_options' in input_parameters_read: input_parameters['species_options'].update(input_parameters_read['species_options'])
    for i in range(2, int(input_parameters['species_options']['nspec'])+1): 
        input_parameters['species_parameters_'+str(i)] = input_parameters['species_parameters_1'].copy() 
        
    # Overwrite the default values if they have been changed in the <path_input_file> 
    for knob in input_parameters_read.keys():
        if knob in input_parameters.keys():
            input_parameters[knob].update(input_parameters_read[knob])

    # Return the input_parameters
    return input_parameters

#===============================================================================
#            Read the input parameters and calculate some extra ones           #
#===============================================================================

def read_inputParameters(path_input_file):
    ''' Read '*.in' file and return dict[knobs][variable].'''
    
    # Read the input file and overwrite the default stella parameters
    input_parameters = read_inputFile(path_input_file)
    
    # Read indirect input parameters
    input_parameters = read_extraInputParameters(input_parameters, path_input_file)
    
    # Return the input_parameters
    return input_parameters 

#===============================================================================
#                   Calculate the internal stella variables                    #
#===============================================================================

def read_extraInputParameters(input_parameters, path_input_file):
    
    # This is only needed when we're reading the input file
    if path_input_file.suffix=='.in':
    
        # Calculate the input variables that stella calculates internally
        input_parameters = calculate_stellaVariables(input_parameters)
        
        # Calculate some extra input variables
        input_parameters = calculate_extraInputParameters(input_parameters)
        
        # Calculate some extra parameters from the VMEC file
        input_parameters = calculate_extraInputParametersFromWout(input_parameters, path_input_file) 
        input_parameters = calculate_extraInputParametersFromNetcdf(input_parameters, path_input_file)
        
    return input_parameters
    
#------------------------------------------------
def calculate_stellaVariables(input_parameters):

    # Calculate <nzgrid>
    if input_parameters['z_grid']['nzgrid'] == 'nzed/2 + (nperiod-1)*nzed':
        nzed = input_parameters['z_grid']['nzed']
        nperiod = input_parameters['z_grid']['nperiod']
        input_parameters['z_grid']['nzgrid'] = nzed/2 + (nperiod-1)*nzed
    else: print('WARNING: <nzgrid> was set by the input file but it should be calculated indirectly through <nzed> and <nperiod>.')

    # Caclulate <nakx> and <naky> 
    if input_parameters['kxky_grid_box']['naky'] == '(ny-1)/3 + 1':
        from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_nakx, calculate_naky
        ny = input_parameters['kxky_grid_box']['ny']
        nx = input_parameters['kxky_grid_box']['nx']
        input_parameters['kxky_grid_box']['naky'] = calculate_naky(ny)
        input_parameters['kxky_grid_box']['nakx'] = calculate_nakx(nx)

    # Calculate <y0> for a full flux surface simulation
    if input_parameters['gyrokinetic_terms']['include_nonlinear'] == True:
        if input_parameters['kxky_grid_box']['y0'] == -1.0:
            if input_parameters['gyrokinetic_terms']['include_full_flux_annulus'] == False:
                print('WARNING: When simulating a flux tube, y0 needs to be set in the input file.')
            if input_parameters['gyrokinetic_terms']['include_full_flux_annulus'] == True:
                input_parameters['kxky_grid_box']['y0'] = '1./(rhostar*geo_surf%rhotor)'
    
    # Return the input parameters
    return input_parameters
    
#===============================================================================
#                        Calculate extra input parameters                      #
#===============================================================================

def calculate_extraInputParameters(input_parameters):
    
    # Fill in the species for the adiabatic electrons: custom knob for the GUI
    if input_parameters['species_options']['nspec'] == 1:
        input_parameters['species_parameters_a']['nine'] = input_parameters['adiabatic_electron_response']['nine']
        input_parameters['species_parameters_a']['tite'] = input_parameters['adiabatic_electron_response']['tite']
        input_parameters['species_parameters_a']['dens'] = round(1/input_parameters['adiabatic_electron_response']['nine'], 4)
        input_parameters['species_parameters_a']['temp'] = round(1/input_parameters['adiabatic_electron_response']['tite'], 4)
    
    # Calculate some extra variables
    input_parameters['geometry_vmec']['rho'] = np.sqrt(input_parameters['geometry_vmec']['torflux'])
    input_parameters['adiabatic_electron_response']['teti'] = 1/input_parameters['adiabatic_electron_response']['tite'] 
        
    # Return the input parameters
    return input_parameters

#-----------------------------------------      
def calculate_extraInputParametersFromWout(input_parameters, path_input_file):
    
    # Step size in real space
    y0 = input_parameters['kxky_grid_box']['y0']
    if input_parameters['gyrokinetic_terms']['include_full_flux_annulus']==True:
        y0 = np.sqrt(input_parameters['geometry_vmec']['torflux'])/input_parameters['physics_inputs']['rhostar']
        input_parameters['kxky_grid_box']['y0'] = y0
        if input_parameters['geometry_options']['geometry_option']!='vmec':
            exit_program('Not made for Miller yet.', calculate_extraInputParametersFromWout, sys._getframe().f_lineno)

    # Abbreviate the needed input parameters
    nzed = input_parameters['z_grid']['nzed'] 
    jtwist = input_parameters['kxky_grid_box']['jtwist']
    svalue = input_parameters['geometry_vmec']['torflux']
    vmec_filename = input_parameters['geometry_vmec']['vmec_filename']
    nfield_periods = input_parameters['geometry_vmec']['nfield_periods']
    if vmec_filename=='wout*.nc': nfield_periods = np.nan
    if vmec_filename=='wout*.nc': svalue = input_parameters['geometry_miller']['rhoc']*input_parameters['geometry_miller']['rhoc']
        
    # Read the VMEC file
    from stellapy.data.paths.load_pathObject import create_dummyPathObject
    from stellapy.data.geometry.read_wout import read_woutFile 
    woutParameters = read_woutFile(create_dummyPathObject(path_input_file, vmec_filename))
    
    # Calculate extra geometric quantities used in stella (depend on rho)
    if 'jtwist' not in woutParameters:
        from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_gridDivisionsAndSize
        woutParameters.update(calculate_gridDivisionsAndSize(y0, nfield_periods, woutParameters, svalue))
        
    # Save the VMEC Parameters 
    input_parameters['kxky_grid_box']['Lx'] = woutParameters['Lx']
    input_parameters['kxky_grid_box']['Ly'] = woutParameters['Ly']
    input_parameters['kxky_grid_box']['dkx'] = woutParameters['dkx']
    input_parameters['kxky_grid_box']['dky'] = woutParameters['dky']
    input_parameters['kxky_grid_box']['shat'] = woutParameters['shat']
    input_parameters['kxky_grid_box']['jtwist'] = woutParameters['jtwist'] if jtwist==-1 else jtwist
    
    # Calculate extra variables 
    poloidal_turns = np.round(nfield_periods/woutParameters['nfp']*abs(woutParameters['iota']),1)
    nzed_per_turn = int(nzed/poloidal_turns) if (not np.isnan(poloidal_turns)) else nzed
    
    # Save the extra variables under the corresponding knobs
    input_parameters['z_grid']['nz'] = nzed_per_turn
    input_parameters['geometry_vmec']['poloidal_turns'] = poloidal_turns
    return input_parameters 
 
#-----------------------------------------      
def calculate_extraInputParametersFromNetcdf(input_parameters, path_input_file): 

    # Initiate
    data = {}
    
    # Netcdf path
    path_netcdf = path_input_file.with_suffix('.out.nc')
    if '_dummy' in str(path_input_file):
        from stellapy.data.input.write_listOfMatchingInputFiles import read_inputFilesInDummyInputFile
        input_files_in_dummy_file = read_inputFilesInDummyInputFile(path_input_file)
        path_netcdf = input_files_in_dummy_file[0].with_suffix('.out.nc')
            
    # Read from the *.out.nc file       
    if os.path.isfile(path_netcdf): 
        netcdf_file = read_outputFile(path_netcdf)
        data['vec_kx'] = read_netcdfVariables('vec_kx', netcdf_file)
        data['vec_ky'] = read_netcdfVariables('vec_ky', netcdf_file)
        data['vec_mu'] = read_netcdfVariables('vec_mu', netcdf_file)
        data['vec_vpa'] = read_netcdfVariables('vec_vpa', netcdf_file)
        netcdf_file.close()
    
    # Read from the *.out.h5 file
    elif os.path.isfile(path_input_file.with_suffix('.out.h5')):
        try: 
            with h5py.File(path_input_file.with_suffix('.out.h5'), 'r') as h5_file: 
                if 'vec_kx' in h5_file.keys():  data['vec_kx'] = h5_file['vec_kx'][:]
                if 'vec_ky' in h5_file.keys():  data['vec_ky'] = h5_file['vec_ky'][:]
                if 'vec_mu' in h5_file.keys():  data['vec_mu'] = h5_file['vec_mu'][:]
                if 'vec_vpa' in h5_file.keys(): data['vec_vpa'] = h5_file['vec_vpa'][:]
        except:
            print('Something went wrong when reading the h5 file for:')
            print('     '+str(path_input_file.with_suffix('.out.h5')))
            sys.exit()
        
    # If both don't exist, return nan
    else:
        print('The netcdf file can not be found for:')
        print('      '+str(path_input_file))  
        return input_parameters
        
    # For a linear simulation, create artificial (kx,ky) vectors
    if len(data['vec_kx'])==1: data['vec_kx'] = np.array([0, data['vec_kx'][0]])
    if len(data['vec_ky'])==1: data['vec_ky'] = np.array([0, data['vec_ky'][0]])
        
    # Calculate the (kx,ky) dimensions   
    input_parameters['kxky_grid_box']['Lx'] = 2*np.pi/(data['vec_kx'][1]-data['vec_kx'][0]) if data['vec_kx'][0]!=data['vec_kx'][1] else np.Inf
    input_parameters['kxky_grid_box']['Ly'] = 2*np.pi/(data['vec_ky'][1]-data['vec_ky'][0]) if data['vec_ky'][0]!=data['vec_ky'][1] else np.Inf
    input_parameters['kxky_grid_box']['kx max'] = np.max(data['vec_kx'])
    input_parameters['kxky_grid_box']['ky max'] = np.max(data['vec_ky'])
    input_parameters['kxky_grid_box']['dkx'] = np.min(np.abs(data['vec_kx'][np.nonzero(data['vec_kx'])])) if len(data['vec_kx'][np.nonzero(data['vec_kx'])])!=0 else 0
    input_parameters['kxky_grid_box']['dky'] = np.min(np.abs(data['vec_ky'][np.nonzero(data['vec_ky'])])) if len(data['vec_ky'][np.nonzero(data['vec_ky'])])!=0 else 0 
    
    # Calculate the (mu,vpa) dimensions
    input_parameters['velocity_grids']['mu max'] = np.max(data['vec_mu'])
    input_parameters['velocity_grids']['vpa max'] = np.max(data['vec_vpa'])
    input_parameters['velocity_grids']['dmu'] = np.min(np.abs(data['vec_mu'][np.nonzero(data['vec_mu'])]))
    input_parameters['velocity_grids']['dvpa'] = np.min(np.abs(data['vec_vpa'][np.nonzero(data['vec_vpa'])]))
    return input_parameters 

        
#===============================================================================
#                 Read specific parameters from the input file                 #
#===============================================================================
 
def read_modeFromInputFile(path_input_file): 
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters('kxky_grid_range', input_parameters, input_parameters_read)
    kx = input_parameters['kxky_grid_range']['akx_min']
    ky = input_parameters['kxky_grid_range']['aky_min']
    return kx, ky
 
def read_numberOfModesFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters()
    input_parameters = update_input_parameters('kxky_grid_option', input_parameters, input_parameters_read)
    input_parameters = update_input_parameters('kxky_grid_box', input_parameters, input_parameters_read)
    input_parameters = update_input_parameters('kxky_grid_range', input_parameters, input_parameters_read)
    if input_parameters['kxky_grid_option']['grid_option'] in {'default', 'range'}:
        nakx = input_parameters['kxky_grid_range']['nakx']
        naky = input_parameters['kxky_grid_range']['naky']
        return nakx, naky
    if input_parameters['kxky_grid_option']['grid_option'] in {'box', 'annulus', 'nonlinear'}:
        from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_nakx, calculate_naky
        nx = input_parameters['kxky_grid_box']['nx']
        ny = input_parameters['kxky_grid_box']['ny']
        nakx = calculate_nakx(nx)
        naky = calculate_naky(ny)
        return nakx, naky
 
def read_vecKxKyFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters()
    input_parameters = update_input_parameters('geometry_options', input_parameters, input_parameters_read)
    input_parameters = update_input_parameters('physics_inputs', input_parameters, input_parameters_read)
    input_parameters = update_input_parameters('kxky_grid_option', input_parameters, input_parameters_read)
    input_parameters = update_input_parameters('kxky_grid_box', input_parameters, input_parameters_read)
    input_parameters = update_input_parameters('kxky_grid_range', input_parameters, input_parameters_read)
    if input_parameters['geometry_options']['geometry_option'] != 'vmec':
        exit_program('Miller is not implemented yet.', read_vecKxKyFromInputFile, sys._getframe().f_lineno)
    if input_parameters['kxky_grid_option']['grid_option'] in {'default', 'range'}:
        nakx = input_parameters['kxky_grid_range']['nakx']
        naky = input_parameters['kxky_grid_range']['naky']
        akx_min = input_parameters['kxky_grid_range']['akx_min']
        akx_max = input_parameters['kxky_grid_range']['akx_max'] #@UnusedVariable
        aky_min = input_parameters['kxky_grid_range']['aky_min']
        aky_max = input_parameters['kxky_grid_range']['aky_max']
        if nakx>1: exit_program('Not implemented nakx>1 yet.', read_vecKxKyFromInputFile, sys._getframe().f_lineno)
        dky = (aky_max - aky_min)/(naky - 1) if naky>1 else 0
        vec_ky = [ aky_min + dky*i for i in range(naky) ]
        vec_kx = [ akx_min ]
        return vec_kx, vec_ky
    if input_parameters['kxky_grid_option']['grid_option'] in {'box', 'annulus', 'nonlinear'}:
        from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_nakx, calculate_naky 
        nx = input_parameters['kxky_grid_box']['nx']
        ny = input_parameters['kxky_grid_box']['ny']
        rhostar = input_parameters['physics_inputs']['rhostar']
        nakx = calculate_nakx(nx); naky = calculate_naky(ny)
        torflux = read_svalueFromInputFile(path_input_file)
        rho = np.sqrt(torflux); y0 = rho/rhostar; dky = 1/y0
        vec_ky = [ dky*i for i in range(naky) ]
        if nakx>1: exit_program('nx>1 is not implemented yet.', read_vecKxKyFromInputFile, sys._getframe().f_lineno)
        vec_kx = [0]
        return vec_kx, vec_ky
 
def read_fullFluxSurfaceFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters('gyrokinetic_terms', input_parameters, input_parameters_read)
    return input_parameters['gyrokinetic_terms']['include_full_flux_annulus']
 
def read_linearNonlinearFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters('gyrokinetic_terms', input_parameters, input_parameters_read)
    return not input_parameters['gyrokinetic_terms']['include_nonlinear'], input_parameters['gyrokinetic_terms']['include_nonlinear']
 
def read_nonlinearFullFluxSurfaceFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters('gyrokinetic_terms', input_parameters, input_parameters_read)
    return input_parameters['gyrokinetic_terms']['include_nonlinear'], input_parameters['gyrokinetic_terms']['include_full_flux_annulus']
 
def read_vmecFileNameFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters('geometry_vmec', input_parameters, input_parameters_read)
    return input_parameters['geometry_vmec']['vmec_filename']
 
def read_nspecFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters('species_options', input_parameters, input_parameters_read)
    return input_parameters['species_options']['nspec']
 
def read_deltFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters('time_step', input_parameters, input_parameters_read)
    return input_parameters['time_step']['delt']
 
def read_svalueFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters('geometry_options', input_parameters, input_parameters_read)
    input_parameters = update_input_parameters('geometry_vmec', input_parameters, input_parameters_read)
    input_parameters = update_input_parameters('geometry_miller', input_parameters, input_parameters_read)
    if input_parameters['geometry_options']['geometry_option']!='vmec':
        svalue = input_parameters['geometry_miller']['rhoc']**2
    if input_parameters['geometry_options']['geometry_option']=='vmec':
        svalue = input_parameters['geometry_vmec']['torflux']
    return svalue
 
def read_tendFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters('time_trace_options', input_parameters, input_parameters_read)
    if 'tend' in input_parameters['time_trace_options']: return input_parameters['time_trace_options']['tend'] 
    return 

#===============================================================================
#              Read specific parameters from the input parameters              #
#===============================================================================
 
def read_nonlinearFromInputParameters(input_parameters):
    try: return input_parameters['gyrokinetic_terms']['include_nonlinear']
    except: return False
    
def read_deltFromInputParameters(input_parameters):
    try: return input_parameters['time_step']['delt']
    except: return 0.1
    
def read_nspecFromInputParameters(input_parameters):
    try: return input_parameters['species_options']['nspec']
    except: return 2
    
def read_vmecfilenameFromInputParameters(input_parameters):
    try: return input_parameters['geometry_vmec']['vmec_filename']
    except: return 'wout*.nc'
    
def read_geoOptionFromInputParameters(input_parameters):
    try: return input_parameters['geometry_options']['geometry_option']
    except: return 'local'
    
def read_nzedFromInputParameters(input_parameters):
    try: return input_parameters['z_grid']['nzed']
    except: return 24
    
def read_fiprimFromInputParameters(input_parameters):
    try: return input_parameters['species_parameters_1']['fprim']
    except: exit_program('Fprim should really be defined in the input file', read_fprimFromInputParameters, sys._getframe().f_lineno)
    
def read_fprimFromInputParameters(input_parameters):
    try: return input_parameters['species_parameters_1']['fprim']
    except: exit_program('Fprim should really be defined in the input file', read_fprimFromInputParameters, sys._getframe().f_lineno)
    
def read_tiprimFromInputParameters(input_parameters):
    try: return input_parameters['species_parameters_1']['tprim']
    except: exit_program('Tiprim should really be defined in the input file', read_tiprimFromInputParameters, sys._getframe().f_lineno)
    
def read_teprimFromInputParameters(input_parameters):
    try: return input_parameters['species_parameters_2']['tprim']
    except: exit_program('Teprim should really be defined in the input file', read_teprimFromInputParameters, sys._getframe().f_lineno)
    
def read_tzprimFromInputParameters(input_parameters):
    try: return input_parameters['species_parameters_3']['tprim']
    except: exit_program('Tzprim should really be defined in the input file', read_tzprimFromInputParameters, sys._getframe().f_lineno)
    
def read_ichargeFromInputParameters(input_parameters):
    try: return input_parameters['species_parameters_1']['z']
    except: exit_program('Main species charge should really be defined in the input file', read_ichargeFromInputParameters, sys._getframe().f_lineno)
    
def read_svalueFromInputParameters(input_parameters):
    try: geometry_option = input_parameters['geometry_options']['geometry_option']
    except: geometry_option = 'local'
    if geometry_option=='vmec':
        try: return input_parameters['geometry_vmec']['torflux']
        except: return 0.6354167
    if geometry_option!='vmec':
        try: return input_parameters['geometry_miller']['rhoc']**2
        except: exit_program('Rhoc must be defined.', read_svalueFromInputFile, sys._getframe().f_lineno)
        
def read_resolutionFromInputParameters(input_parameters):
    try: y0 = input_parameters['kxky_grid_box']['y0']
    except: y0 = -1
    try: nx = input_parameters['kxky_grid_box']['nx']
    except: nx = 1
    try: ny = input_parameters['kxky_grid_box']['ny']
    except: ny = 1
    try: nmu = input_parameters['velocity_grids']['nmu']
    except: nmu = 12
    try: nvgrid = input_parameters['velocity_grids']['nvgrid']
    except: nvgrid = 24
    try: nzed = input_parameters['z_grid']['nzed']
    except: nzed = 24
    return nzed, nmu, nvgrid, nx, ny, y0 

def read_switchesFromInputParameters(input_parameters):
    try: hyper_dissipation = input_parameters['dissipation_and_collisions_options']['hyper_dissipation']
    except: hyper_dissipation = False
    try: flip_flop = input_parameters['numerical_algorithms']['flip_flop']
    except: flip_flop = False
    try: mirror_implicit = input_parameters['numerical_algorithms']['mirror_implicit']
    except: mirror_implicit = True
    try: stream_implicit = input_parameters['numerical_algorithms']['stream_implicit']
    except: stream_implicit = True
    return hyper_dissipation, flip_flop, mirror_implicit, stream_implicit 

def read_explicitOptionFromInputParameters(input_parameters):
    try: explicit_option = input_parameters['numerical_algorithms']['explicit_algorithm']
    except: explicit_option = 'rk3'
    return explicit_option

#===============================================================================
#                Read specific data types from the input file                  #
#===============================================================================

def read_parameterFromIniFile(section, parameter, default_value):
    ''' Read the variables if they exist and convert them to the correct datatype. '''
    if default_value==None: return float(section[parameter])
    if isinstance(default_value, str): return section[parameter]
    if isinstance(default_value, bool): return True if section[parameter]=='True' else False
    if isinstance(default_value, int): return int(section[parameter])
    if isinstance(default_value, float): return float(section[parameter])
    return 
 
#--------------------------------------------
def update_input_parameters(knob, input_parameters, input_parameters_read):
    try: input_parameters[knob].update(input_parameters_read[knob])
    except: pass # The knob was not in the input file
    return input_parameters
