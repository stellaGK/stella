""" 

#===============================================================================
#                              Read the input file                             #
#===============================================================================

Read the input file located at <path_input_file>. Additionally, calculate some 
of the interal stella variables and some extra parameters such as {Lx, Ly}.
    
Hanne Thienpondt
30/01/2023

"""

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
    """ Read the input file and overwrite the default stella parameters. """ 
            
    # Read the input file and overwrite the default stella parameters.   
    if os.path.isfile(path_input_file): 
        if path_input_file.suffix==".in":  input_parameters = read_inFile(path_input_file)
        if path_input_file.suffix==".ini": input_parameters = read_iniFile(path_input_file)
        input_parameters["parameters"]["teti"] = 1/input_parameters["parameters"]["tite"]
        return input_parameters
    if os.path.isfile(path_input_file.with_suffix(".in")): 
        input_parameters = read_inFile(path_input_file.with_suffix(".in"))
        input_parameters["parameters"]["teti"] = 1/input_parameters["parameters"]["tite"]
        return input_parameters

    # Critical error if we didn't find any data   
    exit_reason = "The input data can not be found. The following file does not exist:\n"    
    exit_reason += "    "+str(path_input_file)+"\n\n"
    exit_reason += "You can try to write it through the command:\n"    
    exit_reason += "    write_dataFiles -s ini \n\n"
    exit_reason += "If the data was written before October 2022 try:\n"    
    exit_reason += "    write_listOfMatchingInputFilesForOldStellapy \n" 
    exit_program(exit_reason, read_inputFile, sys._getframe().f_lineno)   
    return

#------------------------- 
def read_inFile(path_input_file):
    
    # Read the input namelist and parse it to a nested dict
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))  
    if 'reinit_knobs' in input_parameters_read: del input_parameters_read['reinit_knobs'] 

    # Initiate the dictionary: load the default stella parameters
    input_parameters = load_defaultInputParameters() 
    input_parameters["species_knobs"].update(input_parameters_read["species_knobs"])

    # Add more default species if nspec>2
    for i in range(2, input_parameters["species_knobs"]["nspec"]+1): 
        input_parameters["species_parameters_"+str(i)] = input_parameters["species_parameters_1"].copy() 

    # Overwrite the default values if they have been changed in the <path_input_file> 
    for knob in input_parameters_read.keys():    
        input_parameters[knob].update(input_parameters_read[knob])

    # Close the input file 
    return input_parameters

#-------------------------
def read_iniFile(path_input_file):
    
    # Read the "*.ini" file  
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
                if parameter=='t_end': 
                    input_parameters_read[knob]['tend'] = input_parameters_read[knob][parameter]
                    del input_parameters_read[knob][parameter]; parameter = 'tend'; 
                if parameter=='write_fluxes_kxky': 
                    input_parameters_read[knob]['write_fluxes_kxkyz'] = input_parameters_read[knob][parameter]
                    del input_parameters_read[knob][parameter]; parameter = 'write_fluxes_kxkyz'; 
                input_parameters_read[knob][parameter] = read_parameterFromIniFile(input_parameters_read[knob], parameter, input_parameters[knob][parameter])
                    
    # Add more default species if nspec>2
    if "species_knobs" in input_parameters_read: input_parameters["species_knobs"].update(input_parameters_read["species_knobs"])
    for i in range(2, input_parameters["species_knobs"]["nspec"]+1): 
        input_parameters["species_parameters_"+str(i)] = input_parameters["species_parameters_1"].copy() 
        
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
    ''' Read "*.in" file and return dict[knobs][variable].  '''    
    
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
    if path_input_file.suffix==".in":
    
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
     
    # Calculate <irad_min> and <irad_max>
    if input_parameters["sfincs_input"]["irad_min"] == "-nradii/2":
        nradii = input_parameters["neoclassical_input"]["nradii"]
        input_parameters["sfincs_input"]["irad_min"] = -nradii/2 
        input_parameters["sfincs_input"]["irad_max"] = nradii/2  
    else: print("WARNING: <irad_min> was set by the input file but it should be calculated indirectly through <nradii>.")

    # Calculate <nzgrid>
    if input_parameters["zgrid_parameters"]["nzgrid"] == "nzed/2 + (nperiod-1)*nzed":
        nzed = input_parameters["zgrid_parameters"]["nzed"]
        nperiod = input_parameters["zgrid_parameters"]["nperiod"]
        input_parameters["zgrid_parameters"]["nzgrid"] = nzed/2 + (nperiod-1)*nzed
    else: print("WARNING: <nzgrid> was set by the input file but it should be calculated indirectly through <nzed> and <nperiod>.")
 
    # Calculate <zgrid_scalefac> and <zgrid_refinement_factor>
    if input_parameters["vmec_parameters"]["zgrid_scalefac"] == "2.0 if zed_equal_arc else 1.0":
        zed_equal_arc = input_parameters["zgrid_parameters"]["zed_equal_arc"]
        input_parameters["vmec_parameters"]["zgrid_scalefac"]           = 2.0 if zed_equal_arc else 1.0
        input_parameters["vmec_parameters"]["zgrid_refinement_factor"]  = 4 if zed_equal_arc else 1
    else: print("WARNING: <zgrid_scalefac> was set by the input file but it should be calculated indirectly through <zed_equal_arc>.")
    
    # Caclulate <nakx> and <naky> 
    if input_parameters["kt_grids_box_parameters"]["naky"] == "(ny-1)/3 + 1":
        from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_nakx, calculate_naky
        ny = input_parameters["kt_grids_box_parameters"]["ny"]
        nx = input_parameters["kt_grids_box_parameters"]["nx"]
        input_parameters["kt_grids_box_parameters"]["naky"] = calculate_naky(ny)
        input_parameters["kt_grids_box_parameters"]["nakx"] = calculate_nakx(nx)

    # Calculate <y0> for a full flux surface simulation
    if input_parameters["physics_flags"]["nonlinear"] == True: 
        if input_parameters["kt_grids_box_parameters"]["y0"] == -1.0:
            if input_parameters["physics_flags"]["full_flux_surface"] == False:
                print("WARNING: When simulating a flux tube, y0 needs to be set in the input file.")
            if input_parameters["physics_flags"]["full_flux_surface"] == True:
                input_parameters["kt_grids_box_parameters"]["y0"] = "1./(rhostar*geo_surf%rhotor)"  
    
    # Return the input parameters
    return input_parameters
    
#===============================================================================
#                        Calculate extra input parameters                      #
#===============================================================================

def calculate_extraInputParameters(input_parameters):
    
    # Fill in the species for the adiabatic electrons: custom knob for the GUI
    if input_parameters["species_knobs"]["nspec"] == 1:
        input_parameters["species_parameters_a"]["nine"] = input_parameters["parameters"]["nine"]
        input_parameters["species_parameters_a"]["tite"] = input_parameters["parameters"]["tite"]
        input_parameters["species_parameters_a"]["dens"] = round(1/input_parameters["parameters"]["nine"], 4)
        input_parameters["species_parameters_a"]["temp"] = round(1/input_parameters["parameters"]["tite"], 4) 
    
    # Calculate some extra variables 
    input_parameters["vmec_parameters"]["rho"] = np.sqrt(input_parameters["vmec_parameters"]["torflux"])  
    input_parameters["parameters"]["teti"] = 1/input_parameters["parameters"]["tite"] 
        
    # Return the input parameters
    return input_parameters

#-----------------------------------------      
def calculate_extraInputParametersFromWout(input_parameters, path_input_file):
    
    # Step size in real space 
    y0 = input_parameters["kt_grids_box_parameters"]["y0"] 
    if input_parameters["physics_flags"]["full_flux_surface"]==True:
        y0 = np.sqrt(input_parameters["vmec_parameters"]["torflux"])/input_parameters["parameters"]["rhostar"]
        input_parameters["kt_grids_box_parameters"]["y0"] = y0
        if input_parameters["geo_knobs"]["geo_option"]!="vmec":
            exit_program("Not made for Miller yet.", calculate_extraInputParametersFromWout, sys._getframe().f_lineno)

    # Abbreviate the needed input parameters
    nzed = input_parameters["zgrid_parameters"]["nzed"] 
    jtwist = input_parameters["kt_grids_box_parameters"]["jtwist"]
    svalue = input_parameters["vmec_parameters"]["torflux"] 
    vmec_filename = input_parameters["vmec_parameters"]["vmec_filename"] 
    nfield_periods = input_parameters["vmec_parameters"]["nfield_periods"] 
    if vmec_filename=='wout*.nc': nfield_periods = np.nan
    if vmec_filename=='wout*.nc': svalue = input_parameters["millergeo_parameters"]["rhoc"]*input_parameters["millergeo_parameters"]["rhoc"]
        
    # Read the VMEC file
    from stellapy.data.paths.load_pathObject import create_dummyPathObject
    from stellapy.data.geometry.read_wout import read_woutFile 
    woutParameters = read_woutFile(create_dummyPathObject(path_input_file, vmec_filename))    
    
    # Calculate extra geometric quantities used in stella (depend on rho)
    if "jtwist" not in woutParameters:
        from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_gridDivisionsAndSize
        woutParameters.update(calculate_gridDivisionsAndSize(y0, nfield_periods, woutParameters, svalue))
        
    # Save the VMEC Parameters 
    input_parameters["kt_grids_box_parameters"]["Lx"] = woutParameters["Lx"] 
    input_parameters["kt_grids_box_parameters"]["Ly"] = woutParameters["Ly"] 
    input_parameters["kt_grids_box_parameters"]["dkx"] = woutParameters["dkx"] 
    input_parameters["kt_grids_box_parameters"]["dky"] = woutParameters["dky"] 
    input_parameters["kt_grids_box_parameters"]["shat"] = woutParameters["shat"]
    input_parameters["kt_grids_box_parameters"]["jtwist"] = woutParameters["jtwist"] if jtwist==-1 else jtwist
    
    # Calculate extra variables 
    poloidal_turns = np.round(nfield_periods/woutParameters['nfp']*abs(woutParameters['iota']),1) 
    nzed_per_turn = int(nzed/poloidal_turns) if (not np.isnan(poloidal_turns)) else nzed
    
    # Save the extra variables under the corresponding knobs
    input_parameters["zgrid_parameters"]["nz"] = nzed_per_turn
    input_parameters["vmec_parameters"]["poloidal_turns"] = poloidal_turns 
    return input_parameters 
 
#-----------------------------------------      
def calculate_extraInputParametersFromNetcdf(input_parameters, path_input_file): 

    # Initiate
    data = {}
    
    # Netcdf path
    path_netcdf = path_input_file.with_suffix('.out.nc')
    if "_dummy" in str(path_input_file):
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
            print("Something went wrong when reading the h5 file for:")
            print("     "+str(path_input_file.with_suffix('.out.h5')))
            sys.exit()
        
    # If both don't exist, return nan
    else:
        print("The netcdf file can not be found for:")
        print("      "+str(path_input_file))  
        return input_parameters
        
    # For a linear simulation, create artificial (kx,ky) vectors
    if len(data['vec_kx'])==1: data['vec_kx'] = np.array([0, data['vec_kx'][0]])  
    if len(data['vec_ky'])==1: data['vec_ky'] = np.array([0, data['vec_ky'][0]])  
        
    # Calculate the (kx,ky) dimensions   
    input_parameters["kt_grids_box_parameters"]["Lx"] = 2*np.pi/(data['vec_kx'][1]-data['vec_kx'][0]) if data['vec_kx'][0]!=data['vec_kx'][1] else np.Inf
    input_parameters["kt_grids_box_parameters"]["Ly"] = 2*np.pi/(data['vec_ky'][1]-data['vec_ky'][0]) if data['vec_ky'][0]!=data['vec_ky'][1] else np.Inf
    input_parameters["kt_grids_box_parameters"]["kx max"] = np.max(data['vec_kx'])
    input_parameters["kt_grids_box_parameters"]["ky max"] = np.max(data['vec_ky'])
    input_parameters["kt_grids_box_parameters"]["dkx"] = np.min(np.abs(data['vec_kx'][np.nonzero(data['vec_kx'])])) if len(data['vec_kx'][np.nonzero(data['vec_kx'])])!=0 else 0
    input_parameters["kt_grids_box_parameters"]["dky"] = np.min(np.abs(data['vec_ky'][np.nonzero(data['vec_ky'])])) if len(data['vec_ky'][np.nonzero(data['vec_ky'])])!=0 else 0 
    
    # Calculate the (mu,vpa) dimensions
    input_parameters["vpamu_grids_parameters"]["mu max"] = np.max(data['vec_mu'])
    input_parameters["vpamu_grids_parameters"]["vpa max"] = np.max(data['vec_vpa'])
    input_parameters["vpamu_grids_parameters"]["dmu"] = np.min(np.abs(data['vec_mu'][np.nonzero(data['vec_mu'])]))
    input_parameters["vpamu_grids_parameters"]["dvpa"] = np.min(np.abs(data['vec_vpa'][np.nonzero(data['vec_vpa'])]))     
    return input_parameters 

        
#===============================================================================
#                 Read specific parameters from the input file                 #
#===============================================================================
 
def read_modeFromInputFile(path_input_file): 
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters("kt_grids_range_parameters", input_parameters, input_parameters_read)  
    kx = input_parameters["kt_grids_range_parameters"]["akx_min"]   
    ky = input_parameters["kt_grids_range_parameters"]["aky_min"]    
    return kx, ky
 
def read_numberOfModesFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters("kt_grids_knobs", input_parameters, input_parameters_read)  
    input_parameters = update_input_parameters("kt_grids_box_parameters", input_parameters, input_parameters_read)  
    input_parameters = update_input_parameters("kt_grids_range_parameters", input_parameters, input_parameters_read)   
    if input_parameters["kt_grids_knobs"]["grid_option"] in {"default", "range"}:
        nakx = input_parameters["kt_grids_range_parameters"]["nakx"]     
        naky = input_parameters["kt_grids_range_parameters"]["naky"]  
        return nakx, naky
    if input_parameters["kt_grids_knobs"]["grid_option"] in {"box", "annulus", "nonlinear"}:
        from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_nakx, calculate_naky 
        nx = input_parameters["kt_grids_box_parameters"]["nx"]     
        ny = input_parameters["kt_grids_box_parameters"]["ny"] 
        nakx = calculate_nakx(nx)
        naky = calculate_naky(ny)  
        return nakx, naky
 
def read_vecKxKyFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters("geo_knobs", input_parameters, input_parameters_read)   
    input_parameters = update_input_parameters("parameters", input_parameters, input_parameters_read)  
    input_parameters = update_input_parameters("kt_grids_knobs", input_parameters, input_parameters_read)  
    input_parameters = update_input_parameters("kt_grids_box_parameters", input_parameters, input_parameters_read)  
    input_parameters = update_input_parameters("kt_grids_range_parameters", input_parameters, input_parameters_read)  
    if input_parameters["geo_knobs"]["geo_option"] != "vmec":
        exit_program("Miller is not implemented yet.", read_vecKxKyFromInputFile, sys._getframe().f_lineno)
    if input_parameters["kt_grids_knobs"]["grid_option"] in {"default", "range"}:
        nakx = input_parameters["kt_grids_range_parameters"]["nakx"]     
        naky = input_parameters["kt_grids_range_parameters"]["naky"] 
        akx_min = input_parameters["kt_grids_range_parameters"]["akx_min"]  
        akx_max = input_parameters["kt_grids_range_parameters"]["akx_max"] #@UnusedVariable 
        aky_min = input_parameters["kt_grids_range_parameters"]["aky_min"]  
        aky_max = input_parameters["kt_grids_range_parameters"]["aky_max"]       
        if nakx>1: exit_program("Not implemented nakx>1 yet.", read_vecKxKyFromInputFile, sys._getframe().f_lineno) 
        dky = (aky_max - aky_min)/(naky - 1) if naky>1 else 0 
        vec_ky = [ aky_min + dky*i for i in range(naky) ]
        vec_kx = [ akx_min ]
        return vec_kx, vec_ky
    if input_parameters["kt_grids_knobs"]["grid_option"] in {"box", "annulus", "nonlinear"}:
        from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_nakx, calculate_naky 
        nx = input_parameters["kt_grids_box_parameters"]["nx"]     
        ny = input_parameters["kt_grids_box_parameters"]["ny"] 
        rhostar = input_parameters["parameters"]["rhostar"]
        nakx = calculate_nakx(nx); naky = calculate_naky(ny)  
        torflux = read_svalueFromInputFile(path_input_file) 
        rho = np.sqrt(torflux); y0 = rho/rhostar; dky = 1/y0
        vec_ky = [ dky*i for i in range(naky) ]
        if nakx>1: exit_program("nx>1 is not implemented yet.", read_vecKxKyFromInputFile, sys._getframe().f_lineno) 
        vec_kx = [0]
        return vec_kx, vec_ky
 
def read_fullFluxSurfaceFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters("physics_flags", input_parameters, input_parameters_read)   
    return input_parameters["physics_flags"]["full_flux_surface"]
 
def read_linearNonlinearFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters("physics_flags", input_parameters, input_parameters_read)    
    return not input_parameters["physics_flags"]["nonlinear"], input_parameters["physics_flags"]["nonlinear"]
 
def read_nonlinearFullFluxSurfaceFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters("physics_flags", input_parameters, input_parameters_read)    
    return input_parameters["physics_flags"]["nonlinear"], input_parameters["physics_flags"]["full_flux_surface"]
 
def read_vmecFileNameFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters("vmec_parameters", input_parameters, input_parameters_read)   
    return input_parameters["vmec_parameters"]["vmec_filename"]
 
def read_nspecFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters("species_knobs", input_parameters, input_parameters_read)    
    return input_parameters["species_knobs"]["nspec"]
 
def read_deltFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters("knobs", input_parameters, input_parameters_read)    
    return input_parameters["knobs"]["delt"]
 
def read_svalueFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters("geo_knobs", input_parameters, input_parameters_read)    
    input_parameters = update_input_parameters("vmec_parameters", input_parameters, input_parameters_read)    
    input_parameters = update_input_parameters("millergeo_parameters", input_parameters, input_parameters_read)      
    if input_parameters["geo_knobs"]["geo_option"]!="vmec": 
        svalue = input_parameters["millergeo_parameters"]["rhoc"]**2
    if input_parameters["geo_knobs"]["geo_option"]=="vmec": 
        svalue = input_parameters["vmec_parameters"]["torflux"]  
    return svalue
 
def read_tendFromInputFile(path_input_file):
    input_parameters_read = f90nml.read(path_input_file)
    input_parameters_read = json.loads(json.dumps(input_parameters_read))
    input_parameters = load_defaultInputParameters() 
    input_parameters = update_input_parameters("knobs", input_parameters, input_parameters_read)    
    if "tend" in input_parameters["knobs"]: return input_parameters["knobs"]["tend"]
    if "t_end" in input_parameters["knobs"]: return input_parameters["knobs"]["t_end"] 
    return 

#===============================================================================
#              Read specific parameters from the input parameters              #
#===============================================================================
 
def read_nonlinearFromInputParameters(input_parameters):
    try: return input_parameters['physics_flags']['nonlinear']
    except: return False
    
def read_deltFromInputParameters(input_parameters):
    try: return input_parameters['knobs']['delt']
    except: return 0.1
    
def read_nspecFromInputParameters(input_parameters):
    try: return input_parameters['species_knobs']['nspec']
    except: return 2
    
def read_vmecfilenameFromInputParameters(input_parameters):
    try: return input_parameters['vmec_parameters']['vmec_filename']
    except: return 'wout*.nc'
    
def read_geoOptionFromInputParameters(input_parameters):
    try: return input_parameters['geo_knobs']['geo_option']
    except: return 'local'
    
def read_nzedFromInputParameters(input_parameters):
    try: return input_parameters['zgrid_parameters']['nzed']
    except: return 24
    
def read_fiprimFromInputParameters(input_parameters):
    try: return input_parameters['species_parameters_1']['fprim']
    except: exit_program("Fprim should really be defined in the input file", read_fprimFromInputParameters, sys._getframe().f_lineno)
    
def read_fprimFromInputParameters(input_parameters):
    try: return input_parameters['species_parameters_1']['fprim']
    except: exit_program("Fprim should really be defined in the input file", read_fprimFromInputParameters, sys._getframe().f_lineno)
    
def read_tiprimFromInputParameters(input_parameters):
    try: return input_parameters['species_parameters_1']['tprim']
    except: exit_program("Tiprim should really be defined in the input file", read_tiprimFromInputParameters, sys._getframe().f_lineno)
    
def read_teprimFromInputParameters(input_parameters):
    try: return input_parameters['species_parameters_2']['tprim']
    except: exit_program("Teprim should really be defined in the input file", read_teprimFromInputParameters, sys._getframe().f_lineno)
    
def read_tzprimFromInputParameters(input_parameters):
    try: return input_parameters['species_parameters_3']['tprim']
    except: exit_program("Tzprim should really be defined in the input file", read_tzprimFromInputParameters, sys._getframe().f_lineno)
    
def read_ichargeFromInputParameters(input_parameters):
    try: return input_parameters['species_parameters_1']['z']
    except: exit_program("Main species charge should really be defined in the input file", read_ichargeFromInputParameters, sys._getframe().f_lineno)
    
def read_svalueFromInputParameters(input_parameters):
    try: geo_option = input_parameters['geo_knobs']['geo_option']
    except: geo_option = 'local'
    if geo_option=='vmec':
        try: return input_parameters['vmec_parameters']['torflux'] 
        except: return 0.6354167
    if geo_option!='vmec':
        try: return input_parameters['millergeo_parameters']['rhoc']**2
        except: exit_program('Rhoc must be defined.', read_svalueFromInputFile, sys._getframe().f_lineno)
        
def read_resolutionFromInputParameters(input_parameters):
    try: y0 = input_parameters["kt_grids_box_parameters"]["y0"]
    except: y0 = -1
    try: nx = input_parameters["kt_grids_box_parameters"]["nx"]
    except: nx = 1
    try: ny = input_parameters["kt_grids_box_parameters"]["ny"]
    except: ny = 1
    try: nmu = input_parameters["vpamu_grids_parameters"]["nmu"]
    except: nmu = 12
    try: nvgrid = input_parameters["vpamu_grids_parameters"]["nvgrid"]
    except: nvgrid = 24
    try: nzed = input_parameters["zgrid_parameters"]["nzed"]
    except: nzed = 24
    return nzed, nmu, nvgrid, nx, ny, y0 

def read_switchesFromInputParameters(input_parameters):
    try: hyper_dissipation = input_parameters["dissipation"]["hyper_dissipation"]
    except: hyper_dissipation = False
    try: flip_flop = input_parameters["time_advance_knobs"]["flip_flop"]
    except: flip_flop = False
    try: mirror_implicit = input_parameters["knobs"]["mirror_implicit"]
    except: mirror_implicit = True
    try: stream_implicit = input_parameters["knobs"]["stream_implicit"]
    except: stream_implicit = True
    return hyper_dissipation, flip_flop, mirror_implicit, stream_implicit 

def read_explicitOptionFromInputParameters(input_parameters):
    try: explicit_option = input_parameters["time_advance_knobs"]["explicit_option"]
    except: explicit_option = 'rk3'
    return explicit_option

#===============================================================================
#                Read specific data types from the input file                  #
#===============================================================================

def read_parameterFromIniFile(section, parameter, default_value):
    ''' Read the variables if they exist and convert them to the correct datatype. ''' 
    if default_value==None: return float(section[parameter])
    if isinstance(default_value, str): return section[parameter]
    if isinstance(default_value, bool): return True if section[parameter]=="True" else False
    if isinstance(default_value, int): return int(section[parameter])
    if isinstance(default_value, float): return float(section[parameter])
    return 
 
#--------------------------------------------
def update_input_parameters(knob, input_parameters, input_parameters_read):
    try: input_parameters[knob].update(input_parameters_read[knob]) 
    except: pass # The knob was not in the input file
    return input_parameters

#--------------------------------------------
if __name__ == "__main__":
    
    # Names
    input_file = pathlib.Path("/home/hanne/CIEMAT/RUNS/input.in")
    new_input_file = pathlib.Path("/home/hanne/CIEMAT/RUNS/input_new.in") 
     
    # Read and write a name list
    input_parameters = f90nml.read(input_file)
    input_parameters.indent = "  "
    input_parameters.write(new_input_file, force=True)
     
    # Remove the spacings between knobs
    with open(new_input_file, 'r') as file :
        filedata = file.read()
    filedata = filedata.replace('/\n', '/')
    with open(new_input_file, 'w') as file:
        file.write(filedata) 
         
    # Get the same dictionary
    input_parameters = read_inputFile(input_file) 
    print(input_parameters)