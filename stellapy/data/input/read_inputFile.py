
import h5py
import os, sys
import numpy as np
import configparser 
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.utils.decorators.printWhichFileWeRead import printWhichFileWeRead
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables
from stellapy.data.stella.load_defaultInputParameters import load_defaultInputParameters
# TODO: replace --> grep -r "read_inputParameters" by read_inFile or simulation.input!
    
#===============================================================================
#                            READ THE INPUT FILE
#===============================================================================

def read_inputFile(input_file):
    """ Read the input file and overwrite the default stella parameters. """ 
            
    # Read the input file and overwrite the default stella parameters.   
    if os.path.isfile(input_file): 
        if input_file.suffix==".in":  input_parameters = read_inFile(input_file)
        if input_file.suffix==".ini": input_parameters = read_iniFile(input_file)
        input_parameters["parameters"]["teti"] = 1/input_parameters["parameters"]["tite"]
        return input_parameters

    # Critical error if we didn't find any data  
    exit_reason = "The input data can not be found.\n"   
    exit_reason += "    "+str(input_file)
    exit_program(exit_reason, read_inputFile, sys._getframe().f_lineno)   
    return

#------------------------- 
def read_inFile(input_file):
    printWhichFileWeRead("READING THE INPUT PARAMETERS FROM THE TEXT FILE") 
    
    # Read the "*.in" file  
    input_data = open(input_file, 'r')
    input_text = input_data.read().replace(' ', '')

    # Initiate the dictionary: load the default stella parameters
    input_parameters = load_defaultInputParameters() 

    # Add more default species if nspec>2
    nspec = read_integerInput(input_text, 'nspec') 
    for i in range(2, nspec+1): 
        input_parameters["species_parameters_"+str(i)] = input_parameters["species_parameters_1"].copy() 

    # Overwrite the default values if they have been changed in the <input_file> 
    for knob in input_parameters.keys():  
        if ("&"+knob) in input_text:
            input_text_knob = input_text.split("&"+knob)[1].split("&")[0]   
            for parameter in input_parameters[knob].keys(): 
                input_parameters[knob][parameter] = read_parameterFromInputFile(input_text_knob, parameter, input_parameters[knob][parameter])
    
    # Close the input file
    input_data.close() 
    return input_parameters

#-------------------------
def read_iniFile(input_file):
    printWhichFileWeRead("READING THE INPUT PARAMETERS FROM THE CONFIGURATION FILE") 
    
    # Read the "*.ini" file  
    ini_data = configparser.ConfigParser()
    ini_data.read(input_file)

    # Initiate the dictionary: load the default stella parameters
    input_parameters = load_defaultInputParameters() 

    # Add more default species if nspec>2
    if ini_data.has_option("species_knobs", "nspec"): nspec = int(ini_data["species_knobs"]["nspec"])
    else: nspec = int(input_parameters["species_knobs"]["nspec"])
    for i in range(2, nspec+1): 
        input_parameters["species_parameters_"+str(i)] = input_parameters["species_parameters_1"].copy() 

    # Overwrite the default values if they have been changed in the <input_file> 
    for knob in input_parameters.keys(): 
        if knob in ini_data:  
            for parameter in input_parameters[knob].keys(): 
                if parameter in ini_data[knob]:  
                    input_parameters[knob][parameter] = read_parameterFromIniFile(ini_data[knob], parameter, input_parameters[knob][parameter])
    
    # Return the input_parameters 
    return input_parameters

#===============================================================================
#                  CREATE THE INPUT PARAMETERS FROM THE INPUT FILE
#===============================================================================

def read_inputParameters(input_file):
    ''' Read "*.in" file and return dict[knobs][variable].  '''    
    
    # Read the input file and overwrite the default stella parameters 
    input_parameters = read_inputFile(input_file)
    
    # Read indirect input parameters
    input_parameters = read_extraInputParameters(input_parameters, input_file)
    
    # Return the input_parameters   
    return input_parameters 

#===============================================================================
#                   CALCULATE THE INTERNAL STELLA VARIABLES
#===============================================================================

def read_extraInputParameters(input_parameters, input_file):
    
    # This is only needed when we're reading the input file 
    if input_file.suffix==".in":
    
        # Calculate the input variables that stella calculates internally
        input_parameters = calculate_stellaVariables(input_parameters)
        
        # Calculate some extra input variables
        input_parameters = calculate_extraInputParameters(input_parameters) 
        
        # Calculate some extra parameters from the VMEC file
        input_parameters = calculate_extraInputParametersFromWout(input_parameters, input_file) 
        input_parameters = calculate_extraInputParametersFromNetcdf(input_parameters, input_file)
        
    return input_parameters
    
    
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
    else: print("WARNING: <nakx> was set by the input file but it should be calculated indirectly through <nx>.")
    
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
#                     CALCULATE EXTRA INPUT VARIABLES
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
    #input_parameters["parameters"]["teti"] = 1/input_parameters["parameters"]["tite"] 
        
    # Return the input parameters
    return input_parameters

#-----------------------------------------      
def calculate_extraInputParametersFromWout(input_parameters, input_file, geometry_path=None):
    printWhichFileWeRead("READING THE VMEC FILE FOR EXTRA INPUT PARAMETERS") 
    
    # Abbreviate the needed input parameters
    y0 = input_parameters["kt_grids_box_parameters"]["y0"] 
    nzed = input_parameters["zgrid_parameters"]["nzed"] 
    jtwist = input_parameters["kt_grids_box_parameters"]["jtwist"]
    svalue = input_parameters["vmec_parameters"]["torflux"] 
    nonlinear = input_parameters["physics_flags"]["nonlinear"]  
    vmec_filename = input_parameters["vmec_parameters"]["vmec_filename"] 
    nfield_periods = input_parameters["vmec_parameters"]["nfield_periods"] 
    if vmec_filename=='wout*.nc': nfield_periods = np.nan
    if vmec_filename=='wout*.nc': svalue = input_parameters["millergeo_parameters"]["rhoc"]*input_parameters["millergeo_parameters"]["rhoc"]
        
    # Read the VMEC file
    from stellapy.data.paths.load_pathObject import create_dummyPathObject
    from stellapy.data.geometry.read_wout import read_woutFile 
    woutParameters = read_woutFile(create_dummyPathObject(input_file, vmec_filename, nonlinear))    
    
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
def calculate_extraInputParametersFromNetcdf(input_parameters, input_file, data={}):
    printWhichFileWeRead("READING THE NETCDF FILE FOR EXTRA INPUT PARAMETERS")
    
    # Read from the *.out.h5 file
    if os.path.isfile(input_file.with_suffix('.out.h5')):  
        try: 
            with h5py.File(input_file.with_suffix('.out.h5'), 'r') as h5_file: 
                if 'vec_kx' in h5_file.keys():  data['vec_kx'] = h5_file['vec_kx'][:]
                if 'vec_ky' in h5_file.keys():  data['vec_ky'] = h5_file['vec_ky'][:]
                if 'vec_mu' in h5_file.keys():  data['vec_mu'] = h5_file['vec_mu'][:]
                if 'vec_vpa' in h5_file.keys(): data['vec_vpa'] = h5_file['vec_vpa'][:] 
        except:
            print("Something went wrong when reading the h5 file for:")
            print("     "+str(input_file.with_suffix('.out.h5')))
            sys.exit()
            
    # Read from the *.out.nc file       
    elif os.path.isfile(input_file.with_suffix('.out.nc')): 
        netcdf_file = read_outputFile(input_file.with_suffix('.out.nc'))
        data['vec_kx'] = read_netcdfVariables('vec_kx', netcdf_file)
        data['vec_ky'] = read_netcdfVariables('vec_ky', netcdf_file)
        data['vec_mu'] = read_netcdfVariables('vec_mu', netcdf_file)
        data['vec_vpa'] = read_netcdfVariables('vec_vpa', netcdf_file)
        netcdf_file.close()
        
    # If both don't exist, return nan
    else:
        print("The netcdf file can not be found for:")
        print("      "+str(input_file))  
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
#                READ SPECIFIC PARAMETERS FROM THE INPUT FILE
#===============================================================================

def read_modeFromInputFile(input_file):
    with open(input_file, 'r') as input_data:
        input_text = input_data.read().replace(' ', '')
        input_text_knob = input_text.split("&kt_grids_range_parameters")[1].split("/")[0] 
        kx = read_floatInput(input_text_knob, 'akx_min')  
        ky = read_floatInput(input_text_knob, 'aky_min')  
    return kx, ky

def read_numberOfModesFromInputFile(input_file):
    with open(input_file, 'r') as input_data:
        input_text = input_data.read().replace(' ', '') 
        input_text_knob = input_text.split("&kt_grids_range_parameters")[1].split("/")[0] 
        naky = read_integerInput(input_text_knob, 'naky')  
        nakx = read_integerInput(input_text_knob, 'nakx')  
    return nakx, naky

def read_vecKxKyFromInputFile(input_file):
    with open(input_file, 'r') as input_data:
        input_text = input_data.read().replace(' ', '') 
        input_text_knob = input_text.split("&kt_grids_range_parameters")[1].split("/")[0] 
        naky = read_integerInput(input_text_knob, 'naky')  
        nakx = read_integerInput(input_text_knob, 'nakx')  
        kx_min = read_floatInput(input_text_knob, 'akx_min')  
        ky_min = read_floatInput(input_text_knob, 'aky_min')
        ky_max = read_floatInput(input_text_knob, 'aky_max')
        if nakx>1: exit_program("Not implemented nakx>1 yet.", read_vecKxKyFromInputFile, sys._getframe().f_lineno) 
        dky = (ky_max - ky_min)/(naky - 1) if naky>1 else 0 
        vec_ky = [ ky_min + dky*i for i in range(naky) ]
        vec_kx = [ kx_min ]
    return vec_kx, vec_ky

def read_linearNonlinearFromInputFile(input_file):
    with open(input_file, 'r') as input_data:
        input_text = input_data.read().replace(' ', '')
        input_text_knob = input_text.split("&physics_flags")[1].split("/")[0] 
        nonlinear = read_booleanInput(input_text_knob, 'nonlinear')  
    return not nonlinear, nonlinear

def read_vmecFileNameFromInputFile(input_file):
    with open(input_file, 'r') as input_data:
        input_text = input_data.read().replace(' ', '')
        if "&vmec_parameters" in input_text:
            input_text_knob = input_text.split("&vmec_parameters")[1].split("/")[0] 
            vmec_filename = read_stringInput(input_text_knob, 'vmec_filename')  
        else:
            vmec_filename = 'wout*.nc'
    return vmec_filename

def read_nspecFromInputFile(input_file):
    with open(input_file, 'r') as input_data:
        input_text = input_data.read().replace(' ', '')
        input_text_knob = input_text.split("&species_knobs")[1].split("/")[0] 
        nspec = read_integerInput(input_text_knob, 'nspec')  
    return nspec

def read_deltFromInputFile(input_file):
    with open(input_file, 'r') as input_data:
        input_text = input_data.read().replace(' ', '')
        input_text_knob = input_text.split("&knobs")[1].split("/")[0] 
        delt = read_floatInput(input_text_knob, 'delt')  
    return delt

def read_svalueFromInputFile(input_file):
    with open(input_file, 'r') as input_data:
        input_text = input_data.read().replace(' ', '')
        if "&vmec_parameters" in input_text:
            input_text_knob = input_text.split("&vmec_parameters")[1].split("/")[0] 
            svalue = read_floatInput(input_text_knob, 'torflux') 
        elif "&millergeo_parameters" in input_text:
            input_text_knob = input_text.split("&millergeo_parameters")[1].split("/")[0] 
            svalue = read_floatInput(input_text_knob, 'rhoc')**2  
    return svalue

def read_tendFromInputFile(input_file):
    with open(input_file, 'r') as input_data:
        input_text = input_data.read().replace(' ', '')
        if 't_end' in input_text:
            input_text_knob = input_text.split("&knobs")[1].split("/")[0] 
            tend = read_floatInput(input_text_knob, 't_end')  
        if 'tend' in input_text:
            input_text_knob = input_text.split("&knobs")[1].split("/")[0] 
            tend = read_floatInput(input_text_knob, 'tend')   
    return tend
    
#===============================================================================
#                READ SPECIFIC DATA TYPES FROM THE INPUT FILE
#===============================================================================

def read_parameterFromIniFile(section, parameter, default_value):
    ''' Read the variables if they exist and convert them to the correct datatype. ''' 
    if default_value == None: return float(section[parameter])
    if isinstance(default_value, str): return section[parameter]
    if isinstance(default_value, bool): return True if section[parameter]=="True" else False
    if isinstance(default_value, int): return int(section[parameter])
    if isinstance(default_value, float): return float(section[parameter])
    return 

#--------------------------------------------
def read_parameterFromInputFile(input_text, parameter, default_value):
    ''' Read the variables if they exist and convert them to the correct datatype. ''' 

    if default_value == None:
        input_value = read_floatInput(input_text, parameter) 
        return input_value if input_value != "USE DEFAULT" else default_value        

    if isinstance(default_value, bool):
        input_value = read_booleanInput(input_text, parameter)
        return input_value if input_value != "USE DEFAULT" else default_value

    if isinstance(default_value, int):
        input_value = read_integerInput(input_text, parameter)
        return input_value if input_value != "USE DEFAULT" else default_value

    if isinstance(default_value, float):
        input_value = read_floatInput(input_text, parameter)
        return input_value if input_value != "USE DEFAULT" else default_value

    if isinstance(default_value, str):
        input_value = read_stringInput(input_text, parameter)
        return input_value if input_value != "USE DEFAULT" else default_value
    return

#--------------------------------------------
def read_integerInput(input_text, variable):
    ''' Read the [variable] from the "*.in" file, if it doesn't exist, return NaN. '''
    if variable not in input_text:
        return "USE DEFAULT"
    try:    
        return int(float(input_text.split(variable)[1].split('=')[1].split('\n')[0]))
    except: 
        print("WARNING: ", variable, "is not an integer so failed to read input file.")
        return np.NaN
    return

#--------------------------------------------
def read_floatInput(input_text, variable):
    ''' Read the [variable] from the "*.in" file, if it doesn't exist, return NaN. '''
    if variable not in input_text: 
        if variable=="delt": print("USE DEFAULT")
        return "USE DEFAULT"  
    try:
        if variable=="delt": # Make sure we don't read delt_max =...
            dummy_text = input_text.split(variable)  
            dummy_text = [i for i in dummy_text if "_max" not in i]
            return float(dummy_text[1].split('=')[1].split('\n')[0].split('!')[0])
        return float(input_text.split(variable)[1].split('=')[1].split('\n')[0].split('!')[0])
    except:  
        if variable!="delt": 
            print("WARNING: ", variable, "is not a float so failed to read input file.")
        return np.NaN
    return

#--------------------------------------------
def read_booleanInput(input_text, variable):
    ''' Read the [variable] from the "*.in" file and convert .false. to False '''
    if variable not in input_text:
        return "USE DEFAULT"
    try:    
        value = input_text.split(variable)[1].split('=')[1].split('\n')[0]
    except: 
        print("WARNING: ", variable, "is not a boolean so failed to read input file.")
        value = np.NaN
    if value == ".false.": value = False
    if value == ".true." : value = True
    if value == "F": value = False
    if value == "T" : value = True
    return value

#--------------------------------------------
def read_stringInput(input_text, variable):
    ''' Read the [variable] from the "*.in" file and encode it to a np.string_ or b'' string '''
    if variable not in input_text:
        return "USE DEFAULT"
    try:    
        value = str(input_text.split(variable)[1].split('=')[1].split('\n')[0])
    except: 
        print("WARNING: ", variable, "is not a string so failed to read input file.")
        value = str(np.NaN)
    value = value.replace(" ","")
    value = value.replace("'","")
    value = value.replace('"','') 
    return value

#--------------------------------------------
if __name__ == "__main__": 
    inputParameters = read_inputParameters("/home/hanne/CIEMAT/RUNS/JUMPERS/LHD_fprim6tprim0/input_ky3.3125.in")
    size = sys.getsizeof(inputParameters)
    for knob in inputParameters.keys():
        size += sys.getsizeof(inputParameters[knob])
        size += sum(map(sys.getsizeof, inputParameters[knob].values())) + sum(map(sys.getsizeof, inputParameters[knob].keys()))
    print(size)
    print(sys.getsizeof(np.ones(1000, np.int8)))
    print(sys.getsizeof(np.ones(1000, np.int16)))
    print(sys.getsizeof(np.ones(1000, np.float16)))
    print(sys.getsizeof(np.ones(1000, np.float32)))
    print(sys.getsizeof(np.ones(1000, np.float64)))
