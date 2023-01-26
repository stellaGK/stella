 
import os, sys
import pathlib, configparser
from stellapy.data.input.read_inputFile import  read_extraInputParameters
from stellapy.data.input.get_gridParameters import get_gridParameters
from stellapy.data.input.get_modeParameters import get_modeParameters
from stellapy.data.input.get_inputParameters import get_inputParameters
from stellapy.data.input.get_basicParameters import get_basicParameters 
from stellapy.data.stella.load_defaultInputParameters import load_defaultInputParameters 
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime
from stellapy.utils.decorators.exit_program import exit_program

#===============================================================================
#                        CREATE THE INPUT OBJECT
#===============================================================================

class Input:
    
    # Copy the data from <simulation> that is needed to construct <input>
    def __init__(self, simulation):
        
        # Remember the path object of <simulation> or <mode>  
        self.path = simulation.path
        return
        
    # Load all the data for the input object
    def load_data(self):
        get_inputParameters(self)
        get_basicParameters(self)
        get_gridParameters(self)
        get_modeParameters(self)
        return 
    
    # Read all input parameters
    @calculate_attributeWhenReadFirstTime
    def inputParameters(self):  get_inputParameters(self);    return self.inputParameters 
    
    # Read basic input parameters
    @calculate_attributeWhenReadFirstTime
    def nspec(self):            get_basicParameters(self);    return self.nspec
    @calculate_attributeWhenReadFirstTime
    def linear(self):           get_basicParameters(self);    return self.linear
    @calculate_attributeWhenReadFirstTime
    def nonlinear(self):        get_basicParameters(self);    return self.nonlinear
         
    # Read basic input parameters for the equilibrium
    @calculate_attributeWhenReadFirstTime
    def y0(self):               get_basicParameters(self);    return self.y0
    @calculate_attributeWhenReadFirstTime
    def rho(self):              get_basicParameters(self);    return self.rho
    @calculate_attributeWhenReadFirstTime
    def svalue(self):           get_basicParameters(self);    return self.svalue
    @calculate_attributeWhenReadFirstTime
    def vmec(self):             get_basicParameters(self);    return self.vmec 
    @calculate_attributeWhenReadFirstTime
    def miller(self):           get_basicParameters(self);    return self.miller
    @calculate_attributeWhenReadFirstTime
    def vmec_filename(self):    get_basicParameters(self);    return self.vmec_filename
    @calculate_attributeWhenReadFirstTime
    def nperiod(self):          get_basicParameters(self);    return self.nperiod
    @calculate_attributeWhenReadFirstTime
    def nfield_periods(self):   get_basicParameters(self);    return self.nfield_periods
    
    # Read the grid parameters
    @calculate_attributeWhenReadFirstTime
    def nzed(self):             get_gridParameters(self);     return self.nzed
    @calculate_attributeWhenReadFirstTime
    def nzgrid(self):           get_gridParameters(self);     return self.nzgrid
    @calculate_attributeWhenReadFirstTime
    def nvgrid(self):           get_gridParameters(self);     return self.nvgrid
    @calculate_attributeWhenReadFirstTime
    def nmu(self):              get_gridParameters(self);     return self.nmu  
    
    # Read the nonlinear mode parameters
    @calculate_attributeWhenReadFirstTime
    def nx(self):               get_modeParameters(self);     return self.nx
    @calculate_attributeWhenReadFirstTime
    def ny(self):               get_modeParameters(self);     return self.ny
    @calculate_attributeWhenReadFirstTime
    def nakx(self):             get_modeParameters(self);     return self.nakx
    @calculate_attributeWhenReadFirstTime
    def naky(self):             get_modeParameters(self);     return self.naky  
    
#     # Read the linear mode parameters
#     @calculate_attributeWhenReadFirstTime 
#     def kx(self):               get_modeParameters(self);     return self.kx
#     @calculate_attributeWhenReadFirstTime
#     def ky(self):               get_modeParameters(self);     return self.ky   

#===============================================================================
#                                  LOAD INPUT                                  #
#===============================================================================

def load_inputObject(self):     
    
    # Add the input for a simulation
    if self.object=="Simulation":
        self.input = Input(self)
        return
        
    # For a linear simulation, add the input per mode. Since most inputs are
    # identical, only read each unique input once (self==mode)
    if self.object=="Mode":
        load_onlyReferenceInputs(self)
        return 

#----------------------------
def load_onlyReferenceInputs(self):
    """ Load an input object for each unique input. """
    
    # Initialize 
    loaded_inputs={}
    
    # Only read one input per unique input
    for imode, mode in enumerate(self.simulation.modes): 
         
        # If the input is already read, create a reference to the existing input
        if mode.path.input in self.simulation.path.loaded_inputs: 
            mode.input = self.simulation.modes[loaded_inputs[mode.path.input]].input
        if mode.path.input not in self.simulation.path.loaded_inputs: 
            self.simulation.path.loaded_inputs.append(mode.path.input) 
            loaded_inputs[mode.path.input] = imode
            mode.input = Input(mode) 
    return  

#-----------------------------------
def read_uniqueInputs(path):  
    
    # Make sure we have the overview of unique inputs  
    if not os.path.isfile(path.folder / (path.name + ".list.inputs.ini")): 
        from stellapy.data.input.write_iniFileForInputs import write_iniFileForInputs 
        write_iniFileForInputs(path.folder)
        
    # Read the overview of unique inputs
    list_of_reference_inputs = configparser.ConfigParser()
    list_of_reference_inputs.read(path.folder / (path.name + ".list.inputs.ini")) 

    # Get the identifier of each reference input: "Input X"
    ids_of_reference_inputs = sorted([ i for i in list_of_reference_inputs.keys() if "Input " in i and "Unique" not in i])   
    return list_of_reference_inputs, ids_of_reference_inputs

#-----------------------------------
def find_referenceInput(path, mode, list_of_reference_inputs, ids_of_reference_inputs):
    
    # Initiate
    reference_input_path = None
    
    # If it is written on marconi, translate to local paths
    for i in ids_of_reference_inputs:
        for key in list_of_reference_inputs[i].keys():
            if 'marconi' in list_of_reference_inputs[i][key]:  
                list_of_reference_inputs[i][key] = str(path.folder)+"/"+(list_of_reference_inputs[i][key].split("/"+path.folder.name+"/")[-1])
            elif '/mnt/lustre/' in list_of_reference_inputs[i][key]:   
                list_of_reference_inputs[i][key] = str(path.folder)+"/"+(list_of_reference_inputs[i][key].split("/"+path.folder.name+"/")[-1])
            else:
                break  
            
    # Find the reference input for <mode> as well as (kx,ky) 
    for reference_input_id in ids_of_reference_inputs:   
        if str(mode.input_file) in list_of_reference_inputs[reference_input_id].values():  
            mode_name = list(list_of_reference_inputs[reference_input_id].values()).index(str(mode.input_file))
            mode_name = list(list_of_reference_inputs[reference_input_id].keys())[mode_name]
            #mode.kx = float(mode_name.split("(")[-1].split(",")[0])
            #mode.ky = float(mode_name.split(",")[-1].split(")")[0])   
            reference_input_path = path.folder / (path.name + ".unique.input" + reference_input_id.split("Input ")[-1]+".ini") 
            break

    # Return the reference input and the path to the file  
    exit_reason = "Couldn't find "+str(mode.input_file)+" in the list of input files. \n"
    exit_reason += "List of input files: "+str(path.folder / (path.name + ".list.inputs.ini"))
    if reference_input_path==None: exit_program(exit_reason, find_referenceInput, sys._getframe().f_lineno)
    return reference_input_path

#===============================================================================
#                             SAVE THE INPUT FILE                              #
#===============================================================================

def save_inputFile(mode):
     
    # Define the configuration file
    ini_path = mode.path.input 
    ini_data = configparser.ConfigParser()
    
    # Make sure we have all the input parameters
    mode.inputParameters = read_extraInputParameters(mode.inputParameters, mode.input_file) 
    
    # Get the input parameters
    input_parameters = mode.inputParameters 
    default_parameters = load_defaultInputParameters() 

    # Look at the knobs in the following order
    knobs = list(input_parameters.keys()) 
    knobs = ["="*15+" SPECIES "+"="*15] + sorted([s for s in knobs if "species_" in s]) + ["parameters"] 
    knobs += ["="*15+" GEOMETRY "+"="*15] + ["geo_knobs", "millergeo_parameters", "vmec_parameters"]
    knobs += ["="*15+" GRIDS "+"="*15] + ["kt_grids_box_parameters", "zgrid_parameters", "vpamu_grids_parameters"]
    knobs += ["="*15+" SIMULATION "+"="*15] + ["knobs", "stella_diagnostics_knobs"]
    knobs += ["="*15+" OTHER "+"="*15] + ["time_advance_knobs", "kt_grids_knobs", "dist_fn_knobs", "physics_flags", "init_g_knobs", "layouts_knobs", "neoclassical_input", "sfincs_input"]
    
    # Only save the variables that differ from the default parameters 
    for knob in knobs: 
        if "====" in knob: ini_data.add_section(knob) 
        if "====" not in knob:
            for key in sorted(list(input_parameters[knob].keys())):
                if key!="teti":
                    if default_parameters[knob][key]!=input_parameters[knob][key]: 
                        if knob not in ini_data: ini_data.add_section(knob) 
                        ini_data[knob][key] = str(input_parameters[knob][key])
                
    # Save the configuration file
    with open(ini_path, 'w') as ini_file:
        ini_data.write(ini_file) 
        print("    ----> Saved the input file as " + ini_path.name)
    
    return ini_path


################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    
    from stellapy.simulations.Simulation import create_simulations
    import copy, timeit; start = timeit.timeit()
    
    folder = pathlib.Path("/home/hanne/CIEMAT/PREVIOUSRUNS/LINEARMAPS/W7Xstandard_rho0.7_aLTe0/ResolutionScan/fprim4tprim4_ky1.5/nzed") 
    folder = pathlib.Path("/home/hanne/CIEMAT/PREVIOUSRUNS/LINEARMAPS/W7Xstandard_rho0.7_aLTe0/ResolutionScan/fprim4tprim4_ky1.5/nmu") 
    folder = pathlib.Path("/home/hanne/CIEMAT/PREVIOUSRUNS/LINEARMAPS/W7Xstandard_rho0.7_aLTe0/LinearMap/fprim4tprim4")   
    simulations = create_simulations(folders=folder, input_files=None, ignore_resolution=True, number_variedVariables=5) 
    print("\nWe have "+str(len(simulations))+" simulations. Test whether we have") 
    print("references to the data instead of copies.")
    for simulation in simulations:
        print("\nSimulation:", simulation.id)  
        input_test = copy.deepcopy(simulation.modes[0].input) 
        para_test = copy.deepcopy(simulation.modes[0].input.inputParameters) 
        nfp_test = copy.deepcopy(simulation.modes[0].input.nfield_periods) 
        for mode in simulation.modes:  
            inputt = mode.input  
            para = mode.input.inputParameters
            nfp = mode.input.nfield_periods 
            reference_test1f = inputt is input_test 
            reference_test2f = para is para_test 
            reference_test3f = nfp is nfp_test 
            reference_test1 = inputt is simulation.modes[0].input
            reference_test2 = para is simulation.modes[0].input.inputParameters
            reference_test3 = nfp is simulation.modes[0].input.nfield_periods
            name = "("+str(mode.kx)+", "+str(mode.ky)+")" 
            turns = mode.input.inputParameters["vmec_parameters"]["poloidal_turns"]
            delt = mode.input.inputParameters["knobs"]["delt"]
            print("{:<15}".format(name), delt, turns, " References:", reference_test1, reference_test2, reference_test3,\
                                             "   TEST:", reference_test1f, reference_test2f, reference_test3f) 


