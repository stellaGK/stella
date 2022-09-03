
# External modules 
import sys, pathlib 
from stellapy.utils.decorators.exit_program import exit_program 
from stellapy.utils.commandprompt.print_progressbar import print_progressbar

# Functions to create the simulations 
from stellapy.utils.files import get_filesInFolder 
from stellapy.utils.files.remove_simulationsWithoutOutncFile import remove_simulationsWithoutOutncFile
from stellapy.simulations.utils.get_simulationIdentifier import get_simulationsIdentifier 

# Input file
from stellapy.data.input.read_inputFile import read_inputFile 

# Plotting style
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import load_labelsLinesMarkers

# Modes
from stellapy.simulations.Modes import create_modes, finish_initializationOfModes

# Data
from stellapy.data.potential import load_potentialObject
from stellapy.data.time.load_timeObject import load_timeObject
from stellapy.data.paths.load_pathObject import load_pathObject
from stellapy.data.omega.load_omegaObject import load_omegaObject
from stellapy.data.input.load_inputObject import load_inputObject
from stellapy.data.fluxes.load_fluxesObject import load_fluxesObject
from stellapy.data.moments.load_momentsObject import load_momentsObject
from stellapy.data.geometry.load_geometryObject import load_geometryObject 
from stellapy.data.dimensions.load_vectorsObject import load_vectorsObject
from stellapy.data.saturated.load_saturatedObject import load_saturatedObject
from stellapy.data.dimensions.load_dimensionsObject import load_dimensionsObject
from stellapy.data.lineardata.load_lineardataObject import load_linearDataObject
from stellapy.data.referenceunits.load_referenceObject import load_referenceObject
from stellapy.data.distribution.load_distributionObject import load_distributionObject
from stellapy.simulations.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime

 
#################################################################
#                        CREATE SIMULATIONS
#################################################################

def create_simulations(folders=None, input_files=None, ignore_resolution=True,
        number_variedVariables=10, folderIsExperiment=False, write_uniqueFiles=False,
        # Print a progress bar on the terminal
        progress_start=None, progress_stop=None):
    """ Will take multiple folders and input files and will group the input files
    together that only differ in (kx,ky) since they are in essence the same simulation. """ 

    # Make sure we have lists of PosixPath objects
    if isinstance(folders, str): folders = [folders]
    if isinstance(input_files, str): input_files = [input_files]
    if isinstance(folders, pathlib.PurePath): folders = [folders]
    if isinstance(input_files, pathlib.PurePath): input_files = [input_files]
    if input_files is not None and not isinstance(input_files[0], pathlib.PurePath):
        input_files = [pathlib.Path(i) for i in input_files]
    if folders is not None and not isinstance(folders[0], pathlib.PurePath):
        folders = [pathlib.Path(f) for f in folders]
           
    # If <input_files>=None read all input_files in <folders>
    if input_files==None:   
        input_files = get_filesInFolder(folders, end=".in") 
             
    # If <folders> and <input_files> are given they match
    elif folders!=None and input_files!=None:
        input_files = [folders[i] / input_files[i] for i in range(len(folders))]
    
    # Only look at input files that have a netcdf output file   
    input_files = remove_simulationsWithoutOutncFile(input_files)  
    if input_files == []: exit_program("No netcdf files were found.", create_simulations, sys._getframe().f_lineno)   
    
    # Sort the input files
    input_files = sorted(input_files)
    
    # If we ignore the resolution, don't look at inputfiles in an "OLD" folder
    # since these are low resolution simulations that have been corrected
    if ignore_resolution: input_files = [i for i in input_files if "OLD" not in str(i)]
    
    # With <restartModes> a file "_ky1.0.in" gets restarted as "_ky1.0_kx0.0.in"
    if not write_uniqueFiles: input_files = [i for i in input_files if pathlib.Path(str(i).replace(".in", "_kx0.0.in")) not in input_files]
        
    # Group the simulations together and save them as simulations[id] = [input_files] 
    simulation_ids = {} 
    
    # Go through the PosixPaths of the input files 
    for i, input_file in enumerate(input_files):  
         
        # Print the progress
        if progress_start!=None: print_progressbar(progress_start+(progress_stop-progress_start)*i/len(input_files)*1/2, 100, prefix = '   Progress:', suffix = 'Creating research object: input files.', length=50)
        
        # Get the simulation identifier and the input parameters
        simulation_id   = get_simulationsIdentifier(input_file) 
        inputParameters = read_inputFile(input_file) 
        
        # The (kx,ky) values are allowed to be different in the "same" simulation
        del inputParameters['kt_grids_range_parameters']
        
        # Remove other parameters that are allowed to differ
        del inputParameters['knobs']['mat_gen']
        del inputParameters['knobs']['nstep']
        del inputParameters['knobs']['tend']
        del inputParameters['knobs']['delt_option']
        del inputParameters['knobs']['delt_max']
        del inputParameters['knobs']['delt_min']
        del inputParameters['init_g_knobs']['phiinit']
        del inputParameters['init_g_knobs']['restart_file']
        del inputParameters['init_g_knobs']['restart_dir']
        del inputParameters['init_g_knobs']['ginit_option'] 
        del inputParameters['stella_diagnostics_knobs']
        if inputParameters['physics_flags']['nonlinear']==True: del inputParameters['knobs']['delt']
         
        # If ignore_resolution = False: the resolution can differ because convergence studies are performed
        # If ignore_resolution = True: the resolution can differ to make sure each mode is converged 
        if not ignore_resolution:
            del inputParameters['zgrid_parameters']['nz']
        if ignore_resolution: 
            del inputParameters['zgrid_parameters']['nzed']
            del inputParameters['vmec_parameters']['nfield_periods']
            del inputParameters['vpamu_grids_parameters']['nvgrid']
            del inputParameters['vpamu_grids_parameters']['nmu'] 
            if inputParameters['physics_flags']['nonlinear']==False: del inputParameters['knobs']['delt']
        
        # Assume we have a new simulation and it is in a new folder
        newInputFileIsUnique = True
        newFolder = {}
        for simulation in simulation_ids.keys():
            newFolder[simulation]=False
        
        # We can sort by "1 folder = 1 experiment"
        if number_variedVariables<0 or folderIsExperiment==True:
            for simulation in simulation_ids.keys():
                newFolder[simulation]=True
            parent_newInput = input_file.parent
            if ("run" in parent_newInput.name) and parent_newInput.name.split("run")[-1].isdigit():
                parent_newInput = parent_newInput.parent
            for simulation in simulation_ids.keys():
                old_input = simulation_ids[simulation]['input files'][0]
                parent_oldInput = old_input.parent
                if ("run" in parent_oldInput.name) and parent_oldInput.name.split("run")[-1].isdigit():
                    parent_oldInput = parent_oldInput.parent
                if parent_newInput == parent_oldInput:
                    newFolder[simulation] = False
            
        # Add the input file to the simulation if it has the same input parameters.
        # If the input parameters are different, create a new simulation for this input file
        if len(list(simulation_ids.keys()))!=0:
            for simulation in simulation_ids.keys():
                if newFolder[simulation]==False and simulation_ids[simulation]['input parameters'] == inputParameters:  
                    newInputFileIsUnique = False
                    simulation_ids[simulation]['input files'].append(input_file)
        if newInputFileIsUnique==True:
            simulation_ids[simulation_id] = {'input files' : [input_file], 'input parameters' : inputParameters}
     
    # For each simulation create a simulation object
    simulations = []; length = len(simulation_ids.keys())
    for i, simulation_id in enumerate(simulation_ids.keys()):  
        if progress_start!=None: print_progressbar(progress_start+(progress_stop-progress_start)*(1+i/length)*1/2, 100, prefix = '   Progress:', suffix = 'Creating research object: simulations.', length=50)
        input_files = simulation_ids[simulation_id]['input files']
        inputParameters = simulation_ids[simulation_id]['input parameters']  
        simulations.append(Simulation(simulation_id, input_files, inputParameters, write_uniqueFiles)) 
   
    # Let a <simulation> see the overall <simulations>
    for simulation in simulations: simulation.simulations = simulations 
    return simulations 

    # Dont collapse header on the next line
    if True: return
 
################################################################################
#                               CLASS SIMULATION                               #
################################################################################
 
class Simulation:
    ''' Saves all the information from a simulation run with stella.
    
    Notes
    -----
    A simulation is a combination of input_files that are identical but have different modes, 
    or that are identical but are restarted. These files need to be in the same folder!
        
    '''
     
    def __init__(self, simulation_id, input_files, inputParameters, write_uniqueFiles, Progress=None):
         
        # Save the information that defines the simulation
        self.id = simulation_id                
        self.object = "Simulation"                    
        self.inputParameters = inputParameters  
        self.formula_quasiLinearFluxes = 'per'
        
        # The data is read per <simulation> or <mode> depending on linear/nonlinear
        self.linear = not self.inputParameters['physics_flags']['nonlinear'] 
        self.nonlinear = self.inputParameters['physics_flags']['nonlinear'] 
        
        # Print information on the GUI (set from the GUI itself) 
        self.Progress = Progress
        
        # For a nonlinear simulation, we can only have one input file 
        if self.nonlinear and len(input_files)>1:  
            exit_reason = "\nNonlinear simulation with multiple input files:" 
            for i in input_files: exit_reason+="\n     "+str(i)  
            exit_program(exit_reason, Simulation, sys._getframe().f_lineno) 
        
        # For a nonlinear simulation use <simulation> to access the data 
        if self.nonlinear: 
            self.input_file = input_files[0]   
        
        # For a linear simulation use <mode> to access the data
        if self.linear:  
            self.modes = create_modes(input_files, self, write_uniqueFiles) 
            finish_initializationOfModes(self)   

        # Get some default values for the markers and lines 
        self.set_labelsLinesMarkers()  

        # Dont collapse header on the next line 
        if True: return 

#===============================================================================
#                         READ THE DATA THROUGH CLASSES                        #
#=============================================================================== 
 
    @calculate_attributeWhenReadFirstTime
    def dim(self):          load_dimensionsObject(self);    return self.dim  
    @calculate_attributeWhenReadFirstTime
    def vec(self):          load_vectorsObject(self);       return self.vec  
    @calculate_attributeWhenReadFirstTime
    def path(self):         load_pathObject(self);          return self.path
    @calculate_attributeWhenReadFirstTime
    def time(self):         load_timeObject(self);          return self.time
    @calculate_attributeWhenReadFirstTime
    def input(self):        load_inputObject(self);         return self.input  
    @calculate_attributeWhenReadFirstTime
    def omega(self):        load_omegaObject(self);         return self.omega
    @calculate_attributeWhenReadFirstTime
    def fluxes(self):       load_fluxesObject(self);        return self.fluxes
    @calculate_attributeWhenReadFirstTime
    def moments(self):      load_momentsObject(self);        return self.moments
    @calculate_attributeWhenReadFirstTime 
    def geometry(self):     load_geometryObject(self);      return self.geometry  
    @calculate_attributeWhenReadFirstTime 
    def saturated(self):    load_saturatedObject(self);     return self.saturated  
    @calculate_attributeWhenReadFirstTime 
    def potential(self):    load_potentialObject(self);     return self.potential  
    @calculate_attributeWhenReadFirstTime 
    def lineardata(self):   load_linearDataObject(self);    return self.lineardata  
    @calculate_attributeWhenReadFirstTime 
    def distribution(self): load_distributionObject(self);  return self.distribution   
    @calculate_attributeWhenReadFirstTime 
    def referenceunits(self):load_referenceObject(self);    return self.referenceunits   
 
#===============================================================================
#                   LINES/MARKERS STYLE FOR EACH SIMULATION                    #
#=============================================================================== 
 
    def set_labelsLinesMarkers(self, line_label=None, marker_label=None, line_style=None, line_color=None, marker_style=None, marker_color=None):
         
        # If the styles are given apply them
        if line_label is not None:
            self.line_label = line_label
            self.marker_label = marker_label
            self.line_style = line_style
            self.line_color = line_color
            self.marker_style = marker_style
            self.marker_color = marker_color
         
        # Otherise load some basic styles
        else: 
            self.line_label,\
            self.marker_label,\
            self.line_style,\
            self.line_color,\
            self.marker_style,\
            self.marker_color = load_labelsLinesMarkers([self.id], [self.id]) 
  
    