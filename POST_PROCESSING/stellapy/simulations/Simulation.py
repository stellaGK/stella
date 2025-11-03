
# External modules 
import sys, pathlib 
from stellapy.utils.decorators.exit_program import exit_program  

# Functions to create the simulations 
from stellapy.utils.files import get_filesInFolder 
from stellapy.utils.files.keep_simulationsWithOutputFiles import keep_simulationsWithOutputFiles
from stellapy.simulations.utils.get_simulationIdentifier import get_simulationsIdentifier 

# Input file
from stellapy.data.input.read_inputFile import read_inputFile, read_numberOfModesFromInputFile 

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
from stellapy.data.input.write_listOfMatchingInputFiles import read_inputFilesInDummyInputFiles, read_inputFilesInDummyInputFile
from stellapy.simulations.utils.get_inputFilesWhoOnlyDifferInResolution import get_inputFilesWhoOnlyDifferInResolution
 
#===============================================================================
#                             CREATE SIMULATIONS
#===============================================================================

def create_simulations(folders=None, input_files=None,  
        ignore_resolution=True, include_knobs=[], ignore_oldFiles=True): 
    
    # Initiate a vector of simulations
    simulations = []
    
    # Get the input files defined by <folders> or <input_files>
    input_files = get_inputFiles(folders, input_files, ignore_oldFiles) 
    
    # With <restartModes> a file '_ky1.0.in' gets restarted as '_ky1.0_kx0.0.in'
    input_files = [i for i in input_files if pathlib.Path(str(i).replace('.in', '_kx0.0.in')) not in input_files]
    
    # For linear flux tube simulations, read the dummy input file instead 
    # Remove the input files in the dummy inputs from the list of input files  
    dummy_input_files = get_filesInFolder(folders, end='_dummy.in') if folders!=None else [i for i in input_files if '_dummy.in' in str(i)]
    input_files_in_dummy_files = read_inputFilesInDummyInputFiles(dummy_input_files)
    input_files = [i for i in input_files if i not in input_files_in_dummy_files] 
    input_files += dummy_input_files 
    
    # Read the input parameters of each input files
    input_parameters = [read_inputFile(input_file) for input_file in input_files]
    
    # Combine input files from linear simulations who are identical except for their resolution
    # This allows to create spectra gamma(ky) with modes that have a smaller delta t at small ky, and 
    # larger nmu for large ky, which is generally necesarry to obtain a converged spectrum 
    input_files, input_parameters, similar_input_files = get_inputFilesWhoOnlyDifferInResolution(input_files, input_parameters, ignore_resolution=ignore_resolution, include_knobs=include_knobs)

    # For each input file create a simulation object 
    for i in range(len(input_files)):    
        simulations.append(Simulation(input_files[i], input_parameters[i], similar_input_files[i]))  
    return simulations  
 
#===============================================================================
#                               CLASS SIMULATION                               #
#===============================================================================
 
class Simulation:
    ''' Saves all the information from a simulation run with stella. '''
     
    def __init__(self, input_file, input_parameters, similar_input_files):
         
        # Give a unique ID to the simulation
        self.id = get_simulationsIdentifier(input_file)
        
        # Save some relevant info 
        self.input_file = input_file 
        self.input_parameters = self.inputParameters = input_parameters 
        self.linear = not self.input_parameters['gyrokinetic_terms']['include_nonlinear']
        self.nonlinear = self.input_parameters['gyrokinetic_terms']['include_nonlinear']
        self.full_flux_surface = self.input_parameters['gyrokinetic_terms']['include_full_flux_annulus']
        self.nakxnaky = self.input_parameters['kxky_grid_range']['nakx']*self.input_parameters['kxky_grid_range']['naky'] if self.input_parameters['kxky_grid_option']['grid_option']=='range' else self.input_parameters['kxky_grid_box']['nx']*self.input_parameters['kxky_grid_box']['ny']
        
        # Remove some input parameters that may vary between simulations without affecting the simulation (much)
        del self.input_parameters['time_step']['delt_option']
        del self.input_parameters['time_step']['delt_max']
        del self.input_parameters['time_step']['delt_min']
        del self.input_parameters['parallelisation']['mat_gen']
        del self.input_parameters['time_trace_options']['nstep']
        del self.input_parameters['time_trace_options']['tend']
        del self.input_parameters['time_trace_options']['avail_cpu_time']
        del self.input_parameters['restart_options']['restart_file']
        del self.input_parameters['restart_options']['restart_dir']
        del self.input_parameters['initialise_distribution']['phiinit']
        del self.input_parameters['initialise_distribution']['initialise_distribution_option']
        del self.input_parameters['diagnostics']
        del self.input_parameters['diagnostics_potential']
        del self.input_parameters['diagnostics_omega']
        del self.input_parameters['diagnostics_distribution']
        del self.input_parameters['diagnostics_fluxes']
        del self.input_parameters['diagnostics_moments']
        del self.input_parameters['kxky_grid_range']
        if self.input_parameters['gyrokinetic_terms']['include_nonlinear']==True: del self.input_parameters['time_step']['delt'] 
        
        # For linear flux tube simulations, add the input files in the dummy input or merge similar simulations
        self.multiple_input_files = False
        if '_dummy.in' in str(self.input_file):
            self.input_files = read_inputFilesInDummyInputFile(self.input_file)
            self.input_files += [self.input_file]
            self.multiple_input_files = True  
        if isinstance(similar_input_files, list):
            if len(similar_input_files)>1: 
                self.input_files = similar_input_files 
                self.multiple_input_files = True 
                for input_file in self.input_files: 
                    if '_dummy.in' in str(input_file):
                        self.input_files += read_inputFilesInDummyInputFile(input_file) 
                
        # If we have multiple input files, count the number of modes
        if self.multiple_input_files:
            self.nakxnaky = 0
            for input_file in self.input_files: 
                if '_dummy.in' not in str(input_file):
                    nakx, naky = read_numberOfModesFromInputFile(input_file) 
                    self.nakxnaky += nakx*naky 
           
        # Get some default values for the markers and lines  
        self.set_labelsLinesMarkers()  
        
        # The GUI will asign Progress objects
        self.Progress = None 

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
         
        # Prevent the loading of the plotting module if we don't need it
        from stellapy.plot.utils.style.get_styleForLinesAndMarkers import load_labelsLinesMarkers
         
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
            
          
#===============================================================================
#                              GET THE INPUT FILES                             #
#===============================================================================
        
def get_inputFiles(folders, input_files, ignore_oldFiles):

    # Make sure we have lists of PosixPath objects
    if isinstance(folders, str): folders = [folders]
    if isinstance(input_files, str): input_files = [input_files]
    if isinstance(folders, pathlib.PurePath): folders = [folders]
    if isinstance(input_files, pathlib.PurePath): input_files = [input_files]
    if input_files!=None and not isinstance(input_files[0], pathlib.PurePath):
        input_files = [pathlib.Path(i) for i in input_files]
    if folders!=None and not isinstance(folders[0], pathlib.PurePath):
        folders = [pathlib.Path(f) for f in folders]
           
    # If <input_files>=None read all input_files in <folders>
    if input_files==None:   
        input_files = get_filesInFolder(folders, end='.in')
             
    # If <folders> and <input_files> are given they match
    elif folders!=None and input_files!=None:
        input_files = [folders[i] / input_files[i] for i in range(len(folders))]
    
    # Only look at input files that have a output files
    input_files = keep_simulationsWithOutputFiles(input_files)
    if input_files == []: 
        exit_reason = 'No netcdf files were found inside the following folders:\n'
        for folder in folders: exit_reason += str(folder)
        exit_program(exit_reason, create_simulations, sys._getframe().f_lineno)
    
    # Ignore simulations in and 'OLD' folder
    if ignore_oldFiles: input_files = [i for i in input_files if 'OLD' not in str(i)]
    
    # With <restartModes> a file '_ky1.0.in' gets mistakingly restarted as '_ky1.0_kx0.0.in'
    input_files = [i for i in input_files if pathlib.Path(str(i).replace('.in', '_kx0.0.in')) not in input_files]
    
    # Sort the input files
    input_files = sorted(input_files)
    return input_files
  
    
