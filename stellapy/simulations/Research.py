import time
import os, pathlib 
from stellapy.simulations.Experiment import create_experiments 
from stellapy.utils.files.sort_listByNumbers import sort_listByLabel 
from stellapy.utils.files.sort_listByNumbers import sort_listByNumbers 
from stellapy.plot.utils.devices.recognize_device import recognize_device 
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import get_plottingOption, load_labelsLinesMarkers
from stellapy.plot.utils.labels.change_stellaParametersForLabels import change_stellaParametersForLabels

################################################################################
#                              CREATE RESEARCH                                 #
################################################################################

def create_research(\
    # To create the simulations we need their location  
    folders=None, input_files=None, 
    # To group by experiments and simulations we need the following data
    variables=10, knob1="-", key1="-", knob2="-", key2="-", knob3="-", key3="-",
    ignoreResolution=True, folderIsExperiment=False, resolutionScan=False, verbose=False, 
    # Print a progress bar on the terminal
    progress_start=None, progress_stop=None):
    ''' Group multiple experiments so they can be saved as a pickly object and handed to plotting functions. '''

    # Add timer 
    start = time.time()
    
    # Collect the details that define how the research is constructed
    creationDetails = CreationDetails(variables, knob1, key1, knob2, key2, knob3, key3, ignoreResolution, folderIsExperiment, resolutionScan, progress_start, progress_stop)

    # Create the experiments 
    experiments = create_experiments(folders, input_files, creationDetails)
    
    # Look for the total amount of varied values throughout the experiments
    variedValues=[]
    for experiment in experiments: 
        variedValues += experiment.variedValues 
    variedValues = sort_listByNumbers(list(set(variedValues)))  
    
    # Create a research object  
    research = Research(experiments, variedValues, creationDetails) 
    
    # Print how much time this took
    if verbose: print("    ---> It took "+str(time.time() - start)+" seconds to create the research object.")
    return research


#-------------------
class CreationDetails:

    def __init__(self, variables, knob1, key1, knob2, key2, knob3, key3, ignoreResolution, folderIsExperiment, resolutionScan, progress_start, progress_stop):
        self.key1 = key1
        self.key2 = key2
        self.key3 = key3
        self.knob1 = knob1
        self.knob2 = knob2
        self.knob3 = knob3
        self.variables = variables
        self.resolutionScan = resolutionScan
        self.ignoreResolution = ignoreResolution
        self.folderIsExperiment = folderIsExperiment
        self.progress_start = progress_start
        self.progress_stop = progress_stop
        self.progress_stop1 = progress_start + (progress_stop-progress_start)*1/2 if progress_start!=None else None 
        return 

################################################################################
#                             CLASS RESEARCH                                   #
################################################################################
    
class Research:
    ''' Make a research object which groups multiple experiment objects. '''
    
    def __init__(self, experiments, variedValues, creationDetails):
        
        # Save the experiments, the varied values and the creation details
        self.experiments = experiments 
        self.variedValues = variedValues 
        self.creationDetails = creationDetails
    
        # Update the experiment labels and get the input_files, folders and file names for the GUI
        self.update_experimentLabels()       
        self.get_foldersAndNamesForGui()
        self.save_numberOfExperimentsAndSimulations()
        
        # Sort the experiments
        self.experiments = sort_listByNumbers(self.experiments, source="id") if not self.creationDetails.resolutionScan else sort_listByLabel(self.experiments)
     
        
        # Get some default values for the markers and lines for each experiment/simulation
        self.set_labelsLinesMarkers()       
        self.plotting_option = get_plottingOption(self.experiments)  
        if True: return
        
################################################################################
#                                METHODS                                       #
################################################################################

    def save_numberOfExperimentsAndSimulations(self):
        self.numberOfExperiments = 0
        self.numberOfSimulations = 0
        for experiment in self.experiments:
            self.numberOfExperiments += 1
            for _ in experiment.simulations: 
                self.numberOfSimulations += 1
    
    #-----------------------------------   
    def get_foldersAndNamesForGui(self):
        
        # Get the input files
        self.input_files = []
        for experiment in self.experiments:
            for simulation in experiment.simulations: 
                if simulation.nonlinear: self.input_files += [simulation.input_file]
                if simulation.linear: self.input_files += [mode.input_file for mode in simulation.modes]
        self.input_files = sorted(self.input_files)
        
        # Split the input files in folder/input_file pairs
        self.files   = [ i.name   for i in self.input_files ]
        self.folders = [ i.parent for i in self.input_files ] 
        self.e_id    = [ i.name   for i in self.input_files ]
        self.s_id    = [ i.name   for i in self.input_files ]
        
        # Now move the "run01" folders to the file name
        for folder in self.folders:
            if "run" in folder.name and folder.name.replace("run","").isdigit():
                index = self.folders.index(folder)
                self.files[index]   = folder.name+"/"+self.files[index]
                self.folders[index] = folder.parent
                
        # Make sure the folders are strings and reduce folder to the path behind
        if len(self.folders) != 0:
            common_prefix = pathlib.Path(os.path.commonprefix(self.folders))
            common_prefix = str(self.folders[0].parent)+"/" if common_prefix.name in str(self.folders[0]) else common_prefix
            self.folders = [ str(f).replace(str(common_prefix), "") for f in self.folders]
            
            # Add the experiment and simulation id to the files
            for experiment in self.experiments:
                for simulation in experiment.simulations:
                    if simulation.nonlinear:  
                        index_i = self.input_files.index(simulation.input_file)
                        self.e_id[index_i] = experiment.id
                        self.s_id[index_i] = simulation.marker_label.replace("$","").replace("\,"," ")
                    if simulation.linear: 
                        for mode in simulation.modes:
                            index_i = self.input_files.index(mode.input_file)
                            self.e_id[index_i] = experiment.id
                            self.s_id[index_i] = simulation.marker_label.replace("$","").replace("\,"," ")
                        
        # Now get the unique folders
        self.unique_folders = list(set(self.folders))
        return 
                  
    #-----------------------------------   
    def update_experimentLabels(self):
        
        # Initiate the new experiment id's
        experiment_ids = []

        # Name the experiments based on the selected stella knob and stella key
        if not self.creationDetails.resolutionScan and self.creationDetails.knob1!="-":
            for experiment in self.experiments: 
                
                # Get the first simulation and its input parameters
                simulation = experiment.simulations[0]
                if simulation.linear: inputobject = simulation.modes[0].input
                if simulation.nonlinear: inputobject = simulation.input
                    
                # Update the experiment key and value to a math notation for the line labels for each experiment 
                if self.creationDetails.key1!='vmec_filename':
                    value = inputobject.inputParameters[self.creationDetails.knob1][self.creationDetails.key1]  
                    variable, value = change_stellaParametersForLabels(simulation, self.creationDetails.knob1, self.creationDetails.key1, value)
                    experiment_id = variable + "$\,=\,$" + str(value)
                if self.creationDetails.key1=='vmec_filename': 
                    experiment_id = str(recognize_device(inputobject.vmec_filename))  
                
                # If two experiment_knobs were selected, add the second experiment_key and value
                if self.creationDetails.knob2!="-":
                    if self.creationDetails.key2!='vmec_filename':
                        value = inputobject.inputParameters[self.creationDetails.knob2][self.creationDetails.key2]   
                        variable, value = change_stellaParametersForLabels(simulation, self.creationDetails.knob2, self.creationDetails.key2, value)
                        experiment_id += "; " + variable + "$\,=\,$" + str(value) 
                    if self.creationDetails.key2=='vmec_filename': 
                        experiment_id += "; " + str(recognize_device(inputobject.vmec_filename))  
                        
                # Set the experiment id to "variable = value" or "variable1 = value1; variable2 = value2"
                experiment_ids.append(experiment_id)

        # When performing a resolution scan, show the resolution parameters
        if self.creationDetails.resolutionScan:
            for experiment in self.experiments: 
                experiment_ids.append(experiment.variedValues[0]+" ") 
            
        # Only use the previously created experiment ids if they are unique for each experiment  
        if len(experiment_ids) != 0 and len(experiment_ids)==len(list(set(experiment_ids))): 
            for experiment in self.experiments:
                experiment.id = experiment_ids[self.experiments.index(experiment)]
                
        # If the folder is the experiment use the folder name
        else:
            experiment_ids = []
            for experiment in self.experiments:
                if experiment.simulations[0].nonlinear:
                    experiment_ids.append(str(experiment.simulations[0].input_file).split("/")[-2]) 
                if experiment.simulations[0].linear:
                    experiment_ids.append(str(experiment.simulations[0].modes[0].input_file).split("/")[-2]) 
            if len(experiment_ids) != 0 and len(experiment_ids)==len(list(set(experiment_ids))): 
                for experiment in self.experiments:
                    experiment.id = experiment_ids[self.experiments.index(experiment)]
        return 
    
    #-----------------------------------   
    def print_research(self):
        
        # Print some information
        print("")
        print("     ","".center(110,"#")) 
        print("     "," RESEARCH ".center(110,"#"))
        print("     ","".center(110,"#")) 
        
        for experiment in self.experiments:
            print("")
            print("      Experiment "+str(self.experiments.index(experiment)+1)+": "+experiment.id)
            print("      --------------------------------")
            if self.creationDetails.variables==-1: print("        Each simulation is assumed to be its own experiment at the following radial position:  ")
            if self.creationDetails.variables==1:  print("        There is 1 varied parameter with the following values:  ")
            if self.creationDetails.variables>1:   print("        There are "+str(self.creationDetails.variables)+" varied parameters with the following values:  ")
            for varied_value in experiment.variedValues:
                print("               "+varied_value)
            print("        The simulations associated to this experiment are:  ")
            for simulation in experiment.simulations:
                print("               "+'"'+simulation.id)
            print("")
             
        print("     ","".center(110,"#")) 
        if True: return

#===============================================================================
#                   LINES/MARKERS STYLE FOR EACH SIMULATION                    #
#=============================================================================== 

    def set_labelsLinesMarkers(self):
        
        # Get the ids of the experiments
        experiment_ids = [experiment.id for experiment in self.experiments]
               
        # Load the default style ased on the experiment ids and varied values
        line_label, marker_label, line_style, line_color, marker_style, marker_color = load_labelsLinesMarkers(experiment_ids, self.variedValues)
        
        # Change the lines/markers of the underlying experiment objects based on self.variedValues
        for experiment in self.experiments: experiment.set_labelsLinesMarkers(line_label, marker_label, line_style, line_color, marker_style, marker_color)
        if True: return

    
    