
import numpy as np 
from itertools import takewhile 
from stellapy.simulations.Simulation import create_simulations
from stellapy.utils.files.sort_listByNumbers import sort_listByNumbers
from stellapy.utils.commandprompt.print_progressbar import print_progressbar
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import load_labelsLinesMarkers
from stellapy.simulations.utils.get_differenceInDictionaries import get_differenceInDictionaries 
from stellapy.plot.utils.labels.change_stellaParametersForLabels import change_stellaParametersForLabels
from stellapy.simulations.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime

################################################################################
#                             CREATE EXPERIMENTS                               #
################################################################################

def create_experiments(folders, input_files, creationDetails):
    ''' Divide the simulations in experiments, based on the number of varied values
    or based on a parameter that differs for each experiment. For example, when we 
    scan multiple radial positions of multiple shots, then each shot (=experiment) 
    is defined by its unique VMEC, while we have multiple simulations for each shot, 
    which correspond to the different radial position. '''
    
    # First create the simulations 
    simulations = create_simulations(folders, input_files, creationDetails.ignoreResolution, creationDetails.variables, creationDetails.folderIsExperiment, progress_start=creationDetails.progress_start, progress_stop=creationDetails.progress_stop1)
 
    # Create the experiments based on how many varied values there are
    if not creationDetails.resolutionScan: 
        experiments = create_experimentsBasedOnVariedValuesOrKnobKey(simulations, creationDetails) 
        experiments = create_listOfVariedValues(experiments, creationDetails.knob1, creationDetails.key1, creationDetails.knob2, creationDetails.key2)
    
    # Create experiments for a resolution scan
    if creationDetails.resolutionScan: 
        experiments = create_experimentsForAResolutionScan(simulations, creationDetails)
    
    # Make sure we don't have the same label for each simulation
    for experiment in experiments:
        remove_identicalSimulationLabels(experiment)
                    
    # Sort the simulations and varied values and set the labels and markers
    for experiment in experiments: 
        sort_listByNumbers(experiment, source="variedValues")   
        experiment.set_labelsLinesMarkers()   
        
    # Return the experiments
    return experiments

#------------------------------  
def create_experimentsBasedOnVariedValuesOrKnobKey(simulations, creationDetails):
    """ Iterate through the simulations in sort them in different experiments. """ 

    # Create an empty list to hold the experiment objects
    experiments = []
    
    # Extract the needed creation details
    key1 = creationDetails.key1 if creationDetails.key1!="rho" else "torflux"
    key2 = creationDetails.key2 if creationDetails.key2!="rho" else "torflux"
    key3 = creationDetails.key3 if creationDetails.key3!="rho" else "torflux"
    knob1 = creationDetails.knob1 
    knob2 = creationDetails.knob2
    knob3 = creationDetails.knob3 
    folderIsExperiment = creationDetails.folderIsExperiment
    number_variedVariables = creationDetails.variables 
    
    # Iterate over the simulations objects 
    for i, simulation in enumerate(simulations):  
        
        # Progress
        if creationDetails.progress_start!=None: print_progressbar(creationDetails.progress_stop1+(creationDetails.progress_stop-creationDetails.progress_stop1)*i/len(simulations), 100, prefix = '   Progress:', suffix = 'Creating research object: experiments.', length=50) 
        
        # Assume we have a new experiment, if not add the simulation to an existing experiment
        create_newExperiment = True
        
        # Go through all the current experiments and see if the simulation belongs to one of these 
        for experiment in experiments:   
            
            # If the simulations are in a different folder and <folderIsExperiment>=True, its a new experiment 
            if not (folderIsExperiment==True and simulation.path.input.parent != experiment.simulations[0].path.input.parent):
                 
                # Look at the difference in the input parameters of the experiment and the simulation
                dict_difference, numberOfDifferences = get_differenceInDictionaries(experiment.inputParameters, simulation.inputParameters)
                
                # Sometimes we rerun the exact same simulation with different modes: we have the same experiment
                # However if we want this simulation to be a new simulation, force this with number_variedValues=-1
                # If number_variedVariables==-2 then we want to sort by "1 folder = 1 experiment"
                # If number_variedVariables==-1 then we want to sort by "1 folder = 1 simulation" 
                if numberOfDifferences == 0: 
                    if number_variedVariables==-2:
                        create_newExperiment = True
                    if number_variedVariables!=-2:
                        create_newExperiment = False
                        experiment.simulations.append(simulation) 

                # The simulation belongs to the current experiment if inputParameters[knob1][key1] is the same.  
                if create_newExperiment==True and numberOfDifferences!=0: 
                    check1 = (knob1 in experiment.inputParameters) and (knob1 in simulation.inputParameters) 
                    check2 = (knob2 in experiment.inputParameters) and (knob2 in simulation.inputParameters) 
                    check3 = (knob3 in experiment.inputParameters) and (knob3 in simulation.inputParameters) 
                    value_exp1 = experiment.inputParameters[knob1][key1] if check1 else 0
                    value_sim1 = simulation.inputParameters[knob1][key1] if check1 else 1 if (knob1=="-") else 0
                    value_exp2 = experiment.inputParameters[knob2][key2] if check2 else 0
                    value_sim2 = simulation.inputParameters[knob2][key2] if check2 else 0 if (knob2=="-") else 1    
                    value_exp3 = experiment.inputParameters[knob3][key3] if check3 else 0
                    value_sim3 = simulation.inputParameters[knob3][key3] if check3 else 0 if (knob3=="-") else 1    
                    if (value_exp1==value_sim1) and (value_sim2==value_exp2) and (value_sim3==value_exp3):    
                        create_newExperiment = False
                        experiment.simulations.append(simulation)
                        if (value_exp1==value_sim1):  
                            varied_variable = knob1 + ": " + key1
                            if varied_variable not in experiment.variedVariables:
                                experiment.variedVariables.append(varied_variable)
                        if (value_sim2==value_exp2):  
                            varied_variable2 = knob2 + ": " + key2
                            if varied_variable2 not in experiment.variedVariables:
                                experiment.variedVariables.append(varied_variable2) 
                        if (value_sim3==value_exp3):  
                            varied_variable3 = knob3 + ": " + key3
                            if varied_variable3 not in experiment.variedVariables:
                                experiment.variedVariables.append(varied_variable3) 
                        for knob in list(dict_difference.keys()):
                            for key in dict_difference[knob]:
                                varied_variable = knob + ": " + key
                                if varied_variable not in experiment.variedVariables:
                                    experiment.variedVariables.append(varied_variable)
                            
                # The simulation is part of the experiment if up to <number_variedVariables> input variables are different  
                if create_newExperiment==True and numberOfDifferences!=0:
                    newExperimentBasedOnSelectedKeys = False   
                    if knob1 in dict_difference:
                        if key1 in dict_difference[knob1]:
                            newExperimentBasedOnSelectedKeys = True
                    if knob2 in dict_difference:
                        if key2 in dict_difference[knob2]: 
                            newExperimentBasedOnSelectedKeys = True            
                    if numberOfDifferences<=number_variedVariables and not newExperimentBasedOnSelectedKeys: 
                        create_newExperiment = False 
                        experiment.simulations.append(simulation)
                        for knob in list(dict_difference.keys()):
                            for key in dict_difference[knob]:
                                varied_variable = knob + ": " + key
                                if varied_variable not in experiment.variedVariables:
                                    experiment.variedVariables.append(varied_variable)
        
        # If more than <number_variedVariables> variables are different we have a new experiment
        # If inputParameters[knob1][key1] is different we have a new experiment
        # If the experiments list is empty we have a new experiment
        # In this case create an <Experiment> object and add it to the list 
        if create_newExperiment==True:    
            experiments.append(Experiment(simulation))

    return experiments

#--------------------------------------------- 
def create_listOfVariedValues(experiments, knob1, key1, knob2, key2, variedVariables=[]):

    # If experiment.variedVariables==[], try to take the value of the other experiments 
    for experiment in experiments: 
        variedVariables += experiment.variedVariables 
    variedVariables=sorted([v for v in set(variedVariables) if len(v)>0], key=variedVariables.count, reverse=True)
    for experiment in experiments: 
        if experiment.variedVariables==[] and variedVariables!=[]:
            experiment.variedVariables = [variedVariables[0]] 
    
    # For the variables that are varied, get their exact values for the simulation labels in the graph
    # For example if experiment.variedVariables = ["vmec_parameter: torflux", "knobs: delt"]
    # then experiment.variedValues = ["rho=0.5; delt=0.1", "rho=0.5; delt=0.01"; "rho=0.6; delt=0.1"]
    for experiment in experiments:
        for simulation in experiment.simulations:
            
            # Initiate the varied values
            varied_values = []  
            
            # Iterate through the varied variables and if it is not the chosen 
            #(knob, key) then add this (key, value) as the label for the simulation
            for variable in experiment.variedVariables:
                knob = variable.split(":")[0]
                variable = variable.split(": ")[-1] 
                if knob!="-":
                    if (knob!=knob1 or variable!=key1) and (knob!=knob2 or variable!=key2):  
                        value = simulation.inputParameters[knob][variable]
                        variable, value = change_stellaParametersForLabels(simulation, knob, variable, value)
                        varied_values += [variable + "$\,=\,$" + str(value)]
             
            # If a/Lne = a/Lni replace it by a/Ln
            if any("$a/L_{ni}$" in v for v in varied_values) and any("$a/L_{ne}$" in v for v in varied_values):
                _elements = [e for e in varied_values if e!='']
                _ionDensityGrad = float([e for e in _elements if "$a/L_{ni}$" in e][0].split("$\,=\,$")[-1])
                _others = [e for e in _elements if ("$a/L_{ni}$" not in e) and ("$a/L_{ne}$" not in e)] 
                if len(_others) != 0: varied_values = _others + ["$a/L_{n}$" + "$\,=\,$" + str(_ionDensityGrad)]
                if len(_others) == 0: varied_values = ["$a/L_{n}$" + "$\,=\,$" + str(_ionDensityGrad)]

            # If kymax and ny are both in varied_values then remove ny
            if any("ky max" in v for v in varied_values) and any("ny" in v for v in varied_values):
                _elements = [e for e in varied_values if e!='']
                _kymax  = [e for e in _elements if "ky max" in e][0]
                _others = [e for e in _elements if ("ky max" not in e) and ("ny" not in e)] 
                if len(_others) != 0: varied_values = _others + [_kymax]
                if len(_others) == 0: varied_values = [_kymax]
                          
            # If nothing was different because we only have one experiment, use the radial position as the label 
            if varied_values==[]: 
                if simulation.linear: varied_values = ["$\\rho$" + "$\,=\,$" +  str(simulation.modes[0].input.rho)] 
                if simulation.nonlinear: varied_values = ["$\\rho$" + "$\,=\,$" +  str(simulation.input.rho)] 
             
            # Sort the varied values and then join them in a string 
            varied_values.sort()
            varied_values = "; ".join(varied_values)
            
            # Add the string of varied values to the dictionary  
            experiment.variedValues.append(varied_values) 

    return experiments

#---------------------------------------------
def create_experimentsForAResolutionScan(simulations, creationDetails):
    
    # Create an empty list to hold the experiment objects
    experiments = []
    
    # For each simulation, remember its resolution
    resolution = {
        "nzed"   : { "values" : [None]*len(simulations), "knob" : "zgrid_parameters", "key" : "nzed"},\
        "nmu"    : { "values" : [None]*len(simulations), "knob" : "vpamu_grids_parameters", "key" : "nmu"},\
        "nvgrid" : { "values" : [None]*len(simulations), "knob" : "vpamu_grids_parameters", "key" : "nvgrid"},\
        "y0"     : { "values" : [None]*len(simulations), "knob" : "kt_grids_box_parameters", "key" : "y0"},\
        "kx max" : { "values" : [None]*len(simulations), "knob" : "kt_grids_box_parameters", "key" : "kx max"},\
        "ky max" : { "values" : [None]*len(simulations), "knob" : "kt_grids_box_parameters", "key" : "ky max"}}
     
    # Iterate over the simulations objects and collect the resolution parameters
    for i, simulation in enumerate(simulations):
        for key in resolution.keys():
            if resolution[key]["key"]!="kx max": 
                resolution[key]["values"][i] = simulation.input.inputParameters[resolution[key]["knob"]][resolution[key]["key"]]
            if resolution[key]["key"]=="kx max":  
                resolution[key]["values"][i] = np.round(simulation.input.inputParameters[resolution[key]["knob"]][resolution[key]["key"]],1)

    # Determine the low resolution choices 
    low_resolution  =  {
        "nzed"    :  np.min(resolution["nzed"]["values"]),\
        "nmu"     :  np.min(resolution["nmu"]["values"]),\
        "nvgrid"  :  np.min(resolution["nvgrid"]["values"]),\
        "y0"      :  np.min(resolution["y0"]["values"]),\
        "kx max"  :  np.min(resolution["kx max"]["values"]),\
        "ky max"  :  np.min(resolution["ky max"]["values"])}
    
    # Make experiments based in the resolution
    for i, simulation in enumerate(simulations):
        
        # Progress
        if creationDetails.progress_start!=None: print_progressbar(creationDetails.progress_stop1+(creationDetails.progress_stop-creationDetails.progress_stop1)*i/len(simulations), 100, prefix = '   Progress:', suffix = 'Creating research object: experiments.', length=50) 
        
        # Get the low resolution experiment
        if all(resolution[key]["values"][i] == low_resolution[key] for key in resolution.keys()):
            
            # Assume we need to make a new experiment
            make_newExperiment = True
                    
            # First try to add it to the old experiment
            for experiment in experiments: 
                if "Low resolution" in experiment.variedValues[0]:
                    experiment.simulations.append(simulation)
                    experiment.variedVariables.append("vmec_parameter: torflux")
                    experiment.variedValues.append("Low resolution"+" ("+str(len(experiment.variedValues))+")")
                    make_newExperiment = False
                    break

            # Otherwise, make a new experiment
            if make_newExperiment:                
                experiment = Experiment(simulation)
                experiment.variedVariables.append("vmec_parameter: torflux")
                experiment.variedValues.append("Low resolution")
                experiments.append(experiment)

        # Get the higher resolution experiments
        else:
            for key in resolution.keys():
                if resolution[key]["values"][i] != low_resolution[key]:
                    
                    # Assume we need to make a new experiment
                    make_newExperiment = True
                    
                    # First try to add it to the old experiments
                    for experiment in experiments: 
                        if key in experiment.variedValues[0]:
                            varied_variable = resolution[key]["knob"]+": "+resolution[key]["key"]
                            varied_value = resolution[key]["key"]+": "+str(low_resolution[key])+" $\\rightarrow$ "+str(resolution[key]["values"][i])
                            if varied_value==experiment.variedValues[0]:
                                experiment.simulations.append(simulation)
                                experiment.variedVariables.append(varied_variable)
                                experiment.variedValues.append(varied_value+" ("+str(len(experiment.variedValues))+")")
                                make_newExperiment = False
                                break
                            
                    # Otherwise, make a new experiment
                    if make_newExperiment:
                        varied_variable = resolution[key]["knob"]+": "+resolution[key]["key"]
                        varied_value = resolution[key]["key"]+": "+str(low_resolution[key])+" $\\rightarrow$ "+str(resolution[key]["values"][i])
                        experiment = Experiment(simulation)
                        experiment.variedVariables.append(varied_variable)
                        experiment.variedValues.append(varied_value)
                        experiments.append(experiment)

                    # Only add it once per keys
                    break
                    
    return experiments

#---------------------------------------------
def remove_identicalSimulationLabels(experiment):
    if len(set(experiment.variedValues))==1 and len(experiment.variedValues)!=1:
        simulation_ids = [simulation.id for simulation in experiment.simulations] 
        common_prefix = ''.join(c[0] for c in takewhile(lambda x:  all(x[0] == y for y in x), zip(*simulation_ids)))
        common_prefix = common_prefix if common_prefix.endswith("_") else '_'.join(common_prefix.split('_')[0:-1])+"_"
        simulation_ids2 = [simulation.id.split("__")[0] for simulation in experiment.simulations]
        common_prefix2 = ''.join(c[0] for c in takewhile(lambda x:  all(x[0] == y for y in x), zip(*simulation_ids2)))
        common_prefix2 = common_prefix if common_prefix2.endswith("_") else '_'.join(common_prefix2.split('_')[0:-1])+"_"
        for i in range(len(experiment.variedValues)): 
            if "__" in common_prefix: 
                simulation_ids = [s[::-1] for s in simulation_ids]
                common_suffix = ''.join(c[0] for c in takewhile(lambda x:  all(x[0] == y for y in x), zip(*simulation_ids)))
                simulation_ids = [s[::-1] for s in simulation_ids]; common_suffix = common_suffix[::-1]
                if common_suffix=="": common_suffix = "__"
                experiment.variedValues[i] = str(experiment.simulations[i].id).split(common_prefix)[-1].split(common_suffix)[0]
            if "__" not in common_prefix:  
                try: experiment.variedValues[i] = str(experiment.simulations[i].id).split("__")[0].split(common_prefix2)[1]
                except: experiment.variedValues[i] = str(experiment.simulations[i].id).split("__")[0] 
    return 

################################################################################
#                               CLASS EXPERIMENT                               #
################################################################################

class Experiment:
    ''' Make an experiment object which groups multiple simulation objects. '''
    
    def __init__(self, simulation):
        
        # Save the simulation objects and its input parameters
        # Note that the input parameters of the simulation can differ in (kx,ky)
        # And also in the resolution if ignore_resolution = True
        self.id = simulation.id
        self.simulations = [simulation]     
        self.inputParameters = simulation.inputParameters 
       
        # Save whichs [knob][parameter] is different in this experiment compared to the
        # other experiments: e.g. variedVariables = ["vmec_parameter: torflux", "knobs: delt"]
        # For each simulation add the exact value of this parameter that is different:
        # e.g. variedValues = ["rho=0.5; delt=0.1", "rho=0.5; delt=0.01"; "rho=0.6; delt=0.1"]
        self.variedVariables = []
        self.variedValues = []
        
        # Print information on the GUI (set from the GUI itself)
        self.Progress = None 
        if True: return 
 
#===============================================================================
#                GET THE TOTAL AMOUNT OF MODES IN THE EXPERIMENT               #
#=============================================================================== 

    def get_modes(self, kx=[], ky=[]): 
        for simulation in self.simulations:          
            kx += list(simulation.vec.kx)
            ky += list(simulation.vec.ky)
        kx = sorted(list(set(kx)))
        ky = sorted(list(set(ky)))
        self.vec = type('Vector', (object,), {'kx' : kx, 'ky' : ky}) 
        self.dim = type('Vector', (object,), {'kx' : len(kx), 'ky' : len(ky)}) 
        return
        
    @calculate_attributeWhenReadFirstTime 
    def dim(self):      self.get_modes();    return self.dim
    @calculate_attributeWhenReadFirstTime 
    def vec(self):      self.get_modes();    return self.vec
    

#===============================================================================
#                   LINES/MARKERS STYLE FOR EACH EXPERIMENT                    #
#=============================================================================== 

    def set_labelsLinesMarkers(self, line_label=None, marker_label=None, line_style=None, line_color=None, marker_style=None, marker_color=None):
        
        # If the styles are given, save them 
        if line_label is not None:
            self.line_label   = line_label[self.id]
            self.marker_label = marker_label[self.id]
            self.line_style   = line_style[self.id]
            self.line_color   = line_color[self.id]
            self.marker_style = marker_style[self.id]
            self.marker_color = marker_color[self.id] 
        
        # Otherise load some basic styles
        else: 
            line_label, marker_label, line_style, line_color, marker_style, marker_color = load_labelsLinesMarkers([self.id], self.variedValues)
            self.line_label = line_label[self.id] 
            self.marker_label = marker_label[self.id]
            self.line_style = line_style[self.id]
            self.line_color = line_color[self.id]
            self.marker_style = marker_style[self.id]
            self.marker_color = marker_color[self.id]

        # Change the lines/markers of the underlying simulation objects based on self.variedValues
        for simulation in self.simulations:
            variedValue = self.variedValues[self.simulations.index(simulation)] 
            simulation.set_labelsLinesMarkers(line_label[variedValue], marker_label[variedValue], line_style[variedValue], \
                                              line_color[variedValue], marker_style[variedValue], marker_color[variedValue]) 
        
 



