
import copy
import numpy as np
import os, configparser
from stellapy.utils.decorators.verbose import noverbose   
from stellapy.data.input.read_inputFile import read_inputFile
from stellapy.data.geometry.load_geometryObject import save_geometryObject
from stellapy.data.geometry.load_geometryObject import load_geometryObjectDirectly    
from stellapy.utils.files.get_filesInFolder import get_filesInFolder
from stellapy.utils.files import remove_simulationsWithoutOutncFile  
 
#===============================================================================
#                          WRITE THE GEOMETRY FILE
#===============================================================================
 
@noverbose
def write_h5FileForGeometry(folder, path=None, mode=None):  
    """ We give the path to the geometry file to avoid a recursion error, 
    when we ask for path.geometry it initially doesn't exist, and to find it
    the <path> module will look for the file and if it doesn't exist, it will
    call <write_h5FileForGeometry> and it will give the correct path. """    

#     # For each <mode> save the geometry to an *rho70.geometry file  
#     if mode!=None and path!=None:  
#         if not os.path.isfile(path):
#             save_geometryObject(mode)  
#             return

    # Get the input_files
    input_files = get_filesInFolder(folder, end=".in")
    
    # Only look at the input files in each parent folder
    for folder in list(set([i.parent if "OLD" not in i.parent.name else i.parent.parent for i in input_files])):    
        
        # Initialize
        modes = []
        
        # Get the files in this folder
        inputfiles = [i for i in input_files if (i.parent==folder) or (i.parent.parent==folder and "OLD" in i.parent.name)]
        inputfiles = remove_simulationsWithoutOutncFile(inputfiles)   
        if inputfiles==[]: continue
         
        # Create a simulation for this folder
        from stellapy.simulations.Simulation import create_simulations 
        simulations = create_simulations(None, inputfiles, write_uniqueFiles=True)
        
        # Collect the modes 
        if simulations[0].linear:  
            for simulation in simulations: 
                modes += simulation.modes  
                    
        # For each <mode> save the geometry to an *.unique.geometryX file 
        # Do this first to update the <input.list.geometries.ini> file
        if simulations[0].linear:  
            for unique_mode in get_uniqueGeometries(modes): 
                if not os.path.isfile(unique_mode.path.geometry):
                    save_geometryObject(unique_mode) 
                    
#         # For each <mode> save the geometry to an *rho70.geometry file  
#         if simulations[0].linear:   
#             for mode in modes: 
#                 if os.path.isfile(mode.path.output_stella):
#                     if not os.path.isfile(mode.path.geometry):
#                         save_geometryObject(mode)   
            
        # For each <simulation> save the geometry to an *rho70.geometry file 
        if simulations[0].nonlinear:  
            for simulation in simulations:  
                if path==None: path = simulation.path.geometry
                simulation.path.geometry = path 
                simulation.input.inputParameters = read_inputFile(simulation.input_file) 
                if not os.path.isfile(simulation.path.geometry):  
                    save_geometryObject(simulation) 
                path = None

    # Finish script  
    return

#-------------------------------
def get_uniqueGeometries(modes): 
    
    # Save each mode in unique_geometries[i] = [mode1, ..., modeX] 
    unique_geometries = {}
      
    # Only keep the knobs that are able to influence the geometry
    for mode in modes:  
        mode.inputParameters = read_inputFile(mode.input_file)
        mode.inputParameters = {"geo_knobs" : mode.inputParameters["geo_knobs"],
            "vpamu_grids_parameters" : mode.inputParameters["vpamu_grids_parameters"],
            "vmec_parameters" : mode.inputParameters["vmec_parameters"],
            "zgrid_parameters" : mode.inputParameters["zgrid_parameters"],
            "millergeo_parameters" : mode.inputParameters["millergeo_parameters"]}
    
    # First read the existing file and get the modes for each reference geometry
    path_of_reference_geometries = modes[0].path.folder / (modes[0].path.name + ".list.geometries.ini")
    if os.path.isfile(path_of_reference_geometries):
        
        # Read the reference geometries
        list_of_reference_geometries = configparser.ConfigParser()
        list_of_reference_geometries.read(path_of_reference_geometries)
        
        # Get the identifier of each reference geometry: "Geometry X"
        ids_of_reference_geometries = sorted([ g for g in list_of_reference_geometries.keys() if "Geometry " in g])  
        modes_to_be_sorted = [m for m in modes]   
        for reference_geometry_id in ids_of_reference_geometries:   
            reference_geometry_integer =  int(reference_geometry_id.split("Geometry ")[-1])
            for input_file in list_of_reference_geometries[reference_geometry_id].values():   
                for mode in modes_to_be_sorted: 
                    if str(mode.input_file)==input_file:
                        if reference_geometry_integer in unique_geometries.keys():
                            unique_geometries[reference_geometry_integer].append(mode)
                        if reference_geometry_integer not in unique_geometries.keys():
                            unique_geometries[reference_geometry_integer] = [mode] 
                        modes_to_be_sorted.remove(mode)
                        break
        
        # Check whether the loaded geometries were correct
        for reference_geometry_id in ids_of_reference_geometries: 
            if reference_geometry_id in unique_geometries.keys():
                reference_geometry_integer =  int(reference_geometry_id.split("Geometry ")[-1])
                reference_geometry_path = modes[0].path.folder / (modes[0].path.name + ".unique.geometry" + str(reference_geometry_integer))
                dummy_mode = copy.deepcopy(unique_geometries[reference_geometry_integer][0])
                dummy_mode.path.geometry = reference_geometry_path 
                load_geometryObjectDirectly(dummy_mode) 
                remove_modes = []  
            
                # In the newer version the geometry data knobs the relevant input geometry knobs
                if hasattr(dummy_mode.geometry, 'geo_knobs'):
                    for mode in unique_geometries[reference_geometry_integer]:
                        for knob in ["geo_knobs", "vmec_parameters", "zgrid_parameters", "millergeo_parameters", "vpamu_grids_parameters"]: 
                            for key in mode.inputParameters[knob].keys():
                                if mode.inputParameters[knob][key]!=getattr(getattr(mode.geometry, knob), key):
                                    remove_modes.append(mode); break
                            else: continue
                            break
                    
                # In the older version, check (nmu, nvgrid, nzed) and hope for the best
                if not hasattr(dummy_mode.geometry, 'geo_knobs'):
                    nzgrid, nmu = np.shape(dummy_mode.geometry.mu_weights)
                    nvgrid = np.shape(dummy_mode.geometry.vpa_weights)[0] 
                    for mode in unique_geometries[reference_geometry_integer]: 
                        nperiod = mode.inputParameters["zgrid_parameters"]["nperiod"]
                        nmu_mode = mode.inputParameters["vpamu_grids_parameters"]["nmu"]
                        nzed_mode = mode.inputParameters["zgrid_parameters"]["nzed"]
                        nzgrid_mode = 2*(int(nzed_mode/2) + int((nperiod-1)*nzed_mode)) + 1   
                        nvgrid_mode = 2*mode.inputParameters["vpamu_grids_parameters"]["nvgrid"] 
                        if nmu!=nmu_mode: remove_modes.append(mode)
                        elif nzgrid!=nzgrid_mode: remove_modes.append(mode)
                        elif nvgrid!=nvgrid_mode: remove_modes.append(mode)
                
                # Remove the modes that don't match their unique geometry
                for remove_mode in remove_modes: unique_geometries[reference_geometry_integer].remove(remove_mode)
      
    # Find the modes with an unique geometry
    for mode in modes:
        for ireference in unique_geometries.keys():
            if mode.inputParameters==unique_geometries[ireference][0].inputParameters: 
                if mode not in unique_geometries[ireference]:
                    unique_geometries[ireference].append(mode)
                break
        else:
            unique_geometries[len(list(unique_geometries.keys()))+1] = [mode] 
    
    # Save the unique geometries to a configuration file  
    unique_modes = save_uniqueGeometries(modes, unique_geometries)
    
    # Return a reference input file for each geometry 
    return unique_modes


#-------------------------------
def save_uniqueGeometries(modes, unique_geometries, difference={}, unique_modes=[]): 
 
    # Define the configuration file  
    ini_path = modes[0].path.folder / (modes[0].path.name +".list.geometries.ini") 
    if os.path.isfile(ini_path): os.system("rm "+str(ini_path))
    ini_data = configparser.ConfigParser()
    
    # Give some general information: the amount of unique geometries and the parent folder 
    ini_data.add_section('Unique Geometries')
    ini_data['Unique Geometries']['folder'] = modes[0].path.folder.name
    ini_data['Unique Geometries']['amount'] = count = str(len(list(set(list(unique_geometries.keys())))))
    print("  We found "+count+" unique geometries inside "+modes[0].path.folder.name+" containing "+str(len(modes))+" modes.") 

    # Find the knobs and values that are varied between the geometries
    keys = list(unique_geometries.keys())
    knobs = list(modes[0].inputParameters.keys())
    for i, j in ((i,j) for i in keys for j in keys): 
        iparameters = unique_geometries[i][0].inputParameters
        jparameters = unique_geometries[j][0].inputParameters
        for knob in knobs:
            if iparameters[knob]!=jparameters[knob]:
                if knob not in difference.keys():
                    difference[knob] = {}
                for parameter in iparameters[knob].keys():
                    if iparameters[knob][parameter]!=jparameters[knob][parameter]:
                        if parameter in difference[knob].keys():
                            difference[knob][parameter].append(iparameters[knob][parameter])
                            difference[knob][parameter].append(jparameters[knob][parameter])
                        if parameter not in difference[knob].keys():
                            difference[knob][parameter] = [iparameters[knob][parameter]]
                            difference[knob][parameter].append(jparameters[knob][parameter])
    
    # Save the information of the varied parameters to the configuration file
    ini_data.add_section('Varied parameters')    
    if len(list(difference.keys()))==0:
        ini_data['Varied parameters']['None'] = str([])
    for knob in difference.keys():
        for parameter in difference[knob].keys():
            ini_data['Varied parameters'][parameter] = str(sorted(list(set(difference[knob][parameter]))))
    
    # Save the modes that belong to each reference geometry
    ireferences = sorted(list(unique_geometries.keys()))
    for ireference in ireferences:
        
        # Make a section for each unique geometry
        ini_data.add_section("Geometry "+str(ireference))  
        
        # Save the unique mode and create its file path
        unique_modes_ireference = unique_geometries[ireference]
        unique_modes_ireference_with_output = [i for i in unique_modes_ireference if os.path.isfile(i.input_file.with_suffix('.out.nc'))]
        if len(unique_modes_ireference_with_output)!=0: reference_mode = unique_modes_ireference_with_output[0]
        if len(unique_modes_ireference_with_output)==0: reference_mode = unique_modes_ireference[0]
        unique_modes.append(reference_mode)
        reference_mode.path.geometry = reference_mode.path.folder / (reference_mode.path.name + ".unique.geometry"+str(ireference))
         
        # Sort the modes
        modes = unique_geometries[ireference]
        numbers = [mode.ky for mode in modes] 
        sorted_indexes = list(np.array(numbers).argsort()) 
        modes = [modes[i] for i in sorted_indexes] 
        modes_temp = [mode for mode in modes if "OLD" not in str(mode.input_file)] 
        modes_temp += [mode for mode in modes if "OLD/" in str(mode.input_file)] 
        for i in range(10): modes_temp += [mode for mode in modes if "OLD"+str(i) in str(mode.input_file)] 
        modes = modes_temp
        
        # Add the mode that belong to this geometry 
        for ifile, mode in enumerate(modes):
            ini_data["Geometry "+str(ireference)][str(ifile+1)] = str(mode.input_file)
            mode.path.geometry = reference_mode.path.geometry
            
    # Save the configuration file
    with open(ini_path, 'w') as ini_file:
        ini_data.write(ini_file) 

    # Return the unique input files  
    return unique_modes
 




