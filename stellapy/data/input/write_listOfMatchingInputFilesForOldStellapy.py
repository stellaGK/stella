"""

#===============================================================================
#                    WRITE LIST OF MATCHING INPUT FILES
#===============================================================================
 
The old stellapy loaded each (kx,ky) simulation as a <mode> object. This has 
beenremoved, and the old files are rewritten with this script in order to 
work with the new version of stellapy (September 2022).

If something is broken, restart by doing
    >>  rm */*_dummy.in

Hanne Thienpondt
13/10/2022
 
"""

#!/usr/bin/python3
import h5py
import os, sys
import pathlib
import configparser
import numpy as np

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.input.write_listOfMatchingInputFiles import read_relevantInputParameters, get_namesForDummyInputFiles, save_dummyInputFiles
from stellapy.utils.files.keep_linearFluxTubeSimulations import keep_linearFluxTubeSimulations  
from stellapy.data.input.read_inputFile import read_modeFromInputFile
from stellapy.data.omega.read_omegaFile import read_omegaFromTxtFile
from stellapy.utils.files.get_filesInFolder import get_filesInFolder   
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                    WRITE LIST OF MATCHING INPUT FILES
#===============================================================================

def write_listOfMatchingInputFilesForOldStellapy(folder):
    
    # Read the <list> of input files (part of old stellapy)
    lists_of_input_files = get_filesInFolder(folder, end="list.inputs.ini") 
    if len(lists_of_input_files)==0: return 
    
    # Get the input_files inside <folder>
    # Only keep linear flux tube simulations  
    input_files = get_filesInFolder(folder, end=".in") 
    input_files = keep_linearFluxTubeSimulations(input_files) 
    
    # Collect the parent directories since we will write a list for each directory
    # For linear restarts, the old files are put in an "OLD" folder, include these in the list of OLD.parent
    parent_directories = list(set([i.parent if "OLD" not in i.parent.name else i.parent.parent for i in input_files]))
    
    # Only look at the input files in each parent folder
    for directory in parent_directories:    
        
        # Read the <list> of input files
        list_of_input_files = [i for i in lists_of_input_files if i.parent==directory] 
        if len(list_of_input_files)>1: exit_program("Not implemented yet", write_listOfMatchingInputFilesForOldStellapy, sys._getframe().f_lineno)
        if len(list_of_input_files)==0: continue 
        
        # Open the configuration file 
        ini_path = list_of_input_files[0] 
        ini_data = configparser.ConfigParser() 
        ini_data.read(ini_path) 
        
        # Number of unique input files = number of dummy input files
        number_of_dummy_files = int(ini_data["Unique Input Files"]["amount"])

        # Initialize a dictionary simulation{ID}{list,input_parameters}
        simulations = {}; input_files_per_dummy_file = {}; remove_keys = []
        
        # Iterate over the dummy input files
        for i in range(number_of_dummy_files):
            
            # Read the input files
            input_files_per_dummy_file[i] = []
            for key in ini_data["Input "+str(i+1)].keys():
                input_file_in_ini_file = pathlib.Path(ini_data["Input "+str(i+1)][key])  
                input_files_per_dummy_file[i].append(input_file_in_ini_file)
                
        # The absolute path of the input files might be wrong, the parent should be <directory>
        for dummy_input_file_key in input_files_per_dummy_file.keys():
            input_files_temp = input_files_per_dummy_file[dummy_input_file_key] 
            input_files_temp = [i.parts for i in input_files_temp] 
            commonprefix = pathlib.Path(*os.path.commonprefix(input_files_temp))
            input_files_temp = input_files_per_dummy_file[dummy_input_file_key] 
            input_files_temp = [pathlib.Path(str(i).replace(str(commonprefix), str(directory))) for i in input_files_temp] 
            input_files_temp = [i for i in input_files_temp if os.path.isfile(i)]
            input_files_per_dummy_file[dummy_input_file_key] = input_files_temp
            if len(input_files_temp)==0: remove_keys.append(dummy_input_file_key)
        for key in remove_keys:
            del input_files_per_dummy_file[key]

        # Create a dummy input file for each unique input file 
        for dummy_input_file_key in input_files_per_dummy_file.keys():

            # Read the input files that are identicial besides (kx,ky)
            input_files_in_dummy_file = input_files_per_dummy_file[dummy_input_file_key] 
            
            # Read the input parameters 
            input_parameters = read_relevantInputParameters(input_files_in_dummy_file[0])
                    
            # Create the simulation
            simulations[str(input_files_in_dummy_file[0])] = {}
            simulations[str(input_files_in_dummy_file[0])]["input_parameters"] = input_parameters
            simulations[str(input_files_in_dummy_file[0])]["input_files"] = input_files_in_dummy_file 
            simulations[str(input_files_in_dummy_file[0])]["unique ID"] = dummy_input_file_key+1
            
            # Find an input file which as a dimensions file
            for input_file in input_files_in_dummy_file: 
                if os.path.isfile(input_file.with_suffix(".dimensions")): 
                    simulations[str(input_files_in_dummy_file[0])]["dimensions_file"] = str(input_file).replace(".in", ".dimensions")
                    break
            else:
                exit_program("No dimension file could be found.", write_listOfMatchingInputFilesForOldStellapy, sys._getframe().f_lineno)
            
        # Get the names of the dummy input files, saved to simulations[simulation_id]["input_file_name"]
        get_namesForDummyInputFiles(simulations)
        
        # Save each dummy input file to a txt file
        save_dummyInputFiles(simulations, folder, directory)
        
        # Copy the ini files to the dummy files
        for simulation_id in simulations.keys():  
            old_ini_file = str(directory)+"/input.unique.input"+str(simulations[simulation_id]["unique ID"])+".ini"
            new_ini_file = str(directory)+"/"+simulations[simulation_id]["input_file_name"].replace(".in", ".ini")
            os.system("cp "+old_ini_file+" "+new_ini_file) 
            
        # Copy the geometry files to the dummy files
        for simulation_id in simulations.keys():  
            old_geo_file = read_unqiueGeometryFile(directory, simulations[simulation_id]["input_files"][0], commonprefix) 
            new_geo_file = str(directory)+"/"+simulations[simulation_id]["input_file_name"].replace(".in", ".geo")
            os.system("cp "+old_geo_file+" "+new_geo_file) 
        
        # Copy the dimension files to the dummy files
        for simulation_id in simulations.keys():  
            old_dim_file = str(simulations[simulation_id]["dimensions_file"]) 
            new_dim_file = str(directory)+"/"+simulations[simulation_id]["input_file_name"].replace(".in", ".dimensions")
            os.system("cp "+old_dim_file+" "+new_dim_file) 
            
            # Update <vec_ky> and <vec_kx> in the dimensions file 
            vec_kx = []; vec_ky = []; dim_time = 0
            for input_file in simulations[simulation_id]["input_files"]: 
                if not os.path.isfile(input_file.with_suffix(".dimensions")): 
                    kx, ky = read_modeFromInputFile(input_file)
                    vec_kx += [kx]; vec_ky += [ky]  
                    dim_time = np.max([dim_time, np.nanmax(np.loadtxt(input_file.with_suffix(".dt1.omega_vs_t"), dtype='float').reshape(-1, 3)[:,0])])
                if os.path.isfile(input_file.with_suffix(".dimensions")): 
                    with h5py.File(input_file.with_suffix(".dimensions"), 'r') as h5_file:  
                        vec_kx += list(h5_file['vec_kx'][()])
                        vec_ky += list(h5_file['vec_ky'][()])
                        dim_time = np.max([dim_time, len(h5_file['vec_time'][()])]) 
            vec_kx = np.array(sorted(list(set(vec_kx)))); dim_kx = len(vec_kx)
            vec_ky = np.array(sorted(list(set(vec_ky)))); dim_ky = len(vec_ky) 
            
            # Update the dummy dimensions file 
            with h5py.File(new_dim_file, 'r+') as h5_file:
                del h5_file["vec_kx"]; h5_file.create_dataset("vec_kx", data=vec_kx)  
                del h5_file["vec_ky"]; h5_file.create_dataset("vec_ky", data=vec_ky)  
                del h5_file["dim_kx"]; h5_file.create_dataset("dim_kx", data=dim_kx)  
                del h5_file["dim_ky"]; h5_file.create_dataset("dim_ky", data=dim_ky)  
                del h5_file["dim_time"]; h5_file.create_dataset("dim_time", data=dim_time)  
             
    return

#--------------------------
def read_unqiueGeometryFile(directory, input_file, commonprefix):

    # Open the configuration file 
    lists_of_geometry_files = get_filesInFolder(directory, end="list.geometries.ini") 
    if len(lists_of_geometry_files)>1: exit_program("Not implemented yet", read_unqiueGeometryFile, sys._getframe().f_lineno)
    ini_path = lists_of_geometry_files[0]
    ini_data = configparser.ConfigParser() 
    ini_data.read(ini_path) 
    
    # Number of unique input files
    found_unique_geometry = False
    number_of_geometry_files = int(ini_data["Unique Geometries"]["amount"])
    for i in range(number_of_geometry_files):
        for key in ini_data["Geometry "+str(i+1)].keys():
            input_file_read = pathlib.Path(str(ini_data["Geometry "+str(i+1)][key]).replace(str(commonprefix), str(directory)))
            if input_file==input_file_read:
                found_unique_geometry = True
                break  
    if found_unique_geometry==False: exit_program("Couldn't find unique geometry file.", read_unqiueGeometryFile, sys._getframe().f_lineno)
            
    # Get the path of the unique geometry file
    old_geo_file = str(directory) + "/input.unique.geometry"+str(i+1)
    return old_geo_file
         
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface
    bash = Bash(write_listOfMatchingInputFilesForOldStellapy, __doc__)     
    
    # Get the bash arguments and execute the script 
    write_listOfMatchingInputFilesForOldStellapy(**bash.get_arguments())   
