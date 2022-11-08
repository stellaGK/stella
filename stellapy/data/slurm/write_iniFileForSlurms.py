
import time
import os, sys 
import numpy as np
import configparser
from datetime import datetime, timedelta
from stellapy.utils.decorators.verbose import noverbose   
from stellapy.utils.decorators.exit_program import exit_program 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder
 
#===============================================================================
#                 WRITE A FILE WITH AN OVERVIEW OF THE SLURMS
#===============================================================================
# [input files]
# 5830814 = input_ky10.0
# 
# [start simulation]
# 5830814 = 2 oct 2018
# 
# [end simulation]
# 5830814 = 2 oct 2018

@noverbose
def write_iniFileForSlurms(folder):   
    
    # Only excute this command on the supercomputer 
    if ("marconi" not in folder.parts[0]) and ("marconi" not in folder.parts[1]):
        print("Can only determine simulation times on the supercomputer.")
        return
     
    # Check whether we have slurm files
    if get_filesInFolder(folder, start="slurm-"):
     
        # Define the configuration file
        input_files = get_filesInFolder(folder, end=".in") 
        input_files = [f for f in input_files if "OLD" not in input_files] + [f for f in input_files if "OLD" in input_files] 
        input_file = input_files[0].name if "_ky" not in input_files[0].name else input_files[0].name.split("_ky")[0]
        ini_path = folder / (input_file +".list.slurms.ini") 
            
        # Get the slurm files
        slurm_inputfiles = sorted(get_filesInFolder(folder, start="slurm_"))   
        slurm_statusfiles = sorted(get_filesInFolder(folder, start="slurm-", end=".out")) 
        slurm_outputfiles = sorted(get_filesInFolder(folder, end=".out")) 
        slurm_outputfiles = [ f for f in slurm_outputfiles if f.name[0:7].isdigit()] 
            
        # If the file exists, check whether it was made after the simulation had finished
        if os.path.isfile(ini_path): 
            times = [os.path.getctime(path) for path in slurm_inputfiles]
            times += [os.path.getctime(path) for path in slurm_statusfiles]
            times += [os.path.getctime(path) for path in slurm_outputfiles]
            times = [datetime.fromtimestamp(time) for time in times]
            times = [time - times[0] for time in times] 
            if np.all([np.abs(times) < timedelta(minutes=5)]):
                print("A touchall was performed and we can not calculate the simulation times.") 
            else:
                os.system("rm "+str(ini_path))  
        
        # If the file doesn't exit, create it
        if not os.path.isfile(ini_path): 
            
            # Store the data in a nested dictionary
            slurmdata = {}; inputfilesdata = {}
            
            # Get all the slurm ids
            ids = [f.name[0:7] for f in slurm_outputfiles]
            for identifier in ids: slurmdata[identifier] = {}
            
            # Link the slurm ids to input files
            for slurm in slurm_statusfiles:
                input_data = open(slurm, 'r')
                input_text = input_data.read() 
                input_file = slurm.parent / input_text.split("INPUT FILE:")[-1].split("\n")[0].replace(" ", "") 
                identifier = slurm.name.split("-")[-1].split(".out")[0]
                starttime = os.path.getctime(slurm) 
                slurmdata[identifier]['input_file'] = input_file
                slurmdata[identifier]['starttime'] = starttime
                if str(input_file) in inputfilesdata: inputfilesdata[str(input_file)] += [identifier]
                if str(input_file) not in inputfilesdata: inputfilesdata[str(input_file)] = [identifier]
                    
            # If there are multiple slurms for the same input file, sort them from oldest to newest 
            for i in inputfilesdata.keys():
                if len(inputfilesdata[i])>1:
                    exit_program('Implement this.', write_iniFileForSlurms, sys._getframe().f_lineno)
            
            # Save the requested simulation time and the number of nodes and the setup time
            for slurm in slurm_inputfiles:
                input_data = open(slurm, 'r')
                input_text = input_data.read() 
                input_file = input_text.split("INPUT FILE:")[-1].split("\n")[0].replace(" ", "").replace('"','')
                wall_time  = input_text.split("#SBATCH --time")[-1].split("\n")[0].replace(" ", "") 
                account = input_text.split("#SBATCH -A")[-1].split("\n")[0].replace(" ", "")
                nodes = int(input_text.split("#SBATCH -N ")[-1].split("\n")[0].replace(" ", ""))
                cores = str(nodes*48); nodes = str(nodes)
                setup = os.path.getctime(slurm) 
                identifier = inputfilesdata[str(folder / input_file)][0]
                slurmdata[identifier]['wall_time'] = wall_time
                slurmdata[identifier]['account'] = account
                slurmdata[identifier]['nodes'] = nodes
                slurmdata[identifier]['cores'] = cores
                slurmdata[identifier]['setup'] = setup
                
            # Try the find the total simulation time from the output file 
            for slurm in slurm_outputfiles:
                input_data = open(slurm, 'r')
                input_text = input_data.read() 
                if "total:" in input_text: simulation_time = float(input_text.split("total:")[-1].split("min")[0])
                if "total:" not in input_text: simulation_time = "/"
                endtime = os.path.getctime(slurm)
                identifier = slurm.name[0:7]
                slurmdata[identifier]['totaltime'] =simulation_time
                slurmdata[identifier]['endtime'] = endtime
            
            # Write the data to a configuration file
            ini_data = configparser.ConfigParser() 
            ini_data.add_section('General')
            ini_data['General']['Account'] = account
            ini_data['General']['Nodes'] = nodes
                
            ini_data.add_section('Input files')
            ini_data.add_section('Total times')
            ini_data.add_section('Wall times')
            ini_data.add_section('Setup times')
            ini_data.add_section('Start times')
            ini_data.add_section('Stop times')
            ini_data.add_section('Cores')
            
            for identifier in slurmdata.keys():   
                ini_data['Input files'][identifier] = str(slurmdata[identifier]['input_file'])
                ini_data['Total times'][identifier] = str(timedelta(minutes=slurmdata[identifier]['totaltime'])) if slurmdata[identifier]['totaltime']!="/" else "/"
                ini_data['Wall times'][identifier]  = slurmdata[identifier]['wall_time']
                ini_data['Setup times'][identifier] = time.ctime(slurmdata[identifier]['setup'])
                ini_data['Start times'][identifier] = time.ctime(slurmdata[identifier]['starttime'])
                ini_data['Stop times'][identifier]  = time.ctime(slurmdata[identifier]['endtime'])
                ini_data['Cores'][identifier] = slurmdata[identifier]['cores']
    
            with open(ini_path, 'w') as ini_file:
                ini_data.write(ini_file) 
                print("Wrote the slurm data to:", ini_file.name)
            
    # Finish script  
    return
    
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    import pathlib
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")
    write_iniFileForSlurms(folder)

