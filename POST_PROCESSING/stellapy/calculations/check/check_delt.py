#!/usr/bin/python3
import sys, os, pathlib  
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)
from stellapy.utils import get_filesInFolder 
    
#======================
# Read the out files
#======================

# Only show a selected number of modes
ky_max = 2

# Count the simulation times and the files per folder
resolution = {"delt" : {}, "tend" : {}} 
last_times = {}
        
# Get the folder
directory = pathlib.Path(os.getcwd())

# Go throught the subfolders
#folders = [pathlib.Path(x[0]) for x in os.walk(directory)]

# Only go through the main folders
folders = sorted([directory / name for name in os.listdir(directory) if os.path.isdir(name)])
folders = [f for f in folders if f.name!="mat"]
folders_old = [f for f in folders if "OLD" in str(f)]
if folders_old != []: folders += [directory]
if folders == []: folders = [directory]

# Get the stella output files
for folder in folders:
    
    # Check whether we have simulations
    if get_filesInFolder(folder, end=".in"):
        
        # Get the files
        output_files  = get_filesInFolder(folder, end=".out") 
        output_files = [f for f in output_files if "input" not in str(f)] 
        output_files = [f for f in output_files if "test" not in str(f)] 
        output_files = [f for f in output_files if "slurm" not in str(f)] 
            
        # Go through the files
        for i in range(len(output_files)):
            
            # Get actual end time in the output file
            file   = output_files[i]
            data   = open(file, 'r')
            text   = data.read()   
            
            # Get the time step
            if "Decreasing code_dt to cfl_dt*cfl_cushion/delt_adjust" not in text:
                last_line = text.split("\n") 
                last_line = [l for l in last_line if l!=""][-1]
                last_line = last_line.split(" ")
                last_line = [l for l in last_line if l!="" and l!="\n"]  
                tend = float(last_line[1])
                delt = float(last_line[2])
            else:
                try:
                    delt = text.split("Decreasing code_dt to cfl_dt*cfl_cushion/delt_adjust")[-1].split("\n")[0]
                    test = text.split(str(delt)) 
                    tend = text.split(str(delt))[-2].split("\n")[1].split(" ")
                    tend = [t for t in tend if t!="" and t!="\n"][1] 
                    try: delt = float(delt)
                    except: delt = float(delt.replace("-", "e-"))
                    tend = float(tend)
                except:
                    try:
                        delt = text.split("Decreasing code_dt to cfl_dt*cfl_cushion/delt_adjust")[-2].split("\n")[0]
                        test = text.split(str(delt)) 
                        tend = text.split(str(delt))[-2].split(" ")
                        tend = [t for t in tend if t!="" and t!="\n"][1] 
                        try: delt = float(delt)
                        except: delt = float(delt.replace("-", "e-"))
                        tend = float(tend)
                    except:
                        continue
        
            # Save the simulation time per parent folder 
            if str(folder)+"/"+file.name in resolution["tend"].keys():
                resolution["delt"][str(folder)+"/"+file.name] += [delt]
                resolution["tend"][str(folder)+"/"+file.name] += [tend]
            if str(folder)+"/"+file.name not in resolution["tend"].keys():
                resolution["delt"][str(folder)+"/"+file.name] = [delt]
                resolution["tend"][str(folder)+"/"+file.name] = [tend]                
                        
        
# Print the data
print()
print("############################################################################")
print("                      SIMULATION TIME PER FILE                              ")
print("############################################################################")
print()
print("{0:^40}".format("FILE") + \
      "{0:^15}".format("DELTA T") + \
      "{0:^20}".format("END TIME"))
print("{0:^40}".format("----") + \
      "{0:^15}".format("------") + \
      "{0:^20}".format("--------"))

# Go through all the folders
keys = list(resolution["tend"].keys())
commonprefix = os.path.commonprefix(keys)
for key in list(resolution["tend"].keys()):
    
        # Get the data
        file = key.split(str(commonprefix))[-1];    
        name = "      "+pathlib.Path(key).parent.name
        folder = "/".join(key.split("/")[:-1]);       
        delt = "{:.2E}".format(resolution["delt"][key][0])
        tend = int(resolution["tend"][key][0])
    
        # Make sure the file name isn't too long
        if len(file) > 50: file = file[0:50]
    
        # Print the simulation time
        print("{0:40}".format(name) + \
              "{0:^15}".format(delt) + \
              "{0:^20}".format(tend))
print()

