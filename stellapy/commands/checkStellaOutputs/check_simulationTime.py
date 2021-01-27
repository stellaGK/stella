#!/usr/bin/python3


import sys, os, pathlib  
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.utils import get_filesInFolder 

#===========================================
# Convert minutes to hours:minutes:seconds
#===========================================

def get_hours(time_inMinutes):
    ''' Convert the minutes to hours:minutes:seconds. '''
    hours   = int(time_inMinutes/60)
    minutes = int(time_inMinutes % 60)
    seconds = int((time_inMinutes*60) % 60)
    return "%d:%02d:%02d" % (hours, minutes, seconds)
 
    
#======================
# Read the out files
#======================

# Count the simulation times and the files
simulations_times = {}
cores = {}
        
# Get the folder
directory = pathlib.Path(os.getcwd())

# Go throught the subfolders
folders = [pathlib.Path(x[0]) for x in os.walk(directory)]

# Get the stella output files
for folder in folders:
    
    if get_filesInFolder(folder, end="stella.out"):
        
        # Read the eu.slurm to get the number of cores used
        input_data   = open(folder / 'eu.slurm', 'r')
        input_text   = input_data.read() 
        cores[str(folder)] = float(input_text.split("-genv -n")[-1].split("$EXEFILE")[0]) 
        
        # Get the stella.out files
        files = sorted(get_filesInFolder(folder, end="stella.out"))
        
        # Go through the files
        for file in files:
            
            # Read the file
            data  = open(file, 'r')
            text  = data.read()   
            
            # Get the simulation time
            if "total:" in text:
                
                # Extract the simulation time
                simulation_time = float(text.split("total:")[-1].split("min")[0])
        
                # Save the simulation time per parent folder
                if str(file.parent) in simulations_times.keys():
                    simulations_times[str(file.parent)] += [simulation_time]
                if str(file.parent) not in simulations_times.keys():
                    simulations_times[str(file.parent)] = [simulation_time]
        
# Print the data
print()
print("############################################################################")
print("                     SIMULATION TIME PER FOLDER                             ")
print("############################################################################")
print()
print("{0:^30}".format("FOLDERS") + \
      "{0:^10}".format("FILES") + \
      "{0:^15}".format("TOTAL TIME") + \
      "{0:^15}".format("AVERAGE TIME") + \
      "{0:^25}".format("TOTAL TIME x CORES") + \
      "{0:^25}".format("AVERAGE TIME x CORES"))
print("{0:^30}".format("-------") + \
      "{0:^10}".format("-----") + \
      "{0:^15}".format("----------") + \
      "{0:^15}".format("------------")+ \
      "{0:^25}".format("------------------") + \
      "{0:^25}".format("--------------------"))

# Keep track of the total amount
total_folders = 0
total_files = 0
total_time = 0 
total_time_xCores = 0 

# Go through all the folders
for key in simulations_times.keys():
    
    # Get the data
    folder = key.split("/")[-1];                   total_folders+=1
    files  = len(simulations_times[key]);          total_files+=files
    total_minutes = sum(simulations_times[key]);   total_time+=total_minutes
    average_minutes = total_minutes/files;         total_time_xCores+=total_minutes*cores[key]

    # Make sure the folder name isn't too long
    if len(folder) > 30: folder = folder[0:30]

    # Print the simulation time
    print("{0:30}".format(folder) + \
          "{0:^10}".format(str(files)) + \
          "{0:^15}".format(get_hours(total_minutes)) + \
          "{0:^15}".format(get_hours(average_minutes)) + \
          "{0:^25}".format(get_hours(total_minutes*cores[key])) + \
          "{0:^25}".format(get_hours(average_minutes*cores[key])))
  


print()
print("############################################################################")
print("                        TOTAL SIMULATION TIME                              ")
print("############################################################################")
print()
print("{0:^30}".format("FOLDERS") + \
      "{0:^10}".format("FILES") + \
      "{0:^15}".format("TOTAL TIME") + \
      "{0:^15}".format("AVERAGE TIME") + \
      "{0:^25}".format("TOTAL TIME x CORES") + \
      "{0:^25}".format("AVERAGE TIME x CORES"))
print("{0:^30}".format("-------") + \
      "{0:^10}".format("-----") + \
      "{0:^15}".format("----------") + \
      "{0:^15}".format("------------") + \
      "{0:^25}".format("------------------") + \
      "{0:^25}".format("--------------------"))
              
# Print the total simulation time
print("{0:^30}".format(str(total_folders)) + \
      "{0:^10}".format(str(total_files)) + \
      "{0:^15}".format(get_hours(total_time)) + \
      "{0:^15}".format(get_hours(total_time/total_files)) + \
      "{0:^25}".format(get_hours(total_time_xCores)) + \
      "{0:^25}".format(get_hours(total_time_xCores/total_files)) + "\n")

    