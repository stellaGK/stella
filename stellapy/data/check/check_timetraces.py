#!/usr/bin/python3
import sys, os 
import pathlib
import numpy as np
import netCDF4 as nc4  
 
# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)
from stellapy.utils.commandprompt.print_progressbar import print_progressbar
from stellapy.data.output.read_outFile import read_outFile
from stellapy.utils.files import get_filesInFolder 
from stellapy.utils.commandprompt.bash import Bash
 
#===============================================================================
#                               CHECK TIME TRACES
#===============================================================================
 
def check_timetraces(folder):
     
    # Count the simulation times and the files per folder
    time_step = {}
    last_times = {}
    output_files = {} 
    folder_size = {}
         
    # Check whether we have simulations
    if get_filesInFolder(folder, end=".in"): 
         
        # Go through the input files
        input_files = get_filesInFolder(folder, end=".in"); count = 1; length = len(input_files); print();  
        print_progressbar(0, length, prefix = '   Progress:', suffix = 'Complete', length = 50)
        for input_file in input_files:
             
            # Define the possible output files
            netcdf_file   = input_file.with_suffix(".out.nc")
            output_file   = input_file.with_suffix(".out")
            stellapy_file = input_file.with_suffix(".dt1.phi_vs_t") 
            
            # Get the memory usage of the folder
            directory = input_file.parent; directory2 = input_file.parent/"restart" 
            folder_size[input_file] = sum(os.path.getsize(directory/f) for f in os.listdir(directory) if os.path.isfile(directory/f))
            if os.path.isdir(directory2):
                folder_size[input_file] += sum(os.path.getsize(directory2/f) for f in os.listdir(directory2) if os.path.isfile(directory2/f))
             
            # Read tend and dt from the netcdf file
            if os.path.isfile(netcdf_file):
                netcdf_file = nc4.Dataset(netcdf_file)   
                tend = netcdf_file.variables["t"][-1]
                tprevious = netcdf_file.variables["t"][-2]
                time_step[input_file] = tend - tprevious
                last_times[input_file] = tend
                output_files[input_file] = "netcdf"
                
            # Read tend and dt from the output file
            elif os.path.isfile(output_file):  
                time_data = read_outFile(input_file)["time"]
                tend = time_data[-1]
                tprevious = time_data[-2]
                time_step[input_file] = tend - tprevious
                last_times[input_file] = tend
                output_files[input_file] = "out"
                
            # Read tend and dt from the stellapy output file
            elif os.path.isfile(stellapy_file):  
                data = np.loadtxt(stellapy_file,skiprows=1,dtype='float').reshape(-1, 6)
                time_data = data[:,0]   
                tend = time_data[-1]
                tprevious = time_data[-2]
                time_step[input_file] = tend - tprevious
                last_times[input_file] = tend
                output_files[input_file] = "phi_vs_t"
                
            # Show the progress
            print_progressbar(count, length, prefix = '   Progress:', suffix = 'Complete', length = 50); count += 1                         
             
    # Print the data
    print() 
    print("".center(80,"="))
    print("Time Traces".center(80))
    print("".center(80,"="))
    print()
    print("{0:^10}".format("DELTA") + \
          "{0:^10}".format("TEND") + \
          "{0:^15}".format("DATA FILE") + \
          "{0:^15}".format("SIZE (Gb)") + \
          "{0:^5}".format("") + \
          "{0:^10}".format("FILE"))
    print("{0:^10}".format("------") + \
          "{0:^10}".format("------") + \
          "{0:^15}".format("----------") + \
          "{0:^15}".format("---------") + \
          "{0:^5}".format("") + \
          "{0:^10}".format("----------"))
     
    # Collect the input files 
    input_files = list(last_times.keys())
    commonprefix = os.path.commonprefix(input_files)
    input_files = sorted(input_files)

    # Iterate over the input files
    for input_file in input_files: 
         
        # Get the data 
        file = str(input_file).split(str(pathlib.Path(commonprefix).parent)+"/")[-1] if commonprefix!=str(input_file) else input_file.name  
        delt = "{:.2f}".format(np.round(time_step[input_file],2))
        t_end = int(last_times[input_file])
        data_file = output_files[input_file]
        memory_usage = "{:.1f}".format(np.round(folder_size[input_file]/10**9,1)) 
     
        # Make sure the file name isn't too long
        if len(file) > 50: file = file[0:50]
     
        # Print the simulation time
        print("{0:^10}".format(delt) + \
              "{0:^10}".format(t_end) + \
              "{0:^15}".format(data_file) + \
              "{0:^15}".format(memory_usage) + \
              "{0:^5}".format("") + \
              "{0:50}".format(file))
        
    print()
    return
           
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================

if __name__ == "__main__":
    bash = Bash(check_timetraces, __doc__)  
    check_timetraces(**bash.get_arguments())   
