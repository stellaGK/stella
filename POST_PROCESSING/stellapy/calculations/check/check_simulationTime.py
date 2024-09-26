#!/usr/bin/python3
 
import sys, os, pathlib  
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 

class initiate_nesteddictt(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def get_filesInFolder(folders=None, name_inputFile=None, start=None, end=None):
    ''' Returns all files inside <folder> that starts with <start> and end with <end>. 

    If <input_file> is given, look for corresponding output files.
    Return the file names if full_path = False, otherwise return the full path names.
    '''
    
    # Initiate the list
    files_inside = [] 
    
    # Make sure we have a list of folders
    if isinstance(folders, pathlib.PurePath):
        folders = [folders]
    
    # Go through the folders to find the required files
    for folder in folders:
        
        # Files inside this folder
        files_inside_folder = []
    
        # If <name_inputFile> is given, look for the corresponding output file
        if name_inputFile and end:
            if not end.startswith('.'): end = '.' + end
            file_inside = folder / name_inputFile.split('.in')[0] + end
            if os.path.isfile(file_inside):
                return file_inside
    
        # Read the files in <folder> that start with <start> and end with <end> 
        for file_name in os.listdir(folder):
            
            if not file_name.startswith('.'):
                
                # If the <file_name> starts with <start> and ends with <end>, add it to files_inside.
                if end and start:
                    if file_name.endswith(end) and file_name.startswith(start):
                        files_inside_folder.append(folder / file_name)    
                elif end:
                    if file_name.endswith(end):
                        files_inside_folder.append(folder / file_name)  
                elif start: 
                    if file_name.startswith(start) : 
                        files_inside_folder.append(folder / file_name)  
                else:
                    files_inside_folder.append(folder / file_name)  
        
                # If <file_name> starts with "run" then this refers to a sub_folder containing a simulation.
                if file_name.startswith("run") and os.path.isdir(file_name):
                    for file_name_run in os.listdir(folder / file_name):
                        if not file_name_run.startswith("."):
                            # If the <file_name> starts with <start> and ends with <end>, add it to files_inside.
                            if end and start:
                                if file_name_run.endswith(end) and file_name_run.startswith(start):
                                    files_inside_folder.append(folder / file_name / file_name_run)    
                            elif end:
                                if file_name_run.endswith(end):
                                    files_inside_folder.append(folder / file_name / file_name_run)  
                            elif start:
                                if file_name_run.startswith(start) :
                                    files_inside_folder.append(folder / file_name / file_name_run) 
                                    
                # If their are more folders, look into these as well
                if os.path.isdir(folder / file_name):
                    files = get_filesInFolder(folder / file_name, name_inputFile, start, end)
                    if files:
                        files_inside_folder += files

        # Add these files to the big list
        files_inside += files_inside_folder 
    
    # If no files were identified, return None
    if files_inside == []: return None
    else:                  return files_inside 

    
def get_hoursFromSeconds(time_inSeconds):
    ''' Convert the minutes to hours:minutes:seconds. '''
    hours   = int(time_inSeconds/3600)
    minutes = int((time_inSeconds/60) % 60)
    seconds = int((time_inSeconds) % 60)
    return "%d:%02d:%02d" % (hours, minutes, seconds)
  
def get_hoursFromMinutes(time_inMinutes):
    ''' Convert the minutes to hours:minutes:seconds. '''
    hours   = int(time_inMinutes/60)
    minutes = int(time_inMinutes % 60)
    seconds = int((time_inMinutes*60) % 60)
    return "%d:%02d:%02d" % (hours, minutes, seconds)

def get_hoursFromHours(time_inHours):
    ''' Convert the minutes to hours:minutes:seconds. '''
    hours   = int(time_inHours)
    minutes = int((time_inHours*60) % 60)
    seconds = int((time_inHours*3600) % 60) 
    return "%d:%02d:%02d" % (hours, minutes, seconds)
 
# Nested dictt
dictt = initiate_nesteddictt()
processed_files = 0
 
# Get the folder
directory = pathlib.Path(os.getcwd())
 
# Get the stella.out files which look like id_name.out
files = sorted(get_filesInFolder(directory, end=".out")) 
files = [ f for f in files if f.name[0:7].isdigit()] 
files = [ f for f in files if not os.path.isfile(f.with_suffix(".dt1.out"))]  
totalsimulations = len(files)
 
if len(files)<500:
    print("\n\n")
    print("  |"+"{0:^8}".format("DEVICE") + \
          "{0:^8}".format("FPRIM") + \
          "{0:^8}".format("TPRIM") + \
          "{0:^10}".format("RESOLUTION") + \
          "  |"+"{0:^10}".format("T_LAST") + \
          "{0:^10}".format("DT") + \
          "{0:^10}".format("RUN TIME") + \
          "{0:^10}".format("STEPS") + \
          "{0:^10}".format("TIME/STEP") + \
          "  |"+"{0:^8}".format("CORES") + \
          "{0:^10}".format("CPU TIME")+ \
          "{0:^17}".format("RUNTIME t=1000")+ \
          "{0:^17}".format("CPUTIME t=1000")+"|")
    print("  |"+"{0:^8}".format("-----") + \
          "{0:^8}".format("-----") + \
          "{0:^8}".format("-----") + \
          "{0:^10}".format("----------") + \
          "  |"+"{0:^10}".format("------") + \
          "{0:^10}".format("--") + \
          "{0:^10}".format("--------") + \
          "{0:^10}".format("-----") + \
          "{0:^10}".format("---------") + \
          "  |"+"{0:^8}".format("-----") + \
          "{0:^10}".format("-------")+ \
          "{0:^17}".format("--------------")+ \
          "{0:^17}".format("--------------")+"|") 
        
# Go through the files
for file in files:
    
    #try:    
    
        # Read the file 
        data = open(file, 'r') 
        slurmid = str(file.name[0:7])
        lines = data.readlines()
        if file.name.endswith("dt1.out"):
            cores_lines = "".join(lines[17:30])
            timer_lines = "".join(lines[-22:])
            laststep_lines = "".join(lines[-30:-25])  
        else:
            cores_lines = "".join(lines[22:27])
            timer_lines = "".join(lines[-22:])
            laststep_lines = "".join(lines[-30:-24]) 
        
        # Only look at files that contain the time information 
        if "total:" in timer_lines: 
 
            # If we have less than 100 lines then nstep<actual step
            indexstart = None; indexstop = None
            for i in range(len(lines)): 
                if "OVERVIEW OF THE SIMULATION" in lines[i]: indexstart = i
                if "ELAPSED TIME" in lines[-i]: indexstop = i
                if indexstart and indexstop:  
                    break 
            
            if len(lines)-indexstop-indexstart > 100: 
                  
                # Get the number of processors from stella.out
                cores = int(int(cores_lines.split("Running on")[-1].split("processors.")[0]))
                  
                # Get the simulation time from stella.out if it finished by stella itself
                runtime_min = float(timer_lines.split("total:")[-1].split("min")[0]) 
                 
                # Get the total number of steps and the reached time in the simulation  
                laststep_lines = laststep_lines.split("\n")[-2].split(" ")
                laststep_lines = [i for i in laststep_lines if i!='']
                steps = int(laststep_lines[0])
                simulationtime = float(laststep_lines[1])
                dt = float(laststep_lines[2])  
                 
                # Get the device 
                if "TJII" in str(file): device = "TJII"
                elif "W7X" in str(file): device = "W7X"
                elif "LHD" in str(file): device = "LHD"
                elif "CBC" in str(file): device = "CBC"
                elif "NCSX" in str(file): device = "NCSX"
                elif "ASDEX" in str(file): device = "ASDEX"
                 
                # Get the gradients 
                fprim = -1; tprim = -1 
                folder = file.parent.name if ("OLD" not in file.parent.name) else file.parent.parent.name
                if "tprim" in folder:
                    fprim = float(folder.split("fprim")[-1].split("tprim")[0])
                    tprim = float(folder.split("tprim")[-1].split("/")[0].split("_")[0])
                if "tiprim" in file.parent.name:
                    fprim = float(folder.split("fprim")[-1].split("tiprim")[0])
                    tprim = float(folder.split("tiprim")[-1].split("/")[0].split("_")[0]) 
                     
                # Resolution scan
                resolution = "/"
                if "double" in folder:
                    resolution =  folder.split("double")[-1].split("_")[0]
                    
                # Print the data
                if len(files)<500: 
                    
                    timeperstep_s = runtime_min*60/steps
                    cputime_hours = runtime_min/60*cores     
                    stepstoreacht1000 = 1000/dt 
                    runtimefor1000 = timeperstep_s*stepstoreacht1000/3600
                    cputimefor1000 = runtimefor1000*96
                      
                    strdt = str(round(dt*1000,1))+"ms"
                    strtimeperstep = str(round(timeperstep_s,4))+"s"
                    strruntime = get_hoursFromMinutes(runtime_min)
                    strcputime = get_hoursFromHours(cputime_hours)
                    strruntimefor1000 = get_hoursFromHours(runtimefor1000)  
                    strcputimefor1000 = get_hoursFromHours(cputimefor1000)  
                    
                    print("  |"+"{0:^8}".format(device) + \
                      "{0:^8}".format(fprim) + \
                      "{0:^8}".format(tprim) + \
                      "{0:^10}".format(resolution) + \
                      "  |"+"{0:^10}".format(simulationtime) + \
                      "{0:^10}".format(strdt) + \
                      "{0:^10}".format(strruntime) + \
                      "{0:^10}".format(steps) + \
                      "{0:^10}".format(strtimeperstep) + \
                      "  |"+"{0:^8}".format(cores) + \
                      "{0:^10}".format(strcputime) + \
                      "{0:^17}".format(strruntimefor1000) + \
                      "{0:^17}".format(strcputimefor1000)+"|")
                
                # Save the data
                processed_files += 1
                if "dt" not in dictt[device][resolution][tprim][fprim]:
                    dictt[device][resolution][tprim][fprim]["dt"] = 0
                    dictt[device][resolution][tprim][fprim]["sims"] = 0
                    dictt[device][resolution][tprim][fprim]["cores"] = 0
                    dictt[device][resolution][tprim][fprim]["steps"] = 0
                    dictt[device][resolution][tprim][fprim]["tlast"] = 0
                    dictt[device][resolution][tprim][fprim]["runtime_min"] = 0
                dictt[device][resolution][tprim][fprim]["dt"] += dt
                dictt[device][resolution][tprim][fprim]["sims"] += 1
                dictt[device][resolution][tprim][fprim]["cores"] += cores
                dictt[device][resolution][tprim][fprim]["steps"] += steps
                dictt[device][resolution][tprim][fprim]["tlast"] += simulationtime
                dictt[device][resolution][tprim][fprim]["runtime_min"] += runtime_min 
    #except:
        #print("SOMETHING WENT WRONG FOR "+str(file))
     
print("\n\n")
print("     ##########################################################################################################################################")
print("                                                    AVERAGE PER DEVICE; FPRIM; TPRIM; RESOLUTION")
print("     ##########################################################################################################################################")
print("                    ---> Examined "+str(processed_files)+" of the "+str(totalsimulations)+" simulations.")
print("\n")
print("  |"+"{0:^4}".format("SIMS") + \
      "{0:^8}".format("DEVICE") + \
      "{0:^8}".format("FPRIM") + \
      "{0:^8}".format("TPRIM") + \
      "{0:^10}".format("RESOLUTION") + \
      "  |"+"{0:^10}".format("T_LAST") + \
      "{0:^10}".format("DT") + \
      "{0:^10}".format("RUN TIME") + \
      "{0:^10}".format("STEPS") + \
      "{0:^10}".format("TIME/STEP") + \
      "  |"+"{0:^8}".format("CORES") + \
      "{0:^10}".format("CPU TIME")+ \
      "{0:^17}".format("RUNTIME t=1000")+ \
      "{0:^17}".format("CPUTIME t=1000")+"|")
print("  |"+"{0:^4}".format("----") + \
      "{0:^8}".format("-----") + \
      "{0:^8}".format("-----") + \
      "{0:^8}".format("-----") + \
      "{0:^10}".format("----------") + \
      "  |"+"{0:^10}".format("------") + \
      "{0:^10}".format("--") + \
      "{0:^10}".format("--------") + \
      "{0:^10}".format("-----") + \
      "{0:^10}".format("---------") + \
      "  |"+"{0:^8}".format("-----") + \
      "{0:^10}".format("-------")+ \
      "{0:^17}".format("--------------")+ \
      "{0:^17}".format("--------------")+"|")

for device in sorted(dictt.keys()):
    for resolution in sorted(dictt[device].keys()): 
        for tprim in sorted(dictt[device][resolution].keys()):
            for fprim in sorted(dictt[device][resolution][tprim].keys()): 
                sims = dictt[device][resolution][tprim][fprim]["sims"]
                dt = dictt[device][resolution][tprim][fprim]["dt"]/sims
                cores = int(dictt[device][resolution][tprim][fprim]["cores"]/sims) 
                steps = int(dictt[device][resolution][tprim][fprim]["steps"]/sims) 
                simulationtime = int(dictt[device][resolution][tprim][fprim]["tlast"]/sims)
                runtime_min = dictt[device][resolution][tprim][fprim]["runtime_min"]/sims
 
                # Print the data
                timeperstep_s = runtime_min*60/steps
                cputime_hours = runtime_min/60*cores     
                stepstoreacht1000 = 1000/dt 
                runtimefor1000 = timeperstep_s*stepstoreacht1000/3600
                cputimefor1000 = runtimefor1000*96
                  
                strdt = str(round(dt*1000,1))+"ms"
                strtimeperstep = str(round(timeperstep_s,4))+"s"
                strruntime = get_hoursFromMinutes(runtime_min)
                strcputime = get_hoursFromHours(cputime_hours)
                strruntimefor1000 = get_hoursFromHours(runtimefor1000)  
                strcputimefor1000 = get_hoursFromHours(cputimefor1000)  
                        
                print("  |"+"{0:^4}".format(sims) + \
                  "{0:^8}".format(device) + \
                  "{0:^8}".format(fprim) + \
                  "{0:^8}".format(tprim) + \
                  "{0:^10}".format(resolution) + \
                  "  |"+"{0:^10}".format(simulationtime) + \
                  "{0:^10}".format(strdt) + \
                  "{0:^10}".format(strruntime) + \
                  "{0:^10}".format(steps) + \
                  "{0:^10}".format(strtimeperstep) + \
                  "  |"+"{0:^8}".format(cores) + \
                  "{0:^10}".format(strcputime) + \
                  "{0:^17}".format(strruntimefor1000) + \
                  "{0:^17}".format(strcputimefor1000)+"|")
  
print("\n\n")
print("     ##########################################################################################################################################")
print("                                                         AVERAGE PER DEVICE AND RESOLUTION")
print("     ##########################################################################################################################################")
print("                    ---> Examined "+str(processed_files)+" of the "+str(totalsimulations)+" simulations.")


for device in sorted(dictt.keys()):
    print("\n")
    print("  |"+"{0:^4}".format("SIMS") + \
          "{0:^8}".format("DEVICE") + \
          "{0:^10}".format("RESOLUTION") + \
          "  |"+"{0:^10}".format("T_LAST") + \
          "{0:^10}".format("DT") + \
          "{0:^10}".format("RUN TIME") + \
          "{0:^10}".format("STEPS") + \
          "{0:^10}".format("TIME/STEP") + \
          "  |"+"{0:^8}".format("CORES") + \
          "{0:^10}".format("CPU TIME")+ \
          "{0:^17}".format("RUNTIME t=1000")+ \
          "{0:^17}".format("CPUTIME t=1000")+"|")
    print("  |"+"{0:^4}".format("----") + \
          "{0:^8}".format("-----") + \
          "{0:^10}".format("----------") + \
          "  |"+"{0:^10}".format("------") + \
          "{0:^10}".format("--") + \
          "{0:^10}".format("--------") + \
          "{0:^10}".format("-----") + \
          "{0:^10}".format("---------") + \
          "  |"+"{0:^8}".format("-----") + \
          "{0:^10}".format("-------")+ \
          "{0:^17}".format("--------------")+ \
          "{0:^17}".format("--------------")+"|")

    for resolution in sorted(dictt[device].keys()):
        sims = 0; dt = 0; cores = 0; steps = 0; simulationtime = 0; runtime_min = 0
        for tprim in sorted(dictt[device][resolution].keys()):
            for fprim in sorted(dictt[device][resolution][tprim].keys()): 
                sims += dictt[device][resolution][tprim][fprim]["sims"]
                dt += dictt[device][resolution][tprim][fprim]["dt"]
                cores += dictt[device][resolution][tprim][fprim]["cores"] 
                steps += dictt[device][resolution][tprim][fprim]["steps"] 
                simulationtime += dictt[device][resolution][tprim][fprim]["tlast"] 
                runtime_min += dictt[device][resolution][tprim][fprim]["runtime_min"]
        dt = dt/sims
        cores = int(cores/sims)
        steps = int(steps/sims)
        simulationtime = int(simulationtime/sims)
        runtime_min = runtime_min/sims
 
        # Print the data
        timeperstep_s = runtime_min*60/steps
        cputime_hours = runtime_min/60*cores     
        stepstoreacht1000 = 1000/dt 
        runtimefor1000 = timeperstep_s*stepstoreacht1000/3600
        cputimefor1000 = runtimefor1000*96
          
        strdt = str(round(dt*1000,1))+"ms"
        strtimeperstep = str(round(timeperstep_s,4))+"s"
        strruntime = get_hoursFromMinutes(runtime_min)
        strcputime = get_hoursFromHours(cputime_hours)
        strruntimefor1000 = get_hoursFromHours(runtimefor1000)  
        strcputimefor1000 = get_hoursFromHours(cputimefor1000)  
        print("  |"+"{0:^4}".format(sims) + \
          "{0:^8}".format(device) + \
          "{0:^10}".format(resolution) + \
          "  |"+"{0:^10}".format(simulationtime) + \
          "{0:^10}".format(strdt) + \
          "{0:^10}".format(strruntime) + \
          "{0:^10}".format(steps) + \
          "{0:^10}".format(strtimeperstep) + \
          "  |"+"{0:^8}".format(cores) + \
          "{0:^10}".format(strcputime) + \
          "{0:^17}".format(strruntimefor1000) + \
          "{0:^17}".format(strcputimefor1000)+"|")

  
print("\n\n")
print("     ##########################################################################################################################################")
print("                                                         AVERAGE OF ALL SIMULATIONS")
print("     ##########################################################################################################################################")
print("                    ---> Examined "+str(processed_files)+" of the "+str(totalsimulations)+" simulations.")
print("\n")
print("  |"+"{0:^4}".format("SIMS") + \
      "  |"+"{0:^10}".format("T_LAST") + \
      "{0:^10}".format("DT") + \
      "{0:^10}".format("RUN TIME") + \
      "{0:^10}".format("STEPS") + \
      "{0:^10}".format("TIME/STEP") + \
      "  |"+"{0:^8}".format("CORES") + \
      "{0:^10}".format("CPU TIME")+ \
      "{0:^17}".format("RUNTIME t=1000")+ \
      "{0:^17}".format("CPUTIME t=1000")+"|")
print("  |"+"{0:^4}".format("----") + \
      "  |"+"{0:^10}".format("------") + \
      "{0:^10}".format("--") + \
      "{0:^10}".format("--------") + \
      "{0:^10}".format("-----") + \
      "{0:^10}".format("---------") + \
      "  |"+"{0:^8}".format("-----") + \
      "{0:^10}".format("-------")+ \
      "{0:^17}".format("--------------")+ \
      "{0:^17}".format("--------------")+"|")


sims = 0; dt = 0; cores = 0; steps = 0; simulationtime = 0; runtime_min = 0
for device in sorted(dictt.keys()):
    for resolution in sorted(dictt[device].keys()):
        for tprim in sorted(dictt[device][resolution].keys()):
            for fprim in sorted(dictt[device][resolution][tprim].keys()): 
                sims += dictt[device][resolution][tprim][fprim]["sims"]
                dt += dictt[device][resolution][tprim][fprim]["dt"]
                cores += dictt[device][resolution][tprim][fprim]["cores"] 
                steps += dictt[device][resolution][tprim][fprim]["steps"] 
                simulationtime += dictt[device][resolution][tprim][fprim]["tlast"] 
                runtime_min += dictt[device][resolution][tprim][fprim]["runtime_min"]
                
dt = dt/sims
cores = int(cores/sims)
steps = int(steps/sims)
simulationtime = int(simulationtime/sims)
runtime_min = runtime_min/sims

# Print the data
timeperstep_s = runtime_min*60/steps
cputime_hours = runtime_min/60*cores     
stepstoreacht1000 = 1000/dt 
runtimefor1000 = timeperstep_s*stepstoreacht1000/3600
cputimefor1000 = runtimefor1000*96
  
strdt = str(round(dt*1000,1))+"ms"
strtimeperstep = str(round(timeperstep_s,4))+"s"
strruntime = get_hoursFromMinutes(runtime_min)
strcputime = get_hoursFromHours(cputime_hours)
strruntimefor1000 = get_hoursFromHours(runtimefor1000)  
strcputimefor1000 = get_hoursFromHours(cputimefor1000)  
print("  |"+"{0:^4}".format(sims) + \
  "  |"+"{0:^10}".format(simulationtime) + \
  "{0:^10}".format(strdt) + \
  "{0:^10}".format(strruntime) + \
  "{0:^10}".format(steps) + \
  "{0:^10}".format(strtimeperstep) + \
  "  |"+"{0:^8}".format(cores) + \
  "{0:^10}".format(strcputime) + \
  "{0:^17}".format(strruntimefor1000) + \
  "{0:^17}".format(strcputimefor1000)+"|")

    
print("\n\n")



    







