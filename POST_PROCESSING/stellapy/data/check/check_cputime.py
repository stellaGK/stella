#!/usr/bin/python3
import glob
import sys, os 
import pathlib
import subprocess
 
# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)
from stellapy.utils.commandprompt.print_progressbar import print_progressbar
from stellapy.utils.commandprompt.bash import Bash
 
#===============================================================================
#                                CHECK CPU TIME                                #
#===============================================================================
# sacct -u hthienpo --format "Account,JobID%15,JobName%20,State,ExitCode,Submit,CPUTime" -S 2022-02-10
# sacct -j 10960378 --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss,MaxVMSize,nnodes,ncpus,nodelist
# sacct -j 10960378 --format=User,JobID,Jobname,state,time,start,end,elapsed,nnodes,ncpus,nodelist
# sacct --helpformat
# TODO: if the simulation was running, rewrite the sacct!
 
def check_cputime(folder, verbose=False):

    # Get the subfolders
    subfolders = get_subdirectories(folder) 
    slurm_infos = {}; input_files_and_slurm_ids = {}
    
    # Show reading progress 
    if verbose: count = print_progressbar(0, len(subfolders), prefix='   Reading slurm info:', suffix='Complete', length=50)+1

    # Iterate over the subfolders
    for subfolder in subfolders:
        
        # Show the progress
        if verbose: print_progressbar(count, len(subfolders), prefix='   Reading slurm info:', suffix=subfolder.name, length=50, endWithNewLine=False); count += 1      

        # Find the slurm id numbers
        slurm_files = [pathlib.Path(i) for i in glob.glob(str(subfolder/'slurm-*.out'))]
        slurm_files += [pathlib.Path(i) for i in glob.glob(str(subfolder/'slurm-*.info'))]
        slurm_ids = list(set([int(i.name.replace('slurm-','').replace('.out','').replace('.info','')) for i in slurm_files]))
        if len(slurm_files)==0: continue 
        
        # Get the matching input file
        input_file = pathlib.Path(glob.glob(str(subfolder/'*.in'))[0])
        input_files_and_slurm_ids[input_file] = slurm_ids
        
        # Get the elasped time
        for slurm_id in slurm_ids: 
            
            # Read the slurm info   
            path_slurm_info = subfolder/('slurm-'+str(slurm_id)+'.info')
            if os.path.isfile(path_slurm_info):
                with open(path_slurm_info) as f:
                    data = f.readlines() 
                    data = [line.replace('\n','') for line in data]
            if not os.path.isfile(path_slurm_info):
                data = write_slurmInfo(path_slurm_info, slurm_id)
                
            # Save the slurm info
            for line in data:
                line = [i for i in line.split(' ') if i!=''] 
                if line[1]==str(slurm_id):  
                    slurm_infos[str(slurm_id)] = {
                        'directory' : subfolder,
                        'job_name' : line[2],
                        'status' : line[3],
                        'wall_time' : line[4], 
                        'elapsed_time' : line[7], 
                        'nodes' : line[8], 
                        'cpus' : line[9]}
                    if slurm_infos[str(slurm_id)]['status']=='RUNNING':
                        try: data = write_slurmInfo(path_slurm_info, slurm_id)
                        except: pass
                    break   
    
    # Finish reading progress 
    if verbose: print_progressbar(len(subfolders), len(subfolders), prefix='   Reading slurm info:', suffix='Complete', length=50)

    # Print the data 
    if verbose:
        print() 
        print("".center(80,"="))
        print("Slurm info".center(80))
        print("".center(80,"="))
        print()
        print("     "+
              "{0:^15}".format("WALL TIME") + \
              "{0:^15}".format("ELASPED TIME") + \
              "{0:^10}".format("STATUS") + \
              "{0:^10}".format("NODES") + \
              "{0:^10}".format("CPUs")+ \
              "{0:^10}".format("COST [kh]")+\
              "    "+"{0:<35}".format("DIRECTORY"))
        print("     "+
              "{0:^15}".format("---------") + #WALL TIME
              "{0:^15}".format("------------") + #ELASPED TIME
              "{0:^10}".format("-----") + #STATUS
              "{0:^10}".format("-----") + #NODES
              "{0:^10}".format("----") + #CPUs
              "{0:^10}".format("---------") + #COST [kh]
              "    "+"{0:<35}".format("---------")) #DIRECTORY
        running_simulations = 0
        total_cpu_time = 0
        
        # Slurm id keys
        slurm_ids = list(slurm_infos.keys()) 
        directories = [slurm_infos[slurm_id]['directory'].name for slurm_id in slurm_ids]
        directories = [directories[i]+str(i) for i in range(len(directories))]
        sorted_indices = [directories.index(i) for i in sorted(directories)]
        slurm_ids = [slurm_ids[i] for i in sorted_indices]
        
        # Iterate over the slurm ids
        for slurm_id in slurm_ids:  
         
            # Make sure the file name isn't too long
            directory = slurm_infos[slurm_id]['directory'].name
            directory = directory+'/'+slurm_infos[slurm_id]['job_name']
            if len(directory) > 60: directory = directory[0:60] 
            
            # Cost of the simulation
            elapsed_time = [int(i) for i in slurm_infos[slurm_id]['elapsed_time'].split(':')] 
            cpu_time = (elapsed_time[0] + elapsed_time[1]*1/60 + elapsed_time[2]*1/3600)*int(slurm_infos[slurm_id]['cpus'])
            if slurm_infos[slurm_id]['status']=='RUNNING': running_simulations += 1
            total_cpu_time += cpu_time 
            
            
            # Print the simulation time
            print("     "+\
                  "{0:^15}".format(slurm_infos[slurm_id]['wall_time']) + \
                  "{0:^15}".format(slurm_infos[slurm_id]['elapsed_time']) + \
                  "{0:^10}".format(slurm_infos[slurm_id]['status']) + \
                  "{0:^10}".format(slurm_infos[slurm_id]['nodes']) + \
                  "{0:^10}".format(slurm_infos[slurm_id]['cpus'])+ \
                  "{0:^10.2f}".format(cpu_time/1000) +\
                  "    "+"{0:<60}".format(directory))
        
        if len(slurm_ids)>1:
            print("     "+
                  "{0:^15}".format("---------") + #WALL TIME
                  "{0:^15}".format("------------") + #ELASPED TIME
                  "{0:^10}".format("-----") + #STATUS
                  "{0:^10}".format("-----") + #NODES
                  "{0:^10}".format("----") + #CPUs
                  "{0:^10}".format("---------") + #COST [kh]
                  "    "+"{0:<50}".format("---------"))  #DIRECTORY
            print("     "+
                  "{0:^15}".format('') + \
                  "{0:^15}".format('') + \
                  "{0:^10}".format(running_simulations) + \
                  "{0:^10}".format('') + \
                  "{0:^10}".format('')+ \
                  "{0:^10.2f}".format(total_cpu_time/1000) + \
                  "    "+"{0:<50}".format(''))
        print()
        
    # Return the slurm info
    return input_files_and_slurm_ids, slurm_infos

#===============================================================================
#                               WRITE SLURM INFO
#===============================================================================

def write_slurmInfo(path_slurm_info, slurm_id):
    args = ['sacct', '-j', str(slurm_id), '--format=User,JobID,Jobname,state,time,start,end,elapsed,nnodes,ncpus']
    data = subprocess.Popen(args, stdout=subprocess.PIPE).communicate()[0]
    try: data.decode('utf-8')
    except: data.decode('latin-1')
    with open(path_slurm_info, 'w') as f:
        f.write(str(data, 'utf-8'))
    data = str(data, 'utf-8').split('\\n')
    return data
    
    
#===============================================================================
#                              GET SUB DIRECTORIES
#===============================================================================

def get_subdirectories(directory):
    subdirectories = [x[0] for x in os.walk(directory, followlinks=True)] 
    subdirectories = [pathlib.Path(f) for f in subdirectories]
    subdirectories = [f for f in subdirectories if f.name!="restart"] 
    return subdirectories

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================

if __name__ == "__main__":
    bash = Bash(check_cputime, __doc__)  
    args = bash.get_arguments()
    args['verbose'] = True
    check_cputime(**args)   
