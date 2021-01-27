#!/usr/bin/python3
''' When launching a simulation, check some variables in the input file and create a custom slurm file. '''

# Tell python where to find the personal modules since this script is run from the terminal
import sys, os, pathlib, time
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.utils import get_filesInFolder 
from stellapy.commands.utils.get_bashArguments import get_bashArguments
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.utils.decorators import print_functionInformationOnCommandPrompt

#################################################################
#                     MAIN PROGRAM
#################################################################

def launch(
        # Execute folder
        folder,
        # Simulation
        nodes=3,\
        github=False,\
        wall_time="22:00:00",\
        # Supercomputer
        account="FUA34_KINCIEMA",\
        partition="skl_fua_prod",\
        # Personal 
        email="hanne.thienpondt@outlook.com"):
    
    # Check whether the correct files are present
    input_file, default_text = check_presenceFiles(folder)
         
    # Check whether it is a linear or nonlinear simulation
    text_nonlinear = default_text.split("nonlinear")[-1].split("\n")[0]
    nonlinear = True if (".true." in text_nonlinear) else False
    
    # If we have a nonlinear simulation, make sure there is a restart folder
    if nonlinear:
        restart_folder = None; count=1
        while (restart_folder==None):
            if "restart"+str(count) not in os.listdir(folder):
                restart_folder = "restart"+str(count)
                os.system("mkdir restart"+str(count))
            count += 1
            
    # Edit the input file a bit
    create_inputFile(folder, input_file, default_text, nonlinear, restart_folder)
    
    # Create the eu.slurm file
    create_euSlurmFile(nonlinear, nodes, folder, input_file, github, wall_time, account, partition, email) 
 
    # Now launch the simulation
    print()
    os.system("sbatch -D " + str(folder) +" eu.slurm") 
    print()
    return

#################################################################
#                     CREATE THE INPUT FILE
#################################################################

def create_inputFile(folder, input_file, default_text, nonlinear, restart_folder):
    
    # If ginit_option==noise then change the restart folder:
    if nonlinear:
        
        # Create the new input file
        new_file = open(input_file, "w" )
        
        # Replace ginit_option
        text_before = default_text.split("ginit_option")[0]
        text_after  = default_text.split("ginit_option")[-1].split("\n")[1:]
        text_after  = '\n'.join(text_after)
        new_text    = text_before + 'ginit_option = "noise"' + "\n" + text_after
        
        # Replace restart directory
        text_before = new_text.split("restart_dir")[0]
        text_after  = new_text.split("restart_dir")[-1].split("\n")[1:]
        text_after  = '\n'.join(text_after)
        new_text    = text_before + 'restart_dir = "'+restart_folder+'"' + "\n" + text_after
        
        # Replace delt_option with delt
        if "delt_option" in new_text: 
            if "delt_option=" in new_text:   delt = "delt_option="
            if "delt_option =" in new_text:  delt = "delt_option ="
            if "delt_option  =" in new_text: delt = "delt_option  ="
            text_before = new_text.split(delt)[0]
            text_after  = new_text.split(delt)[-1].split("\n")[1:]
            text_after  = '\n'.join(text_after)
            new_text    = text_before + 'delt = 0.1' + "\n" + text_after
        
        # Make sure write moments is off
        if "write_moments" in new_text:
            text_before = new_text.split("write_moments")[0]
            text_after  = new_text.split("write_moments")[-1].split("\n")[1:]
            text_after  = '\n'.join(text_after)
            new_text    = text_before + 'write_moments = .false.' + "\n" + text_after
            
        # Make sure write moments is off
        if "write_omega" in new_text:
            text_before = new_text.split("write_omega")[0]
            text_after  = new_text.split("write_omega")[-1].split("\n")[1:]
            text_after  = '\n'.join(text_after)
            new_text    = text_before + 'write_omega = .false.' + "\n" + text_after

        # Write the new input file
        new_file.write(new_text) 
        new_file.close() 

    return
                    
#################################################################
#          CHECK WHETHER THE CORRECT FILES ARE PRESENT
#################################################################

def check_presenceFiles(folder):
    
    # Select the default input file
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [f for f in input_files if "_dt" not in str(f)]
    
    # If there are multiple default input files: exit
    if len(input_files) > 1:
        print(" Multiple default input files were found inside "+str(folder)+":")
        for f in input_files:
            print("     ", f)
        print(" Please make sure there is only one default input file.\n")
        return
    
    # If there are no input_files: exit
    if len(input_files) == 0:
        print(" No input file found in "+str(folder)+".\n")
        return
    
    # Read the magnetic field file name and make sure the file is in the folder
    input_data = open(input_files[0], 'r')
    input_text = input_data.read() 
    vmec_file  = str(input_text.split("vmec_filename")[-1].split("\n")[0])
    vmec_file  = vmec_file.replace(" ", "").replace("=", "").replace("'", "")
    if not os.path.isfile(folder / vmec_file) and not os.path.islink(folder / vmec_file):
        print(" The vmec file "+vmec_file+" was not found in "+str(folder)+".\n")
        return
    
    # Make sure we have a stella executable
    if not os.path.isfile(folder / "stella") and not os.path.islink(folder / "stella"):
        print(" No stella executable found in "+str(folder)+".\n")
        return
    
    # Read the default input file
    input_file   = input_files[0]
    default_file = open(input_file, "r" )
    default_text = default_file.read()
    default_file.close()
    
    # If we made it true, return the input file text
    return input_file, default_text
    
#################################################################
#                     EU.SLURM FILE
#################################################################

def create_euSlurmFile(\
        # Simulation
        nonlinear=False,\
        nodes=3,\
        folder="$RUNS",\
        input_file="input.in",\
        github=False,\
        wall_time="22:00:00",\
        # Supercomputer
        account="FUA34_KINCIEMA",\
        partition="skl_fua_prod",\
        # Personal 
        email="hanne.thienpondt@outlook.com"): 
    
    # Decide the wall time based on whether it is a linear or nonlinear simulation
    if wall_time=="22:00:00":
        wall_time = "22:00:00" if nonlinear else "1:00:00"
    
    # Open a new eu.slurm file 
    new_file = open(folder / "eu.slurm", "w" )
    
    # Calculate the necesarry variables
    nodes = int(nodes)
    cores = str(int(nodes*48))
    stella_link = folder / "stella"
    try: compilation_data = time.ctime(os.path.getmtime(stella_link))
    except: compilation_data = "Symbolic link is broken"
    stella_file = stella_link
    while os.path.islink(stella_file):
        stella_file = os.readlink(stella_file)
    
    # Set the sbatch settings
    sbatch_settings = '\n'.join((   
        '#!/bin/bash',
        '#SBATCH -N '+str(nodes),
        '#SBATCH -A '+account,
        '#SBATCH -p '+partition,
        '#SBATCH --time '+wall_time,
        '#SBATCH --job-name=stella',
        '#SBATCH --mail-type=ALL',
        '#SBATCH --mail-user='+email,
        ))
    
    path_settings = '\n'.join((   
        '\n',
        '# Run the simulation inside SLURM_SUBMIT_DIR',
        'SLURM_SUBMIT_DIR='+str(folder),
        'cd ${SLURM_SUBMIT_DIR}\n', 
        '# Set the path variables',
        'EXEFILE=./stella',
        'OUTFILE="'+input_file.stem+'.stella."$SLURM_JOBID',
        'ERRFILE="'+input_file.stem+'.err."$SLURM_JOBID',
        ))
    
    information = '\n'.join((   
        '\n',
        'echo " "',
        'echo "#########################################################"',
        'echo "                   STELLA SIMULATION                     "',
        'echo "#########################################################"',
        'echo "RUN DATE:          "$(date)',
        'echo "COMPILATION DATE:  '+compilation_data+'"',
        'echo "STELLA DIRECTORY:  '+stella_file+'"',
        'echo "RUN DIRECTORY:     '+str(folder)+'"',
        'echo "INPUT FILE:        '+input_file.name+'"',
        'echo "NODES:             '+str(nodes)+'"',
        'echo "CORES:             '+cores+'"',
        'echo "PATH:              "$PATH',
        'echo " "',
        ))

    if github==False: 
        load_modules = '\n'.join((  
            '\n',
            'echo "#########################################################"',
            'echo "                         MODULES                         "',
            'echo "#########################################################"',
            'echo "module purge"',
            'module purge',
            'echo "load intel/pe-xe-2018--binary"',
            'module load intel/pe-xe-2018--binary',
            'echo "load intelmpi/2018--binary"',
            'module load intelmpi/2018--binary',
            'echo "load zlib/1.2.8--gnu--6.1.0"',
            'module load zlib/1.2.8--gnu--6.1.0',
            'echo "load szip/2.1--gnu--6.1.0"',
            'module load szip/2.1--gnu--6.1.0',
            'echo "load hdf5/1.10.4--intel--pe-xe-2018--binary"',
            'module load hdf5/1.10.4--intel--pe-xe-2018--binary',
            'echo "load netcdf/4.6.1--intel--pe-xe-2018--binary"',
            'module load netcdf/4.6.1--intel--pe-xe-2018--binary',
            'echo "load netcdff/4.4.4--intel--pe-xe-2018--binary"',
            'module load netcdff/4.4.4--intel--pe-xe-2018--binary',
            'echo "load fftw/3.3.7--intelmpi--2018--binary"',
            'module load fftw/3.3.7--intelmpi--2018--binary',
            'echo "load lapack"',
            'module load lapack',
            'echo "load blas"',
            'module load blas', 
            'echo " "',
        ))

    if github==True:
        load_modules = '\n'.join((  
            '\n',
            'echo "#########################################################"',
            'echo "                         MODULES                         "',
            'echo "#########################################################"',
            'echo "module purge"',
            'module purge',
            'echo "load env-skl"',
            'module load env-skl',
            'echo "load intel/pe-xe-2018--binary"',
            'module load intel/pe-xe-2018--binary',
            'echo "load intelmpi/2018--binary"',
            'module load intelmpi/2018--binary',
            'echo "load szip/2.1--gnu--6.1.0"',
            'module load szip/2.1--gnu--6.1.0',
            'echo "load zlib/1.2.8--gnu--6.1.0"',
            'module load zlib/1.2.8--gnu--6.1.0',
            'echo "load hdf5/1.10.4--intel--pe-xe-2018--binary"',
            'module load hdf5/1.10.4--intel--pe-xe-2018--binary',
            'echo "load netcdf/4.6.1--intel--pe-xe-2018--binary"',
            'module load netcdf/4.6.1--intel--pe-xe-2018--binary',
            'echo "load netcdff/4.4.4--intel--pe-xe-2018--binary"',
            'module load netcdff/4.4.4--intel--pe-xe-2018--binary',
            'echo "load fftw/3.3.7--intelmpi--2018--binary"',
            'module load fftw/3.3.7--intelmpi--2018--binary',
            'echo "load mkl/2018--binary "',
            'module load mkl/2018--binary ',
            'echo "load scalapack"',
            'module load scalapack',
            'echo "load blas"',
            'module load blas', 
            'echo " "',
        ))
    
    execute_simulation = '\n'.join((   
        '\n',
        'echo "#########################################################"',
        'echo "                 EXECUTE SIMULATION                     "',
        'echo "#########################################################"',
        'mpirun -errfile-pattern $ERRFILE -outfile-pattern $OUTFILE'+\
        ' -envall -genv -n '+cores+' $EXEFILE '+input_file.name,
        '\n',
        '\n',
    ))
    
    # Combine the texts
    euslurm_text = sbatch_settings+path_settings+information+load_modules+execute_simulation

    # Write the eu.slurm file
    new_file.write(euslurm_text)

#################################################################
#                     BASH INTERFACE
#################################################################

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        function = "runSimulations.launch"
        default_arguments = {\
                'folder': None,\
                'github' : False,\
                'nodes' : 3,\
                'wall_time' : "22:00:00"}
        arguments = get_arguments(function, default_arguments) 
        launch(**arguments) 

# Explain the program
def print_help(function, default_arguments):
    print_decoratorOnCommandPrompt("begin")
    print_functionInformationOnCommandPrompt(function, default_arguments, None)
    print("\n The arguments are:") 
    print("      -g         Load modules for the new code on github" )
    print("      -m         Wall_time in minutes" )  
    print("      -h         Wall_time in hours" )  
    print("      -n         Number of nodes" )  
    print_decoratorOnCommandPrompt("end")

# Get the arguments from the terminal
def get_arguments(function, default_arguments):

    # Make a new dictionary which will hold the new arguments
    arguments = default_arguments.copy()

    # Initiate mandetory variables
    arguments['folder'] = pathlib.Path(os.getcwd())
    
    # Read the options and arguments from the terminal
    options = "gm:h:n:"
    long_options = ["folder=","help"]
    opts_args = get_bashArguments(options, long_options, print_help, function, default_arguments)

    # Asign the options and arguments to variables
    for opt, arg in opts_args:
        if opt == '--folder':
            arguments['folder'] = str(arg)
        elif opt == '-g':
            arguments['github'] = True
        elif opt == '-n':
            arguments['nodes'] = int(arg)
        elif opt == '-m':
            if int(arg)<10:  arguments['wall_time'] = "00:0"+str(int(arg))+":00" 
            if int(arg)>=10: arguments['wall_time'] = "00:"+str(int(arg))+":00" 
        elif opt == '-h':
            if int(arg)<10:  arguments['wall_time'] = "0"+str(int(arg))+":00:00" 
            if int(arg)>=10: arguments['wall_time'] = str(int(arg))+":00:00" 
        elif opt in ("--help"):
            print_help(function, default_arguments)
            sys.exit()
            
    return arguments

# Execute the python script
main_program()
