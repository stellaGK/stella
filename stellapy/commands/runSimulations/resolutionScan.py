#!/usr/bin/python3
''' Run kyscan(identifier) as a bash command on the command prompt.'''

# Tell python where to find the personal modules since this script is run from the terminal
import sys, os, configparser, pathlib, time
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.utils import get_filesInFolder 
from stellapy.commands.utils.get_bashArguments import get_bashArguments
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.utils.decorators import print_functionInformationOnCommandPrompt

#################################################################
#                     MAIN PROGRAM
#################################################################


def resolutionScan(
        # Execute folder
        folder,
        # File containing the ky of the modes
        ky_path=None,\
        # Different time steps and gradients
        vec_mu=[6, 12, 24, 36, 48],\
        vec_vgrid=[6, 12, 24, 36, 48, 60],\
        vec_nzed=[32, 64, 128, 256, 512, 1024],\
        fprims=[1.0, 0.0, 4.0, 4.0, 8.0],\
        tiprims=[1.0, 4.0, 0.0, 4.0, 8.0],\
        end_time=500,\
        # Simulation
        nodes=3,\
        wall_time="1:00:00",\
        run_directory="$RUNS",\
        input_file="input.in",\
        # Supercomputer
        account="FUA34_KINCIEMA",\
        partition="skl_fua_prod",\
        # Personal 
        email="hanne.thienpondt@outlook.com"):
    ''' Scan delta t for the most unstable mode and different values of (fprim, tiprim). '''
    
    # Check whether the correct files are present
    input_file, default_text, vmec_file = check_presenceFiles(folder)
    
    # Count the number of simulations
    count = 0
        
    # For each (fprim, tiprim) create a folder, add the files and run the simulations
    for i in range(len(fprims)):
        fprim = fprims[i]
        tiprim = tiprims[i]
        
        # Scan the three paramaters
        for parameter in ["nmu", "nvgrid","nzed"]:
            
            # Read the ky of the most unstable mode at (fprim, tiprim) based on the "linearmap.ini" file
            ky_file = configparser.ConfigParser() 
            ky_path = get_filesInFolder(folder, start="linearmap")[0] if ky_path==None else ky_path
            ky_file.read(ky_path)
            parameters = '(' + str(fprim) + ", " + str(tiprim) + ')'
            if parameters in ky_file["Ky of the most unstable mode"].keys():
                ky = ky_file["Ky of the most unstable mode"][parameters] 
                 
            # Make a folder for this simulation
            run_directory = parameter+"scan_fprim"+str(int(fprim))+"tprim"+str(int(tiprim))
            os.system("mkdir "+run_directory)
            print("\n------------------------------------------")
            print("\nNEW FOLDER: "+run_directory)
     
            # Add symbolic links to the stella code and the magnetic field file inside this subfolder
            if not os.path.islink(folder / run_directory / "stella"):
                os.symlink(folder / "stella", folder / str(run_directory+"/stella"))
                os.symlink(folder / vmec_file, folder / str(run_directory+"/"+vmec_file))
             
            # For each time step, create an input file and execute it
            if parameter=="nmu":    variables=vec_mu
            if parameter=="nvgrid": variables=vec_vgrid
            if parameter=="nzed":   variables=vec_nzed
            for variable in variables:
                 
                # Create the input file: replace the ky and the dt value
                file_name = create_inputFile(folder, run_directory, input_file, default_text, variable, ky, fprim, tiprim, end_time, parameter)
             
                # Create the eu.slurm file
                create_euSlurmFile(nodes, wall_time, folder, run_directory, file_name, account, partition, email) 
             
                # Now launch the simulation
                print("sbatch -D " + str(folder)+"/"+run_directory +" eu.slurm    --->     "+file_name)
                os.system("sbatch -D " + str(folder)+"/"+run_directory +" eu.slurm"); count += 1

    print("\nFinished launching the "+str(count)+" modes.\n")
    return

#################################################################
#                     CREATE THE INPUT FILE
#################################################################

def create_inputFile(folder, run_directory, input_file, default_text, variable, ky, fprim, tiprim, end_time, parameter):
    # Create the new input file
    file_name = str(input_file.stem)+"_"+parameter+str(variable)+".in"
    new_file = open(input_file.parent / run_directory / file_name, "w" )
    # Replace aky_min
    text_before = default_text.split("aky_min")[0]
    text_after  = default_text.split("aky_min")[-1].split("\n")[1:]
    text_after  = '\n'.join(text_after)
    new_text    = text_before + "aky_min = " + ky + "\n" + text_after
    # Replace aky_max
    text_before = new_text.split("aky_max")[0]
    text_after  = new_text.split("aky_max")[-1].split("\n")[1:]
    text_after  = '\n'.join(text_after)
    new_text    = text_before + "aky_max = " + ky + "\n" + text_after    
    # Replace tprim of species 1
    text_before = new_text.split("species_parameters_1")[0] + "species_parameters_1"
    text_specie = new_text.split("species_parameters_1")[-1].split("/")[0] + '/' 
    text_after  = new_text.split("species_parameters_1")[-1].split("/")[1:]  
    text_after  = '/'.join(text_after)
    specie_bef  = text_specie.split("tprim")[0]
    specie_aft  = text_specie.split("tprim")[-1].split("\n")[1:]
    specie_aft  = '\n'.join(specie_aft)
    new_text    = text_before + specie_bef + "tprim = " + str(tiprim) + "\n" + specie_aft + text_after
    # Replace fprim of species 1
    text_before = new_text.split("species_parameters_1")[0] + "species_parameters_1"
    text_specie = new_text.split("species_parameters_1")[-1].split("/")[0] + '/' 
    text_after  = new_text.split("species_parameters_1")[-1].split("/")[1:]  
    text_after  = '/'.join(text_after)
    specie_bef  = text_specie.split("fprim")[0]
    specie_aft  = text_specie.split("fprim")[-1].split("\n")[1:]
    specie_aft  = '\n'.join(specie_aft)
    new_text    = text_before + specie_bef + "fprim = " + str(fprim) + "\n" + specie_aft + text_after
    # Replace fprim of species 2
    text_before = new_text.split("species_parameters_2")[0] + "species_parameters_2"
    text_specie = new_text.split("species_parameters_2")[-1].split("/")[0] + '/' 
    text_after  = new_text.split("species_parameters_2")[-1].split("/")[1:] 
    text_after  = '/'.join(text_after)
    specie_bef  = text_specie.split("fprim")[0]
    specie_aft  = text_specie.split("fprim")[-1].split("\n")[1:]
    specie_aft  = '\n'.join(specie_aft)
    new_text    = text_before + specie_bef + "fprim = " + str(fprim) + "\n" + specie_aft + text_after
    # Replace delt
    if parameter+"=" in new_text:   delt = parameter+"="
    if parameter+" =" in new_text:  delt = parameter+" ="
    if parameter+"  =" in new_text: delt = parameter+"  ="
    text_before = new_text.split(delt)[0]
    text_after  = new_text.split(delt)[-1].split("\n")[1:]
    text_after  = '\n'.join(text_after)
    new_text    = text_before + parameter+" = " + str(variable) + "\n" + text_after
    # Replace nstep
    text_before = new_text.split("nstep")[0]
    text_after  = new_text.split("nstep")[-1].split("\n")[1:]
    text_after  = '\n'.join(text_after)
    new_text    = text_before + "nstep = " + str(int(end_time/0.0075)) + "\n" + text_after
    new_file.write(new_text) 
    new_file.close() 
    return file_name
                    
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
    
    # Make sure n_aky and n_akx are set to 1
    nakx = int(default_text.split("nakx")[-1].split("\n")[0].replace("=", ""))
    naky = int(default_text.split("naky")[-1].split("\n")[0].replace("=", ""))
    if nakx!=1 or naky!=1:
        print(" Please set nakx and naky to 1 if you want to run one mode per input file.\n")
        return
    
    # If we made it true, return the input file text
    return input_file, default_text, vmec_file
    
#################################################################
#                     EU.SLURM FILE
#################################################################

def create_euSlurmFile(\
        # Simulation
        nodes=3,\
        wall_time="1:00:00",\
        folder="$RUNS",\
        run_directory="Linearmap",\
        input_file="input.in",\
        # Supercomputer
        account="FUA34_KINCIEMA",\
        partition="skl_fua_prod",\
        # Personal 
        email="hanne.thienpondt@outlook.com"): 
    
    # Open a new eu.slurm file 
    new_file = open(folder / str(run_directory+"/eu.slurm"), "w" )
    
    # Calculate the necesarry variables
    nodes = int(nodes)
    cores = str(int(nodes*48))
    run_directory = str(run_directory)
    stella_file = folder / str(run_directory+"/stella")
    try: compilation_data = time.ctime(os.path.getmtime(stella_file))
    except: compilation_data = "Symbolic link is broken"
    
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
        'SLURM_SUBMIT_DIR='+str(folder/run_directory),
        'cd ${SLURM_SUBMIT_DIR}\n', 
        '# Set the path variables',
        'EXEFILE=./stella',
        'OUTFILE=$SLURM_JOBID"_"$SLURM_JOB_NAME".out"',
        'ERRFILE=$SLURM_JOBID"_"$SLURM_JOB_NAME".err"',
        ))
    
    information = '\n'.join((   
        '\n',
        'echo " "',
        'echo "#########################################################"',
        'echo "                   STELLA SIMULATION                     "',
        'echo "#########################################################"',
        'echo "RUN DATE:         "$(date)',
        'echo "COMPILATION DATE: '+compilation_data+'"',
        'echo "RUN DIRECTORY:    '+run_directory+'"',
        'echo "INPUT FILE:       '+input_file+'"',
        'echo "NODES:            '+str(nodes)+'"',
        'echo "CORES:            '+cores+'"',
        'echo "PATH:             "$PATH',
        'echo " "',
        ))
    
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
    

    execute_simulation = '\n'.join((   
        '\n',
        'echo "#########################################################"',
        'echo "                 EXECUTE SIMULATION                     "',
        'echo "#########################################################"',
        'mpirun -errfile-pattern $ERRFILE -outfile-pattern $OUTFILE'+\
        ' -envall -genv -n '+cores+' $EXEFILE '+input_file,
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
        function = "runSimulations.kyscan"
        default_arguments = {\
                'folder': None}
        arguments = get_arguments(function, default_arguments) 
        resolutionScan(**arguments) 

# Explain the program
def print_help(function, default_arguments):
    print_decoratorOnCommandPrompt("begin")
    print_functionInformationOnCommandPrompt(function, default_arguments, None)
    print("\n The arguments are:") 
    #print("      -i         <identifier>" )
    #print("      -c         Sets confirmation=True" )  
    print_decoratorOnCommandPrompt("end")

# Get the arguments from the terminal
def get_arguments(function, default_arguments):

    # Make a new dictionary which will hold the new arguments
    arguments = default_arguments.copy()

    # Initiate mandetory variables
    arguments['folder'] = pathlib.Path(os.getcwd())
    
    # Read the options and arguments from the terminal
    options = "h"
    long_options = ["folder=","help"]
    opts_args = get_bashArguments(options, long_options, print_help, function, default_arguments)

    # Asign the options and arguments to variables
    for opt, arg in opts_args:
        if opt in ("-h", "--help"):
            print_help(function, default_arguments)
            sys.exit()
        elif opt == '--folder':
            arguments['folder'] = str(arg)
            
    # Show the parameters
    print_functionInformationOnCommandPrompt(function, default_arguments, arguments)
    return arguments

# Execute the python script
main_program()
