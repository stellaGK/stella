"""============================================================================================================================================================================= """
""" --------------------------------------------------------------------- A python script for scanning nmu. -------------------------------------------------------------------- """
""" ============================================================================================================================================================================ """

import os
import f90nml
import subprocess
import shutil
from tqdm import tqdm 
import numpy as np 

# Define values to scan. 

nmu = [4, 8, 16, 24, 32, 40, 48, 56, 64, 72]

# Read the default CBC.in file. 
base_nml = f90nml.read('CBC.in')

# Define a function for copying a file from the source directory into each run directory.

def copy_file(filename): 
    src_file = filename # File to copy.
    dst_dir = run_dir   # Destination Directory. 
    shutil.copy(src_file, dst_dir) # Copy file into the directory. 
  
for i in range(len(nmu)):
    # Output directory for all runs.
    output_dir = f'NEO_stella_CBC_aky_2.0_beta_0.015_nmu_scan_no_neoclassics'
    os.makedirs(output_dir, exist_ok=True)
    run_name = f'nmu_{nmu[i]}'
    run_dir = os.path.join(output_dir, run_name)
    os.makedirs(run_dir, exist_ok=True)
    
    # copy_file(f'/users/rjs659/NEO_stella/large_grid_CBC_NEO_data/RHO_STAR_0.005455594781168515/out.neo.equil')  
    # copy_file(f'/users/rjs659/NEO_stella/large_grid_CBC_NEO_data/RHO_STAR_0.005455594781168515/out.neo.grid')
    # copy_file(f'/users/rjs659/NEO_stella/large_grid_CBC_NEO_data/RHO_STAR_0.005455594781168515/out.neo.species')
    # copy_file(f'/users/rjs659/NEO_stella/large_grid_CBC_NEO_data/RHO_STAR_0.005455594781168515/out.neo.version') # For metadata. 

    # copy_file(f'/users/rjs659/NEO_stella/large_grid_CBC_NEO_data/RHO_STAR_0.005455594781168515/out.neo.f.left') # Left flux surface.
    # copy_file(f'/users/rjs659/NEO_stella/large_grid_CBC_NEO_data/RHO_STAR_0.005455594781168515/out.neo.phi.left')
   
    # copy_file(f'/users/rjs659/NEO_stella/large_grid_CBC_NEO_data/RHO_STAR_0.005455594781168515/out.neo.f') # Central flux surface.
    # copy_file(f'/users/rjs659/NEO_stella/large_grid_CBC_NEO_data/RHO_STAR_0.005455594781168515/out.neo.phi')
   
    # copy_file(f'/users/rjs659/NEO_stella/large_grid_CBC_NEO_data/RHO_STAR_0.005455594781168515/out.neo.f.right') # Right flux surface.
    # copy_file(f'/users/rjs659/NEO_stella/large_grid_CBC_NEO_data/RHO_STAR_0.005455594781168515/out.neo.phi.right')

    # Modify the namelist.
    mod_nml = base_nml.copy()
    mod_nml['velocity_grids']['nmu'] = nmu[i]
    
    # Write new input.in
    input_path = os.path.join(run_dir, 'CBC.in')
    mod_nml.write(input_path)
    
    # Create a jobscript.
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name=neo_stella                          # Job name
#SBATCH --partition=nodes                              # What partition the job should run on
#SBATCH --time=0-01:00:00                              # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=96                                    # Number of MPI tasks to request
#SBATCH --cpus-per-task=1                              # Number of CPU cores per MPI task
#SBATCH --account=pet-gspt-2019                        # Project account to use
#SBATCH --mail-type=NONE                               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rjs659@york.ac.uk                  # Where to send mail
#SBATCH --output=%x-%j.log                             # Standard output log
#SBATCH --error=%x-%j.err                              # Standard error log

# Purge any previously loaded modules
module purge

# Load modules
module load gompi/2022b OpenMPI/4.1.4-GCC-12.2.0 netCDF-Fortran/4.6.0-gompi-2022b FFTW/3.3.10-GCC-12.2.0 OpenBLAS/0.3.21-GCC-12.2.0 Python/3.10.8-GCCcore-12.2.0

export PSM2_CUDA=0

export GK_SYSTEM='viking'
export MAKEFLAGS='-IMakefiles'
export HDF5_USE_FILE_LOCKING=FALSE

ulimit -s unlimited

######################## ABove is in bashrc anyway

export OMP_NUM_THREADS=1

EXE=/users/rjs659/NEO_stella/stella
INPUT=CBC.in
OUTPUT=OUTPUT

srun --hint=nomultithread --distribution=block:block -n 96 $EXE $INPUT | tee $OUTPUT
"""

    script_path = os.path.join(run_dir, 'jobscript.job')
    with open(script_path, 'w') as f:
        f.write(slurm_script)

    # Submit the job.
    os.chdir(run_dir)
    os.system("sbatch jobscript.job")
    os.chdir("../../")
