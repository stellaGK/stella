 
def read_cpuTimeFromSlurmFile(simulation):
    
    # Read the slurm file
    slurm_data  = open(simulation.input_file.with_suffix('.slurm'), 'r')
    slurm_text  = slurm_data.read()
    
    # Read the number of nodes
    nodes = int(slurm_text.split('#SBATCH -N')[-1].split('\n')[0].replace(' ',''))
    
    # Read the set wall time
    wall_time = [int(i) for i in slurm_text.split('#SBATCH --time')[-1].split('\n')[0].replace(' ','').split(':')]
    wall_time = wall_time[0] + wall_time[1]*1/60 + wall_time[2]*1/3600
    
    # CPU hours is wall_time*nodes*48cpus for Marconi
    cpu_time = wall_time*nodes*48 
    return cpu_time