
import os
import numpy as np   
from stellapy.utils.files.get_filesInFolder import get_filesInFolder
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile,\
    read_numberOfModesFromInputFile 

#===============================================================================
#                      REPLACE *.OMEGA BY *.DT1.OMEGA
#===============================================================================
# Make sure to run this bash command on marconi after every simulation.
#     > reduce_sizeOmega -t 1
#===============================================================================

def write_txtFileForOmegaVsTime(folder, dt=1):
    """ Reduce the size of the *.omega file by reducing the time axis. """
    
    # Time step
    dt = int(dt) if int(dt)==dt else dt   
    
    # Suffix of the new file
    suffix = ".dt"+str(dt)+".omega_vs_t"
    
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.omega'))]
    if input_files==[]: return 
    
    # Go through the input files  
    for input_file in input_files: 
        
        # Processing status
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
        
        # Check whether the reduced file has already been created
        already_written = False
        if os.path.isfile(input_file.with_suffix(suffix)): 
            time_h5File = input_file.with_suffix(suffix).stat().st_mtime
            time_netcdfFile = input_file.with_suffix(".omega").stat().st_mtime
            if time_h5File>time_netcdfFile:
                already_written = True
                
        # Check whether we have a linear or nonlinear simulation   
        nonlinear = read_linearNonlinearFromInputFile(input_file)[1]
        
        # Check how many modes are written 
        dim_kx, dim_ky = read_numberOfModesFromInputFile(input_file)
        multiple_modes_per_file = True if (dim_kx*dim_ky>1) else False

        # Only write it if we have to
        if already_written==False and nonlinear==False:

            # Read the omega file and return omega_data[kx,ky,time,{0:kx,1:ky,2:time,3:omega,4:gamma}]
            omega_data = read_omega(input_file) 
            
            # Find the indices at which to keep the data (time axis)
            indices = [0]; current_time = omega_data[0,0,0,2]; time_dim = len(omega_data[0,0,:,2])
            for it in range(time_dim):
                if omega_data[0,0,it,2] >= current_time+dt:
                    current_time = omega_data[0,0,it,2]
                    indices.append(it) 

            # Data and header if there is one mode in the file
            if not multiple_modes_per_file: 
                header = "      #time              omega            gamma  \n" 
                omega_data = omega_data[0,0,indices,2:]
            
            # Data and header if there are multiple modes in the file
            if multiple_modes_per_file: 
                header = "      #time               omega              gamma               kx                ky  \n"
                header = "      # kx                  ky               time               omega              gamma               \n"
                omega_data = omega_data[:,:,indices,:].reshape(-1, 5)

            # Create the new file and add the header
            omega_file = open(input_file.with_suffix(suffix),'w') 
            omega_file.write(header)
            np.savetxt(omega_file, omega_data, fmt='%16.8e', delimiter="   ") 
            omega_file.close()
            
            # Write a message
            print(status+"The new file is", str(input_file.with_suffix(suffix)))
            print(status+"   ---> Reduce the time dimension of omega(t) and gamma(t) from", np.shape(omega_data)[0], "to", len(indices), "points.")

def read_omega(input_file):
    ''' Read the *.omega file: [time,  ky,  kx,  Re(om),  Im(om),  Re(omavg), Im(omavg)]
    and assume we only have one mode (kx, ky) for this reduction!'''
    
    # Find the omega file corresponding to the input file
    omega_file = input_file.with_suffix(".omega")
    
    # Read the data in the omega file
    dim_kx, dim_ky = read_numberOfModesFromInputFile(input_file)
    omega_data = np.loadtxt(omega_file, dtype='float').reshape(-1, dim_kx, dim_ky, 7)[:,:,:,:]  
    
    # In the fourth axis we have [time, ky, kx, omega, gamma, omega_avg, gamma_avg], trow away the last two columns
    omega_data = omega_data[:,:,:,:5]
    
    # Move the data so that in the fourth columns we have [kx, ky, time, omega, gamma]
    omega_data = np.concatenate((omega_data[:,:,:,np.newaxis,2], omega_data[:,:,:,np.newaxis,1], omega_data[:,:,:,np.newaxis,0], omega_data[:,:,:,np.newaxis,3], omega_data[:,:,:,np.newaxis,4]), axis=3)
    omega_data = np.swapaxes(omega_data,0,2)
    omega_data = np.swapaxes(omega_data,0,1)

    # Return the data
    return omega_data

