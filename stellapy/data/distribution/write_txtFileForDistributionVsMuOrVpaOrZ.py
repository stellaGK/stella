

import os
import pathlib
import numpy as np
from datetime import datetime, timedelta
from stellapy.utils.files.get_filesInFolder import get_filesInFolder  
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile 
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables

#===============================================================================
#                       DISTRIBUTION OF THE GUIDING CENTERS
#===============================================================================

def write_txtFileForDistributionVsMuOrVpaOrZ(folder):   
    
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 

    # Iterate through the input files 
    for input_file in input_files: 
        
        # Only write the distribution at the final time step for linear simulations  
        nonlinear = read_linearNonlinearFromInputFile(input_file)[1]
        if nonlinear: continue
        
        # Processing status
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
        
        # Path of the new file  
        file_pathz = input_file.with_suffix(".g_vs_z")  
        file_pathmu = input_file.with_suffix(".g_vs_mu")
        file_pathvpa = input_file.with_suffix(".g_vs_vpa")
        netcdf_file = input_file.with_suffix(".out.nc")        
        output_file = input_file.with_suffix(".out") if os.path.isfile(input_file.with_suffix(".out")) else input_file.with_suffix(".out.nc")
        
        # If the file doesn't exist, check whether gvmus and gzvs were written
        if not os.path.isfile(file_pathz): 
            netcdf_data  = read_outputFile(netcdf_file)  
            if "gvmus" not in netcdf_data.variables.keys() or "gzvs" not in netcdf_data.variables.keys():
                continue
            
        # Check whether the txt file is older than the simulation
        simulationWasntFinished = False
        if os.path.isfile(file_pathz): 
            time_txtFile = datetime.fromtimestamp(file_pathz.stat().st_mtime)
            time_simulation = datetime.fromtimestamp(output_file.stat().st_mtime) 
            if time_simulation>time_txtFile+timedelta(minutes=5):
                simulationWasntFinished = True 
        
        # Create text files for g(mu) and g(vpa) 
        if not os.path.isfile(file_pathz) or simulationWasntFinished:   
            
            # Read the geometry data in the output file
            vmec_filename = read_vmecFileNameFromInputFile(input_file)
            path = create_dummyPathObject(input_file, vmec_filename)
            geometry = read_outputFileForGeometry(path) 
            vpa_weights = geometry["vpa_weights"] 
            mu_weights = geometry["mu_weights"] 
            dl_over_B = geometry["dl_over_B"] 
            header="        "
            
            # Get the distribution data from the output file
            netcdf_data = read_outputFile(netcdf_file)  
            g_vs_tsvpaz = read_netcdfVariables('g_vs_tsvpaz', netcdf_data) 
            g_vs_tsmuvpa = read_netcdfVariables('g_vs_tsmuvpa', netcdf_data) 
            netcdf_data.close()
            
            # Get the (t, species, mu, vpa) dimensions  
            dim_time, dim_species, dim_mu, dim_vpa = np.shape(g_vs_tsmuvpa)
            dim_time, dim_species, dim_vpa, dim_z = np.shape(g_vs_tsvpaz) 
            
            # Get the last time step without nans 
            for it in range(dim_time-1, -1, -1):
                if not np.any(np.isnan(g_vs_tsmuvpa[it,:,:,:])):
                    g_vs_smuvpa = g_vs_tsmuvpa[it,:,:,:] 
                    break   
    
            # Calculate the final g(z), g(mu) and g(vpa)
            g_vs_sz = np.zeros((dim_species, dim_z)) 
            g_vs_smu = np.zeros((dim_species, dim_mu)) 
            g_vs_svpa = np.zeros((dim_species, dim_vpa)) 
            
            # Sum over (vpa,mu,z) in g_vs_smuvpa(s,mu,vpa) 
            for vpa in range(dim_vpa): 
                for mu in range(dim_mu):
                    for z in range(len(dl_over_B)): 
                        product = g_vs_smuvpa[:,mu,vpa]*mu_weights[z,mu]*vpa_weights[vpa]*dl_over_B[z] 
                        g_vs_sz[:,z] += product
                        g_vs_smu[:,mu] += product
                        g_vs_svpa[:,vpa] += product 
                        
            # Header for the files 
            for s in range(dim_species): header += "s="+str(s)+"                "  
        
            # Write the g(z) data to a text file 
            gz_file = open(file_pathz,'w') 
            gz_file.write(header+"\n")
            for i in range(dim_z):  
                np.savetxt(gz_file, [g_vs_sz[:,i].transpose()], fmt='%16.8e', delimiter="   ") 
            gz_file.close()
            
            # Write the g(mu) data to a text file 
            gmu_file = open(file_pathmu,'w') 
            gmu_file.write(header+"\n")
            for i in range(dim_mu):  
                np.savetxt(gmu_file, [g_vs_smu[:,i].transpose()], fmt='%16.8e', delimiter="   ") 
            gmu_file.close()
 
            # Write the g(vpa) data to a text file 
            gvpa_file = open(file_pathvpa,'w') 
            gvpa_file.write(header+"\n")
            for i in range(dim_vpa):
                np.savetxt(gvpa_file, [g_vs_svpa[:,i].transpose()], fmt='%16.8e', delimiter="   ") 
            gvpa_file.close()
                
            # Notify that we finished creating the file
            print(status+"    ---> The g(z), g(mu) and g(vpa) files are saved as " + file_pathz.parent.name+"/"+file_pathz.name)    
         
        # Notify that the file already existed   
        else:
            print(status+"The g(z), g(mu) and g(vpa) files already exist:", file_pathz.parent.name+"/"+file_pathz.name)
            
    return  
 
