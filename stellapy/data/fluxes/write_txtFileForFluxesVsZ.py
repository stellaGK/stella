
#!/usr/bin/python3  
import sys, os    
import pathlib
import numpy as np  

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)    
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile 
from stellapy.data.output.read_outputFile import read_netcdfVariables
from stellapy.utils.files.get_filesInFolder import get_filesInFolder  
from stellapy.data.output.read_outputFile import read_outputFile

#===============================================================================
#                               WRITE FLUXES(Z)
#=============================================================================== 

def write_txtFileForFluxesVsZ(folder, verbose=False):
    
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    input_files = [i for i in input_files if read_linearNonlinearFromInputFile(i)[0]==True]
    if input_files==[]: return 
    
    # Go through the input files
    for input_file in input_files:  
        
        try:
        
            # Processing status
            status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   " 
    
            # Path of the new file 
            fluxes_path = input_file.with_suffix(".fluxes_vs_z") 
            netcdf_path = input_file.with_suffix(".out.nc")   
            
            # If the file doesn't exist, check whether the fluxes were written to the netcdf file
            if not os.path.isfile(fluxes_path): 
                netcdf_data = read_outputFile(netcdf_path)   
                if "qflx_kxky" not in netcdf_data.variables.keys():
                    continue  
                
            # Check whether the fluxes file is older than the simulation
            outputFileHasChanged = False
            if os.path.isfile(fluxes_path):  
                if netcdf_path.stat().st_mtime > fluxes_path.stat().st_mtime:
                    outputFileHasChanged = True 
                    
            # Notify that the file already existed   
            if os.path.isfile(fluxes_path) and not outputFileHasChanged:
                if verbose: print(status+"The fluxes(z) file already exists:", fluxes_path.parent.name+"/"+fluxes_path.name)
                continue 
                    
            # Check whether we have <qflux_vs_tskxky> or <qflux_vs_tszkxky>
            netcdf_data = read_outputFile(netcdf_path) 
            qflux_vs_tskxky = read_netcdfVariables('qflux_vs_tskxky', netcdf_data)  
            if len(np.shape(qflux_vs_tskxky))!=6: z_dimension = False
            if len(np.shape(qflux_vs_tskxky))==6: z_dimension = True  
            dim_species = np.shape(qflux_vs_tskxky)[1]
            netcdf_data.close()
            if z_dimension==False: continue
    
            # Read the fluxes and the potential squared from the netcdf file  
            qflux_vs_szkxky, pflux_vs_szkxky, vflux_vs_szkxky, phi2_vs_zkxky = read_fluxesAndPhi2(netcdf_path, status) 
            
            # Header for the file 
            header = "       phi2            "
            for s in range(dim_species): header += "pflux (s="+str(s)+")        "    
            for s in range(dim_species): header += "vflux (s="+str(s)+")        "    
            for s in range(dim_species): header += "qflux (s="+str(s)+")        "    
            header += "\n"
            
            # Create the new file and add the header
            fluxes_file = open(fluxes_path,'w') 
            fluxes_file.write(header) 
            data = np.append(phi2_vs_zkxky[np.newaxis,:,0,0], pflux_vs_szkxky[:,:,0,0], axis=0)   
            data = np.append(data, vflux_vs_szkxky[:,:,0,0], axis=0)   
            data = np.append(data, qflux_vs_szkxky[:,:,0,0], axis=0)    
            np.savetxt(fluxes_file, np.transpose(data), fmt='%16.7e', delimiter="   ")
                    
            # Close the file
            fluxes_file.close()
            
            # Write a message
            print(status+"   ---> The fluxes(z) file is saved as " + fluxes_path.name)  
            
        except:
            print(status+"Something went wrong for:", fluxes_path.parent.parent.name+"/"+fluxes_path.parent.name+"/"+fluxes_path.name)
            sys.exit()
    return


#---------------------------------------------
def read_fluxesAndPhi2(netcdf_file, status): 
        
    # Read the fluxes from the output file at the final time step
    netcdf_data = read_outputFile(netcdf_file)   
    qflux_vs_szkxky = read_netcdfVariables('qflux_vs_tszkxky', netcdf_data, time_indices=[-1])[0,:,:,:,:]
    pflux_vs_szkxky = read_netcdfVariables('pflux_vs_tszkxky', netcdf_data, time_indices=[-1])[0,:,:,:,:]
    vflux_vs_szkxky = read_netcdfVariables('vflux_vs_tszkxky', netcdf_data, time_indices=[-1])[0,:,:,:,:]
    phi_vs_zkxkyri = read_netcdfVariables('phi_vs_tzkxkyri', netcdf_data, time_indices=[-1])[0,:,:,:,:]  
    netcdf_data.close()    
        
    # If the value is NaN, take the last time step without nans
    if np.any(~np.isfinite(qflux_vs_szkxky)) or np.any(~np.isfinite(pflux_vs_szkxky)) or np.any(~np.isfinite(vflux_vs_szkxky)) or np.any(~np.isfinite(phi_vs_zkxkyri)):
        print(status+"   ---> qflux(tlast) = NaN, finding the last non-NaN time step.")   
        netcdf_data = read_outputFile(netcdf_file)   
        qflux_vs_tszkxky = read_netcdfVariables('qflux_vs_tszkxky', netcdf_data)
        pflux_vs_tszkxky = read_netcdfVariables('pflux_vs_tszkxky', netcdf_data)
        vflux_vs_tszkxky = read_netcdfVariables('vflux_vs_tszkxky', netcdf_data)
        phi_vs_tzkxkyri = read_netcdfVariables('phi_vs_tzkxkyri', netcdf_data)
        netcdf_data.close() 
        dim_time = np.shape(qflux_vs_tszkxky)[0] 
        for it in range(dim_time-1, -1, -1): 
            if (not np.any(~np.isfinite(qflux_vs_tszkxky[it,:,:,:,:]))) and \
            (not np.any(~np.isfinite(pflux_vs_tszkxky[it,:,:,:,:]))) and \
            (not np.any(~np.isfinite(vflux_vs_tszkxky[it,:,:,:,:]))) and \
            (not np.any(~np.isfinite(phi_vs_tzkxkyri[it,:,:,:,:]))):
                qflux_vs_szkxky = qflux_vs_tszkxky[it,:,:,:,:] 
                pflux_vs_szkxky = pflux_vs_tszkxky[it,:,:,:,:] 
                vflux_vs_szkxky = vflux_vs_tszkxky[it,:,:,:,:] 
                phi_vs_zkxkyri = phi_vs_tzkxkyri[it,:,:,:,:] 
                break      
    
    # Calculate the potential squared
    phi_vs_zkxky = phi_vs_zkxkyri[:,:,:,0] + 1j*phi_vs_zkxkyri[:,:,:,1]
    phi2_vs_zkxky = np.abs(phi_vs_zkxky)**2   
    return qflux_vs_szkxky, pflux_vs_szkxky, vflux_vs_szkxky, phi2_vs_zkxky
