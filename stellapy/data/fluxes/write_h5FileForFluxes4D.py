  
import numpy as np
import os, h5py, pathlib 
from stellapy.utils.decorators.verbose import noverbose 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder   
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep  
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables 
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry  
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile  
from stellapy.data.paths.load_pathObject import create_dummyPathObject

################################################################################
#                       WRITE THE FLUXES TO AN H5 FILE
################################################################################

@noverbose
def write_h5FileForFluxes4D(folder, dt=10):  
    
    # Time step
    dt = int(dt) if int(dt)==dt else dt   
     
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 
 
    # Iterate through the input files  
    for input_file in input_files:
        
        # Processing status
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
        
        # Only write the fluxes(t,kx,ky) file for nonlinear simulations 
        if read_linearNonlinearFromInputFile(input_file)[1]:    
            
            # Path of the new file 
            fluxes_file = input_file.with_suffix(".dt"+str(dt)+".fluxes4D") 
            netcdf_file = input_file.with_suffix(".out.nc") if not os.path.isfile(input_file.with_suffix(".out.nc.fluxes")) else input_file.with_suffix(".out.nc.fluxes")  
            
            # If the file doesn't exist, check whether the fluxes were written to the netcdf file
            if not os.path.isfile(fluxes_file): 
                netcdf_data  = read_outputFile(netcdf_file)   
                if "qflx_kxky" not in netcdf_data.variables.keys():
                    continue  
                
            # Check whether the fluxes file is older than the simulation
            outputFileHasChanged = False
            if os.path.isfile(fluxes_file):  
                if netcdf_file.stat().st_mtime > fluxes_file.stat().st_mtime:
                    outputFileHasChanged = True 
                    
            # If the file does exist, check whether we have <qflux_vs_tskxky> or <qflux_vs_tszkxky>
            netcdf_data  = read_outputFile(netcdf_file) 
            qflux_vs_tskxky = read_netcdfVariables('qflux_vs_tskxky', netcdf_data)  
            if len(np.shape(qflux_vs_tskxky))!=6: z_dimension = False
            if len(np.shape(qflux_vs_tskxky))==6: z_dimension = True
            
            # Make sure to add this data to already existing files (REMOVE THIS IN FUTURE)
            if os.path.isfile(fluxes_file):
                with h5py.File(fluxes_file, 'r') as f:
                    if ('qflux_vs_tsz' not in f.keys()) and z_dimension:
                        outputFileHasChanged = True
                    
            # Notify that the file already existed   
            if os.path.isfile(fluxes_file) and not outputFileHasChanged:
                print(status+"The 4D fluxes file already exists:", fluxes_file.parent.name+"/"+fluxes_file.name)
                continue
                 
            # If the output file changed, then append to the h5 file
            elif os.path.isfile(fluxes_file) and outputFileHasChanged:
                
                # Check whether the output file (*.out.nc) contains extra time points
                netcdf_data  = read_outputFile(netcdf_file)   
                vec_time     = read_netcdfVariables('vec_time', netcdf_data) 
                netcdf_data.close()  
                
                # Edit the time vector in the output file
                indices  = get_indicesAtFixedStep(vec_time, dt)
                vec_time = vec_time[indices]
                vec_time = [round(n, 8) for n in vec_time]
                tlast_outputfile = vec_time[-1]  
        
                # Read the data in the h5 file, if the flux(t,s,z) data is missing
                # Either remove and rewrite the file if we can add it, otherwise add nans           
                with h5py.File(fluxes_file, 'r') as f: 
                    vec_time_h5file = f['vec_time'][()]
                    qflux_vs_tskxky_h5file = f['qflux_vs_tskxky'][()]
                    pflux_vs_tskxky_h5file = f['pflux_vs_tskxky'][()]
                    vflux_vs_tskxky_h5file = f['vflux_vs_tskxky'][()] 
                    tlast_h5file = vec_time_h5file[-1]  
                
                # If the output file (*.out.nc) contains extra time steps, append to the h5 file 
                if tlast_outputfile > tlast_h5file: 
                    index_tlast = np.argwhere(tlast_h5file < vec_time)[0][0] 
                    vec_time, qflux_vs_tskxky, pflux_vs_tskxky, vflux_vs_tskxky = read_fluxes4D(dt, netcdf_file, input_file, z_dimension) 
                    vec_time = np.append(vec_time_h5file, vec_time[index_tlast:], axis=0)  
                    qflux_vs_tskxky = np.append(qflux_vs_tskxky_h5file, qflux_vs_tskxky[index_tlast:,:,:,:], axis=0)  
                    pflux_vs_tskxky = np.append(pflux_vs_tskxky_h5file, pflux_vs_tskxky[index_tlast:,:,:,:], axis=0)  
                    vflux_vs_tskxky = np.append(vflux_vs_tskxky_h5file, vflux_vs_tskxky[index_tlast:,:,:,:], axis=0) 

                # If the output file has the same time steps, touch the text file
                elif tlast_outputfile <= tlast_h5file:
                    print(status+"The 4D fluxes file is up to date:", fluxes_file.parent.name+"/"+fluxes_file.name)  
                    os.system("touch "+str(fluxes_file))
                    continue
                
            # Create the text file for dphiz(t)
            elif not os.path.isfile(fluxes_file):    
                
                # Read the fluxes versus time from the output file
                vec_time, qflux_vs_tskxky, pflux_vs_tskxky, vflux_vs_tskxky = read_fluxes4D(dt, netcdf_file, input_file, z_dimension) 
                         
            # Create the new h5 file 
            with h5py.File(fluxes_file, 'w') as h5_file: 
                h5_file.create_dataset("vec_time", data=vec_time) 
                h5_file.create_dataset("qflux_vs_tskxky", data=qflux_vs_tskxky) 
                h5_file.create_dataset("pflux_vs_tskxky", data=pflux_vs_tskxky) 
                h5_file.create_dataset("vflux_vs_tskxky", data=vflux_vs_tskxky)  
                                     
            # Notify that we finished creating the file
            if outputFileHasChanged: print(status+"   ---> The 4D fluxes file is updated as " +  fluxes_file.parent.name+"/"+fluxes_file.name)   
            if not outputFileHasChanged:  print(status+"   ---> The 4D fluxes file is saved as " +  fluxes_file.parent.name+"/"+fluxes_file.name)  
        
    return 



#---------------------------------------------
def read_fluxes4D(dt, netcdf_file, input_file, z_dimension):
    
    # If we have fluxes(t,s,z,kx,ky) then write fluxes(t,s,z); fluxes(t,s,kx) and fluxes(t,s,ky)
    if z_dimension:
        
        # Read the fluxes from the output file
        netcdf_data = read_outputFile(netcdf_file)  
        vec_time = read_netcdfVariables('vec_time', netcdf_data) 
        qflux_vs_tszkxky = read_netcdfVariables('qflux_vs_tszkxky', netcdf_data) 
        pflux_vs_tszkxky = read_netcdfVariables('pflux_vs_tszkxky', netcdf_data) 
        vflux_vs_tszkxky = read_netcdfVariables('vflux_vs_tszkxky', netcdf_data)  
        netcdf_data.close()  
            
        # Get the data at every <dt> timestep 
        indices = get_indicesAtFixedStep(vec_time, dt)
        qflux_vs_tszkxky = qflux_vs_tszkxky[indices,:,:,:,:] 
        pflux_vs_tszkxky = pflux_vs_tszkxky[indices,:,:,:,:]
        vflux_vs_tszkxky = vflux_vs_tszkxky[indices,:,:,:,:]
        vec_time = vec_time[indices]
        
        # Read the geometry data in the output file
        vmec_filename = read_vmecFileNameFromInputFile(input_file)
        path = create_dummyPathObject(input_file, vmec_filename)
        geometry = read_outputFileForGeometry(path)  
        dl_over_B = geometry["dl_over_B"]  
        
        # Remove the z-dimension 
        qflux_vs_tskxky = np.sum(qflux_vs_tszkxky[:,:,:,:,:]*dl_over_B[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis], axis=2)
        pflux_vs_tskxky = np.sum(pflux_vs_tszkxky[:,:,:,:,:]*dl_over_B[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis], axis=2)
        vflux_vs_tskxky = np.sum(vflux_vs_tszkxky[:,:,:,:,:]*dl_over_B[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis], axis=2)
        return vec_time, qflux_vs_tskxky, pflux_vs_tskxky, vflux_vs_tskxky
        
    # If we only have fluxes(t,s,kx,ky) then write fluxes(t,s,kx) and fluxes(t,s,ky)
    if not z_dimension:
        
        # Read the fluxes from the output file
        netcdf_data = read_outputFile(netcdf_file)  
        vec_time = read_netcdfVariables('vec_time', netcdf_data) 
        qflux_vs_tskxky = read_netcdfVariables('qflux_vs_tskxky', netcdf_data) 
        pflux_vs_tskxky = read_netcdfVariables('pflux_vs_tskxky', netcdf_data) 
        vflux_vs_tskxky = read_netcdfVariables('vflux_vs_tskxky', netcdf_data)  
        netcdf_data.close()         
            
        # Get the data at every <dt> timestep 
        indices = get_indicesAtFixedStep(vec_time, dt)
        qflux_vs_tskxky = qflux_vs_tskxky[indices,:,:,:] 
        pflux_vs_tskxky = pflux_vs_tskxky[indices,:,:,:]
        vflux_vs_tskxky = vflux_vs_tskxky[indices,:,:,:]
        vec_time = vec_time[indices]
        return vec_time, qflux_vs_tskxky, pflux_vs_tskxky, vflux_vs_tskxky
             
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")   
    write_h5FileForFluxes4D(folder) 
    
