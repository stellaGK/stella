  
import numpy as np
import netCDF4 as nc4 
import os, h5py, pathlib 
from stellapy.utils.decorators.verbose import noverbose 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder   
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep  
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables 
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry  
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile, read_fullFluxSurfaceFromInputFile
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile  
from stellapy.data.paths.load_pathObject import create_dummyPathObject

################################################################################
#                       WRITE THE FLUXES TO AN H5 FILE
################################################################################

@noverbose
def write_h5FileForFluxes3D(folder, dt=10, automatic=False):  
    
    # Time step
    dt = int(dt) if int(dt)==dt else dt   
     
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 
 
    # Iterate through the input files  
    for input_file in input_files:
        
        try:
        
            # Processing status
            status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
            nonlinear = read_linearNonlinearFromInputFile(input_file)[1] 
            if (not nonlinear) and automatic: continue  
                
            # Path of the new file 
            fluxes_path = input_file.with_suffix(".dt"+str(dt)+".fluxes3D") 
            netcdf_path = input_file.with_suffix(".out.nc")  
            
            # If the file doesn't exist, check whether the fluxes were written to the netcdf file
            if not os.path.isfile(fluxes_path): 
                netcdf_data  = read_outputFile(netcdf_path)   
                if "qflx_kxky" not in netcdf_data.variables.keys():
                    continue  
                
            # Check whether the fluxes file is older than the simulation
            outputFileHasChanged = False
            if os.path.isfile(fluxes_path):  
                if netcdf_path.stat().st_mtime > fluxes_path.stat().st_mtime:
                    outputFileHasChanged = True 
                    
            # If the file does exist, check whether we have <qflux_vs_tskxky> or <qflux_vs_tszkxky>
            with nc4.Dataset(netcdf_path) as netcdf_data:  
                name = netcdf_data.variables["qflx_kxky"]
                shape = str(name).split("shape = ")[1].split('\n')[0].replace(" ", "")  
            shape = shape.replace("(","").replace(")","").split(",") 
            if len(shape)!=6: z_dimension = False
            if len(shape)==6: z_dimension = True
            if (not z_dimension) and (not nonlinear): continue 
            
            # Make sure to add this data to already existing files (REMOVE THIS IN FUTURE)
            if os.path.isfile(fluxes_path):
                with h5py.File(fluxes_path, 'r') as f:
                    if ('qflux_vs_tsz' not in f.keys()) and z_dimension:
                        outputFileHasChanged = True
                    
            # Notify that the file already existed   
            if os.path.isfile(fluxes_path) and not outputFileHasChanged:
                print(status+"The 3D fluxes file already exists:", fluxes_path.parent.name+"/"+fluxes_path.name)
                continue
                 
            # If the output file changed, then append to the h5 file
            elif os.path.isfile(fluxes_path) and outputFileHasChanged:
                
                # Check whether the output file (*.out.nc) contains extra time points
                netcdf_data  = read_outputFile(netcdf_path)   
                vec_time     = read_netcdfVariables('vec_time', netcdf_data) 
                netcdf_data.close()  
                
                # Edit the time vector in the output file
                indices  = get_indicesAtFixedStep(vec_time, dt)
                vec_time = vec_time[indices]
                vec_time = [round(n, 8) for n in vec_time]
                tlast_outputfile = vec_time[-1]  
        
                # Read the data in the h5 file, if the flux(t,s,z) data is missing
                # Either remove and rewrite the file if we can add it, otherwise add nans           
                with h5py.File(fluxes_path, 'r') as f: 
                    vec_time_h5file = f['vec_time'][()]
                    tlast_h5file = vec_time_h5file[-1] 
                    if 'qflux_vs_tsz' in f.keys():
                        qflux_vs_tsz_h5file = f['qflux_vs_tsz'][()]
                        pflux_vs_tsz_h5file = f['pflux_vs_tsz'][()]
                        vflux_vs_tsz_h5file = f['vflux_vs_tsz'][()] 
                    if nonlinear:
                        qflux_vs_tskx_h5file = f['qflux_vs_tskx'][()]
                        pflux_vs_tskx_h5file = f['pflux_vs_tskx'][()]
                        vflux_vs_tskx_h5file = f['vflux_vs_tskx'][()]
                        qflux_vs_tsky_h5file = f['qflux_vs_tsky'][()]
                        pflux_vs_tsky_h5file = f['pflux_vs_tsky'][()]
                        vflux_vs_tsky_h5file = f['vflux_vs_tsky'][()]
                    else: 
                        if z_dimension:
                            print(status+"Data was missing in the already existing fluxes 3D file: ")
                            print(status+"     "+fluxes_path.parent.parent.name+"/"+fluxes_path.parent.name+"/"+fluxes_path.name)
                            print(status+"          ---> Removed this file, PLEASE RERUN <write_dataFiles>.")
                            os.system("rm "+str(fluxes_path))
                            continue 
                
                # If the output file (*.out.nc) contains extra time steps, append to the h5 file 
                if tlast_outputfile > tlast_h5file: 
                    index_tlast = np.argwhere(tlast_h5file < vec_time)[0][0] 
                    dl_over_B = read_dlOverB(input_file)
                    vec_time, qflux_vs_tskx, qflux_vs_tsky, qflux_vs_tsz = read_fluxes3D("qflux_vs_tszkxky", dt, netcdf_path, dl_over_B, z_dimension) 
                    vec_time, pflux_vs_tskx, pflux_vs_tsky, pflux_vs_tsz = read_fluxes3D("pflux_vs_tszkxky", dt, netcdf_path, dl_over_B, z_dimension) 
                    vec_time, vflux_vs_tskx, vflux_vs_tsky, vflux_vs_tsz = read_fluxes3D("vflux_vs_tszkxky", dt, netcdf_path, dl_over_B, z_dimension) 
                    vec_time = np.append(vec_time_h5file, vec_time[index_tlast:], axis=0)  
                    qflux_vs_tskx = np.append(qflux_vs_tskx_h5file, qflux_vs_tskx[index_tlast:,:,:], axis=0)  
                    pflux_vs_tskx = np.append(pflux_vs_tskx_h5file, pflux_vs_tskx[index_tlast:,:,:], axis=0)  
                    vflux_vs_tskx = np.append(vflux_vs_tskx_h5file, vflux_vs_tskx[index_tlast:,:,:], axis=0) 
                    qflux_vs_tsky = np.append(qflux_vs_tsky_h5file, qflux_vs_tsky[index_tlast:,:,:], axis=0)  
                    pflux_vs_tsky = np.append(pflux_vs_tsky_h5file, pflux_vs_tsky[index_tlast:,:,:], axis=0)  
                    vflux_vs_tsky = np.append(vflux_vs_tsky_h5file, vflux_vs_tsky[index_tlast:,:,:], axis=0)  
                    if z_dimension:
                        qflux_vs_tsz = np.append(qflux_vs_tsz_h5file, qflux_vs_tsz[index_tlast:,:,:], axis=0)  
                        pflux_vs_tsz = np.append(pflux_vs_tsz_h5file, pflux_vs_tsz[index_tlast:,:,:], axis=0)  
                        vflux_vs_tsz = np.append(vflux_vs_tsz_h5file, vflux_vs_tsz[index_tlast:,:,:], axis=0) 
                    
                # If the output file has the same time steps, touch the text file
                elif tlast_outputfile <= tlast_h5file:
                    print(status+"The 3D fluxes file is up to date:", fluxes_path.parent.name+"/"+fluxes_path.name)  
                    os.system("touch "+str(fluxes_path))
                    continue
                
            # Create the text file for dphiz(t)
            elif not os.path.isfile(fluxes_path):    
                
                # Read the fluxes versus time from the output file
                dl_over_B = read_dlOverB(input_file)
                vec_time, qflux_vs_tskx, qflux_vs_tsky, qflux_vs_tsz = read_fluxes3D("qflux_vs_tszkxky", dt, netcdf_path, dl_over_B, z_dimension) 
                vec_time, pflux_vs_tskx, pflux_vs_tsky, pflux_vs_tsz = read_fluxes3D("pflux_vs_tszkxky", dt, netcdf_path, dl_over_B, z_dimension) 
                vec_time, vflux_vs_tskx, vflux_vs_tsky, vflux_vs_tsz = read_fluxes3D("vflux_vs_tszkxky", dt, netcdf_path, dl_over_B, z_dimension) 
                 
            # Create the new h5 file 
            with h5py.File(fluxes_path, 'w') as h5_file: 
                h5_file.create_dataset("vec_time", data=vec_time) 
                if nonlinear:
                    h5_file.create_dataset("qflux_vs_tskx", data=qflux_vs_tskx) 
                    h5_file.create_dataset("pflux_vs_tskx", data=pflux_vs_tskx) 
                    h5_file.create_dataset("vflux_vs_tskx", data=vflux_vs_tskx) 
                    h5_file.create_dataset("qflux_vs_tsky", data=qflux_vs_tsky) 
                    h5_file.create_dataset("pflux_vs_tsky", data=pflux_vs_tsky) 
                    h5_file.create_dataset("vflux_vs_tsky", data=vflux_vs_tsky)
                if z_dimension:
                    h5_file.create_dataset("qflux_vs_tsz", data=qflux_vs_tsz) 
                    h5_file.create_dataset("pflux_vs_tsz", data=pflux_vs_tsz) 
                    h5_file.create_dataset("vflux_vs_tsz", data=vflux_vs_tsz) 
                                     
            # Notify that we finished creating the file
            if outputFileHasChanged: print(status+"   ---> The 3D fluxes file is updated as " +  fluxes_path.parent.name+"/"+fluxes_path.name)   
            if not outputFileHasChanged:  print(status+"   ---> The 3D fluxes file is saved as " +  fluxes_path.parent.name+"/"+fluxes_path.name)  
        
        except:
            print(status+"Something went wrong for:", input_file.parent.parent.name+"/"+input_file.parent.name+"/"+input_file.name)
            import sys; sys.exit()
        
    return 




#---------------------------------------------
def read_dlOverB(input_file): 
    vmec_filename = read_vmecFileNameFromInputFile(input_file) 
    path = create_dummyPathObject(input_file, vmec_filename)
    geometry = read_outputFileForGeometry(path)  
    dl_over_B = geometry["dl_over_B"]  
    return dl_over_B

#---------------------------------------------
def read_fluxes3D(y_quantity, dt, netcdf_path, dl_over_B, z_dimension):
    """ qflux_vs_tszkxky, pflux_vs_tszkxky, pflux_vs_tszkxky """
    
    # If we have fluxes(t,s,z,kx,ky) then write fluxes(t,s,z); fluxes(t,s,kx) and fluxes(t,s,ky)
    if z_dimension:
        
        # Read the fluxes from the output file
        netcdf_data = read_outputFile(netcdf_path)  
        vec_time = read_netcdfVariables('vec_time', netcdf_data) 
        indices = get_indicesAtFixedStep(vec_time, dt)
        flux_vs_tszkxky = read_netcdfVariables(y_quantity, netcdf_data, indices) 
        vec_time = vec_time[indices]
        netcdf_data.close()   
        
        # Get flux(t,z) by summing over the (kx,ky) components
        flux_vs_tsz = np.sum(flux_vs_tszkxky[:,:,:,:,0] + 2*np.sum(flux_vs_tszkxky[:,:,:,:,1:],axis=4),axis=3)   
        
        # Remove the z-dimension 
        flux_vs_tskxky = np.sum(flux_vs_tszkxky[:,:,:,:,:]*dl_over_B[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis], axis=2)

        # Get flux(t,kx)  
        flux_vs_tskx = flux_vs_tskxky[:,:,:,0] + 2*np.sum(flux_vs_tskxky[:,:,:,1:],axis=3) 
        flux_vs_tsky = np.sum(flux_vs_tskxky[:,:,:,:],axis=2)
        return vec_time, flux_vs_tskx, flux_vs_tsky, flux_vs_tsz
        
    # If we only have fluxes(t,s,kx,ky) then write fluxes(t,s,kx) and fluxes(t,s,ky)
    if not z_dimension:
        
        # Read the fluxes from the output file
        netcdf_data = read_outputFile(netcdf_path)  
        vec_time = read_netcdfVariables('vec_time', netcdf_data) 
        indices = get_indicesAtFixedStep(vec_time, dt)
        flux_vs_tskxky = read_netcdfVariables(y_quantity, netcdf_data, indices) 
        netcdf_data.close()         
            
        # Get the data at every <dt> timestep 
        vec_time = vec_time[indices]
        
        # Get flux(t,kx)  
        flux_vs_tskx = flux_vs_tskxky[:,:,:,0] + 2*np.sum(flux_vs_tskxky[:,:,:,1:],axis=3) 
        flux_vs_tsky = np.sum(flux_vs_tskxky[:,:,:,:],axis=2)
        return vec_time, flux_vs_tskx, flux_vs_tsky, None
                     
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")   
    write_h5FileForFluxes3D(folder) 
    
