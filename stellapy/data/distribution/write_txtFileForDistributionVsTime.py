

import os
import numpy as np 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder  
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep

################################################################################
#                       DISTRIBUTION OF THE GUIDING CENTERS
################################################################################

def write_txtFileForDistributionVsTime(folder, dt=1):   
    
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
        
        # Path of the new file  
        txt_file = input_file.with_suffix(".dt"+str(dt)+".g_vs_t") 
        netcdf_file = input_file.with_suffix(".out.nc") 
        
        # If the file doesn't exist, check whether gvmus and gzvs were written
        if not os.path.isfile(txt_file): 
            netcdf_data  = read_outputFile(netcdf_file)  
            if "gvmus" not in netcdf_data.variables.keys() or "gzvs" not in netcdf_data.variables.keys():
                continue
        
        # Check whether the txt file is older than the simulation
        outputFileHasChanged = False
        if os.path.isfile(txt_file):  
            if netcdf_file.stat().st_mtime > txt_file.stat().st_mtime:
                outputFileHasChanged = True 
                
        # Notify that the file already existed   
        if os.path.isfile(txt_file) and not outputFileHasChanged:
            print(status+"The g(t) file already exists:", txt_file.parent.name+"/"+txt_file.name)
            continue
                
        # If the output file changed, then append to the txt file
        elif os.path.isfile(txt_file) and outputFileHasChanged:
            
            # Check whether the output file contains extra time points
            netcdf_data  = read_outputFile(netcdf_file)  
            dim_species  = len(read_netcdfVariables('type_of_species', netcdf_data))
            vec_time     = read_netcdfVariables('vec_time', netcdf_data) 
            netcdf_data.close()  
            
            # Edit the time vector in the output file
            indices  = get_indicesAtFixedStep(vec_time, dt)
            vec_time = vec_time[indices]
            vec_time = [round(n, 8) for n in vec_time]
            tlast_outputfile = vec_time[-1] 
    
            # Read the distribution versus time from the txt file
            data = np.loadtxt(txt_file,skiprows=1, dtype='float').reshape(-1, dim_species+1)
            g_vs_ts_txtfile = data[:,1:]
            vec_time_txtfile = data[:,0]
            tlast_txtfile = vec_time_txtfile[-1]
            
            # If the output file contains extra time steps, append to the text file
            if tlast_outputfile > tlast_txtfile: 
                index_tlast = np.argwhere(tlast_txtfile < vec_time)[0][0]
                g_vs_ts_outputFile, vec_time_outputFile = read_distributionVsTime(dt, input_file, netcdf_file)
                dim_time, dim_species = np.shape(g_vs_ts_outputFile)  
                vec_time = np.append(vec_time_txtfile, vec_time_outputFile[index_tlast:], axis=0)  
                g_vs_ts = np.append(g_vs_ts_txtfile, g_vs_ts_outputFile[index_tlast:,:], axis=0)  
                
            # If the output file has the same time steps, touch the text file
            elif tlast_outputfile <= tlast_txtfile:
                print(status+"The g(t) file is up to date:", txt_file.parent.name+"/"+txt_file.name)  
                os.system("touch "+str(txt_file))
                continue
            
        # Create the text file for g(t)
        elif not os.path.isfile(txt_file):   
            
            # Read the distribution versus time from the output file
            g_vs_ts, vec_time = read_distributionVsTime(dt, input_file, netcdf_file)
            dim_time, dim_species = np.shape(g_vs_ts) 
                        
        # Header for the text file
        header="       time                "
        for s in range(dim_species): header += "s="+str(s)+"                "   
        
        # Write the g(t) data to a text file 
        gt_file = open(txt_file,'w') 
        gt_file.write(header+"\n") 
        for i in range(len(vec_time)):    
            np.savetxt(gt_file, [[vec_time[i], *g_vs_ts[i,:]]], fmt='%16.8e', delimiter="   ") 
        gt_file.close() 
            
        # Notify that we finished creating the file
        if outputFileHasChanged: print(status+"   ---> The g(t) file is updated as " + txt_file.parent.name+"/"+txt_file.name)   
        if not outputFileHasChanged:  print(status+"   ---> The g(t) file is saved as " + txt_file.parent.name+"/"+txt_file.name)   
    return  

#---------------------------------------------
def read_distributionVsTime(dt, input_file, netcdf_file):
    ''' If something goes wrong in here: rewrite geometry files.
     >> *.unique*; *.ini '''

    # Read the geometry data in the output file 
    vmec_filename = read_vmecFileNameFromInputFile(input_file)
    nonlinear = read_linearNonlinearFromInputFile(input_file)[1]
    path = create_dummyPathObject(input_file, vmec_filename, nonlinear)
    geometry = read_outputFileForGeometry(path)  
    try: geometry["vpa_weights"] 
    except: 
        print(geometry['source'])
        print("vpa_weights could not be found in the geometry data (write_txtFileForDistributionVsTime)"); 
        print(input_file); import sys; sys.exit()
    vpa_weights = geometry["vpa_weights"] 
    mu_weights = geometry["mu_weights"] 
    dl_over_B = geometry["dl_over_B"]  
    
    # Get the distribution data from the output file
    netcdf_data  = read_outputFile(netcdf_file)  
    g_vs_tsmuvpa = read_netcdfVariables('g_vs_tsmuvpa', netcdf_data) 
    vec_time     = read_netcdfVariables('vec_time', netcdf_data) 
    netcdf_data.close() 
    
    # Check the dimension
    if np.shape(mu_weights)[1] != np.shape(g_vs_tsmuvpa)[2]: 
        print('\n          Something is wrong with the unique geometry files!', np.shape(mu_weights)[0], np.shape(g_vs_tsmuvpa)[2])
        print('          Input file:', input_file)
        print('          Geometry file:', path.geometry)
        print('          Geometry source:', geometry["source"]) 
        print('          dl_over_B:', np.shape(dl_over_B))
        print('          mu_weights:', np.shape(mu_weights))
        print('          vpa_weights:', np.shape(vpa_weights))
        print('          g_vs_tsmuvpa:', np.shape(g_vs_tsmuvpa))
        print('          Please run >> rm *unique*; rm *.ini; write_dataFiles\n')
        import sys; sys.exit()
    
    # Get the data at every <dt> timestep
    indices = get_indicesAtFixedStep(vec_time, dt)
    g_vs_tsmuvpa = g_vs_tsmuvpa[indices] 
    vec_time = vec_time[indices]
    
    # Get the (t, species, mu, vpa) dimensions  
    dim_time, dim_species, dim_mu, dim_vpa = np.shape(g_vs_tsmuvpa)  
    
    # Calculate the final g(vpa) and g(mu)
    g_vs_ts = np.zeros((dim_time, dim_species))  
    
    # Sum over (vpa,mu,z) in g_vs_smuvpa(s,mu,vpa) 
    for vpa in range(dim_vpa): 
        for mu in range(dim_mu):
            for z in range(len(dl_over_B)): 
                product = g_vs_tsmuvpa[:,:,mu,vpa]*mu_weights[z,mu]*vpa_weights[vpa]*dl_over_B[z] 
                g_vs_ts[:,:] += product 
                
    # Return the data 
    return g_vs_ts, vec_time 

################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    import pathlib
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI/test")
    write_txtFileForDistributionVsTime(folder)
    
