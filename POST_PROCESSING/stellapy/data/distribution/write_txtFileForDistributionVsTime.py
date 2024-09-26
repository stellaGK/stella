
#!/usr/bin/python3 
import os, sys
import pathlib
import numpy as np 

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile 
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile 
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.utils.files.get_filesInFolder import get_filesInFolder  

#===============================================================================
#                  DISTRIBUTION SQUARED OF THE GUIDING CENTERS
#===============================================================================

def write_txtFileForDistributionVsTime(folder, dt=1, verbose=False):   
    
    # Time step
    dt = int(dt) if int(dt)==dt else dt   
    
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 
    
    # Iterate through the input files 
    for input_file in input_files: 
        
        # Check progress
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
    
        # Check whether we have a linear or nonlinear simulation   
        nonlinear = read_linearNonlinearFromInputFile(input_files[0])[1]
        
        # Depending on whether the simulation is linear or nonlinear write different data
        if not nonlinear: write_txtFileForDistributionVsTimeLinearSimulations(input_file, dt, status, verbose)
        if nonlinear: write_txtFileForDistributionVsTimeNonlinearSimulations(input_file, dt, status) 
        
    return 

#===============================================================================
#                            NONLINEAR SIMULATIONS                             #
#===============================================================================
 
def write_txtFileForDistributionVsTimeNonlinearSimulations(input_file, dt, status):  
    """ Append to the txt file if it already exist. """  
    
    try:
                
        # Path of the new file  
        distribution_path = input_file.with_suffix(".dt"+str(dt)+".g2_vs_t") 
        netcdf_path = input_file.with_suffix(".out.nc") 
        
        # If the file doesn't exist, check whether gvmus and gzvs were written
        if not os.path.isfile(distribution_path): 
            netcdf_data  = read_outputFile(netcdf_path)  
            if "gvmus" not in netcdf_data.variables.keys() or "gzvs" not in netcdf_data.variables.keys():
                return
        
        # Check whether the txt file is older than the simulation
        outputFileHasChanged = False
        if os.path.isfile(distribution_path):  
            if netcdf_path.stat().st_mtime > distribution_path.stat().st_mtime:
                outputFileHasChanged = True 
                
        # Notify that the file already existed   
        if os.path.isfile(distribution_path) and not outputFileHasChanged:
            print(status+"The g(t) file already exists:", distribution_path.parent.name+"/"+distribution_path.name)
            return
                
        # If the output file changed, then append to the txt file
        elif os.path.isfile(distribution_path) and outputFileHasChanged:
            
            # Check whether the output file contains extra time points
            netcdf_data  = read_outputFile(netcdf_path)  
            dim_species  = len(read_netcdfVariables('type_of_species', netcdf_data))
            vec_time     = read_netcdfVariables('vec_time', netcdf_data) 
            netcdf_data.close()  
            
            # Edit the time vector in the output file
            indices  = get_indicesAtFixedStep(vec_time, dt)
            vec_time = vec_time[indices]
            vec_time = [round(n, 8) for n in vec_time]
            tlast_outputfile = vec_time[-1] 
    
            # Read the distribution versus time from the txt file
            data = np.loadtxt(distribution_path,skiprows=1, dtype='float').reshape(-1, dim_species+1)
            g2_vs_ts_txtfile = data[:,1:]
            vec_time_txtfile = data[:,0]
            tlast_txtfile = vec_time_txtfile[-1]
            
            # If the output file contains extra time steps, append to the text file
            if tlast_outputfile > tlast_txtfile: 
                index_tlast = np.argwhere(tlast_txtfile < vec_time)[0][0]
                g2_vs_ts_outputFile, vec_time_outputFile = read_distributionVsTime(dt, input_file, netcdf_path)
                _, dim_species = np.shape(g2_vs_ts_outputFile)  
                vec_time = np.append(vec_time_txtfile, vec_time_outputFile[index_tlast:], axis=0)  
                g2_vs_ts = np.append(g2_vs_ts_txtfile, g2_vs_ts_outputFile[index_tlast:,:], axis=0)  
                
            # If the output file has the same time steps, touch the text file
            elif tlast_outputfile <= tlast_txtfile:
                print(status+"The g(t) file is up to date:", distribution_path.parent.name+"/"+distribution_path.name)  
                os.system("touch "+str(distribution_path))
                return
            
        # Create the data if the file doesn't exist yet
        elif not os.path.isfile(distribution_path):   
            
            # Read the distribution versus time from the output file
            g2_vs_ts, vec_time = read_distributionVsTime(dt, input_file, netcdf_path)
            _, dim_species = np.shape(g2_vs_ts) 
                        
        # Header for the text file
        header="       time                "
        for s in range(dim_species): header += "s="+str(s)+"                "   
        
        # Write the g(t) data to a text file 
        gt_file = open(distribution_path,'w') 
        gt_file.write(header+"\n") 
        for i in range(len(vec_time)):    
            np.savetxt(gt_file, [[vec_time[i], *g2_vs_ts[i,:]]], fmt='%16.8e', delimiter="   ") 
        gt_file.close() 
            
        # Notify that we finished creating the file
        if outputFileHasChanged: print(status+"   ---> The g(t) file is updated as " + distribution_path.parent.name+"/"+distribution_path.name)   
        if not outputFileHasChanged:  print(status+"   ---> The g(t) file is saved as " + distribution_path.parent.name+"/"+distribution_path.name)   

    except:
        print(status+"   ---> Something went wrong writing g(t) for " +  distribution_path.parent.name+"/"+distribution_path.name)   
    return  

#===============================================================================
#                             LINEAR SIMULATIONS                               #
#===============================================================================
 
def write_txtFileForDistributionVsTimeLinearSimulations(input_file, dt, status, verbose):  
    """ Always create a new txt file. """  
    
    try:
                
        # Path of the new file  
        distribution_path = input_file.with_suffix(".dt"+str(dt)+".g2_vs_t") 
        netcdf_path = input_file.with_suffix(".out.nc") 
        
        # If the file doesn't exist, check whether gvmus and gzvs were written
        if not os.path.isfile(distribution_path): 
            netcdf_data  = read_outputFile(netcdf_path)  
            if "gvmus" not in netcdf_data.variables.keys() or "gzvs" not in netcdf_data.variables.keys():
                return
        
        # Check whether the txt file is older than the simulation
        outputFileHasChanged = False
        if os.path.isfile(distribution_path):  
            if netcdf_path.stat().st_mtime > distribution_path.stat().st_mtime:
                outputFileHasChanged = True 
                
        # Notify that the file already existed   
        if os.path.isfile(distribution_path) and not outputFileHasChanged:
            if verbose: print(status+"The g(t) file already exists:", distribution_path.parent.name+"/"+distribution_path.name)
            return 
            
        # Read the distribution versus time from the output file
        g2_vs_ts, vec_time = read_distributionVsTime(dt, input_file, netcdf_path)
        _, dim_species = np.shape(g2_vs_ts)  
                   
        # Header for the text file
        header="       time                "
        for s in range(dim_species): header += "s="+str(s)+"                "   
        
        # Write the g(t) data to a text file 
        gt_file = open(distribution_path,'w') 
        gt_file.write(header+"\n") 
        for i in range(len(vec_time)):    
            np.savetxt(gt_file, [[vec_time[i], *g2_vs_ts[i,:]]], fmt='%16.8e', delimiter="   ") 
        gt_file.close() 
                
        # Notify that we finished creating the file
        if outputFileHasChanged: print(status+"   ---> The g(t) file is updated as " + distribution_path.parent.name+"/"+distribution_path.name)   
        if not outputFileHasChanged:  print(status+"   ---> The g(t) file is saved as " + distribution_path.parent.name+"/"+distribution_path.name)   

    except:
        print(status+"   ---> Something went wrong writing g(t) for " +  distribution_path.parent.name+"/"+distribution_path.name)   
    return  

#---------------------------------------------
def read_distributionVsTime(dt, input_file, netcdf_path):
    ''' If something goes wrong in here: rewrite geometry files.
     >> *.unique*; *.ini '''

    # Read the geometry data in the output file 
    vmec_filename = read_vmecFileNameFromInputFile(input_file)
    path = create_dummyPathObject(input_file, vmec_filename)
    geometry = read_outputFileForGeometry(path)  
    try: geometry["vpa_weights"] 
    except: 
        print(geometry['source'])
        print("vpa_weights could not be found in the geometry data (write_txtFileForDistributionVsTime)"); 
        print(input_file); sys.exit()
    vpa_weights = geometry["vpa_weights"]  
    dl_over_B = geometry["dl_over_B"]  
    
    # Get the distribution data from the output file
    netcdf_data = read_outputFile(netcdf_path)  
    vec_time = read_netcdfVariables('vec_time', netcdf_data) 
    indices = get_indicesAtFixedStep(vec_time, dt)
    g2_vs_tsvpaz = read_netcdfVariables('g2_vs_tsvpaz', netcdf_data, indices) 
    vec_time = vec_time[indices]
    netcdf_data.close()       
    
    # Integrate over (vpa,mu,z) 
    g2_vs_ts = np.sum(g2_vs_tsvpaz[:,:,:,:]*vpa_weights[np.newaxis,np.newaxis,:,np.newaxis]*dl_over_B[np.newaxis,np.newaxis,np.newaxis,:], axis=(2,3))
                
    # Return the data 
    return g2_vs_ts, vec_time 
 
