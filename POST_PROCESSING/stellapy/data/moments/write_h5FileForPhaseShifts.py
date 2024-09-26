  
# WARNING: The phase shift file can be incorrect when a netcdf has been reduced and then
# continued, since it won't overwrite the old file if the timestep is bigger. 
# Implement checks when reading the phaseshift file to make sure vec_time is correct!

import copy
import os, h5py
import numpy as np
import time as timer   
import netCDF4 as nc4  
from scipy.io import netcdf as scnetcdf  
from datetime import datetime, timedelta 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder     
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry  
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables 
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile, read_nspecFromInputFile  
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.stella.load_stellaVariables import stella_variables
  
################################################################################
#                       WRITE THE DIMENSIONS TO AN H5 FILE
################################################################################
  
def write_h5FileForPhaseShifts(folder, dt=50, skip=1):    
       
    # Time step
    dt = int(dt) if int(dt)==dt else dt     
       
    # Get the input files for which a netcdf file exists
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 
   
    # Iterate through the input files  
    for input_file in input_files:  
              
        # Only write the moments file for nonlinear simulations 
        if read_linearNonlinearFromInputFile(input_file)[1]:
              
            # Start timer
            start = timer.time()
          
            # Processing status
            status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
              
            # Path of the new file
            netcdf_path = input_file.with_suffix(".out.nc")
            phaseshifts_path = input_file.with_suffix(".dt"+str(dt)+".phaseshifts")   
            nspec = read_nspecFromInputFile(input_file)
                      
            # Check whether the phase shifts file is older than the simulation
            outputFileHasChanged = False
            if os.path.isfile(phaseshifts_path):  
                if datetime.fromtimestamp(netcdf_path.stat().st_mtime) > (datetime.fromtimestamp(phaseshifts_path.stat().st_mtime)+timedelta(minutes=5)):
                    outputFileHasChanged = True
                if not outputFileHasChanged:
                    with h5py.File(phaseshifts_path, 'r') as f: 
                        if nspec==1 and ("Finished" not in f.keys() or "n0_and_T0_vs_tangleky" not in f.keys()):
                                outputFileHasChanged = True
                        if nspec==2 and ("Finished" not in f.keys() or "T0_and_T1_vs_tangleky" not in f.keys()):
                                outputFileHasChanged = True
                        if nspec==3 and ("Finished" not in f.keys() or "phi_and_T2_vs_tangleky" not in f.keys()):
                                outputFileHasChanged = True
                        if "T0_and_T1_vs_tangleky" in f.keys():
                            if np.isnan(f["phi_and_n0_vs_tanglekx"][-1,0,0]):
                                print("Delete the existing "+phaseshifts_path.parent.name+"/"+phaseshifts_path.name+" file since the kx-data is broken.")
                                os.system("rm "+str(phaseshifts_path))
                                outputFileHasChanged = True
                        variables = ["phi_and_n0", "phi_and_T0", "phi_and_n1", "phi_and_T1", "n0_and_T0", "n1_and_T1", "n0_and_n1", "T0_and_T1"]
                        for variable in variables:
                            if variable+"_skip" in f.keys(): 
                                if f[variable+"_skip"][()]>skip: outputFileHasChanged = True
                      
            # Notify that the file already existed   
            if os.path.isfile(phaseshifts_path) and not outputFileHasChanged:
                print(status+"The phase shifts file already exists:", phaseshifts_path.parent.name+"/"+phaseshifts_path.name)
                continue
                  
            # Otherwise read the 3D moments data
            elif not os.path.isfile(phaseshifts_path) or outputFileHasChanged: 
                  
                # Read the time axis and variables in the netcdf file
                netcdf_data = read_outputFile(netcdf_path)  
                vec_time = read_netcdfVariables('vec_time', netcdf_data)   
                vec_kx = read_netcdfVariables('vec_kx', netcdf_data)   
                vec_ky = read_netcdfVariables('vec_ky', netcdf_data)   
                variables = list(netcdf_data.variables.keys())
                netcdf_data.close()   
                  
                # Check whether the moments are written
                if "density" not in variables:
                    print(status+"The moments were not written to the netcdf file.")
                    continue    
                  
                # Get the integration weigths and the time indices
                dl_over_B = get_integrationWeightsAlongZ(input_file) 
      
                # Get the data at every <dt> timestep 
                time_indices = get_indicesAtFixedStep(vec_time, dt)
                time_old = copy.deepcopy(vec_time)
                time_new = vec_time[time_indices]  
                
                # Get the indicies of each time piece
                tstarts, tends, itstarts, itends = get_time_indices(time_new, time_old)
                
                # If the file already exist, check whether the time vector matches, otherwise create new file
                create_new_file = False; add_original_time_axis = False
                if os.path.isfile(phaseshifts_path):  
                    with h5py.File(phaseshifts_path, 'r') as f: 
                        tends_old = f["tends"][()]
                        if "vec_time_original" not in f:
                            add_original_time_axis = True
                        if "vec_time_original" in f:
                            vec_time_original = f["vec_time_original"][()]
                    if add_original_time_axis:
                        with h5py.File(phaseshifts_path, 'r+') as h5_file:    
                            potential_file = input_file.with_suffix(".dt1.phi_vs_t")
                            if os.path.isfile(potential_file):  
                                data = np.loadtxt(potential_file,skiprows=1,dtype='float').reshape(-1, 6)
                                vec_time_original = data[:,0]
                                h5_file.create_dataset("vec_time_original", data=vec_time_original)  
                    dt_old_end = vec_time_original[-2] - vec_time_original[-3]
                    dt_old_start = vec_time_original[1] - vec_time_original[0]
                    dt_end = time_old[-2] - time_old[-3]
                    dt_start = time_old[1] - time_old[0]
                    if np.any(np.round(tends_old,2)!=np.round(tends,2)):
                        if dt_old_end<dt_end or dt_old_start<dt_start:
                            print("\n     The netcdf file has a bigger time step (dt="+"{:.2f}".format(dt_start)+") then the phaseshifts file (dt="+"{:.2f}".format(dt_old_start)+").")
                            print("     The netcdf file has probably been reduced in the past to a bigger dt.")
                            print("     Please backup the file or remove the file and run the script again.\n")
                            os.system("touch "+str(phaseshifts_path))
                            continue
                        else: 
                            print("Delete the existing "+phaseshifts_path.parent.name+"/"+phaseshifts_path.name+" file since the time vector does not match.")
                            os.system("rm "+str(phaseshifts_path))
                            create_new_file = True
                    if "T0_and_T1_vs_tangleky" in f.keys():
                        print(f["phi_and_n0_vs_tanglekx"][-1,0,0])
                        if np.isnan(f["phi_and_n0_vs_tanglekx"][-1,0,0]):
                            print("Delete the existing "+phaseshifts_path.parent.name+"/"+phaseshifts_path.name+" file since the kx-data is broken.")
                            os.system("rm "+str(phaseshifts_path))
                            create_new_file = True
                if (not os.path.isfile(phaseshifts_path)) or create_new_file:  
                    with h5py.File(phaseshifts_path, 'a') as h5_file:   
                        h5_file.create_dataset("tends", data=np.array(tends)) 
                        h5_file.create_dataset("tstarts", data=np.array(tstarts)) 
                        h5_file.create_dataset("vec_time", data=np.array(time_new)) 
                        h5_file.create_dataset("vec_time_original", data=np.array(time_old)) 
                    
                # Collect the data needed to calculate the phase shifts
                args = {"vec_kx" : vec_kx, "vec_ky" : vec_ky, "itstarts" : itstarts, "itends" : itends, "dl_over_B" : dl_over_B}      
                args = {"skip" : skip, "nspec" : nspec, "netcdf_path" : netcdf_path, "phaseshifts_path" : phaseshifts_path, **args}
                
                # Calculate the phase shifts vs kx and ky 
                calculate_phaseshifts("phi", "n0", **args)
                calculate_phaseshifts("phi", "T0", **args)
                calculate_phaseshifts("phi", "n1", **args)
                calculate_phaseshifts("phi", "T1", **args)
                calculate_phaseshifts("n0", "T0", **args)
                calculate_phaseshifts("n1", "T1", **args)
                calculate_phaseshifts("n0", "n1", **args)
                calculate_phaseshifts("T0", "T1", **args) 
                calculate_phaseshifts("phi", "n2", **args)
                calculate_phaseshifts("phi", "T2", **args)
                
                # Remember that we finished the script
                with h5py.File(phaseshifts_path, 'r+') as h5_file:   
                    if "Finished" not in h5_file:
                        h5_file.create_dataset("Finished", data=1)  
                      
            # Notify that we finished creating the file
            end = timer.time(); elapsed_time = str(timedelta(seconds=round(end-start,0))) 
            print(status+"---> The phase shifts file is saved as " +  phaseshifts_path.parent.name+"/"+phaseshifts_path.name+" ("+elapsed_time+")")   
   
    return 

#--------------------------------------
def get_time_indices(time_new, time_old):
    
    # Initiate
    tstarts = []; tends = []; itstarts = []; itends = [] 
    
    # Iterate over the new time axis
    for itnew in range(len(time_new)):
        if itnew==len(time_new)-1: continue
         
        # Get the piece of time
        tstart = time_new[itnew]; itstarts.append(np.where(time_old==tstart)[0][0])
        tend = time_new[itnew+1]; itends.append(np.where(time_old==tend)[0][0])
        tstarts.append(tstart); tends.append(tend)
        
    return tstarts, tends, itstarts, itends

#--------------------------------------
def calculate_phaseshifts(y1_variable, y2_variable, nspec=2, skip=1, netcdf_path=None, phaseshifts_path=None,
    vec_kx=None, vec_ky=None, itstarts=None, itends=None, dl_over_B=None):
    
    # Check whether we can write the data
    if (nspec==1) and (("1" in y1_variable) or ("1" in y2_variable)): return 
    if (nspec==1) and (("2" in y1_variable) or ("2" in y2_variable)): return 
    if (nspec==2) and (("2" in y1_variable) or ("2" in y2_variable)): return 
    
    # Check whether this variable is already written
    variable_name = y1_variable + "_and_" + y2_variable
    with h5py.File(phaseshifts_path, 'r+') as f: 
        if (variable_name+"_weights") in f.keys() and (variable_name+"_skip") in f.keys():
            if skip >= f[variable_name+"_skip"][()]: 
                print("              ---> Phase shift between "+y1_variable+" and "+y2_variable+" is already written with skip = "+str(f[variable_name+"_skip"][()])+".")
                return
            else:
                print("              ---> Calculate phase shift between "+y1_variable+" and "+y2_variable+" with smaller skip = "+str(skip)+".")
        elif (variable_name+"_weights") in f.keys() and (variable_name+"_skip") not in f.keys():   
            f.create_dataset(variable_name+"_skip", data=1)  
            print("              ---> Phase shift between "+y1_variable+" and "+y2_variable+" is already written with skip = "+str(f[variable_name+"_skip"][()])+".")
            return
        elif (variable_name+"_weights") not in f.keys():
            print("              ---> Calculate phase shift between "+y1_variable+" and "+y2_variable+" with skip = "+str(skip)+".")

    # Initiate
    phase_shifts_vs_tangleky = np.ones((len(itstarts), 361, len(vec_ky)))*np.nan
    phase_shifts_vs_tanglekx = np.ones((len(itstarts), 361, len(vec_kx)))*np.nan
    weights_vs_t = []
    
    # Put the angles in a distibution between -180° and 180°
    bins = np.linspace(-180, 180, 361)
      
    # Iterate over the modes
    for it in range(len(itstarts)):    
        
        # Read the potential and moments from the output file  
        if y1_variable=="phi":          y1_vs_tzkxky = read_phi(netcdf_path, itstarts[it], itends[it], skip) 
        if y2_variable=="phi":          y2_vs_tzkxky = read_phi(netcdf_path, itstarts[it], itends[it], skip) 
        if y1_variable in ["n0","n1","n2"]:  y1_vs_tzkxky = read_moments(netcdf_path, 'dens_vs_tszkxkyri', itstarts[it], itends[it], skip, specie=int(y1_variable[1]))
        if y2_variable in ["n0","n1","n2"]:  y2_vs_tzkxky = read_moments(netcdf_path, 'dens_vs_tszkxkyri', itstarts[it], itends[it], skip, specie=int(y2_variable[1]))
        if y1_variable in ["T0","T1","T2"]:  y1_vs_tzkxky = read_moments(netcdf_path, 'temp_vs_tszkxkyri', itstarts[it], itends[it], skip, specie=int(y1_variable[1]))
        if y2_variable in ["T0","T1","T2"]:  y2_vs_tzkxky = read_moments(netcdf_path, 'temp_vs_tszkxkyri', itstarts[it], itends[it], skip, specie=int(y2_variable[1]))
                
        # Unmask the arrays
        if isinstance(y1_vs_tzkxky, np.ma.core.MaskedArray): y1_vs_tzkxky = np.ma.filled(y1_vs_tzkxky, np.nan)   
        if isinstance(y2_vs_tzkxky, np.ma.core.MaskedArray): y2_vs_tzkxky = np.ma.filled(y2_vs_tzkxky, np.nan)          
                
        # The weight of each angle of point (kx,ky,z,t) is |y1(z,ky)|*|y2(z,ky)|*dell/B 
        weights_vs_tzkxky = np.abs(y1_vs_tzkxky)*np.abs(y2_vs_tzkxky)*dl_over_B[np.newaxis,:,np.newaxis,np.newaxis]
        weights_vs_t.append(np.sum(weights_vs_tzkxky, axis=(0,1,2,3)))
        weights_vs_tzkxky = weights_vs_tzkxky/weights_vs_t[-1]
                        
        # Iterate over the ky-modes                 
        for iky in range(len(vec_ky)): 
        
            # Get the data versus (t,z,kx) for a specific ky
            y1_vs_tzkx = copy.deepcopy(y1_vs_tzkxky[:,:,:,iky])
            y2_vs_tzkx = copy.deepcopy(y2_vs_tzkxky[:,:,:,iky]) 
            
            # Get the phase shift between (y1) and (y2)
            phase_shift_vs_tzkx = get_phase_shift(y1_vs_tzkx, y2_vs_tzkx) 
            
            # To each angle in <vec_phaseshifts_ni_phi> asign its bin <i> in <digitized>
            digitized = np.digitize(phase_shift_vs_tzkx, bins)
            
            # For each angle in <vec_phaseshifts_ni_phi>, add its weight to the angle in <bins>
            weighted_data_in_bins = [np.sum(weights_vs_tzkxky[:,:,:,iky][digitized==i+1]) for i in range(len(bins))]

            # Save the sum of the weights of each angle in <bins> np.ma.filled(marr, np.nan)
            phase_shifts_vs_tangleky[it,:,iky] = np.array(weighted_data_in_bins)
                        
        # Iterate over the ky-modes                 
        for ikx in range(len(vec_kx)): 
         
            # Get the data versus (t,z,kx) for a specific ky
            y1_vs_tzky = copy.deepcopy(y1_vs_tzkxky[:,:,ikx,:])
            y2_vs_tzky = copy.deepcopy(y2_vs_tzkxky[:,:,ikx,:]) 
             
            # Get the phase shift between (y1) and (y2)
            phase_shift_vs_tzky = get_phase_shift(y1_vs_tzky, y2_vs_tzky) 
             
            # To each angle in <vec_phaseshifts_ni_phi> asign its bin <i> in <digitized>
            digitized = np.digitize(phase_shift_vs_tzky, bins)
             
            # For each angle in <vec_phaseshifts_ni_phi>, add its weight to the angle in <bins>
            weighted_data_in_bins = [np.sum(weights_vs_tzkxky[:,:,ikx,:][digitized==i+1]) for i in range(len(bins))]

            # Save the sum of the weights of each angle in <bins>
            phase_shifts_vs_tanglekx[it,:,ikx] = np.array(weighted_data_in_bins)
        
        # Bin 0 and bin -1 will we filled with angles before -180 and after 180 degrees
        phase_shifts_vs_tangleky[:,-1,:] = 0 
        phase_shifts_vs_tanglekx[:,-1,:] = 0 

    # Add the data to the h5 file
    with h5py.File(phaseshifts_path, 'r+') as h5_file:   
        
        # First remove the data if it already existed
        if variable_name+"_vs_tangleky" in h5_file.keys(): del h5_file[variable_name+"_vs_tangleky"]
        if variable_name+"_vs_tanglekx" in h5_file.keys(): del h5_file[variable_name+"_vs_tanglekx"]
        if variable_name+"_weights" in h5_file.keys(): del h5_file[variable_name+"_weights"]
        if variable_name+"_skip" in h5_file.keys(): del h5_file[variable_name+"_skip"]
        
        # Then write the data
        h5_file.create_dataset(variable_name+"_vs_tangleky", data=phase_shifts_vs_tangleky)   
        h5_file.create_dataset(variable_name+"_vs_tanglekx", data=phase_shifts_vs_tanglekx)   
        h5_file.create_dataset(variable_name+"_weights",     data=weights_vs_t)  
        h5_file.create_dataset(variable_name+"_skip",        data=skip) 
        
    # Clean up the memory 
    del phase_shifts_vs_tangleky
    del phase_shifts_vs_tanglekx
    del weights_vs_t
    del y1_vs_tzkxky
    del y2_vs_tzkxky
    return

#---------------------------
def get_phase_shift(y1, y2):
    vec_phaseshifts = np.angle(y1) - np.angle(y2)
    vec_phaseshifts[vec_phaseshifts>np.pi] = vec_phaseshifts[vec_phaseshifts>np.pi]-2*np.pi
    vec_phaseshifts[vec_phaseshifts<-np.pi] = vec_phaseshifts[vec_phaseshifts<-np.pi]+2*np.pi
    return np.array(vec_phaseshifts)/np.pi*180

#---------------------------------------------
def get_integrationWeightsAlongZ(input_file):
    vmec_filename = read_vmecFileNameFromInputFile(input_file) 
    path = create_dummyPathObject(input_file, vmec_filename)
    geometry = read_outputFileForGeometry(path)  
    dl_over_B = geometry["dl_over_B"]   
    return dl_over_B
  
#---------------------------------------------
def read_phi(netcdf_path, itstart=0, itend=-1, skip=1): 
      
    # First get the potential ["t", "tube", "z", "kx", "ky", "ri"]
    try: netcdf_file = nc4.Dataset(netcdf_path)  
    except: netcdf_file = scnetcdf.netcdf_file(netcdf_path,'r') 
    phi_vs_tzkxkyri = copy.deepcopy(netcdf_file.variables["phi_vs_t"][itstart:itend:skip,0,:,:,:,:]) 
    netcdf_file.close()  
      
    # Combine the real and imaginary part to a complex number
    phi_vs_tzkxky = phi_vs_tzkxkyri[:,:,:,:,0] + 1j*phi_vs_tzkxkyri[:,:,:,:,1]
      
    # Return the data 
    return phi_vs_tzkxky
  
#---------------------------------------------
def read_moments(netcdf_path, variable, itstart=0, itend=-1, skip=1, specie=0): 
    
    # Get the dimensions and the stella key of the variable 
    key = stella_variables[variable][0]
    
    # Open the netcdf file
    try: netcdf_file = nc4.Dataset(netcdf_path)  
    except: netcdf_file = scnetcdf.netcdf_file(netcdf_path,'r') 
    variable_vs_tzkxkyri = copy.deepcopy(netcdf_file.variables[key][itstart:itend:skip,specie,0,:,:,:,:]) 
    netcdf_file.close()    
     
    # Combine the real and imaginary part to a complex number
    variable_vs_tzkxky = variable_vs_tzkxkyri[:,:,:,:,0] + 1j*variable_vs_tzkxkyri[:,:,:,:,1] 
      
    # Return the data 
    return variable_vs_tzkxky


