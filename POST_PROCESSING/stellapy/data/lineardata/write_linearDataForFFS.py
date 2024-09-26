 
import sys
import numpy as np
import configparser  
from stellapy.data.lineardata.write_linearDataPerMode import write_linearDataPerMode
from stellapy.utils.decorators.exit_program import exit_program
    
#===============================================================================
#           WRITE THE LINEAR DATA FOR A FULL FLUX SURFACE SIMULATION           #
#=============================================================================== 

def write_linearDataForFFS(omega, dim=None, vec=None,
        input_file=None, full_flux_surface=None, file=None, percentage=0.85, self=None):
    
    # Extract the data
    if self!=None: 
        full_flux_surface = self.input.full_flux_surface
        input_file = self.input_file
        percentage = self.percentage
        dim = self.dim; vec = self.vec
        omega = self.omega
        
    # For a full flux surface simulation, save the overall gamma
    if full_flux_surface==True:

        # Read the linear data file
        path = input_file.with_suffix(".lineardata.ini") 
        if file==None: file = configparser.ConfigParser(); file.read(path)  
        
        # Write if it doesn't exists
        if "full flux surface" not in file:
            
            # Store (gamma, omega) in a matrix
            gamma_avg = np.ones((dim.kx, dim.ky)) * np.nan
            omega_avg = np.ones((dim.kx, dim.ky)) * np.nan
            gamma_min = np.ones((dim.kx, dim.ky)) * np.nan
            gamma_max = np.ones((dim.kx, dim.ky)) * np.nan
            omega_min = np.ones((dim.kx, dim.ky)) * np.nan
            omega_max = np.ones((dim.kx, dim.ky)) * np.nan
            
            # Read the linear data per mode
            for ikx in range(dim.kx):
                for iky in range(dim.ky):
                    
                    # Get the mode identifier
                    kx = vec.kx[ikx]; ky = vec.ky[iky]
                    mode_identifier = '('+str(kx)+', '+str(ky)+')'  
                    
                    # Make sure it it written 
                    if mode_identifier not in file:  
                        write_linearDataPerMode(omega, dim, vec, input_file, file, percentage)
                    gamma_avg[ikx, iky] = float(file[mode_identifier]["gamma_avg"])
                    omega_avg[ikx, iky] = float(file[mode_identifier]["omega_avg"])
                    gamma_min[ikx, iky] = float(file[mode_identifier]["gamma_min"])
                    gamma_max[ikx, iky] = float(file[mode_identifier]["gamma_max"])
                    omega_min[ikx, iky] = float(file[mode_identifier]["omega_min"])
                    omega_max[ikx, iky] = float(file[mode_identifier]["omega_max"])
                
            # Get the relative error
            relative_error_min = np.array(gamma_min)/np.array(gamma_avg) 
            relative_error_max = np.array(gamma_max)/np.array(gamma_avg)  
            
            # The following code doesn't work for dim.kx!=0
            if dim.kx > 1:
                exit_program("Not implemented for dim.kx>1. ", write_linearDataForFFS, sys._getframe().f_lineno)
                
            # Find the converged modes (defined as error < error_max)
            indices = [i for i in range(len(gamma_avg)) if (relative_error_min[i]<self.error_max and relative_error_max[i]<self.error_max and gamma_min[i]<self.error_max and gamma_max[i]<self.error_max)]
            if len(indices)==0: indices = range(len(gamma_avg)) 
            
            # Grab all modes with ky bigger than the ky of the first converged mode
            indices = [i for i in range(len(gamma_avg)) if i>=indices[0]] 
            omega_avg = np.array(omega_avg)[indices]
            gamma_avg = np.array(gamma_avg)[indices]
            gamma_min = np.array(gamma_min)[indices] 
            gamma_max = np.array(gamma_max)[indices] 
            omega_min = np.array(omega_min)[indices] 
            omega_max = np.array(omega_max)[indices]
            vec_ky = np.array(vec.ky)[indices]
            
            # Get the average gamma by averaging right of the first converged mode
            average_gamma = np.nanmean(gamma_avg)  
            average_omega = np.nanmean(omega_avg)  
            average_gamma_min = np.nanmean(gamma_min)
            average_gamma_max = np.nanmean(gamma_max)
            average_omega_min = np.nanmean(omega_min)
            average_omega_max = np.nanmean(omega_max)
            
            # Write the data to the "*.linear_data" file 
            file["full flux surface"] = {
                'gamma'  : average_gamma,\
                'omega'  : average_omega,\
                'gamma_min' : average_gamma_min,\
                'gamma_max' : average_gamma_max,\
                'omega_min' : average_omega_min,\
                'omega_max' : average_omega_max,\
                'kymin'  : vec_ky[0] } 
        
        # Write the "*.linear_data" file 
        file.write(open(path, 'w'))   
        
        # If the data is written, read it  
        self.ffs_gamma = float(file["full flux surface"]["gamma"]) 
        self.ffs_omega[ikx] = float(file["full flux surface"]["omega"]) 
        self.ffs_gamma_min[ikx] = float(file["full flux surface"]["gamma_min"]) 
        self.ffs_gamma_max[ikx] = float(file["full flux surface"]["gamma_max"]) 
        self.ffs_omega_min[ikx] = float(file["full flux surface"]["omega_min"]) 
        self.ffs_omega_max[ikx] = float(file["full flux surface"]["omega_max"]) 
        self.ffs_kymin[ikx] = float(file["full flux surface"]["kymin"]) 
        return 
    
    # For a flux tube simulation, put nans
    if self.input.full_flux_surface==False:
        self.ffs_kymin = [np.nan]
        self.ffs_gamma = [np.nan]
        self.ffs_omega = [np.nan]
        self.ffs_gamma_min = [np.nan]
        self.ffs_gamma_max = [np.nan]
        self.ffs_omega_min = [np.nan]
        self.ffs_omega_max = [np.nan]
        return         
    