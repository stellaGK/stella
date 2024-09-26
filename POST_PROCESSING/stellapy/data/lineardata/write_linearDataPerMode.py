
import numpy as np
import configparser 
from stellapy.data.lineardata.get_modeIdentifier import get_modeIdentifier

#===============================================================================
#                             WRITE THE LINEAR DATA                            #
#=============================================================================== 

def write_linearDataPerMode(omega=None, dim=None, vec=None,
        input_file=None, file=None, percentage=0.85, self=None, verbose=False):
    ''' Calculate the linear data for each mode (kx,ky). ''' 
    
    # Extract the data
    if self!=None:
        input_file = self.input_file
        percentage = self.percentage
        omega = self.omega 
        dim = self.dim
        vec = self.vec 
    omega_data = omega

    # Write the linear data to "*.linear_data" so open the file 
    path = input_file.with_suffix(".lineardata.ini")
    if file==None: file = configparser.ConfigParser(); file.read(path) 
    
    # Iterate over the modes
    for ikx in range(dim.kx):
        for iky in range(dim.ky): 
            
            # Get the mode identifier 
            mode_identifier, kx, ky = get_modeIdentifier(vec.kx[ikx], vec.ky[iky])
            
            # Get the data for this mode
            omega = omega_data.omega_vs_tkxky.omega[:,ikx,iky]
            gamma = omega_data.gamma_vs_tkxky.gamma[:,ikx,iky] 
    
            # Gamma and omega can be nans 
            if len(omega[np.isfinite(omega)])==0 or(np.isclose(kx,0) and np.isclose(ky,0)):       
                file[mode_identifier] = {
                    'stable'     : np.nan,\
                    'unstable'   : np.nan,\
                    'converged'  : np.nan,\
                    'gamma_max'  : np.nan,\
                    'omega_last' : np.nan,\
                    'gamma_last' : np.nan,\
                    'omega_avg'  : np.nan,\
                    'gamma_avg'  : np.nan,\
                    'omega_min'  : np.nan,\
                    'omega_max'  : np.nan,\
                    'gamma_min'  : np.nan,\
                    'gamma_max'  : np.nan} 
                
            # If no wave is examined then omega=0 and gamma=0 and we divide by 0
            else:  
                
                # Remove the infinite and NaN values 
                gamma = gamma[np.isfinite(omega)]
                omega = omega[np.isfinite(omega)]   
                
                # Grab omega and gamma at the last time value where omega is finite. 
                omega_last = omega[-1] 
                gamma_last = gamma[-1]
                
                # Reduce (omega; gamma) to the last x% of their data to get the stable state 
                omega_percentage = omega[int(np.size(omega)*percentage):]
                gamma_percentage = gamma[int(np.size(gamma)*percentage):] 
                
                # Average over the last x% since before it converges it osicillates around its convergent point.
                omega_avg = np.average(omega_percentage)
                gamma_avg = np.average(gamma_percentage)

                # Calculate the error bars over the last x%: 
                omega_min = np.abs(np.nanmin(omega_percentage)-omega_avg)
                omega_max = np.abs(np.nanmax(omega_percentage)-omega_avg)
                gamma_min = np.abs(np.nanmin(gamma_percentage)-gamma_avg)
                gamma_max = np.abs(np.nanmax(gamma_percentage)-gamma_avg) 
                
                # Intiate each mode to be stable and not converged
                stable = True
                unstable = False
                converged = False
                     
                # A mode is unstable when gamma does not change sign in the last 20% of the time trace 
                last_20percent  = int(np.size(gamma)*8/10) 
                gamma_20percent = gamma[last_20percent:]
                sign_changes    = np.where(np.diff(np.sign(gamma_20percent)))[0] 
                if len(sign_changes)==0 and gamma[-1]>0:  
                    stable = False
                    unstable = True 
                    
                # For the unstable modes, check whether the simulation has converged (constant gamma)
                # Consider the mode converged if gamma/omega doesn't change more than 2% over the last 10% of time 
                converged = np.all( np.isclose(omega[int(len(omega)*0.90):], omega[-1], rtol=0.02) )
                converged = np.all( np.isclose(gamma[int(len(gamma)*0.90):], gamma[-1], rtol=0.02) ) & converged
                
                # Add some information to the command prompt if the mode is not converged
                if converged == False and verbose: 
                    fluctuation_o = (max(omega[int(len(omega)*0.90):])-min(omega[int(len(omega)*0.90):]))/omega[-1]*100
                    fluctuation_g = (max(gamma[int(len(gamma)*0.90):])-min(gamma[int(len(gamma)*0.90):]))/gamma[-1]*100
                    print("\nThe simulation for (kx,kx) = ("+str(kx)+", "+str(ky)+") has not converged yet:")
                    print("    The relative fluctuation of omega between (0.9*t, t) is:  "+str(fluctuation_o)+"%.")
                    print("    The relative fluctuation of gamma between (0.9*t, t) is:  "+str(fluctuation_g)+"%.")
         
                # Write the data to the "*.linear_data" file 
                file[mode_identifier] = {
                    'stable'     : stable,\
                    'unstable'   : unstable,\
                    'converged'  : converged,\
                    'gamma_max'  : gamma_max,\
                    'omega_last' : omega_last,\
                    'gamma_last' : gamma_last,\
                    'omega_avg'  : omega_avg,\
                    'gamma_avg'  : gamma_avg,\
                    'omega_min'  : omega_min,\
                    'omega_max'  : omega_max,\
                    'gamma_min'  : gamma_min,\
                    'gamma_max'  : gamma_max} 
                    
                # Write the "*.linear_data" file 
                file.write(open(path, 'w'))   
    return    


