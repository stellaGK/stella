
import os, sys  
import numpy as np
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program  
from stellapy.data.saturated.write_netcdfFileOverSaturatedPhase import write_netcdfFileOverSaturatedPhase 
from stellapy.data.output.read_outputFile import read_netcdfVariables, read_outputFile 

#===============================================================================
#                        ATTACH THE POTENTIAL DATA
#===============================================================================

def read_saturatedNetcdf(path): 
    ''' Read the potential from *.potential4D or *.out.nc ''' 
               
    # Make sure that the *.out.nc.t500-1000 file exists
    if not os.path.isfile(path.saturated) and os.path.isfile(path.output_stella): 
        write_netcdfFileOverSaturatedPhase(path.folder)
        get_saturatedData(path)  
                              
    # Read from the *.out.nc.t500-1000 file 
    if os.path.isfile(path.saturated):
        return read_fromSaturatedNetcdf(path.saturated)    

    # Critical error if we didn't find any data
    exit_reason = "The saturated data could not be found for:\n"
    exit_reason += "     "+str(path.input_file)
    exit_program(exit_reason, read_saturatedNetcdf, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromSaturatedNetcdf(path):  
    
    # Intitiate
    if "tlast" not in path.name:
        tend = float(path.name.split("-")[-1])
        tstart = float(path.name.split("out.nc.t")[-1].split("-")[0])
        saturated = {"trange" : [tstart, tend]}
    if "tlast" in path.name:
        saturated = {"trange" : [-1, -1]}
    
    # Read the netcdf file
    netcdf_data = read_outputFile(path)   
    
    # Read the complex data
    try:
        phi_vs_zkxkyri = read_netcdfVariables('phi_vs_tzkxkyri', netcdf_data)[0,:,:,:,:]
        saturated["phi_vs_zkxky"] =  phi_vs_zkxkyri[:,:,:,0] + 1j*phi_vs_zkxkyri[:,:,:,1]
    except:
        print("REWRITE "+str(path))
        print("SOME VARIABLES ARE MISSING")
        sys.exit()
        
    # Read the real data
    try:
        saturated["g2_vs_svpaz"]  = read_netcdfVariables('g2_vs_tsvpaz', netcdf_data)[0,:,:,:]
        saturated["g2_vs_smuvpa"] = read_netcdfVariables('g2_vs_tsmuvpa', netcdf_data)[0,:,:,:]
        saturated["phi2_vs_kxky"] = read_netcdfVariables('phi2_vs_tkxky', netcdf_data)[0,:,:]
    except:
        print("REWRITE "+str(path))
        print("SOME VARIABLES ARE MISSING")
        sys.exit()
        
    # Read the moments
    if "density" in netcdf_data.variables.keys():
        dens_vs_szkxkyri = read_netcdfVariables('dens_vs_tszkxkyri', netcdf_data)[0,:,:,:,:]
        upar_vs_szkxkyri = read_netcdfVariables('upar_vs_tszkxkyri', netcdf_data)[0,:,:,:,:]
        temp_vs_szkxkyri = read_netcdfVariables('temp_vs_tszkxkyri', netcdf_data)[0,:,:,:,:]
        saturated["dens_vs_szkxky"] =  dens_vs_szkxkyri[:,:,:,:,0] + 1j*dens_vs_szkxkyri[:,:,:,:,1] 
        saturated["upar_vs_szkxky"] =  upar_vs_szkxkyri[:,:,:,:,0] + 1j*upar_vs_szkxkyri[:,:,:,:,1] 
        saturated["temp_vs_szkxky"] =  temp_vs_szkxkyri[:,:,:,:,0] + 1j*temp_vs_szkxkyri[:,:,:,:,1]   
    
    # Read the fluxes
    if "pflx_kxky" in netcdf_data.variables.keys():
        saturated["pflux_vs_szkxky"]  = read_netcdfVariables('pflux_vs_tszkxky', netcdf_data)[0,:,:,:,:]
        saturated["vflux_vs_szkxky"]  = read_netcdfVariables('vflux_vs_tszkxky', netcdf_data)[0,:,:,:,:]
        saturated["qflux_vs_szkxky"]  = read_netcdfVariables('qflux_vs_tszkxky', netcdf_data)[0,:,:,:,:] 
    return saturated 
            
#===============================================================================
#                  ATTACH THE POTENTIAL TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_saturatedData(self): 
    
    # Read the potential
    saturated = read_saturatedNetcdf(self.path)  
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
    sorted_indexes = np.argsort(self.vec.kx_stella)  
    saturated["phi2_vs_kxky"]    = saturated["phi2_vs_kxky"][sorted_indexes,:]
    saturated["phi_vs_zkxky"]    = saturated["phi_vs_zkxky"][:,sorted_indexes,:] 
    if "pflux_vs_szkxky" in saturated:
        saturated["pflux_vs_szkxky"] = saturated["pflux_vs_szkxky"][:,:,sorted_indexes,:] 
        saturated["vflux_vs_szkxky"] = saturated["vflux_vs_szkxky"][:,:,sorted_indexes,:] 
        saturated["qflux_vs_szkxky"] = saturated["qflux_vs_szkxky"][:,:,sorted_indexes,:] 
    if "dens_vs_szkxky" in saturated:
        saturated["dens_vs_szkxky"] = saturated["dens_vs_szkxky"][:,:,sorted_indexes,:] 
        saturated["upar_vs_szkxky"] = saturated["upar_vs_szkxky"][:,:,sorted_indexes,:] 
        saturated["temp_vs_szkxky"] = saturated["temp_vs_szkxky"][:,:,sorted_indexes,:] 
    
    # Save the potential
    self.trange = saturated["trange"]
    self.phi2_vs_kxky = Data(["phi2","kx","ky"], saturated["phi2_vs_kxky"], self.vec.kx, self.vec.ky) 
    self.phi_vs_zkxky = Data(["phi","z","kx","ky"], saturated["phi_vs_zkxky"], self.vec.z, self.vec.kx, self.vec.ky) 
    self.g2_vs_svpaz = Data(["g2","s","vpa","z"], saturated["g2_vs_svpaz"], self.vec.species, self.vec.vpa, self.vec.z) 
    self.g2_vs_smuvpa = Data(["g2","s","mu","vpa"], saturated["g2_vs_smuvpa"], self.vec.species, self.vec.mu, self.vec.vpa) 
    
    # Save the fluxes
    if "pflux_vs_szkxky" in saturated:
        self.pflux_vs_szkxky = Data(["pflux","s","z","kx","ky"], saturated["pflux_vs_szkxky"], self.vec.species, self.vec.z, self.vec.kx, self.vec.ky) 
        self.vflux_vs_szkxky = Data(["vflux","s","z","kx","ky"], saturated["vflux_vs_szkxky"], self.vec.species, self.vec.z, self.vec.kx, self.vec.ky)
        self.qflux_vs_szkxky = Data(["qflux","s","z","kx","ky"], saturated["qflux_vs_szkxky"], self.vec.species, self.vec.z, self.vec.kx, self.vec.ky) 
    else:
        self.pflux_vs_szkxky = None
        self.vflux_vs_szkxky = None
        self.qflux_vs_szkxky = None
        
    # Save the moments
    if "dens_vs_szkxky" in saturated:
        self.dens_vs_szkxky = Data(["dens","s","z","kx","ky"], saturated["dens_vs_szkxky"], self.vec.species, self.vec.z, self.vec.kx, self.vec.ky)
        self.upar_vs_szkxky = Data(["upar","s","z","kx","ky"], saturated["upar_vs_szkxky"], self.vec.species, self.vec.z, self.vec.kx, self.vec.ky)
        self.temp_vs_szkxky = Data(["temp","s","z","kx","ky"], saturated["temp_vs_szkxky"], self.vec.species, self.vec.z, self.vec.kx, self.vec.ky)
    else:
        self.dens_vs_szkxky = None
        self.upar_vs_szkxky = None
        self.temp_vs_szkxky = None
    return 
