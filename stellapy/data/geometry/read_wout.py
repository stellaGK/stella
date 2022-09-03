
import sys 
import numpy as np
import os, h5py, netCDF4, pathlib
from stellapy.utils.config import CONFIG 
from stellapy.utils.decorators.exit_program import exit_program
 
#===============================================================================
#                             READ THE WOUT FILE
#===============================================================================
 
def read_woutFile(path):
    ''' Read the geometry file or the VMEC file. '''
    
    # Read from the *.geo.h5 or *.unique.geometryX.h5 file  
    if sys._getframe().f_back.f_back!=None:
        if sys._getframe().f_back.f_back.f_back!=None:
            if sys._getframe().f_back.f_back.f_back.f_back!=None:
                if sys._getframe().f_back.f_back.f_back.f_back.f_code.co_name!="write_iniFileForInputs":
                    if os.path.isfile(path.geometry):
                        return read_woutFromH5File(path.geometry)  
                    
    # Read from the VMEC file   
    if os.path.isfile(path.vmec):
        return read_woutFromNcFile(path.vmec)  
          
    # Miller coordinates are not implemented yet (can't determine iota)
    if path.vmec_filename=='wout*.nc':   
        return read_millerParameters(path)

    # Critical error if we didn't find any data  
    exit_reason = "The wout data can not be found.\n"
    exit_reason += "The vmec file: '"+path.vmec_filename+"' was not found in "+CONFIG['PATHS']['VMEC']+".\n"
    exit_reason += "Please make sure that the vmec file is available."
    exit_program(exit_reason, read_woutFile, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_millerParameters(path): 
    from stellapy.data.input.read_inputFile import read_inFile
    input_parameters = read_inFile(path.input_file)
    wout_variables = {"source" : "Miller coordinates"} 
    wout_variables["nfp"] = 1 
    wout_variables["shat"] = input_parameters["millergeo_parameters"]["shat"]
    wout_variables["iota"] = 1/input_parameters["millergeo_parameters"]["qinp"]
    return wout_variables
    
#-------------------------------------
def read_woutFromH5File(path): 
    
    # Recall what file we are reading from
    wout_variables = {"source" : path.name}
    
    # Read the wout data from the "*.geo.h5" file
    with h5py.File(path, 'r') as f: 
        for key in ["iota", "nfp", "sign_B", "b0", "aminor", "rmajor", "volume", "vmec_path"]:
            if key in f.keys(): wout_variables[key] = f[key][()]
        for key in ["P", "Lx", "Ly", "dkx", "dky", "shat", "iota", "jtwist", "diotaOverds", "twist_and_shift_geo_fac"]:
            if key in f.keys(): wout_variables[key] = f[key][()]
            
    # Return the wout variables
    return wout_variables

#-------------------------------------
def read_woutFromNcFile(wout_file_path):
 
    # Read the netcdf file  
    dataset  = netCDF4.Dataset(wout_file_path, mode='r')
     
    # Safe the data to a dictionary
    wout_variables              = {} 
    wout_variables['nfp']       = int(dataset.variables['nfp'][:])
    wout_variables['ns']        = int(dataset.variables['ns'][:])       # number of divisions along s=(r/a)^2=rho^2 during the VMEC run
    wout_variables['b0']        = dataset.variables['b0'][:]            # Magnetic field averaged along the axis
    wout_variables['phi']       = dataset.variables['phi'][:]           # Toroidal flux on full mesh 
    wout_variables['aminor']    = dataset.variables['Aminor_p'][:]
    wout_variables['rmajor']    = dataset.variables['Rmajor_p'][:]
    wout_variables['volume']    = dataset.variables['volume_p'][:]
    wout_variables['iotas']     = dataset.variables['iotas'][:]
    wout_variables['iotaf']     = dataset.variables['iotaf'][:]         # Iota on full mesh 
    wout_variables['vmec_path'] = wout_file_path
    wout_variables['source']    = wout_file_path.name
    return wout_variables

#-------------------------------------
def read_allVariablesFromWoutFile(wout_file_path):
 
    # Read the netcdf file  
    dataset  = netCDF4.Dataset(wout_file_path, mode='r')
     
    # Safe the data to a dictionary
    wout_variables              = {} 
    wout_variables['b0']        = dataset.variables['b0'][:]            # Magnetic field averaged along the axis
    wout_variables['phi']       = dataset.variables['phi'][:]           # Toroidal flux on full mesh 
    wout_variables['aminor']    = dataset.variables['Aminor_p'][:]
    wout_variables['rmajor']    = dataset.variables['Rmajor_p'][:]
    wout_variables['volume']    = dataset.variables['volume_p'][:]
    wout_variables['iotas']     = dataset.variables['iotas'][:]
    wout_variables['iotaf']     = dataset.variables['iotaf'][:]         # Iota on full mesh 
    
    # Poloidal and toroidal mode numbers
    wout_variables['xm'] = dataset.variables['xm'][:]             # Poloidal mode numbers  xm(mn_mode)   
    wout_variables['xn'] = dataset.variables['xn'][:]             # Toroidal mode numbers  xn(mn_mode)  
    wout_variables['xm_nyq'] = dataset.variables['xm_nyq'][:]     # Poloidal mode numbers (Nyquist) xm_nyq(mn_mode_nyq)
    wout_variables['xn_nyq'] = dataset.variables['xn_nyq'][:]     # Toroidal mode numbers (Nyquist) xn_nyq(mn_mode_nyq)
    
    # Dimensions
    wout_variables['ns']  = int(dataset.variables['ns'][:])               
    wout_variables['nfp'] = int(dataset.variables['nfp'][:]) 
    wout_variables['mnmax'] = int(dataset.variables['mnmax'][:]) 
    wout_variables['Rmin'] = dataset.variables['rmin_surf'][:] 
    wout_variables['Rmax'] = dataset.variables['rmax_surf'][:] 
    wout_variables['Zmax'] = dataset.variables['zmax_surf'][:] 
    wout_variables['lasym'] = dataset.variables['lasym__logical__'][:]  
     
    # Read rmnc, zmns and bmnc versus (radius,mn_mode)
    wout_variables['rmnc'] = dataset.variables['rmnc'][:]     # cosmn component of cylindrical R, full mesh  
    wout_variables['zmns'] = dataset.variables['zmns'][:]     # sinmn component of cylindrical Z, full mesh 
    wout_variables['bmnc'] = dataset.variables['bmnc'][:]     # cosmn component of mod-B, half mesh   
    wout_variables['lmns'] = dataset.variables['lmns'][:]     # cosmn component of lambda, half mesh   
    wout_variables['gmnc'] = dataset.variables['gmnc'][:]     # cosmn component of jacobian, half mesh   
    return wout_variables
 
 
#===============================================================================
#              ATTACH THE WOUT DATA TO THE SIMULATION OBJECT
#===============================================================================
 
def get_woutData(self):
    """ Attach the data inside the VMEC file to the simulation object. """
    
    # Read the wout parameters
    woutParameters = read_woutFile(self.path)
    
    # For miller, make sure we can set sign(B)
    if woutParameters["source"]=="Miller coordinates": woutParameters["b0"] = 1; woutParameters["nfp"] = 1 # 
    if "nfp" not in woutParameters: woutParameters["b0"] = 1; woutParameters["nfp"] = 1 # 
        
    # Vmec data
    self.nfp = woutParameters['nfp']  
    self.sign_B = np.sign(woutParameters['b0'])  
    return
 
 
def get_woutDataExtra(self):
    """ Attach the data inside the VMEC file to the simulation object. """
    
    # Read the wout parameters
    woutParameters = read_woutFile(self.path)

    # Calculate extra geometric quantities used in stella (depend on rho) 
    if "jtwist" not in woutParameters:
        from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_gridDivisionsAndSize 
        woutParameters.update(calculate_gridDivisionsAndSize(self.y0, self.nfield_periods, woutParameters, self.svalue))
    
    # For miller, set the parameters to nan
    if woutParameters["source"]=="Miller coordinates":
        woutParameters['vmec_path'] = "Miller coordinates"
        for key in ["b0", "aminor", "rmajor", "volume"]: woutParameters[key] = np.nan
        woutParameters["b0"] = 1 # This is not correct, but b0 is needed for sign(B)
      
    # Vmec data  
    self.b0 = woutParameters['b0']
    self.nfp = woutParameters['nfp'] 
    self.aminor = woutParameters['aminor']
    self.rmajor = woutParameters['rmajor']
    self.volume = woutParameters['volume'] 
    self.sign_B = np.sign(woutParameters['b0'])  
    self.vmec_path = str(woutParameters['vmec_path'])
    
    # Stella woutParameters
    self.P    = woutParameters['P']
    self.Lx   = woutParameters['Lx']
    self.Ly   = woutParameters['Ly']
    self.dkx  = woutParameters['dkx']
    self.dky  = woutParameters['dky']
    self.shat = woutParameters['shat']
    self.iota = woutParameters['iota']
    self.jtwist = woutParameters['jtwist']
    self.diotaOverds = woutParameters['diotaOverds']
    self.twist_and_shift_geo_fac = woutParameters['twist_and_shift_geo_fac'] 
    return
 
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__": 
    import timeit; start = timeit.timeit()
    from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile
    input_file = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_GEOMETRY_OBJECT/input_ky3.3125.in")
    vmec_filename = read_vmecFileNameFromInputFile(input_file) 
    wout_variables = read_woutFile(input_file, vmec_filename)
    print("\n"+"  "+"-"*30+"\n"+" "*12+"WOUT DATA"+"\n"+"  "+"-"*30+"\n")
    for key, value in wout_variables.items():
        print('{:15}'.format("  "+key), value)
    print("\n"+"  Elapsed time: "+str((start-timeit.timeit())*1000)+" ms")

            

