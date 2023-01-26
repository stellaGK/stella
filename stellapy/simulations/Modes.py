""" 

#===============================================================================
#                  Create a <Mode> for each linear (kx,ky)                     #
#===============================================================================      

Hanne Thienpondt
08/09/2022

"""

#!/usr/bin/python3 
import h5py
import sys, os
import numpy as np
from scipy.io import netcdf as scnetcdf   

# Stellapy package
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])  
from stellapy.data.output.read_outputFile import read_outputFile 
from stellapy.data.output.read_outputFile import read_netcdfVariables
from stellapy.data.input.read_inputFile import read_vecKxKyFromInputFile
from stellapy.data.input.read_inputFile import read_numberOfModesFromInputFile
from stellapy.data.potential import load_potentialObject
from stellapy.data.paths.load_pathObject import load_pathObject
from stellapy.data.input.load_inputObject import load_inputObject
from stellapy.data.omega.load_omegaObject import load_omegaObject
from stellapy.data.fluxes.load_fluxesObject import load_fluxesObject
from stellapy.data.moments.load_momentsObject import load_momentsObject
from stellapy.data.geometry.load_geometryObject import load_geometryObject 
from stellapy.data.dimensions.load_vectorsObject import load_vectorsObject
from stellapy.data.saturated.load_saturatedObject import load_saturatedObject
from stellapy.data.dimensions.load_dimensionsObject import load_dimensionsObject
from stellapy.data.lineardata.load_lineardataObject import load_linearDataObject
from stellapy.data.referenceunits.load_referenceObject import load_referenceObject
from stellapy.data.distribution.load_distributionObject import load_distributionObject
from stellapy.simulations.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime
 
################################################################################
#                                CREATE MODES                                  #
################################################################################

def create_modes(input_files, simulation, write_uniqueFiles):
    """ Create a mode object for each mode (kx,ky) in a certain <simulatiVectorsson>.
    This only works for linear simulations since each mode is simulated seperatly. 
    A simulation containing multiple modes can have different resolutions in 
    (nzgrid, nvrgrid, nmu, dt) to ensure converge, therefore, make sure we read
    these unique dimensions, without using too much memory or time. """
    
    # Initiate the vector of modes 
    modes = []
    
    # Find the number of modes in each input file
    simulation.modes_per_input_file = {}
    for input_file in input_files: 
        nakx, naky = read_numberOfModesFromInputFile(input_file) 
        simulation.modes_per_input_file[input_file] = [nakx, naky]
         
    # Create a <mode> object for each (kx,ky)
    for input_file in input_files:
        nakx, naky = simulation.modes_per_input_file[input_file]
        for ikx in range(nakx):
            for iky in range(naky):
                modes.append(Mode(input_file, simulation, ikx, iky))
                
    # Attach the vector of modes to the simulation object
    simulation.modes = modes
    
    # Now that the <mode> objects are made, assign the correct (kx,ky) 
    if not write_uniqueFiles: 
        for mode in modes: 
            vec_kx, vec_ky = read_vecKxKyFromInputFile(mode.input_file) 
            mode.kx = vec_kx[mode.ikx]
            mode.ky = vec_ky[mode.iky] 
            
    # If we are writing the list of unique input files, then don't use <mode.vec>
    if write_uniqueFiles:
        vectors_per_input_file = {}
        for mode in modes:    
            if input_file in vectors_per_input_file:
                vec_kx = vectors_per_input_file[input_file][0]
                vec_ky = vectors_per_input_file[input_file][1] 
            elif os.path.isfile(mode.input_file.with_suffix(".dimensions")):   
                with h5py.File(mode.input_file.with_suffix(".dimensions")) as f:  
                    vec_kx = f["vec_kx"][()] 
                    vec_ky = f["vec_ky"][()] 
                    vectors_per_input_file[input_file] = [vec_kx, vec_ky]
            elif os.path.isfile(mode.input_file.with_suffix(".out.nc")):
                netcdf_file = read_outputFile(mode.input_file.with_suffix(".out.nc"))   
                vec_kx = read_netcdfVariables('vec_kx', netcdf_file)
                vec_ky = read_netcdfVariables('vec_ky', netcdf_file) 
                vectors_per_input_file[input_file] = [vec_kx, vec_ky]
                netcdf_file.close()
            mode.kx = vec_kx[mode.ikx]
            mode.ky = vec_ky[mode.iky]   
            
    # Replace kx=-0 by kx=0
    for mode in modes:
        if mode.kx==-0: mode.kx = 0
                  
    # Return the list of modes
    return modes
 
################################################################################
#                                 CLASS MODE                                   #
################################################################################
 
class Mode: 
     
    def __init__(self, input_file, simulation, ikx, iky): 
        
        # Remember this is a mode
        self.object = "Mode"
        
        # Save the input file
        self.input_file = input_file
        
        # Remember the parent simulation  
        self.Progress = simulation.Progress
        self.simulation = simulation 
        self.linear = simulation.linear
        self.nonlinear = simulation.nonlinear
        
        # Remember the (ikx,iky)
        self.ikx = ikx
        self.iky = iky
        return 

    # Read the data in the various files when it's asked for
    @calculate_attributeWhenReadFirstTime
    def dim(self):          load_dimensionsObject(self);    return self.dim  
    @calculate_attributeWhenReadFirstTime
    def vec(self):          load_vectorsObject(self);       return self.vec  
    @calculate_attributeWhenReadFirstTime
    def path(self):         load_pathObject(self);          return self.path
    @calculate_attributeWhenReadFirstTime
    def input(self):        load_inputObject(self);         return self.input  
    @calculate_attributeWhenReadFirstTime
    def omega(self):        load_omegaObject(self);         return self.omega
    @calculate_attributeWhenReadFirstTime
    def fluxes(self):       load_fluxesObject(self);        return self.fluxes
    @calculate_attributeWhenReadFirstTime
    def moments(self):      load_momentsObject(self);        return self.moments
    @calculate_attributeWhenReadFirstTime 
    def geometry(self):     load_geometryObject(self);      return self.geometry  
    @calculate_attributeWhenReadFirstTime 
    def saturated(self):    load_saturatedObject(self);     return self.saturated  
    @calculate_attributeWhenReadFirstTime 
    def potential(self):    load_potentialObject(self);     return self.potential  
    @calculate_attributeWhenReadFirstTime 
    def lineardata(self):   load_linearDataObject(self);    return self.lineardata  
    @calculate_attributeWhenReadFirstTime 
    def distribution(self): load_distributionObject(self);  return self.distribution   
    @calculate_attributeWhenReadFirstTime 
    def referenceunits(self):load_referenceObject(self);    return self.referenceunits   

#------------------------------------
def finish_initializationOfModes(simulation): 
    
    # Only finish the simulations if we're not writing the *.ini files      
    if sys._getframe().f_back.f_back.f_back.f_code.co_name!="write_iniFileForInputs" \
    and sys._getframe().f_back.f_back.f_back.f_code.co_name!="write_h5FileForGeometry":
    
        # First make sure we can access the input data
        for mode in simulation.modes: mode.input
        
        # Make sure we don't have a mode with kx==0 and ky==0
        simulation.modes = [m for m in simulation.modes if not (m.kx==0 and m.ky==0)]
         
        # Sort the modes based on (kx,ky) 
        vec_kx_ky = np.array([1000*mode.ky+mode.kx for mode in simulation.modes])  
        simulation.modes = [simulation.modes[i] for i in vec_kx_ky.argsort()]
         
        # Remove the modes that are removed in <simulation.ini> 
        remove_modes = []
        for mode in simulation.modes:  
            if mode.lineardata.removeMode: remove_modes.append(mode) 
        for mode in remove_modes:
            simulation.modes.remove(mode)  
                
    return  
