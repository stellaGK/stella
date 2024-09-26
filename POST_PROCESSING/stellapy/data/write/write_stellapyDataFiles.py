"""
#===============================================================================
#             Write the data files needed for the stellapy package             #
#===============================================================================

In order to speed-up the processing time and reduce the memory requirements, the 
stellapy package will write its own data files based on the stella data files. 
Nonetheless, all plotting commands should work directly with the stella data, but
it is recommended to write the stellapy data files instead.  
    \033[3m\033[1m>> write_stellapyDataFiles\033[0m

A specific set of data files is written by specifying the <dimension>
    >> write_stellapyDataFiles --small
    
A specific data file can be rewritten with a different time step (e.g. dt=1) through
    >> write_stellapyDataFiles -s pot5D -t 1
    
To speed up the writing of the phaseshifts file, we can skip time points
    >> write_stellapyDataFiles -s phases --skip 10
    
Arguments
---------
    dt : float
    dimensions : {small, medium, big, large, all, 1D, 2D, 3D, 4D, 5D}
    specific : {ini, geo, dim, omega, jump, flux, flux3D, pot, pot2D, pot3D, pot4D, pot5D, g, g2D, g3D, g4D, mom3D, mom4D, mom5D, sat, phases, freq}

Dimensions
----------
Based on <dimensions> different sets of data files will be written:

        if dimensions in ["1D", "2D", "small", "quick", "all"]:
            input.ini    
            input.geo      
            input.dimensions 
            input_v1_dummy.in                                   (LINEAR) 
        
        if dimensions in ["2D", "small", "quick", "all"]:     
            input.phi_vs_z               (z)                    (LINEAR) 
            input.g_vs_vpa               (vpa)                  (LINEAR) 
            input.g_vs_mu                (mu)                   (LINEAR) 
            input.g_vs_z                 (z)                    (LINEAR) 
            input.dt1.omega_vs_t         (t)                    (LINEAR)
            input.dt1.dphiz_vs_t         (t)                    (LINEAR)         
            input.dt1.fluxes_vs_t        (t)                    (NONLINEAR)
            input.dt1.phi_vs_t           (t)
            input.dt1.g_vs_t             (t) 
            input.dt1.moments2D          (t,s)
            XXXXXXXX_slurm.dt1.out
        
        if dimensions in ["3D", "small", "quick", "all"]: 
            input.dt10.fluxes3D           (t,s,(kx,ky,z))
            input.dt10.moments3D          (t,s,(kx,ky,z)) 
            input.dt10.potential3D        (t,(kx,ky,z))
            input.dt10.distribution3D     (t,s,(mu,vpa,z))   
            input.out.nc.tlast           (final time step)       (LINEAR)
            input.out.nc.t500-1000       (averaged over time)    (NONLINEAR) 
            
        if dimensions in ["3D", "small", "all"]: 
            input.dt50.phaseshifts       (t,angle,(kx,ky))       (NONLINEAR)  
        
        if dimensions in ["4D", "medium", "all"]: 
            input.dt10.fluxes4D           (t,s,kx,ky)
            input.dt10.moments4D          (t,s,kx,ky)
            input.dt10.potential4D        (t,kx,ky)   
            input.dt10.distribution4D     (t,s,vpa,(z,mu)) 
            
        if dimensions in ["5D", "big", "large", "all"]:
            input.dt100.moments5D          (t,s,z,kx,ky) 
            input.dt100.potential5D        (t,z,kx,ky) 
    

Hanne Thienpondt
20/01/2023

"""

#!/usr/bin/python3  
import sys, os
import pathlib

# Personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.input.write_iniFileForInputs import write_iniFileForInputs 
from stellapy.data.geometry.write_h5FileForKperp2 import write_h5FileForKperp2
from stellapy.data.fluxes.write_h5FileForFluxes3D import write_h5FileForFluxes3D
from stellapy.data.fluxes.write_h5FileForFluxes4D import write_h5FileForFluxes4D
from stellapy.data.geometry.write_h5FileForGeometry import write_h5FileForGeometry 
from stellapy.data.moments.write_h5FileForMoments2D import write_h5FileForMoments2D
from stellapy.data.moments.write_h5FileForMoments3D import write_h5FileForMoments3D
from stellapy.data.moments.write_h5FileForMoments4D import write_h5FileForMoments4D
from stellapy.data.moments.write_h5FileForMoments5D import write_h5FileForMoments5D
from stellapy.data.fluxes.write_txtFileForFluxesVsZ import write_txtFileForFluxesVsZ
from stellapy.data.moments.write_h5FileForPhaseShifts import write_h5FileForPhaseShifts
from stellapy.data.omega.write_txtFileForOmegaVsTime import write_txtFileForOmegaVsTime 
from stellapy.data.dimensions.write_h5FileForDimensions import write_h5FileForDimensions
from stellapy.data.potential.write_h5FileForPotential3D import write_h5FileForPotential3D
from stellapy.data.potential.write_h5FileForPotential4D import write_h5FileForPotential4D
from stellapy.data.potential.write_h5FileForPotential5D import write_h5FileForPotential5D
from stellapy.data.output.write_txtFileForOutputVsTime import write_txtFileForOutputVsTime
from stellapy.data.fluxes.write_txtFileForFluxesVsTime import write_txtFileForFluxesVsTime
from stellapy.data.potential.write_txtFileForDPhiZVsTime import write_txtFileForDPhiZVsTime
from stellapy.data.potential.write_txtFileForPotentialVsZ import write_txtFileForPotentialVsZ
from stellapy.data.input.write_listOfMatchingInputFiles import write_listOfMatchingInputFiles
from stellapy.data.distribution.write_h5FileForDistribution3D import write_h5FileForDistribution3D
from stellapy.data.distribution.write_h5FileForDistribution4D import write_h5FileForDistribution4D
from stellapy.data.potential.write_txtFileForPotentialVsTime import write_txtFileForPotentialVsTime 
from stellapy.data.saturated.write_netcdfFileOverSaturatedPhase import write_netcdfFileOverSaturatedPhase
from stellapy.data.distribution.write_txtFileForDistributionVsTime import write_txtFileForDistributionVsTime
from stellapy.data.distribution.write_txtFileForDistributionVsMuOrVpaOrZ import write_txtFileForDistributionVsMuOrVpaOrZ
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#             Write the data files needed for the stellapy package             #
#===============================================================================
 
def write_stellapyDataFiles(folder, dt=None, specific=None, skip=None, dimensions="small"):
    print()
         
    # Write all files if no specific file was selected
    if specific==None:
        
        # Only write 1D files
        if dimensions in ["1D", "2D", "all", "small", "quick", "restart"]:
         
            # Write data files for the input, geometry and slurm files
            write_listOfMatchingInputFiles(folder)
            write_iniFileForInputs(folder)       
            write_h5FileForGeometry(folder)      
             
            # Save the dimensions to an h5 file
            write_h5FileForDimensions(folder)    
        
        # Only write 2D files
        if dimensions in ["2D", "all", "small", "quick", "restart"]: 
            
            # Reduce the size of the output, omega and fluxes file 
            write_txtFileForOutputVsTime(folder)       
            write_txtFileForOmegaVsTime(folder)  
            write_txtFileForFluxesVsTime(folder)   
            write_txtFileForFluxesVsZ(folder) 
            
            # Write data files for the potential
            write_txtFileForPotentialVsTime(folder) 
            write_txtFileForDPhiZVsTime(folder)     
            write_txtFileForPotentialVsZ(folder)  
        
        # Only write 2D files
        if dimensions in ["2D", "all", "small", "quick"]: 
             
            # Write data files for the distribution
            write_txtFileForDistributionVsMuOrVpaOrZ(folder)
            write_txtFileForDistributionVsTime(folder)       
            
            # Write data files for the moments
            write_h5FileForMoments2D(folder)
        
        # Only write 3D files
        if dimensions in ["3D", "all", "small", "quick"]: 
            write_h5FileForFluxes3D(folder, automatic=True) 
            write_h5FileForPotential3D(folder, automatic=True)  
            write_h5FileForDistribution3D(folder)
            write_h5FileForMoments3D(folder, automatic=True)  
            write_netcdfFileOverSaturatedPhase(folder)  
        if dimensions in ["3D", "all", "small"]: 
            write_h5FileForPhaseShifts(folder)
        
        # Only write 4D files
        if dimensions in ["4D", "all", "medium"]: 
            write_h5FileForDistribution4D(folder)
            write_h5FileForPotential4D(folder)   
            write_h5FileForFluxes4D(folder)
            write_h5FileForMoments4D(folder)    
            
        # Only write 5D files
        if dimensions in ["5D", "all", "big", "large"]:
            write_h5FileForPotential5D(folder)   
            write_h5FileForMoments5D(folder)  
             
    # Write specific files
    if specific!=None:
        
        # Get the argiments
        if dt==None: args = {'folder' : folder}
        if dt!=None: args = {'folder' : folder, 'dt' : dt}
        if skip!=None: args = {'skip' : skip, **args}
        
        # Write the data files
        if specific=="ini":         write_listOfMatchingInputFiles(folder); write_iniFileForInputs(**args)
        elif specific=="geo":       write_h5FileForGeometry(**args)
        elif specific=="dim":       write_h5FileForDimensions(**args)
        elif specific=="omega":     write_txtFileForOmegaVsTime(**args)
        elif specific=="jump":      write_txtFileForDPhiZVsTime(**args)
        elif specific=="mom2D":     write_h5FileForMoments2D(**args)
        elif specific=="mom3D":     write_h5FileForMoments3D(**args)
        elif specific=="mom4D":     write_h5FileForMoments4D(**args)
        elif specific=="mom5D":     write_h5FileForMoments5D(**args)
        elif specific=="QL":        write_txtFileForQuasiLinearFluxes(**args)
        elif specific=="flux":      write_txtFileForFluxesVsTime(**args)
        elif specific=="flux3D":    write_h5FileForFluxes3D(**args)
        elif specific=="flux4D":    write_h5FileForFluxes4D(**args)
        elif specific=="pot":       write_txtFileForPotentialVsTime(**args)
        elif specific=="pot2D":     write_txtFileForPotentialVsZ(**args)
        elif specific=="pot3D":     write_h5FileForPotential3D(**args)
        elif specific=="pot4D":     write_h5FileForPotential4D(**args)
        elif specific=="pot5D":     write_h5FileForPotential5D(**args)
        elif specific=="g":         write_txtFileForDistributionVsTime(**args)
        elif specific=="g2D":       write_txtFileForDistributionVsMuOrVpaOrZ(**args)
        elif specific=="g3D":       write_h5FileForDistribution3D(**args)
        elif specific=="g4D":       write_h5FileForDistribution4D(**args)
        elif specific=="kperp":     write_h5FileForKperp2(folder=folder) 
        elif specific=="kperp2":    write_h5FileForKperp2(folder=folder) 
        elif specific=="sat":       write_netcdfFileOverSaturatedPhase(**args) 
        elif specific=="phases":    write_h5FileForPhaseShifts(**args)
        elif specific=="freq":      write_h5FileForPotential3D(folder, dt=0.1)
        else:                       print("The specific quantity '"+str(specific)+"' does not exist. ")
     
    print()       
    return

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    bash = Bash(write_stellapyDataFiles, __doc__)  
    bash.add_option('specific', 'str', 's', '', '{ini, geo, dim, omega, jump, QL, flux, flux3D, flux4D, pot, \npot2D, pot3D, pot4D, pot5D, g, g2D, g3D, g4D\nmom3D, mom4D, mom5D, sat, phases, freq}') 
    bash.add_option('dt', 'float', 't', '', 'Time step dt to write data.')  
    bash.add_option('skip', 'int', 'i', '', 'Skip every <skip> time steps.')  
    bash.add_toggle('dimensions', 'restart', '', '', 'Only write quantities needed to restart unconverged modess.')  
    bash.add_toggle('dimensions', '1D', '', '', 'Only write 1D quantities.')  
    bash.add_toggle('dimensions', '2D', '', '', 'Only write 2D quantities.')  
    bash.add_toggle('dimensions', '3D', '', '', 'Only write 3D quantities.')  
    bash.add_toggle('dimensions', '4D', '', '', 'Only write 4D quantities.')  
    bash.add_toggle('dimensions', '5D', '', '', 'Only write 5D quantities.')  
    bash.add_toggle('dimensions', 'all', '', '', 'Write all quantities.')  
    bash.add_toggle('dimensions', 'small', '', '', 'Only write 1D, 2D and 3D quantities.') 
    bash.add_toggle('dimensions', 'medium', '', '', 'Only write 1D, 2D, 3D and 4D quantities.') 
    bash.add_toggle('dimensions', 'big', '', '', 'Only write 1D, 2D, 3D, 4D and 5D quantities.') 
    bash.add_toggle('dimensions', 'large', '', '', 'Only write 1D, 2D, 3D, 4D and 5D quantities.') 
    bash.add_toggle('dimensions', 'quick', '', '', 'Only write 1D, 2D and 3D quantities, without the phase shifts.') 
    write_stellapyDataFiles(**bash.get_arguments())   




