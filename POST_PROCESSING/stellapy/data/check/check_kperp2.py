
#!/usr/bin/python3  
import sys, os   

# Stellapy package
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])  
from stellapy.utils.files.keep_simulationsWithNetcdfFile import keep_simulationsWithNetcdfFile
from stellapy.data.output.read_outputFile import read_netcdfVariables, read_outputFile 
from stellapy.data.input.read_inputFile import read_modeFromInputFile
from stellapy.utils.files.get_filesInFolder import get_filesInFolder 
from stellapy.utils.commandprompt.bash import Bash 

#===============================================================================
#                                 CHECK KPERP2                                 #
#===============================================================================

def check_kperp2(folder):
    
    # Get the input_files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = keep_simulationsWithNetcdfFile(input_files)   
    for input_file in input_files:
    
        # Read kperp2(z,kx,ky) from the netcdf file
        netcdf_file = read_outputFile(input_file.with_suffix(".out.nc"))  
        kperp2 = read_netcdfVariables('kperp2', netcdf_file) 
        netcdf_file.close()
        
        # Print kperp2
        kx, ky = read_modeFromInputFile(input_file)
        print(f"\n|k_perp(z)|^2 for (kx,ky) = ({kx}, {ky}):")
        print(kperp2[:,0,0])
        
    return    
                
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__": 
    bash = Bash(check_kperp2, __doc__)    
    check_kperp2(**bash.get_arguments())  