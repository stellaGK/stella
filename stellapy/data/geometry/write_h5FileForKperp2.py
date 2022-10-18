  
import os, h5py
from stellapy.utils.decorators.verbose import noverbose      
from stellapy.utils.files.get_filesInFolder import get_filesInFolder
from stellapy.utils.files import remove_simulationsWithoutOutncFile
from stellapy.data.output.read_outputFile import read_netcdfVariables, read_outputFile 
 
#===============================================================================
#                          WRITE THE KPERP2 FILE
#===============================================================================
 
@noverbose
def write_h5FileForKperp2(folder):  
    """ We do not put this in the geometry file, since it is (kx*ky) times bigger 
    than any other geometry quantity. Only write it when needed. """

    # Get the input_files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = remove_simulationsWithoutOutncFile(input_files)   
    for input_file in input_files:
        
        # Processing status
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
        
        # Kperp2 file 
        kperp2_file = input_file.with_suffix(".kperp2") 
        if not os.path.isfile(kperp2_file): 
        
            # Read kperp2 
            netcdf_file = read_outputFile(input_file.with_suffix(".out.nc"))  
            kperp2  = read_netcdfVariables('kperp2', netcdf_file)  
             
            # Write the h5 file
            with h5py.File(kperp2_file, 'w') as h5_file: 
                h5_file.create_dataset("kperp2", data=kperp2) 
                
        # Notify
            print(status+"   ---> The kperp2 file is saved as " +  kperp2_file.parent.name+"/"+kperp2_file.name)  
        else:
            print(status+"The kperp2 file already exists:", kperp2_file.parent.name+"/"+kperp2_file.name) 

#--------------------
def read_kperp2(input_file):  
    kperp2_file = input_file.with_suffix(".kperp2") 
    netcdf_file = input_file.with_suffix(".out.nc") 
    if not os.path.isfile(kperp2_file) and os.path.isfile(netcdf_file): 
        netcdf_file = read_outputFile(input_file.with_suffix(".out.nc"))  
        kperp2  = read_netcdfVariables('kperp2', netcdf_file)  
        return kperp2
    elif os.path.isfile(kperp2_file): 
        with h5py.File(kperp2_file, 'r') as f:  
            kperp2 = f["kperp2"][()]    
        return kperp2
    else:
        print("The netcdf file nor the kperp file are present, can not read kperp(z,kx,ky).")
        return None
    