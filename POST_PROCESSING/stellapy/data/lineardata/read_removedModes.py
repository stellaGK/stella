 
import configparser  

#===============================================================================
#                               READ REMOVED MODE                              #
#===============================================================================    

def read_removedModes(self): 
    
    # Read the removed modes from the "removedModes" file 
    file = configparser.ConfigParser() 
    path = self.path.folder / "removedModes.ini"
    file.read(path)
    
    # Initiate the file if it doesn't exist
    if "Removed modes" not in file: 
        file["Simulated modes"] = {"kx" : "\n", "ky" : "\n"}
        file["Removed modes"] = {"kx" : "\n", "ky" : "\n"}

    # Update the present modes
    kx_modes = [str(kx) for kx in self.vec.kx]
    ky_modes = [str(ky) for ky in self.vec.ky] 
    file["Simulated modes"]["kx"] = "\n" + "\n".join(kx_modes) + "\n"
    file["Simulated modes"]["ky"] = "\n" + "\n".join(ky_modes) + "\n"
    
    # Read the removed modes
    self.removedKxModes = file["Removed modes"]["kx"]
    self.removedKyModes = file["Removed modes"]["ky"]
    self.removedKxModes = self.removedKxModes.split("\n")  
    self.removedKyModes = self.removedKyModes.split("\n")   
    self.removedKxModes = [kx for kx in self.removedKxModes if kx!=''] 
    self.removedKyModes = [ky for ky in self.removedKyModes if ky!=''] 
    self.removedKxModes = [ float(kx) for kx in self.removedKxModes ] 
    self.removedKyModes = [ float(ky) for ky in self.removedKyModes ] 

    # Write the "removedModes" file 
    file.write(open(path, 'w'))
    return    


