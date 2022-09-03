 
import configparser  

#===============================================================================
#                               READ REMOVED MODE                              #
#===============================================================================    

def get_removedModes(self):
    
    # If we have a mode, remember whether we need to remove this mode
    if self.object=="Mode":
        self.removeMode = False
        if self.kx in self.simulation.lineardata.removedKxModes: self.removeMode = True
        if self.ky in self.simulation.lineardata.removedKyModes: self.removeMode = True 
        return 
    
    # If we have a simulation, read the modes that need to be removed
    if self.object=="Simulation":
    
        # Read the removed modes from the "removedModes" file 
        file = configparser.ConfigParser() 
        path = self.path.folder / "removedModes.ini"
        file.read(path)
        
        # Initiate the file if it doesn't exist
        if "Removed modes" not in file: 
            file["Simulated modes"] = {"kx" : "\n", "ky" : "\n"}
            file["Removed modes"] = {"kx" : "\n", "ky" : "\n"}
    
        # Update the present modes
        kx_modes = sorted(list(set([mode.kx for mode in self.modes]))); kx_modes = [str(kx) for kx in kx_modes]
        ky_modes = sorted(list(set([mode.ky for mode in self.modes]))); ky_modes = [str(ky) for ky in ky_modes] 
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


