  
from stellapy.data.input.read_inputFile import read_inputParameters

#===============================================================================
#                          GET THE INPUT PARAMETERS
#===============================================================================

def get_inputParameters(self):  
    self.inputParameters = read_inputParameters(self.path.input) 
    #self.inputParameters["parameters"]["teti"] = 1/self.inputParameters["parameters"]["tite"] 
    return