  
from stellapy.data.input.read_inputFile import read_inputParameters

#===============================================================================
#                          GET THE INPUT PARAMETERS
#===============================================================================

def get_inputParameters(self):   
    self.inputParameters = self.input_parameters = read_inputParameters(self.path.input)   
    return