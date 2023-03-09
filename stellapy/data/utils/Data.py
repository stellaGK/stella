
import copy
import numpy as np

#===============================================================================
#                     STANDARD CLASS FOR EACH PIECE OF DATA                    #
#===============================================================================

class Data():
    
    def __init__(self, dimensions, dim0=None, dim1=None, dim2=None, dim3=None, dim4=None, dim5=None):
        
        # Save the dimensions
        self.quantity = dimensions[0]
        self.dimensions = dimensions[1:]
        
        # Save the data as attributes
        if len(self.dimensions)>=0: setattr(self, self.quantity, np.array(dim0))
        if len(self.dimensions)>=1: setattr(self, self.dimensions[0], np.array(dim1))
        if len(self.dimensions)>=2: setattr(self, self.dimensions[1], np.array(dim2))
        if len(self.dimensions)>=3: setattr(self, self.dimensions[2], np.array(dim3))
        if len(self.dimensions)>=4: setattr(self, self.dimensions[3], np.array(dim4)) 
        if len(self.dimensions)>=5: setattr(self, self.dimensions[4], np.array(dim5))  
        return
    
    #--------------------
    def get_labels(self):
        if len(self.dimensions)==2:
            from stellapy.plot.utils.labels import standardLabels
            xlabel = standardLabels[self.dimensions[0]]
            ylabel = standardLabels[self.quantity]
            return xlabel, ylabel
        
    #--------------------
    def get_xyData(self): 
        if len(self.dimensions)==1:
            x = getattr(self, self.dimensions[0])
            y = getattr(self, self.quantity)
            return x, y
        if len(self.dimensions)==2 and "s" in self.dimensions: 
            dimensions = copy.deepcopy(self.dimensions); dimensions.remove("s") 
            x = getattr(self, dimensions[0])
            y = getattr(self, self.quantity)
            return x, y
        print("SOMETHING IS WRONG: stellapy.data.get_dataForPlots.get_xyData")
  

