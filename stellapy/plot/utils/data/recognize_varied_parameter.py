#!/usr/bin/python3  
import sys, os
import numpy as np

# Stellapy package
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0]) 
from stellapy.data.check.check_differencesInStellaInputs import check_differencesInStellaInputs
from stellapy.plot.utils.labels.standardParameters import standardParameters

#===============================================================================
#                       RECOGNIZE THE SCANNED PARAMETER                        #
#===============================================================================

def recognize_varied_parameter(folder, parameter=None):
    """ Try to recognize the changed parameter. """
    
    # Find the varied parameter
    if parameter==None:
        differences = check_differencesInStellaInputs(folder, verbose=False)
        keys = [i for i in list(differences.keys()) if i not in ["tend", "ginit_option", "delt_option"]]
        if len(keys)!=0:  
            lengths = [len(differences[key]) for key in keys]
            parameter = keys[np.argmax(lengths)]
            
    # Test whether it exists in <standardParameters> 
    try: standardParameters[parameter]["knob"]
    except: parameter = "tiprim"
    return parameter