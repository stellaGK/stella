#!/usr/bin/python3  
import sys, os
import pathlib
import numpy as np

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.data.input.get_differencesInStellaInputs import get_differencesInStellaInputs
from stellapy.plot.utils.labels.standardParameters import standardParameters

#===============================================================================
#                       RECOGNIZE THE SCANNED PARAMETER                        #
#===============================================================================

def recognize_varied_parameter(folder, parameter=None):
    """ Try to recognize the changed parameter. """
    
    # Keys that are allowed to differ
    remove_keys = ["aky_min", "aky_max", "akx_min", "akx_max", "nakx", "naky", 
        "tend", "ginit_option", "delt_option", "nsave", "nwrite", "delt", "nmu"]
    
    # Find the varied parameter
    if parameter==None: 
        differences = get_differencesInStellaInputs(folder) 
        keys = [i for i in list(differences.keys()) if i not in remove_keys]
        for key in remove_keys:
            if key in differences:
                del differences[key] 
        if len(keys)!=0:  
            lengths = [len(differences[key]) for key in keys]
            parameter = keys[np.argmax(lengths)]

    # Test whether it exists in <standardParameters> 
    try: standardParameters[parameter]["knob"]
    except: parameter = "tiprim"
    return parameter