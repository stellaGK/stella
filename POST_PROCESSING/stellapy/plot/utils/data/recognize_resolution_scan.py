#!/usr/bin/python3   
import os, pathlib 

#===============================================================================
#                       RECOGNIZE THE SCANNED PARAMETER                        #
#===============================================================================

def recognize_resolution_scan(folder, resolutionScan=None):
    """ Try to recognize if we are doing a resolution scan. """
    if resolutionScan==None: 
        resolutionScan = False
        subfolders = [pathlib.Path(x[0]) for x in os.walk(folder) if pathlib.Path(x[0])!=folder]
        if len(subfolders)>0:
            for subfolder in subfolders:
                if "double" in subfolder.name: resolutionScan = True
    return resolutionScan