#===============================================================================
#          Load the stellapy package in the python interavtive prompt          #
#===============================================================================
# The python interactive prompt is opened by typing in the command prompt:     #
#      >>> python3                                                             #
# This script is used to automatically load the stellapy package.              # 
#===============================================================================

# Print this line as a check to know whether the stellapy package is loaded.
print("Loaded the stellapy package.")

# Python imports work by searching the directories listed in sys.path.
# Manually add the path where the modules for stella are stored.
import sys, os 
sys.path.append(os.environ['STELLA'])  

# Import the modules for stella to the python environment
# This way of reading allows to use help(stella_dirs) in python
import stellapy


        
