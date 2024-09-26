
# Global variables
CONFIG = {} 
NAME_CONFIGURATIONFILE = "ShouldBeSetALongTheWay"

# Make every function avaliable as package.func instead of package.mod.func
from stellapy import enforce_one_function_one_file
enforce_one_function_one_file(__path__, globals(), locals())
 
