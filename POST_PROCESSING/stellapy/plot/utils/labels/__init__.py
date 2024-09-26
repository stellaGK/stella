 
# Global variables
standardHeaders = {}
standardLabels = {}
standardTitles = {}
standardNames = {}
standardNamesScientific = {}
standardParameters = {}
standardParametersInOrder = {} 

from stellapy import enforce_one_function_one_file
enforce_one_function_one_file(__path__, globals(), locals())
