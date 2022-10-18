
import functools 
from .print_decoratorOnCommandPrompt import print_decoratorOnCommandPrompt

# Keep track of the indentation of the subfunctions
indent = ""
wrapper = True 

def verbose(function):
    @functools.wraps(function)                                                  # Make sure function remembers its __name__ and so on
    def wrapper_verbose(*args, **kwargs):
        try:
            global indent
            module = function.__module__.split('.')
            name = module[0] + '.' + module[1] + '.' + function.__name__            # include stellapy.*in the header @unusedvariable
            name = module[1] + '.' + function.__name__                              # don't include stellapy.* in the header
            
            # Add a header to the function with the function name    
            if wrapper==True: 
                indent = print_decoratorOnCommandPrompt(name)                                           
                return_arg = function(*args, **kwargs)
                indent = print_decoratorOnCommandPrompt("end")
                return return_arg
            else:
                return function(*args, **kwargs)
        except: function(*args, **kwargs)

    return wrapper_verbose


def noverbose(function):
    @functools.wraps(function)
    def wrapper_verbose(*args, **kwargs):
        global wrapper 
        wrapper = False 
        return function(*args, **kwargs)
    return wrapper_verbose


def turn_offVerboseWrapper():
    global wrapper 
    wrapper = False 

