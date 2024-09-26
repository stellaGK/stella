
 
def printv(message):
    
    # Get the indentation for the text printed to the command prompt
    from stellapy.utils.decorators.verbose import indent
    from stellapy.utils.decorators.verbose import wrapper
    
    # Only print if verbose==True in the configuration file.  
    if wrapper==True:
        print(indent, message)
        return
    else:
        return

