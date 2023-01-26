
 
def printWhichFileWeRead(message, verbose=False):
    # Only print if verbose==True 
    if verbose==True:
        from stellapy.utils.decorators.verbose import indent
        print(indent, message)
        return
    else:
        return

