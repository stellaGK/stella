
def get_bashArguments(args):
    """ When redirecting the bash arguments to the next script, make sure the
    arguments are passed on correctly. """
    
    # Cut off the first argument which is the execution folder
    args = args[1:]
    
    # The arguments are read in as a list, convert them back to a string
    args = ' "'+'" "'.join(args)+'"'
    
    # Return the string of bash arguments
    return args
    
     
    