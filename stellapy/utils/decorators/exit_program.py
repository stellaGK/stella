
import sys

def exit_program(reason, function, code_line):
    ''' Exit the program and print the reason to the command prompt.
    
    Example
    -------
    exit_program(reason, function_name, sys._getframe().f_lineno) 
     '''

    from stellapy.utils.decorators.verbose import indent
 
    reasons = reason.split("\n") 
    
    width = 60
    print()
    print(indent, "".center(width,"="))
    print(indent, "EXIT PROGRAM".center(width," ")); 
    print(indent, "".center(width,"="))
    for i in range(len(reasons)):
        print(indent, "", reasons[i]); 
    print(indent, "     --->  ", function.__module__);
    print(indent, "     --->  ", function.__name__);
    print(indent, "     --->   line", code_line);
    print(indent, "".center(width,"="))
    print()
    sys.exit(2)
