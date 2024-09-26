
import sys

def exit_program(reason, function, code_line):
    ''' Exit the program and print the reason to the command prompt.
    
    Example
    -------
    exit_program(reason, function_name, sys._getframe().f_lineno) 
     ''' 
 
    reasons = reason.split("\n") 
    
    width = 80
    print()
    print("", "".center(width,"="))
    print("", "EXIT PROGRAM".center(width," ")); 
    print("", "".center(width,"="))
    for i in range(len(reasons)):
        print(" ", reasons[i]); 
    print("", "".center(width,"="))
    print("", " script:    ", function.__module__);
    print("", " function:  ", function.__name__);
    print("", " line:      ", code_line);
    print("", "".center(width,"="))
    print()
    sys.exit(2)
