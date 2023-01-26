
import os, sys 

def restart_program():
    """Restarts the current program.
    Note: this function does not return. Any cleanup action (like
    saving data) must be done before calling this function."""
    
    # Restart the program 
    python = sys.executable
    os.execl(python, python, *sys.argv)


