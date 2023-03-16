"""
#===============================================================================
#                          LIST ALL THE STELLAPY COMMANDS                      #
#===============================================================================

List all the stellapy commands. 

"""

#!/usr/bin/python3  
import sys, os
import pathlib

# Personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                          LIST ALL THE STELLAPY COMMANDS                      #
#===============================================================================

def list_commands(print_only=None, width=90): 
    
    # Start 
    if print_only==None:
        print("\n"+"".center(width,"="))
        print("\033[1m"+"STELLAPY 6.0".center(width," ")+"\033[0m")
        print("".center(width,"="))   
        print("""
        In order to use the stellapy scripts and functions one can call the functions 
        directly from the python3 interactive prompt, or one can use the bash commands 
        defined in stellapy/source.sh directly from the command prompt. A list of possible
        stellapy functions that work as bash commands is shown through the command:
            >> stellapy
        
        These commands can be used in the same way as bash commands, a list of options 
        is shown for each commands through:
            >> command --help
        
        For example the <plot/linear/gamma_vs_time.py> script is called through:
            plot_gamma_vs_time --unstable --kymin 1 --kymax 2
            
        However, for more complex scripts, with many arguments, this is not feasible.
        If one wants to frequently update the stellapy package, it is not recommended to 
        make edits to the scripts of stellapy directly, since they will be lost when
        one updates the package. Therefore, if one wants to save the arguments that are
        passed onto the stellapy scripts, it is recommended to create a duplicate of the
        script in stellapy/User, and relink the original bash command in source.sh to the 
        duplicate script in User, for example the source scripts contain:
            stellapy/source.sh: alias plot_gamma_vs_time="python3 $STELLAPY/plot/linear/gamma_vs_time.py" 
            stellapy/User/source.sh: alias plot_gamma_vs_time="python3 $STELLAPY/User/plot/linear/gamma_vs_time.py" 
        
        It is evident that in the .alias file, or in the .bashrc file we have:
            source $STELLAPY/source.sh
            source $STELLAPY/User/source.sh
        
        This way one can freely manipulate the general scripts, overwrite commands, and 
        save the used variables. 
        
        Developed by Hanne Thienpondt. 
        09/03/2022
        """)
    
    # Start 
    print("\n"+"".center(width,"="))
    print("\033[1m"+"STELLAPY BASH COMMANDS".center(width," ")+"\033[0m")
    print("".center(width,"="), "\n") 
    print("  Overview of the stellapy commands which work like Bash commands. ")
    print("  Developed by Hanne Thienpondt. ")
    print("  09/03/2022")
    print()
    print()
    
    # GUI (Graphical User Interface) 
    if print_only==None or print_only=="GUI":
        print("\033[1m GUI (Graphical User Interface) \033[0m")
        print("\033[1m ----------------------------- \033[0m", "\n")
        print(" THE GUI IS NOT WORKING RIGHT NOW!")
        print("    ", " >> stellaplotter")
        print("    ", " >> stellaplotter_linear")
        print("    ", " >> stellaplotter_nonlinear")
        print()
    
    # Memory management
    if print_only==None or print_only=="memory":
        print("\033[1m MEMORY MANAGEMENT \033[0m")
        print("\033[1m ----------------- \033[0m", "\n")
        print("    ", " >> reduce_sizeNetcdf")  
        print("    ", " >> replace_netcdfFile")  
        print() 
    
    # Data checking
    if print_only==None or print_only=="check":
        print("\033[1m DATA CHECKING \033[0m")
        print("\033[1m ------------- \033[0m", "\n")
        print(" THESE COMMANDS NEED TO BE CLEANED!")
        print("    ", " >> check_cputime")  
        print("    ", " >> check_timetraces")  
        print("    ", " >> check_netcdf")  
        print("    ", " >> check_dimensions")  
        print("    ", " >> check_netcdfVariables")  
        print("    ", " >> check_differencesInStellaInputs")  
        print("    ", " >> print_geometry") 
        print() 
    
    # Data processing
    if print_only==None or print_only=="write":
        print("\033[1m DATA PROCESSING \033[0m")
        print("\033[1m --------------- \033[0m", "\n")
        print("    ", " >> write_stellapyDataFiles") 
        print("          ", " >> write_stellapyDataFiles -s ini") 
        print("          ", " >> write_stellapyDataFiles -s pot3D -t 1") 
        print("          ", " >> write_stellapyDataFiles -s phases --skip 10")
        print("          ", " >> write_stellapyDataFiles --restart") 
        print("          ", " >> write_stellapyDataFiles --2D") 
        print() 
    
    # Plot linear simulations
    if print_only==None or print_only=="linear" or print_only=="plot":
        print("\033[1m PLOT LINEAR SIMULATIONS \033[0m")
        print("\033[1m ---------------------- \033[0m", "\n")
        print("   ", "\033[1m TIME EVOLUTION \033[0m")
        print("   ", "\033[1m -------------- \033[0m") 
        print("      ", " >> plot_dphiz_vs_time")  
        print("      ", " >> plot_gamma_vs_time")   
        print("      ", " >> plot_omega_vs_time")   
        print("      ", " >> plot_potential_vs_time")     
        print() 
        print("   ", "\033[1m FOURIER SPACE \033[0m")
        print("   ", "\033[1m ------------- \033[0m")
        print("      ", " >> plot_gamma_vs_kx     >> plot_gamma_vs_ky")   
        print("      ", " >> plot_omega_vs_kx     >> plot_omega_vs_ky")  
        print() 
        print("   ", "\033[1m VELOCITY SPACE \033[0m")
        print("   ", "\033[1m -------------- \033[0m")
        print("      ", " >> plot_distribution_vs_vpa") 
        print("      ", " >> plot_distribution_vs_mu")  
        print("      ", " >> plot_distribution_vs_velocity")    
        print()
        print("   ", "\033[1m FIELD LINE \033[0m")
        print("   ", "\033[1m ---------- \033[0m")
        print("      ", " >> plot_geometry_vs_z")  
        print("      ", " >> plot_potential_vs_z")  
        print("      ", " >> plot_distribution_vs_z") 
        print() 
        print("   ", "\033[1m PARAMETER INFLUENCE \033[0m")
        print("   ", "\033[1m ------------------- \033[0m")
        print("      ", " >> plot_gamma_vs_parameter")  
        print() 
        print("   ", "\033[1m COMPLEX ANALYSIS \033[0m")
        print("   ", "\033[1m ---------------- \033[0m")
        print("      ", " >> plot_gamma      a.k.a     plot_spectra  ")   
        print("      ", " >> plot_scans      a.k.a     plot_parameter_influence")   
        print() 
    
    # Plot nonlinear simulations
    if print_only==None or print_only=="nonlinear" or print_only=="plot":
        print("\033[1m PLOT NONLINEAR SIMULATIONS \033[0m")
        print("\033[1m -------------------------- \033[0m", "\n")
        print("   ", "\033[1m TIME EVOLUTION \033[0m")
        print("   ", "\033[1m -------------- \033[0m")   
        print("      ", " >> plot_qflux_vs_time")  
        print("      ", " >> plot_pflux_vs_time")     
        print("      ", " >> plot_vflux_vs_time")  
        print("      ", " >> plot_phi2_vs_time")    
        print("      ", " >> plot_potential_vs_time")           
        print() 
        print("   ", "\033[1m SPECTRA \033[0m")
        print("   ", "\033[1m ------- \033[0m") 
        print("      ", " >> plot_qflux_spectra        >> plot_qflux_vs_kx        >> plot_qflux_vs_ky")   
        print("      ", " >> plot_pflux_spectra        >> plot_pflux_vs_kx        >> plot_pflux_vs_ky") 
        print("      ", " >> plot_vflux_spectra        >> plot_vflux_vs_kx        >> plot_vflux_vs_ky") 
        print("      ", " >> plot_phi2_spectra         >> plot_phi2_vs_kx         >> plot_phi2_vs_ky")  
        print("      ", " >> plot_potential_spectra    >> plot_potential_vs_kx    >> plot_potential_vs_ky") 
        print() 
        print("   ", "\033[1m FIELD LINE \033[0m")
        print("   ", "\033[1m ---------- \033[0m")
        print("      ", " >> plot_geometry_vs_z")  
        print("      ", " >> plot_potential_vs_z")  
        print("      ", " >> plot_moment_vs_z")  
        print() 
        print("   ", "\033[1m PARAMETER \033[0m")
        print("   ", "\033[1m --------- \033[0m")
        print("      ", " >> plot_qflux_vs_parameter")  
        print("      ", " >> plot_pflux_vs_parameter")  
        print("      ", " >> plot_vflux_vs_parameter")   
        print() 
      
    return

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash interface
    bash = Bash(list_commands, __doc__) 
    
    # Data toggles to print specific parts of the commands
    bash.add_toggle('print_only', 'plot', '', '', 'Show the plotting commands.')  
    bash.add_toggle('print_only', 'linear', '', '', 'Show only the linear plotting commands.')  
    bash.add_toggle('print_only', 'nonlinear', '', '', 'Show only the nonlinear plotting commands.')  
    bash.add_toggle('print_only', 'GUI', '', '', 'Show only the GUI commands.')  
    bash.add_toggle('print_only', 'memory', '', '', 'Show only the memory commands.')  
    bash.add_toggle('print_only', 'check', '', '', 'Show only the check commands.')   
    bash.add_toggle('print_only', 'write', '', '', 'Show only the data writing commands.')  
    
    # Get the arguments and launch the script
    args = bash.get_arguments()
    del args['folder']
    list_commands(**args)   
    
