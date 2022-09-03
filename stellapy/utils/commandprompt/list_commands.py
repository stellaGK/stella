"""
#===============================================================================
#                          LIST ALL THE STELLAPY COMMANDS                      #
#===============================================================================

List all the stellapy commands. 

"""

#!/usr/bin/python3  
import sys, os

# Personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])  
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                          LIST ALL THE STELLAPY COMMANDS                      #
#===============================================================================

def list_commands(width=90):
    
    # Start 
    print("\n"+"".center(width,"="))
    print("\033[1m"+"STELLAPY".center(width," ")+"\033[0m")
    print("".center(width,"="))   
    print("""
    In order to use the stellapy scripts and functions one can call the functions 
    directly from the python3 interactive prompt, or one can use the bash commands 
    defined in stellapy/source.sh directly from the command prompt. A list of possible
    stellapy functions that work as bash commands is shown through the command:
        >> stellapy
    
    These commands can be used in the same way as bash commands, a list of options 
    is shown for each commands by performing:
        >> command -h
    
    For example the <launchLinearKyScan()> script can easily be used directly:
        launchLinearKyScan -t 50 -d 0.05 -i scan_upto_5
        
    However, for more complex scripts, with many arguments, this is not feasible.
    If one wants to frequently update the stellapy package, it is not recommended to 
    make edits to the scripts of stellapy directly, since they will be lost when
    one updates the package. Therefore, if one wants to save the arguments that are
    passed onto the stellapy scripts, it is recommended to create a duplicate of the
    script in stellapy/User, and relink the original bash command in source.sh to the 
    duplicate script in User, for example the source scripts contain:
        stellapy/source.sh: alias launch_linearGradientScan="python3 $STELLAPY/launch/massLinearLaunchingScripts/launch_linearGradientScan.py" 
        stellapy/User/source.sh: alias launch_linearGradientScan="python3 $STELLAPY/User/launch/massLinearLaunchingScripts/launch_linearGradientScan.py" 
    
    It is evident that in the .alias file, or in the .bashrc file we have:
        source $STELLAPY/source.sh
        source $STELLAPY/User/source.sh
    
    This way one can freely manipulate the general scripts, overwrite commands, and 
    save the used variables. 
    
    Developed by Hanne Thienpondt. 
    30/08/2022
    """)
    
    # Start 
    print("\n"+"".center(width,"="))
    print("\033[1m"+"STELLAPY BASH COMMANDS".center(width," ")+"\033[0m")
    print("".center(width,"="), "\n") 
    print("  Overview of the stellapy commands which work like Bash commands. ")
    print("  Developed by Hanne Thienpondt. ")
    print("  23/05/2022")
    print()
    print()
    
    # GUI (Graphical User Interface) 
    print("\033[1m GUI (Graphical User Interface) \033[0m")
    print("\033[1m ----------------------------- \033[0m", "\n")
    print("    ", " >> stellaplotter")
    print("    ", " >> stellaplotter_linear")
    print("    ", " >> stellaplotter_nonlinear")
    print()
    
    # Supercomputer
    print("\033[1m MARCONI SUPERCOMPUTER \033[0m")
    print("\033[1m --------------------- \033[0m", "\n")
    print("    ", " >> check_runs")
    print("    ", " >> remove_stellaFiles")
    print("    ", " >> load_modulesToCompileStella")
    print()
    print("    ", " >> sync_runsFromMarconi      ($MARCONI_RUNS --> $RUNS)")
    print("    ", " >> sync_newRunsToMarconi     ($NEWRUNS --> $MARCONI_RUNS")
    print("    ", " >> sync_stellaToMarconi      ($STELLA --> $MARCONI_STELLA")
    print("    ", " >> sync_stellapyToCommon     ($STELLAPY --> $MARCONI_COMMON_STELLAPY)")
    print("    ", " >> sync_stellapyFromCommon   ($MARCONI_COMMON_STELLAPY --> $STELLAPY")
    print()
    
    # Memory management
    print("\033[1m MEMORY MANAGEMENT \033[0m")
    print("\033[1m ----------------- \033[0m", "\n")
    print("    ", " >> reduce_sizeNetcdf")  
    print("    ", " >> replace_netcdfFile")  
    print() 
    
    # Data checking
    print("\033[1m DATA CHECKING \033[0m")
    print("\033[1m ------------- \033[0m", "\n")
    print("    ", " >> check_netcdf")  
    print("    ", " >> check_dimensions")  
    print("    ", " >> check_netcdfVariables")  
    print("    ", " >> check_differencesInStellaInputs")  
    print() 
    
    # Data processing
    print("\033[1m DATA PROCESSING \033[0m")
    print("\033[1m --------------- \033[0m", "\n")
    print("    ", " >> write_dataFiles") 
    print("          ", " >> write_dataFiles -s ini") 
    print("          ", " >> write_dataFiles -s pot3D -t 1") 
    print("          ", " >> write_dataFiles -s phases --skip 10")
    print() 
    
    # Launching
    print("\033[1m LAUNCH SIMULATIONS \033[0m")
    print("\033[1m ------------------ \033[0m", "\n") 
    print("   ", "\033[1m GENERAL \033[0m")
    print("   ", "\033[1m ------- \033[0m") 
    print("      ", " >> launch -n 4 -h 23") 
    print("      ", " >> restart")  
    print() 
    print("   ", "\033[1m LINEAR \033[0m")
    print("   ", "\033[1m ------- \033[0m") 
    print("      ", " >> launchLinearKyScan -i scan_upto_5 -d 0.1 -t 50 ")  
    print("      ", " >> launchLinearGradientScan")  
    print("      ", " >> launchLinearScansOfImpurities")  
    print("      ", " >> launch_linearMapOfFprimVsTiprim")  
    print() 
    print("   ", "\033[1m NONLINEAR \033[0m")
    print("   ", "\033[1m --------- \033[0m")  
    print("      ", " >> launchNonlinearScansOfImpurities")
    print() 
    print("   ", "\033[1m CONVERGENCE \033[0m")
    print("   ", "\033[1m ----------- \033[0m") 
    print("      ", " >> check_convergenceModes")  
    print("      ", " >> restart_unconvergedModes")  
    print("      ", " >> check_tend (TODO)")  
    print("      ", " >> restart_simulations -t 1000 (TODO)")  
    print()
    
    # Plot linear simulations
    print("\033[1m PLOT LINEAR SIMULATIONS \033[0m")
    print("\033[1m ---------------------- \033[0m", "\n")
    print("   ", "\033[1m TIME EVOLUTION \033[0m")
    print("   ", "\033[1m -------------- \033[0m") 
    print("      ", " >> plot_dphiz_vs_time")  
    print("      ", " >> plot_gamma_vs_time")   
    print("      ", " >> plot_omega_vs_time")     
    print() 
    print("   ", "\033[1m SPECTRA \033[0m")
    print("   ", "\033[1m ------- \033[0m")
    print("      ", " >> plot_gamma_vs_ky")   
    print("      ", " >> plot_omega_vs_ky") 
    print("      ", " >> plot_gamma_vs_kx")   
    print("      ", " >> plot_omega_vs_kx")  
    print() 
    print("   ", "\033[1m FIELD LINE \033[0m")
    print("   ", "\033[1m ---------- \033[0m")
    print("      ", " >> plot_phi_vs_z")  
    print() 
    print("   ", "\033[1m COMPLEX ANALYSIS \033[0m")
    print("   ", "\033[1m ---------------- \033[0m")
    print("      ", " >> plot_gamma")   
    print() 
    
    # Plot nonlinear simulations
    print("\033[1m PLOT NONLINEAR SIMULATIONS \033[0m")
    print("\033[1m -------------------------- \033[0m", "\n")
    print("   ", "\033[1m TIME EVOLUTION \033[0m")
    print("   ", "\033[1m -------------- \033[0m")  
    print("      ", " >> plot_flux_vs_time --qflux")  
    print("      ", " >> plot_qflux_vs_time")  
    print("      ", " >> plot_pflux_vs_time")     
    print("      ", " >> plot_vflux_vs_time")        
    print() 
    print("   ", "\033[1m SPECTRA \033[0m")
    print("   ", "\033[1m ------- \033[0m")
    print("      ", " >> plot_flux_vs_ky")   
    print("      ", " >> plot_qflux_vs_ky")   
    print("      ", " >> plot_pflux_vs_ky") 
    print("      ", " >> plot_vflux_vs_ky --log") 
    print() 
    print("   ", "\033[1m FIELD LINE \033[0m")
    print("   ", "\033[1m ---------- \033[0m")
    print("      ", " >> plot_phi_vs_z")  
    print() 
    print("   ", "\033[1m PARAMETER \033[0m")
    print("   ", "\033[1m --------- \033[0m")
    print("      ", " >> plot_saturatedflux_vs_parameter")  
    print() 
    print("   ", "\033[1m COMPLEX ANALYSIS \033[0m")
    print("   ", "\033[1m ---------------- \033[0m")
    print("      ", " >> plot_saturatedflux -p fprim") 
    print("      ", " >> ...")   
    print() 
    
    print()
    return

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    bash = Bash(list_commands, __doc__)  
    args = bash.get_arguments()
    del args['folder']
    list_commands(**args)   
    