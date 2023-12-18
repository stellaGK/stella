"""
#===============================================================================
#                          LIST ALL THE STELLAPY COMMANDS                      #
#===============================================================================

List all the stellapy commands, including some extra launching/syncing commands.

"""

#!/usr/bin/python3  
import sys, os
import pathlib

# Personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.utils.commandprompt.list_commands import list_commands
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                          LIST ALL THE STELLAPY COMMANDS                      #
#===============================================================================

def list_commands_extra(print_only=None, width=90):
    
    # Call the general script 
    list_commands(print_only, width=width) 
    
    # Supercomputer
    if print_only==None or print_only=="hpc" or print_only=="sync":
        print("\033[1m MARCONI SUPERCOMPUTER \033[0m")
        print("\033[1m --------------------- \033[0m", "\n")
        print("    ", " >> check_runs")
        print("    ", " >> remove_stellaFiles")
        print("    ", " >> remove_stellapyFiles")
        print("    ", " >> load_modulesToCompileStella")
        print()
        print("    ", " >> sync_runsFromMarconi      ($MARCONI_RUNS --> $RUNS)")
        print("    ", " >> sync_newRunsToMarconi     ($NEWRUNS --> $MARCONI_RUNS)")
        print("    ", " >> sync_stellaToMarconi      ($STELLA --> $MARCONI_STELLA)")
        print("    ", " >> sync_stellapyToCommon     ($STELLAPY --> $MARCONI_COMMON_STELLAPY)")
        print("    ", " >> sync_stellapyFromCommon   ($MARCONI_COMMON_STELLAPY --> $STELLAPY)")
        print()
    
    # Launching
    if print_only==None or print_only=="launch":
        print("\033[1m LAUNCH SIMULATIONS \033[0m")
        print("\033[1m ------------------ \033[0m", "\n") 
        print("   ", "\033[1m GENERAL \033[0m")
        print("   ", "\033[1m ------- \033[0m") 
        print("      ", " >> launch -n 4 -h 23") 
        print("      ", " >> launch -n 1 -m 10 --debug ") 
        print("      ", " >> restart -t 1000 -h 23")  
        print() 
        print("   ", "\033[1m LINEAR \033[0m")
        print("   ", "\033[1m ------- \033[0m") 
        print("      ", " >> launch_linearKxKyScan -d 0.1 -t 50 --kx 0 --ky scan_box")  
        print("      ", " >> launch_linearVariableScan -d 0.1 -t 50 --kx 0 --ky scan_box --poloidalturns '[1,2,3]'")   
        print("      ", " >> launch_linearGradientScan -d 0.1 -t 50 --kx 0 --ky scan_box --fprims '[0,1,2]' --tiprims '[3]' --teprims '[0]'")  
        print("      ", " >> launch_linearMapOfFprimVsTiprim -d 0.1 -t 50 --kx 0 --ky scan_box --fprims '[0,1,2]' --tiprims '[0,1,2]' --teprim 0 --rows '[0,1,2]'") 
        print("      ", " >> launch_linearResolutionscan -d 0.1 -t 50 --kx 0 --ky '(1,2,1)' --fprims '[3]' --tiprims '[3]' --teprim '[0]' --nmu '[12,20]'")   
        print("      ", " >> launch_linearScansOfImpurities -t 50 --kx 0 --ky scan_box --tzprim 3 --zeffs '[2.0]' --fprims '[1.0]' --charge_and_mass '[(6.0, 12.011)]'")    
        print()   
        print("   ", "\033[1m LINEAR FFS \033[0m")
        print("   ", "\033[1m ---------- \033[0m") 
        print("      ", " >> launch_rhostarScan (TODO)")  
        print() 
        print("   ", "\033[1m NONLINEAR \033[0m")
        print("   ", "\033[1m --------- \033[0m")   
        print("      ", " >> launch_nonlinearVariableScan -d 0.1 -t 1000 --fprim '[1,2,3,4]'")
        print("      ", " >> launch_nonlinearGradientScan -d 0.1 -t 1000 --fprims '[0,1,2]' --tiprims '[3]' --teprims '[0]'")
        print("      ", " >> launch_nonlinearScansOfImpurities -d 0.1 -t 1000 --tzprim 3 --zeffs '[2.0]' --fprims '[1.0]' --charge_and_mass '[(6.0, 12.011)]'")
        print("      ", " >> launch_nonlinearRadialScanOfExperimentalProfiles -t 500 --rhos '[0.5, 0.6, 0.7]'")
        print("      ", " >> launch_nonlinearDoubleTheResolutionOneByOne -t 500 --all")
        print("      ", " >> launch_nonlinearDefaultResolutionScan -t 1000")
        print() 
        print("   ", "\033[1m CONVERGENCE \033[0m")
        print("   ", "\033[1m ----------- \033[0m") 
        print("      ", " >> check_timetraces")  
        print("      ", " >> check_convergenceModes")  
        print("      ", " >> restart_unconvergedModes")    
        print()
     
    print()
    return

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash interface
    bash = Bash(list_commands_extra, __doc__) 
    
    # Data toggles to print specific parts of the commands
    bash.add_toggle('print_only', 'plot', '', '', 'Show the plotting commands.')  
    bash.add_toggle('print_only', 'linear', '', '', 'Show only the linear plotting commands.')  
    bash.add_toggle('print_only', 'nonlinear', '', '', 'Show only the nonlinear plotting commands.')  
    bash.add_toggle('print_only', 'GUI', '', '', 'Show only the GUI commands.')  
    bash.add_toggle('print_only', 'hpc', '', '', 'Show only the supercomputer (or synchronization) commands.')  
    bash.add_toggle('print_only', 'sync', '', '', 'Show only the supercomputer (or synchronization) commands.')  
    bash.add_toggle('print_only', 'memory', '', '', 'Show only the memory commands.')  
    bash.add_toggle('print_only', 'check', '', '', 'Show only the check commands.')   
    bash.add_toggle('print_only', 'write', '', '', 'Show only the data writing commands.')  
    bash.add_toggle('print_only', 'launch', '', '', 'Show only the launching commands.')  
    
    # Get the arguments and launch the script
    args = bash.get_arguments()
    del args['folder']
    list_commands_extra(**args)   
    
