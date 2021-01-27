#!/bin/bash

#########################################################
# Location of the bashscripts
#########################################################

BASHSCRIPTS="$STELLAPY/commands"

#########################################################
# Sync runs and code between local pc and super pc
#########################################################

alias sync_newRunsToMarconi="$BASHSCRIPTS/runSimulations/sync_newRunsToMarconi.sh"
alias sync_runsFromMarconi="$BASHSCRIPTS/supercomputer_marconi/sync_runsFromMarconi.sh"
alias sync_runsToPermanentMemory="$BASHSCRIPTS/supercomputer_marconi/syns_runsToPermanentMemory.sh"
alias sync_stellaToMarconi="$BASHSCRIPTS/supercomputer_marconi/sync_stellaToMarconi.sh"
alias sync_stellaToXula="$BASHSCRIPTS/supercomputer_xula/sync_stellaToXula.sh"

#########################################################
# Commands for simulations
#########################################################

SIMULATIONS="$BASHSCRIPTS/runSimulations"
alias launch="python3 $SIMULATIONS/launch.py"
alias restart="python3 $SIMULATIONS/restart.py"
alias kyscan="python3 $SIMULATIONS/kyscan.py"
alias dtscan="python3 $SIMULATIONS/dtscan.py"
alias cpuscan="python3 $SIMULATIONS/cpuscan.py"
alias resolutionscan="python3 $SIMULATIONS/resolutionScan.py"
alias resolutionscannonlinear="python3 $SIMULATIONS/resolutionScanNonlinear.py"
alias create_newSimulationFolder="python3 $SIMULATIONS/create_newSimulationFolder.py"

#########################################################
# Commands for python: info
#########################################################

alias help_stellapy="python3 $BASHSCRIPTS/help_stellapy.py"

#########################################################
# Commands for python: STELLAPY
#########################################################

# Script directories
REDUCE="$BASHSCRIPTS/reduceSizeFiles"
CHECK="$BASHSCRIPTS/checkStellaOutputs"
PLOT="$BASHSCRIPTS/plots"

# REDUCE SIZE FILES
alias reduce_sizeFiles="$REDUCE/reduce_sizeFiles.sh"
alias reduce_sizeNetcdf="python3 $REDUCE/reduce_sizeNetcdf.py"
alias reduce_sizeWout="python3 $REDUCE/reduce_sizeWout.py"

# CHECK STELLA OUTPUT FILES
alias check_netcdf="python3 $CHECK/check_netcdf.py"
alias check_inputs="python3 $CHECK/check_inputs.py"
alias check_resolution="python3 $CHECK/check_resolution.py"
alias check_wout="python3 $CHECK/check_wout.py" 
alias check_vmecgeo="python3 $CHECK/check_vmecgeo.py"
alias check_referenceUnits="python3 $CHECK/check_referenceUnits.py"
alias check_simulationTime="python3 $CHECK/check_simulationTime.py"

# FREQUENCY AND GROWTHRATE 

    # For linear modes
alias plot_frequencyVsTime="python3 -W ignore $PLOT/plot_frequencyVsTime.py"
alias plot_growthrateVsTime="python3 -W ignore $PLOT/plot_growthrateVsTime.py"
alias plot_frequencyVsModes="python3 $PLOT/plot_frequencyVsModes.py"
alias plot_growthrateVsModes="python3 $PLOT/plot_growthrateVsModes.py"
alias plot_frequencyVsModes_growthrateVsModes="python3 $PLOT/plot_frequencyVsModes_growthrateVsModes.py"

    # Linear comparison of multiple radial position
alias plot_frequencyVsRho="python3 $PLOT/plot_frequencyVsRho.py"
alias plot_growthrateVsRho="python3 $PLOT/plot_growthrateVsRho.py"

# POTENTIAL
alias plot_potentialVsZ="python3 $POTENTIAL/plot_potentialVsZ.py"
 






