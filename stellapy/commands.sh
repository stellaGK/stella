#!/bin/bash

# List all the commands
alias stellapy="python3 $STELLAPY/utils/commandprompt/list_commands.py"   


#===============================================================================
#                           Graphical User Interface                           #
#===============================================================================
alias stellaplotter="python3 $STELLAPY/GUI/stella_GUI.py"   
alias stellaplotter_linear="python3 $STELLAPY/GUI/stella_GUI_linear.py"   
alias stellaplotter_nonlinear="python3 $STELLAPY/GUI/stella_GUI_nonlinear.py"   


#===============================================================================
#                               Write data files                               #
#===============================================================================

# Process the data and write small files
alias write_dataFiles="python3 $STELLAPY/data/write/write_dataFiles.py"   

# Make an indentical copy of the *.out.nc file but with less time points
alias reduce_sizeNetcdf="python3 $STELLAPY/data/write/reduce_sizeNetcdf.py"
alias replace_netcdfFile="python3 $STELLAPY/data/write/replace_netcdfFile.py"


#===============================================================================
#                           Plot linear simulations                            #
#===============================================================================

# Plot more complex analyses 
alias plot_gamma="python3 $STELLAPY/plot/linear/gamma_overview_spectrum.py"
alias plot_parameter_influence="python3 $STELLAPY/plot/linear/gamma_overview_spectra.py"
alias plot_spectrum="python3 $STELLAPY/plot/linear/gamma_overview_spectrum.py"
alias plot_spectra="python3 $STELLAPY/plot/linear/gamma_overview_spectra.py"

# Plot spectra
alias plot_gamma_vs_kx="python3 $STELLAPY/plot/linear/gamma_vs_ky.py --kx"
alias plot_gamma_vs_ky="python3 $STELLAPY/plot/linear/gamma_vs_ky.py"
alias plot_omega_vs_kx="python3 $STELLAPY/plot/linear/gamma_vs_ky.py --kx --omega"
alias plot_omega_vs_ky="python3 $STELLAPY/plot/linear/gamma_vs_ky.py --omega"

# Plot the influence of a parameter on {gamma, omega, ky}
alias plot_gamma_vs_parameter="python3 $STELLAPY/plot/linear/gamma_vs_parameter.py"

# Plot time evolution 
alias plot_gamma_vs_time="python3 $STELLAPY/plot/linear/gamma_vs_time.py"  
alias plot_omega_vs_time="python3 $STELLAPY/plot/linear/gamma_vs_time.py --omega"  
alias plot_dphiz_vs_time="python3 $STELLAPY/plot/linear/gamma_vs_time.py --omega"


#===============================================================================
#                          Plot nonlinear simulations                          #
#===============================================================================

# Plot time evolution 
alias plot_flux_vs_time="python3 $STELLAPY/plot/nonlinear/flux_vs_time.py"  
alias plot_qflux_vs_time="python3 $STELLAPY/plot/nonlinear/flux_vs_time.py --qflux"  
alias plot_pflux_vs_time="python3 $STELLAPY/plot/nonlinear/flux_vs_time.py --pflux"  
alias plot_vflux_vs_time="python3 $STELLAPY/plot/nonlinear/flux_vs_time.py --vflux"   


 
