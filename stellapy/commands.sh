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
alias write_dataFiles="python3 $STELLAPY/data/write/write_stellapyDataFiles.py"   
alias write_stellapyDataFiles="python3 $STELLAPY/data/write/write_stellapyDataFiles.py"   

# Make an indentical copy of the *.out.nc file but with less time points
alias reduce_sizeNetcdf="python3 $STELLAPY/data/write/reduce_sizeNetcdf.py"
alias replace_netcdfFile="python3 $STELLAPY/data/write/replace_netcdfFile.py"

# Convert old to new stellapy (October 2022)
alias write_listOfMatchingInputFilesForOldStellapy="python3 $STELLAPY/data/input/write_listOfMatchingInputFilesForOldStellapy.py"


#===============================================================================
#                               Check data files                               #
#===============================================================================

# Check the data files
alias check_netcdfVariables="python3 $STELLAPY/data/check/check_netcdfVariables.py"  
alias check_dimensions="python3 $STELLAPY/data/check/check_dimensions.py"  
alias check_netcdf="python3 $STELLAPY/data/check/check_netcdf.py"  
alias check_timetraces="python3 $STELLAPY/data/check/check_timetraces.py"
alias check_inputs="python3 $STELLAPY/data/check/check_inputs.py"
alias check_cputime="python3 $STELLAPY/data/check/check_cputime.py"
alias check_kperp2="python3 $STELLAPY/data/check/check_kperp2.py"

# Check the differences in files
alias check_differencesInStellaInputs="python3 $STELLAPY/data/check/check_differencesInStellaInputs.py" 

# Check the VMEC
alias print_geometry="python3 $STELLAPY/data/check/print_geometry.py" 
alias check_wout="python3 $STELLAPY/data/check/print_geometry.py" 


#===============================================================================
#                                Plot Geometry                                 #
#===============================================================================

# Plot geometry versus z 
alias plot_geometry_vs_z="python3 $STELLAPY/plot/geometry/geometry_vs_z.py"


#===============================================================================
#                           Plot linear simulations                            #
#===============================================================================

# Plot more complex analyses 
alias plot_gamma="python3 $STELLAPY/plot/linear/gamma_overview_spectra.py"
alias plot_spectra="python3 $STELLAPY/plot/linear/gamma_overview_spectra.py"
alias plot_scans="python3 $STELLAPY/plot/linear/gamma_overview_scans.py" 
alias plot_parameter_influence="python3 $STELLAPY/plot/linear/gamma_overview_scans.py" 

# Plot spectra
alias plot_gamma_vs_kx="python3 $STELLAPY/plot/linear/gamma_vs_wavenumber.py --kx"
alias plot_gamma_vs_ky="python3 $STELLAPY/plot/linear/gamma_vs_wavenumber.py --ky"
alias plot_omega_vs_kx="python3 $STELLAPY/plot/linear/gamma_vs_wavenumber.py --kx --omega"
alias plot_omega_vs_ky="python3 $STELLAPY/plot/linear/gamma_vs_wavenumber.py --ky --omega"

# Plot the influence of a parameter on {gamma, omega, ky}
alias plot_gamma_vs_parameter="python3 $STELLAPY/plot/linear/gamma_vs_parameter.py"

# Plot time evolution 
alias plot_gamma_vs_time="python3 $STELLAPY/plot/linear/gamma_vs_time.py"  
alias plot_omega_vs_time="python3 $STELLAPY/plot/linear/gamma_vs_time.py --omega"  
alias plot_dphiz_vs_time="python3 $STELLAPY/plot/linear/dphiz_vs_time.py"
alias plot_potential_vs_time="python3 $STELLAPY/plot/nonlinear/potential_vs_time.py"   

# Plot velocity space
alias plot_distribution_vs_mu="python3 $STELLAPY/plot/linear/distribution_vs_velocity.py"  
alias plot_distribution_vs_vpa="python3 $STELLAPY/plot/linear/distribution_vs_velocity.py"  
alias plot_distribution_vs_velocity="python3 $STELLAPY/plot/linear/distribution_vs_velocity.py"  

# Parallel mode structure
alias plot_potential_vs_z="python3 $STELLAPY/plot/linear/potential_vs_z.py"  
alias plot_distribution_vs_z="python3 $STELLAPY/plot/linear/distribution_vs_z.py"  


#===============================================================================
#                          Plot nonlinear simulations                          #
#===============================================================================

# Plot time evolution 
alias plot_flux_vs_time="python3 $STELLAPY/plot/nonlinear/flux_vs_time.py"  
alias plot_qflux_vs_time="python3 $STELLAPY/plot/nonlinear/flux_vs_time.py --qflux"  
alias plot_pflux_vs_time="python3 $STELLAPY/plot/nonlinear/flux_vs_time.py --pflux"  
alias plot_vflux_vs_time="python3 $STELLAPY/plot/nonlinear/flux_vs_time.py --vflux"  
alias plot_potential_vs_time="python3 $STELLAPY/plot/nonlinear/potential_vs_time.py"  
alias plot_phi2_vs_time="python3 $STELLAPY/plot/nonlinear/potential_vs_time.py"  

# Plot flux spectra  
alias plot_qflux_vs_wavenumber="python3 $STELLAPY/plot/nonlinear/flux_vs_wavenumber.py --qflux"  
alias plot_pflux_vs_wavenumber="python3 $STELLAPY/plot/nonlinear/flux_vs_wavenumber.py --pflux"  
alias plot_vflux_vs_wavenumber="python3 $STELLAPY/plot/nonlinear/flux_vs_wavenumber.py --vflux"  
alias plot_qflux_spectra="python3 $STELLAPY/plot/nonlinear/flux_vs_wavenumber.py --qflux"  
alias plot_pflux_spectra="python3 $STELLAPY/plot/nonlinear/flux_vs_wavenumber.py --pflux"  
alias plot_vflux_spectra="python3 $STELLAPY/plot/nonlinear/flux_vs_wavenumber.py --vflux"  
alias plot_phi2_spectra="python3 $STELLAPY/plot/nonlinear/potential_vs_wavenumber.py"  
alias plot_potential_spectra="python3 $STELLAPY/plot/nonlinear/potential_vs_wavenumber.py"  
alias plot_qflux_vs_kx="python3 $STELLAPY/plot/nonlinear/flux_vs_wavenumber.py --qflux --kx"  
alias plot_pflux_vs_kx="python3 $STELLAPY/plot/nonlinear/flux_vs_wavenumber.py --pflux --kx"  
alias plot_vflux_vs_kx="python3 $STELLAPY/plot/nonlinear/flux_vs_wavenumber.py --vflux --kx"
alias plot_qflux_vs_ky="python3 $STELLAPY/plot/nonlinear/flux_vs_wavenumber.py --qflux --ky"  
alias plot_pflux_vs_ky="python3 $STELLAPY/plot/nonlinear/flux_vs_wavenumber.py --pflux --ky"  
alias plot_vflux_vs_ky="python3 $STELLAPY/plot/nonlinear/flux_vs_wavenumber.py --vflux --ky" 
alias plot_phi2_vs_kx="python3 $STELLAPY/plot/nonlinear/potential_vs_wavenumber.py --kx"  
alias plot_phi2_vs_ky="python3 $STELLAPY/plot/nonlinear/potential_vs_wavenumber.py --ky" 
alias plot_potential_vs_kx="python3 $STELLAPY/plot/nonlinear/potential_vs_wavenumber.py --kx"  
alias plot_potential_vs_ky="python3 $STELLAPY/plot/nonlinear/potential_vs_wavenumber.py --ky"   

# Parallel mode structure
alias plot_potential_vs_z="python3 $STELLAPY/plot/nonlinear/potential_vs_z.py"   
alias plot_moment_vs_z="python3 $STELLAPY/plot/nonlinear/moment_vs_z.py"   

# Parameter influence
alias plot_qflux_vs_parameter="python3 $STELLAPY/plot/nonlinear/flux_vs_parameter.py --qflux"  
alias plot_pflux_vs_parameter="python3 $STELLAPY/plot/nonlinear/flux_vs_parameter.py --pflux"  
alias plot_vflux_vs_parameter="python3 $STELLAPY/plot/nonlinear/flux_vs_parameter.py --vflux"  
 
