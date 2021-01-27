
import textwrap

print(textwrap.dedent('''

#########################################################
# Sync runs and code between local pc and super pc
#########################################################

>>>  sync_newRunsToMarconi
>>>  sync_runsFromMarconi
>>>  sync_stellaToMarconi
>>>  sync_stellaToXula

#########################################################
# Commands for simulations
#########################################################

>>>  create_newSimulationFolder
>>>  reduce_sizeNetcdf
>>>  reduce_sizeWout
>>>  reduce_sizeFiles             (applies the above 2 commands to all folders in $RUNS)
>>>  move_to_subfolder

#########################################################
# Commands for python: STELLAPY
#########################################################

# STELLA FILES
>>>  print_netcdf
>>>  check_inputs
>>>  read_netcdf
>>>  read_resolution

# CONVERGENCE
>>>  plot_heatfluxVsT_compareSaturatedHeatflux
>>>  plot_heatfluxVsT_compareGrowthRate
>>>  plot_frequencyVsModes

# FREQUENCY AND GROWTHRATE 

 # Read and write linear data
>>>  read_lineardata
>>>  write_lineardata
>>>  write_growthrateVsRho

 # For linear modes
>>>  plot_frequencyVsTime
>>>  plot_growthrateVsTime
>>>  plot_frequencyVsModes --SI --noSI --stable
>>>  plot_growthrateVsModes --SI --noSI --stable
>>>  plot_frequencyVsModes_growthrateVsModes

 # Linear comparison of multiple radial position
>>>  plot_frequencyVsRho
>>>  plot_growthrateVsRho

 # Linear instability decomposition
>>>  plot_frequencyVsModes_instabilityDecomposition
>>>  plot_growthrateVsModes_instabilityDecompositio
>>>  plot_frequencyVsModes_growthratesVsModes_instabilityDecomposition

 # Nonlinear modes
>>>  plot_frequencyVsRho_nonlinear
>>>  plot_growthrateVsRho_nonlinear

# POTENTIAL
>>>  plot_PotentialVsZ
>>>  plot_potentialVsTime
>>>  plot_potentialVsModesVsZ
>>>  plot_potentialVsKxVsKy

 # Movies
>>>  animate_potentialVsKxVsKy

 # Investigate the confinement of the potential along the flux tube
>>>  write_lengthPotential
>>>  plot_lengthPotentialVsModes

# FLUXES
>>>  plot_fluxesVsTime
>>>  plot_heatfluxVsTime_growthrateVsRho
>>>  plot_heatfluxVsTime_linearAndLog
>>>  plot_saturatedfluxVsRho_heatfluxVsTime
>>>  plot_saturatedfluxVsRho_growthrateVsRho

# PROFILES
>>>  write_profiles
>>>  add_krangeToProfile
>>>  plot_profilesVsRho_gradientsVsRho_loggradVsRho

# MAGNETIC FIELD
>>>  print_wout
>>>  read_vmecgeo
>>>  read_wout
>>>  plot_parametersVsRho
>>>  plot_comparisonGeometry

'''))




