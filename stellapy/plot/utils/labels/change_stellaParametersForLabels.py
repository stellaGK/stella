
import numpy as np  

#===============================================================================
#                 CHANGE THE STELLA (NAME, VALUE) FOR THE PLOTS                #
#===============================================================================  
    
def change_stellaParametersForLabels(simulation, knob, variable, value):
    ''' The labels for the lines/markers have been defined previously based on 
    the stella knobs and keys and can be modified here. ''' 
    
    # Remember the original combination
    variable_temp = variable
    
    # Change the notation of the following stella knob/keys/values
    if variable=="delt":            variable = "$\\Delta t$"
    if variable=="y0":              variable = "$L_y/2\\pi$" 
    if variable=="aky_min":         variable = "$k_y$" 
    if variable=="explicit_option": variable = "Runge kutta order" 
    
    # Show the number of poloidal turns 
    if variable=="nfield_periods":
        if simulation.linear: geometry = simulation.modes[0].geometry
        if simulation.nonlinear: geometry = simulation.geometry
        value = np.round(value/geometry.nfp*abs(geometry.iota),1)
        variable = "Poloidal turns" 

    # Round step sizes to 3 or 4 digits
    if variable=="kx max":  value = str(round(value,1))
    if variable=="ky max":  value = str(round(value,1))
    if variable=="dkx":     value = str(round(value,3))
    if variable=="dkx":     variable = "$\\Delta k_x$" 
    if variable=="dky":     value = str(round(value,3))
    if variable=="dky":     variable = "$\\Delta k_y$" 
    if variable=="dmu":     value = str(round(value,4))
    if variable=="dmu":     variable = "$\\mu_\\textrm{min}$" 
    if variable=="dvpa":    value = str(round(value,4))
    if variable=="dvpa":    variable = "$\\Delta v_{_{\\scriptscriptstyle \\mathbin{\\!/\\mkern-5mu/\\!}}}$" 
    
    # Show rho instead of torflux 
    if variable=="torflux": value = str(round(np.sqrt(float(value)),2))
    if variable=="torflux": variable = "$\\rho$"
    if variable=="rho":     variable = "$\\rho$" 
    
    # Update the species variables
    if knob=='species_parameters_1':
        if variable=="tprim": variable="$a/L_{T_i}$"
        if variable=="fprim": variable="$a/L_{ne}$"
    if knob=='species_parameters_2':
        if variable=="tprim": variable="$a/L_{T_e}$"
        if variable=="fprim": variable="$a/L_{ni}$"
        
    # If we didn't write the variable with a math notation, put it in latex texttt font
    if variable==variable_temp: variable = "\\texttt{"+variable.replace("_", " ")+"}"
    
    # Return the updated variables and values
    return variable, value
