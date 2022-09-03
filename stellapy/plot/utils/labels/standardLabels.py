''' LABELS DICTIONARY

For each quantity (qflux, gamma, phi, time, kx, ky, dt, ...) that we are interested 
in, we define labels for the (x,y) axis in <standardLabels> in normalized units
and in SI units, a name for the (x,y) quantities for the optionswindow in <standardNames>,
a name in latex math notation to complement matplotlib titles in <standardNamesScientific>
a title for the matplotlib plot in <standardTitles>, a header for the tkinter windows
in <standardHeaders> and interesting stella parameters with their corresponding stella
 <knob> and <key> in <standardParameters>.

Note that we do not represent the normalized labels of phi correctly, 
technically they are $\\phi (\\rho_r/a) (T_e/e)$ but this is too much.
'''

# Get the standard labels for the (x,y) axis
standardLabels = {

########################################
# STANDARD LABELS IN NORMALIZED UNITS
########################################

    "normalized" : {
        # TIME AXIS:
        "t"             : "$t\, v_{\mathrm{th},r}/a$",\
        # Z AXIS:
        "x"             : "$x/\\rho_i$",\
        "y"             : "$y/\\rho_i$",\
        "z"             : "$z/a$",\
        "zeta"          : "$\\zeta$",\
        "pol"           : "Poloidal turns",\
        "tor"           : "Toroidal turns",\
        # KX OR KY AXIS:
        "kx"            : "$k_{x}\\rho_i$",\
        "ky"            : "$k_{y}\\rho_i$",\
        "ky**2"         : "$(k_{y}\\rho_i)^{2}$",\
        "1/ky**2"       : "$(k_{y}\\rho_i)^{-2}$",\
        # VPA OR MU AXIS:
        "vpa"           : "$v_\\parallel/v_{\mathrm{th},{s}}$",\
        "mu"            : "$\\mu B_r/v^2_{\mathrm{th},{s}}$",\
        # GROWTHRATE AND FREQUENCY
        "omega"         : "$\\omega a/v_{\\mathrm{th},i}$",\
        "gamma"         : "$\\gamma a/v_{\\mathrm{th},i}$",\
        "gamma/ky"      : "$\\gamma a/(v_{\\mathrm{th},i} * k_{y}\\rho_i)$",\
        "gamma/ky^2"   : "$\\gamma a/(v_{\\mathrm{th},i} * (k_{y}\\rho_i)^2)$",\
        "gamma/ky**2"   : "$\\gamma a/(v_{\\mathrm{th},i} * (k_{y}\\rho_i)^2)$",\
        "gamma/ky**3"   : "$\\gamma a/(v_{\\mathrm{th},i} * (k_{y}\\rho_i)^3)$",\
        "gamma_integrated" : "$\\sum_\\gamma (\\gamma a/v_{\\mathrm{th},i}) * (\\Delta k_{y}\\rho_i)$",\
        # FLUXES
        "qflux"         : "$Q_{s}/Q_{\\text{gB},i}$",\
        "pflux"         : "$\\Gamma_{s}/\\Gamma_{\\text{gB},i}$",\
        "vflux"         : "$\\Pi_{s}/\\Pi_{\\text{gB},i}$",\
        # POTENTIAL
        "phi"           : "$\\langle\\hat\\phi\\rangle$",\
        "phi_real"      : "Re$(\\langle\\hat\\phi\\rangle)$",\
        "phi_imag"      : "Im$(\\langle\\hat\\phi\\rangle)$",\
        "phiRealImag"   : "Re$(\\langle\\hat\\phi\\rangle)$ and Im$(\\langle\\hat\\phi\\rangle)$",\
        # POTENTIAL SQUARED
        "phi2"          : "$|\\langle\\hat\\phi\\rangle^2|$",\
        "phi2_allModes" : "$\\sum_{k_x,k_y} |\\langle\\hat\\phi\\rangle|^2$",\
        "phi2_zonal"    : "$\\sum_{k_x,k_y=0} |\\langle\\hat\\phi\\rangle|^2$",\
        "phi2_nozonal"  : "$\\sum_{k_x,k_y \\neq 0} |\\langle\\hat\\phi\\rangle|^2$",\
        # DISTRIBUTION FUNCTION
        "g"             : "$|g_{s}|^2$",\
        "ge"            : "$|g_{e}|^2$",\
        "gi"            : "$|g_{i}|^2$",\
        "distribution"  : "$|g_{s}|^2$",\
        "gvmus"         : "$|g_{s}|^2$",\
        "gzvs"          : "$|g_{s}|^2$",\
        # MOMENTS 
        "n"             : "$n$/normalization",\
        "T"             : "$T$/normalization",\
        "v"             : "$v_\\parallel$/normalization",\
        # PARAMETERS:
        "parameter"     : "Stella parameter",\
        "rho"           : "$\\rho$",\
        "tprim"         : "$a/L_{T_i}$",\
        "tiprim"        : "$a/L_{T_i}$",\
        "teprim"        : "$a/L_{T_e}$",\
        "fprim"         : "$a/L_{n}$",\
        "tite"          : "$T_i/T_e$",\
        "teti"          : "$T_e/T_i$",\
        "delt"          : "$\\Delta t$",\
        "delta t"       : "$\\Delta t$",\
        "rk"            : "Order of the Runge Kutta scheme",\
        "nz"            : "$N_z$",\
        "nzed"          : "$N_z$",\
        "nzgrid"        : "$N_z$",\
        "nfield"        : "Poloidal turns",\
        "nperiod"       : "nperiod",\
        "nmu"           : "$N_{\\mu}$",\
        "nvgrid"        : "$N_{v_\\parallel}$",\
        "dmu"           : "$\\mu_\\textrm{min}$",\
        "dvpa"          : "$\\Delta v_{_{\\scriptscriptstyle \\mathbin{\\!/\\mkern-5mu/\\!}}}$",\
        "y0"            : "$L_y/2\\pi$",\
        "nx"            : "$N_x$",\
        "ny"            : "$N_y$",\
        "kx max"        : r"$k_{x, {\textrm{\large{max}}}}$",\
        "ky max"        : r"$k_{y, {\textrm{\large{max}}}}$",\
        "dkx"           : "$\\Delta k_x$",\
        "dky"           : "$\\Delta k_y$",\
        "Lx"            : "$L_x/\\rho_i$",\
        "Ly"            : "$L_y/\\rho_i$",\
        "explicit_option" : "Runge kutta scheme",\
        "nfield_periods" : "Number of field periods",\
        "D_hyper"       : "Hyper-dissipation",\
        "boundary_option" : "Boundary option",\
        # EMPTY
        None            : "",\
        "-----"         : "",\
        "-"             : "",\
        },\
    
    
################################
# STANDARD LABELS IN SI UNITS
################################

    "SI" : { 
        # TIME AXIS
        "t"             : "$t$ [ms]",\
        "time"          : "$t$ [ms]",\
        # Z AXIS:
        "z"             : "$z$ [m]",\
        "zeta"          : "$\\zeta$",\
        "pol"           : "Poloidal turns",\
        "tor"           : "Toroidal turns",\
        # KX OR KY AXIS:
        "kx"            : "$k_{x} [cm$^{-1}$]$",\
        "ky"            : "$k_{y} [cm$^{-1}$]$",\
        # VPA OR MU AXIS:
        "vpa"           : "$v_\\parallel [m/s]$",\
        "mu"            : "$\\mu [$J/($kg$*T)$]$",\
        # GROWTHRATE AND FREQUENCY
        "gamma"         : "$\\gamma$  [s$^{-1}$]",\
        "omega"         : "$\\omega$ [s$^{-1}$]",\
        "gamma/ky**2"   : "$\\gamma/k_y^2$ [s$^{-1}$cm$^{-2}$]",\
        "gamma_integrated" : "$\\sum_\\gamma (\\gamma a/v_{\\mathrm{th},i}) * (\\Delta k_{y}\\rho_i)$",\
        # FLUXES
        "qflux"         : "$Q_{s}$ [W/m$^2$]",\
        "pflux"         : "$\\Gamma$ [$1/($m$^2s)$]",\
        "vflux"         : "$\\Pi_{s}$ [N/m]",\
        # POTENTIAL
        "phi"           : "$\\langle\\hat\\phi\\rangle$ [J]",\
        "phi_real"      : "Re$(\\langle\\hat\\phi\\rangle)$ [J]",\
        "phi_imag"      : "Im$(\\langle\\hat\\phi\\rangle)$ [J]",\
        "phiRealImag"   : "Re$(\\langle\\hat\\phi\\rangle)$ [J] and Im$(\\langle\\hat\\phi\\rangle)$ [J]",\
        # POTENTIAL
        "phi2"          : "$|\\langle\\hat\\phi\\rangle^2|$  [J]",\
        "phi2_allModes" : "$\\sum_{k_x,k_y} |\\langle\\hat\\phi\\rangle|^2$  [J]",\
        "phi2_zonal"    : "$\\sum_{k_x,k_y=0} |\\langle\\hat\\phi\\rangle|^2$  [J]",\
        "phi2_nozonal"  : "$\\sum_{k_x,k_y \\neq 0} |\\langle\\hat\\phi\\rangle|^2$  [J]",\
        # DISTRIBUTION FUNCTION
        "g"             : "$g$ [s$^3$/m$^6$]",\
        "distribution"  : "$g$ [s$^3$/m$^6$]",\
        "gvmus"         : "$g_{s}$ [s$^3$/m$^6$]",\
        "gzvs"          : "$g_{s}$ [s$^3$/m$^6$]",\
        # MOMENTS 
        "n"             : "$n$ units",\
        "T"             : "$T$ units",\
        "v"             : "$v_\\parallel$ units",\
        # PARAMETERS:
        "parameter"     : "Stella parameter",\
        "rho"           : "$\\rho$",\
        # EMPTY
        None            : "",\
        },\


#####################################################
# STANDARD LABELS IN NORMALIZED UNITS AND RESCALED
#####################################################

    "rescaled" : {
        # TIME AXIS:
        "t"         : "$t/t_{peak}$",\
        "time"      : "$t/t_{peak}$",\
        # FLUXES
        "qflux"     : "$Q_{s}/Q_{sat}$",\
        "pflux"     : "$\\Gamma_{s}/\\Gamma_{sat}$",\
        "vflux"     : "$\\Pi_{s}/\\Pi_{sat}$",\
        # POTENTIAL
        "phi2"      : "$|\\phi^2|$",\
        # DISTRIBUTION FUNCTION
        "gvmus"     : "$\\sum_{v_\\parallel, v_\\perp} g$",\
        "gzvs"      : "$\\sum_{v_\\parallel, z} g$",\
        # EMPTY
        None        : "",\
        },\
    }



###########################################
# RESCALED LABEL FOR EACH NORMALIZED LABEL
###########################################
for quantity in list(standardLabels["normalized"].keys()):  
    if quantity not in list(standardLabels["rescaled"].keys()): 
        standardLabels["rescaled"][quantity] = standardLabels["normalized"][quantity] 
        
        

#####################################################
# CREATE LABELS WHEN WE HAVE SUMMED AWAY A DIMENSION
#####################################################

for units in ["normalized", "SI"]:
    for quant in ['qflux', 'pflux', 'vflux', 'phi', 'phi_real', 'phi_imag', 'phiRealImag']:
        standardLabels[units][quant+" vs kx"] = "$\\sum_{k_y}\,$"+standardLabels[units][quant]
        standardLabels[units][quant+" vs ky"] = "$\\sum_{k_x}\,$"+standardLabels[units][quant] 
    for quant in ['phi2' ]:
        standardLabels[units][quant+" vs kx"] = "$\\sum_{k_y}\,$"+standardLabels[units][quant]
        standardLabels[units][quant+" vs ky"] = "$\\sum_{k_x}\,$"+standardLabels[units][quant] 
        standardLabels[units][quant+" vs z"]  = "$|$"+standardLabels[units][quant]+"/"+standardLabels[units][quant]+"$_{max}|$"
    for quant in ['phi2_zonal']:
        standardLabels[units][quant+" vs kx"] = "$\\sum_{k_y=0}\,$"+standardLabels[units]["phi2"]
        standardLabels[units][quant+" vs ky"] = standardLabels[units]["phi2"] 
        standardLabels[units][quant+" vs z"]  = standardLabels[units][quant].replace(standardLabels[units]["phi2"], standardLabels[units]["phi2 vs z"])
    for quant in ['phi2_nozonal']:
        standardLabels[units][quant+" vs kx"] = "$\\sum_{k_y\\neq 0}\,$"+standardLabels[units]["phi2"]
        standardLabels[units][quant+" vs ky"] = "$\\sum_{k_x}\,$"+standardLabels[units]["phi2"]
        standardLabels[units][quant+" vs z"]  = standardLabels[units][quant].replace(standardLabels[units]["phi2"], standardLabels[units]["phi2 vs z"])
    for quant in ['gvmus']:
        standardLabels[units][quant+" vs vpa"] = "$\\sum_{\\mu}\,$"+standardLabels[units][quant].replace("_{s}", "_{s}(v_\\parallel, \\mu)")
        standardLabels[units][quant+" vs mu"]  = "$\\sum_{v_\\parallel}\,$"+standardLabels[units][quant].replace("_{s}", "_{s}(v_\\parallel, \\mu)")  
        standardLabels[units][quant+" vs z"]   = "$\\sum_{v_\\parallel}\,$"+standardLabels[units][quant].replace("_{s}", "_{s}(v_\\parallel, \\mu)")
    for quant in ['gzvs']:
        standardLabels[units][quant+" vs vpa"] = "$\\sum_{z}$"+standardLabels[units][quant].replace("_{s}", "_{s}(v_\\parallel, z)")
        standardLabels[units][quant+" vs mu"]  = "$\\sum_{v_\\parallel}\,$"+standardLabels[units][quant].replace("_{s}", "_{s}(v_\\parallel, \\mu)")  
        standardLabels[units][quant+" vs z"]   = "$\\sum_{v_\\parallel}\,$"+standardLabels[units][quant].replace("_{s}", "_{s}(v_\\parallel, \\mu)")
    for quant in ['phi', 'phi_real', 'phi_imag', 'phiRealImag']:
        standardLabels[units][quant+" vs z"] = "$|$"+standardLabels[units][quant]+"/"+standardLabels[units][quant]+"$_{max}|$"
    
    
    
###############################
# SYNONYMS FOR THE LABEL KEYS
################################
for units in ["normalized", "SI"]:
    standardLabels[units]["phi2zonal"] = standardLabels[units]["phi2_zonal"]
    standardLabels[units]["phi2nozonal"] = standardLabels[units]["phi2_nozonal"]
    standardLabels[units]["vec_kx"] = standardLabels[units]["kx"]
    standardLabels[units]["vec_ky"] = standardLabels[units]["ky"]  
    standardLabels[units]["phi2_split"] = standardLabels[units]["phi2"]
    standardLabels[units]["phiRealImag"] = standardLabels[units]["phi"]
    standardLabels[units]["phi2AndQFlux"] = standardLabels[units]["phi2"]+" and "+standardLabels[units]["qflux"]
    for i in ["gamma", "omega"]:
        for j in ["_avg", "_last"]:
            standardLabels[units][i+j] = standardLabels[units][i]
            
