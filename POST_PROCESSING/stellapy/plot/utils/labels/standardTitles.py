

################################################################################
#                        TITLES FOR THE DIFFERENT PLOTS
################################################################################

# Get standard titles for plots, the first key is the y-axis and the second the x-axis 
standardTitles = {
    # FREQUENCY AND GROWTH
    "omega" : {
        "t"   : "Frequency of the modes",\
        "kx"  : "Frequency spectrum",\
        "ky"  : "Frequency spectrum",\
        "para": "Frequency versus parameter",\
        },\
    "gamma" : {
        "t"   : "Growth rate of the modes",\
        "kx"  : "Growth rate spectrum",\
        "ky"  : "Growth rate spectrum",\
        "para": "Growth rate versus parameter",\
        },\
    "gamma/ky**2" : { 
        "kx"  : "Growth rate on $k_y^2$ spectrum",\
        "ky"  : "Growth rate on $k_y^2$ spectrum",\
        "para": "Growth rate on $k_y^2$ versus parameter",\
        },\
    "gamma_integrated" : {
        "t"   : "Integrated growth rate of the modes",\
        "kx"  : "Integrated growth rate spectrum",\
        "ky"  : "Integrated growth rate spectrum",\
        "para": "Integrated growth rate versus parameter",\
        },\
    # POTENTIAL
    "phi" : {
        "t"   : "Time evolution of the electrostatic fluctuations $\\langle \\phi \\rangle$",\
        "z"   : "Parallel mode structure of the electrostatic fluctuations $\\langle \\phi \\rangle$",\
        },\
    "phiReal" : {
        "t"   : "Time evolution of the real part of the electrostatic fluctuations Re$(\\langle \\phi \\rangle)$",\
        "z"   : "Parallel mode structure of the real part of the electrostatic fluctuations Re$(\\langle \\phi \\rangle)$",\
        },\
    "phiImag" : {
        "t"   : "Time evolution of the imaginary part of the electrostatic fluctuations Im$(\\langle \\phi \\rangle)$",\
        "z"   : "Parallel mode structure of the imaginary part of the electrostatic fluctuations Im$(\\langle \\phi \\rangle)$",\
        },\
    "phiRealImag" : {
        "t"   : "Time evolution of the real and imaginary part of the electrostatic fluctuations",\
        "z"   : "Parallel mode structure of the real and imaginary part of the electrostatic fluctuations",\
        },\
    "phi2" : {
        "t"   : "Time evolution of the electrostatic fluctuations squared $\\langle \\phi^2 \\rangle$",\
        "kx"  : "Spectrum along $k_x$ of the\n potential squared $\\sum_{k_y} |\\langle\\hat\\phi\\rangle|^2$",\
        "ky"  : "Spectrum along $k_y$ of the\n potential squared $\\sum_{k_x} |\\langle\\hat\\phi\\rangle|^2$",\
        "z"   : "Parallel mode structure of the electrostatic fluctuations squared",\
        },\
    "phi2_nozonal" : {
        "t"   : "",\
        "kx"  : "Spectrum along $k_x$ of the\n potential squared $\\sum_{k_y\\neq 0} |\\langle\\hat\\phi\\rangle|^2$",\
        "ky"  : "Spectrum along $k_y$ of the\n potential squared $\\sum_{k_x} |\\langle\\hat\\phi\\rangle|^2$",\
        "z"   : "Parallel mode structure of the electrostatic fluctuations squared without zonal modes",\
        } ,\
    "phi2_split" : {
        "t"   : "Time evolution of the electrostatic fluctuations squared, split in zonal contributions",\
        "z"   : "Parallel mode structure of the electrostatic fluctuations squared, split in zonal contributions",\
        },\
    "phi2AndQFlux" : {
        "t"   : "Time evolution of the real and imaginary part of the electrostatic fluctuations",\
        "z"   : "Parallel mode structure of the real and imaginary part of the electrostatic fluctuations",\
        },\
    "phi2_range": {
        "kx"  : "Height of the $\\phi^2(k_x)$ spectrum",\
        "ky"  : "Height of the $\\phi^2(k_y)$ spectrum",\
        "kxNZ": "Height of the $\\phi^2(k_x)_{k_y \\neq 0}$ spectrum",\
        "kyNZ": "Height of the $\\phi^2(k_y)_{k_y \\neq 0}$ spectrum",\
        },\
    # FLUXES
    "qflux" : {
        "t"   : "",\
        "kx"  : "Spectrum along $k_x$ of the\n Heat flux $\\sum_{k_y} Q_{k_x,k_y}$",\
        "ky"  : "Spectrum along $k_y$ of the\n Heat flux $\\sum_{k_x} Q_{k_x,k_y}$",\
        "para": "Saturated heat flux",\
        },\
    "pflux" : {
        "t"   : "",\
        "kx"  : "Spectrum along $k_x$ of the\n Particle flux $\\sum_{k_y} \\Gamma_{k_x,k_y}$",\
        "ky"  : "Spectrum along $k_y$ of the\n Particle flux $\\sum_{k_x} \\Gamma_{k_x,k_y}$",\
        "para": "Saturated particle flux",\
        },\
    "vflux" : {
        "t"   : "",\
        "kx"  : "Spectrum along $k_x$ of the\n Momentum flux $\\sum_{k_y} \\Pi_{k_x,k_y}$",\
        "ky"  : "Spectrum along $k_y$ of the\n Momentum flux $\\sum_{k_x} \\Pi_{k_x,k_y}$",\
        "para": "Saturated momentum flux",\
        },\
    "qflux_range": {
        "kx"  : "Height of the $Q(k_x)$ spectrum",\
        "ky"  : "Height of the $Q(k_y)$ spectrum",\
        "kxNZ": "Height of the $Q(k_x)$ spectrum",\
        "kyNZ": "Height of the $Q(k_y)$ spectrum",\
        },\
    # DISTRIBUTION
    "gzvs" : {
        "-"   : "",\
        "vpa" : "Guiding centre distribution\n along the parallel velocity",\
        "mu"  : "Guiding centre distribution\n along the perpendicular velocity",\
        } ,\
    "gvmus" : {
        "t"   : "",\
        "vpa" : "Guiding centre distribution\n along the parallel velocity",\
        "mu"  : "Guiding centre distribution\n along the perpendicular velocity",\
        },\
    # MOMENTS 
    "n" : {
        "t"   : "Time evolution of the density $n$",\
        } ,\
    "T" : {
        "t"   : "Time evolution of the temperature $T$",\
        } ,\
    "v" : {
        "t"   : "Time evolution of the parallel velocity $v_\\parallel$",\
        } ,\
    # DISTRIBUTION
    "g" : {
        "t"   : "Distribution function of the guiding centres",\
        "vpa" : "Distribution function of the guiding centres versus $v_\\parallel$",\
        "mu"  : "Distribution function of the guiding centres versus $\\mu$",\
        },\
    }

# Load the other dictionaries
from stellapy.plot.utils.labels.standardNames import standardNames
from stellapy.plot.utils.labels.standardParameters import standardParameters

# For each title[quantity][para] replace [para] with all the options <standardParameters>
for quantity in ["omega", "gamma", "gamma/ky**2", "gamma_integrated"]: 
    if "para" in standardTitles[quantity].keys():
        for para_quantity in standardParameters.keys():
            if standardNames[para_quantity]!=None:
                replacement = standardNames[para_quantity][0].lower() + standardNames[para_quantity][1:]
                standardTitles[quantity][para_quantity] = standardTitles[quantity]["para"].replace("parameter", replacement)
            if standardNames[para_quantity]==None:
                standardTitles[quantity][para_quantity] = None 