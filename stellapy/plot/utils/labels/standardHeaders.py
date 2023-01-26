

################################################################################
#                        GUI HEADERS FOR THE DIFFERENT PLOTS
################################################################################

# Get standard headers for windows, the first key is the y-axis and the second the x-axis 
standardHeaders = {
    # FREQUENCY AND GROWTH
    "omega" : {
        "t"   : "Time evolution of the frequency of the modes",\
        "kx"  : "Frequency spectrum",\
        "ky"  : "Frequency spectrum",\
        "para": "Frequency versus parameter",\
        },\
    "gamma" : {
        "t"   : "Time evolution of the growth rate of the modes",\
        "kx"  : "Growth rate spectrum",\
        "ky"  : "Growth rate spectrum",\
        "para": "Growth rate versus parameter",\
        },\
    "gamma/ky**2" : { 
        "kx"  : "Growth rate on ky**2 spectrum",\
        "ky"  : "Growth rate on ky**2 spectrum",\
        "para": "Growth rate on ky**2 versus parameter",\
        },\
    "gamma_integrated" : {
        "t"   : "Integrated growth rate of the modes",\
        "kx"  : "Integrated growth rate spectrum",\
        "ky"  : "Integrated growth rate spectrum",\
        "para": "Integrated growth rate versus parameter",\
        },\
    # POTENTIAL
    "phi" : {
        "t"   : "Time evolution of the electrostatic fluctuations",\
        "z"   : "Parallel mode structure of the electrostatic fluctuations",\
        "kx"  : "Contributions along kx of the electrostatic fluctuations",\
        "ky"  : "Contributions along ky of the electrostatic fluctuations",\
        },\
    "phiReal" : {
        "t"   : "Time evolution of the real part of the electrostatic fluctuations",\
        "z"   : "Parallel mode structure of the real part of the electrostatic fluctuations",\
        "kx"  : "Contributions along kx of the real part of the electrostatic fluctuations",\
        "ky"  : "Contributions along ky of the real part of the electrostatic fluctuations",\
        },\
    "phiImag" : {
        "t"   : "Time evolution of the imaginary part of the electrostatic fluctuations",\
        "z"   : "Parallel mode structure of the imaginary part of the electrostatic fluctuations",\
        "kx"  : "Contributions along kx of the imaginary part of the electrostatic fluctuations",\
        "ky"  : "Contributions along ky of the imaginary part of the electrostatic fluctuations",\
        },\
    "phiRealImag" : {
        "t"   : "Time evolution of the real and imaginary part of the electrostatic fluctuations",\
        "z"   : "Parallel mode structure of the real and imaginary part of the electrostatic fluctuations",\
        "kx"  : "Contributions along kx of the real and imaginary part of the electrostatic fluctuationq",\
        "ky"  : "Contributions along ky of the real and imaginary part of the electrostatic fluctuations",\
        },\
    "phi2" : {
        "t"   : "Time evolution of the electrostatic fluctuations squared",\
        "z"   : "Parallel mode structure of the electrostatic fluctuations squared",\
        "kx"  : "Contributions along kx of the electrostatic fluctuations squared",\
        "ky"  : "Contributions along ky of the electrostatic fluctuations squared",\
        },\
    "phi2_nozonal" : {
        "t"   : "Time evolution of the electrostatic fluctuations squared without zonal modes",\
        "z"   : "Parallel mode structure of the electrostatic fluctuations squared without zonal modes",\
        "kx"  : "Contributions along kx of the electrostatic fluctuations squared without zonal modes",\
        "ky"  : "Contributions along ky of the electrostatic fluctuations squared without zonal modes",\
        } ,\
    "phi2_split" : {
        "t"   : "Time evolution of the electrostatic fluctuations squared, split in zonal contributions",\
        "z"   : "Parallel mode structure of the electrostatic fluctuations squared, split in zonal contributions",\
        "kx"  : "Contributions along kx of the electrostatic fluctuations squared split in zonal contributions",\
        "ky"  : "Contributions along ky of the electrostatic fluctuations squared split in zonal contributions",\
        },\
    "phi2AndQFlux" : {
        "t"   : "Time evolution of the electrostatic fluctuations squared and the heat flux",\
        "z"   : "Parallel mode structure of the the electrostatic fluctuations squared and heat flux",\
        "kx"  : "Contributions along kx of the electrostatic fluctuations squared and heat flux",\
        "ky"  : "Contributions along ky of the electrostatic fluctuations squared and heat flux",\
        },\
    "phi2_range": {
        "kx"  : "Height of the Phi²(kx) spectrum",\
        "ky"  : "Height of the Phi²(ky) spectrum",\
        "kxNZ": "Height of the Phi²(kx) spectrum without zonal modes",\
        "kyNZ": "Height of the Phi²(ky) spectrum without zonal modes",\
        "para": "Height of the Phi²(kx) or Phi²(ky) spectrum",\
        },\
    # FLUXES
    "qflux" : {
        "t"   : "Time evolution of the heat flux",\
        "z"   : "Parallel mode structure of the heat flux",\
        "para": "Saturated heat flux",\
        "kx"  : "Heat flux",\
        "ky"  : "Heat flux",\
        },\
    "pflux" : {
        "t"   : "Time evolution of the particle flux",\
        "z"   : "Parallel mode structure of the particle flux",\
        "para": "Saturated particle flux",\
        },\
    "vflux" : {
        "t"   : "Time evolution of the momentum flux",\
        "z"   : "Parallel mode structure of the momentum flux",\
        "para": "Saturated momentum flux",\
        },\
    "qflux_range": {
        "kx"  : "Height of the Q(kx) spectrum",\
        "ky"  : "Height of the Q(ky) spectrum",\
        "kxNZ": "Height of the Q(kx) spectrum",\
        "kyNZ": "Height of the Q(ky) spectrum",\
        "para": "Height of the Q(ky) spectrum",\
        },\
    # MOMENTS 
    "n" : {
        "t"   : "Time evolution of the density",\
        } ,\
    "T" : {
        "t"   : "Time evolution of the temperature",\
        } ,\
    "v" : {
        "t"   : "Time evolution of the parallel velocity",\
        } ,\
    # DISTRIBUTION
    "g" : {
        "t"   : "Distribution function of the guiding centres",\
        "z"   : "Distribution function of the guiding centres",\
        "kx"  : "Distribution function of the guiding centres",\
        "ky"  : "Distribution function of the guiding centres",\
        "vpa" : "Distribution function of the guiding centres along the parallel velocity",\
        "mu"  : "Distribution function of the guiding centres along the perpendicular velocity",\
        },\
    }


# Load the other dictionaries
from stellapy.plot.utils.labels.standardNames import standardNames
from stellapy.plot.utils.labels.standardParameters import standardParameters

# Create some synonyms (t, time); (z, zeta, pol, tor); (g, distribution, gzvs, gvmus 
for y_quantity in standardHeaders.keys(): 
    if "t" in standardHeaders[y_quantity].keys(): 
        standardHeaders[y_quantity]["time"] = standardHeaders[y_quantity]["t"]
for y_quantity in standardHeaders.keys(): 
    if "z" in standardHeaders[y_quantity].keys():
        for z_quantity in ["zeta", "pol", "tor"]:
            standardHeaders[y_quantity][z_quantity] = standardHeaders[y_quantity]["z"]
if "g" in list(standardHeaders.keys()):  
    standardHeaders["distribution"] = {}
    standardHeaders["gzvs"] = {}
    standardHeaders["gvmus"] = {}
    for x_quantity in standardHeaders["g"].keys(): 
        standardHeaders["distribution"][x_quantity] = standardHeaders["g"][x_quantity]
        standardHeaders["gzvs"][x_quantity] = standardHeaders["g"][x_quantity]
        standardHeaders["gvmus"][x_quantity] = standardHeaders["g"][x_quantity]
            
# For each header[quantity][para] replace [para] with all the options <standardParameters>
for y_quantity in standardHeaders.keys(): 
    if "para" in standardHeaders[y_quantity].keys(): 
        standardHeaders[y_quantity]["parameter"] = standardHeaders[y_quantity]["para"]
for y_quantity in standardHeaders.keys(): 
    if "para" in standardHeaders[y_quantity].keys():
        for para_quantity in standardParameters.keys():
            if standardNames[para_quantity]!=None:
                replacement = standardNames[para_quantity][0].lower() + standardNames[para_quantity][1:]
                standardHeaders[y_quantity][para_quantity] = standardHeaders[y_quantity]["para"].replace("parameter", replacement)
            if standardNames[para_quantity]==None:
                standardHeaders[y_quantity][para_quantity] = None









