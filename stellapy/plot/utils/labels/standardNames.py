################################################################################
#                NAMES FOR THE AXIS LABEL IN THE OPTIONS WINDOW
################################################################################

# Get the names for the options window
standardNames = {
    # TIME AXIS:
    "t"             : "Time (t)",\
    "time"          : "Time (t)",\
    # Z AXIS:
    "z"             : "Parallel coordinate (z)",\
    "zeta"          : "Parallel coordinate (zeta)",\
    "pol"           : "Parallel coordinate (poloidal turns)",\
    "tor"           : "Parallel coordinate (toroidal turns)",\
    # KX OR KY AXIS:
    "kx"            : "Wavenumber (kx)",\
    "ky"            : "Wavenumber (ky)",\
    # VPA OR MU AXIS:
    "vpa"           : "Parallel velocity (vpa)",\
    "mu"            : "Magnetic moment (mu)",\
    # GROWTHRATE AND FREQUENCY
    "omega"         : "Frequency (omega)",\
    "gamma"         : "Growth rate (gamma)",\
    "gamma/ky**2"   : "Growth rate (gamma/ky**2)",\
    "gamma_integrated" : "Integrated growth rate (gamma)",\
    # GROWTHRATE AND FREQUENCY
    "omega_avg"     : "Frequency (omega)",\
    "gamma_avg"     : "Growth rate (gamma)",\
    "gamma_avg/ky**2" : "Growth rate (gamma/ky**2)",\
    # FLUXES
    "qflux"         : "Heat flux (Q)",\
    "pflux"         : "Particle flux (Gamma)",\
    "vflux"         : "Momentum flux (Pi)",\
    "qflux_range"   : "Height of the Q(kx) or Q(ky) spectrum",\
    # POTENTIAL
    "phi"           : "Electrostatic potential (phi)",\
    "phi_real"      : "Electrostatic potential (Re(phi))",\
    "phi_imag"      : "Electrostatic potential (Im(phi))",\
    "phiRealImag"   : "Electrostatic potential (Re(phi) and Im(phi))",\
    # POTENTIAL SQUARED
    "phi2"          : "Electrostatic potential squared (phi**2)",\
    "phi2_allModes" : "Electrostatic potential squared (phi**2)",\
    "phi2_zonal"    : "Electrostatic potential squared (phi**2)",\
    "phi2_nozonal"  : "Electrostatic potential squared (phi**2)",\
    "phi2_split"    : "Electrostatic potential squared (phi**2)",\
    "phi2_range"    : "Height of the (phi**2) spectrum",\
    "phi2AndQFlux"  : "Electrostatic potential squared (phi**2) and Heat flux (Q)",\
    # DISTRIBUTION FUNCTION
    "g"             : "Distribution function of the guiding centers (g)",\
    "distribution"  : "Distribution function of the guiding centers (g)",\
    "gvmus"         : "Distribution function of the guiding centers (g)",\
    "gzvs"          : "Distribution function of the guiding centers (g)",\
    # MOMENTS 
    "n"             : "Density (n)",\
    "T"             : "Temperature (T)",\
    "v"             : "Parallel velocity (vpa)",\
    # PARAMETERS:
    "parameter"     : "Stella parameter",\
    "rho"           : "Radial coordinate (rho)",\
    "tprim"         : "Normalized ion temperature gradient length (a/LTi)",\
    "tiprim"        : "Normalized ion temperature gradient length (a/LTi)",\
    "teprim"        : "Normalized electron temperature gradient length (a/LTe)",\
    "fprim"         : "Normalized density gradient length (a/Ln)",\
    "tite"          : "Ion to electron temperature ratio",\
    "teti"          : "Electron to ion temperature ratio",\
    "delt"          : "Time step (dt)",\
    "delta t"       : "Time step (dt)",\
    "rk"            : "Order of the Runge Kutta scheme",\
    "nz"            : "Number of grid divisions along z (nzed)",\
    "nzed"          : "Number of grid divisions along z (nzed)",\
    "nzgrid"        : "Number of grid divisions along z (nvgrid)",\
    "nfield"        : "Number of field periods along z (nfield)",\
    "nperiod"       : "Number of field periods along z (nperiod)",\
    "pol.turns"     : "Number of poloidal turns",\
    "nmu"           : "Number of grid divisions along mu (nmu)",\
    "nvgrid"        : "Number of grid divisions along vpa (nvgrid)",\
    "dmu"           : "Smallest step along mu (dmu)",\
    "dvpa"          : "Step size along vpa (dvpa)",\
    "y0"            : "Box size along x (Ly/2pi)",\
    "nx"            : "Number of grid divisions along x (nx)",\
    "ny"            : "Number of grid divisions along y (ny)",\
    "kx max"        : "Maximum value of kx (kx max)",\
    "ky max"        : "Maximum value of ky (ky max",\
    "dkx"           : "Step size along kx (dkx)",\
    "dky"           : "Step size along ky (dky)",\
    "Lx"            : "Box size along x (Lx)",\
    "Ly"            : "Box size along y (Ly)",\
    "D_hyper"       : "Hyper-dissipation",\
    "alpha0"        : "Alpha",\
    "kappa"         : "Kappa",\
    "tri"           : "Triangularity",\
    "boundary_option" : "Boundary option",\
    "d_hyper"       : "D hyper",\
    "explicit_option" : "Runge Kutta",\
    "cfl_cushion"   : "CFL cushion",\
    "cfl_cushion_upper"   : "CFL cushion up",\
    "cfl_cushion_lower"   : "CFL cushion down",\
    # EMPTY
    None            : None,\
    "-----"         : None,\
    "-"             : None,\
    }


# Rewrite in scientific notation to be used in matplotlib titles
replacements = [ 
    ["(t)", "$t$"],\
    ["(z)", "$z$"],\
    ["(n)", "$n$"],\
    ["(T)", "$T$"],\
    ["(vpa)", "$v_\\parallel$"],\
    ["(zeta)", "$\\zeta$"],\
    ["(kx)", "$k_x$"],\
    ["(ky)", "$k_y$"],\
    ["(vpa)", "$v_\\parallel$"],\
    ["(mu)", "$\\mu$"],\
    ["(omega)", "$\\omega"],\
    ["(gamma)", "$\\gamma$"],\
    ["(gamma/ky**2)", "$\\gamma/k_y**2$"],\
    ["(Q)", "$Q$"],\
    ["(Gamma)", "$\\Gamma$"],\
    ["(Pi)", "$\\Pi$"],\
    ["(phi)", "$\\phi$"],\
    ["(Re(phi))", "Re$(\\phi)$"],\
    ["(Im(phi))", "Im$(\\phi)$"],\
    ["(Re(phi) and Im(phi))", "Re$(\\phi)$ and Im$(\\phi)$"],\
    ["(phi**2)", "$|\\phi^2|"],\
    ["(g)", "$g$"],\
    ["(rho)", "$\\rho$"],\
    ["(a/LTi)", "$a/L_{T_i}$"],\
    ["(a/LTe)", "$a/L_{T_e}$"],\
    ["(a/Ln)", "$a/L_n$"],\
    ["(dt)", "$\Delta t$"],\
    ["z (nzed)", "$z$ $(N_z)$"],\
    ["z (nzgrid)", "$z$ $(N_z)$"],\
    ["z (nfield)", "$z$ $(N_{fp})$"],\
    ["mu (nmu)", "$\\mu$ $(N_{\\mu})$"],\
    ["vpa (nvgrid)", "$v_\\parallel$ $(N_{v_\\parallel})$"],\
    ["mu (dmu)", "$\\mu$ $(\\mu_{min})$"],\
    ["vpa (dvpa)", "$v_\\parallel$ $(\Delta v_\\parallel)$"],\
    ["x (Ly/2pi)", "$x$ $(L_y/2*\\pi)"],\
    ["x (nx)", "$x$ $(N_x)$"],\
    ["y (ny)", "$y$ $(N_y)$"],\
    ["kx (kx max)", "$k_x$ $(k_{x,max})$"],\
    ["ky (ky max)", "$k_y$ $(k_{y,max})$"],\
    ["kx (dkx)", "$k_x$ $(\Delta k_x)$"],\
    ["ky (dky)", "$k_y$ $(\Delta k_y)$"],\
    ["x (Lx)", "$x$ $L_x$"],\
    ["y (Ly)", "$y$ $L_y$"]] 


# Rewrite in scientific notation to be used in matplotlib titles
standardNamesScientific = {}
for key in standardNames.keys():
    standardNamesScientific[key] = standardNames[key]
    if standardNamesScientific[key]!=None:
        for replacement in replacements:
            standardNamesScientific[key] = standardNamesScientific[key].replace(replacement[0], replacement[1])

