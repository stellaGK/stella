'''

Link stellapy names with stella names.

'''

# Save the netcdf name, and its dimensions
stella_variables = {
    
    # DIMENSIONS WITH Vectors (1 dim)
    "vec_kx"   : ["kx",         ["kx"]],
    "vec_ky"   : ["ky",         ["ky"]],
    "vec_z"    : ["zed",        ["z"]],
    "vec_time" : ["t",          ["t"]],
    "vec_vpa"  : ["vpa",        ["vpa"]],
    "vec_mu"   : ["mu",         ["mu"]],
    "theta0"   : ["theta0",     ["kx"]],
    "vec_species" : ["type_of_species", ["s"]],
    
    # DIMENSIONS WITH ONLY A SIZE
    "species"  : ["species",    ["s"]],
    "tube"     : ["tube",       ["tube"]],
    "alpha"    : ["alpha",      ["alpha"]],
    "ri"       : ["ri",         ["ri"]],
    
    # OTHER DIMENSIONS
    "char10"   : ["char10",     ["-"]],
    "char200"  : ["char200",    ["-"]],
    "nlines"   : ["nlines",     ["-"]],
    
    # OTHER VARIABLES
    "code_info"   : ["code_info",    ["-"]],
    "input_file"  : ["input_file",   ["-"]],
    
    
    # Values: Dimensions
    "nkx"      : ["nkx",        ["-"]],
    "nky"      : ["nky",        ["-"]],
    "nmu"      : ["nmu",        ["-"]],
    "nvpa_tot" : ["nvpa_tot",   ["-"]],
    "nzed_tot" : ["nzed_tot",   ["-"]],
    "nspecies" : ["nspecies",   ["-"]],
    "ntubes"   : ["ntubes",     ["-"]],
    
    # Values: Parallelisation
    "nproc"    : ["nproc",      ["-"]],
    "nmesh"    : ["nmesh",      ["-"]],
    
    # Values: Geometry
    "q"        : ["q",          ["-"]],
    "beta"     : ["beta",       ["-"]],
    "shat"     : ["shat",       ["-"]],
    "jtwist"   : ["jtwist",     ["-"]],
    "drhodpsi" : ["drhodpsi",   ["-"]],
    
    
    # Vectors: dimensions
    "vec_time" : ["t",          ["t"]],
    "vec_kx"   : ["kx",         ["kx"]],
    "vec_ky"   : ["ky",         ["ky"]],
    "vec_mu"   : ["mu",         ["mu"]],
    "vec_vpa"  : ["vpa",        ["vpa"]],
    "vec_zed"  : ["zed",        ["z"]],
    
    # Vectors versus (species)
    "charge"   : ["charge",     ["s"]],
    "mass"     : ["mass",       ["s"]],
    "dens"     : ["dens",       ["s"]],
    "temp"     : ["temp",       ["s"]],
    "tprim"    : ["tprim",      ["s"]],
    "fprim"    : ["fprim",      ["s"]],
    "vnew"     : ["vnew",       ["s"]],
    "type_of_species" : ["type_of_species", ["s"]],
    
    #---------------------------------------------------------------------------
    #                                  Geometry                                 
    #---------------------------------------------------------------------------
    
    # Vectors versus (z)
    "gradpar"  : ["gradpar",    ["z"]],
    
    # Vectors versus (z,alpha) with old geometry names
    "bmag"   : ["bmag",    ["z", "alpha"]],
    "theta0" : ["theta0",  ["z", "alpha"]],
    "grho"   : ["grho",    ["z", "alpha"]],
    "jacob"  : ["jacob",   ["z", "alpha"]],
    "b_dot_gradz:"    : ["b_dot_gradz:",    ["z", "alpha"]],
    "grady_dot_grady" : ["grady_dot_grady", ["z", "alpha"]],
    "gradx_dot_grady" : ["gradx_dot_grady", ["z", "alpha"]],
    "gradx_dot_gradx" : ["gradx_dot_gradx", ["z", "alpha"]],
    "B_times_gradB_dot_grady" : ["B_times_gradB_dot_grady",   ["z", "alpha"]],
    "B_times_gradB_dot_gradx" : ["B_times_gradB_dot_gradx",   ["z", "alpha"]],
    "B_times_kappa_dot_grady" : ["B_times_kappa_dot_grady",   ["z", "alpha"]],
    "B_times_kappa_dot_gradx" : ["B_times_kappa_dot_gradx",   ["z", "alpha"]],
    
    # Vectors versus (zed,alpha,theta0,ky)
    "kperp2"   : ["kperp2",          ["z", "tube", "kx", "ky"]],
    
    #---------------------------------------------------------------------------
    #                                  Potential                                
    #---------------------------------------------------------------------------
    
    # Oscillation frequency and growth rate
    "omega_vs_tkxky"  : ["omega",           ["t", "kx", "ky", "ri"]],
    
    # Perturbed electrostatic potential
    "phi2_vs_t"       : ["phi2",            ["t"]],
    "phi2_vs_tkxky"   : ["phi2_vs_kxky",    ["t", "kx", "ky"]],
    "phi_vs_tzkxkyri" : ["phi_vs_t",        ["t", "tube", "z", "kx", "ky", "ri"]],
    "phi_vs_t"        : ["phi_vs_t",        ["t", "tube", "z", "kx", "ky", "ri"]],
    
    #---------------------------------------------------------------------------
    #                                  Moments                                  
    #---------------------------------------------------------------------------
    
    # Moments
    "dens_vs_tszkxkyri" : ["density",       ["t", "s", "tube", "z", "kx", "ky", "ri"]],
    "upar_vs_tszkxkyri" : ["upar",          ["t", "s", "tube", "z", "kx", "ky", "ri"]],
    "temp_vs_tszkxkyri" : ["temperature",   ["t", "s", "tube", "z", "kx", "ky", "ri"]],
    
    #---------------------------------------------------------------------------
    #                                  Fluxes                                   
    #---------------------------------------------------------------------------
    
    # Fluxes
    "pflux_vs_ts" : ["pflux_vs_s",  ["t", "s"]],
    "vflux_vs_ts" : ["vflux_vs_s",  ["t", "s"]],
    "qflux_vs_ts" : ["qflux_vs_s",  ["t", "s"]],
    "pflux_vs_tskxky" : ["pflux_vs_kxkys",  ["t", "s", "kx", "ky"]],
    "vflux_vs_tskxky" : ["vflux_vs_kxkys",  ["t", "s", "kx", "ky"]],
    "qflux_vs_tskxky" : ["qflux_vs_kxkys",  ["t", "s", "kx", "ky"]],
    "pflux_vs_tszkxky" : ["pflux_vs_kxkyzs",  ["t", "s", "tube", "z", "kx", "ky"]],
    "vflux_vs_tszkxky" : ["vflux_vs_kxkyzs",  ["t", "s", "tube", "z", "kx", "ky"]],
    "qflux_vs_tszkxky" : ["qflux_vs_kxkyzs",  ["t", "s", "tube", "z", "kx", "ky"]],
    
    # Only for radial variation runs
    "pflux_vs_tskx" : ["pflux_x",  ["t", "s", "kx"]],
    "vflux_vs_tskx" : ["vflux_x",  ["t", "s", "kx"]],
    "qflux_vs_tskx" : ["qflux_x",  ["t", "s", "kx"]],
    
    #---------------------------------------------------------------------------
    #                           Distribution function                           
    #---------------------------------------------------------------------------
    
    # Vectors versus (t,species,mu,vpa)
    "g2_vs_tsmuvpa"   : ["g2_vs_vpamus",   ["t", "s", "mu", "vpa"]],
    "h2_vs_tsmuvpa"   : ["h2_vs_vpamus",   ["t", "s", "mu", "vpa"]],
    "f2_vs_tsmuvpa"   : ["f2_vs_vpamus",   ["t", "s", "mu", "vpa"]],
    "g2nozonal_vs_tsmuvpa"   : ["g2nozonal_vs_vpamus",  ["t", "s", "mu", "vpa"]],
    "h2nozonal_vs_tsmuvpa"   : ["h2nozonal_vs_vpamus",  ["t", "s", "mu", "vpa"]],
    "f2nozonal_vs_tsmuvpa"   : ["f2nozonal_vs_vpamus",  ["t", "s", "mu", "vpa"]],
    
    # Vectors versus (t,species,vpa,zed,tube)
    "g_vs_tsvpaz"     : ["g2_vs_zvpas",    ["t", "s", "vpa", "z", "tube"]],
    "g2_vs_tsvpaz"    : ["g2_vs_zvpas",    ["t", "s", "vpa", "z", "tube"]],
    "h2_vs_tsvpaz"    : ["h2_vs_zvpas",    ["t", "s", "vpa", "z", "tube"]],
    "f2_vs_tsvpaz"    : ["f2_vs_zvpas",    ["t", "s", "vpa", "z", "tube"]],
    "g2nozonal_vs_tsvpaz"    : ["g2nozonal_vs_zvpas",    ["t", "s", "vpa", "z", "tube"]],
    "h2nozonal_vs_tsvpaz"    : ["h2nozonal_vs_zvpas",    ["t", "s", "vpa", "z", "tube"]],
    "f2nozonal_vs_tsvpaz"    : ["f2nozonal_vs_zvpas",    ["t", "s", "vpa", "z", "tube"]],
    
    # Vectors versus (t,species,mu,zed,tube)
    "g_vs_tsmuz"     : ["g2_vs_zmus",      ["t", "s", "mu", "z", "tube"]],
    "g2_vs_tsmuz"    : ["g2_vs_zmus",      ["t", "s", "mu", "z", "tube"]],
    "h2_vs_tsmuz"    : ["h2_vs_zmus",      ["t", "s", "mu", "z", "tube"]],
    "f2_vs_tsmuz"    : ["f2_vs_zmus",      ["t", "s", "mu", "z", "tube"]],
    "g2nozonal_vs_tsmuz"    : ["g2nozonal_vs_zmus",      ["t", "s", "mu", "z", "tube"]],
    "h2nozonal_vs_tsmuz"    : ["h2nozonal_vs_zmus",      ["t", "s", "mu", "z", "tube"]],
    "f2nozonal_vs_tsmuz"    : ["f2nozonal_vs_zmus",      ["t", "s", "mu", "z", "tube"]],
    
    # Vectors versus (t,species,kx,ky,z,tube)
    "g2_vs_tskxkyz"   : ["g2_vs_zkykxs",    ["t", "s", "kx", "ky", "z", "tube"]],
    "h2_vs_tskxkyz"   : ["h2_vs_zkykxs",    ["t", "s", "kx", "ky", "z", "tube"]],
    "f2_vs_tskxkyz"   : ["f2_vs_zkykxs",    ["t", "s", "kx", "ky", "z", "tube"]],
    
    # Vectors versus (t,species,mu,vpa,tube,z)
    "g2_vs_tsmuvpaz"  : ["g2_vs_zvpamus",   ["t", "s", "mu", "vpa", "tube", "z"]],
    "h2_vs_tsmuvpaz"  : ["h2_vs_zvpamus",   ["t", "s", "mu", "vpa", "tube", "z"]],
    "f2_vs_tsmuvpaz"  : ["f2_vs_zvpamus",   ["t", "s", "mu", "vpa", "tube", "z"]],
    "g2nozonal_vs_tsmuvpaz"  : ["g2nozonal_vs_zvpamus",   ["t", "s", "mu", "vpa", "tube", "z"]],
    "h2nozonal_vs_tsmuvpaz"  : ["h2nozonal_vs_zvpamus",   ["t", "s", "mu", "vpa", "tube", "z"]],
    "f2nozonal_vs_tsmuvpaz"  : ["f2nozonal_vs_zvpamus",   ["t", "s", "mu", "vpa", "tube", "z"]],
    
    }

