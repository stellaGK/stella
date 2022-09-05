
# Save the netcdf name, and its dimensions
stella_variables = {
    
    # SYNONYMS OF OLDER STELLA VERSIONS
    "phi_vs_t" : ["phi_vs_t",        ["t", "tube", "z", "kx", "ky", "ri"]], 
    
    # DIMENSIONS WITH VECTORS (1 dim)
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
    
    
    # VALUES: Dimensions
    "nkx"      : ["nkx",        ["-"]], 
    "nky"      : ["nky",        ["-"]], 
    "nmu"      : ["nmu",        ["-"]], 
    "nvpa_tot" : ["nvpa_tot",   ["-"]], 
    "nzed_tot" : ["nzed_tot",   ["-"]], 
    "nspecies" : ["nspecies",   ["-"]], 
    "ntubes"   : ["ntubes",     ["-"]], 
    # VALUES: Simulations
    "nproc"    : ["nproc",      ["-"]], 
    "nmesh"    : ["nmesh",      ["-"]], 
    # VALUES: Geometry
    "q"        : ["q",          ["-"]], 
    "beta"     : ["beta",       ["-"]], 
    "shat"     : ["shat",       ["-"]], 
    "jtwist"   : ["jtwist",     ["-"]], 
    "drhodpsi" : ["drhodpsi",   ["-"]], 
    
    
    # VECTORS: dimensions
    "vec_time" : ["t",          ["t"]], 
    "vec_kx"   : ["kx",         ["kx"]], 
    "vec_ky"   : ["ky",         ["ky"]], 
    "vec_mu"   : ["mu",         ["mu"]], 
    "vec_vpa"  : ["vpa",        ["vpa"]], 
    "vec_zed"  : ["zed",        ["z"]], 
    # VECTORS versus (species)
    "charge"   : ["charge",     ["s"]], 
    "mass"     : ["mass",       ["s"]], 
    "dens"     : ["dens",       ["s"]], 
    "temp"     : ["temp",       ["s"]], 
    "tprim"    : ["tprim",      ["s"]], 
    "fprim"    : ["fprim",      ["s"]], 
    "vnew"     : ["vnew",       ["s"]], 
    "type_of_species" : ["type_of_species", ["s"]], 
    # VECTORS versus (t)
    "phi2_vs_t": ["phi2",       ["t"]], 
    # VECTORS versus (z)
    "gradpar"  : ["gradpar",    ["z"]],     
    
    
    # VECTORS versus (z,alpha)
    "bmag"     : ["bmag",       ["z", "alpha"]], 
    "theta0"   : ["theta0",     ["z", "alpha"]], 
    "gbdrift"  : ["gbdrift",    ["z", "alpha"]], 
    "gbdrift0" : ["gbdrift0",   ["z", "alpha"]], 
    "cvdrift"  : ["cvdrift",    ["z", "alpha"]], 
    "cvdrift0" : ["cvdrift0",   ["z", "alpha"]], 
    "gds2"     : ["gds2",       ["z", "alpha"]], 
    "gds21"    : ["gds21",      ["z", "alpha"]], 
    "gds22"    : ["gds22",      ["z", "alpha"]], 
    "grho"     : ["grho",       ["z", "alpha"]], 
    "jacob"    : ["jacob",      ["z", "alpha"]], 
    
    
    # VECTORS versus (t,theta0,ky)
    "phi2_vs_tkxky"   : ["phi2_vs_kxky",    ["t", "kx", "ky"]], 
    # VECTORS versus (zed,alpha,theta0,ky)
    "kperp2"          : ["kperp2",          ["z", "tube", "kx", "ky"]], 
    # VECTORS versus (t,species,mu,vpa)
    "g_vs_tsmuvpa"    : ["gvmus",           ["t", "s", "mu", "vpa"]], 
    "g2_vs_tsmuvpa"    : ["gvmus",          ["t", "s", "mu", "vpa"]], 
    # VECTORS versus (t,species,theta0,ky)
    "pflux_vs_tskxky" : ["pflx_kxky",       ["t", "s", "kx", "ky"]], 
    "vflux_vs_tskxky" : ["vflx_kxky",       ["t", "s", "kx", "ky"]], 
    "qflux_vs_tskxky" : ["qflx_kxky",       ["t", "s", "kx", "ky"]], 
    # VECTORS versus (t,species,vpa,zed,tube)
    "g_vs_tsvpaz"     : ["gzvs",            ["t", "s", "vpa", "z", "tube"]], 
    "g2_vs_tsvpaz"    : ["gzvs",            ["t", "s", "vpa", "z", "tube"]], 
    # VECTORS versus (t,tube,zed,theta0,ky,ri)
    "phi_vs_tzkxkyri" : ["phi_vs_t",        ["t", "tube", "z", "kx", "ky", "ri"]], 
    # VECTORS versus (t,species,tube,z,theta0,ky)
    "pflux_vs_tszkxky" : ["pflx_kxky",      ["t", "s", "tube", "z", "kx", "ky"]], 
    "vflux_vs_tszkxky" : ["vflx_kxky",      ["t", "s", "tube", "z", "kx", "ky"]], 
    "qflux_vs_tszkxky" : ["qflx_kxky",      ["t", "s", "tube", "z", "kx", "ky"]], 
    "trappedw_vs_tszkxky" : ["trappedw_vs_kxkyz", ["t", "s", "tube", "z", "kx", "ky"]],
    # VECTORS versus (t,species,tube,zed,theta0,ky,ri)
    "dens_vs_tszkxkyri" : ["density",       ["t", "s", "tube", "z", "kx", "ky", "ri"]], 
    "upar_vs_tszkxkyri" : ["upar",          ["t", "s", "tube", "z", "kx", "ky", "ri"]],
    "temp_vs_tszkxkyri" : ["temperature",   ["t", "s", "tube", "z", "kx", "ky", "ri"]],
    
    
    # OTHER VARIABLES
    "code_info"   : ["code_info",    ["-"]], 
    "input_file"  : ["input_file",   ["-"]], 
    
    
    # OLD STELLAPY VERSION KEYS
    }


# Save the netcdf name, and its dimensions
oldstella_variables = {
    "pflux_vs_tskxky" : "pflx_kxky",
    "vflux_vs_tskxky" : "vflx_kxky",    
    "qflux_vs_tskxky" : "qflx_kxky",      
    "g_vs_tsmuvpa"    : "gvmus",           
    "g_vs_tsvpaz"     : "gzvs",  
    "phi2_vs_t"       : "phi2_vs_tkxky",
    "phi2_vs_tkxky"   : "phi2_vs_kxky"
    }