################################################################################
#               PARAMETERS THAT CAN BE CHOSEN IN THE GUI PLOTS
################################################################################

# Get the stella <parameter>, <key> and <knob> for the interesting parameters
standardParameters = {}
standardParameters["boundary_option"]   = {"knob" : "zgrid_parameters", "key" : "boundary_option"} 
standardParameters["tite"]   = {"knob" : "parameters",              "key" : "tite"} 
standardParameters["teti"]   = {"knob" : "parameters",              "key" : "teti"} 
standardParameters["rho"]    = {"knob" : "vmec_parameters",         "key" : "rho"} 
standardParameters["tprim"]  = {"knob" : "species_parameters_1",    "key" : "tprim"} 
standardParameters["tiprim"] = {"knob" : "species_parameters_1",    "key" : "tprim"}  
standardParameters["teprim"] = {"knob" : "species_parameters_2",    "key" : "tprim"}   
standardParameters["fprim"]  = {"knob" : "species_parameters_1",    "key" : "fprim"}   
standardParameters["delt"]   = {"knob" : "knobs",                   "key" : "delt"}  
standardParameters["delta t"]= {"knob" : "knobs",                   "key" : "delt"}  
standardParameters["nmu"]    = {"knob" : "vpamu_grids_parameters",  "key" : "nmu"} 
standardParameters["nvgrid"] = {"knob" : "vpamu_grids_parameters",  "key" : "nvgrid"} 
standardParameters["dvpa"]   = {"knob" : "vpamu_grids_parameters",  "key" : "dvpa"} 
standardParameters["dmu"]    = {"knob" : "vpamu_grids_parameters",  "key" : "dmu"} 
standardParameters["nz"]     = {"knob" : "zgrid_parameters",        "key" : "nz"} 
standardParameters["nzed"]   = {"knob" : "zgrid_parameters",        "key" : "nzed"} 
standardParameters["nzgrid"] = {"knob" : "zgrid_parameters",        "key" : "nzgrid"} 
standardParameters["nx"]     = {"knob" : "kt_grids_box_parameters", "key" : "nx"}
standardParameters["ny"]     = {"knob" : "kt_grids_box_parameters", "key" : "ny"}
standardParameters["y0"]     = {"knob" : "kt_grids_box_parameters", "key" : "y0"}
standardParameters["kx max"] = {"knob" : "kt_grids_box_parameters", "key" : "kx max"}
standardParameters["ky max"] = {"knob" : "kt_grids_box_parameters", "key" : "ky max"}
standardParameters["dkx"]    = {"knob" : "kt_grids_box_parameters", "key" : "dkx"}
standardParameters["dky"]    = {"knob" : "kt_grids_box_parameters", "key" : "dky"}
standardParameters["Lx"]     = {"knob" : "kt_grids_box_parameters", "key" : "Lx"}
standardParameters["Ly"]     = {"knob" : "kt_grids_box_parameters", "key" : "Ly"}
standardParameters["rk"]     = {"knob" : "time_advance_knobs",      "key" : "explicit_option"}
standardParameters["nfield"] = {"knob" : "vmec_parameters",         "key" : "nfield_periods"}
standardParameters["pol.turns"] = {"knob" : "vmec_parameters",      "key" : "poloidal_turns"}
standardParameters["nperiod"]= {"knob" : "zgrid_parameters",        "key" : "nperiod"}
standardParameters["D_hyper"]= {"knob" : "dissipation",             "key" : "D_hyper"}
standardParameters["-"]      = {"knob" : "-",                       "key" : "-"}
standardParameters["-----"]  = {"knob" : "-",                       "key" : "-"} 

# Define the order for the dropdown menu in the GUI
standardParametersInOrder = [
    "-", "tiprim", "teprim", "fprim", "-----",\
    "delta t", "rk", "D_hyper", "-----",\
    "rho", "nz", "nzed", "nzgrid", "nfield", "nperiod",  "-----",\
    "nmu", "nvgrid", "dmu", "dvpa", "-----",\
    "y0", "nx", "ny", "kx max", "ky max", "dkx", "dky", "Lx", "Ly"]  