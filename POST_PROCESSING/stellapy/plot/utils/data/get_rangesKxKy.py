 

#===============================================================================
#                                GET (KX,KY) RANGES                            #
#===============================================================================

def get_rangesKxKy(args):
    """ Convert {kxmin, kxmax, kymin, kymax} to {kx_range, ky_range} """

    # Get the kx-range
    if "kxmin" in args and "kxmax" in args: args["kx_range"] = [args["kxmin"],args["kxmax"]]; del args["kxmin"]; del args["kxmax"]
    elif "kxmin" in args: args["kx_range"] = [args["kxmin"],99999]; del args["kxmin"]
    elif "kxmax" in args: args["kx_range"] = [-99999,args["kxmax"]]; del args["kxmax"]

    # Get the ky-range
    if "kymin" in args and "kymax" in args: args["ky_range"] = [args["kymin"],args["kymax"]]; del args["kymin"]; del args["kymax"]
    elif "kymin" in args: args["ky_range"] = [args["kymin"],99999]; del args["kymin"]
    elif "kymax" in args: args["ky_range"] = [0,args["kymax"]]; del args["kymax"]
    return args
