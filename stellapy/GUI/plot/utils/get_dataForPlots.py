
import numpy as np 
 
#===============================================================================
#                       GET THE (X,Y) DATA FOR THE PLOTS                       #
#===============================================================================
 
def get_dataForPlots(plot, simulation, quantity, specie=None): 
    
    # Read the data from the correct file/object  
    if "zonal" not in plot.yname:
        if "omega_vs"   in quantity:    data = getattr(simulation.omega, quantity)
        elif "gamma_vs" in quantity:    data = getattr(simulation.omega, quantity) 
        elif "omega_"   in quantity:    data = getattr(simulation.lineardata, quantity)
        elif "gamma_"   in quantity:    data = getattr(simulation.lineardata, quantity)
        elif "flux"     in quantity:    data = getattr(simulation.fluxes, quantity)
        elif "phi_"     in quantity:    data = getattr(simulation.potential, quantity)
        elif "phi2_"    in quantity:    data = getattr(simulation.potential, quantity)
        elif "g_"       in quantity:    data = getattr(simulation.distribution, quantity)
    if "zonal" in plot.yname: 
        if "phi2_"      in quantity:    data = getattr(simulation.potential, quantity+"_"+plot.yname.split("_")[-1])

    # Read the (x,y) or (x,y,z) data  
    if plot.dimensions==2: x, y = data.get_xyData()
    if plot.dimensions==3: x, y, _ = data.get_xyzData() 
    
    # Get a specific species 
    if "s" in data.dimensions and specie!=None and specie!="both":
        axis = data.dimensions.index("s")
        y = y.take(indices=specie, axis=axis) 
    
    # Sometimes we convert the basic dimensions
    if plot.xdim in ["pol", "tor", "zeta"]:
        x = getattr(simulation.vec, plot.xdim)
        
    # Sometimes we convert the basic dimensions
    if plot.xdim in ["pol", "tor", "zeta"]:
        x = getattr(simulation.vec, plot.xdim)
    
    # Get the real or imaginary part of the data 
    if "_real" in plot.yname: y = y.real/np.max(np.abs(y)) 
    if "_imag" in plot.yname: y = y.imag/np.max(np.abs(y))  

    # Rescale the data to SI units or to compare the shape
    if plot.units=="SI": 
        x = x*simulation.referenceunits.ref[data.dimensions[0]]
        y = y*simulation.referenceunits.ref[plot.yname]  
        
    # Only consider the shape of the quantity
    if plot.normalize and plot.units!="SI": 
        if data.dimensions[0]=="t":    x = x/simulation.time.peakTime
        if data.quantity=="qflux":     y = y/simulation.time.saturatedFluxes["qflux"][specie]
        if data.quantity=="pflux":     y = y/simulation.time.saturatedFluxes["pflux"][specie]
        if data.quantity=="vflux":     y = y/simulation.time.saturatedFluxes["vflux"][specie]  
    
    # Put the y-data in logarithmic scales
    if plot.takelogx: x = np.log10(x)
    if plot.takelogy: y = np.log10(y) 
    
    # Python is limited to numbers within (-1.79769313486e+308, 1.79769313486e+308)
    if np.any(y>1.E270) and not np.all(y>1.E270):
        print("WARNING: The data contains numbers that are bigger than 1.E270")
        x = x[x<1.E270]
        y = y[y<1.E270]
                
    # Return the (x,y) data   
    return x, y 

#===============================================================================
#                           GET MULTIPLE QUANTITIES                            #
#===============================================================================
 
def get_yQuantities(quantity):
 
    # We can have combined graphs (phiReal, phiImag) (Phi2, Phi2_nozonal, Phi2_zonal) 
    # and (Phi2_nozonal, qflux) for which we will iterate over certain y-quantities
    if quantity in ["phiRealImag", "phi2_split", "phi2AndQFlux", "distribution"]:
         
        # Define the y-quantities
        if quantity=="phiRealImag":  y_quantities = ["phi_real", "phi_imag"]
        if quantity=="phi2_split":   y_quantities = ["phi2", "phi2_zonal", "phi2_nozonal"]
        if quantity=="phi2AndQFlux": y_quantities = ["phi2_nozonal", "qflux"]
        if quantity=="distribution": y_quantities = ["gvmus", "gzvs"] 
  
    # If we only have 1 quantity, keep the original data
    else: y_quantities = [quantity] 
    return y_quantities






