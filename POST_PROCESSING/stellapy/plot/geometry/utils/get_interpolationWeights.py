 
from stellapy.data.geometry.calculate_geometricQuantitiesVMEC import vmec 
    
def get_interpolationWeightsOnHalfGrid(ns, svalue):
    """" Get the interpolation weights, for when svalue lies between two flux surfaces. """
    VMEC = vmec()  
    VMEC.alpha0 = 0
    VMEC.normalized_toroidal_flux_used = svalue
    VMEC.calculate_normalizedToroidalFluxOnFullAndHalfGrids(ns)
    VMEC.calculate_interpolationWeights(ns)  
    return VMEC.radial_index_half, VMEC.radial_weight_half 