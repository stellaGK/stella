

#===============================================================================
#                              GET MODES IDENTIFIER                            #
#=============================================================================== 
  
def get_modeIdentifier(kx, ky): 
    """ There are issues with kx = 0 being printed as kx = -0.0. """
    if kx==0: kx = 0.0
    if ky==0: ky = 0.0 
    return '('+str(kx)+', '+str(ky)+')', kx, ky
