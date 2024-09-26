"""

#===============================================================================
#                           CALCULATE QUASINEUTRALITY                          #
#===============================================================================

This function returns the values of the equilibrium density as well as its 
gradient from quasineutrality to lowest order: 
    
    Sum(zs * ns) = 0
    Sum(zs * dns * ns / ne) = 0 

It assumes three species in the plasma: 

    main ions : with charge zi and normalized density gradient dnbulk if bulk='i'
    electrons : with charge -1 and normalized density gradient dnbulk if bulk='e'
    impurities: with charge zz and normalized density gradient dnz

It also considers the value of the effective charge zeff. If main ions are set 
as the reference species (bulk='i'), the electron density and its gradient are 
calculated. If electrons are set as the reference species (bulk='e'), the main 
ion density and its gradient are calculated instead.

"""
    
#===============================================================================
#                           CALCULATE QUASINEUTRALITY                          #
#===============================================================================

def calculate_quasiNeutrality(zi=1, zz=6,zeff=1.05, nbulk=1.0, dnbulk=1.0, dnz=0.0, 
        bulk='i', tite=1.0, verbose=False):
      
    nzz_over_ne = (zeff - zi)/zz/(zz - zi)
    nzi_over_ne = (zz - zeff)/zi/(zz - zi)
    ne_over_nzi = 1/nzi_over_ne
    
    # if bulk == 'i': the line below is what we need
    dne_if_dni_fixed = zi*dnbulk*nzi_over_ne + zz*dnz*nzz_over_ne
    
    # elif bulk == 'e':the line below is what we need
    dni_if_dne_fixed = (dnbulk - zz*dnz*nzz_over_ne)*ne_over_nzi/zi

    if bulk == 'i':
        nzi = nbulk
        ne  = ne_over_nzi * nzi
        nzz = nzz_over_ne * ne
        dni = dnbulk
        dne = dne_if_dni_fixed
    elif bulk == 'e':
        ne  = nbulk
        nzi = nzi_over_ne * ne
        nzz = nzz_over_ne * ne
        dni = dni_if_dne_fixed
        dne = dnbulk
    
    if verbose:
        print('\nZ_eff                   = ', "%14.6e" % zeff)
        print('Z_i                     = ', "%14.6e" % zi)
        print('Z_Z                     = ', "%14.6e" % zz)
        print('nine                    = ', "%14.6e" % (nzi/ne))
        print('dens_i                  = ', "%14.6e" % (nzi))
        print('dens_e                  = ', "%14.6e" % (ne))
        print('dens_z                  = ', "%14.6e" % (nzz))
        print('fprim_i                 = ', "%14.6e" % (dni))
        print('fprim_e                 = ', "%14.6e" % (dne))
        print('fprim_i                 = ', "%14.6e" % (dnz))
        print('\n******* Taking n_e as reference **********')
        print('(n_e/n_e, a/L_ne)       = (',"%14.6e" % nbulk, ', ',"%14.6e" % dnbulk,')')
        print('(n_i/n_e, a/L_ni)       = (',"%14.6e" % nzi_over_ne,', ',"%14.6e" % dni_if_dne_fixed,')')
        print('(n_Z/n_e, a/L_nZ)       = (',"%14.6e" % nzz_over_ne,', ',"%14.6e" % dnz,')')
        print('\n******* Taking n_i as reference **********')
        print('(n_e/n_i, a/L_ne)       = (',"%14.6e" % ne_over_nzi, ', ',"%14.6e" % dne_if_dni_fixed,')')
        print('(n_i/n_i, a/L_ni)       = (',"%14.6e" % nbulk,', ',"%14.6e" % dnbulk,')')
        print('(n_Z/n_i, a/L_nZ)       = (',"%14.6e" % (nzz_over_ne*ne_over_nzi) , ', ',"%14.6e" % dnz,')')
        print('\n******* Orderings *****************')
        print('(Z_Z*n_Z)/(n_e)         = ', "%14.6e" % (zz * nzz_over_ne))
        print('(Z_Z^2*n_Z)/(Z_i^2*n_i) = ', "%14.6e" % (zz**2 * nzz_over_ne / zi**2 / nzi_over_ne))
        print('\n********* Threshold for impurity mode: n\'z/n\'i < threshold *******')
        print('(Zi^2*n_i+ne*Ti/Te)/Zz*Zi*n_i = ', "%14.6e" % ((zi**2)*nzi+ne*tite/zz*zi*nzi),'\n ')
        print('\n******* Quasineutrality checks ****')
        print('Sum(z_s * n_s)          = ', "%14.6e" %  (zi * nzi - ne + zz * nzz))
        print('Sum(z_s*dn_s*n_s /n_e)  = ', "%14.6e" % (zi * dni * (nzi/ne) - dne + zz * dnz * (nzz/ne)))
        print('\n')
          
    return nzi, ne, nzz, dni, dne, dnz, nzi/ne


