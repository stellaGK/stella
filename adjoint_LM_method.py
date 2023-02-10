
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 15:50:32 2022

@author: georgiaacton
"""

import numpy as np
import subprocess
import f90nml
import re

###### Calling overal optimisation

def read_miller_parameters ():
    nml = f90nml.read('example.in')
    # miller =  [nml['millergeo_parameters']['rhoc'], \
    #         nml['millergeo_parameters']['rmaj'], \
    #         nml['millergeo_parameters']['rgeo'], \
    #         nml['millergeo_parameters']['shift'], \
    #         nml['millergeo_parameters']['kappa'], \
    #         nml['millergeo_parameters']['kapprim'], \
    #         nml['millergeo_parameters']['qinp'], \
    #         nml['millergeo_parameters']['shat'], \
    #         nml['millergeo_parameters']['tri'], \
    #         nml['millergeo_parameters']['triprim'], \
    #         nml['millergeo_parameters']['betaprim'] ]

    miller = [nml['millergeo_parameters']['kappa'], \
            nml['millergeo_parameters']['tri'] ]

    delta_out = nml['millergeo_parameters']['del']    
    return miller, delta_out

def write_miller_parameters (p, po):
    nml = f90nml.read('example.in')
    nml['millergeo_parameters']['kappa'] = p[0]
    nml['millergeo_parameters']['tri'] = p[1]

    nml.indent = 1
    nml.write('example.in', force=True)
#    nml['millergeo_parameters']['rhoc'] = po[0]
#    nml['millergeo_parameters']['rmaj'] = po[1]
#    nml['millergeo_parameters']['rgeo'] = p[2]
#    nml['millergeo_parameters']['shift'] = p[3]
#    nml['millergeo_parameters']['kappa'] = p[4]
#    nml['millergeo_parameters']['kapprim'] = po[5]
#    nml['millergeo_parameters']['qinp'] = p[6]
#    nml['millergeo_parameters']['shat'] = p[7]
#    nml['millergeo_parameters']['tri'] = p[8]
#    nml['millergeo_parameters']['triprim'] = po[9]
#    nml['millergeo_parameters']['betaprim'] = po[10]
    
def empty_files (): 

    open('adjoint_files/adjoint_ginit.dat', 'w').close()
    open('adjoint_files/adjoint_gend.dat', 'w').close()
    open('adjoint_files/adjoint_omega.dat', 'w').close()
    open('adjoint_files/adjoint_derivatives.dat', 'w').close()
    open('adjoint_files/adjoint_final_time.dat', 'w').close()

def read_from_files ():

    with open('adjoint_files/adjoint_ginit.dat') as f:
        ginit_value = [float(x) for x in f]

    with open('adjoint_files/adjoint_gend.dat') as f:
        gfinal_value = [float(x) for x in f]
 
    with open('adjoint_files/adjoint_omega.dat') as f:
        omega_value = [[float(x) for x in line.split()] for line in f]

    with open('adjoint_files/adjoint_derivatives.dat') as f:
        gdt_value = [[float(x) for x in line.split()] for line in f]

    with open('adjoint_files/adjoint_final_time.dat') as f:
        time_value = [float(x) for x in f]        

    [p_value,delp_value] = read_miller_parameters ()
    ginit_value = np.array(ginit_value)
    gfinal_value = np.array(gfinal_value)
    omega_value = np.array(omega_value[0])
    time_value = np.array(time_value)
    for i in range (len(gdt_value)):
        gdt_value[i] = gdt_value[i][0]

    gdt_value = np.array(gdt_value)
    hess_value = np.multiply.outer(gdt_value, gdt_value)    
    return p_value, delp_value, ginit_value, gfinal_value, omega_value, time_value, \
        gdt_value, hess_value

def gradient_decent (p_in, grdt) :

    p_in = np.array(p_in)
    grdt = np.array(grdt)
    epsilon = 10**(-1)
    # use gradient decent to calculate next valkue of p                                                                                                       
    p_out = p_in - epsilon * grdt

    return p_out

def calling_routine (mu): 
    ## Read output values from stella simulation
    [p_old,delp, ginit, gfinal, omega, time, gdt, hess] = read_from_files()    
    print('p old=', p_old) 
    stop = False
    p_old = np.array(p_old)
    gdt = np.array(gdt)
    hess = np.array(hess)
    
    p_new = gradient_decent (p_old, gdt) 

#     k = 0
#     k_max = 10
#     tau = float(10**(-2))
#     epsilon1 = float(10**(-2))
#     epsilon2 = float(10**(-2))
#     epsilon3 = float(10**(-2))
    
#     l_up = float(11.0)
#     l_down = 9
    
#     # if first_loop:
#     #     p_new = LM_method (p_old, gdt, hess,delp, mu)
#     # else :
#     while (k <= k_max) and not (stop):
#         p_new = LM_method (p_old, gdt,hess, omega, mu)
# #        p_new = gradient_decent (p_old, gdt)
#         p_new = np.array(p_new)
#         print('pnew=', p_new)
#         ## if omega is negative then stop and update p
#         dp = p_new - p_old
#         diff = dp.max() - epsilon2*p_old.max()
#         rho = (0.5*np.dot(dp, np.dot(dp,hess)))/np.dot(dp,gdt)
#         print('rho=', rho)
#         if abs(rho)<epsilon3 :
#             print('rho <', epsilon3)
#             #                p_old = p_new
# #            mu = max(mu/l_down, 10**(-2))
#             stop = True 
#         else :
#             print('rho >', epsilon3)
#             p_new = p_old
#             mu = min(mu*l_up, 10**2)
#             stop = False
#             k = k + 1
#                 #            if diff.any() < 0.0 :
#                 #                stop = True
#                 #            else :

    write_miller_parameters(p_new, p_old)
    return p_new, omega, ginit, gfinal, gdt, time, mu

def LM_method (p_in, gdt_in, hess_in,omega_in,mu_in) : 
    p_in = np.array(p_in)
    dp = np.array(len(p_in))
    gdt_in = np.array(gdt_in)
    hess_in = np.array(hess_in)
    ## Matrix to invert
    A = hess_in + mu_in * np.identity(len(hess_in[0]))
    p_out = np.linalg.solve(A, np.dot(A,p_in)-gdt_in*omega_in)
    dp = np.linalg.solve(A, -gdt_in*omega_in)
    return p_out

def adjoint_loop (): 
   #### Call this file!!!! ####
    res = True
    it_max = 30
    empty_files()
    mu = float(10**(-3)) 
    omega_store = float(100.0)
    
    it = 0
    while(res == True) and it< it_max:
        mu = float(10**(-3))
        print ('it=', it)
        empty_files()
        ## calling stella
        subprocess.run("mpirun ./stella example.in", shell=True)
        #
        [p_updated,omega,gstart,gend,grad,time,mu] = calling_routine(mu)
        # if(time[0] == time[1]) :
        #     print('Omega did not converge in time limit')
        #     if(gend <= 0.1*gstart):
        #         res = False
        #     else:
        #         res = True
        # else :
        if (omega <= 0): 
            if (omega <= 0 and omega_store <=0): # or np.max(grad) < 0.01):
                p_updated = p_store
                res = False
            else:
                #mu = float(10**2)
                res = True 
        else :
            omega_store = omega
            p_store = p_updated
            res = True
            
        it = it + 1
        if(it >= it_max):
            print('Adjoint iterations exceed max')
            res = False

    return (omega, omega_store) 

def main_program (): 
    j_max = 3
    scale1 = float(1.2)
    scale2 = float(1.1)
    change_scale = float(0.5)
    
    j = 1
    while (j <= j_max) :
        nml = f90nml.read('example.in')
        tprim_old = nml['species_parameters_1']['tprim']

        [omega_new, omega_old] = adjoint_loop ()
        
        if omega_new <= 0: 
            tprim_new = tprim_old*scale1
        elif (omega_new < omega_old) :
            if (j == j_max) :
                tprim_new = tprim_old*scale2
            else : 
                stop = True
        else :
            scale2 = (scale2-1)/2 + 1
            tprim_new = tprim_old
            
        print('tprim', tprim_new, tprim_old)
        nml = f90nml.read('example.in')
        nml['species_parameters_1']['tprim'] = tprim_new
        nml.indent = 1
        nml.write('example.in', force=True)
        j = j + 1

if __name__ == "__main__":
#    main_program ()
    debug_code_1 = False
    [omega_new, omega_old] = adjoint_loop () 

    if(debug_code_1): 
        [p_value,delp_value] = read_miller_parameters ()
        p_new = p_value 
        write_miller_parameters(p_new, p_value)
