import matplotlib.pyplot as plt
import h5py

from numpy import *

data = './data.h5'

datafile = h5py.File(data, 'r')

test1   = datafile['compared/test1']
test1_5 = datafile['compared/test1_5']
test2   = datafile['compared/test2']
test3   = datafile['compared/test3']
test4   = datafile['compared/test4']
test5   = datafile['compared/test5']
 
omega_key='$a\\omega/v_{th,i}$'
gamma_key='$a\\gamma/v_{th,i}$'
kx_key='$k_x\\rho_i$'
ky_key='$k_y\\rho_i$'
t_key='$tv_{th,i}/a$'
pot_key='$\\langle\\mathrm{Re}(\\hat{\\varphi}_{\\mathbf{k}_{\\perp}})\\rangle_{z}/\\langle\\mathrm{Re}(\\hat{\\varphi}_{\\mathbf{k}_{\\perp}})\\rangle_{z}(t=0)$'
pot_norm_key='$|\\hat{\\varphi}_{\\mathbf{k}_{\\perp}}|/|\\hat{\\varphi}_{\\mathbf{k}_{\\perp}}|_{\\max}$'
Q_key='$Q_i/Q_{gB,i}$'
Q_kx_key = '$\\Sigma_{k_y}Q_i(k_x,k_y)/Q_{gB,i}$'
Q_ky_key = '$\\Sigma_{k_x}Q_i(k_x,k_y)/Q_{gB,i}$'


######################################################COMPARED DATA#########################################################

######################test1##############################

ky_stella_test1     = test1['stella/ky'][:]
gamma_stella_test1  = test1['stella/gamma'][:] 
omega_stella_test1  = test1['stella/omega'][:]

ky_GENE_test1       = test1['GENE/ky'][:]
gamma_GENE_test1    = test1['GENE/gamma'][:]
omega_GENE_test1  = test1['GENE/omega'][:]

fig=plt.figure(figsize=(15,8)) 
fig.suptitle('Test 1', fontsize=40)

ax_1 = fig.add_axes([0.1,0.15, 0.35, 0.7], ylim=(0,0.2),xlim=(0,4.5), yticks=[0,0.05,0.1,0.15])
ax_1.tick_params(labelsize=30)
ax_1.set_xlabel(ky_key, fontsize=30);  ax_1.set_ylabel(gamma_key, fontsize=30)
ax_1.xaxis.grid(color='grey', linestyle='-', linewidth=0.3); ax_1.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
print(ky_stella_test1)
print(gamma_stella_test1)
ax_1.plot(ky_stella_test1, gamma_stella_test1,linestyle='-',  marker='o', markersize=14, markerfacecolor="white",color='r', linewidth=5, label="stella")
ax_1.plot(ky_GENE_test1, gamma_GENE_test1 ,linestyle='--', marker='^', markersize=14, markerfacecolor="white",color='b', linewidth=5, label="GENE")  
ax_1.legend(loc='best',labelspacing=0.0, prop={'size':30})    

ax_2 = fig.add_axes([0.6,0.15, 0.35, 0.7], ylim=(0,0.6),xlim=(0,4.5), yticks=[0,0.2,0.4,0.6])
ax_2.tick_params(labelsize=30)
ax_2.set_xlabel(ky_key, fontsize=30);  ax_2.set_ylabel(omega_key, fontsize=30)
ax_2.xaxis.grid(color='grey', linestyle='-', linewidth=0.3); ax_2.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_2.plot(ky_stella_test1, omega_stella_test1,linestyle='-',  marker='o', markersize=14, markerfacecolor="white",color='r', linewidth=5, label="stella")
ax_2.plot(ky_GENE_test1, omega_GENE_test1, linestyle='--', marker='^', markersize=14, markerfacecolor="white",color='b', linewidth=5, label="GENE")
ax_2.legend(loc='best',labelspacing=0.0, prop={'size':30})

######################test1_5###########################

kx_stella_test1_5     = test1_5['stella/kx'][:]
gamma_stella_test1_5  = test1_5['stella/gamma'][:] 
omega_stella_test1_5  = test1_5['stella/omega'][:]

kx_GENE_test1_5       = test1_5['GENE/kx'][:]
gamma_GENE_test1_5    = test1_5['GENE/gamma'][:]
omega_GENE_test1_5  = test1_5['GENE/omega'][:]

fig=plt.figure(figsize=(15,8))
fig.suptitle('Test 1', fontsize=40)

ax_1 = fig.add_axes([0.1,0.15, 0.35, 0.7], ylim=(0,0.2),xlim=(0,3), yticks=[0.05,0.1,0.15])
ax_1.tick_params(labelsize=30)
ax_1.set_xlabel(kx_key, fontsize=30);  ax_1.set_ylabel(gamma_key, fontsize=30)
ax_1.xaxis.grid(color='grey', linestyle='-', linewidth=0.3); ax_1.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_1.plot(kx_stella_test1_5, gamma_stella_test1_5,linestyle='-',  marker='o', markersize=14, markerfacecolor="white",color='r', linewidth=5, label="stella")
ax_1.plot(kx_GENE_test1_5, gamma_GENE_test1_5,linestyle='--', marker='^', markersize=14, markerfacecolor="white",color='b', linewidth=5, label="GENE")  
ax_1.legend(loc='best',labelspacing=0.0, prop={'size':30})    

ax_2 = fig.add_axes([0.6,0.15, 0.35, 0.7], ylim=(0.2,0.35),xlim=(0,3), yticks=[0.25,0.3,0.35])
ax_2.tick_params(labelsize=30)
ax_2.set_xlabel(kx_key, fontsize=30);  ax_2.set_ylabel(omega_key, fontsize=30)
ax_2.xaxis.grid(color='grey', linestyle='-', linewidth=0.3); ax_2.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_2.plot(kx_stella_test1_5, omega_stella_test1_5,linestyle='-',  marker='o', markersize=14, markerfacecolor="white",color='r', linewidth=5, label="stella")
ax_2.plot(kx_GENE_test1_5, omega_GENE_test1_5, linestyle='--', marker='^', markersize=14, markerfacecolor="white",color='b', linewidth=5, label="GENE")
ax_2.legend(loc='best',labelspacing=0.0, prop={'size':30})

#####################test2#############################

kx_stella_test2     = test2['stella/kx'][:]
gamma_stella_test2  = test2['stella/gamma'][:] 
omega_stella_test2  = test2['stella/omega'][:]

kx_GENE_test2       = test2['GENE/kx'][:]
gamma_GENE_test2    = test2['GENE/gamma'][:]
omega_GENE_test2  = test2['GENE/omega'][:]

fig=plt.figure(figsize=(15,8)) 
fig.suptitle('Test 2', fontsize=40)

ax_1 = fig.add_axes([0.1,0.15, 0.35, 0.7], ylim=(0.05,0.2),xlim=(0.5,2.2), yticks=[0.05,0.1,0.15,0.2])
ax_1.tick_params(labelsize=30)
ax_1.set_xlabel(kx_key, fontsize=30);  ax_1.set_ylabel(gamma_key, fontsize=30)
ax_1.xaxis.grid(color='grey', linestyle='-', linewidth=0.3); ax_1.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_1.plot(kx_stella_test2, gamma_stella_test2,linestyle='-',  marker='o', markersize=14, markerfacecolor="white",color='r', linewidth=5, label="stella")
ax_1.plot(kx_GENE_test2, gamma_GENE_test2,linestyle='--', marker='^', markersize=14, markerfacecolor="white",color='b', linewidth=5, label="GENE")  
ax_1.legend(loc='best',labelspacing=0.0, prop={'size':30})   

ax_2 = fig.add_axes([0.6,0.15, 0.35, 0.7], ylim=(0.25,0.35),xlim=(0.5,2.2))
ax_2.tick_params(labelsize=30)
ax_2.set_xlabel(kx_key, fontsize=30);  ax_2.set_ylabel(omega_key, fontsize=30)
ax_2.xaxis.grid(color='grey', linestyle='-', linewidth=0.3); ax_2.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_2.plot(kx_stella_test2, omega_stella_test2,linestyle='-',  marker='o', markersize=14, markerfacecolor="white",color='r', linewidth=5, label="stella")
ax_2.plot(kx_GENE_test2, omega_GENE_test2, linestyle='--', marker='^', markersize=14, markerfacecolor="white",color='b', linewidth=5, label="GENE")
ax_2.legend(loc='best',labelspacing=0.0, prop={'size':30})

#####################test3#############################

ky_stella_test3     = test3['stella/ky'][:]
gamma_stella_test3  = test3['stella/gamma'][:] 
omega_stella_test3  = test3['stella/omega'][:]
z1_stella_test3     = test3['stella/branch_1/z'][:]
z2_stella_test3     = test3['stella/branch_2/z'][:]
z3_stella_test3     = test3['stella/branch_3/z'][:]
pot1_stella_test3   = test3['stella/branch_1/pot'][:]
pot2_stella_test3   = test3['stella/branch_2/pot'][:]
pot3_stella_test3   = test3['stella/branch_3/pot'][:]

ky_GENE_test3       = test3['GENE/ky'][:]
gamma_GENE_test3    = test3['GENE/gamma'][:]
omega_GENE_test3    = test3['GENE/omega'][:]
z1_GENE_test3       = test3['GENE/branch_1/z'][:]
z2_GENE_test3       = test3['GENE/branch_2/z'][:]
z3_GENE_test3       = test3['GENE/branch_3/z'][:]
pot1_GENE_test3     = test3['GENE/branch_1/pot'][:]
pot2_GENE_test3     = test3['GENE/branch_2/pot'][:]
pot3_GENE_test3     = test3['GENE/branch_3/pot'][:]

fig=plt.figure(figsize=(15,8)) 
fig.suptitle('Test 3', fontsize=40)

ax_1 = fig.add_axes([0.1,0.15, 0.35, 0.7], ylim=(0,0.4),xlim=(0,10), yticks=[0,0.1,0.2,0.3,0.4])
ax_1.tick_params(labelsize=30)
ax_1.set_xlabel(ky_key, fontsize=30);  ax_1.set_ylabel(gamma_key, fontsize=30)
ax_1.xaxis.grid(color='grey', linestyle='-', linewidth=0.3); ax_1.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_1.plot(ky_stella_test3, gamma_stella_test3,linestyle='-',  marker='o', markersize=14, markerfacecolor="white",color='r', linewidth=5, label="stella")
ax_1.plot(ky_GENE_test3, gamma_GENE_test3,linestyle='--', marker='^', markersize=14, markerfacecolor="white",color='b', linewidth=5, label="GENE")  
ax_1.legend(loc='best',labelspacing=0.0, prop={'size':30})       
   
ax_2 = fig.add_axes([0.6,0.15, 0.35, 0.7], ylim=(-1.5,0.5),xlim=(0,10), yticks=[-1.5,-1,-0.5,0,0.5])
ax_2.tick_params(labelsize=30)
ax_2.set_xlabel(ky_key, fontsize=30);  ax_2.set_ylabel(omega_key, fontsize=30)
ax_2.xaxis.grid(color='grey', linestyle='-', linewidth=0.3); ax_2.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_2.plot(ky_stella_test3, omega_stella_test3,linestyle='-',  marker='o', markersize=14, markerfacecolor="white",color='r', linewidth=5, label="stella")
ax_2.plot(ky_GENE_test3, omega_GENE_test3, linestyle='--', marker='^', markersize=14, markerfacecolor="white",color='b', linewidth=5, label="GENE")
ax_2.legend(loc='best',labelspacing=0.0, prop={'size':30})

fig_pot=plt.figure(figsize=(7,9))
fig_pot.suptitle('Test 3', fontsize=40)

ax_pot_1 = fig_pot.add_axes([0.25,0.68, 0.65, 0.22], ylim=(0,1),xlim=(-8*pi,8*pi), yticks=[0,0.5,1])
ax_pot_1.tick_params(labelsize=30)
ax_pot_1.set_xticks([-8*pi,-4*pi,0,4*pi,8*pi]); ax_pot_1.set_xticklabels(['$-8\\pi$','$-4\\pi$','$z$','$4\\pi$','$8\\pi$'])
ax_pot_1.set_ylabel(pot_norm_key, fontsize=30)
ax_pot_1.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_pot_1.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_pot_1.plot(z1_stella_test3, pot1_stella_test3, linestyle='-', color='r', linewidth=3, label='stella')
ax_pot_1.plot(z1_GENE_test3, pot1_GENE_test3, linestyle='--', color='b', linewidth=3, label='GENE')
ax_pot_1.legend(loc='best',labelspacing=0.0, prop={'size':17})

ax_pot_2 = fig_pot.add_axes([0.25,0.38, 0.65, 0.22], ylim=(0,1),xlim=(-4*pi,4*pi), yticks=[0,0.5,1])
ax_pot_2.tick_params(labelsize=30)
ax_pot_2.set_xticks([-4*pi,-2*pi,0,2*pi,4*pi]); ax_pot_2.set_xticklabels(['$-4\\pi$','$-2\\pi$','$z$','$2\\pi$','$4\\pi$'])
ax_pot_2.set_ylabel(pot_norm_key, fontsize=30)
ax_pot_2.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_pot_2.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_pot_2.plot(z2_stella_test3, pot2_stella_test3, linestyle='-', color='r', linewidth=3, label='stella')
ax_pot_2.plot(z2_GENE_test3, pot2_GENE_test3, linestyle='--', color='b', linewidth=3, label='GENE')
ax_pot_2.legend(loc='best',labelspacing=0.0, prop={'size':17})

ax_pot_3 = fig_pot.add_axes([0.25,0.08, 0.65, 0.22], ylim=(0,1),xlim=(-2*pi,2*pi), yticks=[0,0.5,1])
ax_pot_3.tick_params(labelsize=30)
ax_pot_3.set_xticks([-2*pi,-pi,0,pi,2*pi]); ax_pot_3.set_xticklabels(['$-2\\pi$','$-\\pi$','$z$','$\\pi$','$2\\pi$'])
ax_pot_3.set_ylabel(pot_norm_key, fontsize=30)
ax_pot_3.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_pot_3.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_pot_3.plot(z3_stella_test3, pot3_stella_test3, linestyle='-', color='r', linewidth=3, label='stella')
ax_pot_3.plot(z3_GENE_test3, pot3_GENE_test3, linestyle='--', color='b', linewidth=3, label='GENE')
ax_pot_3.legend(loc='best',labelspacing=0.0, prop={'size':17})

#########################test4###################################

t_stella_test4_case_1        = test4['case1/stella/t'][:]
pot_norm_stella_test4_case_1 = test4['case1/stella/phi'][:] 
t_stella_test4_case_2        = test4['case2/stella/t'][:]
pot_norm_stella_test4_case_2 = test4['case2/stella/phi'][:] 
t_stella_test4_case_3        = test4['case3/stella/t'][:]
pot_norm_stella_test4_case_3 = test4['case3/stella/phi'][:] 
t_stella_test4_case_4        = test4['case4/stella/t'][:]
pot_norm_stella_test4_case_4 = test4['case4/stella/phi'][:] 

t_GENE_test4_case_1        = test4['case1/GENE/t'][:]
pot_norm_GENE_test4_case_1 = test4['case1/GENE/phi'][:] 
t_GENE_test4_case_2        = test4['case2/GENE/t'][:]
pot_norm_GENE_test4_case_2 = test4['case2/GENE/phi'][:] 
t_GENE_test4_case_3        = test4['case3/GENE/t'][:]
pot_norm_GENE_test4_case_3 = test4['case3/GENE/phi'][:] 
t_GENE_test4_case_4        = test4['case4/GENE/t'][:]
pot_norm_GENE_test4_case_4 = test4['case4/GENE/phi'][:]

fig_pot=plt.figure(figsize=(15,8))
fig_pot.suptitle('Test 4', fontsize=40)

ax_1 = fig_pot.add_axes([0.1,0.53, 0.35, 0.35], ylim=(-0.5,1),xlim=(0,3500), yticks=[-0.4,-0.2,0,0.2,0.4,0.6,0.8,1])
ax_1.tick_params(labelsize=20)
ax_1.set_xticks([0,500,1500,2500,3500])
ax_1.set_xlabel('$tv_{th,i}/a$', fontsize=20);  ax_1.set_ylabel(pot_key, fontsize=18)
ax_1.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_1.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_1.plot(t_stella_test4_case_1, pot_norm_stella_test4_case_1, linestyle='-', color='r', linewidth=3, label="stella")
ax_1.plot(t_GENE_test4_case_1, pot_norm_GENE_test4_case_1, linestyle='--', color='b', linewidth=3, label="GENE")
ax_1.legend(loc=4,labelspacing=0.0, prop={'size':19})

ax_2 = fig_pot.add_axes([0.55,0.53, 0.35, 0.35], ylim=(-0.5,1),xlim=(0,2000), yticks=[-0.4,-0.2,0,0.2,0.4,0.6,0.8,1])
ax_2.tick_params(labelsize=20) 
ax_2.set_xticks([0,500,1000,1500,2000])
ax_2.set_xlabel('$tv_{th,i}/a$', fontsize=20);  ax_2.set_ylabel(pot_key, fontsize=18)
ax_2.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_2.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_2.plot(t_stella_test4_case_2, pot_norm_stella_test4_case_2, linestyle='-', color='r', linewidth=3, label="stella")
ax_2.plot(t_GENE_test4_case_2, pot_norm_GENE_test4_case_2, linestyle='--', color='b', linewidth=3, label="GENE")
ax_2.legend(loc=4,labelspacing=0.0, prop={'size':19})

ax_3 = fig_pot.add_axes([0.1,0.08, 0.35, 0.35], ylim=(-0.5,1),xlim=(0,2000), yticks=[-0.4,-0.2,0,0.2,0.4,0.6,0.8,1])
ax_3.tick_params(labelsize=20)
ax_3.set_xticks([0,500,1000,1500,2000])
ax_3.set_xlabel('$tv_{th,i}/a$', fontsize=20);  ax_3.set_ylabel(pot_key, fontsize=18)
ax_3.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_3.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_3.plot(t_stella_test4_case_3, pot_norm_stella_test4_case_3, linestyle='-', color='r', linewidth=3, label="stella")
ax_3.plot(t_GENE_test4_case_3, pot_norm_GENE_test4_case_3, linestyle='--', color='b', linewidth=3, label="GENE")
ax_3.legend(loc=4,labelspacing=0.0, prop={'size':19})
    
ax_4 = fig_pot.add_axes([0.55,0.08, 0.35, 0.35], ylim=(-0.5,1),xlim=(0,2000), yticks=[-0.4,-0.2,0,0.2,0.4,0.6,0.8,1])
ax_4.tick_params(labelsize=20)
ax_4.set_xticks([0,500,1000,1500,2000])
ax_4.set_xlabel('$tv_{th,i}/a$', fontsize=20);  ax_4.set_ylabel(pot_key, fontsize=18)
ax_4.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_4.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_4.plot(t_stella_test4_case_4, pot_norm_stella_test4_case_4, linestyle='-', color='r', linewidth=3, label="stella")
ax_4.plot(t_GENE_test4_case_4, pot_norm_GENE_test4_case_4, linestyle='--', color='b', linewidth=3, label="GENE")
ax_4.legend(loc=4,labelspacing=0.0, prop={'size':19})

#####################test5########################################

t_stella_test5          = test5['stella/t'][:]
Q_stella_test5          = test5['stella/Q'][:]
kx_stella_test5         = test5['stella/kx'][:]
ky_stella_test5         = test5['stella/ky'][:]
Q_vs_kx_stella_test5    = test5['stella/Q_vs_kx'][:]
Q_vs_ky_stella_test5    = test5['stella/Q_vs_ky'][:]

t_GENE_test5          = test5['GENE/t'][:]
Q_GENE_test5          = test5['GENE/Q'][:]
kx_GENE_test5         = test5['GENE/kx'][:]
ky_GENE_test5         = test5['GENE/ky'][:]
Q_vs_kx_GENE_test5    = test5['GENE/Q_vs_kx'][:]
Q_vs_ky_GENE_test5    = test5['GENE/Q_vs_ky'][:]
    

fig=plt.figure(figsize=(15,8))
fig.suptitle('Test 5', fontsize=40)

ax=fig.add_axes([0.1,0.15, 0.85, 0.77], xlim=(0,1900), ylim=(0,8), yticks=[0,1,2,3,4,5,6,7,8])
ax.tick_params(labelsize=30)
ax.set_xlabel(t_key, fontsize=30);  ax.set_ylabel(Q_key, fontsize=30)
ax.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax.plot(t_stella_test5, Q_stella_test5, linestyle='-',color='r', linewidth=5, label='stella')
ax.plot(t_GENE_test5, Q_GENE_test5, linestyle='--',color='b', linewidth=5, label='GENE')
ax.legend(loc='best',labelspacing=0.0, prop={'size':30})

fig=plt.figure(figsize=(15,8)) 
fig.suptitle('Test 5', fontsize=40)


ax_1 = fig.add_axes([0.1,0.15, 0.35, 0.7], ylim=(0,0.2),xlim=(0,6))
ax_1.tick_params(labelsize=30)
ax_1.set_xlabel(ky_key, fontsize=30); ax_1.set_ylabel(Q_ky_key, fontsize=30)
ax_1.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_1.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_1.plot(ky_stella_test5, Q_vs_ky_stella_test5, linestyle='-', marker='o', markersize=14, markerfacecolor='white', color='r', linewidth=4, label='stella')
ax_1.plot(ky_GENE_test5, Q_vs_ky_GENE_test5, linestyle='--', marker='^', markersize=14, markerfacecolor='white', color='b', linewidth=4, label='GENE')      
ax_1.legend(loc='best',labelspacing=0.0, prop={'size':20})

ax_2 = fig.add_axes([0.6,0.15, 0.35, 0.7], ylim=(0,0.2),xlim=(-2.2,2.2))
ax_2.tick_params(labelsize=30)
ax_2.set_xlabel(kx_key, fontsize=30); ax_2.set_ylabel(Q_kx_key, fontsize=30)
ax_2.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_2.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax_2.plot(kx_stella_test5, Q_vs_kx_stella_test5, linestyle='-', marker='o', markersize=14, markerfacecolor='white', color='r', linewidth=4, label='stella')
ax_2.plot(kx_GENE_test5, Q_vs_kx_GENE_test5, linestyle='--', marker='^', markersize=14, markerfacecolor='white', color='b', linewidth=4, label='GENE')      
ax_2.legend(loc='best',labelspacing=0.0, prop={'size':30})


######################################################NON COMPARED DATA#########################################################

######################################test1##################

test1_nc        = datafile['non_compared/test1']

z1_test1        = test1_nc['z_first_branch'][:]
z2_test1        = test1_nc['z_second_branch'][:]
pot1_test1      = test1_nc['pot_first_branch'][:]
pot2_test1      = test1_nc['pot_second_branch'][:]


fig=plt.figure(figsize=(15,8))
fig.suptitle('Test 1 (not compared)', fontsize=40)

ax=fig.add_axes([0.1,0.15, 0.85, 0.77], xlim=(-12*pi/5,12*pi/5), ylim=(0,1), yticks=[0,0.5,1])
ax.set_xticks([-2*pi,-pi,0,pi,2*pi]); ax.set_xticklabels(['$-2\\pi$','$-\\pi$', '$0$', '$\\pi$', '$2\\pi$'])
ax.tick_params(labelsize=30)
ax.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax.set_xlabel('$z$', fontsize=30); ax.set_ylabel(pot_norm_key, fontsize=30)
ax.plot(z1_test1,pot1_test1,linestyle='-', color='limegreen', linewidth=5, label='$(k_x\\rho_i, k_y\\rho_i)=(0, 1.3)$')
ax.plot(z2_test1,pot2_test1,linestyle='-', color='k', linewidth=5, label='$(k_x\\rho_i, k_y\\rho_i)=(0, 2.1)$')
ax.legend(loc=4,labelspacing=0.0, prop={'size':30})

######################################test1_5##################

test1_5_nc        = datafile['non_compared/test1_5']

z1_test1_5        = test1_5_nc['z_first_branch'][:]
z2_test1_5        = test1_5_nc['z_second_branch'][:]
pot1_test1_5     = test1_5_nc['pot_first_branch'][:]
pot2_test1_5     = test1_5_nc['pot_second_branch'][:]


fig=plt.figure(figsize=(15,8))
fig.suptitle('Test 1 (not compared)', fontsize=40)

ax=fig.add_axes([0.1,0.15, 0.85, 0.77], xlim=(-6*pi,6*pi), ylim=(0,1), yticks=[0,0.5,1])
ax.set_xticks([-6*pi,-3*pi,0,3*pi,6*pi]); ax.set_xticklabels(['$-6\\pi$','$-3\\pi$', '$0$', '$3\\pi$', '$6\\pi$'])
ax.tick_params(labelsize=30)
ax.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax.set_xlabel('$z$', fontsize=30); ax.set_ylabel(pot_norm_key, fontsize=30)
ax.plot(z1_test1_5,pot1_test1_5,linestyle='-', color='limegreen', linewidth=5, label='$(k_x\\rho_i, k_y\\rho_i)=(0.2, 2.1)$')
ax.plot(z2_test1_5,pot2_test1_5,linestyle='-', color='k', linewidth=5, label='$(k_x\\rho_i, k_y\\rho_i)=(1.7, 2.1)$')
ax.legend(loc='best',labelspacing=0.0, prop={'size':30})

######################################test2##################

test2_nc        = datafile['non_compared/test2']

z_test2        = test2_nc['z'][:]
pot_test2     = test2_nc['pot'][:]


fig=plt.figure(figsize=(15,8))
fig.suptitle('Test 2 (not compared)', fontsize=40)

ax=fig.add_axes([0.1,0.15, 0.85, 0.77], xlim=(-3*pi,3*pi), ylim=(0,1), yticks=[0,0.5,1])
ax.set_xticks([-3*pi,-2*pi,-pi,0,pi,2*pi,3*pi]); ax.set_xticklabels(['$-3\\pi$','$-2\\pi$','$-\\pi$', '$0$', '$\\pi$', '$2\\pi$', '$3\\pi$'])
ax.tick_params(labelsize=30)
ax.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)
ax.set_xlabel('$z$', fontsize=30); ax.set_ylabel(pot_norm_key, fontsize=30)
ax.plot(z_test2, pot_test2,linestyle='-', color='limegreen', linewidth=5, label='$(k_x\\rho_i, k_y\\rho_i)=(1.2, 2.1)$')
ax.legend(loc='best',labelspacing=0.0, prop={'size':30})






plt.show()
