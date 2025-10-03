
################################################################################
#                   Plot the electrostatic potential squared                   #
################################################################################
#
# Plot the electrostatic potential squared |phi|^2 of a stella simulation
#
# Hanne Thienpondt
# 03/10/2025
################################################################################

# Python modules
import glob
import h5py
import os, sys
import pathlib 
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Plot the electrostatic potential squared
def plot_phi(folder, fontsize=20, figsize=(18,9)):

    # Find the input files with output files
    input_files = [pathlib.Path(i) for i in glob.glob(str(folder/'*.in'))]
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix(".out"))]
    
    # Iterate over the input files
    for input_file in input_files:

        # Read the output file
        out_file = input_file.with_suffix(".out")
        data = np.loadtxt(out_file,skiprows=2,dtype='float').reshape(-1, 5)
        vec_time = data[:,1]
        phi2_vs_time = data[:,2]
        
        # Normalize the potential to the first time step
        phi2_vs_time = phi2_vs_time/phi2_vs_time[0]

        # Pretty figuers
        params = {'text.usetex' : True, 'font.size' : fontsize, 'font.family' : 'lmodern'}
        plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
        plt.rcParams['axes.unicode_minus'] = True
        plt.rcParams.update(params)

        # Create a figure
        fig = plt.figure(figsize=figsize)
        fig.grid_specifications = gridspec.GridSpec(1, 1, figure=fig)
        fig.grid_specifications.update(top=0.9, left=0.1, right=0.95, bottom=0.15)
        ax = plt.subplot(fig.grid_specifications[0]) 
        fig.set_tight_layout(False)
        ax.fig = fig; fig.ax = ax

        # Plot the time trace of the potential
        ax.plot(vec_time, phi2_vs_time, color='black', label="current stella data")
        
        # Add benchmark data
        plot_phi_from_paper(ax, case=1)
        
        # Finish the figure
        #ax.set_xlim([0, np.max(vec_time)])
        #ax.set_ylim([0, np.max(phi2_vs_time)])
        ax.set_xlabel("$t\\, v_{\\mathrm{th},r}/a$")
        ax.set_ylabel("$\\hat\\varphi$")
        ax.legend(loc=4,labelspacing=0.0, prop={'size':fontsize})

        # Show the figure
        figure_name = input_file.parent/input_file.name.replace('.in','_potential_vs_time.png')
        plt.savefig(figure_name, format='png', dpi=300)
        print(f'\n  ---> Saved {figure_name.name}\n')
        
    return
    
#-------------------------------------------------------------------------------
def plot_phi_from_paper(ax, case=1):
    ''' The benchmark data can be found in the attachments of [2022 - Gonzalez - Electrostatic 
    gyrokinetic simulations in Wendelstein 7-X geometry benchmark between the codes stella and GENE]'''
    
    directory = pathlib.Path(os.path.dirname(os.path.realpath(__file__))) 
    data = directory / 'published_data.h5'
    datafile = h5py.File(data, 'r')
    test4   = datafile['compared/test4']

    # Load published stella data
    if case==1:
        vec_time = test4['case1/stella/t'][:]
        phi2_vs_time = test4['case1/stella/phi'][:] 
    if case==2:
        vec_time = test4['case2/stella/t'][:]
        phi2_vs_time = test4['case2/stella/phi'][:] 
    if case==3:
        vec_time = test4['case3/stella/t'][:]
        phi2_vs_time = test4['case3/stella/phi'][:] 
    if case==4:
        vec_time = test4['case4/stella/t'][:]
        phi2_vs_time = test4['case4/stella/phi'][:] 

    # Load published GENE data
    if case==1:
        vec_time_GENE = test4['case1/GENE/t'][:]
        phi2_vs_time_GENE = test4['case1/GENE/phi'][:] 
    if case==2:
        vec_time_GENE = test4['case2/GENE/t'][:]
        phi2_vs_time_GENE = test4['case2/GENE/phi'][:] 
    if case==3:
        vec_time_GENE = test4['case3/GENE/t'][:]
        phi2_vs_time_GENE = test4['case3/GENE/phi'][:] 
    if case==4:
        vec_time_GENE = test4['case4/GENE/t'][:]
        phi2_vs_time_GENE = test4['case4/GENE/phi'][:]
        
    # Change from GENE to stella coordinates
    # vec_time = vec_time/np.sqrt(2)
    # vec_time_GENE = vec_time_GENE/np.sqrt(2)

    # Plot the data
    ax.plot(vec_time, phi2_vs_time, linestyle='-', color='r', linewidth=3, label="published stella data")
    ax.plot(vec_time_GENE, phi2_vs_time_GENE, linestyle='--', color='b', linewidth=3, label="GENE")

    # Finish the figure
    tmax = 3400 if case==1 else 2000
    if case==1: ax.set_xticks([0,500,1500,2500,3500])
    if case!=1: ax.set_xticks([0,500,1000,1500,2000])
    ax.set_xlim([0, tmax])
    ax.set_ylim([-0.5, 1])
    return

#===============================================================================
#                         RUN FROM THE COMMAND PROMPT                          #
#===============================================================================
 
if __name__ == "__main__":
   folder = pathlib.Path(os.getcwd())
   plot_phi(folder)
