
################################################################################
#               Check all the terms in the gyrokinetic equation                #
################################################################################
# These test are designed to test the terms in the gyrokinetic equation one 
# by one. We test the following terms:
#     - nonlinear term                 (nonlinear)
#     - parallel streaming explicit    (include_parallel_streaming, stream_implicit)
#     - parallel streaming implicit    (include_parallel_streaming, stream_implicit)
#     - mirror explicit                (include_mirror, mirror_implicit)
#     - mirror implicit                (include_mirror, mirror_implicit)
#     - magnetic drifts explicit       (xdriftknob and ydriftknob)
#     - diamagnetic drifts explicit    (wstarknob)
#     - drifts implicit                (xdriftknob, ydriftknob, wstarknob, include_parallel_streaming, drifts_implicit)
# 
# The relevant parameters in the input file are:
#       &physics_flags 
#         nonlinear = .false.
#         include_parallel_streaming = .false.
#         include_mirror = .false.
#       /
#       &time_advance_knobs
#         xdriftknob = 0
#         ydriftknob = 0
#         wstarknob = 0 
#       /
#       &dissipation
#         hyper_dissipation = .false.
#       /  
#       &knobs 
#         stream_implicit = .false.
#         mirror_implicit = .false.
#         drifts_implicit = .false.
#       /
################################################################################

# Python modules
import os, sys
import pathlib 
import numpy as np
import xarray as xr  

# Package to run stella 
module_path1 = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
module_path2 = str(pathlib.Path(__file__).parent.parent.parent / 'create_figure.py')
with open(module_path1, 'r') as file: exec(file.read())
with open(module_path2, 'r') as file: exec(file.read())

#-------------------------------------------------------------------------------
#                                NONLINEAR TERM                                #
#-------------------------------------------------------------------------------
def test_rosenbluth_hinton():

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('Rosenbluth-Hinton.in', tmp_path)
    return
    
    
    
#-------------------------------------------------------------------------------
#                              LOCAL PYTHON SCRIPT                             #
#-------------------------------------------------------------------------------
if __name__ == '__main__': 

	# Import modules 
	import matplotlib.pyplot as plt
	import h5py
	
	# Input file
	input_file = 'Rosenbluth-Hinton.in'
	
	# Create directory for a local run
	current_directory = pathlib.Path(__file__).parent  
	local_run_directory = current_directory / 'local_run'
	os.makedirs(local_run_directory, exist_ok=True)
	
	# Check whether the simulation has run
	if not os.path.isfile(local_run_directory / 'Rosenbluth-Hinton.final_fields'): 
		run_local_stella_simulation(input_file, local_run_directory)
		

	# Read benchmark data from DOI:10.1017/S0022377822000393
	published_benchmark_data_path = current_directory.parent / 'DOI_10_1017_S0022377822000393.h5'  
	published_benchmark_data = h5py.File(published_benchmark_data_path, 'r')
	test4 = published_benchmark_data['compared/test4']
	data = {'case1' : {}, 'case2' : {}, 'case3' : {}, 'case4' : {}}
	for case in ['case1', 'case2', 'case3', 'case4']:
		data[case] = {'stella' : {}, 'GENE' : {}}
		data[case]['stella'] = {'t' : test4[f'{case}/stella/t'][:], 'phi' : test4[f'{case}/stella/phi'][:]}
		data[case]['GENE'] = {'t' : test4[f'{case}/GENE/t'][:], 'phi' : test4[f'{case}/GENE/phi'][:]}
		
	# Local stella run
	final_fields_path = local_run_directory / input_file.replace('.in', '.final_fields')
	final_fields = np.loadtxt(final_fields_path, dtype='float', skiprows=2).reshape(-1, 13) 
	phi_real = final_fields[:, 4]
	phi_imag = final_fields[:, 5]
	phi = np.abs(phi_real + 1j*phi_imag)
	z = final_fields[:, 0]
	
    # Plot 
	#ax = create_figure(xlabel='$z$', ylabel='$\\varphi/$max$(\\varphi)$')
	#ax.plot(z, phi/np.max(phi), lw=3, color='navy')
	#ax.set_xlim([np.min(z), np.max(z)])
	#ax.set_ylim(top=1)
	#plt.show() 
		
	# Create a new figure if no existing axis was given   
	font_size = read_fontsize(); set_latex_font()
	fig = plt.figure(figsize=read_figsize())
	fig.grid_specifications = gridspec.GridSpec(2, 2, figure=fig)
	fig.grid_specifications.update(top=0.9, left=0.1, right=0.95, bottom=0.1, hspace=0.4, wspace=0.3)
	ax1 = plt.subplot(fig.grid_specifications[0]) 
	ax2 = plt.subplot(fig.grid_specifications[1]) 
	ax3 = plt.subplot(fig.grid_specifications[2]) 
	ax4 = plt.subplot(fig.grid_specifications[3]) 
	axes = [ax1, ax2, ax3, ax4] 
	fig.suptitle('Test 4 in DOI:10.1017/S0022377822000393: Rosenbluth Hinton', fontsize=40)
    
	# Labels
	time_label = '$tv_{th,i}/a$'
	potential_label = '$\\langle\\mathrm{Re}(\\hat{\\varphi}_{\\mathbf{k}_{\\perp}})\\rangle_{z}'
	potential_label += '/\\langle\\mathrm{Re}(\\hat{\\varphi}_{\\mathbf{k}_{\\perp}})\\rangle_{z}(t=0)$'

	# Compare the benchmark data with the local stella run
	for i, case in enumerate(['case1', 'case2', 'case3', 'case4']):
		axes[i].set_xlabel(time_label)
		axes[i].set_ylabel(potential_label, labelpad=15)
		axes[i].plot(data[case]['stella']['t'], data[case]['stella']['phi'], linestyle='-', color='r', linewidth=3, label="stella (benchmark)")
		axes[i].plot(data[case]['GENE']['t'], data[case]['GENE']['phi'], linestyle='--', color='b', linewidth=3, label="GENE (benchmark)")
		axes[i].legend(loc='best', labelspacing=0.0, prop={'size':font_size})
		axes[i].set_xlim([0,3500])
		
	# Mode of each test
	axes[0].set_title('$k_x\\rho_i = 0.05, k_y\\rho_i = 0.0$')
	axes[1].set_title('$k_x\\rho_i = 0.07, k_y\\rho_i = 0.0$')
	axes[2].set_title('$k_x\\rho_i = 0.10, k_y\\rho_i = 0.0$')
	axes[3].set_title('$k_x\\rho_i = 0.30, k_y\\rho_i = 0.0$')
	update_figure_style(fig=fig, axes=axes, font_size=font_size, scientific_axis=False)
	

	# Show figure
	plt.show() 
    
    
    
    
    
    
    
    
    
    
   
