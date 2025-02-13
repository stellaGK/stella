
# Variables in the namelist
text = 'include_parallel_streaming, include_mirror, nonlinear, adiabatic_option, prp_shear_enabled, hammett_flow_shear, include_pressure_variation, include_geometric_variation, include_parallel_nonlinearity, suppress_zonal_interaction, full_flux_surface, include_apar, include_bpar, radial_variation, xdriftknob, ydriftknob, wstarknob, beta, zeff, tite, nine, rhostar, vnew_ref, g_exb, g_exbfac, omprimfac, irhostar'
text = 'fphi, delt, nstep, tend, delt_option, lu_option, autostop, avail_cpu_time, delt_max, delt_min, cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower, stream_implicit, mirror_implicit, drifts_implicit, use_deltaphi_for_response_matrix, maxwellian_normalization, stream_matrix_inversion, maxwellian_inside_zed_derivative, mirror_semi_lagrange, mirror_linear_interp, zed_upwind, vpa_upwind, time_upwind, &fields_kxkyz, mat_gen, mat_read, rng_seed, ky_solve_radial, ky_solve_real, nitt, print_extra_info_to_terminal'
text = 'nwrite, navg, nsave, save_for_restart, write_phi_vs_time, write_gvmus, write_gzvs, write_omega, write_kspectra, write_moments, write_radial_fluxes, write_radial_moments, write_fluxes_kxkyz, flux_norm, nc_mult, write_apar_vs_time, write_bpar_vs_time'
text = 'naky, nakx, aky_min, aky_max, theta0_min, theta0_max, akx_min, akx_max, kyspacing_option'
text = 'nx, ny, jtwist, jtwistfac, x0, y0, centered_in_rho, periodic_variation, randomize_phase_shift, phase_shift_angle'

# Get the names of the variables
variables = text.split(' ')
variables = [i.replace(',','') for i in variables]
variables = [i for i in variables if i!= '']
variables = [i.replace(',','') for i in variables]
variables_str = ", ".join(variables)
print("\nVariables:")
print(variables_str)
print()

# List of new variables
new_variables = [i+'_new' for i in variables]
new_variables_str = ", ".join(new_variables)
print("\nNew variables:")
print(new_variables_str)
print()

# Save the variables
save_variables = [i+'_new = '+i for i in variables]
save_variables_str = "\n".join(save_variables)
print("\nSave variables:")
print(save_variables_str)
print()

# Check the variables
check_variables = [f'if ({i}_new /= {i}) double_definitions = .true.' for i in variables]
check_variables_str = "\n".join(check_variables)
print("\Check variables:")
print(check_variables_str)
print()


            
