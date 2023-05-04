
def load_stellaKnobsAndKeys():
    ''' List all the possible stella knobs and the corresponding keys. '''
    
    stellaKnobsAndKeys = {
        "-" : \
            ['-'],\
        "zgrid_parameters" : \
            ['nz', 'nzed', 'nperiod', 'ntubes', 'boundary_option', \
             'zed_equal_arc', 'shat_zero', 'nzgrid'],\
        "geo_knobs" : \
            ['geo_option', 'overwrite_bmag', 'overwrite_gradpar', 'overwrite_gds2', \
             'overwrite_gds21', 'overwrite_gds22', 'overwrite_gds23', 'overwrite_gds24', \
             'overwrite_gbdrift', 'overwrite_cvdrift', 'overwrite_gbdrift0', 'geo_file'],\
        "vmec_parameters" : \
            ['vmec_filename', 'alpha0', 'zeta_center', 'nfield_periods', 'torflux', \
             'surface_option', 'verbose', 'zgrid_scalefac', 'zgrid_refinement_factor'],\
        "parameters" : \
            ['beta', 'vnew_ref', 'rhostar', 'zeff', 'tite', 'nine' ],\
        "vpamu_grids_parameters" : 
            ['nvgrid', 'vpa_max', 'nmu', 'vperp_max', 'equally_spaced_mu_grid', 'dmu', 'dvpa'],\
        "dist_fn_knobs" : 
            ['adiabatic_option'],\
        "time_advance_knobs" : 
            ['explicit_option', 'xdriftknob', 'ydriftknob', 'wstarknob', 'flip_flo'],\
        "kt_grids_knobs" : 
            ['grid_option'],\
        "kt_grids_box_parameters" : 
            ['kx max', 'ky max', 'nx', 'ny', 'dkx', 'dky', 'jtwist', 'y0', 'naky', 'nakx'],\
        "kt_grids_range_parameters" : 
            ['nalpha', 'naky', 'nakx', 'aky_min', 'aky_max', 'akx_min', \
             'akx_max', 'theta0_min', 'theta0_max'],\
        "physics_flags" : 
            ['full_flux_surface', 'include_mirror', 'nonlinear', 'include_apar', \
             'include_parallel_nonlinearity', 'include_parallel_streaming'],\
        "init_g_knobs" : 
            ['tstart', 'scale', 'ginit_option', 'width0', 'refac', 'imfac', \
             'den0', 'par0', 'tperp0', 'den1', 'upar1', 'tpar1', 'tperp1', \
             'den2', 'upar2', 'tpar2', 'tperp2', 'phiinit', 'zf_init', 'chop_side',\
             'left', 'even', 'restart_file', 'restart_dir', 'read_many'],\
        "knobs" : 
            ['t_end', 'nstep', 'delt', 'delt_option', 'zed_upwind', \
             'vpa_upwind', 'time_upwind', 'avail_cpu_time', 'cfl_cushion',\
             'delt_adjust', 'mat_gen', 'mat_read', 'fields_kxkyz', 'stream_implicit', \
             'mirror_implicit', 'driftkinetic_implicit',  'mirror_semi_lagrange', \
             'mirror_linear_interp', 'maxwellian_inside_zed_derivative', 'stream_matrix_inversion'],\
        "species_knobs" : 
            ['nspec', 'species_option'],\
        "species_parameters_1" : 
            ['z', 'mass', 'dens', 'tprim', 'fprim', 'd2ndr2', 'd2Tdr2', 'type' ],\
        "species_parameters_2" : 
            ['z', 'mass', 'dens', 'tprim', 'fprim', 'd2ndr2', 'd2Tdr2', 'type' ],\
        "species_parameters_3" : 
            ['z', 'mass', 'dens', 'tprim', 'fprim', 'd2ndr2', 'd2Tdr2', 'type' ],\
        "species_parameters_4" : 
            ['z', 'mass', 'dens', 'tprim', 'fprim', 'd2ndr2', 'd2Tdr2', 'type' ],\
        "stella_diagnostics_knobs" : 
            ['nwrite', 'navg', 'nmovie', 'nsave', 'save_for_restart', 'write_omega', \
             'write_phi_vs_time', 'write_gvmus', 'write_gzvs', 'write_kspectra', \
             'write_moments', 'flux_norm', 'write_fluxes_kxky' ],\
        "millergeo_parameters" : 
            ['rhoc', 'rmaj', 'shift', 'qinp', 'shat', 'kappa', 'kapprim', 'tri',\
             'triprim', 'rgeo', 'betaprim', 'betadbprim', 'd2qdr2', 'd2psidr2', \
             'nzed_local', 'read_profile_variation', 'write_profile_variation'],\
        "layouts_knobs" : 
            ['xyzs_layout', 'vms_layout'],\
        "neoclassical_input" : 
            ['include_neoclassical_terms', 'nradii', 'drho', 'neo_option'],\
        "sfincs_input" : \
            ['read_sfincs_output_from_file', 'nproc_sfincs', 'irad_min', 'irad_max', 'calculate_radial_electric_field',\
            'includeXDotTerm', 'includeElectricFieldTermInXiDot', 'magneticDriftScheme', 'includePhi1',\
            'includePhi1InKineticEquation', 'geometryScheme', 'VMECRadialOption', 'equilibriumFile',\
            'coordinateSystem', 'inputRadialCoordinate', 'inputRadialCoordinateForGradients', 'aHat',\
            'psiAHat', 'Delta', 'nu_n', 'dPhiHatdrN', 'Er_window', 'nxi', 'nx', 'Ntheta', 'Nzeta'],\
        "dissipation" : \
           ['include_collisions', 'include_krook_operator', 'collisions_implicit',\
           'collision_model', 'momentum_conservation', 'energy_conservation',\
           'vpa_operator', 'mu_operator', 'hyper_dissipation', 'remove_zero_projection',\
           'D_hyper', 'nu_krook', 'delay_krook', 'ikxmax_source', 'krook_odd', 'cfac']}
    return stellaKnobsAndKeys
        
