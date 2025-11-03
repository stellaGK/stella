################################################################################
#               PARAMETERS THAT CAN BE CHOSEN IN THE GUI PLOTS
################################################################################

# Get the stella <parameter>, <key> and <knob> for the interesting parameters
standardParameters = {}
standardParameters['boundary_option']   = {'knob' : 'z_boundary_condition', 'key' : 'boundary_option'}
standardParameters['tite']   = {'knob' : 'adiabatic_electron_response',              'key' : 'tite'}
standardParameters['teti']   = {'knob' : 'adiabatic_electron_response',              'key' : 'teti'}
standardParameters['rho']    = {'knob' : 'geometry_vmec',         'key' : 'rho'}
standardParameters['tprim']  = {'knob' : 'species_parameters_1',    'key' : 'tprim'}
standardParameters['tiprim'] = {'knob' : 'species_parameters_1',    'key' : 'tprim'}
standardParameters['teprim'] = {'knob' : 'species_parameters_2',    'key' : 'tprim'}
standardParameters['fprim']  = {'knob' : 'species_parameters_1',    'key' : 'fprim'}
standardParameters['delt']   = {'knob' : 'time_step',                   'key' : 'delt'}
standardParameters['delta t']= {'knob' : 'time_step',                   'key' : 'delt'}
standardParameters['nmu']    = {'knob' : 'velocity_grids',  'key' : 'nmu'}
standardParameters['nvgrid'] = {'knob' : 'velocity_grids',  'key' : 'nvgrid'}
standardParameters['dvpa']   = {'knob' : 'velocity_grids',  'key' : 'dvpa'}
standardParameters['dmu']    = {'knob' : 'velocity_grids',  'key' : 'dmu'}
standardParameters['nz']     = {'knob' : 'z_grid',        'key' : 'nz'}
standardParameters['nzed']   = {'knob' : 'z_grid',        'key' : 'nzed'}
standardParameters['nzgrid'] = {'knob' : 'z_grid',        'key' : 'nzgrid'}
standardParameters['nx']     = {'knob' : 'kxky_grid_box', 'key' : 'nx'}
standardParameters['ny']     = {'knob' : 'kxky_grid_box', 'key' : 'ny'}
standardParameters['y0']     = {'knob' : 'kxky_grid_box', 'key' : 'y0'}
standardParameters['kx max'] = {'knob' : 'kxky_grid_box', 'key' : 'kx max'}
standardParameters['ky max'] = {'knob' : 'kxky_grid_box', 'key' : 'ky max'}
standardParameters['dkx']    = {'knob' : 'kxky_grid_box', 'key' : 'dkx'}
standardParameters['dky']    = {'knob' : 'kxky_grid_box', 'key' : 'dky'}
standardParameters['Lx']     = {'knob' : 'kxky_grid_box', 'key' : 'Lx'}
standardParameters['Ly']     = {'knob' : 'kxky_grid_box', 'key' : 'Ly'}
standardParameters['rk']     = {'knob' : 'numerical_algorithms',      'key' : 'explicit_algorithm'}
standardParameters['nfield'] = {'knob' : 'geometry_vmec',         'key' : 'nfield_periods'}
standardParameters['pol.turns'] = {'knob' : 'geometry_vmec',      'key' : 'poloidal_turns'}
standardParameters['nperiod']= {'knob' : 'z_grid',        'key' : 'nperiod'}
standardParameters['d_hyper']= {'knob' : 'hyper_dissipation',                   'key' : 'd_hyper'}
standardParameters['D_hyper']= {'knob' : 'hyper_dissipation',             'key' : 'D_hyper'}
standardParameters['alpha0'] = {'knob' : 'geometry_vmec',         'key' : 'alpha0'}
standardParameters['tri']    = {'knob' : 'geometry_miller',    'key' : 'tri'}
standardParameters['kappa']  = {'knob' : 'geometry_miller',    'key' : 'kappa'}
standardParameters['explicit_option']= {'knob' : 'time_advance_knobs', 'key' : 'explicit_option'}
standardParameters['cfl_cushion'] = {'knob' : 'time_step',              'key' : 'cfl_cushion'}
standardParameters['cfl_cushion_upper'] = {'knob' : 'time_step',        'key' : 'cfl_cushion_upper'}
standardParameters['cfl_cushion_lower'] = {'knob' : 'time_step',        'key' : 'cfl_cushion_lower'}
standardParameters['-']      = {'knob' : '-',                       'key' : '-'}
standardParameters['-----']  = {'knob' : '-',                       'key' : '-'} 

