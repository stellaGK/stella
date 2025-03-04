!###############################################################################
!############################## UPDATE INPUT FILE ##############################
!###############################################################################
! Overview of namelists defined in stella:
! 
! Geometry
!    geo_knobs                 contains <geo_option>
!       vmec_parameters            if <geo_option> = 'vmec'
!       millergeo_parameters       if <geo_option> = 'miller'
! 
! Parameters
!    parameters_numerical
!    parameters_physics
! 
! Diagnostics
!    diagnostics
! 
! kxky-grids
!    kt_grids_knobs
!    kt_grids_range_parameters
!    kxky_grids_box
!    kt_grids_range_parameters
!    kt_grids_box_parameters
! 
! Grids
!    zgrid_parameters
!    vpamu_grids_parameters
!    
! Species
!    species_knobs
!    species_parameters
!    
!    init_g_knobs
!    layouts_knobs
!    multibox_parameters
!    sources
!    debug_flags
!    verbose
! 
! Dissipation:
!    dissipation
!    collisions_fp
!    collisions_dougherty
!    hyper
! 
! Neoclassical:
!    sfincs_input
!    euterpe_parameters
!    neoclassical_input
! 
! Deprecated namelists
!    knobs
!    time_advance_knobs
!    physics_flags
!    parameters
! 
! Note that there is a python script 'namelist_to_code.py' to turn the variables
! of a namelist into the pieces of code needed in this script. 
!###############################################################################

module input_file

   implicit none

   public :: update_input_file
   
   private
   
   ! Print verbose to output file
   logical :: print_extra_info_to_terminal
   logical :: print_input_file_warnings
   logical :: debug = .false.

contains

   subroutine update_input_file()

      use file_utils, only: init_file_utils, run_name
      use file_utils, only: open_other_input_file
      use file_utils, only: get_unused_unit

      implicit none
      
      ! Fortran uses the unit number to access the file with later read and write statements  
      ! This routine will write a new file called <run_name>_with_defaults.in
      integer :: input_with_defaults_unit_number
      character(500) :: input_file_with_defaults_name
      
      ! Mandatory argument to use <init_file_utils> even though we do not use it
      logical :: list 
      
      ! For <species_parameters_index> we need <nspec>
      integer :: nspec_local, is
      
      !----------------------------------------------------------------------------
      
      if (debug) write(*, *) 'update_input_file::start'

      ! Allow for --version and --help flags for the executable
      call parse_command_line()

      ! Use stella's routine to read and clean the input file
      call init_file_utils(list)

      ! Get a unit number for the <run_name>_with_defaults.in file we will write
      input_file_with_defaults_name = trim(run_name)//"_with_defaults.in"
      call get_unused_unit(input_with_defaults_unit_number)
      open(unit=input_with_defaults_unit_number, file=input_file_with_defaults_name)
      
      ! First read the <verbose> namelist to see if we defined 
      ! <print_extra_info_to_terminal> or <print_input_file_warnings>
      call read_verbose_namelist(input_with_defaults_unit_number)
      
      ! Read the <run_name>.in file and write it's values + stella default
      ! values to <run_name>_with_defaults.in. Moreover we allow for
      ! backwards compatibility with older stella versions (old namelists). 
      call write_geometry(input_with_defaults_unit_number, .true.) 
      call write_vmec_geometry(input_with_defaults_unit_number, .true.) 
      call write_miller_geometry(input_with_defaults_unit_number, .true.) 
      call write_parameters_physics(input_with_defaults_unit_number, .true.)
      call write_vpamu_grids_parameters(input_with_defaults_unit_number, .true.)
      call write_zgrid_parameters(input_with_defaults_unit_number, .true.)
      call write_species_knobs(input_with_defaults_unit_number, .true., nspec_local)
      do is = 1, nspec_local
         call write_species_parameters(input_with_defaults_unit_number, .true., is)
      end do
      call write_kxky_grids(input_with_defaults_unit_number, .true.)
      call write_kxky_grids_range(input_with_defaults_unit_number, .true.)
      call write_kxky_grids_box(input_with_defaults_unit_number, .true.)
      call write_diagnostics(input_with_defaults_unit_number, .true.)
      call write_init_g_knobs(input_with_defaults_unit_number, .true.)
      call write_dissipation(input_with_defaults_unit_number, .true.)
      call write_collisions_dougherty(input_with_defaults_unit_number, .true.)
      call write_collisions_fp(input_with_defaults_unit_number, .true.)
      call write_hyper(input_with_defaults_unit_number, .true.)
      call write_parameters_numerical(input_with_defaults_unit_number, .true.)
      call write_neoclassical_input(input_with_defaults_unit_number, .true.)
      call write_euterpe_parameters(input_with_defaults_unit_number, .true.)
      call write_sources(input_with_defaults_unit_number, .true.)
      call write_multibox_parameters(input_with_defaults_unit_number, .true.)
      call write_layouts_knobs(input_with_defaults_unit_number, .true.)
      call write_debug_flags(input_with_defaults_unit_number, .true.)

      ! Close <run_name>_with_defaults.in
      close(unit=input_with_defaults_unit_number)
      
      ! For the main stella code, the namelists are read from <run_name.in>,
      ! now switch it to read from <run_name>_with_defaults.in instead. 
      ! This is necessary to use the backwards compatibility.
      call open_other_input_file(input_file_with_defaults_name) 
      
      ! Perform some checks 
      if (debug) write(*, *) 'update_input_file::end'

   end subroutine update_input_file
   
!###############################################################################
!################################# PARAMETERS ##################################
!###############################################################################

   !----------------------------------------------------------------------------
   !                             <verbose> namelist                             
   !----------------------------------------------------------------------------
   subroutine read_verbose_namelist(unit_number)
   
      use file_utils, only: input_unit_exist
   
      implicit none

      integer, intent(in) :: unit_number
      integer :: unit_number_temp
      logical :: nml_exist
      
      namelist /verbose/ print_extra_info_to_terminal, print_input_file_warnings
      
      ! Defaults
      print_extra_info_to_terminal = .true.
      print_input_file_warnings = .true.
      
      ! Read the <verbose> namelist
      unit_number_temp = input_unit_exist("verbose", nml_exist) 
      if (nml_exist) read (unit=unit_number_temp, nml=verbose)
       
      ! Write the namelist to <run_name>_with_defaults.in
      write(unit=unit_number, nml=verbose)
      
   end subroutine read_verbose_namelist
   
   !----------------------------------------------------------------------------
   !                       <parameters_physics> namelist                        
   !----------------------------------------------------------------------------
   subroutine write_parameters_physics(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use parameters_physics, only: set_default_parameters_physics => set_default_parameters
      use parameters_physics, only: prp_shear_enabled
      use parameters_physics, only: hammett_flow_shear, include_pressure_variation, include_geometric_variation
      use parameters_physics, only: full_flux_surface
      use parameters_physics, only: include_apar, include_bpar, radial_variation
      use parameters_physics, only: beta, zeff, rhostar, vnew_ref
      use parameters_physics, only: g_exb, g_exbfac, omprimfac, irhostar
      use namelist_adiabatic_electrons, only: read_adiabatic_electrons_namelist => read_namelist
      use namelist_gyrokinetic_equation, only: read_gyrokinetic_equation_namelist => read_namelist
      use namelist_modify_gyrokinetic_equation, only: read_modify_gyrokinetic_equation_namelist => read_namelist
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Define the variables in the namelist <adiabatic_electrons>
      integer :: adiabatic_option_switch, adiabatic_option_periodic, adiabatic_option_fieldlineavg
      character(30) :: adiabatic_option
      real :: tite, nine
      
      ! Define the variables in the namelist <gyrokinetic_equation>
      logical :: include_parallel_nonlinearity
      logical :: include_parallel_streaming
      logical :: include_mirror, nonlinear
      real :: xdriftknob, ydriftknob, wstarknob
      
      ! Define the variables in the namelist <modify_gyrokinetic_equation>
      logical :: suppress_zonal_interaction
      
      
      ! Current namelist
      namelist /parameters_physics/ prp_shear_enabled, include_apar, include_bpar, radial_variation, &
        hammett_flow_shear, include_pressure_variation, include_geometric_variation, &
        full_flux_surface, g_exb, g_exbfac, omprimfac, irhostar, beta, zeff, rhostar, vnew_ref
      namelist /adiabatic_electrons/ adiabatic_option, tite, nine
      namelist /modify_gyrokinetic_equation/ suppress_zonal_interaction
      namelist /gyrokinetic_equation/ xdriftknob, ydriftknob, wstarknob, nonlinear, &
            include_mirror, include_parallel_streaming, include_parallel_nonlinearity
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_parameters_physics'
   
      ! Read in the default parameters set in stella 
      call set_default_parameters_physics()
      
      ! Read the user-specified input parameters in <run_name>.in
      ! And include backwards compatibility for old namelists
      if (write_input_file) then
         unit_number_temp = input_unit_exist("parameters_physics", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=parameters_physics)
         call backwards_compatibility_parameters()
         call backwards_compatibility_physics_flags()
         call backwards_compatibility_time_advance_knobs()
      end if  
      
      ! TODO - Finish this
      call read_modify_gyrokinetic_equation_namelist(suppress_zonal_interaction)
      call read_gyrokinetic_equation_namelist(xdriftknob, ydriftknob, wstarknob, nonlinear, &
      include_mirror, include_parallel_streaming, include_parallel_nonlinearity)
      call read_adiabatic_electrons_namelist(adiabatic_option_switch, nine, tite, &
         adiabatic_option_fieldlineavg, adiabatic_option_periodic)
      if (adiabatic_option_switch==adiabatic_option_periodic) adiabatic_option = 'no-field-line-average-term'
      if (adiabatic_option_switch==adiabatic_option_fieldlineavg) adiabatic_option = 'field-line-average-term'

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=parameters_physics)
      write(unit=unit_number, nml=adiabatic_electrons)
      write(unit=unit_number, nml=modify_gyrokinetic_equation)
      write(unit=unit_number, nml=gyrokinetic_equation)
      
   contains
   
      !-------------------------------------------------------------------------
      !                          <parameters> namelist                          
      !-------------------------------------------------------------------------
      subroutine backwards_compatibility_parameters()
      
         implicit none 
      
         ! Save parameters which could potentially be defined twice
         real :: beta_new, zeff_new, tite_new, nine_new, rhostar_new, vnew_ref_new
         real :: g_exb_new, g_exbfac_new, omprimfac_new, irhostar_new

         ! Local parameters for this subroutine
         logical :: double_definitions = .false.
         logical :: old_nml_exist
         
         ! Deprecated namelist <parameters>
         namelist /parameters/ beta, zeff, tite, nine, rhostar, vnew_ref, &
            g_exb, g_exbfac, omprimfac, irhostar
            
         !----------------------------------------------------------------------

         ! Check whether the old knob exists
         unit_number_temp = input_unit_exist("parameters", old_nml_exist)
         
         ! If only the old namelist is present in <run_name>.in, and not the new
         ! namelist, simply read it in the old namelist
         if (old_nml_exist .and. (.not. new_nml_exist)) then
            read (unit=unit_number_temp, nml=parameters)
         end if 
         
         ! If both namelists exist, make sure the same variable is not mentioned twice
         if (old_nml_exist .and. new_nml_exist) then
         
            ! Print warning, because the input file should be updated.
            ! We do not want to abort, since we want stella to run the automatic
            ! tests on all versions of stella, including old versions.
            if (print_input_file_warnings) then 
               write(*,*) '   '
                  write(*,*) 'WARNING: The <parameters> namelist has been deprecated. Please move the '
                  write(*,*) 'variables [beta, zeff, tite, nine, rhostar, vnew_ref, g_exb, g_exbfac, '
                  write(*,*) 'omprimfac, irhostar] to the <parameters_physics> namelist.'
               write(*,*) '   '
            end if
            
            ! Save parameters defined in the new namelist <parameters_physics>
            beta_new = beta
            zeff_new = zeff
            tite_new = tite
            nine_new = nine
            rhostar_new = rhostar
            vnew_ref_new = vnew_ref
            g_exb_new = g_exb
            g_exbfac_new = g_exbfac
            omprimfac_new = omprimfac
            irhostar_new = irhostar
            
            ! Read parameters defined in the old namelist <parameters>
            read (unit=unit_number_temp, nml=parameters)
            
            ! Check if any variable from the old namelist <parameters> has changed
            if (beta_new /= beta) double_definitions = .true.
            if (zeff_new /= zeff) double_definitions = .true.
            if (tite_new /= tite) double_definitions = .true.
            if (nine_new /= nine) double_definitions = .true.
            if (rhostar_new /= rhostar) double_definitions = .true.
            if (vnew_ref_new /= vnew_ref) double_definitions = .true.
            if (g_exb_new /= g_exb) double_definitions = .true.
            if (g_exbfac_new /= g_exbfac) double_definitions = .true.
            if (omprimfac_new /= omprimfac) double_definitions = .true.
            if (irhostar_new /= irhostar) double_definitions = .true.
            
            ! If we have double defintions, abort stella
            if (double_definitions) call abort_because_of_double_definitions('parameters_physics', 'parameters')
            
         end if 
         
      end subroutine backwards_compatibility_parameters
      
      !-------------------------------------------------------------------------
      !                        <physics_flags> namelist                         
      !------------------------------------------------------------------------- 
      subroutine backwards_compatibility_physics_flags()
      
         implicit none 
      
         ! Save parameters which could potentially be defined twice
         logical :: include_parallel_streaming_new, include_mirror_new, nonlinear_new
         logical :: include_pressure_variation_new, include_geometric_variation_new
         logical :: include_parallel_nonlinearity_new, suppress_zonal_interaction_new
         logical :: prp_shear_enabled_new, hammett_flow_shear_new
         logical :: full_flux_surface_new, include_apar_new, include_bpar_new, radial_variation_new
         character(30) :: adiabatic_option_new
      
         ! Parameters in old namelists which are not used in <parameters_physics>
         logical :: const_alpha_geo          ! Moved to <debug_flags> 
            
         ! Local parameters for this subroutine
         logical :: double_definitions = .false.
         logical :: old_nml_exist
         
         ! Deprecated namelist <physics_flags>
         namelist /physics_flags/ full_flux_surface, radial_variation, &
            include_parallel_nonlinearity, include_parallel_streaming, &
            include_mirror, include_apar, include_bpar, nonlinear, &
            include_pressure_variation, include_geometric_variation, &
            adiabatic_option, const_alpha_geo, suppress_zonal_interaction

         ! Check whether the old knob exists
         unit_number_temp = input_unit_exist("physics_flags", old_nml_exist)
         
         ! If only the old one exists, and not the new one, simply read it in
         if (old_nml_exist .and. (.not. new_nml_exist)) then
            read (unit=unit_number_temp, nml=physics_flags)
         end if 
         
         ! If both namelists exist, make sure the same variable isn't mentioned twice
         if (old_nml_exist .and. new_nml_exist) then
         
            ! Print warning, because the input file should be updated.
            ! We do not want to abort, since we want stella to run the automatic
            ! tests on all versions of stella, including old versions.
            if (print_input_file_warnings) then 
               write(*,*) '   '
                  write(*,*) 'WARNING: The <physics_flags> namelist has been deprecated. Please move the '
                  write(*,*) 'variables [full_flux_surface, radial_variation, include_parallel_nonlinearity, '
                  write(*,*) 'include_parallel_streaming, include_mirror, include_apar, include_bpar, nonlinear, '
                  write(*,*) 'include_pressure_variation, include_geometric_variation, adiabatic_option, '
                  write(*,*) 'const_alpha_geo, suppress_zonal_interaction] to the <parameters_physics> namelist.'
               write(*,*) '   '
            end if
            
            ! Initialize
            double_definitions = .false.
            
            ! Save parameters defined in the new namelist <parameters_physics>
            include_parallel_streaming_new = include_parallel_streaming
            include_mirror_new = include_mirror
            nonlinear_new = nonlinear
            adiabatic_option_new = adiabatic_option
            prp_shear_enabled_new = prp_shear_enabled
            hammett_flow_shear_new = hammett_flow_shear
            include_pressure_variation_new = include_pressure_variation
            include_geometric_variation_new = include_geometric_variation
            include_parallel_nonlinearity_new = include_parallel_nonlinearity
            suppress_zonal_interaction_new = suppress_zonal_interaction
            full_flux_surface_new = full_flux_surface
            include_apar_new = include_apar
            include_bpar_new = include_bpar
            radial_variation_new = radial_variation
            
            ! Read parameters defined in the old namelist <physics_flags>
            read (unit=unit_number_temp, nml=physics_flags)

            ! Check if any variable from the old namelist <physics_flags> has changed
            if (include_parallel_streaming_new .neqv. include_parallel_streaming) double_definitions = .true.
            if (include_mirror_new .neqv. include_mirror) double_definitions = .true.
            if (nonlinear_new .neqv. nonlinear) double_definitions = .true.
            if (prp_shear_enabled_new .neqv. prp_shear_enabled) double_definitions = .true.
            if (hammett_flow_shear_new .neqv. hammett_flow_shear) double_definitions = .true.
            if (include_pressure_variation_new .neqv. include_pressure_variation) double_definitions = .true.
            if (include_geometric_variation_new .neqv. include_geometric_variation) double_definitions = .true.
            if (include_parallel_nonlinearity_new .neqv. include_parallel_nonlinearity) double_definitions = .true.
            if (suppress_zonal_interaction_new .neqv. suppress_zonal_interaction) double_definitions = .true.
            if (full_flux_surface_new .neqv. full_flux_surface) double_definitions = .true.
            if (include_apar_new .neqv. include_apar) double_definitions = .true.
            if (include_bpar_new .neqv. include_bpar) double_definitions = .true.
            if (radial_variation_new .neqv. radial_variation) double_definitions = .true.
            if (adiabatic_option_new /= adiabatic_option) double_definitions = .true.
            
            ! If we have double defintions, abort stella
            if (double_definitions) call abort_because_of_double_definitions('parameters_physics', 'physics_flags')  
            
         end if 
         
      end subroutine backwards_compatibility_physics_flags
      
      !-------------------------------------------------------------------------
      !                     <time_advance_knobs> namelist                       
      !------------------------------------------------------------------------- 
      subroutine backwards_compatibility_time_advance_knobs()
      
         implicit none 
      
         ! Save parameters which could potentially be defined twice
         real :: xdriftknob_new, ydriftknob_new, wstarknob_new 
      
         ! Parameters in old namelists which are not used in <parameters_physics>
         character(10) :: explicit_option    ! Moved to <parameters_numerical>
         logical :: flip_flop                ! Moved to <parameters_numerical>
            
         ! Local parameters for this subroutine
         logical :: old_nml_exist
         logical :: double_definitions
         
         ! Old namelists 
         namelist /time_advance_knobs/ xdriftknob, ydriftknob, wstarknob, explicit_option, flip_flop 

         ! Check whether the old knob exists
         unit_number_temp = input_unit_exist("time_advance_knobs", old_nml_exist)
         
         ! If only the old one exists, and not the new one, simply read it in
         if (old_nml_exist .and. (.not. new_nml_exist)) then
            read (unit=unit_number_temp, nml=time_advance_knobs)
         end if 
         
         ! If both namelists exist, make sure the same variable isn't mentioned twice
         if (old_nml_exist .and. new_nml_exist) then
         
            ! Print warning, because the input file should be updated.
            ! We do not want to abort, since we want stella to run the automatic
            ! tests on all versions of stella, including old versions.
            if (print_input_file_warnings) then 
               write(*,*) '   '
                  write(*,*) 'WARNING: The <time_advance_knobs> namelist has been deprecated. Please move the '
                  write(*,*) 'variables  [xdriftknob, ydriftknob, wstarknob, explicit_option, flip_flop ] to '
                  write(*,*) 'the <parameters_physics> namelist.'
               write(*,*) '   '
            end if
            
            ! Initialize
            double_definitions = .false.
            
            ! Save parameters defined in the new namelist <parameters_physics>
            xdriftknob_new = xdriftknob
            ydriftknob_new = ydriftknob
            wstarknob_new = wstarknob
            
            ! Read parameters defined in the old namelist <parameters>
            read (unit=unit_number_temp, nml=time_advance_knobs)

            ! Check if any variable from the old namelist <time_advance_knobs> has changed
            if (xdriftknob_new /= xdriftknob) double_definitions = .true.
            if (ydriftknob_new /= ydriftknob) double_definitions = .true.
            if (wstarknob_new /= wstarknob) double_definitions = .true.
            
            ! If we have double defintions, abort stella
            if (double_definitions) call abort_because_of_double_definitions('parameters_physics', 'time_advance_knobs') 
            
         end if 
         
      end subroutine backwards_compatibility_time_advance_knobs

   end subroutine write_parameters_physics
   
   !----------------------------------------------------------------------------
   !                      <parameters_numerical> namelist                       
   !----------------------------------------------------------------------------

   subroutine write_parameters_numerical(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use parameters_numerical, only: set_default_parameters_numerical => set_default_parameters
      use parameters_numerical, only: stream_implicit, stream_iterative_implicit, stream_matrix_inversion
      use parameters_numerical, only: driftkinetic_implicit, mirror_implicit, mirror_semi_lagrange, mirror_linear_interp
      use parameters_numerical, only: drifts_implicit, fully_implicit, fully_explicit, fields_kxkyz, mat_gen, mat_read
      use parameters_numerical, only: maxwellian_inside_zed_derivative, use_deltaphi_for_response_matrix
      use parameters_numerical, only: maxwellian_normalization, zed_upwind, vpa_upwind, time_upwind
      use parameters_numerical, only: fphi, nstep, delt, tend, delt_option, lu_option, avail_cpu_time
      use parameters_numerical, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower, delt_max, delt_min 
      use parameters_numerical, only: ky_solve_radial, ky_solve_real, nitt, print_extra_info_to_terminal
      use parameters_numerical, only: explicit_option, flip_flop, rng_seed, autostop
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Current namelist
      namelist /parameters_numerical/ stream_implicit, stream_iterative_implicit, stream_matrix_inversion, &
            driftkinetic_implicit, mirror_implicit, mirror_semi_lagrange, mirror_linear_interp, &
            drifts_implicit, fully_implicit, fully_explicit, & 
            maxwellian_inside_zed_derivative, use_deltaphi_for_response_matrix, &
            maxwellian_normalization, zed_upwind, vpa_upwind, time_upwind, &
            fphi, nstep, delt, tend, delt_option, lu_option, avail_cpu_time, & 
            cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower, delt_max, delt_min, &
            fields_kxkyz, mat_gen, mat_read, &
            ky_solve_radial, ky_solve_real, nitt, print_extra_info_to_terminal, &
            explicit_option, flip_flop, rng_seed, autostop
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_parameters_numerical'
   
      ! Read in the default parameters set in stella 
      call set_default_parameters_numerical()
      
      ! Read the user-specified input parameters in <run_name>.in
      ! And include backwards compatibility for old namelists
      if (write_input_file) then
         unit_number_temp = input_unit_exist("parameters_numerical", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=parameters_numerical) 
         call backwards_compatibility_time_advance_knobs()
         call backwards_compatibility_knobs()
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=parameters_numerical)
      
   contains
      
      !-------------------------------------------------------------------------
      !                     <time_advance_knobs> namelist                       
      !------------------------------------------------------------------------- 
      subroutine backwards_compatibility_time_advance_knobs()
      
         implicit none 
      
         ! Save parameters which could potentially be defined twice
         logical :: flip_flop_new
         character(30) :: explicit_option_new
              
         ! Parameters in old namelists which are not used in <parameters_numerical>
         real :: xdriftknob, ydriftknob, wstarknob ! Moved to <parameters_physics>
            
         ! Local parameters for this subroutine
         logical :: old_nml_exist
         logical :: double_definitions
         
         ! Old namelists 
         namelist /time_advance_knobs/ xdriftknob, ydriftknob, wstarknob, explicit_option, flip_flop 

         ! Check whether the old knob exists
         unit_number_temp = input_unit_exist("time_advance_knobs", old_nml_exist)
         
         ! If only the old one exists, and not the new one, simply read it in
         if (old_nml_exist .and. (.not. new_nml_exist)) then
            read (unit=unit_number_temp, nml=time_advance_knobs)
         end if 
         
         ! If both namelists exist, make sure the same variable isn't mentioned twice
         if (old_nml_exist .and. new_nml_exist) then
         
            ! Print warning, because the input file should be updated.
            ! We do not want to abort, since we want stella to run the automatic
            ! tests on all versions of stella, including old versions.
            if (print_input_file_warnings) then 
               write(*,*) '   '
                  write(*,*) 'WARNING: The <time_advance_knobs> namelist has been deprecated. Please move the '
                  write(*,*) 'variables  [xdriftknob, ydriftknob, wstarknob, explicit_option, flip_flop ] to '
                  write(*,*) 'the <parameters_physics> namelist.'
               write(*,*) '   '
            end if
            
            ! Initialize
            double_definitions = .false.
            
            ! Save parameters defined in the new namelist <parameters_physics>
            explicit_option_new = explicit_option
            flip_flop_new = flip_flop
            
            ! Read parameters defined in the old namelist <parameters>
            read (unit=unit_number_temp, nml=time_advance_knobs)

            ! Check if any variable from the old namelist <time_advance_knobs> has changed
            if (explicit_option_new /= explicit_option) double_definitions = .true.
            if (flip_flop_new .neqv. flip_flop) double_definitions = .true.
            
            ! If we have double defintions, abort stella
            if (double_definitions) call abort_because_of_double_definitions('parameters_numerical', 'time_advance_knobs') 
            
         end if 
         
      end subroutine backwards_compatibility_time_advance_knobs
      
      !-------------------------------------------------------------------------
      !                            <knobs> namelist                             
      !------------------------------------------------------------------------- 
      subroutine backwards_compatibility_knobs()
      
         implicit none 
      
         ! Save parameters which could potentially be defined twice 
         real :: fphi_new, delt_new, tend_new, avail_cpu_time_new
         real :: delt_max_new, delt_min_new, cfl_cushion_upper_new, cfl_cushion_middle_new, cfl_cushion_lower_new
         real :: zed_upwind_new, vpa_upwind_new, time_upwind_new
         integer :: nstep_new, rng_seed_new, ky_solve_radial_new, nitt_new
         logical :: autostop_new, stream_implicit_new, mirror_implicit_new, drifts_implicit_new
         logical :: use_deltaphi_for_response_matrix_new, maxwellian_normalization_new, stream_matrix_inversion_new
         logical :: maxwellian_inside_zed_derivative_new, mirror_semi_lagrange_new, mirror_linear_interp_new
         logical :: fields_kxkyz_new, mat_gen_new, mat_read_new, ky_solve_real_new, print_extra_info_to_terminal_new
         character(20) :: delt_option_new, lu_option_new
         
         ! Old parameters which have been deprecated for a while
         real :: fapar, fbpar
            
         ! Local parameters for this subroutine
         logical :: old_nml_exist
         logical :: double_definitions
         
         ! Old namelists 
         namelist /knobs/ fphi, fapar, fbpar, delt, nstep, tend, delt_option, lu_option, autostop, &
              avail_cpu_time, delt_max, delt_min, cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower, &
              stream_implicit, mirror_implicit, drifts_implicit, use_deltaphi_for_response_matrix, &
              maxwellian_normalization, stream_matrix_inversion, maxwellian_inside_zed_derivative, &
              mirror_semi_lagrange, mirror_linear_interp, zed_upwind, vpa_upwind, time_upwind, &
              fields_kxkyz, mat_gen, mat_read, rng_seed, ky_solve_radial, ky_solve_real, nitt, &
              print_extra_info_to_terminal

         ! Check whether the old knob exists
         unit_number_temp = input_unit_exist("knobs", old_nml_exist)
         
         ! If only the old one exists, and not the new one, simply read it in
         if (old_nml_exist .and. (.not. new_nml_exist)) then
            read (unit=unit_number_temp, nml=knobs)
         end if 
         
         ! If both namelists exist, make sure the same variable isn't mentioned twice
         if (old_nml_exist .and. new_nml_exist) then
         
            ! Print warning, because the input file should be updated.
            ! We do not want to abort, since we want stella to run the automatic
            ! tests on all versions of stella, including old versions.
            if (print_input_file_warnings) then 
               write(*,*) '   '
                  write(*,*) 'WARNING: The <knobs> namelist has been deprecated. Please move the '
                  write(*,*) 'variables  [fphi, delt, nstep, tend, delt_option, lu_option, autostop, '
                  write(*,*) 'avail_cpu_time, delt_max, delt_min, cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower, '
                  write(*,*) 'stream_implicit, mirror_implicit, drifts_implicit, use_deltaphi_for_response_matrix, '
                  write(*,*) 'maxwellian_normalization, stream_matrix_inversion, maxwellian_inside_zed_derivative, '
                  write(*,*) 'mirror_semi_lagrange, mirror_linear_interp, zed_upwind, vpa_upwind, time_upwind, '
                  write(*,*) 'fields_kxkyz, mat_gen, mat_read, rng_seed, ky_solve_radial, ky_solve_real, nitt] to '
                  write(*,*) 'the <parameters_physics> namelist.'
               write(*,*) '   '
            end if
            
            ! Initialize
            double_definitions = .false.
            
            ! Save parameters defined in the new namelist <parameters_physics>
            fphi_new = fphi
            delt_new = delt
            nstep_new = nstep
            tend_new = tend
            delt_option_new = delt_option
            lu_option_new = lu_option
            autostop_new = autostop
            avail_cpu_time_new = avail_cpu_time
            delt_max_new = delt_max
            delt_min_new = delt_min
            cfl_cushion_upper_new = cfl_cushion_upper
            cfl_cushion_middle_new = cfl_cushion_middle
            cfl_cushion_lower_new = cfl_cushion_lower
            stream_implicit_new = stream_implicit
            mirror_implicit_new = mirror_implicit
            drifts_implicit_new = drifts_implicit
            use_deltaphi_for_response_matrix_new = use_deltaphi_for_response_matrix
            maxwellian_normalization_new = maxwellian_normalization
            stream_matrix_inversion_new = stream_matrix_inversion
            maxwellian_inside_zed_derivative_new = maxwellian_inside_zed_derivative
            mirror_semi_lagrange_new = mirror_semi_lagrange
            mirror_linear_interp_new = mirror_linear_interp
            zed_upwind_new = zed_upwind
            vpa_upwind_new = vpa_upwind
            time_upwind_new = time_upwind
            fields_kxkyz_new = fields_kxkyz
            mat_gen_new = mat_gen
            mat_read_new = mat_read
            rng_seed_new = rng_seed
            ky_solve_radial_new = ky_solve_radial
            ky_solve_real_new = ky_solve_real
            nitt_new = nitt
            print_extra_info_to_terminal_new = print_extra_info_to_terminal
            
            ! Read parameters defined in the old namelist <parameters>
            read (unit=unit_number_temp, nml=knobs)

            ! Check if any variable from the old namelist <time_advanknobsce_knobs> has changed
            if (fphi_new /= fphi) double_definitions = .true.
            if (delt_new /= delt) double_definitions = .true.
            if (nstep_new /= nstep) double_definitions = .true.
            if (tend_new /= tend) double_definitions = .true.
            if (delt_option_new /= delt_option) double_definitions = .true.
            if (lu_option_new /= lu_option) double_definitions = .true.
            if (autostop_new .neqv. autostop) double_definitions = .true.
            if (avail_cpu_time_new /= avail_cpu_time) double_definitions = .true.
            if (delt_max_new /= delt_max) double_definitions = .true.
            if (delt_min_new /= delt_min) double_definitions = .true.
            if (cfl_cushion_upper_new /= cfl_cushion_upper) double_definitions = .true.
            if (cfl_cushion_middle_new /= cfl_cushion_middle) double_definitions = .true.
            if (cfl_cushion_lower_new /= cfl_cushion_lower) double_definitions = .true.
            if (stream_implicit_new .neqv. stream_implicit) double_definitions = .true.
            if (mirror_implicit_new .neqv. mirror_implicit) double_definitions = .true.
            if (drifts_implicit_new .neqv.drifts_implicit) double_definitions = .true.
            if (use_deltaphi_for_response_matrix_new .neqv. use_deltaphi_for_response_matrix) double_definitions = .true.
            if (maxwellian_normalization_new .neqv. maxwellian_normalization) double_definitions = .true.
            if (stream_matrix_inversion_new .neqv. stream_matrix_inversion) double_definitions = .true.
            if (maxwellian_inside_zed_derivative_new .neqv. maxwellian_inside_zed_derivative) double_definitions = .true.
            if (mirror_semi_lagrange_new .neqv. mirror_semi_lagrange) double_definitions = .true.
            if (mirror_linear_interp_new .neqv. mirror_linear_interp) double_definitions = .true.
            if (zed_upwind_new /= zed_upwind) double_definitions = .true.
            if (vpa_upwind_new /= vpa_upwind) double_definitions = .true.
            if (time_upwind_new /= time_upwind) double_definitions = .true.
            if (fields_kxkyz_new .neqv. fields_kxkyz) double_definitions = .true.
            if (mat_gen_new .neqv. mat_gen) double_definitions = .true.
            if (mat_read_new .neqv. mat_read) double_definitions = .true.
            if (rng_seed_new /= rng_seed) double_definitions = .true.
            if (ky_solve_radial_new /= ky_solve_radial) double_definitions = .true.
            if (ky_solve_real_new .neqv. ky_solve_real) double_definitions = .true.
            if (nitt_new /= nitt) double_definitions = .true.
            if (print_extra_info_to_terminal_new .neqv. print_extra_info_to_terminal) double_definitions = .true.
            
            ! If we have double defintions, abort stella
            if (double_definitions) call abort_because_of_double_definitions('parameters_numerical', 'knobs') 
            
         end if 
         
      end subroutine backwards_compatibility_knobs
   
   end subroutine write_parameters_numerical
   
   
!###############################################################################
!################################## GEOMETRY ###################################
!###############################################################################
   
   !----------------------------------------------------------------------------
   !                            <geo_knobs> namelist                            
   !----------------------------------------------------------------------------
   subroutine write_geometry(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use geometry, only: set_defaults_for_geo_knobs => init_defaults_geometry 
      use geometry, only: geo_option, geo_file, overwrite_bmag, overwrite_b_dot_grad_zeta
      use geometry, only: overwrite_gds2, overwrite_gds21, overwrite_gds22, overwrite_gds23, overwrite_gds24
      use geometry, only: overwrite_gbdrift, overwrite_cvdrift, overwrite_gbdrift0, q_as_x, set_bmag_const
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /geo_knobs/ geo_option, geo_file, overwrite_bmag, overwrite_b_dot_grad_zeta, &
         overwrite_gds2, overwrite_gds21, overwrite_gds22, overwrite_gds23, overwrite_gds24, &
         overwrite_gbdrift, overwrite_cvdrift, overwrite_gbdrift0, q_as_x, set_bmag_const
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_geometry'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_geo_knobs()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("geo_knobs", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=geo_knobs) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=geo_knobs)
      
   end subroutine write_geometry
   
   !----------------------------------------------------------------------------
   !                            <millergeo_parameters> namelist                            
   !----------------------------------------------------------------------------
   subroutine write_miller_geometry(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use geometry_miller, only: set_defaults_for_millergeo_parameters => init_local_defaults
      use geometry_miller, only: rhoc, rmaj, shift, qinp, shat, kappa, kapprim
      use geometry_miller, only: tri, triprim, rgeo, betaprim, betadbprim, d2qdr2, d2psidr2
      use geometry_miller, only: nzed_local, read_profile_variation, write_profile_variation
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /millergeo_parameters/ rhoc, rmaj, shift, qinp, shat, &
         kappa, kapprim, tri, triprim, rgeo, betaprim, betadbprim, d2qdr2, d2psidr2, &
         nzed_local, read_profile_variation, write_profile_variation
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_miller_geometry'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_millergeo_parameters()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("millergeo_parameters", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=millergeo_parameters) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=millergeo_parameters)
      
   end subroutine write_miller_geometry
   
   !----------------------------------------------------------------------------
   !                         <vmec_parameters> namelist                         
   !----------------------------------------------------------------------------
   subroutine write_vmec_geometry(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use vmec_geometry, only: set_defaults_for_vmec_parameters => init_vmec_defaults
      use vmec_geometry, only: alpha0, zeta_center, rectangular_cross_section, nfield_periods
      use vmec_geometry, only: torflux, zgrid_refinement_factor, surface_option, radial_coordinate
      use vmec_geometry, only: verbose, vmec_filename, n_tolerated_test_arrays_inconsistencies
      use vmec_geometry, only: zgrid_scalefac ! Backwards compatibility for old stella code
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /vmec_parameters/ alpha0, zeta_center, rectangular_cross_section, nfield_periods, &
         torflux, zgrid_refinement_factor, surface_option, radial_coordinate, &
         verbose, vmec_filename, n_tolerated_test_arrays_inconsistencies, &
         ! Backwards compatibility for old stella code, when we change the name of
         ! <vmec_parameters> we can get rid of this parameter, now old input files
         ! will break if we remove this
         zgrid_scalefac
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_vmec_geometry'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_vmec_parameters()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("vmec_parameters", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=vmec_parameters) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=vmec_parameters)
      
   end subroutine write_vmec_geometry
   
!###############################################################################
!################################ DIAGNOSTICS ##################################
!###############################################################################

   !----------------------------------------------------------------------------
   !                     <diagnostics> namelist                    
   !---------------------------------------------------------------------------- 
  
   subroutine write_diagnostics(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use parameters_diagnostics, only: set_defaults_for_diagnostics => set_default_parameters
      use parameters_diagnostics, only: nwrite, navg, nsave, autostop, save_for_restart, flux_norm, nc_mult
      use parameters_diagnostics, only: write_phi2_vs_time, write_apar2_vs_time, write_bpar2_vs_time, write_fluxes_vs_time
      use parameters_diagnostics, only: write_phi_vs_kxkyz, write_g2_vs_vpamus, write_g2_vs_zvpas, write_g2_vs_zmus
      use parameters_diagnostics, only: write_g2_vs_kxkyzs, write_g2_vs_zvpamus, write_distribution_g
      use parameters_diagnostics, only: write_distribution_h, write_distribution_f, write_omega_vs_kxky
      use parameters_diagnostics, only: write_omega_avg_vs_kxky, write_phi2_vs_kxky, write_moments, write_radial_fluxes
      use parameters_diagnostics, only: write_radial_moments, write_fluxes_kxkyz, write_fluxes_kxky, write_all, flux_norm, nc_mult
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /diagnostics/ nwrite, navg, nsave, autostop, save_for_restart, flux_norm, nc_mult, &
            write_phi2_vs_time, write_apar2_vs_time, write_bpar2_vs_time, write_fluxes_vs_time, &
            write_phi_vs_kxkyz, write_g2_vs_vpamus, write_g2_vs_zvpas, write_g2_vs_zmus, &
            write_g2_vs_kxkyzs, write_g2_vs_zvpamus, write_distribution_g, write_distribution_h, write_distribution_f, &
            write_omega_vs_kxky, write_omega_avg_vs_kxky, write_phi2_vs_kxky, write_moments, write_radial_fluxes, &
            write_radial_moments, write_fluxes_kxkyz, write_fluxes_kxky, write_all, flux_norm, nc_mult
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_diagnostics'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_diagnostics()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("diagnostics", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=diagnostics) 
         call write_all_diagnostics()
         call backwards_compatibility_stella_diagnostics_knobs() 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=diagnostics)
      
   contains
   
      ! To allow both <diagnostics> and <stella_diagnostics_knobs> we need to 
      ! to turn on all diagnostics if <write_all> = .true.
      subroutine write_all_diagnostics()
      
         implicit none
         
         if (write_all) then 
            write_phi2_vs_time = .true.
            write_apar2_vs_time = .true.
            write_bpar2_vs_time = .true.
            write_fluxes_vs_time = .true. 
            write_omega_vs_kxky = .true.
            write_omega_avg_vs_kxky = .true.
            write_phi_vs_kxkyz = .true.
            write_phi2_vs_kxky = .true.
            write_moments = .true.
            write_fluxes_kxkyz = .true.
            write_fluxes_kxky = .true.
            write_g2_vs_vpamus = .true.
            write_g2_vs_zvpas = .true.
            write_g2_vs_zmus = .true.
            write_g2_vs_kxkyzs = .true.
            write_g2_vs_zvpamus = .true.
            write_distribution_g = .true.
            write_distribution_h = .true.
            write_distribution_f = .true.
         end if
         
      end subroutine write_all_diagnostics

      !-------------------------------------------------------------------------
      !                     <stella_diagnostics_knobs> namelist                       
      !------------------------------------------------------------------------- 
      subroutine backwards_compatibility_stella_diagnostics_knobs()
      
         implicit none 
      
         ! Save parameters which could potentially be defined twice
         logical :: save_for_restart_new, write_phi2_vs_time_new, write_apar2_vs_time_new
         logical :: write_bpar2_vs_time_new, write_omega_vs_kxky_new, write_omega_avg_vs_kxky_new
         logical :: write_distribution_g_new, write_g2_vs_vpamus_new
         logical :: write_g2_vs_zvpas_new, write_fluxes_kxky_new, flux_norm_new
         integer :: nwrite_new, navg_new, nsave_new, nc_mult_new
         
         ! Old parameters which are not used anymore
         logical :: write_phi_vs_time, write_gvmus, write_gzvs
         logical :: write_omega, write_kspectra, write_apar_vs_time, write_bpar_vs_time
      
         ! Local parameters for this subroutine
         logical :: old_nml_exist
         logical :: double_definitions
          
         ! Old namelists 
         namelist /stella_diagnostics_knobs/ nwrite, navg, nsave, &
            save_for_restart, write_phi_vs_time, write_gvmus, write_gzvs, &
            write_omega, write_kspectra, write_moments, write_radial_fluxes, &
            write_radial_moments, write_fluxes_kxkyz, flux_norm, nc_mult, &
            write_apar_vs_time, write_bpar_vs_time
            
         ! Set default values for the parameters which are not used anymore
         write_phi_vs_time = .false.
         write_gvmus = .false.
         write_gzvs = .false.
         write_omega = .false.
         write_kspectra = .false.
         
         ! Set the following two values to true since they do not exist on all old versions
         write_apar_vs_time = .true.
         write_bpar_vs_time = .true.
         
         ! Check whether the old knob exists
         unit_number_temp = input_unit_exist("stella_diagnostics_knobs", old_nml_exist)
         
         ! If only the old one exists, and not the new one, simply read it in
         if (old_nml_exist .and. (.not. new_nml_exist)) then
            read (unit=unit_number_temp, nml=stella_diagnostics_knobs)
         end if 
         
         ! If both namelists exist, make sure the same variable isn't mentioned twice
         if (old_nml_exist .and. new_nml_exist) then
         
            ! Print warning, because the input file should be updated.
            ! We do not want to abort, since we want stella to run the automatic
            ! tests on all versions of stella, including old versions.
            if (print_input_file_warnings) then 
               write(*,*) '   '
               write(*,*) 'WARNING: The <stella_diagnostics_knobs> namelist has been deprecated. Please move the '
               write(*,*) 'variables  [nwrite, navg, nsave, save_for_restart, write_phi_vs_time, write_gvmus, '
               write(*,*) 'write_gzvs, write_omega, write_kspectra, write_moments, write_radial_fluxes, '
               write(*,*) 'write_radial_moments, write_fluxes_kxkyz, flux_norm, nc_mult, write_apar_vs_time, '
               write(*,*) 'write_bpar_vs_time] to the <diagnostics> namelist.'
               write(*,*) '   '
            end if
            
            ! Initialize
            double_definitions = .false.
            
            ! Save parameters defined in the new namelist <diagnostics>
            nwrite_new = nwrite
            navg_new = navg
            nsave_new = nsave
            nc_mult_new = nc_mult
            save_for_restart_new = save_for_restart
            flux_norm_new = flux_norm
            
            ! New variables names
            write_phi2_vs_time_new = write_phi2_vs_time
            write_apar2_vs_time_new = write_apar2_vs_time
            write_bpar2_vs_time_new = write_bpar2_vs_time
            write_omega_vs_kxky_new = write_omega_vs_kxky
            write_omega_avg_vs_kxky_new = write_omega_avg_vs_kxky
            write_distribution_g_new = write_distribution_g
            write_g2_vs_vpamus_new = write_g2_vs_vpamus
            write_g2_vs_zvpas_new = write_g2_vs_zvpas
            write_fluxes_kxky_new = write_fluxes_kxky
            
            ! Read parameters defined in the old namelist <stella_diagnostics_knobs>
            read (unit=unit_number_temp, nml=stella_diagnostics_knobs)
            
            ! These flags have been removed/renamed
            if (write_phi_vs_time) write_phi2_vs_time = .true.
            if (write_apar_vs_time) write_apar2_vs_time = .true.
            if (write_bpar_vs_time) write_bpar2_vs_time = .true.
            if (write_omega) write_omega_vs_kxky = .true.
            if (write_omega) write_omega_avg_vs_kxky = .true.
            if (write_gvmus) write_distribution_g = .true.
            if (write_gvmus) write_g2_vs_vpamus = .true.
            if (write_gzvs) write_distribution_g = .true.
            if (write_gzvs) write_g2_vs_zvpas = .true.
            if (write_kspectra) write_fluxes_kxky = .true. 

            ! Check if any variable from the old namelist <stella_diagnostics_knobs> has changed
            if (nwrite_new /= nwrite) double_definitions = .true.
            if (navg_new /= navg) double_definitions = .true.
            if (nsave_new /= nsave) double_definitions = .true.
            if (nc_mult_new /= nc_mult) double_definitions = .true.
            if (flux_norm_new .neqv. flux_norm) double_definitions = .true.
            if (save_for_restart_new .neqv. save_for_restart) double_definitions = .true. 
            if (write_phi2_vs_time_new .neqv. write_phi2_vs_time) double_definitions = .true.
            if (write_apar2_vs_time_new .neqv. write_apar2_vs_time) double_definitions = .true.
            if (write_bpar2_vs_time_new .neqv. write_bpar2_vs_time) double_definitions = .true.
            if (write_omega_vs_kxky_new .neqv. write_omega_vs_kxky) double_definitions = .true.
            if (write_omega_avg_vs_kxky_new .neqv. write_omega_avg_vs_kxky) double_definitions = .true.
            if (write_distribution_g_new .neqv. write_distribution_g) double_definitions = .true.
            if (write_g2_vs_vpamus_new .neqv. write_g2_vs_vpamus) double_definitions = .true.
            if (write_g2_vs_zvpas_new .neqv. write_g2_vs_zvpas) double_definitions = .true.
            if (write_fluxes_kxky_new .neqv. write_fluxes_kxky) double_definitions = .true.

            ! If we have double defintions, abort stella
            if (double_definitions) call abort_because_of_double_definitions('diagnostics', 'stella_diagnostics_knobs') 
            
         end if 
         
      end subroutine backwards_compatibility_stella_diagnostics_knobs
      
   end subroutine write_diagnostics

!###############################################################################
!############################## KX KY GRIDS ####################################
!###############################################################################

   !----------------------------------------------------------------------------
   !                          <kt_grids_knobs> namelist                         
   !----------------------------------------------------------------------------
  
   subroutine write_kxky_grids(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use parameters_kxky_grids, only: set_defaults_for_kxky_grids => set_default_parameters
      use parameters_kxky_grids, only: grid_option
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
       namelist /kxky_grids/ grid_option
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_kxky_grids'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_kxky_grids()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("kxky_grids", new_nml_exist)
         if (new_nml_exist) read (unit=unit_number_temp, nml=kxky_grids)
         call backwards_compatibility_kt_grids_knobs()
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=kxky_grids)
       
   contains
   
      !-------------------------------------------------------------------------
      !                  <kt_grids_knobs> namelist                   
      !-------------------------------------------------------------------------
      subroutine backwards_compatibility_kt_grids_knobs()
      
         implicit none 
      
         ! Save parameters which could potentially be defined twice
         character(20) :: grid_option_new

         ! Local parameters for this subroutine
         logical :: double_definitions = .false.
         logical :: old_nml_exist
         
         ! Deprecated namelist <kt_grids_knobs>
         namelist /kt_grids_knobs/ grid_option
            
         !----------------------------------------------------------------------

         ! Check whether the old knob exists
         unit_number_temp = input_unit_exist("kt_grids_knobs", old_nml_exist)
         
         ! If only the old namelist is present in <run_name>.in, and not the new
         ! namelist, simply read it in the old namelist
         if (old_nml_exist .and. (.not. new_nml_exist)) then
            read (unit=unit_number_temp, nml=kt_grids_knobs)
         end if 
         
         ! If both namelists exist, make sure the same variable is not mentioned twice
         if (old_nml_exist .and. new_nml_exist) then
         
            ! Print warning, because the input file should be updated.
            ! We do not want to abort, since we want stella to run the automatic
            ! tests on all versions of stella, including old versions.
            if (print_input_file_warnings) then 
               write(*,*) '   '
               write(*,*) 'WARNING: The <kt_grids_knobs> namelist has been deprecated. Please move the '
               write(*,*) 'variables [nx, ny, jtwist, jtwistfac, x0, y0,  centered_in_rho, periodic_variation, '
               write(*,*) 'randomize_phase_shift, phase_shift_angle] to the <kxky_grids> namelist.'
               write(*,*) '   '
            end if
            
            ! Save parameters defined in the new namelist <kxky_grids>
            grid_option_new = grid_option
            
            ! Read parameters defined in the old namelist <kt_grids_knobs>
            read (unit=unit_number_temp, nml=kt_grids_knobs)
            
            ! Check if any variable from the old namelist <kt_grids_knobs> has changed
            if (grid_option_new /= grid_option) double_definitions = .true. 

            ! If we have double defintions, abort stella
            if (double_definitions) call abort_because_of_double_definitions('kxky_grids', 'kt_grids_knobs')
            
         end if 
         
      end subroutine backwards_compatibility_kt_grids_knobs
      
   end subroutine write_kxky_grids

   !----------------------------------------------------------------------------
   !                     <kxky_grids_box> namelist                   
   !----------------------------------------------------------------------------
   subroutine write_kxky_grids_box(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use parameters_kxky_grids_box, only: set_defaults_for_kxky_grids_box => read_default_box
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! For <parameters_kxky_grids_box> the variables are not made public
      ! because only the variables in <parameters_kxky_grids> should be used
      ! in stella, so we need to define the namelist variables locally 
      integer :: nx, ny, nalpha, jtwist 
      real :: jtwistfac, phase_shift_angle, x0, y0
      logical :: centered_in_rho, periodic_variation
      logical :: randomize_phase_shift
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /kxky_grids_box/ nx, ny, jtwist, jtwistfac, x0, y0, &
         centered_in_rho, periodic_variation, randomize_phase_shift, phase_shift_angle
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_kxky_grids_box'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_kxky_grids_box(nx, ny, nalpha, &
         jtwist, x0, y0, jtwistfac, phase_shift_angle, &
         centered_in_rho, periodic_variation, randomize_phase_shift)
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("kxky_grids_box", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=kxky_grids_box)
         call backwards_compatibility_kt_grids_box_parameters()
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=kxky_grids_box)

   contains
   
      !-------------------------------------------------------------------------
      !                  <kt_grids_box_parameters> namelist                   
      !-------------------------------------------------------------------------
      subroutine backwards_compatibility_kt_grids_box_parameters()
      
         implicit none 
      
         ! Save parameters which could potentially be defined twice
         integer :: nx_new, ny_new, jtwist_new
         real :: x0_new, y0_new, jtwistfac_new, phase_shift_angle_new 
         logical :: centered_in_rho_new, periodic_variation_new, randomize_phase_shift_new

         ! Local parameters for this subroutine
         logical :: double_definitions = .false.
         logical :: old_nml_exist
         
         ! Deprecated namelist <kt_grids_box_parameters>
         namelist /kt_grids_box_parameters/ nx, ny, jtwist, jtwistfac, x0, y0, &
            centered_in_rho, periodic_variation, randomize_phase_shift, phase_shift_angle
            
         !----------------------------------------------------------------------

         ! Check whether the old knob exists
         unit_number_temp = input_unit_exist("kt_grids_box_parameters", old_nml_exist)
         
         ! If only the old namelist is present in <run_name>.in, and not the new
         ! namelist, simply read it in the old namelist
         if (old_nml_exist .and. (.not. new_nml_exist)) then
            read (unit=unit_number_temp, nml=kt_grids_box_parameters)
         end if 
         
         ! If both namelists exist, make sure the same variable is not mentioned twice
         if (old_nml_exist .and. new_nml_exist) then
         
            ! Print warning, because the input file should be updated.
            ! We do not want to abort, since we want stella to run the automatic
            ! tests on all versions of stella, including old versions.
            if (print_input_file_warnings) then 
               write(*,*) '   '
               write(*,*) 'WARNING: The <kt_grids_box_parameters> namelist has been deprecated. Please move the '
               write(*,*) 'variables [nx, ny, jtwist, jtwistfac, x0, y0,  centered_in_rho, periodic_variation, '
               write(*,*) 'randomize_phase_shift, phase_shift_angle] to the <kxky_grids_box> namelist.'
               write(*,*) '   '
            end if
            
            ! Save parameters defined in the new namelist <kxky_grids_box>
            nx_new = nx
            ny_new = ny
            jtwist_new = jtwist
            jtwistfac_new = jtwistfac
            x0_new = x0
            y0_new = y0
            centered_in_rho_new = centered_in_rho
            periodic_variation_new = periodic_variation
            randomize_phase_shift_new = randomize_phase_shift
            phase_shift_angle_new = phase_shift_angle 
            
            ! Read parameters defined in the old namelist <kt_grids_box_parameters>
            read (unit=unit_number_temp, nml=kt_grids_box_parameters)
            
            ! Check if any variable from the old namelist <kt_grids_box_parameters> has changed
            if (nx_new /= nx) double_definitions = .true.
            if (ny_new /= ny) double_definitions = .true.
            if (jtwist_new /= jtwist) double_definitions = .true.
            if (jtwistfac_new /= jtwistfac) double_definitions = .true.
            if (x0_new /= x0) double_definitions = .true.
            if (y0_new /= y0) double_definitions = .true.
            if (centered_in_rho_new .neqv. centered_in_rho) double_definitions = .true.
            if (periodic_variation_new .neqv. periodic_variation) double_definitions = .true.
            if (randomize_phase_shift_new .neqv. randomize_phase_shift) double_definitions = .true.
            if (phase_shift_angle_new /= phase_shift_angle) double_definitions = .true.

            ! If we have double defintions, abort stella
            if (double_definitions) call abort_because_of_double_definitions('kxky_grids_box', 'kt_grids_box_parameters')
            
         end if 
         
      end subroutine backwards_compatibility_kt_grids_box_parameters
      
   end subroutine write_kxky_grids_box
   
   !----------------------------------------------------------------------------
   !                         <kxky_grids_range> namelist                        
   !----------------------------------------------------------------------------
   subroutine write_kxky_grids_range(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use parameters_kxky_grids_range, only: set_defaults_for_kxky_grids_range => read_default_range
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! For <write_kxky_grids_range> the variables are not made public
      ! because only the variables in <parameters_kxky_grids> should be used
      ! in stella, so we need to define the namelist variables locally 
      integer :: naky, nakx
      real :: aky_min, aky_max, akx_min, akx_max, theta0_min, theta0_max
      character(20) :: kyspacing_option
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /kxky_grids_range/ naky, nakx, aky_min, aky_max, &
         theta0_min, theta0_max, akx_min, akx_max, kyspacing_option
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_kxky_grids_range'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_kxky_grids_range(naky, nakx, aky_min, aky_max, &
         akx_min, akx_max, theta0_min, theta0_max, kyspacing_option)
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("kxky_grids_range", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=kxky_grids_range) 
         call backwards_compatibility_kt_grids_range_parameters()
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=kxky_grids_range)
       
   contains
   
      !-------------------------------------------------------------------------
      !                  <kt_grids_range_parameters> namelist                   
      !-------------------------------------------------------------------------
      subroutine backwards_compatibility_kt_grids_range_parameters()
      
         implicit none 
      
         ! Save parameters which could potentially be defined twice
         integer :: naky_new, nakx_new
         real :: aky_min_new, aky_max_new, akx_min_new, akx_max_new
         real :: theta0_min_new, theta0_max_new
         character(20) :: kyspacing_option_new

         ! Local parameters for this subroutine
         logical :: double_definitions = .false.
         logical :: old_nml_exist
         
         ! Deprecated namelist <kt_grids_range_parameters>
         namelist /kt_grids_range_parameters/ naky, nakx, aky_min, aky_max, &
            theta0_min, theta0_max, akx_min, akx_max, kyspacing_option
            
         !----------------------------------------------------------------------

         ! Check whether the old knob exists
         unit_number_temp = input_unit_exist("kt_grids_range_parameters", old_nml_exist)
         
         ! If only the old namelist is present in <run_name>.in, and not the new
         ! namelist, simply read it in the old namelist
         if (old_nml_exist .and. (.not. new_nml_exist)) then
            read (unit=unit_number_temp, nml=kt_grids_range_parameters)
         end if 
         
         ! If both namelists exist, make sure the same variable is not mentioned twice
         if (old_nml_exist .and. new_nml_exist) then
         
            ! Print warning, because the input file should be updated.
            ! We do not want to abort, since we want stella to run the automatic
            ! tests on all versions of stella, including old versions.
            if (print_input_file_warnings) then 
               write(*,*) '   '
               write(*,*) 'WARNING: The <kt_grids_range_parameters> namelist has been deprecated. Please move the '
               write(*,*) 'variables [naky, nakx, aky_min, aky_max, theta0_min, theta0_max, akx_min, akx_max, '
               write(*,*) 'kyspacing_option] to the <kxky_grids_range> namelist.'
               write(*,*) '   '
            end if
            
            ! Save parameters defined in the new namelist <kxky_grids_range>
            naky_new = naky
            nakx_new = nakx
            aky_min_new = aky_min
            aky_max_new = aky_max
            theta0_min_new = theta0_min
            theta0_max_new = theta0_max
            akx_min_new = akx_min
            akx_max_new = akx_max
            kyspacing_option_new = kyspacing_option
            
            ! Read parameters defined in the old namelist <kt_grids_range_parameters>
            read (unit=unit_number_temp, nml=kt_grids_range_parameters)
            
            ! Check if any variable from the old namelist <kt_grids_range_parameters> has changed
            if (naky_new /= naky) double_definitions = .true.
            if (nakx_new /= nakx) double_definitions = .true.
            if (aky_min_new /= aky_min) double_definitions = .true.
            if (aky_max_new /= aky_max) double_definitions = .true.
            if (theta0_min_new /= theta0_min) double_definitions = .true.
            if (theta0_max_new /= theta0_max) double_definitions = .true.
            if (akx_min_new /= akx_min) double_definitions = .true.
            if (akx_max_new /= akx_max) double_definitions = .true.
            if (kyspacing_option_new /= kyspacing_option) double_definitions = .true.
            
            ! If we have double defintions, abort stella
            if (double_definitions) call abort_because_of_double_definitions('kxky_grids_range', 'kt_grids_range_parameters')
            
         end if 
         
      end subroutine backwards_compatibility_kt_grids_range_parameters
      
   end subroutine write_kxky_grids_range

!###############################################################################
!############################### Z VPA MU GRIDS ################################
!###############################################################################

   !----------------------------------------------------------------------------
   !                        <zgrid_parameters> namelist                         
   !----------------------------------------------------------------------------
   subroutine write_zgrid_parameters(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use zgrid, only: set_defaults_for_zgrid_parameters => set_default_parameters
      use zgrid, only: nzed, nperiod, ntubes, shat_zero, dkx_over_dky
      use zgrid, only: boundary_option, zed_equal_arc, grad_x_grad_y_zero
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /zgrid_parameters/ nzed, nperiod, ntubes, shat_zero, &
         boundary_option, zed_equal_arc, grad_x_grad_y_zero, dkx_over_dky 
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_zgrid_parameters'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_zgrid_parameters()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("zgrid_parameters", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=zgrid_parameters) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=zgrid_parameters)
      
   end subroutine write_zgrid_parameters

   !----------------------------------------------------------------------------
   !                      <vpamu_grids_parameters> namelist                     
   !---------------------------------------------------------------------------- 
   subroutine write_vpamu_grids_parameters(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use vpamu_grids, only: set_defaults_for_vpamu_grids_parameters => set_default_parameters
      use vpamu_grids, only: nvgrid, nmu, vpa_max, vperp_max, equally_spaced_mu_grid, conservative_wgts_vpa
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /vpamu_grids_parameters/ nvgrid, nmu, vpa_max, vperp_max, &
         equally_spaced_mu_grid, conservative_wgts_vpa
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_vpamu_grids_parameters'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_vpamu_grids_parameters()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("vpamu_grids_parameters", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=vpamu_grids_parameters) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=vpamu_grids_parameters)
      
   end subroutine write_vpamu_grids_parameters
   

!###############################################################################
!################################## SPECIES ####################################
!###############################################################################

   !----------------------------------------------------------------------------
   !                          <species_knobs> namelist                          
   !----------------------------------------------------------------------------

   subroutine write_species_knobs(unit_number, write_input_file, nspec_local)

      use file_utils, only: input_unit_exist
      use species, only: set_defaults_for_species_knobs => set_default_parameters
      use species, only: nspec, species_option, read_profile_variation
      use species, only: write_profile_variation, ecoll_zeff
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! For <species_parameters_index> we need <nspec>
      integer, intent(out) :: nspec_local
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /species_knobs/ nspec, species_option, read_profile_variation, &
         write_profile_variation, ecoll_zeff
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_species_knobs'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_species_knobs()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("species_knobs", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=species_knobs) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=species_knobs)
      
      ! Save <nspec>
      nspec_local = nspec
      
   end subroutine write_species_knobs
   

   !----------------------------------------------------------------------------
   !                       <species_parameters> namelist                        
   !----------------------------------------------------------------------------

   subroutine write_species_parameters(unit_number, write_input_file, index)

      use file_utils, only: input_unit_exist
      use species, only: set_defaults_for_species_parameters => set_default_parameters_per_specie
      use species, only: z, mass, dens, temp, tprim, fprim, d2ndr2, d2Tdr2, bess_fac, type
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      integer, intent(in) :: index
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /species_parameters_1/ z, mass, dens, temp, &
         tprim, fprim, d2ndr2, d2Tdr2, bess_fac, type
      namelist /species_parameters_2/ z, mass, dens, temp, &
         tprim, fprim, d2ndr2, d2Tdr2, bess_fac, type
      namelist /species_parameters_3/ z, mass, dens, temp, &
         tprim, fprim, d2ndr2, d2Tdr2, bess_fac, type
      namelist /species_parameters_4/ z, mass, dens, temp, &
         tprim, fprim, d2ndr2, d2Tdr2, bess_fac, type
      namelist /species_parameters_5/ z, mass, dens, temp, &
         tprim, fprim, d2ndr2, d2Tdr2, bess_fac, type
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_species_parameters'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_species_parameters()
   
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         if (index==1) then 
            unit_number_temp = input_unit_exist("species_parameters_1", new_nml_exist)  
            if (new_nml_exist) read (unit=unit_number_temp, nml=species_parameters_1)
         end if
         if (index==2) then
            unit_number_temp = input_unit_exist("species_parameters_2", new_nml_exist)  
            if (new_nml_exist) read (unit=unit_number_temp, nml=species_parameters_2)
         end if
         if (index==3) then
            unit_number_temp = input_unit_exist("species_parameters_3", new_nml_exist)  
            if (new_nml_exist) read (unit=unit_number_temp, nml=species_parameters_3)
         end if
         if (index==4) then
            unit_number_temp = input_unit_exist("species_parameters_4", new_nml_exist)  
            if (new_nml_exist) read (unit=unit_number_temp, nml=species_parameters_4)
         end if
         if (index==5) then
            unit_number_temp = input_unit_exist("species_parameters_5", new_nml_exist)  
            if (new_nml_exist) read (unit=unit_number_temp, nml=species_parameters_5)
         end if
         if (index>5) then 
            write(*,*) ' '
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ABORT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
            write(*,*) 'More than 5 species has not been implemented in update_input_file.f90.' 
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ABORT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
            write(*,*) ' '; 
          stop 
         end if  
      end if 
      
      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      if (index==1) write(unit=unit_number, nml=species_parameters_1)
      if (index==2) write(unit=unit_number, nml=species_parameters_2)
      if (index==3) write(unit=unit_number, nml=species_parameters_3)
      if (index==4) write(unit=unit_number, nml=species_parameters_4)
      if (index==5) write(unit=unit_number, nml=species_parameters_5)
      
   end subroutine write_species_parameters
   

!###############################################################################
!################################ DISTRIBUTION #################################
!###############################################################################

   !----------------------------------------------------------------------------
   !                           <init_g_knobs> namelist                          
   !----------------------------------------------------------------------------

   subroutine write_init_g_knobs(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use init_g, only: set_defaults_for_init_g_knobs => set_default_parameters
      use init_g, only: ginit_option, width0, phiinit, chop_side
      use init_g, only: restart_file, restart_dir, left, scale, tstart, zf_init
      use init_g, only: den0, upar0, tpar0, tperp0, imfac, refac
      use init_g, only: den1, upar1, tpar1, tperp1, den2, upar2, tpar2, tperp2
      use init_g, only: kxmax, kxmin, scale_to_phiinit, oddparity
      use stella_save, only: read_many
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /init_g_knobs/ ginit_option, width0, phiinit, chop_side, &
         restart_file, restart_dir, read_many, left, scale, tstart, zf_init, &
         den0, upar0, tpar0, tperp0, imfac, refac, &
         den1, upar1, tpar1, tperp1, den2, upar2, tpar2, tperp2, &
         kxmax, kxmin, scale_to_phiinit, oddparity
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_init_g_knobs'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_init_g_knobs()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("init_g_knobs", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=init_g_knobs) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=init_g_knobs)
      
   end subroutine write_init_g_knobs
   
!###############################################################################
!############################## PARALLELISATION ################################
!###############################################################################
   
   !----------------------------------------------------------------------------
   !                          <layouts_knobs> namelist                          
   !----------------------------------------------------------------------------
   subroutine write_layouts_knobs(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use stella_layouts, only: set_defaults_for_layouts_knobs => set_default_parameters
      use stella_layouts, only: xyzs_layout, vms_layout
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /layouts_knobs/ xyzs_layout, vms_layout
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_layouts_knobs'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_layouts_knobs()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("layouts_knobs", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=layouts_knobs) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=layouts_knobs)
      
   end subroutine write_layouts_knobs
   

!###############################################################################
!############################## RADIAL VARIATION ###############################
!###############################################################################

   !----------------------------------------------------------------------------
   !                       <multibox_parameters> namelist                       
   !----------------------------------------------------------------------------
   
   
   subroutine write_multibox_parameters(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use multibox, only: set_defaults_for_multibox_parameters => set_default_parameters
      use multibox, only: smooth_ZFs, zf_option, LR_debug_option, krook_option, RK_step, nu_krook_mb, mb_debug_step
      use multibox, only: krook_exponent, comm_at_init, phi_bound, phi_pow, krook_efold, use_dirichlet_BC
      use grids_kxky, only: boundary_size, krook_size
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! TODO: <boundary_size> belong to multibox namelist but is a public variable of grids_kxky.f90
      
      ! Namelist
      namelist /multibox_parameters/ boundary_size, krook_size, &
         smooth_ZFs, zf_option, LR_debug_option, krook_option, RK_step, nu_krook_mb, &
         mb_debug_step, krook_exponent, comm_at_init, phi_bound, phi_pow, krook_efold, use_dirichlet_BC
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_multibox_parameters'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_multibox_parameters()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("multibox_parameters", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=multibox_parameters) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=multibox_parameters)
      
   end subroutine write_multibox_parameters

!###############################################################################
!################################# DISSIPATION #################################
!###############################################################################

   !----------------------------------------------------------------------------
   !                           <dissipation> namelist                           
   !----------------------------------------------------------------------------
   subroutine write_dissipation(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use namelist_dissipation, only: read_dissipation_namelist => read_namelist
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist parameters
      logical :: include_collisions, collisions_implicit, hyper_dissipation
      character(30) :: collision_model
      
      ! Namelist
      namelist /dissipation/ include_collisions, collisions_implicit, collision_model, hyper_dissipation
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_dissipation'
   
      ! Read the <dissipation> namelist
      call read_dissipation_namelist(include_collisions, collisions_implicit, collision_model, hyper_dissipation)

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=dissipation)
      
   end subroutine write_dissipation

   !----------------------------------------------------------------------------
   !                          <collisions_fp> namelist                          
   !----------------------------------------------------------------------------
   subroutine write_collisions_fp(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use coll_fokkerplanck, only: set_defaults_for_collisions_fp => set_default_parameters
      use coll_fokkerplanck, only: testpart, fieldpart, lmax, jmax, nvel_local, interspec, intraspec
      use coll_fokkerplanck, only: iiknob, ieknob, eeknob, eiknob, eiediffknob, eideflknob, deflknob, eimassr_approx, advfield_coll, spitzer_problem
      use coll_fokkerplanck, only: density_conservation, density_conservation_field, density_conservation_tp, exact_conservation, exact_conservation_tp
      use coll_fokkerplanck, only: vpa_operator, mu_operator, cfac, cfac2, nuxfac, i1fac, i2fac, no_j1l1, no_j1l2, no_j0l2
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /collisions_fp/ testpart, fieldpart, lmax, jmax, nvel_local, &
         interspec, intraspec, iiknob, ieknob, eeknob, eiknob, eiediffknob, eideflknob, deflknob, eimassr_approx, advfield_coll, spitzer_problem, &
         density_conservation, density_conservation_field, density_conservation_tp, exact_conservation, exact_conservation_tp, &
         vpa_operator, mu_operator, cfac, cfac2, nuxfac, i1fac, i2fac, no_j1l1, no_j1l2, no_j0l2
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_collisions_fp'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_collisions_fp()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("collisions_fp", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=collisions_fp) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=collisions_fp)
      
   end subroutine write_collisions_fp

   !----------------------------------------------------------------------------
   !                       <collisions_dougherty> namelist                      
   !----------------------------------------------------------------------------
   subroutine write_collisions_dougherty(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use coll_dougherty, only: set_defaults_for_collisions_dougherty => set_default_parameters
      use coll_dougherty, only: momentum_conservation, energy_conservation, vpa_operator, mu_operator
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /collisions_dougherty/ momentum_conservation, energy_conservation, vpa_operator, mu_operator
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_collisions_dougherty'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_collisions_dougherty()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("collisions_dougherty", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=collisions_dougherty) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=collisions_dougherty)
      
   end subroutine write_collisions_dougherty

   !----------------------------------------------------------------------------
   !                             <hyper> namelist                               
   !----------------------------------------------------------------------------
   subroutine write_hyper(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use hyper, only: set_defaults_for_hyper => set_default_parameters
      use hyper, only: D_hyper, D_zed, D_vpa, hyp_zed, hyp_vpa, use_physical_ksqr, scale_to_outboard
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /hyper/ D_hyper, D_zed, D_vpa, hyp_zed, hyp_vpa, use_physical_ksqr, scale_to_outboard
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_hyper'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_hyper()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("hyper", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=hyper) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=hyper)
      
   end subroutine write_hyper

!###############################################################################
!################################# NEOCLASSICS #################################
!###############################################################################

   !----------------------------------------------------------------------------
   !                             <sources> namelist                             
   !----------------------------------------------------------------------------
   subroutine write_sources(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use sources, only: set_defaults_for_sources => set_default_parameters
      use sources, only: source_option, nu_krook, tcorr_source
      use sources, only: ikxmax_source, krook_odd, exclude_boundary_regions
      use sources, only: from_zero, conserve_momentum, conserve_density 
      use arrays_fields, only: tcorr_source_qn, exclude_boundary_regions_qn
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist 
      
      ! Namelist
      namelist /sources/ source_option, nu_krook, tcorr_source, &
         ikxmax_source, krook_odd, exclude_boundary_regions, &
         tcorr_source_qn, exclude_boundary_regions_qn, from_zero, &
         conserve_momentum, conserve_density
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_sources'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_sources()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("sources", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=sources) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=sources)
      
   end subroutine write_sources
   


   !----------------------------------------------------------------------------
   !                         <sfincs_input> namelist                            
   !----------------------------------------------------------------------------
  
   ! Sfincs is turned on with the fpp processor, skip this for now

   !----------------------------------------------------------------------------
   !                       <euterpe_parameters> namelist                        
   !----------------------------------------------------------------------------
   subroutine write_euterpe_parameters(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use euterpe_interface, only: set_defaults_for_euterpe_parameters => set_default_parameters
      use euterpe_interface, only: nradii, data_file
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /euterpe_parameters/ nradii, data_file
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_euterpe_parameters'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_euterpe_parameters()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("euterpe_parameters", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=euterpe_parameters) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=euterpe_parameters)
      
   end subroutine write_euterpe_parameters

   !----------------------------------------------------------------------------
   !                       <neoclassical_input> namelist                        
   !----------------------------------------------------------------------------
   subroutine write_neoclassical_input(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use neoclassical_terms, only: set_defaults_for_neoclassical_input => set_default_parameters
      use neoclassical_terms, only: include_neoclassical_terms, neo_option, nradii, drho
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist
      
      ! Namelist
      namelist /neoclassical_input/ include_neoclassical_terms, neo_option, nradii, drho
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_neoclassical_input'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_neoclassical_input()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("neoclassical_input", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=neoclassical_input) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=neoclassical_input)
      
   end subroutine write_neoclassical_input
   
   
   !----------------------------------------------------------------------------
   !                           <debug_flags> namelist                           
   !----------------------------------------------------------------------------
   subroutine write_debug_flags(unit_number, write_input_file)

      use file_utils, only: input_unit_exist
      use debug_flags, only: set_defaults_for_debug_flags => set_default_parameters
      use debug_flags, only: debug_all, stella_debug, ffs_solve_debug, fields_all_debug, fields_debug
      use debug_flags, only: fields_fluxtube_debug, fields_electromagnetic_debug, fields_ffs_debug
      use debug_flags, only: implicit_solve_debug, parallel_streaming_debug, mirror_terms_debug, neoclassical_terms_debug
      use debug_flags, only: response_matrix_debug, time_advance_debug, extended_grid_debug
      use debug_flags, only: diagnostics_all_debug, diagnostics_parameters, diagnostics_fluxes_fluxtube_debug
      use debug_flags, only: diagnostics_omega_debug, diagnostics_debug, dist_fn_debug
      use debug_flags, only: gyro_averages_debug, fluxes_debug, geometry_debug,  const_alpha_geo
   
      implicit none
      
      ! Either we write the real input file, or we write only the stella defaults,
      ! <unit_number> will point to <run_name>_with_defaults.in or default_stella_input.in
      integer, intent(in) :: unit_number
      logical, intent(in) :: write_input_file 
      
      ! Local parameters for this subroutine
      integer :: unit_number_temp
      logical :: new_nml_exist 
      
      ! Namelist
      namelist /debug_flags/ debug_all, stella_debug, ffs_solve_debug, fields_all_debug, fields_debug, &
        fields_fluxtube_debug, fields_electromagnetic_debug, fields_ffs_debug, & 
        implicit_solve_debug, parallel_streaming_debug, mirror_terms_debug, neoclassical_terms_debug, &
        response_matrix_debug, time_advance_debug, extended_grid_debug, &
        diagnostics_all_debug, diagnostics_parameters, diagnostics_fluxes_fluxtube_debug, &
        diagnostics_omega_debug, diagnostics_debug, dist_fn_debug,&
        gyro_averages_debug, fluxes_debug, geometry_debug,  const_alpha_geo
      
      !------------------------------------------------------------------------- 
      
      if (debug) write(*,*) 'update_input_file::write_debug_flags'
   
      ! Read in the default parameters set in stella 
      call set_defaults_for_debug_flags()
      
      ! Read the user-specified input parameters in <run_name>.in
      if (write_input_file) then
         unit_number_temp = input_unit_exist("debug_flags", new_nml_exist) 
         if (new_nml_exist) read (unit=unit_number_temp, nml=debug_flags) 
      end if  

      ! Write the namelist to <run_name>_with_defaults.in or default_stella_input.in
      write(unit=unit_number, nml=debug_flags)
      
   end subroutine write_debug_flags
   
!###############################################################################
!########################### DEFAULT_STELLA_INPUT.IN ###########################
!###############################################################################
   
   ! Write all default variables used in stella to an input file
   subroutine write_default_input_file()

      use file_utils, only: get_unused_unit

      implicit none 
      
      
      ! Local parameters
      character(23) :: input_file_name = 'default_stella_input.in'
      integer :: unit_number, dummy
      
      !-------------------------------------------------------------------------
      
      ! Message
      write(*,*) " "
      write(*,*) "Write all stella default input variables to 'default_stella_input.in'."
      write(*,*) " " 
       
      ! Open the 'default_stella_input.in' file
      call get_unused_unit(unit_number)
      open(unit=unit_number, file=input_file_name) 
         
      ! Write the namelists to "default_stella_input.in" file
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "!-------------------------------"
      write(unit_number, "(A)") "!----------  Geometry ----------"
      write(unit_number, "(A)") "!-------------------------------"
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "! Define geometry parameters"
      write(unit_number, "(A)") "! Depending on <geo_option>, different extra namelists will be read: "
      write(unit_number, "(A)") "!    <geo_option> = ['local', 'default', 'miller'] uses the <millergeo_parameters> namelist "
      write(unit_number, "(A)") "!    <geo_option> = 'vmec' uses the <vmec_parameters> namelist "
      write(unit_number, "(A)") "!    <geo_option> = 'zpinch' does not use any extra namelists " 
      write(unit_number, "(A)") "!    <geo_option> = 'input.profiles' does not use any extra namelists "
      write(unit_number, "(A)") "! Note that q_as_x = <radial_variation>, from the <parameters_physics> namelist"
      call write_geometry(unit_number, .false.) 
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "! Define parameters related to VMEC geometry"
      call write_vmec_geometry(unit_number, .false.) 
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "! Define parameters related to Miller geometry"
      call write_miller_geometry(unit_number, .false.) 
      
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "!-------------------------------"
      write(unit_number, "(A)") "!----------  Physics  ----------"
      write(unit_number, "(A)") "!-------------------------------"
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "! Define parameters related to physics"
      call write_parameters_physics(unit_number, .false.)
      
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "!-------------------------------"
      write(unit_number, "(A)") "!-----------  Grids  -----------"
      write(unit_number, "(A)") "!-------------------------------"
      write(unit_number, "(A)") "   "
      call write_vpamu_grids_parameters(unit_number, .false.)
      call write_zgrid_parameters(unit_number, .false.)
      call write_species_knobs(unit_number, .false., dummy)
      call write_species_parameters(unit_number, .false.,1)
      call write_species_parameters(unit_number, .false.,2)
      call write_kxky_grids(unit_number, .false.)
      call write_kxky_grids_range(unit_number, .false.)
      call write_kxky_grids_box(unit_number, .false.)
      
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "!--------------------------------"
      write(unit_number, "(A)") "!---------  Diagnostics  --------"
      write(unit_number, "(A)") "!--------------------------------"
      write(unit_number, "(A)") "   "
      call write_diagnostics(unit_number, .false.)
      
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "!--------------------------------"
      write(unit_number, "(A)") "!--------  Init fields  ---------"
      write(unit_number, "(A)") "!--------------------------------"
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "! Note that <restart_file> = trim(run_name)//.nc"
      call write_init_g_knobs(unit_number, .false.)
      
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "!-------------------------------"
      write(unit_number, "(A)") "!--------  Dissipation  --------"
      write(unit_number, "(A)") "!-------------------------------"
      write(unit_number, "(A)") "   "
      call write_dissipation(unit_number, .false.)
      call write_collisions_dougherty(unit_number, .false.)
      call write_collisions_fp(unit_number, .false.)
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "! Note that <use_physical_ksqr> = .not. (full_flux_surface .or. radial_variation)"
      call write_hyper(unit_number, .false.)
      
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "!-------------------------------"
      write(unit_number, "(A)") "!----------  Numerics  ---------"
      write(unit_number, "(A)") "!-------------------------------"
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "! Define parameters related to numerics"
      call write_parameters_numerical(unit_number, .false.)
      
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "!-------------------------------"
      write(unit_number, "(A)") "!--------  Neoclassics  --------"
      write(unit_number, "(A)") "!-------------------------------"
      write(unit_number, "(A)") "   "
      call write_neoclassical_input(unit_number, .false.)
      call write_euterpe_parameters(unit_number, .false.)
      call write_sources(unit_number, .false.)
      
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "!--------------------------------"
      write(unit_number, "(A)") "!------  Radial variation  ------"
      write(unit_number, "(A)") "!--------------------------------"
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "! Note that <boundary_size> belongs to grids_kxky?"
      call write_multibox_parameters(unit_number, .false.)
      
      write(unit_number, "(A)") "   "
      write(unit_number, "(A)") "!--------------------------------"
      write(unit_number, "(A)") "!-------  Parallelisation  ------"
      write(unit_number, "(A)") "!--------------------------------"
      write(unit_number, "(A)") "   "
      call write_layouts_knobs(unit_number, .false.)
      call write_debug_flags(unit_number, .false.)

      ! Close the file
      close(unit=unit_number)
      
   end subroutine write_default_input_file
   
!###############################################################################
!################################### ABORT #####################################
!###############################################################################
   
   ! Stop the program when a variable appears in both the old and new namelist
   subroutine abort_because_of_double_definitions(new_namelist, old_namelist)
   
         implicit none
         
         character(*) :: new_namelist
         character(*) :: old_namelist
         
         !---------------------------------------------------------------------- 
         
         ! Write abort message  
         write(*,*) ' '
         write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ABORT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
         write(*,*) 'There is a parameter defined both in the deprecated namelist <'//trim(old_namelist)//'>'
         write(*,*) 'and in the new namelist <'//trim(new_namelist)//'>. Or both namelists exists and the '
         write(*,*) 'variable is only present in the old namelist, and changes the default value of stella.'
         write(*,*) 'The code can not determine which parameter to read, please update the input '
         write(*,*) 'file and remove the deprecated namelist <'//trim(old_namelist)//'>' 
         write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ABORT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write(*,*) ' ';
         
         ! Stop the program
         stop 
         
     end subroutine abort_because_of_double_definitions
   
   
!###############################################################################
!############################# COMMAND LINE FLAGS ##############################
!###############################################################################
! Parse some basic command line arguments. 
! Currently just 'version' and 'help'. 
    
   subroutine parse_command_line()
   
      use git_version, only: get_git_version

      implicit none
      
      integer :: arg_count, arg_n
      integer :: arg_length
      character(len=:), allocatable :: argument
      character(len=*), parameter :: endl = new_line('a')

      ! Get the number of arguments passed to the command prompt
      arg_count = command_argument_count()
      
      ! Read the argument on the command line, if it matches any of the defined
      ! options, then display some text, and stop the program.
      do arg_n = 0, arg_count
         call get_command_argument(1, length=arg_length)
         if (allocated(argument)) deallocate (argument)
         allocate (character(len=arg_length)::argument)
         call get_command_argument(1, argument)
         
         ! Check if the argument matches any of these flags
         if ((argument == "--version") .or. (argument == "-v")) then
            write (*, '("stella version ", a)') get_git_version()
            stop
         else if ((argument == "--default-input") .or. (argument == "-d")) then
            call write_default_input_file()
            stop
         else if ((argument == "--help") .or. (argument == "-h")) then
            write (*, '(a)') "  " //endl// &
               "stella [--version|-v] [--help|-h] [--default-input|-d] [input file]"//endl//endl// &
               "stella is a flux tube gyrokinetic code for micro-stability and turbulence "// &
               "simulations of strongly magnetised plasma"//endl// &
               "For more help, see the documentation at https://stellagk.github.io/stella/"//endl// &
               "or create an issue https://github.com/stellaGK/stella/issues/new"//endl// &
               endl// &
               "  -h, --help           Print this message"//endl// &
               "  -v, --version        Print the stella version"//endl// &
               "  -d, --default-input  Print default stella input file"//endl//" "
            stop
         end if
         
      end do
      
   end subroutine parse_command_line
   
end module input_file




! New suggestion
! namelist /parameters_adiabatic/ adiabatic_option, tite, nine
! namelist /parameters_gyrokinetic_equation/ include_parallel_streaming, include_mirror, nonlinear, &
!    include_parallel_nonlinearity, xdriftknob, ydriftknob, wstarknob, 
! namelist /parameters_electromagnetic/ include_apar, include_bpar, beta
! namelist /parameters_full_flux_surface/ full_flux_surface, rhostar, irhostar
! namelist /parameters_collisions/ zeff, vnew_ref
! namelist /parameters_unsorted/ prp_shear_enabled, hammett_flow_shear, include_pressure_variation
!    include_geometric_variation, suppress_zonal_interaction, radial_variation,
!    g_exb, g_exbfac, omprimfac, 
