!###############################################################################
!###################### STELLA: Delta-f gyrokinetic code #######################
!###############################################################################
! 
! Initialise all modules, arrays and parameters needed for stella.
! 
!###############################################################################
module init_stella

   ! Load the debug flags
   use debug_flags, only: debug => stella_debug
      
   implicit none
   
   ! Make routines available to stella.f90
   public :: initialise_stella
   public :: finish_stella
   public :: parse_command_line
   
   private
   
   ! Make sure we only initialise mpi once
   logical :: mpi_initialised = .false.

contains

!###############################################################################
!############################## INTIALISE STELLA ###############################
!###############################################################################
! Calls the initialisation routines for all the modules.
!###############################################################################

   !****************************************************************************
   !                       OVERVIEW OF INITIALISE STELLA                       !
   !****************************************************************************
   subroutine initialise_stella(istep0, istatus)

      ! Parallelisation
      use mp, only: proc0
      use job_manage, only: time_message
      use timers, only: time_init
      
      ! Time trace
      use grids_time, only: init_tstart
      use initialise_distribution_function, only: tstart
      
      ! Files
      use file_utils, only: flush_output_file
      use file_utils, only: error_unit
      use file_units, only: close_input_file_with_defaults
      use file_units, only: unit_input_file_with_defaults
      
      ! Init modules
      use parallelisation_layouts, only: init_dist_fn_layouts
      use initialise_redistribute, only: init_redistribute
      use diagnostics, only: init_diagnostics
      
      ! Calculations
      use calculations_volume_averages, only: volume_average
      use calculations_transforms, only: init_transforms
      
      ! Grids
      use grids_kxky, only: naky, nakx, ny, nx, nalpha
      use grids_velocity, only: nvgrid, nmu
      use grids_species, only: nspec
      use grids_z, only: nzgrid, ntubes
      
      implicit none

      ! Starting timestep: zero unless the simulation has been restarted
      integer, intent(out) :: istep0
      integer, intent(in out) :: istatus

      ! Local variables
      logical :: restarted
      logical :: fourier_transformations_are_needed

      !-------------------------------------------------------------------------
      
      ! It is important to first initialise the MPI environment, so that we can 
      ! parallelise computations over the <nproc> available processors, and so 
      ! that a specific processor is assigned to <proc0>. Moreover, we need to
      ! start the overall timer, and initialise the file utilities so we can
      ! read and write files. Once we can open files, we want to read the debug
      ! flags. If a list of input files needs to be launched, we divide the number
      ! of jobs (or input files) over the number of processors using job_fork().
      ! Finally, we start more internal timers, get the <run_name> and initialise
      ! the random nunmber generator based on <rng_seed> in the input file
      if (debug) write (*, *) 'stella::init_stella::init_mpi_files_utils_and_timers'
      call init_mpi_files_utils_and_timers

      ! Write message to command prompt with useful info regarding at the start of the simulation
      ! This routine uses <print_extra_info_to_terminal> from the debug_flags namelist
      if (debug) write (*, *) 'stella::init_stella::write_start_message'
      call write_start_message() 
      
      ! Read the input file
      if (debug) write (*, *) 'stella::init_stella::read_parameters_from_input_file'
      call read_parameters_from_input_file
      
      ! The grid points of the distribution function g(kx,ky,z,mu,vpa,species),
      ! are distributed among the various processors according to different layouts.
      if (debug) write (6, *) "stella::init_stella::init_dist_fn_layouts"
      call init_dist_fn_layouts(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny, nx, nalpha)
      
      ! Sets up the mappings between different layouts, needed
      ! to redistribute data when going from one layout to another
      if (debug) write (6, *) "stella::init_stella::init_redistribute"
      call init_redistribute
      
      ! Check whether Fourier transformations are needed (e.g., the nonlinear term
      ! requires it, radial variation requires it, ...) and if so, set up the FFTW plans
      if (debug) write (*, *) "stella::init_stella::init_transforms"
      call are_fourier_transformations_needed(fourier_transformations_are_needed)
      if (fourier_transformations_are_needed) then
         call init_transforms
      end if
      
      ! Initialise the (kx,ky,z,mu,vpa,species) grids as well as the magnetic geometry
      if (debug) write (*, *) "stella::init_stella::init_grids_and_geometry"
      call init_grids_and_geometry
      
      ! Initialise the arrays for the distribution function g(kx,ky,z,i[mu,vpa,species])
      ! as <gold> and <gnew> as well as gvmu(vpa,mu,i[kx,ky,z,species]). Also allocate
      ! arrays for vperp2 and kperp2 and arrays needed for gyro-averaging (j0 and j1).
      ! Moreover, initialise the redistribute, dissipation, soueces, fields, 
      ! distribution and volume_averaged modules
      if (debug) write (*, *) "stella::init_stella::init_arrays_and_stella_modules"
      call init_arrays_and_stella_modules(restarted, istep0)
      
      ! Initialise <delt>m the time advance module, and the response matrix
      if (debug) write (*, *) "stella::init_stella::init_time_step_and_gyrokinetic_equation"
      call init_time_step_and_gyrokinetic_equation(restarted, istatus)
     
      ! The distribution function g(kx,ky,z,mu,vpa,species) has been initialised
      ! in initialise_distribution(). Use the quasineutrality condition to initialise
      ! the electrostatic and electromagnetic fields (phi, apar, bpar).
      if (debug) write (*, *) "stella::init_stella::init_electrostatic_and_magnetic_potential"
      call init_electrostatic_and_magnetic_potential(restarted)

      ! Open ascii output files and initialise the netcdf file with extension .out.nc
      if (debug) write (6, *) 'stella::init_stella::init_diagnostics'
      call init_diagnostics(restarted, tstart, istep0)
      
      ! Initialise the code_time
      if (debug) write (6, *) 'stella::init_stella::init_tstart'
      call init_tstart(tstart)

      ! Make sure the error file and default input file are written
      if (debug) write (6, *) 'stella::init_stella::write_error_file'
      if (proc0) call flush_output_file(error_unit())
      if (proc0) call flush_output_file(unit_input_file_with_defaults)
      if (proc0) call close_input_file_with_defaults()

      ! Add a header to the output file
      if (debug) write (6, *) 'stella::init_stella::print_header'
      call print_header
      
      ! Stop the timing of the Initialisation
      if (proc0) call time_message(.false., time_init, ' Initialisation')

   end subroutine initialise_stella
   
   !****************************************************************************
   !              Intialise MPI environment, file utils and timers              
   !****************************************************************************
   subroutine init_mpi_files_utils_and_timers
   
      ! After we initialised the MPI environment, we can access <proc0> and <broadcast>
      use mp, only: init_mp
      use mp, only: proc0
      use mp, only: broadcast
      
      ! Start timer, so stella exits 5 minutes before <avail_cpu_time>
      use job_manage, only: check_cpu_time

      ! Initialise file utils
      use file_utils, only: init_file_utils
      use file_units, only: init_file_units
      use file_units, only: open_input_file_with_defaults
      
      ! Read debug flags
      use debug_flags, only: read_debug_flags
      
      ! Divide input files over processors
      use job_manage, only: job_fork
      use file_utils, only: runtype_option_switch
      
      ! Start more timers
      use job_manage, only: time_message
      use timers, only: time_total
      use timers, only: time_init
      
      ! Assign the <run_name> to each job
      use file_utils, only: run_name
      use file_utils, only: init_job_name
      
      ! Initialise the random number generator
      use interface_random_number_generator, only: init_random_number_generator
      
      implicit none
      
      character(500), target :: cbuff
      real :: avail_cpu_time_dummy = 10000000000.0
      logical :: list, exit
      
      !----------------------------------------------------------------------

      ! Initialise MPI (Message Passing Interface) which allows us to parallelise
      ! computations over the <nproc> available processors. It will assign a <proc0>.
      if (.not. mpi_initialised) call init_mp
      mpi_initialised = .true.

      ! Initialise timer, so that stella exits 5 minutes before <avail_cpu_time>
      if (debug) write (*, *) 'stella::init_stella::check_time'
      call check_cpu_time(avail_cpu_time_dummy, exit)

      ! Initialise file utilities, this will open the input and error files
      if (proc0) then 
         if (debug) write (*, *) 'stella::init_stella::init_file_utils'
         call init_file_utils(list)
         if (debug) write (*, *) 'stella::init_stella::init_file_units'
         call init_file_units
         if (debug) write (*, *) 'stella::init_stella::open_input_file_with_default'
         call open_input_file_with_defaults
      end if

      ! Read the debug flags, since the input file is now available
      ! Note that any previous debug statements are not printed since debug = False
      ! by default, and it's value isn't changed until we read the debug flags below
      call read_debug_flags

      ! We can launch a list of input files which are divided over the processors
      call broadcast(list)
      call broadcast(runtype_option_switch)
      if (list) call job_fork

      ! Start internal timers, this is part of job_manage and needs to be called after job_fork
      if (proc0) then
         call time_message(.false., time_total, ' Total')
         call time_message(.false., time_init, ' Initialisation')
      end if

      ! Assign a <run_name> to each job, generally this is the name of the input file
      if (proc0) cbuff = trim(run_name)
      call broadcast(cbuff)
      if (.not. proc0) call init_job_name(cbuff)
      
      ! Set up the random number generator (rng) in the ran.fpp module. It will read 
      ! <rng_seed> from the the "initialise_distribution_noise" namelist in the input file
      if (debug) write (6, *) "stella::init_stella::init_random_number_generator"
      call init_random_number_generator
      
   end subroutine init_mpi_files_utils_and_timers

   !****************************************************************************
   !                              Read input file                               
   !****************************************************************************
   subroutine read_parameters_from_input_file
      
      ! Parameters 
      use parameters_physics, only: read_parameters_physics
      use parameters_numerical, only: read_parameters_numerical
      use parameters_diagnostics, only: read_parameters_diagnostics
      use parameters_multibox, only: read_parameters_multibox
      
      ! Parameters from grids
      use grids_kxky, only: read_parameters_kxky_grids
      use grids_species, only: read_parameters_species
      use grids_z, only: read_parameters_z_grid
      use grids_velocity, only: read_parameters_velocity_grids
      
      ! Other parameters
      use initialise_distribution_function, only: read_parameters_distribution_function
      use parallelisation_layouts, only: read_parameters_parallelisation_layouts
      use dissipation_and_collisions, only: read_parameters_dissipation_and_collisions
      use gk_flow_shear, only: read_parameters_flow_shear
      
      ! Parse collision variables to read_parameters_species to avoid circular dependencies
      use dissipation_and_collisions, only: vnew_ref
      
      implicit none
      
      !----------------------------------------------------------------------
   
      ! Read the physics and numerical parameters from the input file
      ! These namelists contain many variables used by other modules, so read it first
      if (debug) write (6, *) "stella::init_stella::read_parameters_physics"
      call read_parameters_physics
      if (debug) write (6, *) "stella::init_stella::read_parameters_numerical"
      call read_parameters_numerical
      
      ! Read the dissipation and collision variables, since <vnew_ref> is needed
      ! for the species parameters read from the EUTERPE code
      if (debug) write (6, *) 'stella::init_stella::read_parameters_dissipation_and_collisions'
      call read_parameters_dissipation_and_collisions
      
      ! Read the namelists related to the (kx,ky,z,mu,vpa,species) grids
      if (debug) write (6, *) "stella::init_stella::read_parameters_z_grid"
      call read_parameters_z_grid
      if (debug) write (6, *) "stella::init_stella::read_species_options"
      call read_parameters_species(vnew_ref)
      if (debug) write (6, *) "stella::init_stella::read_parameters_kxky_grids"
      call read_parameters_kxky_grids
      if (debug) write (6, *) "stella::init_stella::read_velocity_grids_parameters"
      call read_parameters_velocity_grids
      
      ! Read remaining parameters
      if (debug) write (6, *) "stella::init_stella::read_parameters_flow_shear"
      call read_parameters_flow_shear
      if (debug) write (6, *) "stella::init_stella::read_parameters_multibox"
      call read_parameters_multibox
      if (debug) write (6, *) "stella::init_stella::read_parameters_diagnostics"
      call read_parameters_diagnostics
      if (debug) write (6, *) "stella::init_stella::read_parameters_distribution_function"
      call read_parameters_distribution_function
      if (debug) write (6, *) 'stella::init_stella::read_parameters_parallelisation_layouts'
      call read_parameters_parallelisation_layouts
      
   end subroutine read_parameters_from_input_file
   
   !****************************************************************************
   !                   Initialise grids and magnetic geometry                   
   !****************************************************************************
   subroutine init_grids_and_geometry
   
      ! Intialise grids
      use grids_z, only: init_z_grid
      use grids_kxky, only: init_grids_kxky
      use grids_velocity, only: init_velocity_grids
      use grids_species, only: init_species
      use grids_extended_zgrid, only: init_extended_zgrid
   
      ! Initialise geometry
      use geometry, only: init_geometry
      
      ! The geometry module needs <nalpha> and <naky> which have already been read
      ! by read_parameters_kxky_grids() in read_parameters_from_input_file()
      use grids_kxky, only: naky, nalpha
      
      ! Parse collision variables to init_species to avoid circular dependencies
      use dissipation_and_collisions, only: ecoll_zeff
      use dissipation_and_collisions, only: zeff
      use dissipation_and_collisions, only: vnew_ref
   
      implicit none
      
      !----------------------------------------------------------------------
   
      ! Set-up the z-grid, and calculate all of the required geometric coefficients
      ! Note that the geometry namelists will be read from the input file within init_geometry
      if (debug) write (6, *) "stella::init_stella::init_z_grid"
      call init_z_grid
      if (debug) write (6, *) "stella::init_stella::init_geometry"
      call init_geometry(nalpha, naky)
      
      ! The (kx,ky) grids require <shat> and <rhotor> from the geometry module,
      ! so make sure to initialise the geometry before initialising the (kx,ky) grids
      if (debug) write (6, *) 'stella::init_stella::init_grids_kxky'
      call init_grids_kxky
      
      ! Read species_parameters from input file and use the info to, e.g.,
      ! determine if a modified Boltzmann response is to be used
      ! Note that <ecoll_zeff> has already been read in read_parameters_dissipation_and_collisions()
      if (debug) write (6, *) 'stella::init_stella::init_species'
      call init_species(ecoll_zeff, zeff, vnew_ref)
      
      ! Setup the (kx,ky) grids and (x,y) grids, if applicable
      if (debug) write (6, *) 'stella::init_stella::init_multibox_subcalls'
      call init_multibox_subcalls
      
      ! Setup the (vpa,mu) grids and associated integration weights
      if (debug) write (6, *) 'stella::init_stella::init_velocity_grids'
      call init_velocity_grids(vnew_ref)
      
      ! Set up all of the logic needed to do calculations on an extended grid in z.
      ! this extended grid could be due to use of a ballooning angle so that
      ! z goes from -N*pi to N*pi, or it could be due to the coupling of different
      ! kx modes arising from the twist-and-shift boundary condition
      if (debug) write (6, *) 'stella::init_stella::init_extended_zgrid'
      call init_extended_zgrid
      
   end subroutine init_grids_and_geometry

   !****************************************************************************
   !                    Intialise arrays and stella modules                     
   !****************************************************************************
   subroutine init_arrays_and_stella_modules(restarted, istep0)
   
      ! Initialise arrays
      use initialise_arrays, only: init_arrays_vperp_kperp
      use arrays_gyro_averages, only: init_arrays_bessel_functions
   
      ! Initialise other modules
      use initialise_distribution_function, only: init_distribution_function
      use calculations_volume_averages, only: init_volume_averages
      use dissipation_and_collisions, only: init_dissipation
      use gk_sources, only: init_sources
      use field_equations, only: init_field_equations
      
      implicit none
      
      ! Arguments
      integer, intent(in out) :: istep0
      logical, intent(in out) :: restarted
      
      !----------------------------------------------------------------------
      
      ! When doing a volume average using Fourier coefficients, the
      ! ky=0 mode gets a different weighting than finite ky modes, due
      ! to the reality condition being imposed; init_volume_averages accounts for this
      if (debug) write (6, *) 'stella::init_stella::init_volume_averages'
      call init_volume_averages
      
      ! Initialise the arrays for the distribution function g(kx,ky,z,i[mu,vpa,species])
      ! as <gold> and <gnew> as well as gvmu(vpa,mu,i[kx,ky,z,species]).
      ! Aloso allocate vperp2, kperp2 and arrays needed for gyro-averaging (j0 and j1)
      if (debug) write (6, *) "stella::init_stella::init_arrays"
      call init_arrays_vperp_kperp
      call init_arrays_bessel_functions
      
      ! Read dissipation namelist from the input file and print information
      ! about chosen options to stdout
      if (debug) write (6, *) 'stella::init_stella::init_dissipation'
      call init_dissipation
      
      ! Initialise sources
      if (debug) write (6, *) 'stella::init_stella::init_sources'
      call init_sources
      
      ! Allocate and initialise time-independent arrays needed to
      ! solve the field equations; e.g., sum_s (Z_s^2 n_s / T_s)*(1-Gamma0_s)
      if (debug) write (6, *) 'stella::init_stella::init_field_equations'
      call init_field_equations
      
      ! Initialise the guiding-center distribution functions <gvmu>(nvpa, nmu, -kxkyzs-layout-)
      ! <gnew>(kx, ky, z, -vpamus-layout-) and <old>(kx, ky, z, -vpamus-layout-)
      if (debug) write (6, *) "stella::init_stella::init_distribution_function"
      call init_distribution_function(restarted, istep0)
      
   end subroutine init_arrays_and_stella_modules
   
   !****************************************************************************
   !                    Initialise time step and time advance                   
   !****************************************************************************
   subroutine init_time_step_and_gyrokinetic_equation(restarted, istatus)
      
      ! Set the time step
      use parameters_numerical, only: delt_option_switch
      use parameters_numerical, only: delt_option_auto
      use parameters_numerical, only: delt, delt_max, delt_min
      use grids_time, only: init_delt
      use save_stella_for_restart, only: init_dt
      
      ! Initialise parts of the gyrokinetic equation
      use gyrokinetic_equation_initialisation, only: init_gyrokinetic_equation
      use response_matrix, only: init_response_matrix
      use response_matrix, only: read_response_matrix
      
      ! Flags related to the gyrokinetic equation
      use parameters_numerical, only: stream_implicit
      use parameters_numerical, only: driftkinetic_implicit
      use parallelisation_layouts, only: mat_read
      
      implicit none
      
      ! Arguments
      logical, intent(in out) :: restarted
      integer, intent(in out) :: istatus
      
      ! Local variables
      real :: delt_saved
      
      !----------------------------------------------------------------------

      ! If initializing from restart file, set the initial time step size appropriately
      if (restarted .and. delt_option_switch == delt_option_auto) then
         delt_saved = delt
         if (debug) write (6, *) "stella::init_stella::init_dt"
         call init_dt(delt_saved, istatus)
         if (istatus == 0) delt = delt_saved
      end if
      
      ! Set the internal time step size variable code_dt from the input variable delt
      if (debug) write (6, *) "stella::init_stella::init_delt"
      call init_delt(delt, delt_max, delt_min)
      
      ! Allocate and calculate arrays needed for the mirror, parallel streaming,
      ! magnetic drifts, gradient drive, etc. terms during time advance
      if (debug) write (6, *) 'stella::init_stella::init_gyrokinetic_equation'
      call init_gyrokinetic_equation
      if (stream_implicit .or. driftkinetic_implicit) then
         if (mat_read) then
            if (debug) write (6, *) "stella::init_stella::read_response_matrix"
            call read_response_matrix
         else
            if (debug) write (6, *) "stella::init_stella::init_response_matrix"
            call init_response_matrix
         end if
      end if
      
   end subroutine init_time_step_and_gyrokinetic_equation
   
   !****************************************************************************
   !                        Initialise phi, apar and bpar                       
   !****************************************************************************
   subroutine init_electrostatic_and_magnetic_potential(restarted)
      
      ! The fields are phi(kx,ky,z), apar(kx,ky,z) and bpar(kx,ky,z)
      use field_equations, only: advance_fields
      use field_equations, only: fields_updated
      use field_equations, only: rescale_fields
      use initialise_distribution_function, only: phiinit
      use initialise_distribution_function, only: scale_to_phiinit
      
      ! Load the fields and the distribution function g(kx,ky,z,mu,vpa,species)
      use arrays_fields, only: phi, apar, bpar
      use arrays_distribution_function, only: gnew
      
      ! Radial variation runs
      use parameters_physics, only: radial_variation
      use file_utils, only: runtype_option_switch
      use file_utils, only: runtype_multibox
      use parameters_multibox, only: use_dirichlet_BC
      use multibox, only: apply_radial_boundary_conditions
      use multibox, only: multibox_communicate
      use field_equations_radialvariation, only: get_radial_correction
      
      ! Parallelisation
      use mp, only: job
      
      implicit none
      
      ! Arguments
      logical, intent(in out) :: restarted
      
      !----------------------------------------------------------------------

      ! The distribution function g(kx,ky,z,mu,vpa,species) has been initialised
      ! in initialise_distribution(). Use the quasineutrality condition to initialise 
      ! the electrostatic and electromagnetic fields (phi, apar, bpar).
      if (debug) write (6, *) 'stella::init_stella::advance_fields'
      call advance_fields(gnew, phi, apar, bpar, dist='g')
      
      ! Add the radial variation correction to the fields
      if (radial_variation) then
         if (debug) write (6, *) 'stella::init_stella::get_radial_correction'
         call get_radial_correction(gnew, phi, dist='g')
      end if

      ! Fill in the boundary regions using auxilliary simulations if using
      ! multibox, or zero it out if using Dirichlet boundary conditions
      if (runtype_option_switch == runtype_multibox) then
         if (debug) write (6, *) 'stella::init_stella:multibox_communicate'
         call multibox_communicate(gnew)
         if (job == 1) then
            fields_updated = .false.
            call advance_fields(gnew, phi, apar, bpar, dist='g')
         end if
      else if (use_dirichlet_BC) then
         if (debug) write (6, *) 'stella::init_stella:multibox_radial_BC'
         call apply_radial_boundary_conditions(gnew)
         fields_updated = .false.
         call advance_fields(gnew, phi, apar, bpar, dist='g')
      end if

      ! Rescale to phiinit if just beginning a new run
      if (.not. restarted .and. scale_to_phiinit) call rescale_fields(phiinit)
      
   end subroutine init_electrostatic_and_magnetic_potential
   
   !****************************************************************************
   !****************************************************************************
   !****************************************************************************
   ! Call all the multibox communication subroutines to make sure all the jobs 
   ! have the appropriate information
   !****************************************************************************
   subroutine init_multibox_subcalls

      use mp, only: proc0, job
      use grids_species, only: communicate_species_multibox
      use geometry, only: communicate_geo_multibox
      use calculations_kxky, only: communicate_ktgrids_multibox
      use file_utils, only: runtype_option_switch, runtype_multibox
      use parameters_physics, only: radial_variation
      use multibox, only: init_multibox, rhoL, rhoR
      use multibox, only: communicate_multibox_parameters, multibox_communicate

      implicit none
      
      !----------------------------------------------------------------------

      if (debug) write (6, *) 'stella::init_stella::init_multibox'
      call init_multibox
      if (runtype_option_switch == runtype_multibox) then
         if (proc0 .and. (job == 1) .and. radial_variation) then
            if (debug) write (6, *) 'stella::init_stella::init_multibox_geo'
            call communicate_geo_multibox(rhoL, rhoR)
            if (debug) write (6, *) 'stella::init_stella::init_multibox_spec'
            call communicate_species_multibox(rhoL, rhoR)
         end if
         if (job == 1) then
            call communicate_multibox_parameters
         end if
         if (radial_variation) then
            if (debug) write (6, *) 'stella::init_stella::init_multibox_ktgrid'
            call communicate_ktgrids_multibox
         end if
      end if

   end subroutine init_multibox_subcalls

!###############################################################################
!################################# CALCULATIONS ################################
!###############################################################################

   !****************************************************************************
   !                     Are Fourier Transformations needed                     
   !****************************************************************************
   ! Checks the various physics flag choicesto determine if FFTs are needed.
   !****************************************************************************
   subroutine are_fourier_transformations_needed(fourier_transformations_are_needed)

      use file_utils, only: runtype_option_switch, runtype_multibox
      use parameters_physics, only: include_nonlinear, include_parallel_nonlinearity
      use parameters_physics, only: radial_variation, full_flux_surface
      use gk_flow_shear, only: hammett_flow_shear
      use gk_flow_shear, only: g_exb, g_exbfac
      use parameters_diagnostics, only: write_radial_moments, write_radial_fluxes

      implicit none
      
      !----------------------------------------------------------------------

      logical, intent(out) :: fourier_transformations_are_needed

      ! Assume we don't need Fourier transformations
      fourier_transformations_are_needed = .false.
      
      ! If ExB or parallel nonlinearity included in the simulations, need FFTs
      if (include_nonlinear .or. include_parallel_nonlinearity) fourier_transformations_are_needed = .true.
      
      ! If 'global' in radial or bi-normal directions, need FFTs
      if (radial_variation .or. full_flux_surface) fourier_transformations_are_needed = .true.
      
      ! If running in multibox mode, need FFTs
      if (runtype_option_switch == runtype_multibox) fourier_transformations_are_needed = .true.
      
      ! If including flow shear using anything other than wavenumber re-mapping, need FFTs
      if (abs(g_exb * g_exbfac) > epsilon(0.) .and. .not. hammett_flow_shear) fourier_transformations_are_needed = .true.
      
      ! If printing out flux-surface-averaged radial fluxes or moments, need FFTs
      if (write_radial_fluxes .or. write_radial_moments) fourier_transformations_are_needed = .true.

   end subroutine are_fourier_transformations_needed

   !****************************************************************************
   !                       Write start message to screen                        
   !****************************************************************************
   subroutine write_start_message()
   
      use mp, only: proc0, nproc
      use debug_flags, only: print_extra_info_to_terminal
      use git_version, only: get_git_version, get_git_date

      implicit none

      ! Stella version number and release date
      character(len=40) :: git_commit
      character(len=10) :: git_date
      
      !----------------------------------------------------------------------
      
      ! Only print the header on the first processor
      if (.not. proc0) return
      
      ! Get git data
      git_commit = get_git_version()
      git_date = get_git_date()
      
      ! Print the stella header
      if (print_extra_info_to_terminal) then
         write (*, *) ' '
         write (*, *) ' '
         write (*, *) "              I8            ,dPYb, ,dPYb,            "
         write (*, *) "              I8            IP'`Yb IP'`Yb            "
         write (*, *) "           88888888         I8  8I I8  8I            "
         write (*, *) "              I8            I8  8' I8  8'            "
         write (*, *) "     ,g,      I8    ,ggg,   I8 dP  I8 dP    ,gggg,gg "
         write (*, *) "    ,8'8,     I8   i8' '8i  I8dP   I8dP    dP'  'Y8I "
         write (*, *) "   ,8'  Yb   ,I8,  I8, ,8I  I8P    I8P    i8'    ,8I "
         write (*, *) "  ,8'_   8) ,d88b, `YbadP' ,d8b,_ ,d8b,_ ,d8,   ,d8b,"
         write (*, *) '  P` "YY8P8P8P""Y8888P"Y8888P`"Y888P`"Y88P"Y8888P"`Y8'
         write (*, *) ' '
         write (*, *) ' '
         write (*, '(a48)') git_commit
         write (*, *) '                      ', git_date
         write (*, *) ' '
         write (*, *) '                   The stella team' 
         write (*, *) '                 University of Oxford'
         write (*, *) ' '
         write (*, *) ' '
      end if

   end subroutine write_start_message

   !****************************************************************************
   !                                Print header                                
   !****************************************************************************
   subroutine print_header

      use mp, only: proc0
      use debug_flags, only: print_extra_info_to_terminal
      use parameters_physics, only: include_apar, include_bpar
      
      implicit none
      
      !----------------------------------------------------------------------
      
      ! Only print the header on the first processor
      if (.not. proc0) return 

      ! Note that the actual data is written in <diagnostics_potential.f90>
      if (print_extra_info_to_terminal) then
         write (*, '(A)') "############################################################"
         write (*, '(A)') "                OVERVIEW OF THE SIMULATION"
         write (*, '(A)') "############################################################"
      end if
      if (include_apar .and. include_bpar) then
         write (*, '(A)') " "
         write (*, '(A)') "    istep       time          dt          |phi|^2       |apar|^2      |bpar|^2"
         write (*, '(A)') "--------------------------------------------------------------------------------"
      else if (include_apar) then
         write (*, '(A)') " "
         write (*, '(A)') "    istep       time          dt          |phi|^2       |apar|^2"
         write (*, '(A)') "--------------------------------------------------------------------"
      else if (include_bpar) then
         write (*, '(A)') " "
         write (*, '(A)') "    istep       time          dt          |phi|^2       |bpar|^2"
         write (*, '(A)') "--------------------------------------------------------------------"
      else
         write (*, '(A)') " "
         write (*, '(A)') "    istep       time          dt          |phi|^2  "
         write (*, '(A)') "------------------------------------------------------"
      end if

   end subroutine print_header
   
   !****************************************************************************
   !                   GET ARGUMENTS FROM COMMAND PROMPT LINE                   
   !****************************************************************************
   ! Parse some basic command line arguments. Currently just 'version' and 'help'.
   ! This should be called before anything else, but especially before initialising MPI.
   !****************************************************************************
   subroutine parse_command_line
   
      use git_version, only: get_git_version
      integer :: arg_count, arg_n
      integer :: arg_length
      character(len=:), allocatable :: argument
      character(len=*), parameter :: endl = new_line('a')
      
      !----------------------------------------------------------------------

      arg_count = command_argument_count()

      do arg_n = 0, arg_count
         call get_command_argument(1, length=arg_length)
         if (allocated(argument)) deallocate (argument)
         allocate (character(len=arg_length)::argument)
         call get_command_argument(1, argument)

         if ((argument == "--version") .or. (argument == "-v")) then
            write (*, '("stella version ", a)') get_git_version()
            stop
         else if ((argument == "--help") .or. (argument == "-h")) then
            write (*, '(a)') "stella [--version|-v] [--help|-h] [input file]"//endl//endl// &
               "stella is a flux tube gyrokinetic code for micro-stability and turbulence "// &
               "simulations of strongly magnetised plasma"//endl// &
               "For more help, see the documentation at https://stellagk.github.io/stella/"//endl// &
               "or create an issue https://github.com/stellaGK/stella/issues/new"//endl// &
               endl// &
               "  -h, --help     Print this message"//endl// &
               "  -v, --version  Print the stella version"
            stop
         end if
      end do
      
   end subroutine parse_command_line

!###############################################################################
!################################ FINISH STELLA ################################
!###############################################################################
! Finish a simulation, call the finialisation routines of all modules
!###############################################################################

   subroutine finish_stella(istep, last_call)

      ! Parallelisation
      use mp, only: finish_mp
      use mp, only: proc0
      
      ! Time routines
      use job_manage, only: time_message
      use timers, only: time_total
      use timers, only: time_init
      use timers, only: time_gke
      use timers, only: time_mirror
      use timers, only: time_sources
      use timers, only: time_multibox
      use timers, only: time_collisions
      use timers, only: time_field_solve
      use timers, only: time_parallel_nl
      use timers, only: time_implicit_advance
      use timers, only: time_all_diagnostics
      use timers, only: time_parallel_streaming
      use timers, only: time_individual_diagnostics
      use timers, only: time_response_matrix
      use timers, only: time_lu_decomposition
      
      ! Flags
      use file_utils, only: runtype_option_switch
      use file_utils, only: runtype_multibox
      use parallelisation_layouts, only: fields_kxkyz
      use parameters_physics, only: include_parallel_nonlinearity
      use parameters_physics, only: radial_variation
      use parameters_numerical, only: stream_implicit
      use parameters_numerical, only: drifts_implicit
      use dissipation_and_collisions, only: include_collisions
      use debug_flags, only: print_extra_info_to_terminal
      use gk_sources, only: source_option_switch
      use gk_sources, only: source_option_none
      
      ! Finish modules
      use grids_z, only: finish_z_grid
      use geometry, only: finish_geometry
      use multibox, only: finish_multibox
      use gk_sources, only: finish_sources
      use grids_species, only: finish_species
      use file_utils, only: finish_file_utils
      use grids_kxky, only: finish_grids_kxky
      use diagnostics, only: finish_diagnostics
      use grids_velocity, only: finish_velocity_grids
      use field_equations, only: finish_field_equations
      use response_matrix, only: finish_response_matrix
      use grids_extended_zgrid, only: finish_extended_zgrid
      use initialise_arrays, only: finish_arrays_vperp_kperp
      use initialise_redistribute, only: finish_redistribute
      use parameters_physics, only: finish_read_parameters_physics
      use arrays_gyro_averages, only: finish_arrays_bessel_functions
      use calculations_volume_averages, only: finish_volume_averages
      use parameters_numerical, only: finish_read_parameters_numerical
      use initialise_distribution_function, only: finish_distribution_function
      use gyrokinetic_equation_initialisation, only: finish_gyrokinetic_equation

      implicit none

      ! Arguments
      logical, intent(in), optional :: last_call
      integer, intent(in) :: istep
      
      ! Local variables
      real :: sum_timings, sum_redistribute, sum_response_matrix
      logical :: seconds = .true.
      logical :: hours = .false.
      
      !----------------------------------------------------------------------

      ! Make a clean exit of stella
      if (debug) write (*, *) 'stella::finish_stella::finish_diagnostics'
      call finish_diagnostics(istep)
      if (debug) write (*, *) 'stella::finish_stella::finish_response_matrix'
      call finish_response_matrix
      if (debug) write (*, *) 'stella::finish_stella::finish_field_equations'
      call finish_field_equations
      if (debug) write (*, *) 'stella::finish_stella::finish_gyrokinetic_equation'
      call finish_gyrokinetic_equation
      if (debug) write (*, *) 'stella::finish_stella::finish_sources'
      call finish_sources
      if (debug) write (*, *) 'stella::finish_stella::finish_volume_averages'
      call finish_volume_averages
      if (debug) write (*, *) 'stella::finish_stella::finish_extended_zgrid'
      call finish_extended_zgrid
      if (debug) write (*, *) 'stella::finish_stella::finish_multibox'
      call finish_multibox
      if (debug) write (*, *) 'stella::finish_stella::finish_arrays_vperp_kperp'
      call finish_arrays_vperp_kperp
      if (debug) write (*, *) 'stella::finish_stella::finish_arrays_bessel_functions'
      call finish_arrays_bessel_functions
      if (debug) write (*, *) 'stella::finish_stella::finish_redistribute'
      call finish_redistribute
      if (debug) write (*, *) 'stella::finish_stella::finish_distribution_function'
      call finish_distribution_function
      if (debug) write (*, *) 'stella::finish_stella::finish_velocity_grids'
      call finish_velocity_grids
      if (debug) write (*, *) 'stella::finish_stella::finish_grids_kxky'
      call finish_grids_kxky
      if (debug) write (*, *) 'stella::finish_stella::finish_read_parameters_numerical'
      call finish_read_parameters_numerical
      if (debug) write (*, *) 'stella::finish_stella::finish_species'
      call finish_species
      if (debug) write (*, *) 'stella::finish_stella::finish_parameters_physics'
      call finish_read_parameters_physics
      if (debug) write (*, *) 'stella::finish_stella::finish_geometry'
      call finish_geometry
      if (debug) write (*, *) 'stella::finish_stella::finish_z_grid'
      call finish_z_grid
      if (debug) write (*, *) 'stella::finish_stella::finish_file_utils'
      
      ! Print timings
      if (proc0 .and. print_extra_info_to_terminal) then
         
         ! Finish the files
         call finish_file_utils
         
         ! Finish the total timer
         call time_message(.false., time_total, ' Total')
         
         ! If the total time is more than 5 minutes, print the timers in minutes
         ! If the total time is more than 5 hours, print the timers in hours
         if (time_total(1)/60 > 5) seconds = .false.
         if (time_total(1)/60/60 > 5) hours = .true.
         
         ! Write the time message
         write (*, *)
         write (*, '(A)') "############################################################"
         write (*, '(A)') "                        ELAPSED TIME"
         write (*, '(A)') "############################################################"
         
         ! Redstribute
         write (*, fmt='(A)') ' '
         write (*, fmt='(A)') '                          REDISTRIBUTE'
         write (*, fmt='(A)') '                          -----------'
         sum_timings = time_mirror(1, 2) + time_field_solve(1, 2) + time_parallel_nl(1, 2)
         sum_timings = sum_timings + time_collisions(1, 2) + time_sources(1, 2)
         call write_time_message('mirror:', time_mirror(1, 2), sum_timings)
         if (fields_kxkyz) then
            call write_time_message('fields:', time_field_solve(1, 2), sum_timings)
         end if
         if (include_parallel_nonlinearity) then
            call write_time_message('parallel nonlin:', time_parallel_nl(1, 2), sum_timings)
         end if
         if (include_collisions) then
            call write_time_message('collisions:', time_collisions(1, 2), sum_timings)
         end if
         if (source_option_switch /= source_option_none) then
            call write_time_message('sources:', time_sources(1, 2), sum_timings)
         end if
         write (*, *)
         
         ! Response matrix
         write (*, fmt='(A)') ' '
         write (*, fmt='(A)') '                         RESPONSE MATRIX'
         write (*, fmt='(A)') '                         ---------------'
         sum_timings = time_response_matrix(1) + time_lu_decomposition(1)
         call write_time_message('response matrix:', time_response_matrix(1), sum_timings)
         call write_time_message('LU decomposition:', time_lu_decomposition(1), sum_timings)
         call write_time_message('total:', sum_timings, time_total(1))
         write (*, *)
         
         ! Collisions, sources and radial variation
         if (include_collisions .or. (source_option_switch/=source_option_none) .or. (runtype_option_switch==runtype_multibox)) then
            write (*, fmt='(A)') ' '
            write (*, fmt='(A)') '                             OTHER'
            write (*, fmt='(A)') '                             -----'
            sum_timings = time_collisions(1, 1) + time_sources(1, 1)
            sum_timings = sum_timings + time_multibox(1, 1) + time_multibox(1, 2)
            if (include_collisions) then
               call write_time_message('collisions:', time_collisions(1, 1), time_total(1))
            end if
            if (source_option_switch /= source_option_none) then
               call write_time_message('sources:', time_sources(1, 1), time_total(1))
            end if
            if (runtype_option_switch == runtype_multibox) then
               call write_time_message('multibox comm:', time_multibox(1, 1), time_total(1))
               call write_time_message('multibox krook:', time_multibox(1, 2), time_total(1))
            end if
            write (*, *)
         end if
            
         ! Field equations
         write (*, fmt='(A)') ' '
         write (*, fmt='(A)') '                             FIELDS'
         write (*, fmt='(A)') '                             ------' 
         sum_timings = time_field_solve(1, 3) + time_field_solve(1, 4) + time_field_solve(1, 5)
         call write_time_message('int_dv_g:', time_field_solve(1, 3), sum_timings)
         call write_time_message('calculate_phi:', time_field_solve(1, 4), sum_timings)
         call write_time_message('phi_adia_elec:', time_field_solve(1, 5), sum_timings)
         call write_time_message('total:', time_field_solve(1, 1), time_total(1))
         write (*, *)
         
         ! Gyrokinetic equation
         write (*, fmt='(A)') ' '
         write (*, fmt='(A)') '                      GYROKINETIC EQUATION'
         write (*, fmt='(A)') '                      ---------------------' 
         sum_timings = time_mirror(1, 1) + time_implicit_advance(1, 1) + time_parallel_streaming(1, 1)
         sum_timings = sum_timings + time_gke(1, 4) + time_gke(1, 5) + time_gke(1, 6) + time_gke(1, 7) + time_gke(1, 10)
         sum_timings = sum_timings + time_parallel_nl(1, 1)
         call write_time_message('mirror:', time_mirror(1, 1), sum_timings)
         call write_time_message('ExB nonlin:', time_gke(1, 7), sum_timings)
         if (stream_implicit) then
            call write_time_message('stream implicit:', time_implicit_advance(1, 1), sum_timings)
         else
            call write_time_message('stream explicit:', time_parallel_streaming(1, 1), sum_timings)
         end if
         if (.not. drifts_implicit) then
            call write_time_message('magnetic drift x:', time_gke(1, 5), sum_timings)
            call write_time_message('magnetic drift y:', time_gke(1, 4), sum_timings)
            call write_time_message('diamagnetic drift:', time_gke(1, 6), sum_timings)
         end if
         if (include_parallel_nonlinearity) then
            call write_time_message('parallel nonlin:', time_parallel_nl(1, 1), sum_timings)
         end if
         if (radial_variation) then 
            call write_time_message('radial var:', time_gke(1, 10), sum_timings)
         end if 
         call write_time_message('total:', sum_timings, time_total(1))
         write (*, *)
         
         ! Diagnostics
         write (*, fmt='(A)') ' '
         write (*, fmt='(A)') '                          DIAGNOSTICS'
         write (*, fmt='(A)') '                          -----------'
         sum_timings = time_individual_diagnostics(1, 1) + time_individual_diagnostics(1, 2) + time_individual_diagnostics(1, 3)
         sum_timings = sum_timings + time_individual_diagnostics(1, 4) + time_individual_diagnostics(1, 5) + time_individual_diagnostics(1, 6)
         call write_time_message('omega:', time_individual_diagnostics(1, 1), sum_timings)
         call write_time_message('phi:', time_individual_diagnostics(1, 2), sum_timings)
         call write_time_message('omega:', time_individual_diagnostics(1, 3), sum_timings)
         call write_time_message('fluxes:', time_individual_diagnostics(1, 4), sum_timings)
         call write_time_message('moments:', time_individual_diagnostics(1, 5), sum_timings)
         call write_time_message('distribution:', time_individual_diagnostics(1, 6), sum_timings)
         call write_time_message('total:', time_all_diagnostics(1), time_total(1))
         write (*, *)
        
         ! Gyrokinetic equations
         write (*, fmt='(A)') ' '
         write (*, fmt='(A)') '          GYROKINETIC EQUATION AND FIELD EQUATIONS'
         write (*, fmt='(A)') '          ----------------------------------------'
         sum_timings = time_gke(1, 9) + time_gke(1, 8) + time_field_solve(1, 1)
         call write_time_message('implicit GKE:', time_gke(1, 9), sum_timings)
         call write_time_message('explicit GKE:', time_gke(1, 8), sum_timings)
         call write_time_message('fields:', time_field_solve(1, 1), sum_timings)
         call write_time_message('total:', sum_timings, time_total(1))
         write (*, *)

         ! The entire stella code
         write (*, fmt='(A)') ' '
         write (*, fmt='(A)') '                        THE STELLA CODE'
         write (*, fmt='(A)') '                        ---------------'
         sum_timings = time_init(1) + time_all_diagnostics(1) + time_gke(1, 1)
         sum_redistribute = time_mirror(1, 2) + time_field_solve(1, 2) + time_parallel_nl(1, 2)
         sum_redistribute = sum_redistribute + time_collisions(1, 2) + time_sources(1, 2)
         sum_response_matrix = time_response_matrix(1) + time_lu_decomposition(1)
         call write_time_message('initialisation:', time_init(1), time_total(1))
         call write_time_message('diagnostics:', time_all_diagnostics(1), time_total(1))
         call write_time_message('gyrokinetic equations:', time_gke(1, 1), time_total(1))
         call write_time_message('redistribute:', sum_redistribute, time_total(1))
         call write_time_message('sum:', sum_timings, time_total(1))
         write (*, *)
        
      end if

      ! Finish (clean up) mpi message passing
      if (debug) write (*, *) 'stella::finish_stella::finish_mp'
      if (present(last_call)) then
         call finish_mp
         mpi_initialised = .false.
      end if
      
   contains
   
      !-------------------------------------------------------------------------
      subroutine write_time_message(routine, timer, timer_total)
      
         implicit none
         
         real, intent(in) :: timer
         real, intent(in) :: timer_total
         character(*), intent(in) :: routine
         character(len=25) :: fmt_min = "(a24, f9.2, a4, f9.2, a1)"
         character(len=25) :: fmt_sec = "(a24, f9.2, a4 f9.2, a1)"
         character(len=25) :: fmt_hours = "(a24, f9.2, a6, f9.2, a1)"
         
         !----------------------------------------------------------------------
         
         ! Write timers in hours, seconds or minutes
         if (hours) then
            write (*, fmt=fmt_hours) trim(routine), timer / 60 / 60., 'hours', timer/timer_total*100., '%'
         else if (seconds) then
            write (*, fmt=fmt_sec) trim(routine), timer, 'sec', timer/timer_total*100., '%'
         else
            write (*, fmt=fmt_min) trim(routine), timer / 60., 'min', timer/timer_total*100., '%'
         end if
      
      end subroutine write_time_message

   end subroutine finish_stella

end module init_stella
