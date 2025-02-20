program stella

   use redistribute, only: scatter
   use job_manage, only: time_message, checkstop, job_fork
   use job_manage, only: checktime
   use parameters_numerical, only: nstep, tend
   use parameters_numerical, only: avail_cpu_time
   use stella_time, only: update_time, code_time, code_dt, checkcodedt
   use dist_redistribute, only: kxkyz2vmu
   use time_advance, only: advance_stella
   use diagnostics, only: diagnostics_stella
   use stella_save, only: stella_save_for_restart
   use arrays_dist_fn, only: gnew, gvmu
   use file_utils, only: error_unit, flush_output_file
   use git_version, only: get_git_version, get_git_date
   use diagnostics_omega, only: checksaturation

   use debug_flags, only: debug => stella_debug
   ! Input file
   use parameters_diagnostics, only: nsave

   implicit none

   logical :: stop_stella = .false.
   logical :: mpi_initialized = .false.

   integer :: istep0, istep, ierr
   integer :: istatus
   real, dimension(2) :: time_init = 0.
   real, dimension(2) :: time_diagnose_stella = 0.
   real, dimension(2) :: time_total = 0.

   ! Stella version number and release date 
   character(len=40) :: git_commit 
   character(len=10) :: git_date

   call parse_command_line()

   ! Set git data automatically or manually
   git_commit = get_git_version()
   git_date = get_git_date()

   !> Initialize stella
   call init_stella(istep0, git_commit, git_date)

   !> Diagnose stella
   if (debug) write (*, *) 'stella::diagnostics_stella'
   if (istep0 == 0) call diagnostics_stella(istep0)

   !> Advance stella until istep=nstep
   if (debug) write (*, *) 'stella::advance_stella'
   istep = istep0 + 1
   do while ((code_time <= tend .AND. tend > 0) .OR. (istep <= nstep .AND. nstep > 0))
      if (debug) write (*, *) 'istep = ', istep
      if (mod(istep, 10) == 0) then
         call checkstop(stop_stella)
         call checktime(avail_cpu_time, stop_stella)
         call checkcodedt(stop_stella)
         call checksaturation(istep, stop_stella)
      end if
      if (stop_stella) exit
      call advance_stella(istep, stop_stella)
      if (stop_stella) exit
      call update_time
      if (nsave > 0 .and. mod(istep, nsave) == 0) then
         call scatter(kxkyz2vmu, gnew, gvmu)
         call stella_save_for_restart(gvmu, istep, code_time, code_dt, istatus)
      end if
      call time_message(.false., time_diagnose_stella, ' diagnostics') 
      call diagnostics_stella(istep) 
      call time_message(.false., time_diagnose_stella, ' diagnostics')
      ierr = error_unit()
      call flush_output_file(ierr)
      istep = istep + 1
   end do

   !> Finish stella
   if (debug) write (*, *) 'stella::finish_stella'
   call finish_stella(last_call=.true.)

contains

   !> Initialise stella
   !>
   !> Calls the initialisation routines for all the geometry, physics, and
   !> diagnostic modules
   subroutine init_stella(istep0, git_commit, git_date)

      use mp, only: init_mp, broadcast, sum_allreduce
      use mp, only: proc0, job
      use debug_flags, only: read_debug_flags
      use file_utils, only: init_file_utils
      use file_utils, only: runtype_option_switch, runtype_multibox
      use file_utils, only: run_name, init_job_name
      use file_utils, only: flush_output_file, error_unit
      use job_manage, only: checktime, time_message
      use parameters_physics, only: read_parameters_physics
      use parameters_physics, only: radial_variation
      use parameters_numerical, only: read_parameters_numerical
      use parameters_numerical, only: avail_cpu_time, nstep, rng_seed, delt, delt_max, delt_min
      use parameters_numerical, only: stream_implicit, driftkinetic_implicit
      use parameters_numerical, only: delt_option_switch, delt_option_auto
      use parameters_numerical, only: mat_gen, mat_read
      use parameters_kxky_grids, only: read_kxky_grid_parameters
      use species, only: init_species, read_species_knobs
      use species, only: nspec
      use zgrid, only: init_zgrid
      use zgrid, only: nzgrid, ntubes
      use geometry, only: init_geometry
      use geometry, only: finish_init_geometry
      use stella_layouts, only: init_stella_layouts, init_dist_fn_layouts
      use response_matrix, only: init_response_matrix, read_response_matrix
      use init_g, only: ginit, init_init_g, phiinit, scale_to_phiinit
      use init_g, only: tstart
      use fields, only: init_fields, advance_fields, fields_updated
      use fields_radial_variation, only: get_radial_correction
      use fields, only: rescale_fields
      use stella_time, only: init_tstart, init_delt
      use diagnostics, only: init_diagnostics
      use parameters_diagnostics, only: read_diagnostics_knobs
      use arrays_fields, only: phi, apar, bpar
      use arrays_dist_fn, only: gnew
      use dist_fn, only: init_gxyz, init_dist_fn
      use dist_redistribute, only: init_redistribute
      use time_advance, only: init_time_advance
      use extended_zgrid, only: init_extended_zgrid
      use grids_kxky, only: init_grids_kxky
      use parameters_kxky_grids, only: naky, nakx, ny, nx, nalpha
      use vpamu_grids, only: init_vpamu_grids, read_vpamu_grids_parameters
      use vpamu_grids, only: nvgrid, nmu
      use stella_transforms, only: init_transforms
      use stella_save, only: init_dt
      use multibox, only: read_multibox_parameters, init_multibox
      use multibox, only: use_dirichlet_BC, apply_radial_boundary_conditions
      use multibox, only: multibox_communicate
      use ran, only: get_rnd_seed_length, init_ranf
      use dissipation, only: init_dissipation
      use sources, only: init_sources
      use volume_averages, only: init_volume_averages, volume_average

      implicit none

      !> Starting timestep: zero unless the simulation has been restarted
      integer, intent(out) :: istep0

      ! Stella version number and release date 
      character(len=40), intent(in) :: git_commit 
      character(len=10), intent(in) :: git_date

      logical :: exit, list, restarted, needs_transforms
      character(500), target :: cbuff
      integer, dimension(:), allocatable  :: seed
      integer :: i, n, ierr
      real :: delt_saved

      !> initialize mpi message passing
      if (.not. mpi_initialized) call init_mp
      mpi_initialized = .true.

      !> initialize timer
      if (debug) write (*, *) 'stella::init_stella::check_time'
      call checktime(avail_cpu_time, exit)

      if (proc0) then 
         !> initialize file i/o
         if (debug) write (*, *) 'stella::init_stella::init_file_utils'
         call init_file_utils(list)
      end if  

      call read_debug_flags
!      if(stella_debug) debug = .true. 
      debug = debug .and. proc0

      call broadcast(list)
      call broadcast(runtype_option_switch)
      if (list) call job_fork

      if (proc0) then
         call time_message(.false., time_total, ' Total')
         call time_message(.false., time_init, ' Initialization')
      end if

      if (proc0) cbuff = trim(run_name)
      call broadcast(cbuff)
      if (.not. proc0) call init_job_name(cbuff)
      
      !> read the parameters_physics namelist from the input file
      if (debug) write (6, *) "stella::init_stella::read_parameters_physics"
      call read_parameters_physics
      if (debug) write (6, *) "stella::init_stella::read_parameters_numerical"
      call read_parameters_numerical 
      !> write message to screen with useful info regarding start of simulation
      if (debug) write (*, *) 'stella::init_stella::write_start_message'
      call write_start_message(git_commit, git_date) 
      !> read the zgrid_parameters namelist from the input file and setup the z grid
      if (debug) write (6, *) "stella::init_stella::init_zgrid"
      call init_zgrid
      !> read the species_knobs namelist from the input file
      if (debug) write (6, *) "stella::init_stella::read_species_knobs"
      call read_species_knobs
      !> read the grid option from the kt_grids_knobs namelist in the input file;
      !> depending on the grid option chosen, read the corresponding kt_grids_XXXX_parameters
      !> namelist from the input file and allocate some kx and ky arrays
      if (debug) write (6, *) "stella::init_stella::read_kxky_grid_parameters"
      call read_kxky_grid_parameters
      !> read the vpamu_grids_parameters namelist from the input file
      if (debug) write (6, *) "stella::init_stella::read_vpamu_grids_parameters"
      call read_vpamu_grids_parameters
      if (debug) write (6, *) "stella::init_stella::read_multibox_parameters"
      call read_multibox_parameters
      if (debug) write (6, *) "stella::init_stella::read_diagnostics_knobs"
      call read_diagnostics_knobs
      !> setup the various data layouts for the distribution function;
      !> e.g., vmu_lo is the layout in which vpa, mu and species may be distributed
      !> amongst processors, depending on the number of phase space points and processors
      if (debug) write (6, *) "stella::init_stella::init_dist_fn_layouts"
      call init_dist_fn_layouts(nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny, nx, nalpha)
      !> needs_transforms indicates whether or not FFTs will be needed in the simulation
      call check_transforms(needs_transforms)
      !> if FFTs are needed, init_transforms sets up the various FFTW plans
      !> and allocates the necessary arrays
      if (needs_transforms) then
         if (debug) write (*, *) "stella::init_stella::init_transforms"
         call init_transforms
      end if
      !> read in the geometry option and any necessary magnetic geometry info
      !> and use it to calculate all of the required geometric coefficients
      if (debug) write (6, *) "stella::init_stella::init_geometry"
      call init_geometry(nalpha, naky)
      if (debug) write (6, *) 'stella::init_stella::init_grids_kxky'
      call init_grids_kxky
      !> read species_parameters from input file and use the info to, e.g.,
      !> determine if a modified Boltzmann response is to be used
      if (debug) write (6, *) 'stella::init_stella::init_species'
      call init_species
      !> read init_g_knobs namelist from the input file
      !> and prepare for reading in from restart file if requested
      if (debug) write (6, *) "stella::init_stella::init_init_g"
      call init_init_g
      !> read knobs namelist from the input file
      !> and the info to determine the mixture of implicit and explicit time advance
      if (debug) write (6, *) "stella::init_stella::init_run_parameters"
      call read_parameters_physics

      if (debug) write (6, *) "stella::init_stella::init_ranf"
      n = get_rnd_seed_length()
      allocate (seed(n))
      if (rng_seed < 0) then
         call init_ranf(.true., seed, job + 2)
      else
         seed = rng_seed + 37 * (/(i - 1, i=1, n)/)
         call init_ranf(.false., seed, job + 2)
      end if
      deallocate (seed)

      !> read layouts_knobs namelist from the input file,
      !> which determines the order of parallelisation within the different layouts
      if (debug) write (6, *) 'stella::init_stella::init_stella_layouts'
      call init_stella_layouts
      !> setup the (kx,ky) grids and (x,y) grids, if applicable
      if (debug) write (6, *) 'stella::init_stella::init_multibox_subcalls'
      call init_multibox_subcalls
      !> finish_init_geometry deallocates various geometric arrays that
      !> were defined locally within the geometry_miller module when using Miller geometry
      if (debug) write (6, *) 'stella::init_stella::finish_init_geometry'
      call finish_init_geometry
      !> setup the (vpa,mu) grids and associated integration weights
      if (debug) write (6, *) 'stella::init_stella::init_vpamu_grids'
      call init_vpamu_grids
      !> set up all of the logic needed to do calculations on an extended grid in z.
      !> this extended grid could be due to use of a ballooning angle so that
      !> z goes from -N*pi to N*pi, or it could be due to the coupling of different
      !> kx modes arising from the twist-and-shift boundary condition
      if (debug) write (6, *) 'stella::init_stella::init_extended_zgrid'
      call init_extended_zgrid
      !> when doing a volume average using Fourier coefficients, the
      !> ky=0 mode gets a different weighting than finite ky modes, due
      !> to the reality condition being imposed; init_volume_averages accounts for this
      if (debug) write (6, *) 'stella::init_stella::init_volume_averages'
      call init_volume_averages
      !> allocates and initialises kperp2, vperp2 and arrays needed
      !> for gyro-averaging (j0 and j1 or equivalents)
      if (debug) write (6, *) "stella::init_stella::init_dist_fn"
      call init_dist_fn
      !> sets up the mappings between different layouts, needed
      !> to redistribute data when going from one layout to another
      if (debug) write (6, *) "stella::init_stella::init_redistribute"
      call init_redistribute
      !> read dissipation namelist from the input file and print information
      !> about chosen options to stdout
      if (debug) write (6, *) 'stella::init_stella::init_dissipation'
      call init_dissipation
      if (debug) write (6, *) 'stella::init_stella::init_sources'
      call init_sources
      !> allocate and initialise time-independent arrays needed to
      !> solve the field equations; e.g., sum_s (Z_s^2 n_s / T_s)*(1-Gamma0_s)
      if (debug) write (6, *) 'stella::init_stella::init_fields'
      call init_fields
      !> initialise the distribution function in the kxkyz_lo and store in gvmu
      if (debug) write (6, *) "stella::init_stella::ginit"
      call ginit(restarted, istep0)
      !> use mapping from kxkyz_lo to vmu_lo to get a copy of g that has ky, kx and z local to each core;
      !> stored in gnew and copied to gold
      if (debug) write (6, *) "stella::init_stella::init_gxyz"
      call init_gxyz(restarted)

      !> if initializing from restart file, set the initial time step size appropriately
      if (restarted .and. delt_option_switch == delt_option_auto) then
         delt_saved = delt
         if (debug) write (6, *) "stella::init_stella::init_dt"
         call init_dt(delt_saved, istatus)
         if (istatus == 0) delt = delt_saved
      end if
      !> set the internal time step size variable code_dt from the input variable delt
      if (debug) write (6, *) "stella::init_stella::init_delt"
      call init_delt(delt, delt_max, delt_min)
      !> allocate and calculate arrays needed for the mirror, parallel streaming,
      !> magnetic drifts, gradient drive, etc. terms during time advance
      if (debug) write (6, *) 'stella::init_stella::init_time_advance'
      call init_time_advance
      if (stream_implicit .or. driftkinetic_implicit) then
         if (mat_read) then
            if (debug) write (6, *) "stella::init_stella::read_response_matrix"
            call read_response_matrix
         else
            if (debug) write (6, *) "stella::init_stella::init_response_matrix"
            call init_response_matrix
         end if
      end if

      !> get initial field from initial distribution function
      if (debug) write (6, *) 'stella::init_stella::advance_fields'
      call advance_fields(gnew, phi, apar, bpar, dist='g')

      if (radial_variation) then
         if (debug) write (6, *) 'stella::init_stella::get_radial_correction'
         call get_radial_correction(gnew, phi, dist='g')
      end if

      !> fill in the boundary regions using auxilliary simulations if using
      !> multibox, or zero it out if using Dirichlet boundary conditions
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

      !> rescale to phiinit if just beginning a new run
      if (.not. restarted .and. scale_to_phiinit) call rescale_fields(phiinit)

      !> read diagnostics_knob namelist from the input file,
      !> open ascii output files and initialise the neetcdf file with extension .out.nc
      if (debug) write (6, *) 'stella::init_stella::init_diagnostics'
      call init_diagnostics(restarted, tstart, git_commit, git_date)
      !> initialise the code_time
      if (debug) write (6, *) 'stella::init_stella::init_tstart'
      call init_tstart(tstart)

      ierr = error_unit()
      if (proc0) call flush_output_file(ierr)

      !> Add a header to the output file
      call print_header
      !> stop the timing of the initialization
      if (proc0) call time_message(.false., time_init, ' Initialization')

   end subroutine init_stella

   !> call all the multibox communication subroutines to make sure all the jobs have
   !> the appropriate information
   subroutine init_multibox_subcalls

      use mp, only: proc0, job
      use species, only: communicate_species_multibox
      use geometry, only: communicate_geo_multibox
      use calculations_kxky, only: communicate_ktgrids_multibox
      use file_utils, only: runtype_option_switch, runtype_multibox
      use parameters_physics, only: radial_variation
      use multibox, only: init_multibox, rhoL, rhoR
      use multibox, only: communicate_multibox_parameters, multibox_communicate

      implicit none

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

   !> check_transforms checks the various physics flag choices
   !> to determine if FFTs are needed for the simulation
   subroutine check_transforms(needs_transforms)

      use file_utils, only: runtype_option_switch, runtype_multibox
      use parameters_physics, only: nonlinear, include_parallel_nonlinearity
      use parameters_physics, only: radial_variation, full_flux_surface
      use parameters_physics, only: hammett_flow_shear
      use parameters_physics, only: g_exb, g_exbfac 
      
      ! Input file
      use parameters_diagnostics, only: write_radial_moments, write_radial_fluxes

      implicit none

      logical, intent(out) :: needs_transforms

      needs_transforms = .false.
      !> if ExB or parallel nonlinearity included in the simulations, need FFTs
      if (nonlinear .or. include_parallel_nonlinearity) needs_transforms = .true.
      !> if 'global' in radial or bi-normal directions, need FFTs
      if (radial_variation .or. full_flux_surface) needs_transforms = .true.
      !> if running in multibox mode, need FFTs
      if (runtype_option_switch == runtype_multibox) needs_transforms = .true.
      !> if including flow shear using anything other than wavenumber re-mapping, need FFTs
      if (abs(g_exb * g_exbfac) > epsilon(0.) .and. .not. hammett_flow_shear) &
         needs_transforms = .true.
      !> if printing out flux-surface-averaged radial fluxes or moments, need FFTs
      if (write_radial_fluxes .or. write_radial_moments) needs_transforms = .true.

   end subroutine check_transforms

   !> Write the start message to screen
   subroutine write_start_message(git_commit, git_date)
   
      use mp, only: proc0, nproc
      use parameters_numerical, only: print_extra_info_to_terminal

      implicit none

      ! Stella version number and release date 
      character(len=40), intent(in) :: git_commit 
      character(len=10), intent(in) :: git_date

      ! Strings to format data
      character(len=23) :: str
      
      ! Only print the header on the first processor
      if (.not. proc0) return 
      
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
         write (*, '(A)') "############################################################"
         write (*, '(A)') "                     PARALLEL COMPUTING"
         write (*, '(A)') "############################################################"
      end if
      if (nproc == 1) then
         write (str, '(I10, A)') nproc, " processor."
         write (*,*) ' '; write (*, '(A,A,A)') " Running on ", adjustl(trim(str))
      else
         write (str, '(I10, A)') nproc, " processors."
         write (*, '(A,A,A)') " Running on ", adjustl(trim(str))
      end if
      write (*, *)

   end subroutine write_start_message

   subroutine print_header

      use mp, only: proc0
      use parameters_numerical, only: print_extra_info_to_terminal
      use parameters_physics, only: include_apar, include_bpar
      implicit none
      
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
      else
         write (*, '(A)') " "
         write (*, '(A)') "    istep       time          dt          |phi|^2  "
         write (*, '(A)') "------------------------------------------------------"
      end if

   end subroutine print_header

   !> Parse some basic command line arguments. Currently just 'version' and 'help'.
   !>
   !> This should be called before anything else, but especially before initialising MPI.
   subroutine parse_command_line()
      use git_version, only: get_git_version
      integer :: arg_count, arg_n
      integer :: arg_length
      character(len=:), allocatable :: argument
      character(len=*), parameter :: endl = new_line('a')

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

   !> Finish a simulation, call the finialisation routines of all modules
   subroutine finish_stella(last_call)

      use mp, only: finish_mp
      use mp, only: proc0
      use file_utils, only: finish_file_utils, runtype_option_switch, runtype_multibox
      use job_manage, only: time_message
      use parameters_physics, only: finish_read_parameters_physics
      use parameters_physics, only: include_parallel_nonlinearity, radial_variation
      use parameters_numerical, only: finish_read_parameters_numerical
      use zgrid, only: finish_zgrid
      use species, only: finish_species
      use time_advance, only: time_gke, time_parallel_nl
      use time_advance, only: finish_time_advance
      use parallel_streaming, only: time_parallel_streaming
      use mirror_terms, only: time_mirror
      use dissipation, only: time_collisions, include_collisions 
      use sources, only: finish_sources, time_sources, source_option_switch, source_option_none
      use init_g, only: finish_init_g
      use dist_fn, only: finish_dist_fn
      use dist_redistribute, only: finish_redistribute
      use fields, only: finish_fields
      use arrays_fields, only: time_field_solve
      use diagnostics, only: finish_diagnostics, time_diagnostics
      use response_matrix, only: finish_response_matrix
      use geometry, only: finish_geometry
      use extended_zgrid, only: finish_extended_zgrid
      use vpamu_grids, only: finish_vpamu_grids
      use grids_kxky, only: finish_grids_kxky
      use volume_averages, only: finish_volume_averages
      use multibox, only: finish_multibox, time_multibox
      use parameters_numerical, only: stream_implicit, drifts_implicit, fields_kxkyz
      use implicit_solve, only: time_implicit_advance
      use parameters_numerical, only: print_extra_info_to_terminal
      use parameters_numerical, only: fields_kxkyz

      implicit none

      logical, intent(in), optional :: last_call
      real :: sum_timings

      if (debug) write (*, *) 'stella::finish_stella::finish_diagnostics'
      call finish_diagnostics(istep)
      if (debug) write (*, *) 'stella::finish_stella::finish_response_matrix'
      call finish_response_matrix
      if (debug) write (*, *) 'stella::finish_stella::finish_fields'
      call finish_fields
      if (debug) write (*, *) 'stella::finish_stella::finish_time_advance'
      call finish_time_advance
      if (debug) write (*, *) 'stella::finish_stella::finish_sources'
      call finish_sources
      if (debug) write (*, *) 'stella::finish_stella::finish_volume_averages'
      call finish_volume_averages
      if (debug) write (*, *) 'stella::finish_stella::finish_extended_zgrid'
      call finish_extended_zgrid
      if (debug) write (*, *) 'stella::finish_stella::finish_multibox'
      call finish_multibox
      if (debug) write (*, *) 'stella::finish_stella::finish_dist_fn'
      call finish_dist_fn
      if (debug) write (*, *) 'stella::finish_stella::finish_redistribute'
      call finish_redistribute
      if (debug) write (*, *) 'stella::finish_stella::finish_init_g'
      call finish_init_g
      if (debug) write (*, *) 'stella::finish_stella::finish_vpamu_grids'
      call finish_vpamu_grids
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
      if (debug) write (*, *) 'stella::finish_stella::finish_zgrid'
      call finish_zgrid
      if (debug) write (*, *) 'stella::finish_stella::finish_file_utils'
      if (proc0 .and. print_extra_info_to_terminal) then
         call finish_file_utils
         call time_message(.false., time_total, ' Total')
         write (*, *)
         write (*, '(A)') "############################################################"
         write (*, '(A)') "                        ELAPSED TIME"
         write (*, '(A)') "############################################################"

         sum_timings = time_mirror(1, 1) + time_mirror(1, 2) 
         sum_timings = sum_timings + time_implicit_advance(1, 1) + time_parallel_streaming(1, 1)  
         sum_timings = sum_timings + time_gke(1, 4) + time_gke(1, 5) + time_gke(1, 6) + time_gke(1, 7) + time_gke(1, 10)
         sum_timings = sum_timings + time_parallel_nl(1, 1) + time_parallel_nl(1, 2)
         write (*, fmt='(A)') ' '
         write (*, fmt='(A)') '                    GYROKINETIC EQUATION'
         write (*, fmt='(A)') '                    ---------------------' 
         write (*, fmt=101) 'mirror:', time_mirror(1, 1) / 60., 'min', time_mirror(1, 1)/sum_timings*100., '%'
         write (*, fmt=101) '(redistribute):', time_mirror(1, 2) / 60., 'min', time_mirror(1, 2)/sum_timings*100., '%'
         write (*, fmt=101) 'ExB nonlin:', time_gke(1, 7) / 60., 'min', time_gke(1, 7)/sum_timings*100., '%'
         if (stream_implicit) then
            write (*, fmt=101) 'implicit advance:', time_implicit_advance(1, 1) / 60., 'min', & 
               time_implicit_advance(1, 1)/sum_timings*100., '%' 
         else
            write (*, fmt=101) 'stream:', time_parallel_streaming(1, 1) / 60., 'min', & 
               time_parallel_streaming(1, 1)/sum_timings*100., '%'
         end if
         if (.not. drifts_implicit) then
            write (*, fmt=101) 'dgdx:', time_gke(1, 5) / 60., 'min', time_gke(1, 5)/sum_timings*100., '%'
            write (*, fmt=101) 'dgdy:', time_gke(1, 4) / 60., 'min', time_gke(1, 4)/sum_timings*100., '%'
            write (*, fmt=101) 'wstar:', time_gke(1, 6) / 60., 'min', time_gke(1, 6)/sum_timings*100., '%'
         end if
         if (include_parallel_nonlinearity) then  
            write (*, fmt=101) 'parallel nonlin:', time_parallel_nl(1, 1) / 60., 'min'
            write (*, fmt=101) '(redistribute):', time_parallel_nl(1, 2) / 60., 'min'
         end if
         if (radial_variation) then 
            write (*, fmt=101) 'radial var:', time_gke(1, 10) / 60., 'min'
         end if 
         write (*, fmt=101) 'sum:', sum_timings / 60., 'min', sum_timings/time_total(1)*100., '%'

         sum_timings = time_field_solve(1, 1) + time_field_solve(1, 2) + time_field_solve(1, 3) & 
                        + time_field_solve(1, 4) + time_field_solve(1, 5)
         write (*, fmt='(A)') ' '
         write (*, fmt='(A)') '                           FIELDS'
         write (*, fmt='(A)') '                           ------' 
         write (*, fmt=101) 'fields:', time_field_solve(1, 1) / 60., 'min', time_field_solve(1, 1)/sum_timings*100., '%'
         write (*, fmt=101) '(int_dv_g):', time_field_solve(1, 3) / 60., 'min', time_field_solve(1, 3)/sum_timings*100., '%'
         write (*, fmt=101) '(get_phi):', time_field_solve(1, 4) / 60., 'min', time_field_solve(1, 4)/sum_timings*100., '%'
         write (*, fmt=101) '(phi_adia_elec):', time_field_solve(1, 5) / 60., 'min', time_field_solve(1, 5)/sum_timings*100., '%'
         if (fields_kxkyz) then
            write (*, fmt=101) '(redistribute):', time_field_solve(1, 2) / 60., 'min', sum_timings/time_total(1)*100., '%'
         end if
         write (*, fmt=101) 'sum:', sum_timings / 60., 'min', sum_timings/time_total(1)*100., '%'

         sum_timings = time_diagnostics(1, 1) + time_diagnostics(1, 2) + time_diagnostics(1, 3)  
         sum_timings = sum_timings + time_diagnostics(1, 4) + time_diagnostics(1, 5) + time_diagnostics(1, 6)  
         write (*, fmt='(A)') ' '
         write (*, fmt='(A)') '                        DIAGNOSTICS'
         write (*, fmt='(A)') '                        -----------' 
         write (*, fmt=101) 'calculate omega:', time_diagnostics(1, 1) / 60., 'min', time_diagnostics(1, 1)/sum_timings*100., '%'
         write (*, fmt=101) 'write phi:', time_diagnostics(1, 2) / 60., 'min', time_diagnostics(1, 3)/sum_timings*100., '%'
         write (*, fmt=101) 'write omega:', time_diagnostics(1, 3) / 60., 'min', time_diagnostics(1, 2)/sum_timings*100., '%'
         write (*, fmt=101) 'write fluxes:', time_diagnostics(1, 4) / 60., 'min', time_diagnostics(1, 4)/sum_timings*100., '%'
         write (*, fmt=101) 'write moments:', time_diagnostics(1, 5) / 60., 'min', time_diagnostics(1, 5)/sum_timings*100., '%'
         write (*, fmt=101) 'write distribution:', time_diagnostics(1, 6) / 60., 'min', time_diagnostics(1, 6)/sum_timings*100., '%'
         write (*, fmt=101) 'sum:', sum_timings / 60., 'min', sum_timings/time_total(1)*100., '%'

         sum_timings = time_collisions(1, 1) +time_collisions(1, 2)
         sum_timings = sum_timings + time_sources(1, 1) + time_sources(1, 2)
         sum_timings = sum_timings + time_sources(1, 1) + time_sources(1, 2)
         sum_timings = sum_timings + time_multibox(1, 1) + time_multibox(1, 2)
         write (*, fmt='(A)') ' '
         write (*, fmt='(A)') '                           OTHER'
         write (*, fmt='(A)') '                           -----' 
         if (include_collisions) then
            write (*, fmt=101) 'collisions:', time_collisions(1, 1) / 60., 'min'
            write (*, fmt=101) '(redistribute):', time_collisions(1, 2) / 60., 'min'
         end if
         if (source_option_switch /= source_option_none) then
            write (*, fmt=101) 'sources:', time_sources(1, 1) / 60., 'min'
            write (*, fmt=101) '(redistribute):', time_sources(1, 2) / 60., 'min'
         end if
         if (runtype_option_switch == runtype_multibox) then    
            write (*, fmt=101) 'multibox comm:', time_multibox(1, 1) / 60., 'min'
            write (*, fmt=101) 'multibox krook:', time_multibox(1, 2) / 60., 'min'
         end if
         write (*, fmt=101) 'sum:', sum_timings / 60., 'min', sum_timings/time_total(1)*100., '%'
 
         sum_timings = time_init(1) + time_diagnose_stella(1) + time_gke(1, 9) + time_gke(1, 8) 
         write (*, fmt='(A)') ' '
         write (*, fmt='(A)') '                          TOTALS'
         write (*, fmt='(A)') '                          ------'
         write (*, fmt=101) 'initialization:', time_init(1) / 60., 'min', time_init(1)/time_total(1)*100., '%'
         write (*, fmt=101) 'diagnostics:', time_diagnose_stella(1) / 60., 'min', time_diagnose_stella(1)/time_total(1)*100., '%' 
         write (*, fmt=101) 'total implicit:', time_gke(1, 9) / 60., 'min', time_gke(1, 9)/time_total(1)*100., '%'
         write (*, fmt=101) 'total explicit:', time_gke(1, 8) / 60., 'min', time_gke(1, 8)/time_total(1)*100., '%'
         write (*, fmt=101) 'sum:', sum_timings / 60., 'min', sum_timings/time_total(1)*100., '%'
         write (*, fmt=101) 'total:', time_total(1) / 60., 'min', time_total(1)/time_total(1)*100., '%'
         write (*, *)

      end if
101   format(a20, f9.2, a4, f9.2, a1)

      if (debug) write (*, *) 'stella::finish_stella::finish_mp'
      ! finish (clean up) mpi message passing
      if (present(last_call)) then
         call finish_mp
         mpi_initialized = .false.
      end if

   end subroutine finish_stella

   ! subroutine test_redistribute

   !   use stella_layouts, only: kxyz_lo, vmu_lo
   !   use zgrid, only: nzgrid, ntubes
   !   use vpamu_grids, only: nvpa, nmu
   !   use kt_grids, only: ny, ikx_max
   !   use dist_redistribute, only: kxyz2vmu
   !   use redistribute, only: scatter

   !   implicit none

   !   complex, dimension (:,:,:), allocatable :: g_kxyz_lo
   !   complex, dimension (:,:,:,:,:), allocatable :: g_vmu_lo

   !   allocate (g_kxyz_lo(nvpa,nmu,kxyz_lo%llim_proc:kxyz_lo%ulim_alloc))
   !   allocate (g_vmu_lo(ny,ikx_max,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

   !   g_kxyz_lo = 1.0
   !   g_vmu_lo = 2.0

   !   call scatter (kxyz2vmu, g_vmu_lo, g_kxyz_lo)

   !   write (*,*) 'g_vmu_lo', maxval(cabs(g_vmu_lo)), minval(cabs(g_vmu_lo))
   !   write (*,*) 'g_kxyz_lo', maxval(cabs(g_kxyz_lo)), minval(cabs(g_kxyz_lo))

   !   deallocate (g_vmu_lo, g_kxyz_lo)

   ! end subroutine test_redistribute

end program stella
