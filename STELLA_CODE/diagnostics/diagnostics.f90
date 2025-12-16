!###############################################################################
!                                  DIAGNOSTICS                                  
!###############################################################################
! 
! Routines for calculating and writing various physical diagnostics.
! 
!---------------------------------- Input file ---------------------------------
! 
! The diagnostics are writen to text files after every <nwrite> time steps. Moreover, 
! they are written to the NetCDF file at every <nwrite>*<nc_mult> time steps. At
! every <nsave> time steps, the distribution function is saved to a NetCDF file, 
! which allows the user to restart a simulation in the future, if <save_for_restart>
! if toggled to True. The frequency and growth rates are calculated for linear 
! simulations at every time step, and are also averaged over <navg> time steps.
! 
!&diagnostics
!   nwrite = 50.0
!   navg = 50.0
!   nsave = -1.0
!   nc_mult = 1.0
!   save_for_restart = .false.
!   write_all = .false.
!   write_all_time_traces = .true.
!   write_all_spectra_kxkyz = .false.
!   write_all_spectra_kxky = .false.
!   write_all_velocity_space = .false.
!   write_all_potential = .false.
!   write_all_omega = .false.
!   write_all_distribution = .false.
!   write_all_fluxes = .false.
!   write_all_moments = .false.
!/
!&diagnostics_potential
!   write_all_potential_time_traces = .false.
!   write_all_potential_spectra = .false.
!   write_phi2_vs_time = .true.
!   write_apar2_vs_time = .true.
!   write_bpar2_vs_time = .true.
!   write_phi_vs_kxkyz = .false.
!   write_apar_vs_kxkyz = .false.
!   write_bpar_vs_kxkyz = .false.
!   write_phi2_vs_kxky = .false.
!   write_apar2_vs_kxky = .false.
!   write_bpar2_vs_kxky = .false.
!/
!&diagnostics_omega
!   write_omega_vs_kxky = .not. include_nonlinear
!   write_omega_avg_vs_kxky = .not. include_nonlinear
!/
!&diagnostics_distribution
!   write_g2_vs_vpamus = .false.
!   write_g2_vs_zvpas = .false.
!   write_g2_vs_zmus = .false.
!   write_g2_vs_kxkyzs = .false.
!   write_g2_vs_zvpamus = .false.
!   write_distribution_g = .true.
!   write_distribution_h = .false.
!   write_distribution_f = .false.
!/
!&diagnostics_fluxes
!   flux_norm = .true.
!   write_fluxes_vs_time = .true.
!   write_radial_fluxes = radial_variation
!   write_fluxes_kxkyz = .false.
!   write_fluxes_kxky = .false.
!/
!&diagnostics_moments
!   write_moments = .false.
!   write_radial_moments = .false.
!/
! 
!---------------------------------- Time traces --------------------------------
! 
! By default, the following time traces are always calculated:
!  - |phi|^2(t)
!  - |apar|^2(t)
!  - |bpar|^2(t)
! 
! Moreover, for linear simulations the frequency and growth rate are calculated by default:
!  - omega(t)
!  - gamma(t)
! 
!------------------------------------ Fluxes -----------------------------------
! 
! The following fluxes can be calculated within stella:
!   - Particle flux (pflux)
!   - Heat flux (qflux)
!   - Momentum flux (vflux)
! 
! If <write_fluxes_vs_time> = .true. the following time traces are calculated:
!   - pflux(t, s)  -->  pflux_vs_s in the NetCDF file, and written to *.fluxes
!   - qflux(t, s)  -->  qflux_vs_s in the NetCDF file, and written to *.fluxes
!   - vflux(t, s)  -->  vflux_vs_s in the NetCDF file, and written to *.fluxes
! 
! If <write_fluxes_kxky> = .true. the flux spectra are calculated, note that these
! do not represent the Fourier components of the fluxes, but rather the contribution 
! of each (kx,ky) mode to the flux. To obtain the total flux, simply sum over all 
! contributions, since the reality condition has already been taking into account.
!   - pflux(t, kx, ky, s)  -->  pflux_vs_kxkys in the NetCDF file
!   - qflux(t, kx, ky, s)  -->  qflux_vs_kxkys in the NetCDF file
!   - vflux(t, kx, ky, s)  -->  vflux_vs_kxkys in the NetCDF file
!   
! If <write_fluxes_kxkyz> = .true. the following quantities are calculated:
!   - pflux(t, kx, ky, z, s)  -->  pflux_vs_kxkyzs in the NetCDF file
!   - qflux(t, kx, ky, z, s)  -->  qflux_vs_kxkyzs in the NetCDF file
!   - vflux(t, kx, ky, z, s)  -->  vflux_vs_kxkyzs in the NetCDF file
!   
! If <write_radial_fluxes> = .true. the following quantities are calculated:
!   - pflux(t, kx, s)  -->  pflux_x in the NetCDF file
!   - qflux(t, kx, s)  -->  qflux_x in the NetCDF file
!   - vflux(t, kx, s)  -->  vflux_x in the NetCDF file
! 
!###############################################################################
module diagnostics

   ! Debug Flags
   use debug_flags, only: debug => diagnostics_debug

   implicit none

   ! Make routines available to other modules
   public :: diagnose_distribution_function_and_fields
   public :: init_diagnostics
   public :: finish_diagnostics

   private

   ! Current maximum index of the time dimension in the netCDF file
   integer :: nout = 1

   ! Only initialise once
   logical :: initialised_diagnostics = .false.

contains

!###############################################################################
!############################ WRITE DIAGNOSTICS ################################
!###############################################################################

   ! Calculate and write diagnostics.
   subroutine diagnose_distribution_function_and_fields(istep)
   
      ! Parallelisation
      use mp, only: proc0
      use job_manage, only: time_message
      use timers, only: time_all_diagnostics
      use timers, only: time_individual_diagnostics

      ! Fields and distribution function
      use arrays_fields, only: phi, apar, bpar
      use arrays_distribution_function, only: gnew
      
      ! If required, we will update the fields
      use field_equations, only: advance_fields
      use field_equations, only: fields_updated

      ! Flags
      use parameters_physics, only: radial_variation
      use parameters_diagnostics, only: nc_mult, nwrite

      ! Write data
      use diagnostics_omega, only: write_omega_to_netcdf_file, calculate_omega
      use diagnostics_potential, only: write_potential_to_netcdf_file 
      use diagnostics_fluxes, only: write_fluxes_to_netcdf_file
      use diagnostics_moments, only: write_moments_to_netcdf_file
      use diagnostics_distribution, only: write_distribution_to_netcdf_file
      use write_diagnostics_to_netcdf, only: sync_nc

      implicit none

      ! The current time step
      integer, intent(in) :: istep

      ! Writing variables
      logical :: write_to_ascii_files, write_to_netcdf_file

      !-------------------------------------------------------------------------
      
      ! Start the timer
      call time_message(.false., time_all_diagnostics, ' diagnostics')

      ! We only write data at every <nwrite> or every <nwrite>*<nc_mult> time steps
      write_to_ascii_files = (mod(istep, nwrite) == 0)
      write_to_netcdf_file = (mod(istep, nwrite * nc_mult) == 0)
      
      !**********************************************************************
      !                 RUNNING AVERAGES AT EVERY TIME STEP                 !
      !**********************************************************************

      ! Calculate Omega from <phi> = exp(-i*<0mega>*t) at every time step
      call calculate_omega(istep, time_individual_diagnostics(:, 1))    
      
      !**********************************************************************
      !                 WRITE TO ASCII FILES EVERY <NWRITE>                 !
      !**********************************************************************

      ! 0nly write data to the ascii and netcdf files every <nwrite> time steps
      if (.not. write_to_ascii_files) then
         call time_message(.false., time_all_diagnostics, ' diagnostics')
         return
      end if

      ! Get the updated fields <phi>(ky,kx,z,tube) corresponding to <gnew>(ky,kx,z,tube,i[vpa,mu,s])
      if (radial_variation) fields_updated = .false.
      call advance_fields(gnew, phi, apar, bpar, dist='g')

      ! First write data that also has ascii files (do potential first since it will update the fields)
      call write_potential_to_netcdf_file(istep, nout, time_individual_diagnostics(:, 2), write_to_netcdf_file)
      call write_omega_to_netcdf_file(istep, nout, time_individual_diagnostics(:, 3), write_to_netcdf_file)  
      call write_fluxes_to_netcdf_file(nout, time_individual_diagnostics(:, 4), write_to_netcdf_file) 
      
      !**********************************************************************
      !             WRITE TO NETCDF FILES EVERY <NWRITE*NC_MULT>            !
      !**********************************************************************

      ! The ascii files are finished, the netcdf files are written every <nwrite*nc_mult> time steps
      if (.not. write_to_netcdf_file) then
         call time_message(.false., time_all_diagnostics, ' diagnostics')
         return
      end if
 
      ! Write data to the netcdf files
      call write_moments_to_netcdf_file(nout, time_individual_diagnostics(:, 5))
      call write_distribution_to_netcdf_file(nout, time_individual_diagnostics(:, 6))

      ! Synchronize the disk copy of a netCDF dataset with in-memory buffers
      if (proc0) call sync_nc

      ! Keep track of the netcdf pointer
      nout = nout + 1
      
      ! End the timer
      call time_message(.false., time_all_diagnostics, ' diagnostics')

   end subroutine diagnose_distribution_function_and_fields


!###############################################################################
!############################ INITALIZE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !======================== INITALIZE THE DIAGNOSTICS =========================
   !============================================================================
   ! Initialize the <diagnostics> module. Make sure that the other modules
   ! are initialised (zgrid, kt_grids, ...). Open/append the netcdf file with
   ! extension '.out.nc'. Open/append the ascii files ('.out'; '.fluxes'; '.omega').
   ! Gets called in the <init_stella> subroutine in the <stella> module. 
   !============================================================================
   subroutine init_diagnostics(restart, tstart, istep0)

      ! Parallelisation
      use mp, only: broadcast
      use mp, only: proc0

      ! Make sure all grids have been initialised
      use parameters_physics, only: read_parameters_physics
      use parameters_numerical, only: read_parameters_numerical
      use grids_z, only: init_z_grid
      use grids_z, only: read_parameters_z_grid
      use grids_kxky, only: init_grids_kxky
      use grids_kxky, only: read_parameters_kxky_grids
      use grids_species, only: init_species
      use grids_species, only: read_parameters_species
      use initialise_distribution_function, only: read_parameters_distribution_function
      use initialise_distribution_function, only: init_distribution_function
      use initialise_arrays, only: init_arrays_vperp_kperp
      use arrays_gyro_averages, only: init_arrays_bessel_functions
      use parameters_diagnostics, only: read_parameters_diagnostics
      use diagnostics_omega, only: init_diagnostics_omega
      use diagnostics_fluxes, only: init_diagnostics_fluxes
      use diagnostics_potential, only: init_diagnostics_potential
      use dissipation_and_collisions, only: ecoll_zeff
      use dissipation_and_collisions, only: vnew_ref
      use dissipation_and_collisions, only: zeff
      use parallelisation_layouts, only: fields_kxkyz
      
      ! Netcdf output file
      use git_version, only: get_git_version, get_git_date
      use write_diagnostics_to_netcdf, only: init_write_diagnostics_to_netcdf
      use write_diagnostics_to_netcdf, only: get_nout

      implicit none

      ! Has this simulation been restarted?
      logical, intent(in out) :: restart

      ! Current simulation time (in case of a restart)
      real, intent(in) :: tstart
      integer, intent(in out) :: istep0

      ! Stella version number and release date
      character(len=40) :: git_commit
      character(len=10) :: git_date
      
      !-------------------------------------------------------------------------

      ! Only initialise the diagnostics once
      if (initialised_diagnostics) return
      initialised_diagnostics = .true.

      ! Get git data
      git_commit = get_git_version()
      git_date = get_git_date()

      ! Should have been taken care off in the <init_stella> subroutine in the <stella> module.
      ! Nonetheless, make sure that the other routines are intialised.
      call read_parameters_physics
      call read_parameters_numerical(fields_kxkyz)
      call read_parameters_z_grid
      call read_parameters_species(vnew_ref)
      call read_parameters_kxky_grids
      call read_parameters_diagnostics
      call init_z_grid
      call init_grids_kxky
      call init_species(ecoll_zeff, zeff, vnew_ref)
      call read_parameters_distribution_function
      call init_distribution_function(restart, istep0)
      call init_arrays_vperp_kperp
      call init_arrays_bessel_functions

      ! Initialize the submodules
      call init_diagnostics_omega(restart)
      call init_diagnostics_fluxes(restart)
      call init_diagnostics_potential(restart)

      ! Open the netcdf file with extension '.out.nc'
      call init_write_diagnostics_to_netcdf(restart, git_commit, git_date)

      ! Get the final position <nout> of the time axis in the netcdf file
      if (proc0) call get_nout(tstart, nout)
      call broadcast(nout)

   end subroutine init_diagnostics

   !============================================================================
   !========================= FINALIZE THE DIAGNOSTICS =========================
   !============================================================================ 
   subroutine finish_diagnostics(istep)
      
      ! Grids 
      use grids_time, only: code_dt
      use grids_time, only: code_time
      
      ! Finish modules
      use write_diagnostics_to_netcdf, only: finish_write_diagnostics_to_netcdf
      use diagnostics_omega, only: finish_diagnostics_omega
      use diagnostics_fluxes, only: finish_diagnostics_fluxes
      use diagnostics_potential, only: finish_diagnostics_potential
      
      ! Save the stella data so the simulation can be restarted
      use save_stella_for_restart, only: save_stella_data_for_restart
      use parameters_diagnostics, only: save_for_restart

      implicit none

      integer :: istatus
      integer, intent(in) :: istep
      
      !-------------------------------------------------------------------------

      ! Save stella to be restarted if stella is being exited cleanly (hence exit_in = .true.)
      if (save_for_restart) then
         call save_stella_data_for_restart(istep, code_time, code_dt, istatus, .true.)
      end if

      ! Close the netcdf file
      call finish_write_diagnostics_to_netcdf

      ! Finish submodules
      call finish_diagnostics_omega
      call finish_diagnostics_fluxes
      call finish_diagnostics_potential

      nout = 1
      initialised_diagnostics = .false.

   end subroutine finish_diagnostics

end module diagnostics
