! Routines for calculating and writing various physical diagnostics.
module diagnostics

   ! Debug Flags
   use debug_flags, only: debug => diagnostics_debug

   implicit none

   public :: diagnostics_stella, init_diagnostics, finish_diagnostics 
   public :: time_diagnostics

   private

   ! Current maximum index of the time dimension in the netCDF file
   integer :: nout = 1

   ! Has this module been initialised?
   logical :: diagnostics_initialised = .false.

   ! Needed for timing various pieces of the diagnostics
   real, dimension(2, 6) :: time_diagnostics = 0.

contains

!###############################################################################
!############################ WRITE DIAGNOSTICS ################################
!###############################################################################

   ! Calculate and write diagnostics.
   subroutine diagnostics_stella(istep)

      ! Data 
      use arrays_store_fields, only: phi, apar, bpar
      use arrays_store_distribution_fn, only: gnew 
      use fields, only: advance_fields 
      use constants, only: zi 

      ! Flags  
      use parameters_physics, only: radial_variation
      use fields, only: fields_updated    
      use parameters_diagnostics, only: nc_mult, nwrite

      ! Write data 
      use diagnostics_omega, only: write_omega_to_netcdf_file, calculate_omega
      use diagnostics_potential, only: write_potential_to_netcdf_file 
      use diagnostics_fluxes, only: write_fluxes_to_netcdf_file
      use diagnostics_moments, only: write_moments_to_netcdf_file
      use diagnostics_distribution, only: write_distribution_to_netcdf_file
      use stella_io, only: sync_nc
   
      ! Routines
      use job_manage, only: time_message
      use mp, only: proc0

      implicit none

      ! The current time step
      integer, intent(in) :: istep

      ! Writing variables
      logical :: write_to_ascii_files, write_to_netcdf_file

      !---------------------------------------------------------------------- 

      ! We only write data at every <nwrite> or every <nwrite>*<nc_mult> time steps
      if (debug) write (*, *) 'COOKIE diagnostics::diagnostics_stella::start - 1'
      if (debug) write (*, *) 'COOKIE istep', istep
      if (debug) write (*, *) 'COOKIE nwrite', nwrite
      if (debug) write (*, *) 'COOKIE nc_mult', nc_mult
      write_to_ascii_files = (mod(istep, nwrite) == 0)
      write_to_netcdf_file = (mod(istep, nwrite * nc_mult) == 0)
      if (debug) write (*, *) 'COOKIE diagnostics::diagnostics_stella::start - 2'
      
      !**********************************************************************
      !                 RUNNING AVERAGES AT EVERY TIME STEP                 !
      !**********************************************************************

      ! Calculate Omega from <phi> = exp(-i*<0mega>*t) at every time step
      if (debug) write (*, *) 'COOKIE diagnostics::diagnostics_stella::start - 3'
      call calculate_omega(istep, time_diagnostics(:, 1))    
      if (debug) write (*, *) 'COOKIE diagnostics::diagnostics_stella::start - 4'
      
      !**********************************************************************
      !                 WRITE TO ASCII FILES EVERY <NWRITE>                 !
      !**********************************************************************

      ! 0nly write data to the ascii and netcdf files every <nwrite> time steps
      if (.not. write_to_ascii_files) return 
      if (debug) write (*, *) 'diagnostics::diagnostics_stella::txt_files'

      ! Get the updated fields <phi>(ky,kx,z,tube) corresponding to <gnew>(ky,kx,z,tube,i[vpa,mu,s])
      if (radial_variation) fields_updated = .false. 
      call advance_fields(gnew, phi, apar, bpar, dist='g')

      ! First write data that also has ascii files (do potential first since it will update the fields)
      call write_potential_to_netcdf_file(istep, nout, time_diagnostics(:, 2), write_to_netcdf_file)
      call write_omega_to_netcdf_file(istep, nout, time_diagnostics(:, 3), write_to_netcdf_file)  
      call write_fluxes_to_netcdf_file(nout, time_diagnostics(:, 4), write_to_netcdf_file) 
      
      !**********************************************************************
      !             WRITE TO NETCDF FILES EVERY <NWRITE*NC_MULT>            !
      !**********************************************************************

      ! The ascii files are finished, the netcdf files are written every <nwrite*nc_mult> time steps
      if (.not. write_to_netcdf_file) return
      if (debug) write (*, *) 'diagnostics::diagnostics_stella::netcdf_files' 
 
      ! Write data to the netcdf files
      if (debug) write (*, *) 'diagnostics::diagnostics_stella::netcdf_files_moments' 
      call write_moments_to_netcdf_file(nout, time_diagnostics(:, 5))
      if (debug) write (*, *) 'diagnostics::diagnostics_stella::netcdf_files_distribution' 
      call write_distribution_to_netcdf_file(nout, time_diagnostics(:, 6))

      ! Synchronize the disk copy of a netCDF dataset with in-memory buffers    
      if (proc0) call sync_nc

      ! Keep track of the netcdf pointer  
      nout = nout + 1       

   end subroutine diagnostics_stella


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
   subroutine init_diagnostics(restart, tstart)

      ! Parallelisation
      use mp, only: broadcast, proc0

      ! Make sure all grids have been initialised
      use parameters_physics, only: read_parameters_physics
      use parameters_numerical, only: read_parameters_numerical
      use grids_z, only: init_z_grid
      use grids_z, only: read_parameters_z_grid
      use grids_kxky, only: init_grids_kxky
      use grids_kxky, only: read_parameters_kxky_grids
      use grids_species, only: init_species
      use grids_species, only: read_parameters_species
      use initialise_distribution_fn, only: read_parameters_init_distribution
      use arrays_distribution_fn, only: init_arrays_distribution_fn
      use arrays_constants, only: init_arrays_vperp_kperp
      use parameters_diagnostics, only: read_parameters_diagnostics
      use diagnostics_omega, only: init_diagnostics_omega
      use diagnostics_fluxes, only: init_diagnostics_fluxes
      use diagnostics_potential, only: init_diagnostics_potential
      use dissipation_and_collisions, only: ecoll_zeff
      
      ! Netcdf output file
      use git_version, only: get_git_version, get_git_date
      use stella_io, only: init_stella_io
      use stella_io, only: get_nout

      implicit none

      ! Has this simulation been restarted?
      logical, intent(in) :: restart

      ! Current simulation time (in case of a restart)
      real, intent(in) :: tstart

      ! Stella version number and release date
      character(len=40) :: git_commit
      character(len=10) :: git_date
      
      !-------------------------------------------------------------------------

      ! Only initialise the diagnostics once
      if (diagnostics_initialised) return
      diagnostics_initialised = .true.

      ! Get git data
      git_commit = get_git_version()
      git_date = get_git_date()

      ! Should have been taken care off in the <init_stella> subroutine in the <stella> module.
      ! Nonetheless, make sure that the other routines are intialised.
      call read_parameters_physics
      call read_parameters_numerical
      call read_parameters_z_grid
      call read_parameters_species
      call read_parameters_kxky_grids
      call read_parameters_diagnostics
      call init_z_grid
      call init_grids_kxky
      call init_species(ecoll_zeff)
      call read_parameters_init_distribution
      call init_arrays_distribution_fn
      call init_arrays_vperp_kperp

      ! Initialize the submodules
      call init_diagnostics_omega(restart)
      call init_diagnostics_fluxes(restart)
      call init_diagnostics_potential(restart)

      ! Open the netcdf file with extension '.out.nc'
      call init_stella_io(restart, git_commit, git_date)

      ! Get the final position <nout> of the time axis in the netcdf file
      if (proc0) call get_nout(tstart, nout)
      call broadcast(nout)

   end subroutine init_diagnostics

   !============================================================================
   !========================= FINALIZE THE DIAGNOSTICS =========================
   !============================================================================ 
   subroutine finish_diagnostics(istep)
 
      use redistribute, only: scatter
      use stella_io, only: finish_stella_io
      use stella_time, only: code_dt, code_time
      use stella_save, only: stella_save_for_restart
      use calculations_redistribute, only: kxkyz2vmu
      use arrays_store_distribution_fn, only: gnew, gvmu
      use diagnostics_omega, only: finish_diagnostics_omega
      use diagnostics_fluxes, only: finish_diagnostics_fluxes 
      use diagnostics_potential, only: finish_diagnostics_potential
      use parameters_diagnostics, only: save_for_restart 

      implicit none

      integer :: istatus
      integer, intent(in) :: istep

      ! Save stella to be restarted
      if (save_for_restart) then
         call scatter(kxkyz2vmu, gnew, gvmu)
         call stella_save_for_restart(gvmu, istep, code_time, code_dt, istatus, .true.)
      end if

      ! Close the netcdf file
      call finish_stella_io

      ! Finish submodules
      call finish_diagnostics_omega    
      call finish_diagnostics_fluxes    
      call finish_diagnostics_potential   

      nout = 1
      diagnostics_initialised = .false.

   end subroutine finish_diagnostics

end module diagnostics
