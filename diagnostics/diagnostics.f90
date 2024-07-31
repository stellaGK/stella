! Routines for calculating and writing various physical diagnostics.
module diagnostics

   implicit none

   public :: diagnostics_stella, init_diagnostics, finish_diagnostics 
   public :: time_diagnostics

   private

   ! Current maximum index of the time dimension in the netCDF file
   integer :: nout = 1

   ! Has this module been initialised?
   logical :: diagnostics_initialized = .false.

   ! Needed for timing various pieces of the diagnostics
   real, dimension(2, 6) :: time_diagnostics = 0.

contains

!###############################################################################
!############################ WRITE DIAGNOSTICS ################################
!###############################################################################

   ! Calculate and write diagnostics.
   subroutine diagnostics_stella(istep)

      ! Data 
      use fields_arrays, only: phi, apar, bpar
      use dist_fn_arrays, only: gnew 
      use fields, only: advance_fields 
      use constants, only: zi 

      ! Flags  
      use physics_flags, only: radial_variation
      use fields, only: fields_updated    
      use parameters_diagnostics, only: nc_mult, nwrite, debug

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
      write_to_ascii_files = (mod(istep, nwrite) == 0)
      write_to_netcdf_file = (mod(istep, nwrite * nc_mult) == 0)
      
      !**********************************************************************
      !                 RUNNING AVERAGES AT EVERY TIME STEP                 !
      !**********************************************************************

      ! Calculate Omega from <phi> = exp(-i*<0mega>*t) at every time step
      call calculate_omega(istep, time_diagnostics(:, 1))    
      
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
   ! are initialized (zgrid, kt_grids, ...). Open/append the netcdf file with
   ! extension '.out.nc'. Open/append the ascii files ('.out'; '.fluxes'; '.omega').
   ! Gets called in the <init_stella> subroutine in the <stella> module. 
   subroutine init_diagnostics(restart, tstart, git_commit, git_date)

      use zgrid, only: init_zgrid
      use kt_grids, only: init_kt_grids
      use physics_parameters, only: init_physics_parameters
      use run_parameters, only: init_run_parameters
      use species, only: init_species
      use dist_fn, only: init_dist_fn
      use init_g, only: init_init_g
      use stella_io, only: init_stella_io, get_nout
      use diagnostics_omega, only: init_diagnostics_omega
      use diagnostics_fluxes, only: init_diagnostics_fluxes 
      use diagnostics_potential, only: init_diagnostics_potential 
      use mp, only: broadcast, proc0

      implicit none

      ! Has this simulation been restarted?
      logical, intent(in) :: restart

      ! Current simulation time (in case of a restart)
      real, intent(in) :: tstart

      ! Print git information to netcdf file
      character(len=40), intent(in) :: git_commit 
      character(len=10), intent(in) :: git_date

      ! Only initialize the diagnostics once
      if (diagnostics_initialized) return
      diagnostics_initialized = .true.

      ! Should have been taken care off in the <init_stella> subroutine in the <stella> module. 
      ! Nonetheless, make sure that the other routines are intialized.
      call init_zgrid
      call init_physics_parameters
      call init_kt_grids
      call init_run_parameters
      call init_species
      call init_init_g
      call init_dist_fn

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
      use dist_redistribute, only: kxkyz2vmu
      use dist_fn_arrays, only: gnew, gvmu
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
      diagnostics_initialized = .false.

   end subroutine finish_diagnostics

end module diagnostics
