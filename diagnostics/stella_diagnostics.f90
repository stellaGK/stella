! Routines for calculating and writing various physical diagnostics.
module stella_diagnostics

   implicit none

   public :: init_stella_diagnostics, finish_stella_diagnostics
   public :: diagnose_stella, read_stella_diagnostics_knobs
   public :: write_radial_fluxes, write_radial_moments
   public :: nsave, time_diagnostics, debug

   private

   ! Variables used to write diagnostics
   integer :: nwrite, nsave, navg, nc_mult
   logical :: save_for_restart, autostop

   ! Write potential in <diagnose_potential>
   logical :: write_phi_vs_kxkyz
   logical :: write_phi2_vs_kxky 
   logical :: write_apar_vs_time
   logical :: write_bpar_vs_time

   ! Write omega in <diagnose_omega>
   logical :: write_omega

   ! Write fluxes in <diagnose_fluxes>
   logical :: write_fluxes_kxkyz
   logical :: write_radial_fluxes
   logical :: flux_norm 

   ! Write distribution in <diagnose_distribution>
   logical :: write_g2_vs_vpamus
   logical :: write_g2_vs_zvpas
   logical :: write_g2_vs_zmus 
   logical :: write_g2_vs_kxkyzs 
   logical :: write_g2_vs_zvpamus 
   logical :: write_distribution_g
   logical :: write_distribution_h
   logical :: write_distribution_f

   ! Write moments in <diagnose_moments>
   logical :: write_radial_moments
   logical :: write_moments  

   ! Current maximum index of the time dimension in the netCDF file
   integer :: nout = 1

   ! Has this module been initialised?
   logical :: diagnostics_initialized = .false.

   ! Needed for timing various pieces of the diagnostics
   real, dimension(2, 6) :: time_diagnostics = 0.

   ! Debugging
   logical :: debug = .false.

contains

!###############################################################################
!############################ WRITE DIAGNOSTICS ################################
!###############################################################################

   ! Calculate and write diagnostics.
   subroutine diagnose_stella(istep)

      ! Data 
      use fields_arrays, only: phi, apar
      use dist_fn_arrays, only: gnew 
      use fields, only: advance_fields 
      use constants, only: zi 

      ! Flags  
      use physics_flags, only: radial_variation
      use fields, only: fields_updated    

      ! Write data 
      use diagnose_omega, only: write_omega_to_netcdf_file, calculate_omega
      use diagnose_potential, only: write_potential_to_netcdf_file 
      use diagnose_fluxes, only: write_fluxes_to_netcdf_file
      use diagnose_moments, only: write_moments_to_netcdf_file
      use diagnose_distribution, only: write_distribution_to_netcdf_file
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

      ! Calculate Omega from <phi> = exp(-i*<0mega>*t) at every time step
      call calculate_omega(istep, time_diagnostics(:, 1))    

      ! 0nly write data to the ascii and netcdf files every <nwrite> time steps
      if (.not. write_to_ascii_files) return

      ! Get the updated fields <phi>(ky,kx,z,tube) corresponding to <gnew>(ky,kx,z,tube,i[vpa,mu,s])
      if (radial_variation) fields_updated = .false. 
      call advance_fields(gnew, phi, apar, bpar, dist='g')

      ! First write data that also has ascii files (do potential first since it will update the fields)
      call write_potential_to_netcdf_file(istep, nout, time_diagnostics(:, 2), write_to_netcdf_file)
      call write_omega_to_netcdf_file(istep, nout, time_diagnostics(:, 3), write_to_netcdf_file)  
      call write_fluxes_to_netcdf_file(nout, time_diagnostics(:, 4), write_to_netcdf_file) 

      ! The ascii files are finished, the netcdf files are written every <nwrite*nc_mult> time steps
      if (.not. write_to_netcdf_file) return
 
      ! Write data to the netcdf files
      call write_moments_to_netcdf_file(nout, time_diagnostics(:, 5))
      call write_distribution_to_netcdf_file(nout, time_diagnostics(:, 6))

      ! Synchronize the disk copy of a netCDF dataset with in-memory buffers    
      if (proc0) call sync_nc

      ! Keep track of the netcdf pointer  
      nout = nout + 1       

   end subroutine diagnose_stella


!###############################################################################
!############################ INITALIZE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !======================== INITALIZE THE DIAGNOSTICS =========================
   !============================================================================ 
   ! Initialize the <stella_diagnostics> module. Make sure that the other modules
   ! are initialized (zgrid, kt_grids, ...). Open/append the netcdf file with
   ! extension '.out.nc'. Open/append the ascii files ('.out'; '.fluxes'; '.omega').
   ! Gets called in the <init_stella> subroutine in the <stella> module. 
   subroutine init_stella_diagnostics(restart, tstart, git_branch, git_commit, git_date)

      use zgrid, only: init_zgrid
      use kt_grids, only: init_kt_grids
      use physics_parameters, only: init_physics_parameters
      use run_parameters, only: init_run_parameters
      use species, only: init_species
      use dist_fn, only: init_dist_fn
      use init_g, only: init_init_g
      use stella_io, only: init_stella_io, get_nout
      use diagnose_omega, only: init_diagnose_omega
      use diagnose_fluxes, only: init_diagnose_fluxes
      use diagnose_moments, only: init_diagnose_moments
      use diagnose_potential, only: init_diagnose_potential
      use diagnose_distribution, only: init_diagnose_distribution
      use mp, only: broadcast, proc0

      implicit none

      ! Has this simulation been restarted?
      logical, intent(in) :: restart

      ! Current simulation time (in case of a restart)
      real, intent(in) :: tstart

      ! Print git information to netcdf file
      character(len=50), intent(in) :: git_branch
      character(len=40), intent(in) :: git_commit 
      character(len=10), intent(in) :: git_date

      ! Only initialize the diagnostics once
      if (diagnostics_initialized) return
      diagnostics_initialized = .true.

      ! Only debug on the first processor
      debug = debug .and. proc0

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
      call init_diagnose_omega(write_omega, navg, autostop, restart) 
      call init_diagnose_fluxes(write_fluxes_kxkyz, write_radial_fluxes, flux_norm, restart)  
      call init_diagnose_potential(restart, write_phi_vs_kxkyz, write_phi2_vs_kxky) 
      call init_diagnose_moments(write_moments, write_radial_moments) 
      call init_diagnose_distribution(write_g2_vs_vpamus, write_g2_vs_zvpas, write_g2_vs_zmus, & 
               write_g2_vs_kxkyzs, write_g2_vs_zvpamus, write_distribution_g, write_distribution_h, write_distribution_f) 

      ! Open the netcdf file with extension '.out.nc'
      call init_stella_io(restart, git_branch, git_commit, git_date) 

      ! Get the final position <nout> of the time axis in the netcdf file
      if (proc0) call get_nout(tstart, nout)
      call broadcast(nout)

   end subroutine init_stella_diagnostics

   !============================================================================
   !========================= FINALIZE THE DIAGNOSTICS =========================
   !============================================================================ 
   subroutine finish_stella_diagnostics(istep)
 
      use redistribute, only: scatter
      use stella_io, only: finish_stella_io
      use stella_time, only: code_dt, code_time
      use stella_save, only: stella_save_for_restart
      use dist_redistribute, only: kxkyz2vmu
      use dist_fn_arrays, only: gnew, gvmu
      use diagnose_omega, only: finish_diagnose_omega
      use diagnose_fluxes, only: finish_diagnose_fluxes 
      use diagnose_potential, only: finish_diagnose_potential

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
      call finish_diagnose_omega    
      call finish_diagnose_fluxes    
      call finish_diagnose_potential   

      nout = 1
      diagnostics_initialized = .false.

   end subroutine finish_stella_diagnostics

   !============================================================================
   !====================== READ AND BROADCAST INPUT FILE =======================
   !============================================================================ 
   ! Read-in the parameters from the namelist "stella_diagnostics_knobs" in the 
   ! input file and broadcast the parameters to all the processors.
   ! Gets called in the <init_stella> subroutine in the <stella> module. 
   subroutine read_stella_diagnostics_knobs

      use mp, only: broadcast

      implicit none

      ! Read the namelist "stella_diagnostics_knobs" in the input file
      call read_parameters

      ! Broadcast the variables to all processors
      call broadcast(nwrite)
      call broadcast(navg)
      call broadcast(nsave)
      call broadcast(nc_mult)
      call broadcast(autostop) 
      call broadcast(save_for_restart)
      call broadcast(write_omega)
      call broadcast(write_phi2_vs_kxky)
      call broadcast(write_moments)
      call broadcast(write_phi_vs_kxkyz)
      call broadcast(write_g2_vs_vpamus)
      call broadcast(write_g2_vs_zvpas)
      call broadcast(write_g2_vs_zmus)
      call broadcast(write_g2_vs_kxkyzs)
      call broadcast(write_g2_vs_zvpamus)
      call broadcast(write_distribution_g)
      call broadcast(write_distribution_f)
      call broadcast(write_distribution_h)
      call broadcast(write_radial_fluxes)
      call broadcast(write_radial_moments)
      call broadcast(write_fluxes_kxkyz)
      call broadcast(write_apar_vs_time)
      call broadcast(write_bpar_vs_time)
      call broadcast(flux_norm)

   end subroutine read_stella_diagnostics_knobs

   !============================================================================
   !=================== READ DIAGNOSTICS KNOB IN INPUT FILE ====================
   !============================================================================ 
   ! Define default parameters for the <stella_diagnostics> module and overwrite them
   ! with those defined in the namelist "stella_diagnostics_knobs" in the input file.
   ! Gets called in the <read_stella_diagnostics_knobs> subroutine above. 
   subroutine read_parameters

      use mp, only: proc0
      use file_utils, only: input_unit_exist
      use physics_flags, only: radial_variation, nonlinear

      implicit none

      logical :: exist
      integer :: in_file

      ! Define the namelist "stella_diagnostics_knobs" in the input file.
      namelist /stella_diagnostics_knobs/ nwrite, navg, nsave, autostop, &
         save_for_restart, write_phi_vs_kxkyz, write_g2_vs_vpamus, write_g2_vs_zvpas, write_g2_vs_zmus, &
         write_g2_vs_kxkyzs, write_g2_vs_zvpamus, write_distribution_g, write_distribution_h, write_distribution_f, &
         write_omega, write_phi2_vs_kxky, write_moments, write_radial_fluxes, &
         write_apar_vs_time, write_bpar_vs_time, &
         write_radial_moments, write_fluxes_kxkyz, flux_norm, nc_mult

      ! Read the namelist "stella_diagnostics_knobs" only with the first processor.
      if (proc0) then

         ! Stop linear simulations when gamma is constant (careful since we won't catch jumpers!)
         ! It will check gamma over <navg> time steps
         autostop = .true.

         ! Write data to the ascii files at every <nwrite> time steps.
         ! Write data to the netcdf file every <nwrite*nc_mult> time steps.
         nwrite = 50 
         nc_mult = 1

         ! Save gvmu(vpa, nmu, i[s,kx,ky,z]) at every <nsave> time steps so that we can restart the simulation.
         ! From <gvmu> we can calculate <phi> through quasi-neutrality, so <gvmu> is all that we need to save.
         save_for_restart = .false.
         nsave = -1

         ! We calculate running averages over <navg> time points.
         navg = 50

         ! By default, do not write any data (to save memory).
         write_omega = .false.
         write_phi_vs_kxkyz = .false.
         write_apar_vs_time = .false.
         write_bpar_vs_time = .false.
         write_phi2_vs_kxky = .true.
         write_moments = .false.
         write_fluxes_kxkyz = .true.      ! Note if you turn this off the code will break
         write_g2_vs_vpamus = .false.
         write_g2_vs_zvpas = .false.
         write_g2_vs_zmus = .false.
         write_g2_vs_kxkyzs = .false.
         write_g2_vs_zvpamus = .false.
         write_distribution_g = .false.
         write_distribution_h = .false.
         write_distribution_f = .false. 

         ! Flux definition with an extra factor 1/<nabla rho> in front.
         flux_norm = .true.

         ! If <radial_variation> = True, automatically write the corresponding data
         write_radial_fluxes = radial_variation
         write_radial_moments = radial_variation

         ! Read the namelist "stella_diagnostics_knobs"in the input file and overwrite the default variables
         in_file = input_unit_exist("stella_diagnostics_knobs", exist)
         if (exist) read (unit=in_file, nml=stella_diagnostics_knobs)

         ! If <save_for_restart> = False then we need <nsave> = -1
         if (.not. save_for_restart) nsave = -1

         ! For nonlinear simulations, don't stop automatically 
         if (nonlinear) autostop = .false.

      end if

   end subroutine read_parameters 

end module stella_diagnostics
