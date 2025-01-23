!> Routines for calculating and writing various physical diagnostics
module stella_diagnostics

   implicit none

   public :: init_stella_diagnostics, finish_stella_diagnostics
   public :: diagnose_stella, read_stella_diagnostics_knobs
   public :: write_radial_fluxes, write_radial_moments
   public :: nsave

   private

   integer :: ntg_out
   integer :: nwrite, nsave, navg, nc_mult
   integer :: stdout_unit, fluxes_unit, omega_unit
   logical :: save_for_restart
   logical :: write_omega
   logical :: write_moments
   logical :: write_phi_vs_time
   logical :: write_apar_vs_time
   logical :: write_bpar_vs_time
   logical :: write_gvmus
   logical :: write_gzvs
   logical :: write_kspectra
   logical :: write_radial_fluxes
   logical :: write_radial_moments
   logical :: write_fluxes_kxkyz
   logical :: flux_norm

   !> Arrays needed for averaging in x,y,z
   real, dimension(:), allocatable :: pflux_avg, vflux_avg, qflux_avg, heat_avg
   real, dimension(:, :, :), allocatable :: pflux, vflux, qflux, exchange

   !> Needed for calculating growth rates and frequencies
   complex, dimension(:, :, :), allocatable :: omega_vs_time

   !> Current maximum index of the time dimension in the netCDF file
   integer :: nout = 1
   !> Has this module been initialised?
   logical :: diagnostics_initialized = .false.

   !> Debugging
   logical :: debug = .false.

   interface get_one_flux_vmulo
      module procedure get_one_flux_vmulo_int
      module procedure get_one_flux_vmulo_kxkyz
   end interface

contains
   !> Read-in the parameters for the [[stella_diagnostics]] module
   !>
   !> Broadcast the parameters from the namelist "stella_diagnostics_knobs"
   subroutine read_stella_diagnostics_knobs

      use mp, only: broadcast

      implicit none

      !> Read the namelist "stella_diagnostics_knobs" in the input file
      call read_parameters

      !> Broadcast the variables to all processors
      call broadcast(nwrite)
      call broadcast(navg)
      call broadcast(nsave)
      call broadcast(nc_mult)
      call broadcast(save_for_restart)
      call broadcast(write_omega)
      call broadcast(write_kspectra)
      call broadcast(write_moments)
      call broadcast(write_phi_vs_time)
      call broadcast(write_apar_vs_time)
      call broadcast(write_bpar_vs_time)
      call broadcast(write_gvmus)
      call broadcast(write_gzvs)
      call broadcast(write_radial_fluxes)
      call broadcast(write_radial_moments)
      call broadcast(write_fluxes_kxkyz)
      call broadcast(flux_norm)

   end subroutine read_stella_diagnostics_knobs

   !> Initialise the [[stella_diagnostics]] module
   !>
   !> Broadcast the parameters from the namelist "stella_diagnostics_knobs"
   !> and open/append the netcdf file and the ascii files.
   subroutine init_stella_diagnostics(restart, tstart)

      use zgrid, only: init_zgrid, nperiod, nzed
      use kt_grids, only: init_kt_grids
      use physics_parameters, only: init_physics_parameters
      use run_parameters, only: init_run_parameters
      use species, only: init_species
      use dist_fn, only: init_dist_fn
      use init_g, only: init_init_g
      use stella_io, only: init_stella_io, get_nout
      use mp, only: broadcast, proc0

      implicit none

      !> Has this simulation been restarted?
      logical, intent(in) :: restart
      !> Current simulation time
      real, intent(in) :: tstart

      ! Only initialize the diagnostics once
      if (diagnostics_initialized) return
      diagnostics_initialized = .true.

      ! Only debug on the first processor
      debug = debug .and. proc0

      !> Make sure the other routines are intialized
      call init_zgrid
      call init_physics_parameters
      call init_kt_grids
      call init_run_parameters
      call init_species
      call init_init_g
      call init_dist_fn

      !> Allocate the necessary arrays
      call allocate_arrays

      ntg_out = nzed / 2 + (nperiod - 1) * nzed

      !> Initiate the netcdf file with extension '.out.nc'
      call init_stella_io(restart)

      !> Open the '.out', '.fluxes' and '.omega' file
      if (proc0) call open_loop_ascii_files(restart)

      !> Get the final position [[nout]] of the time axis in the netcdf file
      if (proc0) call get_nout(tstart, nout)
      call broadcast(nout)

   end subroutine init_stella_diagnostics

   !> Read the diagnostic input parameters from the input file
   !>
   !> Namelist: `stella_diagnostics_knobs`
   subroutine read_parameters

      use mp, only: proc0
      use file_utils, only: input_unit_exist
      use physics_flags, only: radial_variation

      implicit none

      logical :: exist
      integer :: in_file

      namelist /stella_diagnostics_knobs/ nwrite, navg, nsave, &
         save_for_restart, write_phi_vs_time, write_apar_vs_time, &
         write_bpar_vs_time, write_gvmus, write_gzvs, &
         write_omega, write_kspectra, write_moments, write_radial_fluxes, &
         write_radial_moments, write_fluxes_kxkyz, flux_norm, nc_mult

      if (proc0) then
         nwrite = 50
         navg = 50
         nsave = -1
         save_for_restart = .false.
         write_omega = .false.
         write_phi_vs_time = .false.
         write_apar_vs_time = .false.
         write_bpar_vs_time = .false.
         write_gvmus = .false.
         write_gzvs = .false.
         write_kspectra = .false.
         write_moments = .false.
         write_radial_fluxes = radial_variation
         write_radial_moments = radial_variation
         write_fluxes_kxkyz = .false.
         nc_mult = 1
         flux_norm = .true.

         in_file = input_unit_exist("stella_diagnostics_knobs", exist)
         if (exist) read (unit=in_file, nml=stella_diagnostics_knobs)

         if (.not. save_for_restart) nsave = -1
      end if

   end subroutine read_parameters

   !> Allocate the module-level arrays
   subroutine allocate_arrays

      use species, only: nspec
      use kt_grids, only: nakx, naky

      implicit none

      if (.not. allocated(pflux)) allocate (pflux(nakx, naky, nspec)); pflux = 0.
      if (.not. allocated(qflux)) allocate (qflux(nakx, naky, nspec)); qflux = 0.
      if (.not. allocated(vflux)) allocate (vflux(nakx, naky, nspec)); vflux = 0.
      if (.not. allocated(exchange)) allocate (exchange(nakx, naky, nspec)); exchange = 0.
      if (.not. allocated(pflux_avg)) allocate (pflux_avg(nspec)); pflux_avg = 0.
      if (.not. allocated(qflux_avg)) allocate (qflux_avg(nspec)); qflux_avg = 0.
      if (.not. allocated(vflux_avg)) allocate (vflux_avg(nspec)); vflux_avg = 0.
      if (.not. allocated(heat_avg)) allocate (heat_avg(nspec)); heat_avg = 0.
      if (.not. allocated(omega_vs_time)) then
         if (write_omega) then
            allocate (omega_vs_time(navg, naky, nakx))
            omega_vs_time = 0.
         else
            allocate (omega_vs_time(1, 1, 1))
            navg = 1
         end if
      end if

   end subroutine allocate_arrays

   !> Open the '.out' and the '.fluxes' file.
   !>
   !> When running a new simulation, create a new file or replace an old file.
   !> When restarting a simulation, append the old files.
   subroutine open_loop_ascii_files(restart)

      use file_utils, only: open_output_file
      use species, only: nspec

      implicit none

      logical, intent(in) :: restart
      character(3) :: nspec_str
      character(100) :: str
      logical :: overwrite

      ! Do not overwrite, but append files, when we restart the simulation.
      overwrite = .not. restart

      ! Open the '.out' and the '.fluxes' files.
      call open_output_file(stdout_unit, '.out', overwrite)
      call open_output_file(fluxes_unit, '.fluxes', overwrite)

      ! Create the header for the .fluxes file.
      ! Every column is made up of 12 spaces, so make sure the headers
      ! are placed correctly since we have the following columns for nspec=3:
      ! #time pflx1 pflx2 pflx3 vflx1 vflx2 vflx32 qflx1 qflx2 qflx3
      if (.not. restart) then
         write (nspec_str, '(i3)') nspec * 12
         str = trim('(2a12,2a'//trim(nspec_str)//')')
         write (fluxes_unit, str) '#time', 'pflx', 'vflx', 'qflx'
      end if

      ! Open the '.omega' file and create its header.
      if (write_omega) then
         call open_output_file(omega_unit, '.omega', overwrite)
         if (.not. restart) then
            write (omega_unit, '(7a12)') '#time', 'ky', 'kx', &
               'Re[om]', 'Im[om]', 'Re[omavg]', 'Im[omavg]'
         end if
      end if

   end subroutine open_loop_ascii_files

   !> Close the text files opened by [[open_loop_ascii_files]]
   subroutine close_loop_ascii_files

      use file_utils, only: close_output_file

      implicit none

      call close_output_file(stdout_unit)
      call close_output_file(fluxes_unit)
      if (write_omega) call close_output_file(omega_unit)

   end subroutine close_loop_ascii_files

   !> Calculate and write diagnostics
   subroutine diagnose_stella(istep)

      use mp, only: proc0
      use constants, only: zi
      use redistribute, only: scatter
      use fields_arrays, only: phi, apar, bpar
      use fields_arrays, only: phi_old, phi_corr_QN, apar_old
      use fields, only: fields_updated, advance_fields
      use dist_fn_arrays, only: gvmu, gnew
      use g_tofrom_h, only: g_to_h
      use stella_io, only: write_time_nc
      use stella_io, only: write_phi2_nc, write_apar2_nc, write_bpar2_nc
      use stella_io, only: write_phi_nc, write_apar_nc, write_bpar_nc
      use stella_io, only: write_gvmus_nc
      use stella_io, only: write_gzvs_nc
      use stella_io, only: write_kspectra_nc, write_kspectra_species_nc
      use stella_io, only: write_moments_nc
      use stella_io, only: write_omega_nc
      use stella_io, only: write_fluxes_nc
      use stella_io, only: write_radial_fluxes_nc
      use stella_io, only: write_radial_moments_nc
      use stella_io, only: write_fluxes_kxkyz_nc
      use stella_io, only: sync_nc
      use stella_time, only: code_time, code_dt, cfl_dt_ExB
      use zgrid, only: nztot, nzgrid, ntubes
      use vpamu_grids, only: nmu, nvpa
      use species, only: nspec
      use kt_grids, only: naky, nakx, ikx_max, ny
      use dist_redistribute, only: kxkyz2vmu
      use physics_flags, only: radial_variation, full_flux_surface
      use physics_flags, only: include_apar, include_bpar
      use volume_averages, only: volume_average, fieldline_average
      use run_parameters, only: fphi

      !> TODO-GA: add phi2 diagnostic to different folder for ffs
      use stella_geometry, only: dl_over_b
      use kt_grids, only: swap_kxky, swap_kxky_back
      use kt_grids, only: naky_all
      use stella_transforms, only: transform_ky2y, transform_y2ky
      use volume_averages, only: flux_surface_average_ffs

      use stella_geometry, only: jacob, grad_x, geo_surf
      use kt_grids, only: dy
      use volume_averages, only: mode_fac
      use zgrid, only: delzed

      implicit none

      !> The current timestep
      integer, intent(in) :: istep
      integer :: is
      real :: phi2, apar2, bpar2
      real :: zero
      real, dimension(:, :, :), allocatable :: gvmus
      real, dimension(:, :, :, :), allocatable :: gzvs
!    real, dimension (:,:,:), allocatable :: pflx_zvpa, vflx_zvpa, qflx_zvpa
      real, dimension(:), allocatable :: part_flux, mom_flux, heat_flux
      real, dimension(:, :), allocatable :: part_flux_x, mom_flux_x, heat_flux_x
      real, dimension(:, :), allocatable :: dens_x, upar_x, temp_x
      real, dimension(:, :), allocatable :: phi2_vs_kxky, apar2_vs_kxky, bpar2_vs_kxky
      real, dimension(:, :, :), allocatable :: pflx_vs_kxky, vflx_vs_kxky, qflx_vs_kxky
      real, dimension(:, :, :, :, :), allocatable :: pflx_kxkyz, vflx_kxkyz, qflx_kxkyz
      complex, dimension(:, :, :, :, :), allocatable :: density, upar, temperature, spitzer2

      complex, dimension(:, :), allocatable :: omega_avg
      complex, dimension(:, :), allocatable :: phiavg, phioldavg, aparavg, aparoldavg
      complex, dimension(:, :, :, :), allocatable :: phi_out, apar_out, bpar_out

      !> needed when simulating a full flux surface
      complex, dimension(:, :, :, :), allocatable :: dens_ffs, upar_ffs, pres_ffs

      complex, dimension(:, :, :, :), allocatable :: phi2_kxkyz
      complex, dimension(:, :), allocatable :: phi_swap
      complex, dimension(:, :, :, :), allocatable :: phi2_y, phi2_mod
      real, dimension(:, :), allocatable :: flxfac, norm

      integer :: ikx, iz, it, iky, iy
      real :: area, area1
      ! calculation of omega requires computation of omega more
      ! frequently than every nwrite time steps
      if (write_omega .and. proc0) then
         zero = 100.*epsilon(0.)
         if (istep > 0) then
            allocate (phiavg(naky, nakx))
            allocate (phioldavg(naky, nakx))
            allocate (aparavg(naky, nakx))
            allocate (aparoldavg(naky, nakx))
            call fieldline_average(phi, phiavg)
            call fieldline_average(phi_old, phioldavg)
            if (include_apar) then 
                call fieldline_average(apar, aparavg)
                call fieldline_average(apar_old, aparoldavg)
                ! add <apar> to <phi> in the case <phi> = 0 because the mode has tearing parity
                ! if the mode is a purely growing mode then 
                ! (<phi^n+1> + <apar^n+1> )/(<phi^n> + <apar^n>) = exp (-i delta t omega)  
                phiavg = phiavg + aparavg
                phioldavg = phioldavg + aparoldavg
            end if
            where (abs(phiavg) < zero .or. abs(phioldavg) < zero)
               omega_vs_time(mod(istep, navg) + 1, :, :) = 0.0
            elsewhere
               omega_vs_time(mod(istep, navg) + 1, :, :) = log(phiavg / phioldavg) * zi / code_dt
            end where
            deallocate (phiavg, phioldavg)
            deallocate (aparavg, aparoldavg)
         end if
      end if

      !> only write data to file every nwrite time steps
      if (mod(istep, nwrite) /= 0) return

      if (radial_variation) fields_updated = .false.
      !> get the updated fields corresponding to gnew
      call advance_fields(gnew, phi, apar, bpar, dist='g')

      allocate (phi_out(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (apar_out(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (bpar_out(naky, nakx, -nzgrid:nzgrid, ntubes))
      phi_out = phi
      if (radial_variation) then
         phi_out = phi_out + phi_corr_QN
      end if
      apar_out = apar
      bpar_out = bpar
      allocate (part_flux(nspec))
      allocate (mom_flux(nspec))
      allocate (heat_flux(nspec))

      if (write_radial_fluxes) then
         allocate (part_flux_x(nakx, nspec))
         allocate (mom_flux_x(nakx, nspec))
         allocate (heat_flux_x(nakx, nspec))
      end if

      if (full_flux_surface) then
         allocate (dens_ffs(ny, ikx_max, -nzgrid:nzgrid, nspec))
         allocate (upar_ffs(ny, ikx_max, -nzgrid:nzgrid, nspec))
         allocate (pres_ffs(ny, ikx_max, -nzgrid:nzgrid, nspec))
         allocate (pflx_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         allocate (vflx_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         allocate (qflx_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      else
         allocate (pflx_kxkyz(naky, nakx, nztot, ntubes, nspec))
         allocate (vflx_kxkyz(naky, nakx, nztot, ntubes, nspec))
         allocate (qflx_kxkyz(naky, nakx, nztot, ntubes, nspec))
      end if

      !> obtain turbulent fluxes
      if (radial_variation .or. write_radial_fluxes) then
         !> handle g_to_h in get_fluxes_vmulo to eliminate x^2 terms
         !> call g_to_h (gnew, phi, fphi, phi_corr_QN)
         if (debug) write (*, *) 'stella_diagnostics::write_fluxes_vmulo'
         call get_fluxes_vmulo(gnew, phi_out, part_flux, mom_flux, heat_flux, &
                               part_flux_x, mom_flux_x, heat_flux_x, &
                               pflx_kxkyz, vflx_kxkyz, qflx_kxkyz)
!     call g_to_h (gnew, phi, -fphi, phi_corr_QN)
      else if (full_flux_surface) then
         !> calculate the particle density, parallel flow and pressure in (y,kx,z) space
         !> for all species
         if (debug) write (*, *) 'stella_diagnostics::write_fluxes_ffs'
         call get_moments_ffs(gnew, dens_ffs, upar_ffs, pres_ffs)
         !> calculate the (ky,kx) contributions to the particle, parallel momentum and energy fluxes
         call get_fluxes_ffs(dens_ffs, upar_ffs, pres_ffs, part_flux, mom_flux, heat_flux, &
                             pflx_kxkyz, vflx_kxkyz, qflx_kxkyz)
      else
         if (debug) write (*, *) 'stella_diagnostics::write_fluxes'
         !> redistribute data so that data for each vpa and mu are guaranteed to be on each processor
         call scatter(kxkyz2vmu, gnew, gvmu)
         !> get_fluxes assumes the non-Boltzmann part of the distribution, h, is passed in;
         !> convert from <delta f> = g to h
         call g_to_h(gvmu, phi, bpar, fphi)
         !> compute the fluxes
         call get_fluxes(gvmu, part_flux, mom_flux, heat_flux, &
                         pflx_kxkyz, vflx_kxkyz, qflx_kxkyz)
         !> convert back from h to g
         call g_to_h(gvmu, phi, bpar, -fphi)
      end if

      if (proc0) then
         if (write_omega) then
            allocate (omega_avg(naky, nakx))
            omega_avg = sum(omega_vs_time, dim=1) / real(navg)
         else
            allocate (omega_avg(1, 1))
         end if
         ! if (full_flux_surface) then
         !    it = 1
            
         !    allocate (phi2_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes))
         !    allocate (phi_swap(naky_all, ikx_max))
         !    allocate (phi2_y(ny, ikx_max, -nzgrid:nzgrid, ntubes))
         !    allocate (flxfac(ny, -nzgrid:nzgrid))
         !    allocate (phi2_mod(naky, nakx, -nzgrid:nzgrid, ntubes))

         !    phi2 = 0.0
         !    phi2_kxkyz = real(phi_out * conjg(phi_out))

         !    do iz = -nzgrid, nzgrid
         !       call swap_kxky(phi2_kxkyz(:, :, iz, it), phi_swap(:, :))
         !       call transform_ky2y(phi_swap(:, :), phi2_y(:, :, iz, it))
         !    end do

         !    flxfac = spread(delzed * dy, 1, ny) * jacob
         !    flxfac(:, -nzgrid) = 0.5 * flxfac(:, -nzgrid)
         !    flxfac(:, nzgrid) = 0.5 * flxfac(:, nzgrid)

         !    !! Area of flux annulus
         !    area = sum(flxfac) / ny 
         !    !! Ny * Area of fluxtube for ia = 1 
         !    area1 = sum(flxfac(1,:))

         !    flxfac = flxfac * area / area1**2

         !    call get_modified_fourier_coefficient(phi2_y, phi2_mod, flxfac)

         !    do iz = -nzgrid, nzgrid
         !       do ikx = 1, nakx
         !          do iky = 1, naky
         !             phi2 = phi2 + mode_fac(iky) * phi2_mod(iky, ikx, iz, it)
         !          end do
         !       end do
         !    end do
            
         !    apar2 = 0.0 
         !    deallocate (phi2_kxkyz, phi_swap, phi2_y, flxfac, phi2_mod)
         ! else

         call volume_average(phi_out, phi2)
         call volume_average(apar_out, apar2)
         call volume_average(bpar_out, bpar2)
         ! Print information to stella.out, the header is printed in stella.f90
         write (*, '(A2,I7,A2,ES12.4,A2,ES12.4,A2,ES12.4,A2,ES12.4)') &
            " ", istep, " ", code_time, " ", code_dt, " ", cfl_dt_ExB, " ", phi2
         call write_loop_ascii_files(istep, phi2, apar2, bpar2, part_flux, mom_flux, heat_flux, &
                                     omega_vs_time(mod(istep, navg) + 1, :, :), omega_avg)

         ! do not need omega_avg again this time step
         deallocate (omega_avg)
      end if

      if (mod(istep, nwrite * nc_mult) == 0) then
         if (proc0) then
            if (debug) write (*, *) 'stella_diagnostics::write_time_nc'
            call write_time_nc(nout, code_time)
            if (write_omega) call write_omega_nc(nout, omega_vs_time(mod(istep, navg) + 1, :, :))
            call write_phi2_nc(nout, phi2)
            call write_apar2_nc(nout, apar2)
            call write_bpar2_nc(nout, bpar2)
            call write_fluxes_nc(nout, part_flux, mom_flux, heat_flux)
            if (write_phi_vs_time) then
               if (debug) write (*, *) 'stella_diagnostics::diagnose_stella::write_phi_nc'
               call write_phi_nc(nout, phi_out)
            end if
            if (write_apar_vs_time) then
               if (debug) write (*, *) 'stella_diagnostics::diagnose_stella::write_apar_nc'
               call write_apar_nc(nout, apar_out)
            end if
            if (write_bpar_vs_time) then
               if (debug) write (*, *) 'stella_diagnostics::diagnose_stella::write_bpar_nc'
               call write_bpar_nc(nout, bpar_out)
            end if
            if (write_kspectra) then
               if (debug) write (*, *) 'stella_diagnostics::diagnose_stella::write_kspectra'
               allocate (phi2_vs_kxky(naky, nakx))
               call fieldline_average(real(phi_out * conjg(phi_out)), phi2_vs_kxky)
               call write_kspectra_nc(nout, phi2_vs_kxky, "phi2_vs_kxky", "electrostatic potential")
               deallocate (phi2_vs_kxky)
               if (include_apar) then
                  allocate (apar2_vs_kxky(naky, nakx))
                  call fieldline_average(real(apar_out * conjg(apar_out)), apar2_vs_kxky)
                  call write_kspectra_nc(nout, apar2_vs_kxky, "apar2_vs_kxky", "parallel vector potential")
                  deallocate (apar2_vs_kxky)
               end if
               if (include_bpar) then
                  allocate (bpar2_vs_kxky(naky, nakx))
                  call fieldline_average(real(bpar_out * conjg(bpar_out)), bpar2_vs_kxky)
                  call write_kspectra_nc(nout, bpar2_vs_kxky, "bpar2_vs_kxky", "parallel magnetic field fluctuation")
                  deallocate (bpar2_vs_kxky)
               end if
               !> here write out the spectrum of contributions to the fluxes, after averaging over zed
               allocate (pflx_vs_kxky(naky, nakx, nspec))
               do is = 1, nspec
                  !call fieldline_average(pflx_kxkyz(:,:,:,:,is), pflx_vs_kxky(:,:,is))
                  pflx_vs_kxky(:,:,is) = sum(sum(pflx_kxkyz(:,:,:,:,is),dim=4),dim=3)
               end do
               call write_kspectra_species_nc(nout, pflx_vs_kxky, "pflx_vs_kxky", "particle flux contributions by (kx,ky)")
               deallocate (pflx_vs_kxky)
               allocate (vflx_vs_kxky(naky, nakx, nspec))
               do is = 1, nspec
                  !call fieldline_average(vflx_kxkyz(:,:,:,:,is), vflx_vs_kxky(:,:,is))
                  vflx_vs_kxky(:,:,is) = sum(sum(vflx_kxkyz(:,:,:,:,is),dim=4),dim=3)
               end do
               call write_kspectra_species_nc(nout, vflx_vs_kxky, "vflx_vs_kxky", "momentum flux contributions by (kx,ky)")
               deallocate (vflx_vs_kxky)
               allocate (qflx_vs_kxky(naky, nakx, nspec))
               do is = 1, nspec
                  !call fieldline_average(qflx_kxkyz(:,:,:,:,is), qflx_vs_kxky(:,:,is))
                  qflx_vs_kxky(:,:,is) = sum(sum(qflx_kxkyz(:,:,:,:,is),dim=4),dim=3)
               end do
               call write_kspectra_species_nc(nout, qflx_vs_kxky, "qflx_vs_kxky", "heat flux contributions by (kx,ky)")
               deallocate (qflx_vs_kxky)
            end if
            if (write_radial_fluxes) then
               call write_radial_fluxes_nc(nout, part_flux_x, mom_flux_x, heat_flux_x)
            end if
         end if
         if (write_moments .or. write_radial_moments) then
            if (debug) write (*, *) 'stella_diagnostics::diagnose_stella::write_moments'
            allocate (density(naky, nakx, nztot, ntubes, nspec))
            allocate (upar(naky, nakx, nztot, ntubes, nspec))
            allocate (temperature(naky, nakx, nztot, ntubes, nspec))
            allocate (spitzer2(naky, nakx, nztot, ntubes, nspec))
            if (write_radial_moments) then
               allocate (dens_x(nakx, nspec))
               allocate (upar_x(nakx, nspec))
               allocate (temp_x(nakx, nspec))
            end if
            call get_moments(gnew, density, upar, temperature, dens_x, upar_x, temp_x, spitzer2)
            if (proc0 .and. write_moments) call write_moments_nc(nout, density, upar, temperature, spitzer2)
            if (proc0 .and. write_radial_moments) call write_radial_moments_nc(nout, dens_x, upar_x, temp_x)
            deallocate (density, upar, temperature, spitzer2)
            if (allocated(dens_x)) deallocate (dens_x)
            if (allocated(upar_x)) deallocate (upar_x)
            if (allocated(temp_x)) deallocate (temp_x)
         end if
         if (write_fluxes_kxkyz) then
            if (debug) write (*, *) 'stella_diagnostics::diagnose_stella::write_fluxes_kxkyz'
            if (proc0) call write_fluxes_kxkyz_nc(nout, pflx_kxkyz, vflx_kxkyz, qflx_kxkyz)
         end if
         if (write_gvmus) then
            allocate (gvmus(nvpa, nmu, nspec))
            if (debug) write (*, *) 'stella_diagnostics::diagnose_stella::get_gvmus'
            ! note that gvmus is h at this point
            call get_gvmus(gvmu, gvmus)
            if (debug) write (*, *) 'stella_diagnostics::diagnose_stella::write_gvmus_nc'
            if (proc0) call write_gvmus_nc(nout, gvmus)
            deallocate (gvmus)
         end if
         if (write_gzvs) then
            allocate (gzvs(ntubes, nztot, nvpa, nspec))
            if (debug) write (*, *) 'stella_diagnostics::diagnose_stella::get_gzvs'
            call get_gzvs(gnew, gzvs)
            if (debug) write (*, *) 'stella_diagnostics::diagnose_stella::write_gzvs_nc'
            if (proc0) call write_gzvs_nc(nout, gzvs)
            deallocate (gzvs)
         end if
         if (proc0) call sync_nc

         nout = nout + 1
      end if

      deallocate (part_flux, mom_flux, heat_flux)
      deallocate (pflx_kxkyz, vflx_kxkyz, qflx_kxkyz)
      deallocate (phi_out)
      deallocate (apar_out)
      deallocate (bpar_out)
      if (allocated(part_flux_x)) deallocate (part_flux_x)
      if (allocated(mom_flux_x)) deallocate (mom_flux_x)
      if (allocated(heat_flux_x)) deallocate (heat_flux_x)
      if (allocated(dens_ffs)) deallocate (dens_ffs)
      if (allocated(upar_ffs)) deallocate (upar_ffs)
      if (allocated(pres_ffs)) deallocate (pres_ffs)

!    nout = nout + 1

   end subroutine diagnose_stella

   !> Calculate fluxes
   !>
   !> Assumes that the non-Boltzmann part of df is passed in (aka h)
   subroutine get_fluxes(g, pflx, vflx, qflx, &
                         pflx_vs_kxkyz, vflx_vs_kxkyz, qflx_vs_kxkyz)

      use mp, only: sum_reduce
      use constants, only: zi
      use fields_arrays, only: phi, apar, bpar
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use species, only: spec, nspec
      use stella_geometry, only: jacob, grho, bmag, btor
      use stella_geometry, only: gds21, gds22
      use stella_geometry, only: geo_surf
      use zgrid, only: delzed, nzgrid, ntubes
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: vperp2, vpa, mu
      use physics_flags, only: include_apar, include_bpar
      use run_parameters, only: fphi
      use kt_grids, only: aky, theta0
      use gyro_averages, only: gyro_average, gyro_average_j1
      
      use vpamu_grids, only: maxwell_vpa, maxwell_mu   
      !> For momentum flux -- to allow for non-axisymmetry                                                                                                                                                                                                                                                                                                                                                          
      use stella_geometry, only: gds2
      use stella_geometry, only: gradzeta_grady, gradzeta_gradx, gradpar_zeta
      use kt_grids, only: akx
      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      real, dimension(:), intent(out) :: pflx, vflx, qflx
      real, dimension(:, :, -nzgrid:, :, :), intent(out) :: pflx_vs_kxkyz, vflx_vs_kxkyz, qflx_vs_kxkyz

      integer :: ikxkyz, iky, ikx, iz, it, is, ia
      real, dimension(:), allocatable :: flx_norm
      real :: flx_norm_partial
      complex, dimension(:, :), allocatable :: gtmp1, gtmp2, gtmp3

      allocate (flx_norm(-nzgrid:nzgrid))
      allocate (gtmp1(nvpa, nmu), gtmp2(nvpa, nmu), gtmp3(nvpa, nmu))

      pflx = 0.; vflx = 0.; qflx = 0.
      pflx_vs_kxkyz = 0.; vflx_vs_kxkyz = 0.; qflx_vs_kxkyz = 0.

      flx_norm = jacob(1, :) * delzed
      flx_norm(-nzgrid) = 0.5 * flx_norm(-nzgrid)
      flx_norm(nzgrid) = 0.5 * flx_norm(nzgrid)

      if (flux_norm) then
         ! Flux definition with an extra factor 1/<\nabla\rho> in front.
         flx_norm_partial = sum(flx_norm) / sum(flx_norm * grho(1, :))
         flx_norm = flx_norm / sum(flx_norm * grho(1, :))
      else
         ! Flux definition without the extra factor.
         flx_norm_partial = 1.0
         flx_norm = flx_norm / sum(flx_norm)
      end if

      ia = 1
      ! get electrostatic contributions to fluxes
      if (fphi > epsilon(0.0)) then
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! get particle flux
            call gyro_average(g(:, :, ikxkyz), ikxkyz, gtmp1)   
            !> TODO-GA: use g to f routine rather than g to h -- but check first 
            if(.not. include_apar .and. .not. include_bpar) then 
               gtmp1 = gtmp1 - spec(is)%zt * phi(iky,ikx,iz,it) * spread(maxwell_vpa(:, is),2,nmu) &
                    * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) 
            end if
            call get_one_flux(iky, iz, flx_norm(iz), gtmp1, phi(iky, ikx, iz, it), pflx(is))
            call get_one_flux(iky, iz, flx_norm_partial, gtmp1, phi(iky, ikx, iz, it), pflx_vs_kxkyz(iky, ikx, iz, it, is))

            ! get heat flux
            ! NEEDS TO BE MODIFIED TO TREAT ENERGY = ENERGY(ALPHA)
            gtmp1 = gtmp1 * (spread(vpa**2, 2, nmu) + spread(vperp2(1, iz, :), 1, nvpa))
            call get_one_flux(iky, iz, flx_norm(iz), gtmp1, phi(iky, ikx, iz, it), qflx(is))
            call get_one_flux(iky, iz, flx_norm_partial, gtmp1, phi(iky, ikx, iz, it), qflx_vs_kxkyz(iky, ikx, iz, it, is))

            !> Get momentum flux - allowing for non-axisymmetry
            !> Parallel component
            call gyro_average(g(:, :, ikxkyz), ikxkyz, gtmp1)
            gtmp2 = gtmp1 * spread(vpa, 2, nmu) * gradpar_zeta(ia, iz)
            !> Perpendicular component
            call gyro_average_j1(g(:, :, ikxkyz), ikxkyz, gtmp1)
            gtmp3 = -gtmp1 * spread(vperp2(ia, iz, :), 1, nvpa) * spec(is)%smz &
                    * zi * aky(iky) * (gradzeta_grady(ia, iz) * (gds21(ia, iz) + theta0(iky, ikx) * gds22(ia, iz)) / geo_surf%shat &
                                       - gradzeta_gradx(ia, iz) * (theta0(iky, ikx) * gds21(ia, iz) + gds2(ia, iz)))
            !> Sum parallel and perpendicular components together
            gtmp1 = gtmp2 + gtmp3

            call get_one_flux(iky, iz, flx_norm(iz), gtmp1, phi(iky, ikx, iz, it), vflx(is))
            call get_one_flux(iky, iz, flx_norm_partial, gtmp1, phi(iky, ikx, iz, it), vflx_vs_kxkyz(iky, ikx, iz, it, is))
         end do
      end if

      if (include_apar) then
         ! particle flux
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! Apar contribution to particle flux
            gtmp1 = -2.0*g(:, :, ikxkyz) * spec(is)%stm * spread(vpa, 2, nmu)
            call gyro_average(gtmp1, ikxkyz, gtmp2)
            call get_one_flux(iky, iz, flx_norm(iz), gtmp2, apar(iky, ikx, iz, it), pflx(is))
            !> TODO-GA: check use of flx_norm_partial below
            call get_one_flux(iky, iz, flx_norm_partial, gtmp2, apar(iky, ikx, iz, it), pflx_vs_kxkyz(iky, ikx, iz, it, is))

            ! Apar contribution to heat flux
            gtmp2 = gtmp2 * (spread(vpa**2, 2, nmu) + spread(vperp2(ia, iz, :), 1, nvpa))
            call get_one_flux(iky, iz, flx_norm(iz), gtmp2, apar(iky, ikx, iz, it), qflx(is))
            !> TODO-GA: check use of flx_norm_partial below
            call get_one_flux(iky, iz, flx_norm_partial, gtmp2, apar(iky, ikx, iz, it), qflx_vs_kxkyz(iky, ikx, iz, it, is))

            ! Apar contribution to momentum flux
            ! parallel component
            gtmp1 = -2.0*spread(vpa**2, 2, nmu) * spec(is)%stm * g(:, :, ikxkyz) &
                    * geo_surf%rmaj * btor(iz) / bmag(1, iz)
            call gyro_average(gtmp1, ikxkyz, gtmp2)
            ! perp component
            gtmp1 = 2.0*spread(vpa, 2, nmu) * spec(is)%stm * g(:, :, ikxkyz) &
                    * zi * aky(iky) * spread(vperp2(ia, iz, :), 1, nvpa) * geo_surf%rhoc &
                    * (gds21(ia, iz) + theta0(iky, ikx) * gds22(ia, iz)) * spec(is)%smz &
                    / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2)
            call gyro_average_j1(gtmp1, ikxkyz, gtmp3)
            ! FLAG -- NEED TO ADD IN CONTRIBUTION FROM BOLTZMANN PIECE !!

            gtmp1 = gtmp2 + gtmp3
            !> TODO-GA: check use of flx_norm_partial below
            call get_one_flux(iky, iz, flx_norm(iz), gtmp1, apar(iky, ikx, iz, it), vflx(is))
            call get_one_flux(iky, iz, flx_norm_partial, gtmp1, apar(iky, ikx, iz, it), vflx_vs_kxkyz(iky, ikx, iz, it, is))
         end do
      end if
      
      if (include_bpar) then
         ! particle flux
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! Bpar contribution to particle flux
            gtmp1 = 4.0 * g(:, :, ikxkyz) * spec(is)%tz * spread(mu, 1, nvpa)
            call gyro_average_j1(gtmp1, ikxkyz, gtmp2)
            call get_one_flux(iky, iz, flx_norm(iz), gtmp2, bpar(iky, ikx, iz, it), pflx(is))
            !> TODO-GA: check use of flx_norm_partial below
            call get_one_flux(iky, iz, flx_norm(iz), gtmp2, bpar(iky, ikx, iz, it), pflx_vs_kxkyz(iky, ikx, iz, it, is))

            ! Bpar contribution to heat flux
            gtmp2 = gtmp2 * (spread(vpa**2, 2, nmu) + spread(vperp2(ia, iz, :), 1, nvpa))
            call get_one_flux(iky, iz, flx_norm(iz), gtmp2, bpar(iky, ikx, iz, it), qflx(is))
            !> TODO-GA: check use of flx_norm_partial below
            call get_one_flux(iky, iz, flx_norm(iz), gtmp2, bpar(iky, ikx, iz, it), qflx_vs_kxkyz(iky, ikx, iz, it, is))

            ! Bpar contribution to momentum flux
            ! parallel component
            gtmp1 = 4.0 * spec(is)%tz * spread(mu, 1, nvpa) * g(:, :, ikxkyz) &
                    * spread(vpa, 2, nmu) * geo_surf%rmaj * btor(iz) / bmag(ia, iz)
            call gyro_average_j1(gtmp1, ikxkyz, gtmp2)
            ! perp component
            gtmp3 = 0.0
            ! NOT SUPPORTED, REQUIRES d J1(x)/ d x 
            gtmp1 = gtmp2 + gtmp3

            call get_one_flux(iky, iz, flx_norm(iz), gtmp1, bpar(iky, ikx, iz, it), vflx(is))
            !> TODO-GA: check use of flx_norm_partial below
            call get_one_flux(iky, iz, flx_norm(iz), gtmp1, bpar(iky, ikx, iz, it), vflx_vs_kxkyz(iky, ikx, iz, it, is))
         end do
      end if

      call sum_reduce(pflx, 0); pflx = pflx * spec%dens_psi0
      call sum_reduce(qflx, 0); qflx = qflx * spec%dens_psi0 * spec%temp_psi0
      call sum_reduce(vflx, 0); vflx = vflx * spec%dens_psi0 * sqrt(spec%mass * spec%temp_psi0)

      ! normalise to account for contributions from multiple flux tubes
      ! in flux tube train
      pflx = pflx / real(ntubes)
      qflx = qflx / real(ntubes)
      vflx = vflx / real(ntubes)

      call sum_reduce(pflx_vs_kxkyz, 0)
      call sum_reduce(qflx_vs_kxkyz, 0)
      call sum_reduce(vflx_vs_kxkyz, 0)

      do is = 1, nspec
         pflx_vs_kxkyz(:, :, :, :, is) = pflx_vs_kxkyz(:, :, :, :, is) * spec(is)%dens_psi0
         qflx_vs_kxkyz(:, :, :, :, is) = qflx_vs_kxkyz(:, :, :, :, is) * spec(is)%dens_psi0 * spec(is)%temp_psi0
         vflx_vs_kxkyz(:, :, :, :, is) = vflx_vs_kxkyz(:, :, :, :, is) * spec(is)%dens_psi0 * sqrt(spec(is)%mass * spec(is)%temp_psi0)
      end do

      deallocate (gtmp1, gtmp2, gtmp3)
      deallocate (flx_norm)

   end subroutine get_fluxes

   !==============================================
   !============ GET FLUXES VMULO ================
   !==============================================
   subroutine get_fluxes_vmulo(g, phi, pflx, vflx, qflx, &
                               pflx_x, vflx_x, qflx_x, &
                               pflx_vs_kxkyz, vflx_vs_kxkyz, qflx_vs_kxkyz)

      use mp, only: sum_reduce
      use constants, only: zi
      use dist_fn_arrays, only: g1, g2, kperp2, dkperp2dr
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use species, only: spec
      use stella_geometry, only: grho_norm, bmag, btor
      use stella_geometry, only: drhodpsi
      use stella_geometry, only: gds21, gds22
      use stella_geometry, only: dgds21dr, dgds22dr
      use stella_geometry, only: geo_surf
      use stella_geometry, only: dBdrho, dIdrho
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: vperp2, vpa, mu
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use run_parameters, only: fphi
      use run_parameters, only: maxwellian_normalization
      use kt_grids, only: aky, theta0, naky, nakx, multiply_by_rho
      use physics_flags, only: radial_variation
      use gyro_averages, only: gyro_average, gyro_average_j1, aj0x, aj1x

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      real, dimension(:), intent(out) :: pflx, vflx, qflx
      real, dimension(:, :), intent(out) :: pflx_x, vflx_x, qflx_x
      real, dimension(:, :, -nzgrid:, :, :), intent(out) :: pflx_vs_kxkyz, vflx_vs_kxkyz, qflx_vs_kxkyz

      integer :: ivmu, imu, iv, iz, it, is, ia
      real :: flx_norm
      complex, dimension(:, :), allocatable :: g0k, g1k

      pflx = 0.; vflx = 0.; qflx = 0.
      pflx_x = 0.; vflx_x = 0.; qflx_x = 0.
      pflx_vs_kxkyz = 0.; vflx_vs_kxkyz = 0.; qflx_vs_kxkyz = 0.

      ia = 1
      if (flux_norm) then
         flx_norm = 1./grho_norm
      else
         flx_norm = 1.
      end if

      allocate (g0k(naky, nakx))
      allocate (g1k(naky, nakx))

      ! FLAG - electrostatic for now
      ! get electrostatic contributions to fluxes

      if (fphi > epsilon(0.0)) then
         ia = 1

         !get particle flux
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            call gyro_average(g(:, :, :, :, ivmu), ivmu, g1(:, :, :, :, ivmu))

            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  if (radial_variation) then
                     g0k = g1(:, :, iz, it, ivmu) &
                           * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                              * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                              * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                              + dBdrho(iz) / bmag(ia, iz))

                     call multiply_by_rho(g0k)
                     g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) + g0k
                  end if

                  !subtract adiabatic contribution part of g
                  g0k = spec(is)%zt * fphi * phi(:, :, iz, it) * aj0x(:, :, iz, ivmu)**2
                  if (.not. maxwellian_normalization) then
                     g0k = g0k * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                  end if
                  if (radial_variation) then
                     g1k = g0k * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                                  - spec(is)%fprim - 2.0 * dBdrho(iz) * mu(imu) &
                                  - aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                                  * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                  * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                  + dBdrho(iz) / bmag(ia, iz))

                     call multiply_by_rho(g1k)

                     g0k = g0k + g1k
                  end if
                  g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) + g0k

               end do
            end do
         end do
         call get_one_flux_vmulo(flx_norm * spec%dens_psi0, g1, phi, pflx)
         call get_one_flux_vmulo(flx_norm * spec%dens_psi0, g1, phi, pflx_vs_kxkyz)

         if (write_radial_fluxes) then
            call get_one_flux_radial(flx_norm * spec%dens_psi0, g1, phi, pflx_x)
         end if

         !get heat flux
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            call gyro_average(g(:, :, :, :, ivmu), ivmu, g1(:, :, :, :, ivmu))
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid

                  g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) * (vpa(iv)**2 + vperp2(ia, iz, imu))

                  if (radial_variation) then
                     g0k = g1(:, :, iz, it, ivmu) &
                           * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                              * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                              * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                              + dBdrho(iz) / bmag(ia, iz) &
                              + 2.0 * mu(imu) * dBdrho(iz) / (vpa(iv)**2 + vperp2(ia, iz, imu)))

                     call multiply_by_rho(g0k)

                     g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) + g0k

                  end if

                  !subtract adiabatic contribution part of g
                  g0k = spec(is)%zt * fphi * phi(:, :, iz, it) * aj0x(:, :, iz, ivmu)**2 &
                        * (vpa(iv)**2 + vperp2(ia, iz, imu))
                  if (.not. maxwellian_normalization) then
                     g0k = g0k * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                  end if
                  if (radial_variation) then
                     g1k = g0k * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                                  - spec(is)%fprim - 2.0 * dBdrho(iz) * mu(imu) &
                                  - aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                                  * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                  * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                  + dBdrho(iz) / bmag(ia, iz) &
                                  + 2.0 * mu(imu) * dBdrho(iz) / (vpa(iv)**2 + vperp2(ia, iz, imu)))

                     call multiply_by_rho(g1k)

                     g0k = g0k + g1k
                  end if
                  g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) + g0k
               end do
            end do
         end do
         call get_one_flux_vmulo(flx_norm * spec%dens_psi0 * spec%temp_psi0, g1, phi, qflx)
         call get_one_flux_vmulo(flx_norm * spec%dens_psi0 * spec%temp_psi0, g1, phi, qflx_vs_kxkyz)

         if (write_radial_fluxes) then
            call get_one_flux_radial(flx_norm * spec%dens_psi0 * spec%temp_psi0, g1, phi, qflx_x)
         end if

         ! get momentum flux
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  ! parallel component
                  g0k = g(:, :, iz, it, ivmu) * vpa(iv) * geo_surf%rmaj * btor(iz) / bmag(ia, iz)
                  call gyro_average(g0k, iz, ivmu, g1(:, :, iz, it, ivmu))

                  if (radial_variation) then
                     g0k = g1(:, :, iz, it, ivmu) &
                           * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                              * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                              * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                              + dIdrho / (geo_surf%rmaj * btor(iz)))

                     call multiply_by_rho(g0k)

                     g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) + g0k

                  end if
                  !subtract adiabatic contribution part of g
                  g0k = spec(is)%zt * fphi * phi(:, :, iz, it) * aj0x(:, :, iz, ivmu)**2 &
                        * vpa(iv) * geo_surf%rmaj * btor(iz) / bmag(ia, iz)
                  if (.not. maxwellian_normalization) then
                     g0k = g0k * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                  end if
                  if (radial_variation) then
                     g1k = g0k * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                                  - spec(is)%fprim - 2.0 * dBdrho(iz) * mu(imu) &
                                  - aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                                  * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                  * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                  + dIdrho / (geo_surf%rmaj * btor(iz)))

                     call multiply_by_rho(g1k)

                     g0k = g0k + g1k
                  end if
                  g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) + g0k

                  ! perpendicular component
                  g0k = -g(:, :, iz, it, ivmu) * zi * spread(aky, 2, nakx) * vperp2(ia, iz, imu) * geo_surf%rhoc &
                        * (gds21(ia, iz) + theta0 * gds22(ia, iz)) * spec(is)%smz &
                        / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2)

                  call gyro_average_j1(g0k, iz, ivmu, g2(:, :, iz, it, ivmu))
                  if (radial_variation) then
                     g0k = -g(:, :, iz, it, ivmu) * zi * spread(aky, 2, nakx) * vperp2(ia, iz, imu) * geo_surf%rhoc &
                           * (dgds21dr(ia, iz) + theta0 * dgds22dr(ia, iz)) * aj1x(:, :, iz, ivmu) * spec(is)%smz &
                           / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2)

                     g0k = g0k - g(:, :, iz, it, ivmu) * zi * spread(aky, 2, nakx) * vperp2(ia, iz, imu) * geo_surf%rhoc &
                           * (gds21(ia, iz) + theta0 * gds22(ia, iz)) * spec(is)%smz &
                           / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2) &
                           * (0.5 * aj0x(:, :, iz, ivmu) - aj1x(:, :, iz, ivmu)) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz))

                     g0k = g0k + g2(:, :, iz, it, ivmu) &
                           * (-geo_surf%d2qdr2 * geo_surf%rhoc / (geo_surf%shat * geo_surf%qinp) &
                              - geo_surf%d2psidr2 * drhodpsi)

                     call multiply_by_rho(g0k)

                     g2(:, :, iz, it, ivmu) = g2(:, :, iz, it, ivmu) + g0k
                  end if

                  !subtract adiabatic contribution part of g
                  g0k = -spec(is)%zt * fphi * phi(:, :, iz, it) * aj0x(:, :, iz, ivmu) * aj1x(:, :, iz, ivmu) &
                        * zi * spread(aky, 2, nakx) * vperp2(ia, iz, imu) * geo_surf%rhoc &
                        * (gds21(ia, iz) + theta0 * gds22(ia, iz)) * spec(is)%smz &
                        / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2)
                  if (.not. maxwellian_normalization) then
                     g0k = g0k * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                  end if

                  if (radial_variation) then
                     g1k = -spec(is)%zt * fphi * phi(:, :, iz, it) * aj0x(:, :, iz, ivmu) * aj1x(:, :, iz, ivmu) &
                           * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is) &
                           * zi * spread(aky, 2, nakx) * vperp2(ia, iz, imu) * geo_surf%rhoc &
                           * (dgds21dr(ia, iz) + theta0 * dgds22dr(ia, iz)) * spec(is)%smz &
                           / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2)

                     g1k = g1k - spec(is)%zt * fphi * phi(:, :, iz, it) * aj0x(:, :, iz, ivmu) &
                           * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is) &
                           * zi * spread(aky, 2, nakx) * vperp2(ia, iz, imu) * geo_surf%rhoc &
                           * (gds21(ia, iz) + theta0 * gds22(ia, iz)) * spec(is)%smz &
                           / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2) &
                           * (0.5 * aj0x(:, :, iz, ivmu) - aj1x(:, :, iz, ivmu)) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz))

                     g1k = g1k + &
                           g0k * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                                  - spec(is)%fprim - 2.0 * dBdrho(iz) * mu(imu) &
                                  - 0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                                  * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                  * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                  - geo_surf%d2qdr2 * geo_surf%rhoc / (geo_surf%shat * geo_surf%qinp) &
                                  - geo_surf%d2psidr2 * drhodpsi)

                     call multiply_by_rho(g1k)

                     g0k = g0k + g1k
                  end if
                  g2(:, :, iz, it, ivmu) = g2(:, :, iz, it, ivmu) + g0k
               end do
            end do
         end do

         g1 = g1 + g2
         call get_one_flux_vmulo(flx_norm * spec%dens_psi0 * sqrt(spec%mass * spec%temp_psi0), g1, phi, vflx)
         call get_one_flux_vmulo(flx_norm * spec%dens_psi0 * sqrt(spec%mass * spec%temp_psi0), g1, phi, vflx_vs_kxkyz)

         if (write_radial_fluxes) then
            call get_one_flux_radial(flx_norm * spec%dens_psi0 * sqrt(spec%mass * spec%temp_psi0), g1, phi, vflx_x)
         end if

      end if

      if (allocated(g0k)) deallocate (g0k)
      if (allocated(g1k)) deallocate (g1k)

   end subroutine get_fluxes_vmulo

   !> Calculate the total particle, momentum and heat fluxes (pflx, vflx, qflx)
   !> and the contributions from a given (kx,ky,z) location (pflx_kxkyz, vflx_kxkyz, qflx_kxkyz)
   !> inputs are the particle density (dens), parallel flow (upar) and pressure (pres)
   subroutine get_fluxes_ffs(dens, upar, pres, pflx, vflx, qflx, &
                             pflx_vs_kxkyz, vflx_vs_kxkyz, qflx_vs_kxkyz)

      use constants, only: zi
      use zgrid, only: nzgrid, delzed
      use kt_grids, only: naky, nakx, ny
      use kt_grids, only: aky, dy
      use fields_arrays, only: phi
      use stella_geometry, only: grad_x, jacob

      use physics_flags, only: const_alpha_geo

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: dens, upar, pres
      real, dimension(:), intent(out) :: pflx, vflx, qflx
      real, dimension(:, :, -nzgrid:, :, :), intent(out) :: pflx_vs_kxkyz, vflx_vs_kxkyz, qflx_vs_kxkyz

      integer :: iky, it
      real, dimension(:, :), allocatable :: flxfac
      complex, dimension(:, :, :), allocatable :: dphidy

      integer :: iy, iz

      !> assume a single flux annulus
      it = 1

      !> Set all fluxes to zero just in case
      pflx = 0.; vflx = 0.; qflx = 0.
      pflx_vs_kxkyz = 0.; vflx_vs_kxkyz = 0.; qflx_vs_kxkyz = 0.

      allocate (dphidy(naky, nakx, -nzgrid:nzgrid))

      !> Obtain the y-component of the electric field that appears as a factor
      !> in the flux expression due to the radial component of the ExB velocity
      do iky = 1, naky
	 !dphidy(iky, :, :) = zi * aky(iky) * conjg(phi(iky, :, :, it))
         dphidy(iky, :, :) = phi(iky, :, :, it) * aky(iky)
      end do
      
      !> Calculate Jacobian for spacial integral 
      allocate (flxfac(ny, -nzgrid:nzgrid))

      flxfac = spread(delzed * dy, 1, ny) * jacob
      flxfac(:, -nzgrid) = 0.5 * flxfac(:, -nzgrid)
      flxfac(:, nzgrid) = 0.5 * flxfac(:, nzgrid)
      
      flxfac = flxfac/(sum(flxfac * grad_x) )

      if(const_alpha_geo) flxfac = flxfac / ny
      
      call get_one_flux_ffs(dens, dphidy, flxfac, pflx, pflx_vs_kxkyz(:, :, :, it, :))
      call get_one_flux_ffs(pres, dphidy, flxfac, qflx, qflx_vs_kxkyz(:, :, :, it, :))
      call get_one_flux_ffs(upar, dphidy, flxfac, vflx, vflx_vs_kxkyz(:, :, :, it, :))

      deallocate (dphidy, flxfac)

   end subroutine get_fluxes_ffs

   subroutine get_one_flux_ffs(mom, dphidy, flxfac, flx, flx_vs_kxkyz)

      use species, only: nspec
      use zgrid, only: nzgrid, delzed
      use kt_grids, only: naky, nakx
      use volume_averages, only: mode_fac

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: mom
      complex, dimension(:, :, -nzgrid:), intent(in) :: dphidy
      real, dimension(:, -nzgrid:), intent(in) :: flxfac
      real, dimension(:), intent(out) :: flx
      real, dimension(:, :, -nzgrid:, :), intent(out) :: flx_vs_kxkyz

      integer :: iky, ikx, iz, is
      complex, dimension(:, :, :, :), allocatable :: mom_ky
      complex, dimension (:,:), allocatable :: flx_x

      allocate (mom_ky(naky, nakx, -nzgrid:nzgrid, nspec))

      flx = 0.
      flx_vs_kxkyz = 0.

      !> Multiply input moment (e.g. density, momentum, or energy) but Jacobian (called flxfac)
      !> Need to do this because Jacobian has y-dependence 
      !> then Fourier transform back to (ky,kx)-space
      
      call get_modified_fourier_coefficient(mom, mom_ky, flxfac)
      do is = 1, nspec
         !> pflx_vs_kxkyz is the particle flux before summing over (kx,ky) and integrating over z
!flx_vs_kxkyz(:, :, :, is) = flxfac * aimag(mom_ky(:, :, :, is) * dphidy)
         flx_vs_kxkyz(:, :, :, is) = aimag(mom_ky(:, :, :, is) * conjg(dphidy(:, :, :)))
         !> calculate the volume average of the particle flux
         !> note that the factor of 1/B that appears in the Jacobian has already been taken into account
         !> in the numerator of the flux surface average
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               do iky = 1, naky
                  !if(iky == ikx == 1) then
                  !   cycle 
                  !else
                  flx(is) = flx(is) + 0.5 * mode_fac(iky) * flx_vs_kxkyz(iky, ikx, iz, is)
                  !end if
               end do
            end do
         end do
      end do

      deallocate (mom_ky)

   end subroutine get_one_flux_ffs

   subroutine get_modified_fourier_coefficient(moment, moment_ky, flxfac)

      use species, only: nspec
      use zgrid, only: nzgrid
      use kt_grids, only: ikx_max, naky_all, ny
      use kt_grids, only: swap_kxky_back
      use stella_geometry, only: bmag
      use stella_transforms, only: transform_y2ky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: moment
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: moment_ky

      integer :: ikx, iz, is
      complex, dimension(:, :), allocatable :: tmp_kykx
      complex, dimension(:, :), allocatable :: tmp_ykx
      real, dimension(:, -nzgrid:), intent(in) :: flxfac

      integer :: is_end

      allocate (tmp_kykx(naky_all, ikx_max))
      allocate (tmp_ykx(ny, ikx_max))
      !> is_end is either nspec for fluxes, or ntubes when calcualting |phi|^2
      is_end = size(moment, dim=4)

      do is = 1, is_end
         do iz = -nzgrid, nzgrid
            do ikx = 1, ikx_max
               !> divide the input moment by the magnetic field strength
               !> to account for Jacobian in flux-surface average
               tmp_ykx(:, ikx) = moment(:, ikx, iz, is) * flxfac(:, iz)
            end do
            !> transform the B-modified input moment from y to ky space
            call transform_y2ky(tmp_ykx, tmp_kykx)
            !> swap from all ky and kx >= 0 to all kx and ky >= 0
            call swap_kxky_back(tmp_kykx, moment_ky(:, :, iz, is))
         end do
      end do
      deallocate (tmp_kykx, tmp_ykx)

   end subroutine get_modified_fourier_coefficient

   !============================================================================
   !======================== CALCULATE ONE FLUX(IKY,IZ) ========================
   !============================================================================
   ! For example the particle flux is defined as,
   !      Re[Γ̃_{s,k}] = FAC*int(Im[conj(phi)*VELOCITY_INTEGRAL]*J*dz)/int(J*dz)
   ! With FAC = -sgn(psi_t)*(ñ_s/2)*(k̃_y/<|\tilde{∇}ρ>_ζ)
   ! With VELOCITY_INTEGRAL = (2B̃/√π) int(int(g̃*J_0) dṽparallel) dμ̃)
   ! With <norm>(iz) = J(iz)*dz / <|\tilde{∇}ρ>_ζ and <|\tilde{∇}ρ0>_ζ = sum(<grho>(iz)*J(iz)*dz)
   ! With <grho> = a|∇ρ0| = a (dρ0/dψ) ∇ψ = 1/ρ0 * ∇ψ/(a*Bref) = sqrt(|grad_psi_grad_psi|)/ρ0
   subroutine get_one_flux(iky, iz, norm, gin, fld, flxout)

      use vpamu_grids, only: integrate_vmu
      use volume_averages, only: mode_fac
      use stella_geometry, only: sign_torflux
      use kt_grids, only: aky

      implicit none

      integer, intent(in) :: iky, iz
      real, intent(in) :: norm
      complex, dimension(:, :), intent(in) :: gin
      complex, intent(in) :: fld
      real, intent(in out) :: flxout

      complex :: flx

      call integrate_vmu(gin, iz, flx)
      flxout = flxout - sign_torflux * 0.5 * mode_fac(iky) * aky(iky) * aimag(flx * conjg(fld)) * norm

   end subroutine get_one_flux

   !==============================================
   !============ GET ONE FLUX VMULO ==============
   !==============================================
   subroutine get_one_flux_vmulo_int(weights, gin, fld, flxout)

      use vpamu_grids, only: integrate_vmu
      use stella_layouts, only: vmu_lo
      use kt_grids, only: aky, nakx, naky, boundary_size
      use zgrid, only: nzgrid, ntubes
      use species, only: nspec
      use volume_averages, only: mode_fac
      use stella_geometry, only: dVolume, sign_torflux
      use stella_transforms, only: transform_kx2x_unpadded
      use physics_flags, only: radial_variation

      implicit none

      real, dimension(:), intent(in) :: weights
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: fld
      real, dimension(:), intent(in out) :: flxout

      complex, dimension(:, :, :, :, :), allocatable :: totals
      complex, dimension(:, :), allocatable :: g0x, g1x

      integer :: ia, is, it, iz, ikx
      real, dimension(nspec) :: flux_sum
      real :: volume, factor

      allocate (totals(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))

      ia = 1
      flux_sum = 0.

      ! The factor in front of all the flux definitions is <factor> = -sgn(psi_t)/2
      ! The constants which differ for each flux are gathered in <weights> and added in <integrate_vmu()>
      ! For the particle flux <weights> = ñ_s/<|\tilde{∇}ρ>_ζ = <flx_norm> * <spec%dens_psi0>
      ! For the heat flux <weights> = ñ_s*T̃_s/<|\tilde{∇}ρ>_ζ = <flx_norm> * <spec%dens_psi0> * <spec%temp_psi0>
      ! For the momentum flux <weights> = ñ_s*sqrt(m̃*T̃_s)/<|\tilde{∇}ρ>_ζ = <flx_norm> * <spec%dens_psi0> * sqrt(<spec%mass> * <spec%temp_psi0>)
      factor = -sign_torflux * 0.5

      call integrate_vmu(gin, weights, totals)
      if (radial_variation) then !do it in real-space
         allocate (g0x(naky, nakx))
         allocate (g1x(naky, nakx))
         do is = 1, nspec
            volume = 0.
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call transform_kx2x_unpadded(totals(:, :, iz, it, is), g0x)
                  call transform_kx2x_unpadded(fld(:, :, iz, it), g1x)
                  do ikx = boundary_size + 1, nakx - boundary_size
                     flux_sum(is) = flux_sum(is) + &
                                    sum(factor * mode_fac * aky * aimag(g0x(:, ikx) * conjg(g1x(:, ikx))) * dVolume(ia, ikx, iz))
                     volume = volume + dVolume(ia, ikx, iz)
                  end do
               end do
            end do
         end do
         deallocate (g0x, g1x)
      else
         do is = 1, nspec
            volume = 0.
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do ikx = 1, nakx
                     flux_sum(is) = flux_sum(is) + &
                                    sum(factor * mode_fac * aky * aimag(totals(:, ikx, iz, it, is) * conjg(fld(:, ikx, iz, it))) * dVolume(ia, 1, iz))
                  end do
                  volume = volume + dVolume(ia, 1, iz)
               end do
            end do
         end do
      end if

      flxout = flxout + flux_sum / volume

      deallocate (totals)

   end subroutine get_one_flux_vmulo_int

   subroutine get_one_flux_vmulo_kxkyz(weights, gin, fld, flxout)

      use vpamu_grids, only: integrate_vmu
      use stella_layouts, only: vmu_lo
      use kt_grids, only: aky, nakx, naky
      use zgrid, only: nzgrid, ntubes
      use species, only: nspec
      use volume_averages, only: mode_fac
      use stella_geometry, only: sign_torflux

      implicit none

      real, dimension(:), intent(in) :: weights
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: fld
      real, dimension(:, :, -nzgrid:, :, :), intent(in out) :: flxout

      complex, dimension(:, :, :, :, :), allocatable :: totals

      integer :: ia, is, it, iz, ikx
      real :: factor

      allocate (totals(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))

      ia = 1
      factor = -sign_torflux * 0.5
      call integrate_vmu(gin, weights, totals)
      do is = 1, nspec
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  flxout(:, ikx, iz, it, is) = factor * mode_fac * aky &
                                               * aimag(totals(:, ikx, iz, it, is) * conjg(fld(:, ikx, iz, it)))
               end do
            end do
         end do
      end do

      deallocate (totals)

   end subroutine get_one_flux_vmulo_kxkyz

   !==============================================
   !=========== GET ONE FLUX RADIAL ==============
   !==============================================
   subroutine get_one_flux_radial(weights, gin, fld, flxout)

      use vpamu_grids, only: integrate_vmu
      use stella_geometry, only: dVolume
      use stella_layouts, only: vmu_lo
      use kt_grids, only: aky, nakx, naky
      use zgrid, only: nzgrid, ntubes
      use species, only: nspec
      use volume_averages, only: mode_fac
      use stella_transforms, only: transform_kx2x_unpadded
      use stella_geometry, only: sign_torflux

      implicit none

      real, dimension(:), intent(in) :: weights
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: fld
      real, dimension(:, :), intent(in out) :: flxout

      real, dimension(:), allocatable :: dV_rad
      complex, dimension(:, :, :, :, :), allocatable :: totals

      complex, dimension(:, :), allocatable :: g0x, g1x

      integer :: ia, is, it, iz, ikx
      real :: factor

      allocate (dV_rad(nakx))
      allocate (g0x(naky, nakx))
      allocate (g1x(naky, nakx))
      allocate (totals(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))

      ia = 1

      dV_rad = sum(sum(dVolume, 3), 1) * ntubes

      ! NB: this returns the flux-surface-averaged radial fluxes. Keep in mind that the
      !     volume element in a flux-surface, dV, may not be uniform across surfaces, so
      !     one cannot simply sum across the radius here to get the total flux; rather, one
      !     would have to multiply by dV/V across the radius first
      call integrate_vmu(gin, weights, totals)
      factor = -sign_torflux * 0.5
      do is = 1, nspec
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               call transform_kx2x_unpadded(totals(:, :, iz, it, is), g0x)
               call transform_kx2x_unpadded(fld(:, :, iz, it), g1x)
               do ikx = 1, nakx
                  flxout(ikx, is) = flxout(ikx, is) &
                                    + sum(factor * mode_fac * aky * aimag(g0x(:, ikx) * conjg(g1x(:, ikx))) * dVolume(ia, ikx, iz) / dV_rad(ikx))
               end do
            end do
         end do
      end do

      deallocate (dV_rad, g0x, g1x, totals)

   end subroutine get_one_flux_radial

   !==============================================
   !=============== GET MOMENTS ==================
   !==============================================
   subroutine get_moments(g, dens, upar, temp, dens_x, upar_x, temp_x, spitzer2)

      use zgrid, only: nzgrid, ntubes
      use species, only: spec, nspec
      use vpamu_grids, only: integrate_vmu
      use vpamu_grids, only: vpa, vperp2, mu
      use vpamu_grids, only: maxwell_mu, ztmax, maxwell_fac, maxwell_vpa
      use kt_grids, only: ny, naky, nakx, multiply_by_rho, rho_d_clamped
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use dist_fn_arrays, only: g1, g2, kperp2, dkperp2dr
      use stella_geometry, only: bmag, dBdrho
      use stella_geometry, only: dl_over_b, d_dl_over_b_drho
      use gyro_averages, only: aj0x, aj1x, gyro_average
      use fields_arrays, only: phi, phi_corr_QN, phi_proj
      use run_parameters, only: fphi
      use run_parameters, only: maxwellian_normalization
      use physics_flags, only: radial_variation
      use stella_transforms, only: transform_kx2x_unpadded

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: dens, upar, temp, spitzer2
      real, dimension(:, :), intent(out) :: dens_x, upar_x, temp_x

      complex, dimension(:, :), allocatable :: g0k, g1k, g1x
      real :: zero

      integer :: ivmu, iv, imu, is, ia
      integer :: iz, it

      if (radial_variation) then
         allocate (g0k(naky, nakx))
      end if
      if (write_radial_moments) then
         allocate (g1k(1, nakx))
         allocate (g1x(1, nakx))
      end if

      ! Hack below. Works since J0^2 - 1 and its derivative are zero at the origin
      zero = 100.*epsilon(0.)

      ! h is gyrophase independent, but is in gyrocenter coordinates,
      ! so requires a J_0 to get to particle coordinates
      ! <f>_r = h J_0 - Ze*phi/T * F0
      ! g     = h     - Ze*<phi>_R/T * F0
      ! <f>_r = g J_0 + Ze*(J_0<phi>_R-phi)/T * F0

      ! calculate the integrand appearing in the integral for the density
      ia = 1
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         ! obtain the gyro-average of g that appears in the density integral
         call gyro_average(g(:, :, :, :, ivmu), ivmu, g1(:, :, :, :, ivmu))
         ! FLAG -- AJ0X NEEDS DEALING WITH BELOW
         g2(:, :, :, :, ivmu) = spread(aj0x(:, :, :, ivmu)**2 - 1.0, 4, ntubes) * spec(is)%zt * fphi * phi
         if (.not. maxwellian_normalization) then
            g2(:, :, :, :, ivmu) = g2(:, :, :, :, ivmu) * maxwell_vpa(iv, is) * &
                                   spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx) * maxwell_fac(is), 4, ntubes)
         end if
         g2(:, :, :, :, ivmu) = g2(:, :, :, :, ivmu) + g1(:, :, :, :, ivmu)
         ! g2(:, :, :, :, ivmu) = g1(:, :, :, :, ivmu) + ztmax(iv, is) &
         !                        * spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx) &
         !                                 * maxwell_fac(is) * (aj0x(:, :, :, ivmu)**2 - 1.0), 4, ntubes) * fphi * phi

         if (radial_variation) then
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  !phi
                  g0k = ztmax(iv, is) * maxwell_mu(ia, iz, imu, is) &
                        * maxwell_fac(is) * (aj0x(:, :, iz, ivmu)**2 - 1.0) * fphi * phi(:, :, iz, it) &
                        * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                           - spec(is)%fprim + (dBdrho(iz) / bmag(ia, iz)) * (1.0 - 2.0 * mu(imu) * bmag(ia, iz)) &
                           - aj1x(:, :, iz, ivmu) * aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                           * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                           / (aj0x(:, :, iz, ivmu)**2 - 1.0 + zero))

                  !g
                  g0k = g0k + g1(:, :, iz, it, ivmu) &
                        * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                           * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                           + dBdrho(iz) / bmag(ia, iz))

                  call multiply_by_rho(g0k)

                  ! g0k(1,1) = 0.0

                  !phi QN
                  g0k = g0k + ztmax(iv, is) * maxwell_mu(ia, iz, imu, is) * fphi &
                        * maxwell_fac(is) * (aj0x(:, :, iz, ivmu)**2 - 1.0) * phi_corr_QN(:, :, iz, it)

                  g2(:, :, iz, it, ivmu) = g2(:, :, iz, it, ivmu) + g0k
               end do
            end do
         end if
      end do
      call integrate_vmu(g2, spec%dens_psi0, dens)

      if (write_radial_moments) then
         dens_x = 0.0
         do is = 1, nspec
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g1k(1, :) = dens(1, :, iz, it, is) - phi_proj(:, 1, it)
                  call transform_kx2x_unpadded(g1k, g1x)
                  dens_x(:, is) = dens_x(:, is) &
                                  + real(g1x(1, :) * (dl_over_b(ia, iz) + rho_d_clamped * d_dl_over_b_drho(ia, iz)))
               end do
            end do
         end do
         dens_x = dens_x / ntubes
      end if

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         if (maxwellian_normalization) then
            g2(:, :, :, :, ivmu) = (g1(:, :, :, :, ivmu) + spec(is)%zt &
                                    * spread(aj0x(:, :, :, ivmu)**2 - 1.0, 4, ntubes) * phi * fphi) &
                                   * (vpa(iv)**2 + spread(spread(spread(vperp2(1, :, imu), 1, naky), 2, nakx), 4, ntubes) - 1.5) / 1.5
         else
            g2(:, :, :, :, ivmu) = (g1(:, :, :, :, ivmu) + ztmax(iv, is) &
                                    * spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx) &
                                             * maxwell_fac(is) * (aj0x(:, :, :, ivmu)**2 - 1.0), 4, ntubes) * phi * fphi) &
                                   * (vpa(iv)**2 + spread(spread(spread(vperp2(1, :, imu), 1, naky), 2, nakx), 4, ntubes) - 1.5) / 1.5
         end if
         if (radial_variation) then
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  !phi
                  g0k = ztmax(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is) &
                        * (vpa(iv)**2 + vperp2(ia, iz, imu) - 1.5) / 1.5 &
                        * (aj0x(:, :, iz, ivmu)**2 - 1.0) * phi(:, :, iz, it) * fphi &
                        * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                           - spec(is)%fprim + (dBdrho(iz) / bmag(ia, iz)) * (1.0 - 2.0 * mu(imu) * bmag(ia, iz)) &
                           - aj1x(:, :, iz, ivmu) * aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                           * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                           / (aj0x(:, :, iz, ivmu)**2 - 1.0 + zero) &
                           + 2.0 * mu(imu) * dBdrho(iz) / (vpa(iv)**2 + vperp2(ia, iz, imu) - 1.5))

                  !g
                  g0k = g0k + g1(:, :, iz, it, ivmu) * (vpa(iv)**2 + vperp2(ia, iz, imu) - 1.5) / 1.5 &
                        * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                           * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                           + dBdrho(iz) / bmag(ia, iz) &
                           + 2.0 * mu(imu) * dBdrho(iz) / (vpa(iv)**2 + vperp2(ia, iz, imu) - 1.5))

                  call multiply_by_rho(g0k)

                  !phi QN
                  g0k = g0k + fphi * ztmax(iv, is) * maxwell_mu(ia, iz, imu, is) &
                        * (vpa(iv)**2 + vperp2(ia, iz, imu) - 1.5) / 1.5 &
                        * maxwell_fac(is) * (aj0x(:, :, iz, ivmu)**2 - 1.0) * phi_corr_QN(:, :, iz, it)

                  g2(:, :, iz, it, ivmu) = g2(:, :, iz, it, ivmu) + g0k
               end do
            end do
         end if
      end do
      ! integrate to get dTs/Tr
      !    call integrate_vmu (g2, spec%temp, temp)
      call integrate_vmu(g2, spec%temp_psi0 * spec%dens_psi0, temp)

      if (write_radial_moments) then
         temp_x = 0.0
         do is = 1, nspec
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g1k(1, :) = temp(1, :, iz, it, is)
                  call transform_kx2x_unpadded(g1k, g1x)
                  temp_x(:, is) = temp_x(:, is) &
                                  + real(g1x(1, :) * (dl_over_b(ia, iz) + rho_d_clamped * d_dl_over_b_drho(ia, iz)))
               end do
            end do
         end do
         temp_x = temp_x / ntubes
      end if

      ! for Spitzer problem tests of the collision operator
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         g2(:,:,:,:,ivmu) = g(:,:,:,:,ivmu)*( vpa(iv)*(vpa(iv)**2+spread(spread(spread(vperp2(1,:,imu),1,naky),2,nakx),4,ntubes)) - 5./2.*vpa(iv) )
      end do
      call integrate_vmu(g2, spec%stm, spitzer2) ! AVB: stm is the thermal speed

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         g2(:, :, :, :, ivmu) = vpa(iv) * g1(:, :, :, :, ivmu)
         if (radial_variation) then
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  !g
                  g0k = vpa(iv) * g1(:, :, iz, it, ivmu) &
                        * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                           * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                           + dBdrho(iz) / bmag(ia, iz))

                  call multiply_by_rho(g0k)

                  g2(:, :, iz, it, ivmu) = g2(:, :, iz, it, ivmu) + g0k
               end do
            end do
         end if
      end do
      call integrate_vmu(g2, spec%stm_psi0, upar)
      if (write_radial_moments) then
         upar_x = 0.0
         do is = 1, nspec
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g1k(1, :) = upar(1, :, iz, it, is)
                  call transform_kx2x_unpadded(g1k, g1x)
                  upar_x(:, is) = upar_x(:, is) &
                                  + real(g1x(1, :) * (dl_over_b(ia, iz) + rho_d_clamped * d_dl_over_b_drho(ia, iz)))
               end do
            end do
         end do
         upar_x = upar_x / ntubes
      end if

      if (allocated(g0k)) deallocate (g0k)
      if (allocated(g1k)) deallocate (g1k)
      if (allocated(g1x)) deallocate (g1x)

   end subroutine get_moments

   subroutine get_moments_ffs(g, dens, upar, pres)

      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use species, only: spec, nspec
      use zgrid, only: nzgrid
      use vpamu_grids, only: integrate_vmu_ffs
      use vpamu_grids, only: vpa, vperp2
      use kt_grids, only: naky_all, ikx_max, ny
      use kt_grids, only: swap_kxky
      use dist_fn_arrays, only: g0, g1, g2
      use gyro_averages, only: gyro_average, j0_ffs
      use fields_arrays, only: phi
      use stella_transforms, only: transform_ky2y
      use kt_grids, only: akx, aky
      use kt_grids, only: nakx, theta0, naky
      use constants, only: zi, pi
      use zgrid, only: ntubes

      use run_parameters, only: fphi
      use g_tofrom_h, only: g_to_h

      !> For momentum flux
      use gyro_averages, only: j1_ffs
      use stella_geometry, only: gds21, gds22, gds2
      use stella_geometry, only: geo_surf
      use stella_geometry, only: gradzeta_grady, gradzeta_gradx, gradpar_zeta

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dens, upar, pres

      real, dimension(:), allocatable :: dens_wgts, upar_wgts, pres_wgts
      !> f_swap will contain delta f(ky,kx) on a grid with all kys and kx >= 0
      complex, dimension(:, :), allocatable :: f_swap
      !> fy will contain delta f(y,kx) on a grid with kx >= 0
      complex, dimension(:, :, :), allocatable :: fy, f2y, f3y
      !> integrand will contain the integrand in the velocity moment integrals
      complex, dimension(:), allocatable :: integrand

      integer :: iy, ikx, iz, it
      integer :: ivmu, iv, imu, is
      real :: fac1, fac2

      !> species-dependent factor by which velocity moments must be multiplied
      !> to get density, pressure, etc.
      allocate (dens_wgts(nspec))
      allocate (upar_wgts(nspec))
      allocate (pres_wgts(nspec))

      dens = 0.; upar = 0.; pres = 0.
      !> set species-dependent factors needed for density, parallel flow and pressure
      dens_wgts = spec%dens
      pres_wgts = spec%dens * spec%temp
      upar_wgts = spec%stm !spec%dens * sqrt(spec%mass * spec%temp)

      !> Already allocated arrays (allocated in time_advance.f90). Set to zero just in case
      g0 = 0. ; g1 = 0. ; g2 = 0.

      g0 = g
      !> This returns f = g + Z/T * (<phi> - phi)
      !> NB: g_to_f0 TODO-GA: check 
      call g_to_f1(g, phi, g0)
      !> g1 = J0 * f
      call gyro_average(g0, g1, j0_ffs)
      !> g2 = J1 * f
      call gyro_average(g0, g2, j1_ffs)
      
      allocate (f_swap(naky_all, ikx_max))
      allocate (integrand(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      allocate (fy(ny, ikx_max, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); fy = 0.
      allocate (f2y(ny, ikx_max, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); f2y = 0.
      allocate (f3y(ny, ikx_max, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); f3y = 0.

      !> assume only a single flux annulus
      it = 1
      do iz = -nzgrid, nzgrid
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            !> for every (z,vpa,mu,spec) point
            !> switch from ky >= 0 and kx = [-kxmax, kxmax]
            !> to ky = [-kymax, kymax] and kx >= 0
            call swap_kxky(g1(:, :, iz, it, ivmu), f_swap)
            !> for every (z,vpa,mu,spec) point, Fourier tranform from ky to y space to get
            !> the kx component of <f(y,x)>_r
            call transform_ky2y(f_swap, fy(:, :, ivmu))

            !! The following are only needed for the momentum flux because we need J1*f 
            !! J1* zi * ky * f
            g2(:, :, iz, it, ivmu) = zi * g2(:, :, iz, it, ivmu) * spread(aky, 2, nakx)
            call swap_kxky(g2(:, :, iz, it, ivmu), f_swap)
            call transform_ky2y(f_swap, f2y(:, :, ivmu))
            !! J1 * zi * kx * f
            g2(:, :, iz, it, ivmu) = g2(:, :, iz, it, ivmu) * theta0(:, :)
            call swap_kxky(g2(:, :, iz, it, ivmu), f_swap)
            call transform_ky2y(f_swap, f3y(:, :, ivmu))
         end do

         do ikx = 1, ikx_max
            do iy = 1, ny
               !> the integrand for the density moment is the distribution function
               integrand = fy(iy, ikx, :)
               !> integrate over v-space to get the density, normalised by the reference density.
               call integrate_vmu_ffs(integrand, dens_wgts, iy, iz, dens(iy, ikx, iz, :))

               !> the integrand for the pressure moment is the energy-weighted distribution function
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                  iv = iv_idx(vmu_lo, ivmu)
                  imu = imu_idx(vmu_lo, ivmu)
                  integrand(ivmu) = fy(iy, ikx, ivmu) * (vpa(iv)**2 + vperp2(iy, iz, imu))
               end do
               !> integrate over v-space to get the pressure, normalised by the reference pressure.
               call integrate_vmu_ffs(integrand, pres_wgts, iy, iz, pres(iy, ikx, iz, :))

               fac1 = gradzeta_grady(iy, iz) * gds21(iy, iz) / geo_surf%shat - gradzeta_gradx(iy, iz) * gds2(iy, iz)
               fac2 = gradzeta_grady(iy, iz) * gds22(iy, iz) / geo_surf%shat - gradzeta_gradx(iy, iz) * gds21(iy, iz)
               !> the integrand for the parallel flow moment is the parallel velocity
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                  iv = iv_idx(vmu_lo, ivmu)
                  is = is_idx(vmu_lo, ivmu)
                  imu = imu_idx(vmu_lo, ivmu)
                  integrand(ivmu) = fy(iy, ikx, ivmu) * gradpar_zeta(iy, iz) * vpa(iv) &
                       - vperp2(iy, iz, imu) * spec(is)%smz * (f2y(iy, ikx, ivmu) * fac1 + f3y(iy, ikx, ivmu) * fac2)
               end do
               !> integrate over v-space to get the parallel flow, normalised by the reference thermal speed.
               call integrate_vmu_ffs(integrand, upar_wgts, iy, iz, upar(iy, ikx, iz, :))
            end do
         end do
      end do

      deallocate (dens_wgts, upar_wgts, pres_wgts)
      deallocate (f_swap)
      deallocate (integrand)
      deallocate (fy, f2y, f3y)

   end subroutine get_moments_ffs

   !> the Fourier components of the guiding centre distribution function
   !> normalized by the equilibrium Maxwellian is passed in as g,
   !> along with the Fourier components of the electrostatic potential, phi.
   !> g_to_f calculates the Maxwellian-normalized distribution function f,
   !> which is related to g via
   !> f = g + (Ze/T)*(<phi>_R - phi)
   subroutine g_to_f0(g, phi, f)

      use stella_layouts, only: vmu_lo, is_idx
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use gyro_averages, only: gyro_average, j0_ffs

      use stella_transforms, only: transform_ky2y, transform_y2ky
      use kt_grids, only: swap_kxky, swap_kxky_back
      use vpamu_grids, only: maxwell_vpa, maxwell_mu
      use kt_grids, only: naky, naky_all, nakx, ikx_max, ny
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: f
      complex, dimension(:, :), allocatable :: phi_swap
      complex, dimension(:, :), allocatable :: phiy
      complex, dimension(:, :, :), allocatable :: adjust

      complex, dimension (:,:,:,:), allocatable :: g_store
      integer :: ivmu, is, it
      integer :: iz, iv, imu, ia

      allocate (phi_swap(naky_all, ikx_max))
      allocate (phiy(ny, ikx_max))
      allocate (adjust(naky, nakx, -nzgrid:nzgrid))

      allocate(g_store(naky, nakx, -nzgrid:nzgrid, ntubes)) ; g_store = 0.0

      it = 1

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
	      !> compute phi * maxwellian in real space and transform back to k-space
         do iz = -nzgrid, nzgrid
!!!            call swap_kxky(adjust(:, :, iz), phi_swap)
            call swap_kxky(phi(:, :, iz,1), phi_swap)
            call transform_ky2y(phi_swap, phiy(:, :))

            !> phiy = maxwellian * (<phi> - phi) 
            do ia = 1, ny
               phiy(ia, :) = phiy(ia, :) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is)
            end do

            !> Put this back into k-space for future calculations
            call transform_y2ky(phiy(:, :), phi_swap)
            call swap_kxky_back(phi_swap, adjust(:, :, iz))
         end do

         !> calculate the normalized f, given phi and <phi>_R (temporarily stored in f)
         !! f = g + Z/T * (<phi> - phi) * maxwellian 
         is = is_idx(vmu_lo, ivmu)
         f(:, :, :, :, ivmu) = g(:, :, :, :, ivmu) + spec(is)%zt * spread(adjust(:, :, :), 4, ntubes)
         call gyro_average(f(:, :, :, :, ivmu), g_store(:,:,:,:), j0_ffs(:, :, :, ivmu)) 

         f(:, :, :, :, ivmu) = g_store - spec(is)%zt * spread(adjust(:, :, :), 4, ntubes)
      end do

      deallocate(adjust, phi_swap, phiy)
      deallocate(g_store) 
   end subroutine g_to_f0

   subroutine g_to_f1(g, phi, f)

      use stella_layouts, only: vmu_lo, is_idx
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use gyro_averages, only: gyro_average, j0_ffs

      use stella_transforms, only: transform_ky2y, transform_y2ky
      use kt_grids, only: swap_kxky, swap_kxky_back
      use vpamu_grids, only: maxwell_vpa, maxwell_mu
      use kt_grids, only: naky, naky_all, nakx, ikx_max, ny
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: f
      complex, dimension(:, :), allocatable :: phi_swap
      complex, dimension(:, :), allocatable :: phiy
      complex, dimension(:, :, :), allocatable :: adjust

      complex, dimension (:,:,:,:), allocatable :: g_store
      integer :: ivmu, is, it
      integer :: iz, iv, imu, ia

      allocate (phi_swap(naky_all, ikx_max))
      allocate (phiy(ny, ikx_max))
      allocate (adjust(naky, nakx, -nzgrid:nzgrid))

      allocate(g_store(naky, nakx, -nzgrid:nzgrid, ntubes)) ; g_store = 0.0

      it = 1

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         call gyro_average(phi, g_store, j0_ffs(:, :, :, ivmu)) 
         adjust = g_store(:, :, :, 1) - phi(:, :, :, 1)
	      !> compute phi * maxwellian in real space and transform back to k-space
         do iz = -nzgrid, nzgrid
            call swap_kxky(adjust(:, :, iz), phi_swap)
            !! call swap_kxky(phi(:, :, iz,1), phi_swap)
            call transform_ky2y(phi_swap, phiy(:, :))

            !> phiy = maxwellian * phi 
            do ia = 1, ny
               phiy(ia, :) = phiy(ia, :) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is)
            end do

            !> Put this back into k-space for future calculations
            call transform_y2ky(phiy(:, :), phi_swap)
            call swap_kxky_back(phi_swap, adjust(:, :, iz))
         end do
         !> Note that adjust here is (<phi>_R - phi) * maxwellian 
         !> calculate the normalized f, given (<phi>_R - phi) 
         !> f = g + Z/T * (<phi> - phi) * maxwellian 
         f(:, :, :, :, ivmu) = g(:, :, :, :, ivmu) + spec(is)%zt * spread(adjust(:, :, :), 4, ntubes)
         is = is_idx(vmu_lo, ivmu)
         !> f = g + J0 * phi * maxwellian
      end do

      deallocate(adjust, phi_swap, phiy)
      deallocate(g_store) 
   end subroutine g_to_f1


   !==============================================
   !================ GET GVMUS ===================
   !==============================================
   ! get_gvmus takes g(kx,ky,z) and returns average over z of int dxdy g(x,y,z)^2
   ! SHOULD MODIFY TO TAKE ADVANTAGE OF FACT THAT G(KY,KX,Z) LOCAL IS AVAILABLE
   subroutine get_gvmus(g, gv)

      use mp, only: nproc, sum_reduce
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: is_idx, iky_idx, iz_idx
      use zgrid, only: ntubes
      use vpamu_grids, only: nvpa, nmu
      use stella_geometry, only: dl_over_b
      use volume_averages, only: mode_fac

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      real, dimension(:, :, :), intent(out) :: gv

      integer :: ikxkyz, iv, is, imu, iz, iky, ia!, ivp

      ia = 1

      ! when doing volume averages, note the following:
      ! int dxdy g(x,y)^2 = sum_ky |g(ky=0,kx)|^2 + 2 * sum_{kx,ky} |g(ky>0,kx)|^2
      ! factor of 2 accounted for in fac

      gv = 0.
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         is = is_idx(kxkyz_lo, ikxkyz)
         iky = iky_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            do iv = 1, nvpa
               gv(iv, imu, is) = gv(iv, imu, is) + real(g(iv, imu, ikxkyz) * conjg(g(iv, imu, ikxkyz))) * mode_fac(iky) * dl_over_b(ia, iz)
            end do
         end do
      end do
      gv = gv / real(ntubes)

      if (nproc > 1) call sum_reduce(gv, 0)

   end subroutine get_gvmus

   !==============================================
   !================= GET GVZS ===================
   !==============================================
   ! get_gzvs takes g(kx,ky,z,vpa,mu,s) and returns int dmudxdy g(x,y,z,vpa,mu,s)^2
   subroutine get_gzvs(g, gz)

      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: integrate_mu
      use kt_grids, only: nakx, naky
      use volume_averages, only: mode_fac

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(:, :, :, :), intent(out) :: gz

      integer :: ivmu, iz, it, ikx, iky, izp

      real, dimension(:, :, :), allocatable :: gtmp

      allocate (gtmp(-nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      ! when doing volume averages, note the following:
      ! int dxdy g(x,y)^2 = sum_kx |g(kx,ky=0)|^2 + 2 * sum_{kx,ky} |g(kx,ky>0)|^2
      ! factor of 2 accounted for in mode_fac

      gtmp = 0.
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do ikx = 1, nakx
            do iky = 1, naky
               gtmp(:, :, ivmu) = gtmp(:, :, ivmu) + real(g(iky, ikx, :, :, ivmu) * conjg(g(iky, ikx, :, :, ivmu))) * mode_fac(iky)
            end do
         end do
      end do

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            izp = iz + nzgrid + 1
            call integrate_mu(iz, gtmp(iz, it, :), gz(it, izp, :, :))
         end do
      end do

      deallocate (gtmp)

   end subroutine get_gzvs

   !==============================================
   !======== FINISH STELLA DIAGNOSTIC ============
   !==============================================
   subroutine finish_stella_diagnostics(istep)

      use mp, only: proc0
      use redistribute, only: scatter
      use stella_io, only: finish_stella_io
      use stella_time, only: code_dt, code_time
      use stella_save, only: stella_save_for_restart
      use dist_redistribute, only: kxkyz2vmu
      use dist_fn_arrays, only: gnew, gvmu

      implicit none

      integer :: istatus
      integer, intent(in) :: istep

      if (proc0) then
         call write_final_ascii_files
         call write_final_ascii_files_extended
         call close_loop_ascii_files
      end if
      if (save_for_restart) then
         call scatter(kxkyz2vmu, gnew, gvmu)
         call stella_save_for_restart(gvmu, istep, code_time, code_dt, istatus, .true.)
      end if
      call finish_stella_io
      call deallocate_arrays

      nout = 1
      diagnostics_initialized = .false.

   end subroutine finish_stella_diagnostics

   !==============================================
   !========= WRITE LOOP ASCII FILES =============
   !==============================================
   subroutine write_loop_ascii_files(istep, phi2, apar2, bpar2, pflx, vflx, qflx, om, om_avg)

      use stella_time, only: code_time
      use species, only: nspec
      use kt_grids, only: naky, nakx
      use kt_grids, only: aky, akx

      implicit none

      integer, intent(in) :: istep
      real, intent(in) :: phi2, apar2, bpar2
      real, dimension(:), intent(in) :: pflx, vflx, qflx
      complex, dimension(:, :), intent(in) :: om, om_avg

      character(3) :: nspec_str
      character(100) :: str
      integer :: ikx, iky

      write (stdout_unit, '(a7,i7,a6,e12.4,a10,e12.4,a11,e12.4,a12,e12.4)') 'istep=', istep, &
         'time=', code_time, '|phi|^2=', phi2, '|apar|^2= ', apar2, '|bpar|^2= ', bpar2

      call flush (stdout_unit)

      write (nspec_str, '(i3)') 3 * nspec + 1
      str = trim('('//trim(nspec_str)//'es15.4e3)')
      write (fluxes_unit, str) code_time, pflx, vflx, qflx

      call flush (stdout_unit)
      call flush (fluxes_unit)

      if (write_omega .and. istep > 0) then
         do iky = 1, naky
            do ikx = 1, nakx
               write (omega_unit, '(7e16.8)') code_time, aky(iky), akx(ikx), &
                  real(om(iky, ikx)), aimag(om(iky, ikx)), &
                  real(om_avg(iky, ikx)), aimag(om_avg(iky, ikx))
            end do
            if (nakx > 1) write (omega_unit, *)
         end do
         if (naky > 1) write (omega_unit, *)
         call flush (omega_unit)
      end if

   end subroutine write_loop_ascii_files

   !==============================================
   !========= WRITE FINAL ASCII FILES ============
   !==============================================
   subroutine write_final_ascii_files

      use file_utils, only: open_output_file, close_output_file
      use fields_arrays, only: phi, apar, bpar
      use zgrid, only: nzgrid, ntubes
      use zgrid, only: zed
      use kt_grids, only: naky, nakx
      use kt_grids, only: aky, akx, zed0
      use stella_geometry, only: zed_eqarc
      USE dist_fn_arrays, ONLY: kperp2

      implicit none

      integer :: tmpunit
      integer :: iky, ikx, iz, it

      call open_output_file(tmpunit, '.final_fields')
      write (tmpunit, '(10a14)') '# z', 'z-zed0', 'aky', 'akx', &
         'real(phi)', 'imag(phi)', 'real(apar)', 'imag(apar)', &
         'real(bpar)', 'imag(bpar)', 'z_eqarc-zed0', 'kperp2'
      do iky = 1, naky
         do ikx = 1, nakx
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  write (tmpunit, '(12es15.4e3,i3)') zed(iz), zed(iz) - zed0(iky, ikx), aky(iky), akx(ikx), &
                     real(phi(iky, ikx, iz, it)), aimag(phi(iky, ikx, iz, it)), &
                     real(apar(iky, ikx, iz, it)), aimag(apar(iky, ikx, iz, it)), &
                     real(bpar(iky, ikx, iz, it)), aimag(bpar(iky, ikx, iz, it)), &
                     zed_eqarc(iz) - zed0(iky, ikx), kperp2(iky, ikx, it, iz), it
               end do
               write (tmpunit, *)
            end do
         end do
      end do
      call close_output_file(tmpunit)

   end subroutine write_final_ascii_files

   subroutine write_final_ascii_files_extended

      use file_utils, only: open_output_file, close_output_file
      use fields_arrays, only: phi, apar, bpar
      use zgrid, only: nzgrid, ntubes
      use zgrid, only: zed
      use kt_grids, only: naky, nakx
      use kt_grids, only: aky, akx, zed0
      use stella_geometry, only: zed_eqarc
      USE dist_fn_arrays, ONLY: kperp2
      use extended_zgrid, only: map_to_extended_zgrid, map_to_iz_ikx_from_izext
      use extended_zgrid, only: nsegments, nzed_segment
      use extended_zgrid, only: neigen
      implicit none

      integer :: tmpunit

      complex, dimension (:,:), allocatable :: phiext
      integer, dimension(:), allocatable :: iz_from_izext, ikx_from_izext
      integer :: ie, it, iky, ivmu, izext
      integer :: ulim
      integer :: nz_ext
      
      call open_output_file(tmpunit, '.final_fields_extended')
      write (tmpunit, '(10a14)') '# iz','aky', &
           'real(phi)', 'imag(phi)'
                
      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               nz_ext = nsegments(ie, iky) * nzed_segment + 1
               allocate (phiext(naky, nz_ext))
               allocate (iz_from_izext(nz_ext))
               allocate (ikx_from_izext(nz_ext))
               call map_to_extended_zgrid(it, ie, iky, phi(iky, :, :, :), phiext(iky,:), ulim)
               call map_to_iz_ikx_from_izext(iky, ie, iz_from_izext, ikx_from_izext)
               
               do izext = 1, nz_ext
                  write (tmpunit, '(i3,4es15.4e3)') izext , zed(iz_from_izext(izext)) - zed0(iky, ikx_from_izext(izext)), &
                        aky(iky), real(phiext(iky, izext)), aimag(phiext(iky, izext))
               end do
               write (tmpunit, *)
               deallocate (phiext)
               deallocate (iz_from_izext, ikx_from_izext)
             end do
          end do
       end do
       
       call close_output_file(tmpunit)

    end subroutine write_final_ascii_files_extended
   
   !==============================================
   !============ DEALLCOATE ARRAYS ===============
   !==============================================
   subroutine deallocate_arrays

      implicit none

      if (allocated(pflux)) deallocate (pflux)
      if (allocated(qflux)) deallocate (qflux)
      if (allocated(vflux)) deallocate (vflux)
      if (allocated(exchange)) deallocate (exchange)
      if (allocated(pflux_avg)) deallocate (pflux_avg)
      if (allocated(qflux_avg)) deallocate (qflux_avg)
      if (allocated(vflux_avg)) deallocate (vflux_avg)
      if (allocated(heat_avg)) deallocate (heat_avg)
      if (allocated(omega_vs_time)) deallocate (omega_vs_time)

   end subroutine deallocate_arrays

end module stella_diagnostics
