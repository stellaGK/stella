!###############################################################################
!############################ DIAGNOSE DISTRIBUTION ############################
!###############################################################################
! 
! Routines for calculating and writing the distribution.
! 
! The following distribution functions are evolved within stella:
!  - The perturbed distribution function (f)
!  - The perturbed gyroaveraged distribution function (g)
!  - The non-adiabatic part of the perturbed distribution function (h)
! 
!---------------------------------- Input file ---------------------------------
! 
! The diagnostics are writen to text files after every <nwrite> time steps. Moreover,
! they are written to the NetCDF file at every <nwrite>*<nc_mult> time steps. 
! 
! It is important to toggle which dimensions you wish to write, e.g., choose
! <write_g2_vs_vpamus> = True, and the distribution function, <write_distribution_h> = True.
! 
!&diagnostics
!   nwrite = 50.0
!   nc_mult = 1.0
!   navg = 50.0
!   write_all = .false.
!   write_all_spectra_kxkyz = .false.
!   write_all_spectra_kxky = .false.
!   write_all_velocity_space = .false.
!   write_all_distribution = .false.
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
! 
! 
!---------------------------------- Diagnostics --------------------------------
! 
! If <write_g2_vs_vpamus> = .true. the following quantities are calculated:
!   - |g|^2(t, mu, vpa, s)           -->  g2_vs_vpamus in the NetCDF file
!   - |g_nozonal|^2(t, mu, vpa, s)   -->  g2nozonal_vs_vpamus in the NetCDF file
! 
! If <write_g2_vs_zvpas> = .true. the following quantities are calculated:
!   - |g|^2(t, z, vpa, s)            -->  g2_vs_zvpas in the NetCDF file
!   - |g_nozonal|^2(t, z, vpa, s)    -->  g2nozonal_vs_zvpas in the NetCDF file
! 
! If <write_g2_vs_zmus> = .true. the following quantities are calculated:
!   - |g|^2(t, z, mu, s)             -->  g2_vs_zmus in the NetCDF file
!   - |g_nozonal|^2(t, z, mu, s)     -->  g2nozonal_vs_zmus in the NetCDF file
! 
! If <write_g2_vs_kxkyzs> = .true. the following quantities are calculated:
!   - |g|^2(t, kx, ky, z, s)          -->  g2_vs_zkykxs in the NetCDF file
!   - |g_nozonal|^2(t, kx, ky, z, s)  -->  g2nozonal_vs_zkykxs in the NetCDF file
! 
! If <write_g2_vs_zvpamus> = .true. the following quantities are calculated:
!   - |g|^2(t, z, vpa, mu, s)          -->  g2_vs_zvpamus in the NetCDF file
!   - |g_nozonal|^2(t, z, vpa, mu, s)  -->  g2nozonal_vs_zvpamus in the NetCDF file
! 
! If <write_distribution_h> = .true., and depending on the previous flags, we calculate:
!   - |h|^2(t, mu, vpa, s)             -->  h2_vs_vpamus in the NetCDF file
!   - |h|^2(t, z, vpa, s)              -->  h2_vs_zvpas in the NetCDF file
!   - |h|^2(t, z, mu, s)               -->  h2_vs_zmus in the NetCDF file
!   - |h|^2(t, kx, ky, z, s)           -->  h2_vs_zkykxs in the NetCDF file
!   - |h|^2(t, z, vpa, mu, s)          -->  h2_vs_zvpamus in the NetCDF file
!   - |h_nozonal|^2(t, mu, vpa, s)     -->  h2nozonal_vs_vpamus in the NetCDF file
!   - |h_nozonal|^2(t, z, vpa, s)      -->  h2nozonal_vs_zvpas in the NetCDF file
!   - |h_nozonal|^2(t, z, mu, s)       -->  h2nozonal_vs_zmus in the NetCDF file
!   - |h_nozonal|^2(t, kx, ky, z, s)   -->  h2nozonal_vs_zkykxs in the NetCDF file
!   - |h_nozonal|^2(t, z, vpa, mu, s)  -->  h2nozonal_vs_zvpamus in the NetCDF file
! 
! If <write_distribution_f> = .true., and depending on the previous flags, we calculate:
!   - |f|^2(t, mu, vpa, s)             -->  f2_vs_vpamus in the NetCDF file
!   - |f|^2(t, z, vpa, s)              -->  f2_vs_zvpas in the NetCDF file
!   - |f|^2(t, z, mu, s)               -->  f2_vs_zmus in the NetCDF file
!   - |f|^2(t, kx, ky, z, s)           -->  f2_vs_zkykxs in the NetCDF file
!   - |f|^2(t, z, vpa, mu, s)          -->  f2_vs_zvpamus in the NetCDF file
!   - |f_nozonal|^2(t, mu, vpa, s)     -->  f2nozonal_vs_vpamus in the NetCDF file
!   - |f_nozonal|^2(t, z, vpa, s)      -->  f2nozonal_vs_zvpas in the NetCDF file
!   - |f_nozonal|^2(t, z, mu, s)       -->  f2nozonal_vs_zmus in the NetCDF file
!   - |f_nozonal|^2(t, kx, ky, z, s)   -->  f2nozonal_vs_zkykxs in the NetCDF file
!   - |f_nozonal|^2(t, z, vpa, mu, s)  -->  f2nozonal_vs_zvpamus in the NetCDF file
! 
!###############################################################################

module diagnostics_distribution

   ! Load debug flags
   use debug_flags, only: debug => diagnostics_debug

   implicit none
 
   ! Make routines available to other modules
   public :: init_diagnostics_distribution
   public :: write_distribution_to_netcdf_file

   private

contains

!###############################################################################
!############################## WRITE DISTRIBUTION #############################
!###############################################################################

   !============================================================================
   !============== CALCULATE AND WRITE DISTRIBUTION TO NETCDF FILE =============
   !============================================================================
   subroutine write_distribution_to_netcdf_file(nout, timer)

      ! Redistribute data from  i[vpa,mu,s] to i[kx,ky,z,s]
      use redistribute, only: scatter
      use initialise_redistribute, only: kxkyz2vmu
      use job_manage, only: time_message
      use mp, only: proc0

      ! Distribution function
      use arrays_distribution_function, only: gnew
      use arrays_distribution_function, only: gvmu
      
      ! Fields
      use arrays_fields, only: phi, bpar
      use parameters_physics, only: fphi

      ! Dimensions
      use grids_kxky, only: nakx, naky
      use grids_velocity, only: nvpa, nmu
      use grids_z, only: nztot, ntubes
      use grids_species, only: nspec

      ! Write to netcdf file
      use write_diagnostics_to_netcdf, only: write_g2_vs_vpamus_nc
      use write_diagnostics_to_netcdf, only: write_g2_vs_zvpas_nc
      use write_diagnostics_to_netcdf, only: write_g2_vs_zmus_nc
      use write_diagnostics_to_netcdf, only: write_g2_vs_zkykxs_nc
      use write_diagnostics_to_netcdf, only: write_g2_vs_zvpamus_nc
      use write_diagnostics_to_netcdf, only: write_g2nozonal_vs_vpamus_nc
      use write_diagnostics_to_netcdf, only: write_g2nozonal_vs_zvpas_nc
      use write_diagnostics_to_netcdf, only: write_g2nozonal_vs_zmus_nc
      use write_diagnostics_to_netcdf, only: write_g2nozonal_vs_zvpamus_nc
      
      use write_diagnostics_to_netcdf, only: write_h2_vs_vpamus_nc
      use write_diagnostics_to_netcdf, only: write_h2_vs_zvpas_nc
      use write_diagnostics_to_netcdf, only: write_h2_vs_zmus_nc
      use write_diagnostics_to_netcdf, only: write_h2_vs_zkykxs_nc
      use write_diagnostics_to_netcdf, only: write_h2_vs_zvpamus_nc
      use write_diagnostics_to_netcdf, only: write_h2nozonal_vs_vpamus_nc
      use write_diagnostics_to_netcdf, only: write_h2nozonal_vs_zvpas_nc
      use write_diagnostics_to_netcdf, only: write_h2nozonal_vs_zmus_nc
      use write_diagnostics_to_netcdf, only: write_h2nozonal_vs_zvpamus_nc
      use write_diagnostics_to_netcdf, only: write_f2_vs_vpamus_nc
      use write_diagnostics_to_netcdf, only: write_f2_vs_zvpas_nc
      use write_diagnostics_to_netcdf, only: write_f2_vs_zmus_nc
      use write_diagnostics_to_netcdf, only: write_f2_vs_zkykxs_nc
      use write_diagnostics_to_netcdf, only: write_f2_vs_zvpamus_nc
      use write_diagnostics_to_netcdf, only: write_f2nozonal_vs_vpamus_nc
      use write_diagnostics_to_netcdf, only: write_f2nozonal_vs_zvpas_nc
      use write_diagnostics_to_netcdf, only: write_f2nozonal_vs_zmus_nc
      use write_diagnostics_to_netcdf, only: write_f2nozonal_vs_zvpamus_nc

      use write_diagnostics_to_netcdf, only: write_g_vs_zvpas_zonal_nc
      use write_diagnostics_to_netcdf, only: write_free_energy_nc
      
      ! Input file
      use parameters_diagnostics, only: write_distribution_g
      use parameters_diagnostics, only: write_distribution_h
      use parameters_diagnostics, only: write_distribution_f
      use parameters_diagnostics, only: write_g2_vs_vpamus
      use parameters_diagnostics, only: write_g2_vs_zvpas
      use parameters_diagnostics, only: write_g2_vs_zmus
      use parameters_diagnostics, only: write_g2_vs_kxkyzs
      use parameters_diagnostics, only: write_g2_vs_zvpamus

      use parameters_diagnostics, only: write_g_vs_zvpas_zonal
      use parameters_diagnostics, only: write_free_energy
      use parameters_diagnostics, only: number_zonals_kxs, zonal_iks
      
      ! Calculations
      use calculations_tofrom_ghf, only: g_to_f, g_to_h

      implicit none

      ! The pointer in the netcdf file and a timer
      real, dimension(:), intent(in out) :: timer
      integer, intent(in) :: nout

      ! We will write g(vpa,mu,s), g(tube,z,vpa,s), g(tube,z,mu,s)
      ! And we will do this for g, h and delta f
      real, dimension(:, :, :), allocatable :: g2_vs_vpamus, g2nozonal_vs_vpamus
      real, dimension(:, :, :, :), allocatable :: g2_vs_zvpas, g2_vs_zmus, g2nozonal_vs_zvpas, g2nozonal_vs_zmus

      ! Zonals
      real, dimension(:, :, :, :, :), allocatable :: g_vs_zvpas_zonal
      real, dimension(:), allocatable :: free_energy_vs_kx

      real, dimension(:, :, :, :, :), allocatable :: g2_vs_zkykxs, g2_vs_zvpamus, g2nozonal_vs_zvpamus

      !------------------------------------------------------------------------- 

      ! Only continue if we write data 
      if ((.not. write_distribution_g) .and. (.not. write_distribution_h) .and. (.not. write_distribution_f)) return
      if ((.not. write_g2_vs_vpamus) .and. (.not. write_g2_vs_zvpas) .and. (.not. write_g2_vs_zmus) &
            .and. (.not. write_g2_vs_kxkyzs) .and. (.not. write_g2_vs_zvpamus) .and. (.not. write_g_vs_zvpas_zonal) &
            .and. (.not. write_free_energy)) return

      ! Start timer
      if (proc0) call time_message(.false., timer(:), 'Write distribution')
      if (debug) write (*, *) 'diagnostics::diagnostics_distribution::write_distribution_to_netcdf_file'

      ! Allocate arrays
      if (write_g2_vs_vpamus) allocate (g2_vs_vpamus(nvpa, nmu, nspec))
      if (write_g2_vs_vpamus) allocate (g2nozonal_vs_vpamus(nvpa, nmu, nspec))
      if (write_g2_vs_zvpas) allocate (g2_vs_zvpas(ntubes, nztot, nvpa, nspec))
      if (write_g2_vs_zvpas) allocate (g2nozonal_vs_zvpas(ntubes, nztot, nvpa, nspec))
      if (write_g2_vs_zmus) allocate (g2_vs_zmus(ntubes, nztot, nmu, nspec))
      if (write_g2_vs_zmus) allocate (g2nozonal_vs_zmus(ntubes, nztot, nmu, nspec))
      if (write_g2_vs_zvpamus) allocate (g2_vs_zvpamus(nztot, ntubes, nvpa, nmu, nspec))
      if (write_g2_vs_zvpamus) allocate (g2nozonal_vs_zvpamus(nztot, ntubes, nvpa, nmu, nspec))
      if (write_g2_vs_kxkyzs) allocate (g2_vs_zkykxs(ntubes, nztot, naky, nakx, nspec))

      ! Zonal diagnostics
      if (write_g_vs_zvpas_zonal) allocate (g_vs_zvpas_zonal(2, number_zonals_kxs, nztot, nvpa, nspec))
      if (write_free_energy) allocate (free_energy_vs_kx(nakx)) 

      ! Redistribute the data from <gnew>(ky,kx,z,tube,i[vpa,mu,s]) to <gvmu>(vpa,mu,i[kx,ky,z,s]) for |g|^2(z,kx,ky,s)
      if (write_g2_vs_kxkyzs) call scatter(kxkyz2vmu, gnew, gvmu)

      ! Write distribution data for g
      if (debug) write (*, *) 'diagnostics::diagnostics_distribution::write_distribution_g'
      if (write_distribution_g) then

         ! Use gnew(ky, kx, z, tube, ivmus) to calculate |g|^2(z, vpa, s), |g|^2(z, mu, s) and |g|^2(vpa, mu, s)
         call calculate_distribution(gnew, gvmu, g2_vs_zmus, g2_vs_zvpas, g2_vs_vpamus, g2_vs_zkykxs, g2_vs_zvpamus, &
                  g2nozonal_vs_zmus, g2nozonal_vs_zvpas, g2nozonal_vs_vpamus, g2nozonal_vs_zvpamus, g_vs_zvpas_zonal) 

         ! Write the distribution data to the netcdf file
         if (write_g2_vs_vpamus .and. proc0) call write_g2_vs_vpamus_nc(nout, g2_vs_vpamus)
         if (write_g2_vs_zvpas .and. proc0) call write_g2_vs_zvpas_nc(nout, g2_vs_zvpas)
         if (write_g2_vs_zmus .and. proc0) call write_g2_vs_zmus_nc(nout, g2_vs_zmus)
         if (write_g2_vs_kxkyzs .and. proc0) call write_g2_vs_zkykxs_nc(nout, g2_vs_zkykxs)
         if (write_g2_vs_zvpamus .and. proc0) call write_g2_vs_zvpamus_nc(nout, g2_vs_zvpamus)
         if (write_g2_vs_vpamus .and. proc0) call write_g2nozonal_vs_vpamus_nc(nout, g2nozonal_vs_vpamus)
         if (write_g2_vs_zvpas .and. proc0) call write_g2nozonal_vs_zvpas_nc(nout, g2nozonal_vs_zvpas)
         if (write_g2_vs_zmus .and. proc0) call write_g2nozonal_vs_zmus_nc(nout, g2nozonal_vs_zmus)
         if (write_g2_vs_zvpamus .and. proc0) call write_g2nozonal_vs_zvpamus_nc(nout, g2nozonal_vs_zvpamus)
         if (write_g_vs_zvpas_zonal .and. proc0) call write_g_vs_zvpas_zonal_nc(nout, g_vs_zvpas_zonal)
         if (write_free_energy) then 
            call g_to_h(gnew, phi, bpar, fphi)
            call calculate_free_energy(gnew, free_energy_vs_kx)
            if (proc0) call write_free_energy_nc(nout, free_energy_vs_kx)
            call g_to_h(gnew, phi, bpar, -fphi)
         end if

      end if

      ! Write distribution data for h 
      if (debug) write (*, *) 'diagnostics::diagnostics_distribution::write_distribution_h'
      if (write_distribution_h) then

         ! Switch to h
         if (write_g2_vs_kxkyzs) call g_to_h(gvmu, phi, bpar, fphi)
         call g_to_h(gnew, phi, bpar, fphi)

         ! Use gnew(ky, kx, z, tube, ivmus) to calculate |g|^2(z, vpa, s), |g|^2(z, mu, s) and |g|^2(vpa, mu, s)
         call calculate_distribution(gnew, gvmu, g2_vs_zmus, g2_vs_zvpas, g2_vs_vpamus, g2_vs_zkykxs, g2_vs_zvpamus, &
                  g2nozonal_vs_zmus, g2nozonal_vs_zvpas, g2nozonal_vs_vpamus, g2nozonal_vs_zvpamus, g_vs_zvpas_zonal)
 
         ! Write the distribution data to the netcdf file
         if (write_g2_vs_vpamus .and. proc0) call write_h2_vs_vpamus_nc(nout, g2_vs_vpamus)
         if (write_g2_vs_zvpas .and. proc0) call write_h2_vs_zvpas_nc(nout, g2_vs_zvpas)
         if (write_g2_vs_zmus .and. proc0) call write_h2_vs_zmus_nc(nout, g2_vs_zmus)
         if (write_g2_vs_kxkyzs .and. proc0) call write_h2_vs_zkykxs_nc(nout, g2_vs_zkykxs)
         if (write_g2_vs_zvpamus .and. proc0) call write_h2_vs_zvpamus_nc(nout, g2_vs_zvpamus)
         if (write_g2_vs_vpamus .and. proc0) call write_h2nozonal_vs_vpamus_nc(nout, g2nozonal_vs_vpamus)
         if (write_g2_vs_zvpas .and. proc0) call write_h2nozonal_vs_zvpas_nc(nout, g2nozonal_vs_zvpas)
         if (write_g2_vs_zmus .and. proc0) call write_h2nozonal_vs_zmus_nc(nout, g2nozonal_vs_zmus)
         if (write_g2_vs_zvpamus .and. proc0) call write_h2nozonal_vs_zvpamus_nc(nout, g2nozonal_vs_zvpamus)

         call calculate_free_energy(gnew, free_energy_vs_kx)
         if (proc0) call write_free_energy_nc(nout, free_energy_vs_kx)

         ! Switch back to g
         if (write_g2_vs_kxkyzs) call g_to_h(gvmu, phi, bpar, -fphi)
         call g_to_h(gnew, phi, bpar, -fphi)

      end if

      ! Write distribution data for g
      if (debug) write (*, *) 'diagnostics::diagnostics_distribution::write_distribution_f'
      if (write_distribution_f) then

         ! Switch to f
         if (write_g2_vs_kxkyzs) call g_to_h(gvmu, phi, bpar, fphi)
         call g_to_f(gnew, phi, fphi)

         ! Use gnew(ky, kx, z, tube, ivmus) to calculate |g|^2(z, vpa, s), |g|^2(z, mu, s) and |g|^2(vpa, mu, s)
         call calculate_distribution(gnew, gvmu, g2_vs_zmus, g2_vs_zvpas, g2_vs_vpamus, g2_vs_zkykxs, g2_vs_zvpamus, &
                  g2nozonal_vs_zmus, g2nozonal_vs_zvpas, g2nozonal_vs_vpamus, g2nozonal_vs_zvpamus, g_vs_zvpas_zonal)

         ! Write the distribution data to the netcdf file
         if (write_g2_vs_vpamus .and. proc0) call write_f2_vs_vpamus_nc(nout, g2_vs_vpamus)
         if (write_g2_vs_zvpas .and. proc0) call write_f2_vs_zvpas_nc(nout, g2_vs_zvpas)
         if (write_g2_vs_zmus .and. proc0) call write_f2_vs_zmus_nc(nout, g2_vs_zmus)
         if (write_g2_vs_kxkyzs .and. proc0) call write_f2_vs_zkykxs_nc(nout, g2_vs_zkykxs)
         if (write_g2_vs_zvpamus .and. proc0) call write_f2_vs_zvpamus_nc(nout, g2_vs_zvpamus)
         if (write_g2_vs_vpamus .and. proc0) call write_f2nozonal_vs_vpamus_nc(nout, g2nozonal_vs_vpamus)
         if (write_g2_vs_zvpas .and. proc0) call write_f2nozonal_vs_zvpas_nc(nout, g2nozonal_vs_zvpas)
         if (write_g2_vs_zmus .and. proc0) call write_f2nozonal_vs_zmus_nc(nout, g2nozonal_vs_zmus)
         if (write_g2_vs_zvpamus .and. proc0) call write_f2nozonal_vs_zvpamus_nc(nout, g2nozonal_vs_zvpamus)
         ! Switch back to g
         if (write_g2_vs_kxkyzs) call g_to_h(gvmu, phi, bpar, -fphi)
         call g_to_f(gnew, phi, -fphi)

      end if

      ! Deallocate arrays
      if (write_g2_vs_vpamus) deallocate (g2_vs_vpamus)
      if (write_g2_vs_zvpas) deallocate (g2_vs_zvpas)
      if (write_g2_vs_zmus) deallocate (g2_vs_zmus)
      if (write_g2_vs_kxkyzs) deallocate (g2_vs_zkykxs)
      if (write_g2_vs_zvpamus) deallocate (g2_vs_zvpamus)
      if (write_g_vs_zvpas_zonal) deallocate (g_vs_zvpas_zonal)
      if (write_free_energy) deallocate (free_energy_vs_kx)
      ! End timer
      if (proc0) call time_message(.false., timer(:), 'Write distribution')

   end subroutine write_distribution_to_netcdf_file

   !============================================================================
   !======================== Calculate |g|^2(z,mu,s) ===========================
   !============================================================================
   ! Calculate int dxdydmu  |g(x,y,z,mu,vpa,s)|^2 = |g(vpa,z,s)|^2
   ! We use Parsevals theorem: int dxdy |g(x,y)|^2 = sum_{kx,ky} |g(kx,ky)|^2 
   ! And the mirror condition: int dxdy |g(x,y)|^2 = sum_ky |g(ky=0,kx)|^2 + 2 * sum_{kx,ky} |g(ky>0,kx)|^2
   !============================================================================
   subroutine calculate_distribution(g_vs_kykxztube, g_vs_vpamuikxkyzs, g2_vs_tzmus, &
      g2_vs_zvpas, g2_vs_vpamus, g2_vs_zkykxs, g2_vs_zvpamus, &
      g2nozonal_vs_tzmus, g2nozonal_vs_zvpas, g2nozonal_vs_vpamus, g2nozonal_vs_zvpamus, &
      g_vs_zvpas_zonal)

      ! Geometry
      use geometry, only: dl_over_b

      ! Dimensions
      use grids_kxky, only: zonal_mode
      use calculations_volume_averages, only: mode_fac
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: nvpa, nmu
      use grids_kxky, only: nakx, naky, akx
      use grids_velocity, only: maxwell_vpa, maxwell_mu
      use grids_species, only: nspec   

      use arrays_fields, only: phi
      use grids_species, only: spec
      use grids_species, only: tite, nine
      use grids_species, only: ion_species
      use grids_velocity, only: wgts_mu, wgts_vpa

      ! Calculations
      use calculations_velocity_integrals, only: integrate_vpa, integrate_mu, integrate_vmu

      ! Routines
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx, kxkyz_lo
      use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use mp, only: nproc, sum_reduce
      
      ! Input file
      use parameters_diagnostics, only: write_g2_vs_vpamus
      use parameters_diagnostics, only: write_g2_vs_zvpas
      use parameters_diagnostics, only: write_g2_vs_zmus
      use parameters_diagnostics, only: write_g2_vs_kxkyzs
      use parameters_diagnostics, only: write_g2_vs_zvpamus

      use parameters_diagnostics, only: write_g_vs_zvpas_zonal
      use parameters_diagnostics, only: number_zonals_kxs, zonal_iks

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g_vs_kykxztube
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g_vs_vpamuikxkyzs
      real, dimension(:, :, :, :, :), intent(out) :: g2_vs_zkykxs, g2_vs_zvpamus, g2nozonal_vs_zvpamus
      real, dimension(:, :, :, :), intent(out) :: g2_vs_tzmus, g2_vs_zvpas, g2nozonal_vs_tzmus, g2nozonal_vs_zvpas
      real, dimension(:, :, :), intent(out) :: g2_vs_vpamus, g2nozonal_vs_vpamus
      
      ! Zonals
      real, dimension(:, :,  :, :, :), intent(out) :: g_vs_zvpas_zonal

      ! Arrays needed to perform calculations
      real, dimension(:, :, :), allocatable :: g2_vs_ztubeivmus, g2nozonal_vs_ztubeivmus 
      real, dimension(:, :), allocatable :: g2_vs_ztube, g2_vs_vpamu

      ! Zonals
      real, dimension(:, :, :, :), allocatable :: g_vs_zivmus_zonal

      integer :: ivmus, ia, iz, it, ikx, iky, izp, is, iv, imu, ikxkyzs
      integer :: j, ikx_plot
      !---------------------------------------------------------------------- 

      ! Allocate local arrays
      if (debug) write (*, *) 'diagnostics::diagnostics_distribution::calculate_distribution'
      if (write_g2_vs_zvpas .or. write_g2_vs_zmus) then
         allocate (g2_vs_ztubeivmus(-nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); g2_vs_ztubeivmus = 0.
         allocate (g2nozonal_vs_ztubeivmus(-nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); g2nozonal_vs_ztubeivmus = 0.
      end if
      if (write_g2_vs_kxkyzs) then; allocate (g2_vs_vpamu(nvpa, nmu)); g2_vs_vpamu = 0.; end if
         allocate (g2_vs_ztube(-nzgrid:nzgrid, ntubes)); g2_vs_ztube = 0.

      if (write_g_vs_zvpas_zonal) allocate (g_vs_zivmus_zonal(2, number_zonals_kxs, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); g_vs_zivmus_zonal = 0.0

      ! We add contributions to <g2_vs_tzmus, g2_vs_zvpas, g2_vs_vpamus> so make sure they are zero first 
      if (write_g2_vs_zvpas) then; g2_vs_zvpas = 0.; g2nozonal_vs_zvpas = 0.; end if
      if (write_g2_vs_zmus) then; g2_vs_tzmus = 0.; g2nozonal_vs_tzmus = 0.; end if
      if (write_g2_vs_vpamus) then; g2_vs_vpamus = 0.; g2nozonal_vs_vpamus = 0.; end if
      if (write_g2_vs_zvpamus) then; g2_vs_zvpamus = 0.; g2nozonal_vs_zvpamus = 0.; end if
      if (write_g2_vs_kxkyzs) then; g2_vs_zkykxs = 0.; end if
      ! Assume we only have one flux tube (assume <radial_variation> = False)
      ia = 1

      ! Construct g(z, tube, i[vpa, mu, s]) by using Parseval's theorem to perform sum_{kx,ky} |g^2(kx,ky)| 
      ! int dxdy |g(x,y)|^2 = sum_kx |g(kx,ky=0)|^2 + 2 * sum_{kx,ky} |g(kx,ky>0)|^2 (factor of 2 is accounted for in mode_fac)
      if (debug) write (*, *) 'diagnostics::diagnostics_distribution::calculate_distribution::start_calculations' 
      do ivmus = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmus)
         imu = imu_idx(vmu_lo, ivmus)
         is = is_idx(vmu_lo, ivmus)
         do ikx = 1, nakx
            do iky = 1, naky

               ! Calculate int dxdy |g(x,y)|^2 for each z-point
               g2_vs_ztube = real(g_vs_kykxztube(iky, ikx, :, :, ivmus) * conjg(g_vs_kykxztube(iky, ikx, :, :, ivmus))) * mode_fac(iky)

               if (write_g2_vs_zvpamus) g2_vs_zvpamus(:, :, iv, imu, is) = g2_vs_zvpamus(:, :, iv, imu, is) + g2_vs_ztube(:, :)
               if (write_g2_vs_zvpas .or. write_g2_vs_zmus) g2_vs_ztubeivmus(:, :, ivmus) = g2_vs_ztubeivmus(:, :, ivmus) + g2_vs_ztube(:, :)
               
               if (.not. zonal_mode(iky)) then
                  if (write_g2_vs_zvpamus) g2nozonal_vs_zvpamus(:, :, iv, imu, is) = g2nozonal_vs_zvpamus(:, :, iv, imu, is) + g2_vs_ztube(:, :)
                  if (write_g2_vs_zvpas .or. write_g2_vs_zmus) g2nozonal_vs_ztubeivmus(:, :, ivmus) = g2nozonal_vs_ztubeivmus(:, :, ivmus) + g2_vs_ztube(:, :)
               end if

               ! Now take the field line average
               if (write_g2_vs_vpamus) then
                  do it = 1, ntubes
                     g2_vs_vpamus(iv, imu, is) = g2_vs_vpamus(iv, imu, is) + sum(g2_vs_ztube(:, it) * dl_over_b(ia, :))
                     if (.not. zonal_mode(iky)) then
                        g2nozonal_vs_vpamus(iv, imu, is) = g2nozonal_vs_vpamus(iv, imu, is) + sum(g2_vs_ztube(:, it) * dl_over_b(ia, :))
                     end if
                  end do
               end if
               
            end do
         end do

         ! Zonals
         if (write_g_vs_zvpas_zonal) then
            if (number_zonals_kxs == 1) then 
               if (abs(akx(1)) < epsilon(0.)) then
                  ikx_plot = 2
               else
                  ikx_plot = 1
               end if
               g_vs_zivmus_zonal(1, 1, :, ivmus) = real(g_vs_kykxztube(1, ikx_plot, :, 1, ivmus))
               g_vs_zivmus_zonal(2, 1, :, ivmus) = aimag(g_vs_kykxztube(1, ikx_plot, :, 1, ivmus))

            else 

               ikx_plot = 0

               do ikx = 1, nakx
                  if (any(zonal_iks == ikx)) then
                     ikx_plot = ikx_plot + 1

                     g_vs_zivmus_zonal(1, ikx_plot, :, ivmus) = real(g_vs_kykxztube(1, ikx, :, 1, ivmus))
                     g_vs_zivmus_zonal(2, ikx_plot, :, ivmus) = aimag(g_vs_kykxztube(1, ikx, :, 1, ivmus))

                  end if
               end do
            end if 
         end if 

      end do

      ! Perform the velocity integration for each (kx,ky,z)-point
      if (debug) write (*, *) 'diagnostics::diagnostics_distribution::calculate_distribution::write_g2_vs_kxkyzs'
      if (write_g2_vs_kxkyzs) then
         do ikxkyzs = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc

            ! Get the (kx,ky,z,tube,s) indices on this processor
            iky = iky_idx(kxkyz_lo, ikxkyzs)
            ikx = ikx_idx(kxkyz_lo, ikxkyzs)
            iz = iz_idx(kxkyz_lo, ikxkyzs)
            it = it_idx(kxkyz_lo, ikxkyzs)
            is = is_idx(kxkyz_lo, ikxkyzs)
            izp = iz + nzgrid + 1

            ! Velocity integration allocate (g2_vs_zkykxs(ntubes, nztot, nakx, naky, nspec))
            g2_vs_vpamu = real(g_vs_vpamuikxkyzs(:, :, ikxkyzs) * conjg(g_vs_vpamuikxkyzs(:, :, ikxkyzs))) * mode_fac(iky)
            call integrate_vmu(g2_vs_vpamu, iz, g2_vs_zkykxs(it, izp, iky, ikx, is))

         end do
      end if
 
      ! Next integrate over the perpendicular/parallel velocity vpa for each z-point
      ! The velocity perpendicular integration takes |g|^2(ivmus) and returns int dmu |g|^2(ivmus) = |g|^2(vpa,s)
      ! The velocity parallel integration takes |g|^2(ivmus) and returns int dvpa |g|^2(ivmus) = |g|^2(mu,s)
      ! There is a sum_reduce() inside of the velocity integration, so <g2_vs_tzmus> holds the total sum
      if (debug) write (*, *) 'diagnostics::diagnostics_distribution::calculate_distribution::write_g2_vs_zvpas .or. write_g2_vs_zmus'
      if (write_g2_vs_zvpas .or. write_g2_vs_zmus) then
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               izp = iz + nzgrid + 1
               if (write_g2_vs_zvpas) call integrate_mu(iz, g2_vs_ztubeivmus(iz, it, :), g2_vs_zvpas(it, izp, :, :))
               do j = 1, 2
                  do ikx = 1, number_zonals_kxs
                     if (write_g_vs_zvpas_zonal) call integrate_mu(iz, g_vs_zivmus_zonal(j, ikx, iz, :), g_vs_zvpas_zonal(j, ikx, izp, :, :))
                  end do 
               end do
               if (write_g2_vs_zmus) call integrate_vpa(g2_vs_ztubeivmus(iz, it, :), g2_vs_tzmus(it, izp, :, :))
               if (write_g2_vs_zvpas) call integrate_mu(iz, g2nozonal_vs_ztubeivmus(iz, it, :), g2nozonal_vs_zvpas(it, izp, :, :))
               if (write_g2_vs_zmus) call integrate_vpa(g2nozonal_vs_ztubeivmus(iz, it, :), g2nozonal_vs_tzmus(it, izp, :, :))
            end do
         end do
      end if

      ! For the field line average normalise to account for contributions from multiple flux tubes
      ! in a flux tube train, and sum the values on all processors and to send them to <proc0>
      if (debug) write (*, *) 'diagnostics::diagnostics_distribution::calculate_distribution::sum_all_reduce'
      if (write_g2_vs_vpamus) then
         g2_vs_vpamus = g2_vs_vpamus / real(ntubes)
         if (nproc > 1) call sum_reduce(g2_vs_vpamus, 0)
         if (nproc > 1) call sum_reduce(g2nozonal_vs_vpamus, 0)
      end if
      if (write_g2_vs_zvpamus) then
         if (nproc > 1) call sum_reduce(g2_vs_zvpamus, 0)
         if (nproc > 1) call sum_reduce(g2nozonal_vs_zvpamus, 0)
      end if
      if (write_g2_vs_kxkyzs) then
         if (nproc > 1) call sum_reduce(g2_vs_zkykxs, 0)
      end if

      ! Deallocate local arrays
      if (write_g2_vs_zvpas) deallocate (g2nozonal_vs_ztubeivmus)
      if (write_g2_vs_zvpas) deallocate (g2_vs_ztubeivmus)
      if (write_g2_vs_kxkyzs) deallocate (g2_vs_vpamu)
      if (write_g_vs_zvpas_zonal) deallocate (g_vs_zivmus_zonal)
      deallocate (g2_vs_ztube)

   end subroutine calculate_distribution

   subroutine calculate_free_energy(g_vs_kykxztube, free_energy_vs_kx)

      ! Geometry
      use geometry, only: dl_over_b

      ! Dimensions
      use grids_kxky, only: zonal_mode
      use calculations_volume_averages, only: mode_fac
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: nvpa, nmu
      use grids_kxky, only: nakx, naky, akx
      use grids_velocity, only: maxwell_vpa, maxwell_mu
      use grids_species, only: nspec 

      use arrays_fields, only: phi
      use grids_species, only: spec
      use grids_species, only: tite, nine
      use grids_species, only: ion_species
      use grids_velocity, only: wgts_mu, wgts_vpa

      ! Calculations
      use calculations_velocity_integrals, only: integrate_vmu

      ! Routines
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx, kxkyz_lo
      use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use mp, only: nproc, sum_reduce

      ! use parameters_diagnostics, only: write_free_energy

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g_vs_kykxztube
      ! Zonals
      real, dimension(:), intent(out) :: free_energy_vs_kx

      real, dimension(:, :), allocatable :: g2_vs_ztube

      ! Arrays needed to perform calculations
      real, dimension(:, :, :, :), allocatable :: free_energy_vs_kxkyz_vmu
      real, dimension(:, :, :, :), allocatable ::  free_energy_vs_kxkyz

      integer :: ivmus, ia, iz, it, ikx, iky, is, iv, imu, ikxkyzs
      !---------------------------------------------------------------------- 

      ! Allocate local arrays
      allocate (g2_vs_ztube(-nzgrid:nzgrid, ntubes)); g2_vs_ztube = 0.
      allocate (free_energy_vs_kxkyz_vmu(naky, nakx, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); free_energy_vs_kxkyz_vmu = 0.0
      allocate (free_energy_vs_kxkyz(naky, nakx, -nzgrid:nzgrid, nspec)); free_energy_vs_kxkyz = 0.0


      ! Assume we only have one flux tube (assume <radial_variation> = False)
      ia = 1
      ! Construct g(z, tube, i[vpa, mu, s]) by using Parseval's theorem to perform sum_{kx,ky} |g^2(kx,ky)| 
      ! int dxdy |g(x,y)|^2 = sum_kx |g(kx,ky=0)|^2 + 2 * sum_{kx,ky} |g(kx,ky>0)|^2 (factor of 2 is accounted for in mode_fac)
      if (debug) write (*, *) 'diagnostics::diagnostics_distribution::calculate_distribution::start_calculations' 
      do ivmus = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmus)
         imu = imu_idx(vmu_lo, ivmus)
         is = is_idx(vmu_lo, ivmus)
         do ikx = 1, nakx
            do iky = 1, naky
               ! Calculate int dxdy |g(x,y)|^2 for each z-point
               g2_vs_ztube = real(g_vs_kykxztube(iky, ikx, :, :, ivmus) * conjg(g_vs_kykxztube(iky, ikx, :, :, ivmus))) 
                  
               it = 1
               free_energy_vs_kxkyz_vmu(iky, ikx, :, ivmus) = dl_over_b(ia, :) * & 
                  spec(ion_species)%temp * g2_vs_ztube(:, it) / (maxwell_vpa(iv, ion_species) * maxwell_mu(ia, :, imu, ion_species) )! &
                  ! + spec(ion_species)%dens * spec(ion_species)%temp * (conjg(phi(iky, ikx, :, it)) * mode_fac(iky) &
                  ! + tite / nine * (conjg(phi(iky, ikx, :, it)) * mode_fac(iky) - conjg(phi(0, ikx, :, it)) * mode_fac(0)) ) &
                  ! * real(phi(iky, ikx, :, it)) * conjg(phi(iky, ikx, :, it)) * mode_fac(iky))
            end do 
         end do
      end do

      it = 1 
      do iz = -nzgrid, nzgrid
         do ikx = 1, nakx
            do iky = 1, naky
               call integrate_vmu(free_energy_vs_kxkyz_vmu(iky,ikx,iz,:), ia, iz, free_energy_vs_kxkyz(iky,ikx,iz,:) )
               free_energy_vs_kxkyz(iky, ikx, iz, ion_species) = free_energy_vs_kxkyz(iky, ikx, iz, ion_species) &
                     + spec(ion_species)%dens * spec(ion_species)%temp * (conjg(phi(iky, ikx, iz, it)) * mode_fac(iky) &
                     + tite / nine * (conjg(phi(iky, ikx, iz, it)) * mode_fac(iky) - conjg(phi(0, ikx, iz, it)) * mode_fac(0)) ) &
                     * real(phi(iky, ikx, iz, it)) * conjg(phi(iky, ikx, iz, it)) * mode_fac(iky) * dl_over_b(ia, iz) 
            end do
         end do 
      end do 

      free_energy_vs_kx(:) = sum(free_energy_vs_kxkyz(0,:,:,ion_species), dim=2)/ sum(free_energy_vs_kxkyz(:,:,:,ion_species))

      deallocate (free_energy_vs_kxkyz, free_energy_vs_kxkyz_vmu)
      deallocate (g2_vs_ztube)

   end subroutine calculate_free_energy
!###############################################################################
!############################ INITALIZE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !======================== INITALIZE THE DIAGNOSTICS =========================
   !============================================================================
   subroutine init_diagnostics_distribution()

      implicit none

   end subroutine init_diagnostics_distribution

end module diagnostics_distribution

