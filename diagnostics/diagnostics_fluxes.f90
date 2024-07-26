
!###############################################################################
!############################### DIAGNOSE FLUXES ###############################
!###############################################################################
! 
! Routines for calculating and writing the turbulent fluxes. 
! 
! The particle flux is denoted by pflux.
! The momentum flux is denoted by vflux.
! The heat flux is denoted by qflux.
! 
!###############################################################################
 
module diagnostics_fluxes

   implicit none
  
   public :: init_diagnostics_fluxes
   public :: finish_diagnostics_fluxes
   public :: write_fluxes_to_ascii_file 
   public :: write_fluxes_to_netcdf_file 

   private 

   ! The <units> are used to identify the external ascii files
   integer :: fluxes_unit 

   ! When writing the netcdf data, remember the fluxes versus time for the ascii file
   real, dimension(:), allocatable :: pflux_vs_s, qflux_vs_s, vflux_vs_s  

   ! Debugging
   logical :: debug = .false.

contains

!###############################################################################
!################################# WRITE FLUXES ################################
!###############################################################################

   !============================================================================
   !================= CALCULATE AND WRITE FLUXES TO NETCDF FILE ================
   !============================================================================
   subroutine write_fluxes_to_netcdf_file(nout, timer, write_to_netcdf_file)

      ! Dimensions
      use kt_grids, only: naky, nakx
      use zgrid, only: nztot, ntubes
      use species, only: nspec

      ! Flags 
      use physics_flags, only: radial_variation
      use physics_flags, only: full_flux_surface
      use parameters_diagnostics, only: write_fluxes_vs_time
      use parameters_diagnostics, only: write_fluxes_kxkyz
      use parameters_diagnostics, only: write_radial_fluxes 

      ! Routines
      use job_manage, only: time_message
      use mp, only: proc0
      
      ! Write to netCDF file 
      use stella_io, only: write_fluxes_kxkyz_nc
      use stella_io, only: write_fluxes_vs_time_nc

      implicit none 

      integer, intent(in) :: nout   ! The pointer in the netcdf file
      logical, intent(in) :: write_to_netcdf_file    
      real, dimension(:), intent(in out) :: timer    

      ! We want to write flux(ky,kx,z,tube,s) to the netcdf file
      real, dimension(:, :, :, :, :), allocatable :: pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts

      !---------------------------------------------------------------------- 

      ! Start timer
      if (proc0) call time_message(.false., timer(:), 'Write fluxes')

      ! Allocate the arrays that we want to write to the netcdf file 
      if (write_fluxes_kxkyz) then; allocate (pflux_kxkyzts(naky, nakx, nztot, ntubes, nspec)); pflux_kxkyzts = 0.0; end if
      if (write_fluxes_kxkyz) then; allocate (vflux_kxkyzts(naky, nakx, nztot, ntubes, nspec)); vflux_kxkyzts = 0.0; end if
      if (write_fluxes_kxkyz) then; allocate (qflux_kxkyzts(naky, nakx, nztot, ntubes, nspec)); qflux_kxkyzts = 0.0; end if
      
      !**********************************************************************
      !                          WRITE TO TXT FILE                          !
      !**********************************************************************
      ! Note that to obtain <pflux_vs_s>, <vflux_vs_s>, <qflux_vs_s> to write
      ! the fluxes to the txt files, we already calculate fluxes(kx,ky,z,s)
      !**********************************************************************

      ! Calculate the fluxes if <radial_variation> = True
      if (radial_variation .or. write_radial_fluxes) then
         if (debug) write (*, *) 'diagnostics::diagnostics_fluxes::write_fluxes_for_fluxtube_radialvariation'
         call write_fluxes_for_fluxtube_radialvariation(nout, pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts, write_to_netcdf_file)

      ! Calculate the fluxes if <full_flux_surface> = True
      else if (full_flux_surface) then 
         if (debug) write (*, *) 'diagnostics::diagnostics_fluxes::write_fluxes_for_fullfluxsurface'
         call write_fluxes_for_fullfluxsurface(pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts)

      ! Calculate the fluxes for a flux tube simulation
      else
         if (debug) write (*, *) 'diagnostics::diagnostics_fluxes::write_fluxes_for_fluxtube' 
         call write_fluxes_for_fluxtube(pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts)
      end if

      ! Write fluxes to the ascii files (these variables have been set along the way)
      if (debug) write (*, *) 'diagnostics::diagnostics_fluxes::write_fluxes_to_ascii_file'
      if (proc0) call write_fluxes_to_ascii_file(pflux_vs_s, vflux_vs_s, qflux_vs_s)
      
      !**********************************************************************
      !                         WRITE TO NETCDF FILE                        !
      !**********************************************************************
      ! We already calculate fluxes(kx,ky,z,s) when we were calculating 
      ! <pflux_vs_s>, <vflux_vs_s>, <qflux_vs_s> for the txt files so now
      ! we simply need to write the fluxes to the netcdf file
      !**********************************************************************
      
      ! Do not continue if we do not wish to write to the netCDF file right now
      if (.not. write_to_netcdf_file) return  

      ! Write fluxes(s) to the netcdf file
      if (write_fluxes_vs_time) then 
         if (debug) write (*, *) 'diagnostics::diagnostics_fluxes::write_fluxes_nc'
         if (proc0) call write_fluxes_vs_time_nc(nout, pflux_vs_s, vflux_vs_s, qflux_vs_s)  
      end if

      ! Write fluxes(kx,ky,z,s) to the netcdf file
      if (write_fluxes_kxkyz) then 
         if (debug) write (*, *) 'diagnostics::diagnostics_fluxes::write_fluxes_kxkyz_nc'
         if (proc0) call write_fluxes_kxkyz_nc(nout, pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts)  
      end if

      ! Deallocate the arrays 
      if (write_fluxes_kxkyz) deallocate (pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts)

      ! End timer
      if (proc0) call time_message(.false., timer(:), 'Write fluxes')

   end subroutine write_fluxes_to_netcdf_file

   !============================================================================
   !================================ FLUX TUBE =================================
   !============================================================================
   subroutine write_fluxes_for_fluxtube(pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts)
   
      ! Flags
      use physics_flags, only: include_apar, include_bpar

      ! Load data 
      use dist_fn_arrays, only: gnew, gvmu
      use fields_arrays, only: phi, bpar
      use run_parameters, only: fphi

      ! Redistribute data from  i[vpa,mu,s] to i[kx,ky,z,s] 
      use redistribute, only: scatter
      use dist_redistribute, only: kxkyz2vmu

      ! Calculations 
      use diagnostics_fluxes_fluxtube, only: calculate_fluxes_fluxtube
      use g_tofrom_h, only: g_to_h, g_to_f

      implicit none    

      ! We want to write flux(ky,kx,z,tube,s) to the netcdf file
      real, dimension(:, :, :, :, :), intent(out) :: pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts

      !---------------------------------------------------------------------- 

      ! Redistribute the data from <gnew>(ky,kx,z,tube,i[vpa,mu,s]) to <gvmu>(vpa,mu,i[kx,ky,z,s]),
      ! to ensure that the velocity data is available on each processor.
      call scatter(kxkyz2vmu, gnew, gvmu) 

      ! The <get_fluxes> subroutine assumes that the trubulent component of the distribution function, δf, is passed in.
      ! We know that δf = h - Zs*e/Ts * phi * F0, assuming now that h = <h>_gyroaverage and defining g = <delta f>_gyroaverage 
      ! we have g = h - Zs*e/Ts * <phi>_gyroaverage * F0 and δf = g + Zs*e/Ts * (<phi>_gyroaverage - phi) * F0   
      ! It's been tested numerically, and whether we give <g>, <h> or <δf> does not make a difference for 
      ! <qflux> or <pflux>, but it does matter for <vflux>! Only <δf> is the correct options for <vflux> TODO is it?
      ! TODO-GA for electromagnetic stella the equations are written for f, for electromagnetic stella the equations are written for h
      if (include_apar .or. include_bpar) then  
         call g_to_h(gvmu, phi, bpar, fphi)
      else if (.not. include_apar .and. .not. include_bpar) then 
         call g_to_f(gvmu, phi, fphi)
      end if 

      ! Now calculate the fluxes explicitly
      call calculate_fluxes_fluxtube(gvmu, pflux_vs_s, vflux_vs_s, qflux_vs_s, pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts)

      ! Convert <δf> back to <g> since it will be used by other routines 
      ! TODO-GA for electromagnetic stella the equations are written for f, for electromagnetic stella the equations are written for h
      if (include_apar .or. include_bpar) then  
         call g_to_h(gvmu, phi, bpar, -fphi)
      else if (.not. include_apar .and. .not. include_bpar) then 
         call g_to_f(gvmu, phi, -fphi)
      end if 

   end subroutine write_fluxes_for_fluxtube

   !============================================================================
   !============================= RADIAL VARIATION =============================
   !============================================================================
   subroutine write_fluxes_for_fluxtube_radialvariation(nout, pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts, write_to_netcdf_file)
   
      ! Flags 
      use parameters_diagnostics, only: write_radial_fluxes 

      ! Data 
      use fields_arrays, only: phi, phi_corr_QN
      use physics_flags, only: radial_variation
      use dist_fn_arrays, only: gnew
   
      ! Dimensions
      use kt_grids, only: nakx, naky
      use zgrid, only: nzgrid, ntubes
      use species, only: nspec
   
      ! Write data 
      use stella_io, only: write_radial_fluxes_nc
      use stella_io, only: write_radial_fluxes_nc

      ! Calculations 
      use diagnostics_fluxes_radialvariation, only: calculate_fluxes_radialvariation

      ! Routines
      use job_manage, only: time_message
      use mp, only: proc0

      implicit none 

      ! When writing to the ascii files we don't always write to the netcdf file
      logical, intent(in) :: write_to_netcdf_file  
      integer, intent(in) :: nout    

      ! We want to write flux(ky,kx,z,tube,s) to the netcdf file
      real, dimension(:, :, :, :, :), intent(out) :: pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts

      ! Variables needed to write and calculate diagnostics  
      real, dimension(:, :), allocatable :: pflux_vs_kxs, vflux_vs_kxs, qflux_vs_kxs 
      complex, dimension(:, :, :, :), allocatable :: phi_out

      !---------------------------------------------------------------------- 

      ! Allocate the local arrays
      allocate (phi_out(naky, nakx, -nzgrid:nzgrid, ntubes)); phi_out = 0.0
      allocate (pflux_vs_kxs(nakx, nspec)); pflux_vs_kxs = 0.0
      allocate (vflux_vs_kxs(nakx, nspec)); vflux_vs_kxs = 0.0
      allocate (qflux_vs_kxs(nakx, nspec)); qflux_vs_kxs = 0.0

      ! Get <phi_out>(ky,kx,z,tube)  
      phi_out = phi
      if (radial_variation) then
         phi_out = phi_out + phi_corr_QN
      end if

      ! Calculate the fluxes explicitly
      call calculate_fluxes_radialvariation(gnew, phi_out, pflux_vs_s, vflux_vs_s, qflux_vs_s,  pflux_vs_kxs, vflux_vs_kxs, & 
            qflux_vs_kxs, pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts)

      ! Write fluxes to the netcdf file
      if (write_radial_fluxes .and. write_to_netcdf_file) then 
         if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_radial_fluxes_nc'
         if (proc0) call write_radial_fluxes_nc(nout, pflux_vs_kxs, vflux_vs_kxs, qflux_vs_kxs)
      end if
      
      ! Deallocate the local arrays
      deallocate (pflux_vs_kxs, vflux_vs_kxs, qflux_vs_kxs, phi_out)    

   end subroutine write_fluxes_for_fluxtube_radialvariation

   !============================================================================
   !============================= FULL FLUX SURFACE ============================
   !============================================================================
   subroutine write_fluxes_for_fullfluxsurface(pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts)

      ! Data 
      use dist_fn_arrays, only: gnew

      ! Dimensions 
      use kt_grids, only: ny, ikx_max
      use species, only: nspec
      use zgrid, only: nzgrid

      ! Calculations
      use diagnostics_fluxes_fullfluxsurface, only: calculate_moments_fullfluxsurface
      use diagnostics_fluxes_fullfluxsurface, only: calculate_fluxes_fullfluxsurface

      ! Routines 
      use job_manage, only: time_message 

      implicit none 

      ! We want to write flux(ky,kx,z,tube,s) to the netcdf file
      real, dimension(:, :, :, :, :), intent(out) :: pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts

      ! Variables needed to write and calculate diagnostics   
      complex, dimension(:, :, :, :), allocatable :: dens_vs_ykxzs, upar_vs_ykxzs, pres_vs_ykxzs 

      !---------------------------------------------------------------------- 

      ! Allocate the arrays 
      allocate (dens_vs_ykxzs(ny, ikx_max, -nzgrid:nzgrid, nspec)); dens_vs_ykxzs = 0.0
      allocate (upar_vs_ykxzs(ny, ikx_max, -nzgrid:nzgrid, nspec)); upar_vs_ykxzs = 0.0
      allocate (pres_vs_ykxzs(ny, ikx_max, -nzgrid:nzgrid, nspec)); pres_vs_ykxzs = 0.0

      ! Calculate the particle density, parallel flow and pressure in (y,kx,z) space for all species
      call calculate_moments_fullfluxsurface(gnew, dens_vs_ykxzs, upar_vs_ykxzs, pres_vs_ykxzs) 

      ! Calculate the (ky,kx) contributions to the particle, parallel momentum and energy fluxes
      call calculate_fluxes_fullfluxsurface(dens_vs_ykxzs, upar_vs_ykxzs, pres_vs_ykxzs, & 
            pflux_vs_s, vflux_vs_s, qflux_vs_s, pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts)

      ! Deallocate the arrays for the fluxes
      deallocate (dens_vs_ykxzs, upar_vs_ykxzs, pres_vs_ykxzs) 

   end subroutine write_fluxes_for_fullfluxsurface 

!###############################################################################
!############################### WRITE FILES ###################################
!###############################################################################

   !=========================================================================
   !====================== WRITE FLUXES TO ASCII FILE =======================
   !=========================================================================  
   subroutine write_fluxes_to_ascii_file(pflux_vs_s, vflux_vs_s, qflux_vs_s)

      use stella_time, only: code_time  
      use species, only: nspec
      use mp, only: proc0

      implicit none

      real, dimension(:), intent(in) :: pflux_vs_s, vflux_vs_s, qflux_vs_s
 
      ! Strings to define the format specifier 
      character(3) :: nspec_str
      character(100) :: str 

      !----------------------------------------------------------------------

      ! We only write to the ascii file on the first processor
      if (.not. proc0) return 

      ! For <nspec>=2 the format specifier is <str>='(10es15.4e3)'.
      ! We print 10 columns in scientific form using a total of 15 spaces, hence (10es15), e.g., '4.6942E-011'.
      ! There are 4 spaces after the comma and 3 spaces for the exponential (4e3) and hence 4 spaces between the columns.
      write (nspec_str, '(i3)') 3 * nspec + 1
      str = trim('('//trim(nspec_str)//'es15.4e3)')
      write (fluxes_unit, str) code_time, pflux_vs_s, vflux_vs_s, qflux_vs_s 

      ! Flush the data from the buffer to the actual ascii file
      call flush(fluxes_unit) 

   end subroutine write_fluxes_to_ascii_file 

   !============================================================================
   !========================== OPEN FLUXES ASCII FILE ==========================
   !============================================================================ 
   ! Open the '.fluxes' ascii files. When running a new simulation, create a new file
   ! or replace an old file. When restarting a simulation, append to the old files.
   subroutine open_fluxes_ascii_file(restart)

      use file_utils, only: open_output_file
      use species, only: nspec
      use mp, only: proc0

      implicit none

      logical, intent(in) :: restart

      character(3) :: nspec_str
      character(100) :: str
      logical :: overwrite

      !----------------------------------------------------------------------

      ! We only open the ascii file on the first processor
      if (.not. proc0) return   
 
      ! For a new simulation <overwrite> = True since we wish to create a new ascii file.   
      ! For a restart <overwrite> = False since we wish to append to the existing file. 
      overwrite = .not. restart

      ! Open the '.fluxes' files.
      call open_output_file(fluxes_unit, '.fluxes', overwrite)

      ! Write the header for the '.fluxes' file.
      ! Every column is made up of 15 (12+3) spaces, we have the following columns for nspec=3:
      !   #time   pflux1   pflux2   pflux3   vflux1   vflux2   vflux3   qflux1   qflux2   qflux3
      ! Here <str> is the format string, for 2 species we have <str> = (a12,a14,a30,a30)
      if (.not. restart) then
         write (nspec_str, '(i3)') nspec * 15
         str = trim('(a12,a14,a'//trim(nspec_str)//',a'//trim(nspec_str)//')')
         write (fluxes_unit, str) '#time', 'pflux', 'vflux', 'qflux' 
      end if

   end subroutine open_fluxes_ascii_file

!###############################################################################
!############################ INITALIZE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !======================== INITALIZE THE DIAGNOSTICS =========================
   !============================================================================  
   subroutine init_diagnostics_fluxes(restart)
  
      use species, only: nspec 
      use mp, only: proc0

      implicit none 
 
      logical, intent(in) :: restart 

      !----------------------------------------------------------------------

      ! Only debug on the first processor
      debug = debug .and. proc0

      ! Allocate the arrays for the fluxes
      ! These are needed on all processors since <get_one_flux> will add data to it from each processor
      allocate (qflux_vs_s(nspec)); qflux_vs_s = 0.
      allocate (pflux_vs_s(nspec)); pflux_vs_s = 0.
      allocate (vflux_vs_s(nspec)); vflux_vs_s = 0.

      ! We only open the ascii file on the first processor 
      if (.not. proc0) return         

      ! Open the '.fluxes' ascii file  
      call open_fluxes_ascii_file(restart) 

   end subroutine init_diagnostics_fluxes

   !============================================================================
   !======================== FINALIZE THE DIAGNOSTICS =========================
   !============================================================================  
   subroutine finish_diagnostics_fluxes  

      use file_utils, only: close_output_file
      use mp, only: proc0 

      implicit none

      ! We only have the module arrays on the first processor 
      if (.not. proc0) return    

      ! Deallocate the arrays  
      if (allocated(qflux_vs_s)) deallocate (qflux_vs_s)
      if (allocated(pflux_vs_s)) deallocate (pflux_vs_s)
      if (allocated(vflux_vs_s)) deallocate (vflux_vs_s)   

      ! Close the ascii file
      call close_output_file(fluxes_unit)

   end subroutine finish_diagnostics_fluxes 

end module diagnostics_fluxes

