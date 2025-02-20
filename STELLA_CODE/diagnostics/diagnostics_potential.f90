
!###############################################################################
!############################# DIAGNOSE POTENTIAL ##############################
!###############################################################################
! 
! Routines for calculating and writing the potential.  
! 
!###############################################################################
 
module diagnostics_potential

   implicit none
 
   public :: init_diagnostics_potential
   public :: finish_diagnostics_potential
   public :: write_potential_to_netcdf_file  

   private 

   ! The <units> are used to identify the external ascii files
   integer :: stdout_unit     

   ! Debugging
   logical :: debug = .false. 

contains

!###############################################################################
!############################### WRITE POTENTIAL ###############################
!###############################################################################

   !============================================================================
   !=============== CALCULATE AND WRITE POTENTIAL TO NETCDF FILE ===============
   !============================================================================
   subroutine write_potential_to_netcdf_file(istep, nout, timer, write_to_netcdf_file)

      ! Data 
      use arrays_fields, only: phi, apar, bpar, phi_corr_QN

      ! Dimensions
      use parameters_kxky_grids, only: naky, nakx
      use zgrid, only: ntubes, nzgrid 

      ! Flags 
      use parameters_physics, only: radial_variation
      use parameters_physics, only: include_apar, include_bpar

      ! Calculations 
      use volume_averages, only: volume_average, fieldline_average
      use fields, only: advance_fields

      ! Write to netcdf file 
      use stella_io, only: write_time_nc
      use stella_io, only: write_phi2_nc
      use stella_io, only: write_apar2_nc
      use stella_io, only: write_bpar2_nc
      use stella_io, only: write_phi_nc
      use stella_io, only: write_apar_nc
      use stella_io, only: write_bpar_nc
      use stella_io, only: write_kspectra_nc
      use stella_time, only: code_dt, code_time

      ! Routines
      use job_manage, only: time_message
      use mp, only: proc0
      
      ! Input file
      use parameters_diagnostics, only: write_phi_vs_kxkyz
      use parameters_diagnostics, only: write_apar_vs_kxkyz
      use parameters_diagnostics, only: write_bpar_vs_kxkyz
      use parameters_diagnostics, only: write_phi2_vs_kxky
      use parameters_diagnostics, only: write_apar2_vs_kxky
      use parameters_diagnostics, only: write_bpar2_vs_kxky

      implicit none 

      integer, intent(in) :: istep  ! The current time step  
      integer, intent(in) :: nout   ! The pointer in the netcdf file
      logical, intent(in) :: write_to_netcdf_file    
      real, dimension(:), intent(in out) :: timer   

      ! Variables needed to write and calculate diagnostics
      complex, dimension(:, :, :, :), allocatable :: phi_vs_kykxzt, apar_vs_kykxzt, bpar_vs_kykxzt
      real, dimension(:, :), allocatable :: phi2_vs_kxky, apar2_vs_kxky, bpar2_vs_kxky
      real :: phi2, apar2, bpar2 

      !---------------------------------------------------------------------- 

      ! Start timer
      if (proc0) call time_message(.false., timer(:), 'Write phi')

      ! Allocate arrays 
      allocate (phi_vs_kykxzt(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (apar_vs_kykxzt(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (bpar_vs_kykxzt(naky, nakx, -nzgrid:nzgrid, ntubes))

      ! Get <phi_vs_kykxzt>(ky,kx,z,tube)  
      phi_vs_kykxzt = phi
      apar_vs_kykxzt = apar
      bpar_vs_kykxzt = bpar
      if (radial_variation) then
         phi_vs_kykxzt = phi_vs_kykxzt + phi_corr_QN
      end if
     
      ! Write the potential squared to the ascii file and to the standard output file with unit=*.
      ! The header of the standard output file is printed in the <stella.f90> routine.
      if (proc0) then 
         call volume_average(phi_vs_kykxzt, phi2)
         call volume_average(apar_vs_kykxzt, apar2) 
         call volume_average(bpar_vs_kykxzt, bpar2)
         call write_potential_to_ascii_file(istep, phi2, apar2, bpar2) 
         if (include_apar .and. include_bpar) then 
            write (*, '(A2,I7,A2,ES12.4,A2,ES12.4,A2,ES12.4,A2,ES12.4,A2,ES12.4)') &
            " ", istep, " ", code_time, " ", code_dt, " ", phi2, " ", apar2, " ", bpar2
         else if (include_apar) then 
            write (*, '(A2,I7,A2,ES12.4,A2,ES12.4,A2,ES12.4,A2,ES12.4)') &
            " ", istep, " ", code_time, " ", code_dt, " ", phi2, " ", apar2
         else
            write (*, '(A2,I7,A2,ES12.4,A2,ES12.4,A2,ES12.4)') &
            " ", istep, " ", code_time, " ", code_dt, " ", phi2
         end if
      end if

      ! Write the potential to the netcdf file 
      if (proc0 .and. write_to_netcdf_file) then

         ! Write the time axis to the netcdf file
         if (debug) write (*, *) 'diagnostics::write_time_nc'
         call write_time_nc(nout, code_time)

         ! Write the phi2(t) to the netcdf file (always on)
         if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_phi2_nc'
         call write_phi2_nc(nout, phi2)
         call write_apar2_nc(nout, apar2)
         call write_bpar2_nc(nout, bpar2)

         ! Write phi(t,ky,kx,z,tube), apar(t,ky,kx,z,tube) and bpar(t,ky,kx,z,tube) to the netcdf file
         if (write_phi_vs_kxkyz) then
            if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_phi_nc'
            call write_phi_nc(nout, phi_vs_kykxzt)
         end if
         if (write_apar_vs_kxkyz .and. include_apar) then
            if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_apar_nc'
            call write_apar_nc(nout, apar_vs_kykxzt)
         end if
         if (write_bpar_vs_kxkyz .and. include_bpar) then
            if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_bpar_nc'
            call write_bpar_nc(nout, bpar_vs_kykxzt)
         end if

         ! Write phi2(t,ky,kx) to the netcdf file
         if (write_phi2_vs_kxky) then
            if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_kspectra_nc' 
            allocate (phi2_vs_kxky(naky, nakx)) 
            call fieldline_average(real(phi_vs_kykxzt * conjg(phi_vs_kykxzt)), phi2_vs_kxky)
            call write_kspectra_nc(nout, phi2_vs_kxky, "phi2_vs_kxky", "electrostatic potential")
            deallocate (phi2_vs_kxky)
         end if
         if (write_apar2_vs_kxky .and. include_apar) then 
            allocate (apar2_vs_kxky(naky, nakx))
            call fieldline_average(real(apar_vs_kykxzt * conjg(apar_vs_kykxzt)), apar2_vs_kxky)
            call write_kspectra_nc(nout, apar2_vs_kxky, "apar2_vs_kxky", "parallel vector potential")
            deallocate (apar2_vs_kxky)
         end if 
         if (write_bpar2_vs_kxky .and. include_bpar) then
            allocate (bpar2_vs_kxky(naky, nakx))
            call fieldline_average(real(bpar_vs_kykxzt * conjg(bpar_vs_kykxzt)), bpar2_vs_kxky)
            call write_kspectra_nc(nout, bpar2_vs_kxky, "bpar2_vs_kxky", "parallel magnetic field fluctuation")
            deallocate (bpar2_vs_kxky)
         end if 
         
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

      ! Deallocate arrays 
      deallocate (phi_vs_kykxzt)
      deallocate (apar_vs_kykxzt)
      deallocate (bpar_vs_kykxzt)

      ! End timer
      if (proc0) call time_message(.false., timer(:), 'Write phi')

   end subroutine write_potential_to_netcdf_file

!###############################################################################
!############################### WRITE FILES ###################################
!###############################################################################

   !=========================================================================
   !===================== WRITE POTENTIAL TO ASCII FILE =====================
   !=========================================================================  
   subroutine write_potential_to_ascii_file(istep, phi2, apar2, bpar2)

      use stella_time, only: code_time 

      implicit none

      integer, intent(in) :: istep
      real, intent(in) :: phi2, apar2, bpar2

      write (stdout_unit, '(i10,es14.4e3,es14.4e3,es14.4e3,es14.4e3)') istep, code_time, phi2, apar2, bpar2
      call flush (stdout_unit) 

   end subroutine write_potential_to_ascii_file 

   !=========================================================================
   !======== WRITE POTENTIAL(Z) AT THE FINAL TIME STEP TO ASCII FILE ========
   !=========================================================================  
   subroutine write_potential_to_ascii_file_atfinaltimestep

      ! Data 
      use arrays_fields, only: phi, apar, bpar

      ! Geometry 
      USE arrays_dist_fn, only: kperp2
      use geometry, only: zed_eqarc

      ! Dimensions
      use parameters_kxky_grids, only: naky, nakx
      use grids_kxky, only: aky, akx, zed0
      use zgrid, only: nzgrid, ntubes, zed

      ! Routines 
      use file_utils, only: open_output_file, close_output_file

      implicit none

      integer :: tmpunit
      integer :: iky, ikx, iz, it

      !----------------------------------------------------------------------

      ! Open a new ascii file
      call open_output_file(tmpunit, '.final_fields')

      ! Write the header
      write (tmpunit, '(10a14)') '# z', 'z-zed0', 'aky', 'akx', &
         'real(phi)', 'imag(phi)', 'real(apar)', 'imag(apar)', &
         'real(bpar)', 'imag(bpar)', 'z_eqarc-zed0', 'kperp2'

      ! Write the final fields versus z
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

      ! Close the ascii file
      call close_output_file(tmpunit)

   end subroutine write_potential_to_ascii_file_atfinaltimestep

   !============================================================================
   !========================== OPEN FLUXES ASCII FILE ==========================
   !============================================================================ 
   ! Open the '.fluxes' ascii files. When running a new simulation, create a new file
   ! or replace an old file. When restarting a simulation, append to the old files.
   subroutine open_potential_ascii_file(restart)

      use file_utils, only: open_output_file 
      use mp, only: proc0

      implicit none

      logical, intent(in) :: restart
      logical :: overwrite

      !----------------------------------------------------------------------

      ! We only open the ascii file on the first processor
      if (.not. proc0) return   
 
      ! For a new simulation <overwrite> = True since we wish to create a new ascii file.   
      ! For a restart <overwrite> = False since we wish to append to the existing file. 
      overwrite = .not. restart

      ! Open the '.out' files.
      call open_output_file(stdout_unit, '.out', overwrite)

      ! Write the header for the '.out' file.
      if (.not. restart) then
         write (stdout_unit, '(a10,a11,a15,a15,a15)') 'istep', 'time', '|phi|^2', '|apar|^2', '|bpar|^2' 
         write (stdout_unit, '(a)') '---------------------------------------------------------------------' 
      end if

   end subroutine open_potential_ascii_file

!###############################################################################
!############################ INITALIZE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !======================== INITALIZE THE DIAGNOSTICS =========================
   !============================================================================  
   subroutine init_diagnostics_potential(restart)  
  
      use mp, only: proc0

      implicit none 
 
      logical, intent(in) :: restart 

      !---------------------------------------------------------------------

      ! We only need the module arrays on the first processor 
      if (.not. proc0) return          

      ! Open the '.out' ascii file  
      call open_potential_ascii_file(restart) 

   end subroutine init_diagnostics_potential

   !============================================================================
   !======================== FINALIZE THE DIAGNOSTICS =========================
   !============================================================================  
   subroutine finish_diagnostics_potential

      use file_utils, only: close_output_file
      use mp, only: proc0 

      implicit none

      ! We only have the module arrays on the first processor 
      if (.not. proc0) return     

      ! Close the ascii file
      call close_output_file(stdout_unit)

      ! Write potential(z) at the last time step to an ascii file
      call write_potential_to_ascii_file_atfinaltimestep

   end subroutine finish_diagnostics_potential

end module diagnostics_potential
