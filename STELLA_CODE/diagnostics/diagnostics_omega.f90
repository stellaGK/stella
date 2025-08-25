
!###############################################################################
!############################## DIAGNOSE OMEGA #################################
!###############################################################################
! 
! Routines for calculating and writing the growth rate and frequency data.
!  
! The potential can be written as <phi> = exp(-i*<0mega>*t) = exp(-i*(<omega>+i<gamma>)*t)
!     Thus <Omega> = log(dphi)/(-i*dt) = i*log(dphi)/dt 
!     with <omega> = real(<Omega>) and <gamma> = imag(<Omega>)  
! We take the running average of <Omega> over the last <navg> time points. 
! 
!###############################################################################

module diagnostics_omega

   use debug_flags, only: debug => diagnostics_omega_debug
  
   implicit none 

   ! These routines are called from diagnostics.f90
   public :: calculate_omega, write_omega_to_netcdf_file
   public :: init_diagnostics_omega, finish_diagnostics_omega 
   
   ! This routine is called from stella.f90 to automatically stop
   ! stella if <omega> has saturated during linear simulations
   public :: checksaturation

   private    

   ! The <units> are used to identify the external ascii files
   integer :: omega_unit

   ! Keep track of <omega>(t, ky, kx) to calculate the running average   
   ! These arrays are only allocated on the first processor
   complex, dimension(:, :, :), allocatable :: omega_vs_tkykx
   complex, dimension(:, :), allocatable :: omega_vs_kykx
   
   ! Write omega diagnostics
   logical :: write_omega

contains
 
   !=========================================================================
   !================== CALCULATE OMEGA AT EVERY TIME STEP ===================
   !========================================================================= 
   subroutine calculate_omega(istep, timer)
      
      ! Multiprocessing
      use job_manage, only: time_message
      use mp, only: proc0
      
      ! Physics flags
      use parameters_physics, only: include_apar
      
      ! Fields
      use store_arrays_fields, only: phi, phi_old
      use store_arrays_fields, only: apar, apar_old
      
      ! Grids
      use parameters_kxky_grid, only: nakx, naky 
      use stella_time, only: code_dt
      use constants, only: zi
      
      ! Calculations 
      use calculations_volume_averages, only: fieldline_average
      
      ! Input file
      use parameters_diagnostics, only: navg
      use parameters_diagnostics, only: write_omega_vs_kxky
      use parameters_diagnostics, only: write_omega_avg_vs_kxky

      implicit none 

      integer, intent(in) :: istep  ! The current time step   
      real, dimension(:), intent(in out) :: timer  

      ! Temporary arrays 
      complex, dimension(:, :), allocatable :: phi_vs_kykx, phiold_vs_kykx, aparavg, aparoldavg
      integer :: it_runningaverage
      real :: zero

      !----------------------------------------------------------------------

      debug = debug .and. proc0
      
      ! We only calculate omega on the first processor and if <write_omega> = True
      if (.not. write_omega) return  

      ! Start the timer
      if (proc0) call time_message(.false., timer(:), 'calculate omega')
 
      ! Allocate temporary arrays and define a <zero>
      allocate (phi_vs_kykx(naky, nakx)); phi_vs_kykx = 0.0
      allocate (phiold_vs_kykx(naky, nakx)); phiold_vs_kykx = 0.0
      allocate (aparavg(naky, nakx)); aparavg = 0.0
      allocate (aparoldavg(naky, nakx)); aparoldavg = 0.0
      zero = 100.*epsilon(0.) 

      ! Field line average <phi>(ky,kx,z,tube) to obtain <phi>(ky,kx)
      call fieldline_average(phi, phi_vs_kykx)
      call fieldline_average(phi_old, phiold_vs_kykx)
      
      ! If we include electromagnetic terms, field line average <apar> and <par_old>
      if (include_apar) then 
          call fieldline_average(apar, aparavg)
          call fieldline_average(apar_old, aparoldavg)
          ! add <apar> to <phi> in the case <phi> = 0 because the mode has tearing parity
          ! if the mode is a purely growing mode then 
          ! (<phi^n+1> + <apar^n+1> )/(<phi^n> + <apar^n>) = exp (-i delta t omega)  
          phi_vs_kykx = phi_vs_kykx + aparavg
          phiold_vs_kykx = phiold_vs_kykx + aparoldavg
      end if

      ! We save the running average in one of the <navg> time points of <omega_vs_tkykx>
      it_runningaverage = mod(istep, navg) + 1

      ! The potential can be written as <phi> = exp(-i*<0mega>*t) = exp(-i*(<omega>+i<gamma>)*t)
      ! Thus <Omega> = log(dphi)/(-i*dt) = i*log(dphi)/dt with <omega> = real(<Omega>) and <gamma> = imag(<Omega>)  
      where (abs(phi_vs_kykx) < zero .or. abs(phiold_vs_kykx) < zero)
         omega_vs_tkykx(it_runningaverage, :, :) = 0.0
      elsewhere
         omega_vs_tkykx(it_runningaverage, :, :) = log(phi_vs_kykx / phiold_vs_kykx) * zi / code_dt
      end where

      ! Deallocate temporary arrays
      deallocate (phi_vs_kykx, phiold_vs_kykx) 
      deallocate (aparavg, aparoldavg)

      ! End the timer
      if (proc0) call time_message(.false., timer(:), 'calculate omega')

   end subroutine calculate_omega  

   !=========================================================================
   !=============== WRITE OMEGA AT EVERY <NWRITE> TIME STEPS ================
   !========================================================================= 
   subroutine write_omega_to_netcdf_file(istep, nout, timer, write_to_netcdf_file)
 
      ! Input file
      use parameters_diagnostics, only: navg
      use parameters_diagnostics, only: write_omega_vs_kxky
      use parameters_diagnostics, only: write_omega_avg_vs_kxky
      
      ! Write to netCDF file
      use stella_io, only: write_omega_nc
      
      ! Parallelisation
      use job_manage, only: time_message
      use mp, only: proc0

      implicit none 

      integer, intent(in) :: istep  ! The current time step  
      integer, intent(in) :: nout   ! The pointer in the netcdf file
      logical, intent(in) :: write_to_netcdf_file    
      real, dimension(:), intent(in out) :: timer  

      integer :: it_runningaverage

      !----------------------------------------------------------------------

      ! We only calculate omega on the first processor
      if (.not. write_omega) return     
      
      ! Start timer
      if (proc0) call time_message(.false., timer(:), 'write omega')

      ! Get the index of the current time point in <omega_vs_tkykx>
      it_runningaverage = mod(istep, navg) + 1

      ! Calculate the running average of <omega_vs_tkykx>  
      omega_vs_kykx = sum(omega_vs_tkykx, dim=1) / real(navg) 

      ! Write omega and the running average of omega to the ascii file
      call write_omega_to_ascii_file(istep, omega_vs_tkykx(it_runningaverage, :, :), omega_vs_kykx) 

      ! Write omega to the netcdf file  
      if (write_to_netcdf_file) call write_omega_nc(nout, omega_vs_tkykx(it_runningaverage, :, :)) 

      ! Deallocate temporary arrays
      deallocate (omega_vs_kykx) 

     ! End timer 
      if (proc0) call time_message(.false., timer(:), 'write omega') 

   end subroutine write_omega_to_netcdf_file   

   !=========================================================================
   !=========================== CHECK SATURATION ============================
   !========================================================================= 

   subroutine checksaturation(istep, stop_stella)
   
      use parameters_numerical, only: autostop
      use parameters_diagnostics, only: navg
      use mp, only: proc0, broadcast

      integer, intent(in) :: istep  
      logical, intent(in out) :: stop_stella

      logical :: equal
      real :: max_difference
      integer :: i

      !----------------------------------------------------------------------

      ! Only check if gamma is saturated if we want stella to stop automatically
      if (.not. autostop) return 

      ! Check whether (omega, gamma) has saturated
      if (proc0) then 
         if (istep > navg+1) then

            ! Check whether all elements in <omega_vs_tkykx> are the same
            equal = .true. 
            do i = 1, navg
               max_difference = maxval(abs(omega_vs_tkykx(i,:,:) - omega_vs_tkykx(1,:,:))) 
               if (max_difference > 0.000001) then 
                  equal = .false.; exit 
               end if
            end do  

            ! If gamma has saturated, stop stella
            if (equal) then 
               write (*, *)
               write (*, '(A, I0, A)') 'EXITING STELLA BECAUSE (OMEGA, GAMMA) HAS SATURATED OVER ', navg, ' TIMESTEPS' 
               stop_stella = .true.
            end if
      
         end if 
      end if 

      call broadcast(stop_stella)

   end subroutine checksaturation

!###############################################################################
!############################# INITALIZE & FINALIZE ############################
!###############################################################################

   !============================================================================
   !======================== INITALIZE THE DIAGNOSTICS =========================
   !============================================================================  
   subroutine init_diagnostics_omega(restart) 
   
      ! Input file
      use parameters_diagnostics, only: write_omega_avg_vs_kxky
      use parameters_diagnostics, only: write_omega_vs_kxky
      use parameters_diagnostics, only: navg
      
      ! Grids
      use parameters_kxky_grid, only: nakx, naky
      
      ! Multiprocessing
      use mp, only: proc0

      implicit none 
 
      logical, intent(in) :: restart 

      !----------------------------------------------------------------------
      
      ! Determine whether we want to write any omega diagnostics. Only write
      ! on the first processor. Here <write_omega> is a global variable
      write_omega = write_omega_vs_kxky .or. write_omega_avg_vs_kxky
      if (.not. proc0) write_omega = .false.

      ! If we don't write omega, we don't need to open the ASCII or allocate arrays
      if (.not. write_omega) return     

      ! Allocate omega versus (<navg>, ky, kx) to calculate the running average 
      allocate (omega_vs_tkykx(navg, naky, nakx)); omega_vs_tkykx = 0.0
      allocate (omega_vs_kykx(naky, nakx)); omega_vs_kykx = 0.0

      ! Open the '.omega' ascii file  
      call open_omega_ascii_file(restart)

   end subroutine init_diagnostics_omega  

   !============================================================================
   !======================== FINALIZE THE DIAGNOSTICS =========================
   !============================================================================  
   subroutine finish_diagnostics_omega  

      use parameters_diagnostics, only: write_omega_vs_kxky
      use file_utils, only: close_output_file

      implicit none  

      !----------------------------------------------------------------------

      if (.not. write_omega) return     
      if (allocated(omega_vs_tkykx)) deallocate (omega_vs_tkykx)   
      if (allocated(omega_vs_kykx)) deallocate (omega_vs_kykx)  
      call close_output_file(omega_unit)

   end subroutine finish_diagnostics_omega

!###############################################################################
!############################### ASCII FILES ###################################
!###############################################################################

   !============================================================================
   !============================= OPEN ASCII FILES =============================
   !============================================================================ 
   ! Open the '.omega' ascii files. When running a new simulation, create a new file
   ! or replace an old file. When restarting a simulation, append to the old files.
   subroutine open_omega_ascii_file(restart)
 
      ! Input file
      use parameters_diagnostics, only: write_omega_avg_vs_kxky
      use parameters_diagnostics, only: write_omega_vs_kxky
      
      ! ASCII files
      use file_utils, only: open_output_file  

      implicit none

      logical, intent(in) :: restart 
      logical :: overwrite 

      !---------------------------------------------------------------------- 

      ! For a new simulation <overwrite> = True since we wish to create a new ascii file.   
      ! For a restart <overwrite> = False since we wish to append to the existing file. 
      overwrite = .not. restart 

      ! Open the '.omega' file
      call open_output_file(omega_unit, '.omega', overwrite)

      ! Write the header to the '.omega' file
      if (.not. restart) then
         if (write_omega_avg_vs_kxky .and. write_omega_vs_kxky) then 
            write (omega_unit, '(a12,a14,a15,a20,a15,a18,a16)') '#time', 'ky', 'kx', 'Re[om]', 'Im[om]', 'Re[omavg]', 'Im[omavg]'
         else if (write_omega_avg_vs_kxky) then 
            write (omega_unit, '(a12,a14,a15,a20,a15)') '#time', 'ky', 'kx', 'Re[omavg]', 'Im[omavg]'
         else if (write_omega_vs_kxky) then 
            write (omega_unit, '(a12,a14,a15,a20,a17)') '#time', 'ky', 'kx', 'frequency', 'growth rate'
         end if 
      end if 

   end subroutine open_omega_ascii_file

   !=========================================================================
   !======================== WRITE LOOP ASCII FILES =========================
   !=========================================================================  
   subroutine write_omega_to_ascii_file(istep, omega_vs_kykx, omega_runningavg_vs_kykx)
 
      use parameters_diagnostics, only: write_omega_avg_vs_kxky
      use parameters_diagnostics, only: write_omega_vs_kxky
      use stella_time, only: code_time 
      use parameters_kxky_grid, only: naky, nakx
      use grids_kxky, only: aky, akx 

      implicit none

      integer, intent(in) :: istep ! The current time step   
      complex, dimension(:, :), intent(in) :: omega_vs_kykx, omega_runningavg_vs_kykx
      integer :: ikx, iky

      !--------------------------------------------------------------------
      
      ! Dont write omega for istep = 0
      if (istep<=0) return
 
      ! For each mode (kx,ky) write <omega>(ky,kx) to the ascii file
      do iky = 1, naky
         do ikx = 1, nakx
            if (write_omega_avg_vs_kxky .and. write_omega_vs_kxky) then 
               write (omega_unit, '(7e16.8)') code_time, aky(iky), akx(ikx), &
                  real(omega_vs_kykx(iky, ikx)), aimag(omega_vs_kykx(iky, ikx)), &
                  real(omega_runningavg_vs_kykx(iky, ikx)), aimag(omega_runningavg_vs_kykx(iky, ikx))
            else if (write_omega_avg_vs_kxky) then 
               write (omega_unit, '(5e16.8)') code_time, aky(iky), akx(ikx), &
                  real(omega_runningavg_vs_kykx(iky, ikx)), aimag(omega_runningavg_vs_kykx(iky, ikx))
            else if (write_omega_vs_kxky) then 
               write (omega_unit, '(5e16.8)') code_time, aky(iky), akx(ikx), &
                  real(omega_vs_kykx(iky, ikx)), aimag(omega_vs_kykx(iky, ikx))
            end if 
         end do
         if (nakx > 1) write (omega_unit, *)
      end do
      if (naky > 1) write (omega_unit, *)

      ! Flush the data from the buffer to the actual ascii file
      call flush(omega_unit) 

   end subroutine write_omega_to_ascii_file

end module diagnostics_omega


