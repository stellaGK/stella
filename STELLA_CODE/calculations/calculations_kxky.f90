!###############################################################################
!                                                                               
!###############################################################################
! This module ...
!###############################################################################
module calculations_kxky

   implicit none
   
   ! Make routines available to other modules
   public :: swap_kxky, swap_kxky_back
   public :: swap_kxky_ordered, swap_kxky_back_ordered
   public :: multiply_by_rho
   public :: communicate_ktgrids_multibox
   
   private
   
   interface swap_kxky
      module procedure swap_kxky_real
      module procedure swap_kxky_complex
   end interface swap_kxky
   
   interface swap_kxky_ordered
      module procedure swap_kxky_ordered_real
      module procedure swap_kxky_ordered_complex
   end interface swap_kxky_ordered
   
contains

!###############################################################################
!################################# CALCULATIONS ################################
!###############################################################################

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   ! Take an array with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
   ! and uses reality condition to return array
   ! with kx >= 0 and all ky (ordered like 0, ..., kymax, -kymax, ..., -dky)
   !****************************************************************************
   subroutine swap_kxky_complex(gin, gout)
      
      use grids_kxky, only: naky, naky_all, ikx_max, nakx
      
      implicit none
      
      ! Arguments
      complex, dimension(:, :), intent(in) :: gin
      complex, dimension(:, :), intent(out) :: gout
      
      ! Local variables
      integer :: ikx, ikxneg
      integer :: iky, ikyneg

      !----------------------------------------------------------------------
      
      ! First set arrays equal for ky >= 0 and kx >= 0
      gout(:naky, :) = gin(:, :ikx_max)
      
      ! Next fill in ky < 0, kx >= 0 elements of array using reality
      ikx = 1
      ikxneg = ikx
      
      do iky = naky + 1, naky_all
         ! This is the ky index corresponding to +ky in original array
         ikyneg = naky_all - iky + 2
         gout(iky, ikx) = conjg(gin(ikyneg, ikxneg))
      end do
      
      do ikx = 2, ikx_max
         ikxneg = nakx - ikx + 2
         do iky = naky + 1, naky_all
            ! This is the ky index corresponding to +ky in original array
            ikyneg = naky_all - iky + 2
            gout(iky, ikx) = conjg(gin(ikyneg, ikxneg))
         end do
      end do
      
   end subroutine swap_kxky_complex

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   ! Take an array with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
   ! and uses reality condition to return array
   ! with kx >= 0 and all ky (ordered like 0, ..., kymax, -kymax, ..., -dky)
   !****************************************************************************
   subroutine swap_kxky_real(gin, gout)
      
      use grids_kxky, only: naky, naky_all, ikx_max, nakx
      
      implicit none
      
      ! Arguments
      real, dimension(:, :), intent(in) :: gin
      real, dimension(:, :), intent(out) :: gout
      
      ! Local variables
      integer :: ikx, ikxneg
      integer :: iky, ikyneg

      !----------------------------------------------------------------------
      
      ! First set arrays equal for ky >= 0 and kx >= 0
      gout(:naky, :) = gin(:, :ikx_max)
      
      ! Next fill in ky < 0, kx >= 0 elements of array using reality
      ikx = 1
      ikxneg = ikx
      
      do iky = naky + 1, naky_all
         ! This is the ky index corresponding to +ky in original array
         ikyneg = naky_all - iky + 2
         gout(iky, ikx) = gin(ikyneg, ikxneg)
      end do
      
      do ikx = 2, ikx_max
         ikxneg = nakx - ikx + 2
         do iky = naky + 1, naky_all
            ! This is the ky index corresponding to +ky in original array
            ikyneg = naky_all - iky + 2
            gout(iky, ikx) = gin(ikyneg, ikxneg)
         end do
      end do
      
   end subroutine swap_kxky_real
   
   !****************************************************************************
   !                                      Title
   !****************************************************************************
   ! Take an array with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
   ! and uses reality condition to return array
   ! with kx >= 0 and all ky (ordered like -kymax, ..., 0, ..., kymax)
   !****************************************************************************
   subroutine swap_kxky_ordered_real(gin, gout)
      
      use grids_kxky, only: ikx_max, naky, nakx
      
      implicit none
      
      ! Arguments
      real, dimension(:, :), intent(in) :: gin
      real, dimension(:, :), intent(out) :: gout
      
      ! Local variables
      integer :: ikx, ikxneg
      integer :: iky, ikyneg

      !----------------------------------------------------------------------
      
      ! First set arrays equal for ky >= 0 and kx >= 0
      gout(naky:, :) = gin(:, :ikx_max)
      
      ! Next fill in ky < 0, kx >= 0 elements of array using reality
      ikx = 1
      ikxneg = ikx
      
      do iky = 1, naky - 1
         ! this is the ky index corresponding to +ky in original array
         ikyneg = naky - iky + 1
         gout(iky, ikx) = gin(ikyneg, ikxneg)
      end do
      
      do ikx = 2, ikx_max
         ikxneg = nakx - ikx + 2
         do iky = 1, naky - 1
            ! This is the ky index corresponding to +ky in original array
            ikyneg = naky - iky + 1
            gout(iky, ikx) = gin(ikyneg, ikxneg)
         end do
      end do
      
   end subroutine swap_kxky_ordered_real
   
   !****************************************************************************
   !                                      Title
   !****************************************************************************
   ! take an array with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
   ! and uses reality condition to return array
   ! with kx >= 0 and all ky (ordered like -kymax, ..., 0, ..., kymax)
   !****************************************************************************
   subroutine swap_kxky_ordered_complex(gin, gout)
      
      use grids_kxky, only: naky, nakx, ikx_max
      
      implicit none
         
      ! Arguments
      complex, dimension(:, :), intent(in) :: gin
      complex, dimension(:, :), intent(out) :: gout
      
      ! Local variables
      integer :: ikx, ikxneg
      integer :: iky, ikyneg

      !----------------------------------------------------------------------
      
      ! First set arrays equal for ky >= 0 and kx >= 0
      gout(naky:, :) = gin(:, :ikx_max)
      
      ! Next fill in ky < 0, kx >= 0 elements of array using reality
      ikx = 1
      ikxneg = ikx
      
      do iky = 1, naky - 1
         ! This is the ky index corresponding to +ky in original array
         ikyneg = naky - iky + 1
         gout(iky, ikx) = conjg(gin(ikyneg, ikxneg))
      end do
      
      do ikx = 2, ikx_max
         ikxneg = nakx - ikx + 2
         do iky = 1, naky - 1
            ! This is the ky index corresponding to +ky in original array
            ikyneg = naky - iky + 1
            gout(iky, ikx) = conjg(gin(ikyneg, ikxneg))
         end do
      end do
      
   end subroutine swap_kxky_ordered_complex
   
   !****************************************************************************
   !                                      Title
   !****************************************************************************
   ! Take an array with kx >= 0 and all ky (ordered like 0, ..., kymax, -kymax, ..., -dky)
   ! and uses reality condition to return array
   ! with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
   !****************************************************************************
   subroutine swap_kxky_back(gin, gout)
      
      use grids_kxky, only: naky, nakx, naky_all, ikx_max
      
      implicit none
      
      ! Arguments
      complex, dimension(:, :), intent(in) :: gin
      complex, dimension(:, :), intent(out) :: gout
      
      ! Local variables
      integer :: ikx, ikxneg
      integer :: iky, ikyneg

      !----------------------------------------------------------------------
      
      ! First set arrays equal for ky >= 0 and kx >= 0
      gout(:, :ikx_max) = gin(:naky, :)
      
      ! Next fill in kx < 0, ky >= 0 elements of array using reality
      do ikx = ikx_max + 1, nakx
         ikxneg = nakx - ikx + 2
         iky = 1
         ikyneg = iky
         gout(iky, ikx) = conjg(gin(ikyneg, ikxneg))
         do iky = 2, naky
            ikyneg = naky_all - iky + 2
            gout(iky, ikx) = conjg(gin(ikyneg, ikxneg))
         end do
      end do
      
   end subroutine swap_kxky_back
   
   !****************************************************************************
   !                                      Title
   !****************************************************************************
   ! Take an array with kx >= 0 and all ky (ordered like -kymax, ..., 0, ..., kymax)
   ! and uses reality condition to return array
   ! with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
   !****************************************************************************
   subroutine swap_kxky_back_ordered(gin, gout)
      
      use grids_kxky, only: ikx_max, naky, nakx
      
      implicit none
      
      ! Arguments
      complex, dimension(:, :), intent(in) :: gin
      complex, dimension(:, :), intent(out) :: gout
      
      ! Local variables
      integer :: ikx, ikxneg
      integer :: iky, ikyneg

      !----------------------------------------------------------------------
      
      ! First set arrays equal for ky >= 0 and kx >= 0
      gout(:, :ikx_max) = gin(naky:, :)
      
      ! Next fill in kx < 0, ky >= 0 elements of array using reality
      do ikx = ikx_max + 1, nakx
         ikxneg = nakx - ikx + 2
         do iky = 1, naky
            ikyneg = naky - iky + 1
            gout(iky, ikx) = conjg(gin(ikyneg, ikxneg))
         end do
      end do
      
   end subroutine swap_kxky_back_ordered
   
   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine communicate_ktgrids_multibox
   
      use job_manage, only: njobs
      use mp, only: job, scope, crossdomprocs, subprocs, send, receive
      use grids_kxky, only: phase_shift_angle
      
      implicit none

      !----------------------------------------------------------------------
      
      call scope(crossdomprocs)
      
      if (job == 1) then
         call send(phase_shift_angle, 0, 120)
         call send(phase_shift_angle, njobs - 1, 130)
      elseif (job == 0) then
         call receive(phase_shift_angle, 1, 120)
      elseif (job == njobs - 1) then
         call receive(phase_shift_angle, 1, 130)
      end if
      
      call scope(subprocs)
      
   end subroutine communicate_ktgrids_multibox
   
   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine multiply_by_rho(gin)
      
      use calculations_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use grids_kxky, only: rho_d_clamped, zonal_mode, g0x
      use grids_kxky, only: nakx, naky
      
      implicit none
      
      ! Arguments
      complex, dimension(:, :), intent(inout) :: gin

      !----------------------------------------------------------------------
      
      if (.not. allocated(g0x)) allocate (g0x(naky, nakx))
      
      call transform_kx2x_unpadded(gin, g0x)
      g0x = spread(rho_d_clamped, 1, naky) * g0x
      if (zonal_mode(1)) g0x(1, :) = real(g0x(1, :))
      call transform_x2kx_unpadded(g0x, gin)
      
   end subroutine multiply_by_rho

   end module calculations_kxky
