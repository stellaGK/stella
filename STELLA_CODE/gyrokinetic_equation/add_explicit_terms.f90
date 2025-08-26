!###############################################################################
!###############################################################################
!###############################################################################
! 
! This module is used in the following subroutines, to add factor*field
! to the right-hand-side of the gyrokinetic equation, looping over ivmu.
!     - advance_wdriftx_explicit
!     - advance_wdrifty_explicit
!     - advance_wstar_explicit
! 
!###############################################################################
module add_explicit_terms

   implicit none
   
   ! Make these routine available
   public :: add_explicit_term
   public :: add_explicit_term_ffs

   private

contains

   subroutine add_explicit_term(g, pre_factor, src)

      use stella_layouts, only: vmu_lo
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(-nzgrid:, vmu_lo%llim_proc:), intent(in) :: pre_factor
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: src

      integer :: ivmu
      integer :: iky, ikx, iz, it

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
               do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  do iky = 1, naky
                     src(iky, ikx, iz, it, ivmu) = src(iky, ikx, iz, it, ivmu) + pre_factor(iz, ivmu) * g(iky, ikx, iz, it, ivmu)
                  end do
               end do
               end do
         end do
      end do

   end subroutine add_explicit_term

   ! add vM . grad y d<phi>/dy or vM . grad x d<phi>/dx (or equivalents with g) or omega_* * d<phi>/dy term to RHS of GK equation
   subroutine add_explicit_term_ffs(g, pre_factor, src)

      use stella_layouts, only: vmu_lo
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: ikx_max, nalpha

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(in) :: pre_factor
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: src

      integer :: ivmu
      integer :: ia, ikx, iz, it

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
               do iz = -nzgrid, nzgrid
               do ikx = 1, ikx_max
                  do ia = 1, nalpha
                     src(ia, ikx, iz, it, ivmu) = src(ia, ikx, iz, it, ivmu) + pre_factor(ia, iz, ivmu) * g(ia, ikx, iz, it, ivmu)
                  end do
               end do
               end do
         end do
      end do

   end subroutine add_explicit_term_ffs

end module add_explicit_terms
