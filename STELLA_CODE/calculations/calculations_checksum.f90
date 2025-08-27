!###############################################################################
!                                                                               
!###############################################################################
! This module ...
!###############################################################################
module calculations_checksum

   use debug_flags, only: debug => calculations_debug

   implicit none

   public :: checksum

   private

   interface checksum
      module procedure checksum_field
      module procedure checksum_dist
   end interface

contains 

!###############################################################################
!################################# CALCULATIONS ################################
!############################################################################### 

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine checksum_field(field, total)

      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen, nsegments, ikxmod
      use grids_extended_zgrid, only: iz_low, iz_up

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: field
      real, intent(out) :: total

      integer :: it, iky, ie, iseg
      integer :: ikx

      !----------------------------------------------------------------------

      total = 0.

      do iky = 1, naky
         do it = 1, ntubes
               do ie = 1, neigen(iky)
               iseg = 1
               ikx = ikxmod(iseg, ie, iky)
               total = total + sum(cabs(field(iky, ikx, iz_low(iseg):iz_up(iseg), it)))
               if (nsegments(ie, iky) > 1) then
                  do iseg = 2, nsegments(ie, iky)
                     ikx = ikxmod(iseg, ie, iky)
                     total = total + sum(cabs(field(iky, ikx, iz_low(iseg) + 1:iz_up(iseg), it)))
                  end do
               end if
               end do
         end do
      end do

   end subroutine checksum_field

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine checksum_dist(dist, total, norm)

      use mp, only: sum_allreduce
      use grids_z, only: nzgrid, ntubes
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use grids_kxky, only: naky, nakx
      use grids_velocity, only: maxwell_vpa, maxwell_mu

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: dist
      real, intent(out) :: total
      logical, intent(in), optional :: norm

      integer :: ivmu, iv, imu, is
      integer :: iky, ikx, it
      real :: subtotal

      complex, dimension(:, :, :, :), allocatable :: dist_single

      !----------------------------------------------------------------------

      total = 0.

      allocate (dist_single(naky, nakx, -nzgrid:nzgrid, ntubes))
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         dist_single = dist(:, :, :, :, ivmu)
         if (present(norm)) then
               if (norm) then
               iv = iv_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               is = is_idx(vmu_lo, ivmu)
               do it = 1, ntubes
                  do ikx = 1, nakx
                     do iky = 1, naky
                           dist_single(iky, ikx, :, it) = dist_single(iky, ikx, :, it) * maxwell_vpa(iv, is) * maxwell_mu(1, :, imu, is)
                     end do
                  end do
               end do
               else
               end if
         end if
         call checksum(dist_single, subtotal)
         total = total + subtotal
      end do
      deallocate (dist_single)

      call sum_allreduce(total)

   end subroutine checksum_dist

end module calculations_checksum
