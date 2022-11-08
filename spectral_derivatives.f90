
module spectral_derivatives

   public :: get_dgdy, get_dgdx

   private

   interface get_dgdy
      module procedure get_dgdy_2d
      module procedure get_dgdy_3d
      module procedure get_dgdy_4d
   end interface

   interface get_dgdx
      module procedure get_dgdx_2d
      module procedure get_dgdx_3d
      module procedure get_dgdx_4d
   end interface

contains
   !> compute dg/dy in k-space
   !> accepts g(ky,kx)
   subroutine get_dgdy_2d(g, dgdy)

      use constants, only: zi
      use kt_grids, only: nakx, aky

      implicit none

      complex, dimension(:, :), intent(in) :: g
      complex, dimension(:, :), intent(out) :: dgdy

      dgdy = zi * spread(aky, 2, nakx) * g

   end subroutine get_dgdy_2d

   !> compute dg/dy in k-space
   !> accepts g(ky,kx,z,tube)
   subroutine get_dgdy_3d(g, dgdy)

      use constants, only: zi
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, aky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdy

      integer :: it, iz, ikx

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               dgdy(:, ikx, iz, it) = zi * aky(:) * g(:, ikx, iz, it)
            end do
         end do
      end do

   end subroutine get_dgdy_3d

   !> compute dg/dy in k-space
   !> accepts g(ky,kx,z,tube,(vpa,mu,spec))
   subroutine get_dgdy_4d(g, dgdy)

      use constants, only: zi
      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, aky

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dgdy

      integer :: ivmu, ikx, iz, it

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  dgdy(:, ikx, iz, it, ivmu) = zi * aky(:) * g(:, ikx, iz, it, ivmu)
               end do
            end do
         end do
      end do

   end subroutine get_dgdy_4d

   !> compute dg/dx in k-space
   !> accepts g(ky,kx)
   subroutine get_dgdx_2d(g, dgdx)

      use constants, only: zi
      use kt_grids, only: naky, akx

      implicit none

      complex, dimension(:, :), intent(in) :: g
      complex, dimension(:, :), intent(out) :: dgdx

      dgdx = zi * spread(akx, 1, naky) * g

   end subroutine get_dgdx_2d

   !> compute dg/dx in k-space
   !> accepts g(ky,kx,z,tube)
   subroutine get_dgdx_3d(g, dgdx)

      use constants, only: zi
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: akx, nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdx

      integer :: ikx, iz, it

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               dgdx(:, ikx, iz, it) = zi * akx(ikx) * g(:, ikx, iz, it)
            end do
         end do
      end do

   end subroutine get_dgdx_3d

   !> compute dg/dx in k-space
   !> accepts g(ky,kx,z,tube,(vpa,mu,spec))
   subroutine get_dgdx_4d(g, dgdx)

      use constants, only: zi
      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: akx, nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dgdx

      integer :: ivmu, ikx, iz, it

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  dgdx(:, ikx, iz, it, ivmu) = zi * akx(ikx) * g(:, ikx, iz, it, ivmu)
               end do
            end do
         end do
      end do

   end subroutine get_dgdx_4d

end module spectral_derivatives
