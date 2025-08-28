!###############################################################################
!                             (KX,KY) DERIVATIVES                               
!###############################################################################
! 
! This module computes derivatives in Fourier space.
! 
!                                  MATHEMATICS                                  
! 
! The derivates in real space are calulated in Fourier space as,
!        Fourier [ dg/dy ] = i ky g
!        Fourier [ dg/dx ] = i kx g
! 
! Compute d<chi>_theta/dy and d<chi>/dx in (ky,kx) space where <.>_theta is a gyroaverage
!    Fourier [ d<chi>_theta/dy ] = i * ky * J0 * chi
!    Fourier [ d<chi>_theta/dx ] = i * kx * J0 * chi
! 
! There are different routines depending on the size of the input array.
! 
!###############################################################################
module calculations_kxky_derivatives

   implicit none

   ! Make routines available to other modules
   public :: get_dgdy, get_dgdx
   public :: get_dchidy, get_dchidx

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

   interface get_dchidy
      module procedure get_dchidy_4d
      module procedure get_dchidy_2d
   end interface get_dchidy

contains

!###############################################################################
!################ DERIVATIVES OF DISTRIBUTION FUNCTION W.R.T. Y ################
!###############################################################################

   !****************************************************************************
   !                      Fourier [ dg/dy ] = i ky g(ky,kx)                      
   !****************************************************************************
   subroutine get_dgdy_2d(g, dgdy)

      use constants, only: zi
      use grids_kxky, only: nakx
      use grids_kxky, only: aky

      implicit none

      complex, dimension(:, :), intent(in) :: g
      complex, dimension(:, :), intent(out) :: dgdy

      !----------------------------------------------------------------------

      dgdy = zi * spread(aky, 2, nakx) * g

   end subroutine get_dgdy_2d

   !****************************************************************************
   !                  Fourier [ dg/dy ] = i ky g(ky,kx,z,tube)                  
   !****************************************************************************
   subroutine get_dgdy_3d(g, dgdy)

      use constants, only: zi
      use grids_kxky, only: nakx
      use grids_kxky, only: aky
      use grids_z, only: nzgrid, ntubes
      
      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdy

      integer :: it, iz, ikx

      !----------------------------------------------------------------------

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
               dgdy(:, ikx, iz, it) = zi * aky(:) * g(:, ikx, iz, it)
               end do
         end do
      end do

   end subroutine get_dgdy_3d

   !****************************************************************************
   !               Fourier [ dg/dy ] = i ky g(ky,kx,z,tube,ivpamus)             
   !****************************************************************************
   subroutine get_dgdy_4d(g, dgdy)

      use constants, only: zi
      use stella_layouts, only: vmu_lo
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: nakx
      use grids_kxky, only: aky

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dgdy

      integer :: ivmu, ikx, iz, it

      !----------------------------------------------------------------------

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
   
!###############################################################################
!################ DERIVATIVES OF DISTRIBUTION FUNCTION W.R.T. X ################
!###############################################################################

   !****************************************************************************
   !                      Fourier [ dg/dx ] = i kx g(ky,kx)                     
   !****************************************************************************
   subroutine get_dgdx_2d(g, dgdx)

      use constants, only: zi
      use grids_kxky, only: naky
      use grids_kxky, only: akx

      implicit none

      complex, dimension(:, :), intent(in) :: g
      complex, dimension(:, :), intent(out) :: dgdx

      dgdx = zi * spread(akx, 1, naky) * g

   end subroutine get_dgdx_2d

   !****************************************************************************
   !                  Fourier [ dg/dx ] = i kx g(ky,kx,z,tube)                  
   !****************************************************************************
   subroutine get_dgdx_3d(g, dgdx)

      use constants, only: zi
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: nakx
      use grids_kxky, only: akx

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdx

      integer :: ikx, iz, it

      !----------------------------------------------------------------------

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               dgdx(:, ikx, iz, it) = zi * akx(ikx) * g(:, ikx, iz, it)
            end do
         end do
      end do

   end subroutine get_dgdx_3d

   !****************************************************************************
   !               Fourier [ dg/dx ] = i kx g(ky,kx,z,tube,ivpamus)             
   !****************************************************************************
   subroutine get_dgdx_4d(g, dgdx)

      use constants, only: zi
      use stella_layouts, only: vmu_lo
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: nakx
      use grids_kxky, only: akx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dgdx

      integer :: ivmu, ikx, iz, it

      !----------------------------------------------------------------------

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


!###############################################################################
!############## DERIVATIVES OF FIELDS (PHI, APAR, BPAR) W.R.T. Y ###############
!###############################################################################
   
   !****************************************************************************
   !            Fourier [ d<phi>_theta/dy ] = i ky J0 phi(ky,kx,z,tube)         
   !****************************************************************************
   subroutine get_dchidy_4d(phi, apar, bpar, dchidy)

      ! Constants
      use constants, only: zi
      
      ! Parallelisation
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx, imu_idx
      
      ! Flags
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: full_flux_surface
      use parameters_physics, only: fphi
      
      ! Grids
      use grids_species, only: spec
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: vpa, mu
      use grids_kxky, only: nakx, naky, aky
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average
      use calculations_gyro_averages, only: gyro_average_j1
      use arrays_gyro_averages, only: j0_ffs

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dchidy

      ! Local variables
      integer :: ivmu, iv, is, iky, imu
      complex, dimension(:, :, :, :), allocatable :: field, gyro_tmp
      
      !-------------------------------------------------------------------------
      
      ! Allocate temporary arrays
      allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (gyro_tmp(naky, nakx, -nzgrid:nzgrid, ntubes))

      ! Iterate over the (mu,vpa,s) points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         is = is_idx(vmu_lo, ivmu)
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         
         ! Calculate [i ky phi(ky,kx)]
         field = fphi * phi
         if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
         do iky = 1, naky
            field(iky, :, :, :) = zi * aky(iky) * field(iky, :, :, :)
         end do
         
         ! Take the gyro-average, i.e., add the factor J_0
         if (full_flux_surface) then
            call gyro_average(field, dchidy(:, :, :, :, ivmu), j0_ffs(:, :, :, ivmu))
         else
            call gyro_average(field, ivmu, dchidy(:, :, :, :, ivmu))
         end if
         
         ! Add bpar contribution
         if (include_bpar) then
            field = 4.0 * mu(imu) * (spec(is)%tz) * bpar
            do iky = 1, naky
               field(iky, :, :, :) = zi * aky(iky) * field(iky, :, :, :)
            end do
            call gyro_average_j1(field, ivmu, gyro_tmp)
            dchidy(:, :, :, :, ivmu) = dchidy(:, :, :, :, ivmu) + gyro_tmp
         end if
         
      end do

      ! Deallocate temporary arrays
      deallocate (field)
      deallocate (gyro_tmp)

   end subroutine get_dchidy_4d

   !****************************************************************************
   !                Fourier [ d<phi>_theta/dy ] = i ky J0 phi(ky,kx)            
   !**************************************************************************** 
   subroutine get_dchidy_2d(iz, ivmu, phi, apar, bpar, dchidy)

      ! Constants
      use constants, only: zi
      
      ! Parallelisation
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx, imu_idx
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average
      use calculations_gyro_averages, only: gyro_average_j1
      use arrays_gyro_averages, only: j0_ffs
      
      ! Flags
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: full_flux_surface
      use parameters_physics, only: fphi
      
      ! Grids
      use grids_species, only: spec
      use grids_velocity, only: vpa, mu
      use grids_kxky, only: nakx, naky, aky

      implicit none

      ! Arguments
      integer, intent(in) :: ivmu, iz
      complex, dimension(:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :), intent(out) :: dchidy

      ! Local variables
      integer :: iv, is, imu
      complex, dimension(:, :), allocatable :: field, gyro_tmp
      
      !-------------------------------------------------------------------------
      
      ! Allocate temporary arrays
      allocate (field(naky, nakx))
      allocate (gyro_tmp(naky, nakx))

      ! Get the (mu,vpa,s) point
      is = is_idx(vmu_lo, ivmu)
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      
      ! Calculate [i ky phi(ky,kx)]
      field = fphi * phi
      if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
      field = zi * spread(aky, 2, nakx) * field

      ! Take the gyro-average, i.e., add the factor J_0
      if (full_flux_surface) then
         call gyro_average(field, dchidy, j0_ffs(:, :, iz, ivmu))
      else
         call gyro_average(field, iz, ivmu, dchidy)
      end if

      ! Add bpar contribution
      if (include_bpar) then
         field = 4.0 * mu(imu) * (spec(is)%tz) * bpar
         field = zi * spread(aky, 2, nakx) * field
         call gyro_average_j1(field, iz, ivmu, gyro_tmp)
         dchidy = dchidy + gyro_tmp
      end if
      
      ! Deallocate temporary arrays
      deallocate (field)
      deallocate (gyro_tmp)

   end subroutine get_dchidy_2d


!###############################################################################
!############## DERIVATIVES OF FIELDS (PHI, APAR, BPAR) W.R.T. X ###############
!###############################################################################

   !****************************************************************************
   !                Fourier [ d<phi>_theta/dx] = i kx J0 phi(ky,kx)            
   !****************************************************************************
   subroutine get_dchidx(iz, ivmu, phi, apar, bpar, dchidx)

      ! Constants
      use constants, only: zi
      
      ! Parallelisation
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx, imu_idx
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average
      use calculations_gyro_averages, only: gyro_average_j1
      use arrays_gyro_averages, only: j0_ffs
      
      ! Flags
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: full_flux_surface
      use parameters_physics, only: fphi
      
      ! Grids
      use grids_species, only: spec
      use grids_velocity, only: vpa, mu
      use grids_kxky, only: naky, nakx, akx

      implicit none

      ! Arguments
      integer, intent(in) :: ivmu, iz
      complex, dimension(:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :), intent(out) :: dchidx

      ! Local variables
      integer :: iv, is, imu
      complex, dimension(:, :), allocatable :: field, gyro_tmp
      
      !-------------------------------------------------------------------------
      
      ! Allocate temporary arrays
      allocate (field(naky, nakx))
      allocate (gyro_tmp(naky, nakx))

      ! Get the (mu,vpa,s) point
      is = is_idx(vmu_lo, ivmu)
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      
      ! Calculate [i kx phi(ky,kx)]
      field = fphi * phi
      if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
      field = zi * spread(akx, 1, naky) * field

      ! Take the gyro-average, i.e., add the factor J_0
      if (full_flux_surface) then
         call gyro_average(field, dchidx, j0_ffs(:, :, iz, ivmu))
      else
         call gyro_average(field, iz, ivmu, dchidx)
      end if

      ! Add bpar contribution
      if (include_bpar) then
         field = 4 * mu(imu) * (spec(is)%tz) * bpar
         field = zi * spread(akx, 1, naky) * field
         call gyro_average_j1(field, iz, ivmu, gyro_tmp)
         dchidx = dchidx + gyro_tmp
      end if
      
      ! Deallocate temporary arrays
      deallocate (field)
      deallocate (gyro_tmp)

   end subroutine get_dchidx

end module calculations_kxky_derivatives
