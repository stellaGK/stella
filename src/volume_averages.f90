module volume_averages

   public :: init_volume_averages, finish_volume_averages
   public :: fieldline_average
   public :: volume_average
   public :: flux_surface_average_ffs

   public :: mode_fac
   
   public :: alpha_average_ffs_realspace

   private

   interface fieldline_average
      module procedure fieldline_average_real
      module procedure fieldline_average_complex
   end interface

!   interface alpha_average_ffs
!      module procedure alpha_average_ffs_nonlocal
!      module procedure alpha_average_ffs_local
!      module procedure alpha_average_ffs_realspace
!   end interface alpha_average_ffs
      
   real, dimension(:), allocatable :: mode_fac
   !> Fourier coefficients in y of the Jacobian;
   !> needed for full flux surface simulations
   complex, dimension(:, :), allocatable :: jacobian_ky

contains

   subroutine init_volume_averages

      use zgrid, only: nzgrid, nztot, delzed
      use kt_grids, only: nalpha, aky, nakx, naky, rho_d_clamped
      use stella_geometry, only: geo_surf, drhodpsi
      use stella_geometry, only: geo_surf, jacob, djacdrho, q_as_x, dVolume
      use physics_flags, only: full_flux_surface, radial_variation

      implicit none

      real :: dqdrho

      if (.not. allocated(mode_fac)) then
         allocate (mode_fac(naky)); mode_fac = 2.0
         if (aky(1) < epsilon(0.)) mode_fac(1) = 1.0
      end if
      if (full_flux_surface) then
         call init_flux_surface_average_ffs
      end if

      dqdrho = geo_surf%shat * geo_surf%qinp / geo_surf%rhoc
      if (.not. allocated(dVolume)) allocate (dVolume(nalpha, nakx, -nzgrid:nzgrid))

      !dVolume contains the volume element jacob, which may vary with x or alpha
      ! NB: dVolume does not contain the factor dx, as this should always be uniform
      dVolume = spread(jacob * spread(delzed, 1, nalpha), 2, nakx)
      if (q_as_x) then
         dVolume = dVolume / (dqdrho * drhodpsi)
      end if

      if (radial_variation) then
         if (q_as_x) then
            dVolume = dVolume * (1.+spread(spread(rho_d_clamped, 1, nalpha), 3, nztot) &
                                 * (spread(djacdrho / jacob, 2, nakx) - geo_surf%d2qdr2 / dqdrho &
                                    + geo_surf%d2psidr2 * drhodpsi))
         else
            dVolume = dVolume * (1.+spread(spread(rho_d_clamped, 1, nalpha), 3, nztot) &
                                 * spread(djacdrho / jacob, 2, nakx))
         end if
      end if

      if (full_flux_surface) then
         !something should go here
      end if

      !avoid the double counting at the zed boundaries
      dVolume(:, :, -nzgrid) = 0.5 * dVolume(:, :, -nzgrid)
      dVolume(:, :, nzgrid) = 0.5 * dVolume(:, :, nzgrid)

   end subroutine init_volume_averages

   subroutine finish_volume_averages

      use stella_geometry, only: dVolume
      use physics_flags, only: full_flux_surface

      implicit none

      if (allocated(mode_fac)) deallocate (mode_fac)
      if (allocated(dVolume)) deallocate (dVolume)
      if (full_flux_surface) then
         if (allocated(jacobian_ky)) deallocate (jacobian_ky)
      end if

   end subroutine finish_volume_averages

   !==============================================
   !============ FIELD LINE AVERAGE ==============
   !==============================================
   subroutine fieldline_average_real(unavg, avg)

      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, naky
      use stella_geometry, only: dl_over_b

      implicit none

      real, dimension(:, :, -nzgrid:, :), intent(in) :: unavg
      real, dimension(:, :), intent(out) :: avg

      integer :: it, ia

      ia = 1

      avg = 0.0
      do it = 1, ntubes
         avg = avg + sum(spread(spread(dl_over_b(ia, :), 1, naky), 2, nakx) * unavg(:, :, :, it), dim=3)
      end do
      avg = avg / real(ntubes)

   end subroutine fieldline_average_real

   subroutine fieldline_average_complex(unavg, avg)

      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, naky
      use stella_geometry, only: dl_over_b

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: unavg
      complex, dimension(:, :), intent(out) :: avg

      integer :: it, ia

      ia = 1

      avg = 0.0
      do it = 1, ntubes
         avg = avg + sum(spread(spread(dl_over_b(ia, :), 1, naky), 2, nakx) * unavg(:, :, :, it), dim=3)
      end do
      avg = avg / real(ntubes)

   end subroutine fieldline_average_complex

   !==============================================
   !============== VOLUME AVERAGE ================
   !==============================================
   subroutine volume_average(unavg, avg)

      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use stella_geometry, only: dl_over_b

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: unavg
      real, intent(out) :: avg

      integer :: iky, ikx, iz, it, ia

      ia = 1

      avg = 0.
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               do iky = 1, naky
                  avg = avg + real(unavg(iky, ikx, iz, it) * conjg(unavg(iky, ikx, iz, it))) * mode_fac(iky) * dl_over_b(ia, iz)
               end do
            end do
         end do
      end do
      avg = avg / real(ntubes)

   end subroutine volume_average

   subroutine init_flux_surface_average_ffs

      use zgrid, only: nzgrid
      use kt_grids, only: naky
      use stella_geometry, only: jacob
      use stella_transforms, only: transform_alpha2kalpha

      implicit none

      integer :: iz

      if (.not. allocated(jacobian_ky)) allocate (jacobian_ky(naky, -nzgrid:nzgrid))

      !> calculate the Fourier coefficients in y of the Jacobian
      !> this is needed in the computation of the flux surface average of phi
      do iz = -nzgrid, nzgrid
         call transform_alpha2kalpha(jacob(:, iz), jacobian_ky(:, iz))
      end do

   end subroutine init_flux_surface_average_ffs

   subroutine flux_surface_average_ffs(no_fsa, fsa)

      use zgrid, only: nzgrid, delzed
      use stella_geometry, only: jacob
      use kt_grids, only: naky, naky_all, nalpha
      use kt_grids, only: dy

      implicit none

      complex, dimension(:, -nzgrid:), intent(in) :: no_fsa
      complex, intent(out) :: fsa

      integer :: iky, ikymod
      real :: area

      ! the normalising factor int dy dz Jacobian
      area = sum(spread(delzed * dy, 1, nalpha) * jacob)

      fsa = 0.0
      ! get contribution from negative ky values
      ! for no_fsa, iky=1 corresponds to -kymax, and iky=naky-1 to -dky
      do iky = 1, naky - 1
         ! jacobian_ky only defined for positive ky values
         ! use reality of the jacobian to fill in negative ky values
         ! i.e., jacobian_ky(-ky) = conjg(jacobian_ky(ky))
         ! ikymod runs from naky down to 2, which corresponds
         ! to ky values in jacobian_ky from kymax down to dky
         ikymod = naky - iky + 1
         ! for each ky, add the integral over zed
         fsa = fsa + sum(delzed * no_fsa(iky, :) * jacobian_ky(ikymod, :))
      end do
      ! get contribution from zero and positive ky values
      ! iky = naky correspond to ky=0 for no_fsa and iky=naky_all to ky=kymax
      do iky = naky, naky_all
         ! ikymod runs from 1 to naky
         ! ikymod = 1 corresponds to ky=0 for jacobian_ky and ikymod=naky to ky=kymax
         ikymod = iky - naky + 1
         ! for each ky, add the integral over zed
         fsa = fsa + sum(delzed * no_fsa(iky, :) * conjg(jacobian_ky(ikymod, :)))
      end do
      ! normalise by the flux surface area
      fsa = fsa / area

   end subroutine flux_surface_average_ffs

   ! subroutine alpha_average_ffs_nonlocal (no_avg, avg) 

   !   use zgrid, only: nzgrid, ntubes
   !   use kt_grids, only: naky, naky_all, nalpha
   !   use stella_geometry, only: jacob
   !   use kt_grids, only: dy, nakx

   !   implicit none 

   !   complex, dimension(:,:,-nzgrid:), intent(in) :: no_avg
   !   complex, dimension(:,-nzgrid:), intent(out) :: avg
   !   real, dimension (-nzgrid:nzgrid) :: norm
   !   integer :: iky, ikymod
   !   integer :: ia 

   !   avg = 0.0
   !   norm = sum(jacob* dy , dim = 1)
     
   !   do iky = 1, naky - 1
   !      ikymod = naky - iky + 1
   !      avg(:,:) = avg(:,:) + no_avg(iky, :, :) &
   !           *spread(jacobian_ky(ikymod, :),1,nakx)
   !   end do

   !   do iky = naky, naky_all
   !      ikymod = iky - naky + 1
   !      avg(:,:) = avg(:,:) + no_avg(iky, :, :) &
   !           *spread(conjg(jacobian_ky(ikymod, :)),1,nakx)
   !   end do
     
   !   avg = avg/spread(norm,1,nakx)
     
   ! end subroutine alpha_average_ffs_nonlocal

   ! subroutine alpha_average_ffs_local (no_avg, avg, iz)

   !   use zgrid, only: nzgrid, ntubes
   !   use kt_grids, only: naky, naky_all, nalpha
   !   use stella_geometry, only: jacob
   !   use kt_grids, only: dy, nakx
   !   use stella_transforms , only: transform_kalpha2alpha

   !   implicit none
     
   !   complex, dimension(naky), intent(in) :: no_avg
   !   complex, intent(out) :: avg
   !   real, dimension (:), allocatable :: no_avg_alpha
   !   real :: norm
   !   integer :: iky, ikymod
   !   integer, intent (in) :: iz
   !   integer :: ia

   !   allocate(no_avg_alpha(nalpha)) ; no_avg_alpha = 0.0

   !   avg = 0.0
   !   norm = sum(jacob(:,iz)* dy)
     
   !   call transform_kalpha2alpha (no_avg, no_avg_alpha)
     
   !   do ia = 1, nalpha
   !      avg = avg + no_avg_alpha(ia)* jacob(ia,iz)
   !   end do

   !   ! do iky = 1, naky - 1
   !   !    ikymod = naky - iky + 1
   !   !    avg = avg + no_avg(iky) &
   !   !         *jacobian_ky(ikymod, iz)
   !   ! end do

   !   ! do iky = naky, naky_all
   !   !    ikymod = iky - naky + 1
   !   !    avg = avg + no_avg(iky) &
   !   !         *conjg(jacobian_ky(ikymod, :))
   !   ! end do

   !   avg = avg/norm

   !   deallocate (no_avg_alpha) 

   ! end subroutine alpha_average_ffs_local

   subroutine alpha_average_ffs_realspace (no_avg, avg, iz) 

     use zgrid, only: nzgrid, ntubes
     use kt_grids, only: naky, naky_all, nalpha
     use stella_geometry, only: jacob
     use kt_grids, only: dy, nakx
     
     use mp, only: proc0

     implicit none

     real, dimension(nalpha), intent(in) :: no_avg
     integer, intent (in) :: iz
     real, intent(out) :: avg
     real :: norm
     integer :: ia
     
     avg = 0.0
     norm = sum(jacob(:,iz)*dy) 

     do ia = 1, nalpha 
        avg = avg + no_avg(ia)* jacob(ia,iz) *dy 
     end do

     avg = avg/norm

   end subroutine alpha_average_ffs_realspace

end module volume_averages
