!###############################################################################
!                            ARRAYS CONTAINING CONSTANTS                        
!###############################################################################
! 
! This module will fill the following arrays:
!        - kperp2(ky, kx, alpha, z)
!        - dkperp2dr(ky, kx, alpha, z)
!        - vperp2(alpha, z, mu)
! 
!###############################################################################
module arrays_constants

   ! Load debug flags
   use debug_flags, only: debug => dist_fn_debug
  
   implicit none

   ! Make routines available to other modules
   public :: init_arrays_vperp_kperp
   public :: finish_arrays_vperp_kperp

   private

   ! Only initialize the arrays once
   logical :: initialised_constant_arrays = .false.
   logical :: initialised_kperp2 = .false.
   logical :: initialised_dkperp2dr = .false.
   logical :: initialised_vperp2 = .false.

contains

!###############################################################################
!######################### INITIALISE CONSTANT ARRAYS ##########################
!###############################################################################
   subroutine init_arrays_vperp_kperp

      use mp, only: proc0
      use parameters_physics, only: radial_variation
      use stella_layouts, only: init_dist_fn_layouts

      implicit none
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_constant_arrays) return
      initialised_constant_arrays = .true.

      ! allocate and initialise kperp2 and dkperp2dr
      if (debug) write (*, *) 'dist_fn::init_dist_fn::init_array_kperp2'
      call init_array_kperp2
      if (radial_variation) call init_dkperp2dr

      ! allocate and initialise vperp2
      if (debug) write (*, *) 'dist_fn::init_dist_fn::init_array_vperp2'
      call init_array_vperp2

   end subroutine init_arrays_vperp_kperp

   ! init_array_kperp2 allocates and initialises the kperp2 array
   subroutine init_array_kperp2

      use arrays_store_useful, only: kperp2
      use geometry, only: gds2, gds21, gds22
      use geometry, only: geo_surf, q_as_x
      use grids_z, only: nzgrid
      use grids_kxky, only: naky, nakx, nalpha
      use grids_kxky, only: akx, aky, theta0
      use grids_kxky, only: zonal_mode

      implicit none

      integer :: iky, ikx

      if (initialised_kperp2) return
      initialised_kperp2 = .true.

      ! allocate the kperp2 array to contain |k_perp|^2
      allocate (kperp2(naky, nakx, nalpha, -nzgrid:nzgrid))

      do iky = 1, naky
         if (zonal_mode(iky)) then
            do ikx = 1, nakx
               if (q_as_x) then
                  kperp2(iky, ikx, :, :) = akx(ikx) * akx(ikx) * gds22
               else
                  kperp2(iky, ikx, :, :) = akx(ikx) * akx(ikx) * gds22 / (geo_surf%shat**2)
               end if
            end do
         else
            do ikx = 1, nakx
               kperp2(iky, ikx, :, :) = aky(iky) * aky(iky) &
                                        * (gds2 + 2.0 * theta0(iky, ikx) * gds21 &
                                           + theta0(iky, ikx) * theta0(iky, ikx) * gds22)
            end do
         end if
      end do

      ! NB: should really avoid this by using higher resolution when reading in VMEC geometry and then
      ! NB: course-graining if necessary to map onto lower-resolution stella grid
      ! ensure kperp2 is positive everywhere (only might go negative if using full-flux-surface due to interpolation)
      where (kperp2 < 0.0)
         kperp2 = 0.0
      end where

      call enforce_single_valued_kperp2

   end subroutine init_array_kperp2

   ! init_dkperp2dr allocates and initialises the dkperp2dr array, needed for radial variation
   subroutine init_dkperp2dr

      use arrays_store_useful, only: kperp2, dkperp2dr
      use geometry, only: dgds2dr, dgds21dr, dgds22dr
      use geometry, only: geo_surf, q_as_x
      use grids_z, only: nzgrid
      use grids_kxky, only: naky, nakx, nalpha
      use grids_kxky, only: akx, aky, theta0
      use grids_kxky, only: zonal_mode

      implicit none

      integer :: iky, ikx

      if (initialised_dkperp2dr) return
      initialised_dkperp2dr = .true.

      allocate (dkperp2dr(naky, nakx, nalpha, -nzgrid:nzgrid))
      do iky = 1, naky
         if (zonal_mode(iky)) then
            do ikx = 1, nakx
               if (q_as_x) then
                  where (kperp2(iky, ikx, :, :) > epsilon(0.0))
                     dkperp2dr(iky, ikx, :, :) = akx(ikx) * akx(ikx) * dgds22dr / kperp2(iky, ikx, :, :)
                  elsewhere
                     dkperp2dr(iky, ikx, :, :) = 0.0
                  end where
               else
                  where (kperp2(iky, ikx, :, :) > epsilon(0.0))
                     dkperp2dr(iky, ikx, :, :) = akx(ikx) * akx(ikx) * dgds22dr / (geo_surf%shat**2 * kperp2(iky, ikx, :, :))
                  elsewhere
                     dkperp2dr(iky, ikx, :, :) = 0.0
                  end where
               end if
            end do
         else
            do ikx = 1, nakx
               dkperp2dr(iky, ikx, :, :) = aky(iky) * aky(iky) &
                                           * (dgds2dr + 2.0 * theta0(iky, ikx) * dgds21dr &
                                              + theta0(iky, ikx) * theta0(iky, ikx) * dgds22dr)
               dkperp2dr(iky, ikx, :, :) = dkperp2dr(iky, ikx, :, :) / kperp2(iky, ikx, :, :)
               if (any(kperp2(iky, ikx, :, :) < epsilon(0.))) dkperp2dr(iky, ikx, :, :) = 0.
            end do
         end if
      end do

   end subroutine init_dkperp2dr

   subroutine enforce_single_valued_kperp2

      use arrays_store_useful, only: kperp2
      use grids_kxky, only: naky, nalpha
      use grids_z, only: nzgrid
      use grids_extended_zgrid, only: neigen, nsegments, ikxmod

      implicit none

      integer :: iky, ie, iseg
      real, dimension(:), allocatable :: tmp

      allocate (tmp(nalpha)); tmp = 0.0

      do iky = 1, naky
         do ie = 1, neigen(iky)
            if (nsegments(ie, iky) > 1) then
               do iseg = 2, nsegments(ie, iky)
                  tmp = 0.5 * (kperp2(iky, ikxmod(iseg - 1, ie, iky), :, nzgrid) + kperp2(iky, ikxmod(iseg, ie, iky), :, -nzgrid))
                  kperp2(iky, ikxmod(iseg, ie, iky), :, -nzgrid) = tmp
                  kperp2(iky, ikxmod(iseg - 1, ie, iky), :, nzgrid) = tmp
               end do
            end if
         end do
      end do

      deallocate (tmp)

   end subroutine enforce_single_valued_kperp2

   subroutine init_array_vperp2

      use geometry, only: bmag
      use grids_z, only: nzgrid
      use grids_velocity, only: vperp2
      use grids_velocity, only: nmu, mu
      use grids_kxky, only: nalpha

      implicit none

      integer :: imu

      if (initialised_vperp2) return
      initialised_vperp2 = .true.

      if (.not. allocated(vperp2)) allocate (vperp2(nalpha, -nzgrid:nzgrid, nmu)); vperp2 = 0.

      do imu = 1, nmu
         vperp2(:, :, imu) = 2.0 * mu(imu) * bmag
      end do

   end subroutine init_array_vperp2

   subroutine finish_arrays_vperp_kperp

      implicit none

      call finish_kperp2
      call finish_vperp2

   end subroutine finish_arrays_vperp_kperp

   subroutine finish_kperp2

      use arrays_store_useful, only: kperp2, dkperp2dr

      implicit none

      if (allocated(kperp2)) deallocate (kperp2)
      if (allocated(dkperp2dr)) deallocate (dkperp2dr)

      initialised_kperp2 = .false.
      initialised_dkperp2dr = .false.

   end subroutine finish_kperp2

   subroutine finish_vperp2

      use grids_velocity, only: vperp2

      implicit none

      if (allocated(vperp2)) deallocate (vperp2)

      initialised_vperp2 = .false.

   end subroutine finish_vperp2

end module arrays_constants
