!###############################################################################
!                            ARRAYS CONTAINING CONSTANTS                        
!###############################################################################
! 
! This module will fill the following arrays:
!        - kperp2(ky, kx, alpha, z) = |k_perp|^2
!        - dkperp2dr(ky, kx, alpha, z) 
!        - vperp2(alpha, z, mu) = |v_perp|^2
! 
!###############################################################################
module initialise_arrays

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

      use parameters_physics, only: radial_variation

      implicit none
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_constant_arrays) return
      initialised_constant_arrays = .true.

      ! Allocate and initialise kperp2
      if (debug) write (*, *) 'dist_fn::init_dist_fn::init_array_kperp2'
      call init_array_kperp2
      
      ! Allocate and initialise dkperp2dr
      if (debug) write (*, *) 'dist_fn::init_dist_fn::init_array_dkperp2dr'
      if (radial_variation) call init_array_dkperp2dr

      ! Allocate and initialise vperp2
      if (debug) write (*, *) 'dist_fn::init_dist_fn::init_array_vperp2'
      call init_array_vperp2

   end subroutine init_arrays_vperp_kperp

   !****************************************************************************
   !                            Calculate |k_perp|^2  
   !****************************************************************************
   subroutine init_array_kperp2

      use arrays, only: kperp2
      use geometry, only: grady_dot_grady, gradx_dot_grady, gradx_dot_gradx
      use geometry, only: geo_surf, q_as_x
      use grids_z, only: nzgrid
      use grids_kxky, only: naky, nakx, nalpha
      use grids_kxky, only: akx, aky

      implicit none

      integer :: iky, ikx
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_kperp2) return
      initialised_kperp2 = .true.

      ! Allocate the array
      allocate (kperp2(naky, nakx, nalpha, -nzgrid:nzgrid))

      ! Calculate |k_perp|^2
      do iky = 1, naky

         ! Calculate |k_perp|^2 = kx² |∇x|² + 2 kx ky ∇x . ∇y + ky² |∇y|²
         if (.not. q_as_x) then
            do ikx = 1, nakx
               kperp2(iky, ikx, :, :) = akx(ikx) * akx(ikx) * gradx_dot_gradx &
                   + 2.0 * akx(ikx) * aky(iky) * gradx_dot_grady &
                   + aky(iky) * aky(iky) * grady_dot_grady
            end do
         end if
         
         ! Calculate |k_perp|^2 = q'² kr² |∇r|² + 2 q' kr kα ∇r . ∇α + kα² |∇α|²
         ! See equation A.28 in [2022 - St-Onge - A novel approach to radially global gyrokinetic simulations]
         if (q_as_x) then
            do ikx = 1, nakx
               kperp2(iky, ikx, :, :) = akx(ikx) * akx(ikx) * gradx_dot_gradx * (geo_surf%shat**2) &
                   + 2.0 * akx(ikx) * aky(iky) * gradx_dot_grady * (geo_surf%shat) &
                   + aky(iky) * aky(iky) * grady_dot_grady
            end do
         end if
         
      end do

      ! Ensure kperp2 is positive everywhere (only might go negative if using full-flux-surface due to interpolation)
      ! NB: should really avoid this by using higher resolution when reading in VMEC geometry and then
      ! NB: course-graining if necessary to map onto lower-resolution stella grid
      where (kperp2 < 0.0)
         kperp2 = 0.0
      end where

      ! Make sure |k_perp^2| has the same value on z-points where (kx,ky) modes
      ! are connected to each other through the boundary conditions along z
      call enforce_single_valued_kperp2

   end subroutine init_array_kperp2

   !****************************************************************************
   !                            Calculate dkperp2dr
   !****************************************************************************
   ! This array is only needed for radial variation
   !****************************************************************************
   subroutine init_array_dkperp2dr

      use arrays, only: kperp2, dkperp2dr
      use geometry, only: dgds2dr, dgds21dr, dgds22dr
      use geometry, only: geo_surf, q_as_x
      use grids_z, only: nzgrid
      use grids_kxky, only: naky, nakx, nalpha
      use grids_kxky, only: akx, aky

      implicit none

      integer :: iky, ikx
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_dkperp2dr) return
      initialised_dkperp2dr = .true.

      ! Allocate the array
      allocate (dkperp2dr(naky, nakx, nalpha, -nzgrid:nzgrid))
      
      ! Calculate dkperp2dr
      do iky = 1, naky

         if (q_as_x) then
             do ikx = 1, nakx
                where (kperp2(iky, ikx, :, :) > epsilon(0.0))
                   dkperp2dr(iky, ikx, :, :) = aky(iky) * aky(iky) * dgds2dr &
                        + 2.0 * aky(iky) * akx(ikx) * dgds21dr &
                        + akx(ikx) * akx(ikx) * dgds22dr
                   dkperp2dr(iky, ikx, :, :) = dkperp2dr(iky, ikx, :, :) / kperp2(iky, ikx, :, :)
                  elsewhere
                     dkperp2dr(iky, ikx, :, :) = 0.0
                  end where
             end do
         else
             do ikx = 1, nakx
                where (kperp2(iky, ikx, :, :) > epsilon(0.0))
                   dkperp2dr(iky, ikx, :, :) = aky(iky) * aky(iky) * dgds2dr &
                       + 2.0 * aky(iky) * akx(ikx) * dgds21dr / geo_surf%shat &
                       + akx(ikx) * akx(ikx) * dgds22dr / (geo_surf%shat)**2
                   dkperp2dr(iky, ikx, :, :) = dkperp2dr(iky, ikx, :, :) / kperp2(iky, ikx, :, :)
                  elsewhere
                     dkperp2dr(iky, ikx, :, :) = 0.0
                  end where
             end do
         end if
         
      end do

   end subroutine init_array_dkperp2dr

   !****************************************************************************
   !                            Calculate dkperp2dr
   !****************************************************************************
   ! Make sure |k_perp^2| has the same value on z-points where (kx,ky) modes
   ! are connected to each other through the boundary conditions along z
   !****************************************************************************
   subroutine enforce_single_valued_kperp2

      use arrays, only: kperp2
      use grids_kxky, only: naky, nalpha
      use grids_z, only: nzgrid
      use grids_extended_zgrid, only: neigen, nsegments, ikxmod

      implicit none

      integer :: iky, ie, iseg
      real, dimension(:), allocatable :: tmp
      
      !-------------------------------------------------------------------------

      ! Allocate temporary array
      allocate (tmp(nalpha)); tmp = 0.0

      ! Set connected z-points to their average |k_perp^2| value
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

      ! Deallocate temporary array
      deallocate (tmp)

   end subroutine enforce_single_valued_kperp2

   !****************************************************************************
   !                            Calculate |v_perp|^2
   !****************************************************************************
   subroutine init_array_vperp2

      use geometry, only: bmag
      use grids_z, only: nzgrid
      use grids_velocity, only: vperp2
      use grids_velocity, only: nmu, mu
      use grids_kxky, only: nalpha

      implicit none

      integer :: imu
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_vperp2) return
      initialised_vperp2 = .true.

      ! Allocate the array
      if (.not. allocated(vperp2)) allocate (vperp2(nalpha, -nzgrid:nzgrid, nmu)); vperp2 = 0.

      ! Calculate |v_perp|^2
      do imu = 1, nmu
         vperp2(:, :, imu) = 2.0 * mu(imu) * bmag
      end do

   end subroutine init_array_vperp2

!###############################################################################
!########################## FINALISE CONSTANT ARRAYS ###########################
!###############################################################################

   subroutine finish_arrays_vperp_kperp

      implicit none

      call finish_kperp2
      call finish_vperp2

   end subroutine finish_arrays_vperp_kperp

   !------------------------ Finish kperp2 and dkperp2dr -----------------------
   subroutine finish_kperp2

      use arrays, only: kperp2, dkperp2dr

      implicit none

      if (allocated(kperp2)) deallocate (kperp2)
      if (allocated(dkperp2dr)) deallocate (dkperp2dr)
      initialised_kperp2 = .false.
      initialised_dkperp2dr = .false.

   end subroutine finish_kperp2

   !------------------------------ Finish vperp2 -------------------------------
   subroutine finish_vperp2

      use grids_velocity, only: vperp2

      implicit none

      if (allocated(vperp2)) deallocate (vperp2)
      initialised_vperp2 = .false.

   end subroutine finish_vperp2

end module initialise_arrays
