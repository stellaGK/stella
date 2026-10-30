!###############################################################################
!################## IMPLICIT TERMS OF THE GYROKINETIC EQUATION #################
!###############################################################################
! 
! This module adds the implicit gyrokinetic terms to the right-hand-side of the 
! gyrokinetic equation, in order to advance the distribution function in time.
! 
! First advance the distribution function <g> in time using the gyrokinetic equation.
! Next, adance the fields (electrostatic potential <phi>, as well as the electromagnetic
! fields <apar> and <bpar>) in time using the quasi-neutrality condition.
! 
!###############################################################################
module gyrokinetic_equation_implicit

   ! Load debug flags
   use debug_flags, only: debug => time_advance_debug

   implicit none
   
   ! Make routines available to stella.f90
   public :: advance_distribution_function_using_implicit_gyrokinetic_terms

   private

contains

   !****************************************************************************
   !                      IMPLICIT TIME ADVANCE SUBROUTINE                      
   !****************************************************************************
   subroutine advance_distribution_function_using_implicit_gyrokinetic_terms(istep, phi, apar, bpar, g)

      ! Parallelisation
      use mp, only: proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo
      use timers, only: time_gke
      
      ! Fields
      use field_equations, only: advance_fields
      use field_equations, only: fields_updated

      ! Grids
      use grids_z, only: nzgrid

      ! Physics flags
      use parameters_physics, only: include_parallel_streaming
      use parameters_physics, only: include_mirror
      use gk_flow_shear, only: prp_shear_enabled
      
      ! Numerical flags
      use parameters_numerical, only: stream_implicit
      use parameters_numerical, only: mirror_implicit
      use parameters_numerical, only: flip_flop
      
      ! Radial variation
      use parameters_multibox, only: rk_step
      use gk_radial_variation, only: mb_communicate

      ! Collisions
      use dissipation_and_collisions, only: hyper_dissipation
      use dissipation_and_collisions, only: include_collisions
      use dissipation_and_collisions, only: collisions_implicit
      use dissipation_and_collisions, only: advance_collisions_implicit
      
      ! Advance implicit terms of the gyrokinetic equation
      use dissipation_hyper, only: advance_hyper_dissipation
      use gk_implicit_terms, only: advance_implicit_terms
      use gk_mirror, only: advance_mirror_implicit
      use gk_flow_shear, only: advance_perp_flow_shear

      implicit none

      ! Arguments
      integer, intent(in) :: istep
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      
      !-------------------------------------------------------------------------

      ! Start the timer for the implicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 9), ' implicit')

      ! Reverse the order of operations every time step
      ! as part of alternating direction operator splitting
      ! this is needed to ensure 2nd order accuracy in time
      ! if (mod(istep,2)==0) then
      ! g^{*} (coming from explicit solve) is input
      ! get g^{**}, with g^{**}-g^{*} due to mirror term

      if (rk_step) call mb_communicate(g)

      if (mod(istep, 2) == 1 .or. .not. flip_flop) then

         if (prp_shear_enabled) then
            call remap_flow_shear(g, phi, apar, bpar)
            fields_updated = .false.
         end if

         if (hyper_dissipation) then
            call advance_hyper_dissipation(g)
            fields_updated = .false.
         end if

         if (collisions_implicit .and. include_collisions) then
            call advance_fields(g, phi, apar, bpar, dist='g')
            call advance_collisions_implicit(mirror_implicit, phi, apar, bpar, g)
            fields_updated = .false.
         end if

         if (mirror_implicit .and. include_mirror) then
            call advance_mirror_implicit(collisions_implicit, g, apar)
            fields_updated = .false.
         end if

         ! If the distribution function has been updated due to the addition
         ! of any implicit terms of the gyrokinetic equation, we need to use the 
         ! quasi-neutrality condition to update the fields <phi>, <apar> and <bpar>
         call advance_fields(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 
         
         ! g^{**} is input
         ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
         if (stream_implicit .and. include_parallel_streaming) then
            call advance_implicit_terms(g, phi, apar, bpar)
            fields_updated = .false.
         end if

         ! Update the fields if not already updated
         call advance_fields(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 
         
      else

         ! Get updated fields corresponding to advanced g
         ! note that hyper-dissipation and mirror advances
         ! depended only on g and so did not need field update
         call advance_fields(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 

         ! g^{**} is input
         ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
         if (stream_implicit .and. include_parallel_streaming) then
            call advance_implicit_terms(g, phi, apar, bpar)
            fields_updated = .false.            
         end if

         if (mirror_implicit .and. include_mirror) then
            call advance_mirror_implicit(collisions_implicit, g, apar)
            fields_updated = .false.
         end if

         if (collisions_implicit .and. include_collisions) then
            call advance_fields(g, phi, apar, bpar, dist='g')
            call advance_collisions_implicit(mirror_implicit, phi, apar, bpar, g)
            fields_updated = .false.
         end if

         if (hyper_dissipation) then
            call advance_hyper_dissipation(g)
            fields_updated = .false.
         end if

         if (prp_shear_enabled) then
            call remap_flow_shear(g, phi, apar, bpar)
            fields_updated = .false.
         end if

      end if

      ! Stop the timer for the implict part of the solve
      if (proc0) call time_message(.false., time_gke(:, 9), ' implicit')

   end subroutine advance_distribution_function_using_implicit_gyrokinetic_terms

   !****************************************************************************
   !                  Hammett remap in electromagnetic runs (remap_em)
   !****************************************************************************
   ! The Hammett remap moves each k_y by one k_x cell once its drifted k_x has moved half a cell.  With
   ! c = (Z F0/T) 2 v_par v_th and gbar = g + c J0 A_par:
   !   remap_em = 0 (default): shift g.  A_par is then solved from g with the new k_perp, A_par = j/k_perp^2, so it
   !     jumps by k_perp,old^2/k_perp,new^2 (about 2.6 at the smallest k_y and k_x = 0): an impulsive inductive
   !     E_par at every remap.
   !   remap_em = 1: shift gbar, built at the old cell.  A_par is then solved from gbar,
   !     A_par = j(gbar)/(k_perp^2 + skin term), and stays nearly continuous; g is rebuilt from gbar.
   !     See DOCUMENTATION/flow_shear_gbar_remap.md.
   ! Only the k_y shifted on this call go through g -> gbar -> g.  Every other k_y (including the zonal ones and
   ! k_x = k_y = 0, where A_par from g is 0 but A_par from gbar is j/skin) is restored unchanged afterwards.  The
   ! field solve is local in (k_y, k_x), so this is exact.
   ! advance_fields returns at once if fields_updated is set, so the flag is cleared before each solve.
   ! Electrostatic runs (include_apar = .false.) have gbar = g and always take the plain shift.
   !****************************************************************************
   subroutine remap_flow_shear(g, phi, apar, bpar)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid
      use grids_kxky, only: naky
      use parameters_physics, only: include_apar
      use field_equations, only: advance_fields, fields_updated
      use calculations_tofrom_ghf, only: gbar_to_g
      use gk_flow_shear, only: advance_perp_flow_shear, perp_shift_pending, perp_shift_mask, remap_em

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
      complex, dimension(:, :, :, :, :), allocatable :: gkeep
      logical, dimension(:), allocatable :: shifted
      integer :: iky

      if (remap_em == 0 .or. .not. include_apar .or. .not. perp_shift_pending()) then
         call advance_perp_flow_shear(g)
         return
      end if

      ! Remember which k_y move this call, and keep the others untouched
      allocate (shifted(naky))
      call perp_shift_mask(shifted)
      allocate (gkeep, source=g)

      ! Build gbar from g and the A_par solved from g
      fields_updated = .false.
      call advance_fields(g, phi, apar, bpar, dist='g')
      call gbar_to_g(g, apar, -1.0)

      ! Shift gbar (this also advances shift_state)
      call advance_perp_flow_shear(g)

      ! Solve the fields from the shifted gbar and rebuild g
      fields_updated = .false.
      call advance_fields(g, phi, apar, bpar, dist='gbar')
      call gbar_to_g(g, apar, 1.0)

      ! Restore every k_y that was not shifted (the remap leaves them unchanged; the round trip would not)
      do iky = 1, naky
         if (.not. shifted(iky)) g(iky, :, :, :, :) = gkeep(iky, :, :, :, :)
      end do
      deallocate (gkeep, shifted)
      fields_updated = .false.

   end subroutine remap_flow_shear

end module gyrokinetic_equation_implicit
