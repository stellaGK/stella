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
      use timers, only: region_start, region_end
      use parallelisation_layouts, only: vmu_lo
      
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
      use parameters_physics, only: include_apar
      
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

      ! Local variables
      logical :: apply_collisions, apply_mirror, collisions_with_induction

      !-------------------------------------------------------------------------

      ! Start the timer for the implicit part of the solve
      call region_start('implicit_gke')

      ! Reverse the order of operations every time step
      ! as part of alternating direction operator splitting
      ! this is needed to ensure 2nd order accuracy in time
      ! if (mod(istep,2)==0) then
      ! g^{*} (coming from explicit solve) is input
      ! get g^{**}, with g^{**}-g^{*} due to mirror term

      if (rk_step) call mb_communicate(g)

      ! The implicit collisions and the implicit mirror term both work on the velocity-local
      ! distribution <gvmu>. When the collisions are directly followed by the mirror term, the
      ! collisions leave their result in <gvmu> and the mirror term continues from it, which
      ! saves two redistributions. This hand-over is only valid in that order: when the order is
      ! reversed (flip_flop), the mirror term has to scatter <g> itself, and the collisions have
      ! to gather their result back into <g>; otherwise the mirror term acts on a stale <gvmu>
      ! and the collision update is lost on every other time step.
      apply_collisions = collisions_implicit .and. include_collisions
      apply_mirror = mirror_implicit .and. include_mirror
      ! With A_par, the collisions change the parallel current, so they change A_par, and the inductive
      ! electric field -dA_par/dt acts back on the distribution (see add_collisional_induction). This needs
      ! <g> itself after the collisions, so there is no hand-over of <gvmu> to the mirror term.
      collisions_with_induction = apply_collisions .and. include_apar

      if (mod(istep, 2) == 1 .or. .not. flip_flop) then

         if (prp_shear_enabled) then
            call advance_perp_flow_shear(g)
            fields_updated = .false.
         end if

         if (hyper_dissipation) then
            call advance_hyper_dissipation(g)
            fields_updated = .false.
         end if

         ! Collisions followed by the mirror term: hand over <gvmu> (see above)
         if (apply_collisions) then
            call advance_fields(g, phi, apar, bpar, dist='g')
            if (collisions_with_induction) then
               call add_collisional_induction(g, phi, apar, bpar)
            else
               call advance_collisions_implicit(apply_mirror, phi, apar, bpar, g)
            end if
            fields_updated = .false.
         end if

         if (apply_mirror) then
            call advance_mirror_implicit(apply_collisions .and. .not. collisions_with_induction, g, apar)
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

         ! Mirror term followed by the collisions: no hand-over, so the mirror term scatters <g>
         ! itself and the collisions gather their result back into <g>
         if (apply_mirror) then
            call advance_mirror_implicit(.false., g, apar)
            fields_updated = .false.
         end if

         if (apply_collisions) then
            call advance_fields(g, phi, apar, bpar, dist='g')
            if (collisions_with_induction) then
               call add_collisional_induction(g, phi, apar, bpar)
            else
               call advance_collisions_implicit(.false., phi, apar, bpar, g)
            end if
            fields_updated = .false.
         end if

         if (hyper_dissipation) then
            call advance_hyper_dissipation(g)
            fields_updated = .false.
         end if

         if (prp_shear_enabled) then
            call advance_perp_flow_shear(g)
            fields_updated = .false.
         end if

      end if

      ! Stop the timer for the implict part of the solve
      call region_end('implicit_gke')

   end subroutine advance_distribution_function_using_implicit_gyrokinetic_terms

   !****************************************************************************
   !               IMPLICIT COLLISIONS WITH THE INDUCTIVE ELECTRIC FIELD
   !****************************************************************************
   ! With A_par, the time derivative in the gyrokinetic equation acts on
   !     gbar = g + 2 (Z_s/T_s) J_0 vpa A_par F_s,
   ! so the collisions advance gbar: gbar^{n+1} = gbar^{n} + Delta t C[h], and A_par^{n+1} follows from
   ! Ampere's law written for gbar. This is what the explicit time advance does. The implicit collision
   ! operators advance g with A_par fixed instead; A_par was then recomputed from g, so the change of A_par
   ! caused by the collisional change of the parallel current had no inductive electric field. For the
   ! electrons this lets the collisions damp the current with nu_ei, instead of resistively, with
   ! nu_ei (k_perp d_e)^2 / (1 + (k_perp d_e)^2), where d_e is the electron skin depth.
   ! Here the change of g in the implicit collision step is added to gbar (with the A_par of the start of
   ! the step, so C[h] is evaluated with A_par lagged), and A_par^{n+1} is obtained from gbar.
   !****************************************************************************
   subroutine add_collisional_induction(g, phi, apar, bpar)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid
      use field_equations, only: advance_fields
      use field_equations, only: fields_updated
      use calculations_tofrom_ghf, only: gbar_to_g
      use dissipation_and_collisions, only: advance_collisions_implicit

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g

      complex, dimension(:, :, :, :), allocatable :: apar_old

      !-------------------------------------------------------------------------

      ! A_par of the start of the collision step, consistent with the incoming g
      allocate (apar_old, source=apar)

      ! g^{n+1} = g^{n} + Delta t C[h^{n+1}] with A_par fixed
      call advance_collisions_implicit(.false., phi, apar, bpar, g)

      ! gbar^{n+1} = gbar^{n} + (g^{n+1} - g^{n}) = g^{n+1} + 2 (Z_s/T_s) J_0 vpa A_par^{n} F_s
      call gbar_to_g(g, apar_old, -1.0)

      ! A_par^{n+1} (and phi, B_par) from gbar^{n+1}, then back to g^{n+1}
      fields_updated = .false.
      call advance_fields(g, phi, apar, bpar, dist='gbar')
      call gbar_to_g(g, apar, 1.0)
      fields_updated = .false.

      deallocate (apar_old)

   end subroutine add_collisional_induction

end module gyrokinetic_equation_implicit
