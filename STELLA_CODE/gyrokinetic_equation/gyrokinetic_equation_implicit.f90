!###############################################################################
!###############################################################################
!###############################################################################
! 
! This module ...
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
   !                           IMPLICIT TIME ADVANCE SUBROUTINE
   !****************************************************************************
   subroutine advance_distribution_function_using_implicit_gyrokinetic_terms(istep, phi, apar, bpar, g)

      ! Parallelisation
      use mp, only: proc0
      use job_manage, only: time_message
      use stella_layouts, only: vmu_lo
      use arrays, only: time_gke
      
      ! Fields
      use quasineutrality_equation, only: advance_fields_using_quasineutrality_equation
      use quasineutrality_equation, only: fields_updated

      ! Grids
      use grids_z, only: nzgrid

      ! Physics flags
      use parameters_physics, only: include_parallel_streaming
      use parameters_physics, only: include_mirror
      use parameters_physics, only: radial_variation
      use parameters_physics, only: full_flux_surface
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
            call advance_perp_flow_shear(g)
            fields_updated = .false.
         end if

         if (hyper_dissipation) then
            call advance_hyper_dissipation(g)
            fields_updated = .false.
         end if

         if (collisions_implicit .and. include_collisions) then
            call advance_fields_using_quasineutrality_equation(g, phi, apar, bpar, dist='g')
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
         call advance_fields_using_quasineutrality_equation(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 
         
         ! g^{**} is input
         ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
         if (stream_implicit .and. include_parallel_streaming) then
            call advance_implicit_terms(g, phi, apar, bpar)
            if (radial_variation .or. full_flux_surface) fields_updated = .false.
         end if

         ! Update the fields if not already updated
         call advance_fields_using_quasineutrality_equation(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 
         
      else

         ! Get updated fields corresponding to advanced g
         ! note that hyper-dissipation and mirror advances
         ! depended only on g and so did not need field update
         call advance_fields_using_quasineutrality_equation(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 

         ! g^{**} is input
         ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
         if (stream_implicit .and. include_parallel_streaming) then
            call advance_implicit_terms(g, phi, apar, bpar)
            if (radial_variation .or. full_flux_surface) fields_updated = .false.
         end if

         if (mirror_implicit .and. include_mirror) then
            call advance_mirror_implicit(collisions_implicit, g, apar)
            fields_updated = .false.
         end if

         if (collisions_implicit .and. include_collisions) then
            call advance_fields_using_quasineutrality_equation(g, phi, apar, bpar, dist='g')
            call advance_collisions_implicit(mirror_implicit, phi, apar, bpar, g)
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
      if (proc0) call time_message(.false., time_gke(:, 9), ' implicit')

   end subroutine advance_distribution_function_using_implicit_gyrokinetic_terms

end module gyrokinetic_equation_implicit
