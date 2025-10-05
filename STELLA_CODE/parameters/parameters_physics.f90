!###############################################################################
!############################ READ PHYSICS PARAMETES ###########################
!###############################################################################
! These are different logicals and parameters that adjust the physics in the 
! problem. These will allow you to toggles whether you want to include different
! terms in the gyrokinetic equation, as well as allowing you to include different
! large scale effects such as whether the system allows for electromanetic 
! effect, full flux surface effect, or radially global effects.

! Here, the following namelists are read: 
! &gyrokinetic_terms
! &scale_gyrokinetic_terms
! &electromagnetic
! &flow_shear
! &extra - to be changed 
!###############################################################################

module parameters_physics

   ! Read the parameters for <simulation_domain_switch> from namelist_parameters_physics.f90
   use namelist_parameters_physics, only: simulation_domain_fluxtube
   use namelist_parameters_physics, only: simulation_domain_multibox
   use namelist_parameters_physics, only: simulation_domain_flux_annulus

   implicit none

   ! Public subroutines that are read by the main stella routine.
   public :: read_parameters_physics
   public :: finish_read_parameters_physics

   ! Standard gyrokinetic terms that can be turned on/off
   ! with the following toggles.
   public :: include_parallel_streaming
   public :: include_mirror
   public :: include_xdrift
   public :: include_ydrift
   public :: include_drive
   public :: include_nonlinear
   public :: include_parallel_nonlinearity
   public :: include_electromagnetic
   public :: include_flow_shear
   public :: full_flux_annulus
   public :: radial_variation

   ! Scaling options
   public :: xdriftknob, ydriftknob, wstarknob
   public :: fphi, suppress_zonal_interaction
   
   ! Electromagnetic effects
   public :: include_apar
   public :: include_bpar
   public :: beta 

   ! Full flux annulus effects
   public :: rhostar

   ! Only initialise once
   public :: initialised_parameters_physics

   private

   ! Gyrokinetic terms
   integer :: simulation_domain_switch
   logical :: include_parallel_streaming
   logical :: include_mirror
   logical :: include_xdrift
   logical :: include_ydrift
   logical :: include_drive
   logical :: include_nonlinear
   logical :: include_parallel_nonlinearity
   logical :: include_electromagnetic
   logical :: include_flow_shear
   logical :: full_flux_annulus
   logical :: radial_variation

   ! Scaling options
   real :: xdriftknob, ydriftknob, wstarknob
   real :: fphi
   logical :: suppress_zonal_interaction

   ! Electromagnetic effects
   logical :: include_apar
   logical :: include_bpar
   real :: beta

   ! Full flux annulus effects
   real :: rhostar
   
   ! Only initialise once
   logical :: initialised_parameters_physics = .false.

contains

   !*************************************************************************
   !                         Read physics parameters                         
   !*************************************************************************
   subroutine read_parameters_physics

      use namelist_parameters_physics, only: read_namelist_gyrokinetic_terms
      use namelist_parameters_physics, only: read_namelist_scale_gyrokinetic_terms
      use namelist_parameters_physics, only: read_namelist_electromagnetic
      use namelist_parameters_physics, only: read_namelist_physics_inputs

      implicit none
            
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_parameters_physics) return
      initialised_parameters_physics = .true.

      ! Read the physics namelists in the input file
      call read_namelist_gyrokinetic_terms (simulation_domain_switch, & 
         include_parallel_streaming, include_mirror, &
         include_xdrift, include_ydrift, include_drive, include_nonlinear, &
         include_parallel_nonlinearity, include_electromagnetic, include_flow_shear, &
         full_flux_annulus, radial_variation)
      call read_namelist_scale_gyrokinetic_terms(include_xdrift, include_ydrift, include_drive, & 
         xdriftknob, ydriftknob, wstarknob, fphi, suppress_zonal_interaction)
      call read_namelist_electromagnetic(include_electromagnetic, include_apar, include_bpar, beta) 
      call read_namelist_physics_inputs(rhostar)

      ! Broadcast the input parameters from proc0 to all processors
      call broadcast_parameters

      ! Set the simulation domain
      if (simulation_domain_switch == simulation_domain_fluxtube) then
         full_flux_annulus = .false.
         radial_variation = .false.
      else if (simulation_domain_switch == simulation_domain_multibox) then
         full_flux_annulus = .false.
         radial_variation = .true. ! Full flux is not compatible with multibox
      else if (simulation_domain_switch == simulation_domain_flux_annulus) then
         full_flux_annulus = .true.
         radial_variation = .false.
      else
         write(*,*) "Error: simulation_domain must be 'fluxtube', 'multibox', or 'full_flux_annulus'."
         stop
      end if

   contains
   
      !*************************************************************************
      !                           BROADCAST OPTIONS                            !
      !*************************************************************************
      ! Broadcast these parameters to all the processors - necessary because
      ! the above was only done for the first processor (proc0).
      !*************************************************************************
      subroutine broadcast_parameters

         use mp, only: broadcast

         implicit none
            
         !----------------------------------------------------------------------

         ! Gyrokinetic terms
         call broadcast(simulation_domain_switch)
         call broadcast(include_parallel_streaming)
         call broadcast(include_mirror)
         call broadcast(include_nonlinear)
         call broadcast(include_xdrift)
         call broadcast(include_ydrift)
         call broadcast(include_drive)
         call broadcast(include_parallel_nonlinearity)
         call broadcast(include_electromagnetic)
         call broadcast(include_flow_shear)
         call broadcast(full_flux_annulus)
         call broadcast(radial_variation)

         ! Scaling options
         call broadcast(xdriftknob)
         call broadcast(ydriftknob)
         call broadcast(wstarknob)
         call broadcast(fphi)
         call broadcast(suppress_zonal_interaction)

         ! Electromagnetic effects
         call broadcast(include_apar)
         call broadcast(include_bpar)
         call broadcast(beta)

         ! Full flux annulus effects
         call broadcast(rhostar)

      end subroutine broadcast_parameters

   end subroutine read_parameters_physics

   !****************************************************************************
   !                          FINISH READ PARAMETERS                           !
   !****************************************************************************
   ! Set the initialised flag to be false such that we do not initialise twice.
   !****************************************************************************
   subroutine finish_read_parameters_physics
      implicit none
      initialised_parameters_physics = .false.
   end subroutine finish_read_parameters_physics

end module parameters_physics
