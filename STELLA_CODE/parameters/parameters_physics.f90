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
! &adiabatic_electron_response
! &electromagnetic
! &flow_shear
! &extra - to be changed 
!###############################################################################

module parameters_physics

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
   public :: full_flux_surface
   public :: radial_variation

   ! Scaling options
   public :: xdriftknob, ydriftknob, wstarknob
   public :: fphi, suppress_zonal_interaction
   
   ! Adiabatic options: This is used when nspec = 1. The non-kinetic
   ! species (usually electrons) is set to have an adiabatic response.
   ! This can be either the classic adiabatic option, or the modified
   ! adiabatic option (i.e. modified Boltzmann electrons).
   public :: adiabatic_option_switch, adiabatic_option_fieldlineavg
   public :: tite, nine

   ! Flow shear physics effects
   public :: prp_shear_enabled
   public :: hammett_flow_shear
   public :: g_exb, g_exbfac, omprimfac 
   
   ! Electromagnetic effects
   public :: include_apar
   public :: include_bpar
   public :: beta 

   ! Full flux annulus effects
   public :: rhostar

   !!!! NEED TO MOVE 
   ! public :: include_pressure_variation 
   ! public :: include_geometric_variation  
   public :: zeff
   public :: vnew_ref   

   private

   ! Gyrokinetic terms
   integer :: simulation_domain_switch
   integer, parameter :: simulation_domain_fluxtube = 1, &
                           simulation_domain_multibox = 2, & 
                           simulation_domain_flux_annulus = 3
   logical :: include_parallel_streaming
   logical :: include_mirror
   logical :: include_xdrift
   logical :: include_ydrift
   logical :: include_drive
   logical :: include_nonlinear
   logical :: include_parallel_nonlinearity
   logical :: include_electromagnetic
   logical :: include_flow_shear
   logical :: full_flux_surface
   logical :: radial_variation

   ! Scaling options
   real :: xdriftknob, ydriftknob, wstarknob
   real :: fphi
   logical :: suppress_zonal_interaction
 
   ! Adiabatic options
   integer :: adiabatic_option_switch
   integer, parameter :: adiabatic_option_periodic = 1, &
                       adiabatic_option_zero = 2, &
                       adiabatic_option_fieldlineavg = 3
   real :: tite, nine
   
   ! Flow shear physics effects
   logical :: prp_shear_enabled
   logical :: hammett_flow_shear 
   real :: g_exb, g_exbfac, omprimfac

   ! Electromagnetic effects
   logical :: include_apar
   logical :: include_bpar
   real :: beta

   ! Full flux annulus effects
   real :: rhostar 

   !!!! NEED TO MOVE ??
   real :: zeff, vnew_ref
   
   logical :: initialised = .false.

contains

  !======================================================================
  !====================== READ PHYSICS PARAMETERS =======================
  !======================================================================
  subroutine read_parameters_physics

   use mp, only: proc0
   use namelist_parameters_physics, only: read_namelist_gyrokinetic_terms, &
      read_namelist_scale_gyrokinetic_terms, read_namelist_adiabatic_electron_response, &
      read_namelist_electromagnetic, read_namelist_flow_shear, read_namelist_physics_inputs

   implicit none

   if (initialised) return

   if (proc0) call read_namelist_gyrokinetic_terms (simulation_domain_switch, & 
      include_parallel_streaming, include_mirror, &
      include_xdrift, include_ydrift, include_drive, include_nonlinear, &
      include_parallel_nonlinearity, include_electromagnetic, include_flow_shear, &
      full_flux_surface, radial_variation)

   if (proc0) call read_namelist_scale_gyrokinetic_terms(include_xdrift, include_ydrift, include_drive, & 
      xdriftknob, ydriftknob, wstarknob, fphi, suppress_zonal_interaction)

   if (proc0) call read_namelist_adiabatic_electron_response(adiabatic_option_switch, tite, nine)

   if (proc0) call read_namelist_flow_shear(prp_shear_enabled, hammett_flow_shear, g_exb, g_exbfac, omprimfac)

   if (proc0) call read_namelist_electromagnetic(include_electromagnetic, include_apar, include_bpar, beta) 

   if (proc0) call read_namelist_physics_inputs(rhostar, zeff, vnew_ref)

   call broadcast_parameters
   
   if (simulation_domain_switch == simulation_domain_fluxtube) then
      full_flux_surface = .false.
      radial_variation = .false.
   else if (simulation_domain_switch == simulation_domain_multibox) then
      ! Full flux is not compatible with multibox
      full_flux_surface = .false.  
      radial_variation = .true.
   else if (simulation_domain_switch == simulation_domain_flux_annulus) then
      full_flux_surface = .true. 
      radial_variation = .false.
   else
      write(*,*) "Error: simulation_domain must be 'fluxtube', 'multibox', or 'full_flux_surface'."
      stop
   end if

   initialised = .true.

 contains 
   
   !**********************************************************************
   !                         BROADCAST OPTIONS                           !
   !**********************************************************************
   ! Broadcast these parameters to all the processors - necessary because
   ! the above was only done for the first processor (proc0).
   !**********************************************************************
   subroutine broadcast_parameters

      use mp, only: broadcast

      implicit none 

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
      call broadcast(full_flux_surface)
      call broadcast(radial_variation)

      ! Scaling options
      call broadcast(xdriftknob)
      call broadcast(ydriftknob)
      call broadcast(wstarknob)
      call broadcast(fphi)
      call broadcast(suppress_zonal_interaction)

      ! Adiabatic options
      call broadcast(adiabatic_option_switch)
      call broadcast(tite)
      call broadcast(nine)

      ! Flow shear physics effects
      call broadcast(prp_shear_enabled)
      call broadcast(hammett_flow_shear) 
      call broadcast(g_exb)
      call broadcast(g_exbfac)
      call broadcast(omprimfac)

      ! Electromagnetic effects
      call broadcast(include_apar)
      call broadcast(include_bpar)
      call broadcast(beta)

      ! Full flux annulus effects
      call broadcast(rhostar)

      ! EXTRA - NEED TO MOVE 
      ! call broadcast(include_pressure_variation)
      ! call broadcast(include_geometric_variation)
      call broadcast(vnew_ref)
      call broadcast(zeff)

   end subroutine broadcast_parameters

 end subroutine read_parameters_physics

 !**********************************************************************
 !                      FINISH READ PARAMETERS                         !
 !**********************************************************************
 ! Set the initialised flag to be false such that we do not initialise
 ! twice.
 !**********************************************************************
 subroutine finish_read_parameters_physics
   implicit none
   initialised = .false.
 end subroutine finish_read_parameters_physics

end module parameters_physics
