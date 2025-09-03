!###############################################################################
!################ READ STELLA NAMELISTS FOR PHYSICS PARAMETERS #################
!###############################################################################
! 
! This module will read the namelists associated with physics parameters:
! 
!   gyrokinetic_terms
!     simulation_domain = 'fluxtube'
!     include_parallel_streaming = .true.
!     include_mirror = .true.
!     include_xdrift = .true.
!     include_ydrift = .true.
!     include_drive = .true. !! wstarknob = 1.0 or 0.0
!     include_nonlinear = .false.
!     include_parallel_nonlinearity = .false.
!     include_electromagnetic = .false.
!     include_full_flux_annulus = .false.
!     include_radial_variation = .false.
!   
!   scale_gyrokinetic_terms
!     xdriftknob = 1.0
!     ydriftknob = 1.0
!     wstarknob = 1.0
!     fphi = 1.0
!     suppress_zonal_interaction = .false.
! 
!   electromagnetic
!     include_apar = .false.
!     include_bpar = .false.
!     beta = 0.0
!   
!   flow_shear
!     prp_shear_enabled = .false.
!     hammett_flow_shear = .true.
!     g_exb = 0.0
!     g_exbfac = 1.0
!     omprimfac = 1.0
!   
!   physics_inputs
!     rhostar = - 1.0
! 
! Text options for <simulation_domain>:
!    - Flux tube: {default, fluxtube, ft}
!    - Flux flux annulus: {full_flux_annulus, ffa}
!    - Multibox from radial variation: {multibox, radial}
! 
! For each namelists two (or three) routines exist:
!    - set_default_parameters_<namelist>
!    - read_namelist_<namelist>
!    - check_inputs_<namelist>
! 
! First the default input parameters are set, then the default options are
! overwritten with those specified in the input file. Optionally, it is
! checked whether any input variables are clashing with each other.
! 
!###############################################################################
module namelist_parameters_physics

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_gyrokinetic_terms
   public :: read_namelist_scale_gyrokinetic_terms
   public :: read_namelist_electromagnetic
   public :: read_namelist_flow_shear
   public :: read_namelist_physics_inputs

   ! Parameters need to be public (simulation_domain)
   public :: simulation_domain_fluxtube, simulation_domain_multibox, simulation_domain_flux_annulus

   private 

   ! Create parameters for <simulation_domain>
   integer, parameter :: simulation_domain_fluxtube = 1
   integer, parameter :: simulation_domain_multibox = 2
   integer, parameter :: simulation_domain_flux_annulus = 3
   
   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                             GYROKINETIC TERMS                             !
   !****************************************************************************
   subroutine read_namelist_gyrokinetic_terms(simulation_domain_switch, include_parallel_streaming, & 
      include_mirror, include_xdrift, include_ydrift, include_drive, include_nonlinear, &
      include_parallel_nonlinearity, include_electromagnetic, include_flow_shear, &
      include_full_flux_annulus, include_radial_variation)

      use mp, only: proc0

      implicit none
      
      ! Variables that are read from the input file
      integer, intent(out) :: simulation_domain_switch
      logical, intent(out) :: include_parallel_streaming, include_mirror
      logical, intent(out) :: include_xdrift, include_ydrift, include_drive
      logical, intent(out) :: include_nonlinear, include_parallel_nonlinearity
      logical, intent(out) :: include_electromagnetic, include_flow_shear
      logical, intent (out) :: include_full_flux_annulus, include_radial_variation

      ! Local variable to set <simulation_domain_switch>
      character(30) :: simulation_domain
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_gyrokinetic_terms
      call read_input_file_gyrokinetic_terms
      call check_inputs_gyrokinetic_terms

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_gyrokinetic_terms

         implicit none

         ! By default the domain is a flux tube
         simulation_domain = 'fluxtube'

         ! By default all the gyrokinetic terms are included
         include_parallel_streaming = .true.
         include_mirror = .true.
         include_xdrift = .true.
         include_ydrift = .true.
         include_drive = .true.
         include_nonlinear = .true.
         
         ! By default parallel nonlinearity is not included
         include_parallel_nonlinearity = .false.
         
         ! By default electromagnetic and flow shear effects are not included
         include_electromagnetic = .false.
         include_flow_shear = .false.
         include_full_flux_annulus = .false.
         include_radial_variation = .false.

      end subroutine set_default_parameters_gyrokinetic_terms

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_gyrokinetic_terms

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr
         
         ! Link text options for <simulation_domain> to an integer value
         type(text_option), dimension(7), parameter :: simulation_domain_options = &
               (/text_option('default', simulation_domain_fluxtube), &
               text_option('fluxtube', simulation_domain_fluxtube), &
               text_option('ft', simulation_domain_fluxtube), &
               text_option('multibox', simulation_domain_multibox), &
               text_option('radial', simulation_domain_multibox), &
               text_option('full_flux_annulus', simulation_domain_flux_annulus), & 
               text_option('ffa', simulation_domain_flux_annulus)/)
               
         ! Variables in the gyrokinetic_terms namelist
         namelist /gyrokinetic_terms/ simulation_domain, include_parallel_streaming, &
            include_mirror, include_xdrift, include_ydrift, include_drive, &
            include_nonlinear, include_parallel_nonlinearity, include_electromagnetic, &
            include_flow_shear, include_full_flux_annulus, include_radial_variation

         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('gyrokinetic_terms', dexist)
         if (dexist) read (unit=in_file, nml=gyrokinetic_terms)

         ! Read the text option in <simulation_domain> and store it in <simulation_domain_switch>
         ierr = error_unit()
         call get_option_value(simulation_domain, simulation_domain_options, simulation_domain_switch, &
            ierr, 'simulation_domain in namelist_parameters.f90')

      end subroutine read_input_file_gyrokinetic_terms

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_gyrokinetic_terms

         implicit none
         
         if (.not. include_nonlinear) include_parallel_nonlinearity = .false.

      end subroutine check_inputs_gyrokinetic_terms

   end subroutine read_namelist_gyrokinetic_terms

   !****************************************************************************
   !                          SCALE GYROKINETIC TERMS                          !
   !****************************************************************************
   subroutine read_namelist_scale_gyrokinetic_terms(include_xdrift, include_ydrift, include_drive, &
      xdriftknob, ydriftknob, wstarknob, fphi, suppress_zonal_interaction)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      logical, intent (in) :: include_xdrift, include_ydrift, include_drive
      logical, intent (out) :: suppress_zonal_interaction
      real, intent (out) :: xdriftknob, ydriftknob, wstarknob 
      real, intent (out) :: fphi
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_scale_gyrokinetic_terms
      call read_input_file_scale_gyrokinetic_terms
      call check_inputs_scale_gyrokinetic_terms

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_scale_gyrokinetic_terms

         implicit none
         
         ! Scale the electrostatic potential
         fphi = 1.0

         ! Scale the omega_{d,k,s} and omega_{*,k,s} terms in the gyrokinetic equation
         xdriftknob = 1.0
         ydriftknob = 1.0
         wstarknob = 1.0

         ! The zonal modes can be set to zero at every time step to eliminate their effect
         suppress_zonal_interaction = .false.

      end subroutine set_default_parameters_scale_gyrokinetic_terms

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_scale_gyrokinetic_terms

         use file_utils, only: input_unit_exist

         implicit none

         namelist /scale_gyrokinetic_terms/ xdriftknob, ydriftknob, wstarknob, fphi, &
            suppress_zonal_interaction

         in_file = input_unit_exist('scale_gyrokinetic_terms', dexist)
         if (dexist) read (unit=in_file, nml=scale_gyrokinetic_terms)

      end subroutine read_input_file_scale_gyrokinetic_terms

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_scale_gyrokinetic_terms

         implicit none
         
         ! If the corresponding <include> variable is false, the knob is set to 0.0
         if (.not. include_xdrift) xdriftknob = 0.0
         if (.not. include_ydrift) ydriftknob = 0.0
         if (.not. include_drive) wstarknob = 0.0

      end subroutine check_inputs_scale_gyrokinetic_terms

   end subroutine read_namelist_scale_gyrokinetic_terms

   !****************************************************************************
   !                            ELECTROMAGNETIC TERMS                          !
   !****************************************************************************
   subroutine read_namelist_electromagnetic(include_electromagnetic, include_apar, include_bpar, beta)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      logical, intent(in out) :: include_electromagnetic
      logical, intent(out) :: include_apar, include_bpar
      real, intent(out) :: beta
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_electromagnetic
      call read_input_file_electromagnetic
      call check_inputs_electromagnetic

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_electromagnetic

         implicit none
         
         ! If include_electromagnetic is true, apar and bpar are both included
         ! and beta is set to 0.0 (i.e. no electromagnetic effects)
         include_apar = include_electromagnetic
         include_bpar = include_electromagnetic
         beta = 0.0

      end subroutine set_default_parameters_electromagnetic

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_electromagnetic

         use file_utils, only: input_unit_exist

         implicit none

         namelist /electromagnetic/ include_apar, include_bpar, beta
         in_file = input_unit_exist('electromagnetic', dexist)
         if (dexist) read (unit=in_file, nml=electromagnetic)

      end subroutine read_input_file_electromagnetic

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_electromagnetic

         implicit none
         
         ! If both include_apar and include_bpar are false, then include_electromagnetic is set to false
         if (.not. (include_apar .and. include_bpar)) include_electromagnetic = .false.

      end subroutine check_inputs_electromagnetic


   end subroutine read_namelist_electromagnetic

   !****************************************************************************
   !                               FLOW SHEAR TERMS                            !
   !****************************************************************************
   subroutine read_namelist_flow_shear(prp_shear_enabled, hammett_flow_shear, &
      g_exb, g_exbfac, omprimfac)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      logical, intent(out) :: prp_shear_enabled, hammett_flow_shear 
      real :: g_exb, g_exbfac, omprimfac
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_flow_shear
      call read_input_file_flow_shear

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_flow_shear

         implicit none

         prp_shear_enabled = .false.
         hammett_flow_shear = .true.
         g_exb = 0.0
         g_exbfac = 1.0
         omprimfac = 1.0

      end subroutine set_default_parameters_flow_shear

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_flow_shear

         use file_utils, only: input_unit_exist
         implicit none

         namelist /flow_shear/ prp_shear_enabled, hammett_flow_shear, g_exb, g_exbfac, omprimfac
         in_file = input_unit_exist('flow_shear', dexist)
         if (dexist) read (unit=in_file, nml=flow_shear)

      end subroutine read_input_file_flow_shear

   end subroutine read_namelist_flow_shear

   !****************************************************************************
   !                                PHYSICAL INPUTS                            !
   !****************************************************************************
   subroutine read_namelist_physics_inputs(rhostar)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      real, intent (out) :: rhostar
      
      !-------------------------------------------------------------------------
      
      if (.not. proc0) return
      call set_default_parameters_physics_inputs
      call read_input_file_physics_inputs

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_physics_inputs

         implicit none

         rhostar = -1.0

      end subroutine set_default_parameters_physics_inputs

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_physics_inputs

         use file_utils, only: input_unit_exist
         implicit none

         namelist /physics_inputs/ rhostar
         in_file = input_unit_exist('physics_inputs', dexist)
         if (dexist) read (unit=in_file, nml=physics_inputs)

      end subroutine read_input_file_physics_inputs

   end subroutine read_namelist_physics_inputs

end module namelist_parameters_physics
