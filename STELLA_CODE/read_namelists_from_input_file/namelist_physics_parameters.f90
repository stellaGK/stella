module namelist_parameters_physics

   implicit none

   public :: read_namelist_gyrokinetic_terms
   public :: read_namelist_scale_gyrokinetic_terms
   public :: read_namelist_adiabatic_electron_response
   public :: read_namelist_electromagnetic
   public :: read_namelist_flow_shear
   public :: read_namelist_physics_inputs

   public :: simulation_domain_fluxtube, simulation_domain_multibox, simulation_domain_flux_annulus
   public :: adiabatic_option_periodic, adiabatic_option_zero, adiabatic_option_fieldlineavg

   private 

   integer, parameter :: simulation_domain_fluxtube = 1
   integer, parameter :: simulation_domain_multibox = 2
   integer, parameter :: simulation_domain_flux_annulus = 3

   integer, parameter :: adiabatic_option_periodic = 1, &
                        adiabatic_option_zero = 2, &
                        adiabatic_option_fieldlineavg = 3
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
      
      integer, intent(out) :: simulation_domain_switch
      logical, intent(out) :: include_parallel_streaming, &
         include_mirror, include_xdrift, include_ydrift, include_drive, &
         include_nonlinear, include_parallel_nonlinearity, & 
         include_electromagnetic, include_flow_shear
      logical, intent (out) :: include_full_flux_annulus, include_radial_variation

      character(30) :: simulation_domain

      if (.not. proc0) return
      call set_default_parameters_gyrokinetic_terms
      call read_input_file_gyrokinetic_terms
      call check_inputs_gyrokinetic_terms

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_gyrokinetic_terms

         implicit none

         ! By default the domain is a flux tube
         ! This is the default for all simulations, but can be changed in the input file
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

         integer :: ierr 
         type(text_option), dimension(7), parameter :: simulation_domain_options = &
               (/text_option('default', simulation_domain_fluxtube), &
               text_option('fluxtube', simulation_domain_fluxtube), &
               text_option('ft', simulation_domain_fluxtube), &
               text_option('multibox', simulation_domain_multibox), &
               text_option('radial', simulation_domain_multibox), &
               text_option('full_flux_annulus', simulation_domain_flux_annulus), & 
               text_option('ffa', simulation_domain_flux_annulus)/)

         !----------------------------------------------------------------------
         ! Variables in the gyrokinetic_terms namelist
         namelist /gyrokinetic_terms/ simulation_domain, include_parallel_streaming, &
            include_mirror, include_xdrift, include_ydrift, include_drive, &
            include_nonlinear, include_parallel_nonlinearity, include_electromagnetic, &
            include_flow_shear, include_full_flux_annulus, include_radial_variation

         in_file = input_unit_exist("gyrokinetic_terms", dexist)
         if (dexist) read (unit=in_file, nml=gyrokinetic_terms)

         ! Read the text option in <initialise_distribution> and store it in <init_distribution_switch>
         ierr = error_unit()
         call get_option_value(simulation_domain, simulation_domain_options, simulation_domain_switch, &
            ierr, "simulation_domain in namelist_parameters.f90")

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

      logical, intent (in) :: include_xdrift, include_ydrift, include_drive

      logical, intent (out) :: suppress_zonal_interaction
      real, intent (out) :: xdriftknob, ydriftknob, wstarknob 
      real, intent (out) :: fphi

      if (.not. proc0) return
      call set_default_parameters_scale_gyrokinetic_terms 
      call read_input_file_scale_gyrokinetic_terms

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_scale_gyrokinetic_terms

         implicit none

         ! By default all scaling knobs are set to 1.0 (i.e. no scaling) if they are included
         ! If the corresponding include variable is false, the knob is set to 0.0
         xdriftknob = merge(1.0, 0.0, include_xdrift)
         ydriftknob = merge(1.0, 0.0, include_ydrift)
         wstarknob = merge(1.0, 0.0, include_drive)
         
         fphi = 1.0

         ! By default, suppress zonal interaction is false
         suppress_zonal_interaction = .false.

      end subroutine set_default_parameters_scale_gyrokinetic_terms

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_scale_gyrokinetic_terms

         use file_utils, only: input_unit_exist

         implicit none

         namelist /scale_gyrokinetic_terms/ xdriftknob, ydriftknob, wstarknob, fphi, &
            suppress_zonal_interaction

         in_file = input_unit_exist("scale_gyrokinetic_terms", dexist)
         if (dexist) read (unit=in_file, nml=scale_gyrokinetic_terms)

      end subroutine read_input_file_scale_gyrokinetic_terms

   end subroutine read_namelist_scale_gyrokinetic_terms

   !****************************************************************************
   !                        ADIABATIC ELECTRON RESPONSE                        !
   !****************************************************************************
   subroutine read_namelist_adiabatic_electron_response(adiabatic_option_switch, tite, nine)

      use mp, only: proc0

      implicit none

      integer, intent(out) :: adiabatic_option_switch
      real, intent(out) :: tite, nine

      ! Local option 
      character(30):: adiabatic_option

      if (.not. proc0) return
      call set_default_parameters_adiabatic_electron_response
      call read_input_file_adiabatic_electron_response

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_adiabatic_electron_response

         implicit none

         ! By default if nspec = 1 the adiabatic electron response is modified adiabatic
         adiabatic_option = 'field-line-average-term'
         tite = 1.0
         nine = 1.0

      end subroutine set_default_parameters_adiabatic_electron_response

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_adiabatic_electron_response

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         integer :: ierr

         type(text_option), dimension(6), parameter :: adiabaticopts = &
            (/text_option('default', adiabatic_option_fieldlineavg), &
            text_option('no-field-line-average-term', adiabatic_option_periodic), &
            text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
            text_option('iphi00=0', adiabatic_option_periodic), &
            text_option('iphi00=1', adiabatic_option_periodic), &
            text_option('iphi00=2', adiabatic_option_fieldlineavg)/)

         !----------------------------------------------------------------------
         namelist /adiabatic_electron_response/ adiabatic_option, tite, nine
         in_file = input_unit_exist("adiabatic_electron_response", dexist)
         if (dexist) read (unit=in_file, nml=adiabatic_electron_response)

         ierr = error_unit()
         call get_option_value &
            (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
            ierr, "adiabatic_option in namelist_parameters.f90")

      end subroutine read_input_file_adiabatic_electron_response

   end subroutine read_namelist_adiabatic_electron_response

   !****************************************************************************
   !                            ELECTROMAGNETIC TERMS                          !
   !****************************************************************************
   subroutine read_namelist_electromagnetic(include_electromagnetic, include_apar, include_bpar, beta)

      use mp, only: proc0

      implicit none

      logical, intent(in out) :: include_electromagnetic
      logical, intent(out) :: include_apar, include_bpar
      real, intent(out) :: beta

      if (.not. proc0) return
      call set_default_parameters_electromagnetic
      call read_input_file_electromagnetic
      call check_inputs_electromagnetic

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_electromagnetic

         implicit none
         ! By default, if include_electromagnetic is true then both apar and bpar are included
         ! and beta is set to 0.0 (i.e. no electromagnetic effects)
         ! If include_electromagnetic is true, apar and bpar are both included by default
         include_apar = include_electromagnetic
         include_bpar = include_electromagnetic
         beta = 0.0

      end subroutine set_default_parameters_electromagnetic

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_electromagnetic

         use file_utils, only: input_unit_exist

         implicit none

         namelist /electromagnetic/ include_apar, include_bpar, beta
         in_file = input_unit_exist("electromagnetic", dexist)
         if (dexist) read (unit=in_file, nml=electromagnetic)

      end subroutine read_input_file_electromagnetic

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_electromagnetic

         implicit none
         ! If both include_apar and include_bpar are false, then include_electromagnetic is set to false
         if (.not. (include_apar .and. include_bpar)) include_electromagnetic = .false.
         ! If both include_apar and include_bpar are true, then include_electromagnetic is set to true
         if (include_apar .or. include_bpar) include_electromagnetic = .true.

      end subroutine check_inputs_electromagnetic


   end subroutine read_namelist_electromagnetic

   !****************************************************************************
   !                               FLOW SHEAR TERMS                            !
   !****************************************************************************
   subroutine read_namelist_flow_shear(prp_shear_enabled, hammett_flow_shear, &
   g_exb, g_exbfac, omprimfac)

      use mp, only: proc0

      implicit none

      logical, intent(out) :: prp_shear_enabled, hammett_flow_shear 
      real :: g_exb, g_exbfac, omprimfac

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

         namelist /flow_shear/ prp_shear_enabled, hammett_flow_shear, &
            g_exb, g_exbfac, omprimfac
         in_file = input_unit_exist("flow_shear", dexist)
         if (dexist) read (unit=in_file, nml=flow_shear)

      end subroutine read_input_file_flow_shear

   end subroutine read_namelist_flow_shear

   !****************************************************************************
   !                                PHYSICAL INPUTS                            !
   !****************************************************************************
   subroutine read_namelist_physics_inputs(rhostar, zeff, vnew_ref)

      use mp, only: proc0

      implicit none

      real, intent (out) :: rhostar
      real, intent(out) :: zeff, vnew_ref
      
      if (.not. proc0) return
      call set_default_parameters_physics_inputs
      call read_input_file_physics_inputs

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_physics_inputs

         implicit none
         ! By default
         rhostar = -1.0 ! = m_ref * vt_ref / (e * B_ref * a_ref), with refs in SI
         zeff = 1.0
         vnew_ref = -1.0 ! various input options will override this value if it is negative

      end subroutine set_default_parameters_physics_inputs

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_physics_inputs

         use file_utils, only: input_unit_exist
         implicit none

         namelist /physics_inputs/ rhostar, zeff, vnew_ref
         !, nitt  !, field_tol, itt_tol

         in_file = input_unit_exist("physics_inputs", dexist)
         if (dexist) read (unit=in_file, nml=physics_inputs)

      end subroutine read_input_file_physics_inputs

   end subroutine read_namelist_physics_inputs

end module namelist_parameters_physics