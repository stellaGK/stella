!###############################################################################
!############################ READ PHYSICS PARAMETES ###########################
!###############################################################################
! Namelist: &parameters_physics
! These are different logicals and parameters that adjust the physics in the 
! problem. These will allow you to toggles whether you want to include different
! terms in the gyrokinetic equation, as well as allowing you to include different
! large scale effects such as whether the system allows for electromanetic 
! effect, full flux surface effect, or radially global effects.
!###############################################################################

module parameters_physics

   implicit none

   !> Public subroutines that are read by the main stella routine.
   public :: read_parameters_physics
   public :: finish_read_parameters_physics

   !> Available physics options: These are standard gyrokinetic terms that
   !> can be turned on/off with the following toggles.
   public :: include_parallel_streaming
   public :: include_mirror
   public :: nonlinear
   public :: xdriftknob, ydriftknob, wstarknob
   
   !> Adiabatic options: This is used when nspec = 1. The non-kinetic
   !> species (usually electrons) is set to have an adiabatic response.
   !> This can be either the classic adiabatic option, or the modified
   !> adiabatic option (i.e. modified Boltzmann electrons).
   public :: adiabatic_option_switch, adiabatic_option_fieldlineavg
   
   !> Additional physics effects
   public :: prp_shear_enabled
   public :: hammett_flow_shear
   public :: include_pressure_variation
   public :: include_geometric_variation
   public :: include_parallel_nonlinearity
   public :: suppress_zonal_interaction
   
   !> Large scale physics options of the system - e.g. whether we have full flux effects, 
   !> electromagnetic effects, or radially global effects.
   public :: full_flux_surface
   public :: include_apar
   public :: include_bpar
   public :: radial_variation
   
   public :: beta, zeff, tite, nine, rhostar, vnew_ref
   public :: g_exb, g_exbfac, omprimfac

   public :: flip_flop_old
   public :: explicit_option_old
   public :: time_advance_knob_exists
   
   private

   logical :: include_parallel_streaming
   logical :: include_mirror
   logical :: nonlinear
   real :: xdriftknob, ydriftknob, wstarknob
 
   integer :: adiabatic_option_switch
   integer, parameter :: adiabatic_option_periodic = 1, &
                       adiabatic_option_zero = 2, &
                       adiabatic_option_fieldlineavg = 3
 
 
   logical :: prp_shear_enabled
   logical :: hammett_flow_shear 
   logical :: include_pressure_variation 
   logical :: include_geometric_variation
   logical :: include_parallel_nonlinearity
   logical :: suppress_zonal_interaction
   
   logical :: full_flux_surface
   logical :: include_apar
   logical :: include_bpar
   logical :: radial_variation

   real :: beta, zeff, tite, nine, rhostar, irhostar, vnew_ref
   real :: g_exb, g_exbfac, omprimfac
   logical :: initialised = .false.

   !!> Need to fix for the warning messages
   logical :: debug = .false.

   !!> TODO-GA: need to fix somehow
   logical :: flip_flop_old
   character(10) :: explicit_option_old
   logical :: time_advance_knob_exists    
contains

  !======================================================================
  !====================== READ PHYSICS PARAMETERS =======================
  !======================================================================
  subroutine read_parameters_physics

   use mp, only: proc0
   use text_options, only: text_option, get_option_value
   use file_utils, only: input_unit, error_unit, input_unit_exist

   implicit none
   
   character(30) :: adiabatic_option

   if (initialised) return

   if (proc0) call set_default_parameters
   if (proc0) call read_input_file
   call broadcast_parameters

   initialised = .true.

 contains 
   
   !**********************************************************************
   !                        SET DEFAULT PARAMETERS                       !
   !**********************************************************************
   ! If not specified in the input file these are the default options that 
   ! will be set for all parameters under the namelist 
   ! &parameters_physics'.
   !**********************************************************************
   subroutine set_default_parameters

      implicit none 

      !> Standard gyrokinetic terms
      include_parallel_streaming = .true.
      include_mirror = .true.
      nonlinear = .false.
      xdriftknob = 1.0
      ydriftknob = 1.0
      wstarknob = 1.0

      !> If not chose we set adiabatic option to be adiabatic electrons (no modified Boltzmann response)
      adiabatic_option = 'field-line-average-term'

      !> Additional effects that can be included but are not by default
      prp_shear_enabled = .false.
      hammett_flow_shear = .true.
      include_pressure_variation = .false.
      include_geometric_variation = .true.
      include_parallel_nonlinearity = .false.
      suppress_zonal_interaction = .false.
      
      full_flux_surface = .false.
      include_apar = .false.
      include_bpar = .false.
      radial_variation = .false.

      beta = 0.0 ! beta = 8 * pi * p_ref / B_ref^2
      zeff = 1.0
      tite = 1.0
      nine = 1.0
      rhostar = -1.0 ! = m_ref * vt_ref / (e * B_ref * a_ref), with refs in SI
      vnew_ref = -1.0 ! various input options will override this value if it is negative

      !> Zonal flow options -> TODO-HT: how to turn on/off
      g_exb = 0.0
      g_exbfac = 1.0
      omprimfac = 1.0
      irhostar = -1.0

      flip_flop_old = .false.
      explicit_option_old = 'default'
      time_advance_knob_exists = .false. 
   end subroutine set_default_parameters

   !**********************************************************************
   !                         READ INPUT OPTIONS                          !
   !**********************************************************************
   ! Overwrite any default options with those specified in the input file. 
   ! Then change the other parameters consistently.
   !**********************************************************************
   subroutine read_input_file

      use file_utils, only: input_unit_exist, error_unit

      implicit none

      type(text_option), dimension(6), parameter :: adiabaticopts = &
      (/text_option('default', adiabatic_option_fieldlineavg), &
      !> TODO-HT or TODO-GA: sed: adiabatic_option_default -> adiabatic_option_periodic
      text_option('no-field-line-average-term', adiabatic_option_periodic), &
      text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
      text_option('iphi00=0', adiabatic_option_periodic), &
      text_option('iphi00=1', adiabatic_option_periodic), &
      text_option('iphi00=2', adiabatic_option_fieldlineavg)/)

      integer :: ierr, in_file
      logical :: nml_exist

      namelist /parameters_physics/ include_parallel_streaming, include_mirror, nonlinear, &
        xdriftknob, ydriftknob, wstarknob, adiabatic_option, prp_shear_enabled, &
        hammett_flow_shear, include_pressure_variation, include_geometric_variation, &
        include_parallel_nonlinearity, suppress_zonal_interaction, full_flux_surface, &
        include_apar, include_bpar, radial_variation, &
        beta, zeff, tite, nine, rhostar, vnew_ref, &
        g_exb, g_exbfac, omprimfac, irhostar
        
     !> Overwrite the default options with any that are explicitly given in the input file
     !> under the heading '&parameters_physics'
     in_file = input_unit_exist("parameters_physics", nml_exist)
     if (nml_exist) read (unit=in_file, nml=parameters_physics)

     call check_backwards_compatability

     if (irhostar > 0) rhostar = 1./irhostar
     !> Don't allow people to set rhostar when its not full flux                                                                                                                                        
     !> Otherwise phase_shift_angle will be changed in grids_kxky.f90
     if (.not. full_flux_surface) rhostar = 0

     ierr = error_unit()
     call get_option_value &
       (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
         ierr, "adiabatic_option in parameters_physics")

   end subroutine

   !**********************************************************************
   !                    CHECK BACKWARDS COMPATIBILITY                    !
   !**********************************************************************
   ! Make sure stella either runs or aborts old names for variables or
   ! namelists are used
   !**********************************************************************
   subroutine check_backwards_compatability

      use mp, only: mp_abort, broadcast
      use debug_flags, only: const_alpha_geo
      implicit none

      logical :: old_nml_exist
      integer :: in_file

      character(10) :: explicit_option
      logical :: flip_flop
      
      namelist /physics_flags/ full_flux_surface, radial_variation, &
         include_parallel_nonlinearity, include_parallel_streaming, &
         include_mirror, include_apar, include_bpar, nonlinear, &
         include_pressure_variation, include_geometric_variation, &
         adiabatic_option, const_alpha_geo, suppress_zonal_interaction

      namelist /parameters/ beta, zeff, tite, nine, rhostar, vnew_ref, &
         g_exb, g_exbfac, omprimfac, irhostar

      namelist /time_advance_knobs/ xdriftknob, ydriftknob, wstarknob, explicit_option, flip_flop
      
      in_file = input_unit_exist("physics_flags", old_nml_exist)
      if (old_nml_exist) then
         read (unit=in_file, nml=physics_flags)
         if(debug) then 
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(*,*) 'Please change the namelist <phyiscs_flags> in the input file'
            write(*,*) 'to <parameters_physics>. You can inlclude the old flags under'
            write(*,*) 'this new namelist'
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         end if
!         call broadcast(const_alpha_geo) 
         !         write(*,*) "Aborting in parameters_physics.f90.&
         !              The namelist <physics_flags> does not exist. &
         !      Please replace this with the title <parameters_physics>"
         ! call mp_abort('Aborting in parameters_physics.f90.& 
         !      The namelist <physics_flags> does not exist. &
         !      Please replace this with the title <parameters_physics>')
      end if
      in_file = input_unit_exist("parameters", old_nml_exist)
      if (old_nml_exist) then
         read (unit=in_file, nml=parameters)
         if(debug) then
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(*,*) 'Please change the namelist <parameters> in the input file'
            write(*,*) 'to <parameters_physics>. You can inlclude the old flags under'
            write(*,*) 'this new namelist'
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         end if
         ! write(*,*) "Aborting in parameters_physics.f90.&
         !      The namelist <physics_parameters> does not exist. &
         !      Please replace this with the title <parameters_physics>"
         ! call mp_abort("Aborting in parameters_physics.f90.& 
         !      The namelist <physics_parameters> does not exist. &
         !      Please replace this with the title <parameters_physics>")
      end if

      in_file = input_unit_exist("time_advance_knobs", old_nml_exist)
      if (old_nml_exist) then
         time_advance_knob_exists = .true. 
         read(unit=in_file, nml=time_advance_knobs)
         explicit_option_old = explicit_option
         write(*,*) 'explicit_option_old', explicit_option_old
         flip_flop_old = flip_flop
      end if


   end subroutine check_backwards_compatability
    
   !**********************************************************************
   !                         BROADCAST OPTIONS                           !
   !**********************************************************************
   ! Broadcast these parameters to all the processors - necessary because
   ! the above was only done for the first processor (proc0).
   !**********************************************************************
   subroutine broadcast_parameters

     use mp, only: broadcast

     implicit none 

     call broadcast(include_parallel_streaming)
     call broadcast(include_mirror)
     call broadcast(nonlinear)
     call broadcast(xdriftknob)
     call broadcast(ydriftknob)
     call broadcast(wstarknob)

     call broadcast(adiabatic_option_switch)

     call broadcast(prp_shear_enabled)
     call broadcast(hammett_flow_shear) 
     call broadcast(include_pressure_variation)
     call broadcast(include_geometric_variation)
     call broadcast(include_parallel_nonlinearity)
     call broadcast(suppress_zonal_interaction)
     
     call broadcast(full_flux_surface)
     call broadcast(include_apar)
     call broadcast(include_bpar)
     call broadcast(radial_variation)

     call broadcast(beta)
     call broadcast(vnew_ref)
     call broadcast(zeff)
     call broadcast(rhostar)
     call broadcast(tite)
     call broadcast(nine)
     call broadcast(g_exb)
     call broadcast(g_exbfac)
     call broadcast(omprimfac)


     !> TODO-GA: FIX PLEASE
     call broadcast(flip_flop_old)
     call broadcast(explicit_option_old)
     call broadcast(time_advance_knob_exists) 
   end subroutine broadcast_parameters

 end subroutine read_parameters_physics

 !**********************************************************************
 !                      FINISH READ PARAMETERS                         !
 !**********************************************************************
 ! Set the initialised flag to be false such that we do not initialise
 ! twice.
 !> TODO-HT or TODO-GA: sed: initialised -> initialized
 !**********************************************************************
 subroutine finish_read_parameters_physics
   implicit none
   initialised = .false.
 end subroutine finish_read_parameters_physics

end module parameters_physics
