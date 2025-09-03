!###############################################################################
!############ READ STELLA NAMELISTS FOR RADIAL VARIATION (MULTIBOX) ############
!###############################################################################
! 
! This module will read the namelists associated with radial variation:
! 
!   sources
!     source_option = 'none'
!     nu_krook = 0.05
!     tcorr_source = 0.02
!     ikxmax_source = 1.0
!     krook_odd = .true.
!     exclude_boundary_regions = .false.
!     tcorr_source_qn = 0.0
!     exclude_boundary_regions_qn = .false.
!     from_zero = .true.
!     conserve_momentum = .false.
!     conserve_density = .false.
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
module namelist_sources

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_sources

   ! Parameters need to be public
   public :: source_option_none, source_option_krook, source_option_projection
   
   private
   
   ! Create parameters for <source_option>
   integer, parameter :: source_option_none = 1
   integer, parameter :: source_option_krook = 2
   integer, parameter :: source_option_projection = 3

   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                                  SOURCES                                  !
   !****************************************************************************
   subroutine read_namelist_sources(source_option_switch, nu_krook, tcorr_source, &
      ikxmax_source, krook_odd, exclude_boundary_regions, &
      tcorr_source_qn, exclude_boundary_regions_qn, from_zero, &
      conserve_momentum, conserve_density)

      use mp, only: proc0, mp_abort
      use parameters_physics, only: initialised_parameters_physics
      use grids_kxky, only: initialised_grids_kxky

      implicit none

      ! Variables that are read from the input file
      integer, intent (out) :: source_option_switch, ikxmax_source
      logical, intent (out) :: krook_odd, exclude_boundary_regions
      logical, intent (out) :: exclude_boundary_regions_qn, from_zero
      logical, intent (out) :: conserve_momentum, conserve_density
      real, intent (out) :: nu_krook, tcorr_source
      real, intent (out) :: tcorr_source_qn
      
      ! Local variable to set <source_option>
      character(30) :: source_option
      
      !-------------------------------------------------------------------------
      
      ! Some sources flags are based on <radial_variation> and <periodic_variation>
      ! Therefore, we need to read the physics parameters and (kx,ky) grids first
      if (.not. initialised_parameters_physics) then
         call mp_abort('Initialise physics parameters before reading sources namelists. Aborting.')
      end if
      if (.not. initialised_grids_kxky) then
         call mp_abort('Initialise (kx,ky) grids before reading sources namelists. Aborting.')
      end if

      ! Set the default input parameters and read the input file
      if (.not. proc0) return
      call set_default_parameters_sources
      call read_input_file_sources
      call check_inputs_sources
      call write_parameters_to_input_file

   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_sources

         use parameters_physics, only: radial_variation
         use grids_kxky, only: periodic_variation

         implicit none

         ! Default parameters for the sources
         source_option = 'none'
         nu_krook = 0.05
         tcorr_source = 0.02
         tcorr_source_qn = 0.0
         from_zero = .true.
         conserve_momentum = .false.
         conserve_density = .false.
         exclude_boundary_regions = radial_variation .and. .not. periodic_variation
         
         ! Damp only the odd mode that can affect profiles
         krook_odd = .true.
         
         ! For periodic variation we use kx=0 and kx=1
         ikxmax_source = 1.0
         if (periodic_variation) ikxmax_source = 2 

      end subroutine set_default_parameters_sources

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_sources

         use file_utils, only: input_unit_exist
         use file_units, only: unit_error_file
         use text_options, only: text_option, get_option_value

         implicit none
         
         ! Link text options for <source_option> to an integer value
         type(text_option), dimension(4), parameter :: sourceopts = &
            (/text_option('default', source_option_none), &
              text_option('none', source_option_none), &
              text_option('krook', source_option_krook), &
              text_option('projection', source_option_projection)/)

         ! Variables in the <sources> namelist
         namelist /sources/ source_option, nu_krook, tcorr_source, &
            ikxmax_source, krook_odd, exclude_boundary_regions, &
            tcorr_source_qn, exclude_boundary_regions_qn, from_zero, &
            conserve_momentum, conserve_density
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('sources', dexist)
         if (dexist) read (unit=in_file, nml=sources)

         ! Read the text option in <source_option> and store it in <source_option_switch>
         call get_option_value(source_option, sourceopts, source_option_switch, &
            unit_error_file, 'source_option in sources')

      end subroutine read_input_file_sources

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_sources

         implicit none

         if (tcorr_source_qn < 0) tcorr_source_qn = tcorr_source

      end subroutine check_inputs_sources
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------
         
         ! Don't print values if no sources are selected
         if (source_option_switch == source_option_none) return

         ! Print input parameters
         write (unit, '(A)') '&sources'
         write (unit, '(A, A, A)') '  source_option = "', trim(source_option), '"'
         write (unit, '(A, ES0.4)') '  nu_krook = ', nu_krook
         write (unit, '(A, ES0.4)') '  tcorr_source = ', tcorr_source
         write (unit, '(A, I0)') '  ikxmax_source = ', ikxmax_source
         write (unit, '(A, L0)') '  krook_odd = ', krook_odd
         write (unit, '(A, L0)') '  exclude_boundary_regions = ', exclude_boundary_regions
         write (unit, '(A, ES0.4)') '  tcorr_source_qn = ', tcorr_source_qn
         write (unit, '(A, L0)') '  exclude_boundary_regions_qn = ', exclude_boundary_regions_qn
         write (unit, '(A, L0)') '  from_zero = ', from_zero
         write (unit, '(A, L0)') '  conserve_momentum = ', conserve_momentum
         write (unit, '(A, L0)') '  conserve_density = ', conserve_density
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_sources

end module namelist_sources
