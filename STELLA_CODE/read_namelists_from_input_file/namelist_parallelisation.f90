!###############################################################################
!################## READ STELLA NAMELISTS FOR PARALLELISATION ##################
!###############################################################################
! 
! This module will read the namelists associated with parallelisation:
! 
!   parallelisation
!     xyzs_layout = 'yxzs'
!     vms_layout = 'vms'
!     mat_gen = .false.
!     mat_read = .false.
!     lu_option = 'default'
!     fields_kxkyz = .false.
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
module namelist_parallelisation

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_parallelisation
   
   ! Parameters need to be public
   public :: lu_option_local, lu_option_none, lu_option_global

   private

   ! Create parameters for lu_option
   integer, parameter :: lu_option_none = 1
   integer, parameter :: lu_option_local = 2
   integer, parameter :: lu_option_global = 3
      
   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist
   
contains

   !****************************************************************************
   !                              PARALLELISATION                              !
   !****************************************************************************
   subroutine read_namelist_parallelisation(xyzs_layout, vms_layout, kymus_layout, &
      mat_gen, mat_read, lu_option_switch, fields_kxkyz)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      character(len=4), intent (out) :: xyzs_layout
      character(len=3), intent (out) :: vms_layout
      character(len=5), intent (out) :: kymus_layout
      logical, intent (out)  :: fields_kxkyz, mat_gen, mat_read
      integer, intent (out)  :: lu_option_switch

      ! Local variable to set <lu_option_switch>
      character(30) :: lu_option
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_parallelisation
      call read_input_file_parallelisation
      call check_inputs_parallelisation
      call write_parameters_to_input_file

   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_parallelisation

         implicit none

         xyzs_layout = 'yxzs'
         vms_layout = 'vms'
         kymus_layout = 'kymus'
         fields_kxkyz = .false.
         lu_option = 'default'
         mat_gen = .false.
         mat_read = .false.

      end subroutine set_default_parameters_parallelisation

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_parallelisation

         use file_utils, only: input_unit_exist
         use file_units, only: unit_error_file
         use text_options, only: text_option, get_option_value

         implicit none
         
         ! Link text options for <lu_option> to an integer value
         type(text_option), dimension(3), parameter :: luopts = &
               (/text_option('default', lu_option_none), &
               text_option('local', lu_option_local), &
               text_option('global', lu_option_global)/)

         ! Variables in the <parallelisation> namelist
         namelist /parallelisation/ xyzs_layout, vms_layout, kymus_layout, &
             mat_gen, mat_read, lu_option, fields_kxkyz
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('parallelisation', dexist)
         if (dexist) read (unit=in_file, nml=parallelisation)

         ! Read the text option in <lu_option> and store it in <lu_option_switch>
         call get_option_value(lu_option, luopts, lu_option_switch, &
             unit_error_file, 'lu_option in parameters_numerical')

      end subroutine read_input_file_parallelisation

      !---------------------------- Check inputs --------------------------------
      subroutine check_inputs_parallelisation

         use mp, only: mp_abort

         implicit none

         if (xyzs_layout /= 'xyzs' .and. &
               xyzs_layout /= 'xzys' .and. &
               xyzs_layout /= 'yxzs' .and. &
               xyzs_layout /= 'yzxs' .and. &
               xyzs_layout /= 'zxys' .and. &
               xyzs_layout /= 'zyxs') then
               call mp_abort('parallelisation_layouts: read_parameters finds illegal xyzs_layout. Aborting.')
         end if
         if (vms_layout /= 'vms' .and. &
               vms_layout /= 'mvs') then
               call mp_abort('parallelisation_layouts: read_parameters finds illegal vms_layout. Aborting.')
         end if
         if (kymus_layout /= 'kymus' .and. &
               kymus_layout /= 'mukys') then
               call mp_abort('parallelisation_layouts: read_parameters finds illegal kymus_layout. Aborting.')
         end if

      end subroutine check_inputs_parallelisation
      
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&parallelisation'
         write (unit, '(A, A, A)') '  lu_option = "', trim(lu_option), '"'
         write (unit, '(A, A, A)') '  xyzs_layout = "', trim(xyzs_layout), '"'
         write (unit, '(A, A, A)') '  vms_layout = "', trim(vms_layout), '"'
         write (unit, '(A, A, A)') '  kymus_layout = "', trim(kymus_layout), '"'
         write (unit, '(A, L0)') '  fields_kxkyz = ', fields_kxkyz
         write (unit, '(A, L0)') '  mat_gen = ', mat_gen
         write (unit, '(A, L0)') '  mat_read = ', mat_read
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_parallelisation

end module namelist_parallelisation
