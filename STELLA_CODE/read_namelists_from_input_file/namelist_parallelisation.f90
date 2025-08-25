module namelist_parallelisation

   implicit none

   public :: read_namelist_parallelisation
   
   public :: lu_option_local, lu_option_none, lu_option_global

   private

   integer, parameter :: lu_option_none = 1, lu_option_local = 2, &
      lu_option_global = 3
      
   ! Internal variables
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

      character(len=4), intent (out) :: xyzs_layout
      character(len=3), intent (out) :: vms_layout
      character(len=5), intent (out) :: kymus_layout
      

      logical, intent (out)  :: fields_kxkyz, mat_gen, mat_read
      integer, intent (out)  :: lu_option_switch

      character(30) :: lu_option

      if (.not. proc0) return
      call set_default_parameters_parallelisation
      call read_input_file_parallelisation
      call check_inputs_parallelisation

   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_parallelisation

         implicit none

         xyzs_layout = 'yxzs'
         vms_layout = 'vms'
         kymus_layout = 'kymus'

         mat_gen = .false.
         mat_read = .false.
         lu_option = "default"
         fields_kxkyz = .false.

      end subroutine set_default_parameters_parallelisation

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_parallelisation

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         integer :: ierr 
         type(text_option), dimension(3), parameter :: luopts = &
               (/text_option('default', lu_option_none), &
               text_option('local', lu_option_local), &
               text_option('global', lu_option_global)/)

         namelist /parallelisation/ xyzs_layout, vms_layout, kymus_layout, &
               mat_gen, mat_read, lu_option, fields_kxkyz

         in_file = input_unit_exist("parallelisation", dexist)
         if (dexist) read (unit=in_file, nml=parallelisation)

         ierr = error_unit()
         call get_option_value(lu_option, luopts, lu_option_switch, ierr, &
                                 "lu_option in parameters_numerical")

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
               call mp_abort('stella_layouts: read_parameters finds illegal xyzs_layout. aborting')
         end if
         if (vms_layout /= 'vms' .and. &
               vms_layout /= 'mvs') then
               call mp_abort('stella_layouts: read_parameters finds illegal vms_layout. aborting')
         end if
         if (kymus_layout /= 'kymus' .and. &
               kymus_layout /= 'mukys') then
               call mp_abort('stella_layouts: read_parameters finds illegal kymus_layout. aborting')
         end if

      end subroutine check_inputs_parallelisation

   end subroutine read_namelist_parallelisation

end module namelist_parallelisation