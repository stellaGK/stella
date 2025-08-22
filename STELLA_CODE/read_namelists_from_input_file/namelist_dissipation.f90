!###############################################################################
!#################### READ STELLA NAMELISTS FOR DISSIPATION ####################
!###############################################################################
! 
! This module will set the default input input parameters for each name list,
! and it will read the stella input file per namelist.
! 
! For each namelists two (or three) routines will exist:
!    - set_default_parameters_<namelist>
!    - read_namelist_<namelist>
!    - check_inputs_<namelist>
! 
! First we will set the default input parameters, and then we will overwrite
! any default options with those specified in the input file. Optionally
! we can check if any input variables are clashing with each other.
! 
! Overview of stella namelists for collisions and dissipation:
!   dissipation_and_collisions_options (renamed from &dissipation)
!   collisions_dougherty
!   collisions_fokker_planck (renamed from &collisions_fp)
!   hyper_dissipation
! 
!###############################################################################
module namelist_dissipation

   implicit none

   public :: read_namelist_dissipation
   
   private
   
   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                                DISSIPATION                                !
   !****************************************************************************
   subroutine read_namelist_dissipation(include_collisions, collisions_implicit, collision_model, hyper_dissipation)

      use mp, only: proc0
      
      implicit none

      logical, intent(out) :: include_collisions, collisions_implicit, hyper_dissipation
      character(30), intent(out) :: collision_model

      if (.not. proc0) return
      call set_default_parameters_dissipation
      call read_input_file_dissipation
      call check_inputs_dissipation

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_dissipation

         implicit none

         ! By default we do not include collisions nor dissipation
         include_collisions = .false.
         collisions_implicit = .true.
         hyper_dissipation = .false.

         ! Options: dougherty or fokker-planck
         collision_model = "dougherty"

      end subroutine set_default_parameters_dissipation

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_dissipation

         use file_utils, only: input_unit_exist
         implicit none

         namelist /dissipation/ include_collisions, collisions_implicit, collision_model, hyper_dissipation
         in_file = input_unit_exist("dissipation", dexist)
         if (dexist) read (unit=in_file, nml=dissipation)

      end subroutine read_input_file_dissipation

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_dissipation

         implicit none

         if (.not. include_collisions) collisions_implicit = .false.

      end subroutine check_inputs_dissipation

   end subroutine read_namelist_dissipation

end module namelist_dissipation

