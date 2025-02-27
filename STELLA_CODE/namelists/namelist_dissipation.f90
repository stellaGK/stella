!###############################################################################
!########################## READ DISSIPATION NAMELIST ##########################
!###############################################################################
! Read the <dissipation> namelist in the input file
!###############################################################################
module namelist_dissipation

   implicit none

   public :: read_namelist
   
   private

contains

   subroutine read_namelist(include_collisions, collisions_implicit, collision_model, hyper_dissipation)

      use mp, only: proc0

      implicit none

      logical, intent(out) :: include_collisions, collisions_implicit, hyper_dissipation
      character(30), intent(out) :: collision_model

      if (.not. proc0) return
      call set_default_parameters
      call read_input_file
      call check_inputs

   contains

      !**********************************************************************
      !                          READ INPUT FILE                            !
      !**********************************************************************
      ! Overwrite any default options with those specified in the input file. 
      ! Then change the other parameters consistently.
      !**********************************************************************
      subroutine read_input_file

         use file_utils, only: input_unit_exist

         implicit none

         integer :: in_file
         logical :: dexist

         ! Variables in the <dissipation> namelist
         namelist /dissipation/ include_collisions, collisions_implicit, collision_model, hyper_dissipation

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist("dissipation", dexist)
         if (dexist) read (unit=in_file, nml=dissipation)

      end subroutine read_input_file


      !**********************************************************************
      !                           CHECK INPUTS                             !
      !**********************************************************************
      ! Make sure that the input variables are set correctly.
      !**********************************************************************
      subroutine check_inputs

         implicit none

         if (.not. include_collisions) collisions_implicit = .false.

      end subroutine check_inputs


      !**********************************************************************
      !                         DEFAULT PARAMETERS                          !
      !**********************************************************************
      ! Set the default input parameters.
      !**********************************************************************
      subroutine set_default_parameters()

         implicit none

         include_collisions = .false.
         collisions_implicit = .true.
         hyper_dissipation = .false.

         ! dougherty or fokker-planck
         collision_model = "dougherty" 

      end subroutine set_default_parameters

   end subroutine read_namelist

end module namelist_dissipation

