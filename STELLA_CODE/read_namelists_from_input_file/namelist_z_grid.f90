!###############################################################################
!####################### READ STELLA NAMELISTS FOR Z-GRID ######################
!###############################################################################
! 
! This module will read the namelists associated with the z-grid:
! 
!   z_grid
!     nzed = 24.0
!     nperiod = 1.0
!     ntubes = 1.0
!     zed_equal_arc = .false.
!   
!   z_boundary_condition
!     boundary_option = 'default'
!     shat_zero = 1e-05
!     grad_x_grad_y_zero = 1e-05
!     dkx_over_dky = -1.0
! 
! Text options for <boundary_option>:
!    - Unconnected (kx,ky) modes: {default, zero, unconnected}
!    - Periodic boundary conditions: {self-periodic, periodic}
!    - Conventional twist-and-shift boundary conditions: {linked}
!    - Stellarator-symmetric twist-and-shift boundary conditions: {stellarator}
! 
! For each namelists two (or three) routines exist:
!    - set_default_parameters_<namelist>
!    - read_namelist_<namelist>
!    - check_inputs_<namelist>
! 
! First the default input parameters are set, then the default options are
! overwritten with those specified in the input file. Optionally, it is
! checked whether any input variables are clashing with each other.
!###############################################################################
module namelist_z_grid

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_z_grid
   public :: read_namelist_z_boundary_condition

   ! Parameters need to be public
   public :: boundary_option_zero, boundary_option_self_periodic
   public :: boundary_option_linked, boundary_option_linked_stellarator

   private

   ! Create parameters for <boundary_option>
   integer, parameter :: boundary_option_zero = 1
   integer, parameter :: boundary_option_self_periodic = 2
   integer, parameter :: boundary_option_linked = 3
   integer, parameter :: boundary_option_linked_stellarator = 4
   
   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                                   Z GRID                                  !
   !****************************************************************************
   subroutine read_namelist_z_grid(nzed, nzgrid, nperiod, ntubes, zed_equal_arc)

      use mp, only: proc0, mp_abort
      use parameters_physics, only: initialised_parameters_physics

      implicit none

      ! Variables that are read from the input file
      integer, intent (out) :: nzed, nzgrid
      integer, intent (out) :: nperiod, ntubes
      logical, intent (out) :: zed_equal_arc
      
      !-------------------------------------------------------------------------
      
      ! Only read parameters on the first processor
      if (.not. proc0) return
      
      ! Some z-grid flags are based on <full_flux_surface>
      ! Therefore, we need to read the physics parameters first
      if (.not. initialised_parameters_physics) then
         call mp_abort('Initialise physics parameters before reading z-grid namelists. Aborting.')
      end if

      ! Set the default input parameters and read the input file
      call set_default_z_grid
      call read_input_file_z_grid
      call check_inputs_z_grid
      call calculate_z_grid_parameters

   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine set_default_z_grid

         implicit none 

         ! Default z-grid variables
         nzed = 24
         nperiod = 1
         ntubes = 1
         
         ! If <zed_equal_arc> = True, then zed is chosen to be arc lengthd
         ! If <zed_equal_arc> = False, then 
         !     - zed is the poloidal angle when using Miller geometry (axisymmetric)
         !     - zed is the toroidal angle when using VMEC geometry
         zed_equal_arc = .false.
         
         ! The following variable is calculated based on <nzed>
         nzgrid = -1

      end subroutine set_default_z_grid

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_z_grid

         use file_utils, only: input_unit_exist

         implicit none

         namelist /z_grid/ nzed, nperiod, ntubes, zed_equal_arc
         in_file = input_unit_exist('z_grid', dexist)
         if (dexist) read (unit=in_file, nml=z_grid)

      end subroutine read_input_file_z_grid
      
      !-------------------------- Calculate parameters -------------------------
      subroutine calculate_z_grid_parameters

         implicit none
         
         ! While <nzed> specifies the grid points left and right of z=0 in a 
         ! single segment, <nzgrid> is the total number of grid points
         nzgrid = nzed / 2 + (nperiod - 1) * nzed
      
      end subroutine calculate_z_grid_parameters

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_z_grid

         use parameters_physics, only: full_flux_surface

         implicit none

         ! Make sure <nzed> is an even integer, otherwise the potential explodes
         if (MOD(nzed,2) .eq. 1) nzed = nzed + 1

         ! In full-flux-surface simulations, force the use of z = arc_length to
         ! ensure that gradpar is alpha-independent, which is necessary to obtain
         ! efficient numerical solution of parallel streaming
         if (full_flux_surface) zed_equal_arc = .true.
      
      end subroutine check_inputs_z_grid

   end subroutine read_namelist_z_grid

   !****************************************************************************
   !                           Z GRID BOUNDARY OPTIONS                         !
   !****************************************************************************
   subroutine read_namelist_z_boundary_condition (boundary_option_switch, shat_zero, &
      grad_x_grad_y_zero, dkx_over_dky) ! GA: put rhostar here?

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      integer, intent (out) :: boundary_option_switch
      real, intent (out) :: shat_zero, grad_x_grad_y_zero, dkx_over_dky

      ! Local variable to set <boundary_option_switch>
      character (20) :: boundary_option
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_z_boundary_condition
      call read_input_file_z_boundary_condition

   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine set_default_z_boundary_condition

         implicit none
         boundary_option = 'default'
         
         ! Set minimum shat value, below this value we assume periodic BC
         shat_zero = 1.e-5
      
         ! Set the ideal ratio between dkx and dky, assuming jtwist = 1, for stellarator-symmetric BC
         ! If it is < 0, the code will just use the nfield_periods in the input file
         dkx_over_dky = -1
         
         ! Set the minimum <nabla x . nabla y> value at the end of the FT, below 
         ! this value we assume periodic BC instead of stellarator symmetric BC
         grad_x_grad_y_zero = 1.e-5

      end subroutine set_default_z_boundary_condition

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_z_boundary_condition
      
         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none 

         ! Variables needed to read the input file
         integer :: ierr
         
         ! Link text options for <boundary_option> to an integer value
         type(text_option), dimension(7), parameter :: boundaryopts = &
            (/text_option('default', boundary_option_zero), &
              text_option('zero', boundary_option_zero), &
              text_option('unconnected', boundary_option_zero), &
              text_option('self-periodic', boundary_option_self_periodic), &
              text_option('periodic', boundary_option_self_periodic), &
              text_option('linked', boundary_option_linked), &
              text_option('stellarator', boundary_option_linked_stellarator)/)

         ! Variables in the <z_boundary_condition> namelist
         namelist /z_boundary_condition/ boundary_option, shat_zero, &
               grad_x_grad_y_zero, dkx_over_dky
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('z_boundary_condition', dexist)
         if (dexist) read (unit=in_file, nml=z_boundary_condition)

         ! Read the text option in <boundary_option> and store it in <boundary_option_switch>
         ierr = error_unit()
         call get_option_value(boundary_option, boundaryopts, boundary_option_switch, &
            ierr, 'boundary_option in namelist_z_grid')

      end subroutine read_input_file_z_boundary_condition

   end subroutine read_namelist_z_boundary_condition

end module namelist_z_grid
