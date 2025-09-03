!###############################################################################
!############### READ STELLA NAMELISTS FOR NUMERICAL PARAMETERS ################
!###############################################################################
! 
! This module will read the namelists associated with numerical parameters:
! 
!   time_trace_options
!     nstep = -1.0
!     tend = -1.0
!     autostop = .true.
!     avail_cpu_time = 10000000000.0
!   
!   time_step
!     delt = 0.03
!     delt_option = 'default'
!     delt_max = -1.0
!     delt_min = 1e-10
!     cfl_cushion_upper = 0.5
!     cfl_cushion_middle = 0.25
!     cfl_cushion_lower = 1e-05
!   
!   numerical_algorithms
!     explicit_algorithm = 'rk3'
!     flip_flop = .false.  
!     stream_implicit = .true.
!     stream_iterative_implicit = .false.
!     stream_matrix_inversion = .false.
!     driftkinetic_implicit = .false.
!     mirror_implicit = .true.
!     mirror_semi_lagrange = .true.
!     mirror_linear_interp = .false.
!     drifts_implicit = .false.
!     fully_implicit = .false.
!     fully_explicit = .false.
!     split_parallel_dynamics = .false.
!     maxwellian_inside_zed_derivative = .false.
!     use_deltaphi_for_response_matrix = .false.
!     maxwellian_normalization = .false.
!   
!   numerical_upwinding_for_derivatives
!     zed_upwind = 0.02
!     vpa_upwind = 0.02
!     time_upwind = 0.02
! 
!   flux_annulus
!     nitt = 1
! 
! Text options for <delt_option>:
!    - Automatic: {check_restart, default}
!    - Manual: {set_by_hand}
! 
! Text options for <explicit_algorithm>:
!    - Runge-Kutta 2: {rk2}
!    - Runge-Kutta 3: {rk3, default}
!    - Runge-Kutta 4: {rk4}
!    - Euler: {euler}
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
module namelist_parameters_numerical

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_time_trace_options
   public :: read_namelist_time_step
   public :: read_namelist_numerical_algorithms
   public :: read_namelist_numerical_upwinding_for_derivatives
   public :: read_namelist_flux_annulus

   ! Parameters need to be public (delt_option)
   public :: delt_option_hand, delt_option_auto
   
   ! Parameters need to be public (explicit_algorithm)
   public :: explicit_algorithm_rk3, explicit_algorithm_rk2
   public :: explicit_algorithm_rk4, explicit_algorithm_euler

   private

   ! Create parameters for <delt_option>
   integer, parameter :: delt_option_hand = 1
   integer, parameter :: delt_option_auto = 2
   
   ! Create parameters for <explicit_algorithm>
   integer, parameter :: explicit_algorithm_rk3 = 1
   integer, parameter :: explicit_algorithm_rk2 = 2
   integer, parameter :: explicit_algorithm_rk4 = 3
   integer, parameter :: explicit_algorithm_euler = 4

   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                             TIME TRACE OPTIONS                            !
   !****************************************************************************
   subroutine read_namelist_time_trace_options(nstep, tend, autostop, avail_cpu_time)

      use mp, only: proc0

      implicit none

      integer, intent (out) :: nstep
      logical, intent (out) :: autostop
      real, intent (out) :: tend
      real, intent (out) :: avail_cpu_time
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_time_trace_options
      call read_input_file_time_trace_options
      call write_parameters_to_input_file

   contains

      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_time_trace_options

         implicit none

         nstep = -1.0
         tend = -1.0
         autostop = .true.
         avail_cpu_time = 10000000000.0

      end subroutine set_default_parameters_time_trace_options

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_time_trace_options

         use file_utils, only: input_unit_exist
         implicit none

         namelist /time_trace_options/ nstep, tend, autostop, avail_cpu_time
         in_file = input_unit_exist('time_trace_options', dexist)
         if (dexist) read (unit=in_file, nml=time_trace_options)

      end subroutine read_input_file_time_trace_options
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&time_trace_options'
         write (unit, '(A, I0)') '  nstep = ', nstep
         write (unit, '(A, F0.4)') '  tend = ', tend
         write (unit, '(A, L0)') '  autostop = ', autostop
         write (unit, '(A, ES0.4)') '  avail_cpu_time = ', avail_cpu_time
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_time_trace_options

   !****************************************************************************
   !                             TIME STEP OPTIONS                             !
   !****************************************************************************
   subroutine read_namelist_time_step(delt, delt_option_switch, &
      delt_max, delt_min, cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      real, intent (out) :: delt
      integer, intent (out) :: delt_option_switch
      real, intent (out) :: delt_max, delt_min
      real, intent (out) :: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower

      ! Local variable to set <delt_option_switch>
      character(30) :: delt_option
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_time_step
      call read_input_file_time_step
      call write_parameters_to_input_file

   contains

      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_time_step

         implicit none

         delt = 0.03
         delt_option = 'default'
         delt_max = -1.0
         delt_min = 1e-10
         cfl_cushion_upper = 0.5
         cfl_cushion_middle = 0.25
         cfl_cushion_lower = 1e-05

      end subroutine set_default_parameters_time_step

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_time_step

         use file_utils, only: input_unit_exist
         use file_units, only: unit_error_file
         use text_options, only: text_option, get_option_value

         implicit none

         ! Link text options for <delt_option> to an integer value
         type(text_option), dimension(3), parameter :: deltopts = &
             (/text_option('default', delt_option_auto), &
             text_option('set_by_hand', delt_option_hand), &
             text_option('check_restart', delt_option_auto)/)
             
         ! Variables in the <time_step> namelist
         namelist /time_step/ delt, delt_option, delt_max, delt_min, &
             cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower

         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('time_step', dexist)
         if (dexist) read (unit=in_file, nml=time_step)

         ! Read the text option in <delt_option> and store it in <delt_option_switch>
         call get_option_value(delt_option, deltopts, delt_option_switch, &
             unit_error_file, 'delt_option in namelist_parameters_numerical.f90')

      end subroutine read_input_file_time_step
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&time_step'
         write (unit, '(A, ES0.4)') '  delt = ', delt
         write (unit, '(A, A, A)') '  delt_option = "', trim(delt_option), '"'
         write (unit, '(A, ES0.4)') '  delt_max = ', delt_max
         write (unit, '(A, ES0.4)') '  delt_min = ', delt_min
         write (unit, '(A, ES0.4)') '  cfl_cushion_upper = ', cfl_cushion_upper
         write (unit, '(A, ES0.4)') '  cfl_cushion_middle = ', cfl_cushion_middle
         write (unit, '(A, ES0.4)') '  cfl_cushion_lower = ', cfl_cushion_lower
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_time_step

   !****************************************************************************
   !                           NUMERICAL ALGORITHMS                            !
   !****************************************************************************
   subroutine read_namelist_numerical_algorithms(explicit_algorithm_switch, flip_flop, &
      stream_implicit, stream_iterative_implicit, stream_matrix_inversion, driftkinetic_implicit, &
      mirror_implicit, mirror_semi_lagrange, mirror_linear_interp, drifts_implicit, &
      fully_implicit, fully_explicit, maxwellian_inside_zed_derivative, &
      use_deltaphi_for_response_matrix, split_parallel_dynamics, maxwellian_normalization)

      use mp, only: proc0

      implicit none
      
      ! Variables that are read from the input file
      integer, intent(out) :: explicit_algorithm_switch
      logical, intent(out) :: stream_implicit, stream_iterative_implicit, stream_matrix_inversion, driftkinetic_implicit
      logical, intent(out) :: mirror_implicit, mirror_semi_lagrange, mirror_linear_interp
      logical, intent(out) :: drifts_implicit, fully_implicit, fully_explicit, flip_flop
      logical, intent(out) :: maxwellian_inside_zed_derivative, use_deltaphi_for_response_matrix
      logical, intent(out) :: split_parallel_dynamics, maxwellian_normalization

      ! Local variable to set <explicit_algorithm_switch>
      character(30) :: explicit_algorithm
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_numerical_algorithms
      call read_input_file_numerical_algorithms
      call write_parameters_to_input_file

   contains

     !------------------------ Default input parameters -----------------------
     subroutine set_default_parameters_numerical_algorithms

         implicit none

         explicit_algorithm = 'rk3'
         flip_flop = .false.
         stream_implicit = .true.
         stream_iterative_implicit = .false.
         stream_matrix_inversion = .false.
         driftkinetic_implicit = .false.
         mirror_implicit = .true.
         mirror_semi_lagrange = .true.
         mirror_linear_interp = .false.
         drifts_implicit = .false.
         fully_implicit = .false.
         fully_explicit = .false.
         maxwellian_inside_zed_derivative = .false.
         use_deltaphi_for_response_matrix = .false.
         split_parallel_dynamics = .false.
         maxwellian_normalization = .false.

     end subroutine set_default_parameters_numerical_algorithms

     !---------------------------- Read input file ----------------------------
     subroutine read_input_file_numerical_algorithms

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr
         
         ! Link text options for <explicit_algorithm> to an integer value
         type(text_option), dimension(5), parameter :: explicitopts = &
           (/text_option('default', explicit_algorithm_rk3), &
           text_option('rk3', explicit_algorithm_rk3), &
           text_option('rk2', explicit_algorithm_rk2), &
           text_option('rk4', explicit_algorithm_rk4), &
           text_option('euler', explicit_algorithm_euler)/)
           
         ! Variables in the <numerical_algorithms> namelist
         namelist /numerical_algorithms/ explicit_algorithm, flip_flop, &
            stream_implicit, stream_iterative_implicit, stream_matrix_inversion, driftkinetic_implicit, &
            mirror_implicit, mirror_semi_lagrange, mirror_linear_interp, drifts_implicit, &
            fully_implicit, fully_explicit, maxwellian_inside_zed_derivative, &
            use_deltaphi_for_response_matrix, split_parallel_dynamics, maxwellian_normalization
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('numerical_algorithms', dexist)
         if (dexist) read (unit=in_file, nml=numerical_algorithms)

         ! Read the text option in <explicit_algorithm> and store it in <explicit_algorithm_switch>
         call get_option_value(explicit_algorithm, explicitopts, explicit_algorithm_switch, &
             ierr, 'explicit_algorithm in parameters_numerical')

     end subroutine read_input_file_numerical_algorithms

      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&numerical_algorithms'
         write (unit, '(A, A, A)') '  explicit_algorithm = "', trim(explicit_algorithm), '"'
         write (unit, '(A, L0)') '  flip_flop = ', flip_flop
         write (unit, '(A, L0)') '  stream_implicit = ', stream_implicit
         write (unit, '(A, L0)') '  stream_iterative_implicit = ', stream_iterative_implicit
         write (unit, '(A, L0)') '  stream_matrix_inversion = ', stream_matrix_inversion
         write (unit, '(A, L0)') '  driftkinetic_implicit = ', driftkinetic_implicit
         write (unit, '(A, L0)') '  mirror_implicit = ', mirror_implicit
         write (unit, '(A, L0)') '  mirror_semi_lagrange = ', mirror_semi_lagrange
         write (unit, '(A, L0)') '  mirror_linear_interp = ', mirror_linear_interp
         write (unit, '(A, L0)') '  drifts_implicit = ', drifts_implicit
         write (unit, '(A, L0)') '  fully_implicit = ', fully_implicit
         write (unit, '(A, L0)') '  fully_explicit = ', fully_explicit
         write (unit, '(A, L0)') '  maxwellian_inside_zed_derivative = ', maxwellian_inside_zed_derivative
         write (unit, '(A, L0)') '  use_deltaphi_for_response_matrix = ', use_deltaphi_for_response_matrix
         write (unit, '(A, L0)') '  split_parallel_dynamics = ', split_parallel_dynamics
         write (unit, '(A, L0)') '  maxwellian_normalization = ', maxwellian_normalization
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_numerical_algorithms

   !****************************************************************************
   !                NUMERICAL UPWINDING FOR DERIVATIVES OPTIONS                !
   !****************************************************************************
   subroutine read_namelist_numerical_upwinding_for_derivatives(time_upwind, zed_upwind, vpa_upwind) 

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      real, intent (out) :: time_upwind, zed_upwind, vpa_upwind
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_numerical_upwinding_for_derivatives
      call read_input_file_numerical_upwinding_for_derivatives
      call write_parameters_to_input_file

   contains
        
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_numerical_upwinding_for_derivatives

         implicit none

         time_upwind = 0.02
         zed_upwind = 0.02
         vpa_upwind = 0.02

      end subroutine set_default_parameters_numerical_upwinding_for_derivatives

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_numerical_upwinding_for_derivatives

         use file_utils, only: input_unit_exist
         implicit none

         namelist /numerical_upwinding_for_derivatives/ time_upwind, zed_upwind, vpa_upwind
         in_file = input_unit_exist('numerical_upwinding_for_derivatives', dexist)
         if (dexist) read (unit=in_file, nml=numerical_upwinding_for_derivatives)

      end subroutine read_input_file_numerical_upwinding_for_derivatives

      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&numerical_upwinding_for_derivatives'
         write (unit, '(A, F0.4)') '  time_upwind = ', time_upwind
         write (unit, '(A, F0.4)') '  zed_upwind = ', zed_upwind
         write (unit, '(A, F0.4)') '  vpa_upwind = ', vpa_upwind
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_numerical_upwinding_for_derivatives

   !****************************************************************************
   !                              FULL FLUX ANNULUS                            !
   !****************************************************************************
   subroutine read_namelist_flux_annulus(nitt) !, field_tol, itt_tol)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      integer, intent (out) :: nitt
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_flux_annulus
      call read_input_file_flux_annulus
      call write_parameters_to_input_file

   contains

      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_flux_annulus

         implicit none

         nitt = 1

      end subroutine set_default_parameters_flux_annulus

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_flux_annulus

         use file_utils, only: input_unit_exist

         implicit none

         namelist /flux_annulus/ nitt
         in_file = input_unit_exist('flux_annulus', dexist)
         if (dexist) read (unit=in_file, nml=flux_annulus)

      end subroutine read_input_file_flux_annulus

      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&flux_annulus'
         write (unit, '(A, I0)') '  nitt = ', nitt
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_flux_annulus

end module namelist_parameters_numerical
