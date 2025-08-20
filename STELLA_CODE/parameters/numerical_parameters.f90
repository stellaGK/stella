!###############################################################################
!########################### READ NUMERICAL PARAMETES ##########################
!###############################################################################
! Namelists: 
! &time_trace_options
! &time_step
! &numerical_algorithms 
! &numerical_upwinding_for_derivatives
! &numerical_extra - to be changed
! These flags will allow you to toggle the algorithm choices in stella.
!###############################################################################
module numerical_parameters

   implicit none

   ! Time trace options
   public :: nstep, tend
   public :: autostop
   public :: avail_cpu_time

   ! Time step options
   public :: delt
   public :: delt_option_switch
   public :: delt_option_hand, delt_option_auto
   public :: delt_max, delt_min
   public :: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower

   ! Numerical algorithms
   public :: explicit_algorithm_switch
   public :: explicit_algorithm_rk3, explicit_algorithm_rk2, explicit_algorithm_rk4, explicit_algorithm_euler
   public :: flip_flop
   public :: stream_implicit, stream_iterative_implicit, stream_matrix_inversion, driftkinetic_implicit
   public :: mirror_implicit, mirror_semi_lagrange, mirror_linear_interp
   public :: drifts_implicit
   public :: fully_implicit, fully_explicit
   public :: maxwellian_inside_zed_derivative, use_deltaphi_for_response_matrix
   public :: split_parallel_dynamics
   public :: maxwellian_normalization

   ! numerical_upwinding_for_derivatives
   public :: time_upwind, time_upwind_plus, time_upwind_minus
   public :: zed_upwind, zed_upwind_plus, zed_upwind_minus
   public :: vpa_upwind

   ! extra - need to move
   public :: nitt
   public :: fphi
   public :: ky_solve_radial, ky_solve_real
   public :: fields_kxkyz, mat_gen, mat_read
   public :: lu_option_switch
   public :: lu_option_local, lu_option_none, lu_option_global
   public :: rng_seed
   public :: print_extra_info_to_terminal

   ! Public subroutines that are read by the main stella routine.
   public :: read_numerical_parameters, finish_read_numerical_parameters

   !------------------------------------------------------------------
   private
   
   ! Time trace options
   integer :: nstep 
   real :: tend
   real :: avail_cpu_time
   logical :: autostop

   ! Time step options
   real :: delt
   integer :: delt_option_switch
   integer, parameter :: delt_option_hand = 1, delt_option_auto = 2
   real :: delt_max, delt_min
   real :: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower

   ! Numerical algorithms
   integer :: explicit_algorithm_switch
   integer, parameter :: explicit_algorithm_rk3 = 1, &
        explicit_algorithm_rk2 = 2, &
        explicit_algorithm_rk4 = 3, &
        explicit_algorithm_euler = 4
   logical :: flip_flop
   logical :: stream_implicit, stream_iterative_implicit, stream_matrix_inversion, driftkinetic_implicit
   logical :: mirror_implicit, mirror_semi_lagrange, mirror_linear_interp
   logical :: drifts_implicit, fully_implicit, fully_explicit
   logical :: maxwellian_inside_zed_derivative, use_deltaphi_for_response_matrix
   ! if split_parallel_dynamics = .true. (default), use operator splitting
   ! to treat parallel streaming and mirror term separately
   logical :: split_parallel_dynamics
   ! REMOVE
   logical :: maxwellian_normalization

   ! numerical_upwinding_for_derivatives
   real :: time_upwind, time_upwind_plus, time_upwind_minus
   real :: zed_upwind, zed_upwind_plus, zed_upwind_minus
   real :: vpa_upwind

   ! Extra - need to move
   integer :: nitt
   real :: fphi
   logical :: ky_solve_real
   integer :: ky_solve_radial
   logical :: fields_kxkyz, mat_gen, mat_read
   integer :: lu_option_switch
   integer, parameter :: lu_option_none = 1, &
        lu_option_local = 2, &
        lu_option_global = 3
   integer :: rng_seed
   logical :: print_extra_info_to_terminal

   ! Internal
   logical :: initialised = .false.
   logical :: old_nml_exist
   
contains

  !======================================================================
  !===================== READ NUMERICAL PARAMETERS ======================
  !======================================================================
   subroutine read_numerical_parameters

      use mp, only: proc0, mp_abort
      use text_options, only: text_option, get_option_value
      use file_utils, only: input_unit, error_unit, input_unit_exist
      use input_file_numerical_parameters, only: &
         read_namelist_time_trace_options, read_namelist_time_step, &
         read_namelist_numerical_algorithms, read_namelist_numerical_upwinding_for_derivatives, &
         read_namelist_numerical_extra

      implicit none

      logical :: error = .false.

      if (initialised) return

      if (proc0) call read_namelist_time_trace_options(nstep, tend, autostop, avail_cpu_time)

      if (proc0) call read_namelist_time_step(delt, delt_option_switch, &
                                       delt_max, delt_min, &
                                       cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower)

      if (proc0) call read_namelist_numerical_algorithms(explicit_algorithm_switch, flip_flop, &
                                                  stream_implicit, stream_iterative_implicit, &
                                                  stream_matrix_inversion, driftkinetic_implicit, &
                                                  mirror_implicit, mirror_semi_lagrange, &
                                                  mirror_linear_interp, drifts_implicit, &
                                                  fully_implicit, fully_explicit, &
                                                  maxwellian_inside_zed_derivative, &
                                                  use_deltaphi_for_response_matrix, & 
                                                  split_parallel_dynamics, &
                                                  maxwellian_normalization)

      if (proc0) call read_namelist_numerical_upwinding_for_derivatives(time_upwind, zed_upwind, vpa_upwind)

      if (proc0) call read_namelist_numerical_extra(nitt, fphi, ky_solve_real, &
                                             ky_solve_radial, fields_kxkyz, mat_gen, mat_read, &
                                             lu_option_switch, rng_seed, print_extra_info_to_terminal)

      if (proc0) call check_numerical_inputs 
      call broadcast_parameters

      initialised = .true.

   contains

      !**********************************************************************
      !                         READ INPUT OPTIONS                          !
      !**********************************************************************
      ! Change the other parameters for consistently.
      !**********************************************************************
      subroutine check_numerical_inputs

         use physics_parameters, only: full_flux_surface
         use physics_parameters, only: include_apar, include_bpar
         use physics_parameters, only: include_parallel_streaming
         use physics_parameters, only: include_mirror
         use physics_parameters, only: include_nonlinear
         use physics_parameters, only: rhostar

         implicit none

         integer :: ierr

         ! Abort if neither tend nor nstep are set
         if (tend < 0 .and. nstep < 0) then
            ierr = error_unit()
            write (ierr, *) ''
            write (ierr, *) 'Please specify either <nstep> or <tend> in the <numerical_parameters> namelist.'
            write (ierr, *) 'Aborting.'
            write (*, *) ''
            write (*, *) 'Please specify either <nstep> or <tend> in the <numerical_parameters> namelist.'
            write (*, *) 'Aborting.'
            error = .true.
         end if

         ! Abort if cfl_cushion_lower>cfl_cushion_upper or if cfl_cushion_lower==cfl_cushion_upper
         if ((cfl_cushion_lower > cfl_cushion_upper - 0.001) &
            .or. (cfl_cushion_middle > cfl_cushion_upper - 0.001) &
            .or. (cfl_cushion_middle < cfl_cushion_lower + 0.001)) then
            ierr = error_unit()
            write (ierr, *) ''
            write (ierr, *) 'Please make sure that <cfl_cushion_upper> is bigger than <cfl_cushion_lower>,'
            write (ierr, *) 'and that <cfl_cushion_middle> lies in between <cfl_cushion_upper> and <cfl_cushion_lower>.'
            write (ierr, *) 'Aborting.'
            write (*, *) ''
            write (*, *) 'Please make sure that <cfl_cushion_upper> is bigger than <cfl_cushion_lower>,'
            write (*, *) 'and that <cfl_cushion_middle> lies in between <cfl_cushion_upper> and <cfl_cushion_lower>.'
            write (*, *) 'Aborting.'
            error = .true.
         end if

         ! Semi-lagrange advance of mirror term is not supported for EM simulations
         if (include_apar .and. mirror_semi_lagrange) then
            write (*, *) '!!!WARNING!!!'
            write (*, *) 'mirror_semi_lagrange = .true. is not supported for electromagnetic simulations.'
            write (*, *) 'forcing mirror_semi_lagrange = .false.'
            write (*, *) '!!!WARNING!!!'
            mirror_semi_lagrange = .false.
         end if

         if (drifts_implicit) then
            if (.not. stream_implicit) then
               write (*, *) '!!!WARNING!!!'
               write (*, *) 'drifts_implicit = T requires stream_implicit = T.'
               write (*, *) 'forcing drifts_implicit = F.'
               write (*, *) '!!!WARNING!!!'
               drifts_implicit = .false.
            else if (.not. include_parallel_streaming) then
               write (*, *) '!!!WARNING!!!'
               write (*, *) 'drifts_implicit = T requires include_parallel_streaming = T.'
               write (*, *) 'forcing drifts_implicit = F.'
               write (*, *) '!!!WARNING!!!'
               drifts_implicit = .false.
            end if
            
            if (rhostar > epsilon(0.0)) then
               write (*, *) '!!!WARNING!!!'
               write (*, *) 'drifts_implicit = T, coupled with rhostar > 0, has been observed'
               write (*, *) 'to lead to numerical instability.  unless you know what you are doing,'
               write (*, *) 'it is suggested that you set drifts_implicit = F or rhostar = 0.'
               write (*, *) '!!!WARNING!!!'
            end if
         end if

         if (.not. full_flux_surface) then  
            nitt = 1
         end if

         ! Print warning messages and override inconsistent or unsupported options for full_flux_surface = T
         if (full_flux_surface) then
            if (fields_kxkyz) then
               write (*, *)
               write (*, *) '!!!WARNING!!!'
               write (*, *) 'The option fields_kxkyz=T is not currently supported for full_flux_surface=T.'
               write (*, *) 'Forcing fields_kxkyz=F.'
               write (*, *) '!!!WARNING!!!'
               write (*, *)
               fields_kxkyz = .false.
            end if
            if (mirror_semi_lagrange) then
               write (*, *)
               write (*, *) '!!!WARNING!!!'
               write (*, *) 'The option mirror_semi_lagrange=T is not consistent with full_flux_surface=T.'
               write (*, *) 'Forcing mirror_semi_lagrange=F.'
               write (*, *) '!!!WARNING!!!'
               mirror_semi_lagrange = .false.
            end if
            if (stream_implicit) then
               driftkinetic_implicit = .true.
            end if
            if (maxwellian_normalization) then 
               write (*, *)
               write (*, *) '!!!WARNING!!!'
               write (*, *) 'The option maxwellian_normalisation=T is not consistent with full_flux_surface=T.'
               write (*, *) 'Forcing maxwellian_normalisation=F.'
               write (*, *) '!!!WARNING!!!'
               maxwellian_normalization = .false.
            end if
            
         else
            driftkinetic_implicit = .false.
         end if

         if (fully_explicit) flip_flop = .false.
         
         ! Calculate some useful derived quantities that are used repeatedly across modules
         time_upwind_plus = 0.5 * (1.0 + time_upwind)
         time_upwind_minus = 0.5 * (1.0 - time_upwind)
         zed_upwind_plus = 0.5 * (1.0 + zed_upwind)
         zed_upwind_minus = 0.5 * (1.0 - zed_upwind)
         
         if (.not. include_mirror) mirror_implicit = .false.
         if (.not. include_parallel_streaming) then
            stream_implicit = .false.
            driftkinetic_implicit = .false.
         end if
         
         if (mirror_implicit .or. stream_implicit .or. drifts_implicit) then
            fully_explicit = .false.
         else
            fully_explicit = .true.
         end if
         
         if (mirror_implicit .and. stream_implicit .and. drifts_implicit .and. .not. include_nonlinear) then
            fully_implicit = .true.
         else
            fully_implicit = .false.
         end if

       end subroutine check_numerical_inputs

      !**********************************************************************
      !                         BROADCAST OPTIONS                           !
      !**********************************************************************
      ! Broadcast these parameters to all the processors - necessary because
      ! the above was only done for the first processor (proc0).
      !**********************************************************************
      subroutine broadcast_parameters

         use mp, only: broadcast
         
         implicit none    
         ! Exit stella if we ran into an error
         call broadcast(error)
         if (error) call mp_abort('Aborting in numerical_parameters.f90')

         ! Time trace options
         call broadcast(nstep)
         call broadcast(tend)
         call broadcast(avail_cpu_time)
         call broadcast(autostop)

         ! Time step options
         call broadcast(delt)
         call broadcast(delt_option_switch)
         call broadcast(delt_max)
         call broadcast(delt_min)
         call broadcast(cfl_cushion_upper)
         call broadcast(cfl_cushion_middle)
         call broadcast(cfl_cushion_lower)

         ! Numerical algorithms
         call broadcast(explicit_algorithm_switch)
         call broadcast(flip_flop)
         call broadcast(stream_implicit)
         call broadcast(stream_iterative_implicit)
         call broadcast(stream_matrix_inversion)
         call broadcast(driftkinetic_implicit)
         call broadcast(mirror_implicit)
         call broadcast(mirror_semi_lagrange)
         call broadcast(mirror_linear_interp)
         call broadcast(drifts_implicit)
         call broadcast(fully_explicit)
         call broadcast(fully_implicit)
         call broadcast(maxwellian_inside_zed_derivative)
         call broadcast(use_deltaphi_for_response_matrix)
         call broadcast(split_parallel_dynamics)
         call broadcast(maxwellian_normalization)

         ! numerical_upwinding_for_derivatives
         call broadcast(time_upwind_plus)
         call broadcast(time_upwind_minus)
         call broadcast(zed_upwind_plus)
         call broadcast(zed_upwind_minus)
         call broadcast(zed_upwind)
         call broadcast(vpa_upwind)
         call broadcast(time_upwind)

         ! Extra - need to move
         call broadcast(nitt)
         call broadcast(fphi)
         call broadcast(ky_solve_radial)
         call broadcast(ky_solve_real)
         call broadcast(fields_kxkyz)
         call broadcast(mat_gen)
         call broadcast(mat_read)
         call broadcast(lu_option_switch)
         call broadcast(rng_seed)
         call broadcast(print_extra_info_to_terminal)

      end subroutine broadcast_parameters
    
    end subroutine read_numerical_parameters
  
  subroutine finish_read_numerical_parameters
    
    implicit none
    initialised = .false.
    
  end subroutine finish_read_numerical_parameters

end module numerical_parameters
