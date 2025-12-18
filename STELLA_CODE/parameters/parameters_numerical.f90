!###############################################################################
!########################## READ NUMERICAL PARAMETERS ##########################
!###############################################################################
! Namelists: 
! &time_trace_options
! &time_step
! &numerical_algorithms 
! &numerical_upwinding_for_derivatives
! &flux_annulus - to be changed
! These flags will allow you to toggle the algorithm choices in stella.
!###############################################################################
module parameters_numerical
   
   ! Read the parameters for <delt_option_switch> from namelist_parameters_numerical.f90
   use namelist_parameters_numerical, only: delt_option_hand
   use namelist_parameters_numerical, only: delt_option_auto
   
   ! Read the parameters for <explicit_algorithm_switch> from namelist_parameters_numerical.f90
   use namelist_parameters_numerical, only: explicit_algorithm_rk2
   use namelist_parameters_numerical, only: explicit_algorithm_rk3
   use namelist_parameters_numerical, only: explicit_algorithm_rk4
   use namelist_parameters_numerical, only: explicit_algorithm_euler

   implicit none
   
   ! Although the parameters are available through namelist_species
   ! make them available through grids_species as well
   public :: delt_option_hand, delt_option_auto
   public :: explicit_algorithm_rk3, explicit_algorithm_rk2
   public :: explicit_algorithm_rk4, explicit_algorithm_euler

   ! Time trace options
   public :: nstep, tend
   public :: autostop
   public :: avail_cpu_time

   ! Time step options
   public :: delt
   public :: delt_option_switch
   public :: delt_max, delt_min
   public :: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower

   ! Numerical algorithms
   public :: explicit_algorithm_switch
   public :: flip_flop
   public :: stream_implicit, stream_iterative_implicit, stream_matrix_inversion, driftkinetic_implicit
   public :: mirror_implicit, mirror_semi_lagrange, mirror_linear_interp
   public :: drifts_implicit
   public :: fully_implicit, fully_explicit
   public :: maxwellian_inside_zed_derivative, use_deltaphi_for_response_matrix
   public :: split_parallel_dynamics

   ! numerical_upwinding_for_derivatives
   public :: time_upwind, time_upwind_plus, time_upwind_minus
   public :: zed_upwind, zed_upwind_plus, zed_upwind_minus
   public :: vpa_upwind

   ! Flux annulus options
   public :: nitt
   
   ! Public subroutines that are read by the main stella routine.
   public :: read_parameters_numerical, finish_read_parameters_numerical

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
   real :: delt_max, delt_min
   real :: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower

   ! Numerical algorithms
   integer :: explicit_algorithm_switch
   logical :: flip_flop
   logical :: stream_implicit, stream_iterative_implicit, stream_matrix_inversion, driftkinetic_implicit
   logical :: mirror_implicit, mirror_semi_lagrange, mirror_linear_interp
   logical :: drifts_implicit, fully_implicit, fully_explicit
   logical :: maxwellian_inside_zed_derivative, use_deltaphi_for_response_matrix
   
   ! if split_parallel_dynamics = .true. (default), use operator splitting
   ! to treat parallel streaming and mirror term separately
   logical :: split_parallel_dynamics

   ! numerical_upwinding_for_derivatives
   real :: time_upwind, time_upwind_plus, time_upwind_minus
   real :: zed_upwind, zed_upwind_plus, zed_upwind_minus
   real :: vpa_upwind

   ! Extra - need to move
   integer :: nitt

   ! Internal
   logical :: initialised = .false.
   
contains

  !======================================================================
  !===================== READ NUMERICAL PARAMETERS ======================
  !======================================================================
   subroutine read_parameters_numerical(fields_kxkyz)

      use mp, only: proc0, mp_abort
      use text_options, only: text_option, get_option_value
      use file_utils, only: input_unit, input_unit_exist
      use namelist_parameters_numerical, only: read_namelist_time_trace_options
      use namelist_parameters_numerical, only: read_namelist_time_step
      use namelist_parameters_numerical, only: read_namelist_numerical_algorithms
      use namelist_parameters_numerical, only: read_namelist_numerical_upwinding_for_derivatives
      use namelist_parameters_numerical, only: read_namelist_flux_annulus

      implicit none
         
      logical, intent(in out) :: fields_kxkyz
      logical :: error = .false.

      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised) return
      initialised = .true.

      ! Read the input parameters
      if (proc0) call read_namelist_time_trace_options(nstep, tend, autostop, avail_cpu_time)
      if (proc0) call read_namelist_time_step(delt, delt_option_switch, delt_max, delt_min, &
         cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower)
      if (proc0) call read_namelist_numerical_algorithms(explicit_algorithm_switch, flip_flop, &
         stream_implicit, stream_iterative_implicit, stream_matrix_inversion, driftkinetic_implicit, &
         mirror_implicit, mirror_semi_lagrange, mirror_linear_interp, drifts_implicit, &
         fully_implicit, fully_explicit, maxwellian_inside_zed_derivative, &
         use_deltaphi_for_response_matrix, split_parallel_dynamics)
      if (proc0) call read_namelist_numerical_upwinding_for_derivatives(time_upwind, zed_upwind, vpa_upwind)
      if (proc0) call read_namelist_flux_annulus(nitt)
      
      ! Check the numerical inputs
      if (proc0) call check_numerical_inputs(fields_kxkyz)

      ! Broadcast the input parameters to all processors
      call broadcast_parameters

   contains

      !**********************************************************************
      !                         READ INPUT OPTIONS                          !
      !**********************************************************************
      ! Change the other parameters for consistently.
      !**********************************************************************
      subroutine check_numerical_inputs(fields_kxkyz)

         use parameters_physics, only: full_flux_surface
         use parameters_physics, only: include_apar
         use parameters_physics, only: include_parallel_streaming
         use parameters_physics, only: include_mirror
         use parameters_physics, only: include_nonlinear
         use parameters_physics, only: include_xdrift, include_ydrift
         use parameters_physics, only: rhostar
         use mp, only: mp_abort
         
         implicit none
         
         logical, intent(in out) :: fields_kxkyz

         !----------------------------------------------------------------------

         ! Abort if neither tend nor nstep are set
         if (tend < 0 .and. nstep < 0) then
            call mp_abort('Please specify either <nstep> or <tend> in the <parameters_numerical> namelist. Aborting')
         end if

         ! Abort if cfl_cushion_lower>cfl_cushion_upper or if cfl_cushion_lower==cfl_cushion_upper
         if (cfl_cushion_lower > cfl_cushion_upper - 0.001) then
            call mp_abort('Please make sure that <cfl_cushion_upper> is bigger than <cfl_cushion_lower>. Aborting')
         end if
         if (cfl_cushion_middle > cfl_cushion_upper - 0.001) then
            call mp_abort('Please make sure that <cfl_cushion_upper> is bigger than <cfl_cushion_middle>. Aborting')
         end if
         if (cfl_cushion_middle < cfl_cushion_lower + 0.001) then
            call mp_abort('Please make sure that <cfl_cushion_middle> is bigger than <cfl_cushion_lower>. Aborting')
         end if
         
         ! The flag <split_parallel_dynamics> does not do anything for now
         if (.not. split_parallel_dynamics) then
            call mp_abort('The option split_parallel_dynamics = .false. has not been implemented yet. Aborting.')
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
            
         else
            driftkinetic_implicit = .false.
         end if
         
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

         if ( (.not. include_mirror .or. .not. mirror_implicit ) .and. &
              (.not. include_parallel_streaming .or. .not. stream_implicit ) .and. &
              (.not. (include_xdrift .and. include_ydrift) .or. .not. drifts_implicit ) ) then 
            fully_explicit = .true.
         else
            fully_explicit = .false.
         end if

         if (fully_explicit) flip_flop = .false.
         
         if ( (.not. include_mirror .or. mirror_implicit ) .and. &
              (.not. include_parallel_streaming .or. stream_implicit ) .and. &
              (.not. (include_xdrift .or. include_ydrift) .or. drifts_implicit ) .and. &
               .not. include_nonlinear ) then
            fully_implicit = .true.
         else
            fully_implicit = .false.
         end if

         if (fully_explicit .and. fully_implicit) then
            write(*,*) 'WARNING: NO TERMS BEING INCLUDED - setting fully_explicit = .true.'
            fully_implicit = .false.
            fully_explicit = .true.
          end if
         
         ! Linear simulations can be stopped automatically when the growth rate
         ! becomes saturated. Make sure this flag is off for nonlinear simulations.
         if (include_nonlinear) autostop = .false.

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
         if (error) call mp_abort('Aborting in parameters_numerical.f90')

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
         
      end subroutine broadcast_parameters
    
    end subroutine read_parameters_numerical
  
  subroutine finish_read_parameters_numerical
    
    implicit none
    initialised = .false.
    
  end subroutine finish_read_parameters_numerical

end module parameters_numerical
