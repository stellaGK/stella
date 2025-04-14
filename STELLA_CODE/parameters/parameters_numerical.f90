!###############################################################################
!########################### READ NUMERICAL PARAMETES ##########################
!###############################################################################
! Namelist: &numerical
! These flags will allow you to toggle the algorithm choices in stella.
!###############################################################################
module parameters_numerical

   implicit none

   !> Public subroutines that are read by the main stella routine.
   public :: read_parameters_numerical, finish_read_parameters_numerical
   
   !> Algorithm schemes
   public :: stream_implicit, stream_iterative_implicit, stream_matrix_inversion, driftkinetic_implicit
   public :: mirror_implicit, mirror_semi_lagrange, mirror_linear_interp
   public :: drifts_implicit
   public :: fully_implicit, fully_explicit
   public :: maxwellian_inside_zed_derivative, use_deltaphi_for_response_matrix
   public :: split_parallel_dynamics
   
   public :: maxwellian_normalization
   
   !> Upwinding options 
   public :: time_upwind, time_upwind_plus, time_upwind_minus
   public :: zed_upwind, zed_upwind_plus, zed_upwind_minus
   public :: vpa_upwind
   
   !> Time options -> TODO-GA or TODO-HT need better description
   public :: nstep, delt, tend, nitt
   
   public :: flip_flop
   public :: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
   public :: avail_cpu_time
   public :: delt_max, delt_min
   public :: autostop
   !> TODO-HT or TODO-GA get rid of this??
   public :: fphi
   
   public :: ky_solve_radial, ky_solve_real
   public :: fields_kxkyz, mat_gen, mat_read

   public :: print_extra_info_to_terminal

   !> Explicit time-stepping option 
   public :: explicit_option_switch, explicit_option_rk3, explicit_option_rk2
   public :: explicit_option_rk4, explicit_option_euler

   !> Delt-options
   public :: delt_option_switch
   public :: delt_option_hand, delt_option_auto

   !> LU-option
   public :: lu_option_switch
   public :: lu_option_local, lu_option_none, lu_option_global
   
   !> TODO-GA - REMOVE
   public :: rng_seed
   
   private

   !> Algorithm schemes
   logical :: stream_implicit, stream_iterative_implicit, stream_matrix_inversion, driftkinetic_implicit
   logical :: mirror_implicit, mirror_semi_lagrange, mirror_linear_interp
   logical :: drifts_implicit, fully_implicit, fully_explicit
   logical :: maxwellian_inside_zed_derivative, use_deltaphi_for_response_matrix
   ! if split_parallel_dynamics = .true. (default), use operator splitting
   ! to treat parallel streaming and mirror term separately
   logical :: split_parallel_dynamics
   !> REMOVE
   logical :: maxwellian_normalization

   !> Upwinding options 
   real :: time_upwind, time_upwind_plus, time_upwind_minus
   real :: zed_upwind, zed_upwind_plus, zed_upwind_minus
   real :: vpa_upwind

   integer :: delt_option_switch, lu_option_switch
   integer, parameter :: delt_option_hand = 1, delt_option_auto = 2

   integer, parameter :: lu_option_none = 1, &
        lu_option_local = 2, &
        lu_option_global = 3

   !> For CFL conditions
   integer :: nstep 
   real :: delt, tend
   real :: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
   real :: avail_cpu_time
   real :: delt_max, delt_min
   integer :: explicit_option_switch
   logical :: flip_flop
   integer, parameter :: explicit_option_rk3 = 1, &
        explicit_option_rk2 = 2, &
        explicit_option_rk4 = 3, &
        explicit_option_euler = 4

   logical :: autostop
   real :: fphi
   logical :: fields_kxkyz, mat_gen, mat_read
   logical :: ky_solve_real
   integer :: ky_solve_radial
   integer :: nitt
   logical :: print_extra_info_to_terminal
   
   logical :: initialised = .false.

   !> REMOVE
   integer :: rng_seed

   !> Internal
   logical :: old_nml_exist
   
contains

  !======================================================================
  !===================== READ NUMERICAL PARAMETERS ======================
  !======================================================================
   subroutine read_parameters_numerical

      use mp, only: proc0, mp_abort
      use text_options, only: text_option, get_option_value
      use file_utils, only: input_unit, error_unit, input_unit_exist

      implicit none

      logical :: error = .false.
      character(10) :: explicit_option
      character(20) :: delt_option, lu_option

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
      ! &numerical'.
      !**********************************************************************
      subroutine set_default_parameters

         implicit none

         stream_implicit = .true.
         stream_iterative_implicit = .false.
         mirror_implicit = .true.
         drifts_implicit = .false.

         stream_matrix_inversion = .false.
         mirror_semi_lagrange = .true.
         mirror_linear_interp = .false. 
         maxwellian_inside_zed_derivative = .false. 
         use_deltaphi_for_response_matrix = .false.
         split_parallel_dynamics = .true.
         maxwellian_normalization = .false. 
         zed_upwind = 0.02
         vpa_upwind = 0.02
         time_upwind = 0.02
         fphi = 1.0
         !> Stella runs until t*v_{th,i}/a=tend or until istep=nstep
         nstep = - 1
         delt = 0.03 !> Set some number for the time step - QUESTION: should we set to be a silly number?
         tend = -1.0
         delt_option = 'default'
         !> The response matrix is solved with a none, local or global scheme, local seems to be the most efficient
         lu_option = 'default'
         avail_cpu_time = 1.e10
         !> Code_dt needs to stay within [cfl_dt*cfl_cushion_upper, cfl_dt*cfl_cushion_lower]
         !> code_dt can be increased if cfl_dt increases, however, never increase above delt_max (=delt by default)
         !> Exit stella if code_dt < delt_min (e.g. when the code blows up)
         cfl_cushion_upper = 0.5       !> Stay a factor of 2 under the CFL condition, otherwise it might run out of control
         cfl_cushion_middle = 0.25     !> If code_dt>cfl_dt/2 or code_dt<cfl_dt/100000, set code_dt to cfl_dt/4
         cfl_cushion_lower = 0.00001   !> Default is very low to not trigger it.
         delt_max = -1
         delt_min = 1.e-10
         autostop = .true. 
         
         fields_kxkyz = .false. 
         mat_gen = .false. 
         mat_read = .false.

         ky_solve_radial = 0
         ky_solve_real = .false.
         nitt = 1
         print_extra_info_to_terminal = .true.

         explicit_option = 'default'
         flip_flop = .false.

         !> TODO-GA REMOVE
         rng_seed = -1 !negative values use current time as seed   
      end subroutine

      !**********************************************************************
      !                         READ INPUT OPTIONS                          !
      !**********************************************************************
      ! Overwrite any default options with those specified in the input file. 
      ! Then change the other parameters consistently.
      !**********************************************************************
      subroutine read_input_file

         use parameters_physics, only: full_flux_surface
         use parameters_physics, only: include_apar, include_bpar
         use parameters_physics, only: include_parallel_streaming
         use parameters_physics, only: include_mirror
         use parameters_physics, only: nonlinear
         use parameters_physics, only: rhostar
         !> For FFS - need to delete
         
         implicit none 

         type(text_option), dimension(3), parameter :: deltopts = &
            (/text_option('default', delt_option_auto), &
            text_option('set_by_hand', delt_option_hand), &
            text_option('check_restart', delt_option_auto)/)
         
         type(text_option), dimension(4), parameter :: lu_opts = &
            (/text_option('default', lu_option_local), &
            text_option('none', lu_option_none), &
            text_option('local', lu_option_local), &
            text_option('global', lu_option_global)/)

         type(text_option), dimension(5), parameter :: explicitopts = &
              (/text_option('default', explicit_option_rk3), &
              text_option('rk3', explicit_option_rk3), &
              text_option('rk2', explicit_option_rk2), &
              text_option('rk4', explicit_option_rk4), &
              text_option('euler', explicit_option_euler)/)
         
         integer :: ierr, in_file
         logical :: nml_exist

         namelist /parameters_numerical/ stream_implicit, stream_iterative_implicit, stream_matrix_inversion, &
            driftkinetic_implicit, mirror_implicit, mirror_semi_lagrange, mirror_linear_interp, &
            drifts_implicit, fully_implicit, fully_explicit, & 
            maxwellian_inside_zed_derivative, use_deltaphi_for_response_matrix, &
            maxwellian_normalization, zed_upwind, vpa_upwind, time_upwind, &
            fphi, nstep, delt, tend, delt_option, lu_option, avail_cpu_time, & 
            cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower, delt_max, delt_min, &
            fields_kxkyz, mat_gen, mat_read, &
            ky_solve_radial, ky_solve_real, nitt, print_extra_info_to_terminal, &
            explicit_option, flip_flop, rng_seed, autostop, &
            split_parallel_dynamics
         
         !> Overwrite the default input parameters by those specified in the input file
         !> under the heading '&numerical'
         in_file = input_unit_exist("parameters_numerical", nml_exist)
         if (nml_exist) read (unit=in_file, nml=parameters_numerical)
         
         call check_backwards_compatability

         ierr = error_unit()
         call get_option_value &
            (delt_option, deltopts, delt_option_switch, ierr, &
            "delt_option in parameters_numerical")
         call get_option_value &
            (lu_option, lu_opts, lu_option_switch, ierr, &
            "lu_option in parameters_numerical")
         call get_option_value &
            (explicit_option, explicitopts, explicit_option_switch, &
            ierr, "explicit_option in parameters_numerical")

         !> Abort if neither tend nor nstep are set
         if (tend < 0 .and. nstep < 0) then
            ierr = error_unit()
            write (ierr, *) ''
            write (ierr, *) 'Please specify either <nstep> or <tend> in the <parameters_numerical> namelist.'
            write (ierr, *) 'Aborting.'
            write (*, *) ''
            write (*, *) 'Please specify either <nstep> or <tend> in the <parameters_numerical> namelist.'
            write (*, *) 'Aborting.'
            error = .true.
         end if

         !> Abort if cfl_cushion_lower>cfl_cushion_upper or if cfl_cushion_lower==cfl_cushion_upper
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

         !> Semi-lagrange advance of mirror term is not supported for EM simulations
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

         !> Print warning messages and override inconsistent or unsupported options for full_flux_surface = T
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
         
         !> Calculate some useful derived quantities that are used repeatedly across modules
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
         
         if (mirror_implicit .and. stream_implicit .and. drifts_implicit .and. .not. nonlinear) then
            fully_implicit = .true.
         else
            fully_implicit = .false.
         end if
         

       end subroutine read_input_file

       !**********************************************************************
       !                    CHECK BACKWARDS COMPATIBILITY                    !
       !**********************************************************************
       ! Make sure stella either runs or aborts old names for variables or
       ! namelists are used
       !**********************************************************************
       subroutine check_backwards_compatability
         
         use file_utils, only: input_unit, input_unit_exist
         
         implicit none
         
         integer :: in_file
         
         ! These variables belonged to <time_advance_knobs> and are now read in <parameters_physics>
         ! We define them here so we can read the namelist, but we will not use them.
         real :: xdriftknob, ydriftknob, wstarknob

         namelist /knobs/ fphi, delt, nstep, tend, &
              delt_option, lu_option, autostop, &
              avail_cpu_time, delt_max, delt_min, &
              cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower, &
              stream_implicit, mirror_implicit, &
              drifts_implicit, use_deltaphi_for_response_matrix, &
              maxwellian_normalization, &
              stream_matrix_inversion, maxwellian_inside_zed_derivative, &
              mirror_semi_lagrange, mirror_linear_interp, &
              zed_upwind, vpa_upwind, time_upwind, &
              fields_kxkyz, mat_gen, mat_read, rng_seed, &
              ky_solve_radial, ky_solve_real, nitt, print_extra_info_to_terminal
         
         namelist /time_advance_knobs/ xdriftknob, ydriftknob, wstarknob, explicit_option, flip_flop
         
         in_file = input_unit_exist("knobs", old_nml_exist)
         if (old_nml_exist) then
            read (unit=in_file, nml=knobs)
            if(print_extra_info_to_terminal) then 
               write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(*,*) 'Please replace the namelist <knobs> in the input file with'
               write(*,*) '<parameters_numerical> and include all variable names under'
               write(*,*) 'this new namelist.'
               write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            end if
            ! write(*,*) "Aborting in parameters_numerical.f90. & 
            !      The namelist <knobs> does not exist. & 
            !      Please replace this with the title <parameters_numerical>"
            ! call mp_abort("Aborting in parameters_numerical.f90. &
            !      The namelist <knobs> does not exist. & 
            !      Please replace this with the title <parameters_numerical>")
         end if 

         in_file = input_unit_exist("time_advance_knobs", old_nml_exist)
         if (old_nml_exist) then
            read(unit=in_file, nml=time_advance_knobs) 
            if (print_extra_info_to_terminal) then 
               write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(*,*) 'Please replace the namelist <time_advance_knobs> in the input file.'
               write(*,*) 'Refer to the input paramters text file as to which namelist to use.'
               write(*,*) 'Some of these parameters have been moved to <run_parameters>'
               write(*,*) 'and others have been moves to <physics_parameters>.'
               write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            end if
  
              ! write(*,*) "Aborting in run_parameters.f90.&
            !      The namelist <time_advance_knobs> does not exist.&
            !      Please replace this with the title <numerical>"
            ! call mp_abort("Aborting in run_parameters.f90.&
            !      The namelist <time_advance_knobs> does not exist.&
            !      Please replace this with the title <run_parameters>")
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
         !> Exit stella if we ran into an error
         call broadcast(error)
         if (error) call mp_abort('Aborting in parameters_numerical.f90')

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

         call broadcast(time_upwind_plus)
         call broadcast(time_upwind_minus)
         call broadcast(zed_upwind_plus)
         call broadcast(zed_upwind_minus)
         call broadcast(zed_upwind)
         call broadcast(vpa_upwind)
         call broadcast(time_upwind)

         call broadcast(nstep) 
         call broadcast(delt)
         call broadcast(tend)
         call broadcast(nitt)
         
         call broadcast(explicit_option_switch)
         call broadcast(flip_flop)

         call broadcast(cfl_cushion_upper)
         call broadcast(cfl_cushion_middle)
         call broadcast(cfl_cushion_lower)
         call broadcast(avail_cpu_time)
         
         call broadcast(delt_max)
         call broadcast(delt_min)
         call broadcast(autostop)
         call broadcast(fphi)

         call broadcast(ky_solve_radial)
         call broadcast(ky_solve_real)
         call broadcast(fields_kxkyz)
         
         call broadcast(delt_option_switch)
         call broadcast(lu_option_switch)

         call broadcast(mat_gen)
         call broadcast(mat_read)

         call broadcast(print_extra_info_to_terminal)
         
         !> GA REMOVE
         call broadcast(rng_seed)

   
      end subroutine broadcast_parameters
    
    end subroutine read_parameters_numerical
  
  subroutine finish_read_parameters_numerical
    
    implicit none
    initialised = .false.
    
  end subroutine finish_read_parameters_numerical

end module parameters_numerical
