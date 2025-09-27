!###############################################################################
!###################### STELLA: Delta-f gyrokinetic code #######################
!###############################################################################
! 
! This is the main stella program. The parameters, grids and modules are
! initialised within initialise_stella(). Here, the distribution function and
! fields (phi, apar, bpar) are evolved in time using the gyrokinetic equation
! and the field equations. At every <nwrite> time steps the distribution function
! and fields are examined, and diagnostics are written to an output file.
! 
!###############################################################################
program stella

   ! Load the debug flags
   use debug_flags, only: debug => stella_debug
   
   ! While evolving the distribution function and fields in time, we will
   ! run diagnostics. Moreover we need to initialise and finialise stella.
   use diagnostics, only: diagnose_distribution_function_and_fields
   use init_stella, only: initialise_stella
   use init_stella, only: finish_stella
      
   implicit none
   
   ! Keep track of the time step
   integer :: istep

   ! For a restarted simulation, istep0 > 0
   integer :: istep0
   
   ! Used for restarted simulations
   integer :: istatus
      
   !----------------------------------------------------------------------------

   ! Read options from the command line such as --help and --version
   call parse_command_line()

   ! Initialise stella
   if (debug) write (*, *) 'stella::init_stella'
   call initialise_stella(istep0, istatus)

   ! Diagnose stella for istep=0
   if (debug) write (*, *) 'stella::diagnose_distribution_function_and_fields'
   if (istep0 == 0) call diagnose_distribution_function_and_fields(istep0)

   ! Advance the guiding-center distribution function <g>, as well as the electrostatic
   ! potential <phi> and the electromagnetic fields <apar> and <bpar> in time, using
   ! the gyrokinetic equation and the field equations (e.g. quasineutrality).
   if (debug) write (*, *) 'stella::run_stella'
   call run_stella
   
   ! If we stopped running stella, we did not have to increase the time step counter
   ! This is important for restarted simulations, so they have the correct <istep> value
   istep = istep - 1
   
   ! Finish stella
   if (debug) write (*, *) 'stella::finish_stella'
   call finish_stella(istep, last_call=.true.)

contains

   !****************************************************************************
   !      EVOLVE THE DISTRIBUTION FUNCTION AND FIELDS AND RUN DIAGNOSTICS      !
   !****************************************************************************
   ! The guiding-center distribution function <g>, the electrostatic potential
   ! <phi> and the electromagnetic fields <apar> and <bpar> are evolved in time,
   ! untill <code_time> is bigger than <tend> or <nstep> is bigger than <istep>. 
   ! 
   ! Moreover, every time step, the distribution function and fields are diagnosed.
   ! Every 10 time steps, it is checked whether the program should be stopped, and
   ! every <nsave> time steps, the data necessary to restart a simulation is saved.
   !****************************************************************************
   subroutine run_stella
      
      ! Check whether stella should be stopped
      use grids_time, only: check_code_dt
      use job_manage, only: check_cpu_time
      use job_manage, only: check_stop_file
      use diagnostics_omega, only: check_saturation_omega
      use parameters_numerical, only: avail_cpu_time
      
      ! Time trace
      use grids_time, only: code_dt
      use grids_time, only: code_time
      use grids_time, only: update_time
      use parameters_numerical, only: nstep
      use parameters_numerical, only: tend
      
      ! Write output files
      use file_utils, only: error_unit
      use file_utils, only: flush_output_file
      use save_stella_for_restart, only: save_stella_data_for_restart
      use parameters_diagnostics, only: nsave
      
      implicit none

      ! Local variables
      logical :: stop_stella = .false.
      
      !-------------------------------------------------------------------------

      ! Initialise the time step counter
      istep = istep0 + 1
      
      ! Evolve the gyrokinetic equation in time until <istep>=<nstep> or <code_time>=<tend>
      do while ((code_time <= tend .and. tend > 0) .or. (istep <= nstep .and. nstep > 0))
      
         ! Keep track of the time step when debugging
         if (debug) write (*, *) 'istep = ', istep
         
         ! Every 10 time steps, check if stella should be stopped, so stella can make a clean exit
         !     - stop stella is a *.stop file has appeared
         !     - stop stella if <elapsed_time> is within 5 minutes of <avail_cpu_time>
         !     - stop stella of <code_dt> < <code_dt_min> which happens when |phi| blows up
         !     - stop stella if <autostop> = True and <omega_vs_tkykx> has saturated over <navg> time steps
         if (mod(istep, 10) == 0) then
            call check_stop_file(stop_stella)
            call check_cpu_time(avail_cpu_time, stop_stella)
            call check_code_dt(stop_stella)
            call check_saturation_omega(istep, stop_stella)
         end if
         if (stop_stella) exit
         
         ! Advance the guiding-center distribution function <g>, as well as the fields
         ! (electrostatic potential <phi> and the electromagnetic fields <apar> and <bpar>)
         ! in time, using the gyrokinetic equation and the quasineutrality equation.
         call advance_distribution_function_and_fields(istep, stop_stella)
         
         ! During the time advance, the time step is sometimes changed and the
         ! time advance is retried, if it is changed unsuccesfully more than 5 times,
         ! we give up on evolving the gyrokinetic equation since <stop_stella> = True
         if (stop_stella) exit
         
         ! After succesfully evolving the fields, we can update code_time to code_time + code_dt
         call update_time
         
         ! Every <nsave> time steps, save the data necessary to restart stella
         if (nsave > 0 .and. mod(istep, nsave) == 0) then
            call save_stella_data_for_restart(istep, code_time, code_dt, istatus)
         end if
         
         ! Calculate diagnostics, e.g., turbulent fluxes, growth rates, density fluctuations, ...
         call diagnose_distribution_function_and_fields(istep)
         
         ! Make sure the error file is written
         call flush_output_file(error_unit())
         
         ! Increase the time step counter
         istep = istep + 1
         
      end do
   
   end subroutine run_stella
   
!###############################################################################
!############ ADVANCE THE DISTRIBUTION FUNCTION AND FIELDS IN TIME #############
!###############################################################################

   !****************************************************************************
   !            ADVANCE THE DISTRIBUTION FUNCTION AND FIELDS IN TIME            
   !****************************************************************************
   ! First advance the distribution function <g> in time using the 
   ! gyrokinetic equation. Next, adance the fields (electrostatic potential
   ! <phi>, as well as the electromagnetic fields <apar> and <bpar>) in time
   ! using the quasi-neutrality condition.
   !****************************************************************************
   subroutine advance_distribution_function_and_fields(istep, stop_stella)

      ! Parallelisation
      use mp, only: proc0, broadcast
      use debug_flags, only: debug => time_advance_debug
      
      ! Distribution function
      use arrays_distribution_function, only: gold
      use arrays_distribution_function, only: gnew
      use gyrokinetic_equation_explicit, only: advance_distribution_function_using_explicit_gyrokinetic_terms
      use gyrokinetic_equation_implicit, only: advance_distribution_function_using_implicit_gyrokinetic_terms
      
      ! Fields
      use arrays_fields, only: phi, apar, bpar
      use arrays_fields, only: phi_old, apar_old
      use field_equations, only: advance_fields
      use field_equations, only: fields_updated
      
      ! Physics flags
      use parameters_physics, only: include_apar
      
      ! Numerical schemes
      use parameters_numerical, only: flip_flop
      use parameters_numerical, only: fully_explicit
      use parameters_numerical, only: fully_implicit
      use parameters_multibox, only: rk_step
      
      ! Sources and sinks used in radial variation
      ! TODO - move radial variation stuff to a subroutine
      use gk_sources, only: include_qn_source
      use gk_sources, only: update_quasineutrality_source
      use gk_sources, only: source_option_switch
      use gk_sources, only: source_option_projection
      use gk_sources, only: source_option_krook
      use gk_sources, only: update_tcorr_krook
      use gk_sources, only: project_out_zero
      use gk_radial_variation, only: mb_communicate

      implicit none

      ! Arguments
      integer, intent(in) :: istep
      logical, intent(in out) :: stop_stella

      ! Local variables
      logical :: restart_time_step, time_advance_successful
      integer :: count_restarts

      !-------------------------------------------------------------------------

      ! Unless running in multibox mode, no need to worry about
      ! mb_communicate calls as the subroutine is immediately exited
      ! if not in multibox mode.
      if (.not. rk_step) then
         if (debug) write (*, *) 'time_advance::multibox'
         call mb_communicate(gnew)
      end if

      ! Save value of phi & apar for use in diagnostics (to obtain frequency)
      phi_old = phi
      if (include_apar) apar_old = apar

      ! Flag which is set to true once we've taken a step without needing to
      ! reset dt (which can be done by the nonlinear term(s))
      time_advance_successful = .false.

      ! If cfl_cushion_lower is chosen too close to cfl_cushion_upper, then
      ! we might get stuck restarting the time step over and over, so exit stella
      count_restarts = 1

      ! Attempt the Lie or flip-flop time advance until we've done it without the timestep changing.
      do while (.not. time_advance_successful)

         ! If we've already attempted a time advance then we've updated gnew, so reset it.
         gnew = gold

         ! Ensure fields are consistent with gnew.
         ! Use the quasi-neutrality equation to advance the fields in time
         call advance_fields(gnew, phi, apar, bpar, dist='g')

         ! Keep track whether any routine wants to modify the time step
         restart_time_step = .false.

         ! Advance the distribution equation in time using the gyrokinetic equation.
         ! We need to advance the explicit and implicit parts of the gyrokinetic
         ! equation (GKE). We use operator splitting to separately evolve all terms 
         ! treated implicitly. By default, <flip_flop> = True, and the order of 
         ! operations is reversed every time step as part of alternating direction 
         ! operator splitting. This is needed to ensure 2nd order accuracy in time. 
         ! If <flip_flop> is False instead, we always advance the explicit terms first.
         if (mod(istep, 2) == 1 .or. .not. flip_flop) then
            if (.not. fully_implicit) then
               if (debug) write (*, *) 'time_advance::advance_explicit'
               call advance_distribution_function_using_explicit_gyrokinetic_terms(gnew, restart_time_step, istep)
            end if
            if (.not. restart_time_step .and. .not. fully_explicit) then
               if (debug) write (*, *) 'time_advance::advance_implicit'
               call advance_distribution_function_using_implicit_gyrokinetic_terms(istep, phi, apar, bpar, gnew)
            end if
         else
            if (.not. fully_explicit) then
               if (debug) write (*, *) 'time_advance::advance_implicit'
               call advance_distribution_function_using_implicit_gyrokinetic_terms(istep, phi, apar, bpar, gnew)
            end if
            if (.not. fully_implicit) then
               if (debug) write (*, *) 'time_advance::advance_explicit'
               call advance_distribution_function_using_explicit_gyrokinetic_terms(gnew, restart_time_step, istep)
            end if
         end if

         ! If the time step has not been restarted, the time advance was succesfull.
         ! Otherwise, discard changes to <gnew> and start the time step again. 
         ! In this case, the fields will have to be recalculated
         if (.not. restart_time_step) then
            time_advance_successful = .true.
         else
            count_restarts = count_restarts + 1
            fields_updated = .false.
         end if

         ! At some point, give up on restarting the time step
         if (count_restarts > 5) then
            stop_stella = .true.
            call broadcast(stop_stella)
            gnew = gold
            fields_updated = .false.
            if (proc0) then
               write (*, *)
               write (*, *) 'EXITING STELLA BECAUSE WE ALREADY RESTARTED THE TIME STEP 5 TIMES.'
               write (*, *) 'CHANGE CFL_CUSHION_UPPER AND CFL_CUSHION_LOWER AND RESTART THE SIMULATION.'
               write (*, *)
            end if
            exit
         end if

      end do

      ! Presumably this is to do with the radially global version of the code?
      ! perhaps it could be packaged together with the update_delay_krook code
      ! below and made into a single call where all of this happens so that
      ! users of the flux tube version of the code need not worry about it.
      if (source_option_switch == source_option_projection) then
         call project_out_zero(gold, gnew)
         fields_updated = .false.
      end if

      gold = gnew

      ! Ensure fields are updated so that omega calculation is correct.
      call advance_fields(gnew, phi, apar, bpar, dist='g')

      ! Update the delay parameters for the Krook operator
      if (source_option_switch == source_option_krook) call update_tcorr_krook(gnew)
      if (include_qn_source) call update_quasineutrality_source

   end subroutine advance_distribution_function_and_fields

   !****************************************************************************
   !****************************************************************************
   !****************************************************************************
   ! Parse some basic command line arguments. Currently just 'version' and 'help'.
   ! This should be called before anything else, but especially before initialising MPI.
   !****************************************************************************
   subroutine parse_command_line()
   
      use git_version, only: get_git_version
      integer :: arg_count, arg_n
      integer :: arg_length
      character(len=:), allocatable :: argument
      character(len=*), parameter :: endl = new_line('a')
      
      !----------------------------------------------------------------------

      arg_count = command_argument_count()

      do arg_n = 0, arg_count
         call get_command_argument(1, length=arg_length)
         if (allocated(argument)) deallocate (argument)
         allocate (character(len=arg_length)::argument)
         call get_command_argument(1, argument)

         if ((argument == "--version") .or. (argument == "-v")) then
            write (*, '("stella version ", a)') get_git_version()
            stop
         else if ((argument == "--help") .or. (argument == "-h")) then
            write (*, '(a)') "stella [--version|-v] [--help|-h] [input file]"//endl//endl// &
               "stella is a flux tube gyrokinetic code for micro-stability and turbulence "// &
               "simulations of strongly magnetised plasma"//endl// &
               "For more help, see the documentation at https://stellagk.github.io/stella/"//endl// &
               "or create an issue https://github.com/stellaGK/stella/issues/new"//endl// &
               endl// &
               "  -h, --help     Print this message"//endl// &
               "  -v, --version  Print the stella version"
            stop
         end if
      end do
   end subroutine parse_command_line
   
end program stella
