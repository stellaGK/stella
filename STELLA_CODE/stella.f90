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
   use init_stella, only: parse_command_line

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
      use job_manage, only: time_message
      use timers, only: time_gke
      
      ! Distribution function
      use arrays_distribution_function, only: gold
      use arrays_distribution_function, only: gnew
      
      ! Advance the distribution function using the gyrokinetic equation
      use gyrokinetic_equation_explicit, only: advance_distribution_function_using_explicit_gyrokinetic_terms
      use gyrokinetic_equation_implicit, only: advance_distribution_function_using_implicit_gyrokinetic_terms
      use gk_radial_variation, only: advance_distribution_function_and_fields_radial_variation
      use gk_radial_variation, only: advance_distribution_function_and_fields_radial_variation_init
      
      ! Fields
      use arrays_fields, only: phi, apar, bpar
      use arrays_fields, only: phi_old, apar_old
      
      ! Advance the fields using the field equations (e.g., the quasineutrality condition)
      use field_equations, only: advance_fields
      use field_equations, only: fields_updated
      
      ! Physics flags
      use parameters_physics, only: include_apar
      
      ! Numerical schemes
      use parameters_numerical, only: flip_flop
      use parameters_numerical, only: fully_explicit
      use parameters_numerical, only: fully_implicit

      implicit none

      ! Arguments
      integer, intent(in) :: istep
      logical, intent(in out) :: stop_stella

      ! Local variables
      logical :: restart_time_step, time_advance_successful
      integer :: count_restarts

      !-------------------------------------------------------------------------
      
      ! Start timer
      call time_message(.false., time_gke(:,1), 'gyrokinetic equations')

      ! Advance the distribution function and fields for radial variation runs
      call advance_distribution_function_and_fields_radial_variation_init()

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
         ! treated implicitly. By default, <flip_flop> = False, and we always advance
         ! the explicit terms first (Lie splitting, first order in time). With
         ! <flip_flop> = True, the order of the explicit and implicit parts (and of the
         ! terms inside the implicit part) is reversed every other time step, so that two
         ! steps form a symmetric splitting, which is second order in time for non-stiff
         ! problems (e.g. adiabatic electrons with time_upwind = 0).
         !
         ! Keep <flip_flop> = False by default, in particular for kinetic electrons.
         ! With kinetic electrons, flip_flop = True is unstable at time steps that
         ! flip_flop = False handles (H. Thienpondt, November 2026). The time-centred
         ! implicit scheme (time_upwind = 0.02) barely damps the stiff electron components
         ! of the parallel streaming and the mirror term, and flips their sign every time
         ! step; reversing the order every step, in step with these sign flips, makes the
         ! splitting errors of consecutive steps add up instead of cancelling. In nonlinear
         ! stellarator runs with kinetic electrons a mode at the smallest ky, one to three
         ! grid points wide at a maximum of |B|, grew from t = 0 and gave a burst of the
         ! heat flux up to 370 Q_gB; in a linear run of that mode the largest stable time
         ! step dropped from 0.047 (flip_flop = False) to 0.013 a/v_th, below the CFL time
         ! step of the nonlinear runs (0.045). It needs the implicit mirror term and the
         ! explicit drive, appears with every Runge-Kutta scheme and on master too, and
         ! disappears with time_upwind >= 0.3, which makes the scheme first order again.
         ! Its second order is also lost with kinetic electrons (order reduction), so
         ! flip_flop = True has no advantage where it matters.
         !
         ! Strang splitting, E(dt/2) I(dt) E(dt/2) every time step, was tried as a stable
         ! second-order alternative (November 2026) and removed again. With the half steps
         ! of consecutive steps merged, its repeated step is exactly the step I(dt) E(dt) of
         ! flip_flop = False, plus a half step E(dt/2) at the start and before each output:
         ! its stability, growth rates and frequencies are identical to those of
         ! flip_flop = False, and only the distribution function at the output is more
         ! accurate. It is still first order, because the implicit part is itself split into
         ! the mirror term and the parallel streaming in a fixed order. To make the time
         ! advance second order without alternating any order, one would need:
         !   - a symmetric implicit part, e.g. M(dt/2) S(dt) M(dt/2) with M the mirror and
         !     S the streaming step, with the mirror arrays (semi-Lagrange shifts or the
         !     tridiagonal matrices, and the apar response with electromagnetic fields) for
         !     dt/2, rebuilt in reset_dt; for full second order also half steps of the
         !     implicit collisions, the hyper-dissipation and the flow shear;
         !   - time_upwind = 0, since the theta-scheme is only second order when it is
         !     time-centred;
         !   - Strang splitting of the explicit and implicit parts, with merged half steps.
         ! The mirror step and its redistributions take 20-35% of the time advance, so this
         ! costs roughly 15-30% per step. It was not done, because the time step of stella
         ! is set by the CFL condition and by stability, not by accuracy.
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
      
      ! Advance the distribution function and fields for radial variation runs
      call advance_distribution_function_and_fields_radial_variation()
      
      ! Stop timer
      call time_message(.false., time_gke(:,1), 'gyrokinetic equations')

   end subroutine advance_distribution_function_and_fields
   
end program stella
