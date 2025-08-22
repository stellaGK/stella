
module time_advance

   use debug_flags, only: debug => time_advance_debug

   implicit none 
   
   public :: init_time_advance, finish_time_advance
   public :: advance_stella

   private

   logical :: time_advance_initialized = .false.
   logical :: readinit = .false.

contains

   !****************************************************************************
   !                          INITIALISE TIME ADVANCE                          !
   !****************************************************************************

   subroutine init_time_advance

      use mp, only: proc0

      use parameters_physics, only: radial_variation
      use parameters_physics, only: include_parallel_nonlinearity
      use store_arrays_useful, only: wdriftinit, wstarinit, parnlinit, &
            radialinit, driftimpinit

      use arrays_drifts, only: init_wdrift, init_wstar
      use neoclassical_terms, only: init_neoclassical_terms
      use dissipation, only: init_collisions, include_collisions
      use parallel_streaming, only: init_parallel_streaming
      use mirror_terms, only: init_mirror
      use flow_shear, only: init_flow_shear
      use sources, only: init_quasineutrality_source, init_source_timeaverage

      use radial_variation_time_advance, only: init_radial_variation
      use parallel_nonlinearity, only: init_parallel_nonlinearity

      implicit none

      if (time_advance_initialized) return
      time_advance_initialized = .true.

      debug = debug .and. proc0

      ! Read time_advance_knobs namelist from the input file;
      ! sets the explicit time advance option, as well as allows for scaling of
      ! the x and y components of the magnetic drifts and of the drive term
      ! allocate distribution function sized arrays needed, e.g., for Runge-Kutta time advance
      if (debug) write (6, *) 'time_advance::init_time_advance::allocate_arrays'
      call allocate_arrays
      ! Set up neoclassical corrections to the equilibrium Maxwellian;
      ! only calculated/needed when simulating higher order terms in rhostar for intrinsic rotation
      if (debug) write (6, *) 'time_advance::init_time_advance::init_neoclassical_terms'
      call init_neoclassical_terms
      ! Calculate the term multiplying dg/dvpa in the mirror term
      ! and set up either the semi-Lagrange machinery or the tridiagonal matrix to be inverted
      ! if solving implicitly
      if (debug) write (6, *) 'time_advance::init_time_advance::init_mirror'
      call init_mirror
      ! Calculate the term multiplying dg/dz in the parallel streaming term
      ! and set up the tridiagonal matrix to be inverted if solving implicitly
      if (debug) write (6, *) 'time_advance::init_time_advance::init_parstream'
      call init_parallel_streaming
      ! Allocate and calculate the factors multiplying dg/dx, dg/dy, dphi/dx and dphi/dy
      ! in the magnetic drift terms
      if (debug) write (6, *) 'time_advance::init_time_advance::init_wdrift'
      call init_wdrift 
      ! Allocate and calculate the factor multiplying dphi/dy in the gradient drive term
      if (debug) write (6, *) 'time_advance::init_time_advance::init_wstar'
      call init_wstar 
      if (debug) write (6, *) 'time_advance::init_time_advance::init_flow_shear'
      call init_flow_shear
      if (debug) write (6, *) 'time_advance::init_time_advance::init_parallel_nonlinearity'
      if (include_parallel_nonlinearity) call init_parallel_nonlinearity 
      if (debug) write (6, *) 'time_advance::init_time_advance::init_radial_variation'
      if (radial_variation) call init_radial_variation
      if (include_collisions) then
         if (debug) write (6, *) 'time_advance::init_time_advance::init_collisions'
         call init_collisions
      end if
      if (debug) write (6, *) 'time_advance::init_time_advance::init_cfl'
      call init_cfl

      if (debug) write (6, *) 'time_advance::init_time_advance::init_source_timeaverage'
      call init_source_timeaverage
      if (debug) write (6, *) 'time_advance::init_time_advance::init_quasineutrality_source'
      call init_quasineutrality_source

   end subroutine init_time_advance

   subroutine allocate_arrays

      use stella_layouts, only: vmu_lo
      use z_grid, only: nzgrid, ntubes
      use parameters_kxky_grid, only: naky, nakx
      use store_arrays_distribution_fn, only: g0, g1, g2, g3
      use parameters_numerical, only: explicit_algorithm_switch, explicit_algorithm_rk3, &
           explicit_algorithm_rk2, explicit_algorithm_rk4, explicit_algorithm_euler

      implicit none

      if (.not. allocated(g0)) &
         allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g0 = 0.
      if (.not. allocated(g1)) &
         allocate (g1(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g1 = 0.
      if (.not. allocated(g2)) &
         allocate (g2(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g2 = 0.
      if (.not. allocated(g3) .and. explicit_algorithm_switch == explicit_algorithm_rk4) then
         allocate (g3(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         g3 = 0.
      else
         allocate (g3(1, 1, 1, 1, 1))
      end if

   end subroutine allocate_arrays

   subroutine init_cfl

      use mp, only: proc0, nproc, max_allreduce, min_allreduce
      use mp, only: scope, allprocs, subprocs
      use store_arrays_distribution_fn, only: wdriftx_g, wdrifty_g
      use stella_time, only: code_dt, write_dt, cfl_dt_linear
      use parameters_numerical, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
      use parameters_numerical, only: stream_implicit, mirror_implicit, drifts_implicit
      use parameters_physics, only: radial_variation, prp_shear_enabled
      use parameters_kxky_grid, only: nx
      use z_grid, only: delzed
      use velocity_grids, only: dvpa
      use grids_kxky, only: akx, aky, rho
      use parallel_streaming, only: stream
      use parallel_streaming, only: stream_rad_var1, stream_rad_var2
      use mirror_terms, only: mirror
      use flow_shear, only: prl_shear, shift_times
      use file_utils, only: runtype_option_switch, runtype_multibox
      use dissipation, only: include_collisions, collisions_implicit
      use dissipation, only: cfl_dt_vpadiff, cfl_dt_mudiff
      use debug_flags, only: print_extra_info_to_terminal

       use reset_timestep, only: reset_dt

      implicit none

      real :: cfl_dt_mirror, cfl_dt_stream, cfl_dt_shear
      real :: cfl_dt_wdriftx, cfl_dt_wdrifty
      real :: zero
      real :: wdriftx_max, wdrifty_max

      ! Avoid divide by zero in cfl_dt terms below
      zero = 100.*epsilon(0.)

      ! FLAG -- assuming equal spacing in zed!

      if (cfl_dt_linear < 0) cfl_dt_linear = code_dt / cfl_cushion_upper

      if (.not. drifts_implicit) then
         ! Get the local max value of wdriftx on each processor
         wdriftx_max = maxval(abs(wdriftx_g))
         ! Xompare these max values across processors to get global max
         if (nproc > 1) then
            call max_allreduce(wdriftx_max)
         end if
         ! NB: wdriftx_g has code_dt built-in, which accounts for code_dt factor here
         cfl_dt_wdriftx = abs(code_dt) / max(maxval(abs(akx)) * wdriftx_max, zero)
         cfl_dt_linear = cfl_dt_wdriftx
      end if

      cfl_dt_shear = abs(code_dt) / max(maxval(abs(aky)) * maxval(abs(prl_shear)), zero)
      cfl_dt_linear = min(cfl_dt_linear, cfl_dt_shear)

      if (prp_shear_enabled) then
         cfl_dt_shear = minval(shift_times)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_shear)
      end if

      if (.not. stream_implicit) then
         ! NB: stream has code_dt built-in, which accounts for code_dt factor here
         cfl_dt_stream = abs(code_dt) * delzed(0) / max(maxval(abs(stream)), zero)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_stream)
      end if

      ! TODO:GA- add correct CFL condition 
      ! if (driftkinetic_implicit) then
      !    cfl_dt_stream = abs(code_dt) * delzed(0) / max(maxval(abs(stream_correction)), zero)
      !    cfl_dt_linear = min(cfl_dt_linear, cfl_dt_stream)
      ! end if

      if (.not. mirror_implicit) then
         ! NB: mirror has code_dt built-in, which accounts for code_dt factor here
         cfl_dt_mirror = abs(code_dt) * dvpa / max(maxval(abs(mirror)), zero)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_mirror)
      end if

      if (radial_variation) then
         ! While other quantities should go here, parallel streaming with electrons
         ! is what will limit us
         cfl_dt_stream = abs(code_dt) * delzed(0) / max(maxval(abs(stream_rad_var1)), zero)
         cfl_dt_stream = cfl_dt_stream / abs(rho(nx) + zero)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_stream)

         cfl_dt_stream = abs(code_dt) * delzed(0) / max(maxval(abs(stream_rad_var2)), zero)
         cfl_dt_stream = cfl_dt_stream / abs(rho(nx) + zero)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_stream)

      end if

      if (include_collisions .and. .not. collisions_implicit) then
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_vpadiff)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_mudiff)
      end if

      if (.not. drifts_implicit) then
         ! Get the local max value of wdrifty on each processor
         wdrifty_max = maxval(abs(wdrifty_g))
         ! Compare these max values across processors to get global max
         if (nproc > 1) then
            call max_allreduce(wdrifty_max)
         end if
         ! NB: wdrifty_g has code_dt built-in, which accounts for code_dt factor here
         cfl_dt_wdrifty = abs(code_dt) / max(maxval(abs(aky)) * wdrifty_max, zero)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_wdrifty)
      end if

      if (runtype_option_switch == runtype_multibox) call scope(allprocs)
      call min_allreduce(cfl_dt_linear)
      if (runtype_option_switch == runtype_multibox) call scope(subprocs)

      if (proc0 .and. print_extra_info_to_terminal) then
         write (*, '(A)') "############################################################"
         write (*, '(A)') "                        CFL CONDITION"
         write (*, '(A)') "############################################################"
         write (*, '(A16)') 'LINEAR CFL_DT: '
         if (.not. drifts_implicit) write (*, '(A12,ES12.4)') '   wdriftx: ', cfl_dt_wdriftx
         if (.not. drifts_implicit) write (*, '(A12,ES12.4)') '   wdrifty: ', cfl_dt_wdrifty
         if (.not. stream_implicit) write (*, '(A12,ES12.4)') '   stream: ', cfl_dt_stream
         if (.not. mirror_implicit) write (*, '(A12,ES12.4)') '   mirror: ', cfl_dt_mirror
         write (*, '(A12,ES12.4)') '   total: ', cfl_dt_linear
         write (*, *)
      end if

      if (abs(code_dt) > cfl_dt_linear * cfl_cushion_upper) then
         if (proc0) then
            write (*, *) 'CHANGING TIME STEP:'
            write (*, '(A22, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
            write (*, '(A22, ES10.2E2)') "   cfl_dt_linear:"//REPEAT(' ', 50), cfl_dt_linear
            write (*, '(A22, ES10.2E2)') "   cfl_cushion_upper:"//REPEAT(' ', 50), cfl_cushion_upper
            write (*, '(A22, ES10.2E2)') "   cfl_cushion_middle:"//REPEAT(' ', 50), cfl_cushion_middle
            write (*, '(A22, ES10.2E2)') "   cfl_cushion_lower:"//REPEAT(' ', 50), cfl_cushion_lower
            write (*, '(A70)') '     ==> User-specified delt is larger than cfl_dt*cfl_cushion_upper.'//REPEAT(' ', 50)
            write (*, '(A55,ES12.4)') '     ==> Changing code_dt to cfl_dt*cfl_cushion_upper ='//REPEAT(' ', 50), cfl_dt_linear * cfl_cushion_upper
            write (*, *)
         end if
         code_dt = sign(1.0, code_dt) * cfl_dt_linear * cfl_cushion_upper
         call reset_dt
      else if (proc0 .and. print_extra_info_to_terminal) then
         call write_dt
         write (*, *)
      end if

   end subroutine init_cfl

   !****************************************************************************
   !****************************************************************************
   !                        MAIN TIME ADVANCE OF STELLA                        !
   !****************************************************************************
   !****************************************************************************

   subroutine advance_stella(istep, stop_stella)

      use store_arrays_distribution_fn, only: gold, gnew
      use store_arrays_fields, only: phi, apar, bpar
      use store_arrays_fields, only: phi_old, apar_old
      use fields, only: advance_fields, fields_updated
      use parameters_numerical, only: fully_explicit, fully_implicit
      use parameters_multibox, only: rk_step
      use sources, only: include_qn_source, update_quasineutrality_source
      use sources, only: source_option_switch, source_option_projection
      use sources, only: source_option_krook
      use sources, only: update_tcorr_krook, project_out_zero
      use parameters_physics, only: include_apar
      use mp, only: proc0, broadcast

      use parameters_numerical, only: flip_flop
      use radial_variation_time_advance, only: mb_communicate

      implicit none

      integer, intent(in) :: istep
      logical, intent(in out) :: stop_stella

      logical :: restart_time_step, time_advance_successful
      integer :: count_restarts

      ! Unless running in multibox mode, no need to worry about
      ! mb_communicate calls as the subroutine is immediately exited
      ! if not in multibox mode.
      if (.not. rk_step) then
         if (debug) write (*, *) 'time_advance::multibox'
         call mb_communicate(gnew)
      end if

      ! Save value of phi & apar
      ! for use in diagnostics (to obtain frequency)
      phi_old = phi
      if (include_apar) apar_old = apar

      ! Flag which is set to true once we've taken a step without needing to
      ! reset dt (which can be done by the nonlinear term(s))
      time_advance_successful = .false.

      ! If cfl_cushion_lower is chosen too close to cfl_cushion_upper, then
      ! we might get stuck restarting the time step over and over, so exit stella
      count_restarts = 1

      ! Attempt the Lie or flip-flop time advance until we've done it without the
      ! timestep changing.
      do while (.not. time_advance_successful)

         ! If we've already attempted a time advance then we've updated gnew, so reset it.
         gnew = gold

         ! Ensure fields are consistent with gnew.
         call advance_fields(gnew, phi, apar, bpar, dist='g')

         ! Keep track whether any routine wants to modify the time step
         restart_time_step = .false.

         ! Reverse the order of operations every time step
         ! as part of alternating direction operator splitting
         ! this is needed to ensure 2nd order accuracy in time
         if (mod(istep, 2) == 1 .or. .not. flip_flop) then
            ! Advance the explicit parts of the GKE
            if (debug) write (*, *) 'time_advance::advance_explicit'
            if (.not. fully_implicit) call advance_explicit(gnew, restart_time_step, istep)
            if (debug) write (*, *) 'time_advance::advance_implicit'
            ! Use operator splitting to separately evolve all terms treated implicitly
            if (.not. restart_time_step .and. .not. fully_explicit) call advance_implicit(istep, phi, apar, bpar, gnew)
         else
            ! Advance the explicit parts of the GKE
            if (debug) write (*, *) 'time_advance::advance_implicit'
            if (.not. fully_explicit) call advance_implicit(istep, phi, apar, bpar, gnew)
            if (debug) write (*, *) 'time_advance::advance_explicit'
            if (.not. fully_implicit) call advance_explicit(gnew, restart_time_step, istep)
         end if

         ! If the time step has not been restarted, the time advance was succesfull
         ! Otherwise, discard changes to gnew and start the time step again, fields
         ! will have to be recalculated
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
      ! perhaps it could be packaged together with thee update_delay_krook code
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

   end subroutine advance_stella

   !-------------------------------------------------------------------------------
   !                       EXPLICIT TIME ADVANCE SUBROUTINES
   !-------------------------------------------------------------------------------
   !> advance_explicit takes as input the guiding centre distribution function
   !> in k-space and updates it to account for all of the terms in the GKE that
   !> are advanced explicitly in time
   subroutine advance_explicit(g, restart_time_step, istep)

      use mp, only: proc0
      use job_manage, only: time_message
      use stella_layouts, only: vmu_lo, iv_idx
      use store_arrays_fields, only: phi, apar, bpar
      use store_arrays_useful, only: time_gke
      use z_grid, only: nzgrid
      use extended_zgrid, only: periodic, phase_shift
      use parameters_kxky_grid, only: naky
      use parameters_physics, only: include_apar
      use parallel_streaming, only: stream_sign
      use parameters_numerical, only: explicit_algorithm_switch, explicit_algorithm_rk3, &
           explicit_algorithm_rk2, explicit_algorithm_rk4, explicit_algorithm_euler
      use fields, only: advance_fields
      use g_tofrom_h, only: gbar_to_g

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      integer :: ivmu, iv, sgn, iky

      ! Start the timer for the explicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 8), ' explicit')

      ! Incoming pdf is g = h - (Z F0/T) (J0 phi + 4 mu (T/Z) (J1/gamma) bpar)
      ! if include_apar = T, convert from g to gbar = g + Z F0/T (2J0 vpa vth apar),
      ! as gbar appears in time derivative
      if (include_apar) then
         ! If the fields are not already updated, then update them
         call advance_fields(g, phi, apar, bpar, dist='g')
         call gbar_to_g(g, apar, -1.0)
      end if

      select case (explicit_algorithm_switch)
      case (explicit_algorithm_euler)
         ! Forward Euler
         call advance_explicit_euler(g, restart_time_step, istep)
      case (explicit_algorithm_rk2)
         ! SSP RK2
         call advance_explicit_rk2(g, restart_time_step, istep)
      case (explicit_algorithm_rk3)
         ! SSP RK3 -> This is the default option 
         call advance_explicit_rk3(g, restart_time_step, istep)
      case (explicit_algorithm_rk4)
         ! RK4
         call advance_explicit_rk4(g, restart_time_step, istep)
      end select

      if (include_apar) then
         ! If the fields are not already updated, then update them
         call advance_fields(g, phi, apar, bpar, dist='gbar')
         ! Implicit solve will use g rather than gbar for advance,
         ! so convert from gbar to g
         call gbar_to_g(g, apar, 1.0)
      end if

      !> Enforce periodicity for periodic (including zonal) modes
      do iky = 1, naky
         if (periodic(iky)) then
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               !> stream_sign > 0 corresponds to dz/dt < 0
               sgn = stream_sign(iv)
               g(iky, :, sgn * nzgrid, :, ivmu) = &
                  g(iky, :, -sgn * nzgrid, :, ivmu) * phase_shift(iky)**(-sgn)
            end do
         end if
      end do

      ! Stop the timer for the explicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 8), ' explicit')

   end subroutine advance_explicit

   !-------------------------------------------------------------------------------
   !                        EXPLICIT EULER TIME ADVANCE SUBROUTINE
   !-------------------------------------------------------------------------------
   !> advance_explicit_euler uses forward Euler to advance one time step
   subroutine advance_explicit_euler(g, restart_time_step, istep)

      use store_arrays_distribution_fn, only: g0
      use z_grid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step
      use radial_variation_time_advance, only: mb_communicate

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      ! rk_step only true if running in multibox mode
      if (rk_step) call mb_communicate(g)

      g0 = g

      call solve_gke(g0, g, restart_time_step, istep)

      g = g0 + g

   end subroutine advance_explicit_euler

   !-------------------------------------------------------------------------------
   !                         EXPLICIT RK2 TIME ADVANCE SUBROUTINE
   !-------------------------------------------------------------------------------
   !> advance_expliciit_rk2 uses strong stability-preserving rk2 to advance one time step
   subroutine advance_explicit_rk2(g, restart_time_step, istep)

      use store_arrays_distribution_fn, only: g0, g1
      use z_grid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step
      use radial_variation_time_advance, only: mb_communicate

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      integer :: icnt

      ! rk_step only true if running in multibox mode
      if (rk_step) call mb_communicate(g)

      g0 = g
      icnt = 1

      ! SSP rk2 algorithm to advance explicit part of code
      ! if GK equation written as dg/dt = rhs - vpar . grad h,
      ! solve_gke returns rhs*dt
      do while (icnt <= 2)
         select case (icnt)
         case (1)
            call solve_gke(g0, g1, restart_time_step, istep)
         case (2)
            g1 = g0 + g1
            if (rk_step) call mb_communicate(g1)
            call solve_gke(g1, g, restart_time_step, istep)
         end select
         if (restart_time_step) then
            ! If the code_dt is reset, we need to quit this loop and restart the timestep again
            icnt = 10
         else
            icnt = icnt + 1
         end if
      end do

      ! This is g at intermediate time level
      g = 0.5 * g0 + 0.5 * (g1 + g)

   end subroutine advance_explicit_rk2

   !-------------------------------------------------------------------------------
   !                         EXPLICIT RK3 TIME ADVANCE SUBROUTINE
   !-------------------------------------------------------------------------------
   !> advance_expliciit_rk3 uses strong stability-preserving rk3 to advance one time step
   subroutine advance_explicit_rk3(g, restart_time_step, istep)

      use store_arrays_distribution_fn, only: g0, g1, g2
      use z_grid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step
      use radial_variation_time_advance, only: mb_communicate

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      integer :: icnt

      ! rk_STEP = false unless in multibox mode
      if (rk_step) call mb_communicate(g)

      g0 = g
      icnt = 1

      ! SSP rk3 algorithm to advance explicit part of code
      ! if GK equation written as dg/dt = rhs - vpar . grad h,
      ! solve_gke returns rhs*dt
      do while (icnt <= 3)
         select case (icnt)
         case (1)
            call solve_gke(g0, g1, restart_time_step, istep)
         case (2)
            g1 = g0 + g1
            if (rk_step) call mb_communicate(g1)
            call solve_gke(g1, g2, restart_time_step, istep)
         case (3)
            g2 = g1 + g2
            if (rk_step) call mb_communicate(g2)
            call solve_gke(g2, g, restart_time_step, istep)
         end select
         if (restart_time_step) then
            ! If the code_dt is reset, we need to quit this loop and restart the timestep again
            icnt = 10
         else
            icnt = icnt + 1
         end if
      end do

      ! This is g at intermediate time level
      g = g0 / 3.+0.5 * g1 + (g2 + g) / 6.

   end subroutine advance_explicit_rk3

   !-------------------------------------------------------------------------------
   !                         EXPLICIT RK4 TIME ADVANCE SUBROUTINE
   !-------------------------------------------------------------------------------
   !> advance_expliciit_rk4 uses rk4 to advance one time step
   subroutine advance_explicit_rk4(g, restart_time_step, istep)

      use store_arrays_distribution_fn, only: g0, g1, g2, g3
      use z_grid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step
      use radial_variation_time_advance, only: mb_communicate

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      integer :: icnt

      ! rk_step is false unless in multibox mode
      if (rk_step) call mb_communicate(g)

      g0 = g
      icnt = 1

      ! RK4 algorithm to advance explicit part of code
      ! if GK equation written as dg/dt = rhs - vpar . grad h,
      ! solve_gke returns rhs*dt
      do while (icnt <= 4)
         select case (icnt)
         case (1)
            call solve_gke(g0, g1, restart_time_step, istep)
         case (2)
            ! g1 is h*k1
            g3 = g0 + 0.5 * g1
            if (rk_step) call mb_communicate(g3)
            call solve_gke(g3, g2, restart_time_step, istep)
            g1 = g1 + 2.*g2
         case (3)
            ! g2 is h*k2
            g2 = g0 + 0.5 * g2
            if (rk_step) call mb_communicate(g2)
            call solve_gke(g2, g3, restart_time_step, istep)
            g1 = g1 + 2.*g3
         case (4)
            ! g3 is h*k3
            g3 = g0 + g3
            if (rk_step) call mb_communicate(g3)
            call solve_gke(g3, g, restart_time_step, istep)
            g1 = g1 + g
         end select
         if (restart_time_step) then
            ! If the code_dt is reset, we need to quit this loop and restart the timestep again
            icnt = 10
         else
            icnt = icnt + 1
         end if
      end do

      ! This is g at intermediate time level
      g = g0 + g1 / 6.

   end subroutine advance_explicit_rk4

   !-------------------------------------------------------------------------------
   !                  NEEDED FOR ALL EXPLICIT TIME ADVANCE SUBROUTINES
   !-------------------------------------------------------------------------------
   !> solve_gke accepts as argument pdf, the guiding centre distribution function in k-space,
   !> and returns rhs_ky, the right-hand side of the gyrokinetic equation in k-space;
   !> i.e., if dg/dt = r, then rhs_ky = r*dt;
   !> note that if include_apar = T, then the input pdf is actually gbar = g + (Ze/T)*(vpa/c)*<Apar>*F0
   subroutine solve_gke(pdf, rhs_ky, restart_time_step, istep)

      use job_manage, only: time_message
      use multibox, only: add_multibox_krook

      use stella_layouts, only: vmu_lo
      use stella_transforms, only: transform_y2ky

      use store_arrays_fields, only: phi, apar, bpar
      use store_arrays_distribution_fn, only: g_scratch
      use arrays_gyro_averages, only: j0_ffs

      use calculations_kxky, only: swap_kxky_back
      use gyro_averages, only: gyro_average
      use g_tofrom_h, only: gbar_to_g 

      use parameters_physics, only: include_parallel_nonlinearity
      use parameters_physics, only: include_parallel_streaming
      use parameters_physics, only: include_mirror, include_apar
      use parameters_physics, only: include_nonlinear, include_bpar
      use parameters_physics, only: full_flux_surface, radial_variation
      use parameters_physics, only: xdriftknob, ydriftknob
      use parameters_physics, only: g_exb
      use parameters_kxky_grid, only: ikx_max, ny, naky_all
      use parameters_numerical, only: stream_implicit, mirror_implicit, drifts_implicit
      use parameters_multibox, only: include_multibox_krook

      use z_grid, only: nzgrid, ntubes
      use grids_kxky, only: zonal_mode, akx
      
      use dissipation, only: include_collisions, advance_collisions_explicit, collisions_implicit
      use sources, only: source_option_switch, source_option_krook
      use sources, only: add_krook_operator
      use parallel_streaming, only: advance_parallel_streaming_explicit
      use mirror_terms, only: advance_mirror_explicit
      use flow_shear, only: advance_parallel_flow_shear
      use advance_explicit_drifts, only: advance_wstar_explicit
      use advance_explicit_drifts, only: advance_wdriftx_explicit, advance_wdrifty_explicit
      use parallel_nonlinearity, only: advance_parallel_nonlinearity
      use radial_variation_time_advance, only: advance_radial_variation
      use advance_nonlinearity, only: advance_ExB_nonlinearity
      use dissipation, only: hyper_dissipation

      use fields, only: fields_updated, advance_fields
      use fields_radial_variation, only: get_radial_correction

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: pdf
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out), target :: rhs_ky
      logical, intent(out) :: restart_time_step
      integer, intent(in) :: istep

      complex, dimension(:, :, :, :, :), allocatable, target :: rhs_y
      complex, dimension(:, :, :, :, :), pointer :: rhs
      complex, dimension(:, :), allocatable :: rhs_ky_swap

      integer :: iz, it, ivmu

      rhs_ky = 0.

      ! If full_flux_surface = .true., then initially obtain the RHS of the GKE in alpha-space;
      ! will later inverse Fourier transform to get RHS in k_alpha-space
      if (full_flux_surface) then
         ! rhs_ky will always be needed as the array returned by the subroutine,
         ! but intermediate array rhs_y (RHS of gke in alpha-space) only needed for full_flux_surface = .true.
         allocate (rhs_y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         rhs_y = 0.
         ! rhs is array referred to for both flux tube and full-flux-surface simulations;
         ! for full-flux-surface it should point to rhs_y
         rhs => rhs_y
      else
         ! rhs is array referred to for both flux tube and full-flux-surface simulations;
         ! for flux tube it should point to rhs_ky
         rhs => rhs_ky
      end if

      ! Start with g in k-space and (ky,kx,z) local
      ! obtain fields corresponding to g
      if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_fields'

      ! If advancing apar, then gbar is evolved in time rather than g
      if (include_apar) then
         call advance_fields(pdf, phi, apar, bpar, dist='gbar')

         ! Convert from gbar to g = h - (Z F0/T)( J0 phi + 4 mu (T/Z) (J1/gamma) bpar),
         ! as all terms on RHS of GKE use g rather than gbar
         call gbar_to_g(pdf, apar, 1.0)
      else
         call advance_fields(pdf, phi, apar, bpar, dist='g')
      end if

      if (radial_variation) call get_radial_correction(pdf, phi, dist='gbar')

      ! Obtain the gyro-average of the electrostatic potential phi and store in g_scratch;
      ! this can be a particularly costly operation when simulating a full flux surface
      ! due to the coupling of different k-alphas inherent in the gyro-average;
      ! calculate once here to avoid repeated calculation later
      ! TODO-GA : can this be spec up??
      if (full_flux_surface) call gyro_average(phi, g_scratch, j0_ffs)

      !! INSERT TEST HERE TO SEE IF dg/dy, dg/dx, d<phi>/dy, d<phi>/dx WILL BE NEEDED
      !! IF SO, PRE-COMPUTE ONCE HERE

      ! Default is to continue with same time step size.
      ! if estimated CFL condition for nonlinear terms is violated
      ! then restart_time_step will be set to .true.
      restart_time_step = .false.
      ! Calculate and add ExB nonlinearity to RHS of GK eqn
      ! do this first, as the CFL condition may require a change in time step
      ! and thus recomputation of mirror, wdrift, wstar, and parstream
      if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_ExB_nonlinearity'
      if (include_nonlinear) call advance_ExB_nonlinearity(pdf, rhs, restart_time_step, istep)

      ! Include contribution from the parallel nonlinearity (aka turbulent acceleration)
      if (include_parallel_nonlinearity .and. .not. restart_time_step) &
         call advance_parallel_nonlinearity(pdf, rhs, restart_time_step)

      if (.not. restart_time_step) then

         ! Include contribution from perp flow shear in the parallel component of the toroidal flow
         if ((g_exb**2) > epsilon(0.0)) call advance_parallel_flow_shear(rhs)

         ! Calculate and add mirror term to RHS of GK eqn
         if (include_mirror .and. .not. mirror_implicit) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_mirror_explicit'
            call advance_mirror_explicit(pdf, rhs)
         end if

         if (.not. drifts_implicit) then
            ! Calculate and add alpha-component of magnetic drift term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wdrifty_explicit'
            if (abs(ydriftknob) > epsilon(0.0)) then
               call advance_wdrifty_explicit(pdf, phi, bpar, rhs)
            end if

            ! Calculate and add psi-component of magnetic drift term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wdriftx_explicit'
            if (abs(xdriftknob) > epsilon(0.0)) then
               call advance_wdriftx_explicit(pdf, phi, bpar, rhs)
            end if
            
            ! Calculate and add omega_* term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wstar_explicit'
            call advance_wstar_explicit(phi, rhs)
         end if

         ! Calculate and add contribution from collisions to RHS of GK eqn
         if (include_collisions .and. .not. collisions_implicit) call advance_collisions_explicit(pdf, phi, bpar, rhs)

         ! Calculate and add parallel streaming term to RHS of GK eqn
         if (include_parallel_streaming .and. (.not. stream_implicit)) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_parallel_streaming_explicit'
            call advance_parallel_streaming_explicit(pdf, phi, bpar, rhs)
         end if
         
         if (hyper_dissipation) then
            call advance_hyper_explicit(pdf, rhs)
         end if

         ! If simulating a full flux surface (flux annulus), all terms to this point have been calculated
         ! in real-space in alpha (y); transform to kalpha (ky) space before adding to RHS of GKE.
         ! NB: it may be that for fully explicit calculation, this transform can be eliminated with additional code changes
         if (full_flux_surface) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::transform_y2ky'
            allocate (rhs_ky_swap(naky_all, ikx_max))
            it = 1
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               do iz = -nzgrid, nzgrid
                  call transform_y2ky(rhs_y(:, :, iz, it, ivmu), rhs_ky_swap)
                  call swap_kxky_back(rhs_ky_swap, rhs_ky(:, :, iz, it, ivmu))
               end do
               ! ensure that the kx=ky=0 mode is zeroed out
               if (zonal_mode(1) .and. akx(1) < epsilon(0.)) then
                  rhs_ky(1, 1, :, it, ivmu) = 0.0
               end if
            end do
            deallocate (rhs_ky_swap)
         end if

         if (radial_variation) call advance_radial_variation(pdf, rhs)

         if (source_option_switch == source_option_krook) call add_krook_operator(pdf, rhs)

         if (include_multibox_krook) call add_multibox_krook(pdf, rhs)

      end if

      ! if advancing apar, need to convert input pdf back from g to gbar
      if (include_apar) call gbar_to_g(pdf, apar, -1.0)

      fields_updated = .false.

      if (allocated(rhs_y)) deallocate (rhs_y)
      nullify (rhs)

   end subroutine solve_gke

   subroutine advance_hyper_explicit(gin, gout)

      use stella_layouts, only: vmu_lo
      use z_grid, only: nzgrid, ntubes
      use parameters_kxky_grid, only: naky, nakx
      use hyper, only: advance_hyper_vpa, advance_hyper_zed
      use hyper, only: hyp_zed, hyp_vpa

      implicit none
      
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      complex, dimension(:, :, :, :, :), allocatable :: dg

      allocate (dg(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); dg = 0.0

      if (hyp_zed) then
         call advance_hyper_zed(gin, dg)
         gout = gout + dg
      end if
      if (hyp_vpa) then
         call advance_hyper_vpa(gin, dg)
      end if
      deallocate (dg)
      
   end subroutine advance_hyper_explicit

   !-------------------------------------------------------------------------------
   !                           IMPLICIT TIME ADVANCE SUBROUTINE
   !-------------------------------------------------------------------------------

   subroutine advance_implicit(istep, phi, apar, bpar, g)

      use mp, only: proc0
      use job_manage, only: time_message

      use stella_layouts, only: vmu_lo

      use z_grid, only: nzgrid
      use store_arrays_useful, only: time_gke

      use parameters_physics, only: include_parallel_streaming
      use parameters_physics, only: radial_variation, full_flux_surface
      use parameters_physics, only: include_mirror, prp_shear_enabled
      use parameters_numerical, only: stream_implicit, mirror_implicit, drifts_implicit
      use parameters_multibox, only: rk_step
      use parameters_numerical, only: flip_flop

      use fields, only: advance_fields, fields_updated
      use dissipation, only: hyper_dissipation
      use dissipation, only: collisions_implicit, include_collisions
      use dissipation, only: advance_collisions_implicit
      use hyper, only: advance_hyper_dissipation
      use implicit_solve, only: advance_implicit_terms
      use mirror_terms, only: advance_mirror_implicit
      use flow_shear, only: advance_perp_flow_shear
      use radial_variation_time_advance, only: mb_communicate

      implicit none

      integer, intent(in) :: istep
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g

      ! Start the timer for the implicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 9), ' implicit')

      ! Reverse the order of operations every time step
      ! as part of alternating direction operator splitting
      ! this is needed to ensure 2nd order accuracy in time
      ! if (mod(istep,2)==0) then
      ! g^{*} (coming from explicit solve) is input
      ! get g^{**}, with g^{**}-g^{*} due to mirror term

      if (rk_step) call mb_communicate(g)

      if (mod(istep, 2) == 1 .or. .not. flip_flop) then

         if (prp_shear_enabled) then
            call advance_perp_flow_shear(g)
            fields_updated = .false.
         end if

         if (hyper_dissipation) then
            call advance_hyper_dissipation(g)
            fields_updated = .false.
         end if

         if (collisions_implicit .and. include_collisions) then
            call advance_fields(g, phi, apar, bpar, dist='g')
            call advance_collisions_implicit(mirror_implicit, phi, apar, bpar, g)
            fields_updated = .false.
         end if

         if (mirror_implicit .and. include_mirror) then
            call advance_mirror_implicit(collisions_implicit, g, apar)
            fields_updated = .false.
         end if

         call advance_fields(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 
         ! g^{**} is input
         ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
         if (stream_implicit .and. include_parallel_streaming) then
            call advance_implicit_terms(g, phi, apar, bpar)
            if (radial_variation .or. full_flux_surface) fields_updated = .false.
         end if

         ! Update the fields if not already updated
         call advance_fields(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 
      else

         ! Get updated fields corresponding to advanced g
         ! note that hyper-dissipation and mirror advances
         ! depended only on g and so did not need field update
         call advance_fields(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 

         ! g^{**} is input
         ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
         if (stream_implicit .and. include_parallel_streaming) then
            call advance_implicit_terms(g, phi, apar, bpar)
            if (radial_variation .or. full_flux_surface) fields_updated = .false.
         end if

         if (mirror_implicit .and. include_mirror) then
            call advance_mirror_implicit(collisions_implicit, g, apar)
            fields_updated = .false.
         end if

         if (collisions_implicit .and. include_collisions) then
            call advance_fields(g, phi, apar, bpar, dist='g')
            call advance_collisions_implicit(mirror_implicit, phi, apar, bpar, g)
            fields_updated = .false.
         end if

         if (hyper_dissipation) then
            call advance_hyper_dissipation(g)
            fields_updated = .false.
         end if

         if (prp_shear_enabled) then
            call advance_perp_flow_shear(g)
            fields_updated = .false.
         end if

      end if

      ! Stop the timer for the implict part of the solve
      if (proc0) call time_message(.false., time_gke(:, 9), ' implicit')

   end subroutine advance_implicit

   
   !******************************************************************************
   !                           FINISH TIME ADVANCE SUBROUTINE
   !******************************************************************************
   subroutine finish_time_advance

      use stella_transforms, only: finish_transforms
      use parameters_physics, only: full_flux_surface
      use extended_zgrid, only: finish_extended_zgrid
      use parallel_streaming, only: finish_parallel_streaming
      use mirror_terms, only: finish_mirror
      use flow_shear, only: finish_flow_shear
      use neoclassical_terms, only: finish_neoclassical_terms
      use dissipation, only: finish_dissipation
      use arrays_drifts, only: finish_wstar, finish_wdrift
      use parallel_nonlinearity, only: finish_parallel_nonlinearity
      use radial_variation_time_advance, only: finish_radial_variation

      implicit none

      if (full_flux_surface) call finish_transforms
      call finish_dissipation
      call finish_parallel_nonlinearity 
      call finish_wstar 
      call finish_wdrift 
      call finish_radial_variation
      call finish_parallel_streaming
      call finish_flow_shear
      call finish_mirror
      call finish_neoclassical_terms
      call deallocate_arrays

      time_advance_initialized = .false.
      readinit = .false.

   end subroutine finish_time_advance

   subroutine deallocate_arrays

      use store_arrays_distribution_fn, only: g0, g1, g2, g3

      implicit none

      if (allocated(g0)) deallocate (g0)
      if (allocated(g1)) deallocate (g1)
      if (allocated(g2)) deallocate (g2)
      if (allocated(g3)) deallocate (g3)

   end subroutine deallocate_arrays

end module time_advance
