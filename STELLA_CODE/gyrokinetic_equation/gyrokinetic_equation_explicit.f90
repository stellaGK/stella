!###############################################################################
!################## EXPLICIT TERMS OF THE GYROKINETIC EQUATION #################
!###############################################################################
! 
! This module adds the explicit gyrokinetic terms to the right-hand-side of the 
! gyrokinetic equation, in order to advance the distribution function in time.
! 
! First advance the distribution function <g> in time using the gyrokinetic equation.
! Next, adance the fields (electrostatic potential <phi>, as well as the electromagnetic
! fields <apar> and <bpar>) in time using the quasi-neutrality condition.
! 
!###############################################################################
module gyrokinetic_equation_explicit

   ! Load debug flags
   use debug_flags, only: debug => time_advance_debug

   implicit none
   
   ! Make routines available to stella.f90
   public :: advance_distribution_function_using_explicit_gyrokinetic_terms

   private

contains

   !****************************************************************************
   !                       EXPLICIT TIME ADVANCE SUBROUTINES
   !****************************************************************************
   ! Advance the guiding center distribution equation <g> in k-space in time using 
   ! the gyrokinetic equation. Specifically, the distribution function <g> is updated 
   ! based on all of the terms in the GKE that are advanced explicitly in time.
   !****************************************************************************
   subroutine advance_distribution_function_using_explicit_gyrokinetic_terms(g, restart_time_step, istep)

      ! Parallelisation
      use mp, only: proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo, iv_idx
      use timers, only: time_gke
      
      ! Fields
      use parameters_physics, only: include_apar
      use arrays_fields, only: phi, apar, bpar
      use field_equations, only: advance_fields
      
      ! Grids
      use grids_z, only: nzgrid
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: periodic
      use grids_extended_zgrid, only: phase_shift
      use gk_parallel_streaming, only: stream_sign
      
      ! Calculations
      use calculations_tofrom_ghf, only: gbar_to_g, g_or_gbar_to_gbarneo
      
      ! Numerical time advance schemes
      use parameters_numerical, only: explicit_algorithm_switch
      use parameters_numerical, only: explicit_algorithm_rk3
      use parameters_numerical, only: explicit_algorithm_rk2
      use parameters_numerical, only: explicit_algorithm_rk4
      use parameters_numerical, only: explicit_algorithm_euler

      ! For NEO's neoclassical corrections.
      use neoclassical_terms_neo, only: neoclassical_is_enabled

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      ! Local variables
      integer :: ivmu, iv, sgn, iky

      !-------------------------------------------------------------------------

      ! Start the timer for the explicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 8), ' explicit')
      
      ! If the fields are not already updated, then update them
      if (include_apar) then
         call advance_fields(g, phi, apar, bpar, dist='g')
      end if

      ! If NEO's higher order corrections are included,  convert from g to gbarneo, as gbarneo appears in time derivatives
      if (neoclassical_is_enabled()) then
         call g_or_gbar_to_gbarneo(g, phi, apar, bpar, 1.0)
      end if

      ! Incoming distribution function is g = h - (Z F0/T) (J0 phi + 4 mu (T/Z) (J1/gamma) bpar)
      ! If <include_apar> = T, convert from g to gbar = g + Z F0/T (2J0 vpa vth apar), as gbar appears in time derivatives
      if (include_apar) then
         call gbar_to_g(g, apar, -1.0)
      end if

      ! Use a numerical time advance scheme to advance the distribution function
      ! in time, based on the explicit terms in the gyrokinetic equation. 
      ! Choose between Forward Euler, or 2nd, 3rd or 4rd order Runge-Kutta schemes
      ! The 3rd order Runge-Kutta scheme is the default option.
      select case (explicit_algorithm_switch)
      case (explicit_algorithm_euler)
         call advance_explicit_euler(g, restart_time_step, istep)
      case (explicit_algorithm_rk2)
         call advance_explicit_rk2(g, restart_time_step, istep)
      case (explicit_algorithm_rk3)
         call advance_explicit_rk3(g, restart_time_step, istep)
      case (explicit_algorithm_rk4)
         call advance_explicit_rk4(g, restart_time_step, istep)
      end select

      ! If the fields are not already updated, then update them
      if (include_apar) then
         call advance_fields(g, phi, apar, bpar, dist='gbar')
      end if

      ! Later, the implicit solve will use <g> rather than <gbar> to advance the 
      ! distribution function in time. Therefore, convert <gbar> to <g> again
      if (include_apar) then
         call gbar_to_g(g, apar, 1.0)
      end if

      ! We now switch back to g. 
      if (neoclassical_is_enabled()) then
         call g_or_gbar_to_gbarneo(g, phi, apar, bpar, -1.0)
      end if

      ! Enforce periodicity for periodic (including zonal) modes
      ! <stream_sign> > 0 corresponds to dz/dt < 0
      do iky = 1, naky
         if (periodic(iky)) then
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               sgn = stream_sign(iv)
               g(iky, :, sgn * nzgrid, :, ivmu) = g(iky, :, -sgn * nzgrid, :, ivmu) * phase_shift(iky)**(-sgn)
            end do
         end if
      end do

      ! Stop the timer for the explicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 8), ' explicit')

   end subroutine advance_distribution_function_using_explicit_gyrokinetic_terms
   
   !****************************************************************************
   !                NEEDED FOR ALL EXPLICIT TIME ADVANCE SUBROUTINES
   !****************************************************************************
   ! solve_gyrokinetic_equation_explicit accepts as argument pdf, the guiding center distribution function in k-space,
   ! and returns rhs_ky, the right-hand side of the gyrokinetic equation in k-space;
   ! i.e., if dg/dt = r, then rhs_ky = r*dt;
   ! note that if include_apar = T, then the input pdf is actually gbar = g + (Ze/T)*(vpa/c)*<Apar>*F0
   !****************************************************************************
   subroutine add_explicit_gyrokinetic_terms(pdf, rhs_ky, restart_time_step, istep)

      ! Parallelisation
      use job_manage, only: time_message
      use multibox, only: add_multibox_krook
      use parallelisation_layouts, only: vmu_lo
      use calculations_transforms, only: transform_y2ky

      ! Fields
      use arrays_fields, only: phi, apar, bpar
      
      ! Distribution function
      use arrays_distribution_function, only: phi_gyro
      
      ! Calculations
      use arrays_gyro_averages, only: j0_ffs
      use calculations_gyro_averages, only: gyro_average
      use calculations_kxky, only: swap_kxky_back
      use calculations_tofrom_ghf, only: gbar_to_g

      ! Physics flags
      use parameters_physics, only: include_parallel_nonlinearity
      use parameters_physics, only: include_parallel_streaming
      use parameters_physics, only: include_mirror
      use parameters_physics, only: include_apar
      use parameters_physics, only: include_nonlinear
      use parameters_physics, only: full_flux_surface
      use parameters_physics, only: radial_variation
      use parameters_physics, only: xdriftknob, ydriftknob
      use dissipation_and_collisions, only: include_collisions
      use dissipation_and_collisions, only: hyper_dissipation
      use gk_flow_shear, only: g_exb
      
      ! Krook source
      use parameters_multibox, only: include_multibox_krook
      use gk_sources, only: source_option_switch
      use gk_sources, only: source_option_krook
      use gk_sources, only: add_krook_operator
      
      ! Numerical schemes
      use parameters_numerical, only: stream_implicit
      use parameters_numerical, only: mirror_implicit
      use parameters_numerical, only: drifts_implicit
      use dissipation_and_collisions, only: advance_collisions_explicit
      use dissipation_and_collisions, only: collisions_implicit

      ! Grids
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: zonal_mode, akx
      use grids_kxky, only: ikx_max, ny, naky_all
      
      ! Routines to add terms to the gyrokinetic equation
      use gk_parallel_streaming, only: advance_parallel_streaming_explicit
      use gk_mirror, only: advance_mirror_explicit
      use gk_flow_shear, only: advance_parallel_flow_shear
      use gk_drive, only: advance_wstar_explicit
      use gk_magnetic_drift, only: advance_wdriftx_explicit, advance_wdrifty_explicit
      use gk_nonlinearity, only: advance_parallel_nonlinearity
      use gk_radial_variation, only: advance_radial_variation
      use gk_nonlinearity, only: advance_ExB_nonlinearity
      use field_equations_radialvariation, only: get_radial_correction
      use field_equations, only: advance_fields
      use field_equations, only: fields_updated

      ! For advancing NEO's neoclassical corrections explicitly. 

      use neoclassical_terms_neo, only: neoclassical_is_enabled
      use gk_neo_chi_terms, only: advance_neo_chi_terms_explicit
      use gk_neo_apar_terms, only: advance_neo_apar_terms_explicit
      use gk_neo_dchidz_terms, only: advance_neo_dchidz_terms_explicit
      use gk_neo_drive, only: advance_wstar1_explicit, advance_wpol_explicit
      use gk_neo_drifts, only: advance_neo_mag_drift_explicit, advance_neo_curv_drift_explicit
 
      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: pdf
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out), target :: rhs_ky
      logical, intent(out) :: restart_time_step
      integer, intent(in) :: istep

      ! Local variables
      complex, dimension(:, :, :, :, :), allocatable, target :: rhs_y
      complex, dimension(:, :, :, :, :), pointer :: rhs
      complex, dimension(:, :), allocatable :: rhs_ky_swap
      integer :: iz, it, ivmu


      !-------------------------------------------------------------------------

      ! Initialise the right-hand-side of the gyrokinetic equation to zero
      rhs_ky = 0.

      ! If full_flux_surface = .true., then initially obtain the RHS of the GKE in alpha-space;
      ! We will later inverse Fourier transform to get RHS in k_alpha-space
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
      if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_fields'

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

      ! Obtain the gyro-average of the electrostatic potential phi and store in phi_gyro;
      ! this can be a particularly costly operation when simulating a full flux surface
      ! due to the coupling of different k-alphas inherent in the gyro-average;
      ! calculate once here to avoid repeated calculation later
      ! TODO-GA : can this be spec up??
      if (full_flux_surface) call gyro_average(phi, phi_gyro, j0_ffs)

      !! INSERT TEST HERE TO SEE IF dg/dy, dg/dx, d<phi>/dy, d<phi>/dx WILL BE NEEDED
      !! IF SO, PRE-COMPUTE ONCE HERE

      ! Default is to continue with same time step size.
      ! if estimated CFL condition for nonlinear terms is violated
      ! then restart_time_step will be set to .true.
      restart_time_step = .false.
      
      ! Calculate and add ExB nonlinearity to RHS of GK eqn
      ! do this first, as the CFL condition may require a change in time step
      ! and thus recomputation of mirror, wdrift, wstar, and parstream
      if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_ExB_nonlinearity'
      if (include_nonlinear) call advance_ExB_nonlinearity(pdf, rhs, restart_time_step, istep)

      ! Include contribution from the parallel nonlinearity (aka turbulent acceleration)
      if (include_parallel_nonlinearity .and. .not. restart_time_step) &
         call advance_parallel_nonlinearity(pdf, rhs, restart_time_step)

      if (.not. restart_time_step) then

         ! Include contribution from perp flow shear in the parallel component of the toroidal flow
         if ((g_exb**2) > epsilon(0.0)) call advance_parallel_flow_shear(rhs)

         ! Calculate and add mirror term to RHS of GK eqn
         if (include_mirror .and. .not. mirror_implicit) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_mirror_explicit'
            call advance_mirror_explicit(pdf, rhs)
         end if

         if (.not. drifts_implicit) then
            ! Calculate and add alpha-component of magnetic drift term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_wdrifty_explicit'
            if (abs(ydriftknob) > epsilon(0.0)) then
               call advance_wdrifty_explicit(pdf, phi, bpar, rhs)
            end if

            ! Calculate and add psi-component of magnetic drift term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_wdriftx_explicit'
            if (abs(xdriftknob) > epsilon(0.0)) then
               call advance_wdriftx_explicit(pdf, phi, bpar, rhs)
            end if
            
            ! Calculate and add omega_* term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_wstar_explicit'
            call advance_wstar_explicit(phi, rhs) 

            ! =========================================================================== !
            ! If NEO's corrections are included, then...

            if (neoclassical_is_enabled()) then
                 ! Advance the neoclassical chi terms.                 
                 ! call advance_neo_chi_terms_explicit(phi, rhs)

                 ! If apar is switched on, we must advance the neoclassical apar terms. 
                 if (include_apar) then
                     call advance_neo_apar_terms_explicit(rhs)
                 end if
 
                 ! Advance the neoclassical dchi/dz terms.
                 ! call advance_neo_dchidz_terms_explicit(phi, rhs)
 
                 ! Advance the neoclassical equilibrium gradient drive terms. 
                 call advance_wstar1_explicit(phi, rhs)
                 call advance_wpol_explicit(phi, rhs)

                 ! Advance the neoclassical magnetic and curvature drift terms.
                 ! call advance_neo_mag_drift_explicit(phi, rhs)
                 ! call advance_neo_curv_drift_explicit(phi, rhs)
            end if
            ! ============================================================================ !

         end if
 
         ! Calculate and add contribution from collisions to RHS of GK eqn
         if (include_collisions .and. .not. collisions_implicit) call advance_collisions_explicit(pdf, phi, bpar, rhs)

         ! Calculate and add parallel streaming term to RHS of GK eqn
         if (include_parallel_streaming .and. (.not. stream_implicit)) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_parallel_streaming_explicit'
            call advance_parallel_streaming_explicit(pdf, phi, bpar, rhs)
         end if
         
         if (hyper_dissipation) then
            call advance_hyper_explicit(pdf, rhs)
         end if

         ! If simulating a full flux surface (flux annulus), all terms to this point have been calculated
         ! in real-space in alpha (y); transform to kalpha (ky) space before adding to RHS of GKE.
         ! NB: it may be that for fully explicit calculation, this transform can be eliminated with additional code changes
         if (full_flux_surface) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::transform_y2ky'
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

   end subroutine add_explicit_gyrokinetic_terms

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine advance_hyper_explicit(gin, gout)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, nakx
      use dissipation_hyper, only: advance_hyper_vpa, advance_hyper_zed
      use dissipation_hyper, only: hyp_zed, hyp_vpa

      implicit none
      
      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      ! Local variables
      complex, dimension(:, :, :, :, :), allocatable :: dg

      !-------------------------------------------------------------------------

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
   
   
!###############################################################################
!######################## EXPLICIT TIME ADVANCE SCHEMES ########################
!###############################################################################

   !****************************************************************************
   !                        EXPLICIT EULER TIME ADVANCE SUBROUTINE
   !****************************************************************************
   ! Uses forward Euler to advance one time step.
   !****************************************************************************
   subroutine advance_explicit_euler(g, restart_time_step, istep)
   
      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      use gk_radial_variation, only: mb_communicate
      use parameters_multibox, only: rk_step

      ! Distribution function
      use arrays_distribution_function, only: g0
      
      ! Grids
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      !-------------------------------------------------------------------------

      ! <rk_step> is only true if running in multibox mode
      if (rk_step) call mb_communicate(g)

      ! Store the distribution function
      g0 = g

      ! Calculate the right-hand-side of the gyrokinetic equation and store it in <g>
      call add_explicit_gyrokinetic_terms(g0, g, restart_time_step, istep)

      ! Add the right-hand-side of the gyrokinetic equation to <g0>
      g = g0 + g

   end subroutine advance_explicit_euler

   !****************************************************************************
   !                         EXPLICIT RK2 TIME ADVANCE SUBROUTINE
   !****************************************************************************
   ! Uses strong stability-preserving rk2 to advance one time step.
   !****************************************************************************
   subroutine advance_explicit_rk2(g, restart_time_step, istep)
   
      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step
      use gk_radial_variation, only: mb_communicate
      
      ! Distribution function
      use arrays_distribution_function, only: g0, g1

      ! Grids
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      ! Local variables
      integer :: counter

      !-------------------------------------------------------------------------

      ! <rk_step> is only true if running in multibox mode
      if (rk_step) call mb_communicate(g)

      ! Store the distribution function
      g0 = g
      
      ! Count the iterations
      counter = 1

      ! SSP rk2 algorithm to advance explicit part of code
      ! if GK equation written as dg/dt = rhs - vpar . grad h,
      ! add_explicit_gyrokinetic_terms returns rhs*dt
      do while (counter <= 2)
      
         ! Second order Runge-Kutta Scheme
         select case (counter)
         case (1)
            call add_explicit_gyrokinetic_terms(g0, g1, restart_time_step, istep)
         case (2)
            g1 = g0 + g1
            if (rk_step) call mb_communicate(g1)
            call add_explicit_gyrokinetic_terms(g1, g, restart_time_step, istep)
         end select
         
         ! If the code_dt is reset, we need to quit this loop and restart the timestep again
         ! Otherwise increase the <counter> that keeps track of the iterations
         if (restart_time_step) then 
            counter = 10
         else
            counter = counter + 1
         end if
         
      end do

      ! This is g at intermediate time level
      g = 0.5 * g0 + 0.5 * (g1 + g)

   end subroutine advance_explicit_rk2

   !****************************************************************************
   !                         EXPLICIT RK3 TIME ADVANCE SUBROUTINE
   !****************************************************************************
   ! Uses strong stability-preserving rk3 to advance one time step.
   !****************************************************************************
   subroutine advance_explicit_rk3(g, restart_time_step, istep)
   
      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step
      use gk_radial_variation, only: mb_communicate

      ! Distribution function
      use arrays_distribution_function, only: g0, g1, g2
      
      ! Grids
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      ! Local variables
      integer :: counter

      !-------------------------------------------------------------------------

      ! rk_STEP = false unless in multibox mode
      if (rk_step) call mb_communicate(g)

      ! Store the distribution function
      g0 = g
      
      ! Count the iterations
      counter = 1

      ! SSP rk3 algorithm to advance explicit part of code
      ! if GK equation written as dg/dt = rhs - vpar . grad h,
      ! add_explicit_gyrokinetic_terms returns rhs*dt
      do while (counter <= 3)
      
         ! Third order Runge-Kutta Scheme
         select case (counter)
         case (1)
            call add_explicit_gyrokinetic_terms(g0, g1, restart_time_step, istep)
         case (2)
            g1 = g0 + g1
            if (rk_step) call mb_communicate(g1)
            call add_explicit_gyrokinetic_terms(g1, g2, restart_time_step, istep)
         case (3)
            g2 = g1 + g2
            if (rk_step) call mb_communicate(g2)
            call add_explicit_gyrokinetic_terms(g2, g, restart_time_step, istep)
         end select
         
         ! If the code_dt is reset, we need to quit this loop and restart the timestep again
         ! Otherwise increase the <counter> that keeps track of the iterations
         if (restart_time_step) then
            counter = 10
         else
            counter = counter + 1
         end if
      end do

      ! This is g at intermediate time level
      g = g0 / 3.+0.5 * g1 + (g2 + g) / 6.

   end subroutine advance_explicit_rk3

   !****************************************************************************
   !                         EXPLICIT RK4 TIME ADVANCE SUBROUTINE
   !****************************************************************************
   ! Uses rk4 to advance one time step.
   !****************************************************************************
   subroutine advance_explicit_rk4(g, restart_time_step, istep)
   
      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step
      use gk_radial_variation, only: mb_communicate

      ! Distribution function
      use arrays_distribution_function, only: g0, g1, g2, g3
      
      ! Grids
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      ! Local variables
      integer :: counter

      !-------------------------------------------------------------------------

      ! rk_step is false unless in multibox mode
      if (rk_step) call mb_communicate(g)

      ! Store the distribution function
      g0 = g
      
      ! Count the iterations
      counter = 1

      ! RK4 algorithm to advance explicit part of code
      ! if GK equation written as dg/dt = rhs - vpar . grad h,
      ! add_explicit_gyrokinetic_terms returns rhs*dt
      do while (counter <= 4)
      
         ! Fourth order Runge-Kutta Scheme
         select case (counter)
         case (1)
            call add_explicit_gyrokinetic_terms(g0, g1, restart_time_step, istep)
         case (2)
            g3 = g0 + 0.5 * g1 ! g1 is h*k1
            if (rk_step) call mb_communicate(g3)
            call add_explicit_gyrokinetic_terms(g3, g2, restart_time_step, istep)
            g1 = g1 + 2.*g2
         case (3)
            g2 = g0 + 0.5 * g2 ! g2 is h*k2
            if (rk_step) call mb_communicate(g2)
            call add_explicit_gyrokinetic_terms(g2, g3, restart_time_step, istep)
            g1 = g1 + 2.*g3
         case (4)
            g3 = g0 + g3 ! g3 is h*k3
            if (rk_step) call mb_communicate(g3)
            call add_explicit_gyrokinetic_terms(g3, g, restart_time_step, istep)
            g1 = g1 + g
         end select
         
         ! If the code_dt is reset, we need to quit this loop and restart the timestep again
         ! Otherwise increase the <counter> that keeps track of the iterations
         if (restart_time_step) then
            counter = 10
         else
            counter = counter + 1
         end if
         
      end do

      ! This is g at intermediate time level
      g = g0 + g1 / 6.

   end subroutine advance_explicit_rk4

end module gyrokinetic_equation_explicit
