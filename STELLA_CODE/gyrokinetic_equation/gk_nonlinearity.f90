!###############################################################################
!################# NONLINEAR TERM OF THE GYROKINETIC EQUATION ##################
!###############################################################################
! 
! This module evolves the nonlinear term:
!     (B_r/2) (dy/dalpha) (dx/dpsi) F_k ( F^{-1}_k [ik_y J_0 ϕ_k] F^{-1}_k [ik_x g_{k,s}]
!        - F^{-1}_k [ik_x J_0 ϕ_k] F^{-1}_k [ik_y g_{k,s} ] )
! 
!###############################################################################
module gk_nonlinearity

   ! Load debug flags
   use debug_flags, only: debug => time_advance_debug

   implicit none

   ! Make routines available to other modules
   public :: init_parallel_nonlinearity
   public :: advance_ExB_nonlinearity
   public :: advance_parallel_nonlinearity
   public :: finish_parallel_nonlinearity

   private

   real, dimension(:, :), allocatable :: par_nl_fac, d_par_nl_fac_dr
   real, dimension(:, :), allocatable :: par_nl_curv, d_par_nl_curv_dr
   real, dimension(:), allocatable :: par_nl_driftx, par_nl_drifty
   real, dimension(:), allocatable :: d_par_nl_driftx_dr, d_par_nl_drifty_dr

contains
      
   !****************************************************************************
   !                     INITIALISE PARALLEL NONLINEARITY                      !
   !****************************************************************************
   subroutine init_parallel_nonlinearity

      use parameters_physics, only: rhostar
      use parameters_physics, only: ydriftknob
      use parameters_physics, only: radial_variation

      use grids_species, only: spec, nspec
      use grids_z, only: nztot, nzgrid
      use arrays, only: initialised_parallel_streaming

      use geometry, only: geo_surf, drhodpsi, q_as_x
      use geometry, only: gradpar, dbdzed, bmag
      use geometry, only: B_times_kappa_dot_grady, B_times_kappa_dot_gradx
      use geometry, only: dIdrho, dgradpardrho, dBdrho, d2Bdrdth
      use geometry, only: dcvdriftdrho, dcvdrift0drho
      
      implicit none

      !-------------------------------------------------------------------------

      if (.not. allocated(par_nl_fac)) allocate (par_nl_fac(-nzgrid:nzgrid, nspec))
      ! This is the factor multiplying -dphi/dz * dg/dvpa in the parallel nonlinearity
      par_nl_fac = 0.5 * rhostar * spread(spec%stm_psi0 * spec%zt_psi0, 1, nztot) * spread(gradpar, 2, nspec)

      if (.not. allocated(par_nl_curv)) allocate (par_nl_curv(-nzgrid:nzgrid, nspec))
      ! ydriftknob is here because this term comes from bhat x curvature . grad B
      par_nl_curv = -ydriftknob * rhostar * geo_surf%rgeo * geo_surf%betaprim * drhodpsi &
                     * spread(dbdzed(1, :) * gradpar / bmag(1, :), 2, nspec) / spread(spec%zt_psi0, 1, nztot)

      if (.not. allocated(par_nl_drifty)) allocate (par_nl_drifty(-nzgrid:nzgrid))
      par_nl_drifty = 0.5 * rhostar * B_times_kappa_dot_grady(1, :)
      if (.not. allocated(par_nl_driftx)) allocate (par_nl_driftx(-nzgrid:nzgrid))
      if (q_as_x) then
         par_nl_driftx = 0.25 * rhostar * B_times_kappa_dot_gradx(1, :) * 2. * geo_surf%shat
      else
         par_nl_driftx = 0.25 * rhostar * B_times_kappa_dot_gradx(1, :) * 2.
      end if

      if (radial_variation) then
         if (.not. allocated(d_par_nl_fac_dr)) allocate (d_par_nl_fac_dr(-nzgrid:nzgrid, nspec))
         ! This is the factor multiplying -dphi/dz * dg/dvpa in the parallel nonlinearity
         d_par_nl_fac_dr = 0.5 * rhostar * spread(spec%stm_psi0 * spec%zt_psi0, 1, nztot) * spread(dgradpardrho, 2, nspec)

         if (.not. allocated(d_par_nl_curv_dr)) allocate (d_par_nl_curv_dr(-nzgrid:nzgrid, nspec))
         ! ydriftknob is here because this term comes from bhat x curvature . grad B
         ! Handle terms with no zeroes
         d_par_nl_curv_dr = par_nl_curv * (dIdrho / geo_surf%rgeo - drhodpsi * geo_surf%d2psidr2 &
                                          - spread(dBdrho / bmag(1, :) + dgradpardrho / gradpar, 2, nspec))
         ! Handle terms with possible zeroes
         d_par_nl_curv_dr = d_par_nl_curv_dr &
                              - ((ydriftknob * rhostar * geo_surf%rgeo * drhodpsi * spread(gradpar / bmag(1, :), 2, nspec)) &
                              / spread(spec%zt_psi0, 1, nztot)) &
                              * (geo_surf%betadbprim * spread(dbdzed(1, :), 2, nspec) &
                              + geo_surf%betaprim * spread(d2Bdrdth, 2, nspec))

         if (.not. allocated(d_par_nl_drifty_dr)) allocate (d_par_nl_drifty_dr(-nzgrid:nzgrid))
         d_par_nl_drifty_dr = 0.25 * rhostar * dcvdriftdrho(1, :)
         if (.not. allocated(d_par_nl_drifty_dr)) allocate (d_par_nl_driftx_dr(-nzgrid:nzgrid))
         if (q_as_x) then
               d_par_nl_driftx_dr = 0.25 * rhostar * dcvdrift0drho(1, :)
         else
               d_par_nl_driftx_dr = 0.25 * rhostar * dcvdrift0drho(1, :) / geo_surf%shat
         end if
      end if

      initialised_parallel_streaming = .true.

   end subroutine init_parallel_nonlinearity

   !****************************************************************************
   !                         ADVANCE ExB NONLINEARITY                          !
   !****************************************************************************
   subroutine advance_ExB_nonlinearity(g, gout, restart_time_step, istep)

      use constants, only: pi, zi

      use job_manage, only: time_message
      use mp, only: proc0, min_allreduce
      use mp, only: scope, allprocs, subprocs
      use parallelisation_layouts, only: vmu_lo, imu_idx, is_idx
      use file_utils, only: runtype_option_switch, runtype_multibox

      use calculations_timestep, only: reset_dt
      use calculations_tofrom_ghf, only: g_to_h
      use calculations_gyro_averages, only: gyro_average
      use calculations_kxky, only: swap_kxky, swap_kxky_back
      use calculations_kxky_derivatives, only: get_dgdy, get_dgdx
      use calculations_kxky_derivatives, only: get_dchidx, get_dchidy
      use calculations_transforms, only: transform_y2ky, transform_x2kx
      use calculations_transforms, only: transform_y2ky_xfirst, transform_x2kx_xfirst
      
      use parameters_physics, only: fphi
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: suppress_zonal_interaction
      use parameters_physics, only: full_flux_surface, radial_variation
      use parameters_numerical, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
      
      use timers, only: time_gke
      use arrays, only: shift_state
      use arrays_fields, only: phi, apar, bpar
      use arrays_fields, only: phi_corr_QN, phi_corr_GA
      use arrays_distribution_function, only: g_scratch
      
      use grids_kxky, only: x
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: akx, aky, rho_clamped
      use grids_kxky, only: nakx, ikx_max, naky, naky_all, nx, ny
      use grids_time, only: cfl_dt_ExB, cfl_dt_linear, code_dt, code_dt_max

      use gk_flow_shear, only: prp_shear_enabled, hammett_flow_shear
      use gk_flow_shear, only: g_exb, g_exbfac

      use geometry, only: exb_nonlin_fac, exb_nonlin_fac_p, gfac

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
      logical, intent(out) :: restart_time_step
      integer, intent(in) :: istep

      ! Local variables
      complex, dimension(:, :), allocatable :: g0k, g0a, g0k_swap
      complex, dimension(:, :), allocatable :: g0kxy, g0xky, prefac
      real, dimension(:, :), allocatable :: g0xy, g1xy, bracket
      real :: zero, cfl_dt
      integer :: ivmu, iz, it, imu, is
      logical :: yfirst

      !-------------------------------------------------------------------------

      ! Alpha-component of magnetic drift (requires ky -> y)
      if (proc0) call time_message(.false., time_gke(:, 7), ' ExB nonlinear advance')

      if (debug) write (*, *) 'time_advance::solve_gke::advance_ExB_nonlinearity::get_dgdy'

      ! Avoid divide by zero in cfl_dt terms below
      zero = 100.*epsilon(0.)

      ! Initialize cfl_dt_ExB
      cfl_dt_ExB = 10000000.

      restart_time_step = .false.
      ! This statement seems to imply that flow shear is not compatible with FFS 
      ! need to check
      yfirst = .not. prp_shear_enabled

      allocate (g0k(naky, nakx))
      allocate (g0a(naky, nakx))
      allocate (g0xy(ny, nx))
      allocate (g1xy(ny, nx))
      allocate (bracket(ny, nx))
      allocate (prefac(naky, nx))

      if (yfirst) then
         allocate (g0k_swap(naky_all, ikx_max))
         allocate (g0kxy(ny, ikx_max))
      else
         allocate (g0xky(naky, nx))
      end if

      ! Compute phase factor needed when running with equilibrium flow shear
      prefac = 1.0
      if (prp_shear_enabled .and. hammett_flow_shear) then
         prefac = exp(-zi * g_exb * g_exbfac * spread(x, 1, naky) * spread(aky * shift_state, 2, nx))
      end if

      ! Incoming pdf is g = <f>. 
      ! For EM simulations, the pdf entering the ExB nonlinearity needs to be
      ! the non-Boltzmann part of f (h = f + (Ze/T)*phi*F0)
      if (include_apar .or. include_bpar) call g_to_h(g, phi, bpar, fphi)

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
               do iz = -nzgrid, nzgrid
               ! Compute i*ky*g
               call get_dgdy(g(:, :, iz, it, ivmu), g0k)
               ! Take the FFT to get dg/dy in (y,x) space
               call forward_transform(g0k, g0xy)
               ! Compute i*kx*<chi>
               if (full_flux_surface) then
                  call get_dgdx(g_scratch(:, :, iz, it, ivmu), g0k)
               else
                  call get_dchidx(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k)
               end if
               ! Zero out the zonal contribution to d<chi>/dx if requested
               if (suppress_zonal_interaction) then
                  g0k(1,:) = 0.0
               end if
               ! If running with equilibrium flow shear, make adjustment to
               ! The term multiplying dg/dy
               if (prp_shear_enabled .and. hammett_flow_shear) then
                  call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0a)
                  g0k = g0k - g_exb * g_exbfac * spread(shift_state, 2, nakx) * g0a
               end if
               ! Take the FFT to get d<chi>/dx in (y,x) space
               call forward_transform(g0k, g1xy)
               ! Multiply by the geometric factor appearing in the Poisson bracket;
               ! i.e., (dx/dpsi*dy/dalpha)*0.5
               g1xy = g1xy * exb_nonlin_fac
               ! Compute the contribution to the Poisson bracket from dg/dy*d<chi>/dx
               bracket = g0xy * g1xy

               ! Estimate the CFL dt due to the above contribution
               cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * aky(naky), zero))

               if (radial_variation) then
                  bracket = bracket + gfac * g0xy * g1xy * exb_nonlin_fac_p * spread(rho_clamped, 1, ny)
                  call gyro_average(phi_corr_QN(:, :, iz, it), iz, ivmu, g0a)
                  g0a = fphi * (g0a + phi_corr_GA(:, :, iz, it, ivmu))
                  call get_dgdx(g0a, g0k)
                  call forward_transform(g0k, g1xy)
                  g1xy = g1xy * exb_nonlin_fac
                  bracket = bracket + g0xy * g1xy
                  ! Estimate the CFL dt due to the above contribution
                  cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * aky(naky), zero))
               end if

               ! Compute dg/dx in k-space (= i*kx*g)
               call get_dgdx(g(:, :, iz, it, ivmu), g0k)
               ! Zero out the zonal contribution to dg/dx if requested
               if (suppress_zonal_interaction) then
                  g0k(1,:) = 0.0
               end if
               ! If running with equilibrium flow shear, correct dg/dx term
               if (prp_shear_enabled .and. hammett_flow_shear) then
                  call get_dgdy(g(:, :, iz, it, ivmu), g0a)
                  g0k = g0k - g_exb * g_exbfac * spread(shift_state, 2, nakx) * g0a
               end if
               ! Take the FFT to get dg/dx in (y,x) space
               call forward_transform(g0k, g0xy)
               ! Compute d<chi>/dy in k-space
               if (full_flux_surface) then
                  call get_dgdy(g_scratch(:, :, iz, it, ivmu), g0k)
               else
                  call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k)
               end if
               ! Take the FFT to get d<chi>/dy in (y,x) space
               call forward_transform(g0k, g1xy)
               ! Multiply by the geometric factor appearing in the Poisson bracket;
               ! i.e., (dx/dpsi*dy/dalpha)*0.5
               g1xy = g1xy * exb_nonlin_fac
               ! Compute the contribution to the Poisson bracket from dg/dy*d<chi>/dx
               bracket = bracket - g0xy * g1xy

               ! Estimate the CFL dt due to the above contribution
               cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * akx(ikx_max), zero))

               if (radial_variation) then
                  bracket = bracket - gfac * g0xy * g1xy * exb_nonlin_fac_p * spread(rho_clamped, 1, ny)
                  call gyro_average(phi_corr_QN(:, :, iz, it), iz, ivmu, g0a)
                  g0a = fphi * (g0a + phi_corr_GA(:, :, iz, it, ivmu))
                  call get_dgdy(g0a, g0k)
                  call forward_transform(g0k, g1xy)
                  g1xy = g1xy * exb_nonlin_fac
                  bracket = bracket - g0xy * g1xy
                  ! Estimate the CFL dt due to the above contribution
                  cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * akx(ikx_max), zero))
               end if

               if (yfirst) then
                  call transform_x2kx(bracket, g0kxy)
                  if (full_flux_surface) then
                     gout(:, :, iz, it, ivmu) = g0kxy
                  else
                     call transform_y2ky(g0kxy, g0k_swap)
                     call swap_kxky_back(g0k_swap, gout(:, :, iz, it, ivmu))
                  end if
               else
                  call transform_y2ky_xfirst(bracket, g0xky)
                  g0xky = g0xky / prefac
                  call transform_x2kx_xfirst(g0xky, gout(:, :, iz, it, ivmu))
               end if
               end do
         end do
      end do

      ! Convert back from h to g = <f> (only needed for EM sims)
      if (include_apar .or. include_bpar) call g_to_h(g, phi, bpar, -fphi)

      deallocate (g0k, g0a, g0xy, g1xy, bracket)
      if (allocated(g0k_swap)) deallocate (g0k_swap)
      if (allocated(g0xky)) deallocate (g0xky)
      if (allocated(g0kxy)) deallocate (g0kxy)

      if (runtype_option_switch == runtype_multibox) call scope(allprocs)

      call min_allreduce(cfl_dt_ExB)

      if (runtype_option_switch == runtype_multibox) call scope(subprocs)

      ! Check estimated cfl_dt to see if the time step size needs to be changed
      cfl_dt = min(cfl_dt_ExB, cfl_dt_linear)
      if (code_dt > cfl_dt * cfl_cushion_upper) then
         if (proc0) then
               write (*, *) ' '
               write (*, '(A30,I0,A1)') 'CHANGING TIME STEP: (istep = ', istep, ')'
               write (*, '(A22, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
               write (*, '(A22, ES10.2E2)') "   cfl_dt_ExB:"//REPEAT(' ', 50), cfl_dt_ExB
               write (*, '(A22, ES10.2E2)') "   cfl_dt_linear:"//REPEAT(' ', 50), cfl_dt_linear
               write (*, '(A22, ES10.2E2)') "   cfl_cushion_upper:"//REPEAT(' ', 50), cfl_cushion_upper
               write (*, '(A22, ES10.2E2)') "   cfl_cushion_middle:"//REPEAT(' ', 50), cfl_cushion_middle
               write (*, '(A22, ES10.2E2)') "   cfl_cushion_lower:"//REPEAT(' ', 50), cfl_cushion_lower
               write (*, '(A62)') '     ==> The code_dt is larger than cfl_dt*cfl_cushion_upper.'
               write (*, '(A59,ES11.4)') '      ==> Decreasing code_dt to cfl_dt*cfl_cushion_middle =', cfl_dt * cfl_cushion_middle
               write (*, *) ' '
         end if
         code_dt = cfl_dt * cfl_cushion_middle
         call reset_dt
         restart_time_step = .true.
      else if (code_dt < min(cfl_dt * cfl_cushion_lower, code_dt_max)) then
         if (proc0) then
               write (*, *) ' '
               write (*, '(A30,I0,A1)') 'CHANGING TIME STEP: (istep = ', istep, ')'
               write (*, '(A22, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
               write (*, '(A22, ES10.2E2)') "   cfl_dt_ExB:"//REPEAT(' ', 50), cfl_dt_ExB
               write (*, '(A22, ES10.2E2)') "   cfl_dt_linear:"//REPEAT(' ', 50), cfl_dt_linear
               write (*, '(A22, ES10.2E2)') "   cfl_cushion_upper:"//REPEAT(' ', 50), cfl_cushion_upper
               write (*, '(A22, ES10.2E2)') "   cfl_cushion_middle:"//REPEAT(' ', 50), cfl_cushion_middle
               write (*, '(A22, ES10.2E2)') "   cfl_cushion_lower:"//REPEAT(' ', 50), cfl_cushion_lower
               write (*, '(A63)') '     ==> The code_dt is smaller than cfl_dt*cfl_cushion_lower.'
               write (*, '(A59,ES11.4)') '      ==> Increasing code_dt to cfl_dt*cfl_cushion_middle =', cfl_dt * cfl_cushion_middle
               write (*, *) ' '
         end if
         code_dt = min(cfl_dt * cfl_cushion_middle, code_dt_max)
         call reset_dt
         restart_time_step = .true.
      else
         gout = code_dt * gout
      end if

      if (proc0) call time_message(.false., time_gke(:, 7), ' ExB nonlinear advance')

   contains

      !-------------------------------------------------------------------------
      subroutine forward_transform(gk, gx)

         use calculations_transforms, only: transform_ky2y, transform_kx2x
         use calculations_transforms, only: transform_ky2y_xfirst, transform_kx2x_xfirst

         implicit none
   
         ! Arguments
         complex, dimension(:, :), intent(in) :: gk
         real, dimension(:, :), intent(out) :: gx

         !----------------------------------------------------------------------

         if (yfirst) then
               ! Ee have i*ky*g(kx,ky) for ky >= 0 and all kx.
               ! We want to do 1D complex to complex transform in y, 
               ! which requires i*ky*g(kx,ky) for all ky and kx >= 0 .
               ! Use the reality condition: g(kx,-ky) = conjg(g(-kx,ky))
               ! so i*(-ky)*g(kx,-ky) = -i*ky*conjg(g(-kx,ky)) = conjg(i*ky*g(-kx,ky))
               ! and i*kx*g(kx,-ky) = i*kx*conjg(g(-kx,ky)) = conjg(i*(-kx)*g(-kx,ky))
               ! and i*(-ky)*J0(kx,-ky)*phi(kx,-ky) = conjg(i*ky*J0(-kx,ky)*phi(-kx,ky))
               ! and i*kx*J0(kx,-ky)*phi(kx,-ky) = conjg(i*(-kx)*J0(-kx,ky)*phi(-kx,ky))
               ! i.e., can calculate dg/dx, dg/dy, d<phi>/dx and d<phi>/dy
               ! on stella (kx,ky) grid, then conjugate and flip sign of (kx,ky)
               ! NB: J0(kx,ky) = J0(-kx,-ky)
               ! TODO DSO: coordinate change for shearing
               call swap_kxky(gk, g0k_swap)
               call transform_ky2y(g0k_swap, g0kxy)
               call transform_kx2x(g0kxy, gx)
         else
               call transform_kx2x_xfirst(gk, g0xky)
               g0xky = g0xky * prefac
               call transform_ky2y_xfirst(g0xky, gx)
         end if

      end subroutine forward_transform

   end subroutine advance_ExB_nonlinearity

   !****************************************************************************
   !                       ADVANCE PARALLEL NONLINEARITY                       !
   !****************************************************************************
   subroutine advance_parallel_nonlinearity(g, gout, restart_time_step)

      ! Constants + Paralellisation
      use constants, only: zi
      use initialise_redistribute, only: xyz2vmu
      use job_manage, only: time_message
      use mp, only: scope, allprocs, subprocs
      use mp, only: proc0, min_allreduce, mp_abort
      use redistribute, only: gather, scatter
      use parallelisation_layouts, only: vmu_lo, xyz_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      use file_utils, only: runtype_option_switch, runtype_multibox
      
      ! Parameters
      use parameters_physics, only: full_flux_surface, radial_variation
      use parameters_numerical, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
      
      ! Grids + Arrays
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: dvpa, vpa, mu
      use grids_extended_zgrid, only: periodic
      use grids_z, only: nzgrid, delzed, ntubes
      use grids_kxky, only: akx, aky, rho_clamped
      use grids_kxky, only: nakx, naky, nx, ny, ikx_max
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: fill_zed_ghost_zones
      use grids_extended_zgrid, only: neigen, nsegments, ikxmod
      use grids_time, only: cfl_dt_parallel, cfl_dt_linear, code_dt, code_dt_max
      use timers, only: time_parallel_nl
      use arrays_fields, only: phi, phi_corr_QN, phi_corr_GA
      
      ! Calculations
      use calculations_timestep, only: reset_dt
      use calculations_gyro_averages, only: gyro_average
      use calculations_kxky, only: swap_kxky, swap_kxky_back
      use calculations_finite_differences, only: third_order_upwind
      use calculations_finite_differences, only: second_order_centered_zed
      use calculations_transforms, only: transform_ky2y, transform_y2ky
      use calculations_transforms, only: transform_kx2x, transform_x2kx
      
      ! GK equation
      use gk_parallel_streaming, only: stream_sign

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
      logical, intent(out) :: restart_time_step

      ! Local variables
      integer :: ivmu, ixyz
      integer :: iz, it, iv, imu, is
      integer :: iky, ie, iseg
      integer :: advect_sign
      real :: cfl_dt
      real, dimension(:), allocatable :: dgdv
      real, dimension(:, :, :, :, :), allocatable :: g0xy
      real, dimension(:, :, :), allocatable :: gxy_vmulocal
      real, dimension(:, :), allocatable :: g1xy, advect_speed
      complex, dimension(2) :: gleft, gright
      complex, dimension(:, :, :, :), allocatable :: phi_gyro, dphidz
      complex, dimension(:, :), allocatable :: g0k, g0kxy, g0k_swap
      complex, dimension(:, :), allocatable :: tmp

      !-------------------------------------------------------------------------
      
      ! WARNING this routine will probably break if neigen_max = 0
      ! Which happens when be set grid_option = 'range'

      ! Alpha-component of magnetic drift (requires ky -> y)
      if (proc0) call time_message(.false., time_parallel_nl(:, 1), ' parallel nonlinearity advance')

      ! Initialize cfl_dt_parallel
      cfl_dt_parallel = 10000000.

      restart_time_step = .false.

      ! Overview:
      ! ---------
      ! Need g and d<phi>/dz in (x,y) space in order to upwind dg/dvpa
      ! 
      ! 1) transform d<phi>/dz from (kx,ky) to (x,y). layout: vmu_lo
      ! 2) need sign of parnl advection in xyz_lo (since dg/dvpa
      !    requires vpa local), so d<phi>/dz(vmu_lo) --> d<phi>/dz(xyz_lo)
      ! 3) transform g from (kx,ky) to (x,y). layout: vmu_lo
      ! 4) dg/dvpa requires vpa local, so g(vmu_lo) --> g(xyz_lo)
      ! 5) calculate dg/dvpa
      ! 6) multiply dg/dvpa with d<phi>/dz
      ! 7) product(xyz_lo) --> product(vmu_lo)
      ! 8) inverse transform product(vmu_lo)

      allocate (g0k(naky, nakx))
      allocate (g0xy(ny, nx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      allocate (g0kxy(ny, ikx_max))
      if (radial_variation) allocate (g1xy(ny, nx))
      allocate (phi_gyro(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dphidz(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g0k_swap(2 * naky - 1, ikx_max))
      allocate (tmp(size(gout, 1), size(gout, 2)))

      ! Get d<phi>/dz in vmu_lo
      ! we will need to transform it to real-space
      ! as its sign is needed for upwinding of dg/dvpa
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         ! Construct <phi>
         dphidz = phi
         if (radial_variation) dphidz = dphidz + phi_corr_QN
         call gyro_average(dphidz, ivmu, phi_gyro)
         if (radial_variation) phi_gyro = phi_gyro + phi_corr_GA(:, :, :, :, ivmu)

         do iky = 1, naky
               do it = 1, ntubes
               do ie = 1, neigen(iky)
                  do iseg = 1, nsegments(ie, iky)
                     ! First fill in ghost zones at boundaries in g(z)
                     call fill_zed_ghost_zones(it, iseg, ie, iky, phi_gyro, gleft, gright)
                     ! Now get d<phi>/dz
                     call second_order_centered_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
                        phi_gyro(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                        delzed(0), stream_sign(iv), gleft, gright, periodic(iky), &
                        dphidz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
                  end do
               end do
               end do
         end do

         if (radial_variation) then
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  ! Use reality to swap from ky >= 0, all kx to kx >= 0 , all ky
                  call swap_kxky(dphidz(:, :, iz, it), g0k_swap)
                  ! Transform in y
                  call transform_ky2y(g0k_swap, g0kxy)
                  ! Transform in x
                  call transform_kx2x(g0kxy, g1xy)
                  g0xy(:, :, iz, it, ivmu) = g1xy * (par_nl_fac(iz, is) + d_par_nl_fac_dr(iz, is) * spread(rho_clamped, 1, ny))

                  g0k = zi * spread(aky, 2, nakx) * phi_gyro(:, :, iz, it)
                  call swap_kxky(g0k, g0k_swap)
                  call transform_ky2y(g0k_swap, g0kxy)
                  call transform_kx2x(g0kxy, g1xy)
                  g0xy(:, :, iz, it, ivmu) = g0xy(:, :, iz, it, ivmu) &
                     + vpa(iv) * g1xy * (par_nl_drifty(iz) + d_par_nl_drifty_dr(iz) * spread(rho_clamped, 1, ny))

                  g0k = zi * spread(akx, 1, naky) * phi_gyro(:, :, iz, it)
                  call swap_kxky(g0k, g0k_swap)
                  call transform_ky2y(g0k_swap, g0kxy)
                  call transform_kx2x(g0kxy, g1xy)
                  g0xy(:, :, iz, it, ivmu) = g0xy(:, :, iz, it, ivmu) &
                     + vpa(iv) * g1xy * (par_nl_driftx(iz) + d_par_nl_driftx_dr(iz) * spread(rho_clamped, 1, ny))

                  g0xy(:, :, iz, it, ivmu) = g0xy(:, :, iz, it, ivmu) &
                     + vpa(iv) * mu(imu) * (par_nl_curv(iz, is) + d_par_nl_curv_dr(iz, is) * spread(rho_clamped, 1, ny))

               end do
            end do
         else
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g0k = dphidz(:, :, iz, it) * par_nl_fac(iz, is) + vpa(iv) * mu(imu) * par_nl_curv(iz, is) &
                           + zi * vpa(iv) * phi_gyro(:, :, iz, it) * (spread(akx, 1, naky) * par_nl_driftx(iz) &
                                                                  + spread(aky, 2, nakx) * par_nl_drifty(iz))
                  ! Use reality to swap from ky >= 0, all kx to kx >= 0 , all ky
                  call swap_kxky(g0k, g0k_swap)
                  ! Transform in y
                  call transform_ky2y(g0k_swap, g0kxy)
                  ! Transform in x
                  call transform_kx2x(g0kxy, g0xy(:, :, iz, it, ivmu))
               end do
            end do
         end if
      end do

      ! Do not need phi_gyro or dphidz  again so deallocate
      deallocate (phi_gyro, dphidz)
      deallocate (g0k)
      if (allocated(g1xy)) deallocate (g1xy)

      allocate (gxy_vmulocal(nvpa, nmu, xyz_lo%llim_proc:xyz_lo%ulim_alloc))
      allocate (advect_speed(nmu, xyz_lo%llim_proc:xyz_lo%ulim_alloc))

      ! We now have the advection velocity in vpa in (x,y) space.
      ! Next redistribute it so that (vpa,mu) are local
      if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
      call scatter(xyz2vmu, g0xy, gxy_vmulocal)
      if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
      ! advect_speed does not depend on vpa
      advect_speed = gxy_vmulocal(1, :, :)

      ! Transform g from (kx,ky) to (x,y)
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
               do iz = -nzgrid, nzgrid
               call swap_kxky(g(:, :, iz, it, ivmu), g0k_swap)
               ! Transform in y
               call transform_ky2y(g0k_swap, g0kxy)
               ! Transform in x
               call transform_kx2x(g0kxy, g0xy(:, :, iz, it, ivmu))
               end do
         end do
      end do

      ! Redistribute so that (vpa,mu) local
      if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
      call scatter(xyz2vmu, g0xy, gxy_vmulocal)
      if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')

      allocate (dgdv(nvpa))

      ! We now need to form dg/dvpa and obtain product of dg/dvpa with advection speed
      do ixyz = xyz_lo%llim_proc, xyz_lo%ulim_proc
         do imu = 1, nmu
            ! advect_sign set to +/- 1 depending on sign of the parallel nonlinearity
            ! advection velocity
            ! NB: advect_sign = -1 corresponds to positive advection velocity
            advect_sign = int(sign(1.0, advect_speed(imu, ixyz)))
            call third_order_upwind(1, gxy_vmulocal(:, imu, ixyz), dvpa, advect_sign, dgdv)
            gxy_vmulocal(:, imu, ixyz) = dgdv * advect_speed(imu, ixyz)
            cfl_dt_parallel = min(cfl_dt_parallel, dvpa / abs(advect_speed(imu, ixyz)))
         end do
      end do

      ! Finished with dgdv and advect_speed
      deallocate (dgdv, advect_speed)

      ! Now that we have the full parallel nonlinearity in (x,y)-space.
      ! Need to redistribute so that (x,y) local for transforms
      if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
      call gather(xyz2vmu, gxy_vmulocal, g0xy)
      if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')

      ! Finished with gxy_vmulocal - deallocate
      deallocate (gxy_vmulocal)

      ! g0xy is parallel nonlinearity term with (x,y) on processor
      ! Need to inverse Fourier transform
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
               do iz = -nzgrid, nzgrid
               call transform_x2kx(g0xy(:, :, iz, it, ivmu), g0kxy)
               if (full_flux_surface) then
                  gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + code_dt * g0kxy
               else
                  call transform_y2ky(g0kxy, g0k_swap)
                  call swap_kxky_back(g0k_swap, tmp)
                  gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + code_dt * tmp
               end if
               end do
         end do
      end do
      deallocate (g0k_swap, g0kxy, g0xy)

      if (runtype_option_switch == runtype_multibox) call scope(allprocs)

      call min_allreduce(cfl_dt_parallel)

      if (runtype_option_switch == runtype_multibox) call scope(subprocs)

      ! Check estimated cfl_dt to see if the time step size needs to be changed
      cfl_dt = min(cfl_dt_parallel, cfl_dt_linear)
      if (code_dt > cfl_dt * cfl_cushion_upper) then
         if (proc0) then
               write (*, *) ' '
               write (*, *) 'CHANGING TIME STEP:'
               write (*, '(A22, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
               write (*, '(A22, ES10.2E2)') "   cfl_dt_parallel:"//REPEAT(' ', 50), cfl_dt_parallel
               write (*, '(A22, ES10.2E2)') "   cfl_dt_linear:"//REPEAT(' ', 50), cfl_dt_linear
               write (*, '(A22, ES10.2E2)') "   cfl_cushion_upper:"//REPEAT(' ', 50), cfl_cushion_upper
               write (*, '(A22, ES10.2E2)') "   cfl_cushion_middle:"//REPEAT(' ', 50), cfl_cushion_middle
               write (*, '(A22, ES10.2E2)') "   cfl_cushion_lower:"//REPEAT(' ', 50), cfl_cushion_lower
               write (*, '(A62)') '     ==> The code_dt is larger than cfl_dt*cfl_cushion_upper.'
               write (*, '(A59,ES11.4)') '      ==> Decreasing code_dt to cfl_dt*cfl_cushion_middle =', cfl_dt * cfl_cushion_middle
               write (*, *) ' '
         end if
         code_dt = cfl_dt * cfl_cushion_middle
         call reset_dt
         restart_time_step = .true.
      else if (code_dt < min(cfl_dt * cfl_cushion_lower, code_dt_max)) then
         if (proc0) then
               write (*, *) ' '
               write (*, *) 'CHANGING TIME STEP:'
               write (*, '(A22, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
               write (*, '(A22, ES10.2E2)') "   cfl_dt_parallel:"//REPEAT(' ', 50), cfl_dt_parallel
               write (*, '(A22, ES10.2E2)') "   cfl_dt_linear:"//REPEAT(' ', 50), cfl_dt_linear
               write (*, '(A22, ES10.2E2)') "   cfl_cushion_upper:"//REPEAT(' ', 50), cfl_cushion_upper
               write (*, '(A22, ES10.2E2)') "   cfl_cushion_middle:"//REPEAT(' ', 50), cfl_cushion_middle
               write (*, '(A22, ES10.2E2)') "   cfl_cushion_lower:"//REPEAT(' ', 50), cfl_cushion_lower
               write (*, '(A63)') '     ==> The code_dt is smaller than cfl_dt*cfl_cushion_lower.'
               write (*, '(A59,ES11.4)') '      ==> Increasing code_dt to cfl_dt*cfl_cushion_middle =', cfl_dt * cfl_cushion_middle
               write (*, *) ' '
         end if
         code_dt = min(cfl_dt * cfl_cushion_middle, code_dt_max)
         call reset_dt
         restart_time_step = .true.
      end if

      if (proc0) call time_message(.false., time_parallel_nl(:, 1), ' parallel nonlinearity advance')

   end subroutine advance_parallel_nonlinearity

   !****************************************************************************
   !                      FINALISE PARALLEL NONLINEARITY                       !
   !****************************************************************************
   subroutine finish_parallel_nonlinearity

      use arrays, only: initialised_parallel_streaming

      implicit none

      if (allocated(par_nl_fac)) deallocate (par_nl_fac)
      if (allocated(par_nl_curv)) deallocate (par_nl_curv)
      if (allocated(par_nl_driftx)) deallocate (par_nl_driftx)
      if (allocated(par_nl_drifty)) deallocate (par_nl_drifty)

      initialised_parallel_streaming = .false.

   end subroutine finish_parallel_nonlinearity

end module gk_nonlinearity
