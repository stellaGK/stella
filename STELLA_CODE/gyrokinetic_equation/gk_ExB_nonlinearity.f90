!###############################################################################
!################# NONLINEAR TERM OF THE GYROKINETIC EQUATION ##################
!###############################################################################
! 
! This module evolves the nonlinear term:
!     (B_r/2) (dy/dalpha) (dx/dpsi) F_k ( F^{-1}_k [ik_y J_0 ϕ_k] F^{-1}_k [ik_x g_{k,s}]
!        - F^{-1}_k [ik_x J_0 ϕ_k] F^{-1}_k [ik_y g_{k,s} ] )
! 
!###############################################################################
module gk_ExB_nonlinearity

   ! Load debug flags
   use debug_flags, only: debug => time_advance_debug

   implicit none

   ! Make routines available to other modules
   public :: advance_ExB_nonlinearity

   private

contains

   !****************************************************************************
   !                         ADVANCE ExB NONLINEARITY                          !
   !****************************************************************************
   subroutine advance_ExB_nonlinearity(g, gout, restart_time_step, istep)

      ! Constants
      use constants, only: pi, zi

      ! Parallelisation
      use timers, only: region_start, region_end
      use mp, only: proc0, min_allreduce
      use mp, only: scope, allprocs, subprocs
      use parallelisation_layouts, only: vmu_lo, imu_idx, is_idx, iv_idx
      use file_utils, only: runtype_option_switch, runtype_multibox

      ! Calculations
      use calculations_timestep, only: reset_dt
      use calculations_tofrom_ghf, only: g_to_h
      use calculations_gyro_averages, only: gyro_average
      use calculations_kxky, only: swap_kxky, swap_kxky_back
      use calculations_kxky_derivatives, only: get_dgdy, get_dgdx
      use calculations_kxky_derivatives, only: get_dchidx, get_dchidy
      use calculations_transforms, only: transform_y2ky, transform_x2kx
      use calculations_transforms, only: transform_y2ky_xfirst, transform_x2kx_xfirst
      
      ! Physics
      use parameters_physics, only: fphi
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: suppress_zonal_interaction
      use parameters_physics, only: full_flux_surface, radial_variation
      use parameters_numerical, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
      
      ! Fields
      use arrays, only: shift_state
      use arrays_fields, only: phi, apar, bpar
      use arrays_fields, only: phi_corr_QN, phi_corr_GA
      use arrays_distribution_function, only: phi_gyro
      
      ! Grids
      use grids_kxky, only: x
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: akx, aky, rho_clamped
      use grids_kxky, only: nakx, ikx_max, naky, naky_all, nx, ny
      use grids_kxky, only: non_twisting, delta_kx
      use grids_time, only: cfl_dt_ExB, cfl_dt_linear, code_dt, code_dt_max
      use grids_velocity, only: vpa, nmu
      use grids_species, only: nspec

      ! Flow shear
      use gk_flow_shear, only: prp_shear_enabled
      use gk_flow_shear, only: hammett_flow_shear
      use gk_flow_shear, only: g_exb, g_exbfac

      ! Geometry
      use geometry, only: exb_nonlin_fac, exb_nonlin_fac_p, gfac
      
      ! TImers

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
      logical, intent(out) :: restart_time_step
      integer, intent(in) :: istep

      ! Local variables
      complex, dimension(:, :), allocatable :: apar_kykx, bpar_kykx
      complex, dimension(:, :), allocatable :: g0k, g0a, g0k_swap
      complex, dimension(:, :), allocatable :: g0kxy, g0xky, prefac, prefac_z
      real, dimension(:, :), allocatable :: dkx_over_ky
      complex, dimension(:, :, :), allocatable :: phase_ntft
      real, dimension(:, :), allocatable :: g0xy, g1xy, bracket
      complex, dimension(:, :), allocatable :: tmp
      integer :: ivmu, iz, it
      real :: zero, cfl_dt
      logical :: yfirst

      ! d<chi>/dx and d<chi>/dy in (y,x) space per (mu, species) point <imus> on this processor, see below
      logical :: chi_by_mu
      integer :: imus, nmus, imu, is, iv
      integer, dimension(:), allocatable :: imus_of_ivmu, ivmu_of_imus
      integer, dimension(:, :), allocatable :: imus_index
      real, dimension(:, :, :, :, :), allocatable :: chix, chiy, chix_apar, chiy_apar
      complex, dimension(:, :), allocatable :: dchidx, dchidy, dchidx_apar, dchidy_apar

      !-------------------------------------------------------------------------

      ! Start timer
      call region_start('exb_nonlinearity')
      if (debug) write (*, *) 'time_advance::solve_gke::advance_ExB_nonlinearity::get_dgdy'

      ! Avoid dividing by zero in cfl_dt terms below
      zero = 100.*epsilon(0.)

      ! Initialize cfl_dt_ExB and restart_time_step
      cfl_dt_ExB = 10000000.
      restart_time_step = .false.
      
      ! By default, prp_shear_enabled = .false. and thus yfirst = .true., and we Fourier transform y first
      ! If perpendicular flow shear is included, it is important to Fourier transform x first.
      ! The non-twisting flux tube also needs the x transform first, to multiply by exp(i ΔKx x) in (ky, x)
      yfirst = .not. (prp_shear_enabled .or. non_twisting)

      ! Allocate arrays
      allocate (apar_kykx(naky, nakx))
      allocate (bpar_kykx(naky, nakx))
      allocate (g0k(naky, nakx))
      allocate (g0a(naky, nakx))
      allocate (g0xy(ny, nx))
      allocate (g1xy(ny, nx))
      allocate (bracket(ny, nx))
      allocate (prefac(naky, nx))
      allocate (prefac_z(naky, nx))
      allocate (tmp(size(gout, 1), size(gout, 2)))

      ! Work arrays of <forward_transform>, here naky_all = 2*naky-1. Both sets are allocated, so that
      ! they can always be passed to <forward_transform>, which gets the private copies of the threads
      ! as arguments (through host association it would see the shared originals)
      allocate (g0k_swap(naky_all, ikx_max))
      allocate (g0kxy(ny, ikx_max))
      allocate (g0xky(naky, nx))

      ! Compute phase factor needed when running with equilibrium flow shear
      prefac = 1.0
      if (prp_shear_enabled .and. hammett_flow_shear) then
         prefac = exp(-zi * g_exb * g_exbfac * spread(x, 1, naky) * spread(aky * shift_state, 2, nx))
      end if

      ! Non-twisting flux tube (Ball and Brunner 2021, eq. B.7-B.9): the modes on the kx grid are
      ! exp(i (akx + ΔKx) x + i ky Y) with ΔKx = <delta_kx>(ky, z), and Y the non-twisting binormal
      ! coordinate. The x derivatives are i (akx + ΔKx), and the phase exp(i ΔKx x) is applied in
      ! (ky, x) space, before the transform to Y (and removed after the transform back). The Poisson
      ! bracket has the same form in (x, Y) as in (x, y). i ΔKx <chi> is obtained from i ky <chi>.
      if (non_twisting) then
         allocate (dkx_over_ky(naky, -nzgrid:nzgrid)); dkx_over_ky = 0.0
         where (spread(abs(aky), 2, 2 * nzgrid + 1) > epsilon(0.0)) &
            dkx_over_ky = delta_kx / spread(aky, 2, 2 * nzgrid + 1)
         allocate (phase_ntft(naky, nx, -nzgrid:nzgrid))
         do iz = -nzgrid, nzgrid
            phase_ntft(:, :, iz) = exp(zi * spread(delta_kx(:, iz), 2, nx) * spread(x, 1, naky))
         end do
      end if

      ! Incoming pdf is g = <f>.
      ! For EM simulations, the pdf entering the ExB nonlinearity needs to be
      ! the non-Boltzmann part of f (h = f + (Ze/T)*phi*F0)
      if (include_apar .or. include_bpar) call g_to_h(g, phi, bpar, fphi)

      ! <chi> = J_0 (phi - 2 vpa (stm) apar) + (J_1 term of bpar) depends on vpa only through its apar part, which
      ! is linear in vpa, so the Fourier transforms of d<chi>/dx and d<chi>/dy are done once for every (mu, species)
      ! on this processor instead of for every vpa: d<chi>/dx(y,x) = chix(:, :, iz, it, imus) + vpa * chix_apar(...),
      ! and the same for y, with the corrections for the non-twisting flux tube, flow shear and the zonal modes
      ! applied to both parts as in the loop below. With radial variation or a full flux surface, <chi> depends
      ! on vpa in other ways as well, and is transformed for every vpa below
      chi_by_mu = .not. (radial_variation .or. full_flux_surface)
      if (chi_by_mu) then
         ! The (mu, species) points on this processor, and one vmu point <ivmu_of_imus> of each
         allocate (imus_index(nmu, nspec)); imus_index = 0
         allocate (imus_of_ivmu(vmu_lo%llim_proc:vmu_lo%ulim_alloc)); imus_of_ivmu = 0
         nmus = 0
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            if (imus_index(imu, is) == 0) then
               nmus = nmus + 1
               imus_index(imu, is) = nmus
            end if
            imus_of_ivmu(ivmu) = imus_index(imu, is)
         end do
         allocate (ivmu_of_imus(max(nmus, 1)))
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            ivmu_of_imus(imus_of_ivmu(ivmu)) = ivmu
         end do
         allocate (chix(ny, nx, -nzgrid:nzgrid, ntubes, max(nmus, 1)))
         allocate (chiy(ny, nx, -nzgrid:nzgrid, ntubes, max(nmus, 1)))
         if (include_apar) then
            allocate (chix_apar(ny, nx, -nzgrid:nzgrid, ntubes, max(nmus, 1)))
            allocate (chiy_apar(ny, nx, -nzgrid:nzgrid, ntubes, max(nmus, 1)))
         else
            allocate (chix_apar(1, 1, 1, 1, 1), chiy_apar(1, 1, 1, 1, 1))
         end if
         allocate (dchidx(naky, nakx), dchidy(naky, nakx), dchidx_apar(naky, nakx), dchidy_apar(naky, nakx))

         !$omp parallel default(none) &
         !$omp firstprivate(nmus, ntubes, nzgrid, nakx, shift_state) &
         !$omp private(imus, it, iz, ivmu, prefac_z, apar_kykx, bpar_kykx, g0k, g0k_swap, g0kxy, g0xky, &
         !$omp         dchidx, dchidy, dchidx_apar, dchidy_apar) &
         !$omp shared(ivmu_of_imus, chix, chiy, chix_apar, chiy_apar, phi, apar, bpar, include_apar, include_bpar, &
         !$omp        non_twisting, phase_ntft, prefac, dkx_over_ky, suppress_zonal_interaction, prp_shear_enabled, &
         !$omp        hammett_flow_shear, g_exb, g_exbfac)
         !$omp do collapse(3)
         do imus = 1, nmus
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  ivmu = ivmu_of_imus(imus)
                  if (non_twisting) then
                     prefac_z = prefac * phase_ntft(:, :, iz)
                  else
                     prefac_z = prefac
                  end if
                  apar_kykx = 0.
                  bpar_kykx = 0.
                  if (include_apar) apar_kykx = apar(:, :, iz, it)
                  if (include_bpar) bpar_kykx = bpar(:, :, iz, it)
                  call get_dchi_vpa_parts(iz, ivmu, phi(:, :, iz, it), apar_kykx, bpar_kykx, &
                                          dchidx, dchidy, dchidx_apar, dchidy_apar)
                  call transform_dchidx(dchidx, dchidy, iz, chix(:, :, iz, it, imus), prefac_z, g0k, g0k_swap, g0kxy, g0xky)
                  call forward_transform(dchidy, chiy(:, :, iz, it, imus), prefac_z, g0k_swap, g0kxy, g0xky)
                  if (include_apar) then
                     call transform_dchidx(dchidx_apar, dchidy_apar, iz, chix_apar(:, :, iz, it, imus), prefac_z, g0k, g0k_swap, &
                                          g0kxy, g0xky)
                     call forward_transform(dchidy_apar, chiy_apar(:, :, iz, it, imus), prefac_z, g0k_swap, g0kxy, g0xky)
                  end if
               end do
            end do
         end do
         !$omp end do
         !$omp end parallel
      end if

      ! Iterate over velocity space
      !$omp parallel default(none) &
      !$omp firstprivate(vmu_lo, shift_state, exb_nonlin_fac, naky,nakx, ny ,rho_clamped, exb_nonlin_fac_p, ikx_max, code_dt, tmp) &
      !$omp reduction(min: cfl_dt_ExB) &
      !$omp private(g0k, g0xy, g0a, g1xy, bracket, g0xky, g0k_swap, g0kxy, apar_kykx, bpar_kykx, prefac_z, imus, iv) &
      !$omp shared(chi_by_mu, imus_of_ivmu, chix, chiy, chix_apar, chiy_apar, vpa, g, gout, include_apar,phi_gyro, include_bpar, gfac, apar, phi, bpar, aky, ntubes, full_flux_surface, suppress_zonal_interaction,radial_variation, prp_shear_enabled, hammett_flow_shear, g_exb, g_exbfac, prefac, yfirst, akx, fphi, phi_corr_GA, phi_corr_qn, zero, nzgrid, non_twisting, phase_ntft, dkx_over_ky, delta_kx)
      !$omp do collapse(3)
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid

               ! Phase factor in (ky, x) space: flow shear and/or the non-twisting flux tube
               if (non_twisting) then
                  prefac_z = prefac * phase_ntft(:, :, iz)
               else
                  prefac_z = prefac
               end if
            
               ! Get the fields (to avoid indexing arrays which aren't allocated)
               if (include_apar) apar_kykx = apar(:, :, iz, it)
               if (include_bpar) bpar_kykx = bpar(:, :, iz, it)
            
               ! Compute i*ky*g
               call get_dgdy(g(:, :, iz, it, ivmu), g0k)
               
               ! Take the FFT to get dg/dy in (y,x) space
               call forward_transform(g0k, g0xy, prefac_z, g0k_swap, g0kxy, g0xky)
               
               ! d<chi>/dx in (y,x) space: from the transforms per (mu, species) above, or transformed here
               if (chi_by_mu) then
                  imus = imus_of_ivmu(ivmu)
                  g1xy = chix(:, :, iz, it, imus)
                  if (include_apar) then
                     iv = iv_idx(vmu_lo, ivmu)
                     g1xy = g1xy + vpa(iv) * chix_apar(:, :, iz, it, imus)
                  end if
               else
                  ! Compute i*kx*<chi>
                  if (full_flux_surface) then
                     call get_dgdx(phi_gyro(:, :, iz, it, ivmu), g0k)
                  else
                     call get_dchidx(iz, ivmu, phi(:, :, iz, it), apar_kykx, bpar_kykx, g0k)
                  end if

                  ! Non-twisting flux tube: add i ΔKx <chi> = (ΔKx/ky) i ky <chi>
                  if (non_twisting) then
                     call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar_kykx, bpar_kykx, g0a)
                     g0k = g0k + spread(dkx_over_ky(:, iz), 2, nakx) * g0a
                  end if
               
                  ! Zero out the zonal contribution to d<chi>/dx if requested
                  if (suppress_zonal_interaction) then
                     g0k(1,:) = 0.0
                  end if
               
                  ! If running with equilibrium flow shear, make adjustment to
                  ! The term multiplying dg/dy
                  if (prp_shear_enabled .and. hammett_flow_shear) then
                     call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar_kykx, bpar_kykx, g0a)
                     g0k = g0k - g_exb * g_exbfac * spread(shift_state, 2, nakx) * g0a
                  end if
               
                  ! Take the FFT to get d<chi>/dx in (y,x) space
                  call forward_transform(g0k, g1xy, prefac_z, g0k_swap, g0kxy, g0xky)
               end if
               
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
                  call forward_transform(g0k, g1xy, prefac_z, g0k_swap, g0kxy, g0xky)
                  g1xy = g1xy * exb_nonlin_fac
                  bracket = bracket + g0xy * g1xy
                  ! Estimate the CFL dt due to the above contribution
                  cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * aky(naky), zero))
               end if

               ! Compute dg/dx in k-space (= i*kx*g); for the non-twisting flux tube i (akx + ΔKx) g
               call get_dgdx(g(:, :, iz, it, ivmu), g0k)
               if (non_twisting) g0k = g0k + zi * spread(delta_kx(:, iz), 2, nakx) * g(:, :, iz, it, ivmu)
               
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
               call forward_transform(g0k, g0xy, prefac_z, g0k_swap, g0kxy, g0xky)
               
               ! d<chi>/dy in (y,x) space: from the transforms per (mu, species) above, or transformed here
               if (chi_by_mu) then
                  g1xy = chiy(:, :, iz, it, imus)
                  if (include_apar) g1xy = g1xy + vpa(iv) * chiy_apar(:, :, iz, it, imus)
               else
                  ! Compute d<chi>/dy in k-space
                  if (full_flux_surface) then
                     call get_dgdy(phi_gyro(:, :, iz, it, ivmu), g0k)
                  else
                     call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar_kykx, bpar_kykx, g0k)
                  end if
               
                  ! Take the FFT to get d<chi>/dy in (y,x) space
                  call forward_transform(g0k, g1xy, prefac_z, g0k_swap, g0kxy, g0xky)
               end if
               
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
                  call forward_transform(g0k, g1xy, prefac_z, g0k_swap, g0kxy, g0xky)
                  g1xy = g1xy * exb_nonlin_fac
                  bracket = bracket - g0xy * g1xy
                  ! Estimate the CFL dt due to the above contribution
                  cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * akx(ikx_max), zero))
               end if

               if (yfirst) then
                  call transform_x2kx(bracket, g0kxy)
                  if (full_flux_surface) then
                     gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + code_dt * g0kxy
                  else
                     call transform_y2ky(g0kxy, g0k_swap)
                     call swap_kxky_back(g0k_swap, tmp)
                     gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + code_dt * tmp
                  end if
               else
                  call transform_y2ky_xfirst(bracket, g0xky)
                  g0xky = g0xky / prefac_z
                  call transform_x2kx_xfirst(g0xky, tmp)
                  gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + code_dt * tmp
               end if
            end do
         end do
      end do
      !$omp end do
      !$omp end parallel

      ! Convert back from h to g = <f> (only needed for EM sims)
      if (include_apar .or. include_bpar) call g_to_h(g, phi, bpar, -fphi)

      deallocate (g0k, g0a, g0xy, g1xy, bracket)
      deallocate (prefac, prefac_z)
      if (allocated(dkx_over_ky)) deallocate (dkx_over_ky, phase_ntft)
      if (allocated(g0k_swap)) deallocate (g0k_swap)
      if (allocated(g0xky)) deallocate (g0xky)
      if (allocated(g0kxy)) deallocate (g0kxy)
      if (chi_by_mu) deallocate (imus_index, imus_of_ivmu, ivmu_of_imus, chix, chiy, chix_apar, chiy_apar, &
                                 dchidx, dchidy, dchidx_apar, dchidy_apar)

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
      end if

      call region_end('exb_nonlinearity')

   contains

      !-------------------------------------------------------------------------
      ! d<chi>/dx (or one of its parts) <dx> in (y,x) space, with the same corrections as in the loop over the
      ! vmu points: the non-twisting flux tube, the zonal modes and the flow shear (from d<chi>/dy, <dy>).
      ! Everything that is private in the OpenMP loop is passed as an argument (<iz>, the phase factor and the
      ! work arrays), since host association would give the shared originals instead of the private copies
      subroutine transform_dchidx(dx, dy, iz, gx, phase, gk, gk_swap, gkxy, gxky)

         implicit none

         complex, dimension(:, :), intent(in) :: dx, dy
         integer, intent(in) :: iz
         real, dimension(:, :), intent(out) :: gx
         complex, dimension(:, :), intent(in) :: phase
         complex, dimension(:, :), intent(out) :: gk, gk_swap, gkxy, gxky

         !----------------------------------------------------------------------

         gk = dx
         if (non_twisting) gk = gk + spread(dkx_over_ky(:, iz), 2, nakx) * dy
         if (suppress_zonal_interaction) gk(1, :) = 0.0
         if (prp_shear_enabled .and. hammett_flow_shear) gk = gk - g_exb * g_exbfac * spread(shift_state, 2, nakx) * dy
         call forward_transform(gk, gx, phase, gk_swap, gkxy, gxky)

      end subroutine transform_dchidx

      !-------------------------------------------------------------------------
      subroutine forward_transform(gk, gx, phase, gk_swap, gkxy, gxky)

         use calculations_transforms, only: transform_ky2y, transform_kx2x
         use calculations_transforms, only: transform_ky2y_xfirst, transform_kx2x_xfirst

         implicit none
   
         ! Arguments
         complex, dimension(:, :), intent(in) :: gk
         real, dimension(:, :), intent(out) :: gx
         ! Phase factor in (ky, x) space and work arrays of this thread: they are private in the
         ! OpenMP loop, so they are passed explicitly instead of through host association, and the
         ! transforms of different threads can run at the same time (the work arrays of the FFTs in
         ! calculations_transforms are threadprivate)
         complex, dimension(:, :), intent(in) :: phase
         complex, dimension(:, :), intent(out) :: gk_swap, gkxy, gxky

         !----------------------------------------------------------------------
         
         ! Fourier transform g(kx,ky) to g(kx,y) and then to g(x,y)
         if (yfirst) then
               ! We have i*ky*g(kx,ky) for ky >= 0 and all kx.
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
               call swap_kxky(gk, gk_swap)
               call transform_ky2y(gk_swap, gkxy)
               call transform_kx2x(gkxy, gx)
               
         ! Fourier transform g(kx,ky) to g(x,ky) and then to g(x,y)
         else
               call transform_kx2x_xfirst(gk, gxky)
               gxky = gxky * phase
               call transform_ky2y_xfirst(gxky, gx)
         end if
      end subroutine forward_transform

   end subroutine advance_ExB_nonlinearity

   !****************************************************************************
   !                 PARTS OF d<chi>/dx AND d<chi>/dy WITHOUT vpa
   !****************************************************************************
   ! At the (mu, species) point of <ivmu>, d<chi>/dx = <dchidx> + vpa * <dchidx_apar>, and the same for y,
   ! with d<chi>/dx and d<chi>/dy as in get_dchidx and get_dchidy: <dchidx> holds the phi and bpar parts,
   ! which do not depend on vpa, and <dchidx_apar> the apar part, -2 (stm) J_0 i kx apar, without the vpa
   !****************************************************************************
   subroutine get_dchi_vpa_parts(iz, ivmu, phi, apar, bpar, dchidx, dchidy, dchidx_apar, dchidy_apar)

      use constants, only: zi
      use parallelisation_layouts, only: vmu_lo, is_idx
      use calculations_gyro_averages, only: gyro_average
      use calculations_kxky_derivatives, only: get_dchidx, get_dchidy
      use parameters_physics, only: include_apar
      use grids_species, only: spec
      use grids_kxky, only: naky, nakx, akx, aky

      implicit none

      integer, intent(in) :: iz, ivmu
      complex, dimension(:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :), intent(out) :: dchidx, dchidy, dchidx_apar, dchidy_apar

      complex, dimension(naky, nakx) :: no_apar, field
      integer :: is

      !-------------------------------------------------------------------------

      ! The phi and bpar parts
      no_apar = 0.
      call get_dchidx(iz, ivmu, phi, no_apar, bpar, dchidx)
      call get_dchidy(iz, ivmu, phi, no_apar, bpar, dchidy)

      ! The apar part, without the factor vpa
      if (include_apar) then
         is = is_idx(vmu_lo, ivmu)
         field = -2.0 * spec(is)%stm_psi0 * apar
         call gyro_average(zi * spread(akx, 1, naky) * field, iz, ivmu, dchidx_apar)
         call gyro_average(zi * spread(aky, 2, nakx) * field, iz, ivmu, dchidy_apar)
      else
         dchidx_apar = 0.
         dchidy_apar = 0.
      end if

   end subroutine get_dchi_vpa_parts

end module gk_ExB_nonlinearity
