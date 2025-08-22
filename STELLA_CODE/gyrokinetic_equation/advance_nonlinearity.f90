module advance_nonlinearity

    use debug_flags, only: debug => time_advance_debug

    implicit none

    public :: advance_ExB_nonlinearity

    private

contains

    !****************************************************************************
    !                         ADVANCE ExB NONLINEARITY                          !
    !****************************************************************************
    subroutine advance_ExB_nonlinearity(g, gout, restart_time_step, istep)

        use mp, only: proc0, min_allreduce
        use mp, only: scope, allprocs, subprocs
        use job_manage, only: time_message
        use constants, only: pi, zi
        use file_utils, only: runtype_option_switch, runtype_multibox

        use stella_transforms, only: transform_y2ky, transform_x2kx
        use stella_transforms, only: transform_y2ky_xfirst, transform_x2kx_xfirst
        use stella_time, only: cfl_dt_ExB, cfl_dt_linear, code_dt, code_dt_max
        use stella_layouts, only: vmu_lo, imu_idx, is_idx
        
        use g_tofrom_h, only: g_to_h
        use gyro_averages, only: gyro_average
        use calculations_kxky_derivatives, only: get_dchidx, get_dchidy
        use calculations_kxky_derivatives, only: get_dgdy, get_dgdx
        use calculations_kxky, only: swap_kxky, swap_kxky_back

        use store_arrays_fields, only: phi, apar, bpar, shift_state
        use store_arrays_fields, only: phi_corr_QN, phi_corr_GA
        use store_arrays_distribution_fn, only: g_scratch
        use store_arrays_useful, only: time_gke
        
        use parameters_numerical, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
        use parameters_physics, only: g_exb, g_exbfac, fphi
        use parameters_physics, only: full_flux_surface, radial_variation
        use parameters_physics, only: prp_shear_enabled, hammett_flow_shear
        use parameters_physics, only: include_apar, include_bpar
        use parameters_physics, only: suppress_zonal_interaction
        use parameters_kxky_grid, only: nakx, ikx_max, naky, naky_all, nx, ny

        use z_grid, only: nzgrid, ntubes
        use grids_kxky, only: akx, aky, rho_clamped
        use grids_kxky, only: x
        use geometry, only: exb_nonlin_fac, exb_nonlin_fac_p, gfac

         use reset_timestep, only: reset_dt

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
        logical, intent(out) :: restart_time_step
        integer, intent(in) :: istep

        complex, dimension(:, :), allocatable :: g0k, g0a, g0k_swap
        complex, dimension(:, :), allocatable :: g0kxy, g0xky, prefac
        real, dimension(:, :), allocatable :: g0xy, g1xy, bracket

        real :: zero, cfl_dt
        integer :: ivmu, iz, it, imu, is
        logical :: yfirst

        ! alpha-component of magnetic drift (requires ky -> y)
        if (proc0) call time_message(.false., time_gke(:, 7), ' ExB nonlinear advance')

        if (debug) write (*, *) 'time_advance::solve_gke::advance_ExB_nonlinearity::get_dgdy'

        ! avoid divide by zero in cfl_dt terms below
        zero = 100.*epsilon(0.)

        ! Initialize cfl_dt_ExB
        cfl_dt_ExB = 10000000.

        restart_time_step = .false.
        ! this statement seems to imply that flow shear is not compatible with FFS
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

        !> compute phase factor needed when running with equilibrium flow shear
        prefac = 1.0
        if (prp_shear_enabled .and. hammett_flow_shear) then
            prefac = exp(-zi * g_exb * g_exbfac * spread(x, 1, naky) * spread(aky * shift_state, 2, nx))
        end if

        ! incoming pdf is g = <f>
        ! for EM simulations, the pdf entering the ExB nonlinearity needs to be
        ! the non-Boltzmann part of f (h = f + (Ze/T)*phi*F0)
        if (include_apar .or. include_bpar) call g_to_h(g, phi, bpar, fphi)

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
                do iz = -nzgrid, nzgrid
                !> compute i*ky*g
                call get_dgdy(g(:, :, iz, it, ivmu), g0k)
                !> FFT to get dg/dy in (y,x) space
                call forward_transform(g0k, g0xy)
                !> compute i*kx*<chi>
                if (full_flux_surface) then
                    call get_dgdx(g_scratch(:, :, iz, it, ivmu), g0k)
                else
                    call get_dchidx(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k)
                end if
                !> zero out the zonal contribution to d<chi>/dx if requested
                if (suppress_zonal_interaction) then
                    g0k(1,:) = 0.0
                end if
                !> if running with equilibrium flow shear, make adjustment to
                !> the term multiplying dg/dy
                if (prp_shear_enabled .and. hammett_flow_shear) then
                    call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0a)
                    g0k = g0k - g_exb * g_exbfac * spread(shift_state, 2, nakx) * g0a
                end if
                !> FFT to get d<chi>/dx in (y,x) space
                call forward_transform(g0k, g1xy)
                !> multiply by the geometric factor appearing in the Poisson bracket;
                !> i.e., (dx/dpsi*dy/dalpha)*0.5
                g1xy = g1xy * exb_nonlin_fac
                !> compute the contribution to the Poisson bracket from dg/dy*d<chi>/dx
                bracket = g0xy * g1xy

                !> estimate the CFL dt due to the above contribution
                cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * aky(naky), zero))

                if (radial_variation) then
                    bracket = bracket + gfac * g0xy * g1xy * exb_nonlin_fac_p * spread(rho_clamped, 1, ny)
                    call gyro_average(phi_corr_QN(:, :, iz, it), iz, ivmu, g0a)
                    g0a = fphi * (g0a + phi_corr_GA(:, :, iz, it, ivmu))
                    call get_dgdx(g0a, g0k)
                    call forward_transform(g0k, g1xy)
                    g1xy = g1xy * exb_nonlin_fac
                    bracket = bracket + g0xy * g1xy
                    !> estimate the CFL dt due to the above contribution
                    cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * aky(naky), zero))
                end if

                !> compute dg/dx in k-space (= i*kx*g)
                call get_dgdx(g(:, :, iz, it, ivmu), g0k)
                !> zero out the zonal contribution to dg/dx if requested
                if (suppress_zonal_interaction) then
                    g0k(1,:) = 0.0
                end if
                !> if running with equilibrium flow shear, correct dg/dx term
                if (prp_shear_enabled .and. hammett_flow_shear) then
                    call get_dgdy(g(:, :, iz, it, ivmu), g0a)
                    g0k = g0k - g_exb * g_exbfac * spread(shift_state, 2, nakx) * g0a
                end if
                !> FFT to get dg/dx in (y,x) space
                call forward_transform(g0k, g0xy)
                !> compute d<chi>/dy in k-space
                if (full_flux_surface) then
                    call get_dgdy(g_scratch(:, :, iz, it, ivmu), g0k)
                else
                    call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k)
                end if
                !> FFT to get d<chi>/dy in (y,x) space
                call forward_transform(g0k, g1xy)
                !> multiply by the geometric factor appearing in the Poisson bracket;
                !> i.e., (dx/dpsi*dy/dalpha)*0.5
                g1xy = g1xy * exb_nonlin_fac
                !> compute the contribution to the Poisson bracket from dg/dy*d<chi>/dx
                bracket = bracket - g0xy * g1xy

                !> estimate the CFL dt due to the above contribution
                cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * akx(ikx_max), zero))

                if (radial_variation) then
                    bracket = bracket - gfac * g0xy * g1xy * exb_nonlin_fac_p * spread(rho_clamped, 1, ny)
                    call gyro_average(phi_corr_QN(:, :, iz, it), iz, ivmu, g0a)
                    g0a = fphi * (g0a + phi_corr_GA(:, :, iz, it, ivmu))
                    call get_dgdy(g0a, g0k)
                    call forward_transform(g0k, g1xy)
                    g1xy = g1xy * exb_nonlin_fac
                    bracket = bracket - g0xy * g1xy
                    !> estimate the CFL dt due to the above contribution
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

        ! convert back from h to g = <f> (only needed for EM sims)
        if (include_apar .or. include_bpar) call g_to_h(g, phi, bpar, -fphi)

        deallocate (g0k, g0a, g0xy, g1xy, bracket)
        if (allocated(g0k_swap)) deallocate (g0k_swap)
        if (allocated(g0xky)) deallocate (g0xky)
        if (allocated(g0kxy)) deallocate (g0kxy)

        if (runtype_option_switch == runtype_multibox) call scope(allprocs)

        call min_allreduce(cfl_dt_ExB)

        if (runtype_option_switch == runtype_multibox) call scope(subprocs)

        !> check estimated cfl_dt to see if the time step size needs to be changed
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

        subroutine forward_transform(gk, gx)

            use stella_transforms, only: transform_ky2y, transform_kx2x
            use stella_transforms, only: transform_ky2y_xfirst, transform_kx2x_xfirst

            implicit none

            complex, dimension(:, :), intent(in) :: gk
            real, dimension(:, :), intent(out) :: gx

            if (yfirst) then
                ! we have i*ky*g(kx,ky) for ky >= 0 and all kx
                ! want to do 1D complex to complex transform in y
                ! which requires i*ky*g(kx,ky) for all ky and kx >= 0
                ! use g(kx,-ky) = conjg(g(-kx,ky))
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

end module advance_nonlinearity