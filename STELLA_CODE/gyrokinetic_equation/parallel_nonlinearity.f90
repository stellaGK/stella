module parallel_nonlinearity

    implicit none

    public :: init_parallel_nonlinearity
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
        use species, only: spec, nspec
        use z_grid, only: nztot, nzgrid
        use geometry, only: geo_surf, drhodpsi, q_as_x
        use geometry, only: gradpar, dbdzed, bmag
        use geometry, only: cvdrift, cvdrift0
        use geometry, only: dIdrho, dgradpardrho, dBdrho, d2Bdrdth
        use geometry, only: dcvdriftdrho, dcvdrift0drho
        use parameters_physics, only: radial_variation

        use parameters_physics, only: ydriftknob
        use store_arrays_useful, only: parnlinit

        implicit none

        if (.not. allocated(par_nl_fac)) allocate (par_nl_fac(-nzgrid:nzgrid, nspec))
        ! this is the factor multiplying -dphi/dz * dg/dvpa in the parallel nonlinearity
        par_nl_fac = 0.5 * rhostar * spread(spec%stm_psi0 * spec%zt_psi0, 1, nztot) * spread(gradpar, 2, nspec)

        if (.not. allocated(par_nl_curv)) allocate (par_nl_curv(-nzgrid:nzgrid, nspec))
        ! ydriftknob is here because this term comes from bhat x curvature . grad B
        par_nl_curv = -ydriftknob * rhostar * geo_surf%rgeo * geo_surf%betaprim * drhodpsi &
                        * spread(dbdzed(1, :) * gradpar / bmag(1, :), 2, nspec) / spread(spec%zt_psi0, 1, nztot)

        if (.not. allocated(par_nl_drifty)) allocate (par_nl_drifty(-nzgrid:nzgrid))
        par_nl_drifty = 0.25 * rhostar * cvdrift(1, :)
        if (.not. allocated(par_nl_driftx)) allocate (par_nl_driftx(-nzgrid:nzgrid))
        if (q_as_x) then
            par_nl_driftx = 0.25 * rhostar * cvdrift0(1, :)
        else
            par_nl_driftx = 0.25 * rhostar * cvdrift0(1, :) / geo_surf%shat
        end if

        if (radial_variation) then
            if (.not. allocated(d_par_nl_fac_dr)) allocate (d_par_nl_fac_dr(-nzgrid:nzgrid, nspec))
            ! this is the factor multiplying -dphi/dz * dg/dvpa in the parallel nonlinearity
            d_par_nl_fac_dr = 0.5 * rhostar * spread(spec%stm_psi0 * spec%zt_psi0, 1, nztot) * spread(dgradpardrho, 2, nspec)

            if (.not. allocated(d_par_nl_curv_dr)) allocate (d_par_nl_curv_dr(-nzgrid:nzgrid, nspec))
            ! ydriftknob is here because this term comes from bhat x curvature . grad B
            ! handle terms with no zeroes
            d_par_nl_curv_dr = par_nl_curv * (dIdrho / geo_surf%rgeo - drhodpsi * geo_surf%d2psidr2 &
                                            - spread(dBdrho / bmag(1, :) + dgradpardrho / gradpar, 2, nspec))
            ! handle terms with possible zeroes
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

        parnlinit = .true.

    end subroutine init_parallel_nonlinearity

    !****************************************************************************
    !                       ADVANCE PARALLEL NONLINEARITY                       !
    !****************************************************************************
    subroutine advance_parallel_nonlinearity(g, gout, restart_time_step)

        use constants, only: zi
        use mp, only: proc0, min_allreduce, mp_abort
        use mp, only: scope, allprocs, subprocs
        use stella_layouts, only: vmu_lo, xyz_lo
        use stella_layouts, only: iv_idx, imu_idx, is_idx
        use job_manage, only: time_message
        use finite_differences, only: second_order_centered_zed
        use finite_differences, only: third_order_upwind
        use redistribute, only: gather, scatter
        use store_arrays_fields, only: phi, phi_corr_QN, phi_corr_GA
        use stella_transforms, only: transform_ky2y, transform_y2ky
        use stella_transforms, only: transform_kx2x, transform_x2kx
        use stella_time, only: cfl_dt_parallel, cfl_dt_linear, code_dt, code_dt_max
        use parameters_numerical, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
        use z_grid, only: nzgrid, delzed, ntubes
        use extended_zgrid, only: neigen, nsegments, ikxmod
        use extended_zgrid, only: iz_low, iz_up
        use extended_zgrid, only: periodic
        use parameters_physics, only: full_flux_surface, radial_variation
        use grids_kxky, only: akx, aky, rho_clamped
        use parameters_kxky_grid, only: nakx, naky, nx, ny, ikx_max
        use calculations_kxky, only: swap_kxky, swap_kxky_back
        use velocity_grids, only: nvpa, nmu
        use velocity_grids, only: dvpa, vpa, mu
        use gyro_averages, only: gyro_average
        use gk_parallel_streaming, only: stream_sign
        use dist_redistribute, only: xyz2vmu
        use file_utils, only: runtype_option_switch, runtype_multibox
        use extended_zgrid, only: fill_zed_ghost_zones

        use store_arrays_useful, only: time_parallel_nl
         use reset_timestep, only: reset_dt

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
        logical, intent(out) :: restart_time_step

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
        
        ! WARNING this routine will probably break if neigen_max = 0
        ! which happens when be set grid_option = 'range'

        ! alpha-component of magnetic drift (requires ky -> y)
        if (proc0) call time_message(.false., time_parallel_nl(:, 1), ' parallel nonlinearity advance')

        ! Initialize cfl_dt_parallel
        cfl_dt_parallel = 10000000.

        restart_time_step = .false.

        ! overview:
        ! need g and d<phi>/dz in (x,y) space in
        ! order to upwind dg/dvpa
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

        ! get d<phi>/dz in vmu_lo
        ! we will need to transform it to real-space
        ! as its sign is needed for upwinding of dg/dvpa
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            ! construct <phi>
            dphidz = phi
            if (radial_variation) dphidz = dphidz + phi_corr_QN
            call gyro_average(dphidz, ivmu, phi_gyro)
            if (radial_variation) phi_gyro = phi_gyro + phi_corr_GA(:, :, :, :, ivmu)

            do iky = 1, naky
                do it = 1, ntubes
                do ie = 1, neigen(iky)
                    do iseg = 1, nsegments(ie, iky)
                        ! first fill in ghost zones at boundaries in g(z)
                        call fill_zed_ghost_zones(it, iseg, ie, iky, phi_gyro, gleft, gright)
                        ! now get d<phi>/dz
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
                    ! use reality to swap from ky >= 0, all kx to kx >= 0 , all ky
                    call swap_kxky(dphidz(:, :, iz, it), g0k_swap)
                    ! transform in y
                    call transform_ky2y(g0k_swap, g0kxy)
                    ! transform in x
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
                    ! use reality to swap from ky >= 0, all kx to kx >= 0 , all ky
                    call swap_kxky(g0k, g0k_swap)
                    ! transform in y
                    call transform_ky2y(g0k_swap, g0kxy)
                    ! transform in x
                    call transform_kx2x(g0kxy, g0xy(:, :, iz, it, ivmu))
                end do
                end do
            end if
        end do

        ! do not need phi_gyro or dphidz  again so deallocate
        deallocate (phi_gyro, dphidz)
        deallocate (g0k)
        if (allocated(g1xy)) deallocate (g1xy)

        allocate (gxy_vmulocal(nvpa, nmu, xyz_lo%llim_proc:xyz_lo%ulim_alloc))
        allocate (advect_speed(nmu, xyz_lo%llim_proc:xyz_lo%ulim_alloc))

        ! we now have the advection velocity in vpa in (x,y) space
        ! next redistribute it so that (vpa,mu) are local
        if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
        call scatter(xyz2vmu, g0xy, gxy_vmulocal)
        if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
        ! advect_speed does not depend on vpa
        advect_speed = gxy_vmulocal(1, :, :)

        ! transform g from (kx,ky) to (x,y)
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
                do iz = -nzgrid, nzgrid
                call swap_kxky(g(:, :, iz, it, ivmu), g0k_swap)
                ! transform in y
                call transform_ky2y(g0k_swap, g0kxy)
                ! transform in x
                call transform_kx2x(g0kxy, g0xy(:, :, iz, it, ivmu))
                end do
            end do
        end do

        ! redistribute so that (vpa,mu) local
        if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
        call scatter(xyz2vmu, g0xy, gxy_vmulocal)
        if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')

        allocate (dgdv(nvpa))

        ! we now need to form dg/dvpa and obtain product of dg/dvpa with advection speed
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

        ! finished with dgdv and advect_speed
        deallocate (dgdv, advect_speed)

        ! now that we have the full parallel nonlinearity in (x,y)-space
        ! need to redistribute so that (x,y) local for transforms
        if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
        call gather(xyz2vmu, gxy_vmulocal, g0xy)
        if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')

        ! finished with gxy_vmulocal
        deallocate (gxy_vmulocal)

        ! g0xy is parallel nonlinearity term with (x,y) on processor
        ! need to inverse Fourier transform
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

        !> check estimated cfl_dt to see if the time step size needs to be changed
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
    !    else
    !       gout = code_dt*gout
        end if

        if (proc0) call time_message(.false., time_parallel_nl(:, 1), ' parallel nonlinearity advance')

    end subroutine advance_parallel_nonlinearity

    !****************************************************************************
    !                      FINALISE PARALLEL NONLINEARITY                       !
    !****************************************************************************
    subroutine finish_parallel_nonlinearity

        use store_arrays_useful, only: parnlinit

        implicit none

        if (allocated(par_nl_fac)) deallocate (par_nl_fac)
        if (allocated(par_nl_curv)) deallocate (par_nl_curv)
        if (allocated(par_nl_driftx)) deallocate (par_nl_driftx)
        if (allocated(par_nl_drifty)) deallocate (par_nl_drifty)

        parnlinit = .false.

   end subroutine finish_parallel_nonlinearity

end module parallel_nonlinearity