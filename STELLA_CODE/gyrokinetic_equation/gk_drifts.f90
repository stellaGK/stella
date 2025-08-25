module gk_drifts

    use debug_flags, only: debug => time_advance_debug

    implicit none
    
    public :: init_wdrift
    public :: init_wstar
    public :: finish_wdrift
    public :: finish_wstar

    public :: advance_wdriftx_explicit
    public :: advance_wdrifty_explicit
    public :: advance_wstar_explicit

    private

contains

    !*****************************************************************************
    !                           Initialise explicit drifts
    !*****************************************************************************
    subroutine init_wdrift

        use mp, only: mp_abort
        use store_arrays_distribution_fn, only: wdriftx_g, wdrifty_g
        use store_arrays_distribution_fn, only: wdriftx_phi, wdrifty_phi
        use store_arrays_distribution_fn, only: wdriftx_bpar, wdrifty_bpar
        use stella_layouts, only: vmu_lo
        use stella_layouts, only: iv_idx, imu_idx, is_idx
        use stella_time, only: code_dt
        use species, only: spec
        use z_grid, only: nzgrid
        use parameters_kxky_grid, only: nalpha
        use geometry, only: cvdrift, gbdrift
        use geometry, only: cvdrift0, gbdrift0
        use geometry, only: gds23, gds24
        use geometry, only: geo_surf, q_as_x
        use geometry, only: dxdpsi, drhodpsi, dydalpha
        use velocity_grids, only: vpa, vperp2, mu
        use velocity_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
        use neoclassical_terms, only: include_neoclassical_terms
        use neoclassical_terms, only: dphineo_dzed, dphineo_drho, dphineo_dalpha
        use neoclassical_terms, only: dfneo_dvpa, dfneo_dzed, dfneo_dalpha
        use parameters_numerical, only: maxwellian_normalization

        use parameters_physics, only: xdriftknob, ydriftknob
        use store_arrays_useful, only: wdriftinit
        implicit none

        integer :: ivmu, iv, imu, is
        real :: fac
        real, dimension(:, :), allocatable :: wcvdrifty, wgbdrifty
        real, dimension(:, :), allocatable :: wcvdriftx, wgbdriftx

        if (wdriftinit) return
        wdriftinit = .true.

        ! Allocate wdriftx_phi, the factor multiplying dphi/dx in the magnetic drift term
        if (.not. allocated(wdriftx_phi)) then
            allocate (wdriftx_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            wdriftx_phi = 0.0
        end if
        ! Allocate wdrifty_phi, the factor multiplying dphi/dy in the magnetic drift term
        if (.not. allocated(wdrifty_phi)) then
            allocate (wdrifty_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            wdrifty_phi = 0.0
        end if
        ! Allocate wdriftx_bpar, the factor multiplying dbpar/dx in the magnetic drift term
        if (.not. allocated(wdriftx_bpar)) then
            allocate (wdriftx_bpar(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            wdriftx_bpar = 0.0
        end if
        ! Allocate wdrifty_bpar, the factor multiplying dbpar/dy in the magnetic drift term
        if (.not. allocated(wdrifty_bpar)) then
            allocate (wdrifty_bpar(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            wdrifty_bpar = 0.0
        end if
        ! Allocate wdriftx_g, the factor multiplying dg/dx in the magnetic drift term
        if (.not. allocated(wdriftx_g)) then
            allocate (wdriftx_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            wdriftx_g = 0.0
        end if
        ! Allocate wdrifty_g, the factor multiplying dg/dy in the magnetic drift term
        if (.not. allocated(wdrifty_g)) then
            allocate (wdrifty_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            wdrifty_g = 0.0
        end if

        allocate (wcvdrifty(nalpha, -nzgrid:nzgrid))
        allocate (wgbdrifty(nalpha, -nzgrid:nzgrid))
        allocate (wcvdriftx(nalpha, -nzgrid:nzgrid))
        allocate (wgbdriftx(nalpha, -nzgrid:nzgrid))

        ! FLAG -- need to deal with shat=0 case.  ideally move away from q as x-coordinate
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            fac = -ydriftknob * 0.5 * code_dt * spec(is)%tz_psi0
            ! This is the curvature drift piece of wdrifty with missing factor of vpa
            ! vpa factor is missing to avoid singularity when including
            ! non-Maxwellian corrections to equilibrium
            wcvdrifty = fac * cvdrift * vpa(iv)
            ! This is the grad-B drift piece of wdrifty
            wgbdrifty = fac * gbdrift * 0.5 * vperp2(:, :, imu)
            wdrifty_g(:, :, ivmu) = wcvdrifty * vpa(iv) + wgbdrifty
            ! if including neoclassical correction to equilibrium Maxwellian,
            ! then add in v_E^{nc} . grad y dg/dy coefficient here
            if (include_neoclassical_terms) then
                wdrifty_g(:, :, ivmu) = wdrifty_g(:, :, ivmu) + code_dt * 0.5 * (gds23 * dphineo_dzed + drhodpsi * dydalpha * dphineo_drho)
            end if

            wdrifty_phi(:, :, ivmu) = spec(is)%zt * (wgbdrifty + wcvdrifty * vpa(iv))

            ! if maxwwellian_normalization = .true., evolved distribution function is normalised by a Maxwellian
            ! otherwise, it is not; a Maxwellian weighting factor must thus be included
            if (.not. maxwellian_normalization) then
                wdrifty_phi(:, :, ivmu) = wdrifty_phi(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
            end if
            ! assign wdrifty_bpar, neoclassical terms not supported
            wdrifty_bpar(:,:,ivmu) = 4.0 * mu(imu) * wdrifty_phi(:, :, ivmu) * spec(is)%tz
            ! if including neoclassical corrections to equilibrium,
            ! add in -(Ze/m) * v_curv/vpa . grad y d<phi>/dy * dF^{nc}/dvpa term
            ! and v_E . grad z dF^{nc}/dz (here get the dphi/dy part of v_E)
            if (include_neoclassical_terms) then
                ! NB: the below neoclassical correction needs to be divided by an equilibrium Maxwellian
                ! if maxwellian_normalization = .true.
                if (maxwellian_normalization) then
                call mp_abort("include_neoclassical_terms=T not currently supported for maxwellian_normalization=T.  aborting")
                end if
                wdrifty_phi(:, :, ivmu) = wdrifty_phi(:, :, ivmu) &
                                        - 0.5 * spec(is)%zt * dfneo_dvpa(:, :, ivmu) * wcvdrifty &
                                        - code_dt * 0.5 * dfneo_dzed(:, :, ivmu) * gds23
            end if

            if (q_as_x) then
                fac = -xdriftknob * 0.5 * code_dt * spec(is)%tz_psi0
            else
                fac = -xdriftknob * 0.5 * code_dt * spec(is)%tz_psi0 / geo_surf%shat
            end if
            ! This is the curvature drift piece of wdriftx with missing factor of vpa
            ! vpa factor is missing to avoid singularity when including
            ! non-Maxwellian corrections to equilibrium
            wcvdriftx = fac * cvdrift0 * vpa(iv)
            ! This is the grad-B drift piece of wdriftx
            wgbdriftx = fac * gbdrift0 * 0.5 * vperp2(:, :, imu)
            wdriftx_g(:, :, ivmu) = wcvdriftx * vpa(iv) + wgbdriftx
            ! if including neoclassical correction to equilibrium Maxwellian,
            ! then add in v_E^{nc} . grad x dg/dx coefficient here
            if (include_neoclassical_terms) then
                wdriftx_g(:, :, ivmu) = wdriftx_g(:, :, ivmu) + code_dt * 0.5 * (gds24 * dphineo_dzed - dxdpsi * dphineo_dalpha)
            end if
            wdriftx_phi(:, :, ivmu) = spec(is)%zt * (wgbdriftx + wcvdriftx * vpa(iv))
            ! if maxwellian_normalizatiion = .true., evolved distribution function is normalised by a Maxwellian
            ! otherwise, it is not; a Maxwellian weighting factor must thus be included
            if (.not. maxwellian_normalization) then
                wdriftx_phi(:, :, ivmu) = wdriftx_phi(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
            end if
            ! assign wdriftx_bpar, neoclassical terms not supported
            wdriftx_bpar(:,:,ivmu) = 4.0 * mu(imu) * wdriftx_phi(:, :, ivmu) * spec(is)%tz
            ! if including neoclassical corrections to equilibrium,
            ! add in (Ze/m) * v_curv/vpa . grad x d<phi>/dx * dF^{nc}/dvpa term
            ! and v_E . grad z dF^{nc}/dz (here get the dphi/dx part of v_E)
            ! and v_E . grad alpha dF^{nc}/dalpha (dphi/dx part of v_E)
            if (include_neoclassical_terms) then
                ! NB: the below neoclassical correction needs to be divided by an equilibrium Maxwellian
                ! if running with maxwellian_normalzation = .true.
                if (maxwellian_normalization) then
                call mp_abort("include_neoclassical_terms=T not currently supported for maxwellian_normalization=T.  aborting")
                end if
                wdriftx_phi(:, :, ivmu) = wdriftx_phi(:, :, ivmu) &
                                        - 0.5 * spec(is)%zt * dfneo_dvpa(:, :, ivmu) * wcvdriftx &
                                        + code_dt * 0.5 * (dfneo_dalpha(:, :, ivmu) * dxdpsi - dfneo_dzed(:, :, ivmu) * gds24)
            end if

        end do

        deallocate (wcvdriftx, wgbdriftx, wcvdrifty, wgbdrifty)

    end subroutine init_wdrift

     subroutine init_wstar

        use mp, only: mp_abort
        use stella_layouts, only: vmu_lo
        use stella_layouts, only: iv_idx, imu_idx, is_idx
        use stella_time, only: code_dt
        use species, only: spec
        use z_grid, only: nzgrid
        use parameters_kxky_grid, only: nalpha
        use geometry, only: dydalpha, drhodpsi, clebsch_factor
        use velocity_grids, only: vperp2, vpa
        use velocity_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
        use store_arrays_distribution_fn, only: wstar
        use neoclassical_terms, only: include_neoclassical_terms
        use neoclassical_terms, only: dfneo_drho
        use parameters_numerical, only: maxwellian_normalization

        use parameters_physics, only: wstarknob
        use store_arrays_useful, only: wstarinit

        implicit none

        integer :: is, imu, iv, ivmu
        real, dimension(:, :), allocatable :: energy

        if (wstarinit) return
        wstarinit = .true.

        if (.not. allocated(wstar)) &
            allocate (wstar(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstar = 0.0

        allocate (energy(nalpha, -nzgrid:nzgrid))

        ! We have wstar = - (1/C) (a*Bref) (dy/dalpha) (d/dpsi) ... 
        ! The profile gradients are given with respect to r, i.e., <fprim> = -(a/n)(dn/dr)
        ! wstar = - (1/C) (a*Bref) (dy/dalpha) (1/n) (dn/dpsi) ... 
        !       = - (1/C) (a*Bref) (dy/dalpha) (dr/dpsi) (1/n) (dn/dr) ... 
        !       = - (1/C) * (1/a) (dy/dalpha) * (a*Bref) (dr/dpsi) * (a/n) (dn/dr) ... 
        !       = - (1/<clebsch_factor>) * <dydalpha> * <drhodpsi> * <fprim> ...
        ! Note that for psi=psit we have B = sign_torflux ∇ψ x ∇α 
        ! Note that for psi=psip we have B = - ∇ψ x ∇α  
        !     <dydalpha> = (rhor/a)(d(y/rhor)/dalpha) = (1/a)(dy/dalpha) 
        !     <drhodpsi> = drho/dψ̃ = d(r/a)/d(psi/(a^2*Br)) = (a*Bref) * dr/dpsi
        !     1/<clebsch_factor> = -1 or sign_torflux
        ! Note that for psi=q we have B = - (dpsi_p/dq) ∇ψ x ∇α and,
        !     <drhodpsi> = drho/dq = d(r/a)/dq = (1/a) * dr/dq
        !     <dydalpha> = (rhor/a)(d(y/rhor)/dalpha) = (1/a)(dy/dalpha) 
        !     1/<clebsch_factor> = - dq/d(psip/(a^2*Br)) = - (a^2*Bref) (qd/dpsi_p) 
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            energy = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
            if (include_neoclassical_terms) then
                if (maxwellian_normalization) then
                call mp_abort("include_neoclassical_terms = T not yet supported for maxwellian_normalization = T. Aborting.")
                else
                wstar(:, :, ivmu) = - (1/clebsch_factor) * dydalpha * drhodpsi * wstarknob * 0.5 * code_dt &
                                    * (maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is) &
                                        * (spec(is)%fprim + spec(is)%tprim * (energy - 1.5)) &
                                        - dfneo_drho(:, :, ivmu))
                end if
            else
                wstar(:, :, ivmu) = - (1/clebsch_factor) * dydalpha * drhodpsi * wstarknob * 0.5 * code_dt &
                                    * (spec(is)%fprim + spec(is)%tprim * (energy - 1.5))
            end if
            if (.not. maxwellian_normalization) then
                wstar(:, :, ivmu) = wstar(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
            end if
        end do

        deallocate (energy)

    end subroutine init_wstar

    !*****************************************************************************
    !                           Advance explicit drifts
    !*****************************************************************************
    subroutine advance_wstar_explicit(phi, gout)

        use mp, only: proc0, mp_abort
        use job_manage, only: time_message
        use store_arrays_fields, only: apar, bpar
        use stella_layouts, only: vmu_lo
        use calculations_stella_transforms, only: transform_ky2y
        use z_grid, only: nzgrid, ntubes
        use parameters_kxky_grid, only: naky, naky_all, nakx, ikx_max, ny
        use calculations_kxky, only: swap_kxky
        use parameters_physics, only: full_flux_surface
        use store_arrays_distribution_fn, only: wstar, g_scratch
        use calculations_gyro_averages, only: gyro_average

        use calculations_kxky_derivatives, only: get_dgdy, get_dchidy
        use store_arrays_useful, only: time_gke

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

        complex, dimension(:, :, :, :, :), allocatable :: g0, g0y
        complex, dimension(:, :), allocatable :: g0_swap

        integer :: iz, it, ivmu

        ! start timing the time advance due to the driving gradients
        if (proc0) call time_message(.false., time_gke(:, 6), ' wstar advance')

        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        if (full_flux_surface) then
            ! assume only a single flux surface simulated
            it = 1
            allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            allocate (g0_swap(naky_all, ikx_max))

            ! calculate d<phi>/dy in k-space
            ! Here g_scratch is <phi> in k-space that has been pre-calculated and stored
            call get_dgdy(g_scratch, g0)
            
            ! transform d<phi>/dy from ky-space to y-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                do iz = -nzgrid, nzgrid
                call swap_kxky(g0(:, :, iz, it, ivmu), g0_swap)
                call transform_ky2y(g0_swap, g0y(:, :, iz, it, ivmu))
                end do
            end do
            ! multiply d<chi>/dy with omega_* coefficient and add to source (RHS of GK eqn)
            !       call add_wstar_term_ffs (g0y, gout)
            call add_explicit_term_ffs(g0y, wstar, gout)
            deallocate (g0y, g0_swap)
        else
            ! get d<chi>/dy in k-space
            if (debug) write (*, *) 'time_advance::solve_gke::get_dchidy'
            call get_dchidy(phi, apar, bpar, g0)
            ! omega_* stays in ky,kx,z space with ky,kx,z local
            ! multiply d<chi>/dy with omega_* coefficient and add to source (RHS of GK eqn)
            if (debug) write (*, *) 'time_advance::solve_gke::add_wstar_term'
            !       call add_wstar_term (g0, gout)
            call add_explicit_term(g0, wstar(1, :, :), gout)
        end if
        deallocate (g0)

        ! stop timing the time advance due to the driving gradients
        if (proc0) call time_message(.false., time_gke(:, 6), ' wstar advance')

    end subroutine advance_wstar_explicit

    ! advance_wdrifty_explicit subroutine calculates and adds the y-component of the
    ! magnetic drift term to the RHS of the GK equation
    subroutine advance_wdrifty_explicit(g, phi, bpar, gout)

        use mp, only: proc0
        use stella_layouts, only: vmu_lo
        use job_manage, only: time_message
        use calculations_stella_transforms, only: transform_ky2y
        use z_grid, only: nzgrid, ntubes
        use parameters_kxky_grid, only: nakx, ikx_max, naky, naky_all, ny
        use calculations_kxky, only: swap_kxky
        use parameters_physics, only: full_flux_surface, include_bpar
        use calculations_gyro_averages, only: gyro_average, gyro_average_j1
        use store_arrays_distribution_fn, only: wdrifty_g, wdrifty_phi, wdrifty_bpar
        use store_arrays_distribution_fn, only: g_scratch

        use calculations_kxky_derivatives, only: get_dgdy

        use store_arrays_useful, only: time_gke

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

        integer :: ivmu, iz, it
        complex, dimension(:, :, :, :), allocatable :: dphidy, dbpardy
        complex, dimension(:, :, :, :, :), allocatable :: g0k, g0y
        complex, dimension(:, :), allocatable :: g0k_swap

        ! start the timing of the y component of the magnetic drift advance
        if (proc0) call time_message(.false., time_gke(:, 4), ' dgdy advance')

        allocate (dphidy(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (dbpardy(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (g0k(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wdrifty_explicit::get_dgdy'
        ! calculate dg/dy in (ky,kx) space
        call get_dgdy(g, g0k)
        ! calculate dbpar/dy in (ky,kx) space
        if (include_bpar) call get_dgdy(bpar, dbpardy)

        if (full_flux_surface) then
            ! assume only a single flux surface simulated
            it = 1
            allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            allocate (g0k_swap(naky_all, ikx_max))
            ! transform dg/dy from k-space to y-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                do iz = -nzgrid, nzgrid
                call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
                call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
                end do
            end do

            ! add vM . grad y dg/dy term to equation
            call add_explicit_term_ffs(g0y, wdrifty_g, gout)

            ! > calculate dphi/dy in (ky,kx) space
            ! Here g_scratch is <phi> in k-space that has been pre-calculated and stored
            call get_dgdy(g_scratch, g0k)

            ! transform d<phi>/dy from k-space to y-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                do iz = -nzgrid, nzgrid
                call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
                call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
                end do
            end do

            ! add vM . grad y d<phi>/dy term to equation
            call add_explicit_term_ffs(g0y, wdrifty_phi, gout)

            deallocate (g0y, g0k_swap)
        else
            if (debug) write (*, *) 'time_advance::solve_gke::add_dgdy_term'
            ! add vM . grad y dg/dy term to equation
            call add_explicit_term(g0k, wdrifty_g(1, :, :), gout)

            ! Note that this is here because for FFS te gyro-average is calculated once outside this routine
            ! TODO-GA: can we do something similar for fluxtube to save cpu time?
            ! calculate dphi/dy in (ky,kx) space
            call get_dgdy(phi, dphidy)

            ! get <dphi/dy> in k-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                call gyro_average(dphidy, ivmu, g0k(:, :, :, :, ivmu))
            end do

            ! add vM . grad y d<phi>/dy term to equation
            call add_explicit_term(g0k, wdrifty_phi(1, :, :), gout)
            
            if (include_bpar) then
                ! get <dbpar/dy> in k-space
                do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                call gyro_average_j1(dbpardy, ivmu, g0k(:, :, :, :, ivmu))
                end do
                ! add vM . grad y (4 mu d<bpar>/dy) term to equation
                call add_explicit_term(g0k, wdrifty_bpar(1, :, :), gout)            
            end if
        end if
        deallocate (g0k, dphidy, dbpardy)

        ! stop the timing of the y component of the magnetic drift advance
        if (proc0) call time_message(.false., time_gke(:, 4), ' dgdy advance')

    end subroutine advance_wdrifty_explicit

    ! advance_wdriftx_explicit subroutine calculates and adds the x-component of the
    ! magnetic drift term to the RHS of the GK equation
    subroutine advance_wdriftx_explicit(g, phi, bpar, gout)

        use mp, only: proc0
        use stella_layouts, only: vmu_lo
        use job_manage, only: time_message
        use calculations_stella_transforms, only: transform_ky2y
        use z_grid, only: nzgrid, ntubes
        use parameters_kxky_grid, only: nakx, ikx_max, naky, naky_all, ny
        use grids_kxky, only: akx
        use calculations_kxky, only: swap_kxky
        use parameters_physics, only: full_flux_surface, include_bpar
        use calculations_gyro_averages, only: gyro_average
        use store_arrays_distribution_fn, only: wdriftx_g, wdriftx_phi, wdriftx_bpar
        use store_arrays_distribution_fn, only: g_scratch
        use calculations_kxky_derivatives, only: get_dgdx

        use store_arrays_useful, only: time_gke

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

        integer :: ivmu, iz, it
        complex, dimension(:, :, :, :), allocatable :: dphidx, dbpardx
        complex, dimension(:, :, :, :, :), allocatable :: g0k, g0y
        complex, dimension(:, :), allocatable :: g0k_swap

        ! start the timing of the x component of the magnetic drift advance
        if (proc0) call time_message(.false., time_gke(:, 5), ' dgdx advance')

        ! do not calculate if wdriftx terms are all zero
        if (maxval(abs(akx)) < epsilon(0.)) then
            if (proc0) call time_message(.false., time_gke(:, 5), ' dgdx advance')
            return
        end if

        allocate (dphidx(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (dbpardx(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (g0k(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        if (debug) write (*, *) 'time_advance::solve_gke::get_dgdx'
        ! calculate dg/dx in (ky,kx) space
        call get_dgdx(g, g0k)

        ! calculate dbpar/dx in (ky,kx) space
        if (include_bpar) call get_dgdx(bpar, dbpardx)

        if (full_flux_surface) then
            ! assume a single flux surface is simulated
            it = 1
            allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            allocate (g0k_swap(naky_all, ikx_max))
            ! transform dg/dx from k-space to y-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                do iz = -nzgrid, nzgrid
                call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
                call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
                end do
            end do
            ! add vM . grad x dg/dx term to equation
            call add_explicit_term_ffs(g0y, wdriftx_g, gout)

            ! Here g_scratch is <phi> in k-space that has been pre-calculated and stored
            ! get <dphi/dx> in k-space
            call get_dgdx(g_scratch, g0k)

            ! transform d<phi>/dx from k-space to y-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                do iz = -nzgrid, nzgrid
                call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
                call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
                end do
            end do
            ! add vM . grad x d<phi>/dx term to equation
            call add_explicit_term_ffs(g0y, wdriftx_phi, gout)
            deallocate (g0y, g0k_swap)
        else
            if (debug) write (*, *) 'time_advance::solve_gke::add_dgdx_term'
            ! add vM . grad x dg/dx term to equation
            call add_explicit_term(g0k, wdriftx_g(1, :, :), gout)
            ! calculate dphi/dx in (ky,kx) space
            call get_dgdx(phi, dphidx)
            ! get <dphi/dx> in k-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                call gyro_average(dphidx, ivmu, g0k(:, :, :, :, ivmu))
            end do
            ! add vM . grad x d<phi>/dx term to equation
            call add_explicit_term(g0k, wdriftx_phi(1, :, :), gout)
            if (include_bpar) then
                ! get <dbpar/dx> in k-space
                do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                call gyro_average(dbpardx, ivmu, g0k(:, :, :, :, ivmu))
                end do
                ! add vM . grad x ( 4 mu d<bpar>/dx ) term to equation
                call add_explicit_term(g0k, wdriftx_bpar(1, :, :), gout)
            end if
        end if
        deallocate (g0k, dphidx, dbpardx)

        ! stop the timing of the x component of the magnetic drift advance
        if (proc0) call time_message(.false., time_gke(:, 5), ' dgdx advance')
        
    end subroutine advance_wdriftx_explicit

    subroutine add_explicit_term(g, pre_factor, src)

        use stella_layouts, only: vmu_lo
        use z_grid, only: nzgrid, ntubes
        use parameters_kxky_grid, only: naky, nakx

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
        real, dimension(-nzgrid:, vmu_lo%llim_proc:), intent(in) :: pre_factor
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: src

        integer :: ivmu
        integer :: iky, ikx, iz, it

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
                do iz = -nzgrid, nzgrid
                do ikx = 1, nakx
                    do iky = 1, naky
                        src(iky, ikx, iz, it, ivmu) = src(iky, ikx, iz, it, ivmu) + pre_factor(iz, ivmu) * g(iky, ikx, iz, it, ivmu)
                    end do
                end do
                end do
            end do
        end do

    end subroutine add_explicit_term

    ! add vM . grad y d<phi>/dy or vM . grad x d<phi>/dx (or equivalents with g) or omega_* * d<phi>/dy term to RHS of GK equation
    subroutine add_explicit_term_ffs(g, pre_factor, src)

        use stella_layouts, only: vmu_lo
        use z_grid, only: nzgrid, ntubes
        use parameters_kxky_grid, only: ikx_max, nalpha

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
        real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(in) :: pre_factor
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: src

        integer :: ivmu
        integer :: ia, ikx, iz, it

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
                do iz = -nzgrid, nzgrid
                do ikx = 1, ikx_max
                    do ia = 1, nalpha
                        src(ia, ikx, iz, it, ivmu) = src(ia, ikx, iz, it, ivmu) + pre_factor(ia, iz, ivmu) * g(ia, ikx, iz, it, ivmu)
                    end do
                end do
                end do
            end do
        end do

    end subroutine add_explicit_term_ffs

    !*****************************************************************************
    !                           Finalise explicit drifts
    !*****************************************************************************
    subroutine finish_wdrift

        use store_arrays_distribution_fn, only: wdriftx_g, wdrifty_g
        use store_arrays_distribution_fn, only: wdriftx_phi, wdrifty_phi
        use store_arrays_useful, only: wdriftinit

        implicit none

        if (allocated(wdriftx_g)) deallocate (wdriftx_g)
        if (allocated(wdrifty_g)) deallocate (wdrifty_g)
        if (allocated(wdriftx_phi)) deallocate (wdriftx_phi)
        if (allocated(wdrifty_phi)) deallocate (wdrifty_phi)

        wdriftinit = .false.

    end subroutine finish_wdrift

    subroutine finish_wstar 

        use store_arrays_distribution_fn, only: wstar, wstarp
        use store_arrays_useful, only: wstarinit

        implicit none

        if (allocated(wstar)) deallocate (wstar)

        wstarinit = .false.

    end subroutine finish_wstar
    
end module gk_drifts