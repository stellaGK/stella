module arrays_drifts

    implicit none

    public :: init_wdrift
    public :: init_wstar
    public :: finish_wdrift
    public :: finish_wstar

    private 

    logical :: wdriftinit = .false.
    logical :: wstarinit = .false.

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

        !> allocate wdriftx_phi, the factor multiplying dphi/dx in the magnetic drift term
        if (.not. allocated(wdriftx_phi)) then
            allocate (wdriftx_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            wdriftx_phi = 0.0
        end if
        !> allocate wdrifty_phi, the factor multiplying dphi/dy in the magnetic drift term
        if (.not. allocated(wdrifty_phi)) then
            allocate (wdrifty_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            wdrifty_phi = 0.0
        end if
        !> allocate wdriftx_bpar, the factor multiplying dbpar/dx in the magnetic drift term
        if (.not. allocated(wdriftx_bpar)) then
            allocate (wdriftx_bpar(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            wdriftx_bpar = 0.0
        end if
        !> allocate wdrifty_bpar, the factor multiplying dbpar/dy in the magnetic drift term
        if (.not. allocated(wdrifty_bpar)) then
            allocate (wdrifty_bpar(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            wdrifty_bpar = 0.0
        end if
        !> allocate wdriftx_g, the factor multiplying dg/dx in the magnetic drift term
        if (.not. allocated(wdriftx_g)) then
            allocate (wdriftx_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            wdriftx_g = 0.0
        end if
        !> allocate wdrifty_g, the factor multiplying dg/dy in the magnetic drift term
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
            !> this is the curvature drift piece of wdrifty with missing factor of vpa
            !> vpa factor is missing to avoid singularity when including
            !> non-Maxwellian corrections to equilibrium
            wcvdrifty = fac * cvdrift * vpa(iv)
            !> this is the grad-B drift piece of wdrifty
            wgbdrifty = fac * gbdrift * 0.5 * vperp2(:, :, imu)
            wdrifty_g(:, :, ivmu) = wcvdrifty * vpa(iv) + wgbdrifty
            !> if including neoclassical correction to equilibrium Maxwellian,
            !> then add in v_E^{nc} . grad y dg/dy coefficient here
            if (include_neoclassical_terms) then
                wdrifty_g(:, :, ivmu) = wdrifty_g(:, :, ivmu) + code_dt * 0.5 * (gds23 * dphineo_dzed + drhodpsi * dydalpha * dphineo_drho)
            end if

            wdrifty_phi(:, :, ivmu) = spec(is)%zt * (wgbdrifty + wcvdrifty * vpa(iv))

            !> if maxwwellian_normalization = .true., evolved distribution function is normalised by a Maxwellian
            !> otherwise, it is not; a Maxwellian weighting factor must thus be included
            if (.not. maxwellian_normalization) then
                wdrifty_phi(:, :, ivmu) = wdrifty_phi(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
            end if
            !> assign wdrifty_bpar, neoclassical terms not supported
            wdrifty_bpar(:,:,ivmu) = 4.0 * mu(imu) * wdrifty_phi(:, :, ivmu) * spec(is)%tz
            !> if including neoclassical corrections to equilibrium,
            !> add in -(Ze/m) * v_curv/vpa . grad y d<phi>/dy * dF^{nc}/dvpa term
            !> and v_E . grad z dF^{nc}/dz (here get the dphi/dy part of v_E)
            if (include_neoclassical_terms) then
                !> NB: the below neoclassical correction needs to be divided by an equilibrium Maxwellian
                !> if maxwellian_normalization = .true.
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
            !> this is the curvature drift piece of wdriftx with missing factor of vpa
            !> vpa factor is missing to avoid singularity when including
            !> non-Maxwellian corrections to equilibrium
            wcvdriftx = fac * cvdrift0 * vpa(iv)
            !> this is the grad-B drift piece of wdriftx
            wgbdriftx = fac * gbdrift0 * 0.5 * vperp2(:, :, imu)
            wdriftx_g(:, :, ivmu) = wcvdriftx * vpa(iv) + wgbdriftx
            !> if including neoclassical correction to equilibrium Maxwellian,
            !> then add in v_E^{nc} . grad x dg/dx coefficient here
            if (include_neoclassical_terms) then
                wdriftx_g(:, :, ivmu) = wdriftx_g(:, :, ivmu) + code_dt * 0.5 * (gds24 * dphineo_dzed - dxdpsi * dphineo_dalpha)
            end if
            wdriftx_phi(:, :, ivmu) = spec(is)%zt * (wgbdriftx + wcvdriftx * vpa(iv))
            !> if maxwellian_normalizatiion = .true., evolved distribution function is normalised by a Maxwellian
            !> otherwise, it is not; a Maxwellian weighting factor must thus be included
            if (.not. maxwellian_normalization) then
                wdriftx_phi(:, :, ivmu) = wdriftx_phi(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
            end if
            !> assign wdriftx_bpar, neoclassical terms not supported
            wdriftx_bpar(:,:,ivmu) = 4.0 * mu(imu) * wdriftx_phi(:, :, ivmu) * spec(is)%tz
            !> if including neoclassical corrections to equilibrium,
            !> add in (Ze/m) * v_curv/vpa . grad x d<phi>/dx * dF^{nc}/dvpa term
            !> and v_E . grad z dF^{nc}/dz (here get the dphi/dx part of v_E)
            !> and v_E . grad alpha dF^{nc}/dalpha (dphi/dx part of v_E)
            if (include_neoclassical_terms) then
                !> NB: the below neoclassical correction needs to be divided by an equilibrium Maxwellian
                !> if running with maxwellian_normalzation = .true.
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
    !                           Finalise explicit drifts
    !*****************************************************************************
    subroutine finish_wdrift

        use store_arrays_distribution_fn, only: wdriftx_g, wdrifty_g
        use store_arrays_distribution_fn, only: wdriftx_phi, wdrifty_phi
        use store_arrays_distribution_fn, only: wdriftpx_g, wdriftpy_g
        use store_arrays_distribution_fn, only: wdriftpx_phi, wdriftpy_phi
        use store_arrays_useful, only: wdriftinit
    !use store_arrays_distribution_fn, only: adiabatic_phi

        implicit none

        if (allocated(wdriftx_g)) deallocate (wdriftx_g)
        if (allocated(wdrifty_g)) deallocate (wdrifty_g)
        if (allocated(wdriftx_phi)) deallocate (wdriftx_phi)
        if (allocated(wdrifty_phi)) deallocate (wdrifty_phi)
        if (allocated(wdriftpx_g)) deallocate (wdriftpx_g)
        if (allocated(wdriftpy_g)) deallocate (wdriftpy_g)
        if (allocated(wdriftpx_phi)) deallocate (wdriftpx_phi)
        if (allocated(wdriftpy_phi)) deallocate (wdriftpy_phi)
    !   if (allocated(adiabatic_phi)) deallocate (adiabatic_phi)

        wdriftinit = .false.

    end subroutine finish_wdrift

    subroutine finish_wstar 

        use store_arrays_distribution_fn, only: wstar, wstarp
        use store_arrays_useful, only: wstarinit

        implicit none

        
        if (allocated(wstar)) deallocate (wstar)
        if (allocated(wstarp)) deallocate (wstarp)

        wstarinit = .false.

    end subroutine finish_wstar

end module arrays_drifts