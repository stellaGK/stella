module radial_variation_time_advance

    use debug_flags, only: debug => parallel_streaming_debug 

    implicit none

    public :: init_radial_variation
    public :: advance_radial_variation
    public :: finish_radial_variation
    public :: mb_communicate

    private

contains

    !****************************************************************************
    !                 INITIALISE RADIAL VARIATION TIME ADVANCE                  !
    !****************************************************************************

    subroutine init_radial_variation

        use stella_layouts, only: vmu_lo
        use stella_layouts, only: iv_idx, imu_idx, is_idx
        use stella_time, only: code_dt

        use species, only: spec, pfac
        use z_grid, only: nzgrid

        use parameters_kxky_grid, only: nalpha
        use parameters_physics, only: xdriftknob, ydriftknob, wstarknob

        use geometry, only: drhodpsi, dydalpha, gfac
        use geometry, only: dBdrho, geo_surf, q_as_x
        use geometry, only: dcvdriftdrho, dcvdrift0drho
        use geometry, only: dgbdriftdrho, dgbdrift0drho

        use velocity_grids, only: vperp2, vpa, mu
        use velocity_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac

        use store_arrays_distribution_fn, only: wstarp
        use store_arrays_distribution_fn, only: wdriftx_phi, wdrifty_phi
        use store_arrays_distribution_fn, only: wdriftpx_g, wdriftpy_g
        use store_arrays_distribution_fn, only: wdriftpx_phi, wdriftpy_phi
        use store_arrays_useful, only: radialinit
      
        implicit none

        integer :: is, imu, iv, ivmu
        real :: fac
        real, dimension(:, :), allocatable :: energy

        real, dimension(:, :), allocatable :: wcvdrifty, wgbdrifty
        real, dimension(:, :), allocatable :: wcvdriftx, wgbdriftx

        if (radialinit) return
        radialinit = .true.

        allocate (wcvdrifty(nalpha, -nzgrid:nzgrid))
        allocate (wgbdrifty(nalpha, -nzgrid:nzgrid))
        allocate (wcvdriftx(nalpha, -nzgrid:nzgrid))
        allocate (wgbdriftx(nalpha, -nzgrid:nzgrid))

        if (.not. allocated(wstarp)) &
            allocate (wstarp(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstarp = 0.0
        if (.not. allocated(wdriftpx_phi)) &
            allocate (wdriftpx_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        if (.not. allocated(wdriftpy_phi)) &
            allocate (wdriftpy_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        if (.not. allocated(wdriftpx_g)) &
            allocate (wdriftpx_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        if (.not. allocated(wdriftpy_g)) &
            allocate (wdriftpy_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        allocate (energy(nalpha, -nzgrid:nzgrid))

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            energy = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
            !FLAG DSO - THIS NEEDS TO BE ADDED SOMEDAY!
            !if (include_neoclassical_terms) then
            !   wstarp(:,:,ivmu) = dydalpha*drhodpsi*wstarknob*0.5*code_dt &
            !        * (maxwell_vpa(iv)*maxwell_mu(:,:,imu) &
            !        * (spec(is)%fprim+spec(is)%tprim*(energy-1.5)) &
            !        - dfneo_drho(:,:,ivmu))
            !else
            !recall that fprim = -dn/dr and trpim = -dt/dr

            wstarp(:, :, ivmu) = -wstarknob * 0.5 * code_dt &
                                * dydalpha * drhodpsi * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is) &
                                * (pfac * (spec(is)%d2ndr2 - (spec(is)%fprim)**2 - (spec(is)%tprim)**2 * energy) &
                                    + pfac * (spec(is)%d2Tdr2 - (spec(is)%tprim)**2) * (energy - 1.5) &
                                    - gfac * 2 * spec(is)%tprim * mu(imu) * spread(dBdrho, 1, nalpha) &
                                    + (spec(is)%fprim + spec(is)%tprim * (energy - 1.5)) &
                                    * (pfac * (spec(is)%fprim + spec(is)%tprim * (energy - 1.5)) &
                                        + gfac * 2 * mu(imu) * spread(dBdrho, 1, nalpha) &
                                        + gfac * drhodpsi * geo_surf%d2psidr2))

            !end if

            !wdrift
            fac = -ydriftknob * 0.5 * code_dt * spec(is)%tz_psi0
            ! this is the curvature drift piece of wdrifty with missing factor of vpa
            ! vpa factor is missing to avoid singularity when including
            ! non-Maxwellian corrections to equilibrium
            wcvdrifty = gfac * fac * dcvdriftdrho * vpa(iv)
            ! this is the grad-B drift piece of wdrifty
            wgbdrifty = gfac * fac * dgbdriftdrho * 0.5 * vperp2(:, :, imu)
            wdriftpy_g(:, :, ivmu) = wcvdrifty * vpa(iv) + wgbdrifty

            wdriftpy_phi(:, :, ivmu) = spec(is)%zt * (wgbdrifty + wcvdrifty * vpa(iv)) &
                                        * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is) &
                                        - wdrifty_phi(:, :, ivmu) * (pfac * (spec(is)%fprim + spec(is)%tprim * (energy - 2.5)) &
                                                                    + gfac * 2.*mu(imu) * spread(dBdrho, 1, nalpha))

            if (q_as_x) then
                fac = -xdriftknob * 0.5 * code_dt * spec(is)%tz_psi0
            else
                fac = -xdriftknob * 0.5 * code_dt * spec(is)%tz_psi0 / geo_surf%shat
            end if
            ! this is the curvature drift piece of wdriftx with missing factor of vpa
            ! vpa factor is missing to avoid singularity when including
            ! non-Maxwellian corrections to equilibrium
            wcvdriftx = gfac * fac * dcvdrift0drho * vpa(iv)
            ! this is the grad-B drift piece of wdriftx
            wgbdriftx = gfac * fac * dgbdrift0drho * 0.5 * vperp2(:, :, imu)
            wdriftpx_g(:, :, ivmu) = wgbdriftx + wcvdriftx * vpa(iv)

            wdriftpx_phi(:, :, ivmu) = spec(is)%zt * (wgbdriftx + wcvdriftx * vpa(iv)) &
                                        * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is) &
                                        - wdriftx_phi(:, :, ivmu) * (pfac * (spec(is)%fprim + spec(is)%tprim * (energy - 2.5)) &
                                                                    + gfac * 2.*mu(imu) * spread(dBdrho, 1, nalpha))

    !      !the next piece is everything under the x derivative, as this needs to be
    !      !transformed separately
    !      wdriftpx_phi(:,:,ivmu) = spec(is)%zt*(wgbdriftx + wcvdriftx*vpa(iv))  &
    !           * maxwell_vpa(iv,is)*maxwell_mu(:,:,imu,is)*maxwell_fac(is)

    !      !this is variation in the Maxwellian part of the adiabatic response of phi,
    !      !which needs to be transformed separately before differentiation wrt x
    !      !the gyroaveraging and quasineutrality is already done in fields
    !      adiabatic_phi(:,:,ivmu) = -(pfac*(spec(is)%fprim+spec(is)%tprim*(energy-2.5)) &
    !                                 +gfac*2.*mu(imu)*spread(dBdrho,1,nalpha))

        end do

        deallocate (energy, wcvdriftx, wgbdriftx, wcvdrifty, wgbdrifty)

    end subroutine init_radial_variation

    !******************************************************************************
    !                        TIME ADVANCE FOR RADIAL VARIATION
    !******************************************************************************
    subroutine advance_radial_variation(g, gout)

        use mp, only: mp_abort, proc0
        use job_manage, only: time_message

        use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst

        use stella_layouts, only: vmu_lo
        use stella_layouts, only: iv_idx, imu_idx, is_idx

        use calculations_kxky_derivatives, only: get_dgdy, get_dgdx
        use calculations_kxky_derivatives, only: get_dchidy
        use calculations_kxky, only: multiply_by_rho
        use gyro_averages, only: gyro_average, gyro_average_j1

        use store_arrays_fields, only: phi, apar, bpar
        use store_arrays_fields, only: phi_corr_QN, phi_corr_GA
        
        use z_grid, only: nzgrid, ntubes

        use parameters_kxky_grid, only: nakx, naky
        use parameters_physics, only: full_flux_surface, fphi
        use parameters_physics, only: include_parallel_streaming, include_mirror

        use store_arrays_distribution_fn, only: wdriftx_phi, wdrifty_phi
        use store_arrays_distribution_fn, only: wdriftpx_g, wdriftpy_g
        use store_arrays_distribution_fn, only: wdriftpx_phi, wdriftpy_phi 
        use store_arrays_distribution_fn, only: wstar, wstarp
        use store_arrays_useful, only: time_gke

        use mirror_terms, only: add_mirror_radial_variation
        use gk_flow_shear, only: prl_shear, prl_shear_p
        use gk_parallel_streaming, only: add_parallel_streaming_radial_variation

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

        integer :: ia, ivmu, iv, imu, is, iz, it

        complex, dimension(:, :), allocatable :: g0k, g1k, g0a
        complex, dimension(:, :, :, :, :), allocatable :: g_corr

        allocate (g0k(naky, nakx))

        allocate (g1k(naky, nakx))
        allocate (g0a(naky, nakx))

        if (debug) write (*, *) 'time_advance::solve_gke::advance_radial_variation'

        if (proc0) call time_message(.false., time_gke(:, 10), ' radial variation advance')

        if (include_mirror .or. include_parallel_streaming) then
            allocate (g_corr(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            g_corr = 0.
        end if

        !grab the mirror and parallel streaming corrections here to save on FFTs
        if (include_mirror) then
            call add_mirror_radial_variation(g, g_corr)
        end if
        if (include_parallel_streaming) then
            call add_parallel_streaming_radial_variation(g, g_corr, gout)
        end if

        if (full_flux_surface) then
            ! FLAG -- ADD SOMETHING HERE
            call mp_abort('wstarp term not yet setup for full_flux_surface = .true. aborting.')
        end if

        ia = 1
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
                do iz = -nzgrid, nzgrid
                g0k = 0.

                !wstar variation
                call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0a)
                g0k = g0k + g0a * wstarp(ia, iz, ivmu)

                !radial variation in ExB nonlinearity is handled in advance_ExB_nonlinearity

                !wdrift(x/y) - g

                call get_dgdx(g(:, :, iz, it, ivmu), g0a)
                g0k = g0k + g0a * wdriftpx_g(ia, iz, ivmu)

                call get_dgdy(g(:, :, iz, it, ivmu), g0a)
                g0k = g0k + g0a * wdriftpy_g(ia, iz, ivmu)

                !wdrift - phi
                call get_dgdx(phi(:, :, iz, it), g1k)
                !wdriftx variation
                call gyro_average(g1k, iz, ivmu, g0a)
                g0k = g0k + g0a * wdriftpx_phi(ia, iz, ivmu)

                call get_dgdy(phi(:, :, iz, it), g1k)
                !wdrifty variation
                call gyro_average(g1k, iz, ivmu, g0a)
                g0k = g0k + g0a * wdriftpy_phi(ia, iz, ivmu)

                !prl_shear variation
                g0k = g0k + g0a * prl_shear_p(ia, iz, ivmu)

                !mirror term and/or parallel streaming
                if (include_mirror .or. include_parallel_streaming) then
                    g0k = g0k + g_corr(:, :, iz, it, ivmu)
                end if

                !inverse and forward transforms
                call multiply_by_rho(g0k)

                !quasineutrality/gyroaveraging
                call gyro_average(phi_corr_QN(:, :, iz, it), iz, ivmu, g0a)
                g0a = fphi * (g0a + phi_corr_GA(:, :, iz, it, ivmu))

                !wstar - gyroaverage/quasineutrality variation
                call get_dgdy(g0a, g1k)
                g0k = g0k + g1k * wstar(ia, iz, ivmu)

                !wdrifty gyroaverage/quasineutrality variation
                g0k = g0k + g1k * wdrifty_phi(ia, iz, ivmu)

                !prl_shear gyroaverage/quasineutrality variation
                g0k = g0k + g1k * prl_shear(ia, iz, ivmu)

                !wdriftx gyroaverage/quasineutrality variation
                call get_dgdx(g0a, g1k)
                g0k = g0k + g1k * wdriftx_phi(ia, iz, ivmu)

                gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + g0k
                end do
            end do
        end do

        deallocate (g0k, g1k, g0a)
        if (allocated(g_corr)) deallocate (g_corr)

        if (proc0) call time_message(.false., time_gke(:, 10), ' radial variation advance')

    end subroutine advance_radial_variation

    !******************************************************************************
    !                           MULTIBOX COMMUNICATION SUBROUTINE
    !******************************************************************************

    subroutine mb_communicate(g_in)

        use mp, only: job
        use stella_layouts, only: vmu_lo
        use z_grid, only: nzgrid
        use multibox, only: multibox_communicate, apply_radial_boundary_conditions
        use parameters_multibox, only: use_dirichlet_bc
        use fields, only: fields_updated, advance_fields
        use store_arrays_fields, only: phi, apar, bpar
        use file_utils, only: runtype_option_switch, runtype_multibox

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g_in

        if (runtype_option_switch == runtype_multibox) then
            if (job /= 1) then
                call advance_fields(g_in, phi, apar, bpar, dist='g')
            end if

            call multibox_communicate(g_in)

            if (job == 1) then
                fields_updated = .false.
                call advance_fields(g_in, phi, apar, bpar, dist='g')
            end if
        else if (use_dirichlet_BC) then
            call apply_radial_boundary_conditions(g_in)
            fields_updated = .false.
            call advance_fields(g_in, phi, apar, bpar, dist='g')
        end if

    end subroutine mb_communicate

    subroutine finish_radial_variation

        use store_arrays_useful, only: radialinit
        use store_arrays_distribution_fn, only: wdriftpx_g, wdriftpy_g
        use store_arrays_distribution_fn, only: wdriftpx_phi, wdriftpy_phi
        use store_arrays_distribution_fn, only: wstarp

        implicit none

        if (allocated(wstarp)) deallocate (wstarp)
        if (allocated(wdriftpx_g)) deallocate (wdriftpx_g)
        if (allocated(wdriftpy_g)) deallocate (wdriftpy_g)
        if (allocated(wdriftpx_phi)) deallocate (wdriftpx_phi)
        if (allocated(wdriftpy_phi)) deallocate (wdriftpy_phi)

        radialinit = .false.

    end subroutine finish_radial_variation


end module radial_variation_time_advance