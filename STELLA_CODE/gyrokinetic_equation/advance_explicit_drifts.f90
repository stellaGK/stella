module advance_explicit_drifts

    use debug_flags, only: debug => time_advance_debug
    use store_arrays_useful, only: time_gke

    implicit none

    public :: advance_wdriftx_explicit
    public :: advance_wdrifty_explicit
    public :: advance_wstar_explicit

    private

contains

    !*****************************************************************************
    !                           Advance explicit drifts
    !*****************************************************************************
    subroutine advance_wstar_explicit(phi, gout)

        use mp, only: proc0, mp_abort
        use job_manage, only: time_message
        use store_arrays_fields, only: apar, bpar
        use stella_layouts, only: vmu_lo
        use stella_transforms, only: transform_ky2y
        use z_grid, only: nzgrid, ntubes
        use parameters_kxky_grid, only: naky, naky_all, nakx, ikx_max, ny
        use calculations_kxky, only: swap_kxky
        use parameters_physics, only: full_flux_surface
        use store_arrays_distribution_fn, only: wstar, g_scratch
        use gyro_averages, only: gyro_average

        use calculations_kxky_derivatives, only: get_dgdy, get_dchidy

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

        complex, dimension(:, :, :, :, :), allocatable :: g0, g0y
        complex, dimension(:, :), allocatable :: g0_swap

        integer :: iz, it, ivmu

        !> start timing the time advance due to the driving gradients
        if (proc0) call time_message(.false., time_gke(:, 6), ' wstar advance')

        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        if (full_flux_surface) then
            !> assume only a single flux surface simulated
            it = 1
            allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            allocate (g0_swap(naky_all, ikx_max))

            !> calculate d<phi>/dy in k-space
            !> Here g_scratch is <phi> in k-space that has been pre-calculated and stored
            call get_dgdy(g_scratch, g0)
            
            !> transform d<phi>/dy from ky-space to y-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                do iz = -nzgrid, nzgrid
                call swap_kxky(g0(:, :, iz, it, ivmu), g0_swap)
                call transform_ky2y(g0_swap, g0y(:, :, iz, it, ivmu))
                end do
            end do
            !> multiply d<chi>/dy with omega_* coefficient and add to source (RHS of GK eqn)
            !       call add_wstar_term_ffs (g0y, gout)
            call add_explicit_term_ffs(g0y, wstar, gout)
            deallocate (g0y, g0_swap)
        else
            !> get d<chi>/dy in k-space
            if (debug) write (*, *) 'time_advance::solve_gke::get_dchidy'
            call get_dchidy(phi, apar, bpar, g0)
            !> omega_* stays in ky,kx,z space with ky,kx,z local
            !> multiply d<chi>/dy with omega_* coefficient and add to source (RHS of GK eqn)
            if (debug) write (*, *) 'time_advance::solve_gke::add_wstar_term'
            !       call add_wstar_term (g0, gout)
            call add_explicit_term(g0, wstar(1, :, :), gout)
        end if
        deallocate (g0)

        !> stop timing the time advance due to the driving gradients
        if (proc0) call time_message(.false., time_gke(:, 6), ' wstar advance')

    end subroutine advance_wstar_explicit

    !> advance_wdrifty_explicit subroutine calculates and adds the y-component of the
    !> magnetic drift term to the RHS of the GK equation
    subroutine advance_wdrifty_explicit(g, phi, bpar, gout)

        use mp, only: proc0
        use stella_layouts, only: vmu_lo
        use job_manage, only: time_message
        use stella_transforms, only: transform_ky2y
        use z_grid, only: nzgrid, ntubes
        use parameters_kxky_grid, only: nakx, ikx_max, naky, naky_all, ny
        use calculations_kxky, only: swap_kxky
        use parameters_physics, only: full_flux_surface, include_bpar
        use gyro_averages, only: gyro_average, gyro_average_j1
        use store_arrays_distribution_fn, only: wdrifty_g, wdrifty_phi, wdrifty_bpar
        use store_arrays_distribution_fn, only: g_scratch

        use calculations_kxky_derivatives, only: get_dgdy

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

        integer :: ivmu, iz, it
        complex, dimension(:, :, :, :), allocatable :: dphidy, dbpardy
        complex, dimension(:, :, :, :, :), allocatable :: g0k, g0y
        complex, dimension(:, :), allocatable :: g0k_swap

        !> start the timing of the y component of the magnetic drift advance
        if (proc0) call time_message(.false., time_gke(:, 4), ' dgdy advance')

        allocate (dphidy(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (dbpardy(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (g0k(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wdrifty_explicit::get_dgdy'
        !> calculate dg/dy in (ky,kx) space
        call get_dgdy(g, g0k)
        !> calculate dbpar/dy in (ky,kx) space
        if (include_bpar) call get_dgdy(bpar, dbpardy)

        if (full_flux_surface) then
            !> assume only a single flux surface simulated
            it = 1
            allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            allocate (g0k_swap(naky_all, ikx_max))
            !> transform dg/dy from k-space to y-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                do iz = -nzgrid, nzgrid
                call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
                call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
                end do
            end do

            !> add vM . grad y dg/dy term to equation
            call add_explicit_term_ffs(g0y, wdrifty_g, gout)

            !> > calculate dphi/dy in (ky,kx) space
            !> Here g_scratch is <phi> in k-space that has been pre-calculated and stored
            call get_dgdy(g_scratch, g0k)

            !> transform d<phi>/dy from k-space to y-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                do iz = -nzgrid, nzgrid
                call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
                call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
                end do
            end do

            !> add vM . grad y d<phi>/dy term to equation
            call add_explicit_term_ffs(g0y, wdrifty_phi, gout)

            deallocate (g0y, g0k_swap)
        else
            if (debug) write (*, *) 'time_advance::solve_gke::add_dgdy_term'
            ! add vM . grad y dg/dy term to equation
            call add_explicit_term(g0k, wdrifty_g(1, :, :), gout)

            !> Note that this is here because for FFS te gyro-average is calculated once outside this routine
            !> TODO-GA: can we do something similar for fluxtube to save cpu time?
            !> calculate dphi/dy in (ky,kx) space
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

        !> stop the timing of the y component of the magnetic drift advance
        if (proc0) call time_message(.false., time_gke(:, 4), ' dgdy advance')

    end subroutine advance_wdrifty_explicit

    !> advance_wdriftx_explicit subroutine calculates and adds the x-component of the
    !> magnetic drift term to the RHS of the GK equation
    subroutine advance_wdriftx_explicit(g, phi, bpar, gout)

        use mp, only: proc0
        use stella_layouts, only: vmu_lo
        use job_manage, only: time_message
        use stella_transforms, only: transform_ky2y
        use z_grid, only: nzgrid, ntubes
        use parameters_kxky_grid, only: nakx, ikx_max, naky, naky_all, ny
        use grids_kxky, only: akx
        use calculations_kxky, only: swap_kxky
        use parameters_physics, only: full_flux_surface, include_bpar
        use gyro_averages, only: gyro_average
        use store_arrays_distribution_fn, only: wdriftx_g, wdriftx_phi, wdriftx_bpar
        use store_arrays_distribution_fn, only: g_scratch
        use calculations_kxky_derivatives, only: get_dgdx

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

        integer :: ivmu, iz, it
        complex, dimension(:, :, :, :), allocatable :: dphidx, dbpardx
        complex, dimension(:, :, :, :, :), allocatable :: g0k, g0y
        complex, dimension(:, :), allocatable :: g0k_swap

        !> start the timing of the x component of the magnetic drift advance
        if (proc0) call time_message(.false., time_gke(:, 5), ' dgdx advance')

        !> do not calculate if wdriftx terms are all zero
        if (maxval(abs(akx)) < epsilon(0.)) then
            if (proc0) call time_message(.false., time_gke(:, 5), ' dgdx advance')
            return
        end if

        allocate (dphidx(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (dbpardx(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (g0k(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        if (debug) write (*, *) 'time_advance::solve_gke::get_dgdx'
        !> calculate dg/dx in (ky,kx) space
        call get_dgdx(g, g0k)

        !> calculate dbpar/dx in (ky,kx) space
        if (include_bpar) call get_dgdx(bpar, dbpardx)

        if (full_flux_surface) then
            !> assume a single flux surface is simulated
            it = 1
            allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            allocate (g0k_swap(naky_all, ikx_max))
            !> transform dg/dx from k-space to y-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                do iz = -nzgrid, nzgrid
                call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
                call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
                end do
            end do
            !> add vM . grad x dg/dx term to equation
            call add_explicit_term_ffs(g0y, wdriftx_g, gout)

            !> Here g_scratch is <phi> in k-space that has been pre-calculated and stored
            !> get <dphi/dx> in k-space
            call get_dgdx(g_scratch, g0k)

            !> transform d<phi>/dx from k-space to y-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                do iz = -nzgrid, nzgrid
                call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
                call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
                end do
            end do
            !> add vM . grad x d<phi>/dx term to equation
            call add_explicit_term_ffs(g0y, wdriftx_phi, gout)
            deallocate (g0y, g0k_swap)
        else
            if (debug) write (*, *) 'time_advance::solve_gke::add_dgdx_term'
            !> add vM . grad x dg/dx term to equation
            call add_explicit_term(g0k, wdriftx_g(1, :, :), gout)
            !> calculate dphi/dx in (ky,kx) space
            call get_dgdx(phi, dphidx)
            !> get <dphi/dx> in k-space
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                call gyro_average(dphidx, ivmu, g0k(:, :, :, :, ivmu))
            end do
            !> add vM . grad x d<phi>/dx term to equation
            call add_explicit_term(g0k, wdriftx_phi(1, :, :), gout)
            if (include_bpar) then
                !> get <dbpar/dx> in k-space
                do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                call gyro_average(dbpardx, ivmu, g0k(:, :, :, :, ivmu))
                end do
                !> add vM . grad x ( 4 mu d<bpar>/dx ) term to equation
                call add_explicit_term(g0k, wdriftx_bpar(1, :, :), gout)
            end if
        end if
        deallocate (g0k, dphidx, dbpardx)

        !> stop the timing of the x component of the magnetic drift advance
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

    !> add vM . grad y d<phi>/dy or vM . grad x d<phi>/dx (or equivalents with g) or omega_* * d<phi>/dy term to RHS of GK equation
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
    

end module advance_explicit_drifts