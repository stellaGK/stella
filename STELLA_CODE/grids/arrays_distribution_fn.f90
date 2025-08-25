module arrays_distribution_fn

    use debug_flags, only: debug => dist_fn_debug
  
    implicit none

    public :: init_array_gxyz
    public :: init_arrays_distribution_fn
    public :: finish_arrays_distribution_fn

    private

    logical :: dist_fn_initialised = .false.
    logical :: gxyz_initialised = .false.

contains


    !****************************************************************************
    !                   INITIALISE DISTRIBUTION FUNCTION ARRAYS                 !
    !****************************************************************************
    !> This subroutine initialises the distribution function arrays used in the
    !> STELLA code. It allocates the arrays and sets them to zero.
    subroutine init_arrays_distribution_fn

        use mp, only: proc0
        use stella_layouts, only: kxkyz_lo, vmu_lo, kymus_lo
        use grids_z, only: nzgrid, ntubes
        use parameters_kxky_grid, only: naky, nakx
        use grids_velocity, only: nvpa, nmu
        use store_arrays_distribution_fn, only: gnew, gold, g_scratch
        use store_arrays_distribution_fn, only: gvmu, g_kymus
        use parameters_numerical, only: split_parallel_dynamics
        
        implicit none

        if (dist_fn_initialised) return
        dist_fn_initialised = .true.

        debug = debug .and. proc0

        if (debug) write (*, *) 'dist_fn::init_arrays_distribution_fn::allocate_arrays'
        call allocate_arrays

    contains

        subroutine allocate_arrays

            implicit none

            if (.not. allocated(gnew)) &
                allocate (gnew(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            gnew = 0.
            if (.not. allocated(gold)) &
                allocate (gold(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            gold = 0.
            if (.not. allocated(g_scratch)) &
                allocate (g_scratch(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            g_scratch = 0.
            if (.not. allocated(gvmu)) &
                allocate (gvmu(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
            gvmu = 0.
            if (.not. allocated(g_kymus)) then
                if (.not. split_parallel_dynamics) then
                    allocate (g_kymus(nakx, -nzgrid:nzgrid, ntubes, nvpa, kymus_lo%llim_proc:kymus_lo%ulim_alloc))
                else
                    allocate (g_kymus(1, 1, 1, 1, 1))
                end if
                g_kymus = 0.
            end if

        end subroutine allocate_arrays

    end subroutine init_arrays_distribution_fn

    subroutine init_array_gxyz(restarted)

        use store_arrays_distribution_fn, only: gvmu, gold, gnew
        use redistribute, only: gather, scatter
        use calculations_redistribute, only: kxkyz2vmu
        use parameters_physics, only: radial_variation
        use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
        use calculations_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
        use parameters_kxky_grid, only: nalpha, nakx, naky
        use calculations_kxky, only: multiply_by_rho
        use grids_velocity, only: mu, vpa, vperp2
        use grids_z, only: nzgrid, ntubes
        use grids_species, only: spec, pfac
        use geometry, only: dBdrho, gfac

        implicit none

        real :: corr
        integer :: ivmu, is, imu, iv, it, iz, ia
        real, dimension(:, :), allocatable :: energy
        complex, dimension(:, :), allocatable :: g0k
        logical, intent(in) :: restarted

        if (gxyz_initialised) return
        gxyz_initialised = .false.

        ! get version of g that has ky,kx,z local
        call gather(kxkyz2vmu, gvmu, gnew)

        ia = 1

        !calculate radial corrections to F0 for use in Krook operator, as well as g1 from initialisation
        if (radial_variation) then
            !init_g uses maxwellians, so account for variation in temperature, density, and B

            allocate (energy(nalpha, -nzgrid:nzgrid))
            allocate (g0k(naky, nakx))

            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                is = is_idx(vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)
                iv = iv_idx(vmu_lo, ivmu)
                energy = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
                do it = 1, ntubes
                do iz = -nzgrid, nzgrid

                    corr = -(pfac * (spec(is)%fprim + spec(is)%tprim * (energy(ia, iz) - 1.5)) &
                            + 2 * gfac * mu(imu) * dBdrho(iz))

                    if (.not. restarted) then
                        g0k = corr * gnew(:, :, iz, it, ivmu)
                        call multiply_by_rho(g0k)
                        gnew(:, :, iz, it, ivmu) = gnew(:, :, iz, it, ivmu) + g0k
                    end if
                end do
                end do
            end do
            deallocate (energy, g0k)

            if (.not. restarted) call scatter(kxkyz2vmu, gnew, gvmu)
        end if

        gold = gnew

    end subroutine init_array_gxyz

    !-------------------------------------------------------------------------
    !                 INITIALISE GXYZ DISTRIBUTION FUNCTION ARRAYS           !
    !-------------------------------------------------------------------------
    !> This subroutine initialises the gxyz arrays, which are used to store the
    !> distribution function in the kxkyz layout. It gathers the gvmu array
    !> and calculates the radial corrections to F0 for use in the Krook operator.
    !> It also initialises the gold array with the values of gnew.
    !> If the radial variation is not enabled, it simply copies gnew to gold.
    !> The gxyz arrays are used in the Krook operator and projection method.
    subroutine init_gxyz(restarted)

        use store_arrays_distribution_fn, only: gvmu, gold, gnew
        use redistribute, only: gather, scatter
        use calculations_redistribute, only: kxkyz2vmu
        use parameters_physics, only: radial_variation
        use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
        use calculations_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
        use parameters_kxky_grid, only: nalpha, nakx, naky
        use calculations_kxky, only: multiply_by_rho
        use grids_velocity, only: mu, vpa, vperp2
        use grids_z, only: nzgrid, ntubes
        use grids_species, only: spec, pfac
        use geometry, only: dBdrho, gfac

        implicit none

        real :: corr
        integer :: ivmu, is, imu, iv, it, iz, ia
        real, dimension(:, :), allocatable :: energy
        complex, dimension(:, :), allocatable :: g0k
        logical, intent(in) :: restarted

        if (gxyz_initialised) return
        gxyz_initialised = .false.

        ! get version of g that has ky,kx,z local
        call gather(kxkyz2vmu, gvmu, gnew)

        ia = 1

        ! Calculate radial corrections to F0 for use in Krook operator, as well as g1 from initialisation
        if (radial_variation) then
            !init_g uses maxwellians, so account for variation in temperature, density, and B

            allocate (energy(nalpha, -nzgrid:nzgrid))
            allocate (g0k(naky, nakx))

            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                is = is_idx(vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)
                iv = iv_idx(vmu_lo, ivmu)
                energy = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
                do it = 1, ntubes
                do iz = -nzgrid, nzgrid

                    corr = -(pfac * (spec(is)%fprim + spec(is)%tprim * (energy(ia, iz) - 1.5)) &
                            + 2 * gfac * mu(imu) * dBdrho(iz))

                    if (.not. restarted) then
                        g0k = corr * gnew(:, :, iz, it, ivmu)
                        call multiply_by_rho(g0k)
                        gnew(:, :, iz, it, ivmu) = gnew(:, :, iz, it, ivmu) + g0k
                    end if
                end do
                end do
            end do
            deallocate (energy, g0k)

            if (.not. restarted) call scatter(kxkyz2vmu, gnew, gvmu)
        end if

        gold = gnew

    end subroutine init_gxyz

    !****************************************************************************
    !                   FINALISE DISTRIBUTION FUNCTION ARRAYS                 !
    !****************************************************************************
    !> This subroutine finalises the distribution function arrays, deallocating
    !> them and resetting the initialised flags.
    subroutine finish_arrays_distribution_fn

        use store_arrays_distribution_fn, only: gnew, gold, g_scratch, gvmu, g_kymus

        implicit none

        if (allocated(gnew)) deallocate (gnew)
        if (allocated(gold)) deallocate (gold)
        if (allocated(g_scratch)) deallocate (g_scratch)
        if (allocated(gvmu)) deallocate (gvmu)
        if (allocated(g_kymus)) deallocate (g_kymus)

        dist_fn_initialised = .false.
        gxyz_initialised = .false.

    end subroutine finish_arrays_distribution_fn

end module arrays_distribution_fn