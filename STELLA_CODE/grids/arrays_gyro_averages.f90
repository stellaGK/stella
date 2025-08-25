module arrays_gyro_averages

    use common_types, only: coupled_alpha_type
    use debug_flags, only: debug => gyro_averages_debug
    use debug_flags, only: debug_test_gyro_average

    implicit none
    
    public :: init_bessel, finish_bessel

    public :: aj0x, aj0v, aj1x, aj1v    
    public :: j0_B_ffs, j0_ffs

    public :: find_max_required_kalpha_index
    public :: j1_ffs
    public :: j0_const, j0_B_const
    public :: j0max_const

    private

    ! Local variables for gyro averages
    real, dimension(:, :, :, :), allocatable :: aj0x, aj1x
    ! (naky, nakx, nalpha, -nzgrid:nzgrid, -vmu-layout-)
    real, dimension(:, :), allocatable :: aj0v, aj1v
    ! (nmu, -kxkyz-layout-)

    ! Flux surface gyro averages
    type(coupled_alpha_type), dimension(:, :, :, :), allocatable :: j0_ffs, j0_B_ffs
    type(coupled_alpha_type), dimension(:, :, :, :), allocatable :: j1_ffs
    real, dimension(:, :, :, :), allocatable :: j0_const, j0_B_const, j0max_const

    logical :: bessinit = .false.

contains

    !###############################################################################
    !################## INITIALISE BESSEL FUNCTIONS FOR FLUX TUBE ##################
    !###############################################################################
    subroutine init_bessel

        use store_arrays_useful, only: kperp2
        use parameters_physics, only: full_flux_surface
        use species, only: spec
        use geometry, only: bmag
        use z_grid, only: nzgrid
        use velocity_grids, only: vperp2, nmu
        use parameters_kxky_grid, only: naky, nakx
        use stella_layouts, only: kxkyz_lo, vmu_lo
        use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, imu_idx
        use spfunc, only: j0, j1

        implicit none

        integer :: iz, iky, ikx, imu, is, ia
        integer :: ikxkyz, ivmu
        real :: arg

        if (bessinit) return
        bessinit = .true.

        if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::allocate_aj0v_aj1v'
        if (.not. allocated(aj0v)) then
            allocate (aj0v(nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
            aj0v = 0.
        end if
        if (.not. allocated(aj1v)) then
            allocate (aj1v(nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
            aj1v = 0.
        end if

        if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::calculate_aj0v_aj1v'
        ia = 1
        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            do imu = 1, nmu
                arg = spec(is)%bess_fac * spec(is)%smz_psi0 * sqrt(vperp2(ia, iz, imu) * kperp2(iky, ikx, ia, iz)) / bmag(ia, iz)
                aj0v(imu, ikxkyz) = j0(arg)
                ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)),
                aj1v(imu, ikxkyz) = j1(arg)
            end do
        end do

        if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::full_flux_surface'
        if (full_flux_surface) then
            call init_bessel_ffs 
        else
            if (.not. allocated(aj0x)) then
                allocate (aj0x(naky, nakx, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
                aj0x = 0.
            end if
            
            if (.not. allocated(aj1x)) then
                allocate (aj1x(naky, nakx, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
                aj1x = 0.
            end if

            ia = 1
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                is = is_idx(vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)
                do iz = -nzgrid, nzgrid
                do ikx = 1, nakx
                    do iky = 1, naky
                        arg = spec(is)%bess_fac * spec(is)%smz_psi0 * sqrt(vperp2(ia, iz, imu) * kperp2(iky, ikx, ia, iz)) / bmag(ia, iz)
                        aj0x(iky, ikx, iz, ivmu) = j0(arg)
                        ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)),
                        aj1x(iky, ikx, iz, ivmu) = j1(arg)
                    end do
                end do
                end do
            end do
        end if
        if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::test_gyro_average'

        ! if(debug_test_gyro_average) call test_gyro_average

    end subroutine init_bessel

    !###############################################################################
    !################## INITIALISE BESSEL FUNCTIONS FOR FULL FLUX ##################
    !###############################################################################

    subroutine init_bessel_ffs

        use mp, only: sum_allreduce, proc0
        use spfunc, only: j0
        use stella_layouts, only: vmu_lo
        use stella_layouts, only: iv_idx, imu_idx, is_idx
        use stella_transforms, only: transform_alpha2kalpha
        use species, only: nspec, spec
        use geometry, only: bmag
        use z_grid, only: nzgrid, nztot
        use velocity_grids, only: nmu, nvpa
        use velocity_grids, only: vperp2, maxwell_vpa, maxwell_mu
        use parameters_kxky_grid, only: nalpha, naky, naky_all, ikx_max
        use calculations_kxky, only: swap_kxky_ordered
        use store_arrays_useful, only: kperp2

        use parameters_kxky_grid, only: nakx
        use calculations_kxky, only: swap_kxky_back_ordered
        use spfunc, only: j1

        implicit none

        integer :: j0_ffs_unit, j0_B_ffs_unit
        integer :: iky, ikx, ia, iz
        integer :: ivmu, iv, imu, is
        integer :: ia_max_j0_count, ia_max_j0_B_count
        real :: arg, rtmp
        real :: ia_max_j0_reduction_factor, ia_max_j0_B_reduction_factor
        real, dimension(:), allocatable :: wgts
        real, dimension(:), allocatable :: aj0_alpha, j0_B
        real, dimension(:, :, :), allocatable :: kperp2_swap
        complex, dimension(:), allocatable :: aj0_kalpha, j0_B_kalpha

        real, dimension(:), allocatable :: aj1_alpha
        complex, dimension(:), allocatable :: aj1_kalpha
        integer :: ia_max_j1_count
        real :: ia_max_j1_reduction_factor

        complex, dimension(:, :), allocatable :: j0_const_in_kalpha, j0_B_const_in_kalpha
        complex, dimension(:, :), allocatable :: j0_const_c, j0_B_const_c

        real, dimension(:), allocatable :: j0max
        complex, dimension(:, :), allocatable :: j0max_const_in_kalpha, j0max_const_c

        ! if (debug_test_gyro_average) call open_output_file (j0_ffs_unit, '.j0_ffs')
        ! if (debug_test_gyro_average) call open_output_file (j0_B_ffs_unit, '.j0_over_B_ffs')

        ! wgts are species-dependent factors appearing in Gamma0 factor
        allocate (wgts(nspec))
        wgts = spec%dens * spec%z**2 / spec%temp

        if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::full_flux_surface::allocate_arrays'
        !> aj0_alpha will contain J_0 as a function of k_alpha and alpha
        allocate (aj0_alpha(nalpha))
        allocate (aj0_kalpha(naky))
        !> j0_B will contain J_0*B*exp(-v^2) as a function of k_alpha and alpha
        allocate (j0_B(nalpha))
        allocate (j0_B_kalpha(naky))
        allocate (kperp2_swap(naky_all, ikx_max, nalpha))
        if (.not. allocated(j0_ffs)) then
            allocate (j0_ffs(naky_all, ikx_max, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        end if
        if (.not. allocated(j0_B_ffs)) then
            allocate (j0_B_ffs(naky_all, ikx_max, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        end if

        allocate (aj1_alpha(nalpha))
        allocate (aj1_kalpha(naky))
        if (.not. allocated(j1_ffs)) then
            allocate (j1_ffs(naky_all, ikx_max, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        end if

        allocate (j0_const_in_kalpha(naky_all, ikx_max)); j0_const_in_kalpha = 0.0
        allocate (j0_B_const_in_kalpha(naky_all, ikx_max)); j0_B_const_in_kalpha = 0.0
        allocate (j0_const_c(naky, nakx)); j0_const_c = 0.0
        allocate (j0_B_const_c(naky, nakx)); j0_B_const_c = 0.0
        allocate (j0_const(naky, nakx, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); j0_const = 0.0
        allocate (j0_B_const(naky, nakx, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); j0_B_const = 0.0

        allocate (j0max(nalpha))
        allocate (j0max_const_in_kalpha(naky_all, ikx_max)); j0max_const_in_kalpha = 0.0
        allocate (j0max_const_c(naky, nakx)); j0max_const_c = 0.0
        allocate (j0max_const(naky, nakx, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); j0max_const = 0.0

        ia_max_j0_count = 0; ia_max_j0_B_count = 0
        do iz = -nzgrid, nzgrid
    !         if (proc0) write (*, *) 'calculating Fourier coefficients needed for gyro-averaging with alpha variation; zed index: ', iz
            !> for each value of alpha, take kperp^2 calculated on domain kx = [-kx_max, kx_max] and ky = [0, ky_max]
            !> and use symmetry to obtain kperp^2 on domain kx = [0, kx_max] and ky = [-ky_max, ky_max]
            !> this makes later convolutions involving sums over all ky more straightforward
            if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::full_flux_surface::swap_kxky'
            do ia = 1, nalpha
                call swap_kxky_ordered(kperp2(:, :, ia, iz), kperp2_swap(:, :, ia))
            end do
            if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::full_flux_surface::j0_loop'
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                is = is_idx(vmu_lo, ivmu)
                iv = iv_idx(vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)
                do ikx = 1, ikx_max
                do iky = 1, naky_all
                    do ia = 1, nalpha
                        !> calculate the argument of the Bessel function, which depends on both alpha and k_alpha
                        arg = spec(is)%bess_fac * spec(is)%smz_psi0 * sqrt(vperp2(ia, iz, imu) * kperp2_swap(iky, ikx, ia)) / bmag(ia, iz)
                        ! compute the value of the Bessel function J0 corresponding to argument arg
                        aj0_alpha(ia) = j0(arg)
                        !> compute J_0*B*exp(-v^2), needed when integrating g over v-space in Maxwell's equations,
                        !> due to B in v-space Jacobian
                        j0_B(ia) = aj0_alpha(ia) * bmag(ia, iz)
                        j0max = aj0_alpha(ia) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is)
                        aj1_alpha(ia) = j1(arg)
                    end do

                    !> fourier transform aj0_alpha and j0_B.
                    !> note that fourier coefficients aj0_kalpha and j0_B_kalpha have
                    !> been filtered to avoid aliasing
                    call transform_alpha2kalpha(aj0_alpha, aj0_kalpha)
                    call transform_alpha2kalpha(j0_B, j0_B_kalpha)
                    call transform_alpha2kalpha(aj1_alpha, aj1_kalpha)

                    j0_const_in_kalpha(iky, ikx) = aj0_kalpha(1)
                    j0_B_const_in_kalpha(iky, ikx) = j0_B_kalpha(1)
                    j0max_const_in_kalpha(iky, ikx) = j0max(1)
                        
                    !> given the Fourier coefficients aj0_kalpha, calculate the minimum number of coefficients needed,
                    !> called j0_ffs%max_idx, to ensure that the relative error in the total spectral energy is below a specified tolerance
                    !if (debug) write (*,*) 'arrays_gyro_averages:init_bessel::full_flux_surface::find_max_required_kalpha_index'
                    !                ! TMP FOR TESTING
                    !                j0_ffs(iky,ikx,iz,ivmu)%max_idx = naky
                    call find_max_required_kalpha_index(aj0_kalpha, j0_ffs(iky, ikx, iz, ivmu)%max_idx, imu, iz, is)
                    !> given the Fourier coefficients j0_B_kalpha, calculate the minimum number of coefficients needed,
                    !> called j0_B_ffs%max_idx, to ensure that the relative error in the total spectral energy is below a specified tolerance
                    call find_max_required_kalpha_index(j0_B_kalpha, j0_B_ffs(iky, ikx, iz, ivmu)%max_idx, imu, iz, is)
                    call find_max_required_kalpha_index(aj1_kalpha, j1_ffs(iky, ikx, iz, ivmu)%max_idx, imu, iz, is)
                    !> keep track of the total number of coefficients that must be retained across different phase space points
                    ia_max_j0_count = ia_max_j0_count + j0_ffs(iky, ikx, iz, ivmu)%max_idx
                    !> keep track of the total number of coefficients that must be retained across different phase space points
                    ia_max_j0_B_count = ia_max_j0_B_count + j0_B_ffs(iky, ikx, iz, ivmu)%max_idx
                    ia_max_j1_count = ia_max_j1_count + j1_ffs(iky, ikx, iz, ivmu)%max_idx

                    !> allocate array to hold the reduced number of Fourier coefficients
                    if (.not. associated(j0_ffs(iky, ikx, iz, ivmu)%fourier)) &
                        allocate (j0_ffs(iky, ikx, iz, ivmu)%fourier(j0_ffs(iky, ikx, iz, ivmu)%max_idx))
                    !> fill the array with the requisite coefficients
                    j0_ffs(iky, ikx, iz, ivmu)%fourier = aj0_kalpha(:j0_ffs(iky, ikx, iz, ivmu)%max_idx)
                    ! if (debug_test_gyro_average) call test_ffs_bessel_coefs (j0_ffs(iky,ikx,iz,ivmu)%fourier, aj0_alpha, iky, ikx, iz, j0_ffs_unit, ivmu)
                    if (.not. associated(j0_B_ffs(iky, ikx, iz, ivmu)%fourier)) &
                        allocate (j0_B_ffs(iky, ikx, iz, ivmu)%fourier(j0_B_ffs(iky, ikx, iz, ivmu)%max_idx))
                    !> fill the array with the requisite coefficients
                    j0_B_ffs(iky, ikx, iz, ivmu)%fourier = j0_B_kalpha(:j0_B_ffs(iky, ikx, iz, ivmu)%max_idx)
                    ! if (debug_test_gyro_average) call test_ffs_bessel_coefs (j0_B_ffs(iky,ikx,iz,ivmu)%fourier, j0_B, iky, ikx, iz, j0_B_ffs_unit, ivmu)
                    if (.not. associated(j1_ffs(iky, ikx, iz, ivmu)%fourier)) &
                        allocate (j1_ffs(iky, ikx, iz, ivmu)%fourier(j1_ffs(iky, ikx, iz, ivmu)%max_idx))
                    j1_ffs(iky, ikx, iz, ivmu)%fourier = aj1_kalpha(:j1_ffs(iky, ikx, iz, ivmu)%max_idx)
                end do
                end do
                call swap_kxky_back_ordered(j0_const_in_kalpha, j0_const_c)
                j0_const(:, :, iz, ivmu) = real(j0_const_c)
                call swap_kxky_back_ordered(j0_B_const_in_kalpha, j0_B_const_c)
                j0_B_const(:, :, iz, ivmu) = real(j0_B_const_c)
                call swap_kxky_back_ordered(j0max_const_in_kalpha, j0max_const_c)
                j0max_const(:, :, iz, ivmu) = real(j0max_const_c)
            end do
        end do

        deallocate (j0_B, j0_B_kalpha)
        deallocate (aj0_alpha)

        deallocate (j0_const_in_kalpha, j0_const_c)
        deallocate (j0_B_const_in_kalpha, j0_B_const_c)

        deallocate (j0max_const_in_kalpha, j0max_const_c, j0max)

        !> calculate the reduction factor of Fourier modes
        !> used to represent J0
        !> avoid overflow by converting integers to reals before multiplying
        rtmp = real(naky) * real(naky_all) * real(ikx_max) * real(nztot) * real(nmu) * real(nvpa) * real(nspec)
        call sum_allreduce(ia_max_j0_count)
        ia_max_j0_reduction_factor = real(ia_max_j0_count) / rtmp
        call sum_allreduce(ia_max_j0_B_count)
        ia_max_j0_B_reduction_factor = real(ia_max_j0_B_count) / rtmp

        call sum_allreduce(ia_max_j1_count)
        ia_max_j1_reduction_factor = real(ia_max_j1_count) / rtmp

        if (proc0) then
            write (*, *) 'average number of k-alphas needed to represent J0(kperp(alpha))=', ia_max_j0_reduction_factor * naky, 'out of ', naky
            write (*, *) 'average number of k-alphas needed to represent J0(kperp(alpha))*B(alpha)*exp(-v^2)=', &
                ia_max_j0_B_reduction_factor * naky, 'out of ', naky
            write (*, *)
        end if

        deallocate (wgts)
        deallocate (aj0_kalpha)
        deallocate (kperp2_swap)
        deallocate (aj1_alpha, aj1_kalpha)

        ! if (debug_test_gyro_average) call close_output_file (j0_ffs_unit)
        ! if (debug_test_gyro_average) call close_output_file (j0_B_ffs_unit)
        
    end subroutine init_bessel_ffs

    !> subroutine takes a set of Fourier coefficients (ft)
    !> and returns the minimum number of coeffients that must be retained (idx)
    !> to ensure that the relative error in the total spectral energy is
    !> below a specified tolerance (tol_floor)
    subroutine find_max_required_kalpha_index(ft, idx, imu, iz, is, tol_in)

        use velocity_grids, only: maxwell_mu

        implicit none

        complex, dimension(:), intent(in) :: ft
        integer, intent(out) :: idx
        integer, intent(in), optional :: imu, iz, is
        real, intent(in), optional :: tol_in

        real, parameter :: tol_floor = 1.0e-8
        integer :: i, n
        real :: subtotal, total
        real :: tol
        real, dimension(:), allocatable :: ftmod2

        n = size(ft)

        ! use conservative estimate
        ! when deciding number of modes to retain
        if (present(tol_in)) then
            tol = tol_in
        elseif (present(imu) .and. present(iz) .and. present(is)) then
            !       tol = min(0.1,tol_floor/maxval(maxwell_mu(:,iz,imu,is)))
            tol = min(1.0e-6, tol_floor / maxval(maxwell_mu(:, iz, imu, is)))
        else
            tol = tol_floor
        end if

        allocate (ftmod2(n))
        ! get spectral energy associated with each mode
        ftmod2 = sqrt(real(ft * conjg(ft)))
        ! get total spectral energy
        total = sqrt(sum(ftmod2))
        subtotal = 0.

        ! find minimum spectral index for which
        ! desired percentage of spectral energy contained
        ! in modes with indices at or below it
        if (total > 0.) then
            i = 1
            do while (subtotal < total * (1.0 - tol))
                idx = i
                subtotal = sqrt(sum(ftmod2(:i)))
                i = i + 1
            end do
        else
            idx = 1
        end if

        deallocate (ftmod2)

    end subroutine find_max_required_kalpha_index

    !###############################################################################
    !################### FINALISE BESSEL FUNCTIONS FOR FLUX TUBE ###################
    !###############################################################################

    subroutine finish_bessel

        implicit none

        if (allocated(aj0v)) deallocate (aj0v)
        if (allocated(aj1v)) deallocate (aj1v)
        if (allocated(aj0x)) deallocate (aj0x)
        if (allocated(aj1x)) deallocate (aj1x)
        if (allocated(j0_ffs)) deallocate (j0_ffs)
        if (allocated(j0_B_ffs)) deallocate (j0_B_ffs)

        if (allocated(j0_B_const)) deallocate (j0_B_const)
        if (allocated(j0_const)) deallocate (j0_const)
        if (allocated(j0max_const)) deallocate (j0max_const)

        bessinit = .false.

    end subroutine finish_bessel

end module arrays_gyro_averages