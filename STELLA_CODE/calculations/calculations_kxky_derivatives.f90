module calculations_kxky_derivatives

    implicit none

    public :: get_dgdy, get_dgdx
    public :: get_dchidy, get_dchidx

    private

    interface get_dgdy
      module procedure get_dgdy_2d
      module procedure get_dgdy_3d
      module procedure get_dgdy_4d
   end interface

   interface get_dgdx
      module procedure get_dgdx_2d
      module procedure get_dgdx_3d
      module procedure get_dgdx_4d
   end interface

   interface get_dchidy
      module procedure get_dchidy_4d
      module procedure get_dchidy_2d
   end interface get_dchidy

contains

    !> compute dg/dy in k-space
    !> accepts g(ky,kx)
    subroutine get_dgdy_2d(g, dgdy)

        use constants, only: zi
        use parameters_kxky_grid, only: nakx
        use grids_kxky, only: aky

        implicit none

        complex, dimension(:, :), intent(in) :: g
        complex, dimension(:, :), intent(out) :: dgdy

        dgdy = zi * spread(aky, 2, nakx) * g

    end subroutine get_dgdy_2d

    !> compute dg/dy in k-space
    !> accepts g(ky,kx,z,tube)
    subroutine get_dgdy_3d(g, dgdy)

        use constants, only: zi
        use parameters_kxky_grid, only: nakx
        use grids_kxky, only: aky
        use z_grid, only: nzgrid, ntubes
        
        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdy

        integer :: it, iz, ikx

        do it = 1, ntubes
            do iz = -nzgrid, nzgrid
                do ikx = 1, nakx
                dgdy(:, ikx, iz, it) = zi * aky(:) * g(:, ikx, iz, it)
                end do
            end do
        end do

    end subroutine get_dgdy_3d

    !> compute dg/dy in k-space
    !> accepts g(ky,kx,z,tube,(vpa,mu,spec))
    subroutine get_dgdy_4d(g, dgdy)

        use constants, only: zi
        use stella_layouts, only: vmu_lo
        use z_grid, only: nzgrid, ntubes
        use parameters_kxky_grid, only: nakx
        use grids_kxky, only: aky

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dgdy

        integer :: ivmu, ikx, iz, it

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
                do iz = -nzgrid, nzgrid
                do ikx = 1, nakx
                    dgdy(:, ikx, iz, it, ivmu) = zi * aky(:) * g(:, ikx, iz, it, ivmu)
                end do
                end do
            end do
        end do

    end subroutine get_dgdy_4d

    !> compute dg/dx in k-space
    !> accepts g(ky,kx)
    subroutine get_dgdx_2d(g, dgdx)

        use constants, only: zi
        use parameters_kxky_grid, only: naky
        use grids_kxky, only: akx

        implicit none

        complex, dimension(:, :), intent(in) :: g
        complex, dimension(:, :), intent(out) :: dgdx

        dgdx = zi * spread(akx, 1, naky) * g

    end subroutine get_dgdx_2d

    !> compute dg/dx in k-space
    !> accepts g(ky,kx,z,tube)
    subroutine get_dgdx_3d(g, dgdx)

        use constants, only: zi
        use z_grid, only: nzgrid, ntubes
        use parameters_kxky_grid, only: nakx
        use grids_kxky, only: akx

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdx

        integer :: ikx, iz, it

        do it = 1, ntubes
            do iz = -nzgrid, nzgrid
                do ikx = 1, nakx
                dgdx(:, ikx, iz, it) = zi * akx(ikx) * g(:, ikx, iz, it)
                end do
            end do
        end do

    end subroutine get_dgdx_3d

    !> compute dg/dx in k-space
    !> accepts g(ky,kx,z,tube,(vpa,mu,spec))
    subroutine get_dgdx_4d(g, dgdx)

        use constants, only: zi
        use stella_layouts, only: vmu_lo
        use z_grid, only: nzgrid, ntubes
        use parameters_kxky_grid, only: nakx
        use grids_kxky, only: akx

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dgdx

        integer :: ivmu, ikx, iz, it

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
                do iz = -nzgrid, nzgrid
                do ikx = 1, nakx
                    dgdx(:, ikx, iz, it, ivmu) = zi * akx(ikx) * g(:, ikx, iz, it, ivmu)
                end do
                end do
            end do
        end do

    end subroutine get_dgdx_4d



    !============================================================================
    !============================ FIELDS DERIVATIVES ============================
    !============================================================================
    !> Compute d<chi>/dy and d<chi>/dx in (ky,kx) space where <.> is a gyroaverage
    !>    d<chi>/dy = i * ky * J0 * chi
    !>    d<chi>/dx = i * kx * J0 * chi
    !>    chi = phi - Z/T * vpa * apar
    !> There are different routines depending on the size of the input array
    !============================================================================
    
    !> Compute d<chi>/dy in (ky,kx,z,tube) space
    subroutine get_dchidy_4d(phi, apar, bpar, dchidy)

        use constants, only: zi
        !> Layouts
        use stella_layouts, only: vmu_lo
        use stella_layouts, only: is_idx, iv_idx, imu_idx
        !> Parameters
        use parameters_physics, only: include_apar, include_bpar
        use parameters_physics, only: full_flux_surface
        use parameters_physics, only: fphi
        use parameters_kxky_grid, only: nakx, naky
        !> Grids
        use species, only: spec
        use z_grid, only: nzgrid, ntubes
        use velocity_grids, only: vpa, mu
        use grids_kxky, only: aky
        !> Calculations
        use calculations_gyro_averages, only: gyro_average
        use calculations_gyro_averages, only: gyro_average_j1
        use arrays_gyro_averages, only: j0_ffs

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dchidy

        integer :: ivmu, iv, is, iky, imu
        complex, dimension(:, :, :, :), allocatable :: field, gyro_tmp
        !-------------------------------------------------------------------------
        allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (gyro_tmp(naky, nakx, -nzgrid:nzgrid, ntubes))

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            ! intermediate calculation to get factor involving phi contribution
            field = fphi * phi
            ! add apar contribution if including it
            if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
            ! take spectral y-derivative
            do iky = 1, naky
                field(iky, :, :, :) = zi * aky(iky) * field(iky, :, :, :)
            end do
            if (full_flux_surface) then
                call gyro_average(field, dchidy(:, :, :, :, ivmu), j0_ffs(:, :, :, ivmu))
            else
                call gyro_average(field, ivmu, dchidy(:, :, :, :, ivmu))
            end if
            if (include_bpar) then
                field = 4.0 * mu(imu) * (spec(is)%tz) * bpar
                do iky = 1, naky
                field(iky, :, :, :) = zi * aky(iky) * field(iky, :, :, :)
                end do
                call gyro_average_j1(field, ivmu, gyro_tmp)
                !> include bpar contribution
                dchidy(:, :, :, :, ivmu) = dchidy(:, :, :, :, ivmu) + gyro_tmp
            end if
        end do

        deallocate (field)
        deallocate (gyro_tmp)

    end subroutine get_dchidy_4d

    !> Compute d<chi>/dy in (ky,kx) space
    subroutine get_dchidy_2d(iz, ivmu, phi, apar, bpar, dchidy)

        use constants, only: zi
        !> Layouts
        use stella_layouts, only: vmu_lo
        use stella_layouts, only: is_idx, iv_idx, imu_idx
        !> Parameters
        use parameters_physics, only: include_apar, include_bpar
        use parameters_physics, only: full_flux_surface
        use parameters_physics, only: fphi
        use parameters_kxky_grid, only: nakx, naky
        !> Grids
        use species, only: spec
        use velocity_grids, only: vpa, mu
        use grids_kxky, only: aky
        !> Calculations
        use calculations_gyro_averages, only: gyro_average
        use calculations_gyro_averages, only: gyro_average_j1
        use arrays_gyro_averages, only: j0_ffs

        implicit none

        integer, intent(in) :: ivmu, iz
        complex, dimension(:, :), intent(in) :: phi, apar, bpar
        complex, dimension(:, :), intent(out) :: dchidy

        integer :: iv, is, imu
        complex, dimension(:, :), allocatable :: field, gyro_tmp
        !-------------------------------------------------------------------------
        allocate (field(naky, nakx))
        allocate (gyro_tmp(naky, nakx))

        is = is_idx(vmu_lo, ivmu)
        iv = iv_idx(vmu_lo, ivmu)
        imu = imu_idx(vmu_lo, ivmu)
        field = fphi * phi
        if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
        field = zi * spread(aky, 2, nakx) * field

        if (full_flux_surface) then
            call gyro_average(field, dchidy, j0_ffs(:, :, iz, ivmu))
        else
            call gyro_average(field, iz, ivmu, dchidy)
        end if

        if (include_bpar) then
            field = 4.0 * mu(imu) * (spec(is)%tz) * bpar
            field = zi * spread(aky, 2, nakx) * field
            call gyro_average_j1(field, iz, ivmu, gyro_tmp)
            !> include bpar contribution
            dchidy = dchidy + gyro_tmp
        end if
        deallocate (field)
        deallocate (gyro_tmp)

    end subroutine get_dchidy_2d

    !> Compute d<chi>/dx in (ky,kx) space
    subroutine get_dchidx(iz, ivmu, phi, apar, bpar, dchidx)

        use constants, only: zi
        !> Layouts
        use stella_layouts, only: vmu_lo
        use stella_layouts, only: is_idx, iv_idx, imu_idx
        !> Parameters
        use parameters_kxky_grid, only: naky, nakx
        use parameters_physics, only: include_apar, include_bpar
        use parameters_physics, only: full_flux_surface
        use parameters_physics, only: fphi
        !> Grids
        use species, only: spec
        use velocity_grids, only: vpa, mu
        use grids_kxky, only: akx
        !> Calculations
        use calculations_gyro_averages, only: gyro_average
        use calculations_gyro_averages, only: gyro_average_j1
        use arrays_gyro_averages, only: j0_ffs

        implicit none

        integer, intent(in) :: ivmu, iz
        complex, dimension(:, :), intent(in) :: phi, apar, bpar
        complex, dimension(:, :), intent(out) :: dchidx

        integer :: iv, is, imu
        complex, dimension(:, :), allocatable :: field, gyro_tmp
        !-------------------------------------------------------------------------
        allocate (field(naky, nakx))
        allocate (gyro_tmp(naky, nakx))

        is = is_idx(vmu_lo, ivmu)
        iv = iv_idx(vmu_lo, ivmu)
        imu = imu_idx(vmu_lo, ivmu)
        field = fphi * phi
        if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
        field = zi * spread(akx, 1, naky) * field

        if (full_flux_surface) then
            call gyro_average(field, dchidx, j0_ffs(:, :, iz, ivmu))
        else
            call gyro_average(field, iz, ivmu, dchidx)
        end if

        if (include_bpar) then
            field = 4 * mu(imu) * (spec(is)%tz) * bpar
            field = zi * spread(akx, 1, naky) * field
            call gyro_average_j1(field, iz, ivmu, gyro_tmp)
            !> include bpar contribution
            dchidx = dchidx + gyro_tmp
        end if
        deallocate (field)
        deallocate (gyro_tmp)

    end subroutine get_dchidx

end module calculations_kxky_derivatives