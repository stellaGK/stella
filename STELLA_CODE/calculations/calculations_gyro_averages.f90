module calculations_gyro_averages

   use common_types, only: coupled_alpha_type
   use debug_flags, only: debug => gyro_averages_debug
   use arrays_gyro_averages, only: aj0x, aj1x, aj0v, aj1v
   use arrays_gyro_averages, only: j0_ffs

   implicit none 

   public :: gyro_average
   public :: gyro_average_j1
   public :: band_lu_solve_ffs, band_lu_factorisation_ffs

   private

   interface gyro_average
      module procedure gyro_average_local
      module procedure gyro_average_kxky_local
      module procedure gyro_average_kxkyz_local
      module procedure gyro_average_kxkyzv_local
      module procedure gyro_average_v_local
      module procedure gyro_average_vmu_local
      module procedure gyro_average_vmus_nonlocal
      module procedure gyro_average_ffs_kxky_local
      module procedure gyro_average_ffs_kxkyz_local
      module procedure gyro_average_ffs_field
      module procedure gyro_average_ffs
   end interface gyro_average

   interface gyro_average_j1
      module procedure gyro_average_j1_vmus_nonlocal
      module procedure gyro_average_j1_local
      module procedure gyro_average_j1_kxky_local
      module procedure gyro_average_j1_kxkyz_local
      module procedure gyro_average_j1_v_local
      module procedure gyro_average_j1_vmu_local
      module procedure gyro_average_j1_kxkyzv_local
   end interface

contains

   !> gyro_average_local takes a field at a given ky, kx, z and (vpa, mu, s) value
   !> and returns the gyro-average of that field;
   ! this should never be called for a full flux surface simulation, so no need
   ! to check if flux tube or full flux surface
   subroutine gyro_average_local(field, iky, ikx, iz, ivmu, gyro_field)

      implicit none

      complex, intent(in) :: field
      integer, intent(in) :: iky, ikx, iz, ivmu
      complex, intent(out) :: gyro_field

      gyro_field = aj0x(iky, ikx, iz, ivmu) * field

   end subroutine gyro_average_local

   subroutine gyro_average_kxky_local(field, iz, ivmu, gyro_field)

      use parameters_physics, only: full_flux_surface

      implicit none

      complex, dimension(:, :), intent(in) :: field
      integer, intent(in) :: iz, ivmu
      complex, dimension(:, :), intent(out) :: gyro_field

      if (full_flux_surface) then
         !> if simulating a full flux surface, the alpha dependence present
         !> in kperp makes gyro-averaging non-local in k-space
         call gyro_average(field, gyro_field, j0_ffs(:, :, iz, ivmu))
      else
         !> if simulating a flux tube, a gyro-average is local in k-space
         gyro_field = aj0x(:, :, iz, ivmu) * field
      end if

   end subroutine gyro_average_kxky_local

   subroutine gyro_average_kxkyz_local(field, ivmu, gyro_field)

      use parameters_physics, only: full_flux_surface
      use z_grid, only: nzgrid, ntubes

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: field
      integer, intent(in) :: ivmu
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: gyro_field

      if (full_flux_surface) then
         !> if simulating a full flux surface, the alpha dependence present
         !> in kperp makes gyro-averaging non-local in k-space
         call gyro_average(field, gyro_field, j0_ffs(:, :, :, ivmu))
      else
         !> if simulating a flux tube, a gyro-average is local in k-space
         gyro_field = spread(aj0x(:, :, :, ivmu), 4, ntubes) * field
      end if

   end subroutine gyro_average_kxkyz_local

   subroutine gyro_average_kxkyzv_local(field, gyro_field)

      use parameters_physics, only: full_flux_surface
      use z_grid, only: nzgrid
      use stella_layouts, only: vmu_lo

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: field
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: gyro_field
      integer :: ivmu

      if (full_flux_surface) then
         !> if simulating a full flux surface, the alpha dependence present
         !> in kperp makes gyro-averaging non-local in k-space
         call gyro_average(field, gyro_field, j0_ffs)
      else
         !> if simulating a flux tube, a gyro-average is local in k-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            call gyro_average(field(:, :, :, :, ivmu), ivmu, gyro_field(:, :, :, :, ivmu))
         end do
      end if

   end subroutine gyro_average_kxkyzv_local

   subroutine gyro_average_ffs_kxky_local(field, gyro_field, coefs)

      use parameters_kxky_grid, only: naky
      use parameters_kxky_grid, only: naky_all, ikx_max
      use calculations_kxky, only: swap_kxky_ordered, swap_kxky_back_ordered

      implicit none

      complex, dimension(:, :), intent(in) :: field
      complex, dimension(:, :), intent(out) :: gyro_field
      type(coupled_alpha_type), dimension(:, :), intent(in) :: coefs

      integer :: iky, ikx, ikyp
      integer :: idx
      complex, dimension(:, :), allocatable :: field_kyall, gyro_field_kyall

      ! need to switch from ky>=0 and all kx
      ! to kx>=0 and all ky (using reality condition)
      allocate (field_kyall(naky_all, ikx_max))
      allocate (gyro_field_kyall(naky_all, ikx_max)); gyro_field_kyall = 0.
      call swap_kxky_ordered(field, field_kyall)
      ! NB: J0(kx,ky) = J0(-kx,-ky) and Gamma0(kx,ky) = Gamma0(-kx,-ky)
      do ikx = 1, ikx_max
         do iky = 1, naky_all
            ! account for contributions from less positive ky values (and this ky itself)
            do ikyp = 1, min(naky, iky)
               ! idx is the index corresponding to k_alpha - k_alpha'
               ! runs from iky down to 1
               idx = iky - ikyp + 1
               ! if the Fourier coefficient corresponding to this value of (k_alpha-k_alpha',k_alpha')
               ! is sufficiently small, then it will not have been included in the truncated version
               ! of the coefficients; in this case, it makes no contribution to the gyro-average sum
               if (coefs(idx, ikx)%max_idx >= ikyp) then
                  gyro_field_kyall(iky, ikx) = gyro_field_kyall(iky, ikx) &
                                               + coefs(idx, ikx)%fourier(ikyp) * field_kyall(idx, ikx)
               end if
            end do
            ! if iky = naky_all, then already at max positive ky, so no contributions
            ! from more positive ky value possible
            if (iky == naky_all) cycle
            ! account for contributions from more positive ky values (but not this ky itself,
            ! as already accounted for above
            do ikyp = 2, min(naky, naky_all - iky + 1)
               ! idx is the index corresponding to k_alpha - k_alpha'
               ! runs from iky + 1 up to iky + naky (or until naky_all, if it is reached first)
               idx = iky + ikyp - 1
               ! if the Fourier coefficient corresponding to this value of (k_alpha-k_alpha',k_alpha')
               ! is sufficiently small, then it will not have been included in the truncated version
               ! of the coefficients; in this case, it makes no contribution to the gyro-average sum
               if (coefs(idx, ikx)%max_idx >= ikyp) then
                  ! the k_alpha' values considered in this loop are negative, but only have
                  ! Fourier coefficients for positive ky values;
                  ! must use the reality condition to convert this to the equivalent coefficients for negative ky
                  gyro_field_kyall(iky, ikx) = gyro_field_kyall(iky, ikx) &
                                               + conjg(coefs(idx, ikx)%fourier(ikyp)) * field_kyall(idx, ikx)
               end if
            end do
         end do
      end do

      call swap_kxky_back_ordered(gyro_field_kyall, gyro_field)
      deallocate (field_kyall, gyro_field_kyall)

   end subroutine gyro_average_ffs_kxky_local

   subroutine gyro_average_ffs_kxkyz_local(field, gyro_field, coefs)

      use z_grid, only: nzgrid, ntubes

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: field
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: gyro_field
      type(coupled_alpha_type), dimension(:, :, -nzgrid:), intent(in) :: coefs

      integer :: iz, it

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            call gyro_average(field(:, :, iz, it), gyro_field(:, :, iz, it), coefs(:, :, iz))
         end do
      end do

   end subroutine gyro_average_ffs_kxkyz_local

   subroutine gyro_average_ffs_field(field, gyro_field, coefs)

      use common_types, only: coupled_alpha_type
      use stella_layouts, only: vmu_lo
      use z_grid, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: field
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: gyro_field
      type(coupled_alpha_type), dimension(:, :, -nzgrid:, vmu_lo%llim_proc:), intent(in) :: coefs

      integer :: ivmu

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         call gyro_average(field(:, :, :, :), gyro_field(:, :, :, :, ivmu), coefs(:, :, :, ivmu))
      end do

   end subroutine gyro_average_ffs_field

   subroutine gyro_average_ffs(dist, gyro_dist, coefs)

      use common_types, only: coupled_alpha_type
      use stella_layouts, only: vmu_lo
      use z_grid, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: dist
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: gyro_dist
      type(coupled_alpha_type), dimension(:, :, -nzgrid:, vmu_lo%llim_proc:), intent(in) :: coefs

      integer :: ivmu

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         call gyro_average(dist(:, :, :, :, ivmu), gyro_dist(:, :, :, :, ivmu), coefs(:, :, :, ivmu))
      end do

   end subroutine gyro_average_ffs

   subroutine gyro_average_v_local(distfn, imu, ikxkyz, gyro_distfn)

      implicit none

      complex, dimension(:), intent(in) :: distfn
      integer, intent(in) :: imu, ikxkyz
      complex, dimension(:), intent(out) :: gyro_distfn

      gyro_distfn = aj0v(imu, ikxkyz) * distfn

   end subroutine gyro_average_v_local

   subroutine gyro_average_vmu_local(distfn, ikxkyz, gyro_distfn)

      use velocity_grids, only: nvpa

      implicit none

      complex, dimension(:, :), intent(in) :: distfn
      integer, intent(in) :: ikxkyz
      complex, dimension(:, :), intent(out) :: gyro_distfn

      gyro_distfn = spread(aj0v(:, ikxkyz), 1, nvpa) * distfn

   end subroutine gyro_average_vmu_local

   subroutine gyro_average_vmus_nonlocal(field, iky, ikx, iz, gyro_field)

      use stella_layouts, only: vmu_lo

      implicit none

      complex, dimension(vmu_lo%llim_proc:), intent(in) :: field
      integer, intent(in) :: iky, ikx, iz
      complex, dimension(vmu_lo%llim_proc:), intent(out) :: gyro_field

      gyro_field = aj0x(iky, ikx, iz, :) * field

   end subroutine gyro_average_vmus_nonlocal
   
   subroutine gyro_average_j1_vmus_nonlocal(field, iky, ikx, iz, gyro_field)

      use stella_layouts, only: vmu_lo

      implicit none

      complex, dimension(vmu_lo%llim_proc:), intent(in) :: field
      integer, intent(in) :: iky, ikx, iz
      complex, dimension(vmu_lo%llim_proc:), intent(out) :: gyro_field

      gyro_field = aj1x(iky, ikx, iz, :) * field

   end subroutine gyro_average_j1_vmus_nonlocal
   
   subroutine gyro_average_j1_local(field, iky, ikx, iz, ivmu, gyro_field)

      implicit none

      complex, intent(in) :: field
      integer, intent(in) :: iky, ikx, iz, ivmu
      complex, intent(out) :: gyro_field

      gyro_field = aj1x(iky, ikx, iz, ivmu) * field

   end subroutine gyro_average_j1_local
   
   subroutine gyro_average_j1_kxky_local(field, iz, ivmu, gyro_field)

      implicit none

      complex, dimension(:, :), intent(in) :: field
      integer, intent(in) :: iz, ivmu
      complex, dimension(:, :), intent(out) :: gyro_field

      gyro_field = aj1x(:, :, iz, ivmu) * field

   end subroutine gyro_average_j1_kxky_local

   subroutine gyro_average_j1_kxkyz_local(field, ivmu, gyro_field)

      use z_grid, only: nzgrid, ntubes

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: field
      integer, intent(in) :: ivmu
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: gyro_field

      integer :: iz, it

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            gyro_field(:, :, iz, it) = aj1x(:, :, iz, ivmu) * field(:, :, iz, it)
         end do
      end do

   end subroutine gyro_average_j1_kxkyz_local

   subroutine gyro_average_j1_vmu_local(distfn, ikxkyz, gyro_distfn)

      use velocity_grids, only: nvpa

      implicit none

      complex, dimension(:, :), intent(in) :: distfn
      integer, intent(in) :: ikxkyz
      complex, dimension(:, :), intent(out) :: gyro_distfn

      gyro_distfn = spread(aj1v(:, ikxkyz), 1, nvpa) * distfn

   end subroutine gyro_average_j1_vmu_local

   subroutine gyro_average_j1_v_local(distfn, imu, ikxkyz, gyro_distfn)

      implicit none

      complex, dimension(:), intent(in) :: distfn
      integer, intent(in) :: imu, ikxkyz
      complex, dimension(:), intent(out) :: gyro_distfn

      gyro_distfn = aj1v(imu, ikxkyz) * distfn

   end subroutine gyro_average_j1_v_local

   subroutine gyro_average_j1_kxkyzv_local(field, gyro_field)
      use mp, only: proc0, mp_abort
      use z_grid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use parameters_physics, only: full_flux_surface
      
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: field
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: gyro_field
      integer :: ivmu

      if (full_flux_surface) then
         if (proc0) write (*, *) 'gyro_average_j1_kxkyzv_local does not support full_flux_surface'
         call mp_abort('gyro_average_j1_kxkyzv_local does not support full_flux_surface. aborting')
         return
         !> if simulating a full flux surface, the alpha dependence present
         !> in kperp makes gyro-averaging non-local in k-space
         ! call gyro_average(field, gyro_field, j0_ffs)
      else
         !> if simulating a flux tube, a gyro-average is local in k-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            call gyro_average_j1(field(:, :, :, :, ivmu), ivmu, gyro_field(:, :, :, :, ivmu))
         end do
      end if
   end subroutine 
   
   subroutine band_lu_solve_ffs(lu, solvec)

      use common_types, only: gam0_ffs_type
      use z_grid, only: nzgrid
      use parameters_kxky_grid, only: ikx_max

      implicit none

      type(gam0_ffs_type), dimension(:, -nzgrid:), intent(in) :: lu
      complex, dimension(:, :, -nzgrid:), intent(in out) :: solvec

      integer :: ikx, iz

      do iz = -nzgrid, nzgrid
         do ikx = 1, ikx_max
            call band_lu_solve_ffs_single(lu(ikx, iz), solvec(:, ikx, iz))
         end do
      end do

   end subroutine band_lu_solve_ffs

   subroutine band_lu_solve_ffs_single(lu, solvec)

      use common_types, only: gam0_ffs_type
      use parameters_kxky_grid, only: naky

      implicit none

      type(gam0_ffs_type), intent(in) :: lu
      complex, dimension(:), intent(in out) :: solvec

      integer :: n, nsubdiag, nsupdiag, nrhs
      integer :: info
      complex, dimension(:, :), allocatable :: solmat

      ! n is the order of the matrix for which we have the LU factorisation
      n = size(lu%pivot_index)
      ! nsubdiag and nsupdiag are the number of sub- and super-diagonals within the band of the matrix to be decomposed
      nsubdiag = naky - 1; nsupdiag = naky - 1
      ! nrhs is the number of right-hand sides of the matrix equation lu%matrix * solvec = rhs for which to solve
      nrhs = 1
      ! initialise solmat = rhs, as it will be overwritten by zgbtrs below
      allocate (solmat(size(solvec), 1))
      solmat(:, 1) = solvec

      call zgbtrs('N', n, nsubdiag, nsupdiag, nrhs, lu%matrix, size(lu%matrix, 1), lu%pivot_index, solmat, size(solmat, 1), info)

      solvec = solmat(:, 1)

      deallocate (solmat)

   end subroutine band_lu_solve_ffs_single

   subroutine band_lu_factorisation_ffs(gam0, lu_gam0)

      use common_types, only: coupled_alpha_type, gam0_ffs_type
      use z_grid, only: nzgrid
      use parameters_kxky_grid, only: ikx_max, naky_all, naky

      implicit none

      type(coupled_alpha_type), dimension(:, :, -nzgrid:), intent(in) :: gam0
      type(gam0_ffs_type), dimension(:, -nzgrid:), intent(out) :: lu_gam0

      integer :: iky, ikx, iz

      complex, dimension(:, :), allocatable :: gam0tmp

      allocate (gam0tmp(naky, naky_all))

      do ikx = 1, ikx_max
         do iz = -nzgrid, nzgrid
            ! create array from Fourier coefficients of Gamma_0(ky,y)
            do iky = 1, naky_all
               gam0tmp(:, iky) = gam0(iky, ikx, iz)%fourier
            end do
            call band_lu_factorisation_single(gam0tmp, lu_gam0(ikx, iz))
         end do
      end do

      deallocate (gam0tmp)

   end subroutine band_lu_factorisation_ffs

   subroutine band_lu_factorisation_single(gam0, lu_gam0)

      use common_types, only: gam0_ffs_type
      use parameters_kxky_grid, only: naky, naky_all

      implicit none

      complex, dimension(:, :), intent(in) :: gam0
      type(gam0_ffs_type), intent(out) :: lu_gam0

      integer :: nrows, ncols, nsubdiag, nsupdiag, leading_dim
      integer :: info
      integer :: i, imod

      ! nrows and ncols are the number of rows and columns of the matrix to be LU-decomposed (variant of gam0)
      ! this matrix is naky_all x naky_all
      nrows = naky_all; ncols = naky_all
      ! nsubdiag and nsupdiag are the number of sub- and super-diagonals within the band of the matrix to be decomposed
      nsubdiag = naky - 1; nsupdiag = naky - 1
      ! leading_dim is the 'leading dimension' of the lu_gam0 array
      leading_dim = 2 * nsubdiag + nsupdiag + 1

      ! lu_gam0 is a re-arranged version of gam0 on entry, and on exit contains details of LU factorisation
      ! that can be used by zgbtrs to solve the linear system gam0 * phi = rhs
      if (.not. associated(lu_gam0%matrix)) then
         allocate (lu_gam0%matrix(leading_dim, ncols))
         ! initialise first nsubdiag rows to zero, as they are unused
         lu_gam0%matrix(:nsubdiag, :) = 0.0
         ! fill next supdiag rows using elements from super-diagonals
         ! using reality of gam0 to set fourier(ky < 0) = conjg(fourier(ky > 0))
         do i = 1, nsupdiag
            imod = naky - i + 1
            lu_gam0%matrix(nsubdiag + i, imod:) = conjg(gam0(imod, imod:))
            ! fill unused entries with zero
            lu_gam0%matrix(nsubdiag + i, :imod - 1) = 0.0
         end do
         ! fill next row using main diagonal entries
         lu_gam0%matrix(nsubdiag + nsupdiag + 1, :) = gam0(1, :)
         ! fill remaining nsubdiag rows using elements from sub-diagonals
         do i = 1, nsubdiag
            imod = naky + i - 1
            lu_gam0%matrix(leading_dim - i + 1, :imod) = gam0(naky - i + 1, :imod)
            ! fill unused entries with zeroes
            lu_gam0%matrix(leading_dim - i + 1, imod + 1:) = 0.0
         end do
      end if

      if (.not. associated(lu_gam0%pivot_index)) allocate (lu_gam0%pivot_index(min(nrows, ncols)))
      ! overwrites lu_gam0%matrix with information needed to solve the linear system gam0 * phi = rhs;
      ! also returns pivot_index, with pivot_index(i) giving the row of the matrix interchanged with the ith row,
      ! and info, which should be zero if the LU factorisation is successful
      call zgbtrf(nrows, ncols, nsubdiag, nsupdiag, lu_gam0%matrix, leading_dim, lu_gam0%pivot_index, info)

   end subroutine band_lu_factorisation_single

   ! set up field that varies as x^2 = rho^2 * cos(angle)^2,
        ! with rho the distance from the origin, and 'angle' is the angle made with the horizontal
        ! if considering a particle at x=0, then rho is thee gyro-radius and angle is the gyro-angle
        ! the gyro-average should theen be 1/(2pi) * int_0^2pi dangle rho^2 * cos(angle)^2 = rho^2/2
    !     subroutine test_gyro_average

    !     use constants, only: pi
    !     use parameters_kxky_grid, only: ny, nx, x0, y0
    !     use parameters_kxky_grid, only: nakx, ikx_max, naky, naky_all
    !     use grids_kxky, only: x, y
    !     use calculations_kxky, only: swap_kxky, swap_kxky_back
    !     use stella_transforms, only: transform_x2kx, transform_y2ky
    !     use stella_transforms, only: transform_kx2x, transform_ky2y
    !     use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
    !     use velocity_grids, only: nmu
    !     use species, only: nspec, spec
    !     use geometry, only: alpha, bmag, x_displacement_fac

    !     implicit none

    !     real, dimension (:,:), allocatable :: fld_yx
    !     complex, dimension (:,:), allocatable :: fld_ykx
    !     complex, dimension (:,:), allocatable :: fld_kykx_swapped
    !     complex, dimension (:,:), allocatable :: fld_kykx, gyro_fld
    !     real, dimension (:,:,:,:), allocatable :: gyro_fld_yx
    !     real :: gyroradius
    !     integer :: iy, ix, ivmu, iv, imu, is

    !     integer, parameter :: iz = 0

    !     allocate (fld_yx(ny,nx))
    !     allocate (fld_ykx(ny,ikx_max))
    !     allocate (fld_kykx_swapped(naky_all,ikx_max))
    !     allocate (fld_kykx(naky,nakx))
    !     allocate (gyro_fld(naky,nakx))
    !     allocate (gyro_fld_yx(ny,nx,nmu,nspec))

    !     if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_averages::fld_yx'
    !     ! set up field that varies as x^2 = rho^2 * cos(angle)^2 and is constant in y
    ! !      fld_yx = spread(0.1*(x-pi*x0),1,ny)**2
    !     fld_yx = spread(cos(50.0*(x/x0-pi)),1,ny)
    ! !      fld_yx = spread(exp(-0.1*(x-pi*x0)**2),1,ny)
    !     if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_averages::transform_x2kx'
    !     ! transform from (y,x) to (y,kx), with kx going from 0 to kxmax
    !     call transform_x2kx (fld_yx, fld_ykx)
    !     if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_averages::transform_y2ky'
    !     ! transform from (y,kx) to (ky,kx), with ky going from (0,...,kymax,-kymax,...,-dky)
    !     ! and kx going from 0 to kxmax
    !     call transform_y2ky (fld_ykx, fld_kykx_swapped)
    !     if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_averages::swap_kxky_back'
    !     ! use reality condition to re-arrange array so that ky goes from 0 to kymax
    !     ! and kx goes from (0,...,kxmax,-kxmax,...,-dkx)
    !     call swap_kxky_back (fld_kykx_swapped, fld_kykx)

    !     ! gyro-average the field at z=0 for different values of mu
    !     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
    !         ! get the vpa index
    !         iv = iv_idx(vmu_lo,ivmu)
    !         ! as J0 independent of vpa, pick only one vpa to test
    !         if (iv /= 1) cycle
    !         ! get the species index
    !         is = is_idx(vmu_lo,ivmu)
    !         ! get the mu index
    !         imu = imu_idx(vmu_lo,ivmu)
    ! !         if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_averages::gyro_average'
    !         if (full_flux_surface) then
    !             call gyro_average (fld_kykx, gyro_fld, j0_ffs(:,:,iz,ivmu))
    !         else
    !             call gyro_average (fld_kykx, iz, ivmu, gyro_fld)
    !         end if
    !         ! use reality to re-arrange array entries so that ky goes from (0,...,kymax,-kymax,...,-dky)
    !         ! and kx goes from 0 to kxmax
    !         call swap_kxky (gyro_fld, fld_kykx_swapped)
    !         ! transform from (ky,kx) to (y,kx)
    !         call transform_ky2y (fld_kykx_swapped, fld_ykx)
    !         ! transform from (y,kx) to (y,x)
    !         call transform_kx2x (fld_ykx, gyro_fld_yx(:,:,imu,is))
    !     end do
    !     if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_averages::write_to_screen'
    !     ! NB: this is only set up to work on a single processor at the moment
    !     ! NB: to extend, must move information about gyro_fld onto proc0
    !     do is = 1, nspec
    !         do imu = 1, nmu
    !             do ix = 1, nx
    !             do iy = 1, ny
    !                 ! gyro-radius/reference gyro-radius is v_perp/Omega/rho_ref = (v_perp/vths)*(rho_s/rho_ref)
    !                 ! = vperp * sqrt(T_s/T_ref*m_s/m_ref)(B_ref/Z*B) = vperp / (spec%zstm*bmag)
    !                 gyroradius = sqrt(vperp2(iy,iz,imu))/(spec(is)%zstm*bmag(iy,iz))
    ! !                  gyroradius = sqrt(vperp2(1,iz,imu))/(spec(is)%zstm*bmag(1,iz))
    !                 write (42,*) 'y: ', y(iy), 'x: ', x(ix)-x0*pi, 'gyro_fld: ', gyro_fld_yx(iy,ix,imu,is), 'gyroradius: ', gyroradius, 'spec: ', is, &
    !                     'alpha: ', alpha(iy), 'x_displacement_fac: ', x_displacement_fac(iy,iz)
    !             end do
    !             write (42,*)
    !             end do
    !             write (42,*)
    !         end do
    !         write (42,*)
    !     end do
    !     do iy = 1, 1000
    !         gyroradius = (iy-1)*15.0/999.
    !         !         write (43,*) 'gyroradius: ', gyroradius, 'bes: ', bessi0(0.1*0.5*(gyroradius/x_displacement_fac(1,iz))**2)*exp(-0.1*0.5*(gyroradius/x_displacement_fac(1,iz))**2)
    !         !         write (43,*) 'gyroradius: ', gyroradius, 'analytical: ', 0.5*(0.1*gyroradius/x_displacement_fac(1,iz))**2
    !         write (43,*) 'gyroradius: ', gyroradius, 'analytical: ', j0(50.*gyroradius/(x0*x_displacement_fac(1,iz)))
    !     end do
    !     ! TMP FOR TESTING
    ! !      stop

    !     deallocate (fld_yx, fld_ykx)
    !     deallocate (fld_kykx_swapped, fld_kykx)
    !     deallocate (gyro_fld, gyro_fld_yx)

    !     end subroutine test_gyro_average
    

   ! subroutine test_band_lu_factorisation (gam0, lu_gam0)

   !   use common_types, only: coupled_alpha_type, gam0_ffs_type
   !   use z_grid, only: nzgrid
   !   use kt_grids, only: naky_all, naky

   !   implicit none

   !   type (coupled_alpha_type), dimension (:,:,-nzgrid:), intent (in) :: gam0
   !   type (gam0_ffs_type), dimension (:,-nzgrid:), intent (out) :: lu_gam0

   !   integer :: iky, ikx, ikyp, iz
   !   complex, dimension (naky_all) :: solvec

   !   ikx = 1 ; iz = -nzgrid
   !   do iky = 1, naky_all
   !      do ikyp = 1, naky
   !         gam0(iky,ikx,iz)%fourier(ikyp) = iky-naky + ikyp-1
   !      end do
   !   end do

   !   call band_lu_factorisation_ffs (gam0, lu_gam0)

   !   do iky = 1, naky_all
   !      solvec(iky) = iky
   !   end do

   !   call band_lu_solve_ffs_single (lu_gam0(ikx,iz), solvec)

   !   do iky = 1, naky_all
   !      write (*,*) 'iky: ', iky, 'solution: ', solvec(iky)
   !   end do
   !   stop

   ! end subroutine test_band_lu_factorisation


   !> inverse fourier transform coefs%fourier for several phase space points and compare with
        !> unfiltered version in alpha-space
        ! subroutine test_ffs_bessel_coefs (coefs, f_alpha, iky, ikx, iz, unit, ivmu)

        !     use stella_layouts, only: vmu_lo, iv_idx, is_idx, imu_idx

        !     implicit none

        !     complex, dimension (:), intent (in) :: coefs
        !     real, dimension (:), intent (in) :: f_alpha
        !     integer, intent (in) :: iky, ikx, iz
        !     integer, intent (in) :: unit
        !     integer, intent (in), optional :: ivmu

        !     integer :: iv, imu, is

        !     if (present(ivmu)) then
        !     !> coefficients should all be independent of vpa, so only do comparison for one vpa point
        !     iv = iv_idx(vmu_lo,ivmu)
        !     if (iv == 1) then
        !         !> only sample subset of mu locations
        !         imu = imu_idx(vmu_lo,ivmu)
        !         if (mod(imu-1,nmu/2-1)==0) then
        !             is = is_idx(vmu_lo,ivmu)
        !             call test_ffs_bessel_coefs_subset (coefs, f_alpha, iky, ikx, iz, unit, iv, imu, is)
        !         end if
        !     end if
        !     else
        !     call test_ffs_bessel_coefs_subset (coefs, f_alpha, iky, ikx, iz, unit)
        !     end if

        ! end subroutine test_ffs_bessel_coefs

        ! subroutine test_ffs_bessel_coefs_subset (coefs, f_alpha, iky, ikx, iz, unit, iv, imu, is)

        !     use constants, only: pi
        !     use z_grid, only: nzgrid, zed
        !     use parameters_kxky_grid, only: naky, nalpha
        !     use grids_kxky, only: aky_all_ordered
        !     use velocity_grids, only: mu
        !     use stella_transforms, only: transform_kalpha2alpha
        !     use geometry, only: alpha

        !     implicit none

        !     complex, dimension (:), intent (in) :: coefs
        !     real, dimension (:), intent (in) :: f_alpha
        !     integer, intent (in) :: iky, ikx, iz
        !     integer, intent (in) :: unit
        !     integer, intent (in), optional :: iv, imu, is

        !     complex, dimension (:), allocatable :: coefs_padded
        !     real, dimension (:), allocatable :: f_alpha_approx

        !     integer :: ia
        !     integer :: max_idx
        !     real :: relative_error
        !     real, parameter :: minval = 1.0e-3

        !     ! only sample a subset of z locations
        !     if (mod(iz,nzgrid/2)==0) then
        !     ! consider only kx = 0
        !     if (ikx == 1) then
        !         allocate (coefs_padded(naky))
        !         allocate (f_alpha_approx(nalpha))
        !         ! initialize the padded coefficient array to zero
        !         coefs_padded = 0.0
        !         ! fill in non-zero entries with truncated Fourier coefficients
        !         max_idx = size(coefs)
        !         coefs_padded(:max_idx) = coefs
        !         ! inverse Fourier transform to get alpha-dependent function
        !         call transform_kalpha2alpha (coefs_padded, f_alpha_approx)
        !         if (present(iv)) then
        !             do ia = 1, nalpha
        !                 relative_error = 2.0*abs(f_alpha(ia)-f_alpha_approx(ia))/(abs(f_alpha(ia)) + abs(f_alpha_approx(ia)))
        !                 write (unit,*) alpha(ia), f_alpha(ia), f_alpha_approx(ia), &
        !                     relative_error, aky_all_ordered(iky), ikx, iz, zed(iz), iv, imu, mu(imu), is
        !             end do
        !             ! user 2*pi periodicity in alpha to fill in final point (for visualization purposes)
        !             ia = 1
        !             relative_error = 2.0*abs(f_alpha(ia)-f_alpha_approx(ia))/(abs(f_alpha(ia)) + abs(f_alpha_approx(ia)))
        !             write (unit,*) 2.0*pi, f_alpha(1), f_alpha_approx(1), &
        !                 relative_error, aky_all_ordered(iky), ikx, iz, zed(iz), iv, imu, mu(imu), is
        !         else
        !             do ia = 1, nalpha
        !                 relative_error = 2.0*abs(f_alpha(ia)-f_alpha_approx(ia))/(abs(f_alpha(ia)) + abs(f_alpha_approx(ia)))
        !                 write (unit,*) alpha(ia), f_alpha(ia), f_alpha_approx(ia), &
        !                     relative_error, aky_all_ordered(iky), ikx, iz, zed(iz)
        !             end do
        !             ! user 2*pi periodicity in alpha to fill in final point (for visualization purposes)
        !             ia = 1
        !             relative_error = 2.0*abs(f_alpha(ia)-f_alpha_approx(ia))/(abs(f_alpha(ia)) + abs(f_alpha_approx(ia)))
        !             write (unit,*) 2.0*pi, f_alpha(1), f_alpha_approx(1), &
        !                 relative_error, aky_all_ordered(iky), ikx, iz, zed(iz)
        !         end if
        !         write (unit,*)
        !         deallocate (coefs_padded, f_alpha_approx)
        !     end if
        !     end if

        ! end subroutine test_ffs_bessel_coefs_subset
end module calculations_gyro_averages
