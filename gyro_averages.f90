module gyro_averages

  public :: aj0x, aj0v, aj1v
  public :: init_bessel, finish_bessel
  public :: gyro_average
  public :: gyro_average_j1
  
  private

  interface gyro_average
     module procedure gyro_average_kxky_local
     module procedure gyro_average_kxkyz_local
     module procedure gyro_average_vmu_local
     module procedure gyro_average_vmus_nonlocal
  end interface

  real, dimension (:,:,:,:), allocatable :: aj0x
  ! (naky, nakx, nalpha, -nzgrid:nzgrid, -vmu-layout-)

  real, dimension (:,:), allocatable :: aj0v, aj1v
  ! (nmu, -kxkyz-layout-)

  integer, dimension (:,:,:,:), allocatable :: ia_max_aj0a
  complex, dimension (:,:,:,:,:), allocatable :: aj0a

  logical :: bessinit = .false.

contains

  subroutine init_bessel

    use mp, only: sum_allreduce, proc0
    use dist_fn_arrays, only: kperp2
    use physics_flags, only: full_flux_surface
    use species, only: spec, nspec
    use stella_geometry, only: bmag
    use zgrid, only: nzgrid, nztot
    use vpamu_grids, only: vperp2, nmu, nvpa
    use kt_grids, only: naky, nakx, nalpha
    use stella_layouts, only: kxkyz_lo, vmu_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, imu_idx
    use spfunc, only: j0, j1
    use stella_transforms, only: transform_alpha2kalpha
!    use stella_transforms, only: transform_kalpha2alpha

    implicit none

    integer :: iz, iky, ikx, imu, is, ia
    integer :: ikxkyz, ivmu
    real :: arg
    integer :: ia_max_aj0a_count
    real :: ia_max_aj0a_reduction_factor

    real, dimension (:), allocatable :: aj0_alpha
    complex, dimension (:), allocatable :: aj0_kalpha

    if (bessinit) return
    bessinit = .true.

    if (.not.allocated(aj0v)) then
       allocate (aj0v(nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       aj0v = 0.
    end if
    if (.not.allocated(aj0x)) then
       allocate (aj0x(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       aj0x = 0.
    end if
    if (.not.allocated(aj1v)) then
       allocate (aj1v(nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       aj1v = 0.
    end if
    
    ia = 1
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          arg = spec(is)%smz*sqrt(vperp2(ia,iz,imu)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
          aj0v(imu,ikxkyz) = j0(arg)
          ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
          aj1v(imu,ikxkyz) = j1(arg)
       end do
    end do

    if (full_flux_surface) then
       allocate (aj0_alpha(nalpha))
       allocate (aj0_kalpha(naky))
       if (.not.allocated(ia_max_aj0a)) allocate(ia_max_aj0a(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       if (.not.allocated(aj0a)) then
          allocate(aj0a(naky,naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
          aj0a = 0.
       end if
          
       ia_max_aj0a_count = 0
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          do iz = -nzgrid, nzgrid
             do ikx = 1, nakx
                do iky = 1, naky
                   do ia = 1, nalpha
                      arg = spec(is)%smz*sqrt(vperp2(ia,iz,imu)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
                      aj0x(iky,ikx,iz,ivmu) = j0(arg)
                      aj0_alpha(ia) = aj0x(iky,ikx,iz,ivmu)
                   end do
                   ! fourier transform aj0_alpha
                   ! note that fourier coefficients aj0_kalpha have
                   ! been filter to avoid aliasing
                   call transform_alpha2kalpha (aj0_alpha, aj0_kalpha)
                   call find_max_required_kalpha_index (aj0_kalpha, imu, iz, ia_max_aj0a(iky,ikx,iz,ivmu))
                   ia_max_aj0a_count = ia_max_aj0a_count + ia_max_aj0a(iky,ikx,iz,ivmu)
                   aj0a(:ia_max_aj0a(iky,ikx,iz,ivmu),iky,ikx,iz,ivmu) = aj0_kalpha(:ia_max_aj0a(iky,ikx,iz,ivmu))
                end do
             end do
          end do
       end do
       ! calculate the reduction factor of Fourier modes
       ! used to represent J0
       call sum_allreduce (ia_max_aj0a_count)
       ia_max_aj0a_reduction_factor = real(ia_max_aj0a_count)/real(naky*nakx*nztot*nmu*nvpa*nspec*naky)

       if (proc0) then
          write (*,*) 'average number of k-alphas needed to represent J0(kperp(alpha))=', ia_max_aj0a_reduction_factor*naky, 'out of ', naky
          write (*,*)
       end if

       deallocate (aj0_kalpha, aj0_alpha)
    else
       ia = 1
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          do iz = -nzgrid, nzgrid
             do ikx = 1, nakx
                do iky = 1, naky
                   arg = spec(is)%smz*sqrt(vperp2(ia,iz,imu)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
                   aj0x(iky,ikx,iz,ivmu) = j0(arg)
                end do
             end do
          end do
       end do
    end if

  end subroutine init_bessel

  subroutine find_max_required_kalpha_index (ft, imu, iz, idx)

    use vpamu_grids, only: maxwell_mu

    implicit none

    complex, dimension (:), intent (in) :: ft
    integer, intent (in) :: imu, iz
    integer, intent (out) :: idx

    real, parameter :: tol_floor = 0.01
    integer :: i, n
    real :: subtotal, total
    real :: tol
    real, dimension (:), allocatable :: ftmod2

    n = size(ft)

    ! use conservative estimate
    ! when deciding number of modes to retain
    tol = min(0.1,tol_floor/maxval(maxwell_mu(:,iz,imu)))

    allocate (ftmod2(n))
    ! get spectral energy associated with each mode
    ftmod2 = real(ft*conjg(ft))
    ! get total spectral energy
    total = sum(ftmod2)
    subtotal = 0.

    ! find minimum spectral index for which
    ! desired percentage of spectral energy contained
    ! in modes with indices at or below it
    i = 1
    do while (subtotal < total*(1.0-tol))
       idx = i
       subtotal = sum(ftmod2(:i))
       i = i + 1
    end do
    
    deallocate (ftmod2)

  end subroutine find_max_required_kalpha_index

  subroutine finish_bessel

    implicit none

    if (allocated(aj0v)) deallocate (aj0v)
    if (allocated(aj0x)) deallocate (aj0x)
    if (allocated(aj0a)) deallocate (aj0a)
    if (allocated(ia_max_aj0a)) deallocate (ia_max_aj0a)

    bessinit = .false.

  end subroutine finish_bessel

  subroutine gyro_average_kxky_local (field, iz, ivmu, gyro_field)

    implicit none

    complex, dimension (:,:), intent (in) :: field
    integer, intent (in) :: iz, ivmu
    complex, dimension (:,:), intent (out) :: gyro_field

    gyro_field = aj0x(:,:,iz,ivmu)*field

  end subroutine gyro_average_kxky_local

  subroutine gyro_average_kxkyz_local (field, ivmu, gyro_field)

    use zgrid, only: nzgrid, ntubes

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: field
    integer, intent (in) :: ivmu
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: gyro_field

    gyro_field = spread(aj0x(:,:,:,ivmu),4,ntubes)*field

  end subroutine gyro_average_kxkyz_local

  subroutine gyro_average_vmu_local (distfn, ikxkyz, gyro_distfn)

    use vpamu_grids, only: nvpa

    implicit none

    complex, dimension (:,:), intent (in) :: distfn
    integer, intent (in) :: ikxkyz
    complex, dimension (:,:), intent (out) :: gyro_distfn

    gyro_distfn = spread(aj0v(:,ikxkyz),1,nvpa)*distfn

  end subroutine gyro_average_vmu_local

  subroutine gyro_average_j1 (distfn, ikxkyz, gyro_distfn)

    use vpamu_grids, only: nvpa

    implicit none

    complex, dimension (:,:), intent (in) :: distfn
    integer, intent (in) :: ikxkyz
    complex, dimension (:,:), intent (out) :: gyro_distfn

    gyro_distfn = spread(aj1v(:,ikxkyz),1,nvpa)*distfn

  end subroutine gyro_average_j1

  subroutine gyro_average_vmus_nonlocal (field, iky, ikx, iz, gyro_field)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (vmu_lo%llim_proc:), intent (in) :: field
    integer, intent (in) :: iky, ikx, iz
    complex, dimension (vmu_lo%llim_proc:), intent (out) :: gyro_field

    gyro_field = aj0x(iky,ikx,iz,:)*field

  end subroutine gyro_average_vmus_nonlocal

end module gyro_averages
