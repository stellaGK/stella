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

  logical :: bessinit = .false.

contains

  subroutine init_bessel

    use dist_fn_arrays, only: kperp2
    use species, only: spec, nspec
    use stella_geometry, only: bmag, nalpha
    use zgrid, only: nzgrid
    use vpamu_grids, only: vperp2, nmu
    use kt_grids, only: naky, nakx
    use stella_layouts, only: kxkyz_lo, vmu_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, imu_idx
    use spfunc, only: j0, j1

    implicit none

    integer :: iz, iky, ikx, imu, is, ia
    integer :: ikxkyz, ivmu
    real :: arg

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

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       do iz = -nzgrid, nzgrid
          do ia = 1, nalpha
             do ikx = 1, nakx
                do iky = 1, naky
                   arg = spec(is)%smz*sqrt(vperp2(ia,iz,imu)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
                   aj0x(iky,ikx,iz,ivmu) = j0(arg)
                end do
             end do
          end do
       end do
    end do

  end subroutine init_bessel

  subroutine finish_bessel

    implicit none

    if (allocated(aj0v)) deallocate (aj0v)
    if (allocated(aj0x)) deallocate (aj0x)

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

    use zgrid, only: nzgrid

    implicit none

    complex, dimension (:,:,-nzgrid:), intent (in) :: field
    integer, intent (in) :: ivmu
    complex, dimension (:,:,-nzgrid:), intent (out) :: gyro_field

    gyro_field = aj0x(:,:,:,ivmu)*field

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
