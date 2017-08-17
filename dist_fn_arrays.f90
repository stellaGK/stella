!> A container for the arrays that are used to store the distribution function among other things.
!!  These need to be accessible at a lower dependency level than the dist_fn module itself.
!! These arrays are allocated in the function dist_fn::allocate_arrays. 

module dist_fn_arrays

  public :: gnew, gold
  public :: g1, g2
  public :: gvmu
  public :: aj0x, aj0v
  public :: wstar
  public :: wdriftx, wdrifty
  public :: g_adjust

  ! dist fn
  complex, dimension (:,:,:,:), allocatable :: gnew, gold
  ! (naky, nakx, -ntgrid:ntgrid, -vmu-layout-)

  complex, dimension (:,:,:,:), allocatable :: g1, g2
  ! (naky, nakx, -ntgrid:ntgrid, -vmu-layout-)

  complex, dimension (:,:,:), allocatable :: gvmu
  ! (-nvgrid:nvgrid, nmu, nspec, -kxkyz-layout-)

  real, dimension (:,:), allocatable :: wstar
  ! (-ntgrid:ntgrid, -vmu-layout-)

  real, dimension (:,:,:), allocatable :: wdriftx, wdrifty
  ! (ny_ffs, -ntgrid:ntgrid, -vmu-layout-)

  real, dimension (:,:,:,:), allocatable :: aj0x
  ! (naky, nakx, -ntgrid:ntgrid, -vmu-layout-)

  real, dimension (:,:), allocatable :: aj0v
  ! (nmu, -kxkyz-layout-)

  interface g_adjust
     module procedure g_adjust_kxkyz
     module procedure g_adjust_vmu
  end interface

  private

contains

  subroutine g_adjust_vmu (g, phi, apar, facphi, facapar)

    use species, only: spec
    use theta_grid, only: ntgrid
    use vpamu_grids, only: anon, vpa
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use kt_grids, only: naky, nakx

    implicit none
    complex, dimension (:,:,-ntgrid:,vmu_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-ntgrid:), intent (in) :: phi, apar
    real, intent (in) :: facphi, facapar

    integer :: ivmu, ig, iky, ikx, is, imu, iv
    complex :: adj

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do ig = -ntgrid, ntgrid
          do ikx = 1, nakx
             do iky = 1, naky
                adj = aj0x(iky,ikx,ig,ivmu)*spec(is)%zt*anon(ig,iv,imu) &
                     * ( facphi*phi(iky,ikx,ig) - facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,ig) )
                g(iky,ikx,ig,ivmu) = g(iky,ikx,ig,ivmu) + adj
             end do
          end do
       end do
    end do
  end subroutine g_adjust_vmu

  subroutine g_adjust_kxkyz (g, phi, apar, facphi, facapar)

    use species, only: spec
    use theta_grid, only: ntgrid
    use vpamu_grids, only: anon, vpa
    use vpamu_grids, only: nvgrid, nmu
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, ig_idx, is_idx

    implicit none
    complex, dimension (-nvgrid:,:,kxkyz_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-ntgrid:), intent (in) :: phi, apar
    real, intent (in) :: facphi, facapar

    integer :: ikxkyz, ig, iky, ikx, is, imu, iv
    complex :: adj

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       ig = ig_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iky = iky_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             adj = aj0v(imu,ikxkyz)*spec(is)%zt*anon(ig,iv,imu) &
                  * ( facphi*phi(iky,ikx,ig) - facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,ig) )
             g(iv,imu,ikxkyz) = g(iv,imu,ikxkyz) + adj
          end do
       end do
    end do
  end subroutine g_adjust_kxkyz

end module dist_fn_arrays
