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
  ! (naky, nakx, -ntgrid:ntgrid, -gxyz-layout-)

  complex, dimension (:,:,:,:), allocatable :: g1, g2
  ! (naky, nakx, -ntgrid:ntgrid, -gxyz-layout-)

  complex, dimension (:,:,:), allocatable :: gvmu
  ! (-nvgrid:nvgrid, nmu, nspec, -gvmu-layout-)

  real, dimension (:,:), allocatable :: wstar, wdriftx, wdrifty
  ! (-ntgrid:ntgrid, -gxyz-layout-)

  real, dimension (:,:,:,:), allocatable :: aj0x
  ! (naky, nakx, -ntgrid:ntgrid, -gxyz-layout-)

  real, dimension (:,:), allocatable :: aj0v
  ! (nmu, -gvmus-layout-)

  interface g_adjust
     module procedure g_adjust_vmu
     module procedure g_adjust_xyz
  end interface

  private

contains

  subroutine g_adjust_xyz (g, phi, apar, facphi, facapar)

    use species, only: spec
    use theta_grid, only: ntgrid
    use vpamu_grids, only: anon, vpa
    use stella_layouts, only: gxyz_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use kt_grids, only: naky, nakx

    implicit none
    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-ntgrid:), intent (in) :: phi, apar
    real, intent (in) :: facphi, facapar

    integer :: ixyz, ig, iky, ikx, is, imu, iv
    complex :: adj

    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       iv = iv_idx(gxyz_lo,ixyz)
       imu = imu_idx(gxyz_lo,ixyz)
       is = is_idx(gxyz_lo,ixyz)
       do ig = -ntgrid, ntgrid
          do ikx = 1, nakx
             do iky = 1, naky
                adj = aj0x(iky,ikx,ig,ixyz)*spec(is)%zt*anon(ig,iv,imu) &
                     * ( facphi*phi(iky,ikx,ig) - facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,ig) )
                g(iky,ikx,ig,ixyz) = g(iky,ikx,ig,ixyz) + adj
             end do
          end do
       end do
    end do
  end subroutine g_adjust_xyz

  subroutine g_adjust_vmu (g, phi, apar, facphi, facapar)

    use species, only: spec
    use theta_grid, only: ntgrid
    use vpamu_grids, only: anon, vpa
    use vpamu_grids, only: nvgrid, nmu
    use stella_layouts, only: gvmu_lo
    use stella_layouts, only: iky_idx, ikx_idx, ig_idx, is_idx

    implicit none
    complex, dimension (-nvgrid:,:,gvmu_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-ntgrid:), intent (in) :: phi, apar
    real, intent (in) :: facphi, facapar

    integer :: ivmu, ig, iky, ikx, is, imu, iv
    complex :: adj

    do ivmu = gvmu_lo%llim_proc, gvmu_lo%ulim_proc
       ig = ig_idx(gvmu_lo,ivmu)
       ikx = ikx_idx(gvmu_lo,ivmu)
       iky = iky_idx(gvmu_lo,ivmu)
       is = is_idx(gvmu_lo,ivmu)
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             adj = aj0v(imu,ivmu)*spec(is)%zt*anon(ig,iv,imu) &
                  * ( facphi*phi(iky,ikx,ig) - facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,ig) )
             g(iv,imu,ivmu) = g(iv,imu,ivmu) + adj
          end do
       end do
    end do
  end subroutine g_adjust_vmu

end module dist_fn_arrays
