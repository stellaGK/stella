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
  ! (naky, ntheta0, -ntgrid:ntgrid, -gxyz-layout-)

  complex, dimension (:,:,:,:), allocatable :: g1, g2
  ! (naky, ntheta0, -ntgrid:ntgrid, -gxyz-layout-)

  complex, dimension (:,:,:), allocatable :: gvmu
  ! (-nvgrid:nvgrid, nmu, nspec, -gvmu-layout-)

  real, dimension (:,:), allocatable :: wstar, wdriftx, wdrifty
  ! (-ntgrid:ntgrid, -gxyz-layout-)

  real, dimension (:,:,:,:), allocatable :: aj0x
  ! (naky, ntheta0, -ntgrid:ntgrid, -gxyz-layout-)

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
    use kt_grids, only: naky, ntheta0

    implicit none
    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-ntgrid:), intent (in) :: phi, apar
    real, intent (in) :: facphi, facapar

    integer :: ixyz, ig, ik, it, is, imu, iv
    complex :: adj

    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       iv = iv_idx(gxyz_lo,ixyz)
       imu = imu_idx(gxyz_lo,ixyz)
       is = is_idx(gxyz_lo,ixyz)
       do ig = -ntgrid, ntgrid
          do it = 1, ntheta0
             do ik = 1, naky
                adj = aj0x(ik,it,ig,ixyz)*spec(is)%zt*anon(ig,iv,imu) &
                     * ( facphi*phi(ik,it,ig) - facapar*vpa(iv)*spec(is)%stm*apar(ik,it,ig) )
                g(ik,it,ig,ixyz) = g(ik,it,ig,ixyz) + adj
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
    use stella_layouts, only: ik_idx, it_idx, ig_idx, is_idx

    implicit none
    complex, dimension (-nvgrid:,:,gvmu_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-ntgrid:), intent (in) :: phi, apar
    real, intent (in) :: facphi, facapar

    integer :: ivmu, ig, ik, it, is, imu, iv
    complex :: adj

    do ivmu = gvmu_lo%llim_proc, gvmu_lo%ulim_proc
       ig = ig_idx(gvmu_lo,ivmu)
       it = it_idx(gvmu_lo,ivmu)
       ik = ik_idx(gvmu_lo,ivmu)
       is = is_idx(gvmu_lo,ivmu)
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             adj = aj0v(imu,ivmu)*spec(is)%zt*anon(ig,iv,imu) &
                  * ( facphi*phi(ik,it,ig) - facapar*vpa(iv)*spec(is)%stm*apar(ik,it,ig) )
             g(iv,imu,ivmu) = g(iv,imu,ivmu) + adj
          end do
       end do
    end do
  end subroutine g_adjust_vmu

end module dist_fn_arrays
