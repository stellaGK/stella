!> A container for the arrays that are used to store the distribution function among other things.
!!  These need to be accessible at a lower dependency level than the dist_fn module itself.
!! These arrays are allocated in the function dist_fn::allocate_arrays. 

module dist_fn_arrays

  public :: gnew, gold
  public :: g1, g2, g3
  public :: gvmu
  public :: aj0x
  public :: aj0v, aj1v
  public :: kperp2
  public :: wstar
  public :: wdriftx_g, wdrifty_g
  public :: wdriftx_phi, wdrifty_phi
  public :: gbar_to_g
  public :: gbar_to_h
  public :: gstar_to_g
  public :: g_to_h

  ! dist fn
  complex, dimension (:,:,:,:), allocatable :: gnew, gold
  ! (naky, nakx, -nzgrid:nzgrid, -vmu-layout-)

  complex, dimension (:,:,:,:), allocatable :: g1, g2, g3
  ! (naky, nakx, -nzgrid:nzgrid, -vmu-layout-)

  complex, dimension (:,:,:), allocatable :: gvmu
  ! (nvpa, nmu, nspec, -kxkyz-layout-)

  real, dimension (:,:,:), allocatable :: wstar
  ! (nalpha, -nzgrid:nzgrid, -vmu-layout-)

  real, dimension (:,:,:), allocatable :: wdriftx_g, wdrifty_g
  real, dimension (:,:,:), allocatable :: wdriftx_phi, wdrifty_phi
  ! (ny_ffs, -nzgrid:nzgrid, -vmu-layout-)

  real, dimension (:,:,:,:), allocatable :: aj0x
  ! (naky, nakx, -nzgrid:nzgrid, -vmu-layout-)

  real, dimension (:,:), allocatable :: aj0v, aj1v
  ! (nmu, -kxkyz-layout-)

  real, dimension (:,:,:), allocatable :: kperp2

  interface gbar_to_g
     module procedure gbar_to_g_kxkyz
     module procedure gbar_to_g_vmu
  end interface

  interface gbar_to_h
     module procedure gbar_to_h_kxkyz
     module procedure gbar_to_h_vmu
  end interface

!  interface gbar_to_hhat
!     module procedure gbar_to_hhat_kxkyz
!     module procedure gbar_to_hhat_vmu
!  end interface

  interface g_to_h
     module procedure g_to_h_kxkyz
     module procedure g_to_h_vmu
     module procedure g_to_h_vmu_zext
  end interface

  private

contains

  subroutine gbar_to_h_vmu (g, phi, apar, facphi, facapar)

    use species, only: spec
    use zgrid, only: nzgrid
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use kt_grids, only: naky, nakx

    implicit none
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-nzgrid:), intent (in) :: phi, apar
    real, intent (in) :: facphi, facapar

    integer :: ivmu, iz, iky, ikx, is, imu, iv
    complex :: adj

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do iz = -nzgrid, nzgrid
          do ikx = 1, nakx
             do iky = 1, naky
                adj = aj0x(iky,ikx,iz,ivmu)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
                     * ( facphi*phi(iky,ikx,iz) - facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,iz) )
                g(iky,ikx,iz,ivmu) = g(iky,ikx,iz,ivmu) + adj
             end do
          end do
       end do
    end do

  end subroutine gbar_to_h_vmu

  subroutine gbar_to_h_kxkyz (g, phi, apar, facphi, facapar)

    use species, only: spec
    use zgrid, only: nzgrid
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa
    use vpamu_grids, only: nvpa, nmu
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx

    implicit none
    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-nzgrid:), intent (in) :: phi, apar
    real, intent (in) :: facphi, facapar

    integer :: ikxkyz, iz, iky, ikx, is, imu, iv
    complex :: adj

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iz = iz_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iky = iky_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          do iv = 1, nvpa
             adj = aj0v(imu,ikxkyz)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
                  * ( facphi*phi(iky,ikx,iz) - facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,iz) )
             g(iv,imu,ikxkyz) = g(iv,imu,ikxkyz) + adj
          end do
       end do
    end do

  end subroutine gbar_to_h_kxkyz

  subroutine gbar_to_g_vmu (g, apar, facapar)

    use species, only: spec
    use zgrid, only: nzgrid
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use kt_grids, only: naky, nakx

    implicit none
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-nzgrid:), intent (in) :: apar
    real, intent (in) :: facapar

    integer :: ivmu, iz, iky, ikx, is, imu, iv
    complex :: adj

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do iz = -nzgrid, nzgrid
          do ikx = 1, nakx
             do iky = 1, naky
                adj = -aj0x(iky,ikx,iz,ivmu)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
                     * ( facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,iz) )
                g(iky,ikx,iz,ivmu) = g(iky,ikx,iz,ivmu) + adj
             end do
          end do
       end do
    end do

  end subroutine gbar_to_g_vmu

  subroutine gbar_to_g_kxkyz (g, apar, facapar)

    use species, only: spec
    use zgrid, only: nzgrid
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa
    use vpamu_grids, only: nvpa, nmu
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx

    implicit none
    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-nzgrid:), intent (in) :: apar
    real, intent (in) :: facapar

    integer :: ikxkyz, iz, iky, ikx, is, imu, iv
    complex :: adj

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iz = iz_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iky = iky_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          do iv = 1, nvpa
             adj = -aj0v(imu,ikxkyz)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
                  * ( facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,iz) )
             g(iv,imu,ikxkyz) = g(iv,imu,ikxkyz) + adj
          end do
       end do
    end do

  end subroutine gbar_to_g_kxkyz

  subroutine g_to_h_vmu (g, phi, facphi)

    use species, only: spec
    use zgrid, only: nzgrid
    use vpamu_grids, only: maxwell_vpa, maxwell_mu
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use kt_grids, only: naky, nakx

    implicit none
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-nzgrid:), intent (in) :: phi
    real, intent (in) :: facphi

    integer :: ivmu, iz, iky, ikx, is, imu, iv
    complex :: adj

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do iz = -nzgrid, nzgrid
          do ikx = 1, nakx
             do iky = 1, naky
                adj = aj0x(iky,ikx,iz,ivmu)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
                     * facphi*phi(iky,ikx,iz)
                g(iky,ikx,iz,ivmu) = g(iky,ikx,iz,ivmu) + adj
             end do
          end do
       end do
    end do

  end subroutine g_to_h_vmu

  subroutine g_to_h_vmu_zext (gext, phiext, facphi, iky, ie)

    use species, only: spec
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: nsegments
    use vpamu_grids, only: maxwell_vpa, maxwell_mu
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx

    implicit none
    complex, dimension (:,vmu_lo%llim_proc:), intent (in out) :: gext
    complex, dimension (:), intent (in) :: phiext
    real, intent (in) :: facphi
    integer, intent (in) :: iky, ie

    integer :: ivmu, iseg, iz, ikx, is, imu, iv
    integer :: idx
    complex :: adj

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)

       idx = 0
       iseg = 1
       ikx = ikxmod(iseg,ie,iky)
       do iz = iz_low(iseg), iz_up(iseg)
          idx = idx + 1
          adj = aj0x(iky,ikx,iz,ivmu)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
               * facphi*phiext(idx)
          gext(idx,ivmu) = gext(idx,ivmu) + adj
       end do
       if (nsegments(ie,iky) > 1) then
          do iseg = 2, nsegments(ie,iky)
             do iz = iz_low(iseg)+1, iz_up(iseg)
                adj = aj0x(iky,ikx,iz,ivmu)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
                     * facphi*phiext(idx)
                gext(idx,ivmu) = gext(idx,ivmu) + adj
                idx = idx + 1
             end do
          end do
       end if

    end do

  end subroutine g_to_h_vmu_zext

  subroutine g_to_h_kxkyz (g, phi, facphi)

    use species, only: spec
    use zgrid, only: nzgrid
    use vpamu_grids, only: maxwell_vpa, maxwell_mu
    use vpamu_grids, only: nvpa, nmu
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx

    implicit none
    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-nzgrid:), intent (in) :: phi
    real, intent (in) :: facphi

    integer :: ikxkyz, iz, iky, ikx, is, imu, iv
    complex :: adj

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iz = iz_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iky = iky_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          do iv = 1, nvpa
             adj = aj0v(imu,ikxkyz)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
                  * facphi*phi(iky,ikx,iz)
             g(iv,imu,ikxkyz) = g(iv,imu,ikxkyz) + adj
          end do
       end do
    end do

  end subroutine g_to_h_kxkyz

  subroutine gstar_to_g (g, phi, apar, facphi, facapar)

    use constants, only: zi
    use species, only: spec
    use zgrid, only: nzgrid
    use vpamu_grids, only: vpa
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, is_idx
    use kt_grids, only: naky, nakx
    use kt_grids, only: aky

    implicit none
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-nzgrid:), intent (in) :: phi, apar
    real, intent (in) :: facphi, facapar

    integer :: ivmu, iz, iky, ikx, is, iv
    complex :: adj

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do iz = -nzgrid, nzgrid
          do ikx = 1, nakx
             do iky = 1, naky
                ! BACKWARDS DIFFERENCE FLAG
!                adj = zi*aj0x(iky,ikx,iz,ivmu)*aky(iky)*wstar(iz,ivmu) &
                adj = zi*aj0x(iky,ikx,iz,ivmu)*aky(iky)*2.0*wstar(1,iz,ivmu) &
                     * ( facphi*phi(iky,ikx,iz) - facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,iz) )
                g(iky,ikx,iz,ivmu) = g(iky,ikx,iz,ivmu) + adj
             end do
          end do
       end do
    end do

  end subroutine gstar_to_g

end module dist_fn_arrays
