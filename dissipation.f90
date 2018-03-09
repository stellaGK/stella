module dissipation

  implicit none

  public :: init_dissipation
  public :: include_collisions
  public :: advance_collisions
  public :: time_collisions
  public :: hyper_dissipation
  public :: advance_hyper_dissipation

  private

  logical :: include_collisions
  logical :: conserve_moments
  logical :: hyper_dissipation
  real :: D_hyper

  real, dimension (2,2) :: time_collisions

contains

  subroutine init_dissipation

    implicit none

    call read_parameters

  end subroutine init_dissipation

  subroutine read_parameters

    use file_utils, only: input_unit_exist
    use mp, only: proc0, broadcast

    implicit none

    namelist /dissipation/ hyper_dissipation, D_hyper, &
         include_collisions, conserve_moments

    integer :: in_file
    logical :: dexist

    if (proc0) then
       include_collisions = .false.
       conserve_moments = .false.
       hyper_dissipation = .false.
       D_hyper = 0.05

       in_file = input_unit_exist("dissipation", dexist)
       if (dexist) read (unit=in_file, nml=dissipation)
    end if

    call broadcast (include_collisions)
    call broadcast (conserve_moments)
    call broadcast (hyper_dissipation)
    call broadcast (D_hyper)

  end subroutine read_parameters

  subroutine advance_collisions (g, phi, gke_rhs)

    use mp, only: proc0
    use job_manage, only: time_message
    use redistribute, only: scatter, gather
    use stella_time, only: code_dt
    use zgrid, only: nzgrid
    use species, only: spec
    use run_parameters, only: fphi
    use kt_grids, only: naky, nakx
    use vpamu_grids, only: nvpa, nmu
    use geometry, only: bmag
    use stella_layouts, only: vmu_lo, kxkyz_lo
    use stella_layouts, only: is_idx, iky_idx, ikx_idx, iz_idx
    use dist_redistribute, only: kxkyz2vmu
    use dist_fn_arrays, only: gvmu, g_to_h, kperp2

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:), intent (in) :: phi
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: gke_rhs

    integer :: is, ikxkyz, imu, iv, ivmu, ikx, iky, iz
    complex, dimension (:), allocatable :: mucoll
    complex, dimension (:,:,:), allocatable :: coll
    complex, dimension (:,:,:,:), allocatable :: tmp_vmulo

    if (proc0) call time_message(.false.,time_collisions(:,1),' collisions')

    allocate (tmp_vmulo(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    ! switch from g = <f> to h = f + Z*e*phi/T * F0
    tmp_vmulo = g
    call g_to_h (tmp_vmulo, phi, fphi)

    ! remap so that (vpa,mu) local
    if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')
    call scatter (kxkyz2vmu, tmp_vmulo, gvmu)
    if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')

    deallocate (tmp_vmulo)
    
    allocate (coll(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    allocate (mucoll(nmu))

    ! take vpa derivatives
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          call vpa_differential_operator (gvmu(:,imu,ikxkyz), coll(:,imu,ikxkyz)) 
       end do      
       do iv = 1, nvpa
          call mu_differential_operator (iz, gvmu(iv,:,ikxkyz), mucoll)
          coll(iv,:,ikxkyz) = coll(iv,:,ikxkyz) + mucoll
       end do
       if (conserve_moments) then
          call conserve_momentum (iky, ikx, iz, is, ikxkyz, gvmu(:,:,ikxkyz), coll(:,:,ikxkyz))
          call conserve_energy (iz, ikxkyz, gvmu(:,:,ikxkyz), coll(:,:,ikxkyz))
       end if
       ! save memory by using gvmu and deallocating coll below
       ! before re-allocating tmp_vmulo
       gvmu(:,:,ikxkyz) = coll(:,:,ikxkyz) - 0.5*kperp2(iky,ikx,iz)*(spec(is)%smz/bmag(1,iz))**2*gvmu(:,:,ikxkyz)
    end do
    deallocate (coll, mucoll)
    allocate (tmp_vmulo(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    ! remap so that (ky,kx,z) local
    call gather (kxkyz2vmu, gvmu, tmp_vmulo)

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       gke_rhs(:,:,:,ivmu) =  gke_rhs(:,:,:,ivmu) + code_dt*spec(is)%vnew(is)*tmp_vmulo(:,:,:,ivmu)
    end do

    deallocate (tmp_vmulo)

    if (proc0) call time_message(.false.,time_collisions(:,1),' collisions')

  end subroutine advance_collisions

  subroutine vpa_differential_operator (h, Dh)

    use vpamu_grids, only: nvpa, vpa, dvpa

    implicit none

    complex, dimension (:), intent (in) :: h
    complex, dimension (:), intent (out) :: Dh

    integer :: iv

    ! use h = 0 at ghost cells beyond +/- vpa_max
    iv = 1
    Dh(iv) = (0.5*h(iv+1)*(1.0+vpa(iv+1))-h(iv))/dvpa
    iv = nvpa
    Dh(iv) = (-h(iv)+0.5*h(iv-1)*(1.0-vpa(iv-1)))/dvpa
    do iv = 2, nvpa-1
       Dh(iv) = (0.5*h(iv+1)*(1.0+vpa(iv+1))-h(iv)+0.5*h(iv-1)*(1.0-vpa(iv-1)))/dvpa
    end do

  end subroutine vpa_differential_operator

  subroutine mu_differential_operator (iz, h, Dh)

    use vpamu_grids, only: nmu, mu, dmu
    use finite_differences, only: d2_3pt, fd3pt
    use geometry, only: bmag

    implicit none

    integer, intent (in) :: iz
    complex, dimension (:), intent (in) :: h
    complex, dimension (:), intent (out) :: Dh

    integer :: imu
    complex, dimension (:), allocatable :: h_ghost, Dh_ghost
    real, dimension (:), allocatable :: dmu_ghost

    allocate (h_ghost(nmu+1))
    allocate (Dh_ghost(nmu+1))
    allocate (dmu_ghost(nmu))
    ! pad h_ghost array with ghost cell beyond max(mu) with zero BC
    h_ghost(:nmu) = h ; h_ghost(nmu+1) = 0.
    ! assign extra dmu value at nmu (beyond mu grid)
    ! because it will be accessed (but not later used) 
    ! by generic subroutine d2_3pt
    dmu_ghost(:nmu-1) = dmu(:nmu-1) ; dmu_ghost(nmu) = 1.0

    call d2_3pt (h_ghost, Dh_ghost, dmu_ghost)
    Dh = Dh_ghost(:nmu)*mu/bmag(1,iz)

    ! next add (1/B + 2*mu)*dh/dmu + 2*h
    call fd3pt (h_ghost, Dh_ghost, dmu_ghost)
    Dh = Dh + (1./bmag(1,iz) + 2.*mu)*Dh_ghost(:nmu) + 2.*h

    deallocate (h_ghost, Dh_ghost, dmu_ghost)

  end subroutine mu_differential_operator

  subroutine conserve_momentum (iky, ikx, iz, is, ikxkyz, h, Ch)

    use species, only: spec
    use geometry, only: bmag
    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: vpa, nvpa, nmu, vperp2
    use vpamu_grids, only: maxwell_vpa, maxwell_mu
    use dist_fn_arrays, only: aj0v, aj1v, kperp2
    
    implicit none

    integer, intent (in) :: iky, ikx, iz, is, ikxkyz
    complex, dimension (:,:), intent (in) :: h
    complex, dimension (:,:), intent (in out) :: Ch

    complex, dimension (:,:), allocatable :: u_fac
    complex :: integral

    allocate (u_fac(nvpa,nmu))

    u_fac = spread(aj0v(:,ikxkyz),1,nvpa)*spread(vpa,2,nmu)
    call integrate_vmu (u_fac*h,iz,integral)

    Ch = Ch + 2.0*u_fac*integral*spread(maxwell_mu(1,iz,:),1,nvpa)*spread(maxwell_vpa,2,nmu)

    u_fac = spread(vperp2(1,iz,:)*aj1v(:,ikxkyz),1,nvpa)*sqrt(kperp2(iky,ikx,iz))*spec(is)%smz/bmag(1,iz)
    call integrate_vmu (u_fac*h,iz,integral)
    
    Ch = Ch + 2.0*u_fac*integral*spread(maxwell_mu(1,iz,:),1,nvpa)*spread(maxwell_vpa,2,nmu)

    deallocate (u_fac)

  end subroutine conserve_momentum

  subroutine conserve_energy (iz, ikxkyz, h, Ch)

    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: vpa, nvpa, nmu, vperp2
    use vpamu_grids, only: maxwell_vpa, maxwell_mu
    use dist_fn_arrays, only: aj0v
    
    implicit none

    integer, intent (in) :: iz, ikxkyz
    complex, dimension (:,:), intent (in) :: h
    complex, dimension (:,:), intent (in out) :: Ch

    complex, dimension (:,:), allocatable :: T_fac
    complex :: integral

    allocate (T_fac(nvpa,nmu))

    T_fac = spread(aj0v(:,ikxkyz),1,nvpa)*(spread(vpa**2,2,nmu)+spread(vperp2(1,iz,:),1,nvpa)-1.5)
    call integrate_vmu (T_fac*h,iz,integral)

    Ch = Ch + 2.0*T_fac*integral*spread(maxwell_mu(1,iz,:),1,nvpa)*spread(maxwell_vpa,2,nmu)

    deallocate (T_fac)

  end subroutine conserve_energy

  subroutine advance_hyper_dissipation (g)

    use stella_time, only: code_dt
    use zgrid, only: nzgrid
    use stella_layouts, only: vmu_lo
    use dist_fn_arrays, only: kperp2

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: g

    integer :: ivmu
    real :: k2max

    k2max = maxval(kperp2)

    ! add in hyper-dissipation of form dg/dt = -D*(k/kmax)^4*g
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       g(:,:,:,ivmu) = g(:,:,:,ivmu)/(1.+code_dt*(kperp2/k2max)**2*D_hyper)
    end do

  end subroutine advance_hyper_dissipation

end module dissipation
