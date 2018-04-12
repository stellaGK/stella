module neoclassical_terms

  implicit none

  public :: init_neoclassical_terms
  public :: finish_neoclassical_terms
  public :: include_neoclassical_terms
  public :: dfneo_dzed, dfneo_dvpa, dfneo_drho
  public :: dphineo_dzed, dphineo_drho

  private

  logical :: include_neoclassical_terms
  integer :: nradii
  real :: drho

  integer :: neo_option_switch
  integer, parameter :: neo_option_sfincs = 1

  real, dimension (:,:), allocatable :: dfneo_dzed, dfneo_dvpa, dfneo_drho
  real, dimension (:), allocatable :: dphineo_dzed, dphineo_drho

  logical :: neoinit = .false.
  logical :: debug = .false.

contains

  subroutine init_neoclassical_terms

    use zgrid, only: nzgrid
    use vpamu_grids, only: nvpa, nmu
    use species, only: nspec
    use stella_layouts, only: vmu_lo
    use sfincs_interface, only: get_neo_from_sfincs
    
    implicit none

    real, dimension (:,:,:,:,:), allocatable :: f_neoclassical
    real, dimension (:,:), allocatable :: phi_neoclassical
    
    if (neoinit) return
    neoinit = .true.

    call read_parameters
    if (include_neoclassical_terms) then
       allocate (f_neoclassical(-nzgrid:nzgrid,nvpa,nmu,nspec,-nradii/2:nradii/2))
       allocate (phi_neoclassical(-nzgrid:nzgrid,-nradii/2:nradii/2))
       if (.not.allocated(dfneo_dvpa)) &
            allocate (dfneo_dvpa(-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       if (.not.allocated(dfneo_drho)) &
            allocate (dfneo_drho(-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       if (.not.allocated(dfneo_dzed)) &
            allocate (dfneo_dzed(-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       if (.not.allocated(dphineo_dzed)) &
            allocate (dphineo_dzed(-nzgrid:nzgrid))
       if (.not.allocated(dphineo_drho)) &
            allocate (dphineo_drho(-nzgrid:nzgrid))
       select case (neo_option_switch)
       case (neo_option_sfincs)
          call get_neo_from_sfincs (nradii, drho, f_neoclassical, phi_neoclassical)
       end select
       if (debug) write (6,*) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_dzed'
       call get_dfneo_dzed (f_neoclassical(:,:,:,:,0), dfneo_dzed)
       if (debug) write (6,*) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_dvpa'
       call get_dfneo_dvpa (f_neoclassical(:,:,:,:,0), dfneo_dvpa)
       if (debug) write (6,*) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_drho'
       call get_dfneo_drho (f_neoclassical, dfneo_drho)
       if (debug) write (6,*) 'neoclassical_terms::init_neoclassical_terms::get_dphineo_dzed'
       call get_dphineo_dzed (phi_neoclassical(:,0), dphineo_dzed)
       if (debug) write (6,*) 'neoclassical_terms::init_neoclassical_terms::get_dphineo_drho'
       call get_dphineo_drho (phi_neoclassical, dphineo_drho)
       if (debug) write (6,*) 'neoclassical_terms::init_neoclassical_terms::write_neoclassical'
       call write_neoclassical (f_neoclassical, phi_neoclassical)
       deallocate (f_neoclassical, phi_neoclassical)
    end if

  end subroutine init_neoclassical_terms

  subroutine read_parameters

    use mp, only: proc0, broadcast
    use file_utils, only: error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value

    implicit none

    type (text_option), dimension (2), parameter :: neoopts = (/ &
         text_option('default', neo_option_sfincs), &
         text_option('sfincs', neo_option_sfincs) /)
    character (10) :: neo_option

    namelist /neoclassical_input/ include_neoclassical_terms, &
         neo_option, nradii, drho

    logical :: exist
    integer :: ierr, in_file

    if (proc0) then
       ! set to .true. to include neoclassical terms in GK equation
       include_neoclassical_terms = .false.
       ! number of radial points used for radial derivatives
       ! of neoclassical quantities
       nradii = 5
       ! spacing in rhoc between radial points used for radial derivatives
       drho = 0.01
       ! option for obtaining neoclassical distribution function and potential
       neo_option = 'sfincs'

       in_file = input_unit_exist("neoclassical_input", exist)
       if (exist) read (unit=in_file, nml=neoclassical_input)

       ierr = error_unit()
       call get_option_value &
            (neo_option, neoopts, neo_option_switch, &
            ierr, "neo_option in neoclassical_input")

       if (nradii /= 3 .and. nradii /= 5) then
          write (*,*) 'WARNING: only nradii of 3 or 5 is currently supported in neoclassical_input namelist'
          write (*,*) 'WARNING: forcing nradii=5'
          nradii = 5
       end if
    end if

    call broadcast (include_neoclassical_terms)
    call broadcast (neo_option_switch)
    call broadcast (nradii)
    call broadcast (drho)

  end subroutine read_parameters

  subroutine distribute_vmus_over_procs (local, distributed)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx

    implicit none

    real, dimension (:,:,:), intent (in) :: local
    real, dimension (vmu_lo%llim_proc:), intent (out) :: distributed

    integer :: ivmu, iv, imu, is
    
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       distributed(ivmu) = local(iv,imu,is)
    end do

  end subroutine distribute_vmus_over_procs

  subroutine get_dfneo_dvpa (fneo, dfneo)

    use finite_differences, only: fd5pt
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use vpamu_grids, only: nvpa, nmu, nvpa
    use vpamu_grids, only: dvpa
    use species, only: nspec

    implicit none

    real, dimension (-nzgrid:,:,:,:), intent (in) :: fneo
    real, dimension (-nzgrid:,vmu_lo%llim_proc:), intent (out) :: dfneo

    integer :: iz, imu, is
    real, dimension (:), allocatable :: tmp1, tmp2
    real, dimension (:,:,:,:), allocatable :: dfneo_local

    allocate (tmp1(nvpa), tmp2(nvpa))
    allocate (dfneo_local(-nzgrid:nzgrid,nvpa,nmu,nspec))

    do is = 1, nspec
       do imu = 1, nmu
          do iz = -nzgrid, nzgrid
             ! hack to avoid dealing with negative indices in fd5pt
             tmp1 = fneo(iz,:,imu,is)
             call fd5pt (tmp1, tmp2, dvpa)
             dfneo_local(iz,:,imu,is) = tmp2
          end do
       end do
    end do

    do iz = -nzgrid, nzgrid
       call distribute_vmus_over_procs (dfneo_local(iz,:,:,:), dfneo(iz,:))
    end do

    deallocate (dfneo_local)
    deallocate (tmp1, tmp2)

  end subroutine get_dfneo_dvpa

  subroutine get_dfneo_dzed (fneo, dfneo)

    use finite_differences, only: fd5pt
    use zgrid, only: nztot, nzgrid, delzed
    use vpamu_grids, only: nvpa, nmu
    use species, only: nspec
    use stella_layouts, only: vmu_lo

    implicit none

    real, dimension (-nzgrid:,:,:,:), intent (in) :: fneo
    real, dimension (-nzgrid:,vmu_lo%llim_proc:), intent (out) :: dfneo

    integer :: iv, imu, is, iz
    real, dimension (:), allocatable :: tmp1, tmp2
    real, dimension (:), allocatable :: dfneo_local(:,:,:,:)

    allocate (tmp1(nztot), tmp2(nztot))
    allocate (dfneo_local(-nzgrid:nzgrid,nvpa,nmu,nspec))

    do is = 1, nspec
       do imu = 1, nmu
          do iv = 1, nvpa
             ! hack to avoid dealing with negative indices in fd5pt
!             ! fneo is F_nc * exp(2*mu*B) * ...
!             ! need to get rid of z-dependent exponential before
!             ! taking d/dz
             tmp1 = fneo(:,iv,imu,is)!*exp(-2.0*mu(imu)*bmag(1,:))
             call fd5pt (tmp1, tmp2, delzed(0))
             ! put the z-dependent exponential normalization factor back
             dfneo_local(:,iv,imu,is) = tmp2!*exp(2.0*mu(imu)*bmag(1,:))
          end do
       end do
    end do

    do iz = -nzgrid, nzgrid
       call distribute_vmus_over_procs (dfneo_local(iz,:,:,:), dfneo(iz,:))
    end do

    deallocate (dfneo_local)
    deallocate (tmp1, tmp2)

  end subroutine get_dfneo_dzed

  subroutine get_dfneo_drho (fneo, dfneo)

    use finite_differences, only: fd3pt, fd5pt
    use zgrid, only: nzgrid
    use vpamu_grids, only: nvpa, nmu
    use species, only: nspec
    use stella_layouts, only: vmu_lo

    implicit none

    real, dimension (-nzgrid:,:,:,:,-nradii/2:), intent (in) :: fneo
    real, dimension (-nzgrid:,vmu_lo%llim_proc:), intent (out) :: dfneo

    integer :: iz, iv, imu, is
    real, dimension (:), allocatable :: tmp1, tmp2
    real, dimension (:,:,:,:), allocatable :: dfneo_local

    allocate (tmp1(nradii), tmp2(nradii))
    allocate (dfneo_local(-nzgrid:nzgrid,nvpa,nmu,nspec))

    do is = 1, nspec
       do imu = 1, nmu
          do iv = 1, nvpa
             do iz = -nzgrid, nzgrid
                ! hack to avoid dealing with negative indices in fd5pt
                tmp1 = fneo(iz,iv,imu,is,:)
                if (nradii == 5) then
                   call fd5pt (tmp1, tmp2, drho)
                else
                   call fd3pt (tmp1, tmp2, drho)
                end if
                dfneo_local(iz,iv,imu,is) = tmp2(nradii/2+1)
             end do
          end do
       end do
    end do

    do iz = -nzgrid, nzgrid
       call distribute_vmus_over_procs (dfneo_local(iz,:,:,:), dfneo(iz,:))
    end do

    deallocate (dfneo_local)
    deallocate (tmp1, tmp2)

  end subroutine get_dfneo_drho

  subroutine get_dphineo_dzed (phineo, dphineo)

    use finite_differences, only: fd5pt
    use zgrid, only: nztot, nzgrid, delzed

    implicit none

    real, dimension (-nzgrid:), intent (in) :: phineo
    real, dimension (-nzgrid:), intent (out) :: dphineo

    real, dimension (:), allocatable :: tmp1, tmp2

    allocate (tmp1(nztot), tmp2(nztot))

    ! hack to avoid dealing with negative indices in fd5pt
    tmp1 = phineo
    call fd5pt (tmp1, tmp2, delzed(0))
    dphineo = tmp2

    deallocate (tmp1, tmp2)

  end subroutine get_dphineo_dzed

  subroutine get_dphineo_drho (phineo, dphineo)

    use finite_differences, only: fd3pt, fd5pt
    use zgrid, only: nzgrid

    implicit none

    real, dimension (-nzgrid:,-nradii/2:), intent (in) :: phineo
    real, dimension (-nzgrid:), intent (out) :: dphineo

    integer :: iz
    real, dimension (:), allocatable :: tmp1, tmp2

    allocate (tmp1(nradii), tmp2(nradii))

    do iz = -nzgrid, nzgrid
       ! hack to avoid dealing with negative indices in fd5pt
       tmp1 = phineo(iz,:)
       if (nradii == 5) then
          call fd5pt (tmp1, tmp2, drho)
       else
          call fd3pt (tmp1, tmp2, drho)
       end if
       dphineo(iz) = tmp2(nradii/2+1)
    end do

    deallocate (tmp1, tmp2)

  end subroutine get_dphineo_drho

  subroutine write_neoclassical (fnc, phinc)

    use mp, only: proc0
    use mp, only: send, receive
    use file_utils, only: open_output_file, close_output_file
    use species, only: nspec
    use zgrid, only: nzgrid, zed
    use vpamu_grids, only: vpa, mu, nvpa
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_layouts, only: idx_local, proc_id

    implicit none

    real, dimension (-nzgrid:,:,:,:,-nradii/2:), intent (in) :: fnc
    real, dimension (-nzgrid:,-nradii/2:), intent (in) :: phinc
    
    integer :: neo_unit
    integer :: irad, iz, ivmu, iv, imu, is
    real, dimension (:), allocatable :: dfdv_local, dfdr_local, dfdz_local

    allocate (dfdv_local(-nzgrid:nzgrid))
    allocate (dfdr_local(-nzgrid:nzgrid))
    allocate (dfdz_local(-nzgrid:nzgrid))

    if (proc0) then
       call open_output_file (neo_unit,'.neoclassical')
       write (neo_unit,'(2a8,10a13)') '#1.rad', '2.spec', '3.zed', '4.mu', '5.vpa', &
            '6.f_neo', '7.dfdvpa_neo', '8.dfdrho_neo', '9.dfdzed_neo', '10.phi_neo', &
            '11.dphidrho', '12.dphidzed'
    end if
    do irad = -nradii/2, nradii/2
       do ivmu = vmu_lo%llim_world, vmu_lo%ulim_world
          iv = iv_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          is = is_idx(vmu_lo,ivmu)
!          ! TMP FOR TESTING -- MAB
!          if (is /=1) cycle
          if (idx_local(vmu_lo, iv, imu, is)) then
             if (proc0) then
                dfdv_local = dfneo_dvpa(:,ivmu)
                dfdr_local = dfneo_drho(:,ivmu)
                dfdz_local = dfneo_dzed(:,ivmu)
             else
                call send (dfneo_dvpa(:,ivmu), 0)
                call send (dfneo_drho(:,ivmu), 0)
                call send (dfneo_dzed(:,ivmu), 0)
             end if
          else if (proc0) then
             call receive (dfdv_local, proc_id(vmu_lo,ivmu))
             call receive (dfdr_local, proc_id(vmu_lo,ivmu))
             call receive (dfdz_local, proc_id(vmu_lo,ivmu))
          end if
          if (proc0) then
!             ! TMP FOR TESTING -- MAB
!             iz = 0
             do iz = -nzgrid, nzgrid
                write (neo_unit,'(2i8,10e13.5)') irad, is, zed(iz), mu(imu), vpa(iv), &
                     fnc(iz,iv,imu,is,irad), &!*exp(-2.0*mu(imu)*bmag(1,iz)), &
                     dfdv_local(iz), &!*exp(-2.0*mu(imu)*bmag(1,iz)), &
                     dfdr_local(iz), &!*exp(-2.0*mu(imu)*bmag(1,iz)), &
                     dfdz_local(iz), &!*exp(-2.0*mu(imu)*bmag(1,iz)), &
                     phinc(iz,irad), dphineo_drho(iz), dphineo_dzed(iz)
             end do
          end if
!          ! TMP FOR TESTING -- MAB
!          if (iv == nvpa) write (neo_unit,*)
       end do
       if (proc0) write (neo_unit,*)
    end do
    if (proc0) call close_output_file (neo_unit)

    deallocate (dfdv_local, dfdr_local, dfdz_local)

  end subroutine write_neoclassical

  subroutine finish_neoclassical_terms

    implicit none

    if (allocated(dfneo_dvpa)) deallocate (dfneo_dvpa)
    if (allocated(dfneo_drho)) deallocate (dfneo_drho)
    if (allocated(dfneo_dzed)) deallocate (dfneo_dzed)
    if (allocated(dphineo_dzed)) deallocate (dphineo_dzed)
    if (allocated(dphineo_drho)) deallocate (dphineo_drho)

    neoinit = .false.

  end subroutine finish_neoclassical_terms

end module neoclassical_terms
