module dist_fn

  implicit none

  public :: init_gxyz
  public :: init_dist_fn, finish_dist_fn
  public :: adiabatic_option_switch
  public :: adiabatic_option_fieldlineavg

  private
  
  logical :: dist_fn_initialized = .false.
  logical :: gxyz_initialized = .false.
  logical :: kp2init = .false.
  logical :: vp2init = .false.
  logical :: bessinit = .false.
  logical :: readinit = .false.

  integer :: adiabatic_option_switch
  integer, parameter :: adiabatic_option_default = 1, &
       adiabatic_option_zero = 2, &
       adiabatic_option_fieldlineavg = 3, &
       adiabatic_option_yavg = 4

  logical :: debug = .false.

contains

  subroutine init_gxyz

    use dist_fn_arrays, only: gvmu, gold, gnew
    use redistribute, only: gather
    use dist_redistribute, only: kxkyz2vmu

    implicit none

    if (gxyz_initialized) return
    gxyz_initialized = .false.

    ! get version of g that has ky,kx,z local
    call gather (kxkyz2vmu, gvmu, gnew)
    gold = gnew

  end subroutine init_gxyz

  subroutine init_dist_fn

    use mp, only: proc0
    use stella_layouts, only: init_dist_fn_layouts
    use species, only: nspec
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx, ny, nx
    use vpamu_grids, only: nvgrid, nmu

    implicit none

    if (dist_fn_initialized) return
    dist_fn_initialized = .true.

    debug = debug .and. proc0
    
    if (debug) write (*,*) 'dist_fn::init_dist_fn::read_parameters'
    call read_parameters
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_dist_fn_layouts'
    call init_dist_fn_layouts (nzgrid, naky, nakx, nvgrid, nmu, nspec, ny, nx)
    if (debug) write (*,*) 'dist_fn::init_dist_fn::allocate_arrays'
    call allocate_arrays
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_kperp2'
    call init_kperp2
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_vperp2'
    call init_vperp2
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_bessel'
    call init_bessel

  end subroutine init_dist_fn

  subroutine read_parameters

    use file_utils, only: error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast

    implicit none

    logical :: dfexist

    type (text_option), dimension (7), parameter :: adiabaticopts = &
         (/ text_option('default', adiabatic_option_default), &
            text_option('no-field-line-average-term', adiabatic_option_default), &
            text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
            text_option('iphi00=0', adiabatic_option_default), &
            text_option('iphi00=1', adiabatic_option_default), &
            text_option('iphi00=2', adiabatic_option_fieldlineavg), &
            text_option('iphi00=3', adiabatic_option_yavg)/)
    character(30) :: adiabatic_option
            
    namelist /dist_fn_knobs/ adiabatic_option

    integer :: ierr, in_file

    if (readinit) return
    readinit = .true.

    if (proc0) then
       adiabatic_option = 'default'

       in_file = input_unit_exist("dist_fn_knobs", dfexist)
       if (dfexist) read (unit=in_file, nml=dist_fn_knobs)

       ierr = error_unit()
       call get_option_value &
            (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
            ierr, "adiabatic_option in dist_fn_knobs")
    end if

    call broadcast (adiabatic_option_switch)

  end subroutine read_parameters 

  subroutine init_kperp2

    use dist_fn_arrays, only: kperp2
    use geometry, only: gds2, gds21, gds22
    use geometry, only: geo_surf
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx, theta0
    use kt_grids, only: akx, aky
    use kt_grids, only: zonal_mode

    implicit none

    integer :: iky, ikx

    if (kp2init) return
    kp2init = .true.

    allocate (kperp2(naky,nakx,-nzgrid:nzgrid))
    do iky = 1, naky
       if (zonal_mode(iky)) then
          do ikx = 1, nakx
             kperp2(iky,ikx,:) = akx(ikx)*akx(ikx)*gds22(1,:)/(geo_surf%shat**2)
          end do
       else
          do ikx = 1, nakx
             kperp2(iky,ikx,:) = aky(iky)*aky(iky) &
                  *(gds2(1,:) + 2.0*theta0(iky,ikx)*gds21(1,:) &
                  + theta0(iky,ikx)*theta0(iky,ikx)*gds22(1,:))
          end do
       end if
    end do

  end subroutine init_kperp2

  subroutine allocate_arrays

    use stella_layouts, only: kxkyz_lo, vmu_lo
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx
    use vpamu_grids, only: nvpa, nmu
    use dist_fn_arrays, only: gnew, gold
    use dist_fn_arrays, only: gvmu

    implicit none

    if (.not.allocated(gnew)) &
         allocate (gnew(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    gnew = 0.
    if (.not.allocated(gold)) &
         allocate (gold(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    gold = 0.
    if (.not.allocated(gvmu)) &
         allocate (gvmu(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    gvmu = 0.

  end subroutine allocate_arrays

  subroutine init_vperp2

    use geometry, only: nalpha, bmag
    use zgrid, only: nzgrid
    use vpamu_grids, only: vperp2, mu
    use vpamu_grids, only: nmu

    implicit none

    integer :: imu
    
    if (vp2init) return
    vp2init = .true.

    if (.not.allocated(vperp2)) allocate (vperp2(nalpha,-nzgrid:nzgrid,nmu)) ; vperp2 = 0.
    
    do imu = 1, nmu
       vperp2(:,:,imu) = 2.0*mu(imu)*bmag
    end do

  end subroutine init_vperp2

  subroutine init_bessel

    use dist_fn_arrays, only: aj0v, aj1v
    use dist_fn_arrays, only: aj0x
    use dist_fn_arrays, only: kperp2
    use species, only: spec, nspec
    use geometry, only: bmag
    use zgrid, only: nzgrid
    use vpamu_grids, only: vperp2, nmu
    use kt_grids, only: naky, nakx
    use stella_layouts, only: kxkyz_lo, vmu_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, imu_idx
    use spfunc, only: j0, j1

    implicit none

    integer :: iz, iky, ikx, imu, is
    integer :: ikxkyz, ivmu
    real :: arg

    if (bessinit) return
    bessinit = .true.

    call init_kperp2

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
    
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          arg = spec(is)%smz*sqrt(vperp2(1,iz,imu)*kperp2(iky,ikx,iz))/bmag(1,iz)
          aj0v(imu,ikxkyz) = j0(arg)
          ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
          aj1v(imu,ikxkyz) = j1(arg)
       end do
    end do

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       do iz = -nzgrid, nzgrid
          do ikx = 1, nakx
             do iky = 1, naky
                arg = spec(is)%smz*sqrt(vperp2(1,iz,imu)*kperp2(iky,ikx,iz))/bmag(1,iz)
                aj0x(iky,ikx,iz,ivmu) = j0(arg)
             ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
!             aj1(iz,ivmu) = j1(arg)
             end do
          end do
       end do
    end do

  end subroutine init_bessel

  subroutine finish_dist_fn

    implicit none

    call finish_bessel
    call finish_kperp2
    call finish_vperp2
    call deallocate_arrays

    dist_fn_initialized = .false.
    readinit = .false.
    gxyz_initialized = .false.

  end subroutine finish_dist_fn

  subroutine deallocate_arrays

    use dist_fn_arrays, only: gnew, gold

    implicit none

    if (allocated(gnew)) deallocate (gnew)
    if (allocated(gold)) deallocate (gold)

  end subroutine deallocate_arrays

  subroutine finish_kperp2

    use dist_fn_arrays, only: kperp2

    implicit none

    if (allocated(kperp2)) deallocate (kperp2)

    kp2init = .false.

  end subroutine finish_kperp2

  subroutine finish_vperp2

    use vpamu_grids, only: vperp2

    implicit none

    if (allocated(vperp2)) deallocate (vperp2)

    vp2init = .false.
    
  end subroutine finish_vperp2

  subroutine finish_bessel

    use dist_fn_arrays, only: aj0v, aj0x

    implicit none

    if (allocated(aj0v)) deallocate (aj0v)
    if (allocated(aj0x)) deallocate (aj0x)

    bessinit = .false.

  end subroutine finish_bessel

end module dist_fn
