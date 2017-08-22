module dist_fn

  use redistribute, only: redist_type

  implicit none

  public :: init_get_fields, finish_get_fields
  public :: get_fields
  public :: init_gxyz
  public :: init_dist_fn, finish_dist_fn
  public :: advance_stella
  public :: time_gke

  private

  logical :: get_fields_initialized = .false.
  logical :: dist_fn_initialized = .false.
  logical :: gxyz_initialized = .false.
  logical :: kp2init = .false.
  logical :: wdriftinit = .false.
  logical :: wstarinit = .false.
  logical :: bessinit = .false.
  logical :: mirrorinit = .false.
  logical :: streaminit = .false.
  logical :: redistinit = .false.
  logical :: readinit = .false.

  integer :: boundary_option_switch
  integer, parameter :: boundary_option_zero = 1, &
       boundary_option_self_periodic = 2, &
       boundary_option_linked = 3

  integer :: adiabatic_option_switch
  integer, parameter :: adiabatic_option_default = 1, &
       adiabatic_option_zero = 2, &
       adiabatic_option_fieldlineavg = 3, &
       adiabatic_option_yavg = 4

  real :: xdriftknob, ydriftknob, streamknob, mirrorknob, wstarknob
  real, dimension (:,:,:), allocatable :: gamtot, apar_denom
  real, dimension (:,:), allocatable :: gamtot3

  real, dimension (:,:,:), allocatable :: kperp2

  ! needed for mirror term
  integer, dimension (:,:), allocatable :: mirror_sign
  real, dimension (:,:,:,:), allocatable :: mirror

  ! needed for parallel streaming term
  integer, dimension (:), allocatable :: stream_sign
  real, dimension (:,:,:), allocatable :: stream

  ! these arrays needed to keep track of connections between different
  ! 2pi segments
  integer, dimension (:), allocatable :: neigen
  integer, dimension (:), allocatable :: ig_low, ig_mid, ig_up
  integer, dimension (:,:), allocatable :: nsegments
  integer, dimension (:,:,:), allocatable :: ikxmod
  logical, dimension (:), allocatable :: periodic

  ! needed for timing various pieces of gke solve
  real, dimension (2,6) :: time_gke

  type (redist_type) :: kxkyz2vmu
  type (redist_type) :: kxyz2vmu

  integer :: ny_ffs

  logical :: debug = .false.

contains

  subroutine init_get_fields

    use mp, only: sum_allreduce, proc0
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, onlY: ig_idx, ikx_idx, iky_idx, is_idx
    use dist_fn_arrays, only: aj0v
    use run_parameters, only: fphi, fapar
    use run_parameters, only: tite, nine, beta
    use species, only: spec, has_electron_species
    use geometry, only: dl_over_b
    use zgrid, only: nzgrid
    use vpamu_grids, only: nvpa, nvgrid, nmu
    use vpamu_grids, only: vpa
    use vpamu_grids, only: anon, integrate_vmu
    use species, only: spec
    use kt_grids, only: naky, nakx, aky, akx

    implicit none

    integer :: ikxkyz, ig, ikx, iky, is
    real :: tmp, wgt
    real, dimension (:,:), allocatable :: g0

    if (get_fields_initialized) return
    get_fields_initialized = .true.

    debug = debug .and. proc0

    if (.not.allocated(gamtot)) allocate (gamtot(naky,nakx,-nzgrid:nzgrid)) ; gamtot = 0.
    if (.not.allocated(gamtot3)) then
       if (.not.has_electron_species(spec) &
            .and. adiabatic_option_switch==adiabatic_option_fieldlineavg) then
          allocate (gamtot3(nakx,-nzgrid:nzgrid)) ; gamtot3 = 0.
       else
          allocate (gamtot3(1,1)) ; gamtot3 = 0.
       end if
    end if
    if (.not.allocated(apar_denom)) allocate (apar_denom(naky,nakx,-nzgrid:nzgrid)) ; apar_denom = 0.

    if (fphi > epsilon(0.0)) then
       allocate (g0(-nvgrid:nvgrid,nmu))
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          ig = ig_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          g0 = spread((1.0 - aj0v(:,ikxkyz)**2),1,nvpa)*anon(ig,:,:)
          wgt = spec(is)%z*spec(is)%z*spec(is)%dens/spec(is)%temp
          call integrate_vmu (g0, ig, tmp)
          gamtot(iky,ikx,ig) = gamtot(iky,ikx,ig) + tmp*wgt
       end do
       call sum_allreduce (gamtot)
       ! avoid divide by zero when kx=ky=0
       ! do not evolve this mode, so value is irrelevant
       gamtot(1,1,:) = 1.0

       deallocate (g0)

       if (.not.has_electron_species(spec)) then
          gamtot = gamtot + tite/nine
          if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
             if (abs(aky(1)) < epsilon(0.)) then
                do ikx = 1, nakx
                   ! avoid divide by zero for kx=ky=0 mode,
                   ! which we do not need anyway
                   if (abs(akx(ikx)) < epsilon(0.)) cycle
                   tmp = nine/tite-sum(dl_over_b/gamtot(1,ikx,:))
                   gamtot3(ikx,:) = 1./(gamtot(1,ikx,:)*tmp)
                end do
             end if
          end if
       end if
    end if

    if (fapar > epsilon(0.)) then
       allocate (g0(-nvgrid:nvgrid,nmu))
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          ig = ig_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          g0 = spread(vpa**2,2,nmu)*spread(aj0v(:,ikxkyz)**2,1,nvpa)*anon(ig,:,:)
          wgt = 2.0*beta*spec(is)%z*spec(is)%z*spec(is)%dens/spec(is)%mass
          call integrate_vmu (g0, ig, tmp)
          apar_denom(iky,ikx,ig) = apar_denom(iky,ikx,ig) + tmp*wgt
       end do
       call sum_allreduce (apar_denom)
       apar_denom = apar_denom + kperp2

       deallocate (g0)
    end if

  end subroutine init_get_fields

  subroutine get_fields (g, phi, apar)

    use mp, only: sum_allreduce
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: ig_idx, ikx_idx, iky_idx, is_idx
    use dist_fn_arrays, only: aj0v
    use run_parameters, only: fphi, fapar
    use run_parameters, only: beta
    use geometry, only: dl_over_b
    use zgrid, only: nzgrid
    use vpamu_grids, only: nvgrid, nvpa, nmu
    use vpamu_grids, only: vpa
    use vpamu_grids, only: integrate_vmu
    use kt_grids, only: nakx, aky
    use species, only: spec, has_electron_species

    implicit none
    
    complex, dimension (-nvgrid:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:), intent (out) :: phi, apar

    real :: wgt
    complex, dimension (:,:), allocatable :: g0
    integer :: ikxkyz, ig, ikx, iky, is
    complex :: tmp

    phi = 0.
    if (fphi > epsilon(0.0)) then
       allocate (g0(-nvgrid:nvgrid,nmu))
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          ig = ig_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iky = iky_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          g0 = spread(aj0v(:,ikxkyz),1,nvpa)*g(:,:,ikxkyz)
          wgt = spec(is)%z*spec(is)%dens
          call integrate_vmu (g0, ig, tmp)
          phi(iky,ikx,ig) = phi(iky,ikx,ig) + wgt*tmp
       end do
       call sum_allreduce (phi)
       phi = phi/gamtot

       if (.not.has_electron_species(spec) .and. &
            adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          if (abs(aky(1)) < epsilon(0.)) then
             do ikx = 1, nakx
                tmp = sum(dl_over_b*phi(1,ikx,:))
                phi(1,ikx,:) = phi(1,ikx,:) + tmp*gamtot3(ikx,:)
             end do
          end if
       end if

       deallocate (g0)
    end if

    apar = 0.
    if (fapar > epsilon(0.0)) then
       allocate (g0(-nvgrid:nvgrid,nmu))
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          ig = ig_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iky = iky_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          g0 = spread(aj0v(:,ikxkyz),1,nvpa)*spread(vpa,2,nmu)*g(:,:,ikxkyz)
          wgt = 2.0*beta*spec(is)%z*spec(is)%dens*spec(is)%stm
          call integrate_vmu (g0, ig, tmp)
          apar(iky,ikx,ig) = apar(iky,ikx,ig) + tmp*wgt
       end do
       call sum_allreduce (apar)
       apar = apar/apar_denom
       deallocate (g0)
    end if

  end subroutine get_fields

  subroutine finish_get_fields

    implicit none

    get_fields_initialized = .false.
    if (allocated(gamtot)) deallocate (gamtot)
    if (allocated(gamtot3)) deallocate (gamtot3)

  end subroutine finish_get_fields

  subroutine init_gxyz

    use dist_fn_arrays, only: gvmu, gold, gnew
    use redistribute, only: gather

    implicit none

    if (gxyz_initialized) return
    gxyz_initialized = .false.

    ! get version of g that has ky,kx,z local
    call gather (kxkyz2vmu, gvmu, gnew)
    gold = gnew

  end subroutine init_gxyz

  subroutine init_dist_fn

    use stella_layouts, only: init_stella_layouts, init_dist_fn_layouts
    use stella_transforms, only: init_transforms
    use species, only: init_species
    use species, only: nspec
    use zgrid, only: init_zgrid
    use zgrid, only: nzgrid
    use kt_grids, only: init_kt_grids
    use kt_grids, only: naky, nakx, ny
    use kt_grids, only: alpha_global
    use vpamu_grids, only: init_vpamu_grids
    use vpamu_grids, only: nvgrid, nmu
    use run_parameters, only: init_run_parameters
    use neoclassical_terms, only: init_neoclassical_terms

    implicit none

    if (dist_fn_initialized) return
    dist_fn_initialized = .true.

    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_stella_layouts'
    call init_stella_layouts
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_species'
    call init_species
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_zgrid'
    call init_zgrid
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_kt_grids'
    call init_kt_grids

    ! TMP FOR TESTING
    if (alpha_global) then
       ny_ffs = ny
    else
       ny_ffs = 1
    end if

    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_vpamu_grids'
    call init_vpamu_grids
    if (debug) write (*,*) 'dist_fn::init_dist_fn::read_parameters'
    call read_parameters
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_run_parameters'
    call init_run_parameters
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_dist_fn_layouts'
    call init_dist_fn_layouts (nzgrid, naky, nakx, nvgrid, nmu, nspec, ny)
    if (debug) write (*,*) 'dist_fn::init_dist_fn::allocate_arrays'
    call allocate_arrays
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_kperp2'
    call init_kperp2
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_connections'
    call init_connections
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_vperp2'
    call init_vperp2
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_bessel'
    call init_bessel
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_mirror'
    call init_mirror
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_parstream'
    call init_parstream
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_wdrift'
    call init_wdrift
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_wstar'
    call init_wstar
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_neoclassical_terms'
    call init_neoclassical_terms
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_redistribute'
    call init_redistribute
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_cfl'
    call init_cfl
    if (alpha_global) then
       if (debug) write (*,*) 'dist_fn::init_dist_fn::init_transforms'
       call init_transforms
    end if

  end subroutine init_dist_fn

  subroutine read_parameters

    use file_utils, only: input_unit, error_unit, input_unit_exist
    use geometry, only: shat
    use zgrid, only: shat_zero
    use text_options, only: text_option, get_option_value
    use species, only: nspec
    use mp, only: proc0, broadcast

    implicit none

    logical :: dfexist

    type (text_option), dimension (6), parameter :: boundaryopts = &
         (/ text_option('default', boundary_option_zero), &
            text_option('zero', boundary_option_zero), &
            text_option('unconnected', boundary_option_zero), &
            text_option('self-periodic', boundary_option_self_periodic), &
            text_option('periodic', boundary_option_self_periodic), &
            text_option('linked', boundary_option_linked) /)
    character(20) :: boundary_option

    type (text_option), dimension (7), parameter :: adiabaticopts = &
         (/ text_option('default', adiabatic_option_default), &
            text_option('no-field-line-average-term', adiabatic_option_default), &
            text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
            text_option('iphi00=0', adiabatic_option_default), &
            text_option('iphi00=1', adiabatic_option_default), &
            text_option('iphi00=2', adiabatic_option_fieldlineavg), &
            text_option('iphi00=3', adiabatic_option_yavg)/)
    character(30) :: adiabatic_option
            
    namelist /dist_fn_knobs/ boundary_option, &
         xdriftknob, ydriftknob, streamknob, mirrorknob, wstarknob, &
         adiabatic_option
    integer :: ierr, in_file

    if (readinit) return
    readinit = .true.

    if (proc0) then
       boundary_option = 'default'
       adiabatic_option = 'default'
       xdriftknob = 1.0
       ydriftknob = 1.0
       streamknob = 1.0
       mirrorknob = 1.0
       wstarknob = 1.0

       in_file = input_unit_exist("dist_fn_knobs", dfexist)
       if (dfexist) read (unit=in_file, nml=dist_fn_knobs)

       if(abs(shat) <=  shat_zero) boundary_option = 'periodic'

       ierr = error_unit()
       call get_option_value &
            (boundary_option, boundaryopts, boundary_option_switch, &
            ierr, "boundary_option in dist_fn_knobs")

       call get_option_value &
            (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
            ierr, "adiabatic_option in dist_fn_knobs")

    end if

    call broadcast (boundary_option_switch)
    call broadcast (adiabatic_option_switch)
    call broadcast (xdriftknob)
    call broadcast (ydriftknob)
    call broadcast (streamknob)
    call broadcast (mirrorknob)
    call broadcast (wstarknob)

  end subroutine read_parameters 

  subroutine init_kperp2

    use geometry, only: gds2, gds21, gds22, shat
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx, theta0
    use kt_grids, only: akx, aky

    implicit none

    integer :: iky, ikx

    if (kp2init) return
    kp2init = .true.

    allocate (kperp2(naky,nakx,-nzgrid:nzgrid))
    do iky = 1, naky
       if (abs(aky(iky)) < epsilon(0.)) then
          do ikx = 1, nakx
             kperp2(iky,ikx,:) = akx(ikx)*akx(ikx)*gds22/(shat*shat)
          end do
       else
          do ikx = 1, nakx
             kperp2(iky,ikx,:) = aky(iky)*aky(iky) &
                  *(gds2 + 2.0*theta0(iky,ikx)*gds21 &
                  + theta0(iky,ikx)*theta0(iky,ikx)*gds22)
          end do
       end if
    end do

  end subroutine init_kperp2

  subroutine init_wdrift

    use dist_fn_arrays, only: wdriftx, wdrifty
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use species, only: spec
    use zgrid, only: nzgrid
    use geometry, only: cvdrift, gbdrift
    use geometry, only: cvdrift0, gbdrift0, shat
    use vpamu_grids, only: vpa, vperp2

    implicit none

    integer :: ivmu, iv, imu, is, iy

    if (wdriftinit) return
    wdriftinit = .true.

    if (.not.allocated(wdriftx)) &
         allocate (wdriftx(ny_ffs,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    if (.not.allocated(wdrifty)) &
         allocate (wdrifty(ny_ffs,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    ! FLAG -- need to deal with shat=0 case.  ideally move away from q as x-coordinate
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do iy = 1, ny_ffs
          wdrifty(iy,:,ivmu) = -ydriftknob*0.5*code_dt*spec(is)%tz &
               * (cvdrift*vpa(iv)**2 + gbdrift*0.5*vperp2(:,imu))
          wdriftx(iy,:,ivmu) = -xdriftknob*0.5*code_dt*spec(is)%tz/shat &
               * (cvdrift0*vpa(iv)**2 + gbdrift0*0.5*vperp2(:,imu))
       end do
    end do

  end subroutine init_wdrift

  subroutine init_wstar

    use dist_fn_arrays, only: wstar
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use species, only: spec
    use zgrid, only: nzgrid
    use vpamu_grids, only: energy, anon

    implicit none

    integer :: is, imu, iv, ivmu

    if (wstarinit) return
    wstarinit = .true.

    if (.not.allocated(wstar)) &
         allocate (wstar(-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc)) ; wstar = 0.0

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       iv = iv_idx(vmu_lo,ivmu)
       wstar(:,ivmu) = wstarknob*0.5*code_dt*anon(:,iv,imu) &
            * (spec(is)%fprim+spec(is)%tprim*(energy(:,iv,imu)-1.5))
    end do

  end subroutine init_wstar

  subroutine allocate_arrays

    use stella_layouts, only: kxkyz_lo, vmu_lo
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx
    use vpamu_grids, only: nvgrid, nmu
    use dist_fn_arrays, only: gnew, gold
    use dist_fn_arrays, only: g1, g2
    use dist_fn_arrays, only: gvmu

    implicit none

    if (.not.allocated(gnew)) &
         allocate (gnew(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    gnew = 0.
    if (.not.allocated(gold)) &
         allocate (gold(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    gold = 0.
    if (.not.allocated(g1)) &
         allocate (g1(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    g1 = 0.
    if (.not.allocated(g2)) &
         allocate (g2(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    g2 = 0.
    if (.not.allocated(gvmu)) &
         allocate (gvmu(-nvgrid:nvgrid,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    gvmu = 0.

  end subroutine allocate_arrays

  subroutine init_connections

    use zgrid, only: nperiod, nzgrid, nzed
    use kt_grids, only: ntheta0, nakx, naky
    use kt_grids, only: jtwist_out, aky
    use species, only: nspec
    use vpamu_grids, only: nmu, nvgrid

    implicit none

    integer :: iseg, iky, ie, ntg, ikx, ikxshiftend
    integer :: nseg_max, neigen_max
    integer :: iky_max
    integer, dimension (:), allocatable :: ikx_shift_left_kypos, ikx_shift_left_kyneg
    integer, dimension (:,:), allocatable :: ikx_shift

    ntg = nzed/2
    ! iky_max is the index of the most positive ky
    iky_max = naky/2+1

    if (.not. allocated(neigen)) allocate (neigen(naky))
    if (.not. allocated(periodic)) allocate (periodic(naky)) ; periodic = .false.

    if (boundary_option_switch==boundary_option_self_periodic) then
       periodic = .true.
    else
       where (abs(aky) < epsilon(0.0)) periodic = .true.
    end if

    select case (boundary_option_switch)
    case (boundary_option_linked)

       ! ntheta0 = 2*(nakx-1) + 1
       ! nakx includes kx >= 0
       ! ntheta0 also includes kx < 0
       neigen(1) = ntheta0
       if (naky > 1) then
          do iky = 2, iky_max
             ! must link different kx values at theta = +/- pi
             ! neigen is the number of independent eigenfunctions along the field line
             neigen(iky) = min((iky-1)*jtwist_out,ntheta0)
          end do
          ! number of eigenfunctions for -ky is same as for +ky
          neigen(iky_max+1:) = neigen(iky_max:2:-1)
       end if

       neigen_max = maxval(neigen)

       if (.not. allocated(ikx_shift_left_kypos)) then
          allocate (ikx_shift_left_kypos(neigen_max)) ; ikx_shift_left_kypos = 0
          allocate (ikx_shift_left_kyneg(neigen_max)) ; ikx_shift_left_kyneg = 0
          allocate (ikx_shift(ntheta0,naky)) ; ikx_shift = 0
       end if

       ! figure out how much to shift ikx by to get to
       ! the left-most (theta-theta0) in each set of connected 2pi segments
       ! note that theta0 goes from 0 to theta0_max and then from theta0_min back
       ! to -dtheta0
       do ikx = 1, neigen_max
          ! first ntheta0/2+1 theta0s are 0 and all positive theta0 values
          ! remainder are negative theta0s
          ! note that ntheta0 is always positive for box
          ! theta_0 = kx / ky / shat
          ! if ky > 0, then most positive theta_0 corresponds to most positive kx
          if (ikx <= nakx) then
             ikx_shift_left_kypos(ikx) = nakx-2*ikx+1
          else
             ikx_shift_left_kypos(ikx) = 3*nakx-2*ikx
          end if
          ! if ky < 0, most positive theta_0 corresponds to most negative kx
          if (ikx < nakx) then
             ikx_shift_left_kyneg(ikx) = nakx
          else
             ikx_shift_left_kyneg(ikx) = 1-nakx
          end if
       end do

       do iky = 1, naky
          ! ikx_shift is how much to shift each ikx by to connect
          ! to the next theta0 (from most positive to most negative)

          ! if ky < 0, then going to more negative theta0
          ! corresponds to going to more positive kx
          if (aky(iky) < 0.0) then
             ! first treat kx positive
             ! connect to more positive kx neigen away
             ! but only if not trying to connect to kx
             ! so positive that ikx is not on grid
             if (nakx - neigen(iky) > 0) ikx_shift(:nakx-neigen(iky),iky) = neigen(iky)
             ! next treat kx negative
             ! if kx sufficiently negative, then 
             ! shifting by neigen keeps kx negative
             do ikx = nakx+1, ntheta0
                if (ikx+neigen(iky) <= ntheta0) then
                   ikx_shift(ikx,iky) = neigen(iky)
                   ! if theta0 not sufficiently negative,
                   ! then must shift to postive theta0
                else if (ikx+neigen(iky) <= ntheta0+nakx) then
                   ikx_shift(ikx,iky) = neigen(iky) - ntheta0
                end if
             end do
          else
             ! if ky > 0, then going to more negative theta0
             ! corresponds to going to more negative kx
             do ikx = 1, nakx
                ! if theta0 is sufficiently positive, shifting to more
                ! negative theta0 corresponds to decreasing ikx
                if (ikx-neigen(iky) > 0) then
                   ikx_shift(ikx,iky) = -neigen(iky)
                   ! if a positive theta0 connects to a negative theta0
                   ! must do more complicated mapping of ikx
                else if (ikx-neigen(iky)+ntheta0 >= nakx+1) then
                   ikx_shift(ikx,iky) = ntheta0 - neigen(iky)
                end if
             end do
             ! if theta0 is negative, then shifting to more negative
             ! theta0 corresponds to decreasing ikx
             do ikx = nakx+1, ntheta0
                ! if theta0 is sufficiently negative, it has no
                ! more negative theta0 with which it can connect
                if (ikx-neigen(iky) >= nakx) then
                   ikx_shift(ikx,iky) = -neigen(iky)
                end if
                ! theta0 is positive
             end  do
          end if
       end do

       if (.not. allocated(nsegments)) then
          allocate (nsegments(neigen_max,naky))
       end if

       do iky = 1, naky
          if (neigen(iky) == 0) then
             nsegments(:,iky) = 1
          else
             nsegments(:,iky) = (ntheta0-1)/neigen(iky)

             do ie = 1, mod(ntheta0-1,neigen(iky))+1
                nsegments(ie,iky) = nsegments(ie,iky) + 1
             end do
          end if
       end do

       nseg_max = maxval(nsegments)

       if (.not. allocated(ig_low)) then
          allocate (ig_low(nseg_max)) ; ig_low = -nzgrid
          allocate (ig_mid(nseg_max)) ; ig_mid = 0
          allocate (ig_up(nseg_max)) ; ig_up = nzgrid
       end if
       
    case default
       
       neigen = ntheta0 ; neigen_max = ntheta0
       
       if (.not. allocated(ikx_shift_left_kypos)) then
          allocate (ikx_shift_left_kypos(neigen_max))
          allocate (ikx_shift_left_kyneg(neigen_max))
          allocate (ikx_shift(ntheta0,naky))
       end if
       ikx_shift = 0 ; ikx_shift_left_kypos = 0 ; ikx_shift_left_kyneg = 0
       
       if (.not. allocated(nsegments)) then
          allocate (nsegments(neigen_max,naky))
       end if
       
       ! this is the number of 2pi poloidal segments in the extended theta domain,
       ! which is needed in initializing the reponse matrix and doing the implicit sweep
       nsegments = 2*(nperiod-1) + 1
       
       nseg_max = maxval(nsegments)
       
       if (.not. allocated(ig_low)) then
          allocate (ig_low(nseg_max))
          allocate (ig_mid(nseg_max))
          allocate (ig_up(nseg_max))
       end if

       ! ig_low(j) is the ig index corresponding to the inboard midplane from below (theta=-pi) within the jth segment
       ! ig_mid(j) is the ig index corresponding to the outboard midplane (theta=0) within the jth segment
       do iseg = 1, nseg_max
          ig_low(iseg) = -nzgrid + (iseg-1)*nzed
          ig_mid(iseg) = ig_low(iseg) + nzed/2
          ig_up(iseg) = ig_low(iseg) + nzed
       end do

    end select

    if (.not. allocated(ikxmod)) allocate (ikxmod(nseg_max,neigen_max,naky))
    do iky = 1, naky
       ! only do the following once for each independent set of theta0s
       ! the assumption here is that all kx are on processor and sequential
       do ikx = 1, neigen(iky)
          if (aky(iky) < 0.) then
             ikxshiftend = ikx_shift_left_kyneg(ikx)
          else
             ikxshiftend = ikx_shift_left_kypos(ikx)
          end if
          ! remap to start at theta0 = theta0_max
          ! (so that theta-theta0 is most negative)
          ! for this set of connected theta0s
          iseg = 1
          ikxmod(iseg,ikx,iky) = ikx + ikxshiftend
          if (nsegments(ikx,iky) > 1) then
             do iseg = 2, nsegments(ikx,iky)
                ikxmod(iseg,ikx,iky) = ikxmod(iseg-1,ikx,iky) + ikx_shift(ikxmod(iseg-1,ikx,iky),iky)
             end do
          end if
       end do
    end do

    if (allocated(ikx_shift_left_kypos)) deallocate (ikx_shift_left_kypos)
    if (allocated(ikx_shift_left_kyneg)) deallocate (ikx_shift_left_kyneg)
    if (allocated(ikx_shift)) deallocate (ikx_shift)

  end subroutine init_connections

  subroutine init_vperp2

    use geometry, only: bmag
    use zgrid, only: nzgrid
    use vpamu_grids, only: vperp2, energy, anon
    use vpamu_grids, only: vpa, mu
    use vpamu_grids, only: nmu, nvgrid

    implicit none

    integer :: iv
    
    if (.not.allocated(vperp2)) allocate (vperp2(-nzgrid:nzgrid,nmu)) ; vperp2 = 0.
    if (.not.allocated(energy)) allocate (energy(-nzgrid:nzgrid,-nvgrid:nvgrid,nmu)) ; energy = 0.
    if (.not.allocated(anon)) allocate (anon(-nzgrid:nzgrid,-nvgrid:nvgrid,nmu)) ; anon = 0.

    vperp2 = 2.0*spread(mu,1,2*nzgrid+1)*spread(bmag,2,nmu)

    do iv = -nvgrid, nvgrid
       energy(:,iv,:) = vpa(iv)**2 + 2.0*spread(mu,1,2*nzgrid+1)*spread(bmag,2,nmu)
       anon(:,iv,:) = exp(-energy(:,iv,:))
    end do

  end subroutine init_vperp2

  subroutine init_bessel

    use dist_fn_arrays, only: aj0v, aj0x
    use species, only: spec, nspec
    use geometry, only: bmag
    use zgrid, only: nzgrid
    use vpamu_grids, only: vperp2, nmu
    use kt_grids, only: naky, nakx
    use stella_layouts, only: kxkyz_lo, vmu_lo
    use stella_layouts, only: iky_idx, ikx_idx, ig_idx, is_idx, imu_idx
    use spfunc, only: j0!, j1

    implicit none

    integer :: ig, iky, ikx, imu, is
    integer :: ikxkyz, ivmu
    real :: arg

    if (bessinit) return
    bessinit = .true.

    call init_kperp2

    allocate (aj0v(nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc)) ; aj0v = 0.
    allocate (aj0x(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc)) ; aj0x = 0.
!    allocate (aj1(-nzgrid:nzgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; aj1 = 0.
!    allocate (aj2(-nzgrid:nzgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; aj2 = 0.

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       ig = ig_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          arg = spec(is)%smz*sqrt(vperp2(ig,imu)*kperp2(iky,ikx,ig))/bmag(ig)
          aj0v(imu,ikxkyz) = j0(arg)
             ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
!             aj1(ig,ikxkyz) = j1(arg)
!             aj2(ig,ikxkyz) = 2.0*aj1(ig,ikxkyz)-aj0(ig,ikxkyz)
       end do
    end do

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       do ig = -nzgrid, nzgrid
          do ikx = 1, nakx
             do iky = 1, naky
                arg = spec(is)%smz*sqrt(vperp2(ig,imu)*kperp2(iky,ikx,ig))/bmag(ig)
                aj0x(iky,ikx,ig,ivmu) = j0(arg)
             ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
!             aj1(ig,ivmu) = j1(arg)
!             aj2(ig,ivmu) = 2.0*aj1(ig,ivmu)-aj0(ig,ivmu)
             end do
          end do
       end do
    end do

  end subroutine init_bessel

  subroutine init_mirror

    use stella_time, only: code_dt
    use species, only: spec, nspec
    use vpamu_grids, only: nmu
    use vpamu_grids, only: mu
    use zgrid, only: nzgrid
    use geometry, only: dbdthet, gradpar

    implicit none

    integer :: ig, iy

    if (mirrorinit) return
    mirrorinit = .true.
    
    if (.not.allocated(mirror)) allocate (mirror(ny_ffs,-nzgrid:nzgrid,nmu,nspec)) ; mirror = 0.
    if (.not.allocated(mirror_sign)) allocate (mirror_sign(ny_ffs,-nzgrid:nzgrid)) ; mirror_sign = 0

    do iy = 1, ny_ffs
       mirror(iy,:,:,:) = mirrorknob*code_dt*spread(spread(spec%stm,1,2*nzgrid+1),2,nmu) &
            *spread(spread(mu,1,2*nzgrid+1)*spread(dbdthet*gradpar,2,nmu),3,nspec)
       ! mirror_sign set to +/- 1 depending on the sign of the mirror term.
       ! NB: mirror_sign = -1 corresponds to positive advection velocity
       do ig = -nzgrid, nzgrid
          mirror_sign(iy,ig) = int(sign(1.0,mirror(iy,ig,1,1)))
       end do
    end do

  end subroutine init_mirror

  subroutine init_parstream

    use stella_time, only: code_dt
    use species, only: spec, nspec
    use vpamu_grids, only: nvgrid
    use vpamu_grids, only: vpa, nvpa
    use zgrid, only: nzgrid
    use geometry, only: gradpar

    implicit none

    integer :: iv

    if (streaminit) return
    streaminit = .true.

    if (.not.allocated(stream)) allocate (stream(-nzgrid:nzgrid,-nvgrid:nvgrid,nspec)) ; stream = 0.
    if (.not.allocated(stream_sign)) allocate (stream_sign(-nvgrid:nvgrid)) ; stream_sign = 0

    ! sign of stream corresponds to appearing on RHS of GK equation
    stream = -streamknob*code_dt*spread(spread(spec%stm,1,2*nzgrid+1),2,nvpa) &
         * spread(spread(vpa,1,2*nzgrid+1)*spread(gradpar,2,nvpa),3,nspec)

    ! stream_sign set to +/- 1 depending on the sign of the parallel streaming term.
    ! NB: stream_sign = -1 corresponds to positive advection velocity
    do iv = -nvgrid, nvgrid
       stream_sign(iv) = int(sign(1.0,stream(0,iv,1)))
    end do

  end subroutine init_parstream

  subroutine init_redistribute

    use kt_grids, only: alpha_global

    implicit none

    if (redistinit) return
    redistinit = .true.

    if (debug) write (*,*) 'dist_fn::init_redistribute::init_kxkyz_to_gzkxky_redistribute'
    call init_kxkyz_to_vmu_redistribute
    if (alpha_global) call init_kxyz_to_vmu_redistribute

  end subroutine init_redistribute

  subroutine init_kxkyz_to_vmu_redistribute

    use mp, only: nproc
    use stella_layouts, only: kxkyz_lo, vmu_lo
    use stella_layouts, only: kxkyzidx2vmuidx
    use stella_layouts, only: idx_local, proc_id
    use redistribute, only: index_list_type, init_redist
    use redistribute, only: delete_list, set_redist_character_type
    use vpamu_grids, only: nvgrid, nmu
    use zgrid, only: nzgrid

    implicit none

    type (index_list_type), dimension (0:nproc-1) :: to_list, from_list
    integer, dimension (0:nproc-1) :: nn_to, nn_from
    integer, dimension (3) :: from_low, from_high
    integer, dimension (4) :: to_high, to_low
    integer :: ikxkyz, ivmu
    integer :: iv, imu, iky, ikx, ig
    integer :: ip, n
    logical :: initialized = .false.

    if (initialized) return
    initialized = .true.

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do ikxkyz = kxkyz_lo%llim_world, kxkyz_lo%ulim_world
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             call kxkyzidx2vmuidx (iv, imu, ikxkyz, kxkyz_lo, vmu_lo, iky, ikx, ig, ivmu)
             if (idx_local(kxkyz_lo,ikxkyz)) &
                  nn_from(proc_id(vmu_lo,ivmu)) = nn_from(proc_id(vmu_lo,ivmu)) + 1
             if (idx_local(vmu_lo,ivmu)) &
                  nn_to(proc_id(kxkyz_lo,ikxkyz)) = nn_to(proc_id(kxkyz_lo,ikxkyz)) + 1
          end do
       end do
    end do
    
    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          allocate (to_list(ip)%third(nn_to(ip)))
          allocate (to_list(ip)%fourth(nn_to(ip)))
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0

    ! loop over all vmu indices, find corresponding y indices
    do ikxkyz = kxkyz_lo%llim_world, kxkyz_lo%ulim_world
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             ! obtain corresponding y indices
             call kxkyzidx2vmuidx (iv, imu, ikxkyz, kxkyz_lo, vmu_lo, iky, ikx, ig, ivmu)
             ! if vmu index local, set:
             ! ip = corresponding y processor
             ! from_list%first-third arrays = iv,imu,ikxkyz  (ie vmu indices)
             ! later will send from_list to proc ip
             if (idx_local(kxkyz_lo,ikxkyz)) then
                ip = proc_id(vmu_lo,ivmu)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n) = iv
                from_list(ip)%second(n) = imu
                from_list(ip)%third(n) = ikxkyz
             end if
             ! if y index local, set ip to corresponding vmu processor
             ! set to_list%first,second arrays = iky,iy  (ie y indices)
             ! will receive to_list from ip
             if (idx_local(vmu_lo,ivmu)) then
                ip = proc_id(kxkyz_lo,ikxkyz)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n) = iky
                to_list(ip)%second(n) = ikx
                to_list(ip)%third(n) = ig
                to_list(ip)%fourth(n) = ivmu
             end if
          end do
       end do
    end do

    from_low(1) = -nvgrid
    from_low(2) = 1
    from_low(3) = kxkyz_lo%llim_proc

    from_high(1) = nvgrid
    from_high(2) = nmu
    from_high(3) = kxkyz_lo%ulim_alloc

    to_low(1) = 1
    to_low(2) = 1
    to_low(3) = -nzgrid
    to_low(4) = vmu_lo%llim_proc

    to_high(1) = vmu_lo%naky
    to_high(2) = vmu_lo%nakx
    to_high(3) = vmu_lo%nzed
    to_high(4) = vmu_lo%ulim_alloc

    call set_redist_character_type (kxkyz2vmu, 'kxkyz2vmu')

    call init_redist (kxkyz2vmu, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

  end subroutine init_kxkyz_to_vmu_redistribute

  subroutine init_kxyz_to_vmu_redistribute

    use mp, only: nproc
    use stella_layouts, only: kxyz_lo, vmu_lo
    use stella_layouts, only: kxyzidx2vmuidx
    use stella_layouts, only: idx_local, proc_id
    use redistribute, only: index_list_type, init_redist
    use redistribute, only: delete_list, set_redist_character_type
    use vpamu_grids, only: nvgrid, nmu
    use zgrid, only: nzgrid

    implicit none

    type (index_list_type), dimension (0:nproc-1) :: to_list, from_list
    integer, dimension (0:nproc-1) :: nn_to, nn_from
    integer, dimension (3) :: from_low, from_high
    integer, dimension (4) :: to_high, to_low
    integer :: ikxyz, ivmu
    integer :: iv, imu, iy, ikx, ig
    integer :: ip, n
    logical :: initialized = .false.

    if (initialized) return
    initialized = .true.

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do ikxyz = kxyz_lo%llim_world, kxyz_lo%ulim_world
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             call kxyzidx2vmuidx (iv, imu, ikxyz, kxyz_lo, vmu_lo, iy, ikx, ig, ivmu)
             if (idx_local(kxyz_lo,ikxyz)) &
                  nn_from(proc_id(vmu_lo,ivmu)) = nn_from(proc_id(vmu_lo,ivmu)) + 1
             if (idx_local(vmu_lo,ivmu)) &
                  nn_to(proc_id(kxyz_lo,ikxyz)) = nn_to(proc_id(kxyz_lo,ikxyz)) + 1
          end do
       end do
    end do
    
    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first(nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third(nn_from(ip)))
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first(nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          allocate (to_list(ip)%third(nn_to(ip)))
          allocate (to_list(ip)%fourth(nn_to(ip)))
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0

    ! loop over all vmu indices, find corresponding y indices
    do ikxyz = kxyz_lo%llim_world, kxyz_lo%ulim_world
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             ! obtain corresponding y indices
             call kxyzidx2vmuidx (iv, imu, ikxyz, kxyz_lo, vmu_lo, iy, ikx, ig, ivmu)
             ! if vmu index local, set:
             ! ip = corresponding y processor
             ! from_list%first-third arrays = iv,imu,ikxyz  (ie vmu indices)
             ! later will send from_list to proc ip
             if (idx_local(kxyz_lo,ikxyz)) then
                ip = proc_id(vmu_lo,ivmu)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n) = iv
                from_list(ip)%second(n) = imu
                from_list(ip)%third(n) = ikxyz
             end if
             ! if y index local, set ip to corresponding vmu processor
             ! set to_list%first,second arrays = iy,iy  (ie y indices)
             ! will receive to_list from ip
             if (idx_local(vmu_lo,ivmu)) then
                ip = proc_id(kxyz_lo,ikxyz)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n) = iy
                to_list(ip)%second(n) = ikx
                to_list(ip)%third(n) = ig
                to_list(ip)%fourth(n) = ivmu
             end if
          end do
       end do
    end do

    from_low(1) = -nvgrid
    from_low(2) = 1
    from_low(3) = kxyz_lo%llim_proc

    from_high(1) = nvgrid
    from_high(2) = nmu
    from_high(3) = kxyz_lo%ulim_alloc

    to_low(1) = 1
    to_low(2) = 1
    to_low(3) = -nzgrid
    to_low(4) = vmu_lo%llim_proc

    to_high(1) = vmu_lo%ny
    to_high(2) = vmu_lo%nakx
    to_high(3) = vmu_lo%nzed
    to_high(4) = vmu_lo%ulim_alloc

    call set_redist_character_type (kxyz2vmu, 'kxkyz2vmu')

    call init_redist (kxyz2vmu, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

  end subroutine init_kxyz_to_vmu_redistribute

  subroutine init_cfl
    
    use mp, only: proc0
    use dist_fn_arrays, only: wdriftx, wdrifty
    use stella_time, only: cfl_dt, code_dt, write_dt
    use zgrid, only: delzed
    use vpamu_grids, only: dvpa
    use kt_grids, only: akx, aky

    implicit none
    
    real :: cfl_dt_mirror, cfl_dt_stream
    real :: cfl_dt_wdriftx, cfl_dt_wdrifty
    real :: zero

    ! avoid divide by zero in cfl_dt terms below
    zero = 100.*epsilon(0.)

    cfl_dt_mirror = code_dt*dvpa/max(maxval(abs(mirror)),zero)
    ! FLAG -- assuming equal spacing in zed!
    cfl_dt_stream = code_dt*delzed(0)/max(maxval(abs(stream)),zero)
    cfl_dt_wdriftx = code_dt/max(maxval(akx)*maxval(abs(wdriftx)),zero)
    cfl_dt_wdrifty = code_dt/max(maxval(abs(aky))*maxval(abs(wdrifty)),zero)
    cfl_dt = min(cfl_dt_mirror,cfl_dt_stream,cfl_dt_wdriftx,cfl_dt_wdrifty)
    
    if (proc0) write (*,'(a43,4e12.4)') 'cfl_dt (mirror, stream, wdriftx, wdrifty):', &
         cfl_dt_mirror, cfl_dt_stream, cfl_dt_wdriftx, cfl_dt_wdrifty

    if (code_dt > cfl_dt) then
       if (proc0) then
          write (*,'(a21,e12.4,a32)') 'user-specified delt =', code_dt, 'does not satisfy CFL condition.'
          write (*,'(a32,e12.4)') 'changing code_dt to the CFL dt =', cfl_dt
       end if
       code_dt = cfl_dt
       call reset_dt
    else if (proc0) then
       call write_dt
       write (*,*)
    end if

  end subroutine init_cfl

  subroutine reset_dt

    implicit none

    ! need to recompute mirror and streaming terms
    ! to account for updated code_dt
    mirrorinit = .false.
    streaminit = .false.
    wdriftinit = .false.
    wstarinit = .false.
    call init_wstar
    call init_wdrift
    call init_mirror
    call init_parstream

  end subroutine reset_dt

  subroutine advance_stella

    use dist_fn_arrays, only: gold, gnew
    use dist_fn_arrays, only: g1, g2

    implicit none

    ! SSP rk3 algorithm
    ! if GK equation written as dg/dt = rhs,
    ! solve_gke returns rhs*dt
    call solve_gke (gold, g1)
    g1 = gold + g1
    call solve_gke (g1, g2)
    g2 = g1 + g2
    call solve_gke (g2, gnew)
    gnew = gold/3. + 0.5*g1 + (g2 + gnew)/6.

    gold = gnew

  end subroutine advance_stella

  subroutine solve_gke (gin, rhs_ky)

    use job_manage, only: time_message
    use dist_fn_arrays, only: gvmu
    use dist_fn_arrays, only: g_adjust
    use fields_arrays, only: phi, apar
    use stella_layouts, only: vmu_lo
    use stella_transforms, only: transform_y2ky
    use redistribute, only: gather, scatter
    use run_parameters, only: fphi, fapar
    use zgrid, only: nzgrid
    use vpamu_grids, only: nvgrid, nmu
    use kt_grids, only: nakx, ny
    use kt_grids, only: alpha_global

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: gin
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (out), target :: rhs_ky

    complex, dimension (:,:,:,:), allocatable, target :: rhs_y
    complex, dimension (:,:,:,:), pointer :: rhs

    rhs_ky = 0.

    if (alpha_global) then
       allocate (rhs_y(ny,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_proc))
       rhs_y = 0.
       rhs => rhs_y
    else
       rhs => rhs_ky
    end if

    ! start with g in k-space and (ky,kx,z) local

    ! obtain fields corresponding to g
    call advance_fields (gin)
    
    ! switch from g = h + (Ze/T)*<chi>*F_0 to h = f + (Ze/T)*phi*F_0
    call g_adjust (gin, phi, apar, fphi, fapar)
    call g_adjust (gvmu, phi, apar, fphi, fapar)

    ! calculate and add mirror term to RHS of GK eqn
    call advance_mirror (gin, rhs)
    ! calculate and add alpha-component of magnetic drift term to RHS of GK eqn
    call advance_wdrifty (gin, rhs)
    ! calculate and add psi-component of magnetic drift term to RHS of GK eqn
    call advance_wdriftx (gin, rhs)

    if (alpha_global) then
       call transform_y2ky (rhs_y, rhs_ky)
       deallocate (rhs_y)
    end if

    ! calculate and add parallel streaming term to RHS of GK eqn
    call advance_parallel_streaming (gin, rhs_ky)
    ! calculate and add omega_* term to RHS of GK eqn
    call advance_wstar (rhs_ky)
    ! calculate and add collision term to RHS of GK eqn
!    call advance_collisions

    ! switch from h back to g
    call g_adjust (gin, phi, apar, -fphi, -fapar)

    nullify (rhs)

  end subroutine solve_gke

  subroutine advance_fields (g)

    use mp, only: proc0
    use stella_layouts, only: vmu_lo
    use job_manage, only: time_message
    use redistribute, only: scatter
    use dist_fn_arrays, only: gvmu
    use fields_arrays, only: phi, apar
    use zgrid, only: nzgrid

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g

    ! time the communications + field solve
    if (proc0) call time_message(.false.,time_gke(:,1),' fields')
    ! first gather (vpa,mu) onto processor for v-space operations
    ! v-space operations are field solve, dg/dvpa, and collisions
    if (debug) write (*,*) 'dist_fn::advance_stella::scatter'
    call scatter (kxkyz2vmu, g, gvmu)
    ! given gvmu with vpa and mu local, calculate the corresponding fields
    if (debug) write (*,*) 'dist_fn::advance_stella::get_fields'
    call get_fields (gvmu, phi, apar)
    ! time the communications + field solve
    if (proc0) call time_message(.false.,time_gke(:,1),' fields')

  end subroutine advance_fields

  subroutine advance_parallel_streaming (g, gout)

    use mp, only: proc0
    use stella_layouts, only: vmu_lo
    use job_manage, only: time_message
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: gout

    complex, dimension (:,:,:,:), allocatable :: g0

    allocate (g0(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    ! parallel streaming stays in ky,kx,z space with ky,kx,z local
    if (proc0) call time_message(.false.,time_gke(:,3),' Stream advance')
    ! get dg/dz, with z the parallel coordinate and store in g0_kykxz
    if (debug) write (*,*) 'dist_fn::solve_gke::get_dgdz'
    call get_dgdz (g, g0)
    ! multiply dg/dz with vpa*(b . grad z) and add to source (RHS of GK equation)
    if (debug) write (*,*) 'dist_fn::solve_gke::add_stream_term'
    call add_stream_term (g0, gout)
    if (proc0) call time_message(.false.,time_gke(:,3),' Stream advance')
    deallocate (g0)

  end subroutine advance_parallel_streaming

  subroutine advance_wstar (gout)

    use mp, only: proc0
    use job_manage, only: time_message
    use fields_arrays, only: phi, apar
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: gout

    complex, dimension (:,:,:,:), allocatable :: g0

    allocate (g0(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    ! omega_* stays in ky,kx,z space with ky,kx,z local
    if (proc0) call time_message(.false.,time_gke(:,6),' wstar advance')
    ! get d<chi>/dy
    if (debug) write (*,*) 'dist_fn::solve_gke::get_dchidy'
    call get_dchidy (phi, apar, g0)
    ! multiply with omega_* coefficient and add to source (RHS of GK eqn)
    if (debug) write (*,*) 'dist_fn::solve_gke::add_wstar_term'
    call add_wstar_term (g0, gout)
    if (proc0) call time_message(.false.,time_gke(:,6),' wstar advance')
    deallocate (g0)

  end subroutine advance_wstar

  subroutine advance_mirror (g, gout)

    use mp, only: proc0
    use redistribute, only: gather, scatter
    use dist_fn_arrays, only: gvmu
    use job_manage, only: time_message
    use stella_layouts, only: kxyz_lo, kxkyz_lo, vmu_lo
    use stella_transforms, only: transform_ky2y
    use zgrid, only: nzgrid
    use kt_grids, only: alpha_global
    use kt_grids, only: nakx, naky, ny
    use vpamu_grids, only: nvgrid, nmu

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: gout

    complex, dimension (:,:,:), allocatable :: g0v
    complex, dimension (:,:,:,:), allocatable :: g0x

    if (proc0) call time_message(.false.,time_gke(:,2),' Mirror advance')
    ! the mirror term is most complicated of all when doing full flux surface
    if (alpha_global) then
       allocate (g0v(-nvgrid:nvgrid,nmu,kxyz_lo%llim_proc:kxyz_lo%ulim_alloc))
       allocate (g0x(ny,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

       ! for upwinding, need to evaluate dh/dvpa in y-space
       ! first must take h(ky) and transform to h(y)
       call transform_ky2y (g, g0x)
       ! second, remap h so velocities are local
       call scatter (kxyz2vmu, g0x, g0v)
       ! next, take dh/dvpa and multiply with mirror coefficient in y-space.
       call get_dgdvpa_global (g0v)
       ! then take the results and remap again so y,kx,z local.
       call gather (kxyz2vmu, g0v, g0x)
       ! finally add the mirror term to the RHS of the GK eqn
       call add_mirror_term_global (g0x, gout)
    else
       allocate (g0v(-nvgrid:nvgrid,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       allocate (g0x(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

       ! get dg/dvpa and store in g0_vmu
       if (debug) write (*,*) 'dist_fn::advance_stella::get_dgdvpa'
       call get_dgdvpa (gvmu, g0v)
       if (debug) write (*,*) 'dist_fn::advance_stella::gather'
       ! swap layouts so that (z,kx,ky) are local
       call gather (kxkyz2vmu, g0v, g0x)
       ! get mirror term and add to source
       call add_mirror_term (g0x, gout)
    end if
    deallocate (g0x, g0v)
    if (proc0) call time_message(.false.,time_gke(:,2),' Mirror advance')

  end subroutine advance_mirror

  subroutine advance_wdrifty (g, gout)

    use mp, only: proc0
    use stella_layouts, only: vmu_lo
    use job_manage, only: time_message
    use stella_transforms, only: transform_ky2y
    use zgrid, only: nzgrid
    use kt_grids, only: nakx, naky, ny
    use kt_grids, only: alpha_global

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: gout

    complex, dimension (:,:,:,:), allocatable :: g0k, g0y

    ! alpha-component of magnetic drift (requires ky -> y)
    if (proc0) call time_message(.false.,time_gke(:,4),' dgdy advance')

    allocate (g0k(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    if (debug) write (*,*) 'dist_fn::solve_gke::get_dgdy'
    call get_dgdy (g, g0k)
    if (alpha_global) then
       allocate (g0y(ny,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       call transform_ky2y (g0k, g0y)
       call add_dgdy_term_global (g0y, gout)
       deallocate (g0y)
    else
       if (debug) write (*,*) 'dist_fn::solve_gke::add_dgdy_term'
       call add_dgdy_term (g0k, gout)
    end if
    deallocate (g0k)

    if (proc0) call time_message(.false.,time_gke(:,4),' dgdy advance')

  end subroutine advance_wdrifty

  subroutine advance_wdriftx (g, gout)

    use mp, only: proc0
    use stella_layouts, only: vmu_lo
    use job_manage, only: time_message
    use stella_transforms, only: transform_ky2y
    use zgrid, only: nzgrid
    use kt_grids, only: nakx, naky, ny
    use kt_grids, only: alpha_global

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: gout

    complex, dimension (:,:,:,:), allocatable :: g0k, g0y

    ! psi-component of magnetic drift (requires ky -> y)
    if (proc0) call time_message(.false.,time_gke(:,5),' dgdx advance')

    allocate (g0k(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    if (debug) write (*,*) 'dist_fn::solve_gke::get_dgdx'
    call get_dgdx (g, g0k)
    if (alpha_global) then
       allocate (g0y(ny,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       call transform_ky2y (g0k, g0y)
       call add_dgdx_term_global (g0y, gout)
       deallocate (g0y)
    else
       if (debug) write (*,*) 'dist_fn::solve_gke::add_dgdx_term'
       call add_dgdx_term (g0k, gout)
    end if
    deallocate (g0k)
    if (proc0) call time_message(.false.,time_gke(:,5),' dgdx advance')

  end subroutine advance_wdriftx

  subroutine get_dgdvpa (g, dgdv)

    use finite_differences, only: third_order_upwind
    use stella_layouts, only: kxkyz_lo, ig_idx
    use vpamu_grids, only: nvgrid, nmu, dvpa

    implicit none

    complex, dimension (-nvgrid:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (-nvgrid:,:,kxkyz_lo%llim_proc:), intent (out) :: dgdv

    integer :: ikxkyz, imu, ig

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       ig = ig_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          call third_order_upwind (-nvgrid,g(:,imu,ikxkyz),dvpa,mirror_sign(1,ig),dgdv(:,imu,ikxkyz))
       end do
    end do

  end subroutine get_dgdvpa

  subroutine get_dgdvpa_global (g)

    use finite_differences, only: third_order_upwind
    use stella_layouts, only: kxyz_lo, ig_idx, iy_idx
    use vpamu_grids, only: nvgrid, nmu, dvpa

    implicit none

    complex, dimension (-nvgrid:,:,kxyz_lo%llim_proc:), intent (in out) :: g

    integer :: ikxyz, imu, ig, iy
    complex, dimension (:), allocatable :: tmp

    allocate (tmp(-nvgrid:nvgrid))
    do ikxyz = kxyz_lo%llim_proc, kxyz_lo%ulim_proc
       ig = ig_idx(kxyz_lo,ikxyz)
       iy = iy_idx(kxyz_lo,ikxyz)
       do imu = 1, nmu
          call third_order_upwind (-nvgrid,g(:,imu,ikxyz),dvpa,mirror_sign(iy,ig),tmp)
          g(:,imu,ikxyz) = tmp
       end do
    end do
    deallocate (tmp)

  end subroutine get_dgdvpa_global

  subroutine add_mirror_term (g, src)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: imu_idx, is_idx
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: imu, is, ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       src(:,:,:,ivmu) = src(:,:,:,ivmu) + spread(spread(mirror(1,:,imu,is),1,naky),2,nakx)*g(:,:,:,ivmu)
    end do

  end subroutine add_mirror_term

  subroutine add_mirror_term_global (g, src)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: imu_idx, is_idx
    use zgrid, only: nzgrid
    use kt_grids, only: nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: imu, is, ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       src(:,:,:,ivmu) = src(:,:,:,ivmu) + spread(mirror(:,:,imu,is),2,nakx)*g(:,:,:,ivmu)
    end do

  end subroutine add_mirror_term_global

  subroutine get_dgdz (g, dgdz)

    use finite_differences, only: third_order_upwind_zed
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx
    use zgrid, only: nzgrid, delzed
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (out) :: dgdz

    integer :: ivmu, iseg, ie, iky, iv
    complex, dimension (2) :: gleft, gright

    ! FLAG -- assuming delta zed is equally spaced below!
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       do iky = 1, naky
          do ie = 1, neigen(iky)
             do iseg = 1, nsegments(ie,iky)
                ! if iseg,ie,iky corresponds to negative kx, no need to solve
                ! as it will be constrained by reality condition
                if (ikxmod(iseg,ie,iky) > nakx) cycle
                call fill_zed_ghost_zones (iseg, ie, iky, g(:,:,:,ivmu), gleft, gright)
                call third_order_upwind_zed (ig_low(iseg), iseg, nsegments(ie,iky), &
                     g(iky,ikxmod(iseg,ie,iky),ig_low(iseg):ig_up(iseg),ivmu), &
                     delzed(0), stream_sign(iv), gleft, gright, &
                     dgdz(iky,ikxmod(iseg,ie,iky),ig_low(iseg):ig_up(iseg),ivmu))
             end do
          end do
       end do
    end do

  end subroutine get_dgdz

  subroutine add_stream_term (g, src)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, is_idx
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: iv, is, ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       src(:,:,:,ivmu) = src(:,:,:,ivmu) + spread(spread(stream(:,iv,is),1,naky),2,nakx)*g(:,:,:,ivmu)
    end do

  end subroutine add_stream_term

  subroutine fill_zed_ghost_zones (iseg, ie, iky, g, gleft, gright)

    use zgrid, only: nzgrid
    use kt_grids, only: ntheta0, nakx, naky
    use kt_grids, only: aky

    implicit none

    integer, intent (in) :: iseg, ie, iky
    complex, dimension (:,:,-nzgrid:), intent (in) :: g
    complex, dimension (:), intent (out) :: gleft, gright

    integer :: ikxneg, ikyneg

    ! stream_sign > 0 --> stream speed < 0

    if (iseg == 1) then
       gleft = 0.0
    else
       ! if trying to connect to a segment that 
       ! corresponds to kx negative, use reality
       ! condition to connect to positive kx instead
       ! (with sign of ky flipped)
       if (ikxmod(iseg-1,ie,iky) > nakx) then
          ikxneg = ntheta0-ikxmod(iseg-1,ie,iky)+2
          if (abs(aky(iky)) < epsilon(0.)) then
             ikyneg = iky
          else
             ikyneg = naky-iky+2
          end if
          gleft = conjg(g(ikyneg,ikxneg,ig_up(iseg-1)-2:ig_up(iseg-1)-1))
       else
          gleft = g(iky,ikxmod(iseg-1,ie,iky),ig_up(iseg-1)-2:ig_up(iseg-1)-1)
       end if
    end if
    
    if (nsegments(ie,iky) > iseg) then
       ! if trying to connect to a segment that 
       ! corresponds to kx negative, use reality
       ! condition to connect to positive kx instead
       ! (with sign of ky flipped)
       if (ikxmod(iseg+1,ie,iky) > nakx) then
          ikxneg = ntheta0-ikxmod(iseg+1,ie,iky)+2
          if (abs(aky(iky)) < epsilon(0.)) then
             ikyneg = iky
          else
             ikyneg = naky-iky+2
          end if
          gright = conjg(g(ikyneg,ikxneg,ig_low(iseg+1)+1:ig_low(iseg+1)+2))
       else
          ! connect to segment with larger theta-theta0 (on right)
          gright = g(iky,ikxmod(iseg+1,ie,iky),ig_low(iseg+1)+1:ig_low(iseg+1)+2)
       end if
    else
       gright = 0.0
    end if
    
  end subroutine fill_zed_ghost_zones

  subroutine get_dgdy (g, dgdy)

    use constants, only: zi
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use kt_grids, only: nakx, aky

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (out) :: dgdy

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       dgdy(:,:,:,ivmu) = zi*spread(spread(aky,2,nakx),3,2*nzgrid+1)*g(:,:,:,ivmu)
    end do

  end subroutine get_dgdy

  subroutine add_dgdy_term (g, src)

    use dist_fn_arrays, only: wdrifty
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,ivmu) = src(:,:,:,ivmu) + spread(spread(wdrifty(1,:,ivmu),1,naky),2,nakx)*g(:,:,:,ivmu)
    end do

  end subroutine add_dgdy_term

  subroutine add_dgdy_term_global (g, src)

    use dist_fn_arrays, only: wdrifty
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use kt_grids, only: nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,ivmu) = src(:,:,:,ivmu) + spread(wdrifty(:,:,ivmu),2,nakx)*g(:,:,:,ivmu)
    end do

  end subroutine add_dgdy_term_global

  subroutine get_dgdx (g, dgdx)

    use constants, only: zi
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use kt_grids, only: naky, akx

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (out) :: dgdx

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       dgdx(:,:,:,ivmu) = zi*spread(spread(akx,1,naky),3,2*nzgrid+1)*g(:,:,:,ivmu)
    end do

  end subroutine get_dgdx

  subroutine add_dgdx_term (g, src)

    use dist_fn_arrays, only: wdriftx
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,ivmu) = src(:,:,:,ivmu) + spread(spread(wdriftx(1,:,ivmu),1,naky),2,nakx)*g(:,:,:,ivmu)
    end do

  end subroutine add_dgdx_term

  subroutine add_dgdx_term_global (g, src)

    use dist_fn_arrays, only: wdriftx
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use kt_grids, only: nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,ivmu) = src(:,:,:,ivmu) + spread(wdriftx(:,:,ivmu),2,nakx)*g(:,:,:,ivmu)
    end do

  end subroutine add_dgdx_term_global

  subroutine get_dchidy (phi, apar, dchidy)

    use constants, only: zi
    use dist_fn_arrays, only: aj0x
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: is_idx, iv_idx
    use run_parameters, only: fphi, fapar
    use species, only: spec
    use zgrid, only: nzgrid
    use vpamu_grids, only: vpa
    use kt_grids, only: nakx, aky

    implicit none

    complex, dimension (:,:,-nzgrid:), intent (in) :: phi, apar
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (out) :: dchidy

    integer :: ivmu, iv, is

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       iv = iv_idx(vmu_lo,ivmu)
       dchidy(:,:,:,ivmu) = zi*spread(spread(aky,2,nakx),3,2*nzgrid+1)*aj0x(:,:,:,ivmu) &
            * ( fphi*phi - fapar*vpa(iv)*spec(is)%stm*apar )
    end do

  end subroutine get_dchidy

  subroutine add_wstar_term (g, src)

    use dist_fn_arrays, only: wstar
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,ivmu) = src(:,:,:,ivmu) + spread(spread(wstar(:,ivmu),1,naky),2,nakx)*g(:,:,:,ivmu)
    end do

  end subroutine add_wstar_term

  subroutine finish_dist_fn

    use stella_transforms, only: finish_transforms
    use kt_grids, only: alpha_global
    
    implicit none

    dist_fn_initialized = .false.
    readinit = .false.
    gxyz_initialized = .false.

    if (alpha_global) call finish_transforms
    call finish_redistribute
    call finish_stream
    call finish_mirror
    call finish_bessel
    call finish_wdrift
    call finish_wstar
    call finish_kperp2
    call deallocate_arrays

  end subroutine finish_dist_fn

  subroutine finish_redistribute

    implicit none

    redistinit = .false.

  end subroutine finish_redistribute

  subroutine finish_stream

    implicit none

    if (allocated(stream)) deallocate (stream)
    if (allocated(stream_sign)) deallocate (stream_sign)

    streaminit = .false.

  end subroutine finish_stream

  subroutine finish_mirror

    implicit none

    if (allocated(mirror)) deallocate (mirror)
    if (allocated(mirror_sign)) deallocate (mirror_sign)

    mirrorinit = .false.

  end subroutine finish_mirror

  subroutine finish_wdrift

    use dist_fn_arrays, only: wdriftx, wdrifty

    implicit none

    if (allocated(wdriftx)) deallocate (wdriftx)
    if (allocated(wdrifty)) deallocate (wdrifty)

    wdriftinit = .false.

  end subroutine finish_wdrift

  subroutine finish_wstar

    use dist_fn_arrays, only: wstar

    implicit none

    if (allocated(wstar)) deallocate (wstar)

    wstarinit = .false.

  end subroutine finish_wstar

  subroutine finish_kperp2

    implicit none

    if (allocated(kperp2)) deallocate (kperp2)

    kp2init = .false.

  end subroutine finish_kperp2

  subroutine deallocate_arrays

    use dist_fn_arrays, only: gnew, gold
    use dist_fn_arrays, only: g1, g2

    implicit none

    if (allocated(gnew)) deallocate (gnew)
    if (allocated(gold)) deallocate (gold)
    if (allocated(g1)) deallocate (g1)
    if (allocated(g2)) deallocate (g2)

  end subroutine deallocate_arrays

  subroutine finish_bessel

    use dist_fn_arrays, only: aj0v, aj0x

    implicit none

    if (allocated(aj0v)) deallocate (aj0v)
    if (allocated(aj0x)) deallocate (aj0x)

    bessinit = .false.

  end subroutine finish_bessel

end module dist_fn
