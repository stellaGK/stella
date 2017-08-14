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
       boundary_option_alternate_zero = 3, &
       boundary_option_linked = 4

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
  integer, dimension (:), allocatable :: mirror_sign
  real, dimension (:,:,:), allocatable :: mirror

  ! needed for parallel streaming term
  integer, dimension (:), allocatable :: stream_sign
  real, dimension (:,:,:), allocatable :: stream

  ! these arrays needed to keep track of connections between different
  ! 2pi segments
  integer :: nseg_max, neigen_max
  integer, dimension (:), allocatable :: neigen
  integer, dimension (:), allocatable :: ig_low, ig_mid, ig_up, it_shift_left
  integer, dimension (:,:), allocatable :: nsegments, ir_up, it_shift
  integer, dimension (:,:,:), allocatable :: itmod
  logical, dimension (:), allocatable :: periodic

  ! needed for timing various pieces of gke solve
  real, dimension (2,6) :: time_gke

  type (redist_type) :: vmu2y
  type (redist_type) :: vmu2zkxky

  logical :: debug = .false.

contains

  subroutine init_get_fields

    use mp, only: sum_allreduce, proc0
    use stella_layouts, only: gvmu_lo
    use stella_layouts, onlY: ig_idx, it_idx, ik_idx, is_idx
    use dist_fn_arrays, only: aj0v
    use run_parameters, only: fphi, fapar
    use run_parameters, only: tite, nine, beta
    use species, only: spec, has_electron_species
    use theta_grid, only: ntgrid, dl_over_b
    use vpamu_grids, only: nvpa, nvgrid, nmu
    use vpamu_grids, only: vpa
    use vpamu_grids, only: anon, integrate_vmu
    use species, only: spec
    use kt_grids, only: naky, ntheta0, aky, akx

    implicit none

    integer :: ivmu, ig, it, ik, is
    real :: tmp, wgt
    real, dimension (:,:), allocatable :: g0

    if (get_fields_initialized) return
    get_fields_initialized = .true.

    debug = debug .and. proc0

    if (.not.allocated(gamtot)) allocate (gamtot(naky,ntheta0,-ntgrid:ntgrid)) ; gamtot = 0.
    if (.not.allocated(gamtot3)) then
       if (.not.has_electron_species(spec) &
            .and. adiabatic_option_switch==adiabatic_option_fieldlineavg) then
          allocate (gamtot3(ntheta0,-ntgrid:ntgrid)) ; gamtot3 = 0.
       else
          allocate (gamtot3(1,1)) ; gamtot3 = 0.
       end if
    end if
    if (.not.allocated(apar_denom)) allocate (apar_denom(naky,ntheta0,-ntgrid:ntgrid)) ; apar_denom = 0.

    if (fphi > epsilon(0.0)) then
       allocate (g0(-nvgrid:nvgrid,nmu))
       do ivmu = gvmu_lo%llim_proc, gvmu_lo%ulim_proc
          ik = ik_idx(gvmu_lo,ivmu)
          it = it_idx(gvmu_lo,ivmu)
          ig = ig_idx(gvmu_lo,ivmu)
          is = is_idx(gvmu_lo,ivmu)
          g0 = spread((1.0 - aj0v(:,ivmu)**2),1,nvpa)*anon(ig,:,:)
          wgt = spec(is)%z*spec(is)%z*spec(is)%dens/spec(is)%temp
          call integrate_vmu (g0, ig, tmp)
          gamtot(ik,it,ig) = gamtot(ik,it,ig) + tmp*wgt
       end do
       call sum_allreduce (gamtot)

       deallocate (g0)

       if (.not.has_electron_species(spec)) then
          gamtot = gamtot + tite/nine
          if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
             if (aky(1) < epsilon(0.)) then
                do it = 1, ntheta0
                   ! avoid divide by zero for kx=ky=0 mode,
                   ! which we do not need anyway
                   if (abs(akx(it)) < epsilon(0.)) cycle
                   tmp = nine/tite-sum(dl_over_b/gamtot(1,it,:))
                   gamtot3(it,:) = 1./(gamtot(1,it,:)*tmp)
                end do
             end if
          end if
       end if
    end if

    if (fapar > epsilon(0.)) then
       allocate (g0(-nvgrid:nvgrid,nmu))
       do ivmu = gvmu_lo%llim_proc, gvmu_lo%ulim_proc
          ik = ik_idx(gvmu_lo,ivmu)
          it = it_idx(gvmu_lo,ivmu)
          ig = ig_idx(gvmu_lo,ivmu)
          is = is_idx(gvmu_lo,ivmu)
          g0 = spread(vpa**2,2,nmu)*spread(aj0v(:,ivmu)**2,1,nvpa)*anon(ig,:,:)
          wgt = 2.0*beta*spec(is)%z*spec(is)%z*spec(is)%dens/spec(is)%mass
          call integrate_vmu (g0, ig, tmp)
          apar_denom(ik,it,ig) = apar_denom(ik,it,ig) + tmp*wgt
       end do
       call sum_allreduce (apar_denom)
       apar_denom = apar_denom + kperp2

       deallocate (g0)
    end if

  end subroutine init_get_fields

  subroutine get_fields (g, phi, apar)

    use mp, only: sum_allreduce
    use stella_layouts, only: gvmu_lo
    use stella_layouts, only: ig_idx, it_idx, ik_idx, is_idx
    use dist_fn_arrays, only: aj0v
    use run_parameters, only: fphi, fapar
    use run_parameters, only: beta
    use theta_grid, only: ntgrid, dl_over_b
    use vpamu_grids, only: nvgrid, nvpa, nmu
    use vpamu_grids, only: vpa
    use vpamu_grids, only: integrate_vmu
    use kt_grids, only: ntheta0, aky
    use species, only: spec, has_electron_species

    implicit none
    
    complex, dimension (-nvgrid:,:,gvmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-ntgrid:), intent (out) :: phi, apar

    real :: wgt
    complex, dimension (:,:), allocatable :: g0
    integer :: ivmu, ig, it, ik, is
    complex :: tmp

    phi = 0.
    if (fphi > epsilon(0.0)) then
       allocate (g0(-nvgrid:nvgrid,nmu))
       do ivmu = gvmu_lo%llim_proc, gvmu_lo%ulim_proc
          ig = ig_idx(gvmu_lo,ivmu)
          it = it_idx(gvmu_lo,ivmu)
          ik = ik_idx(gvmu_lo,ivmu)
          is = is_idx(gvmu_lo,ivmu)
          g0 = spread(aj0v(:,ivmu),1,nvpa)*g(:,:,ivmu)
          wgt = spec(is)%z*spec(is)%dens
          call integrate_vmu (g0, ig, tmp)
          phi(ik,it,ig) = phi(ik,it,ig) + wgt*tmp
       end do
       call sum_allreduce (phi)
       phi = phi/gamtot

       if (.not.has_electron_species(spec) .and. &
            adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          if (aky(1) < epsilon(0.)) then
             do it = 1, ntheta0
                tmp = sum(dl_over_b*phi(1,it,:))
                phi(1,it,:) = phi(1,it,:) + tmp*gamtot3(it,:)
             end do
          end if
       end if

       deallocate (g0)
    end if

    apar = 0.
    if (fapar > epsilon(0.0)) then
       allocate (g0(-nvgrid:nvgrid,nmu))
       do ivmu = gvmu_lo%llim_proc, gvmu_lo%ulim_proc
          ig = ig_idx(gvmu_lo,ivmu)
          it = it_idx(gvmu_lo,ivmu)
          ik = ik_idx(gvmu_lo,ivmu)
          is = is_idx(gvmu_lo,ivmu)
          g0 = spread(aj0v(:,ivmu),1,nvpa)*spread(vpa,2,nmu)*g(:,:,ivmu)
          wgt = 2.0*beta*spec(is)%z*spec(is)%dens*spec(is)%stm
          call integrate_vmu (g0, ig, tmp)
          apar(ik,it,ig) = apar(ik,it,ig) + tmp*wgt
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
    call gather (vmu2zkxky, gvmu, gnew)
    gold = gnew

  end subroutine init_gxyz

  subroutine init_dist_fn

    use stella_layouts, only: init_stella_layouts, init_dist_fn_layouts
    use species, only: init_species
    use species, only: nspec
    use theta_grid, only: init_theta_grid
    use theta_grid, only: ntgrid
    use kt_grids, only: init_kt_grids
    use kt_grids, only: naky, ntheta0, ny
    use vpamu_grids, only: init_vpamu_grids
    use vpamu_grids, only: nvgrid, nmu
    use run_parameters, only: init_run_parameters

    implicit none

    if (dist_fn_initialized) return
    dist_fn_initialized = .true.

    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_stella_layouts'
    call init_stella_layouts
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_species'
    call init_species
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_theta_grid'
    call init_theta_grid
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_kt_grids'
    call init_kt_grids
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_vpamu_grids'
    call init_vpamu_grids
    if (debug) write (*,*) 'dist_fn::init_dist_fn::read_parameters'
    call read_parameters
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_run_parameters'
    call init_run_parameters
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_dist_fn_layouts'
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nvgrid, nmu, nspec, ny)
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
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_redistribute'
    call init_redistribute
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_cfl'
    call init_cfl

  end subroutine init_dist_fn

  subroutine read_parameters

    use file_utils, only: input_unit, error_unit, input_unit_exist
    use theta_grid, only: shat
    use text_options, only: text_option, get_option_value
    use species, only: nspec
    use mp, only: proc0, broadcast

    implicit none

    logical :: dfexist

    type (text_option), dimension (8), parameter :: boundaryopts = &
         (/ text_option('default', boundary_option_zero), &
            text_option('zero', boundary_option_zero), &
            text_option('unconnected', boundary_option_zero), &
            text_option('self-periodic', boundary_option_self_periodic), &
            text_option('periodic', boundary_option_self_periodic), &
            text_option('kperiod=1', boundary_option_self_periodic), &
            text_option('linked', boundary_option_linked), &
            text_option('alternate-zero', boundary_option_alternate_zero) /)
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

       if(abs(shat) <=  1.e-5) boundary_option = 'periodic'

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

    use theta_grid, only: ntgrid, gds2, gds21, gds22, shat
    use kt_grids, only: naky, ntheta0, aky, theta0, akx

    implicit none

    integer :: ik, it

    if (kp2init) return
    kp2init = .true.

    allocate (kperp2(naky,ntheta0,-ntgrid:ntgrid))
    do ik = 1, naky
       if (aky(ik) == 0.0) then
          do it = 1, ntheta0
             kperp2(ik,it,:) = akx(it)*akx(it)*gds22/(shat*shat)
          end do
       else
          do it = 1, ntheta0
             kperp2(ik,it,:) = aky(ik)*aky(ik) &
                  *(gds2 + 2.0*theta0(ik,it)*gds21 &
                  + theta0(ik,it)*theta0(ik,it)*gds22)
          end do
       end if
    end do

  end subroutine init_kperp2

  subroutine init_wdrift

    use dist_fn_arrays, only: wdriftx, wdrifty
    use stella_layouts, only: gxyz_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use species, only: spec
    use theta_grid, only: ntgrid
    use theta_grid, only: cvdrift, gbdrift
    use theta_grid, only: cvdrift0, gbdrift0, shat
    use vpamu_grids, only: vpa, vperp2

    implicit none

    integer :: ixyz, iv, imu, is

    if (wdriftinit) return
    wdriftinit = .true.

    if (.not.allocated(wdriftx)) &
         allocate (wdriftx(-ntgrid:ntgrid,gxyz_lo%llim_proc:gxyz_lo%ulim_alloc))
    if (.not.allocated(wdrifty)) &
         allocate (wdrifty(-ntgrid:ntgrid,gxyz_lo%llim_proc:gxyz_lo%ulim_alloc))

    ! FLAG -- need to deal with shat=0 case.  ideally move away from q as x-coordinate
    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       iv = iv_idx(gxyz_lo,ixyz)
       imu = imu_idx(gxyz_lo,ixyz)
       is = is_idx(gxyz_lo,ixyz)
       wdrifty(:,ixyz) = -ydriftknob*0.5*code_dt*spec(is)%tz &
            * (cvdrift*vpa(iv)**2 + gbdrift*0.5*vperp2(:,imu))
       wdriftx(:,ixyz) = -xdriftknob*0.5*code_dt*spec(is)%tz/shat &
            * (cvdrift0*vpa(iv)**2 + gbdrift0*0.5*vperp2(:,imu))
    end do

  end subroutine init_wdrift

  subroutine init_wstar

    use dist_fn_arrays, only: wstar
    use stella_layouts, only: gxyz_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use species, only: spec
    use theta_grid, only: ntgrid
    use vpamu_grids, only: energy, anon

    implicit none

    integer :: is, imu, iv, ixyz

    if (wstarinit) return
    wstarinit = .true.

    if (.not.allocated(wstar)) &
         allocate (wstar(-ntgrid:ntgrid,gxyz_lo%llim_proc:gxyz_lo%ulim_alloc)) ; wstar = 0.0

    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       is = is_idx(gxyz_lo,ixyz)
       imu = imu_idx(gxyz_lo,ixyz)
       iv = iv_idx(gxyz_lo,ixyz)
       wstar(:,ixyz) = wstarknob*0.5*code_dt*anon(:,iv,imu) &
            * (spec(is)%fprim+spec(is)%tprim*(energy(:,iv,imu)-1.5))
    end do

  end subroutine init_wstar

  subroutine allocate_arrays

    use stella_layouts, only: gvmu_lo, gxyz_lo
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0
    use vpamu_grids, only: nvgrid, nmu
    use dist_fn_arrays, only: gnew, gold
    use dist_fn_arrays, only: g1, g2
    use dist_fn_arrays, only: gvmu

    implicit none

    if (.not.allocated(gnew)) &
         allocate (gnew(naky,ntheta0,-ntgrid:ntgrid,gxyz_lo%llim_proc:gxyz_lo%ulim_alloc))
    gnew = 0.
    if (.not.allocated(gold)) &
         allocate (gold(naky,ntheta0,-ntgrid:ntgrid,gxyz_lo%llim_proc:gxyz_lo%ulim_alloc))
    gold = 0.
    if (.not.allocated(g1)) &
         allocate (g1(naky,ntheta0,-ntgrid:ntgrid,gxyz_lo%llim_proc:gxyz_lo%ulim_alloc))
    g1 = 0.
    if (.not.allocated(g2)) &
         allocate (g2(naky,ntheta0,-ntgrid:ntgrid,gxyz_lo%llim_proc:gxyz_lo%ulim_alloc))
    g2 = 0.
    if (.not.allocated(gvmu)) &
         allocate (gvmu(-nvgrid:nvgrid,nmu,gvmu_lo%llim_proc:gvmu_lo%ulim_alloc))
    gvmu = 0.

  end subroutine allocate_arrays

  subroutine init_connections

    use theta_grid, only: nperiod, ntgrid, ntheta
    use kt_grids, only: ntheta0, jtwist_out, naky, aky
    use species, only: nspec
    use vpamu_grids, only: nmu, nvgrid

    implicit none

    integer :: iseg, ik, ie, ntg, it

    ntg = ntheta/2

    if (.not. allocated(neigen)) allocate (neigen(naky))
    if (.not. allocated(periodic)) allocate (periodic(naky)) ; periodic = .false.

    if (boundary_option_switch==boundary_option_self_periodic) then
       periodic = .true.
    else
       do ik = 1, naky
          if (aky(ik) < epsilon(0.0)) periodic(ik) = .true.
       end do
    end if

    select case (boundary_option_switch)
    case (boundary_option_linked)

       ik = 1
       neigen(ik) = ntheta0
       if (naky > 1) then
          do ik = 2, naky
             ! must link different kx values at theta = +/- pi
             ! neigen is the number of independent eigenfunctions along the field line
             neigen(ik) = min((ik-1)*jtwist_out,ntheta0)
          end do
       end if

       neigen_max = maxval(neigen)

       if (.not. allocated(it_shift_left)) then
          allocate (it_shift_left(neigen_max))
          allocate (it_shift(ntheta0,naky)) ; it_shift = 0
       end if

       ! figure out how much to shift it by to get to
       ! the left-most (theta-theta0) in each set of connected 2pi segments
       ! note that theta0 goes from 0 to theta0_max and then from theta0_min back
       ! to -dtheta0
       do it = 1, neigen_max
          ! first ntheta0/2+1 theta0s are 0 and all positive theta0 values
          ! remainder are negative theta0s
          ! note that ntheta0 is always positive for box
          if (it <= ntheta0/2+1) then
             it_shift_left(it) = ntheta0/2-2*it+2
          else
             it_shift_left(it) = 3*(ntheta0/2)-2*it+3
          end if
       end do
       
       do ik = 1, naky
          ! it_shift is how much to shift each it by to connect
          ! to the next theta0 (from most positive to most negative)
          do it = 1, ntheta0
             ! if theta0 is negative, then shifting to more negative
             ! theta0 corresponds to decreasing it
             if (it > ntheta0/2+1) then
                ! if theta0 is sufficiently negative, it has no
                ! more negative theta0 with which it can connect
                if (it-neigen(ik) >= ntheta0/2+1) then
                   it_shift(it,ik) = -neigen(ik)
                end if
                ! theta0 is positive
             else
                ! if theta0 is sufficiently positive, shifting to more
                ! negative theta0 corresponds to decreasing it
                if (it-neigen(ik) > 0) then
                   it_shift(it,ik) = -neigen(ik)
                   ! if a positive theta0 connects to a negative theta0
                   ! must do more complicated mapping of it
                else if (it-neigen(ik)+ntheta0 >= 2+ntheta0/2) then
                   it_shift(it,ik) = ntheta0 - neigen(ik)
                end if
             end if
          end do
       end do

       if (.not. allocated(nsegments)) then
          allocate (nsegments(neigen_max,naky))
       end if

       do ik = 1, naky
          if (neigen(ik) == 0) then
             nsegments(:,ik) = 1
          else
             nsegments(:,ik) = (ntheta0-1)/neigen(ik)

             do ie = 1, mod(ntheta0-1,neigen(ik))+1
                nsegments(ie,ik) = nsegments(ie,ik) + 1
             end do
          end if
       end do

       nseg_max = maxval(nsegments)

       if (.not. allocated(ig_low)) then
          allocate (ig_low(nseg_max)) ; ig_low = -ntgrid
          allocate (ig_mid(nseg_max)) ; ig_mid = 0
          allocate (ig_up(nseg_max)) ; ig_up = ntgrid
       end if
       
    case default
       
       neigen = ntheta0 ; neigen_max = ntheta0
       
       if (.not. allocated(it_shift_left)) then
          allocate (it_shift_left(neigen_max))
          allocate (it_shift(ntheta0,naky))
       end if
       it_shift = 0 ; it_shift_left = 0
       
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
          ig_low(iseg) = -ntgrid + (iseg-1)*ntheta
          ig_mid(iseg) = ig_low(iseg) + ntheta/2
          ig_up(iseg) = ig_low(iseg) + ntheta
       end do

    end select

    if (.not. allocated(ir_up)) allocate (ir_up(neigen_max,naky))
    ir_up = ntg*nsegments+1
    
    if (.not. allocated(itmod)) allocate (itmod(nseg_max,neigen_max,naky))
    do ik = 1, naky
       if (periodic(ik)) ir_up(:,ik) = ir_up(:,ik) + 2*nvgrid
       
       ! only do the following once for each independent set of theta0s
       ! the assumption here is that all kx are on processor and sequential
       do it = 1, neigen(ik)
          
          ! remap to start at theta0 = -theta0_max for this set of connected theta0s
          iseg = 1
          itmod(iseg,it,ik) = it + it_shift_left(it)
          if (nsegments(it,ik) > 1) then
             do iseg = 2, nsegments(it,ik)
                itmod(iseg,it,ik) = itmod(iseg-1,it,ik) + it_shift(itmod(iseg-1,it,ik),ik)
             end do
          end if
       end do
    end do

  end subroutine init_connections

  subroutine init_vperp2

    use theta_grid, only: ntgrid, bmag
    use vpamu_grids, only: vperp2, energy, anon
    use vpamu_grids, only: vpa, mu
    use vpamu_grids, only: nmu, nvgrid

    implicit none

    integer :: iv
    
    if (.not.allocated(vperp2)) allocate (vperp2(-ntgrid:ntgrid,nmu)) ; vperp2 = 0.
    if (.not.allocated(energy)) allocate (energy(-ntgrid:ntgrid,-nvgrid:nvgrid,nmu)) ; energy = 0.
    if (.not.allocated(anon)) allocate (anon(-ntgrid:ntgrid,-nvgrid:nvgrid,nmu)) ; anon = 0.

    vperp2 = 2.0*spread(mu,1,2*ntgrid+1)*spread(bmag,2,nmu)

    do iv = -nvgrid, nvgrid
       energy(:,iv,:) = vpa(iv)**2 + 2.0*spread(mu,1,2*ntgrid+1)*spread(bmag,2,nmu)
       anon(:,iv,:) = exp(-energy(:,iv,:))
    end do

  end subroutine init_vperp2

  subroutine init_bessel

    use dist_fn_arrays, only: aj0v, aj0x
    use species, only: spec, nspec
    use theta_grid, only: bmag, ntgrid
    use vpamu_grids, only: vperp2, nmu
    use kt_grids, only: naky, ntheta0
    use stella_layouts, only: gvmu_lo, gxyz_lo
    use stella_layouts, only: ik_idx, it_idx, ig_idx, is_idx, imu_idx
    use spfunc, only: j0!, j1

    implicit none

    integer :: ig, ik, it, imu, is
    integer :: ivmu, ixyz
    real :: arg

    if (bessinit) return
    bessinit = .true.

    call init_kperp2

    allocate (aj0v(nmu,gvmu_lo%llim_proc:gvmu_lo%ulim_alloc)) ; aj0v = 0.
    allocate (aj0x(naky,ntheta0,-ntgrid:ntgrid,gxyz_lo%llim_proc:gxyz_lo%ulim_alloc)) ; aj0x = 0.
!    allocate (aj1(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; aj1 = 0.
!    allocate (aj2(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc)) ; aj2 = 0.

    do ivmu = gvmu_lo%llim_proc, gvmu_lo%ulim_proc
       ik = ik_idx(gvmu_lo,ivmu)
       it = it_idx(gvmu_lo,ivmu)
       ig = ig_idx(gvmu_lo,ivmu)
       is = is_idx(gvmu_lo,ivmu)
       do imu = 1, nmu
          arg = spec(is)%smz*sqrt(vperp2(ig,imu)*kperp2(ik,it,ig))/bmag(ig)
          aj0v(imu,ivmu) = j0(arg)
             ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
!             aj1(ig,ivmu) = j1(arg)
!             aj2(ig,ivmu) = 2.0*aj1(ig,ivmu)-aj0(ig,ivmu)
       end do
    end do

    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       is = is_idx(gxyz_lo,ixyz)
       imu = imu_idx(gxyz_lo,ixyz)
       do ig = -ntgrid, ntgrid
          do it = 1, ntheta0
             do ik = 1, naky
                arg = spec(is)%smz*sqrt(vperp2(ig,imu)*kperp2(ik,it,ig))/bmag(ig)
                aj0x(ik,it,ig,ixyz) = j0(arg)
             ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
!             aj1(ig,ixyz) = j1(arg)
!             aj2(ig,ixyz) = 2.0*aj1(ig,ixyz)-aj0(ig,ixyz)
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
    use theta_grid, only: ntgrid
    use theta_grid, only: dbdthet, gradpar

    implicit none

    integer :: ig

    if (mirrorinit) return
    mirrorinit = .true.

    if (.not.allocated(mirror)) allocate (mirror(-ntgrid:ntgrid,nmu,nspec)) ; mirror = 0.
    if (.not.allocated(mirror_sign)) allocate (mirror_sign(-ntgrid:ntgrid)) ; mirror_sign = 0

    mirror = mirrorknob*code_dt*spread(spread(spec%stm,1,2*ntgrid+1),2,nmu) &
         *spread(spread(mu,1,2*ntgrid+1)*spread(dbdthet*gradpar,2,nmu),3,nspec)
    ! mirror_sign set to +/- 1 depending on the sign of the mirror term.
    ! NB: mirror_sign = -1 corresponds to positive advection velocity
    do ig = -ntgrid, ntgrid
       mirror_sign(ig) = int(sign(1.0,mirror(ig,1,1)))
    end do

  end subroutine init_mirror

  subroutine init_parstream

    use stella_time, only: code_dt
    use species, only: spec, nspec
    use vpamu_grids, only: nvgrid
    use vpamu_grids, only: vpa, nvpa
    use theta_grid, only: ntgrid
    use theta_grid, only: gradpar

    implicit none

    integer :: iv

    if (streaminit) return
    streaminit = .true.

    if (.not.allocated(stream)) allocate (stream(-ntgrid:ntgrid,-nvgrid:nvgrid,nspec)) ; stream = 0.
    if (.not.allocated(stream_sign)) allocate (stream_sign(-nvgrid:nvgrid)) ; stream_sign = 0

    ! sign of stream corresponds to appearing on RHS of GK equation
    stream = -streamknob*code_dt*spread(spread(spec%stm,1,2*ntgrid+1),2,nvpa) &
         * spread(spread(vpa,1,2*ntgrid+1)*spread(gradpar,2,nvpa),3,nspec)

    ! stream_sign set to +/- 1 depending on the sign of the parallel streaming term.
    ! NB: stream_sign = -1 corresponds to positive advection velocity
    do iv = -nvgrid, nvgrid
       stream_sign(iv) = int(sign(1.0,stream(0,iv,1)))
    end do

  end subroutine init_parstream

  subroutine init_redistribute

    implicit none

    if (redistinit) return
    redistinit = .true.

    if (debug) write (*,*) 'dist_fn::init_redistribute::init_gvmu_to_gy_redistribute'
    call init_gvmu_to_gy_redistribute
    if (debug) write (*,*) 'dist_fn::init_redistribute::init_gvmu_to_gzkxky_redistribute'
    call init_gvmu_to_gzkxky_redistribute

  end subroutine init_redistribute

  subroutine init_gvmu_to_gy_redistribute

    use mp, only: nproc
    use stella_layouts, only: gvmu_lo, gy_lo
    use stella_layouts, only: vmuidx2yidx
    use stella_layouts, only: idx_local, proc_id
    use redistribute, only: index_list_type, init_redist
    use redistribute, only: delete_list, set_redist_character_type
    use vpamu_grids, only: nvgrid, nmu

    implicit none

    type (index_list_type), dimension (0:nproc-1) :: to_list, from_list
    integer, dimension (0:nproc-1) :: nn_to, nn_from
    integer, dimension (3) :: from_low, from_high
    integer, dimension (2) :: to_high, to_low
    integer :: ivmu, iy
    integer :: iv, imu, ik
    integer :: ip, n
    logical :: initialized = .false.

    if (initialized) return
    initialized = .true.

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do ivmu = gvmu_lo%llim_world, gvmu_lo%ulim_world
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             call vmuidx2yidx (iv, imu, ivmu, gvmu_lo, gy_lo, ik, iy)
             if (idx_local(gvmu_lo,ivmu)) &
                  nn_from(proc_id(gy_lo,iy)) = nn_from(proc_id(gy_lo,iy)) + 1
             if (idx_local(gy_lo,iy)) &
                  nn_to(proc_id(gvmu_lo,ivmu)) = nn_to(proc_id(gvmu_lo,ivmu)) + 1
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
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    nn_to = 0
    nn_from = 0

    ! loop over all vmu indices, find corresponding y indices
    do ivmu = gvmu_lo%llim_world, gvmu_lo%ulim_world
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             ! obtain corresponding y indices
             call vmuidx2yidx (iv, imu, ivmu, gvmu_lo, gy_lo, ik, iy)
             ! if vmu index local, set:
             ! ip = corresponding y processor
             ! from_list%first-third arrays = iv,imu,ivmu  (ie vmu indices)
             ! later will send from_list to proc ip
             if (idx_local(gvmu_lo,ivmu)) then
                ip = proc_id(gy_lo,iy)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n) = iv
                from_list(ip)%second(n) = imu
                from_list(ip)%third(n) = ivmu
             end if
             ! if y index local, set ip to corresponding vmu processor
             ! set to_list%first,second arrays = ik,iy  (ie y indices)
             ! will receive to_list from ip
             if (idx_local(gy_lo,iy)) then
                ip = proc_id(gvmu_lo,ivmu)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n) = ik
                to_list(ip)%second(n) = iy
             end if
          end do
       end do
    end do

    from_low(1) = -nvgrid
    from_low(2) = 1
    from_low(3) = gvmu_lo%llim_proc

    from_high(1) = nvgrid
    from_high(2) = nmu
    from_high(3) = gvmu_lo%ulim_alloc

    to_low(1) = 1
    to_low(2) = gy_lo%llim_proc

    to_high(1) = gy_lo%naky
    to_high(2) = gy_lo%ulim_alloc

    call set_redist_character_type (vmu2y, 'vmu2y')

    call init_redist (vmu2y, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

  end subroutine init_gvmu_to_gy_redistribute

  subroutine init_gvmu_to_gzkxky_redistribute

    use mp, only: nproc
    use stella_layouts, only: gvmu_lo, gxyz_lo
    use stella_layouts, only: vmuidx2zkxkyidx
    use stella_layouts, only: idx_local, proc_id
    use redistribute, only: index_list_type, init_redist
    use redistribute, only: delete_list, set_redist_character_type
    use vpamu_grids, only: nvgrid, nmu
    use theta_grid, only: ntgrid

    implicit none

    type (index_list_type), dimension (0:nproc-1) :: to_list, from_list
    integer, dimension (0:nproc-1) :: nn_to, nn_from
    integer, dimension (3) :: from_low, from_high
    integer, dimension (4) :: to_high, to_low
    integer :: ivmu, ixyz
    integer :: iv, imu, ik, it, ig
    integer :: ip, n
    logical :: initialized = .false.

    if (initialized) return
    initialized = .true.

    ! count number of elements to be redistributed to/from each processor
    nn_to = 0
    nn_from = 0
    do ivmu = gvmu_lo%llim_world, gvmu_lo%ulim_world
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             call vmuidx2zkxkyidx (iv, imu, ivmu, gvmu_lo, gxyz_lo, ik, it, ig, ixyz)
             if (idx_local(gvmu_lo,ivmu)) &
                  nn_from(proc_id(gxyz_lo,ixyz)) = nn_from(proc_id(gxyz_lo,ixyz)) + 1
             if (idx_local(gxyz_lo,ixyz)) &
                  nn_to(proc_id(gvmu_lo,ivmu)) = nn_to(proc_id(gvmu_lo,ivmu)) + 1
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
    do ivmu = gvmu_lo%llim_world, gvmu_lo%ulim_world
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             ! obtain corresponding y indices
             call vmuidx2zkxkyidx (iv, imu, ivmu, gvmu_lo, gxyz_lo, ik, it, ig, ixyz)
             ! if vmu index local, set:
             ! ip = corresponding y processor
             ! from_list%first-third arrays = iv,imu,ivmu  (ie vmu indices)
             ! later will send from_list to proc ip
             if (idx_local(gvmu_lo,ivmu)) then
                ip = proc_id(gxyz_lo,ixyz)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n) = iv
                from_list(ip)%second(n) = imu
                from_list(ip)%third(n) = ivmu
             end if
             ! if y index local, set ip to corresponding vmu processor
             ! set to_list%first,second arrays = ik,iy  (ie y indices)
             ! will receive to_list from ip
             if (idx_local(gxyz_lo,ixyz)) then
                ip = proc_id(gvmu_lo,ivmu)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                to_list(ip)%first(n) = ik
                to_list(ip)%second(n) = it
                to_list(ip)%third(n) = ig
                to_list(ip)%fourth(n) = ixyz
             end if
          end do
       end do
    end do

    from_low(1) = -nvgrid
    from_low(2) = 1
    from_low(3) = gvmu_lo%llim_proc

    from_high(1) = nvgrid
    from_high(2) = nmu
    from_high(3) = gvmu_lo%ulim_alloc

    to_low(1) = 1
    to_low(2) = 1
    to_low(3) = -ntgrid
    to_low(4) = gxyz_lo%llim_proc

    to_high(1) = gxyz_lo%naky
    to_high(2) = gxyz_lo%ntheta0
    to_high(3) = gxyz_lo%ntheta
    to_high(4) = gxyz_lo%ulim_alloc

    call set_redist_character_type (vmu2zkxky, 'vmu2zkxky')

    call init_redist (vmu2zkxky, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

  end subroutine init_gvmu_to_gzkxky_redistribute

  subroutine init_cfl
    
    use mp, only: proc0
    use dist_fn_arrays, only: wdriftx, wdrifty
    use stella_time, only: cfl_dt, code_dt, write_dt
    use theta_grid, only: delthet
    use vpamu_grids, only: dvpa
    use kt_grids, only: akx, aky

    implicit none
    
    real :: cfl_dt_mirror, cfl_dt_stream
    real :: cfl_dt_wdriftx, cfl_dt_wdrifty
    real :: zero

    ! avoid divide by zero in cfl_dt terms below
    zero = 100.*epsilon(0.)

    cfl_dt_mirror = code_dt*dvpa/max(maxval(abs(mirror)),zero)
    ! FLAG -- assuming equal spacing in theta!
    cfl_dt_stream = code_dt*delthet(0)/max(maxval(abs(stream)),zero)
    cfl_dt_wdriftx = code_dt/max(maxval(akx)*maxval(abs(wdriftx)),zero)
    cfl_dt_wdrifty = code_dt/max(maxval(aky)*maxval(abs(wdrifty)),zero)
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

  subroutine solve_gke (gin, gout)

    use mp, only: proc0
    use job_manage, only: time_message
    use dist_fn_arrays, only: gvmu
    use dist_fn_arrays, only: g_adjust
    use fields_arrays, only: phi, apar
    use stella_layouts, only: gvmu_lo
    use stella_layouts, only: gxyz_lo
    use redistribute, only: gather, scatter
    use run_parameters, only: fphi, fapar
    use theta_grid, only: ntgrid
    use vpamu_grids, only: nvgrid, nmu
    use kt_grids, only: naky, ntheta0

    implicit none

    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in out) :: gin
    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (out) :: gout

    complex, dimension (:,:,:), allocatable :: g0_vmu
    complex, dimension (:,:,:,:), allocatable :: g0_zkxky

    gout = 0.

    if (.not.allocated(g0_vmu)) &
         allocate (g0_vmu(-nvgrid:nvgrid,nmu,gvmu_lo%llim_proc:gvmu_lo%ulim_alloc))
    if (.not.allocated(g0_zkxky)) &
         allocate (g0_zkxky(naky,ntheta0,-ntgrid:ntgrid,gxyz_lo%llim_proc:gxyz_lo%ulim_alloc))

    ! start with g in k-space and (ky,kx,z) local

    ! first gather (vpa,mu) onto processor for v-space operations
    if (proc0) call time_message(.false.,time_gke(:,1),' fields')
    if (debug) write (*,*) 'dist_fn::advance_stella::scatter'
    call scatter (vmu2zkxky, gin, gvmu)
    if (debug) write (*,*) 'dist_fn::advance_stella::get_fields'
    call get_fields (gvmu, phi, apar)
    if (proc0) call time_message(.false.,time_gke(:,1),' fields')

    ! switch from g to h
    call g_adjust (gin, phi, apar, fphi, fapar)
    call g_adjust (gvmu, phi, apar, fphi, fapar)

    if (proc0) call time_message(.false.,time_gke(:,2),' Mirror advance')
    if (debug) write (*,*) 'dist_fn::advance_stella::get_dgdvpa'
    ! get dg/dvpa and store in g0_vmu
    call get_dgdvpa (gvmu, g0_vmu)
    if (debug) write (*,*) 'dist_fn::advance_stella::gather'
    ! swap layouts so that (z,kx,ky) are local
    call gather (vmu2zkxky, g0_vmu, g0_zkxky)
    ! get mirror term and add to source
    call add_mirror_term (g0_zkxky, gout)
    if (proc0) call time_message(.false.,time_gke(:,2),' Mirror advance')

    ! loop over vpar, mu, species, z, and kx
!    do iglo = gy_lo%llim_proc, gy_lo%ulim_proc
       ! transform dg/dvpa from k_alpha space to alpha space
!       call transform_g (g0_ky(:,iglo), g0_y(:,iglo))
       ! get mirror term and add to source
!       call add_mirror_term (g0_y(:,iglo), source(...))
!    end do

!     ! operate on g(vpa,mu,spec,k_alpha,k_psi,theta) with collision operator
!     call get_collision_term (gold, g0_vmus)
!     ! swap layouts so that k-space is local and (z,v)-space spread over procs
!     call remap_g (vmus_to_ky, g0_vmus, g0_ky)
!     ! loop over vpar, mu, species, z, and kx
!     do iglo = gy_lo%llim_proc, gy_lo%ulim_proc
!        ! transform collision term from k_alpha space to alpha space
!        call transform_g (g0_ky(:,iglo), g0_y(:,iglo))
!        ! add collision term to source
!        call add_collision_term (g0_y(:,iglo), source(...))
!     end do

    if (allocated(g0_vmu)) deallocate (g0_vmu)

    if (proc0) call time_message(.false.,time_gke(:,3),' Stream advance')
    if (debug) write (*,*) 'dist_fn::solve_gke::get_dgdz'
    call get_dgdz (gin, g0_zkxky)
    if (debug) write (*,*) 'dist_fn::solve_gke::add_stream_term'
    call add_stream_term (g0_zkxky, gout)
    if (proc0) call time_message(.false.,time_gke(:,3),' Stream advance')

    if (proc0) call time_message(.false.,time_gke(:,4),' dgdy advance')
    if (debug) write (*,*) 'dist_fn::solve_gke::get_dgdy'
    call get_dgdy (gin, g0_zkxky)
    if (debug) write (*,*) 'dist_fn::solve_gke::add_dgdy_term'
    call add_dgdy_term (g0_zkxky, gout)
    if (proc0) call time_message(.false.,time_gke(:,4),' dgdy advance')

    if (proc0) call time_message(.false.,time_gke(:,5),' dgdx advance')
    if (debug) write (*,*) 'dist_fn::solve_gke::get_dgdx'
    call get_dgdx (gin, g0_zkxky)
    if (debug) write (*,*) 'dist_fn::solve_gke::add_dgdx_term'
    call add_dgdx_term (g0_zkxky, gout)
    if (proc0) call time_message(.false.,time_gke(:,5),' dgdx advance')

    if (proc0) call time_message(.false.,time_gke(:,6),' wstar advance')
    if (debug) write (*,*) 'dist_fn::solve_gke::get_dchidy'
    call get_dchidy (phi, apar, g0_zkxky)
    if (debug) write (*,*) 'dist_fn::solve_gke::add_wstar_term'
    call add_wstar_term (g0_zkxky, gout)
    if (proc0) call time_message(.false.,time_gke(:,6),' wstar advance')

    if (allocated(g0_zkxky)) deallocate (g0_zkxky)

!     ! swap layouts so that k_alpha is local
!     call remap_g (vmus_to_zky, gold, g0_zky)
!     ! loop over vpar, mu, species, z, and kx
!     do iglo = gyz_lo%llim_proc, gyz_lo%ulim_proc
!        ! transform g from k_alpha space to alpha space
!        call transform_g (g0_zky(:,:,iglo), g0_zy(:,:,iglo))
!        ! add parallel streaming and psi-component of magnetic drift to source
!        call add_dgdz_dgdx_terms (g0_y(:,iglo), source(...))
!     end do
!     ! get dg/dalpha in k_alpha space
!     call get_dgdalpha (g0_ky)
!     ! loop over vpar, mu, species, z, and kx
!     do iglo = gy_lo%llim_proc, gy_lo%ulim_proc
!        ! transform g from k_alpha space to alpha space
!        call transform_g (g0_ky(:,iglo), g0_y(:,iglo))
!        ! add alpha-component of magnetic drift to source
!        call add_dgdy_term (g0_y(:,iglo), source(...))
!     end do


!     ! transform phi from k_alpha space to alpha space
!     call transform_phi (g0_zkykx, g0_zykx)
!     ! calculate contributions from source terms associated with parallel streaming
!     ! and psi-component of magnetic drift and add to source
!     call add_dphidz_dphidx_terms (g0_zykx, source)

!     ! get d<phi>/dalpha
!     call get_dphidalpha (g0_zkykx)
!     ! transform dphi/dalpha from k_alpha space to alpha space
!     call transform_phi (g0_zkykx, g0_zykx)
!     ! calculate contribution from omega_star and alpha-component
!     ! of magnetic drift and add to source
!     call add_dphidy_terms (g0_zykx, source)

!     ! 
!     ! multiple dg/dx and dphi/dy to get part of nonlinearity
!     call get_dgdx_dphidy_nonlinearity (g0_zyx, dphidy)
!     g0_zyx = g0_zyx * 

    ! switch from h back to g
    call g_adjust (gin, phi, apar, -fphi, -fapar)

  end subroutine solve_gke

  subroutine get_dgdvpa (g, dgdv)

    use finite_differences, only: third_order_upwind
    use stella_layouts, only: gvmu_lo, ig_idx
    use vpamu_grids, only: nvgrid, nmu, dvpa

    implicit none

    complex, dimension (-nvgrid:,:,gvmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (-nvgrid:,:,gvmu_lo%llim_proc:), intent (out) :: dgdv

    integer :: ivmu, imu, ig

    do ivmu = gvmu_lo%llim_proc, gvmu_lo%ulim_proc
       ig = ig_idx(gvmu_lo,ivmu)
       do imu = 1, nmu
          call third_order_upwind (-nvgrid,g(:,imu,ivmu),dvpa,mirror_sign(ig),dgdv(:,imu,ivmu))
       end do
    end do

  end subroutine get_dgdvpa

  subroutine add_mirror_term (g, src)

    use stella_layouts, only: gxyz_lo
    use stella_layouts, only: imu_idx, is_idx
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0

    implicit none

    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in out) :: src

    integer :: imu, is, ixyz

    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       imu = imu_idx(gxyz_lo,ixyz)
       is = is_idx(gxyz_lo,ixyz)
       src(:,:,:,ixyz) = src(:,:,:,ixyz) + spread(spread(mirror(:,imu,is),1,naky),2,ntheta0)*g(:,:,:,ixyz)
    end do

  end subroutine add_mirror_term

  subroutine get_dgdz (g, dgdz)

    use finite_differences, only: third_order_upwind_theta
    use stella_layouts, only: gxyz_lo
    use stella_layouts, only: iv_idx
    use theta_grid, only: ntgrid, delthet
    use kt_grids, only: naky

    implicit none

    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (out) :: dgdz

    integer :: ixyz, iseg, ie, ik, iv
    complex, dimension (2) :: gleft, gright

    ! FLAG -- assuming delta theta is equally spaced below!
    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       iv = iv_idx(gxyz_lo,ixyz)
       do ik = 1, naky
          do ie = 1, neigen(ik)
             do iseg = 1, nsegments(ie,ik)
                call fill_theta_ghost_zones (iseg, ie, ik, g(ik,:,:,ixyz), gleft, gright)
                call third_order_upwind_theta (ig_low(iseg), iseg, nsegments(ie,ik), &
                     g(ik,itmod(iseg,ie,ik),ig_low(iseg):ig_up(iseg),ixyz), &
                     delthet(0), stream_sign(iv), gleft, gright, &
                     dgdz(ik,itmod(iseg,ie,ik),ig_low(iseg):ig_up(iseg),ixyz))
             end do
          end do
       end do
    end do

  end subroutine get_dgdz

  subroutine add_stream_term (g, src)

    use stella_layouts, only: gxyz_lo
    use stella_layouts, only: iv_idx, is_idx
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0

    implicit none

    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in out) :: src

    integer :: iv, is, ixyz

    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       iv = iv_idx(gxyz_lo,ixyz)
       is = is_idx(gxyz_lo,ixyz)
       src(:,:,:,ixyz) = src(:,:,:,ixyz) + spread(spread(stream(:,iv,is),1,naky),2,ntheta0)*g(:,:,:,ixyz)
    end do

  end subroutine add_stream_term

  subroutine fill_theta_ghost_zones (iseg, ie, ik, g, gleft, gright)

    use theta_grid, only: ntgrid

    implicit none

    integer, intent (in) :: iseg, ie, ik
    complex, dimension (:,-ntgrid:), intent (in) :: g
    complex, dimension (:), intent (out) :: gleft, gright

    ! stream_sign > 0 --> stream speed < 0

    if (iseg == 1) then
       ! zero BC for g for theta < theta_min and vpa > 0
       gleft = 0.0
    else
       ! connect to segment with smaller theta-theta0 (on left)
       gleft = g(itmod(iseg-1,ie,ik),ig_up(iseg-1)-2:ig_up(iseg-1)-1)
    end if
    
    if (nsegments(ie,ik) > iseg) then
       ! connect to segment with larger theta-theta0 (on right)
       gright = g(itmod(iseg+1,ie,ik),ig_low(iseg+1)+1:ig_low(iseg+1)+2)
    else
       ! zero BC for g for theta > theta_max
       gright = 0.0
    end if

  end subroutine fill_theta_ghost_zones

  subroutine get_dgdy (g, dgdy)

    use constants, only: zi
    use stella_layouts, only: gxyz_lo
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, aky

    implicit none

    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (out) :: dgdy

    integer :: ixyz

    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       dgdy(:,:,:,ixyz) = zi*spread(spread(aky,2,ntheta0),3,2*ntgrid+1)*g(:,:,:,ixyz)
    end do

  end subroutine get_dgdy

  subroutine add_dgdy_term (g, src)

    use dist_fn_arrays, only: wdrifty
    use stella_layouts, only: gxyz_lo
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0

    implicit none

    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in out) :: src

    integer :: ixyz

    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       src(:,:,:,ixyz) = src(:,:,:,ixyz) + spread(spread(wdrifty(:,ixyz),1,naky),2,ntheta0)*g(:,:,:,ixyz)
    end do

  end subroutine add_dgdy_term

  subroutine get_dgdx (g, dgdx)

    use constants, only: zi
    use stella_layouts, only: gxyz_lo
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, akx

    implicit none

    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (out) :: dgdx

    integer :: ixyz

    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       dgdx(:,:,:,ixyz) = zi*spread(spread(akx,1,naky),3,2*ntgrid+1)*g(:,:,:,ixyz)
    end do

  end subroutine get_dgdx

  subroutine add_dgdx_term (g, src)

    use dist_fn_arrays, only: wdriftx
    use stella_layouts, only: gxyz_lo
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0

    implicit none

    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in out) :: src

    integer :: ixyz

    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       src(:,:,:,ixyz) = src(:,:,:,ixyz) + spread(spread(wdriftx(:,ixyz),1,naky),2,ntheta0)*g(:,:,:,ixyz)
    end do

  end subroutine add_dgdx_term

  subroutine get_dchidy (phi, apar, dchidy)

    use constants, only: zi
    use dist_fn_arrays, only: aj0x
    use stella_layouts, only: gxyz_lo
    use stella_layouts, only: is_idx, iv_idx
    use run_parameters, only: fphi, fapar
    use species, only: spec
    use theta_grid, only: ntgrid
    use vpamu_grids, only: vpa
    use kt_grids, only: ntheta0, aky

    implicit none

    complex, dimension (:,:,-ntgrid:), intent (in) :: phi, apar
    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (out) :: dchidy

    integer :: ixyz, iv, is

    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       is = is_idx(gxyz_lo,ixyz)
       iv = iv_idx(gxyz_lo,ixyz)
       dchidy(:,:,:,ixyz) = zi*spread(spread(aky,2,ntheta0),3,2*ntgrid+1)*aj0x(:,:,:,ixyz) &
            * ( fphi*phi - fapar*vpa(iv)*spec(is)%stm*apar )
    end do

  end subroutine get_dchidy

  subroutine add_wstar_term (g, src)

    use dist_fn_arrays, only: wstar
    use stella_layouts, only: gxyz_lo
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0

    implicit none

    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-ntgrid:,gxyz_lo%llim_proc:), intent (in out) :: src

    integer :: ixyz

    do ixyz = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
       src(:,:,:,ixyz) = src(:,:,:,ixyz) + spread(spread(wstar(:,ixyz),1,naky),2,ntheta0)*g(:,:,:,ixyz)
    end do

  end subroutine add_wstar_term

  subroutine finish_dist_fn

    implicit none

    dist_fn_initialized = .false.
    readinit = .false.
    gxyz_initialized = .false.

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
