module theta_grid_gridgen
  implicit none

  public :: theta_grid_gridgen_init
  public :: gridgen_get_grids

  private

  ! knobs
  integer :: npadd
  real :: alknob, epsknob, bpknob, extrknob, tension
  real :: thetamax, deltaw, widthw
  logical :: exist

contains

  subroutine theta_grid_gridgen_init
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.
    call read_parameters
  end subroutine theta_grid_gridgen_init

  subroutine read_parameters
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: in_file
    namelist /theta_grid_gridgen_knobs/ &
         npadd, alknob, epsknob, bpknob, extrknob, tension, thetamax, deltaw, widthw

    npadd = 2
    alknob = 0.0
    epsknob = 1e-5
    bpknob = 1.e-8
    extrknob = 0.0
    tension = 1.0
    thetamax = 0.0
    deltaw = 0.0
    widthw = 1.0
    in_file = input_unit_exist("theta_grid_gridgen_knobs", exist)
    if (exist) read (unit=input_unit("theta_grid_gridgen_knobs"), &
         nml=theta_grid_gridgen_knobs)
  end subroutine read_parameters

  subroutine gridgen_get_grids (nperiod, ntheta, ntgrid, nbset, &
       theta, bset, bmag, gradpar, gbdrift, gbdrift0, cvdrift, &
       cvdrift0, cdrift, cdrift0, gbdrift_th, cvdrift_th, &
       gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
       Rplot, Zplot, Rprime, Zprime, aplot, aprime, Bpol)
    use gridgen4mod
    use constants
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (in out) :: theta
    real, dimension (nbset), intent (in out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (in out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
         gbdrift_th, cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
         Rplot, Zplot, Rprime, Zprime, aplot, aprime, Bpol
    integer :: ntheta_old, ntgrid_old, nbset_old
    real, dimension (-ntgrid:ntgrid) :: thetasave
    real, dimension (ntheta+1) :: thetaold, thetanew
    real, dimension (ntheta+1) :: bmagold, bmagnew
    integer :: i
    logical:: debug=.false.
if (debug) write(6,*) 'gridgen_get_grids'


    ntheta_old = ntheta
    ntgrid_old = ntgrid
    nbset_old = nbset

    thetasave = theta
    thetaold = theta(-ntheta/2:ntheta/2)
    bmagold = bmag(-ntheta/2:ntheta/2)

if (debug) write(6,*) 'gridgen_get_grids: call gridgen4_2'
    call gridgen4_2 (1,ntheta_old+1,thetaold,bmagold, npadd, &
         alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
         ntheta,nbset,thetanew,bmagnew,bset)

    if (ntheta_old /= ntheta) then
       write(*,*) 'Error in theta_grid_gridgen?'
       write(*,*) 'ntheta_old = ',ntheta_old
       write(*,*) 'ntheta_new = ',ntheta
       write(*,*) 'Stopping this run would be wise.'
       write(*,*) 'Try again with ntheta = ',ntheta_old + 2
    end if

    ! interpolate to new grid
    ntgrid = ntheta/2 + (nperiod-1)*ntheta

    theta(-ntheta/2:ntheta/2-1) = thetanew(1:ntheta)
    theta(ntheta/2) = thetanew(1) + real(2)*pi
    bmag(-ntheta/2:ntheta/2-1) = bmagnew(1:ntheta)
    bmag(ntheta/2) = bmagnew(1)
    do i = 1, nperiod-1
       theta(-ntheta/2+i*ntheta:ntheta/2-1+i*ntheta) &
            = thetanew(1:ntheta) + real(2*i)*pi
       theta(ntheta/2+i*ntheta) = thetanew(1) + real(2*(i+1))*pi
       theta(-ntheta/2-i*ntheta:ntheta/2-1-i*ntheta) &
            = thetanew(1:ntheta) - real(2*i)*pi
       bmag(-ntheta/2+i*ntheta:ntheta/2-1+i*ntheta) = bmagnew(1:ntheta)
       bmag( ntheta/2+i*ntheta) = bmagnew(1)
       bmag(-ntheta/2-i*ntheta:ntheta/2-1-i*ntheta) = bmagnew(1:ntheta)
    end do

if (debug) write(6,*) 'gridgen_get_grids: call regrid'
    call regrid (ntgrid_old, thetasave, gradpar, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gbdrift, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gbdrift0, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, cvdrift, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, cvdrift0, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, cdrift, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, cdrift0, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gbdrift_th, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, cvdrift_th, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds2, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds21, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds22, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds23, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds24, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, gds24_noq, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, grho, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, Rplot, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, Zplot, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, aplot, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, Rprime, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, Zprime, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, aprime, ntgrid, theta)
    call regrid (ntgrid_old, thetasave, Bpol, ntgrid, theta)

if (debug) write(6,*) 'gridgen_get_grids: end'
  end subroutine gridgen_get_grids

  subroutine regrid (nold, x, y, nnew, xnew)
    use splines
    implicit none
    integer, intent (in) :: nold
    real, dimension (-nold:nold), intent (in) :: x
    real, dimension (-nold:nold), intent (in out) :: y
    integer, intent (in) :: nnew
    real, dimension (-nnew:nnew), intent (in) :: xnew
    type (spline) :: spl
    integer :: i

    call new_spline (2*nold+1, x(-nold:nold), y(-nold:nold), spl)

    do i = -nnew, nnew
       y(i) = splint(xnew(i), spl)
    end do

    call delete_spline (spl)
  end subroutine regrid

end module theta_grid_gridgen

module theta_grid_salpha
  implicit none

  public :: init_theta_grid_salpha
  public :: salpha_get_sizes
  public :: salpha_get_grids

  private

  ! knobs
  real :: alpmhdfac, alpha1

  ! internal variable
  integer :: model_switch
  integer, parameter :: model_salpha = 1, model_alpha1 = 2, &
       model_nocurve = 3, model_ccurv = 4, model_b2 = 5, &
       model_eps = 6, model_normal_only = 7
  

  real :: shift
  logical :: exist


contains

   subroutine init_theta_grid_salpha
     use theta_grid_params, only: init_theta_grid_params, rhoc, eps, epsl
     use geometry, only: rhoc_geo=>rhoc
     implicit none
     logical, save :: initialized = .false.

     if (initialized) return
    initialized = .false.

    call init_theta_grid_params
! make rhoc consistent with eps, epsl, and insert this value into geometry 
    if (epsl > epsilon(0.0)) then
       rhoc = 2.*eps/epsl
    else
       rhoc = 1.
    end if
    rhoc_geo=rhoc

    call read_parameters
  end subroutine init_theta_grid_salpha

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use theta_grid_params, only: shift_in => shift, alpmhd
    use text_options, only: text_option, get_option_value
    implicit none

    character(20) :: model_option
    type (text_option), dimension (8), parameter :: modelopts = &
         (/ text_option('default', model_salpha), &
            text_option('s-alpha', model_salpha), &
            text_option('alpha1', model_alpha1), &
            text_option('rogers', model_eps), &
            text_option('b2', model_b2), &
            text_option('normal_only', model_normal_only), &
            text_option('const-curv', model_ccurv), &
            text_option('no-curvature', model_nocurve) /)

    namelist /theta_grid_salpha_knobs/ alpmhdfac, alpha1, model_option
    integer :: ierr, in_file

    alpmhdfac = 0.0
    alpha1 = 0.0
    model_option = 'default'
    in_file = input_unit_exist("theta_grid_salpha_knobs", exist)
!    if (exist) read (unit=input_unit("theta_grid_salpha_knobs"), &
!         nml=theta_grid_salpha_knobs)
    if (exist) read (unit=in_file,nml=theta_grid_salpha_knobs)

    ierr = error_unit()
    call get_option_value &
         (model_option, modelopts, model_switch, &
         ierr, "model_option in theta_grid_salpha_knobs")

    if (alpmhdfac > epsilon(0.0)) then
       shift = - alpmhd*alpmhdfac
    else
       shift = shift_in
    end if

  end subroutine read_parameters

  subroutine salpha_get_sizes (nthetaout, nperiodout, nbsetout)
    use theta_grid_params, only: ntheta, nperiod
    implicit none
    integer, intent (out) :: nthetaout, nperiodout, nbsetout

    nthetaout = ntheta
    nperiodout = nperiod
    nbsetout = ntheta/2+1 ! upper bound when alpha1 model is used
  end subroutine salpha_get_sizes

  subroutine salpha_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
       gbdrift_th, cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
       Rplot, Zplot, Rprime, Zprime, aplot, aprime, shat, drhodpsi, kxfac, &
       qval, shape, Bpol)
    use constants
    use theta_grid_params, only: eps, epsl, shat_param => shat, pk, rhoc
    use theta_grid_gridgen, only: theta_grid_gridgen_init, gridgen_get_grids
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (out) :: theta
    real, dimension (nbset), intent (out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
         gbdrift_th, cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
         Rplot, Zplot, Rprime, Zprime, aplot, aprime, Bpol
    real, intent (out) :: shat, drhodpsi, kxfac, qval
    character (8), intent(out) :: shape
    integer :: i

    theta = (/ (real(i)*2.0*pi/real(ntheta), i=-ntgrid,ntgrid) /)

! BD: dummy response for graphics in s-alpha mode until I have time to fix this:
    if (abs(epsl) > epsilon(0.)) then
       Rplot = 2./epsl*(1.+eps*cos(theta))  ; Rprime = 0.
    else
       Rplot = 1. ; Rprime = 0.
    end if
    Zplot = 1.  ; Zprime = 0.
    aplot = 1.  ; aprime = 0.

! MB : should look into changing this
    Bpol = 0.

    if (model_switch == model_alpha1) then
       bmag = 1.0-eps*cos(theta)-alpha1*cos(3.0*theta)
    else if (model_switch == model_b2) then
       bmag = 1.0 - eps*cos(theta)
    else
       bmag = 1.0/(1.0 + eps*cos(theta))
    end if

    shat = shat_param
    if (eps > epsilon(0.0)) then
       drhodpsi = 0.5*epsl**2/(pk*eps)
    else
       drhodpsi = 1.0
    end if
    kxfac = 0.5
    if (epsl > epsilon(0.0)) then
       qval = epsl/pk
    else
       qval = 1.
    end if
    select case (model_switch)
    case (model_salpha,model_alpha1,model_b2)
       cvdrift = epsl*(cos(theta) + (shat*theta-shift*sin(theta))*sin(theta))
       cvdrift0 = -epsl*shat*sin(theta)
       gds2 = 1.0 + (shat*theta-shift*sin(theta))**2
       gds21 = -shat*(shat*theta - shift*sin(theta))
       gds22 = shat*shat
       grho = 1.0
       if (model_switch == model_b2) then
          cvdrift = cvdrift/bmag**2
          cvdrift0 = cvdrift0/bmag**2
       end if
       if (epsl < epsilon(0.)) shape = 'slab    '
       gbdrift = cvdrift
       gbdrift0 = cvdrift0
    
    case (model_normal_only)
       cvdrift = epsl*cos(theta)
       cvdrift0 = 0.
       gds2 = 1.0 + (shat*theta-shift*sin(theta))**2
       gds21 = -shat*(shat*theta - shift*sin(theta))
       gds22 = shat*shat
       grho = 1.0
       if (epsl < epsilon(0.)) shape = 'slab    '
       gbdrift = cvdrift
       gbdrift0 = cvdrift0
    
    case (model_eps)
       cvdrift = epsl*(cos(theta) -eps + (shat*theta-shift*sin(theta))*sin(theta))
       cvdrift0 = -epsl*shat*sin(theta)
       gds2 = 1.0 + (shat*theta-shift*sin(theta))**2
       gds21 = -shat*(shat*theta - shift*sin(theta))
       gds22 = shat*shat
       grho = 1.0
       if (epsl < epsilon(0.)) shape = 'slab    '
       gbdrift = cvdrift
       gbdrift0 = cvdrift0
    
    case (model_ccurv,model_nocurve)
       cvdrift = epsl
       cvdrift0 = 0.0

! Some strangeness here to get straight at some point:
!    ccurv == constant curvature should be the case used for cylindrical
!             geometry, but evidently Paolo and Barrett do not like the 
!             gds2 definition there, and have been using the slab
!             option (no_curvature) for their Z-pinch studies.  
!
!    Simply need to look into the shift dependence of gds2
!
       if (model_switch == model_nocurve) then
!dja fix for no curvature
          gbdrift = 0.0
!dja end
          gds2 = 1.0 + (shat*theta)**2
          gds21 = -shat*shat*theta
          shape = 'slab    '
          gbdrift = cvdrift*(1.-shift)
          gbdrift0 = cvdrift0
    
       else
          gds2 = 1.0 + (shat*theta-shift*sin(theta))**2
! probably should be:
!          gds2 = 1.0 + (shat*theta)**2
          gds21 = -shat*shat*theta
          shape = 'cylinder'
          gbdrift = cvdrift*(1.-shift)
          gbdrift0 = cvdrift0
       endif

       gds22 = shat*shat
       grho = 1.0

    end select
    gradpar = pk/2.0

    ! not sure about factor of epsl below...
    cdrift = 2.*epsl*(cos(theta)+shat*theta*sin(theta))
    cdrift0 = -2.*epsl*shat*sin(theta)
    ! MB: gds23 and gds24 are geometrical factors appearing at next order in gk eqn
    ! MB: NEED TO INCLUDE SHIFT IN BELOW EXPRESSIONS
    gds23 = -0.5*epsl*shat*theta*(1.+2.*eps*cos(theta))/eps
    gds24_noq = 0.5*epsl*(1.+eps*cos(theta))/eps
    gds24 = shat*gds24_noq

    ! MB: NEED TO INCLUDE SHIFT BELOW
    cvdrift_th = -0.25*(cos(theta))*epsl**2/eps
    gbdrift_th = cvdrift_th

    if (model_switch /= model_alpha1) then
       bset = bmag(-ntheta/2:0)
    else
       call theta_grid_gridgen_init
       call gridgen_get_grids (nperiod, ntheta, ntgrid, nbset, &
            theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
            gbdrift_th, cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
            Rplot, Zplot, Rprime, Zprime, aplot, aprime, Bpol)
    end if
  end subroutine salpha_get_grids

end module theta_grid_salpha

module theta_grid_eik
  implicit none

  public :: init_theta_grid_eik
  public :: eik_get_sizes
  public :: eik_get_grids

  private
  logical :: exist

contains

  subroutine init_theta_grid_eik

    use geometry, only: init_theta, nperiod_geo => nperiod
    use geometry, only: eikcoefs, itor, delrho, rhoc
    use geometry, only: gen_eq, ppl_eq, transp_eq
    use theta_grid_params, only: init_theta_grid_params, ntheta, nperiod

    implicit none

    real :: rhoc_save
    logical, save :: initialized = .false.
!CMR nov04: adding following debug switch
    logical :: debug=.false.
!CMR

    if (initialized) return
    initialized = .true.
if (debug) write(6,*) "init_theta_grid_eik: call init_theta_grid_params, ntheta=",ntheta
! After this call, would think you have ntheta from input file
! stored in theta_grid_params data structure.
! but when running from numerical equilibrium, this is not right
! Instead, get it stored via the eikcoefs call below.  
    call init_theta_grid_params

if (debug) write(6,*) "init_theta_grid_eik: call read_parameters, ntheta=",ntheta
    call read_parameters
!CMR replace call init_theta(ntheta) with following condition 
!    to avoid inappropriate calls to init_theta (as in geo/et.f90)
    if(.not. gen_eq .and. .not. ppl_eq .and. &
       .not. transp_eq ) then 
       if (debug) write(6,*) "init_theta_grid_eik: call init_theta, ntheta=",ntheta
       call init_theta (ntheta)
    endif
!CMRend
    nperiod_geo = nperiod 
    rhoc_save = rhoc
    if (itor == 0) rhoc = 1.5*delrho
!    print *, 'itor= ',itor, ' rhoc= ',rhoc, 'rhoc_save = ',rhoc_save
if (debug) write(6,*) "init_theta_grid_eik: call eikcoefs, ntheta=",ntheta
    call eikcoefs (ntheta)
if (debug) write(6,*) "init_theta_grid_eik: done, ntheta=",ntheta

!    write (*,*) 'init_theta_grid_eik: ntheta = ',ntheta

    rhoc = rhoc_save

  end subroutine init_theta_grid_eik

  subroutine eik_get_sizes (nthetaout, nperiodout, nbsetout)
    use geometry, only: nperiod
    use theta_grid_params, only: ntheta, nperiod_th => nperiod
    implicit none
    integer, intent (out) :: nthetaout, nperiodout, nbsetout

    nthetaout = ntheta
    nperiodout = nperiod
    nbsetout = ntheta/2+1 ! upper bound
  end subroutine eik_get_sizes

  subroutine eik_get_grids (nperiod, ntheta, ntgrid, nbset, theta, bset, bmag,&
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
            gbdrift_th, cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq, &
            grho, Rplot, Zplot, Rprime, Zprime, aplot, aprime, shat, drhodpsi,&
            kxfac, qval, gb_to_cv, Bpol)
    use theta_grid_gridgen, only: theta_grid_gridgen_init, gridgen_get_grids
    use geometry, only: kxfac_out => kxfac
    use geometry, only: theta_out => theta
    use geometry, only: gradpar_out => gradpar
    use geometry, only: bmag_out => bmag
    use geometry, only: cvdrift_out => cvdrift
    use geometry, only: cvdrift0_out => cvdrift0
    use geometry, only: gbdrift_out => gbdrift
    use geometry, only: gbdrift0_out => gbdrift0
    use geometry, only: cdrift_out => cdrift
    use geometry, only: cdrift0_out => cdrift0
    use geometry, only: gbdrift_th_out => gbdrift_th
    use geometry, only: cvdrift_th_out => cvdrift_th
    use geometry, only: gds2_out => gds2
    use geometry, only: gds21_out => gds21
    use geometry, only: gds22_out => gds22
    use geometry, only: gds23_out => gds23
    use geometry, only: gds24_out => gds24
    use geometry, only: gds24_noq_out => gds24_noq
    use geometry, only: grho_out => grho
    use geometry, only: Rplot_out => Rplot
    use geometry, only: Zplot_out => Zplot
    use geometry, only: aplot_out => aplot
    use geometry, only: Rprime_out => Rprime
    use geometry, only: Zprime_out => Zprime
    use geometry, only: aprime_out => aprime
    use geometry, only: Bpol_out => Bpol
    use geometry, only: qsf
    use geometry, only: s_hat_new, drhodpsin
    implicit none
    integer, intent (in) :: nperiod
    integer, intent (in out) :: ntheta, ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (out) :: theta
    real, dimension (nbset), intent (out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
         gbdrift_th, cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
         Rplot, Zplot, Rprime, Zprime, aplot, aprime, Bpol
    real, intent (out) :: shat, drhodpsi, kxfac, qval
    logical, intent (in) :: gb_to_cv
    integer :: ig
    logical:: debug=.false.
if (debug) write(6,*) 'eik_get_grids: ntgrid=',ntgrid
    do ig=-ntgrid,ntgrid
       theta(ig)     = theta_out(ig)
       gradpar(ig)   = gradpar_out(ig)
       bmag(ig)      = bmag_out(ig)
       cvdrift(ig)   = cvdrift_out(ig)
       cvdrift0(ig)  = cvdrift0_out(ig)
       gbdrift(ig)   = gbdrift_out(ig)
       gbdrift0(ig)  = gbdrift0_out(ig)
       cdrift(ig)    = cdrift_out(ig)
       cdrift0(ig)   = cdrift0_out(ig)
       gbdrift_th(ig)= gbdrift_th_out(ig)
       cvdrift_th(ig)= cvdrift_th_out(ig)
       gds2(ig)      = gds2_out(ig)
       gds21(ig)     = gds21_out(ig)
       gds22(ig)     = gds22_out(ig)
       gds23(ig)     = gds23_out(ig)
       gds24(ig)     = gds24_out(ig)
       gds24_noq(ig) = gds24_noq_out(ig)
       grho(ig)      = grho_out(ig)
       Rplot(ig)     = Rplot_out(ig)
       Zplot(ig)     = Zplot_out(ig)
       aplot(ig)     = aplot_out(ig)
       Rprime(ig)    = Rprime_out(ig)
       Zprime(ig)    = Zprime_out(ig)
       aprime(ig)    = aprime_out(ig)
       Bpol(ig)      = Bpol_out(ig)
    end do

    if (gb_to_cv) then
       do ig=-ntgrid,ntgrid
          gbdrift(ig) = cvdrift_out(ig)
          gbdrift0(ig) = cvdrift0_out(ig)
          gbdrift_th(ig) = cvdrift_th_out(ig)
       end do
    end if

!    do ig=-ntgrid,ntgrid
!       write (*,*) theta(ig), gradpar(ig), bmag(ig), grho(ig), &
!            gbdrift(ig), gbdrift(ig), gds2(ig)
!    end do

if (debug) write(6,*) 'eik_get_grids: call theta_grid_gridgen_init'
    call theta_grid_gridgen_init
if (debug) write(6,*) 'eik_get_grids: call gridgen_get_grids'
    call gridgen_get_grids (nperiod, ntheta, ntgrid, nbset, &
         theta, bset, bmag, &
         gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
         gbdrift_th, cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq, &
         grho, Rplot, Zplot, Rprime, Zprime, aplot, aprime, Bpol)
    shat = s_hat_new
    drhodpsi = drhodpsin
    kxfac = kxfac_out
    qval = qsf

!    write (*,*) 
!    do ig=-ntgrid,ntgrid
!       write (*,*) theta(ig), gradpar(ig), bmag(ig), grho(ig), &
!            gbdrift(ig), gbdrift(ig), gds2(ig)
!    end do

if (debug) write(6,*) 'eik_get_grids: end'
  end subroutine eik_get_grids

  subroutine read_parameters
    use file_utils, only: input_unit, input_unit_exist
    use geometry, only: nperiod
    use geometry, only: rhoc
    use geometry, only: itor, iflux, irho
    use geometry, only: ppl_eq, gen_eq, efit_eq, eqfile, local_eq, dfit_eq, gs2d_eq
    use geometry, only: equal_arc, transp_eq, idfit_eq
    use geometry, only: bishop
    use geometry, only: s_hat_input
    use geometry, only: alpha_input, invLp_input, beta_prime_input, dp_mult
    use geometry, only: rmaj, r_geo
    use geometry, only: shift, qinp, akappa, akappri, tri, tripri, asym, asympri
    use geometry, only: delrho, rmin, rmax
    use geometry, only: isym, in_nt, writelots
    use geometry, only: write_profile_variation, read_profile_variation
    use theta_grid_params, only: nperiod_in => nperiod
    use theta_grid_params, only: rhoc_in => rhoc
    use theta_grid_params, only: rmaj_in => rmaj, r_geo_in => r_geo
    use theta_grid_params, only: qinp_in => qinp, shat_in => shat
    use theta_grid_params, only: alpmhd_in => alpmhd
    use theta_grid_params, only: shift_in => shift
    use theta_grid_params, only: akappa_in => akappa, akappri_in => akappri
    use theta_grid_params, only: tri_in => tri, tripri_in => tripri
    use theta_grid_params, only: asym_in => asym, asympri_in => asympri
    use theta_grid_params, only: betaprim_in => betaprim
    implicit none
    integer :: in_file

    namelist /theta_grid_eik_knobs/ itor, iflux, irho, &
         ppl_eq, gen_eq, efit_eq, eqfile, dfit_eq, &
         equal_arc, bishop, local_eq, idfit_eq, gs2d_eq, transp_eq, &
         s_hat_input, alpha_input, invLp_input, beta_prime_input, dp_mult, &
         delrho, rmin, rmax, isym, writelots, &
         write_profile_variation, read_profile_variation

    nperiod = nperiod_in  
    rhoc = rhoc_in
    s_hat_input = shat_in
    alpha_input = alpmhd_in
    rmaj = rmaj_in
    r_geo = r_geo_in
    shift = shift_in
    qinp = qinp_in
    akappa = akappa_in
    akappri = akappri_in
    tri = tri_in
    tripri = tripri_in
    asym = asym_in
    asympri = asympri_in
    beta_prime_input = betaprim_in

    itor = 1
    iflux = 0
    irho = 2
    equal_arc = .true.
    bishop = 5
    dp_mult = 1.0
    delrho = 1e-3
    rmin = 1e-3
    rmax = 1.0
    isym = 0
    in_nt = .false.
    writelots = .false.
    local_eq = .true.
    write_profile_variation = .false.
    read_profile_variation = .false.
    
    in_file = input_unit_exist("theta_grid_eik_knobs", exist)
    if (exist) read (unit=input_unit("theta_grid_eik_knobs"), nml=theta_grid_eik_knobs)
  end subroutine read_parameters
end module theta_grid_eik

module theta_grid_file
  implicit none

  public :: init_theta_grid_file
  public :: file_get_sizes
  public :: file_get_grids

  private

  character(200) :: gridout_file
  real :: shat_input, drhodpsi_input, kxfac_input, qval_input
  logical :: no_geo_info = .false.
  integer, public :: ntheta, nperiod, ntgrid, nbset
  logical :: exist
contains

  subroutine init_theta_grid_file
    use theta_grid_params, only: init_theta_grid_params
    implicit none
    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call init_theta_grid_params
    call read_parameters
  end subroutine init_theta_grid_file

  subroutine read_parameters
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: in_file
    namelist /theta_grid_file_knobs/ gridout_file, no_geo_info

    gridout_file = "grid.out"
    in_file = input_unit_exist("theta_grid_file_knobs", exist)
    if (exist) read (unit=input_unit("theta_grid_file_knobs"), nml=theta_grid_file_knobs)
  end subroutine read_parameters

  subroutine file_get_sizes
    use file_utils, only: get_unused_unit
    implicit none
    integer :: unit
    character(200) :: line
    integer :: i, ntgrid
    real :: rmaj

    call get_unused_unit (unit)
    open (unit=unit, file=gridout_file, status="old")

    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt=*) nbset
    read (unit=unit, fmt="(a)") line
    do i = 1, nbset
       read (unit=unit, fmt="(a)") line
    end do

    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt=*) ntgrid, nperiod, ntheta, &
         drhodpsi_input, rmaj, shat_input, kxfac_input, qval_input

    close (unit=unit)
  end subroutine file_get_sizes

  subroutine file_get_grids (ntgrid, nbset, theta, bset, &
       bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
       gbdrift_th, cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
       Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
       shat, drhodpsi, kxfac, qval, gb_to_cv)
    use file_utils, only: get_unused_unit
    implicit none
    integer, intent (in out) :: ntgrid, nbset
    real, dimension (-ntgrid:ntgrid), intent (out) :: theta
    real, dimension (nbset), intent (out) :: bset
    real, dimension (-ntgrid:ntgrid), intent (out) :: &
         bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
         gbdrift_th, cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
         Rplot, Zplot, Rprime, Zprime, aplot, aprime
    real, intent (out) :: shat, drhodpsi, kxfac, qval
    logical, intent (in) :: gb_to_cv
    integer :: unit
    character(200) :: line
    integer :: i

    shat = shat_input
    drhodpsi = drhodpsi_input
    kxfac = kxfac_input
    qval = qval_input

    call get_unused_unit (unit)
    open (unit=unit, file=gridout_file, status="old")
    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt="(a)") line
    do i = 1, nbset
       read (unit=unit, fmt=*) bset(i) ! actually alambda
    end do
    bset = 1.0/bset ! switch alambda to bset

    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt="(a)") line

    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) gbdrift(i), gradpar(i), grho(i)
    end do

    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) cvdrift(i), gds2(i), bmag(i), theta(i)
    end do

    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) gds21(i), gds22(i)
    end do

    ! TMP UNTIL WORK OUT HOW TO GET FROM FILE
    gds23 = 0. ; gds24 = 0. ; gds24_noq = 0.
    gbdrift_th = 0. ; cvdrift_th = 0.

    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) cvdrift0(i), gbdrift0(i)
    end do

    if (gb_to_cv) then
       do i =-ntgrid,ntgrid
          gbdrift(i) = cvdrift(i)
          gbdrift0(i) = cvdrift0(i)
          gbdrift_th(i) = cvdrift_th(i)
       end do
    end if

    if (.not. no_geo_info) then

       read (unit=unit, fmt="(a)",err=100) line
       do i = -ntgrid, ntgrid
          read (unit=unit, fmt=*, err=100) Rplot(i), Rprime(i)
       end do

       read (unit=unit, fmt="(a)",err=100) line
       do i = -ntgrid, ntgrid
          read (unit=unit, fmt=*, err=100) Zplot(i), Zprime(i)
       end do

       read (unit=unit, fmt="(a)",err=100) line
       do i = -ntgrid, ntgrid
          read (unit=unit, fmt=*, err=100) aplot(i), aprime(i)
       end do

       close (unit=unit)    
       return
    end if

    ! TMP UNTIL FIGURE OUT HOW TO WORK WITH FILE -- MAB
    ! set coriolis drift to zero
    cdrift = 0. ; cdrift0 = 0.

100 continue

! dummy values for backward compatibility
    Rplot = 1. ; Rprime = 0.
    Zplot = 1. ; Zprime = 0.
    aplot = 1. ; aprime = 0.

    close (unit=unit)

  end subroutine file_get_grids

end module theta_grid_file

module theta_grid

  implicit none

  public :: init_theta_grid
  public :: theta, theta2
  public :: delthet, delthet2
  public :: dbdthet
  public :: gradpar
  public :: bmag
  public :: bset
  public :: dl_over_b
  public :: itor_over_B, IoB
  public :: gbdrift, gbdrift0
  public :: cvdrift, cvdrift0
  public :: cdrift, cdrift0
  public :: gbdrift_th, cvdrift_th
  public :: gds2, gds21, gds22, kxfac, qval
  public :: gds23, gds24, gds24_noq
  public :: grho
  public :: bmin, bmax, eps, shat, drhodpsi, jacob
  public :: ntheta, ntgrid, nperiod, nbset
  public :: Rplot, Zplot, aplot, Rprime, Zprime, aprime, Bpol
  public :: shape, gb_to_cv

  private

  real, dimension (:), allocatable :: theta, theta2, delthet, delthet2
  real, dimension (:), allocatable :: bset
  real, dimension (:), allocatable :: bmag, gradpar, dbdthet
  real, dimension (:), allocatable :: itor_over_B, IoB
  real, dimension (:), allocatable :: gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0
  real, dimension (:), allocatable :: gbdrift_th, cvdrift_th
  real, dimension (:), allocatable :: gds2, gds21, gds22, gds23, gds24, gds24_noq
  real, dimension (:), allocatable :: grho, jacob, dl_over_b
  real, dimension (:), allocatable :: Rplot, Zplot, aplot, Bpol
  real, dimension (:), allocatable :: Rprime, Zprime, aprime
  
  real :: bmin, bmax, eps, shat, drhodpsi, kxfac, qval
  integer :: ntheta, ntgrid, nperiod, nbset
  logical :: gb_to_cv

  ! internal variables
  integer :: eqopt_switch
  integer, parameter :: eqopt_eik = 1, eqopt_salpha = 2, eqopt_file = 3
  character (8) :: shape
  logical :: exist

contains

  subroutine init_theta_grid

    use mp, only: proc0

    implicit none

    logical, save :: initialized = .false.
    logical :: debug=.false.
    if (initialized) return
    initialized = .true.

    if (proc0) then
if (debug) write(6,*) "init_theta_grid: call read_parameters"
       call read_parameters
if (debug) write(6,*) "init_theta_grid: call get_sizes"
       call get_sizes
if (debug) write(6,*) "init_theta_grid: call allocate_arrays"
       call allocate_arrays
if (debug) write(6,*) "init_theta_grid: call get_grids"
       call get_grids
if (debug) write(6,*) "init_theta_grid: call finish_init"
       call finish_init
    end if
    call broadcast_results

    ! get useful quantity for field line averages
    if (.not.allocated(dl_over_b)) then
       allocate (dl_over_b(-ntgrid:ntgrid))
       dl_over_b = delthet*jacob
       dl_over_b = dl_over_b / sum(dl_over_b)
    end if

  end subroutine init_theta_grid

  subroutine broadcast_results

    use mp, only: proc0, broadcast
    use geometry, only: rhoc, d2psidr2

    implicit none

    integer :: i

    call broadcast (bmin)
    call broadcast (bmax)
    call broadcast (eps)
    call broadcast (kxfac)
    call broadcast (rhoc)
    call broadcast (qval)
    call broadcast (ntheta)
    call broadcast (ntgrid)
    call broadcast (nperiod)
    call broadcast (nbset)
    call broadcast (d2psidr2)

    if (.not. proc0) then
       call allocate_arrays
       allocate (dbdthet(-ntgrid:ntgrid))
       allocate (theta2(-ntgrid:ntgrid))
       allocate (delthet(-ntgrid:ntgrid))
       allocate (delthet2(-ntgrid:ntgrid))
    end if
    call broadcast (theta)
    call broadcast (theta2)
    call broadcast (delthet)
    call broadcast (delthet2)
    call broadcast (bset)
    call broadcast (bmag)
    call broadcast (itor_over_B)
    call broadcast (IoB)
    call broadcast (gradpar)
    call broadcast (gbdrift)
    call broadcast (gbdrift0)
    call broadcast (cvdrift)
    call broadcast (cvdrift0)
    call broadcast (cdrift)
    call broadcast (cdrift0)
    call broadcast (gbdrift_th)
    call broadcast (cvdrift_th)
    call broadcast (gds2)
    call broadcast (gds21)
    call broadcast (gds22)
    call broadcast (gds23)
    call broadcast (gds24)
    call broadcast (gds24_noq)
    call broadcast (grho)
    call broadcast (shat)
    call broadcast (jacob)
    call broadcast (Rplot)
    call broadcast (Zplot)
    call broadcast (aplot)
    call broadcast (Rprime)
    call broadcast (Zprime)
    call broadcast (aprime)
    call broadcast (Bpol)
    call broadcast (drhodpsi)
    call broadcast (gb_to_cv)
    call broadcast (dbdthet)
  end subroutine broadcast_results

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options
    implicit none
    type (text_option), dimension (5), parameter :: eqopts = &
         (/ text_option('default', eqopt_eik), &
            text_option('eik', eqopt_eik), &
            text_option('s-alpha', eqopt_salpha), &
            text_option('grid.out', eqopt_file), &
            text_option('file', eqopt_file) /)
    character(20) :: equilibrium_option
    ! 'default' 'eik': call eikcoefs for parameterized equilibrium
    ! 's-alpha': s-alpha
    ! 'grid.out' 'file': read grid from grid.out file generated by rungridgen
    namelist /theta_grid_knobs/ equilibrium_option, gb_to_cv
    integer :: ierr, in_file

    gb_to_cv = .false.
    equilibrium_option = 'default'
    in_file = input_unit_exist("theta_grid_knobs", exist)
!    if (exist) read (unit=input_unit("theta_grid_knobs"), nml=theta_grid_knobs)
    if (exist) read (unit=in_file, nml=theta_grid_knobs)

    ierr = error_unit()
    call get_option_value &
         (equilibrium_option, eqopts, eqopt_switch, &
         ierr, "equilibrium_option in theta_grid_knobs")
  end subroutine read_parameters

  subroutine allocate_arrays
    implicit none
    allocate (theta(-ntgrid:ntgrid))
    allocate (bset(nbset))
    allocate (bmag(-ntgrid:ntgrid))
    allocate (gradpar(-ntgrid:ntgrid))
    allocate (itor_over_B(-ntgrid:ntgrid))
    allocate (IoB(-ntgrid:ntgrid))
    allocate (gbdrift(-ntgrid:ntgrid))
    allocate (gbdrift0(-ntgrid:ntgrid))
    allocate (cvdrift(-ntgrid:ntgrid))
    allocate (cvdrift0(-ntgrid:ntgrid))
    allocate (cdrift(-ntgrid:ntgrid))
    allocate (cdrift0(-ntgrid:ntgrid))
    allocate (gbdrift_th(-ntgrid:ntgrid))
    allocate (cvdrift_th(-ntgrid:ntgrid))
    allocate (gds2(-ntgrid:ntgrid))
    allocate (gds21(-ntgrid:ntgrid))
    allocate (gds22(-ntgrid:ntgrid))
    allocate (gds23(-ntgrid:ntgrid))
    allocate (gds24(-ntgrid:ntgrid))
    allocate (gds24_noq(-ntgrid:ntgrid))
    allocate (grho(-ntgrid:ntgrid))
    allocate (jacob(-ntgrid:ntgrid))
    allocate (Rplot(-ntgrid:ntgrid))
    allocate (Rprime(-ntgrid:ntgrid))
    allocate (Zplot(-ntgrid:ntgrid))
    allocate (Zprime(-ntgrid:ntgrid))
    allocate (aplot(-ntgrid:ntgrid))
    allocate (aprime(-ntgrid:ntgrid))
    allocate (Bpol(-ntgrid:ntgrid))
  end subroutine allocate_arrays

  subroutine finish_init

    implicit none

    real, dimension (nbset) :: bset_save
    real, dimension (-ntgrid:ntgrid) :: eik_save
    
    ! in case nbset changes after gridgen
    if (nbset /= size(bset)) then
       bset_save = bset(:nbset)
       deallocate (bset)
       allocate (bset(nbset))
       bset = bset_save
    end if

    ! in case ntgrid changes after gridgen
    if (ntgrid*2+1 /= size(theta)) then

       eik_save = theta(-ntgrid:ntgrid); deallocate (theta)
       allocate (theta(-ntgrid:ntgrid)); theta = eik_save

       eik_save = bmag(-ntgrid:ntgrid); deallocate (bmag)
       allocate (bmag(-ntgrid:ntgrid)); bmag = eik_save

       eik_save = gradpar(-ntgrid:ntgrid); deallocate (gradpar)
       allocate (gradpar(-ntgrid:ntgrid)); gradpar = eik_save

       eik_save = itor_over_B(-ntgrid:ntgrid); deallocate (itor_over_B)
       allocate (itor_over_B(-ntgrid:ntgrid)); itor_over_B = eik_save

       eik_save = IoB(-ntgrid:ntgrid); deallocate (IoB)
       allocate (IoB(-ntgrid:ntgrid)); IoB = eik_save

       eik_save = gbdrift(-ntgrid:ntgrid); deallocate (gbdrift)
       allocate (gbdrift(-ntgrid:ntgrid)); gbdrift = eik_save

       eik_save = gbdrift0(-ntgrid:ntgrid); deallocate (gbdrift0)
       allocate (gbdrift0(-ntgrid:ntgrid)); gbdrift0 = eik_save

       eik_save = cvdrift(-ntgrid:ntgrid); deallocate (cvdrift)
       allocate (cvdrift(-ntgrid:ntgrid)); cvdrift = eik_save

       eik_save = cvdrift0(-ntgrid:ntgrid); deallocate (cvdrift0)
       allocate (cvdrift0(-ntgrid:ntgrid)); cvdrift0 = eik_save

       eik_save = cdrift(-ntgrid:ntgrid); deallocate (cdrift)
       allocate (cdrift(-ntgrid:ntgrid)); cdrift = eik_save

       eik_save = cdrift0(-ntgrid:ntgrid); deallocate (cdrift0)
       allocate (cdrift0(-ntgrid:ntgrid)); cdrift0 = eik_save

       eik_save = gbdrift_th(-ntgrid:ntgrid); deallocate (gbdrift_th)
       allocate (gbdrift_th(-ntgrid:ntgrid)); gbdrift_th = eik_save

       eik_save = cvdrift_th(-ntgrid:ntgrid); deallocate (cvdrift_th)
       allocate (cvdrift_th(-ntgrid:ntgrid)); cvdrift_th = eik_save

       eik_save = gds2(-ntgrid:ntgrid); deallocate (gds2)
       allocate (gds2(-ntgrid:ntgrid)); gds2 = eik_save

       eik_save = gds21(-ntgrid:ntgrid); deallocate (gds21)
       allocate (gds21(-ntgrid:ntgrid)); gds21 = eik_save

       eik_save = gds22(-ntgrid:ntgrid); deallocate (gds22)
       allocate (gds22(-ntgrid:ntgrid)); gds22 = eik_save

       eik_save = gds23(-ntgrid:ntgrid); deallocate (gds23)
       allocate (gds23(-ntgrid:ntgrid)); gds23 = eik_save

       eik_save = gds24(-ntgrid:ntgrid); deallocate (gds24)
       allocate (gds24(-ntgrid:ntgrid)); gds24 = eik_save

       eik_save = gds24_noq(-ntgrid:ntgrid); deallocate (gds24_noq)
       allocate (gds24_noq(-ntgrid:ntgrid)); gds24_noq = eik_save

       eik_save = grho(-ntgrid:ntgrid); deallocate (grho)
       allocate (grho(-ntgrid:ntgrid)); grho = eik_save

       eik_save = Rplot(-ntgrid:ntgrid); deallocate (Rplot)
       allocate (Rplot(-ntgrid:ntgrid)); Rplot = eik_save

       eik_save = Zplot(-ntgrid:ntgrid); deallocate (Zplot)
       allocate (Zplot(-ntgrid:ntgrid)); Zplot = eik_save

       eik_save = aplot(-ntgrid:ntgrid); deallocate (aplot)
       allocate (aplot(-ntgrid:ntgrid)); aplot = eik_save

       eik_save = Bpol(-ntgrid:ntgrid); deallocate (Bpol)
       allocate (Bpol(-ntgrid:ntgrid)); Bpol = eik_save
    end if

    bmax = maxval(bmag)
    bmin = minval(bmag)
! ?? check Krook collision operator coding which is only place eps is used
! the line with bmin/bmax was the original coding.  Changed in 2002-2004 time 
! frame, now changed back (8.19.04) BD
    eps = 1.0 - sqrt(bmin/bmax)
!    eps = 1.0 - 1.0/bmax

    allocate (dbdthet(-ntgrid:ntgrid))
    allocate (theta2(-ntgrid:ntgrid))
    allocate (delthet(-ntgrid:ntgrid))
    allocate (delthet2(-ntgrid:ntgrid))
    
    theta2 = theta*theta
    delthet(:ntgrid-1) = theta(-ntgrid+1:) - theta(:ntgrid-1)
    delthet(ntgrid) = 0.!delthet(-ntgrid)
    delthet2 = delthet*delthet

    ! first get db/dtheta at grid points
    dbdthet(-ntgrid+1:ntgrid-1) = (bmag(-ntgrid+2:ntgrid)-bmag(-ntgrid:ntgrid-2)) &
         / (delthet(-ntgrid:ntgrid-2)+delthet(-ntgrid+1:ntgrid-1))
    dbdthet(-ntgrid) = (bmag(-ntgrid+1)-bmag(-ntgrid))/delthet(-ntgrid)
    dbdthet(ntgrid) = (bmag(ntgrid)-bmag(ntgrid-1))/delthet(ntgrid-1)

    jacob = 1.0/(drhodpsi*gradpar*bmag)
    
  end subroutine finish_init

  subroutine get_sizes
    use theta_grid_eik, only: eik_get_sizes, init_theta_grid_eik
    use theta_grid_salpha, only: salpha_get_sizes, init_theta_grid_salpha
    use theta_grid_file, only: file_get_sizes, init_theta_grid_file
    use theta_grid_file, only: ntheta_file=>ntheta, nperiod_file=>nperiod
    use theta_grid_file, only: nbset_file=>nbset
    implicit none
    logical:: debug=.false.
if (debug) write(6,*) 'get_sizes: eqopt_switch=',eqopt_switch
    select case (eqopt_switch)
    case (eqopt_eik)
if (debug) write(6,*) 'get_sizes: call init_theta_grid_eik'
       call init_theta_grid_eik
if (debug) write(6,*) 'get_sizes: call eik_get_sizes'
       call eik_get_sizes (ntheta, nperiod, nbset)
    case (eqopt_salpha)
if (debug) write(6,*) 'get_sizes: call init_theta_grid_salpha'
       call init_theta_grid_salpha
if (debug) write(6,*) 'get_sizes: call salpha_get_sizes'
       call salpha_get_sizes (ntheta, nperiod, nbset)
    case (eqopt_file)
if (debug) write(6,*) 'get_sizes: call init_theta_grid_file'
       call init_theta_grid_file
if (debug) write(6,*) 'get_sizes: call file_get_sizes'
       call file_get_sizes
       ntheta=ntheta_file
       nperiod=nperiod_file
       nbset=nbset_file
    end select
    ntgrid = ntheta/2 + (nperiod-1)*ntheta 
if (debug) write(6,*) 'get_sizes: done'
  end subroutine get_sizes

  subroutine get_grids

    use theta_grid_eik, only: eik_get_grids
    use theta_grid_salpha, only: salpha_get_grids
    use theta_grid_file, only: file_get_grids
    use theta_grid_params, only: eps, btor_slab
    use geometry, only: rhoc

    implicit none

    logical:: debug=.false.
    select case (eqopt_switch)
    case (eqopt_eik)
if (debug) write(6,*) 'get_grids: call eik_get_grids'
       call eik_get_grids (nperiod, ntheta, ntgrid, nbset, &
            theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
            gbdrift_th, cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
            Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
            shat, drhodpsi, kxfac, qval, gb_to_cv, Bpol)
       shape = 'torus   '
    case (eqopt_salpha)
if (debug) write(6,*) 'get_grids: call salpha_get_grids'
       call salpha_get_grids (nperiod, ntheta, ntgrid, nbset, &
            theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
            gbdrift_th, cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
            Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
            shat, drhodpsi, kxfac, qval, shape, Bpol)
    case (eqopt_file)
if (debug) write(6,*) 'get_grids: call file_get_grids'
       call file_get_grids (ntgrid, nbset, theta, bset, bmag, &
            gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0, &
            gbdrift_th, cvdrift_th, gds2, gds21, gds22, gds23, gds24, gds24_noq, grho, &
            Rplot, Zplot, Rprime, Zprime, aplot, aprime, &
            shat, drhodpsi, kxfac, qval, gb_to_cv)
       shape = 'torus   '
    end select
    kxfac = abs(kxfac)
    qval = abs(qval)

    itor_over_B=0.
!CMR, 2/2/2011: 
! If using slab geometry, set itor_over_B = btor_slab from "theta_grid_params": 
! cleaner equivalent alternative to using btor_slab in "dist_fn_knobs", and 
! sets geometric paramater itor_over_B in one place for ALL geometries.
!
    if (eqopt_switch .eq. eqopt_salpha .and. eps < epsilon(0.) ) then
       itor_over_B = btor_slab
    else
!CMR, 19/10/10: moved MAB's definition of geometry quantity itor_over_B from 
!               dist_fn.f90 to here.
! Calculate the parallel velocity shear drive factor itor_over_B 
! (which effectively depends on the angle the field lines make with the flow)
! note that the following is only valid in a torus!
! itor_over_B = (q/rho) * Rmaj*Btor/(a*B)
       IoB = sqrt(Rplot**2 - (grho/(bmag*drhodpsi))**2)
       ! working in terms of omprim instead of g_exb
       itor_over_B = IoB
    ! RN> 2011/1/25: fix here avoids dividing by rhoc if rhoc=0
    ! CMR, 2/2/2011: itor_over_B=0 if rhoc=0
    !                this dropping parallel sheared flow source term in GKE
    !                itor_over_B=0 is safer than itor_over_B=NaN!
!       if (abs(rhoc) > epsilon(0.)) itor_over_B = qval / rhoc * IoB
     endif
  end subroutine get_grids

  subroutine finish_theta_grid

    implicit none

    if (allocated(dl_over_b)) deallocate (dl_over_b)

  end subroutine finish_theta_grid

end module theta_grid


