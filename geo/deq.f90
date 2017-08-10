module deq

  implicit none
  private

  integer :: nw, nh, nbbbs, nrho

  real, allocatable, dimension (:) :: psi_bar, fp, qsf, pressure, beta, spsi_bar
  real, allocatable, dimension (:) :: dummy, dfit_R, dfit_Z, sdfit_R, sdfit_Z
  real, allocatable, dimension (:) :: rho_mid, psi_mid
  real, allocatable, dimension (:,:) :: dfit_psi, sdfit_psi
  real, allocatable, dimension (:,:,:) :: dpm, dtm
  real, allocatable, dimension (:) :: rbbbs, zbbbs, thetab, r_bound !boundary of plasma

  real :: dfit_dR, dfit_dZ, psi_N

  real :: R_mag, Z_mag, B_T0, aminor, beta_0

  logical :: init_bound = .true.
  logical :: init_btori = .true.
  logical :: init_dbtori = .true.
  logical :: init_q = .true.
  logical :: init_pressure = .true.
  logical :: init_dpressure = .true.
  logical :: init_beta = .true.
  logical :: init_rho = .true.
  logical :: init_psi = .true.
  logical :: init_R = .true.
  logical :: init_Z = .true.

  public :: dfit_init, dfitin, gradient, eqitem, bgradient

  public :: invR
  public :: Rpos
  public :: Zpos
  public :: btori,    initialize_btori
  public :: dbtori,   initialize_dbtori
  public :: qfun,     initialize_q
  public :: pfun,     initialize_pressure
  public :: dpfun,    initialize_dpressure
  public :: betafun,  initialize_beta
  public :: bound,    initialize_bound
  public :: psi,      initialize_psi
  public :: rhofun,   initialize_rho


contains
  
  subroutine dfitin(eqfile, psi_0, psi_a, rmaj, B_T, amin, initeq, big)
!
!     This subroutine reads an DFIT output file containing 
!     the axisymmetric magnetic field geometry on a rectangular 
!     domain defined by the coordinates (R,Z).
!
!     dfit_R     R grid
!     dfit_Z     Z grid
!     fp    F on psibar grid
!     dfit_psi   psi on (R,Z) grid
!     R_mag        major radius of the magnetic axis
!     Z_mag        elevation of the magnetic axis
!     rwid       total width of the domain
!     rleft      position of leftmost point of domain
!     zhei       total height of domain
!
    use splines, only: inter_cspl
    implicit none

    real :: xdum, p_0
    real :: rwid, rleft, zhei, amin, B_T
    real :: psi_0, psi_a, rmaj, bcentr
    real, dimension(:), allocatable :: tmp1, tmp2, zp, temp, &
             zx1, zxm, zy1, zyn
    real :: zxy11, zxym1, zxy1n, zxymn
    real :: fitp_surf2, delta_R, delta_Z
    
    character*80 :: filename, eqfile
    character char*10
    
    integer :: i, j, init, ndum, initeq, i1, big, nhb, nwb, ierr
    integer :: jmin, jmax
    
    data init /1/
    save init
    
! Need to generalize initialization condition if equilibrium changes

    if(initeq == 0) return
    init=0
    
    i=index(eqfile,' ')-1
    filename = eqfile(1:i)
    open(unit=5,file=filename,status='old',form='formatted')
    
! Read the data

    read(5,*) nw, nh, nrho

!    nwb = nw * big
!    nhb = nh * big

    nwb = nw 
    nhb = nh 

    call alloc_module_arrays(nwb, nwb, nhb, nw, nh, nrho)

    read(5,2020) rwid, zhei
    read(5,2020) R_mag, Z_mag
    read(5,2020) psi_0, psi_a, bcentr
!
! pbar is defined by
! pbar = (psi-psi_0)/(psi_a-psi_0)
! fp and q are functions of pbar
!
    do i=1,nw
       spsi_bar(i) = float(i) / float(nw)
       sdfit_R(i) = rwid * float(i) / float(nw)
    enddo

    do j=1,nh
       sdfit_Z(j) = zhei*(float(j-1)/float(nh-1)-0.5)
    enddo

    do i=1,nwb
       psi_bar(i) = float(i) / float(nwb)
       dfit_R(i) = rwid * float(i) / float(nwb)
    enddo

    do j=1,nhb
       dfit_Z(j) = zhei*(float(j-1)/float(nhb-1)-0.5)
    enddo

    read(5,2020) (dummy(j), j = 1, nw)  ! pressure read
    call inter_cspl(nw, spsi_bar, dummy, nwb, psi_bar, pressure)
    read(5,2020) ((dfit_psi(i,j) , i = 1, nw) , j = 1, nh)
!    read(5,2020) ((sdfit_psi(i,j) , i = 1, nw) , j = 1, nh)

!    allocate(zp(3*nw*nh), temp(nw+2*nh))
!    allocate(zx1(nh), zxm(nh), zy1(nw), zyn(nw))

!    call fitp_surf1(nw, nh, sdfit_R, sdfit_Z, sdfit_psi, &
!         nw, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
!         255, zp, temp, 1., ierr)

    do j = 1, nhb
       do i = 1, nwb
!          dfit_psi(i,j) = fitp_surf2(dfit_R(i), dfit_Z(j), nw, nh, &
!               sdfit_R, sdfit_Z, sdfit_psi, nw, zp, 1.)
       enddo
    enddo
      
!    deallocate(zp, temp)

    read (5,2020) (rho_mid(i), i=1,nrho)
    read (5,2020) (psi_mid(i), i=1,nrho)

    fp = 0.
    qsf = 0.

    nw = nwb
    nh = nhb

    nbbbs = 2*(nw+nh)-5

! stay away from edge of box
    delta_R = 1.02 ; delta_Z = 1.02
    
    allocate(rbbbs(nbbbs), zbbbs(nbbbs), thetab(nbbbs), r_bound(nbbbs))
    jmin = 1
    jmax = nh/2
    do j=jmin, jmax
       rbbbs(j) = (dfit_R(1)-0.5*rwid)/delta_R+0.5*rwid
       zbbbs(j) = dfit_Z(j+nh/2)/delta_Z
    end do

    jmin = jmax + 1
    jmax = jmax + nw - 2
    do j=jmin, jmax
       rbbbs(j) = (dfit_R(j-jmin+2)-0.5*rwid)/delta_R+0.5*rwid
       zbbbs(j) = dfit_Z(nh)/delta_Z
    end do

    jmin = jmax + 1
    jmax = jmax + nh - 1
    do j=jmin, jmax
       rbbbs(j) = (dfit_R(nw)-0.5*rwid)/delta_R+0.5*rwid
       zbbbs(j) = dfit_Z(jmax-j+2)/delta_Z
    end do

    jmin = jmax + 1
    jmax = jmax + nw - 1 
    do j=jmin, jmax
       rbbbs(j) = (dfit_R(jmax-j+2)-0.5*rwid)/delta_R+0.5*rwid
       zbbbs(j) = dfit_Z(1)/delta_Z
    end do

    jmin = jmax + 1
    jmax = jmax + nh/2    
    do j=jmin, jmax
       rbbbs(j) = (dfit_R(1)-0.5*rwid)/delta_R+0.5*rwid
       zbbbs(j) = dfit_Z(j-jmin+2)/delta_Z
    end do

    do j=1,nbbbs-1
       if ((rbbbs(j) == rbbbs(j+1)) .and. (zbbbs(j) == zbbbs(j+1))) then
          write(*,*) 'duplicates: ',rbbbs(j),zbbbs(j)
       end if
    end do

! get r_boundary(theta)

    thetab = atan2 ((zbbbs-Z_mag), (rbbbs-R_mag))
    r_bound = sqrt( (rbbbs - R_mag)**2 + (zbbbs - Z_mag)**2 )

    call sort(thetab, r_bound, zbbbs, rbbbs)

! Allow for duplicated points near +- pi:

    if(thetab(1) == thetab(2)) then
       thetab(1) = thetab(1) + 4.*acos(0.)
       call sort(thetab, r_bound, zbbbs, rbbbs)
    endif

    if(thetab(nbbbs-1) == thetab(nbbbs)) then
       thetab(nbbbs) = thetab(nbbbs) - 4.*acos(0.)
       call sort(thetab, r_bound, zbbbs, rbbbs)
    endif

! It isn't likely that a duplicate point would exist near theta = 0, 
! so I am not allowing this possibility for now.

    do i=1,nbbbs-1
       if(thetab(i) == thetab(i+1)) then
          write(*,*) 'Duplicates near theta = 0 not allowed.'
          write(*,*) i, i+1, ' Stopping.'
          stop
       endif
    enddo
    deallocate (rbbbs, zbbbs)

    aminor=R_mag
    amin=aminor

    close(5)


    r_bound = r_bound / aminor

    R_mag = R_mag / aminor
    Z_mag = Z_mag / aminor
!    rleft = rleft / aminor
!    rwid = rwid / aminor
!    zhei = zhei / aminor
    dfit_R = dfit_R / aminor
    dfit_Z = dfit_Z / aminor
    
! should rmaj be R_mag? use R_mag for now.  

    rmaj = R_mag
    B_T0 = abs(bcentr)
    B_T = B_T0
    psi_a = psi_a / (B_T0*aminor**2)
    psi_0 = psi_0 / (B_T0*aminor**2)
    psi_N = psi_a - psi_0
    
    p_0=pressure(1)

    fp = fp / (B_T0*aminor)

! MKS: beta = 2 mu_0 p / B**2

    beta = 8. * (2. * acos(0.)) * pressure * 1.e-7 / B_T0**2
    beta_0 = beta(1)

    pressure = pressure / p_0
    
    dfit_psi = dfit_psi / (B_T0*aminor**2)
    
    dfit_dR = dfit_R(2) - dfit_R(1)
    dfit_dZ = dfit_Z(2) - dfit_Z(1)

    2020 format (5e16.8)

  end subroutine dfitin

  subroutine dfit_init

    real, dimension(nw, nh) :: eqth 
    integer :: i, j
    
    do i = 1, nw
       do j = 1,nh
          if(dfit_Z(j) == Z_mag .and. dfit_R(i) == R_mag) then
             eqth(i,j) = 0.  ! value should not matter
          else
             eqth(i,j) = atan2( (dfit_Z(j)-Z_mag), (dfit_R(i)-R_mag))
          endif
       enddo
    enddo

    call derm(dfit_psi, dpm)
    call tderm(eqth, dtm)
    
  end subroutine dfit_init

  subroutine tderm(f, dfm)

    implicit none
    integer i, j
    real f(:,:), dfm(:,:,:), pi

    pi = 2.*acos(0.)
    
! DFIT grid is equally spaced in R, Z -- this routine uses that fact and 
! is therefore not completely general.  It is fine for DFIT output.    

    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))/dfit_dR
    
    i=nw
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))/dfit_dR
   
    j=1
    dfm(:,j,2) = -0.5*(3*f(:,j)-4*f(:,j+1)+f(:,j+2))/dfit_dZ
    
    j=nh      
    dfm(:,j,2) = 0.5*(3*f(:,j)-4*f(:,j-1)+f(:,j-2))/dfit_dZ
    
    do i=2,nw-1
       dfm(i,:,1)=0.5*(f(i+1,:)-f(i-1,:))/dfit_dR
    enddo
    
    do j=2,nh-1
       do i = 1,nw
          if(f(i,j+1)-f(i,j-1) > pi) then
             dfm(i,j,2)=0.5*(f(i,j+1)-f(i,j-1)-2.*pi)/dfit_dZ
          else
             dfm(i,j,2)=0.5*(f(i,j+1)-f(i,j-1))/dfit_dZ
          endif
       enddo
    enddo
    
  end subroutine tderm

  subroutine derm(f, dfm)

    implicit none
    integer i, j
    real f(:,:), dfm(:,:,:)
    
! DFIT grid is equally spaced in R, Z -- this routine uses that fact and 
! is therefore not completely general.  It is fine for DFIT output.    

    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))/dfit_dR
    
    i=nw
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))/dfit_dR
   
    j=1
    dfm(:,j,2) = -0.5*(3*f(:,j)-4*f(:,j+1)+f(:,j+2))/dfit_dZ
    
    j=nh      
    dfm(:,j,2) = 0.5*(3*f(:,j)-4*f(:,j-1)+f(:,j-2))/dfit_dZ
    
    do i=2,nw-1
       dfm(i,:,1)=0.5*(f(i+1,:)-f(i-1,:))/dfit_dR
    enddo
    
    do j=2,nh-1
       dfm(:,j,2)=0.5*(f(:,j+1)-f(:,j-1))/dfit_dZ
    enddo
    
  end subroutine derm

  subroutine gradient(rgrid, theta, grad, char, rp, nth, ntm)

    use splines, only: inter_d_cspl
    integer nth, ntm
    character*1 char
    real, dimension(-ntm:), intent(in) :: rgrid, theta
    real, dimension(-ntm:,:), intent(out) :: grad
    real aa(1), daa(1), rp, rpt(1)
    integer i, j
    
    grad = 0.
    do i=-nth, nth
       call eqitem(rgrid(i), theta(i), dpm(:,:,1), grad(i,1))
       call eqitem(rgrid(i), theta(i), dpm(:,:,2), grad(i,2))
    enddo

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nw, psi_bar, pressure, 1, rpt, aa, daa)
       grad = grad*daa(1)*0.5* beta_0/psi_N
    endif

  end subroutine gradient

  subroutine bgradient(rgrid, theta, grad, char, rp, nth_used, ntm)

    use splines, only: inter_d_cspl
    implicit none
    
    integer nth_used, ntm
    character*1 char
    real rgrid(-ntm:), theta(-ntm:), grad(-ntm:,:)
    real tmp(2), aa(1), daa(1), rp, rpt(1)
    real, dimension(nw, nh, 2) ::  dbish
    integer i
    logical :: first = .true.
 
    dbish(:,:,1) = sqrt(dpm(:,:,1)**2 + dpm(:,:,2)**2)
    dbish(:,:,2) = 0.
    
    if(char == 'T') then  ! the order of the next two statements matters
       dbish(:,:,2) = (dtm(:,:,2)*dpm(:,:,1)-dtm(:,:,1)*dpm(:,:,2))/dbish(:,:,1)
       dbish(:,:,1) = (dtm(:,:,1)*dpm(:,:,1)+dtm(:,:,2)*dpm(:,:,2))/dbish(:,:,1)
    endif

    do i=-nth_used,-1
       call eqitem(rgrid(i), theta(i), dbish(:,:,1), tmp(1))
       call eqitem(rgrid(i), theta(i), dbish(:,:,2), tmp(2))
       grad(i,1) = tmp(1)
       grad(i,2) = tmp(2)
    enddo

    do i=0,nth_used
       call eqitem(rgrid(i), theta(i), dbish(:,:,1), tmp(1))
       call eqitem(rgrid(i), theta(i), dbish(:,:,2), tmp(2))
       grad(i,1)=tmp(1)
       grad(i,2)=tmp(2)
    enddo

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nw, psi_bar, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1) * 0.5*beta_0/psi_N
          grad(i,2)=grad(i,2)*daa(1) * 0.5*beta_0/psi_N
       enddo
    endif

  end subroutine bgradient

  subroutine eqitem(r, thetin, f, fstar)
      
    integer :: i, j, istar, jstar
    real, intent (in) :: r, thetin, f(:,:)
    real, intent (out) :: fstar
    real st, dt, sr, dr
    real r_pos, z_pos
    
    r_pos = Rpos(r, thetin)
    z_pos = Zpos(r, thetin)

! find point on R mesh

    if(r_pos >= dfit_R(nw) .or. r_pos <= dfit_R(1)) then
       write(*,*) 'No evaluation of eqitem allowed outside'
       write(*,*) 'or on edge of R domain'
       write(*,*) r, thetin, dfit_R(nw), r_pos
       stop      
    endif

! ensure point is on Z mesh

    if(z_pos >= dfit_Z(nh) .or. z_pos <= dfit_Z(1)) then
       write(*,*) 'No evaluation of eqitem allowed outside'
       write(*,*) 'or on edge of Z domain'
       write(*,*) r, thetin, dfit_Z(1), dfit_Z(nh), z_pos
       stop
    endif
    
    istar=0
    do i=2,nw
       if(istar /= 0) cycle
       if(r_pos < dfit_R(i)) then
          dr = r_pos - dfit_R(i-1)
          sr = dfit_R(i) - r_pos
          istar=i-1
       endif
    enddo
      
! Now do Z direction

    jstar=0
    do j=1,nh
       if(jstar /= 0) cycle
       if(z_pos < dfit_Z(j)) then
          dt = z_pos - dfit_Z(j-1)
          st = dfit_Z(j) - z_pos
          jstar=j-1
       endif
    enddo
      
! use opposite area stencil to interpolate

    fstar=f(istar    , jstar    ) * sr * st &
         +f(istar + 1, jstar    ) * dr * st &
         +f(istar    , jstar + 1) * sr * dt &
         +f(istar + 1, jstar + 1) * dr * dt
    fstar = fstar &
         /abs(dfit_R(istar+1)-dfit_R(istar)) &
         /(dfit_Z(jstar+1)-dfit_Z(jstar))
!     write(*,*) i, dr, dt, sr, st
!     write(*,*) f(istar,jstar+1),f(istar+1,jstar+1)
!     write(*,*) f(istar,jstar),f(istar+1,jstar)
!     write(*,*) dfit_R(istar),dfit_R(istar+1)
!     write(*,*) dfit_Z(jstar),dfit_Z(jstar+1)

  end subroutine eqitem

  function Zpos (r, theta)
   
    real, intent (in) :: r, theta
    real :: Zpos
    
    Zpos = Z_mag + r * sin(theta)

  end function Zpos

  function Rpos (r, theta)
   
    real, intent (in) :: r, theta
    real :: Rpos

    Rpos = R_mag + r * cos(theta)
    
  end function Rpos

  function invR (r, theta)
   
    real, intent (in) :: r, theta
    real :: invR

    invR = 1/(R_mag + r*cos(theta))
    
  end function invR

  function initialize_psi (init) 

    integer :: init, initialize_psi
    
    init_psi = .false.
    if(init == 1) init_psi = .true.
    initialize_psi = 1

  end function initialize_psi

  function psi (r, theta)
   
    real, intent (in) :: r, theta
    real :: f, psi

    call eqitem(r, theta, dfit_psi, f)
    psi = f
    
  end function psi
   
  function initialize_btori (init) 

    integer :: init, initialize_btori
    
    init_btori = .false.
    if(init == 1) init_btori = .true.
    initialize_btori = 1

  end function initialize_btori

  function btori (pbar)
  
    use splines
    real :: pbar, btori
    type (spline), save :: spl

    if(init_btori) call new_spline(nw, psi_bar, fp, spl)

    btori = splint(pbar, spl)

  end function btori

  function initialize_dbtori (init) 

    integer :: init, initialize_dbtori
    
    init_dbtori = .false.
    if(init == 1) init_dbtori = .true.
    initialize_dbtori = 1

  end function initialize_dbtori

  function dbtori (pbar)
  
    use splines
    real :: pbar, dbtori
    type (spline), save :: spl

    if(init_dbtori) call new_spline(nw, psi_bar, fp, spl)

    dbtori = dsplint(pbar, spl)/psi_N

  end function dbtori

  function initialize_q (init) 

    integer :: init, initialize_q
    
    init_q = .false.
    if(init == 1) init_q = .true.
    initialize_q = 1

  end function initialize_q

  function qfun (pbar)
  
    use splines
    real :: pbar, qfun
    type (spline), save :: spl

    if(init_q) then
       call new_spline(nw, psi_bar, qsf, spl)
    endif

    qfun = splint(pbar, spl)

  end function qfun

  function initialize_pressure (init) 

    integer :: init, initialize_pressure
    
    init_pressure = .false.
    if(init == 1) init_pressure = .true.
    initialize_pressure = 1

  end function initialize_pressure

  function pfun (pbar)
  
    use splines
    real :: pbar, pfun
    type (spline), save :: spl

    if(init_pressure) call new_spline(nw, psi_bar, pressure, spl)

    pfun = splint(pbar, spl) * beta_0/2.

  end function pfun

  function initialize_dpressure (init) 

    integer :: init, initialize_dpressure
    
    init_dpressure = .false.
    if(init == 1) init_dpressure = .true.
    initialize_dpressure = 1

  end function initialize_dpressure

  function dpfun (pbar)
  
    use splines
    real :: pbar, dpfun
    type (spline), save :: spl

    if(init_dpressure) call new_spline(nw, psi_bar, pressure, spl)

    dpfun = dsplint(pbar, spl)/psi_N * beta_0/2.

  end function dpfun

  function initialize_beta (init) 

    integer :: init, initialize_beta
    
    init_beta = .false.
    if(init == 1) init_beta = .true.
    initialize_beta = 1

  end function initialize_beta

  function betafun (pbar)
  
    use splines
    real :: pbar, betafun
    type (spline), save :: spl

    if(init_beta) call new_spline(nw, psi_bar, beta, spl)

    betafun = splint(pbar, spl)

  end function betafun

  function initialize_rho (init) 

    integer :: init, initialize_rho
    
    init_rho = .false.
    if(init == 1) init_rho = .true.
    initialize_rho = 1

  end function initialize_rho

  function rhofun (pbar)
  
    use splines
    real :: pbar, rhofun
    type (spline), save :: spl

    if(init_rho) call new_spline(nrho, psi_mid, rho_mid, spl)

    rhofun = splint(pbar, spl)

  end function rhofun

  function initialize_bound (init) 

    integer :: init, initialize_bound
    
    init_bound = .false.
    if(init == 1) init_bound = .true.
    initialize_bound = 1

  end function initialize_bound

  function bound(theta) 

    use splines
    real :: theta, bound
    type (spline), save :: spl
    integer i

    if(init_bound) call new_spline(nbbbs, thetab, r_bound, spl)
    init_bound = .false.
    
    bound = splint(theta, spl)

  end function bound    

  subroutine alloc_module_arrays(np, nw, nh, nws, nhs, nrho)

  integer :: np, nw, nh, nws, nhs, nrho
 
  allocate (rho_mid(nrho), psi_mid(nrho))
  allocate (psi_bar(np), fp(np), qsf(np), pressure(np), beta(np))
  allocate (dummy(nws), dfit_R(nw), dfit_Z(nh))
  allocate (spsi_bar(nws), sdfit_R(nws), sdfit_Z(nhs))
!  allocate (dfit_psi(nw, nh), sdfit_psi(nws, nhs))
  allocate (dfit_psi(nw, nh))
  allocate (dpm(nw, nh, 2), dtm(nw, nh, 2))

  end subroutine alloc_module_arrays

  subroutine sort(a, b, c, d)

    real, dimension(:) :: a, b, c, d
    real tmp
    integer :: i, j, jmax

    jmax = size(a)

    do j=1,jmax
       do i=1,jmax-j
          if(a(i+1) < a(i)) then
             tmp=a(i); a(i)=a(i+1); a(i+1)=tmp             
             tmp=b(i); b(i)=b(i+1); b(i+1)=tmp
             tmp=c(i); c(i)=c(i+1); c(i+1)=tmp
             tmp=d(i); d(i)=d(i+1); d(i+1)=tmp
          endif
       enddo
    enddo
  end subroutine sort

end module deq
