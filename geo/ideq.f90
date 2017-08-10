module ideq

  implicit none
  private

  integer :: nr, nt

  real, allocatable, dimension (:) :: rho_d, eqpsi, psi_bar, fp, qsf, beta, pressure, diam, rc
  real, allocatable, dimension (:) :: rho_mid, psi_mid
  real, allocatable, dimension (:,:)   :: R_psi, Z_psi, B_psi
  real, allocatable, dimension (:,:,:) :: dpm, dtm, drm, dzm, dbm, dbtm
  real, allocatable, dimension (:,:,:) :: dpcart, dbcart, dtcart, dbtcart
  real, allocatable, dimension (:,:,:) :: dpbish, dbbish, dtbish, dbtbish

  real :: psi_N, psi_a, psi_0
  real :: R_mag, Z_mag, B_T0, aminor, beta_0

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
  logical :: init_invR = .true.
  logical :: init_rcenter = .true.
  logical :: init_diameter = .true.

  public :: B_psi
  public :: dfit_init, idfitin, gradient, eqitem, bgradient

  public :: invR
  public :: Rpos
  public :: Zpos
  public :: btori,    initialize_btori
  public :: dbtori,   initialize_dbtori
  public :: qfun,     initialize_q
  public :: pfun,     initialize_pressure
  public :: dpfun,    initialize_dpressure
  public :: betafun,  initialize_beta
  public :: psi,      initialize_psi
  public :: diameter, initialize_diameter
  public :: rcenter,  initialize_rcenter

contains
  
  subroutine idfitin(eqfile, theta, psi_0_out, psi_a_out, rmaj, B_T, amin, initeq)

    use splines
    implicit none

    type :: grid_type
       integer :: nt
       real, dimension (:), pointer :: R => null ()
       real, dimension (:), pointer :: Z => null ()
       real, dimension (:), pointer :: B => null ()
       real, dimension (:), pointer :: theta => null ()
       type (periodic_spline) :: spl 
    end type grid_type

    type :: eq_type
       integer :: nr
       real, dimension (:), pointer :: pbar => null ()
       real, dimension (:), pointer :: psi => null ()
       real, dimension (:), pointer :: pressure => null ()
       real, dimension (:), pointer :: volume => null ()
       type (grid_type), dimension (:), pointer :: grid => null ()
    end type eq_type

    type (eq_type) :: d

    character*80, intent (in) :: eqfile
    real, dimension(:), intent (in) :: theta
    real, intent (out) :: rmaj, amin, psi_a_out, psi_0_out, B_T
    integer, intent (in) :: initeq

    real :: xdum, p_0
    real :: I_ring, twopi
    real :: dpsidr, dpsidz
    
    integer :: i, j, init, ierr
    integer :: idum, jj, nthg
    integer :: jmin, jmax
    
    character*80 :: filename
    character char*10
    character (200) :: line
    
    data init /1/
    save init
    
! Need to generalize initialization condition if equilibrium changes

    if(initeq == 0) return
    init=0

    nt = size(theta)
    twopi = 8.*atan(1.)
    i=index(eqfile,' ')-1
    filename = eqfile(1:i)
    open(unit=5,file=filename,status='old',form='formatted')
    
! Read the data

    read (5, fmt="(a)") line
    read (5, fmt="(a)") line
    read (5,*) R_mag, Z_mag, I_ring
    read (5, fmt="(a)") line
    read (5, *) d%nr
    read (5, fmt="(a)") line
    read (5, fmt="(a)") line

    allocate (d%pbar     (d%nr))
    allocate (d%psi      (d%nr))
    allocate (d%pressure (d%nr))
    allocate (d%volume   (d%nr))
    allocate (d%grid     (d%nr))

    do i=1,d%nr
       read(5,*) idum, d%pbar(i), d%psi(i), d%pressure(i), xdum, d%volume(i), d%grid(i)%nt
       d%grid(i)%nt = d%grid(i)%nt-1
       jmax = d%grid(i)%nt
       allocate (d%grid(i)%R     (jmax))
       allocate (d%grid(i)%Z     (jmax))
       allocate (d%grid(i)%B     (jmax))
       allocate (d%grid(i)%theta (jmax))
       do j=1,jmax
          read (5,*) d%grid(i)%R(j), d%grid(i)%Z(j), dpsidR, dpsidZ
! make B:
          d%grid(i)%B(j) = sqrt(dpsidR**2 + dpsidZ**2)/(d%grid(i)%R(j)*twopi)
          d%grid(i)%theta(j) = atan2 (d%grid(i)%Z(j) - Z_mag, (d%grid(i)%R(j) - R_mag))
       end do
       read (5,*) xdum, xdum, xdum, xdum
    end do
          
! different from tokamak case
! here, psi_0 is the value of psi on the innermost flux surface that is kept in equilibrium
! and psi_a is the last flux surface kept.  
    psi_0 = d%psi(1)
    psi_a = d%psi(d%nr)

    do i=1,d%nr
       call sort (d%grid(i)%theta, d%grid(i)%R, d%grid(i)%Z, d%grid(i)%B)
    end do

! Spline the data onto a regular grid for later use.    

!     R_psi(1:nr, 1:nt)
!     Z_psi(1:nr, 1:nt)
!     B_psi(1:nr, 1:nt)
!
! nt determined from input file by user

!!! need to define our theta grid.  still not right

    call alloc_arrays (d%nr, nt)

    do i = 1, d%nr
       call new_periodic_spline (d%grid(i)%nt, d%grid(i)%theta(:), d%grid(i)%R(:), twopi, d%grid(i)%spl)
       do jj = 1, nt
          R_psi (i, jj) = periodic_splint (theta(jj), d%grid(i)%spl)
       end do
       call delete_periodic_spline (d%grid(i)%spl)
    end do

    do i = 1, d%nr
       call new_periodic_spline (d%grid(i)%nt, d%grid(i)%theta(:), d%grid(i)%Z(:), twopi, d%grid(i)%spl)
       do jj = 1, nt
          Z_psi (i, jj) = periodic_splint (theta(jj), d%grid(i)%spl)
       end do
       call delete_periodic_spline (d%grid(i)%spl)
    end do

    do i = 1, d%nr
       call new_periodic_spline (d%grid(i)%nt, d%grid(i)%theta(:), d%grid(i)%B(:), twopi, d%grid(i)%spl)
       do jj = 1, nt
          B_psi (i, jj) = periodic_splint (theta(jj), d%grid(i)%spl)
       end do
       call delete_periodic_spline (d%grid(i)%spl)
    end do

    qsf = 0.
    rmaj = 1.
! B at center of ring is mu_0 I / (2 a)
    B_T0 = 4.*3.1415*1.e-7*I_ring*0.5/R_mag

    write(*,*) 'B_T0 = ',B_T0
    B_T = B_T0

    psi_N = B_T0 * R_mag**2 
    psi_a_out = psi_a / psi_N
    psi_0_out = psi_0 / psi_N

    R_psi = R_psi/R_mag
    Z_psi = Z_psi/R_mag
    B_psi = B_psi/B_T0

    eqpsi = d%psi / psi_N

    fp = 0.
    nthg = nt

!!! need to define the rho_d grid    

    rho_d = 1.-R_psi(:,1)

! recover some memory: 

    do i = 1, d%nr
       deallocate (d%grid(i)%R)
       deallocate (d%grid(i)%Z)
       deallocate (d%grid(i)%B)
       deallocate (d%grid(i)%theta)
    end do
    deallocate (d%grid)

  end subroutine idfitin

  subroutine alloc_arrays(nr, nt)

    integer :: nr, nt
	
	if(.not.allocated(rho_d)) then
    allocate(rho_d(nr), eqpsi(nr), psi_bar(nr), fp(nr), qsf(nr), beta(nr), pressure(nr), &
         rc(nr), diam(nr))
    allocate(R_psi(nr, nt), Z_psi(nr, nt), B_psi(nr, nt))
    allocate(drm(nr, nt, 2), dzm(nr, nt, 2), dbm(nr, nt, 2), dbtm(nr, nt, 2), &
         dpm(nr, nt, 2), dtm(nr, nt, 2))
    allocate(dpcart(nr, nt, 2), dbcart(nr, nt, 2), dtcart(nr, nt, 2), dbtcart(nr, nt, 2))
    allocate(dpbish(nr, nt, 2), dbbish(nr, nt, 2), dtbish(nr, nt, 2), dbtbish(nr, nt, 2))
	endif
  end subroutine alloc_arrays

  subroutine dfit_init

    real, dimension(nr, nt) :: eqpsi1, eqth 
    integer :: i, j
    real pi

    pi=2*acos(0.)
    do j=1,nt
       do i=1,nr
          eqpsi1(i,j) = eqpsi(i)
          eqth(i,j) = (j-1)*pi/float(nt-1)
       enddo
    enddo
    
    call derm(eqth,   dtm,  'T')
    call derm(R_psi,  drm,  'E')
    call derm(Z_psi,  dzm,  'O')
    call derm(B_psi,  dbm,  'E')
    call derm(eqpsi1, dpm,  'E')
    
! grad(psi) in cartesian form 
    call eqdcart(dpm, dpcart)
! grad(psi) in Bishop form 
    call eqdbish(dpcart, dpbish)

! grad(B) in cartesian form
    call eqdcart(dbm, dbcart)
! grad(B) in Bishop form
    call eqdbish(dbcart, dbbish)

! grad(BT) in cartesian form
    dbtcart = 0.
! grad(BT) in Bishop form
    dbtbish = 0.

! grad(theta) in cartesian form
    call eqdcart(dtm, dtcart)
! grad(theta) in Bishop form
    call eqdbish(dtcart, dtbish)

    


!    do i = 1, nw
!       do j = 1,nh
!          if(dfit_Z(j) == Z_mag .and. dfit_R(i) == R_mag) then
!             eqth(i,j) = 0.  ! value should not matter
!          else
!             eqth(i,j) = atan2( (dfit_Z(j)-Z_mag), (dfit_R(i)-R_mag))
!          endif
!       enddo
!    enddo

!    call derm(dfit_psi, dpm)
!    call tderm(eqth, dtm)
    
  end subroutine dfit_init

  subroutine derm(f, dfm, char)

    implicit none
    integer i, j
    character*1 :: char
    real f(:,:), dfm(:,:,:), pi
    
    pi = 2.*acos(0.)
    
    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))         
    
    i=nr
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))
   
!!! this is a bad assumption
! assume up-down symmetry for now:
 
    select case (char)
    case ('E') 
       j=1
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j+1))
       
       j=nt      
       dfm(:,j,2) = -0.5*(f(:,j-1)-f(:,j-1))
    case ('O')
       j=1
       dfm(:,j,2) = 0.5*(f(:,j+1)+f(:,j+1))
       
       j=nt      
       dfm(:,j,2) = -0.5*(f(:,j-1)+f(:,j-1))
    case ('T')
       j=1
       dfm(:,j,2) = f(:,j+1)
       
       j=nt      
       dfm(:,j,2) = pi - f(:,j-1)
    end select

    do i=2,nr-1
       dfm(i,:,1)=0.5*(f(i+1,:)-f(i-1,:))
    enddo
    
    do j=2,nt-1
       dfm(:,j,2)=0.5*(f(:,j+1)-f(:,j-1))
    enddo
   
  end subroutine derm

  subroutine gradient(rgrid, theta, grad, char, rp, nth_used, ntm)

    use splines, only: inter_d_cspl
    implicit none

    integer, intent (in) :: nth_used, ntm
    character*1, intent (in) :: char
    real, dimension(-ntm:), intent(in) :: rgrid, theta
    real, dimension(-ntm:,:), intent(out) :: grad
    real, dimension (nr, nt, 2) :: dcart
    real, dimension (2) :: tmp
    real, dimension (1) :: aa, daa, rpt
    real :: rp
    integer i
    
    select case(char)
    case('B') 
       dcart = dbcart
    case('D')  ! diagnostic 
       dcart = dbtcart
    case('P') 
       dcart = dpcart
    case('R') 
       dcart = dpcart  ! dpcart is correct for 'R'
    case('T')
       dcart = dtcart
    end select
    
    do i=-nth_used,-1
       call eqitem(rgrid(i), theta(i), dcart(:,:,1), tmp(1), 'R')
       call eqitem(rgrid(i), theta(i), dcart(:,:,2), tmp(2), 'Z')
       if(char == 'T') then
          grad(i,1)=-tmp(1)
          grad(i,2)=-tmp(2)
       else
          grad(i,1)=tmp(1)
          grad(i,2)=tmp(2)
       endif
    enddo

    do i=0,nth_used
       call eqitem(rgrid(i), theta(i), dcart(:,:,1), tmp(1), 'R')
       call eqitem(rgrid(i), theta(i), dcart(:,:,2), tmp(2), 'Z')
       grad(i,1)=tmp(1)
       grad(i,2)=tmp(2)
    enddo

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1) * 0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1) * 0.5*beta_0
       enddo
    endif

  end subroutine gradient

  subroutine bgradient(rgrid, theta, grad, char, rp, nth_used, ntm)

    use splines, only: inter_d_cspl
    implicit none
    
    integer, intent (in) :: nth_used, ntm
    character*1, intent (in) :: char
    real, dimension (-ntm:), intent (in)  ::  rgrid, theta
    real, dimension (-ntm:,:), intent (out) :: grad
    real, dimension (2) :: tmp
    real, dimension (1) :: aa, daa, rpt
    real :: rp
    real, dimension(nr, nt, 2) ::  dbish
    integer i
    logical :: first = .true.
 
    select case(char)
    case('B') 
       dbish = dbbish
    case('D')  ! diagnostic
       dbish = dbtbish
    case('P') 
       dbish = dpbish
    case('R') 
       dbish = dpbish  ! dpcart is correct for 'R'
    case('T')
       dbish = dtbish
    end select

    do i=-nth_used,nth_used
       call eqitem(rgrid(i), theta(i), dbish(:,:,1), grad(i,1), 'R')
       call eqitem(rgrid(i), theta(i), dbish(:,:,2), grad(i,2), 'Z')
    enddo

    if (char == 'T') then
       where (theta(-nth_used:nth_used) < 0.0)
          grad(-nth_used:nth_used,1) = -grad(-nth_used:nth_used,1)
          grad(-nth_used:nth_used,2) = -grad(-nth_used:nth_used,2)
       end where
    end if

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1) * 0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1) * 0.5*beta_0
       enddo
    endif

  end subroutine bgradient

  subroutine eqitem(r, theta_in, f, fstar, char)
      
    integer :: i, j, istar, jstar
    character*1 :: char
    real :: r, thet, fstar, sign, tp, tps, theta_in
    real :: st, dt, sr, dr, pi, rt
    real, dimension(:,:) :: f
    real, dimension(size(f,2)) :: mtheta
    
    pi = 2.*acos(0.)

! check for axis evaluation
      
    if(r == eqpsi(1)) then
       write(*,*) 'no evaluation at axis allowed in eqitem'
       write(*,*) r, theta_in, eqpsi(1)
       stop
    endif
    
! allow psi(r) to be a decreasing function

    sign=1.
    if(eqpsi(2) < eqpsi(1)) sign=-1.
    
    if(r < sign*eqpsi(1)) then
       write(*,*) 'r < Psi_0 in eqitem'
       write(*,*) r,sign,eqpsi(1)
       stop
    endif
      
! find r on psi mesh

! disallow evaluations outside the plasma surface for now

	if(r == eqpsi(nr)) then
		rt = 0.9999999999*r
	if(rt > eqpsi(nr)) rt = 1.0000000001*r
	r = rt
	endif

    if(r >= eqpsi(nr)) then
       write(*,*) 'No evaluation of eqitem allowed outside surface'
       write(*,*) r, theta_in, eqpsi(nr), sign
       stop      
    endif
    
    istar=0
    do i=2,nr
       if(r < sign*eqpsi(i)) then
          dr = r - sign*eqpsi(i-1)
          sr = sign*eqpsi(i) - r
          istar=i-1
          exit
       endif
    enddo
    
! no gradients available at axis, so do not get too close

    if(istar == 1) then
!       write(*,*) 'Too close to axis in eqitem'
!       write(*,*) r, theta_in, eqpsi(1), eqpsi(2)
!       stop
    endif
  
! Now do theta direction

    thet = mod2pi(theta_in)

! assume up-down symmetry

    tp=abs(thet)
    tps=1.
    if(char == 'Z' .and. tp >= 1.e-10) tps=thet/tp
        
! get thet on theta mesh

    do j=1,nt
       mtheta(j)=(j-1)*pi/float(nt-1)
    enddo
  
! note that theta(1)=0 for gen_eq theta 

    jstar=-1
    do j=1,nt
       if(jstar /= -1) cycle
       if(tp < mtheta(j)) then
          dt = tp - mtheta(j-1)
          st = mtheta(j) - tp
          jstar=j-1
       endif
    enddo
      
! treat theta = pi separately
  
    if(jstar == -1) then
       jstar=nt-1
       dt=mtheta(jstar+1)-mtheta(jstar)
       st=0.
    endif

! use opposite area stencil to interpolate

!    if(char == 'R') i=1
!    if(char == 'Z') i=2
    fstar=f(istar    , jstar    ) * sr * st &
         +f(istar + 1, jstar    ) * dr * st &
         +f(istar    , jstar + 1) * sr * dt &
         +f(istar + 1, jstar + 1) * dr * dt
    fstar=fstar*tps &
         /abs(eqpsi(istar+1)-eqpsi(istar)) &
         /(mtheta(jstar+1)-mtheta(jstar))
!     write(*,*) i, dr, dt, sr, st
!     write(*,*) f(istar,jstar+1),f(istar+1,jstar+1)
!     write(*,*) f(istar,jstar),f(istar+1,jstar)
!     write(*,*) eqpsi(istar),eqpsi(istar+1)
!     write(*,*) mtheta(jstar),mtheta(jstar+1)
      

  end subroutine eqitem

  subroutine eqdcart(dfm, dfcart)
      
    implicit none

    real, dimension (:,:,:), intent(in)  :: dfm
    real, dimension (:,:,:), intent(out) :: dfcart
    real, dimension (size(dfm,1),size(dfm,2)) :: denom
    integer :: i, j
      
    denom(:,:) = drm(:,:,1)*dzm(:,:,2) - drm(:,:,2)*dzm(:,:,1)

    dfcart(:,:,1) =   dfm(:,:,1)*dzm(:,:,2) - dzm(:,:,1)*dfm(:,:,2)
    dfcart(:,:,2) = - dfm(:,:,1)*drm(:,:,2) + drm(:,:,1)*dfm(:,:,2)

    do j=1,nt
       do i=2,nr
          dfcart(i,j,:)=dfcart(i,j,:)/denom(i,j)
       enddo
    enddo    

  end subroutine eqdcart

  subroutine eqdbish(dcart, dbish)

    implicit none
    real, dimension(:, :, :), intent (in) :: dcart
    real, dimension(:, :, :), intent(out) :: dbish
    real, dimension(size(dcart,1),size(dcart,2)) :: denom
    integer :: i, j

    denom(:,:) = sqrt(dpcart(:,:,1)**2 + dpcart(:,:,2)**2)

    dbish(:,:,1) = dcart(:,:,1)*dpcart(:,:,1) + dcart(:,:,2)*dpcart(:,:,2)
    dbish(:,:,2) =-dcart(:,:,1)*dpcart(:,:,2) + dcart(:,:,2)*dpcart(:,:,1)
    
    do j=1,nt
       do i=2,nr
          dbish(i,j,:) = dbish(i,j,:)/denom(i,j)
       enddo
    enddo

  end subroutine eqdbish

  function initialize_invR (init) 

    integer :: init, initialize_invR
    
    init_invR = .false.
    if(init == 1) init_invR = .true.
    initialize_invR = 1

  end function initialize_invR

  function invR (r, theta)
   
    real, intent (in) :: r, theta
    real :: f, invR
    real :: th
    
    th = mod2pi( theta)
    
    call eqitem(r, th, R_psi, f, 'R')
    invR=1./f
    
  end function invR

  function Rpos (r, theta)
   
    real, intent (in) :: r, theta
    real :: f, Rpos
    real :: th
    
    th = mod2pi( theta)
    
    call eqitem(r, th, R_psi, f, 'R')
    Rpos=f
    
  end function Rpos

  function Zpos (r, theta)
   
    real, intent (in) :: r, theta
    real :: f, Zpos
    real :: th
    
    th = mod2pi( theta)
    
    call eqitem(r, th, Z_psi, f, 'Z')
    Zpos=f
    
  end function Zpos

  function initialize_psi (init) 

    integer :: init, initialize_psi
    
    init_psi = .false.
    if(init == 1) init_psi = .true.
    initialize_psi = 1

  end function initialize_psi

  function psi (r, theta)
   
    real, intent (in) :: r, theta
    real :: psi

    psi = r
    
  end function psi

  function mod2pi (theta)
    
    real, intent(in) :: theta
    real :: pi, th, mod2pi
    logical :: out
    
    pi=2.*acos(0.)
    
    if(theta <= pi .and. theta >= -pi) then
       mod2pi = theta
       return
    endif
    
    th=theta
    out=.true.
    do while(out)
       if(th > pi) th = th - 2.*pi
       if(th <-pi) th = th + 2.*pi
       if(th <= pi .and. th >= -pi) out=.false.
    enddo
    mod2pi=th
    
  end function mod2pi
   
  function initialize_diameter (init) 

    integer :: init, initialize_diameter
    
    init_diameter = .false.
    if(init == 1) init_diameter = .true.
    initialize_diameter = 1

  end function initialize_diameter

  function diameter (rp)
  
! not really the diameter in this case.  Instead, return the 
! normalized minor radius, measured inside the ring, in the plane of 
! the ring, starting at the ring and going inward.  

    use splines
    real :: rp, diameter
    type (spline), save :: spl

    if(init_diameter) then
       diam = rho_d
       call new_spline(nr, eqpsi, diam, spl)
       init_diameter = .false.
    endif

    diameter = splint(rp, spl)

  end function diameter

  function initialize_rcenter (init) 

    integer :: init, initialize_rcenter
    
    init_rcenter = .false.
    if(init == 1) init_rcenter = .true.
    initialize_rcenter = 1

  end function initialize_rcenter

  function rcenter (rp)
  
    use splines
    real :: rp, rcenter
    type (spline), save :: spl

    if(init_rcenter) then
       rc(:) = 0.5*(R_psi(:,1)+R_psi(:,(nt-1)/2))
       call new_spline(nr, eqpsi, rc, spl)
       init_rcenter = .false.
    endif

    rcenter = splint(rp, spl)

  end function rcenter

  function initialize_dbtori (init) 

    integer :: init, initialize_dbtori
    
    init_dbtori = .false.
    if(init == 1) init_dbtori = .true.
    initialize_dbtori = 1

  end function initialize_dbtori

  function dbtori (pbar)
  
    real :: pbar, dbtori

    dbtori = 0.

  end function dbtori

  function initialize_btori (init) 

    integer :: init, initialize_btori
    
    init_btori = .false.
    if(init == 1) init_btori = .true.
    initialize_btori = 1

  end function initialize_btori

  function btori (pbar)
  
    real :: pbar, btori

    btori = 0.

  end function btori

  function initialize_q (init) 

    integer :: init, initialize_q
    
    init_q = .false.
    if(init == 1) init_q = .true.
    initialize_q = 1

  end function initialize_q

  function qfun (pbar)
  
    real :: pbar, qfun

    qfun = 0.

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

    if(init_pressure) call new_spline(nr, psi_bar, pressure, spl)
    init_pressure = .false.
!
! p_N would be B**2/mu_0 => p = beta/2 in our units
!
    pfun = 0.5*beta_0*splint(pbar, spl)

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
!
! p_N would be B**2/mu_0 => p = beta/2 in our units
!
    if(init_dpressure) then
       call new_spline(nr, psi_bar, pressure, spl)
       init_dpressure = .false.
    endif

    dpfun = dsplint(pbar, spl)/(psi_a-psi_0) * beta_0/2.

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

    if(pbar == 0.) then
       betafun=beta(1)
       return
    endif

    if(init_beta) call new_spline(nr, psi_bar, beta, spl)
    init_beta = .false.

    betafun = splint(pbar, spl)

  end function betafun

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


end module ideq
