! This module contains routines to read from CMR's GS2D equilibrium solver

module gs2d
  real, dimension(:), allocatable :: ps,amin_gs2d,q,f,p,pp
  real, dimension(:), allocatable :: rsep,zsep,rgrid,zgrid
  real, dimension(:,:), allocatable :: psi
  real :: r0,a,rmag,zmag
  real :: psmin,psip,b0,ippsi
  integer :: nsurf,nsep,nr,nz

contains

  subroutine read_gs2d(filename)
    implicit none
!CMR, 2/3/2006: lengthen the string filename to 80
    character(len=80), intent(in) :: filename
    character(len=80) :: line
    integer :: i
    open(1,file=trim(filename),status='unknown')
    do i=1,7
       read(1,fmt='(a80)') line
    enddo
    read(1,*) r0,a,rmag,zmag
    read(1,fmt='(a80)') line ; read(1,*) psmin,psip,b0,ippsi
    read(1,fmt='(a80)') line ; read(1,*) nsurf
    allocate(ps(nsurf),amin_gs2d(nsurf),q(nsurf),f(nsurf),p(nsurf),pp(nsurf))
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') ps
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') amin_gs2d
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') q
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') f
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') p
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') pp
    read(1,fmt='(T22,I6)') nsep
    allocate(rsep(nsep),zsep(nsep))
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') rsep
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') zsep
    read(1,fmt='(a80)') line ; read(1,*) nr,nz
    allocate(rgrid(nr),zgrid(nz),psi(nr,nz))
    read(1,fmt='(a80)') line ; read(1,*) rgrid
    read(1,fmt='(a80)') line ; read(1,*) zgrid
    
    read(1,fmt='(a80)') line ; read(1,fmt='(1p,8e16.8)') psi
    close(1)
  end subroutine read_gs2d
  
  subroutine write_gs2d(filename)
    implicit none
    character(len=40), intent(in) :: filename
    integer :: i,j
    open(1,file=filename,status='unknown')
    write(1,fmt='("GS2 input file",T30,"Produced by GS2D at:",a40)') ' '
    write(1,fmt='(A80/,"GS2D Equilibrium Boundary description:",3(/A80))') repeat('-',80),' ',' ',repeat('-',80)
    write(1,fmt='(T2,"r0",T15,"a",T27,"rmag",T39,"zmag (m)"/1p,4e16.8)') r0,a,rmag,zmag
    write(1,fmt='(T2,"psmin",T15,"psedge (Wb)",T27,"b0 (T)",T39,"ip(A)"/1p,4e16.8)') psmin,psip,b0,ippsi
    write(1,fmt='("nfs"/I6)') nsurf
    write(1,fmt='("Psi on 1d grid (for FS quantities) (T)")')
    write(1,fmt='(1p,8e16.8)') ps
    write(1,fmt='("amin_gs2d (m)")')
    write(1,fmt='(1p,8e16.8)') amin_gs2d
    write(1,fmt='("q")')
    write(1,fmt='(1p,8e16.8)') q
    write(1,fmt='("f =r B_phi (Tm)")')
    write(1,fmt='(1p,8e16.8)') f
    write(1,fmt='("p (Pa)")')
    write(1,fmt='(1p,8e16.8)') p
    write(1,fmt='("dp/dpsi (Pa/Wb)")')
    write(1,fmt='(1p,8e16.8)') pp
    write(1,fmt='("No of points on LCFS=",I6)') nsep
    write(1,fmt='("r(j) (m) on LCFS")')
    write(1,fmt='(1p,8e16.8)') rsep
    write(1,fmt='("z(j) (m) on LCFS")')
    write(1,fmt='(1p,8e16.8)') zsep
    write(1,fmt='("NR",T14,"NZ"/2I6)') NR, NZ
    write(1,fmt='("rgrid (m)")')
    write(1,fmt='(1p,8e16.8)') rgrid
    write(1,fmt='("zgrid (m)")')
    write(1,fmt='(1p,8e16.8)') zgrid
    write(1,fmt='("Psi on grid (Wb) : NB Br=(1/2pi r)*dpsi/dz")')
    write(1,fmt='(1p,8e16.8)') psi
    close(1)
  end subroutine write_gs2d
end module gs2d

module eeq

  implicit none
  private

  integer :: nw, nh, nbbbs, ntime, ntstar

  real, allocatable, dimension (:) :: psi_bar, fp, qsf, pressure, beta, spsi_bar
  real, allocatable, dimension (:) :: dummy, efit_R, efit_Z, sefit_R, sefit_Z, efit_t
  real, allocatable, dimension (:,:) :: dum2, efit_psi, sefit_psi
  real, allocatable, dimension (:,:,:) :: dpm, dtm, dum3
  real, allocatable, dimension (:) :: rbbbs, zbbbs, thetab, r_bound !boundary of plasma

  real :: efit_dR, efit_dZ, psi_N

  real :: R_mag, Z_mag, B_T0, aminor, beta_0

  logical :: init_bound = .true.
  logical :: init_btori = .true.
  logical :: init_dbtori = .true.
  logical :: init_q = .true.
  logical :: init_pressure = .true.
  logical :: init_dpressure = .true.
  logical :: init_beta = .true.
  logical :: init_psi = .true.

  public :: efit_init, efitin, gradient, eqitem, bgradient
  public :: gs2din

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


contains
  
  subroutine gs2din(eqfile, psi_0, psi_a, rmaj, B_T, amin, initeq, big)
!
!     This subroutine reads a GS2D output file containing 
!     the axisymmetric magnetic field geometry on a rectangular 
!     domain defined by the coordinates (R,Z).
!
!     efit_R     R grid
!     efit_Z     Z grid
!     fp    F on psibar grid
!     efit_psi   psi on (R,Z) grid
!     R_mag        major radius of the magnetic axis
!     Z_mag        elevation of the magnetic axis
!     rwid       total width of the domain
!     rleft      position of leftmost point of domain
!     zhei       total height of domain
!
    use splines, only: inter_cspl, fitp_surf2, fitp_surf1
    use gs2d
    implicit none

    real :: xdum, p_0
    real :: rwid, rleft, zhei, amin, B_T
    real :: psi_0, psi_a, rmaj, rcentr, bcentr
    real, dimension(:), allocatable :: zp, temp, zx1, zxm, zy1, zyn
    real :: zxy11, zxym1, zxy1n, zxymn
    
    character*80 :: filename, eqfile
    character char*10
    
    integer :: i, j, init, ndum, initeq, big, nhb, nwb, ierr
    
    data init /1/
    save init

    logical:: debug =.false.
    
! Need to generalize initialization condition if equilibrium changes

if (debug) write(6,*) "gs2din: initeq=",initeq
    if(initeq == 0) return
    init=0
    
    i=index(eqfile,' ')-1
    filename = eqfile(1:i)
if (debug) write(6,*) "gs2din: filename=",filename
! read GS2D datafile
    call read_gs2d(filename)
if (debug) write(6,*) "gs2din: read_gs2d done, "
if (debug) write(6,fmt='(T2,"psmin",T15,"psedge (Wb)",T27,"b0 (T)",T39,"ip(A)"/1p,4e16.8)') psmin,psip,b0,ippsi

    nw=nr ; nh=nz
    nwb = nw * big
    nhb = nh * big
    call alloc_module_arrays(nwb, nwb, nhb, nw, nh)
    rwid=maxval(rgrid)-minval(rgrid) ; zhei=maxval(zgrid)-minval(zgrid)
    rcentr=r0 ; rleft=minval(rgrid)
    R_mag=rmag ; Z_mag=zmag ; psi_0=psmin ; psi_a=psip ; bcentr=b0

    psi_0 = psi_0/(8.0*atan(1.))
    psi_a = psi_a/(8.0*atan(1.))
if (debug) write(6,*) "gs2din: psi_0, psi_a=", psi_0, psi_a

!
! pbar is defined by
! pbar = (psi-psi_0)/(psi_a-psi_0)
! fp and q are functions of pbar
!
    sefit_R = rgrid ; sefit_Z = zgrid
    if (size(ps) /= nw) then
       write(6,*) 'gs2din: size(ps) (',size(ps),') /= nw (',nw,')'
       write(6,*) 'gs2din: => should fix the GS2D output file and try again'
       stop
    endif
    spsi_bar=(ps-minval(ps))/(maxval(ps)-minval(ps))

    do i=1,nwb
       psi_bar(i) = float(i-1) / float(nwb-1)
       efit_R(i) = rleft +rwid * float(i-1) / float(nwb-1)
    enddo

! nb Zgrid is not necessarily UpDown symmetric
    do j=1,nhb
       efit_Z(j) = ((float(j-1) / float(nhb-1) ))*zhei+minval(zgrid)
    enddo

    call inter_cspl(nw, spsi_bar, f, nwb, psi_bar, fp)
    call inter_cspl(nw, spsi_bar, p, nwb, psi_bar, pressure)
! divide GS2D psi by 2pi to get poloidal flux in (Wb/rad)
    sefit_psi=psi/(8.0*atan(1.0))

    allocate(zp(3*nw*nh), temp(nw+2*nh))
    allocate(zx1(nh), zxm(nh), zy1(nw), zyn(nw))

    call fitp_surf1(nw, nh, sefit_R, sefit_Z, sefit_psi, &
         nw, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
         255, zp, temp, 1., ierr)

    do j = 1, nhb
       do i = 1, nwb
          efit_psi(i,j) = fitp_surf2(efit_R(i), efit_Z(j), nw, nh, &
               sefit_R, sefit_Z, sefit_psi, nw, zp, 1.)
       enddo
    enddo
      
    deallocate(zp, temp)

    call inter_cspl(nw, spsi_bar, q, nwb, psi_bar, qsf)

    nw = nwb
    nh = nhb

    nbbbs=size(rsep)
    allocate(rbbbs(nbbbs), zbbbs(nbbbs), thetab(nbbbs), r_bound(nbbbs))
    rbbbs=rsep ; zbbbs=zsep

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
       if(thetab(i+1) == thetab(i)) then
!          write(*,*) 'Duplicates near theta = 0 not allowed.'
!          write(*,*) i, i+1, ' Stopping.'
!          stop
!
! put in kluge for duplicate points, which happens near theta=0:
          thetab(i+1) = thetab(i+1)+1.e-8
       endif
    enddo

if (debug) write(6,*) "gs2din:  rbbbs=",rbbbs
if (debug) write(6,*) "gs2din:  zbbbs=",zbbbs
    call a_minor(rbbbs, zbbbs, Z_mag, amin)
    aminor=amin
if (debug) write(6,*) "gs2din:  aminor=",aminor

    deallocate (rbbbs, zbbbs)

    r_bound = r_bound / aminor

    R_mag = R_mag / aminor
    Z_mag = Z_mag / aminor
!    rleft = rleft / aminor
!    rwid = rwid / aminor
!    zhei = zhei / aminor
    rcentr = rcentr / aminor
    efit_R = efit_R / aminor
    efit_Z = efit_Z / aminor
    
! should rmaj be R_mag or rcentr? use R_mag for now.  

    rmaj = R_mag
    B_T0 = abs(bcentr)
    B_T = B_T0
    psi_a = psi_a / (B_T0*aminor**2)
    psi_0 = psi_0 / (B_T0*aminor**2)
    psi_N = psi_a - psi_0

if (debug) write(6,*) "gs2din: B_T0, aminor, psi_0, psi_a=", B_T0, aminor, psi_0, psi_a

    
    p_0=pressure(1)
!    do i=1,nw
!       psi_bar(i) = float(i-1) / float(nw-1)
!       efit_R(i) = rleft +rwid * float(i-1) / float(nw-1)
!    enddo

    fp = fp / (B_T0*aminor)

! MKS: beta = 2 mu_0 p / B**2

    beta = 8. * (2. * acos(0.)) * pressure * 1.e-7 / B_T0**2
    beta_0 = beta(1)

    pressure = pressure / p_0
    
    efit_psi = efit_psi / (B_T0*aminor**2)
    
!    do j=1,nh
!       efit_Z(j) = ((float(j-1) / float(nh-1) ) - 0.5)*zhei
!    enddo

    efit_dR = efit_R(2) - efit_R(1)
    efit_dZ = efit_Z(2) - efit_Z(1)

    1000 format(5(a10),i2,i4,i4)
    2020 format (5e16.9)
    2022 format (2i5)      

  end subroutine gs2din

  subroutine efitin(eqfile, psi_0, psi_a, rmaj, B_T, amin, initeq, big)
!
!     This subroutine reads an EFIT output file containing 
!     the axisymmetric magnetic field geometry on a rectangular 
!     domain defined by the coordinates (R,Z).
!
!     efit_R     R grid
!     efit_Z     Z grid
!     fp    F on psibar grid
!     efit_psi   psi on (R,Z) grid
!     R_mag        major radius of the magnetic axis
!     Z_mag        elevation of the magnetic axis
!     rwid       total width of the domain
!     rleft      position of leftmost point of domain
!     zhei       total height of domain
!
    use splines, only: inter_cspl, fitp_surf2, fitp_surf1
    implicit none

    real :: xdum, p_0
    real :: rwid, rleft, zhei, amin, B_T
    real :: psi_0, psi_a, rmaj, rcentr, bcentr
    real, dimension(:), allocatable :: zp, temp, zx1, zxm, zy1, zyn
    real :: zxy11, zxym1, zxy1n, zxymn
    
    character*80 :: filename, eqfile
    character char*10
    
    integer :: i, j, init, ndum, initeq, big, nhb, nwb, ierr
    
    data init /1/
    save init
    
! Need to generalize initialization condition if equilibrium changes

    if(initeq == 0) return
    init=0
    
    i=index(eqfile,' ')-1
    filename = eqfile(1:i)
    open(unit=5,file=filename,status='old',form='formatted')
    
! Read the data

   read(5,1000) char, char, char, char, char, i, nw, nh
!   write(*,1000) char, char, char, char, char, i, nw, nh

    nwb = nw * big
    nhb = nh * big

    call alloc_module_arrays(nwb, nwb, nhb, nw, nh)

    read(5,2020) rwid, zhei, rcentr, rleft, xdum      
    read(5,2020) R_mag, Z_mag, psi_0, psi_a, bcentr
!    write(*,2020) rwid, zhei, rcentr, rleft, xdum      
!    write(*,2020) R_mag, Z_mag, psi_0, psi_a, bcentr

!
! pbar is defined by
! pbar = (psi-psi_0)/(psi_a-psi_0)
! fp and q are functions of pbar
!
    do i=1,2
       read(5,2020)xdum,xdum,xdum,xdum,xdum
    enddo
    
    do i=1,nw
       spsi_bar(i) = float(i-1) / float(nw-1)
       sefit_R(i) = rleft +rwid * float(i-1) / float(nw-1)
    enddo

    do j=1,nh
       sefit_Z(j) = ((float(j-1) / float(nh-1) ) - 0.5)*zhei
    enddo

    do i=1,nwb
       psi_bar(i) = float(i-1) / float(nwb-1)
       efit_R(i) = rleft +rwid * float(i-1) / float(nwb-1)
    enddo

    do j=1,nhb
       efit_Z(j) = ((float(j-1) / float(nhb-1) ) - 0.5)*zhei
    enddo

    read(5,2020) (dummy(j),   j = 1, nw)
    call inter_cspl(nw, spsi_bar, dummy, nwb, psi_bar, fp)
    read(5,2020) (dummy(j), j = 1, nw)
    call inter_cspl(nw, spsi_bar, dummy, nwb, psi_bar, pressure)
    read(5,2020) (dummy(j),     j = 1, nw)
    read(5,2020) (dummy(j),     j = 1 ,nw)
    read(5,2020) ((sefit_psi(i,j) , i = 1, nw) , j = 1, nh)


    allocate(zp(3*nw*nh), temp(nw+2*nh))
    allocate(zx1(nh), zxm(nh), zy1(nw), zyn(nw))

    call fitp_surf1(nw, nh, sefit_R, sefit_Z, sefit_psi, &
         nw, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
         255, zp, temp, 1., ierr)

    do j = 1, nhb
       do i = 1, nwb
          efit_psi(i,j) = fitp_surf2(efit_R(i), efit_Z(j), nw, nh, &
               sefit_R, sefit_Z, sefit_psi, nw, zp, 1.)
!          write (9, *) efit_R(i), efit_Z(j), efit_psi(i,j)
       enddo
!       write (9, *)
    enddo
      
    deallocate(zp, temp)



    read(5,2020) (dummy(j) ,   j = 1, nw)
    
!    write(*,2020) (dummy(j) ,   j = 1, nw)
!    write (*,*) 
    call inter_cspl(nw, spsi_bar, dummy, nwb, psi_bar, qsf)

    nw = nwb
    nh = nhb
!    read (5, *) nbbbs, ndum
    read(5,2022) nbbbs, ndum
   allocate(rbbbs(nbbbs), zbbbs(nbbbs), thetab(nbbbs), r_bound(nbbbs))
    read(5,2020) (rbbbs(i), zbbbs(i) , i = 1, nbbbs)
!    write (*,*) (rbbbs(i), i=1,nbbbs)
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
       if(thetab(i+1) == thetab(i)) then
!          write(*,*) 'Duplicates near theta = 0 not allowed.'
!          write(*,*) i, i+1, ' Stopping.'
!          stop
!
! put in kluge for duplicate points, which happens near theta=0:
          thetab(i+1) = thetab(i+1)+1.e-8
       endif
    enddo

    call a_minor(rbbbs, zbbbs, Z_mag, amin)
    aminor=amin

!      read(5,2020)(xlim(i),ylim(i),i=1,limitr)
    close(5)

    deallocate (rbbbs, zbbbs)

    r_bound = r_bound / aminor

    R_mag = R_mag / aminor
    Z_mag = Z_mag / aminor
!    rleft = rleft / aminor
!    rwid = rwid / aminor
!    zhei = zhei / aminor
    rcentr = rcentr / aminor
    efit_R = efit_R / aminor
    efit_Z = efit_Z / aminor
    
! should rmaj be R_mag or rcentr? use R_mag for now.  

    rmaj = R_mag
    B_T0 = abs(bcentr)
    B_T = B_T0
    psi_a = psi_a / (B_T0*aminor**2)
    psi_0 = psi_0 / (B_T0*aminor**2)
    psi_N = psi_a - psi_0
    
    p_0=pressure(1)
!    do i=1,nw
!       psi_bar(i) = float(i-1) / float(nw-1)
!       efit_R(i) = rleft +rwid * float(i-1) / float(nw-1)
!    enddo

    fp = fp / (B_T0*aminor)

! MKS: beta = 2 mu_0 p / B**2

    beta = 8. * (2. * acos(0.)) * pressure * 1.e-7 / B_T0**2
    beta_0 = beta(1)

    pressure = pressure / p_0
    
    efit_psi = efit_psi / (B_T0*aminor**2)
    
!    do j=1,nh
!       efit_Z(j) = ((float(j-1) / float(nh-1) ) - 0.5)*zhei
!    enddo

    efit_dR = efit_R(2) - efit_R(1)
    efit_dZ = efit_Z(2) - efit_Z(1)

    if (.true.) then
      write (*,*) "Finished efitin... imported EFIT equilibrium"
      write (*,*) 'Some important quantities:'
      write (*,*) "aminor", aminor
      write (*,*) 'R_mag', R_mag
      write (*,*) 'B_T0', B_T0
      write (*,*) 'beta', beta_0
    end if

    1000 format(5(a10),i2,i4,i4)
    2020 format (5e16.9)
    2022 format (2i5)      

  end subroutine efitin

  subroutine efit_init

    real, dimension(nw, nh) :: eqth 
    integer :: i, j
!cmr nov04: adding following debug switch
    logical :: debug=.false.
!cmr

if (debug) write(6,*) "efit_init: do i"     
    do i = 1, nw
       do j = 1,nh
          if(efit_Z(j) == Z_mag .and. efit_R(i) == R_mag) then
             eqth(i,j) = 0.  ! value should not matter
          else
             eqth(i,j) = atan2( (efit_Z(j)-Z_mag), (efit_R(i)-R_mag))
          endif
       enddo
    enddo

    call derm(efit_psi, dpm)
    call tderm(eqth, dtm)
    
  end subroutine efit_init

  subroutine tderm(f, dfm)

    implicit none
    integer i, j
    real f(:,:), dfm(:,:,:), pi

    pi = 2.*acos(0.)
    
! EFIT grid is equally spaced in R, Z -- this routine uses that fact and 
! is therefore not completely general.  It is fine for EFIT output.    

    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))/efit_dR
    
    i=nw
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))/efit_dR
   
    j=1
    dfm(:,j,2) = -0.5*(3*f(:,j)-4*f(:,j+1)+f(:,j+2))/efit_dZ
    
    j=nh      
    dfm(:,j,2) = 0.5*(3*f(:,j)-4*f(:,j-1)+f(:,j-2))/efit_dZ
    
    do i=2,nw-1
       dfm(i,:,1)=0.5*(f(i+1,:)-f(i-1,:))/efit_dR
    enddo
    
    do j=2,nh-1
       do i = 1,nw
          if(f(i,j+1)-f(i,j-1) > pi) then
             dfm(i,j,2)=0.5*(f(i,j+1)-f(i,j-1)-2.*pi)/efit_dZ
          else
             dfm(i,j,2)=0.5*(f(i,j+1)-f(i,j-1))/efit_dZ
          endif
       enddo
    enddo
    
  end subroutine tderm

  subroutine derm(f, dfm)

    implicit none
    integer i, j
    real f(:,:), dfm(:,:,:)
    
! EFIT grid is equally spaced in R, Z -- this routine uses that fact and 
! is therefore not completely general.  It is fine for EFIT output.    

    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))/efit_dR
    
    i=nw
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))/efit_dR
   
    j=1
    dfm(:,j,2) = -0.5*(3*f(:,j)-4*f(:,j+1)+f(:,j+2))/efit_dZ
    
    j=nh      
    dfm(:,j,2) = 0.5*(3*f(:,j)-4*f(:,j-1)+f(:,j-2))/efit_dZ
    
    do i=2,nw-1
       dfm(i,:,1)=0.5*(f(i+1,:)-f(i-1,:))/efit_dR
    enddo
    
    do j=2,nh-1
       dfm(:,j,2)=0.5*(f(:,j+1)-f(:,j-1))/efit_dZ
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

    if(r_pos >= efit_R(nw) .or. r_pos <= efit_R(1)) then
       write(*,*) 'No evaluation of eqitem allowed outside'
       write(*,*) 'or on edge of R domain'
       write(*,*) r, thetin, efit_R(nw), r_pos
       stop      
    endif

! ensure point is on Z mesh

    if(z_pos >= efit_Z(nh) .or. z_pos <= efit_Z(1)) then
       write(*,*) 'No evaluation of eqitem allowed outside'
       write(*,*) 'or on edge of Z domain'
       write(*,*) r, thetin, efit_Z(1), efit_Z(nh), z_pos
       stop
    endif
    
    istar=0
    do i=2,nw
       if(istar /= 0) cycle
       if(r_pos < efit_R(i)) then
          dr = r_pos - efit_R(i-1)
          sr = efit_R(i) - r_pos
          istar=i-1
       endif
    enddo
      
! Now do Z direction

    jstar=0
    do j=1,nh
       if(jstar /= 0) cycle
       if(z_pos < efit_Z(j)) then
          dt = z_pos - efit_Z(j-1)
          st = efit_Z(j) - z_pos
          jstar=j-1
       endif
    enddo
      
! use opposite area stencil to interpolate

    fstar=f(istar    , jstar    ) * sr * st &
         +f(istar + 1, jstar    ) * dr * st &
         +f(istar    , jstar + 1) * sr * dt &
         +f(istar + 1, jstar + 1) * dr * dt
    fstar = fstar &
         /abs(efit_R(istar+1)-efit_R(istar)) &
         /(efit_Z(jstar+1)-efit_Z(jstar))
!     write(*,*) i, dr, dt, sr, st
!     write(*,*) f(istar,jstar+1),f(istar+1,jstar+1)
!     write(*,*) f(istar,jstar),f(istar+1,jstar)
!     write(*,*) efit_R(istar),efit_R(istar+1)
!     write(*,*) efit_Z(jstar),efit_Z(jstar+1)

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

    call eqitem(r, theta, efit_psi, f)
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

  subroutine alloc_module_arrays(np, nw, nh, nws, nhs, ntime)

    integer :: np, nw, nh, nws, nhs
    integer, intent (in), optional :: ntime
 
    allocate(psi_bar(np), fp(np), qsf(np), pressure(np), beta(np))
    allocate(dummy(nws), efit_R(nw), efit_Z(nh))
    allocate(spsi_bar(nws), sefit_R(nws), sefit_Z(nhs))
    allocate(efit_psi(nw, nh), sefit_psi(nws, nhs))
    allocate(dpm(nw, nh, 2), dtm(nw, nh, 2))

    if (present(ntime)) then
       allocate (dum2(ntime,nws), dum3(nws, nhs, ntime))
       allocate (efit_t(ntime))
    end if

  end subroutine alloc_module_arrays

  subroutine a_minor(r, z, Z_mag, a)

    use splines
    real, dimension(:), intent (in) :: r, z
    real :: a, Z_mag, r1, r2
    integer, parameter :: nz = 5
    real, dimension(nz) :: rtmp, ztmp
    integer i, j, i1, n
    type (spline) :: spl
!CMR, 28/10/08: add code to avoid duplicate points in 5 point spline 
!               to determine r2 on inboard mid-plane
    integer :: k = 0 
    logical:: debug=.false.
    n = size(r)

    if (debug) write(6,*) "aminor:"
    if (debug) write(6,fmt='("r=",10f8.4)') r
    if (debug) write(6,fmt='("z=",10f8.4)') z

    if(n < nz) then
       write(*,*) 'nbbbs < nz -- very strange.  Stopping.'
       write(*,*) 'Look in eeq.f90.'
!       stop
    endif

    j = 0
    do i = nz/2+1,1,-1
       j = j + 1
       ztmp(j) = z(i)
       rtmp(j) = r(i)
    enddo

    if(r(1) == r(n) .and. z(1) == z(n)) then
       do i = n-1, n-nz/2, -1
          j = j + 1
          ztmp(j) = z(i)
          rtmp(j) = r(i)
       enddo
    else
       do i = n, n-nz/2+1, -1
          j = j + 1
          ztmp(j) = z(i)
          rtmp(j) = r(i)
       enddo
    endif

    if (debug) write(6,fmt='("rtmp=",5f8.4)') rtmp
    if (debug) write(6,fmt='("ztmp=",5f8.4)') ztmp
    if (debug) write(6,fmt='("Z_mag=",f8.4)') Z_mag
    
    call new_spline(nz, ztmp, rtmp, spl)
    r1 = splint(Z_mag, spl)    
    call delete_spline(spl)

! find point near magnetic axis elevation on low field side

    do i = nz, n
       if(z(i)-Z_mag > 0.) then
          i1 = i - 1
          exit
       endif
    enddo

!CMR, 28/10/08: modify this code to avoid duplicate points in 5 point spline 
!               to determine r2 on inboard mid-plane
    do i = 1, nz
       rtmp(i) = r(i1 - nz/2 + i - 1 + k )
       ztmp(i) = z(i1 - nz/2 + i - 1 + k )
       if ( i.gt.1 ) then
          if ((rtmp(i)-rtmp(i-1))**2+(ztmp(i)-ztmp(i-1))**2 .lt. 1.0e-7) then
             k=k+1
             if (debug) write(6,fmt='("a_minor: duplicate pt, set k=",i2)') k
             rtmp(i) = r(i1 - nz/2 + i - 1 + k)
             ztmp(i) = z(i1 - nz/2 + i - 1 + k)
          endif
       endif
    enddo

    if (debug) write(6,fmt='("a_minor: rtmp=",5f8.4)') rtmp
    if (debug) write(6,fmt='("a_minor: ztmp=",5f8.4)') ztmp
    if (debug) write(6,fmt='("a_minor: Z_mag=",f8.4)') Z_mag

    call new_spline(nz, ztmp, rtmp, spl)
    r2 = splint(Z_mag, spl)
    call delete_spline(spl)

    a = (r2 - r1)/2.
    if (debug) write(6,*) "a_minor: r1,r2=",r1,r2
    if (debug) write(6,*) "a_minor: a=",a

  end subroutine a_minor

  subroutine sort(a, b, c, d)

    real, dimension(:) :: a, b, c, d
    real tmp
    integer :: i, j, jmax
    logical :: sorted

    jmax = size(a)

    do 
       sorted = .true.
       do i=1,jmax-1
          if(a(i+1) < a(i)) then
             tmp=a(i); a(i)=a(i+1); a(i+1)=tmp             
             tmp=b(i); b(i)=b(i+1); b(i+1)=tmp
             tmp=c(i); c(i)=c(i+1); c(i+1)=tmp
             tmp=d(i); d(i)=d(i+1); d(i+1)=tmp
             sorted = .false.
          endif
       enddo
       if (sorted) exit
    enddo

  end subroutine sort

! alternative coding
!  subroutine sort(a, b, c, d)
!
!    real, dimension(:) :: a, b, c, d
!    real tmp
!    integer :: i, j, jmax
!
!    jmax = size(a)
!
!    do j=1,jmax
!       do i=1,jmax-j
!          if(a(i+1) < a(i)) then
!             tmp=a(i); a(i)=a(i+1); a(i+1)=tmp             
!             tmp=b(i); b(i)=b(i+1); b(i+1)=tmp
!             tmp=c(i); c(i)=c(i+1); c(i+1)=tmp
!             tmp=d(i); d(i)=d(i+1); d(i+1)=tmp
!          endif
!       enddo
!    enddo
!  end subroutine sort
!
end module eeq
