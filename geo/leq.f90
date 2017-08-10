module leq

  implicit none

  integer :: nr, nt, np
  private

  logical :: write_prof_var, read_prof_var
  
  integer :: ntg, ntg2pi

  real :: dpsidrho, dIdrho, bi, dqdr, d2Idr2
  real, dimension (:), allocatable :: d2Rdrdth, d2Zdrdth, grho, gpsi, bmag, dBdrho, d2Bdrdth
  real, dimension (:), allocatable :: dgradpardrho, gradpar, dgradparBdrho, dBdth, gradparb
  real, dimension (:), allocatable :: cvdrift0, gbdrift0, dcvdrift0drho, dgbdrift0drho, theta, thetatot
  real, dimension (:), allocatable :: varthet, dvarthdr, gradrho_gradthet, cross, d2varthdr2
  real, dimension (:), allocatable :: gradthet2, gradalph_gradthet, gradrho_gradalph, gradalph2
  real, dimension (:), allocatable :: cvdrift, gbdrift, d2Bdr2, d2Rdr2, d2Zdr2, drz, drzdth
  real, dimension (:), allocatable :: d2Rdr2dth, d2Zdr2dth, d2gpsidr2, dcrossdr
  real, dimension (:), allocatable :: dcvdriftdrho, dgbdriftdrho, dgds2dr, dgds21dr, dgds22dr
  real, dimension (:), allocatable :: gds2, gds21, gds22
  real, dimension (:), allocatable :: dgr2dr, dgpsi2dr
  real, dimension (:), allocatable :: dgrgt, dgt2, dgagr, dgagt, dga2
  real, dimension (:,:), allocatable :: Rr, Zr
  real, allocatable, dimension (:)     :: eqpsi, fp, beta, pressure
  real, allocatable, dimension (:,:)   :: R_psi, Z_psi
  real, allocatable, dimension (:,:,:) :: drm, dzm, dbtm, dpm, dtm
  real, allocatable, dimension (:,:,:) :: dpcart, dtcart, dbtcart
  real, allocatable, dimension (:,:,:) :: dpbish, dtbish, dbtbish

  real, dimension (:), allocatable :: jacrho, delthet, djacdrho, djacrdrho
  real, dimension (:), allocatable :: d2jacdr2, dRdrho, dZdrho, dRdth, dZdth

  real, dimension (:), allocatable :: d2R, d2Z

  real :: beta_0
  
  type :: flux_surface
     real :: R_center, R_geo, k, kp, d, dp, r, dr, delp, q, shat, pp, a, ap, &
          betaprim, betadbprim, d2qdr2, d2psidr2
     integer :: nt, np
  end type flux_surface

  type (flux_surface) :: surf

  public :: leq_init, leqin, gradient, eqitem, bgradient, leqcoefs

  public :: invR, Rpos, Zpos, diameter, btori, dbtori,  qfun, pfun, &
       dpfun, betafun, psi, rcenter, dpdrhofun

contains

  subroutine leqin(R0, Ra, k, kp, d, dp, r, dr, s, qq, qs, a, ap, bp, bpp, d2q, &
       d2psidr2, nt_used, np_used, write_profile_variation, read_profile_variation)
        
    real :: R0, Ra, k, kp, d, dp, r, dr, s, qq, qs, a, ap, bp, bpp, d2q, d2psidr2
    integer :: nt_used, np_used

    logical, intent (in) :: write_profile_variation, read_profile_variation
    
    integer :: i
    real :: dum
    character (1000) :: line

    surf%R_center = R0
    surf%R_geo = Ra
    surf%delp = s
    surf%k = k
    surf%kp = kp
    surf%q = qq
    surf%shat = qs
    surf%d = d
    surf%dp = dp
    surf%a = a
    surf%ap = ap
    surf%r = r
    surf%dr = dr
    surf%pp = 0.
    surf%betaprim = bp
    surf%betadbprim = bpp
    surf%d2qdr2 = d2q
    surf%d2psidr2 = d2psidr2

    beta_0 = 1.

    nr = 3
    nt = nt_used
    np = np_used
    ntg2pi = nt-1
    ntg = ntg2pi*(2*np-1)
    if(.not.allocated(beta)) call alloc_arrays(3, nt)
    surf%nt = nt
    surf%np = np
    dqdr = surf%shat*surf%q/surf%r

    write_prof_var = write_profile_variation
    read_prof_var = read_profile_variation

    call leq_init

  end subroutine leqin

  subroutine alloc_arrays(nr, nt)

    integer :: nr, nt

    allocate(eqpsi(nr), fp(nr), beta(nr), pressure(nr))
    ! periodic quantities can be computed on 2*pi grid and replicated
    allocate(R_psi(nr, nt), Z_psi(nr, nt))
    allocate(drm(nr, nt, 2), dzm(nr, nt, 2), dbtm(nr, nt, 2), &
         dpm(nr, nt, 2), dtm(nr, nt, 2))
    allocate(dpcart(nr, nt, 2), dtcart(nr, nt, 2), dbtcart(nr, nt, 2))
    allocate(dpbish(nr, nt, 2), dtbish(nr, nt, 2), dbtbish(nr, nt, 2))
    allocate (Rr(3,-ntg:ntg), Zr(3,-ntg:ntg))
    allocate (jacrho(-ntg:ntg), djacdrho(-ntg:ntg), djacrdrho(-ntg:ntg), d2jacdr2(-ntg:ntg))
    allocate (d2Rdrdth(-ntg:ntg), d2Zdrdth(-ntg:ntg), gpsi(-ntg:ntg))
    allocate (grho(-ntg:ntg), bmag(-ntg:ntg), dBdrho(-ntg:ntg), dgradpardrho(-ntg:ntg), gradpar(-ntg:ntg))
    allocate (d2Bdrdth(-ntg:ntg), dgradparBdrho(-ntg:ntg), dBdth(-ntg:ntg), gradparb(-ntg:ntg))
    allocate (cvdrift0(-ntg:ntg), gbdrift0(-ntg:ntg), dcvdrift0drho(-ntg:ntg), dgbdrift0drho(-ntg:ntg))
    allocate (theta(-ntg:ntg))
    allocate (dRdrho(-ntg:ntg), dZdrho(-ntg:ntg), dRdth(-ntg:ntg), dZdth(-ntg:ntg))
    allocate (gradrho_gradthet(-ntg:ntg), gradthet2(-ntg:ntg), dgr2dr(-ntg:ntg), dgpsi2dr(-ntg:ntg))
    allocate (dgrgt(-ntg:ntg), dgt2(-ntg:ntg), dgagr(-ntg:ntg), dgagt(-ntg:ntg), dga2(-ntg:ntg))
    allocate (d2Rdr2(-ntg:ntg), d2Zdr2(-ntg:ntg), d2Bdr2(-ntg:ntg))
    allocate (drz(-ntg:ntg), drzdth(-ntg:ntg), d2Rdr2dth(-ntg:ntg), d2Zdr2dth(-ntg:ntg))
    allocate (d2gpsidr2(-ntg:ntg))
    allocate (dgds22dr(-ntg:ntg), gds2(-ntg:ntg), gds21(-ntg:ntg), gds22(-ntg:ntg))
    allocate (gradalph_gradthet(-ntg:ntg), gradalph2(-ntg:ntg), gradrho_gradalph(-ntg:ntg))
    allocate (dgds2dr(-ntg:ntg), dgds21dr(-ntg:ntg))
    allocate (dcvdriftdrho(-ntg:ntg), dgbdriftdrho(-ntg:ntg))
    allocate (varthet(-ntg:ntg), dvarthdr(-ntg:ntg), d2varthdr2(-ntg:ntg))
    allocate (cross(-ntg:ntg), cvdrift(-ntg:ntg), gbdrift(-ntg:ntg))
    allocate (dcrossdr(-ntg:ntg))
    allocate (d2R(-ntg:ntg), d2Z(-ntg:ntg))

  end subroutine alloc_arrays

  subroutine leq_init

    implicit none
    real, dimension(nr, nt) :: eqpsi1, eqth, eqbtor
    
    real, dimension (-ntg:ntg) :: d2Rdth2, d2Zdth2

    real dr(3)
    real pi, t, r
    integer i, j

    character (3) :: dum

    pi=2*acos(0.)
    dr(1) = -surf%dr
    dr(2) = 0.
    dr(3) = surf%dr

    ! initialize to zero
    ! will be overwritten if reading in from file
    ! only relevant for profile variation tests
    d2R = 0. ; d2Z = 0.
    
    if (read_prof_var) then
       open (1002,file='RZ.in',status='old')
       do j=-ntg,ntg
          read (1002,'(3e13.5)') theta(j), d2R(j), d2Z(j)
       end do
       close (1002)
    end if
    
    do j=-ntg,ntg
       theta(j) = j*(2*np-1)*pi/real(ntg)
       do i=1,3
          r = surf%r + dr(i)
          Rr(i,j) = Rpos(r,theta(j),j)
          Zr(i,j) = Zpos(r,theta(j),j)
       end do
    end do

    do j=1,nt
       do i=1,nr
          r = surf%r + dr(i)
          t = (j-1)*pi/real(nt-1)
          R_psi(i,j) = Rpos(r, t, j) 
          Z_psi(i,j) = Zpos(r, t, j)
          eqth(i,j) = t
          eqpsi1(i,j) = 1 + dr(i)
          eqbtor(i,j) = surf%r_geo/R_psi(i,j)
       enddo
    enddo
    
    do i=1,nr
       pressure(i) = -dr(i)
    enddo

    eqpsi(:) = eqpsi1(:,1)

    call derm(eqth,   dtm,  'T')
    call derm(R_psi,  drm,  'E')
    call derm(Z_psi,  dzm,  'O')
    call derm(eqbtor, dbtm, 'E')
    call derm(eqpsi1, dpm,  'E')

    allocate (delthet(-ntg:ntg-1))
    ! get delta theta as a function of theta
    delthet = theta(-ntg+1:)-theta(:ntg-1)

    ! get dR/drho and dZ/drho
    call get_drho (Rr, dRdrho)
    call get_drho (Zr, dZdrho)

    ! get dR/dtheta and dZ/dtheta
    call get_dthet(Rr(2,:), dRdth)
    call get_dthet(Zr(2,:), dZdth)

    ! I=Btor*R is a flux function
    ! bi = I/(Btor(psi,theta of Rgeo)*a) = Rgeo/a
    bi = surf%R_geo

    ! get second derivatives of R and Z with respect to theta
    call get_d2dthet2 (Rr(2,:), d2Rdth2)
    call get_d2dthet2 (Zr(2,:), d2Zdth2)
    ! get mixed theta and rho derivatives of R and Z
    call get_dthet (dRdrho, d2Rdrdth)
    call get_dthet (dZdrho, d2Zdrdth)

    ! get the Jacobian of the transformation from (rho,theta,zeta) to (R,Z,zeta)
    ! this is what I call jacr or jacrho in following comments
    ! as opposed to jacobian, which is for tranformation from (psi,theta,zeta) to (R,Z,zeta)
    call get_jacrho

    ! get dpsinorm/drho
    call get_dpsidrho

    ! get |grad rho| and |grad psi|
    call get_gradrho

    ! quantity needed in calculation of dI/drho and djacrho/drho
    drz = (dRdrho*dRdth + dZdrho*dZdth)/jacrho
    call get_dthet (drz, drzdth)

    ! get dI/drho
    call get_dIdrho

    ! get djacobian/drho*dpsi/drho and djacr/drho
    call get_djacdrho

    ! get d2R/drho2 and d2Z/drho2
    call get_d2RZdr2
    
    if (write_prof_var) then
       open (1002,file='RZ.out',status='unknown')
       do j=-ntg,ntg
          write (1002,'(3e13.5)') theta(j), d2Rdr2(j), d2Zdr2(j)
       end do
       close (1002)
    end if
    
    ! get theta derivative of d2R/drho2 and d2Z/drho2
    call get_dthet (d2Rdr2, d2Rdr2dth)
    call get_dthet (d2Zdr2, d2Zdr2dth)

    ! calculate the magnitude of B (normalized by B(psi,theta corresponding to Rgeo))
    ! B/B0 = sqrt(I**2 + |grad psi|**2)/R
    bmag = sqrt(bi**2 + gpsi**2)/Rr(2,:)

    ! get dB/dtheta
    call get_dthet (bmag, dbdth)
    
    ! calculate b . grad theta
    gradpar = dpsidrho/(bmag*jacrho)
    ! b . grad B
    gradparb = gradpar*dBdth

    ! get d|grad rho|^2/drho and d|grad psi|^2/drho
    call get_dgr2dr

    ! get dB/drho and d2B/drho2
    call get_dBdrho

    ! d (b . grad theta) / drho
    dgradpardrho = -gradpar*(dBdrho/bmag + djacdrho/jacrho)

    ! get d/dtheta (dB/drho)
    call get_dthet (dBdrho, d2Bdrdth)
    
    ! d(b . grad B)/drho
    dgradparBdrho = dgradpardrho*dBdth + gradpar*d2Bdrdth

    ! obtain varthet = (I/(q*(dpsi/dr)) * int_0^theta dtheta' jacrho/R^2
    call get_varthet

    ! obtain dvarthet/drho
    call get_dvarthdr

    ! get |grad theta|^2, grad r . grad theta, grad alpha . grad theta, etc.
    call get_graddotgrad

    ! this is (grad alpha x B) . grad theta
    cross = dpsidrho*(gradrho_gradalph*gradalph_gradthet - gradalph2*gradrho_gradthet)

    ! this is bhat/B x (grad B) . grad alpha * 2 * dpsiN/drho
    gbdrift = 2.0*(-dBdrho + cross*dBdth*dpsidrho/bmag**2)
    ! this is bhat/B x (bhat . grad bhat) . grad alpha * 2 * dpsiN/drho
    ! this is assuming betaprim = 4*pi*ptot/B0^2 * (-d ln ptot / drho)
    cvdrift = (gbdrift + 2.0*surf%betaprim/bmag)/bmag

    ! this is 2 *(bhat/B x grad B / B) . (grad q) * dpsiN/drho / (bhat . grad B)
    ! same as usual GS2 definition once bhat . grad B is added in below
    cvdrift0 = -2.*bi*dqdr/bmag**2

    ! this is 2*dpsiN/drho times the rho derivative (bhat/B x grad B / B) . (grad q)
    dcvdrift0drho = cvdrift0*(dgradparbdrho + gradparb*(dIdrho/bi - 2.*dBdrho/bmag - surf%d2psidr2/dpsidrho)) &
         - 2.*bi*gradparb*surf%d2qdr2/bmag**2
    ! this is 2*dpsiN/drho times the rho derivative of (bhat x gradB/B) . (grad q)
    dgbdrift0drho = cvdrift0*bmag*(dgradparbdrho + gradparb*(dIdrho/bi - dBdrho/bmag - surf%d2psidr2/dpsidrho)) &
         - 2.*bi*gradparb*surf%d2qdr2/bmag

    cvdrift0 = cvdrift0*gradparb
    ! this is 2 * dpsiN/drho * (bhat x gradB/B) . (grad q)
    gbdrift0 = cvdrift0*bmag

    ! get d^2I/drho^2 and d^2 Jac / dr^2
    call get_d2Idr2_d2jacdr2

    ! get d^2varhteta/drho^2
    call get_d2varthdr2

    ! get d2B/drho^2
    call get_d2Bdr2

    ! get d/dr [(grad alpha x B) . grad theta]
    call get_dcrossdr

! corrected Jan. 16, 2014
!    dgbdriftdrho = 2.0*(-d2Bdr2 + dpsidrho*(dcrossdr*dBdth+cross*(d2Bdrdth-2.*dBdth))/bmag**2)
!    dcvdriftdrho = (dgbdriftdrho - gbdrift*dBdrho/bmag)/bmag + 2.0*surf%betadbprim/bmag**2 &
!         - 4.0*surf%betaprim/bmag**3

    ! dgbdriftdrho is d/drho (bhat/B x (grad B) . grad alpha) * 2 * dpsiN/drho
    dgbdriftdrho = 2.0*(surf%d2psidr2*dBdrho/dpsidrho - d2Bdr2 &
         + dpsidrho*(dcrossdr*dBdth+cross*(d2Bdrdth-2.*dBdth*dBdrho/bmag))/bmag**2)
!         + dpsidrho*(dcrossdr*dBdth+cross*(d2Bdrdth-2.*dBdth*dBdrho))/bmag**2)
    ! dcvdriftdrho is d/drho (bhat/B x [bhat . grad bhat] . grad alpha) * 2 * dpsiN/drho
    dcvdriftdrho = (dgbdriftdrho - gbdrift*dBdrho/bmag)/bmag &
         + 2.0*surf%betadbprim/bmag**2 - 4.0*surf%betaprim*dBdrho/bmag**3 &
         - 2.0*surf%betaprim*surf%d2psidr2/dpsidrho

    open (1001,file='leq.out',status='unknown')
    write (1001,'(a9,e12.4,a11,e12.4,a11,e12.4)') '#dI/dr: ', dIdrho, 'd2I/dr2: ', d2Idr2, 'dpsi/dr: ', dpsidrho
    write (1001,'(55a13)') '#1.theta', '2.R', '3.dR/dr', '4.d2Rdr2', '5.dR/dth', &
         '6.d2Rdrdth', '7.dZ/dr', '8.d2Zdr2', '9.dZ/dth', '10.d2Zdrdth', &
         '11.bmag', '12.dBdr', '13.d2Bdr2', '14.dB/dth', '15.d2Bdrdth', &
         '16.varthet', '17.dvarthdr', '18.d2varthdr2', '19.jacr', '20.djacrdr', &
         '21.djacdrho', '22.d2jacdr2', '23.grho2', '24.dgr2dr', '25.gthet2', &
         '26.dgt2', '27.grgthet', '28.dgrgt', '29.galphgth', '30.dgagt', &
         '31.grgalph', '32.dgagr', '33.galph2', '34.dga2', '35.cross', &
         '36.dcrossdr', '37.gbdrift0', '38.dgbdrift0', '39.cvdrift0', '40.dcvdrift0', &
         '41.gbdrift', '42.dgbdrift', '43.cvdrift', '44.dcvdrift', '45.drzdth', &
         '46.gradpar', '47.dgpardr', '48.gradparB', '49.dgparBdr', '50.gds2', &
         '51.dgds2dr', '52.gds21', '53.dgds21dr', '54.gds22', '55.dgds22dr'

    do i = -ntg, ntg
       write (1001,'(55e13.4)') theta(i), Rr(2,i),dRdrho(i), d2Rdr2(i), dRdth(i), &
            d2Rdrdth(i), dZdrho(i), d2Zdr2(i), dZdth(i), d2Zdrdth(i), &
            bmag(i), dBdrho(i), d2Bdr2(i), dBdth(i), d2Bdrdth(i), &
            varthet(i), dvarthdr(i), d2varthdr2(i), jacrho(i), djacrdrho(i), &
            djacdrho(i), d2jacdr2(i), grho(i)**2, dgr2dr(i), gradthet2(i), &
            dgt2(i), gradrho_gradthet(i), dgrgt(i), gradalph_gradthet(i), dgagt(i), &
            gradrho_gradalph(i), dgagr(i), gradalph2(i), dga2(i), cross(i), &
            dcrossdr(i), gbdrift0(i), dgbdrift0drho(i), cvdrift0(i), dcvdrift0drho(i), &
            gbdrift(i), dgbdriftdrho(i), cvdrift(i), dcvdriftdrho(i), drzdth(i), &
            gradpar(i), dgradpardrho(i), gradparB(i), dgradparBdrho(i), gds2(i), &
            dgds2dr(i), gds21(i), dgds21dr(i), gds22(i), dgds22dr(i)
    end do
    close (1001)

! below is actually grad(rho) instead of grad(psi),
! and 'cartesian' refers to (R,Z) coordinates -- MAB
! grad(psi) in cartesian form 
    call eqdcart(dpm, dpcart)
! grad(psi) in Bishop form 
    call eqdbish(dpcart, dpbish)

! grad(BT) in cartesian form
    call eqdcart(dbtm, dbtcart)
! grad(BT) in Bishop form
    call eqdbish(dbtcart, dbtbish)

! grad(theta) in cartesian form
    call eqdcart(dtm, dtcart)
! grad(theta) in Bishop form
    call eqdbish(dtcart, dtbish)

  end subroutine leq_init

  ! takes in f(r), with r given at three radial locations
  ! and returns df = df/dr at the middle radius
  subroutine get_drho (f, df)

    implicit none

    real, dimension (:,-ntg:), intent (in) :: f
    real, dimension (-ntg:), intent (out) :: df

    df = 0.5*(f(3,:)-f(1,:))/surf%dr

  end subroutine get_drho

  ! given function f(theta), calculate second derivative
  ! of f with respect to theta
  ! second order accurate, with equal grid spacing assumed
  subroutine get_d2dthet2 (f, d2f)

    implicit none

    real, dimension (-ntg:), intent (in) :: f
    real, dimension (-ntg:), intent (out) :: d2f

    ! assuming equal grid spacing in theta here
    d2f(-ntg+1:ntg-1) = (f(:ntg-2)-2.*f(-ntg+1:ntg-1)+f(-ntg+2:))/delthet(-ntg+1:ntg-1)**2

    ! use periodicity at boundary
    d2f(-ntg) = (f(ntg-1)-2.*f(-ntg)+f(-ntg+1))/delthet(-ntg+1)**2
    d2f(ntg) = d2f(-ntg)

  end subroutine get_d2dthet2

  ! given function f(theta:-pi->pi), calculate theta derivative
  ! second order accurate, with equal grid spacing assumed
  ! assumes periodic in theta -- may need to change this in future
  subroutine get_dthet (f, df)

    implicit none

    real, dimension (-ntg:), intent (in) :: f
    real, dimension (-ntg:), intent (out) :: df

    ! assuming equal grid spacing in theta here
    df(-ntg+1:ntg-1) = (f(-ntg+2:)-f(:ntg-2))/(delthet(:ntg-2)+delthet(-ntg+1:))

    ! use periodicity at boundary
    df(-ntg) = (f(-ntg+1)-f(ntg-1))/(delthet(-ntg)+delthet(ntg-1))
    df(ntg) = df(-ntg)

  end subroutine get_dthet

  subroutine get_jacrho

    implicit none

    ! jacrho = R*(dR/drho * dZ/dtheta - dR/dtheta * dZ/drho)
    jacrho = Rr(2,:)*(dRdrho*dZdth - dRdth*dZdrho)

  end subroutine get_jacrho

  ! get dpsinorm/drho = (I/2*pi*q)*int_0^{2*pi} dthet jacrho/R**2
  subroutine get_dpsidrho

    use constants, only: pi

    implicit none
    
    ! theta_integrate returns integral from 0 -> 2*pi
    call theta_integrate (jacrho(-ntg2pi:ntg2pi)/Rr(2,-ntg2pi:ntg2pi)**2, dpsidrho)

    ! integration done using trapezoidal rule
    dpsidrho = dpsidrho*bi/(2.*pi*surf%q)

  end subroutine get_dpsidrho

  subroutine get_gradrho

    implicit none

    grho = Rr(2,:)*sqrt(dRdth**2 + dZdth**2)/jacrho
    gpsi = grho*dpsidrho

  end subroutine get_gradrho

  subroutine get_dIdrho

    use constants, only: pi

    implicit none

    real :: num1, num2, denom
    real, dimension (:), allocatable :: dum

    allocate (dum(-ntg:ntg)) ; dum = 0.

    dum = jacrho*( 1.0 + (bi/gpsi)**2 ) / Rr(2,:)**2
    call theta_integrate (dum(-ntg2pi:ntg2pi), denom)

    dum = jacrho*( 2.*dRdrho/Rr(2,:) + dqdr/surf%q ) / Rr(2,:)**2
    call theta_integrate (dum(-ntg2pi:ntg2pi), num1)

    ! betaprim below is (4*pi*ptot/B0^2)*(-d ln ptot / drho)
    dum = ( -2.*(dRdth*d2Rdrdth + dZdth*d2Zdrdth)/jacrho &
         + drzdth + surf%betaprim*jacrho/dpsidrho**2 ) / grho**2
    call theta_integrate (dum(-ntg2pi:ntg2pi), num2)

    dIdrho = bi*(num1 + num2)/denom

    deallocate (dum)

  end subroutine get_dIdrho

  subroutine get_djacdrho

    implicit none

    real :: test

    ! this is dpsi/dr * d/dr (jacobian)
    ! betaprim below is (4*pi*ptot/B0^2)*(-d ln ptot / drho)
    djacdrho = (Rr(2,:)/grho)**2*(2.*(dRdth*d2Rdrdth+dZdth*d2Zdrdth)/jacrho &
         - drzdth + jacrho*(bi*dIdrho/Rr(2,:)**2 - surf%betaprim)/dpsidrho**2)

    ! this is d/dr (jacobian_r)
    djacrdrho = djacdrho + jacrho*surf%d2psidr2/dpsidrho

  end subroutine get_djacdrho

  subroutine get_d2RZdr2

    implicit none

    ! get factor common to both d2R/drho2 and d2Z/drho2
    d2Rdr2 = ((djacrdrho-jacrho*dRdrho/Rr(2,:))/Rr(2,:) &
         - dRdrho*d2Zdrdth + dZdrho*d2Rdrdth)/(dRdth**2+dZdth**2)

    d2Zdr2 = -d2Rdr2*dRdth
    d2Rdr2 = d2Rdr2*dZdth

  end subroutine get_d2RZdr2

  subroutine get_dgr2dr

    implicit none

    dgr2dr = 2.*(grho**2*(dRdrho/Rr(2,:)-djacrdrho/jacrho) &
         + (Rr(2,:)/jacrho)**2*(dRdth*d2Rdrdth + d2Zdrdth*dZdth))

    dgpsi2dr = 2.*(gpsi**2*(dRdrho/Rr(2,:)-djacdrho/jacrho) &
         + (Rr(2,:)/jacrho)**2*(dRdth*d2Rdrdth + d2Zdrdth*dZdth)*dpsidrho**2)

  end subroutine get_dgr2dr

  subroutine get_graddotgrad

    implicit none

    ! grad theta . grad theta
    gradthet2 = (Rr(2,:)/jacrho)**2*(dRdrho**2 + dZdrho**2)
    ! grad rho . grad theta
    gradrho_gradthet = -(Rr(2,:)/jacrho)**2*(dRdrho*dRdth+dZdrho*dZdth)

    ! grad alpha . grad theta
    gradalph_gradthet = -(varthet*dqdr + surf%q*dvarthdr)*gradrho_gradthet &
         - bi*jacrho/(dpsidrho*Rr(2,:)**2)*gradthet2
    ! grad rho . grad alpha
    gradrho_gradalph = -(varthet*dqdr + surf%q*dvarthdr)*grho**2 &
         - bi*jacrho/(dpsidrho*Rr(2,:)**2)*gradrho_gradthet
    ! grad alpha . grad alpha
    gradalph2 = (1./Rr(2,:)**2) + ((varthet*dqdr+surf%q*dvarthdr)*grho)**2 &
         + 2.*bi*jacrho*(varthet*dqdr+surf%q*dvarthdr)*gradrho_gradthet/(dpsidrho*Rr(2,:)**2) &
         + (bi*jacrho/(dpsidrho*Rr(2,:)**2))**2*gradthet2

  end subroutine get_graddotgrad

  subroutine get_dBdrho

    implicit none

    integer :: i

    ! dB/drho
    dBdrho = ( bi*dIdrho + 0.5*dgpsi2dr ) / (bmag*Rr(2,:)**2) &
         - bmag*dRdrho/Rr(2,:)

  end subroutine get_dBdrho

  subroutine get_varthet

    implicit none

    call theta_integrate_indef(jacrho/Rr(2,:)**2, varthet)
    varthet = bi*varthet/(dpsidrho*surf%q)

  end subroutine get_varthet

  subroutine get_dvarthdr

    implicit none

    real, dimension (-ntg:ntg) :: dum

    dum = bi*jacrho*( dIdrho/bi - dqdr/surf%q + djacdrho/jacrho &
         - 2.*dRdrho/Rr(2,:) )/Rr(2,:)**2
    call theta_integrate_indef(dum, dvarthdr)
    dvarthdr = dvarthdr/(dpsidrho*surf%q)

  end subroutine get_dvarthdr

  subroutine get_d2Idr2_d2jacdr2

    use constants, only: pi

    implicit none

    integer :: i

    real :: denom, num1, num2, num3, num4
    real, dimension (-ntg:ntg) :: tmp, tmp2

    ! denom is the denominator in the expression for d^2 I / dr^2
    tmp = jacrho/Rr(2,:)**2*(1.0 + (bi/gpsi)**2)
    call theta_integrate (tmp(-ntg2pi:ntg2pi), denom)
    denom = denom/bi

    d2jacdr2 = dIdrho*bi*jacrho/gpsi**2 &
         * (dIdrho/bi + djacrdrho/jacrho - dgpsi2dr/gpsi**2 &
         - 2.*dRdrho/Rr(2,:))
    
    tmp = -d2jacdr2/Rr(2,:)**2 - dIdrho*jacrho/(bi*Rr(2,:)**2) &
         * (djacrdrho/jacrho - dIdrho/bi - 2.*dRdrho/Rr(2,:))
    call theta_integrate (tmp(-ntg2pi:ntg2pi), num1)

    ! tmp = -jacrho/(dpsidrho*Rr(2,:)**2)*(djacdrho/jacrho - 2.*dRdrho/Rr(2,:))
    ! call theta_integrate (tmp(-ntg2pi:ntg2pi), num2)
    ! d2jacdr2 = d2jacdr2 - tmp*Rr(2,:)**2*surf%d2psidr2
    ! num2 = surf%d2psidr2 * (2*pi*surf%q/bi*(dqdr/surf%q - dIdrho/bi) + num2)
    
    tmp = (d2Rdr2*dRdth+dRdrho*d2Rdrdth+d2Zdr2*dZdth+dZdrho*d2Zdrdth)/jacrho &
         - djacrdrho*(dRdrho*dRdth+dZdrho*dZdth)/jacrho**2
    call get_dthet (tmp, tmp2)
    tmp = (tmp2 - 2./jacrho*(-djacrdrho/jacrho*(dRdth*d2Rdrdth+dZdth*d2Zdrdth) &
         + d2Rdrdth**2 + dRdth*d2Rdr2dth + d2Zdrdth**2 + dZdth*d2Zdr2dth))/grho**2 &
         - dgr2dr*(drzdth - 2./jacrho*(dRdth*d2Rdrdth + dZdth*d2Zdrdth))/grho**4
    call theta_integrate (tmp(-ntg2pi:ntg2pi), num2)
    d2jacdr2 = d2jacdr2 - tmp*Rr(2,:)**2

    tmp = jacrho*(surf%betadbprim + surf%betaprim*(djacrdrho/jacrho- dgpsi2dr/gpsi**2))/gpsi**2
    call theta_integrate (tmp(-ntg2pi:ntg2pi), num3)
    d2jacdr2 = d2jacdr2 - tmp*Rr(2,:)**2

    tmp = jacrho/Rr(2,:)**2*(2.*d2Rdr2/Rr(2,:) - 2.*(dRdrho/Rr(2,:))**2 &
         + surf%d2qdr2/surf%q - (dqdr/surf%q)**2 + (2*dRdrho/Rr(2,:) + dqdr/surf%q) &
         * (djacrdrho/jacrho - 2.*dRdrho/Rr(2,:)))
    call theta_integrate (tmp(-ntg2pi:ntg2pi), num4)

    d2Idr2 = (num1+num2+num3+num4)/denom
!    d2jacdr2 = d2jacdr2 + bi*jacrho/(gpsi*Rr(2,:))**2*d2Idr2 + 2.*djacdrho*dRdrho/Rr(2,:)**3
    d2jacdr2 = d2jacdr2 + bi*jacrho/gpsi**2*d2Idr2 + 2.*djacdrho*dRdrho/Rr(2,:)

  end subroutine get_d2Idr2_d2jacdr2

  subroutine get_d2varthdr2

    implicit none

    real, dimension (-ntg:ntg) :: dum

    dum = bi*jacrho/(surf%q*dpsidrho*Rr(2,:)**2)*( (dIdrho/bi - dqdr/surf%q &
!    dum = bi*jacrho/(surf%q*Rr(2,:)**2)*( (dIdrho/bi - dqdr/surf%q &
         + djacdrho/jacrho - 2.*dRdrho/Rr(2,:))**2 &
         + d2Idr2/bi - (dIdrho/bi)**2 - surf%d2qdr2/surf%q &
         + (dqdr/surf%q)**2 + d2jacdr2/jacrho - (djacdrho/jacrho)**2 &
         - djacdrho*surf%d2psidr2/(dpsidrho*jacrho) &
         - 2.*d2Rdr2/Rr(2,:) + 2.*(dRdrho/Rr(2,:))**2 )
    call theta_integrate_indef(dum, d2varthdr2)

  end subroutine get_d2varthdr2

  subroutine get_d2Bdr2

    implicit none

! corrected Jan. 16, 2014
!    d2gpsidr2 = 2.*(dpsidrho*Rr(2,:)/jacrho)**2 &
!         * 2.*(dRdrho/Rr(2,:)-djacdrho/jacrho) &
!         * (dRdrho/Rr(2,:)-djacdrho/jacrho+dRdth*d2Rdrdth+dZdth*d2Zdrdth) &
!         - dRdrho**2/Rr(2,:)**2 + d2Rdr2/Rr(2,:) + djacdrho**2/jacrho**2 &
!         - d2jacdr2/jacrho + d2Rdrdth**2 + dRdth*d2Rdr2dth &
!         + d2Zdrdth**2 + dZdth*d2Zdr2dth
    ! d2gpsidr2 = 2.*( dgr2dr*(dRdrho/Rr(2,:) - djacdrho/jacrho) &
    !      + grho**2*(d2Rdr2/Rr(2,:) - (dRdrho/Rr(2,:))**2 - d2jacdr2/jacrho &
    !      + djacdrho*djacrdrho/jacrho**2) + (Rr(2,:)/jacrho)**2 &
    !      * (dRdth**2 + dRdth*d2Rdr2dth + dZdth**2 + dZdth*d2Zdr2dth &
    !      + 2.*(dRdrho/Rr(2,:) - djacrdrho/jacrho)*(dRdth*d2Rdrdth+dZdth*d2Zdrdth)) )
    d2gpsidr2 = 2.*(dRdrho/Rr(2,:)-djacdrho/jacrho)*dgpsi2dr &
         + 2.*gpsi**2*(d2Rdr2/Rr(2,:)-(dRdrho/Rr(2,:))**2 - d2jacdr2/jacrho + djacdrho*djacrdrho/jacrho**2) &
         + 2.*(Rr(2,:)*gpsi/jacrho)**2*( d2Rdrdth**2 + dRdth*d2Rdr2dth + d2Zdrdth**2 + dZdth*d2Zdr2dth &
         + 2.*(dRdth*d2Rdrdth + dZdth*d2Zdrdth)*(dRdrho/Rr(2,:)-djacdrho/jacrho) )

    ! d2gpsidr2 = 2.*(dpsidrho*Rr(2,:)/jacrho)**2 &
    !      * (2.*(dRdrho/Rr(2,:)-djacdrho/jacrho) &
    !      * ((dRdrho/Rr(2,:)-djacdrho/jacrho)*(dRdth**2+dZdth**2) &
    !      + 2.*(dRdth*d2Rdrdth+dZdth*d2Zdrdth)) &
    !      + (dRdth**2+dZdth**2)*(d2rdr2/Rr(2,:) - (dRdrho/Rr(2,:))**2 &
    !      - d2jacdr2/jacrho + (djacdrho/jacrho)**2) &
    !      + d2Rdrdth**2 + dRdth*d2Rdr2dth + d2Zdrdth**2 + dZdth*d2Zdr2dth) &
    !      + 4.*dpsidrho*surf%d2psidr2*dgr2dr &
    !      + 2.*grho**2*(surf%d2psidr2**2 + dpsidrho*surf%d3psidr3)

    ! get d/drho (dB/drho)
    d2Bdr2 = -dBdrho*dRdrho/Rr(2,:) + bmag*(dRdrho/Rr(2,:))**2 &
         - bmag*d2Rdr2/Rr(2,:) + 0.5*(2.*(dIdrho**2 + bi*d2Idr2) &
         + d2gpsidr2)/(bmag*Rr(2,:)**2) &
         - (dBdrho + bmag*dRdrho/Rr(2,:))*(2.*dRdrho/Rr(2,:)+dBdrho/bmag)

  end subroutine get_d2Bdr2

  subroutine get_dcrossdr

    implicit none

    integer :: i

! corrected Jan. 16, 2014
!    dgr2 = 2.*(Rr(2,:)/jacrho)**2*(dRdrho/Rr(2,:)-djacdrho/jacrho &
!         + dRdth*d2Rdrdth + dZdth*d2Zdrdth)
!    dgrgt = -(Rr(2,:)/jacrho)**2*((2.*dRdrho/Rr(2,:)-2.*djacdrho/jacrho) &
!         + d2Rdr2*dRdth+dRdrho*d2Rdrdth+d2Zdr2*dZdth+dZdrho*d2Zdrdth)
!    dgt2 = (Rr(2,:)/jacrho)**2*(2.*(dRdrho/Rr(2,:)-djacdrho/jacrho) &
!         + 2.*(dRdrho*d2Rdrdth + dZdrho*d2Zdrdth))
!    dga2 = -2*dRdrho/Rr(2,:)**3 + dgr2*(varthet*dqdr+surf%q*dvarthdr) &
!         + (2.0*grho**2*(varthet*dqdr+surf%q*dvarthdr) &
!         + 2.*bi*jacrho*gradrho_gradthet/(dpsidrho*Rr(2,:)**2)) &
!         *(surf%d2qdr2*varthet+2.*dqdr*dvarthdr+surf%q*d2varthdr2) &
!         + 2.*(varthet*dqdr+surf%q*dvarthdr)*bi*jacrho/(dpsidrho*Rr(2,:)**2) &
!         * (dgrgt + gradrho_gradthet*(dIdrho/bi + djacdrho/jacrho - 2.*dRdrho/Rr(2,:))) &
!         + (bi*jacrho/(dpsidrho*Rr(2,:)**2))**2*(dgt2 + gradthet2*(2.*dIdrho/bi + 2.*djacdrho/jacrho &
!         - 4.*dRdrho/Rr(2,:)))
!    dgagr = -grho**2*(dvarthdr*dqdr+2.*varthet*dqdr+surf%q*d2varthdr2) &
!         - dgr2*(varthet*dqdr+surf%q*dvarthdr) - bi*jacrho/(dpsidrho*Rr(2,:)**2) &
!         * (dgrgt + gradrho_gradthet*(dIdrho/bi + djacdrho/jacrho - 2.*dRdrho/Rr(2,:)))

    ! dgr2 = d/drho (|grad rho|^2)
    ! dgr2 = 2.*(Rr(2,:)/jacrho)**2*((dRdrho/Rr(2,:)-djacdrho/jacrho)*(dRdth**2+dZdth**2) &
    !      + dRdth*d2Rdrdth + dZdth*d2Zdrdth)
    ! dgrgt = d/drho (grad rho . grad theta)
!    dgrgt = -(Rr(2,:)/jacrho)**2*(2.*(dRdrho/Rr(2,:)-djacdrho/jacrho)*(dRdrho*dRdth+dZdrho*dZdth) &
!         + d2Rdr2*dRdth+dRdrho*d2Rdrdth+d2Zdr2*dZdth+dZdrho*d2Zdrdth)
    dgrgt = 2.*gradrho_gradthet*(dRdrho/Rr(2,:)-djacrdrho/jacrho) &
         - (Rr(2,:)/jacrho)**2*(d2Rdr2*dRdth + dRdrho*d2Rdrdth + d2Zdr2*dZdth + dZdrho*d2Zdrdth)
    ! dgt2 = d/drho (|grad theta|^2)
    dgt2 = 2.*(Rr(2,:)/jacrho)**2*((dRdrho/Rr(2,:)-djacrdrho/jacrho)*(dRdrho**2+dZdrho**2) &
         + dRdrho*d2Rdr2 + dZdrho*d2Zdr2)
    ! this is d/drho (|grad alph|^2)
    ! will later multiply it by 0.5*dpsidrho**2
    dga2 = -2*dRdrho/Rr(2,:)**3 + dgr2dr*(varthet*dqdr+surf%q*dvarthdr)**2 &
         + (2.0*grho**2*(varthet*dqdr+surf%q*dvarthdr) &
         + 2.*bi*jacrho*gradrho_gradthet/(dpsidrho*Rr(2,:)**2)) &
         * (surf%d2qdr2*varthet+2.*dqdr*dvarthdr+surf%q*d2varthdr2) &
         + 2.*(varthet*dqdr+surf%q*dvarthdr)*bi*jacrho/(dpsidrho*Rr(2,:)**2) &
         * (dgrgt + gradrho_gradthet*(dIdrho/bi + djacdrho/jacrho - 2.*dRdrho/Rr(2,:))) &
         + (bi*jacrho/(dpsidrho*Rr(2,:)**2))**2*(dgt2 + 2.*gradthet2*(dIdrho/bi + djacdrho/jacrho &
         - 2.*dRdrho/Rr(2,:)))

    ! dgagr = d/drho (grad alpha . grad rho)
    dgagr = -grho**2*(2.*dvarthdr*dqdr+varthet*surf%d2qdr2+surf%q*d2varthdr2) &
         - dgr2dr*(varthet*dqdr+surf%q*dvarthdr) - bi*jacrho/(dpsidrho*Rr(2,:)**2) &
         * (dgrgt + gradrho_gradthet*(dIdrho/bi + djacdrho/jacrho - 2.*dRdrho/Rr(2,:)))

    ! dgagt = d/drho (grad alpha . grad theta)
    dgagt = -gradrho_gradthet*(2.*dvarthdr*dqdr+varthet*surf%d2qdr2+surf%q*d2varthdr2) &
         - dgrgt*(varthet*dqdr+surf%q*dvarthdr) - bi*jacrho/(dpsidrho*Rr(2,:)**2) &
         * (dgt2 + gradthet2*(dIdrho/bi + djacdrho/jacrho - 2.*dRdrho/Rr(2,:)))

    ! dcrossdr = d/drho [(grad alpha x B) . grad theta)]
    dcrossdr = dpsidrho*(dgagr*gradalph_gradthet+gradrho_gradalph*dgagt &
         - dga2*gradrho_gradthet - gradalph2*dgrgt) + surf%d2psidr2*cross/dpsidrho

    ! |grad alpha|^2 * (dpsiN/drho)^2 (dpsiN/drho factor accounts for ky normalization)
    gds2 = gradalph2*dpsidrho**2
    ! (grad q . grad alpha) * (dpsiN/drho)^2
    gds21 = gradrho_gradalph*dqdr*dpsidrho**2
    ! |grad q|^2 * (dpsiN/drho)^2
    gds22 = (gpsi*dqdr)**2

    ! this is (dpsi/drho)^2*d|grad alpha|^2/dr
    dgds2dr = dga2*dpsidrho**2
    ! this is (dpsi/drho)^2*d(grad alpha . grad q)/dr
! corrected Jan. 16, 2014
!    dgds21dr = -dgagr*dqdr*dpsidrho**2
    ! note that there will be multiplication by 2 in dist_fn.fpp
    dgds21dr = (dgagr*dqdr+surf%d2qdr2*gradrho_gradalph)*dpsidrho**2
    ! this is (dpsi/drho)^2*d(|grad q|^2)/dr
    dgds22dr = (dqdr**2*dgr2dr + 2.*grho**2*dqdr*surf%d2qdr2)*dpsidrho**2

    ! note that kperp2 = (n0/a)^2*(drho/dpsiN)^2*(gds2 + 2*theta0*gds21 + theta0^2*gds22)
    ! note that dkperp2/dr = (n0/a)^2*(drho/dpsiN)^2*(dgds2dr + 2*theta0*dgds21dr + theta0^2*dgds22dr)

  end subroutine get_dcrossdr

  subroutine theta_integrate (integrand, integral)

    implicit none

    real, dimension (-ntg2pi:), intent (in) :: integrand
    real, intent (out) :: integral

    ! use trapezoidal rule to integrate in theta
    integral = 0.5*sum(delthet(-ntg2pi:ntg2pi-1)*(integrand(-ntg2pi:ntg2pi-1) + integrand(-ntg2pi+1:ntg2pi)))

  end subroutine theta_integrate

  ! get indefinite integral of integrand
  subroutine theta_integrate_indef (integrand, integral)

    implicit none

    real, dimension (-ntg:), intent (in) :: integrand
    real, dimension (-ntg:), intent (out) :: integral

    integer :: i

    ! use trapezoidal rule to integrate in theta
    integral(0) = 0.0
    do i = 1, ntg
       integral(i) = integral(i-1)+0.5*delthet(i-1)*(integrand(i-1)+integrand(i))
    end do
    do i = -1, -ntg, -1
       integral(i) = integral(i+1)-0.5*delthet(i)*(integrand(i+1)+integrand(i))
    end do

  end subroutine theta_integrate_indef

  subroutine derm(f, dfm, char)

    implicit none
    integer i, j
    character(1) :: char
    real f(:,:), dfm(:,:,:), pi

    pi = 2*acos(0.)
    
    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))         
    
    i=nr
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))
   
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
    
    integer nth_used, ntm
    character(1) char
    real rgrid(-ntm:), theta(-ntm:), grad(-ntm:,:)
    real tmp(2), aa(1), daa(1), rp, rpt(1)
    real, dimension(nr,nt,2) :: dcart
    integer i
    
    select case(char)
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
          grad(i,1)=grad(i,1)*daa(1)*0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1)*0.5*beta_0
       enddo
    endif

  end subroutine gradient

  subroutine bgradient(rgrid, theta, grad, char, rp, nth_used, ntm)

    use splines, only: inter_d_cspl
    implicit none
    
    integer nth_used, ntm
    character(1) char
    real rgrid(-ntm:), theta(-ntm:), grad(-ntm:,:)
    real :: aa(1), daa(1), rp, rpt(1)
    real, dimension(nr,nt,2) ::  dbish
    integer i

    select case(char)
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
      
    integer :: j, istar, jstar
    character(1) :: char
    real :: r, thet, fstar, tp, tps, theta_in
    real :: st, dt, pi
    real, dimension(:,:) :: f
    real, dimension(size(f,2)) :: mtheta
    
    pi = 2.*acos(0.)
! find r on psi mesh
    
    istar = 2

! Now do theta direction

    thet = mod2pi(theta_in)

! assume up-down symmetry

    tp=abs(thet)
    tps=1.
    if(char == 'Z' .and. thet /= 0.) tps=thet/abs(thet)
        
! get thet on theta mesh

    mtheta = (/ ( real(j-1)*pi/real(nt-1), j=1,nt) /)
  
! note that theta(1)=0 for local_eq theta 

    jstar=-1    
    do j=1,nt
       if(tp < mtheta(j)) then
          dt = tp - mtheta(j-1)
          st = mtheta(j) - tp
          jstar=j-1
          exit
       endif
       if(jstar /= -1) write(*,*) 'exit error j'
    enddo
      
! treat theta = pi separately
  
    if(jstar == -1) then
       jstar=nt-1
       dt=mtheta(jstar+1)-mtheta(jstar)
       st=0.
    endif

! use opposite area stencil to interpolate

    fstar=f(istar    , jstar    )  * st &
         +f(istar    , jstar + 1)  * dt 
    fstar=fstar*tps &
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

    dfcart = 0.
    
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

  subroutine leqcoefs (dgpdr, dgpbdr, dcvd0dr, dgbd0dr, dcvdr, dgbdr, dg2dr, dg21dr, dg22dr, dbdr, &
       gs2, gs21, gs22, gpar, gbd0, cvd0, gbd, cvd, drdpsi, gr, b, bpol, dkxfacdr)

    implicit none

    real, dimension (-ntg:), intent (out) :: dgpdr, dgpbdr, dcvd0dr, dgbd0dr, dcvdr, dgbdr
    real, dimension (-ntg:), intent (out) :: dg2dr, dg21dr, dg22dr, dbdr
    real, dimension (-ntg:), intent (out) :: gs2, gs21, gs22, gpar, gbd0, cvd0, gbd, cvd
    real, dimension (-ntg:), intent (out) :: gr, b, bpol
    real, intent (out) :: drdpsi, dkxfacdr

    integer :: i

    drdpsi = 1./dpsidrho
    gr = grho
    b = bmag
    bpol = gpsi/Rr(2,:)
    gs2 = gds2
    gs21 = gds21
    gs22 = gds22
    gpar = gradpar
    gbd0 = gbdrift0
    cvd0 = cvdrift0
    gbd = gbdrift
    cvd = cvdrift
    dgpdr = dgradpardrho
    dgpbdr = dgradparbdrho
    dcvd0dr = dcvdrift0drho
    dgbd0dr = dgbdrift0drho
    dcvdr = dcvdriftdrho
    dgbdr = dgbdriftdrho
    dg2dr = dgds2dr
    dg21dr = dgds21dr
    dg22dr = dgds22dr
    dbdr = dBdrho
    ! factor multiplying kx in dkx/dr
    if (dqdr > epsilon(0.)) then 
       dkxfacdr = surf%d2qdr2/dqdr - surf%d2psidr2/dpsidrho
    else
       write (*,*) 'WARNING (if running with profile_variation=T): dq/dr=0 case currently excludes nonlinear piece in dg/dr equation'
       dkxfacdr = 0.
    end if
    ! dgpdr = 0.
    ! dgpbdr = 0.
    ! dcvd0dr = 0.
    ! dgbd0dr = 0.
    ! dcvdr = 0.
    ! dgbdr = 0.
    ! dg2dr = 0.
    ! dg21dr = 0.
    ! dg22dr = 0.
    ! dbdr = 0.

!    do i = -ntg, ntg
!       write (*,'(a9,17e12.4)') 'leqcoefs', theta(i), dgpdr(i), dgpbdr(i), dcvd0dr(i), dgbd0dr(i), dcvdr(i), dgbdr(i), &
!            dg2dr(i), dg21dr(i), dg22dr(i), dbdr(i), cvdrift0(i), dIdrho, dgradparbdrho(i), gradparb(i), dqdr, dpsidrho**2
!    end do

  end subroutine leqcoefs

  function invR (r, theta, i)
   
    integer, intent (in) :: i
    real, intent (in) :: r, theta
    real :: invR
    
    invR=1./Rpos(r, theta, i)
    
  end function invR

  function rcenter(rp)

    real, intent(in) :: rp
    real :: rcenter

    rcenter = surf%R_center
    
  end function rcenter

  function Rpos (r, theta, j)
   
    use constants, only: pi

    integer, intent (in) :: j
    real, intent (in) :: r, theta
    real :: Rpos
    real :: g, gp, dr
    integer :: i
    
    dr = r - surf%r

! For Y Xiao: 
!    g = surf%delp/surf%r + surf%d * sin(theta)**2
!    Rpos = surf%R_center*(1.+r*(cos(theta)-g)-g*dr)

    g = cos(theta + surf%d * sin(theta-0.5*pi*surf%a))
    gp = -sin(theta + surf%d * sin(theta-0.5*pi*surf%a)) &
         *(surf%dp*sin(theta-0.5*surf%a)-surf%d*0.5*pi*surf%ap*cos(theta-0.5*pi*surf%a))

    ! allow for strange specification of R_psi
    if (j==ntg+1) then
       i = -ntg
    else
       i = j
    end if
    
    ! second line here is (1/2)*(r-r0)**2*d2R/dr|_r0
    ! note that d2R=0 unless read_profile_variation = T in input file
    Rpos = surf%R_center + surf%delp*dr + g*surf%r + (g+surf%r*gp)*dr &
         + 0.5*(r-0.5)**2*d2R(i)
    
  end function Rpos

  function Zpos (r, theta, j)
   
    integer, intent (in) :: j
    real, intent (in) :: r, theta
    real :: Zpos, dr
    integer :: i

    ! allow for strange specification of Z_psi
    if (j==ntg+1) then
       i = -ntg
    else
       i = j
    end if

    dr = r - surf%r
    ! note that d2Z=0 unless read_profile_variation=T in input file
    Zpos = surf%k*sin(theta)*surf%r + (surf%r*surf%kp + surf%k)*sin(theta)*dr &
         + 0.5*(r-0.5)**2*d2Z(i)
    
  end function Zpos

  function psi (r, theta)
   
    real, intent (in) :: r, theta
    real :: psi

    psi = r - surf%r
    
  end function psi

  function mod2pi (theta)
    
    real, intent(in) :: theta
    real :: pi, th, mod2pi
    real, parameter :: theta_tol = 1.e-6
    logical :: out
    
    pi=2.*acos(0.)
    
    if(theta <= pi .and. theta >= -pi) then
       mod2pi = theta
       return
    endif
    
    if(theta - theta_tol <= pi .and. theta >= -pi) then
       mod2pi = pi
       return
    endif

    if(theta <= pi .and. theta + theta_tol >= -pi) then
       mod2pi = -pi
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
   
  function diameter (rp)
  
    real :: rp, diameter

    diameter = 2.*rp

  end function diameter

  function dbtori (pbar)
    real :: pbar, dbtori
    dbtori = 1.
  end function dbtori

  function btori (pbar)
    real :: pbar, btori
    btori = surf%r_geo
  end function btori

  function qfun (pbar)
    real :: pbar, qfun
    qfun = surf%q
  end function qfun

  function pfun (pbar)
    real :: pbar, pfun
    pfun = 0.5*beta_0
  end function pfun
  
  function dpfun (pbar)  
    real :: pbar, dpfun    

       dpfun = -1.

  end function dpfun

  function dpdrhofun(rho)

    real :: rho, dpdrhofun

    dpdrhofun = surf%pp

  end function dpdrhofun
  
  function betafun (pbar)  
    real :: pbar, betafun
    betafun = beta_0
  end function betafun

end module leq
