program eiktest

  use geometry
  implicit none

  integer :: i, nth, ntgrid, ntheta
  real :: pi, broot
  real :: q, qq, qplus, qmnus
  real :: bb, bm, bp, profile_fac
  real :: rk, rkm, rkp, rkpri
  real :: rt, rtm, rtp, rtpri
  real :: rd, rdm, rdp, rdpri
  real :: rhocm, rhocp
  real :: bbar

  real :: diffscheme, beta_p1, beta_p2, beta_prime_times, beta_prime_over
  integer :: nbeta

  real :: airat, wgt, wdavg, jac
  logical :: fast, dipole

  namelist/stuff/ntheta,nperiod,rmaj,akappri,akappa,shift,equal_arc, &
       rhoc,rmin,rmax,itor,qinp,iflux,delrho,tri,bishop, &
       irho,isym,tripri,efit_eq,dfit_eq,writelots,R_geo, &
       gen_eq, ppl_eq, eqfile, local_eq, idfit_eq,&
       chs_eq,&
       s_hat_input,p_prime_input,invLp_input,beta_prime_input, &
       diffscheme,nbeta,beta_p1,beta_p2,alpha_input,big, &
       beta_prime_times, beta_prime_over, fast, profile_fac, &
       tstar, shotnum, gs2d_eq, transp_eq, Xanthopoulos

  pi=2.*acos(0.)
     
! set defaults

  gen_eq = .false.
  ppl_eq = .false.
  chs_eq = .false.
  idfit_eq = .false.
  dfit_eq = .false.
  efit_eq = .false.
  gs2d_eq = .false.
  transp_eq = .false.
  Xanthopoulos = .false.

  dipole = .false.
  shotnum=1001109020
  tstar = 1.0

  ntheta=32; nperiod=1
  equal_arc = .false. ; bishop = 0 ; in_nt = .false. ; fast = .true.
  itor=1; iflux=0;  big = 8

  eqfile='dskeq.cdf'
  rhoc=0.5

  rmin=0.05
  rmax=1.0  ! other than unity is an inconsistent normalization

  beta_prime_times = -1.

  profile_fac = 0.5

  rmaj = 3   ;     akappa = 1   ;     tri = 0   ;     qinp = 2
  R_geo = 3  ;     akappri = 0  ;     tripri=0  ;     shat = 1
  shift = 0  
   
  s_hat_input = 1
  beta_prime_input = -1
  p_prime_input = -2
  invLp_input = 5

  irho = 2 ;      isym = 0

  delrho = 0.01

  eqinit = 1       ! Mike K. codes do not have this variable.

      
!     read in data
  open(10,file='eik.in')
  read(10,stuff)
  close(10)
  
  if(bishop >= 7) then
     if(beta_prime_times == -1) then
        write(*,*) 'beta_prime_times should be set for bishop = 7 or 8'
        beta_prime_times = 1.
     endif
     dp_mult = beta_prime_times
  endif

!
! Note that if iflux=1 then always choose itor=1
!
  if(iflux.eq.1 .and. itor.ne.1) then
     itor=1
     write(*,*) 'Forcing itor=1'
  endif
  
  if(iflux.eq.2) then
     write(*,*) 'iflux=2 is not a standalone option.'
     stop
  endif
!
! Catch if no q profile 
!
  if(iflux.ne.1 .and. irho.eq.1) then
     irho=2
     write(*,*) 'Forcing irho=2'
  endif
  
  if(iflux.ne.1) then
     if(gen_eq) write(*,*) 'Forcing gen_eq to be false'
     if(ppl_eq) write(*,*) 'Forcing ppl_eq to be false'
     if(chs_eq) write(*,*) 'Forcing chs_eq to be false'
     if(transp_eq) write(*,*) 'Forcing transp_eq to be false'
     if(gen_eq .or. ppl_eq .or. transp_eq) write(*,*) 'because iflux.ne.1'
     gen_eq=.false.
     ppl_eq=.false.
     chs_eq = .false.
  endif
  
  open(unit=21,file='eik.out',status='unknown')
  open(unit=24,file='shear.out',status='unknown')
  open(unit=25,file='trap.out',status='unknown')
  open(unit=11,file='eik2.out',status='unknown')
  open(unit=99,file='eik3.out',status='unknown')
  open(unit=78,file='eik7.out',status='unknown')
!     
!     compute the theta grid

  if((.not. gen_eq) &
       .and. (.not. transp_eq) .and. (.not. ppl_eq) .and. (.not. chs_eq)) &
       call init_theta(ntheta)
  
  call eikcoefs
! note that ntheta may have changed

  ntgrid = (size(theta)-1)/2
  ntheta = (size(theta)-1)/(2*nperiod-1)

  if (writelots) then
     bbar = 0.
     do i=-ntgrid,ntgrid-1
        bbar = bbar + bmag(i)*(theta(i+1)-theta(i))
     end do
     bbar = bbar/(2.*pi)

     write(*,*) 'B_bar = ',bbar
  end if
           
  if (dipole) then
     shat = max(1.e-7, shat)
     if (shat == 1.e-7) gds22 = shat*shat
  end if

  write(21,*) ' ntgrid  nperiod  ntheta, drhodpsi, rmaj, shat, kxfac, q'
!CMR 14/01/05 change format to prevent line breaks with intel compiler and allow runngridgen to read the output files
  write(21,fmt='(3i8,1p5e12.4)') ntgrid, nperiod, ntheta, drhodpsin, rmaj, shat, kxfac, qsf
!CMR
  
  write(21,*) '    gbdrift     gradpar       grho    tgrid'
  do i= -ntgrid,ntgrid
     write(21,1000) gbdrift(i), gradpar(i),grho(i),theta(i)
  enddo
  
  write(21,*) '    cvdrift            gds2            bmag   tgrid'
  do i= -ntgrid,ntgrid
     write(21,1000) cvdrift(i), gds2(i),bmag(i),theta(i)
  enddo
     
!      write(*,*) 'bmag(0)= ',bmag(0),' bmag(pi)= ',bmag(ntgrid)

  write(21,*) '    gds21             gds22  tgrid'
  do i= -ntgrid,ntgrid
     write(21,1000) gds21(i), gds22(i),theta(i)
  enddo
  
  write(21,*) '    cvdrift0          gbdrift0   tgrid'
  do i= -ntgrid,ntgrid
     write(21,1000) cvdrift0(i),gbdrift0(i),theta(i)
  enddo
     
  call Hahm_Burrell(irho, profile_fac)

  broot=-1.
!      write(*,*) 'rho ',rhoc
!      write(*,*) 'rp ',rpofrho(rhoc)
!      write(*,*) 'r ',rfun(rpofrho(rhoc),0.,broot)
  dipole = dfit_eq .or. idfit_eq

  if (.not. dipole) then
     bb=beta_a_fun(rfun(rpofrho(rhoc),0.,broot),0.)
     write(21,*) '-8*pi/B_0**2 * dp/drho= '
     if(bb.ne.0.) then
        write(21,1000) -dbetadrho
        write(21,1000) -dbetadrho/bb, -0.5*dbetadrho/bb
        write(21,1000) bb, bb/2.
     endif
  end if

  write(11,*) 'dV/drhon= ',dvdrhon
  if (.not. dipole) then
     write (11,*)
     write (11,*) 'The following q and shat values'
     write (11,*) 'were calculated directly'
     write (11,*) 'from the equilibrium.'
     write (11,*) 'They should match the ones above' !EGH
     q=qfun(0.)
     write(11,*) 'q_0= ',q
     q=qfun(pbarofrho(rhoc))
     write(11,*) 'rho= ',rhoc,' q(rho)= ',q
  !      write(*,*) 'rho= ',rhoc,' q(rho)= ',q
     q=qfun(0.95)
     write(11,*) 'q_95= ',q
  
     call geofax(rhoc, rt, rk, rd)
     if(iflux == 1) then
        rhocp=rhoc+delrho
        rhocm=rhoc-delrho
        call geofax(rhocp, rtp, rkp, rdp)
        call geofax(rhocm, rtm, rkm, rdm)
        qplus=qfun(pbarofrho(rhocp))
        qmnus=qfun(pbarofrho(rhocm))
        rkpri=0.5*(rkp-rkm)/delrho
        rtpri=0.5*(rtp-rtm)/delrho
        rdpri=0.5*(rdp-rdm)/delrho
        if (dipole) then
           shat = 0.
        else
           shat = 0.5*rhoc*(qplus-qmnus)/qfun(pbarofrho(rhoc))/delrho
        endif
     else
        rkpri=akappri
        rtpri=tripri
        rdpri=shift
     endif
  
     qq=qfun(pbarofrho(rhoc))
     write(11,*) 'shat= ',shat
     write(11,*) 'qinp= ',qq
     write(11,*) 'shift= ',rdpri
     write(11,*) 'akappa= ',rk
     write(11,*) 'akappri= ',rkpri
     write(11,*) 'tri= ',rt
     write(11,*) 'tripri= ',rtpri
     write(11,*) 'surface area= ',surfarea,' avgrmid**2'
  end if

!      write(*,*) 'mag(Bp) at rho= ',rhoc
!      do i=-nth,nth
!         r=rfun(rpofrho(rhoc),theta(i))
!         write(*,*) r,theta(i),bpmagfun(r,theta(i))
!      enddo

!      write(*,*) '    theta        gbdrift          gradpar        grho'
!      do i= -5,5
!      do i= -ntgrid,ntgrid
!         write(*,1000) theta(i),gbdrift(i), gradpar(i),grho(i)
!      enddo

     
  close(11)
  close(21)
  close(99)

  open (unit=21,file='eik5.out',status='unknown')
  write (21,'(''#theta1 gradpar2 bmag3 grho4 '', &
       &   ''gbdrift5 gbdrift06 cvdrift7 cvdrift08 '', &
       &    ''gds29 gds2110 gds2211'')')
  if (iflux .ne. 1) then
     do i=-ntgrid,ntgrid
        write (21,1000) theta(i),gradpar(i),bmag(i),grho(i), &
             gbdrift(i),gbdrift0(i),cvdrift(i),cvdrift0(i), &
             gds2(i),gds21(i),gds22(i)
     enddo
  else
     do i=-ntgrid,ntgrid
        write (21,1000) theta(i),gradpar(i),bmag(i),grho(i), &
             gbdrift(i),gbdrift0(i),cvdrift(i),cvdrift0(i), &
             gds2(i),gds21(i),gds22(i)
     enddo
  endif
  close (unit=21)

  open (unit=21,file='eik6.out',status='unknown')

  write (21,fmt="('# shape: torus')")
  write (21,fmt="('# q = ',e10.4,' drhodpsi = ',e10.4)") qsf, drhodpsin
  write (21,fmt="('# theta1             R2                  Z3               alpha4      ', &
            &   '       Rprime5              Zprime6           alpha_prime7 ')")
  do i=-ntgrid,ntgrid
     write (21,1000) theta(i),Rplot(i),Zplot(i),aplot(i), &
          Rprime(i),Zprime(i),aprime(i)
  enddo
  close (unit=21)

!  wdavg= 0.
!  wgt = 0.
!  do i=-ntgrid,ntgrid-1
!     jac = abs(1./(drhodpsin*gradpar(i)*bmag(i)))
!     wgt = wgt + (theta(i+1)-theta(i))*jac
!     wdavg = wdavg+(theta(i+1)-theta(i))*jac!*0.5*(cvdrift(i)+cvdrift(i+1))
!  end do
!  wdavg=wdavg!/wgt
!
!  print *, 'Avg curvature = ',wdavg!, wgt
  stop
1000  format(20(1x,1pg18.11))
end program eiktest


