module geometry

  implicit none

  public

  private :: bgrad, root, fluxavg

  real, allocatable, dimension(:) :: grho, theta, bmag, gradpar, &
       cvdrift, cvdrift0, gbdrift, gbdrift0, gds2, gds21, gds22, jacob, &
       Rplot, Zplot, Rprime, Zprime, aplot, aprime, Uk1, Uk2, &
       cdrift, cdrift0, gds23, gds24, gds24_noq, gbdrift_th, cvdrift_th, Bpol, &
       dgradpardrho, dgradparbdrho, dcvdrift0drho, dgbdrift0drho, &
       dcvdriftdrho, dgbdriftdrho, dgds2dr, dgds21dr, dgds22dr, dBdrho
  
  real :: dkxfacdr

  real, allocatable, dimension(:) :: J_X, B_X, g11_X, g12_X, g22_X, &
       K1_X, K2_X, gradpar_X

  real :: Rref_X, Bref_X

  real :: rhoc, rmaj, r_geo, shift, dbetadrho, kxfac
  real :: qinp, shat, akappa, akappri, tri, tripri, dpressdrho, asym, asympri, d2qdr2
  real :: betadbprim, d2psidr2
  real :: delrho, rmin, rmax, qsf, aminor
  
  real :: s_hat_input, p_prime_input, invLp_input, beta_prime_input, &
       alpha_input, dp_mult

  integer :: nperiod
  integer :: itor, iflux, irho
  integer :: shotnum
  real :: tstar
   
  real :: dp_new, di_new
  real, private :: psi_0, psi_a
  real :: B_T0, avgrmid, dvdrhon, surfarea, grho1n, grho2n, drhodpsin
  real :: grhoavg  ! needed for Trinity -- MAB

  real, dimension(3) :: rpval
  real :: rpmin, rpmax
  real :: s_hat_new
  real :: beta_prime_new

  integer :: isym, big
  integer :: eqinit = 1
  
  logical :: gen_eq, efit_eq, ppl_eq, local_eq
  logical :: write_profile_variation, read_profile_variation

  logical :: in_nt, writelots, equal_arc, idfit_eq, gs2d_eq
  logical :: transp_eq, Xanthopoulos

  integer :: bishop

  public :: beta_a_fun, eikcoefs, geofax, iofrho, pbarofrho, &
       qfun, rpofrho, rcenter, init_theta, nth_get, f_trap

  integer, private :: ntgrid, nth, ntheta
  logical :: debug = .false.

  character(800) :: eqfile

  type advanced_parameters_type  
    logical :: equal_arc
    real :: dp_mult
    real :: delrho
    real :: rmin
    real :: rmax
    integer :: isym
    logical :: in_nt
    logical :: writelots
    integer :: itor
  end type advanced_parameters_type

  type miller_parameters_type
     real :: rmaj
     real :: R_geo
     real :: akappa
     real :: akappri
     real :: tri
     real :: tripri
     real :: shift
     real :: qinp
     real :: shat
     real :: asym
     real :: asympri
     real :: d2qdr2
     real :: betadbprim
     real :: d2psidr2
   end type miller_parameters_type

  type(advanced_parameters_type) :: advanced_parameters
  type(miller_parameters_type) :: miller_parameters

  !> These are functions of theta that are output by the 
  !! geometry module
  type coefficients_type

    real :: grho   
    real :: bmag       
    real :: gradpar    
    real :: cvdrift    
    real :: cvdrift0   
    real :: gbdrift    
    real :: gbdrift0   
    real :: cdrift    
    real :: cdrift0    
    real :: gbdrift_th 
    real :: cvdrift_th 
    real :: gds2       
    real :: gds21      
    real :: gds22      
    real :: gds23      
    real :: gds24      
    real :: gds24_noq  
    real :: jacob      
    real :: Rplot      
    real :: Zplot      
    real :: aplot      
    real :: Rprime     
    real :: Zprime     
    real :: aprime     
    real :: Uk1        
    real :: Uk2        
    real :: Bpol       

  end type coefficients_type

  !> These are coefficients that are constant
  !! on a flux surface
  type constant_coefficients_type
    real :: qsf
    real :: rmaj
    real :: shat
    real :: kxfac
    real :: aminor
  end type constant_coefficients_type

  type(coefficients_type) :: output_coefficients

  type(constant_coefficients_type) :: constant_coefficients

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  
!!!  Input: 
!!!     rhoc :: minor radius of surface of interest (see irho below)
!!!     itor :: default one.  Choose zero to get "s-alpha" equilibria (see iflux=0)
!!!     iflux :: chooses mode of operation.  
!!!              case (0)  :: no equilibrium.  Miller parameterization of local equilibrium
!!!                           iflux = 0 forces local_eq = .true. and irho = 2
!!!              case (1)  :: use numerical equilibrium.  5 choices:
!!!                 ppl_eq = .true. :: use PPL style NetCDF equilibrium (Menard)
!!!                 transp_eq = .true. :: use TRXPL style NetCDF equilibrium (Menard)
!!!                 gen_eq = .true. :: use GA style NetCDF equilibrium (TOQ)
!!!     -------------------> [note: if either of the above == .true., 
!!!                          set eqfile = input file name]
!!!                 efit_eq = .true. :: use EFIT  equilibrium
!!!              case (2) :: running inside nt without actual equilibria
!!!                           iflux = 2 forces local_eq = .true. and irho = 2
!!!             case (10) :: used by nt to initialize some variables.  Obsolescent.
!!!     irho :: choose flux surface label
!!!             case (1) :: sqrt(toroidal flux)/sqrt(toroidal flux of LCFS)
!!!             case (2) :: diameter/(diameter of LCFS).  recommended
!!!             case (3) :: poloidal flux/(poloidal flux of LCFS)
!!!     
!!!     equal_arc :: .true. makes grad_par coefficient a constant.
!!!
!!!     !! New feature:  Given an actual equilibrium, allow s_hat and d pressure/d rho
!!!     !! to vary self-consistently (a la Greene and Chance, Bishop).  We use the 
!!!     !! Bishop formalism [C. M. Bishop, et al., NF, Vol. 24, 1579 (1984).]
!!!
!!!     bishop :: uses Bishop relations to generate coefficients.  
!!!            case(0) :: do not use Bishop relations
!!!            case(1) :: use Bishop relations with actual equilibrium shat, p'
!!!            case(3) :: use Bishop relations with new shat, L_p from below
!!!     s_hat_input :: new magnetic shear, rho/q*dq/drho
!!!     invLp_input :: new -1/pressure * d pressure/d rho
!!!            case(4) :: use Bishop relations with new shat, beta' from below
!!!     s_hat_input :: new magnetic shear, rho/q*dq/drho
!!!     beta_prime_input :: new d beta/d rho
!!!            case(5) :: use Bishop relations with new shat, alpha from below
!!!     s_hat_input :: new magnetic shear, rho/q*dq/drho
!!!     alpha_input :: new alpha = -q**2 * rmaj * d pressure/d rho
!!!            case(6) :: use Bishop relations with equilibrium shat, beta' from below
!!!     beta_prime_input :: new d beta/d rho
!!!            case(7) :: use Bishop relations with equilibrium shat, beta' from below
!!!     dp_mult :: d pressure/d rho = equilibrium gradient * dp_mult
!!!             (In each case, the definition of rho is determined by irho above)
!!!
!!!     nperiod :: number of cells of width 2 pi in theta. 
!!!  
!!!     rmaj :: major radius of magnetic axis (in units of aminor) if iflux = 1
!!!          or major radius of center of flux surface of interest (units of aminor) otherwise
!!!     r_geo :: major radius of last closed flux surface
!!!     shift :: d rmaj/drho.  (Typically < 0.)
!!!     dIdrho :: d I/drho where I = R B_T
!!!     qinp :: safety factor q at surface of interest
!!!     shat :: d q/drho at surface of interest
!!!     akappa :: elongation (k) of surface of interest
!!!     akappri :: dk/drho
!!!     tri :: triangularity of surface of interest
!!!     tripri :: d tri/drho
!!!     dpressdrho :: dpressure /drho = 0.5 dbeta /drho at surface of interest.
!!!  
!!!     delrho :: numerical parameter.  Should be "small enough".  Typically 0.001 ok.
!!!     rmin :: minimum minor radius at which items are evaluated
!!!     rmax :: should be equal to unity 
!!!     
!!!     isym :: 1 -> assume up-down symmetric equilibrium.
!!!     in_nt :: .true. if running inside of nt.
!!!     writelots :: .true. for more output
!!!  
!!!  Output: 
!!!     theta :: the theta grid of the results.  The final grid has 
!!!              pi - 2*pi*nperiod <= theta <= 2*pi*nperiod - pi
!!!     bmag :: |B| 
!!!     gradpar :: coefficient of b_hat.grad operator 
!!!     grho :: grad(rho)
!!!     cvdrift  :: coefficient of v_par**2/Omega b_hat x (b_hat.grad b_hat).grad 
!!!     cvdrift0 :: part of total curvature drift operator proportional to theta_0
!!!     gbdrift  :: coefficient of mu/Omega b_hat x (grad B).grad operator
!!!     gbdrift0 :: part of total grad B drift operator proportional to theta_0
!!!     gbdrift_th :: grad B drift . grad theta
!!!     cvdrift_th :: curvature drift . grad theta
!!!     cdrift :: part of coriolis drift operator independent of theta_0
!!!     cdrift0 :: part of coriolis drift operator proportional to theta_0
!!!     gds2 :: part of grad_perp**2 operator independent of theta_0
!!!     gds21 :: part of grad_perp**2 operator proportional to theta_0
!!!     gds22 :: part of grad_perp**2 operator proportional to theta_0**2
!!!     gds23 :: part of v_E . grad theta operator independent of theta_0
!!!     gds24 :: part of v_E . grad theta operator proportional to theta_0
!!!     gds24_noq :: gds24/dqdrp
!!!     jacob :: Jacobian
!!!     
!!!  Use: 
!!!     Include this module, set input variables appropriately and call eikcoefs.  
!!!
!!!     NOTE: If gen_eq and ppl_eq are false, you should call init_theta
!!!           before eikcoefs to set the number of theta gridpoints per 2 pi.
!!!
!!!     I do not recommend calling other routines, but a few are available for 
!!!     diagnostic purposes.  
!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
contains

  subroutine eikcoefs (ntheta_returned)
    
    use  geq, only: eqin, geq_init
    
    use  peq, only: peqin => eqin, teqin, peq_init
    use  eeq, only: efitin, efit_init, gs2din
    use ideq, only: idfitin, ideq_init => dfit_init
    use  leq, only: leqin, dpdrhofun, leqcoefs

    implicit none

    integer, optional, intent (out) :: ntheta_returned

! Local variables:

    real, allocatable, dimension(:)   :: smodel, dsdrp, grho2, seik, seik1, seik2
    real, allocatable, dimension(:)   :: dsdthet, dsdthet1, dsdthet2, grho1
    real, allocatable, dimension(:)   :: gdsdum1, gdsdum2, gdsdum3
    real, allocatable, dimension(:)   :: gdsdum4, gdsdum5  ! MAB
    real, allocatable, dimension(:)   :: th_bish, Rpol, ltheta
    real, allocatable, dimension(:)   :: rgrid, rgrid1, rgrid2, Bpolmag, Bmod, dSdl
    real, allocatable, dimension(:)   :: rmajor, ans, ds, arcl
    real, allocatable, dimension(:)   :: Zoftheta, dZdl
    real, allocatable, dimension(:,:) :: thgrad, rpgrad, crpgrad, grads
    real, allocatable, dimension(:,:) :: bvector, gradrptot, gradstot, gradztot
    real, allocatable, dimension(:,:) :: gradthtot ! MAB

    real, allocatable, dimension(:) :: a_bish, b_bish, c_bish
    real, allocatable, dimension(:) :: da_bish, db_bish, dc_bish, trip

    real :: rho, rp, rp1, rp2, drhodrp, dqdrho, drhodrhod
    real :: dpsidrp, dpsidrho, dqdrp, dsdrptot, drhodpsi
    real :: qval, qval1, qval2, delrp, dum
    real :: pi, bi, Ltot, rinv

    real :: a_b, b_b, c_b
    real :: s_hat, dp, di, pressure, tmp, rbar

    character(1) :: char
    integer :: i, j, k, itot, nthg, n

!     compute the initial constants
    pi=2.*acos(0.)
    
if (debug) write(6,*) "eikcoefs: local_eq=",local_eq
    if (local_eq .and. iflux == 1) then
       write (*,*) 'Forcing iflux = 0'
       iflux = 0
    end if

    if(iflux == 0 .or. iflux == 2) then
       if(.not. local_eq) write(*,*) 'Forcing local_eq = true'
       local_eq = .true.
       if(irho /= 2) write(*,*) 'Forcing irho = 2'
       irho = 2
       if(bishop == 0) then 
          write(*,*) 'Forcing bishop = 1'
          bishop = 1          
       endif
    endif

if (debug) write(6,*) "eikcoefs: call check"
    call check(gen_eq, efit_eq, ppl_eq, &
         local_eq, idfit_eq) 

    if(.not. gen_eq .and. .not. ppl_eq &
    .and. .not. transp_eq &
         .and. .not. allocated(theta)) then
       write(*,*) 'You should call init_theta to specify the '
       write(*,*) 'number of theta grid points per 2 pi '
       write(*,*) 'before calling eikcoefs.'
       write(*,*) 
       write(*,*) 'Proceeding with ntheta = 32.'
       ntheta = 32
       call init_theta(ntheta)
    endif

    if(iflux /= 1 .and. iflux /= 10) then

       call leqin(rmaj, R_geo, akappa, akappri, tri, tripri, rhoc, delrho, shift, &
            qinp, s_hat_input, asym, asympri, beta_prime_input, betadbprim, d2qdr2, &
            d2psidr2, ntheta/2+1, nperiod, write_profile_variation, read_profile_variation)
       if(.not.allocated(gds22)) call alloc_module_arrays(ntgrid)
       call alloc_local_arrays(ntgrid)
       ! get geometry coefficients needed for dgk/drho equation
       call leqcoefs (dgradpardrho, dgradparbdrho, dcvdrift0drho, dgbdrift0drho, &
            dcvdriftdrho, dgbdriftdrho, dgds2dr, dgds21dr, dgds22dr, dBdrho, &
            gds2, gds21, gds22, gradpar, gbdrift0, cvdrift0, gbdrift, cvdrift, &
            drhodpsin, grho, bmag, Bpol, dkxfacdr)
    endif

if (debug) write(6,*) "eikcoefs: iflux=",iflux
    select case (iflux)
       case (0)
          avgrmid=1.
       case (1, 10)
          if (gs2d_eq) efit_eq = .true.
          if(gen_eq) then
             call eqin(eqfile, psi_0, psi_a, rmaj, B_T0, avgrmid, eqinit, in_nt, nthg)
             if(present(ntheta_returned)) then
                call tdef(nthg, ntheta_returned)
             else
                call tdef(nthg)
             endif
!CMR+SSAAR: moved following if clause ahead of ppl_eq to avoid ball problems 
          else if(transp_eq) then
if (debug) write(6,fmt='("eikcoefs: transp_eq, eqfile=",a)') eqfile
             ppl_eq = .true.
             call teqin(eqfile, psi_0, psi_a, rmaj, B_T0, avgrmid, eqinit, in_nt, nthg)
if (debug) write(6,*) 'eikcoefs: transp_eq, called teqin'
if (debug) write(6,fmt='("eikcoefs: teqin returns",1p5e10.2,i6,l2,i6)') psi_0, psi_a, rmaj, B_T0, avgrmid, eqinit, in_nt, nthg
             if(present(ntheta_returned)) then
                call tdef(nthg, ntheta_returned)
             else
                call tdef(nthg)
             endif
if (debug) write(6,*) 'eikcoefs: transp_eq, called tdef'
!CMRend
          else if(ppl_eq) then
             call peqin(eqfile, psi_0, psi_a, rmaj, B_T0, avgrmid, eqinit, in_nt, nthg)
             if(present(ntheta_returned)) then
                call tdef(nthg,ntheta_returned)
             else
                call tdef(nthg)
             endif
          else if(efit_eq) then
             if(big <= 0) big = 8

             if(gs2d_eq) then
if (debug) write(6,*) "eikcoefs: call gs2din eqfile=",eqfile
                call gs2din(eqfile, psi_0, psi_a, rmaj, B_T0, avgrmid, eqinit, big) 
if (debug) write(6,*) "eikcoefs: done gs2din  psi_0,psi_a, rmaj, B_T0, avgrmid=",psi_0,psi_a, rmaj, B_T0, avgrmid
             else
                call efitin(eqfile, psi_0, psi_a, rmaj, B_T0, avgrmid, eqinit, big) 
             endif
          else if(idfit_eq) then
             call idfitin(eqfile, theta, psi_0, psi_a, rmaj, B_T0, avgrmid, eqinit) 
             if(present(ntheta_returned)) then
                call tdef(nthg, ntheta_returned)
             else
                call tdef(nthg)
             endif
          endif
if (debug) write(6,*) 'eikcoefs: iflux=',iflux
          if(iflux == 10) return
if (debug) write(6,*) 'eikcoefs: alloc_module_arrays, ntgrid=',ntgrid
          if(.not.allocated(gds22)) call alloc_module_arrays(ntgrid)
if (debug) write(6,*) 'eikcoefs: alloc_local_arrays'
          call alloc_local_arrays(ntgrid)
       case (2) 
          continue
    end select
if (debug) writelots=.true.
if (debug) write(6,*) "eikcoefs: if (writelots)"
    if(writelots) then
       write(11,*) 'All lengths normalized to the avg midplane minor radius'
       if(irho /= 2) then
          write(11,*) 'except rho, which always goes from 0 to 1.'
          write(11,*) 'The rho normalization is determined by irho.'
       endif
       write(11,*) 'avgrmid = ',avgrmid,' meters'
       if(iflux == 1) write(11,*) 'R_mag = ',rmaj
       if(iflux == 1) write(11,*) 'B_T0 = ',B_T0
    endif

    rho=rhoc

if (debug) write(6,*) "eikcoefs: find rp"
    if(iflux == 1) then
       rpmin = psi_0
       rpmax = psi_a
if (debug) write(6,*) "eikcoefs: rpmin,rpmax=",rpmin,rpmax
       rp = rpofrho(rho)
if (debug) write(6,*) "eikcoefs: rp=",rp
    else
       rpmin = 0.
       rpmax = 1.0
       rp=rho
    endif

    if(iflux == 1) R_geo=rcenter(rpmax)

if (debug) write(6,*) "eikcoefs: find rgrid"    
    if(efit_eq) then
       call rtg(rgrid, rp)
    else
       rgrid = rp
    endif

    if(writelots) then
       write(11,*) 'rp of rho = ',rp
       write(11,*) 'Rcenter(rho) = ',rcenter(rp)
       write(11,*) 'R_geo = ',R_geo,' !or ',btori(rgrid(0), 0.)
    endif

    if(writelots) then
       rbar = 0.
       write(78,*) '# theta_pol,  R,   Z'
       do i= -nth,nth         
          write(99,*) 'i=',i,' thet= ',theta(i),' rgrid= ',rgrid(i)
          write(78,*) theta(i),Rpos(rgrid(i),theta(i),i),Zpos(rgrid(i),theta(i),i)
       enddo

       do i=-nth,nth-1
          rbar = rbar + rgrid(i) * (theta(i+1)-theta(i))
       end do
       rbar = rbar / (2.*pi)
       write (*,*) 'r_bar = ',rbar
       write (*,*)
    endif

    call drho_drp(rp, drhodrp)

    if(nperiod > 1) call periodic_copy(theta, 2.*pi)

    if(writelots) then
       write(11,*) 'drhodrp=',drhodrp
       if (irho == 1) then
          call drho_drhod(rp, drhodrp, drhodrhod)
          write(11,*) 'drho_drhod=',drhodrhod
       end if
    end if

!
! should test whether this is a new equilibrium
!

 if (debug) write(11,*) 'eikcoefs: various logs', gen_eq, ppl_eq, efit_eq, idfit_eq

    if(gen_eq)    call geq_init
    if(ppl_eq)    call peq_init
if (debug) write(6,*) "eikcoefs: call eeq_init"
!    if(efit_eq)  call eeq_init
    if(efit_eq)   call efit_init
if (debug) write(6,*) "eikcoefs: done eeq_init"
    if(idfit_eq)  call ideq_init

if (debug) write(6,*) "eikcoefs: call rmajortgrid"
    call rmajortgrid(rgrid, theta, rmajor)

!     compute gradient of rp
    char='P'

    ! note that rpgrad is given in (R,Z,phi) coordinates if bishop=0
    ! and (rho,l,phi) coordinates if bishop>0
    ! also, rpgrad is grad psi if iflux=1 and grad rho if iflux=0 -- MAB
    if(bishop == 0) then
       call  grad(rgrid, theta, rpgrad, char, dum, nth, ntgrid)
       crpgrad = rpgrad
    else
       call bgrad(rgrid, theta, rpgrad, char, dum, nth, ntgrid)
    endif

    if(bishop /= 0) then       
       call loftheta(rgrid, theta, ltheta)
       call grad(rgrid, theta, crpgrad, char, dum, nth, ntgrid)
       call th_bishop(crpgrad, th_bish, nth)
       ! note that Bpolmag and Bmod are recalculated later if using
       ! Miller (local) equilibrium
       call B_pol(rgrid, theta, crpgrad, Bpolmag, nth)
       call B_mod(rgrid, theta, Bpolmag, Bmod, nth)
       call R_pol(theta, th_bish, ltheta, Rpol, nth)

       if(nperiod > 1) then
          call periodic_copy(rmajor, 0.)
          call periodic_copy(Bpolmag, 0.)
          call periodic_copy(th_bish, 0.)
          call periodic_copy(Rpol, 0.)
          Ltot = ltheta(nth)-ltheta(-nth)
          call periodic_copy(ltheta, Ltot)
       endif
    endif 

!     compute |grad rho|:
    if (.not. local_eq) then
       do k=-nperiod+1,nperiod-1
          do i=-nth,nth
             itot=i+k*ntheta
             grho(itot)=sqrt(rpgrad(i,1)**2+rpgrad(i,2)**2)*abs(drhodrp)
          enddo
       enddo

       if(isym == 1) call sym(grho, 0, ntgrid)
    end if

!     compute  coordinate gradients

    if (debug) write(6,*) "eikcoefs: call thetagrad"

    call thetagrad(rgrid,thgrad)

    if (debug) write(6,*) "eikcoefs: done thetagrad"
     
!     compute eikonal S
    if(iflux == 0 .or. iflux == 2) qval=qfun(pbarofrho(rho))
    if (debug) write(6,*) "eikcoefs: call eikonal"
    call eikonal(rgrid, rpgrad, thgrad, qval, seik, dsdthet, dpsidrp)
    if (debug) write(6,*) "eikcoefs: done eikonal"
    if(writelots) write(11,*) 'q= ',qval
	
    if(iflux == 0 .or. iflux == 2) then
       bi = btori(rgrid(0), theta(0))
       if (.not. local_eq) drhodpsin = drhodrp/dpsidrp
       do i=-nth,nth
          rinv = invRfun(rgrid(i), theta(i), i)
          Bpolmag(i) = sqrt(rpgrad(i,1)**2+rpgrad(i,2)**2)*rinv*abs(dpsidrp)
          Bmod(i) = sqrt( (bi*rinv)**2 + bpolmag(i)**2) 
       enddo       
       if(nperiod > 1) then
          call periodic_copy(Bpolmag, 0.)
          call periodic_copy(Bmod, 0.)
       endif
    endif

!    do i=-nth,nth
!       write(43,*) theta(i), Bpolmag(i), rpos(rgrid(i), theta(i)), zpos(rgrid(i), theta(i)), seik(i)
!    enddo

!
! ******************
! Debug: 
! Set itor=1 and the first two columns should be the 
! same (better agreement results from higher ntheta)
! Set itor=0 and the first and third columns should 
! be equal.  All of this, of course, with iflux=0.
! 
! This amounts to checking the expression for
! nu for axisymmetric, low beta equilibrium:
! nu = - q r/R sin(theta).
!
!      call rmajortgrid(rgrid, theta, rmajor)
!      do i=-nth,nth
!         smodel(i)=-qval*(theta(i)-rgrid(i)*sin(theta(i))/rmajor(i))
!         write(*,*) seik(i),smodel(i),-qval*theta(i)
!      enddo
! ******************

    if(writelots) write(11,*) 'dpsidrp= ',dpsidrp
!     
!     compute derivatives of eikonal

    if(bishop /= 0) then

       drhodpsi = drhodrp / dpsidrp

       n=ntgrid
       if(.not. allocated(a_bish)) then
          allocate(a_bish(-n:n), b_bish(-n:n), c_bish(-n:n), &
               da_bish(-n:n), db_bish(-n:n), dc_bish(-n:n), trip(-n:n))
       endif
          
       bi = btori(rgrid(0),theta(0))
       
       da_bish = 1/rmajor**2 + (bi/(Bpolmag*rmajor**2))**2
       db_bish = bi/(Bpolmag*rmajor)**2
       dc_bish = bi*(sin(th_bish) + rmajor/Rpol)/(-Bpolmag*rmajor**4)
       
       call tripprod2dtgrid(rpgrad, thgrad, rgrid, trip)
       if (nperiod > 1) call periodic_copy(trip, 0.) 

       if(iflux == 0 .or. iflux == 2) trip = trip * dpsidrp

!CMR, June06
if (debug) write(6,*) "eikcoefs:  gradients"
if (debug) write(6,*) 'drhodpsi', drhodpsi
if (debug) write(6,*) 'bi', bi
if (debug) write(6,*) 'rmajor', rmajor(:)
if (debug) write(6,*) 'Bpolmag', Bpolmag(:)
if (debug) write(6,*) rpgrad(-nth:nth,1)
if (debug) write(6,*) rpgrad(-nth:nth,2)
if (debug) write(6,*) thgrad(-nth:nth,1)
if (debug) write(6,*) thgrad(-nth:nth,2)
if (debug) write(6,*) 1.0d0/trip(-nth:nth)
if (debug) write(6,*) -Rpol(-nth:nth)/bpolmag(-nth:nth)
if (debug) write(6,*) "eikcoefs: end gradients"
!CMRend


       da_bish = da_bish/trip
       db_bish = db_bish/trip
       dc_bish = dc_bish/trip
       
       call integrate(da_bish, theta, a_bish, ntgrid)
       call integrate(db_bish, theta, b_bish, ntgrid)
       call integrate(dc_bish, theta, c_bish, ntgrid)

       a_b = a_bish(nth)-a_bish(-nth)
       b_b = b_bish(nth)-b_bish(-nth)
       c_b = c_bish(nth)-c_bish(-nth)
       
       a_bish = - a_bish * Bpolmag * rmajor
       b_bish = - b_bish * Bpolmag * rmajor
       c_bish = - c_bish * Bpolmag * rmajor
       
! find shat, pressure' and report them

       if(iflux == 0 .or. iflux == 2) then
          dp = dpdrhofun(rhoc)*drhodpsi
       else          
          dp = dpfun(rgrid(0),theta(0))
       endif

       di = dbtori(rgrid(0),theta(0))

       if (abs(qval) > epsilon(0.)) then
          tmp = rho/qval/(2.*pi)/drhodpsi
          s_hat = tmp*(a_b*di+ b_b*dp + c_b*2)
!CMR, June06
          if (debug) write(6,*) "a_b (f'), b_b (p'), c_b terms =", a_b*di, b_b*dp, c_b*2
          if (debug) write(6,*) 'qval, rho, drhodpsi, s_hat =', qval, rho, drhodpsi, s_hat 
!CMRend
       else
          tmp = 1.
          s_hat = 0.
       end if

       pressure = pfun(rgrid(0),theta(0))

       write(11,*) 
       write(11,*) 'Quantity, equilibrium value, value used'

       select case (bishop)
       case (2)
          s_hat_new = s_hat_input
          dp_new = p_prime_input
          di_new = (s_hat_new/tmp -2*c_b -dp_new*b_b) / a_b
          write(11,*) 'p_prime = ',dp,', ',dp_new
          if(writelots) write(*,*) 'p_prime = ',dp,', ',dp_new
       case (3)
          s_hat_new = s_hat_input
          dp_new = -invLp_input*pressure*drhodpsi
          di_new = (s_hat_new/tmp -2*c_b -dp_new*b_b) / a_b
          write(11,*) '1/L_p = ',-dp/drhodpsi/pressure,', ',-dp_new/drhodpsi/pressure
          if(writelots) write(*,*) '1/L_p = ',-dp/drhodpsi/pressure,', ',-dp_new/drhodpsi/pressure
       case (4)
          s_hat_new = s_hat_input
          dp_new = 0.5*beta_prime_input*drhodpsi
          di_new = (s_hat_new/tmp -2*c_b -dp_new*b_b) / a_b
          write(11,*) 'd beta/d rho = ',2.*dp/drhodpsi,', ',2.*dp_new/drhodpsi
          if(writelots) write(*,*) 'd beta/d rho = ',2.*dp/drhodpsi,', ',2.*dp_new/drhodpsi
       case (5) 
          s_hat_new = s_hat_input
          dp_new = -alpha_input/qval**2/rmaj*drhodpsi
          di_new = (s_hat_new/tmp -2*c_b -dp_new*b_b) / a_b
          write(11,*) 'alpha = ', -qval**2*rmaj*dp/drhodpsi,', ', &
               -qval**2*rmaj*dp_new/drhodpsi
          if(writelots) write(*,*) 'alpha = ', -qval**2*rmaj*dp/drhodpsi,', ', &
               -qval**2*rmaj*dp_new/drhodpsi
       case (6)
          s_hat_new = s_hat
          dp_new = 0.5*beta_prime_input*drhodpsi
          di_new = (s_hat_new/tmp -2*c_b -dp_new*b_b) / a_b
          write(11,*) 'd beta/d rho = ',2.*dp/drhodpsi,', ',2.*dp_new/drhodpsi
          if(writelots) write(*,*) 'd beta/d rho = ',2.*dp/drhodpsi,', ',2.*dp_new/drhodpsi
       case (7)
          s_hat_new = s_hat
          dp_new = dp * dp_mult
          di_new = (s_hat_new/tmp -2*c_b -dp_new*b_b) / a_b
          beta_prime_new=2.*dp_new/drhodpsi
          write(11,*) 'd beta/d rho = ',2.*dp/drhodpsi,', ',2.*dp_new/drhodpsi
          if(writelots) write(*,*) 'd beta/d rho = ',2.*dp/drhodpsi,', ',2.*dp_new/drhodpsi
       case (8)
          s_hat_new = s_hat_input
          dp_new = 0.5*beta_prime_input*drhodpsi * dp_mult
          di_new = (s_hat_new/tmp -2*c_b -dp_new*b_b) / a_b
          write(11,*) 'd beta/d rho = ',2.*dp/drhodpsi,', ',2.*dp_new/drhodpsi
          if(writelots) write(*,*) 'd beta/d rho = ',2.*dp/drhodpsi,', ',2.*dp_new/drhodpsi
       case (9)
          s_hat_new = s_hat_input
          dp_new = 0.5*beta_prime_input*drhodpsi
          di_new = (s_hat_new/tmp -2*c_b -dp_new*b_b) / a_b
          write(11,*) 'd beta/d rho = ',2.*dp/drhodpsi,', ',2.*dp_new/drhodpsi
          if(writelots) write(*,*) 'd beta/d rho = ',2.*dp/drhodpsi,', ',2.*dp_new/drhodpsi
       case default
          if(iflux /= 1) then
             s_hat_new = s_hat_input
             dp_new = dp
             di_new = (s_hat_new/tmp -2*c_b -dp_new*b_b) / a_b
          else
             s_hat_new = s_hat
             dp_new = dp
             di_new = di
          endif
          write(11,*) 'd beta/d rho = ',2.*dp/drhodpsi,', ',2.*dp_new/drhodpsi
          if(writelots) write(*,*) 'd beta/d rho = ',2.*dp/drhodpsi,', ',2.*dp_new/drhodpsi
      end select

       write(11,*) 's_hat ',s_hat,', ',s_hat_new
       if(writelots) write(*,*) 's_hat ',s_hat,', ',s_hat_new
       write(11,*) 

       shat = s_hat_new
       dqdrp = s_hat_new*qval/rho*drhodrp
       dbetadrho = 2*dp_new/drhodpsi

       call gradl(ltheta, seik, dSdl, 2*pi*qval, nth)
       if(nperiod > 1) call periodic_copy(dSdl, 0.)

       ! here the coordinates are rho and l
       grads(:,1) = a_bish*di_new + b_bish*dp_new + c_bish*2
       grads(:,2) = dSdl

       ! obtain Z(theta)
       call Ztgrid(rgrid, theta, Zoftheta)
       if(nperiod > 1) call periodic_copy(Zoftheta, 0.)
       ! calculate dZ/dl
       call gradl(ltheta, Zoftheta, dZdl, 0., nth)
       if(nperiod > 1) call periodic_copy(dZdl, 0.)

       ! gradztot(:,1) is grad Z . rhohat
       gradztot(:,1) = crpgrad(:,2)/grho
       ! gradztot(:,2) is grad Z . lhat
       gradztot(:,2) = dZdl
       ! gradztot(:,3) is grad Z . phihat
       gradztot(:,3) = 0.0

    else

!     rp derivative
       delrp=delrho/drhodrp
       rp1=rp-delrp
       rp2=rp+delrp
       if (iflux == 0 .or. iflux == 2) then
          dqdrho=dqdrhofun(rho)
          qval1=qval-dqdrho*delrho
          qval2=qval+dqdrho*delrho
       endif
    
       call seikon(rp1, qval1, seik1, dsdthet1, rgrid1, dum)
       call seikon(rp2, qval2, seik2, dsdthet2, rgrid2, dum)
       
       do i=-nth,nth
          dsdrp(i)=(seik2(i)-seik1(i))/(2.*delrp)
          write (*,*) 'dsdrp', dsdrp(i)
       enddo
    
       if(iflux == 1) then
          dqdrho=(qval2-qval1)/(2.*delrho)
          shat=dqdrho*rho/qval
       endif
    
       if(writelots) then
          write(11,*) 'dqdrho=',dqdrho,' qval1,2= ',qval1,qval2
          write(11,*) 's_hat= ',shat
          write (11,*) 'i,rgrid,rgrid1,rgrid2:'
          do i=-nth,nth
             write (11,*) i,rgrid(i),rgrid1(i),rgrid2(i)
          enddo
       endif
 
       s_hat_new = shat

       rpval(1)=rp
       rpval(2)=rp1
       rpval(3)=rp2

       if(iflux == 0 .or. iflux == 2) then
          dbetadrho = 2.*dpdrhofun(rhoc)
       else          
          dbetadrho = 2.*dpfun(rgrid(0),theta(0))/drhodrp
       endif

!
!     compute gradient of S
    
       ! this is really gradient of alpha
       ! note that coordinates are R and Z (not rho and l, as in bishop case)
       dqdrp=dqdrho*drhodrp
       do k=-nperiod+1,nperiod-1
          do j=1,2
             do i=-nth,nth
                itot=i+k*ntheta
                dsdrptot=dsdrp(i)-k*2*pi*dqdrp
                grads(itot,j)=rpgrad(i,j)*dsdrptot+thgrad(i,j)*dsdthet(i)
             enddo
          enddo
       enddo    
       ! gradztot(:,1) is grad Z . Rhat
       gradztot(:,1) = 0.
       ! gradztot(:,2) is grad Z . Zhat
       gradztot(:,2) = 1.
       ! gradztot(:,3) is grad Z . phihat
       gradztot(:,3) = 0.
    endif

    if(writelots) then
       write(11,*) 'i,theta,grads1,2=  '
       do i=-ntgrid,ntgrid
          write(11,*) i,theta(i),grads(i,1),grads(i,2)
       enddo
    endif
!     compute total grad S including toroidal part
    call gradstottgrid(rmajor, grads, gradstot)
     
!     compute the magnitude squared of grad S
    dpsidrho=dpsidrp/drhodrp
    if(iflux == 1 .and. .not. local_eq) drhodpsin = drhodrp/dpsidrp
    call dottgridf(gradstot, gradstot, gdsdum1)

    do k=-nperiod+1,nperiod-1
       do i=-nth,nth
          itot=i+k*ntheta
          gradrptot(itot,1)=rpgrad(i,1)
          gradrptot(itot,2)=rpgrad(i,2)
          gradrptot(itot,3)=0.
       enddo
    enddo
    call dottgridf(gradstot, gradrptot, gdsdum2)
    call dottgridf(gradrptot, gradrptot, gdsdum3)
    ! these quantities are calculated elsewhere if local_eq = T
    if (.not. local_eq) then
       gds2  = gdsdum1            *dpsidrho**2
       gds21 = gdsdum2*dqdrp      *dpsidrho**2
       gds22 = gdsdum3*dqdrp*dqdrp*dpsidrho**2
    end if

! MAB> new geometry terms gds23 and gds24 added because v_E . grad theta
! is needed to simulate low-flow gk equation
! v_E . grad theta = i*ky*phi*(rho_{r}*v_{t,r}/a**2)*(gds23 + gsd24*theta0)
! calculation of gds23 and gds24 assumes rp=rho, which is not necessarily true
! for numerical equilibria
! <MAB

    ! gradthtot is grad theta
    ! since ballooning, theta goes from -ntheta*nperiod -> ntheta*nperiod,
    ! but grad theta is periodic in theta
    do k=-nperiod+1,nperiod-1
       do i=-nth,nth
          itot=i+k*ntheta
          gradthtot(itot,1)=thgrad(i,1)
          gradthtot(itot,2)=thgrad(i,2)
          gradthtot(itot,3)=0.
       enddo
    enddo

    ! calculate grad rho dot grad theta
    call dottgridf(gradrptot, gradthtot, gdsdum4)
    ! calculate grad alpha dot grad theta
    call dottgridf(gradstot, gradthtot, gdsdum5)
    gds23 = (gdsdum4*gds2 - gdsdum5*gdsdum2*dpsidrho**2)/bmod**2
    gds24_noq = (gdsdum4*gdsdum2 - gdsdum5*gdsdum3)*(dpsidrho/bmod)**2
    gds24 = gds24_noq*dqdrp
    gds24_noq = gds24_noq*qval/rhoc

!     compute magnetic field
    call bvectortgrid(rgrid, theta, nth, rpgrad, dpsidrp, bvector)
     
!     compute curvature, grad b, and coriolis drift terms
    call drift(rgrid, rp, bvector, gradstot, gradrptot, gradztot, gradthtot, &
         dqdrp, dpsidrho, drhodrp, Bmod, Bpolmag, Rpol, th_bish, ltheta)
     
!     compute the magnetic field along theta
!    call  bmagtgrid(rgrid, bmagtg)    

    if (.not. local_eq) then
       if(itor /= 0) then
!       call periodic_copy(bmagtg, 0.) 
!       bmag = bmagtg
          if(nperiod > 1 ) call periodic_copy(Bmod, 0.) 
          bmag = Bmod
       else
          do k= -nperiod+1,nperiod-1
             do i=-nth,nth
                itot= i+k*ntheta
                bmag(itot)= 1./(1.+rgrid(i)*cos(theta(i))/rmaj)
             enddo
          enddo
       endif

       if(isym == 1) call sym(bmag, 0, ntgrid)
    end if

!     compute the gradparallel coefficient 
    if (.not. local_eq) then
       if (.not. local_eq) then
          bi = btori(rgrid(0),theta(0))
          do k=-nperiod+1,nperiod-1
             do i=-nth,nth
                itot=i+k*ntheta
                gradpar(itot) = -bi/(rmajor(i)**2*bmag(i)*dsdthet(i))
             enddo
          enddo
       end if

       if(isym == 1) call sym(gradpar, 0, ntgrid)
    end if
      
    if(equal_arc) then
       call arclength (ntheta, nperiod, gradpar, arcl)
       theta = arcl
    endif

!     compute grad rho*Jacobian of the transformation to the new coordinates: 
!     Axisymmetry assumed to get this (nu indep. of phi):
    do i=-nth,nth
       jacob(i)=abs(dpsidrho/(gradpar(i)*bmag(i)))
       ds(i)=grho(i)*jacob(i)
!         write(*,*) theta(i),' Jac_eik= ',jacob(i)
    enddo

!     compute surface area from A=Int(R sqrt(r**2+drdtheta**2) dtheta dphi)
!     does not work with equal_arc grid
!      surfarea=surfareafun(rgrid)
!      write(*,*) 'surface area= ',surfarea,' avgrmid**2'
      
!     compute the surface area from A=Int(J |grad rho| dtheta dalpha)
    call integrate(ds, theta, ans, nth)
    surfarea=2.*pi*(ans(nth)-ans(-nth))
!    write(*,*) 'surface area= ',surfarea,' avgrmid**2'
    if (writelots) write (11,*) 'surfarea=', surfarea
    

!     compute dV/drhon = 2 pi Int(J dtheta)
    call integrate(jacob, theta, ans, nth)
    dvdrhon=2.*pi*(ans(nth)-ans(-nth))
!      write(*,*) 'dV/drhon= ',dvdrhon
!      write(*,*) '< |grad rho| >= ',surfarea/dvdrhon
!      write(*,*) 

!     compute < |grad rho|**2>
    do i=-nth,nth
       grho2(i)=grho(i)**2*jacob(i)
    enddo
    grhoavg=surfarea/dvdrhon
    if(isym == 1) then 
       call sym(gds2, 0, ntgrid)
       call sym(gds21, 1, ntgrid)
       call sym(gds22, 0, ntgrid)
! not sure abou the following yet -- MAB
!       call sym(gds23, 0, ntgrid)
!       call sym(gds24, 1, ntgrid)
    endif

    write(25,*) rhoc, f_trap(bmag(-nth:nth))

    kxfac = 0.5*abs(qval/rhoc/drhodpsin)
    qsf = qval
    dum=qfun(0.5)
    if(eqinit >= 1) eqinit=0
    call plotdata (rgrid, seik, grads, dpsidrho)

    if (.not. local_eq) Bpol = bpolmag

    call dealloc_local_arrays

! Generate metric tensor elements to compare with Xanthopoulos and Jenko    
! Need to change normalization to match XJ

!    if (Xanthopoulos) then
!
!       Rref_phys = (rmaj+rgrid(0)) * avgrmid ! This is R_ref in physical units
!
!       Bref_X = 
!       Rref_X =
!
!       call alloc_Xanth (ntgrid)
!       
!       K1_X = - Rref_X * gbdrift0 ! Must be clear: changing normalization, or unnormalizing?
!       K2_X = - Rref_X * gbdrift  ! Must be clear: changing normalization, or unnormalizing?
!       gradpar_X = gradpar ! dimensionless for both of us, but same normalization?
!       B_X = bmag  ! * Bref ??
!       J_X = jacob ! coordinates are different, so this is trickier
!       g11_X = 
!       g12_X = 
!       g22_X = 
!
!    end if


contains

  subroutine alloc_local_arrays(n)

    integer n

    allocate(smodel (-n:n), &
         dsdrp      (-n:n), &
         grho2      (-n:n), &
         grho1      (-n:n), &
         seik       (-n:n), &
         seik1      (-n:n), &
         seik2      (-n:n), &
         dsdthet    (-n:n), &
         dsdthet1   (-n:n), &
         dsdthet2   (-n:n), &
         gdsdum1    (-n:n), &
         gdsdum2    (-n:n), &
         gdsdum3    (-n:n), &
         gdsdum4    (-n:n), &  ! MAB
         gdsdum5    (-n:n), &  ! MAB
         th_bish    (-n:n), &
         Rpol       (-n:n), &
         ltheta     (-n:n), &
         rgrid      (-n:n), &
         rgrid1     (-n:n), &
         rgrid2     (-n:n), &
         Bpolmag    (-n:n), &
         Bmod       (-n:n), &
         dSdl       (-n:n), &
         rmajor     (-n:n), &
         ans        (-n:n), &
         ds         (-n:n), &
         arcl       (-n:n), &
         Zoftheta   (-n:n), &
         dZdl       (-n:n))

    allocate(thgrad (-n:n,2), &
         rpgrad     (-n:n,2), &
         crpgrad    (-n:n,2), &
         grads      (-n:n,2))

    allocate(bvector(-n:n,3), &
         gradrptot  (-n:n,3), &
         gradstot   (-n:n,3), &
         gradztot   (-n:n,3), &
         gradthtot  (-n:n,3))  ! MAB

  end subroutine alloc_local_arrays

  subroutine dealloc_local_arrays

    deallocate(smodel, &
         dsdrp     , &
         grho2     , &
         grho1     , &
         seik      , &
         seik1     , &
         seik2     , &
         dsdthet   , &
         dsdthet1  , &
         dsdthet2  , &
         gdsdum1   , &
         gdsdum2   , &
         gdsdum3   , &
         gdsdum4   , &  ! MAB
         gdsdum5   , &  ! MAB
         th_bish   , &
         Rpol      , &
         ltheta    , &
         rgrid     , &
         rgrid1    , &
         rgrid2    , &
         Bpolmag   , &
         Bmod      , &
         dSdl      , &
         rmajor    , &
         ans       , &
         ds        , &
         arcl      , &
         Zoftheta  , &
         dZdl)

    deallocate(thgrad, &
         rpgrad    , &
         crpgrad   , &
         grads     )

    deallocate(bvector, &
         gradrptot , &
         gradstot  , &
         gradztot  , &
         gradthtot )  ! MAB

  end subroutine dealloc_local_arrays

end subroutine eikcoefs

  function surfareafun(rgrid)

    real surfareafun
    real, dimension(:) :: rgrid  
    real, dimension(-ntgrid:ntgrid) :: ds, drdth, ans
    real pi
    integer i
    
    write(*,*) 'surfareafun not generalized yet? needs to be checked.'

    pi=2.*acos(0.)

    drdth(-nth)=(rgrid(nth-1)-rgrid(-nth+1))/(theta(nth-1)-theta(-nth+1))

    do i=-nth+1,nth-1
       drdth(i)=(rgrid(i+1)-rgrid(i-1))/(theta(i+1)-theta(i-1))
    enddo

    drdth(nth)=(rgrid(-nth+1)-rgrid(nth-1))/(theta(-nth+1)-theta(nth-1))

    do i=-nth,nth
       ds(i)=2.*pi/invRfun(rgrid(i),theta(i),i)*sqrt(rgrid(i)**2+drdth(i)**2)
    enddo

    call integrate(ds, theta, ans, nth)
    surfareafun=ans(nth)-ans(-nth)
    

  end function surfareafun

  subroutine bvectortgrid(rgrid, theta, nth, rpgrad, dpsidrp, bvector)

    integer, intent (in) :: nth
    real, dimension(-ntgrid:  ), intent (in) :: rgrid, theta
    real, dimension(-ntgrid:,:), intent (in) :: rpgrad
    real, dimension(-ntgrid:,:), intent (out) :: bvector
    real, intent (in) :: dpsidrp
    
    real :: bi
    real, dimension(-ntgrid:ntgrid) :: Rinv

    integer :: i

    bi = btori(rgrid(0), theta(0))
    do i = -nth, nth
       Rinv(i) = invRfun(rgrid(i), theta(i), i)
    enddo

    do i=-nth,nth
       bvector(i,1) =-dpsidrp*rpgrad(i,2)
       bvector(i,2) = dpsidrp*rpgrad(i,1)
       bvector(i,3) = bi
    enddo

    do i=-nth,nth
       bvector(i,1) = bvector(i,1) * Rinv(i)
       bvector(i,2) = bvector(i,2) * Rinv(i)
       bvector(i,3) = bvector(i,3) * Rinv(i)
    enddo

    return
  end subroutine bvectortgrid

  subroutine gradstottgrid (rmajor, grads, gradstot)

    real, dimension(-ntgrid:), intent (in) :: rmajor
    real, dimension(-ntgrid:, :), intent (in) :: grads
    real, dimension(-ntgrid:, :), intent (out) :: gradstot

    integer :: i, k, itot

    do k=-nperiod+1,nperiod-1
       do i=-nth,nth
          itot=i+k*ntheta
          gradstot(itot,1)=grads(itot,1)
          gradstot(itot,2)=grads(itot,2)
          if(itor == 0) then
             gradstot(itot,3)=0.
          else
             gradstot(itot,3)=1./rmajor(i)
          endif
       enddo
    enddo

  end subroutine gradstottgrid

  subroutine crosstgrid(a, b, c)

    real, dimension(-ntgrid:, :), intent (in) ::  a, b
    real, dimension(-ntgrid:, :), intent (out) ::  c
    integer :: i

    do i = -nth, nth
       c(i,1)=a(i,2)*b(i,3)-b(i,2)*a(i,3)
       c(i,2)=a(i,3)*b(i,1)-b(i,3)*a(i,1)
       c(i,3)=a(i,1)*b(i,2)-b(i,1)*a(i,2)
    enddo

  end subroutine crosstgrid

  subroutine dottgrid(a, b, c) 

    real, dimension (-ntgrid:, :) :: a, b
    real, dimension (-ntgrid:) :: c
    integer :: i, k, itot

    do k = -nperiod+1, nperiod-1
       do i = -nth, nth
          itot = i + k*ntheta
          c(itot) = a(i,1)*b(itot,1) + a(i,2)*b(itot,2) + a(i,3)*b(itot,3)
       enddo
    enddo

  end subroutine dottgrid

  subroutine dottgridf(a, b, c)

    real a(-ntgrid:,:),b(-ntgrid:,:),c(-ntgrid:)

    c = a(:,1)*b(:,1) + a(:,2)*b(:,2) + a(:,3)*b(:,3)

  end subroutine dottgridf

  subroutine bmagtgrid(rgrid, bmagtg)

    real, dimension (-ntgrid:), intent (in) :: rgrid
    real, dimension (-ntgrid:), intent (out) :: bmagtg

    integer i      

    do i=-nth,nth
       bmagtg(i)=bmagfun(rgrid(i),theta(i),i)
    enddo

  end subroutine bmagtgrid
      
  real function bmagfun(r, thet, i)

    integer, intent (in) :: i
    real, intent (in) :: r, thet
    real :: bt
    
    if(itor == 0) then
       bmagfun = abs(btori(r, thet)*invRfun(r, thet, i))
    else if(iflux == 1 .and. gen_eq) then
       bmagfun = bmodfun(r, thet)
! needs to be upgraded!!!
    else 
       bt = btori(r, thet)*invRfun(r, thet, i)
       write(*,*) 'error: bmagfun not defined.'
       stop
!       bmagfun = sqrt(bt**2+bpmagfun(r, thet)**2)
    endif
    
  end function bmagfun

  function Rpos(r, thet, i)

    use  geq, only: geq_R => Rpos
    use  peq, only: peq_R => Rpos
    use  leq, only: leq_R => Rpos
    use  eeq, only: eeq_R => Rpos
    use ideq, only: ideq_R => Rpos

    integer, intent (in) :: i
    real, intent (in) :: r, thet
    real :: rpos

    if(gen_eq)  Rpos = geq_R (r, thet)       
    if(ppl_eq)  Rpos = peq_R (r, thet)       
    if(efit_eq) Rpos = eeq_R (r, thet)
    if(idfit_eq)Rpos = ideq_R (r, thet)
    if(local_eq)Rpos = leq_R (r, thet, i)

  end function Rpos

  function Zpos(r, thet, i)

    use  geq, only: geq_Z => Zpos
    use  peq, only: peq_Z => Zpos
    use  leq, only: leq_Z => Zpos
    use  eeq, only: eeq_Z => Zpos
    use ideq, only: ideq_Z => Zpos

    integer, intent (in) :: i
    real, intent (in) :: r, thet
    real :: Zpos

    if(gen_eq)  Zpos =  geq_Z (r, thet)       
    if(ppl_eq)  Zpos =  peq_Z (r, thet)       
    if(efit_eq) Zpos =  eeq_Z (r, thet)
    if(idfit_eq)Zpos = ideq_Z (r, thet)
    if(local_eq)Zpos =  leq_Z (r, thet, i)

  end function Zpos
  
  function invRfun(r, thet, i)

    use  geq, only: geq_invR => invR
    use  peq, only: peq_invR => invR
    use  leq, only: leq_invR => invR
    use  eeq, only: eeq_invR => invR
    use ideq, only: ideq_invR => invR

    integer, intent (in) :: i
    real, intent (in) :: r, thet
    real :: invRfun

    if(itor == 0) then
       invRfun=1./rmaj
    else
       if(gen_eq)  invRfun =  geq_invR (r, thet)
       if(ppl_eq)  invRfun =  peq_invR (r, thet)
       if(efit_eq) invRfun =  eeq_invR (r, thet)
       if(idfit_eq)invRfun = ideq_invR (r, thet)
       if(local_eq)invRfun =  leq_invR (r, thet, i)
    endif
    
  end function invRfun

  subroutine drift(rgrid, rp, bvector, gradstot,  &
       gradrptot, gradztot, gradthtot, dqdrp, dpsidrho, drhodrp, Bmod, &
       Bpolmag, Rpol, th_bish, ltheta)

    real, dimension (-ntgrid:), intent (in) :: rgrid, Bmod, Bpolmag, Rpol, th_bish, ltheta
    real, dimension (-ntgrid:, :), intent (in) :: bvector, gradstot, gradrptot, gradztot, gradthtot
    real, intent (in) :: rp, dqdrp, dpsidrho, drhodrp

    real, dimension (-ntgrid:ntgrid, 3) :: bgradtot, pgradtot, &
         dummy, dummy1, curve
    real, dimension (-ntgrid:ntgrid, 2) :: pgrad, bgrad1
    real, dimension (-ntgrid:ntgrid) :: gbdrift1, gbdrift2, cvdrift1, cvdrift2, cdrift1, cdrift2
    real, dimension (-ntgrid:ntgrid) :: gbdrift3, cvdrift3

    real :: dum
    character(1) char
    integer :: i, k, itot, ndum

    ndum = 2*nth + 1

    if (itor == 0) then
!       if (.not. local_eq) then
          gbdrift = (-sin(theta)*gradstot(:,1)-cos(theta)*gradstot(:,2))/rmaj
          gbdrift=2.*dpsidrho*gbdrift
          
          gbdrift0=(-sin(theta)*gradrptot(:,1)-cos(theta)*gradrptot(:,2))/rmaj
          gbdrift0=2.*dpsidrho*gbdrift0*dqdrp
          
          cvdrift=gbdrift
          cvdrift0=gbdrift0      
!       end if

    else
       char='B'
       if(bishop >= 1) then
          call bishop_gradB(rgrid, Bmod, Bpolmag, Rpol, th_bish, ltheta, bgrad1)
       else
          call grad(rgrid, theta, bgrad1, char, dum, nth, ntgrid)
       endif

       do i=-nth,nth
          bgradtot(i,1)=bgrad1(i,1)
          bgradtot(i,2)=bgrad1(i,2)
          bgradtot(i,3)=0.
       enddo

!     create curvature
       char='R'
       select case (bishop) 
       case (0) 
          call grad(rgrid, theta, pgrad, char, rp, nth, ntgrid)
       case (1) 
          pgradtot = gradrptot*dpsidrho*drhodrp*dpfun(rgrid(0), theta(0))
       case default
          pgradtot = gradrptot*dpsidrho*drhodrp*dp_new
       end select

       if(bishop < 1) then
          do i=-nth,nth
             pgradtot(i,1) = pgrad(i,1)
             pgradtot(i,2) = pgrad(i,2)
             pgradtot(i,3) = 0.
          enddo
       endif

       do i=-nth,nth
          curve(i,1)=bgradtot(i,1)/bmod(i)+pgradtot(i,1)/bmod(i)**2
          curve(i,2)=bgradtot(i,2)/bmod(i)+pgradtot(i,2)/bmod(i)**2
          curve(i,3)=0.
       enddo

!       write(*,*) 'curve 1, 2 ',curve(0,1), curve(0,2)

       call crosstgrid(bvector,bgradtot,dummy1)
       call dottgrid(dummy1,gradstot,gbdrift1)
       call dottgrid(dummy1,gradrptot,gbdrift2)
       call dottgrid(dummy1,gradthtot,gbdrift3)
       call crosstgrid(bvector,curve,dummy)
       call dottgrid(dummy,gradstot,cvdrift1)
       call dottgrid(dummy,gradrptot,cvdrift2)
       call dottgrid(dummy,gradthtot,cvdrift3)
       ! factor of 2 appears below because will later multiply
       ! by wunits, which has a factor of 1/2 built-in
       do k=-nperiod+1,nperiod-1
          do i=-nth,nth
             itot=i+k*ntheta
!             if (.not. local_eq) then
                gbdrift(itot) =2.*dpsidrho*gbdrift1(itot)/bmod(i)**3
                gbdrift0(itot)=2.*dpsidrho*gbdrift2(itot)*dqdrp/bmod(i)**3
                cvdrift(itot) =2.*dpsidrho*cvdrift1(itot)/bmod(i)**2
                cvdrift0(itot)=2.*dpsidrho*cvdrift2(itot)*dqdrp/bmod(i)**2
!             end if
             ! below are (v_M . grad theta) terms
             gbdrift_th(itot)=gbdrift3(itot)/bmod(i)**3
             cvdrift_th(itot)=cvdrift3(itot)/bmod(i)**2 
          enddo
       enddo
    endif

    ! coriolis drift term is -2*vpa*om_phi/Omega * (bhat x (zhat x bhat)) . grad (g+q<phi>F_0/T)
    ! here we calculate 2*(zhat . grad (alpha+q*theta0))/(B/Bref)

    ! this is the Zhat . grad(alpha) part of Zhat . grad(S)
    ! note that a factor of k_y = (n0/a) (drho/dpsiN) has been factored out
    call dottgrid(gradztot,gradstot,cdrift1)
    ! this is the Zhat . grad(q) part of Zhat . grad(S)
    ! note that a factor of k_y * theta0 has been factored out
    call dottgrid(gradztot,gradrptot,cdrift2)

    ! factor of 4 (instead of 2) below because will later mutliply
    ! by wunits, which has a built-in factor of 1/2
    do k=-nperiod+1,nperiod-1
       do i=-nth,nth
          itot=i+k*ntheta
          ! cdrift0 is part of coriolis drift dependent on theta0
          cdrift(itot)  = -4.*dpsidrho*cdrift1(itot)/bmod(i)
          cdrift0(itot) = -4.*dpsidrho*cdrift2(itot)*dqdrp/bmod(i)
       end do
    end do

!    if(isym == 1 .and. .not. local_eq) then
    if(isym == 1) then
       call sym(gbdrift, 0, ntgrid)
       call sym(cvdrift, 0, ntgrid)
       call sym(gbdrift0, 1, ntgrid)
       call sym(cvdrift0, 1, ntgrid)
! should maybe add in cdrift, cdrift0, gbdrift_th, and cvdrift_th here -- MAB
    endif

!    do k=-nperiod+1,nperiod-1
!       do i=-nth,nth
!          itot=i+k*ntheta
!          write (*,'(a6,6e12.4)') 'cveik', theta(itot), gbdrift0(itot), cvdrift0(itot), gbdrift(itot), cvdrift(itot), gds22(itot)
!       end do
!    end do

  end subroutine drift

  subroutine sym(a, isign, ntgrid)

    integer i, isign, ntgrid
    real, dimension(-ntgrid:), intent (in out) :: a
      
    if(isign == 0) then
       do i=1,ntgrid
          a(i)=0.5*(a(i)+a(-i))
          a(-i)=a(i)
       enddo
    else
       do i=1,ntgrid
          a(i)=0.5*(a(i)-a(-i))
          a(i)=-a(-i)
       enddo
       a(0)=0.
    endif

  end subroutine sym
      
  subroutine seikon(rp, qval, seik, dsdthet, rgrid, dpsidrp)
    
    real :: rp, qval
    real, dimension(-ntgrid:) :: seik, dsdthet, rgrid
    real, dimension(-ntgrid:ntgrid, 2) :: thgrad, rpgrad
    real :: dpsidrp, dum
    
    character(1) char    
    
    rgrid = rp
    char='P'
    if (bishop == 0)  then
       call grad(rgrid, theta, rpgrad, char, dum, nth, ntgrid)
    else
       call bgrad(rgrid, theta, rpgrad, char, dum, nth, ntgrid)          
    end if

    call thetagrad(rgrid, thgrad)
    call eikonal(rgrid, rpgrad, thgrad, qval, seik, dsdthet, dpsidrp)
    
  end subroutine seikon

  subroutine eikonal(rgrid, rpgrad, thgrad, qval, seik, dsdthet, dpsidrp)

    real, dimension(-ntgrid:ntgrid) :: trip, seik, dsdthet, rgrid
    real, dimension(-ntgrid:, :) :: rpgrad, thgrad 
    real :: qval, dpsidrp, pi, bi
    integer :: i

    pi=2.*acos(0.)
    call tripprod2dtgrid(rpgrad, thgrad, rgrid, trip)
    
    bi = btori(rgrid(0),theta(0))
    do i=-nth,nth
       dsdthet(i)=-bi/trip(i)*invRfun(rgrid(i),theta(i),i)**2
    enddo
    if(nperiod > 1) call periodic_copy(dsdthet, 0.)

    call integrate(dsdthet, theta, seik, ntgrid)

    if(iflux == 0 .or. iflux == 2) then
       dpsidrp=-(seik(nth)-seik(-nth))/(2.*pi*qval)
       do i=-nth,nth
          seik(i)=seik(i)/dpsidrp
          dsdthet(i)=dsdthet(i)/dpsidrp
       enddo
    else
       dpsidrp=1.
       qval=-(seik(nth)-seik(-nth))/(2.*pi)
    endif

  end subroutine eikonal

  subroutine thetagrad(rgrid, thgrad)

    real, dimension(-ntgrid:), intent (in) :: rgrid
    real, dimension(-ntgrid:, :), intent (out) :: thgrad
    real :: dum
    character(1) char
    
    char='T'
    if(bishop == 0) then
       if (efit_eq) then
          write(*,*) 'error in thetagrad'
          stop
       endif
       call grad(rgrid, theta, thgrad, char, dum, nth, ntgrid)
    else
       call bgrad(rgrid, theta, thgrad, char, dum, nth, ntgrid)
    end if

    if(nperiod > 1) then
       call periodic_copy(thgrad(-ntgrid:ntgrid,1),0.)
       call periodic_copy(thgrad(-ntgrid:ntgrid,2),0.)
    endif

  end subroutine thetagrad

  function diameter(rp)

    use geq, only: geq_diameter => diameter, geq_init_diameter => initialize_diameter
    use peq, only: peq_diameter => diameter, peq_init_diameter => initialize_diameter
    use ideq,only: ideq_diameter => diameter, ideq_init_diameter => initialize_diameter
    use eeq, only: bound
    
    integer :: i, initd = 1

    real :: rp, pi, diameter

    pi=2.*acos(0.)
    if(rp <= rpmin .and. .not. efit_eq) then 
       diameter = 0.
       return
    endif
    
    if(eqinit == 1) initd = 1

    if(iflux /= 1) then
       diameter=2.*rp
    else
       if (gen_eq)  i=geq_init_diameter(initd)
       if (idfit_eq)i=ideq_init_diameter(initd)
       if (ppl_eq)  i=peq_init_diameter(initd)
       initd=0
         
       if(gen_eq)  diameter = geq_diameter(rp)
       if(idfit_eq)diameter = ideq_diameter(rp)
       if(ppl_eq)  diameter = peq_diameter(rp)
       if(efit_eq) diameter = rfun(rp, 0., bound(0.)) + rfun(rp, pi, bound(pi))
    endif
      
  end function diameter

  function psifun (rp)

    real, intent (in) :: rp
    real :: psifun

    psifun=min(1.,max(0.,(rp-psi_0)/(psi_a-psi_0)))
    
  end function psifun
      
  function rpfun(r, thet)

    real, intent (in) :: r, thet
    real :: rpfun

    rpfun=psi(r, thet)

  end function rpfun
      
  function rpofrho(rho)

    real :: rpofrho
    real :: a, b, xerrbi, xerrsec, fval, soln, rho
    integer :: nsolv,ier
    
    a=rpmin
    b=rpmax
    xerrbi=1.e-4
    xerrsec=1.00e-8
    nsolv=1
    if(irho == 1) then
       fval=rho*rho*phi(b)
       call root(phi, fval, a, b, xerrbi, xerrsec, nsolv, ier, soln)
    elseif(irho == 2) then
       fval=rho*diameter(b)
       call root(diameter, fval, a, b, xerrbi, xerrsec, nsolv, ier, soln)
    elseif(irho == 3) then
       fval=rho*psifun(b)
       call root(psifun, fval, a, b, xerrbi, xerrsec, nsolv, ier, soln)
    end if
    rpofrho=soln
    if(ier > 0) write(11,*) 'error in rpofrho,rho=',rho

  end function rpofrho
  
  subroutine tripprod2dtgrid(x, y, rgrid, val)

    real, dimension(-ntgrid:), intent (in) :: rgrid
    real, dimension(-ntgrid:), intent (out) :: val
    real, dimension(-ntgrid:, :), intent(in) :: x, y
    integer :: i

! factor of 1/R from grad zeta

    do i=-nth,nth
       val(i)=(x(i,1)*y(i,2)-x(i,2)*y(i,1))*invRfun(rgrid(i),theta(i),i)
    enddo
    
  end subroutine tripprod2dtgrid

  function btori(r, thet)

    use geq, only: geq_btori => btori , geq_init_btori => initialize_btori
    use peq, only: peq_btori => btori , ppl_init_btori => initialize_btori
    use eeq, only: eeq_btori => btori, efit_init_btori => initialize_btori
    use leq, only: leq_btori => btori

    real :: btori
    real, intent (in) :: r, thet
    real :: pbar, f
    real, save :: r_last, theta_last, I_last
    integer :: i, initb = 1
!
! In the present code, most calls to this routine have the same r, thet, so: 
!
    if(abs(r-r_last) < epsilon(0.) &
         .and. abs(thet-theta_last) < epsilon(0.) &
         .and. eqinit /= 1) then
       btori = I_last
       return
    endif
    
    if(eqinit == 1) initb = 1

    if(iflux == 1) then
       pbar=min(1.,max(0.,(rpfun(r,thet)-psi_0)/(psi_a-psi_0)))

       if (gen_eq) then
          i=geq_init_btori(initb)
          f=geq_btori(pbar)
       elseif (ppl_eq) then
          i=ppl_init_btori(initb)
          f=peq_btori(pbar)
       elseif (efit_eq) then
          i = efit_init_btori(initb)
          f=eeq_btori(pbar)
       endif

       if(initb == 1) initb = 0

       btori=f
    else         
       btori=leq_btori(pbar)
    endif

    r_last = r
    theta_last = thet
    I_last = btori

  end function btori

  function dbtori(r, thet)

 ! returns dI/dpsi

    use geq, only: geq_dbtori => dbtori , geq_init_dbtori => initialize_dbtori
    use peq, only: peq_dbtori => dbtori , ppl_init_dbtori => initialize_dbtori
    use leq, only: leq_dbtori => dbtori
    use eeq, only: eeq_dbtori => dbtori, efit_init_dbtori => initialize_dbtori

    real :: dbtori
    real, intent (in) :: r, thet
    real :: pbar, f
    integer :: i, initdb = 1

    if(eqinit == 1) initdb = 1

    if(iflux == 1) then
       pbar=min(1.,max(0.,(rpfun(r,thet)-psi_0)/(psi_a-psi_0)))

       if (gen_eq) then

          i=geq_init_dbtori(initdb)    ;  f=geq_dbtori(pbar)

       elseif (ppl_eq) then

          i=ppl_init_dbtori(initdb)    ;  f=peq_dbtori(pbar)

       elseif (efit_eq) then        

          i = efit_init_dbtori(initdb) ; f=eeq_dbtori(pbar)

       endif

       if(initdb == 1) initdb = 0

       dbtori=f
    else         
       dbtori=leq_dbtori(pbar)
    endif

  end function dbtori

  function iofrho(rho)
    use geq, only: geq_iofpbar => btori
    use peq, only: peq_iofpbar => btori
    use eeq, only: eeq_iofpbar => btori
    use leq, only: leq_i => btori
    
    real :: iofrho
    real, intent (in) :: rho
            
    if (gen_eq) iofrho = geq_iofpbar(pbarofrho(rho))
    if (ppl_eq) iofrho = peq_iofpbar(pbarofrho(rho))
    if(efit_eq) iofrho = eeq_iofpbar(pbarofrho(rho))
    if(local_eq) iofrho = leq_i(rho)

  end function iofrho

  subroutine integrate(arg, grid, ans, n)

    real, dimension(-ntgrid:), intent (in) :: arg, grid
    real, dimension(-ntgrid:), intent (out) :: ans
    integer :: i, n

    ans=0.
    
    do i=1,n
       ans(i)=ans(i-1)+0.5*(grid(i)-grid(i-1))*(arg(i)+arg(i-1))
    enddo
    do i=-1,-n,-1
       ans(i)=ans(i+1)+0.5*(grid(i)-grid(i+1))*(arg(i)+arg(i+1))
    enddo
       
  end subroutine integrate

  subroutine rmajortgrid(rgrid, theta, rmajor)

    real, dimension(-ntgrid:), intent(in) :: rgrid, theta
    real, dimension(-ntgrid:), intent(out) :: rmajor
    integer :: i

    do i=-nth,nth
       rmajor(i)=Rpos(rgrid(i),theta(i),i)
    enddo
    
    return
  end subroutine rmajortgrid

  subroutine Ztgrid(rgrid, theta, Zoftheta)

    real, dimension(-ntgrid:), intent(in) :: rgrid, theta
    real, dimension(-ntgrid:), intent(out) :: Zoftheta
    integer :: i

    do i=-nth,nth
       Zoftheta(i)=Zpos(rgrid(i),theta(i),i)
    enddo
    
    return
  end subroutine Ztgrid

  function qfun(pbar)
    use geq, only: geq_qfun => qfun, geq_init_q => initialize_q
    use peq, only: peq_qfun => qfun, ppl_init_q => initialize_q
    use eeq, only: eeq_qfun => qfun, efit_init_q => initialize_q
    use leq, only: leq_qfun => qfun

    integer :: i, initq = 1
    real :: pbar, qfun
    
!    if(in_nt .and. iflux == 1) then
!       write(*,*) 'qfun not update in nt.'
!stop
!endif
    
    if(iflux == 0 .or. iflux == 2) then
       qfun=leq_qfun(pbar)
       return
    endif
	    
    if(eqinit ==1 ) initq = 1

    if (gen_eq) then
       i = geq_init_q(initq)
       qfun = geq_qfun(pbar)
    elseif (ppl_eq) then
       i = ppl_init_q(initq)
       qfun = peq_qfun(pbar)
    elseif (efit_eq) then
       i = efit_init_q(initq)
       qfun = eeq_qfun(pbar)
    endif

    if(initq == 1) initq = 0

  end function qfun

  function pfun(r,thet)

    use  geq, only: geq_pfun => pfun, geq_init_pressure => initialize_pressure
    use  peq, only: peq_pfun => pfun, ppl_init_pressure => initialize_pressure
    use  eeq, only: eeq_pfun => pfun, efit_init_pressure => initialize_pressure
    use ideq, only: ideq_pfun => pfun, idfit_init_pressure => initialize_pressure
    use  leq, only: leq_pfun => pfun

    real, intent (in) :: r, thet
    real :: pfun, f
    real pbar
    integer :: i, initp = 1

    if(iflux /= 1) then
       pfun=leq_pfun(0.) 
       return
    endif

    f=0.
    pbar=min(1.,max(0.,(psi(r,thet)-psi_0)/(psi_a-psi_0)))

    if(eqinit == 1) initp = 1
    
    if (gen_eq) then
       i=geq_init_pressure(initp)
       pfun = geq_pfun(pbar)
    elseif (ppl_eq) then
       i=ppl_init_pressure(initp)
       pfun = peq_pfun(pbar)
    elseif (efit_eq) then
       i = efit_init_pressure(initp)
       pfun = eeq_pfun(pbar)
    elseif (idfit_eq) then
       i = idfit_init_pressure(initp)
       pfun = ideq_pfun(pbar)
    endif

    if(initp == 1) initp = 0
    
  end function pfun
      
  function dpfun(r, thet)

    use  geq, only: geq_dpfun => dpfun, geq_init_dpressure => initialize_dpressure
    use  peq, only: peq_dpfun => dpfun, ppl_init_dpressure => initialize_dpressure
    use  eeq, only: eeq_dpfun => dpfun, efit_init_dpressure => initialize_dpressure
    use ideq, only: ideq_dpfun => dpfun, idfit_init_dpressure => initialize_dpressure
    use  leq, only: leq_dpfun => dpfun


    real, intent (in) :: r, thet
    real :: dpfun, f
    real pbar
    integer :: i, initdp = 1

    if(iflux /= 1) then
       dpfun=leq_dpfun(0.) 
       return
    endif

    f=0.
    pbar=min(1.,max(0.,(psi(r,thet)-psi_0)/(psi_a-psi_0)))

    if(eqinit == 1) initdp = 1
    
    if(gen_eq)  i = geq_init_dpressure(initdp)
    if(ppl_eq)  i = ppl_init_dpressure(initdp)
    if(efit_eq) i = efit_init_dpressure(initdp)
    if(idfit_eq)i = idfit_init_dpressure(initdp)
    initdp=0
          
    if(gen_eq)   dpfun = geq_dpfun(pbar)
    if(ppl_eq)   dpfun = peq_dpfun(pbar)
    if(efit_eq)  dpfun = eeq_dpfun(pbar)
    if(idfit_eq) dpfun = ideq_dpfun(pbar)
    
  end function dpfun

  function beta_a_fun(r,thet)
    
    use  geq, only: geq_beta => betafun,  geq_init_beta => initialize_beta
    use  peq, only: peq_beta => betafun,  ppl_init_beta => initialize_beta
    use  eeq, only: eeq_beta => betafun, efit_init_beta => initialize_beta
    use ideq, only: ideq_beta => betafun, idfit_init_beta => initialize_beta

    real, intent (in) :: r, thet
    real :: beta_a_fun
    real :: pbar, f
    integer :: i, initbeta = 1

    if(iflux == 0 .or. iflux == 2) then
       beta_a_fun=0.
       return
    endif
    
    if(in_nt .and. iflux == 1) then
       if(eqinit == 1) initbeta = 1
    endif
    
    f=0.
    pbar=min(1.,max(0.,(psi(r,thet)-psi_0)/(psi_a-psi_0)))
    
    if (gen_eq) then
       i = geq_init_beta(initbeta)
       f = geq_beta(pbar)
    elseif (ppl_eq) then
       i = ppl_init_beta(initbeta)
       f = peq_beta(pbar)
    elseif(efit_eq) then
       i = efit_init_beta(initbeta)
       f = eeq_beta(pbar)
    elseif(idfit_eq) then
       i = idfit_init_beta(initbeta)
       f = ideq_beta(pbar)
    endif

    initbeta = 0
    beta_a_fun=f

  end function beta_a_fun
      
  function pbarofrho(rho)

    real, intent (in) :: rho
    real :: pbarofrho
    
!     punt for now; not used

    if(iflux /= 1) then
!       write(*,*) 'error in pbarofrho!'
       pbarofrho=0.
       return         
    endif

    pbarofrho=min(1.,max(0.,(rpofrho(rho)-psi_0)/(psi_a-psi_0)))
    
  end function pbarofrho

  function dqdrhofun(rho)
    
    real, intent (in) :: rho
    real dqdrhofun
    
    dqdrhofun=shat*qfun(pbarofrho(rho))/(rhoc)

  end function dqdrhofun

  subroutine root(f,fval,a,b,xerrbi,xerrsec,nsolv,ier,soln)

    real, external :: f
    real :: fval,a,b,a1,b1,f1,f2,f3,trial,xerrbi,xerrsec,soln,aold
    integer :: i,ier,nsolv,niter,isolv
    
    ier=0
    a1=a
    b1=b
    f1=f(a1)-fval
    f2=f(b1)-fval
    if(xerrbi < 0.) goto 1000
    if(f1*f2 > 0.) then
       write(11,*) 'f1 and f2 have same sign in bisection routine'
       write(11,*) 'a1,f1,b1,f2=',a1,f1,b1,f2
       write(11,*) 'fval=',fval
       ier=1
       goto 1000
    endif
    niter=int(log(abs(b-a)/xerrbi)/log(2.))+1
    do i=1,niter
       trial=0.5*(a1+b1)
       f3=f(trial)-fval
       if(f3*f1 > 0.) then
          a1=trial
          f1=f3
       else
          b1=trial
          f2=f3
       endif
!      write(11,*) 'i,a1,f1,b1,f2 ',i,a1,f1,b1,f2
    enddo
 1000 continue
    if( abs(f1) > abs(f2) ) then
       f3=f1
       f1=f2
       f2=f3
       aold=a1
       a1=b1
       b1=aold
    endif
!     start secant method
!     write(11,*) 'a1,f1,b1,f2',a1,f1,b1,f2
    isolv=0
    do i=1,10
       aold=a1
       f3=f2-f1
       if (abs(f3) < 1.e-11)  f3=1.e-11
       a1=a1-f1*(b1-a1)/f3
       f2=f1
       b1=aold
!       write(11,*) 'a1,f1,b1,f2',a1,f1,b1,f2
       if(abs(a1-b1) < xerrsec) isolv=isolv+1
       if (isolv >= nsolv) goto 9000
       f1=f(a1)-fval
    enddo
9000 soln=a1

  end subroutine root

  function phi(rp)
      
    integer, parameter :: nimax = 200

    real :: phi
    real, intent (in) :: rp
    real :: pbar, pb(nimax), dpb
    integer :: i, ni
      
    ni=200
    if(ni > nimax) then
       write(*,*) 'Increase nimax in phi function'
       stop
    endif
    
    pbar=(rp-psi_0)/(psi_a-psi_0)
    if(.not.between(pbar,0.,1.)) then
       write(*,*) 'Defining Phi outside of plasma'
       write(*,*) 'rp, psi_0, psi_a = ',rp,psi_0,psi_a
    endif
    
    phi=0.      
    do i=1,ni
       pb(i)=pbar*float(i-1)/float(ni-1)
    enddo
    dpb=pbar/float(ni-1)
    do i=1,ni-1 ! fixed by MB
       phi=phi+0.5*dpb*(qfun(pb(i+1))+qfun(pb(i)))
    enddo
    phi=phi*(psi_a-psi_0)   ! fixed by MB

  end function phi
      
  function psi(r, thet)

    use  geq, only: geq_psi => psi,  geq_init_psi => initialize_psi
    use  peq, only: peq_psi => psi,  ppl_init_psi => initialize_psi
    use  eeq, only: eeq_psi => psi, efit_init_psi => initialize_psi
    use ideq, only: ideq_psi => psi, idfit_init_psi => initialize_psi
    use  leq, only: leq_psi => psi

    real :: psi
    real, intent (in) :: r, thet
    
    integer :: init = 1, i
     
    if(abs(r) < epsilon(0.)) then
       psi=psi_0
       return
    endif
    if(abs(r-1.) < epsilon(0.) .and. abs(thet) < epsilon(0.)) then
       psi=psi_a
       return
    endif
!
! Reinitialize if equilibrium has changed.
!
    if(eqinit == 1) init = 1

    if (gen_eq) then
       i = geq_init_psi(init)
       psi = geq_psi(r, thet)       
    elseif (ppl_eq) then
       i = ppl_init_psi(init)
       psi = peq_psi(r, thet)       
    else if(efit_eq) then
       i = efit_init_psi(init)
       psi = eeq_psi(r, thet)
    else if(idfit_eq) then
       i = idfit_init_psi(init)
       psi = ideq_psi(r, thet)
    else if(local_eq) then
       psi = leq_psi(r, thet)
    endif

    init = 0

  end function psi

  function between(x,x1,x2)

    logical :: between
    real :: x, x1, x2
    if( ((x1 <= x).and.(x <= x2)) .or. ((x1 >= x).and.(x >= x2)) ) then
       between = .true.
    else
       between = .false.
    end if

  end function between

  subroutine geofax(rho, t, e, d)

    real, intent (in) :: rho
    real, intent (out) :: t, e, d
    real :: a, b, aa, bb
    real, dimension(-ntgrid:ntgrid) :: rgrid
    real :: rp, hmaxu, hmaxl, rmaxu, rmaxl
    real :: rmid, Rlo, Rhi, Rinnr, Routr

    integer :: i

    rp=rpofrho(rho)

    if(efit_eq) then
       call rtg(rgrid, rp)
    else
       rgrid = rp
    endif

! do not treat up-down asymmetry completely for now
    hmaxu = -1.e6
    hmaxl = 1.e6
    rmaxu = 0.
    rmaxl = 0.

    do i=0,nth
       a = Rpos(rgrid(i), theta(i), i)
       b = Zpos(rgrid(i), theta(i), i)

       aa = Rpos(rgrid(-i), theta(-i), -i)
       bb = Zpos(rgrid(-i), theta(-i), -i)

       if(hmaxu < b ) then
          hmaxu = b
          rmaxu = a
       endif

       if(hmaxl > bb ) then
          hmaxl = bb
          rmaxl = aa
       endif

       Rlo = min(a,aa)
       Rhi = max(a,aa)
       if (Rinnr > Rlo) then
          Rinnr = Rlo
       end if
       if (Routr < Rhi) then
          Routr = Rhi
       endif
    enddo
    a = (Rpos(rgrid(0), theta(0), 0)-Rpos(rgrid(nth), theta(nth), nth))
    e = (hmaxu-hmaxl)/a
! Ro=(Routr+Rinnr)/2. should = rcenter(rp)    
    rmid=(Routr-Rinnr)/2. ! midplane minor radius of this surface
! Returning only the upper triangularity:
    rmax=rmaxu
! wasn't properly normalized before 29 Feb 08, and don't use asin:
!    t = asin(rcenter(rp)-rmax)
    t = (rcenter(rp)-rmax)/rmid
    d = rcenter(rp) - rcenter(psi_a)
end subroutine geofax
      
  ! function rmagaxis

  !   real :: rmagaxis
    
  !   if(iflux == 1) then
  !      rmagaxis=rmaj
  !   else
  !      rmagaxis=R_geo
  !   endif
    
  ! end function rmagaxis

  function rcenter(rp)

    use geq, only: geq_rcenter => rcenter, geq_init_rc => initialize_rcenter
    use peq, only: peq_rcenter => rcenter, peq_init_rc => initialize_rcenter
    use eeq, only: bound, eeq_init_rc => initialize_bound
    use leq, only: leq_rcenter => rcenter
    use ideq, only: ideq_rcenter => rcenter, ideq_init_rc => initialize_rcenter

    real, intent (in) :: rp
    real :: rcenter
    real :: pi
    integer :: i, init_rc = 1

    pi = 2.*acos(0.)

    if(eqinit == 1) init_rc = 1

    if (gen_eq)  i = geq_init_rc(init_rc) 
    if (ppl_eq)  i = peq_init_rc(init_rc) 
    if (efit_eq) i = eeq_init_rc(init_rc)
    if (idfit_eq)i = ideq_init_rc(init_rc)
    init_rc = 0

    rcenter=rmaj
    if(gen_eq)   rcenter = geq_rcenter(rp) 
    if(idfit_eq) rcenter = ideq_rcenter(rp) 
    if(ppl_eq)   rcenter = peq_rcenter(rp) 

    if(efit_eq)  rcenter = rmaj + 0.5*(rfun(rp, 0., bound(0.))-rfun(rp,pi,bound(pi)))
    if(local_eq) rcenter = leq_rcenter(rp)

  end function rcenter

  function bmodfun(r,thet)

    use geq, only: geqitem => eqitem, eqB_psi => B_psi
    use ideq, only: ideqitem => eqitem, ideqB_psi => B_psi

    real :: bmodfun
    real, intent (in) :: r, thet

    integer :: init = 1
    real :: f

    if(eqinit == 1) init = 1 

    if (gen_eq) then
       
       call geqitem(r, thet, eqB_psi, f, 'R')
       bmodfun=f
       return
    elseif (idfit_eq) then
       
       call ideqitem(r, thet, ideqB_psi, f, 'R')
       bmodfun=f
       return
    else
       
       write(*,*) 'Stopping in bmodfun.'  
       write(*,*) 'You must use gen_eq to call bmodfun.'
       stop
       
    endif
    
! 1000 format(1x,11e16.9)
    
  end function bmodfun

  subroutine arclength (ntheta, nperiod, gpar, arcl)

    integer, intent (in) :: ntheta, nperiod
!    real, dimension(-ntgrid:), intent (in) :: theta
    real, dimension(-ntgrid:), intent (out) :: arcl
    real, dimension(-ntgrid:), intent (in out) :: gpar

    integer :: nth, j, k
    real :: pi, arcfac
    
    nth=ntheta/2
    if(2*nth /= ntheta) write(*,*) 'ntheta should be even. ',nth, ntheta
    pi=2.*acos(0.)
    
    arcl(-nth)=0.
    do j=-nth,nth-1
       arcl(j+1)=arcl(j)+0.5*(theta(j+1)-theta(j))*(1./gpar(j)+1./gpar(j+1))
    enddo
    arcfac=1./arcl(nth)
    
    do j=-nth,nth
       arcl(j)=arcl(j)*2.*pi*arcfac-pi
    enddo
    
    if(nperiod > 1) call periodic_copy(arcl, 2.*pi)
    do k=-nperiod+1,nperiod-1
       do j=-nth,nth
          gpar(j+k*ntheta)=2.*pi*arcfac
       enddo
    enddo
    
  end subroutine arclength

  subroutine loftheta(rgrid, theta, ltheta)

    real, dimension (-ntgrid:) :: rgrid, theta, ltheta
    integer :: i
    
    real :: rold, told, rpos_old, zpos_old, rpos_new, zpos_new
    real :: r, thet  !, L_tot
        
    ltheta(0)=0.
!
! this could be made more efficient
!
    rold=rgrid(0)
    told=theta(0)

    do i=1,nth
       r=rgrid(i)
       thet=theta(i)

       rpos_old = Rpos(rold, told, i-1)
       zpos_old = Zpos(rold, told, i-1)
       rpos_new = Rpos(r, thet, i)
       zpos_new = Zpos(r, thet, i)

       ltheta(i)=ltheta(i-1) &
            +sqrt((Rpos_new-Rpos_old)**2 + (Zpos_new-Zpos_old)**2)
       rold=r
       told=thet
    enddo

    rold=rgrid(0)
    told=theta(0)

    do i=-1,-nth,-1
       r=rgrid(i)
       thet=theta(i)

!       rpos_old = Rpos(rold, told, i-1)
!       zpos_old = Zpos(rold, told, i-1)
       rpos_old = Rpos(rold, told, i+1)
       zpos_old = Zpos(rold, told, i+1)
       rpos_new = Rpos(r, thet, i)
       zpos_new = Zpos(r, thet, i)

       ltheta(i)=ltheta(i+1) &
            -sqrt((Rpos_new-Rpos_old)**2 + (Zpos_new-Zpos_old)**2)
       rold=r
       told=thet
    enddo

  end subroutine loftheta

  subroutine gradl (ltheta, f, dfdl, ext, n)
    real, dimension (-ntgrid:), intent(in) :: ltheta, f
    real, dimension (-ntgrid:), intent(out) :: dfdl
    real, intent(in) :: ext
    integer :: n
    integer :: i

    dfdl(-n) = (f(n-1) + ext - f(-n+1)) &
         /(ltheta(-n) - ltheta(-n+1))/2.

    do i = -n+1, n-1
       dfdl(i) = (f(i+1) - f(i-1))/(ltheta(i+1)-ltheta(i-1))
    enddo
    
    dfdl(n) = dfdl(-n)

  end subroutine gradl

  subroutine th_bishop(rpgrad, th_bish, nth)

    integer, intent (in) :: nth
    
    real, dimension (-ntgrid:, :), intent (in) :: rpgrad
    real, dimension (-ntgrid:), intent (out) :: th_bish

    real, dimension (-ntgrid:ntgrid) :: magrp
    real, dimension (-ntgrid:ntgrid, 2) :: tvec
    real :: pi
    integer :: i

    pi=2*acos(0.)
    do i=-nth,nth
       magrp(i)=sqrt(rpgrad(i,1)**2+rpgrad(i,2)**2)
    enddo
    
! tvec is the tangent vector in the surface of the poloidal plane

    do i=-nth,nth
       tvec(i,1)=rpgrad(i,2)/magrp(i)
       tvec(i,2)=-rpgrad(i,1)/magrp(i)
    enddo
    
! need to find grad(R) and dot it with tvec.
! but grad(R) is simple -- it is just (1,0)
! so simply return acos(tvec(:,1)) 
! with corrections for left-hand quadrants
!
    do i=-nth,nth
       tvec(i,1) = min(1.,max(-1.,tvec(i,1)))
       if(tvec(i,2) > 0.) then
          th_bish(i)=2*pi-acos(max(min(tvec(i,1),1.0),-1.0))
       else               
          th_bish(i)=acos(max(min(tvec(i,1),1.0),-1.0))
       endif
    enddo

  end subroutine th_bishop

  subroutine B_pol(rgrid, theta, rpgrad, Bpolmag, nth)

    integer, intent (in) :: nth
    real, dimension (-ntgrid:), intent (in) :: rgrid, theta
    real, dimension (-ntgrid:, :), intent (in) :: rpgrad
    real, dimension (-ntgrid:)   , intent (out) :: Bpolmag
    
    real, dimension(-ntgrid:ntgrid) :: Rinv
    integer :: i
     
    do i=-nth,nth
       Rinv(i) = invRfun(rgrid(i),theta(i),i)
    enddo

    do i=-nth,nth
       Bpolmag(i)=sqrt(rpgrad(i,1)**2 + rpgrad(i,2)**2)*Rinv(i)
    enddo
    
  end subroutine B_pol

  subroutine B_mod(rgrid, theta, Bpolmag, Bmod, nth)

    integer, intent (in) :: nth
    real, dimension (-ntgrid:), intent (in) :: rgrid, theta, Bpolmag
    real, dimension(-ntgrid:), intent (out) :: Bmod
    real :: bi
    integer :: i

       bi = btori(rgrid(0),theta(0))
       do i=-nth,nth
       Bmod(i)=sqrt((bi*invRfun(rgrid(i),theta(i),i))**2 &
            +Bpolmag(i)**2)
    enddo
    
  end subroutine B_mod

  subroutine R_pol(theta, th_bish, ltheta, Rpol, nth)

    integer, intent (in) :: nth
    real, dimension(-ntgrid:), intent (in) :: th_bish, theta, ltheta
    real, dimension(-ntgrid:), intent (out) :: Rpol
    real, dimension(-ntgrid:ntgrid) :: dthdl
    real :: pi
    integer :: i, is

    pi=2.*acos(0.)
! 
! R = 1/(d theta/dl * d th_bish/d theta)
! 
! the only tricky part is near th_bish = 0.
!
! find this point:

    is = nth-1
    do i=1,nth-1
       if(abs(th_bish(i+1)-th_bish(i)) > pi) then
          is=i
          goto 10
       endif
    enddo
    
10  continue

    do i=1,is-1
       Rpol(i)=(th_bish(i+1)-th_bish(i-1))/(theta(i+1)-theta(i-1))
    enddo

    Rpol(is)=(th_bish(is+1)-2*pi-th_bish(is-1))/(theta(is+1)-theta(is-1))

    Rpol(is+1)=(th_bish(is+2)-th_bish(is)-2*pi)/(theta(is+2)-theta(is))

    do i=is+2,nth-1
       Rpol(i)=(th_bish(i+1)-th_bish(i-1))/(theta(i+1)-theta(i-1))
    enddo

    Rpol(nth)=(th_bish(-nth+1)-th_bish(nth-1))  &
         /(theta(-nth+1)-theta(nth-1)+2*pi)

    Rpol(-nth)=Rpol(nth)
    
    do i=-nth+1,0
       Rpol(i)=(th_bish(i+1)-th_bish(i-1))/(theta(i+1)-theta(i-1))
    enddo
      
    dthdl(-nth) = (theta(nth-1) - theta(-nth+1) -2*pi) &
         /(ltheta(nth-1) - ltheta(-nth+1) - 2*ltheta(nth))

    do i = -nth+1, nth-1
       dthdl(i) = (theta(i+1)-theta(i-1))/(ltheta(i+1)-ltheta(i-1))
    enddo
    
    dthdl(nth) = dthdl(-nth)

    do i=-nth,nth
       Rpol(i) = 1/(Rpol(i)*dthdl(i))
    enddo

  end subroutine R_pol

  subroutine test(rgrid, theta, Bpolmag, Bmod, Rpol, th_bish, bgrad, nth)

    integer, intent (in) :: nth

    real, dimension(-ntgrid:), intent (in) :: rgrid, theta, Bmod, Bpolmag, &
         Rpol, th_bish
    real, dimension(-ntgrid:, :), intent(in) :: bgrad

    real, dimension(-ntgrid:ntgrid) :: rmajor, bbgrad
    real :: pi, dp, bi, di, bp
    integer :: i

    pi = 2.*acos(0.)

    call rmajortgrid(rgrid, theta, rmajor)

    if(iflux /= 1) write(*,*) 'should not be using test.'
    dp = dpfun(rgrid(0),theta(0))
    bi = btori(rgrid(0),theta(0))
    di = dbtori(rgrid(0),theta(0))
    do i = -nth, nth
       bp=-Bpolmag(i)
       bbgrad(i) = bp/Bmod(i) &
            *(bp/Rpol(i) + dp*rmajor(i) &
            - bi**2*sin(th_bish(i))/rmajor(i)**3/bp)
!       bbgrad(i) = -di*bi/rmajor(i)*bp &
!            -bi**2*sin(th_bish(i))/rmajor(i)**3
    enddo

    do i=-nth,nth
       bp=-Bpolmag(i)
       write(*,100) theta(i), &
!            bi/rmajor(i)*bgrad(i,1), &
            bgrad(i,1), &
            bbgrad(i), &
            +Bp**2/Bmod(i)/Rpol(i), &
            +Bp/Bmod(i)*dp*rmajor(i), &
            -1./Bmod(i)*bi**2*sin(th_bish(i))/rmajor(i)**3, &
!            Rpol(i), &
!            th_bish(i), &
!            -di*bi/rmajor(i)*Bp, &
!            -bi**2*sin(th_bish(i))/rmajor(i)**3, &
            bbgrad(i) - bgrad(i,1)
!            Bmod(i), bgrad(i,2)
    enddo

  100 format(20(1x,g13.6))
  end subroutine test

  subroutine tdef(nthg, ntheta_returned)
    
    real :: pi
!RN>
!    integer :: nthg, nthsave, i, ntheta_returned
    integer :: nthg, i
    integer, intent(out), optional :: ntheta_returned
!<RN
!cmr Jun06: adding following debug switch
    !logical :: debug=.false.
!cmr
!    logical :: first = .true.
    
!    if(.not.first) return
!    first = .false.

    pi = 2*acos(0.)
    
!    nthsave=nth
    nth=nthg/2   ! correct, at least for geq
    if (debug) write(6,*) "tdef: nthg,nth=",nthg,nth
!    write(*,*) 'old nth: ',nthsave,' new nth: ',nth
!RN>
!    ntheta_returned=2*nth  ! guarantees even
!    ntheta = ntheta_returned
    ntheta=2*nth  ! guarantees even
    if(present(ntheta_returned)) ntheta_returned = ntheta
!<RN
    ntgrid=(2*nperiod-1)*nth
    
    !     redefine theta grid used by the rest of the code.
    !     note that the gen_eq theta grid has theta(1)=0.
    !     problems avoided by mtheta grid in eqitem.  
    
    if (debug) write(6,*) "tdef: allocated(theta)=",allocated(theta)
    if (debug .and. allocated(theta)) write(6,*) "tdef: size(theta)=",size(theta)
    if(.not.allocated(theta)) allocate(theta(-ntgrid:ntgrid))
    if (debug) write(6,*) "tdef: ntgrid=",ntgrid
    theta(0)=0.
    do i=1,nth
       theta(i)=i*pi/float(nth)
       theta(-i)=-theta(i)
    enddo

  end subroutine tdef

  subroutine alloc_Xanth (n)

    integer n

    if (allocated(g11_X)) return

    allocate (g11_X     (-n:n), &
              g12_X     (-n:n), &
              g22_X     (-n:n), &
              B_X       (-n:n), &
              J_X       (-n:n), &
              K1_X      (-n:n), &
              K2_X      (-n:n), &
              gradpar_X (-n:n))

  end subroutine alloc_Xanth

  subroutine alloc_module_arrays(n)

    integer n
!cmr Jun06: adding following debug switch
    logical :: debug=.false.
!cmr
    if (debug) write(6,*) "alloc_module_arrays: n=",n
!CMR
    allocate(grho   (-n:n), &
         bmag       (-n:n), &
         gradpar    (-n:n), &
         cvdrift    (-n:n), &
         cvdrift0   (-n:n), &
         gbdrift    (-n:n), &
         gbdrift0   (-n:n), &
         cdrift    (-n:n), &
         cdrift0    (-n:n), &
         gbdrift_th (-n:n), &
         cvdrift_th (-n:n), &
         gds2       (-n:n), &
         gds21      (-n:n), &
         gds22      (-n:n), &
         gds23      (-n:n), &
         gds24      (-n:n), &
         gds24_noq  (-n:n), &
         jacob      (-n:n), &
         Rplot      (-n:n), &
         Zplot      (-n:n), &
         aplot      (-n:n), &
         Rprime     (-n:n), &
         Zprime     (-n:n), &
         aprime     (-n:n), &
         Uk1        (-n:n), &
         Uk2        (-n:n), &
         Bpol       (-n:n), &
         dgradpardrho(-n:n), &
         dgradparBdrho(-n:n), &
         dcvdrift0drho(-n:n), &
         dgbdrift0drho(-n:n), &
         dcvdriftdrho(-n:n), &
         dgbdriftdrho(-n:n), &
         dgds2dr(-n:n), &
         dgds21dr(-n:n), &
         dgds22dr(-n:n), &
         dBdrho(-n:n))
    if (debug) write(6,*) "alloc_module_arrays: done"
  end subroutine alloc_module_arrays

  subroutine init_theta(nt)
    
    integer, intent(in) :: nt
    integer i
    real :: pi
!cmr Jun06: adding following debug switch
    logical :: debug=.false.
!cmr
    logical :: first_local = .true.

    pi = 2*acos(0.)
    ntheta=nt
    nth = nt / 2
    ntgrid = (2*nperiod - 1)*nth       
    if (debug) write(6,*) "init_theta: allocated(theta),ntgrid=",allocated(theta),ntgrid
    if (.not. first_local) deallocate (theta)
    allocate(theta(-ntgrid:ntgrid))
    first_local = .false.

    theta = (/ (i*pi/real(nth), i=-ntgrid, ntgrid ) /) 

  end subroutine init_theta

  subroutine periodic_copy(a, ext)

    real, dimension(-ntgrid:ntgrid) :: a
    
    real ext
    integer :: i, k, itot

    if(nperiod == 1) return

    do k=1-nperiod,-1
       do i=-nth,nth
          itot = i + k*ntheta
          a(itot) = a(i) + k*ext
       enddo
    enddo

    do k=1,nperiod-1
       do i=-nth,nth
          itot = i + k*ntheta
          a(itot) = a(i) + k*ext
       enddo
    enddo

  end subroutine periodic_copy

  subroutine bishop_gradB(rgrid, Bmod, Bpolmag, Rpol, th_bish, ltheta, &
       gradB)
    
    real, dimension(-ntgrid:), intent (in) :: rgrid, Bmod, Bpolmag, &
         Rpol, th_bish, ltheta
    real, dimension(-ntgrid:,:), intent(out) :: gradB
    real, dimension(-ntgrid:ntgrid) :: rmajor, dBdl
    real :: bi, bp
    integer :: i 

    call rmajortgrid(rgrid, theta, rmajor)
    call gradl(ltheta, Bmod, dBdl, 0., nth)

    bi = btori(rgrid(0),theta(0))
    do i=-nth, nth
       bp=-Bpolmag(i)
       gradB(i,1)=bp/Bmod(i) *(bp/Rpol(i) + dp_new*rmajor(i) &
            - bi**2*sin(th_bish(i))/rmajor(i)**3/bp)
       gradB(i,2)=dBdl(i)
    enddo

  end subroutine bishop_gradB

  subroutine check(geq, eeq, peq, leq, ideq)
    logical, intent(in) :: geq, eeq, peq, leq, ideq
    
    if(geq .and. eeq) then
       write(*,*) 'Choosing gen_eq = .true. AND efit_eq = .true. is not permitted.'
       write(*,*) 'Stopping.'
       stop
    endif

    if(geq .and. peq) then
       write(*,*) 'Choosing gen_eq = .true. AND ppl_eq = .true. is not permitted.'
       write(*,*) 'Stopping.'
       stop
    endif

    if(geq .and. leq) then
       write(*,*) 'Choosing gen_eq = .true. AND iflux = 0 is not permitted.'
       write(*,*) 'Stopping.'
       stop
    endif

    if(eeq .and. leq) then
       write(*,*) 'Choosing efit_eq = .true. AND iflux = 0 is not permitted.'
       write(*,*) 'Stopping.'
       stop
    endif

    if(eeq .and. peq) then
       write(*,*) 'Choosing efit_eq = .true. AND ppl_eq = .true. is not permitted.'
       write(*,*) 'Stopping.'
       stop
    endif

    if(peq .and. leq) then
       write(*,*) 'Choosing ppl_eq = .true. AND iflux = 0 is not permitted.'
       write(*,*) 'Stopping.'
       stop
    endif

  end subroutine check

  subroutine grad(rgrid, theta, gradf, char, rp, nth, ntgrid)
     
     use  geq, only:  geq_gradient => gradient
     use  peq, only:  peq_gradient => gradient
     use  eeq, only:  eeq_gradient => gradient
     use ideq, only: ideq_gradient => gradient
     use  leq, only:  leq_gradient => gradient

     integer :: nth, ntgrid
     real, dimension(-ntgrid:) :: rgrid, theta
     real, dimension(-ntgrid:,:) :: gradf
     character(1) :: char
     real rp

     if(gen_eq)   call  geq_gradient(rgrid, theta, gradf, char, rp, nth, ntgrid)
     if(ppl_eq)   call  peq_gradient(rgrid, theta, gradf, char, rp, nth, ntgrid)
     if(efit_eq)  call  eeq_gradient(rgrid, theta, gradf, char, rp, nth, ntgrid)
     if(idfit_eq) call ideq_gradient(rgrid, theta, gradf, char, rp, nth, ntgrid)
     if(local_eq) call  leq_gradient(rgrid, theta, gradf, char, rp, nth, ntgrid)

   end subroutine grad

   subroutine bgrad(rgrid, theta, gradf, char, rp, nth, ntgrid)
     
     use  geq, only:  geq_bgradient => bgradient
     use  peq, only:  peq_bgradient => bgradient
     use  eeq, only:  eeq_bgradient => bgradient
     use ideq, only: ideq_bgradient => bgradient
     use  leq, only:  leq_bgradient => bgradient

     integer :: nth, ntgrid
     real, dimension(-ntgrid:) :: rgrid, theta
     real, dimension(-ntgrid:,:) :: gradf
     character(1) :: char
     real rp

     if(gen_eq)   call  geq_bgradient(rgrid, theta, gradf, char, rp, nth, ntgrid)
     if(ppl_eq)   call  peq_bgradient(rgrid, theta, gradf, char, rp, nth, ntgrid)
     if(efit_eq)  call  eeq_bgradient(rgrid, theta, gradf, char, rp, nth, ntgrid)
     if(idfit_eq) call ideq_bgradient(rgrid, theta, gradf, char, rp, nth, ntgrid)
     if(local_eq) call  leq_bgradient(rgrid, theta, gradf, char, rp, nth, ntgrid)

   end subroutine bgrad

  subroutine rtg(rgrid, rp)

    use eeq, only: ebound => bound
    real, dimension(-ntgrid:), intent (out) :: rgrid
    real, intent (in) :: rp
    integer :: i
    
    if(.not. efit_eq) then
       write(*,*) 'error in rtg.  efit_eq = .false.'
       stop          
    endif

    if (efit_eq) then
       do i=-nth,nth
          rgrid(i)=rfun(rp, theta(i), ebound(theta(i)))
       enddo
    endif
       
  end subroutine rtg

  function rfun(rp, thet, broot)

    use eeq, only: ebound => bound
    integer :: i, j, k, imax, jmax, kmax
    real :: rfun, fa, fb, fbroot, bmult, rootval, thetroot
    real :: a, b, xerrsec, rp, thet, broot

    if(.not. efit_eq) then
       rfun = rp
       return
    endif
    
    thetroot=thet
    rootval=rp
    a=0.
    b=broot


    if(broot < 0.) then
       if (efit_eq) broot = ebound(thet)
    endif

    fa    =rpfun(a,thetroot)-rootval
    fbroot=rpfun(b,thetroot)-rootval
    fb=fbroot
    if(fa*fbroot <= 0.) goto 100
    
! need to bracket root.  

    bmult=0.01
    jmax=5
    kmax=5
       
    do k=1,kmax
       bmult=bmult*2.
       imax=10
       do j=1,jmax
          imax=imax*2
          do i=1,imax
             b=broot*(1+bmult*float(i)/float(imax))
             fb=rpfun(b,thetroot)-rootval
             if(fa*fb <= 0.) goto 100
          enddo
          do i=1,imax
             b=broot/(1+bmult*float(i)/float(imax))
             fb=rpfun(b,thetroot)-rootval
             if(fa*fb <= 0.) goto 100
          enddo
       enddo
    enddo
    write(*,*) 'yuk'
 100   continue 
    xerrsec=1.e-9
!         write(*,*) 'a, b = ',a,b
    rfun=zbrent(rpfun, a, b, rootval, thetroot, xerrsec)
!         write(*,1000) a, b, fa, fb, rfun, thet
! 1000 format(1x,11e16.9)

  end function rfun

  function zbrent(func, x1, x2, rootval, thetroot, tol)

    real :: zbrent, func
!    real, parameter :: eps = 3.e-8, eps1 = 2.e-5
    real, parameter :: eps = 3.e-8, eps1 = 2.e-8
    integer, parameter :: itmax = 100
    real, intent (in) :: x1, x2, rootval, thetroot, tol
    real :: a, b, c, fa, fb, fc, d, e, tol1, xm, q, p, r, s
    integer :: iter
    external func

    a=x1
    b=x2
    c=0.
    fa=func(a,thetroot)-rootval
    fb=func(b,thetroot)-rootval
    if(fb*fa > 0.) then
! check to see if fa or fb is close to zero already:
       if(abs(fa) < eps1) then
          zbrent=a
          return
       elseif(abs(fb) < eps1) then
          zbrent=b
          return
       endif
       write(*,*) a,b,fa,fb
       write(*,*) 'root must be bracketed for zbrent.'
       stop
    endif
    fc=fb
    do 11 iter=1,itmax
       if(fb*fc > 0.) then
          c=a
          fc=fa
          d=b-a
          e=d
       endif
       if(abs(fc) < abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       endif
       tol1=2.*eps*abs(b)+0.5*tol
       xm=.5*(c-b)
       if(abs(xm) <= tol1 .or. abs(fb) < epsilon(0.))then
          zbrent=b
          return
       endif
       if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s=fb/fa
          if (abs(a-c) < epsilon(0.)) then
             p=2.*xm*s
             q=1.-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
             q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p > 0.) q=-q
          p=abs(p)
          if(2.*p  <  min(3.*xm*q-abs(tol1*q),abs(e*q))) then
             e=d
             d=p/q
          else
             d=xm
             e=d
          endif
       else
          d=xm
          e=d
       endif
       a=b
       fa=fb
       if(abs(d)  >  tol1) then
          b=b+d
       else
          b=b+sign(tol1,xm)
       endif
       fb=func(b,thetroot)-rootval
11     continue
       write(*,*) 'zbrent exceeding maximum iterations.'
       stop
       zbrent=b
       
  end function zbrent

  subroutine Hahm_Burrell(i, a)   

    use geq, only: geq_Hahm_Burrell => Hahm_Burrell
    use peq, only: peq_Hahm_Burrell => Hahm_Burrell

    integer i
    real a

    if(gen_eq) call geq_Hahm_Burrell(i, a)
    if(ppl_eq) call peq_Hahm_Burrell(i, a)

  end subroutine Hahm_Burrell

  function nth_get()
    
    integer nth_get

    nth_get = nth

  end function nth_get

  subroutine drho_drp(rp, drhodrp)

    real :: rp1, rp2, rp, rho1, rho2
    real, intent (out) :: drhodrp

!     compute d rho / d rp
    rp1=rp*(1.-delrho)
    rp2=rp*(1.+delrho)
    if(irho == 1) then
       rho1=sqrt(abs(phi(rp1)/phi(rpmax)))
       rho2=sqrt(abs(phi(rp2)/phi(rpmax)))
    elseif(irho == 2) then
       rho1=0.5*diameter(rp1)
       rho2=0.5*diameter(rp2)
    elseif(irho == 3) then
       rho1=psifun(rp1)
       rho2=psifun(rp2)
    endif
    drhodrp=(rho2-rho1)/(rp2-rp1)

  end subroutine drho_drp

  subroutine drho_drhod(rp, drhodrp, drhodrhod)

    real :: rp1, rp2, rp, rho1, rho2
    real, intent (in) :: drhodrp
    real, intent (out) :: drhodrhod

!     compute d rho / d rb
    rp1=rp*(1.-delrho)
    rp2=rp*(1.+delrho)
    rho1=0.5*diameter(rp1)
    rho2=0.5*diameter(rp2)

    drhodrhod=drhodrp/((rho2-rho1)/(rp2-rp1))

  end subroutine drho_drhod

  function fluxavg(f)
    real, dimension(-nth:nth), intent(in) :: f
    real, dimension(-nth:nth) :: delth
    real :: fluxavg

    delth(-nth+1:nth-1) = 0.5*(theta(-nth+2:nth)-theta(-nth:nth-2))
    delth(-nth) = 0.5*(theta(-nth+1)-theta(-nth))
    delth(nth) = 0.5*(theta(nth)-theta(nth-1))

    fluxavg = sum(f*jacob(-nth:nth)*delth)/sum(jacob(-nth:nth)*delth)

  end function fluxavg

  function f_trap(b_mag)
    
    real, dimension(-nth:nth), intent (in) :: b_mag
    real :: f_trap, ftu, ftl, havg, h2avg, h(-nth:nth)

! use expressions from Y. R. Lin-Liu and Bob Miller, coding from Q. P. Liu
    
    h = min(b_mag/maxval(b_mag),1.0)
    havg = fluxavg(h)
    h2avg = fluxavg(h**2)
    ftu = 1.0 - h2avg/havg**2*(1.0-sqrt(1.0-havg)*(1.0+0.5*havg))
    ftl = 1.0 - h2avg*fluxavg((1.0-sqrt(1.0-h)*(1.0+0.5*h))/h**2)
    f_trap = 0.75*ftu + 0.25*ftl 
    
  end function f_trap

  subroutine plotdata (rgrid, seik, grads, dpsidrho)

    real, dimension (-ntgrid:), intent (in) :: rgrid, seik
    real, dimension (-ntgrid:,:), intent (in) :: grads
    real, intent (in) :: dpsidrho

    real :: rplus, rminus, dr
    integer :: i
!
! At the present flux surface, calculate these three functions 
! as a function of theta, together with their derivates with 
! respect to rho: 
! 
! R, Z, alpha-phi
!
! Call these variables
!
! Rplot, Zplot, aplot
!
! and their derivatives:
!
! Rprime, Zprime, aprime
!    
! BD: bug??
! Is this correct for EFIT data?  Need to check
!
    do i=-ntgrid,ntgrid
       Rplot(i) = Rpos(rgrid(i), theta(i), i)
       Zplot(i) = Zpos(rgrid(i), theta(i), i)
       aplot(i) = seik(i)
    end do

    do i=-ntgrid,ntgrid
       rplus  = rgrid(i)*(1.+delrho)
       rminus = rgrid(i)*(1.-delrho)
       dr = 2.*delrho*rgrid(i)
       Rprime(i) = (Rpos(rplus, theta(i), i)-Rpos(rminus, theta(i), i))/dr
       Zprime(i) = (Zpos(rplus, theta(i), i)-Zpos(rminus, theta(i), i))/dr
       aprime(i) = grads(i,1)*dpsidrho
    end do

  end subroutine plotdata


end module geometry



  !type(geometry_advanced_parameters_type) :: geometry_advanced_parameters 

!> This subroutine implements a general interface to geometry
!! which can be used outside GS2; in particular, it can be
!! called from C/C++ programs. 
!!
!! Equilibrium types:
!!     1: Local (Miller)
!!     2: EFIT
!!     3: Chease
!!     4: Transp
!!
!! Some equilibria, those that are created from a grid of psi(R,Z),
!! which means EFIT and Chease, change ntheta. 
!! This should then be used as the number of 
!! gridpoints along the fieldline. 
!!
!! Massive note to self, check where iflux and local_eq are set
!! in GS2

!contains
subroutine geometry_set_inputs(equilibrium_type,&
                               eqfile_in, &
                               irho_in, &
                               rhoc_in, &
                               bishop_in, &
                               nperiod_in, &
                               ntheta)
  use geometry

  integer, intent(in) :: equilibrium_type, nperiod_in, irho_in, bishop_in
  real, intent(in) :: rhoc_in
  integer, intent(inout) :: ntheta
  character(len=800), intent(in) :: eqfile_in

  write (*,*) 'Setting inputs'
  write (*,*)

  rhoc = rhoc_in
  nperiod = nperiod_in
  eqfile = eqfile_in
  irho = irho_in

  iflux = 1 ! Numerical unless changed for miller

  local_eq = .false.
  gen_eq = .false.
  ppl_eq = .false.
  transp_eq = .false.
  efit_eq = .false.

  equal_arc = .true.
  bishop = bishop_in
  dp_mult = 1.0
  delrho = 1e-3
  rmin = 1e-3
  rmax = 1.0
  isym = 0
  in_nt = .false.
  writelots = .false.
  !! Advanced use only
  itor = 1

  
  !ntheta_out = -1

  select case (equilibrium_type)
  case (1)  ! Miller
    local_eq = .true.
    iflux = 0
    call init_theta(ntheta)
  case (2)  ! EFIT
    efit_eq = .true.
  case default
    write(*,*) "Whatever geometry you are using hasn't "
    write(*,*) "been tested in a long time!"
    stop 
  end select



end subroutine geometry_set_inputs  


!> This subroutine gets default values for the advanced switches,
!! returning a derived type (advanced_parameters_type) containing
!! all the parameters. It does not modify the geometry variables. 

subroutine geometry_get_default_advanced_parameters(advanced_parameters_out)
  use geometry
  type(advanced_parameters_type), intent(out) :: advanced_parameters_out
 
  advanced_parameters_out%equal_arc = .true.
  advanced_parameters_out%dp_mult = 1.0
  advanced_parameters_out%delrho = 1e-3
  advanced_parameters_out%rmin = 1e-3
  advanced_parameters_out%rmax = 1.0
  advanced_parameters_out%isym = 0
  advanced_parameters_out%in_nt = .false.
  advanced_parameters_out%writelots = .false.
  advanced_parameters_out%itor =  1
end subroutine geometry_get_default_advanced_parameters

!> Returns a derived type containing values of the advanced geometry
!! parameters. 
subroutine geometry_get_advanced_parameters(advanced_parameters_out)
  use geometry
  type(advanced_parameters_type), intent(out) :: advanced_parameters_out

  advanced_parameters%equal_arc = equal_arc
  advanced_parameters%dp_mult = dp_mult
  advanced_parameters%delrho = delrho
  advanced_parameters%rmin = rmin
  advanced_parameters%rmax = rmax
  advanced_parameters%isym = isym
  advanced_parameters%in_nt = in_nt
  advanced_parameters%writelots = writelots
  advanced_parameters%itor = itor

  advanced_parameters_out = advanced_parameters

end subroutine geometry_get_advanced_parameters

!> Accepts a derived type containing values of the advanced geometry
!! parameters and sets the variables in the geometry module accordingly. 
subroutine geometry_set_advanced_parameters(advanced_parameters_in)
  use geometry
  type(advanced_parameters_type), intent(in) :: advanced_parameters_in

  write(*,*) 'Setting advanced parameters'
  advanced_parameters = advanced_parameters_in

  equal_arc = advanced_parameters%equal_arc
  dp_mult = advanced_parameters%dp_mult
  delrho = advanced_parameters%delrho
  rmin = advanced_parameters%rmin
  rmax = advanced_parameters%rmax
  isym = advanced_parameters%isym
  in_nt = advanced_parameters%in_nt
  writelots = advanced_parameters%writelots
  itor = advanced_parameters%itor

  write(*,*) 'delrho was set to', delrho

end subroutine geometry_set_advanced_parameters

subroutine geometry_get_miller_parameters(miller_parameters_out)
  use geometry
  type(miller_parameters_type), intent(out) :: miller_parameters_out


   miller_parameters%rmaj = rmaj
   miller_parameters%r_geo = r_geo
   miller_parameters%akappa = akappa
   miller_parameters%akappri = akappri
   miller_parameters%tri = tri
   miller_parameters%tripri = tripri
   miller_parameters%shift = shift
   miller_parameters%qinp = qinp
   miller_parameters%shat = shat
   miller_parameters%asym = asym
   miller_parameters%asympri = asympri
   miller_parameters%d2qdr2 = d2qdr2
   miller_parameters%betadbprim = betadbprim
   miller_parameters%d2psidr2 = d2psidr2

   miller_parameters_out = miller_parameters

end subroutine geometry_get_miller_parameters

subroutine geometry_set_miller_parameters(miller_parameters_in)
  use geometry
  type(miller_parameters_type), intent(in) :: miller_parameters_in



   miller_parameters = miller_parameters_in

   rmaj = miller_parameters_in%rmaj
   R_geo = miller_parameters_in%R_geo
   akappa = miller_parameters_in%akappa
   akappri = miller_parameters_in%akappri
   tri = miller_parameters_in%tri
   tripri = miller_parameters_in%tripri
   shift = miller_parameters_in%shift
   qinp = miller_parameters_in%qinp
   shat = miller_parameters_in%shat
   asym = miller_parameters_in%asym
   asympri = miller_parameters_in%asympri
   betadbprim = miller_parameters_in%betadbprim

   write(*,*) 's_hat was set to', shat

end subroutine geometry_set_miller_parameters

!> Vary the magnetic shear and pressure gradient 
!! while maintaing a solution to the Grad-Shrafranov eqn,
!! as described variously by Greene & Chance, Bishop and Miller.
!! 
subroutine geometry_vary_s_alpha(s_hat_input_in, beta_prime_input_in)
  use geometry
  real, intent(in) :: s_hat_input_in, beta_prime_input_in
  s_hat_input = s_hat_input_in
  beta_prime_input = beta_prime_input_in
end subroutine geometry_vary_s_alpha

subroutine geometry_calculate_coefficients(grid_size)
  use geometry
  integer, intent(out) :: grid_size
  integer :: ntheta_out
  call eikcoefs(ntheta_out)
  write (*,*) 'ntheta out is ', ntheta_out
  grid_size = (2*nperiod - 1)*ntheta_out + 1 ! = 2*ntgrid + 1
  write (*,*) 'grid_size out is ', grid_size
end subroutine geometry_calculate_coefficients

!> Get the geometric coefficients calculated by the geometry module.
subroutine geometry_get_coefficients(grid_size, coefficients_out)
  use geometry
  integer, intent(in) :: grid_size
  type(coefficients_type), dimension(grid_size) :: coefficients_out
  integer ::ntgrid, i

  !ntgrid = (2*nperiod - 1) * ntheta/2
  ntgrid = (grid_size - 1)/2

   !allocate(coefficients_out(2*ntgrid+1))   
  !write (*,*) 'HERE2'
   !!allocate(coefficients_out%bmag(ntgrid:ntgrid))       
   !allocate(coefficients_out%gradpar(ntgrid:ntgrid))    
   !allocate(coefficients_out%cvdrift(ntgrid:ntgrid))    
   !allocate(coefficients_out%cvdrift0(ntgrid:ntgrid))   
   !allocate(coefficients_out%gbdrift(ntgrid:ntgrid))    
   !allocate(coefficients_out%gbdrift0(ntgrid:ntgrid))   
   !allocate(coefficients_out%cdrift(ntgrid:ntgrid))    
   !allocate(coefficients_out%cdrift0(ntgrid:ntgrid))    
   !allocate(coefficients_out%gbdrift_th(ntgrid:ntgrid)) 
   !allocate(coefficients_out%cvdrift_th(ntgrid:ntgrid)) 
   !allocate(coefficients_out%gds2(ntgrid:ntgrid))       
   !allocate(coefficients_out%gds21(ntgrid:ntgrid))      
   !allocate(coefficients_out%gds22(ntgrid:ntgrid))      
   !allocate(coefficients_out%gds23(ntgrid:ntgrid))      
   !allocate(coefficients_out%gds24(ntgrid:ntgrid))      
   !allocate(coefficients_out%gds24_noq(ntgrid:ntgrid))  
   !allocate(coefficients_out%jacob(ntgrid:ntgrid))      
   !allocate(coefficients_out%Rplot(ntgrid:ntgrid))      
   !allocate(coefficients_out%Zplot(ntgrid:ntgrid))      
   !allocate(coefficients_out%aplot(ntgrid:ntgrid))      
   !allocate(coefficients_out%Rprime(ntgrid:ntgrid))     
   !allocate(coefficients_out%Zprime(ntgrid:ntgrid))     
   !allocate(coefficients_out%aprime(ntgrid:ntgrid))     
   !allocate(coefficients_out%Uk1(ntgrid:ntgrid))        
   !allocate(coefficients_out%Uk2(ntgrid:ntgrid))        
   !allocate(coefficients_out%Bpol(ntgrid:ntgrid))       

  write (*,*) 'HERE'
  write (*,*) 'Grid size should be ', 2*ntgrid + 1
   do i = -ntgrid,ntgrid
   write (*,*) 'i', i
   coefficients_out(i+ntgrid+1)%grho        = grho(i)   
   coefficients_out(i+ntgrid+1)%bmag        = bmag(i)       
   coefficients_out(i+ntgrid+1)%gradpar     = gradpar(i)    
   coefficients_out(i+ntgrid+1)%cvdrift     = cvdrift(i)    
   coefficients_out(i+ntgrid+1)%cvdrift0    = cvdrift0(i)   
   coefficients_out(i+ntgrid+1)%gbdrift     = gbdrift(i)    
   coefficients_out(i+ntgrid+1)%gbdrift0    = gbdrift0(i)   
   coefficients_out(i+ntgrid+1)%cdrift     = cdrift(i)    
   coefficients_out(i+ntgrid+1)%cdrift0     = cdrift0(i)    
   coefficients_out(i+ntgrid+1)%gbdrift_th  = gbdrift_th(i) 
   coefficients_out(i+ntgrid+1)%cvdrift_th  = cvdrift_th(i) 
   coefficients_out(i+ntgrid+1)%gds2        = gds2(i)       
   coefficients_out(i+ntgrid+1)%gds21       = gds21(i)      
   coefficients_out(i+ntgrid+1)%gds22       = gds22(i)      
   coefficients_out(i+ntgrid+1)%gds23       = gds23(i)      
   coefficients_out(i+ntgrid+1)%gds24       = gds24(i)      
   coefficients_out(i+ntgrid+1)%gds24_noq   = gds24_noq(i)  
   coefficients_out(i+ntgrid+1)%jacob       = jacob(i)      
   coefficients_out(i+ntgrid+1)%Rplot       = Rplot(i)      
   coefficients_out(i+ntgrid+1)%Zplot       = Zplot(i)      
   coefficients_out(i+ntgrid+1)%aplot       = aplot(i)      
   coefficients_out(i+ntgrid+1)%Rprime      = Rprime(i)     
   coefficients_out(i+ntgrid+1)%Zprime      = Zprime(i)     
   coefficients_out(i+ntgrid+1)%aprime      = aprime(i)     
   coefficients_out(i+ntgrid+1)%Uk1         = Uk1(i)        
   coefficients_out(i+ntgrid+1)%Uk2         = Uk2(i)        
   coefficients_out(i+ntgrid+1)%Bpol        = Bpol(i)       
   end do

   write (*,*) 'Returning....'

end subroutine geometry_get_coefficients

subroutine geometry_get_constant_coefficients(constant_coefficients_out)
  use geometry
  type(constant_coefficients_type), intent(out) :: constant_coefficients_out


  constant_coefficients%qsf = qsf
  constant_coefficients%rmaj = rmaj
  constant_coefficients%shat = shat
  constant_coefficients%kxfac = kxfac
  constant_coefficients%aminor = aminor
  constant_coefficients_out = constant_coefficients
  
end subroutine geometry_get_constant_coefficients
