module periodic_splines
  use warning_helpers, only: is_not_zero, is_zero, not_exactly_equal, exactly_equal
 
  implicit none

  private

  public :: spline, new_spline, delete_spline
  public :: periodic_spline, new_periodic_spline, delete_periodic_spline
  public :: splint, dsplint, periodic_splint, periodic_dsplint
  public :: inter_d_cspl, inter_cspl
  public :: fitp_curv1, fitp_curvp1, fitp_curv2, fitp_curvp2
  public :: fitp_surf1, fitp_surf2, fitp_curvd

  interface handle_spline_error
     module procedure handle_spline_error_logical
  end interface handle_spline_error

  !> Holds data representing a non-periodic spline. Should be set up
  !> by calling [[new_spline]].
  type :: spline
     private
     !> Length of the data arrays represented by the spline
     integer :: n = 0
     !> Holds the independent and dependent values of the
     !> splined data in `x` and `y`. The second derivative
     !> is held in `y2` and calculated automatically.
     real, dimension (:), allocatable :: x, y, y2
     !> Indicates if the spline corresponding to this data is valid
     !> and can be used with the spline evaluation routines.
     logical, public :: valid = .false.
     !> The tension used in computing the splined data, note this
     !> must be the value used in the initialisation when passed
     !> to the spline evaluation routines.
     real :: tension = 1.0
   contains
     procedure :: interpolate => spline_interp
     procedure :: derivative => spline_deriv
  end type spline

  !> Constructor for spline
  interface spline
     module procedure new_spline
  end interface spline

  !> Holds data representing a periodic spline. Should be set up by
  !> calling [[new_periodic_spline]].
  type :: periodic_spline
     private
     !> Length of the data arrays represented by the spline
     integer :: n = 0
     !> The actual size of the periodic domain
     real :: period = 0
     !> Holds the independent and dependent values of the
     !> splined data in `x` and `y`. The second derivative
     !> is held in `y2` and calculated automatically.
     real, dimension (:), allocatable :: x, y, y2
     !> Indicates if the spline corresponding to this data is valid
     !> and can be used with the spline evaluation routines.
     logical, public :: valid = .false.
     !> The tension used in computing the splined data, note this
     !> must be the value used in the initialisation when passed
     !> to the spline evaluation routines.
     real :: tension = 1.0
   contains
     procedure :: interpolate => periodic_spline_interp
     procedure :: derivative => periodic_spline_deriv
  end type periodic_spline

  !> Constructor for periodic_spline
  interface periodic_spline
     module procedure new_periodic_spline
  end interface periodic_spline

contains

  !> Populates a spline instance `spl` representing the non-periodic
  !> data y(x).
  type(spline) function new_spline (x, y, tension) result(spl)
    use optionals, only: get_option_with_default
    implicit none
    real, dimension (:), intent (in) :: x, y
    real, intent(in), optional :: tension
    real, dimension(:), allocatable :: temp
    integer :: ierr
    spl%valid = .false.
    spl%tension = get_option_with_default(tension, 1.0)
    spl%n = size(x)
    allocate(spl%x, source = x) ; allocate(spl%y, source = y)
    allocate (spl%y2(spl%n), temp(spl%n))
    call fitp_curv1 (spl%n, spl%x, spl%y, 0.0, 0.0, 3, spl%y2, temp, spl%tension, ierr)
    spl%valid = ierr == 0
  end function new_spline

  !> Populates a periodic_spline instance `spl` representing the
  !> periodic data y(x) of length n and periodic on `period`.  Note
  !> that the spline library expects `period > x(n) - x(1)`, which
  !> means the input data shouldn't include the duplicate periodic
  !> point.  As a convenience the user can pass data with the
  !> duplicate point and set `drop_last_point = .true.` to
  !> automatically exclude the duplicate point.
  type(periodic_spline) function new_periodic_spline (x, y, period, &
       drop_last_point, tension) result(spl)
    use optionals, only: get_option_with_default
    implicit none
    real, dimension (:), intent (in) :: x, y
    real, intent (in) :: period
    logical, intent(in), optional :: drop_last_point
    real, intent(in), optional :: tension
    logical :: drop_point
    real, dimension (:), allocatable :: temp
    integer :: ierr
    spl%valid = .false.
    drop_point = get_option_with_default(drop_last_point, .false.)
    spl%tension = get_option_with_default(tension, 1.0)
    spl%n = size(x)
    if (drop_point) spl%n = spl%n - 1
    allocate(spl%x, source = x(:spl%n)) ; allocate(spl%y, source = y(:spl%n))
    allocate (spl%y2(spl%n), temp(2*spl%n))
    spl%period = period
    call fitp_curvp1 (spl%n,spl%x,spl%y,spl%period,spl%y2,temp,spl%tension,ierr)
    spl%valid = ierr == 0
  end function new_periodic_spline

  !> Reset and deallocate variables in passed spline
  subroutine delete_spline (spl)
    implicit none
    type (spline), intent (in out) :: spl
    spl%n = 0
    if (allocated(spl%x)) deallocate (spl%x,spl%y)
    if (allocated(spl%y2)) deallocate (spl%y2)
    spl%valid = .false.
    spl%tension = 1.0
  end subroutine delete_spline

  !> Reset and deallocate variables in passed periodic spline
  subroutine delete_periodic_spline (spl)
    implicit none
    type (periodic_spline), intent (in out) :: spl
    spl%n = 0
    spl%period = 0.0
    if (allocated(spl%x)) deallocate (spl%x,spl%y)
    if (allocated(spl%y2)) deallocate (spl%y2)
    spl%valid = .false.
    spl%tension = 1.0
  end subroutine delete_periodic_spline

  !> Bound wrapper to splint
  real function spline_interp(self, x)
    implicit none
    class (spline), intent(in) :: self
    real, intent(in) :: x
    spline_interp = splint(x, self)
  end function spline_interp

  !> Bound wrapper to dsplint
  real function spline_deriv(self, x)
    implicit none
    class (spline), intent(in) :: self
    real, intent(in) :: x
    spline_deriv = dsplint(x, self)
  end function spline_deriv

  !> FIXME : Add documentation  
  real function splint (x, spl)
    implicit none
    real, intent (in) :: x
    type (spline), intent (in) :: spl
    call handle_spline_error(spl%valid, 'splint')
    splint = fitp_curv2(x, spl%n, spl%x, spl%y, spl%y2, spl%tension)
  end function splint

  !> Bound wrapper to splint
  real function periodic_spline_interp(self, x)
    implicit none
    class (periodic_spline), intent(in) :: self
    real, intent(in) :: x
    periodic_spline_interp = periodic_splint(x, self)
  end function periodic_spline_interp

  !> Bound wrapper to dsplint
  real function periodic_spline_deriv(self, x)
    implicit none
    class (periodic_spline), intent(in) :: self
    real, intent(in) :: x
    periodic_spline_deriv = periodic_dsplint(x, self)
  end function periodic_spline_deriv

  !> FIXME : Add documentation
  real function periodic_splint (x, spl)
    implicit none
    real, intent (in) :: x
    type (periodic_spline), intent (in) :: spl
    call handle_spline_error(spl%valid, 'periodic_splint')
    periodic_splint = fitp_curvp2(x, spl%n, spl%x, spl%y, spl%period, spl%y2, spl%tension)
  end function periodic_splint

  !> FIXME : Add documentation
  real function dsplint (x, spl)
    implicit none
    real, intent (in) :: x
    type (spline), intent (in) :: spl
    call handle_spline_error(spl%valid, 'dsplint')
    dsplint = fitp_curvd(x, spl%n, spl%x, spl%y, spl%y2, spl%tension)
  end function dsplint

  !> FIXME : Add documentation
  real function periodic_dsplint (x, spl)
    implicit none
    real, intent (in) :: x
    type (periodic_spline), intent (in) :: spl
    call handle_spline_error(spl%valid, 'periodic_dsplint')
    periodic_dsplint = fitp_curvpd(x, spl%n, spl%x, spl%y, spl%period, spl%y2, spl%tension)
  end function periodic_dsplint

  !> FIXME : Add documentation
  real function splintint (x0, x1, spl)
    implicit none
    real, intent (in) :: x0, x1
    type (spline), intent (in) :: spl
    call handle_spline_error(spl%valid, 'splintint')
    splintint = fitp_curvi(x0, x1, spl%n, spl%x, spl%y, spl%y2, spl%tension)
  end function splintint

  !> FIXME : Add documentation
  real function periodic_splintint (x0, x1, spl)
    implicit none
    real, intent (in) :: x0, x1
    type (periodic_spline), intent (in) :: spl
    call handle_spline_error(spl%valid, 'periodic_splintint')
    periodic_splintint = fitp_curvpi(x0, x1, spl%n, spl%x, spl%y, spl%period, spl%y2, spl%tension)
  end function periodic_splintint

  !> If not valid abort with error noting invalid spline and which method was invoked
  subroutine handle_spline_error_logical(valid, routine_name)
    use mp, only: mp_abort
    implicit none
    logical, intent(in) :: valid
    character(len = *), intent(in) :: routine_name
    if (.not. valid) call mp_abort('Attempt to use invalid spline in '//routine_name)
  end subroutine handle_spline_error_logical

  !> FIXME : Add documentation
  subroutine inter_d_cspl(r,data,x,dint,ddint)
    implicit none
    real, dimension(:), intent(in) :: r, data
    real, dimension(:), intent(in) :: x
    real, dimension(:), intent(out) :: dint, ddint
    integer :: i, m
    type(spline) :: spl
    real, parameter :: tension = 1.0
    spl = new_spline(r, data, tension)
    m = size(x)
    do i = 1, m
       dint(i) = spl%interpolate(x(i))
       ddint(i)= spl%derivative(x(i))
    end do
  end subroutine inter_d_cspl

  !> FIXME : Add documentation  
  subroutine inter_cspl(r,data,x,dint)
    implicit none
    real, dimension(:), intent(in) :: r, data
    real, dimension(:), intent(in) :: x
    real, dimension(:), intent(out) :: dint
    integer :: i, m
    type(spline) :: spl
    real, parameter :: tension = 1.0
    spl = new_spline(r, data, tension)
    m = size(x)
    do i = 1, m
       dint(i) = spl%interpolate(x(i))
    end do
  end subroutine inter_cspl

! From inet!cs.utexas.edu!cline Tue Oct 31 17:10:31 CST 1989
! Received: from mojave.cs.utexas.edu by cs.utexas.edu (5.59/1.44)
!       id AA29509; Tue, 31 Oct 89 17:11:51 CST
! Posted-Date: Tue, 31 Oct 89 17:10:31 CST
! Message-Id: <8910312310.AA04442@mojave.cs.utexas.edu>
! Received: by mojave.cs.utexas.edu (14.5/1.4-Client)
!       id AA04442; Tue, 31 Oct 89 17:10:34 cst
! Date: Tue, 31 Oct 89 17:10:31 CST
! X-Mailer: Mail User's Shell (6.5 4/17/89)
! From: cline@cs.utexas.edu (Alan Cline)
! To: ehg@research.att.com
! Subject: New FITPACK Subset for netlib
! 
! 
! This new version of FITPACK distributed by netlib is about 20% of 
! the total package in terms of characters, lines of code, and num-
! ber of subprograms. However, these 25 subprograms represent about
! 95% of usages of the package.  What has been omitted are such ca-
! pabilities as:
!   1. Automatic tension determination,
!   2. Derivatives, arclengths, and enclosed areas for planar 
!      curves,
!   3. Three dimensional curves,
!   4. Special surface fitting using equispacing assumptions,
!   5. Surface fitting in annular, wedge, polar, toroidal, lunar,
!      and spherical geometries,
!   6. B-splines in tension generation and usage,
!   7. General surface fitting in three dimensional space.
! 
! (The code previously circulated in netlib is less than 10% of the
! total  package  and is more than a decade old.  Its usage is dis-
! couraged.)
! 
! Please note:  Two versions of the subroutine snhcsh are included.
! Both serve the same purpose:  obtaining approximations to certain
! hyperbolic trigonometric-like functions.  The first is less accu-
! rate (but more efficient) than the second.  Installers should se- 
! lect the one with the precision they desire.
! 
! Interested parties can obtain the entire package on disk or  tape
! from Pleasant  Valley Software, 8603 Altus Cove, Austin TX (USA),
! 78759 at a cost of $495 US. A 340 page manual  is  available  for
! $30  US  per  copy.  The  package  includes  examples and machine
! readable documentation.

  !> FIXME : Add documentation (or tidyup above)
  pure subroutine fitp_curv1 (n,x,y,slp1,slpn,islpsw,yp,temp,sigma,ierr)
    implicit none
    integer, intent(in) :: n, islpsw
    integer, intent(out) :: ierr
    real, dimension(n), intent(in) :: x, y
    real, dimension(n), intent(out) :: yp, temp
    real, intent(in) :: slp1,slpn,sigma
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the parameters necessary to
! compute an interpolatory spline under tension through
! a sequence of functional values. the slopes at the two
! ends of the curve may be specified or omitted.  for actual
! computation of points on the curve it is necessary to call
! the function curv2.
!
! on input--
!
!   n is the number of values to be interpolated (n.ge.2).
!
!   x is an array of the n increasing abscissae of the
!   functional values.
!
!   y is an array of the n ordinates of the values, (i. e.
!   y(k) is the functional value corresponding to x(k) ).
!
!   slp1 and slpn contain the desired values for the first
!   derivative of the curve at x(1) and x(n), respectively.
!   the user may omit values for either or both of these
!   parameters and signal this with islpsw.
!
!   islpsw contains a switch indicating which slope data
!   should be used and which should be estimated by this
!   subroutine,
!          = 0 if slp1 and slpn are to be used,
!          = 1 if slp1 is to be used but not slpn,
!          = 2 if slpn is to be used but not slp1,
!          = 3 if both slp1 and slpn are to be estimated
!              internally.
!
!   yp is an array of length at least n.
!
!   temp is an array of length at least n which is used for
!   scratch storage.
!
! and
!
!   sigma contains the tension factor. this value indicates
!   the curviness desired. if abs(sigma) is nearly zero
!   (e.g. .001) the resulting curve is approximately a
!   cubic spline. if abs(sigma) is large (e.g. 50.) the
!   resulting curve is nearly a polygonal line. if sigma
!   equals zero a cubic spline results.  a standard value
!   for sigma is approximately 1. in absolute value.
!
! on output--
!
!   yp contains the values of the second derivative of the
!   curve at the given nodes.
!
!   ierr contains an error flag,
!        = 0 for normal return,
!        = 1 if n is less than 2,
!        = 2 if x-values are not strictly increasing.
!
! and
!
!   n, x, y, slp1, slpn, islpsw and sigma are unaltered.
!
! this subroutine references package modules ceez, terms,
! and snhcsh.
!
!-----------------------------------------------------------
    integer :: i, ibak, nm1, np1
    real :: sdiag1, diag1, delxnm, dx1, diag, sdiag2, dx2, diag2
    real :: delxn, slpp1, delx1, sigmap, c3, c2, c1, slppn, delx2
    
    nm1 = n-1
    np1 = n+1
    ierr = 0
    if (n .le. 1) go to 8
    if (x(n) .le. x(1)) go to 9
!
! denormalize tension factor
!
    sigmap = abs(sigma) * (n - 1) / (x(n) - x(1))
!
! approximate end slopes
!
    if (islpsw .ge. 2) go to 1
    slpp1 = slp1
    go to 2
1   delx1 = x(2)-x(1)
    delx2 = delx1+delx1
    if (n .gt. 2) delx2 = x(3)-x(1)
    if (delx1 .le. 0. .or. delx2 .le. delx1) go to 9
    call fitp_ceez (delx1,delx2,sigmap,c1,c2,c3,n)
    slpp1 = c1*y(1)+c2*y(2)
    if (n .gt. 2) slpp1 = slpp1+c3*y(3)
2   if (islpsw .eq. 1 .or. islpsw .eq. 3) go to 3
    slppn = slpn
    go to 4
3   delxn = x(n)-x(nm1)
    delxnm = delxn+delxn
    if (n .gt. 2) delxnm = x(n)-x(n-2)
    if (delxn .le. 0. .or. delxnm .le. delxn) go to 9
    call fitp_ceez (-delxn,-delxnm,sigmap,c1,c2,c3,n)
    slppn = c1*y(n)+c2*y(nm1)
    if (n .gt. 2) slppn = slppn+c3*y(n-2)
!
! set up right hand side and tridiagonal system for yp and
! perform forward elimination
!
4   delx1 = x(2)-x(1)
    if (delx1 .le. 0.) go to 9
    dx1 = (y(2)-y(1))/delx1
    call fitp_terms (diag1,sdiag1,sigmap,delx1)
    yp(1) = (dx1-slpp1)/diag1
    temp(1) = sdiag1/diag1
    if (n .eq. 2) go to 6
    do i = 2,nm1
       delx2 = x(i+1)-x(i)
       if (delx2 .le. 0.) go to 9
       dx2 = (y(i+1)-y(i))/delx2
       call fitp_terms (diag2,sdiag2,sigmap,delx2)
       diag = diag1+diag2-sdiag1*temp(i-1)
       yp(i) = (dx2-dx1-sdiag1*yp(i-1))/diag
       temp(i) = sdiag2/diag
       dx1 = dx2
       diag1 = diag2
       sdiag1 = sdiag2
    end do
6   diag = diag1-sdiag1*temp(nm1)
    yp(n) = (slppn-dx1-sdiag1*yp(nm1))/diag
!
! perform back substitution
!
    do i = 2,n
       ibak = np1-i
       yp(ibak) = yp(ibak)-temp(ibak)*yp(ibak+1)
    end do
    return
!
! too few points
!
8   ierr = 1
    return
!
! x-values not strictly increasing
!
9   ierr = 2
    return
  end subroutine fitp_curv1

  !> FIXME : Add documentation
  real function fitp_curv2 (t,n,x,y,yp,sigma)
    implicit none
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: x, y, yp
    real, intent(in) :: t, sigma
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function interpolates a curve at a given point
! using a spline under tension. the subroutine curv1 should
! be called earlier to determine certain necessary
! parameters.
!
! on input--
!
!   t contains a real value to be mapped onto the interpo-
!   lating curve.
!
!   n contains the number of points which were specified to
!   determine the curve.
!
!   x and y are arrays containing the abscissae and
!   ordinates, respectively, of the specified points.
!
!   yp is an array of second derivative values of the curve
!   at the nodes.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, yp, and sigma should be input
! unaltered from the output of curv1.
!
! on output--
!
!   curv2 contains the interpolated value.
!
! none of the input parameters are altered.
!
! this function references package modules intrvl and
! snhcsh.
!
!-----------------------------------------------------------
    integer :: i, im1
    real :: ss, sigdel, s1, s2, sum, sigmap
    real :: del1, del2, dels
!
! determine interval
!
    im1 = fitp_intrvl(t,x,n)
    i = im1+1
!
! denormalize tension factor
!
    sigmap = abs(sigma) * (n - 1) / (x(n) - x(1))
!
! set up and perform interpolation
!
    del1 = t-x(im1)
    del2 = x(i)-t
    dels = x(i)-x(im1)
    sum = (y(i)*del1+y(im1)*del2)/dels
    if (is_zero(sigmap)) then
       fitp_curv2 = sum-del1*del2*(yp(i)*(del1+dels)+yp(im1)*(del2+dels))/(6.*dels)
    else
       sigdel = sigmap*dels
       ss = sinhm_fun(sigdel)
       s1 = sinhm_fun(sigmap * del1)
       s2 = sinhm_fun(sigmap * del2)
       fitp_curv2 = sum+(yp(i)*del1*(s1-ss)+yp(im1)*del2*(s2-ss))/(sigdel*sigmap*(1.+ss))
    end if
  end function fitp_curv2

  !> FIXME : Add documentation  
  real function fitp_curvd (t,n,x,y,yp,sigma)
    implicit none
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: x, y, yp
    real, intent(in) :: t, sigma
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function differentiates a curve at a given point
! using a spline under tension. the subroutine curv1 should
! be called earlier to determine certain necessary
! parameters.
!
! on input--
!
!   t contains a real value at which the derivative is to be
!   determined.
!
!   n contains the number of points which were specified to
!   determine the curve.
!
!   x and y are arrays containing the abscissae and
!   ordinates, respectively, of the specified points.
!
!   yp is an array of second derivative values of the curve
!   at the nodes.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, yp, and sigma should be input
! unaltered from the output of curv1.
!
! on output--
!
!   curvd contains the derivative value.
!
! none of the input parameters are altered.
!
! this function references package modules intrvl and
! snhcsh.
!
!-----------------------------------------------------------
    integer :: i, im1
    real :: ss, sigdel, c1, c2, sum, sigmap
    real :: del1, del2, dels
!
! determine interval
!
    im1 = fitp_intrvl(t,x,n)
    i = im1+1
!
! denormalize tension factor
!
    sigmap = abs(sigma) * (n - 1) / (x(n) - x(1))
!
! set up and perform differentiation
!
    del1 = t-x(im1)
    del2 = x(i)-t
    dels = x(i)-x(im1)
    sum = (y(i)-y(im1))/dels
    if (is_zero(sigmap)) then
       fitp_curvd = sum+(yp(i)*(2.*del1*del1-del2*(del1+dels))- &
            yp(im1)*(2.*del2*del2-del1*(del2+dels))) &
            /(6.*dels)
    else
       sigdel = sigmap*dels
       ss = sinhm_fun(sigdel)
       c1 = coshm_fun(sigmap * del1)
       c2 = coshm_fun(sigmap * del2)
       fitp_curvd = sum+(yp(i)*(c1-ss)-yp(im1)*(c2-ss))/(sigdel*sigmap*(1.+ss))
    end if
  end function fitp_curvd

  !> FIXME : Add documentation  
  real function fitp_curvi (xl,xu,n,x,y,yp,sigma)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: xl, xu, sigma
    real, dimension(n), intent(in) :: x, y, yp
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function integrates a curve specified by a spline
! under tension between two given limits. the subroutine
! curv1 should be called earlier to determine necessary
! parameters.
!
! on input--
!
!   xl and xu contain the upper and lower limits of inte-
!   gration, respectively. (sl need not be less than or
!   equal to xu, curvi (xl,xu,...) .eq. -curvi (xu,xl,...) ).
!
!   n contains the number of points which were specified to
!   determine the curve.
!
!   x and y are arrays containing the abscissae and
!   ordinates, respectively, of the specified points.
!
!   yp is an array from subroutine curv1 containing
!   the values of the second derivatives at the nodes.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, yp, and sigma should be input
! unaltered from the output of curv1.
!
! on output--
!
!   curvi contains the integral value.
!
! none of the input parameters are altered.
!
! this function references package modules intrvl and
! snhcsh.
!
!-----------------------------------------------------------
    integer :: i, ilp1, ilm1, il, ium1, iu
    real :: delu1, delu2, c2, ss, cs, cu2, cl1, cl2, cu1
    real :: dell1, dell2, deli, c1, ssign, sigmap
    real :: xxl, xxu, t1, t2, dels, sum, del1, del2

    ! denormalize tension factor
    sigmap = abs(sigma) * (n - 1) / (x(n) - x(1))

    ! determine actual upper and lower bounds
    if (xl < xu) then
       xxl = xl
       xxu = xu
       ssign = 1.
    else if (xl > xu) then
       xxl = xu
       xxu = xl
       ssign = -1.
    else
       ! return zero if xl .eq. xu
       fitp_curvi = 0.
       return
    end if

    ! search for proper intervals
    ilm1 = fitp_intrvl (xxl,x,n)
    il = ilm1+1
    ium1 = fitp_intrvl (xxu,x,n)
    iu = ium1+1
    if (il == iu) then
       ! integrate from xxl to xxu
       delu1 = xxu-x(ium1)
       delu2 = x(iu)-xxu
       dell1 = xxl-x(ium1)
       dell2 = x(iu)-xxl
       dels = x(iu)-x(ium1)
       deli = xxu-xxl
       t1 = (delu1+dell1)*deli/(2.*dels)
       t2 = (delu2+dell2)*deli/(2.*dels)
       sum = t1*y(iu)+t2*y(ium1)
       if (is_not_zero(sigma)) then
          cu1 = coshmm_fun(sigmap*delu1)
          cu2 = coshmm_fun(sigmap*delu2)
          cl1 = coshmm_fun(sigmap*dell1)
          cl2 = coshmm_fun(sigmap*dell2)
          ss = sinhm_fun(sigmap*dels)

          sum = sum+(yp(iu)*(delu1*delu1*(cu1-ss/2.) &
               -dell1*dell1*(cl1-ss/2.)) &
               +yp(ium1)*(dell2*dell2*(cl2-ss/2.) &
               -delu2*delu2*(cu2-ss/2.)))/ &
               (sigmap*sigmap*dels*(1.+ss))
       else
          sum = sum-t1*(delu2*(dels+delu1)+dell2*(dels+dell1))* &
               yp(iu)/12. &
               -t2*(dell1*(dels+dell2)+delu1*(dels+delu2))* &
               yp(ium1)/12.
       end if
    else
       ! integrate from xxl to x(il)
       sum = 0.
       if (not_exactly_equal(xxl, x(il))) then
          del1 = xxl-x(ilm1)
          del2 = x(il)-xxl
          dels = x(il)-x(ilm1)
          t1 = (del1+dels)*del2/(2.*dels)
          t2 = del2*del2/(2.*dels)
          sum = t1*y(il)+t2*y(ilm1)
          if (is_not_zero(sigma)) then
             c1 = coshmm_fun(sigmap*del1)
             c2 = coshmm_fun(sigmap*del2)
             cs = coshmm_fun(sigmap*dels)
             ss = sinhm_fun(sigmap*dels)
             sum = sum+((dels*dels*(cs-ss/2.)-del1*del1*(c1-ss/2.)) &
                  *yp(il)+del2*del2*(c2-ss/2.)*yp(ilm1))/ &
                  (sigmap*sigmap*dels*(1.+ss))
          else
             sum = sum-t1*t1*dels*yp(il)/6. &
                  -t2*(del1*(del2+dels)+dels*dels)*yp(ilm1)/12.
          end if
       end if

       ! integrate over interior intervals
       if (iu-il /= 1) then
          ilp1 = il+1
          do i = ilp1,ium1
             dels = x(i)-x(i-1)
             sum = sum+(y(i)+y(i-1))*dels/2.
             if (is_not_zero(sigma)) then
                cs = coshmm_fun(sigmap*dels)
                ss = sinhm_fun(sigmap*dels)
                sum = sum+(yp(i)+yp(i-1))*dels*(cs-ss/2.)/(sigmap*sigmap*(1.+ss))
             else
                sum = sum-(yp(i)+yp(i-1))*dels*dels*dels/24.
             end if
          end do
       end if

       ! integrate from x(iu-1) to xxu
       if (not_exactly_equal(xxu, x(ium1))) then
          del1 = xxu-x(ium1)
          del2 = x(iu)-xxu
          dels = x(iu)-x(ium1)
          t1 = del1*del1/(2.*dels)
          t2 = (del2+dels)*del1/(2.*dels)
          sum = sum+t1*y(iu)+t2*y(ium1)
          if (is_not_zero(sigma)) then
             c1 = coshmm_fun(sigmap*del1)
             c2 = coshmm_fun(sigmap*del2)
             cs = coshmm_fun(sigmap*dels)
             ss = sinhm_fun(sigmap*dels)
             sum = sum+(yp(iu)*del1*del1*(c1-ss/2.)+yp(ium1)* &
                  (dels*dels*(cs-ss/2.)-del2*del2*(c2-ss/2.))) &
                  /(sigmap*sigmap*dels*(1.+ss))
          else
             sum = sum-t1*(del2*(del1+dels)+dels*dels)*yp(iu)/12.-t2*t2*dels*yp(ium1)/6.
          end if
       end if
    end if

    ! correct sign and return
    fitp_curvi = ssign*sum
  end function fitp_curvi

  !> FIXME : Add documentation
  pure subroutine fitp_curvp1 (n,x,y,p,yp,temp,sigma,ierr)
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: ierr
    real, intent(in) :: sigma, p
    real, dimension(:), intent(in) :: x, y
    real, dimension(:), intent(out) :: yp, temp
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the parameters necessary to
! compute a periodic interpolatory spline under tension
! through a sequence of functional values. for actual ends
! of the curve may be specified or omitted.  for actual
! computation of points on the curve it is necessary to call
! the function curvp2.
!
! on input--
!
!   n is the number of values to be interpolated (n.ge.2).
!
!   x is an array of the n increasing abscissae of the
!   functional values.
!
!   y is an array of the n ordinates of the values, (i. e.
!   y(k) is the functional value corresponding to x(k) ).
!
!   p is the period (p .gt. x(n)-x(1)).
!
!   yp is an array of length at least n.
!
!   temp is an array of length at least 2*n which is used
!   for scratch storage.
!
! and
!
!   sigma contains the tension factor.  this value indicates
!   the curviness desired. if abs(sigma) is nearly zero
!   (e.g. .001) the resulting curve is approximately a
!   cubic spline. if abs(sigma) is large (e.g. 50.) the
!   resulting curve is nearly a polygonal line. if sigma
!   equals zero a cubic spline results.  a standard value
!   for sigma is approximately 1. in absolute value.
!
! on output--
!
!   yp contains the values of the second derivative of the
!   curve at the given nodes.
!
!   ierr contains an error flag,
!        = 0 for normal return,
!        = 1 if n is less than 2,
!        = 2 if p is less than or equal to x(n)-x(1),
!        = 3 if x-values are not strictly increasing.
!
! and
!
!  n, x, y, and sigma are unaltered.
!
! this subroutine references package modules terms and
! snhcsh.
!
!-----------------------------------------------------------
    integer :: i, npibak, npi, ibak, nm1, np1
    real :: diag, diag2, sdiag2, ypn, dx2, sigmap, delx1
    real :: sdiag1, delx2, dx1, diag1

    nm1 = n-1
    np1 = n+1
    ierr = 0
    if (n .le. 1) go to 6
    if (p .le. x(n)-x(1) .or. p .le. 0.) go to 7
!
! denormalize tension factor
!
    sigmap = abs(sigma) * n / p
!
! set up right hand side and tridiagonal system for yp and
! perform forward elimination
!
    delx1 = p-(x(n)-x(1))
    dx1 = (y(1)-y(n))/delx1
    call fitp_terms (diag1,sdiag1,sigmap,delx1)
    delx2 = x(2)-x(1)
    if (delx2 .le. 0.) go to 8
    dx2 = (y(2)-y(1))/delx2
    call fitp_terms (diag2,sdiag2,sigmap,delx2)
    diag = diag1+diag2
    yp(1) = (dx2-dx1)/diag
    temp(np1) = -sdiag1/diag
    temp(1) = sdiag2/diag
    dx1 = dx2
    diag1 = diag2
    sdiag1 = sdiag2
    if (n .eq. 2) go to 2
    do i = 2,nm1
       npi = n+i
       delx2 = x(i+1)-x(i)
       if (delx2 .le. 0.) go to 8
       dx2 = (y(i+1)-y(i))/delx2
       call fitp_terms (diag2,sdiag2,sigmap,delx2)
       diag = diag1+diag2-sdiag1*temp(i-1)
       yp(i) = (dx2-dx1-sdiag1*yp(i-1))/diag
       temp(npi) = -temp(npi-1)*sdiag1/diag
       temp(i) = sdiag2/diag
       dx1 = dx2
       diag1 = diag2
       sdiag1 = sdiag2
    end do
2   delx2 = p-(x(n)-x(1))
    dx2 = (y(1)-y(n))/delx2
    call fitp_terms (diag2,sdiag2,sigmap,delx2)
    yp(n) = dx2-dx1
    temp(nm1) = temp(2*n-1)-temp(nm1)
    if (n .eq. 2) go to 4
!
! perform first step of back substitution
    !
    do i = 3,n
       ibak = np1-i
       npibak =n+ibak
       yp(ibak) = yp(ibak)-temp(ibak)*yp(ibak+1)
       temp(ibak) =temp(npibak)-temp(ibak)*temp(ibak+1)
    end do
4   yp(n) = (yp(n)-sdiag2*yp(1)-sdiag1*yp(nm1))/ &
         (diag1+diag2+sdiag2*temp(1)+sdiag1*temp(nm1))
!
! perform second step of back substitution
!
    ypn =   yp(n)
    do i = 1,nm1
       yp(i) = yp(i)+temp(i)*ypn
    end do
    return
!
! too few points
!
6   ierr = 1
    return
!
! period too small
!
7   ierr = 2
    return
!
! x-values not strictly increasing
!
8   ierr = 3
    return
  end subroutine fitp_curvp1

  !> FIXME : Add documentation  
  real function fitp_curvp2 (t,n,x,y,p,yp,sigma)
    implicit none
    real, intent(in) :: t, p, sigma
    integer, intent(in) :: n
    real, dimension(:), intent(in) :: x, y, yp
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function interpolates a curve at a given point using
! a periodic spline under tension. the subroutine curvp1
! should be called earlier to determine certain necessary
! parameters.
!
! on input--
!
!   t contains a real value to be mapped onto the interpo-
!   lating curve.
!
!   n contains the number of points which were specified to
!   determine the curve.
!
!   x and y are arrays containing the abscissae and
!   ordinates, respectively, of the specified points.
!
!   p contains the period.
!
!   yp is an array of second derivative values of the curve
!   at the nodes.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, p, yp, and sigma should be input
! unaltered from the output of curvp1.
!
! on output--
!
!   curvp2 contains the interpolated value.
!
! none of the input parameters are altered.
!
! this function references package modules intrvp and
! snhcsh.
!
!-----------------------------------------------------------
    integer :: i, im1
    real :: ss, sigdel, sum, s2, s1
    real :: tp, sigmap, dels, del2, del1

!
! determine interval
!
    tp = fitp_periodic_wrap_value(t, x(1), p)
    im1 = fitp_intrvp (tp, x, n)
    i = im1+1
!
! denormalize tension factor
!
    sigmap = abs(sigma) * n / p
!
! set up and perform interpolation
!
    del1 = tp-x(im1)
    if (im1 == n) then
       i = 1
       del2 = x(1)+p-tp
       dels = p-(x(n)-x(1))
    else
       del2 = x(i)-tp
       dels = x(i)-x(im1)
    end if
    sum = (y(i)*del1+y(im1)*del2)/dels
    if (is_zero(sigmap)) then
       fitp_curvp2 = sum-del1*del2*(yp(i)*(del1+dels)+yp(im1)*(del2+dels))/(6.*dels)
    else
       sigdel = sigmap*dels
       ss = sinhm_fun(sigdel)
       s1 = sinhm_fun(sigmap * del1)
       s2 = sinhm_fun(sigmap * del2)
       fitp_curvp2 = sum+(yp(i)*del1*(s1-ss)+yp(im1)*del2*(s2-ss))/(sigdel*sigmap*(1.+ss))
    end if
  end function fitp_curvp2

  !> FIXME : Add documentation  
  real function fitp_curvpd (t,n,x,y,p,yp,sigma)
    real, intent(in) :: t, p, sigma
    integer, intent(in) :: n
    real, dimension(:), intent(in) :: x, y, yp
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function is the derivative of curvp2
! interpolates a curve at a given point using
! a periodic spline under tension. the subroutine curvp1
! should be called earlier to determine certain necessary
! parameters.
!
! on input--
!
!   t contains a real value to be mapped onto the interpo-
!   lating curve.
!
!   n contains the number of points which were specified to
!   determine the curve.
!
!   x and y are arrays containing the abscissae and
!   ordinates, respectively, of the specified points.
!
!   p contains the period.
!
!   yp is an array of second derivative values of the curve
!   at the nodes.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, p, yp, and sigma should be input
! unaltered from the output of curvp1.
!
! on output--
!
!   curvpd contains the interpolated derivative
!
! none of the input parameters are altered.
!
! this function references package modules intrvp and
! snhcsh.
!
!-----------------------------------------------------------
    integer :: i, im1
    real :: ss, sigdel, sum, c2, c1
    real :: tp, sigmap, dels, del2, del1
!
! determine interval
!
    tp = fitp_periodic_wrap_value(t, x(1), p)
    im1 = fitp_intrvp (tp, x, n)
    i = im1+1
!
! denormalize tension factor
!
    sigmap = abs(sigma) * n / p
!
! set up and perform interpolation
!
    del1 = tp-x(im1)
    if (im1 .eq. n) then
       i = 1
       del2 = x(1)+p-tp
       dels = p-(x(n)-x(1))
    else
       del2 = x(i)-tp
       dels = x(i)-x(im1)
    end if
    ! Here on identical to fitp_curvd
    sum = (y(i)-y(im1))/dels
    if (is_zero(sigmap)) then
       fitp_curvpd = sum+(yp(i)*(2.*del1*del1-del2*(del1+dels))- &
            yp(im1)*(2.*del2*del2-del1*(del2+dels))) &
            /(6.*dels)
    else
       sigdel = sigmap*dels
       ss = sinhm_fun(sigdel)
       c1 = coshm_fun(sigmap * del1)
       c2 = coshm_fun(sigmap * del2)
       fitp_curvpd = sum+(yp(i)*(c1-ss)-yp(im1)*(c2-ss))/(sigdel*sigmap*(1.+ss))
    end if
  end function fitp_curvpd

  !> FIXME : Add documentation
  real function fitp_curvpi (xl,xu,n,x,y,p,yp,sigma)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: xl, xu, p, sigma
    real, dimension(n), intent(in) :: x, y, yp
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function integrates a curve specified by a periodic
! spline under tension between two given limits. the
! subroutine curvp1 should be called earlier to determine
! necessary parameters.
!
! on input--
!
!   xl and xu contain the upper and lower limits of inte-
!   gration, respectively. (sl need not be less than or
!   equal to xu, curvpi (xl,xu,...) .eq. -curvpi (xu,xl,...) ).
!
!   n contains the number of points which were specified to
!   determine the curve.
!
!   x and y are arrays containing the abscissae and
!   ordinates, respectively, of the specified points.
!
!   p contains the period.
!
!   yp is an array from subroutine curvp1 containing
!   the values of the second derivatives at the nodes.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, p, yp, and sigma should be input
! unaltered from the output of curvp1.
!
! on output--
!
!
!   curvpi contains the integral value.
!
! none of the input parameters are altered.
!
! this function references package modules intrvp and
! snhcsh.
!
!--------------------------------------------------------------
    integer :: np1, im1, ii, iup1, ilp1, ideltp
    real :: s7, s6, s3, c2, c1, s5, s4, cl1, cu2, cu1, si, so
    real :: cl2, delu2, delu1, s8, deli, dell2, dell1
    integer :: ium1, il, isave, isign, lper, ilm1, iu, i
    real :: xxu, xsave, x1pp, sigmap, xxl, xil, del1, s2, cs
    real :: t2, t1, del2, s1, xiu, ss, dels

    integer :: uper
    logical :: bdy
!
! denormalize tension factor
!
    sigmap = abs(sigma) * n / p
!
! determine actual upper and lower bounds
!
    x1pp = x(1)+p
    isign = 1
    xxl = fitp_periodic_wrap_value(xl, x(1), p)
    ilm1 = fitp_intrvp (xxl, x, n)
    lper = int((xl-x(1))/p)
    if (xl .lt. x(1)) lper = lper-1
    xxu = fitp_periodic_wrap_value(xu, x(1), p)
    ium1 = fitp_intrvp (xxu, x, n)
    uper = int((xu-x(1))/p)
    if (xu .lt. x(1)) uper = uper-1
    ideltp = uper-lper
    bdy = ideltp * (xxu-xxl) .lt. 0.
    if ((ideltp .eq. 0 .and. xxu .lt. xxl) .or. ideltp .lt. 0) isign = -1
    if (bdy) ideltp = ideltp-isign
    if (xxu .ge. xxl) go to 1
    xsave = xxl
    xxl = xxu
    xxu = xsave
    isave = ilm1
    ilm1 = ium1
    ium1 = isave
1   il = ilm1+1
    if (ilm1 .eq. n) il = 1
    xil = x(il)
    if (ilm1 .eq. n) xil = x1pp
    iu = ium1+1
    if (ium1 .eq. n) iu = 1
    xiu = x(iu)
    if (ium1 .eq. n) xiu = x1pp
    s1 = 0.
    if (ilm1 .eq. 1 .or. (ideltp .eq. 0 .and..not. bdy)) go to 4
!
! integrate from x(1) to x(ilm1), store in s1
!
    do i = 2,ilm1
       dels = x(i)-x(i-1)
       s1 = s1+(y(i)+y(i-1))*dels/2.
       if (is_not_zero(sigma)) then
          ss = sinhm_fun(sigmap * dels)
          cs = coshmm_fun(sigmap * dels)
          s1 = s1+(yp(i)+yp(i-1))*dels*(cs-ss/2.)/(sigmap*sigmap*(1.+ss))
       else
          s1 = s1-(yp(i)+yp(i-1))*dels*dels*dels/24.
       end if
    end do
4   s2 = 0.
    if (x(ilm1) .ge. xxl .or. (ideltp .eq. 0 .and. .not. bdy)) go to 6
!
! integrate from x(ilm1) to xxl, store in s2
!
    del1 = xxl-x(ilm1)
    del2 = xil-xxl
    dels = xil-x(ilm1)
    t1 = del1*del1/(2.*dels)
    t2 = (del2+dels)*del1/(2.*dels)
    s2 = t1*y(il)+t2*y(ilm1)
    if (is_not_zero(sigma)) then
       c1 = coshmm_fun(sigmap * del1)
       c2 = coshmm_fun(sigmap * del2)
       ss = sinhm_fun(sigmap * dels)
       cs = coshmm_fun(sigmap * dels)
       s2 = s2+(yp(il)*del1*del1*(c1-ss/2.)+yp(ilm1)* &
            (dels*dels*(cs-ss/2.)-del2*del2*(c2-ss/2.))) &
            /(sigmap*sigmap*dels*(1.+ss))
    else
       s2 = s2-t1*(del2*(del1+dels) &
            +dels*dels)*yp(il)/12. &
            -t2*t2*dels*yp(ilm1)/6.
    end if

6   s3 = 0.
    if (xxl .ge. xil .or. (ideltp .eq. 0 .and. bdy).or. ilm1 .eq. ium1) go to 8


!
! integrate from xxl to xil, store in s3
!
    del1 = xxl-x(ilm1)
    del2 = xil-xxl
    dels = xil-x(ilm1)
    t1 = (del1+dels)*del2/(2.*dels)
    t2 = del2*del2/(2.*dels)
    s3 = t1*y(il)+t2*y(ilm1)
    if (is_not_zero(sigma)) then
       c1 = coshmm_fun(sigmap * del1)
       c2 = coshmm_fun(sigmap * del2)
       ss = sinhm_fun(sigmap * dels)
       cs = coshmm_fun(sigmap * dels)
       s3 = s3+((dels*dels*(cs-ss/2.)-del1*del1*(c1-ss/2.)) &
            *yp(il)+del2*del2*(c2-ss/2.)*yp(ilm1))/ &
            (sigmap*sigmap*dels*(1.+ss))
    else
       s3 = s3-t1*t1*dels*yp(il)/6. &
            -t2*(del1*(del2+dels)+dels*dels)* &
            yp(ilm1)/12.
    end if
8   s4 = 0.
    if (ilm1 .ge. ium1-1 .or. (ideltp .eq. 0 .and. bdy)) go to 11
!
! integrate from xil to x(ium1), store in s4
!
    ilp1 = il+1
    do i = ilp1,ium1
       dels = x(i)-x(i-1)
       s4 = s4+(y(i)+y(i-1))*dels/2.
       if (is_not_zero(sigma)) then
          ss = sinhm_fun(sigmap * dels)
          cs = coshmm_fun(sigmap * dels)
          s4 = s4+(yp(i)+yp(i-1))*dels*(cs-ss/2.)/(sigmap*sigmap*(1.+ss))
       else
          s4 = s4-(yp(i)+yp(i-1))*dels*dels*dels/24.
       end if
    end do
11  s5 = 0.
    if (x(ium1) .ge. xxu .or. (ideltp .eq. 0 .and. bdy) .or. ilm1 .eq. ium1) go to 13
!
! integrate from x(ium1) to xxu, store in s5
!
    del1 = xxu-x(ium1)
    del2 = xiu-xxu
    dels = xiu-x(ium1)
    t1 = del1*del1/(2.*dels)
    t2 = (del2+dels)*del1/(2.*dels)
    s5 = t1*y(iu)+t2*y(ium1)
    if (is_not_zero(sigma)) then
       c1 = coshmm_fun(sigmap * del1)
       c2 = coshmm_fun(sigmap * del2)
       ss = sinhm_fun(sigmap * dels)
       cs = coshmm_fun(sigmap * dels)
       s5 = s5+(yp(iu)*del1*del1*(c1-ss/2.)+yp(ium1)* &
            (dels*dels*(cs-ss/2.)-del2*del2*(c2-ss/2.))) &
            /(sigmap*sigmap*dels*(1.+ss))
    else
       s5 = s5-t1*(del2*(del1+dels)+dels*dels)*yp(iu)/12.-t2*t2*dels*yp(ium1)/6.
    end if
13  s6 = 0.
    if (xxu .ge. xiu .or. (ideltp .eq. 0 .and. .not. bdy)) go to 15
!
! integrate from xxu to xiu, store in s6
!
    del1 = xxu-x(ium1)
    del2 = xiu-xxu
    dels = xiu-x(ium1)
    t1 = (del1+dels)*del2/(2.*dels)
    t2 = del2*del2/(2.*dels)
    s6 = t1*y(iu)+t2*y(ium1)
    if (is_not_zero(sigma)) then
       c1 = coshmm_fun(sigmap * del1)
       c2 = coshmm_fun(sigmap * del2)
       ss = sinhm_fun(sigmap * dels)
       cs = coshmm_fun(sigmap * dels)
       s6 = s6+((dels*dels*(cs-ss/2.)-del1*del1*(c1-ss/2.)) &
            *yp(iu)+del2*del2*(c2-ss/2.)*yp(ium1))/ &
            (sigmap*sigmap*dels*(1.+ss))
    else
       s6 = s6-t1*t1*dels*yp(iu)/6.-t2*(del1*(del2+dels)+dels*dels)*yp(ium1)/12.
    end if
15  s7 = 0.
    if (iu .eq. 1 .or. (ideltp .eq. 0 .and. .not. bdy)) go to 18
!
! integrate from xiu to x1pp, store in s7
!
    np1 = n+1
    iup1 = iu+1
    do ii = iup1,np1
       im1 = ii-1
       i = ii
       if (i .eq. np1) i=1
       dels = x(i)-x(im1)
       if (dels .le. 0.) dels=dels+p
       s7 = s7+(y(i)+y(im1))*dels/2.
       if (is_not_zero(sigma)) then
          ss = sinhm_fun(sigmap * dels)
          cs = coshmm_fun(sigmap * dels)
          s7 = s7+(yp(i)+yp(im1))*dels*(cs-ss/2.)/(sigmap*sigmap*(1.+ss))
       else
          s7 = s7-(yp(i)+yp(im1))*dels*dels*dels/24.
       end if
    end do
18  s8 = 0.
    if (ilm1 .lt. ium1 .or. (ideltp .eq. 0 .and. bdy))go to 20

!
! integrate from xxl to xxu, store in s8
!
    delu1 = xxu-x(ium1)
    delu2 = xiu-xxu
    dell1 = xxl-x(ium1)
    dell2 = xiu-xxl
    dels = xiu-x(ium1)
    deli = xxu-xxl
    t1 = (delu1+dell1)*deli/(2.*dels)
    t2 = (delu2+dell2)*deli/(2.*dels)
    s8 = t1*y(iu)+t2*y(ium1)
    if (is_not_zero(sigma)) then
       cu1 = coshmm_fun(sigmap * delu1)
       cu2 = coshmm_fun(sigmap * delu2)
       cl1 = coshmm_fun(sigmap * dell1)
       cl2 = coshmm_fun(sigmap * dell2)
       ss = sinhm_fun(sigmap * dels)
       s8 = s8+(yp(iu)*(delu1*delu1*(cu1-ss/2.) &
            -dell1*dell1*(cl1-ss/2.)) &
            +yp(ium1)*(dell2*dell2*(cl2-ss/2.) &
            -delu2*delu2*(cu2-ss/2.)))/ &
            (sigmap*sigmap*dels*(1.+ss))
    else
       s8 = s8-t1*(delu2*(dels+delu1) &
            +dell2*(dels+dell1))*yp(iu)/12. &
            -t2*(dell1*(dels+dell2) &
            +delu1*(dels+delu2))*yp(ium1)/12.
    end if
20  so = s1+s2+s6+s7
    si = s3+s4+s5+s8
    if (bdy) then
      fitp_curvpi = ideltp * (so+si) + isign * so
    else
      fitp_curvpi = ideltp * (so+si) + isign * si
    end if
  end function fitp_curvpi

  !> FIXME : Add documentation
  pure subroutine fitp_surf1 (m,n,x,y,z,iz,zx1,zxm,zy1,zyn,zxy11, &
       zxym1,zxy1n,zxymn,islpsw,zp,temp,sigma,ierr)
    implicit none
    integer, intent(in) :: m, n, iz, islpsw
    real, dimension(m), intent(in) :: x, zy1, zyn
    real, dimension(n), intent(in) :: y, zx1, zxm
    real, dimension(iz,n), intent(in) :: z
    real, intent(in) :: zxy11, zxym1, zxy1n, zxymn, sigma
    real, dimension(m,n,3), intent(out) :: zp
    real, dimension(n+n+m), intent(out) :: temp
    integer, intent(out) :: ierr
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the parameters necessary to
! compute an interpolatory surface passing through a rect-
! angular grid of functional values. the surface determined
! can be represented as the tensor product of splines under
! tension. the x- and y-partial derivatives around the
! boundary and the x-y-partial derivatives at the four
! corners may be specified or omitted. for actual mapping
! of points onto the surface it is necessary to call the
! function surf2.
!
! on input--
!
!   m is the number of grid lines in the x-direction, i. e.
!   lines parallel to the y-axis (m .ge. 2).
!
!   n is the number of grid lines in the y-direction, i. e.
!   lines parallel to the x-axis (n .ge. 2).
!
!   x is an array of the m x-coordinates of the grid lines
!   in the x-direction. these should be strictly increasing.
!
!   y is an array of the n y-coordinates of the grid lines
!   in the y-direction. these should be strictly increasing.
!
!   z is an array of the m * n functional values at the grid
!   points, i. e. z(i,j) contains the functional value at
!   (x(i),y(j)) for i = 1,...,m and j = 1,...,n.
!
!   iz is the row dimension of the matrix z used in the
!   calling program (iz .ge. m).
!
!   zx1 and zxm are arrays of the m x-partial derivatives
!   of the function along the x(1) and x(m) grid lines,
!   respectively. thus zx1(j) and zxm(j) contain the x-part-
!   ial derivatives at the points (x(1),y(j)) and
!   (x(m),y(j)), respectively, for j = 1,...,n. either of
!   these parameters will be ignored (and approximations
!   supplied internally) if islpsw so indicates.
!
!   zy1 and zyn are arrays of the n y-partial derivatives
!   of the function along the y(1) and y(n) grid lines,
!   respectively. thus zy1(i) and zyn(i) contain the y-part-
!   ial derivatives at the points (x(i),y(1)) and
!   (x(i),y(n)), respectively, for i = 1,...,m. either of
!   these parameters will be ignored (and estimations
!   supplied internally) if islpsw so indicates.
!
!   zxy11, zxym1, zxy1n, and zxymn are the x-y-partial
!   derivatives of the function at the four corners,
!   (x(1),y(1)), (x(m),y(1)), (x(1),y(n)), and (x(m),y(n)),
!   respectively. any of the parameters will be ignored (and
!   estimations supplied internally) if islpsw so indicates.
!
!   islpsw contains a switch indicating which boundary
!   derivative information is user-supplied and which
!   should be estimated by this subroutine. to determine
!   islpsw, let
!        i1 = 0 if zx1 is user-supplied (and = 1 otherwise),
!        i2 = 0 if zxm is user-supplied (and = 1 otherwise),
!        i3 = 0 if zy1 is user-supplied (and = 1 otherwise),
!        i4 = 0 if zyn is user-supplied (and = 1 otherwise),
!        i5 = 0 if zxy11 is user-supplied
!                                       (and = 1 otherwise),
!        i6 = 0 if zxym1 is user-supplied
!                                       (and = 1 otherwise),
!        i7 = 0 if zxy1n is user-supplied
!                                       (and = 1 otherwise),
!        i8 = 0 if zxymn is user-supplied
!                                       (and = 1 otherwise),
!   then islpsw = i1 + 2*i2 + 4*i3 + 8*i4 + 16*i5 + 32*i6
!                   + 64*i7 + 128*i8
!   thus islpsw = 0 indicates all derivative information is
!   user-supplied and islpsw = 255 indicates no derivative
!   information is user-supplied. any value between these
!   limits is valid.
!
!   zp is an array of at least 3*m*n locations.
!
!   temp is an array of at least n+n+m locations which is
!   used for scratch storage.
!
! and
!
!   sigma contains the tension factor. this value indicates
!   the curviness desired. if abs(sigma) is nearly zero
!   (e. g. .001) the resulting surface is approximately the
!   tensor product of cubic splines. if abs(sigma) is large
!   (e. g. 50.) the resulting surface is approximately
!   bi-linear. if sigma equals zero tensor products of
!   cubic splines result. a standard value for sigma is
!   approximately 1. in absolute value.
!
! on output--
!
!   zp contains the values of the xx-, yy-, and xxyy-partial
!   derivatives of the surface at the given nodes.
!
!   ierr contains an error flag,
!        = 0 for normal return,
!        = 1 if n is less than 2 or m is less than 2,
!        = 2 if the x-values or y-values are not strictly
!            increasing.
!
! and
!
!   m, n, x, y, z, iz, zx1, zxm, zy1, zyn, zxy11, zxym1,
!   zxy1n, zxymn, islpsw, and sigma are unaltered.
!
! this subroutine references package modules ceez, terms,
! and snhcsh.
!
!-----------------------------------------------------------
    integer :: jbak, jbakp1, jm1, jp1, im1, ibakp1
    integer :: npibak, ip1, ibak, nm1, np1, mm1, mp1, npm, j
    integer :: npmpj, i, npi
    real :: diagi, sdiag1, del2, zxymns, delxmm, del1
    real :: diag1, deli, diag2, diagin, sdiag2, t
    real :: delxm, sigmay, dely1, c1, dely2, c2
    real :: delx1, delx2, zxy1ns, c3, delyn, sigmax, delynm

    mm1 = m-1
    mp1 = m+1
    nm1 = n-1
    np1 = n+1
    npm = n+m
    ierr = 0
    if (n .le. 1 .or. m .le. 1) go to 46
    if (y(n) .le. y(1)) go to 47
!
! denormalize tension factor in y-direction
!
    sigmay = abs(sigma) * (n - 1) / (y(n) - y(1))
!
! obtain y-partial derivatives along y = y(1)
!
    if ((islpsw/8)*2 .ne. (islpsw/4)) go to 2
    do i = 1,m
       zp(i,1,1) = zy1(i)
    end do
    go to 5
2   dely1 = y(2)-y(1)
    dely2 = dely1+dely1
    if (n .gt. 2) dely2 = y(3)-y(1)
    if (dely1 .le. 0. .or. dely2 .le. dely1) go to 47
    call fitp_ceez (dely1,dely2,sigmay,c1,c2,c3,n)
    do i = 1,m
       zp(i,1,1) = c1*z(i,1)+c2*z(i,2)
    end do
    if (n .eq. 2) go to 5
    do i = 1,m
       zp(i,1,1) = zp(i,1,1)+c3*z(i,3)
    end do
!
! obtain y-partial derivatives along y = y(n)
!
5   if ((islpsw/16)*2 .ne. (islpsw/8)) go to 7
    do i = 1,m
       npi = n+i
       temp(npi) = zyn(i)
    end do
    go to 10
7   delyn = y(n)-y(nm1)
    delynm = delyn+delyn
    if (n .gt. 2) delynm = y(n)-y(n-2)
    if (delyn .le. 0. .or. delynm .le. delyn) go to 47
    call fitp_ceez (-delyn,-delynm,sigmay,c1,c2,c3,n)
    do i = 1,m
       npi = n+i
       temp(npi) = c1*z(i,n)+c2*z(i,nm1)
    end do
    if (n .eq. 2) go to 10
    do i = 1,m
       npi = n+i
       temp(npi) = temp(npi)+c3*z(i,n-2)
    end do
10  if (x(m) .le. x(1)) go to 47
!
! denormalize tension factor in x-direction
!
    sigmax = abs(sigma) * (m - 1) / (x(m) - x(1))
!
! obtain x-partial derivatives along x = x(1)
!
    if ((islpsw/2)*2 .ne. islpsw) go to 12
    do j = 1,n
       zp(1,j,2) = zx1(j)
    end do
    if ((islpsw/32)*2 .eq. (islpsw/16) .and. (islpsw/128)*2  .eq. (islpsw/64)) go to 15
12  delx1 = x(2)-x(1)
    delx2 = delx1+delx1
    if (m .gt. 2) delx2 = x(3)-x(1)
    if (delx1 .le. 0. .or. delx2 .le. delx1) go to 47
    call fitp_ceez (delx1,delx2,sigmax,c1,c2,c3,m)
    if ((islpsw/2)*2 .eq. islpsw) go to 15
    do j = 1,n
       zp(1,j,2) = c1*z(1,j)+c2*z(2,j)
    end do
    if (m .eq. 2) go to 15
    do j = 1,n
       zp(1,j,2) = zp(1,j,2)+c3*z(3,j)
    end do
!
! obtain x-y-partial derivative at (x(1),y(1))
!
15  if ((islpsw/32)*2 .ne. (islpsw/16)) go to 16
    zp(1,1,3) = zxy11
    go to 17
16  zp(1,1,3) = c1*zp(1,1,1)+c2*zp(2,1,1)
    if (m .gt. 2) zp(1,1,3) = zp(1,1,3)+c3*zp(3,1,1)
!
! obtain x-y-partial derivative at (x(1),y(n))
!
17  if ((islpsw/128)*2 .ne. (islpsw/64)) go to 18
    zxy1ns = zxy1n
    go to 19
18  zxy1ns = c1*temp(n+1)+c2*temp(n+2)
    if (m .gt. 2) zxy1ns = zxy1ns+c3*temp(n+3)
!
! obtain x-partial derivative along x = x(m)
!
19  if ((islpsw/4)*2 .ne. (islpsw/2)) go to 21
    do j = 1,n
       npmpj = npm+j
       temp(npmpj) = zxm(j)
    end do
    if ((islpsw/64)*2 .eq. (islpsw/32) .and.(islpsw/256)*2 .eq. (islpsw/128)) go to 24
21  delxm = x(m)-x(mm1)
    delxmm = delxm+delxm
    if (m .gt. 2) delxmm = x(m)-x(m-2)
    if (delxm .le. 0. .or. delxmm .le. delxm) go to 47
    call fitp_ceez (-delxm,-delxmm,sigmax,c1,c2,c3,m)
    if ((islpsw/4)*2 .eq. (islpsw/2)) go to 24
    do j = 1,n
       npmpj = npm+j
       temp(npmpj) = c1*z(m,j)+c2*z(mm1,j)
    end do
    if (m .eq. 2) go to 24
    do j = 1,n
       npmpj = npm+j
       temp(npmpj) = temp(npmpj)+c3*z(m-2,j)
    end do
!
! obtain x-y-partial derivative at (x(m),y(1))
!
24  if ((islpsw/64)*2 .ne. (islpsw/32)) go to 25
    zp(m,1,3) = zxym1
    go to 26
25  zp(m,1,3) = c1*zp(m,1,1)+c2*zp(mm1,1,1)
    if (m .gt. 2) zp(m,1,3) = zp(m,1,3)+c3*zp(m-2,1,1)
!
! obtain x-y-partial derivative at (x(m),y(n))
!
26  if ((islpsw/256)*2 .ne. (islpsw/128)) go to 27
    zxymns = zxymn
    go to 28
27  zxymns = c1*temp(npm)+c2*temp(npm-1)
    if (m .gt. 2) zxymns = zxymns+c3*temp(npm-2)
!
! set up right hand sides and tridiagonal system for y-grid
! perform forward elimination
!
28  del1 = y(2)-y(1)
    if (del1 .le. 0.) go to 47
    deli = 1./del1
    do i = 1,m
       zp(i,2,1) = deli*(z(i,2)-z(i,1))
    end do
    zp(1,2,3) = deli*(zp(1,2,2)-zp(1,1,2))
    zp(m,2,3) = deli*(temp(npm+2)-temp(npm+1))
    call fitp_terms (diag1,sdiag1,sigmay,del1)
    diagi = 1./diag1
    do i = 1,m
       zp(i,1,1) = diagi*(zp(i,2,1)-zp(i,1,1))
    end do
    zp(1,1,3) = diagi*(zp(1,2,3)-zp(1,1,3))
    zp(m,1,3) = diagi*(zp(m,2,3)-zp(m,1,3))
    temp(1) = diagi*sdiag1
    if (n .eq. 2) go to 34
    do j = 2,nm1
       jm1 = j-1
       jp1 = j+1
       npmpj = npm+j
       del2 = y(jp1)-y(j)
       if (del2 .le. 0.) go to 47
       deli = 1./del2
       do i = 1,m
          zp(i,jp1,1) = deli*(z(i,jp1)-z(i,j))
       end do
       zp(1,jp1,3) = deli*(zp(1,jp1,2)-zp(1,j,2))
       zp(m,jp1,3) = deli*(temp(npmpj+1)-temp(npmpj))
       call fitp_terms (diag2,sdiag2,sigmay,del2)
       diagin = 1./(diag1+diag2-sdiag1*temp(jm1))
       do i = 1,m
          zp(i,j,1) = diagin*(zp(i,jp1,1)-zp(i,j,1)-sdiag1*zp(i,jm1,1))
       end do
       zp(1,j,3) = diagin*(zp(1,jp1,3)-zp(1,j,3)-sdiag1*zp(1,jm1,3))
       zp(m,j,3) = diagin*(zp(m,jp1,3)-zp(m,j,3)-sdiag1*zp(m,jm1,3))
       temp(j) = diagin*sdiag2
       diag1 = diag2
       sdiag1 = sdiag2
    end do
34  diagin = 1./(diag1-sdiag1*temp(nm1))
    do i = 1,m
       npi = n+i
       zp(i,n,1) = diagin*(temp(npi)-zp(i,n,1)-sdiag1*zp(i,nm1,1))
    end do
    zp(1,n,3) = diagin*(zxy1ns-zp(1,n,3)-sdiag1*zp(1,nm1,3))
    temp(n) = diagin*(zxymns-zp(m,n,3)-sdiag1*zp(m,nm1,3))
!
! perform back substitution
!
    do j = 2,n
       jbak = np1-j
       jbakp1 = jbak+1
       t = temp(jbak)
       do i = 1,m
          zp(i,jbak,1) = zp(i,jbak,1)-t*zp(i,jbakp1,1)
       end do
       zp(1,jbak,3) = zp(1,jbak,3)-t*zp(1,jbakp1,3)
       temp(jbak) = zp(m,jbak,3)-t*temp(jbakp1)
    end do
!
! set up right hand sides and tridiagonal system for x-grid
! perform forward elimination
!
    del1 = x(2)-x(1)
    if (del1 .le. 0.) go to 47
    deli = 1./del1
    do j = 1,n
       zp(2,j,2) = deli*(z(2,j)-z(1,j))
       zp(2,j,3) = deli*(zp(2,j,1)-zp(1,j,1))
    end do
    call fitp_terms (diag1,sdiag1,sigmax,del1)
    diagi = 1./diag1
    do j = 1,n
       zp(1,j,2) = diagi*(zp(2,j,2)-zp(1,j,2))
       zp(1,j,3) = diagi*(zp(2,j,3)-zp(1,j,3))
    end do
    temp(n+1) = diagi*sdiag1
    if (m  .eq. 2) go to 43
    do i = 2,mm1
       im1 = i-1
       ip1 = i+1
       npi = n+i
       del2 = x(ip1)-x(i)
       if (del2 .le. 0.) go to 47
       deli = 1./del2
       do j = 1,n
          zp(ip1,j,2) = deli*(z(ip1,j)-z(i,j))
          zp(ip1,j,3) = deli*(zp(ip1,j,1)-zp(i,j,1))
       end do
       call fitp_terms (diag2,sdiag2,sigmax,del2)
       diagin = 1./(diag1+diag2-sdiag1*temp(npi-1))
       do j = 1,n
          zp(i,j,2) = diagin*(zp(ip1,j,2)-zp(i,j,2)-sdiag1*zp(im1,j,2))
          zp(i,j,3) = diagin*(zp(ip1,j,3)-zp(i,j,3)-sdiag1*zp(im1,j,3))
       end do
       temp(npi) = diagin*sdiag2
       diag1 = diag2
       sdiag1 = sdiag2
    end do
43  diagin = 1./(diag1-sdiag1*temp(npm-1))
    do j = 1,n
       npmpj = npm+j
       zp(m,j,2) = diagin*(temp(npmpj)-zp(m,j,2)-sdiag1*zp(mm1,j,2))
       zp(m,j,3) = diagin*(temp(j)-zp(m,j,3)-sdiag1*zp(mm1,j,3))
    end do
!
! perform back substitution
!
    do i = 2,m
       ibak = mp1-i
       ibakp1 = ibak+1
       npibak = n+ibak
       t = temp(npibak)
       do j = 1,n
          zp(ibak,j,2) = zp(ibak,j,2)-t*zp(ibakp1,j,2)
          zp(ibak,j,3) = zp(ibak,j,3)-t*zp(ibakp1,j,3)
       end do
    end do
    return
!
! too few points
!
46  ierr = 1
    return
!
! points not strictly increasing
!
47  ierr = 2
    return
  end subroutine fitp_surf1

  !> FIXME : Add documentation
  real function fitp_surf2 (xx,yy,m,n,x,y,z,iz,zp,sigma)
    implicit none
    
    real, intent(in) :: xx, yy, sigma
    integer, intent(in) :: m, n, iz
    real, dimension(m), intent(in) :: x
    real, dimension(n), intent(in) :: y
    real, dimension(iz,n), intent(in) :: z
    real, dimension(m,n,3), intent(in) :: zp
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function interpolates a surface at a given coordinate
! pair using a bi-spline under tension. the subroutine surf1
! should be called earlier to determine certain necessary
! parameters.
!
! on input--
!
!   xx and yy contain the x- and y-coordinates of the point
!   to be mapped onto the interpolating surface.
!
!   m and n contain the number of grid lines in the x- and
!   y-directions, respectively, of the rectangular grid
!   which specified the surface.
!
!   x and y are arrays containing the x- and y-grid values,
!   respectively, each in increasing order.
!
!   z is a matrix containing the m * n functional values
!   corresponding to the grid values (i. e. z(i,j) is the
!   surface value at the point (x(i),y(j)) for i = 1,...,m
!   and j = 1,...,n).
!
!   iz contains the row dimension of the array z as declared
!   in the calling program.
!
!   zp is an array of 3*m*n locations stored with the
!   various surface derivative information determined by
!   surf1.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters m, n, x, y, z, iz, zp, and sigma should be
! input unaltered from the output of surf1.
!
! on output--
!
!   surf2 contains the interpolated surface value.
!
! none of the input parameters are altered.
!
! this function references package modules intrvl and
! snhcsh.
!
!-----------------------------------------------------------
    integer :: im1, i, j, jm1
    real :: zxxi, zxxim1, zim1, zi
    real :: sigmax, del1, del2
    real :: sinhms, sinhm2, sinhm1, dels, sigmay
!
! inline one dimensional cubic spline interpolation
!
!/Moved to contained function
!    hermz (f1,f2,fp1,fp2) = (f2*del1+f1*del2)/dels-del1* &
!         del2*(fp2*(del1+dels)+ &
!         fp1*(del2+dels))/(6.*dels)
!
! inline one dimensional spline under tension interpolation
!
!/Moved to contained function
!    hermnz (f1,f2,fp1,fp2,sigmap) = (f2*del1+f1*del2)/dels &
!         +(fp2*del1*(sinhm1-sinhms) &
!         +fp1*del2*(sinhm2-sinhms) &
!         )/(sigmap*sigmap*dels*(1.+sinhms))

!
! denormalize tension factor in x and y direction
!
    sigmax = abs(sigma) * (m - 1) / (x(m) - x(1))
    sigmay = abs(sigma) * (n - 1) / (y(n) - y(1))
!
! determine y interval
!
    jm1 = fitp_intrvl (yy,y,n)
    j = jm1+1
!
! determine x interval
!
    im1 = fitp_intrvl (xx,x,m)
    i = im1+1
    del1 = yy-y(jm1)
    del2 = y(j)-yy
    dels = y(j)-y(jm1)
    if (is_zero(sigmay)) then
       ! perform four interpolations in y-direction
       zim1 = hermz(z(i-1,j-1),z(i-1,j),zp(i-1,j-1,1),zp(i-1,j,1))
       zi = hermz(z(i,j-1),z(i,j),zp(i,j-1,1),zp(i,j,1))
       zxxim1 = hermz(zp(i-1,j-1,2),zp(i-1,j,2),zp(i-1,j-1,3),zp(i-1,j,3))
       zxxi = hermz(zp(i,j-1,2),zp(i,j,2),zp(i,j-1,3),zp(i,j,3))
    else
       sinhm1 = sinhm_fun(sigmay * del1)
       sinhm2 = sinhm_fun(sigmay * del2)
       sinhms = sinhm_fun(sigmay * dels)
       zim1 = hermnz(z(i-1,j-1),z(i-1,j),zp(i-1,j-1,1),zp(i-1,j,1),sigmay)
       zi = hermnz(z(i,j-1),z(i,j),zp(i,j-1,1),zp(i,j,1),sigmay)
       zxxim1 = hermnz(zp(i-1,j-1,2),zp(i-1,j,2),zp(i-1,j-1,3),zp(i-1,j,3),sigmay)
       zxxi = hermnz(zp(i,j-1,2),zp(i,j,2),zp(i,j-1,3),zp(i,j,3),sigmay)
    end if

    ! perform final interpolation in x-direction

    del1 = xx-x(im1)
    del2 = x(i)-xx
    dels = x(i)-x(im1)

    if (is_zero(sigmax)) then
       fitp_surf2 = hermz(zim1,zi,zxxim1,zxxi)
    else
       sinhm1 = sinhm_fun(sigmax * del1)
       sinhm2 = sinhm_fun(sigmax * del2)
       sinhms = sinhm_fun(sigmax * dels)
       fitp_surf2 = hermnz(zim1,zi,zxxim1,zxxi,sigmax)
    end if
  contains
    !> FIXME : Add documentation
    elemental real function hermz(f1,f2,fp1,fp2)
      !Note del1,sinhm1 etc. are known because we are in subroutine scope
      implicit none
      real, intent(in) :: f1,f2,fp1,fp2
      hermz= (f2*del1+f1*del2)/dels-del1* &
           del2*(fp2*(del1+dels)+ &
           fp1*(del2+dels))/(6.*dels)
    end function hermz

    !> FIXME : Add documentation
    elemental real function hermnz(f1,f2,fp1,fp2,sigmap)
      !Note del1,sinhm1 etc. are known because we are in subroutine scope
      implicit none
      real, intent(in) :: f1,f2,fp1,fp2,sigmap
      hermnz= (f2*del1+f1*del2)/dels &
           +(fp2*del1*(sinhm1-sinhms) &
           +fp1*del2*(sinhm2-sinhms) &
           )/(sigmap*sigmap*dels*(1.+sinhms))
    end function hermnz
  end function fitp_surf2

  !> FIXME : Add documentation  
  pure subroutine fitp_ceez (del1,del2,sigma,c1,c2,c3,n)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: del1, del2, sigma
    real, intent(out) :: c1, c2, c3
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the coefficients c1, c2, and c3
! used to determine endpoint slopes. specifically, if
! function values y1, y2, and y3 are given at points x1, x2,
! and x3, respectively, the quantity c1*y1 + c2*y2 + c3*y3
! is the value of the derivative at x1 of a spline under
! tension (with tension factor sigma) passing through the
! three points and having third derivative equal to zero at
! x1. optionally, only two values, c1 and c2 are determined.
!
! on input--
!
!   del1 is x2-x1 (.gt. 0.).
!
!   del2 is x3-x1 (.gt. 0.). if n .eq. 2, this parameter is
!   ignored.
!
!   sigma is the tension factor.
!
! and
!
!   n is a switch indicating the number of coefficients to
!   be returned. if n .eq. 2 only two coefficients are
!   returned. otherwise all three are returned.
!
! on output--
!
!   c1, c2, and c3 contain the coefficients.
!
! none of the input parameters are altered.
!
! this subroutine references package module snhcsh.
!
!-----------------------------------------------------------
    real :: delm, delp, sinhmp, denom, sinhmm, del, coshm2, coshm1

    if (n == 2) then
      ! two coefficients
      c1 = -1./del1
      c2 = -c1
    else if (is_zero(sigma)) then
      ! tension .eq. 0.
      del = del2-del1
      c1 = -(del1+del2)/(del1*del2)
      c2 = del2/(del1*del)
      c3 = -del1/(del2*del)
   else
      ! tension .ne. 0.
      coshm1 = coshm_fun(sigma * del1)
      coshm2 = coshm_fun(sigma * del2)
      delp = sigma*(del2+del1)/2.
      delm = sigma*(del2-del1)/2.
      sinhmp = sinhm_fun(delp)
      sinhmm = sinhm_fun(delm)
      denom = coshm1*(del2-del1)-2.*del1*delp*delm*(1.+sinhmp)*(1.+sinhmm)
      c1 = 2.*delp*delm*(1.+sinhmp)*(1.+sinhmm)/denom
      c2 = -coshm2/denom
      c3 = coshm1/denom
   end if
  end subroutine fitp_ceez

  !> FIXME : Add documentation  
  integer function fitp_intrvl (t,x,n)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: t
    real, dimension(n), intent(in) :: x
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function determines the index of the interval
! (determined by a given increasing sequence) in which
! a given value lies.
!
! on input--
!
!   t is the given value.
!
!   x is a vector of strictly increasing values.
!
! and
!
!   n is the length of x (n .ge. 2).
!
! on output--
!
!   intrvl returns an integer i such that
!
!          i =  1       if         e   t .lt. x(2)  ,
!          i =  n-1     if x(n-1) .le. t            ,
!          otherwise       x(i)  .le. t .le. x(i+1),
!
! none of the input parameters are altered.
!
!-----------------------------------------------------------
    fitp_intrvl = fitp_intrv_helper(t, x, n, .false.)
  end function fitp_intrvl

  !> FIXME : Add documentation  
  integer function fitp_intrvp (t,x,n)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: t
    real, dimension(n), intent(in) :: x
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function determines the index of the interval
! (determined by a given increasing sequence) in which a
! given value lies, after translating the value to within
! the correct period.  it also returns this translated value.
!
! on input--
!
!   t is the given value.
!
!   x is a vector of strictly increasing values.
!
!   n is the length of x (n .ge. 2).
!
! on output--
!
!   intrvl returns an integer i such that
!
!          i = 1       if             tp .lt. x(2)  ,
!          i = n       if   x(n) .le. tp            ,
!          otherwise       x(i)  .le. tp .lt. x(i+1),
!
! none of the input parameters are altered.
!
!-----------------------------------------------------------
    fitp_intrvp = fitp_intrv_helper(t, x, n, .true.)
  end function fitp_intrvp

  !> Given a value t, and stating value of abscissae, x, periodic with period p
  !> map t to the equivalent value in the range of x (i.e. between x(1) and x(1) + p)
  elemental real function fitp_periodic_wrap_value(t, x_start, p) result(wrapped_value)
    real, intent(in) :: t, x_start, p
    integer :: nper
    nper = int((t - x_start) / p)
    wrapped_value = t - nper * p
    if (wrapped_value < x_start) wrapped_value = wrapped_value + p
  end function fitp_periodic_wrap_value

  integer function fitp_intrv_helper(tt, x, n, periodic) result(index)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: tt
    real, dimension(n), intent(in) :: x
    logical, intent(in) :: periodic
    integer, save :: i = 1
    integer :: il, ih, upper
    if (periodic) then
       upper = n
    else
       upper = n - 1
    end if

!
! check for illegal i
!
    if (i .ge. n) i = n / 2
!
! check old interval and extremes
!
    if (tt .lt. x(i)) then
       if (tt .le. x(2)) then
          i = 1
          index = 1
          return
       else
          il = 2
          ih = i
       end if
    else if (tt .le. x(i+1)) then
       index = i
       return
    else if (tt .ge. x(upper)) then
       i = upper
       index = upper
       return
    else
       il = i+1
       ih = upper
    end if

    !
    ! binary search loop
    !
    do while (.true.)
       i = (il + ih) / 2
       if (tt < x(i)) then
          ih = i
       else if (tt > x(i + 1)) then
          il = i + 1
       else
          index = i
          exit
       end if
    end do
  end function fitp_intrv_helper

  !> Calculate sinhm(x) = sinh(x) / x - 1
  elemental real function sinhm_fun(x)
    real, intent(in) :: x
    if (is_zero(x)) then
       sinhm_fun = 0
    else
       sinhm_fun = sinh(x) / x - 1
    end if
  end function sinhm_fun

  !> Calculate coshm(x) = cosh(x) - 1
  elemental real function coshm_fun(x)
    real, intent(in) :: x
    coshm_fun = cosh(x) - 1
  end function coshm_fun

  !> Calculate coshmm(x) = (cosh(x) - 1 - x * x / 2) / (x * x)
  elemental real function coshmm_fun(x)
    real, intent(in) :: x
    if (is_zero(x)) then
       coshmm_fun = 0
    else
       coshmm_fun = (cosh(x) - 1 - x*x/2) / (x*x)
    end if
  end function coshmm_fun

  !> FIXME : Add documentation  
  pure subroutine fitp_terms (diag,sdiag,sigma,del)
    implicit none
    real, intent(in) :: sigma, del
    real, intent(out) :: diag, sdiag
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine computes the diagonal and superdiagonal
! terms of the tridiagonal linear system associated with
! spline under tension interpolation.
!
! on input--
!
!   sigma contains the tension factor.
!
! and
!
!   del contains the step size.
!
! on output--
!
!                sigma*del*cosh(sigma*del) - sinh(sigma*del)
!   diag = del*--------------------------------------------.
!                     (sigma*del)**2 * sinh(sigma*del)
!
!                   sinh(sigma*del) - sigma*del
!   sdiag = del*----------------------------------.
!                (sigma*del)**2 * sinh(sigma*del)
!
! and
!
!   sigma and del are unaltered.
!
! this subroutine references package module snhcsh.
!
!-----------------------------------------------------------
    real :: coshm, denom, sigdel, sinhm

    if (is_zero(sigma)) then
       diag = del/3.
       sdiag = del/6.
    else
       sigdel = sigma*del
       sinhm = sinhm_fun(sigdel)
       coshm = coshm_fun(sigdel)
       denom = sigma*sigdel*(1.+sinhm)
       diag = (coshm-sinhm)/denom
       sdiag = sinhm/denom
       ! Note we could use the cosh and sinh intrinsics to evalute the expression
       ! directly as below. This is found to agree very well with the existing form.
       !diag = del * (sigdel*cosh(sigdel) - sinh(sigdel))/(sigdel*sigdel*sinh(sigdel))
       !sdiag = del * (sinh(sigdel) - sigdel)/(sigdel*sigdel*sinh(sigdel))
    end if
  end subroutine fitp_terms

end module periodic_splines
