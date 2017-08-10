module gridgen4mod

contains

  subroutine gridgen4read (filename,ntheta,nperiod,ntgrid,nlambda,theta,alambda, &
       gbdrift,gradpar,cvdrift,gds2,bmag,gds21,gds22,cvdrift0,gbdrift0)
    implicit none
    character(*), intent (in) :: filename
    integer, intent (in out) :: ntheta, nperiod, ntgrid, nlambda
    real, dimension (-ntgrid:ntgrid) :: theta
    real, dimension (nlambda) :: alambda
    real, dimension (-ntgrid:ntgrid) :: gbdrift, gradpar
    real, dimension (-ntgrid:ntgrid) :: cvdrift, gds2, bmag
    real, dimension (-ntgrid:ntgrid) :: gds21, gds22
    real, dimension (-ntgrid:ntgrid) :: cvdrift0, gbdrift0

    integer :: nthetain, nperiodin, ntgridin, nlambdain
    integer :: i
    integer :: unit
    logical :: od
    character(200) :: line

!CMR, August 2010:
!  (i) modify inquire to use opened instead of read, write and readwrite
!   needed for gfortran compiler: TT adopted same solution in file_utils
!  (ii) terminate with error message if no free LUN found
!CMRend

    unit=0
    do i = 10,100
       inquire (unit=i,opened=od)
       if ( .not. od ) then
          unit=i
          exit
       end if
    end do
    if (unit .eq. 0) then
       write(6,*) "gridgen4read:  no free LUN between 10,100 => force quit"
       stop
    endif

    open (unit=unit, file=filename, status="old")
    
    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt=*) nlambdain
    if (nlambdain > nlambda) then
       print *, "grid.out:nlambda > nlambda: ", nlambdain, nlambda
       stop
    end if
    nlambda = nlambdain
    read (unit=unit, fmt="(a)") line
    do i = 1, nlambda
       read (unit=unit, fmt=*) alambda(i)
    end do
    
    read (unit=unit, fmt="(a)") line
    read (unit=unit, fmt=*) ntgridin, nperiodin, nthetain
    if (ntgridin > ntgrid) then
       print *, "grid.out:ntgrid > ntgrid: ", ntgridin, ntgrid
       stop
    end if
    ntgrid = ntgridin
    nperiod = nperiodin
    ntheta = nthetain
    
    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) gbdrift(i), gradpar(i)
    end do
    
    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) cvdrift(i), gds2(i), bmag(i), theta(i)
    end do
    
    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) gds21(i), gds22(i)
    end do
    
    read (unit=unit, fmt="(a)") line
    do i = -ntgrid, ntgrid
       read (unit=unit, fmt=*) cvdrift0(i), gbdrift0(i)
    end do
    
    close (unit=unit)
  end subroutine gridgen4read
  
  subroutine gridgen4 (n,nbmag,thetain,bmagin, npadd, &
       alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw, &
       ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
    implicit none
    integer, intent (in) :: nbmag, n
    real, dimension (nbmag), intent (in) :: thetain, bmagin
    integer, intent (in) :: npadd
    real, intent (in) :: alknob, epsknob, bpknob, extrknob
    real, intent (in) :: thetamax, deltaw, widthw
    integer, intent (in out) :: ntheta, nlambda
    real, dimension (ntheta+1), intent (out) :: thetagrid, bmaggrid
    real, dimension (nlambda), intent (out) :: alambdagrid
    
    call gridgen4_1 (n,nbmag,thetain,bmagin, npadd, &
         alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,1.0, &
         ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
  end subroutine gridgen4
  
  subroutine gridgen4_1 (n,nbmag,thetain,bmagin, npadd, &
       alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
       ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: nbmag
    real, dimension (nbmag), intent (in) :: thetain, bmagin
    integer, intent (in) :: npadd
    real, intent (in) :: alknob, epsknob, bpknob, extrknob
    real, intent (in) :: thetamax, deltaw, widthw, tension
    integer, intent (in out) :: ntheta, nlambda
    real, dimension (ntheta+1), intent (out) :: thetagrid, bmaggrid
    real, dimension (nlambda), intent (out) :: alambdagrid
    
    call gridgen4_2 (n,nbmag,thetain,bmagin, npadd, &
         alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
         ntheta,nlambda,thetagrid,bmaggrid,alambdagrid)
    alambdagrid(1:nlambda) = 1.0/alambdagrid(1:nlambda)
  end subroutine gridgen4_1
  
  subroutine gridgen4_2 (n,nbmag,thetain,bmagin, npadd, &
       alknob,epsknob,bpknob,extrknob,thetamax,deltaw,widthw,tension, &
       ntheta,nbset,thetagrid,bmaggrid,bset)
    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: nbmag
    real, dimension (nbmag), intent (in) :: thetain, bmagin
    integer, intent (in) :: npadd
    real, intent (in) :: alknob, epsknob, bpknob, extrknob
    real, intent (in) :: thetamax, deltaw, widthw, tension
    integer, intent (in out) :: ntheta, nbset
    real, dimension (ntheta+1), intent (out) :: thetagrid, bmaggrid
    real, dimension (nbset), intent (out) :: bset

    real, parameter :: pi = 3.1415926535897931, twopi = 6.2831853071795862
    real :: npi
    integer, parameter :: maxiter = 100
    logical, parameter :: debug_output = .true.

    real, dimension (:), allocatable :: bmagspl

    integer :: nstart
    real, allocatable, dimension (:) :: thetastart, bmagstart
    logical, allocatable, dimension (:) :: essentialstart
    real, allocatable, dimension (:) :: bmin, bmax, bextr, thetaextr
    real, allocatable, dimension (:) :: bprime

    integer :: nset, nsetset
    real, allocatable, dimension (:) :: thetasetset, bmagset
    integer, allocatable, dimension (:) :: ibmagsetset, icollsetset
    logical, allocatable, dimension (:) :: essentialset
    
    integer, allocatable, dimension (:) :: ithetasort, ibmagsort
    
    real, allocatable, dimension (:) :: thetares, alambdares
    
    integer :: nthetaout, nlambdaout
    
    integer :: debug_unit
    logical :: debug=.false.    
    npi = n * twopi

if (debug) write(6,*) "gridgen4_2: call gg4init"
    call gg4init
if (debug) write(6,*) "gridgen4_2: call gg4debug"
    call gg4debug (nbmag, thetain, bmagin, "input grid")

! 1 Set up starting grid.
if (debug) write(6,*) "gridgen4_2: call gg4start"
    call gg4start (n)
if (debug) write(6,*) "gridgen4_2: call gg4debug"
    call gg4debug (nstart, thetastart, thetastart, "starting grid")

! 2 Collect all bounce points associated with starting grid.
if (debug) write(6,*) "gridgen4_2: call gg4collect"
    call gg4collect
if (debug) write(6,*) "gridgen4_2: call gg4debug"
    call gg4debug (nset, bmagset, bmagset, "bmagset")
if (debug) write(6,*) "gridgen4_2: call gg4debugi"
    call gg4debugi (nsetset,thetasetset,ibmagsetset,icollsetset,"thetasetset")

! 3 Sort collected grids.
if (debug) write(6,*) "gridgen4_2: call gg4sort"
    call gg4sort

! 4 Calculate spacing and resolution metrics.
if (debug) write(6,*) "gridgen4_2: call gg4metrics"
    call gg4metrics

! 5 Remove sets until the number of points is small enough.
if (debug) write(6,*) "gridgen4_2: call gg4metrics"
    call gg4remove

! 6 Build output grids.
if (debug) write(6,*) "gridgen4_2: call gg4results"
    call gg4results (n)
if (debug) write(6,*) "gridgen4_2: call gg4debug"
    call gg4debug (nthetaout, thetagrid, bmaggrid, "output grid")
if (debug) write(6,*) "gridgen4_2: call gg4debugi"
    call gg4debug (nlambdaout, bset, bset, "output 1/lambda")

    ntheta = 2*(nthetaout/2)
    nbset = nlambdaout

! 7 Clean up.
if (debug) write(6,*) "gridgen4_2: call gg4finish"
    call gg4finish

  contains
  
    subroutine gg4init
      use splines, only: fitp_curvp1
      implicit none
      logical :: od
      real, dimension (2*nbmag) :: tmp
      integer :: ierr, i

      ierr = 0
      allocate (bmagspl(nbmag)) ; bmagspl = 0.
      call fitp_curvp1 (nbmag-1,thetain,bmagin,npi,bmagspl,tmp,tension,ierr)
      if (ierr /= 0) then
         print *, "CURVP1: IERR=", ierr
         select case (ierr)
         case (1)
            print *, "N is less than 2"
         case (2)
            print *, "P is less than or equal to X(N)-X(1)"
         case (3)
            print *, "X values are not strictly increasing"
         end select
         write(6,*) 'gg4init: stopping with ierr=',ierr
         stop
      end if

      nstart = 0

      nset = 0
      nsetset = 0

      if (debug_output) then
         debug_unit=0
!CMR, August 2010:
!  (i) modify inquire to use opened instead of read, write and readwrite
!   needed for gfortran compiler: TT adopted same solution in file_utils
!  (ii) terminate with error message if no free LUN found
!CMRend
         do i = 10,100
            inquire (unit=i,opened=od)
            if ( .not. od ) then
               debug_unit=i
               exit
            end if
         end do
         if (debug_unit .eq. 0) then
            write(6,*) "gg4init:  no free LUN between 10,100 => force quit"
            stop
         endif
         open (unit=debug_unit, file="gridgen.200", status="unknown")
         write (unit=debug_unit, fmt=*) "nbmag=", nbmag
         write (unit=debug_unit, fmt=*) "thetain,bmagin="
         do i=1,nbmag
            write (unit=debug_unit, fmt=*) thetain(i), bmagin(i)
         end do
         write (unit=debug_unit, fmt=*) "alknob=", alknob
         write (unit=debug_unit, fmt=*) "epsknob=", epsknob
         write (unit=debug_unit, fmt=*) "bpknob=", bpknob
         write (unit=debug_unit, fmt=*) "extrknob=", extrknob
         write (unit=debug_unit, fmt=*) "thetamax=", thetamax
         write (unit=debug_unit, fmt=*) "deltaw,widthw=", deltaw, widthw
         write (unit=debug_unit, fmt=*) "tension=", tension
         write (unit=debug_unit, fmt=*) "ntheta,nlambda=", ntheta, nbset
      end if
      if (debug) write(6,*) 'gg4init: end'
    end subroutine gg4init

    subroutine gg4finish
      implicit none
      
      if (debug_output) then
         close (unit=debug_unit)
      end if
      
      deallocate (bmagspl)
      deallocate (thetastart)
      deallocate (bmagstart)
      deallocate (bmin,bmax,bextr,thetaextr,bprime)
      deallocate (essentialstart)

      deallocate (thetasetset, bmagset)
      deallocate (ibmagsetset, icollsetset)
      deallocate (essentialset)
      
      deallocate (ithetasort, ibmagsort)

      deallocate (thetares, alambdares)

      if (mod(nthetaout,2) /= 1) then
         print *, "gridgen4_1:gg4results:nthetaout=",nthetaout
         print *, "nthetaout is not odd, so there is a problem with gridgen"
!       stop
      end if
    end subroutine gg4finish

    subroutine old_gg4finish
      implicit none
      
      deallocate (thetastart, essentialstart)

      deallocate (thetasetset, bmagset)
      deallocate (ibmagsetset, icollsetset)
      deallocate (essentialset)
      
      deallocate (ithetasort, ibmagsort)
      
      deallocate (thetares, alambdares)
      
      if (debug_output) then
         close (unit=debug_unit)
      end if
    end subroutine old_gg4finish

    subroutine gg4start (n)
      implicit none
      integer, intent (in) :: n
      integer :: i, j
      real :: thetal, thetar, theta0
      real :: bprime0
      integer :: nextr
      
      allocate (bmin(nbmag-1),bmax(nbmag-1),bextr(nbmag-1),thetaextr(nbmag-1))
      allocate (bprime(nbmag))
      bmin = min(bmagin(1:nbmag-1),bmagin(2:nbmag))
      bmax = max(bmagin(1:nbmag-1),bmagin(2:nbmag))
      bextr = 0.0
      thetaextr = 0.0
      
      do i = 1, nbmag
         bprime(i) = bmagp(thetain(i))
      end do

! Find all extrema.
      nextr = 0
      do i = 2, nbmag-2
         if (abs(bprime(i)) < epsilon(0.)) then
            bextr(i) = bmagin(i)
            thetaextr(i) = thetain(i)
            nextr = nextr + 1
         else if (bprime(i) < 0.0 .and. bprime(i+1) > 0.0) then
            thetal = thetain(i)
            thetar = thetain(i+1)
            do j = 1, maxiter
               theta0 = 0.5*(thetal+thetar)
               bprime0 = bmagp(theta0)
               if (bprime0 < 0.0) then
                  thetal = theta0
               else if (bprime0 > 0.0) then
                  thetar = theta0
               else
                  exit
               end if
            end do
            if (theta0-thetain(i) < epsknob) then
               theta0 = thetain(i)
               bprime(i) = 0.0
            else if (thetain(i+1)-theta0 < epsknob) then
               theta0 = thetain(i+1)
               bprime(i+1) = 0.0
            end if
            thetaextr(i) = theta0
            bextr(i) = bmagint(theta0)
            nextr = nextr + 1
         else if (bprime(i) > 0.0 .and. bprime(i+1) < 0.0) then
            thetal = thetain(i)
            thetar = thetain(i+1)
            do j = 1, maxiter
               theta0 = 0.5*(thetal+thetar)
               bprime0 = bmagp(theta0)
               if (bprime0 > 0.0) then
                  thetal = theta0
               else if (bprime0 < 0.0) then
                  thetar = theta0
               else
                  exit
               end if
            end do
            if (theta0-thetain(i) < epsknob) then
               theta0 = thetain(i)
               bprime(i) = 0.0
            else if (thetain(i+1)-theta0 < epsknob) then
               theta0 = thetain(i+1)
               bprime(i+1) = 0.0
            end if
            thetaextr(i) = theta0
            bextr(i) = bmagint(theta0)
            nextr = nextr + 1
         end if
      end do
    
! Collect -pi, all local extrema, original grid points, points between
! original grid points
      nstart = 1 + nextr + (nbmag-1)*(1+npadd)

      allocate (thetastart(nstart), essentialstart(nstart), bmagstart(nstart))
      essentialstart(:1+nextr) = .true.
      essentialstart(2+nextr:) = .false.
      thetastart(1) = -pi*n
      thetastart(2:nextr+1) = pack(thetaextr,abs(bextr) > epsilon(0.))
      
      nstart = nextr + 1
      do i = 1, nbmag-1
         do j = 0, npadd
            thetastart(nstart+1) &
                 = thetain(i) + real(j)/real(npadd+1)*(thetain(i+1)-thetain(i))
            if (all(abs(thetastart(nstart+1)-thetastart(1:nstart)) > epsknob)) then
               nstart = nstart + 1
            end if
         end do
      end do
      
      do i = 1, nstart
         bmagstart(i) = bmagint(thetastart(i))
      end do
    end subroutine gg4start

    subroutine old_gg4start
      implicit none
      integer :: i, j
      real :: bprimel, bprimer, bprime0
      real :: thetal, thetar, theta0
      real, parameter :: tol=1e-9
      
      allocate (thetastart(nbmag*(2+npadd)))
      allocate (essentialstart(nbmag*(2+npadd)))
      
! 1 Collect essential points: -pi + all local extrema.
! 1.1 Collect -pi.
      call add_start (-pi, .true.)
! 1.2 Collect extrema.
      bprimer = bmagp(thetain(2))
      do i = 2, nbmag-2
         bprimel = bprimer
         bprimer = bmagp(thetain(i+1))
         if (abs(bprimel) < epsilon(0.)) then
            call add_start (thetain(i), .true.)
         else if (bprimel < 0.0 .and. bprimer > 0.0) then
! 1.2.1 Find local minimum.
            thetal = thetain(i)
            thetar = thetain(i+1)
            do j = 1, maxiter
               theta0 = 0.5*(thetal+thetar)
               bprime0 = bmagp(theta0)
               if (bprime0 < 0.0) then
                  thetal = theta0
               else if (bprime0 > 0.0) then
                  thetar = theta0
               else
                  exit
               end if
            end do
            if (theta0-thetain(i) < tol) theta0 = thetain(i)
            if (thetain(i+1)-theta0 < tol) theta0 = thetain(i+1)
            call add_start (theta0, .true.)
         else if (bprimel > 0.0 .and. bprimer < 0.0) then
! 1.2.2 Find local maximum.
            thetal = thetain(i)
            thetar = thetain(i+1)
            do j = 1, maxiter
               theta0 = 0.5*(thetal+thetar)
               bprime0 = bmagp(theta0)
               if (bprime0 > 0.0) then
                  thetal = theta0
               else if (bprime0 < 0.0) then
                  thetar = theta0
               else
                  exit
               end if
            end do
            if (theta0-thetain(i) < tol) then
               theta0 = thetain(i)
               bprime(i) = 0.0
            else if (thetain(i+1)-theta0 < tol) then
               theta0 = thetain(i+1)
               bprime(i+1) = 0.0
            end if
            call add_start (theta0, .true.)
         end if
      end do

! 2 Collect original grid, except for extrema.
      do i = 2, nbmag-1
         if (all(abs(thetain(i)-thetastart(1:nstart)) > tol)) then
            call add_start(thetain(i), .false.)
         end if
      end do

! 3 Collect points between original grid points.
      do i = 1, nbmag-1
         do j = 1, npadd
            theta0 = thetain(i) + real(j)/real(npadd+1)*(thetain(i+1)-thetain(i))
            if (all(abs(theta0-thetastart(1:nstart)) > tol)) then
               call add_start(theta0, .false.)
            end if
         end do
      end do
    end subroutine old_gg4start

    subroutine add_start (theta, essential)
      implicit none
      real, intent (in) :: theta
      logical, intent (in) :: essential

      real, allocatable, dimension (:) :: thetatmp
      logical, allocatable, dimension (:) :: essentialtmp

      if (nstart >= size(thetastart)) then
         allocate (thetatmp(nstart))
         allocate (essentialtmp(nstart))
         thetatmp = thetastart(1:nstart)
         essentialtmp = essentialstart(1:nstart)
         deallocate (thetastart,essentialstart)
         allocate (thetastart(nstart+nbmag))
         allocate (essentialstart(nstart+nbmag))
         thetastart(1:nstart) = thetatmp
         essentialstart(1:nstart) = essentialtmp
         deallocate (thetatmp,essentialtmp)
      end if
      nstart = nstart + 1
      thetastart(nstart) = theta
      essentialstart(nstart) = essential
    end subroutine add_start

    subroutine gg4collect
      implicit none
      integer :: i, iset, j
      integer :: nsetsetmax
      real :: thetai, bmagi
      
! Estimate upper bound on number of points.
!    nsetsetmax = count( &
!     (spread(bmagstart,2,nbmag-1) >= spread(bmin,1,nstart) &
!      .and. spread(bmagstart,2,nbmag-1) <= spread(bmax,1,nstart))) &
!             + 2*count( &
!     (spread(bextr,1,nstart) /= 0 &
!      .and. &
!       ((spread(bmagstart,2,nbmag-1) >= spread(bmax,1,nstart) &
!         .and. spread(bmagstart,2,nbmag-1) <= spread(bextr,1,nstart)) &
!        .or. &
!        (spread(bmagstart,2,nbmag-1) <= spread(bmin,1,nstart) &
!         .and. spread(bmagstart,2,nbmag-1) >= spread(bextr,1,nstart)))))

      nsetsetmax = 0
      do i = 1, nbmag-1
         do j = 1, nstart
            if (bmagstart(j) >= bmin(i) &
                 .and. bmagstart(j) <= bmax(i)) &
                 nsetsetmax = nsetsetmax + 1
            if (abs(bextr(i)) > epsilon(0.) &
                 .and. ((bmagstart(j) >= bmax(i) &
                 .and. bmagstart(j) <= bextr(i)) &
                 .or. (bmagstart(j) <= bmin(i) &
                 .and. bmagstart(j) >= bextr(i)))) &
                 nsetsetmax = nsetsetmax + 2
         end do
      end do
! End of Liu's fix.

! Allocate
      allocate (thetasetset(nsetsetmax))
      allocate (ibmagsetset(nsetsetmax))
      allocate (icollsetset(nsetsetmax))

      nsetset = 0
      starting_points: do iset = 1, nstart
         thetai = thetastart(iset)
         bmagi = bmagstart(iset)
         
! 1 For extrema, check previous sets to attach to.
         if (essentialstart(iset)) then
            do i = 1, iset-1
               if (abs(bmagstart(i)-bmagi) < epsknob) then
! 1.1.1 If the extremum does belong in this set, eliminate points
!       near the extremum from this set.
                  where (ibmagsetset(1:nsetset) == i &
                       .and. abs(thetasetset(1:nsetset)-thetai) < epsknob)
                  ibmagsetset(1:nsetset) = 0
                  end where
! 1.1.2 Attach the extremum to this set.
                  call add_setset (thetai, i, 112)
                  bmagstart(iset) = 0.0
                  cycle starting_points
               end if
            end do
         end if

! 2 Start a new set.
         call add_setset (thetai, iset, 2)

! 3 Check each original grid interval for matching bmag.
         grid_interval: do i = 1, nbmag-1
! 3.0.1 Stoopid problems near -pi.
! 'original' coding:
!            if (iset == 1 .and. i == 1) cycle grid_interval
! modification 5.16.02:
            if (iset == 1 .and. (i == 1 .or. i == nbmag-1)) cycle grid_interval
            if (bmagin(i) > bmagi .and. abs(thetain(i)-thetai) > epsilon(0.)) then
! 3.1 Consider when the left grid point is greater than the target bmag.
!    Then, there are three cases in which there are matching points.
!    (1) The right grid point is equal to the target bmag, and the slope
!        at the right point is positive.
!    (2) The right gridpoint is less than the target bmag.
!    (3) The right grid point is greater than the target bmag,
!        and the interval is concave down, and the minimum in
!        this interval, which is guaranteed to exist, is less than
!        the target bmag.
               if ((abs(bmagin(i+1)-bmagi) < epsilon(0.) &
                    .or. abs(thetain(i+1)-thetai) < epsilon(0.)) &
                    .and. bprime(i+1) > 0.0) then
! 3.1.1 Consider when the right grid point is equal to the target bmag.
                  call add_setset_root (thetain(i+1),thetain(i),bmagi,iset,311)
                  cycle grid_interval
               end if
               if (bmagin(i+1) < bmagi) then
! 3.1.2 Consider when the right grid point is less than the target bmag.
! 3.1.2.1 If this interval bounds the starting theta point, that is what
!         the target point is, and it is already collected.
                  if (thetai >= thetain(i) .and. thetai <= thetain(i+1)) then
                     cycle grid_interval
                  end if
! 3.1.2.2 Otherwise, find and collect the target point.
                  call add_setset_root (thetain(i+1),thetain(i),bmagi,iset,3122)
                  cycle grid_interval
               end if
! 3.1.3 Check if the grid interval is concave down.
! 3.1.3.1 If not, skip to the next interval.
               if (bprime(i) >= 0.0 .or. bprime(i+1) <= 0.0) cycle grid_interval
! 3.1.3.2 Consider the case where the starting theta point is within
!         this interval.
               if (thetai > thetain(i) .and. thetai < thetain(i+1)) then
! 3.1.3.2.1 If the starting point is an extremum, skip to next interval.
                  if (essentialstart(iset)) cycle grid_interval
! 3.1.3.2.2 Otherwise, the other target point is right of the starting
!           point is the slope is negative, and left if positive.
                  if (bprime(i) < 0.0) then
                     call add_setset_root (thetai,thetain(i+1),bmagi,iset,313221)
                  else
                     call add_setset_root (thetai,thetain(i),bmagi,iset,313222)
                  end if
                  cycle grid_interval
               end if
! 3.1.3.3 If the minimum within this interval is less than the target bmag,
!         then there are two target points in this interval, one on each side
!         of the minimum.
! 3.1.3.3.1 If this interval is on the edge, there will not be any minimum.
               if (i == 1 .or. i == nbmag-1) cycle grid_interval
! 3.1.3.3.2 Find the minimum in this interval.
               if (abs(bextr(i)) < epsilon(0.)) then
                  print *, "gridgen4.f90:gg4collect:3.1.3.3.2:"," missing extremum"
                  print *, "iset,i:",iset,i
                  print *, "bmagi:",bmagi
                  print *, "bprimel,bprimer:", bprime(i),bprime(i+1)
                  print *, "thetain(i):", thetain(i)
                  print *, "thetain(i+1):", thetain(i+1)
                  print *, "thetai:", thetai
                  stop
               end if
! 3.1.3.3.2.1 If the minimum is greater than the target bmag, skip to
!             the next interval.
               if (bextr(i) > bmagi) cycle grid_interval
! 3.1.3.3.2.2 Collect the point left of the minimum.
               call add_setset_root (thetaextr(i),thetain(i),bmagi,iset,313322)
! 3.1.3.3.2.3 Collect the point right of the minimum.
               call add_setset_root (thetaextr(i),thetain(i+1),bmagi,iset,313323)
               cycle grid_interval
            else if (bmagin(i) < bmagi .and. abs(thetain(i)-thetai) > epsilon(0.)) then
! 3.2 Consider then the left grid point is less than the target bmag.
!     Then, there are three cases in which there are matching points.
!     (1) The right grid point is equal to the target bmag, and the
!         slope at the right point is negative.
!     (2) The right grid point is greater than the target bmag.
!     (3) The right grid point is less than the target bmag,
!         and the interval is concave up, and the maximum in
!         this interval, which is guaranteed to exist, is greater
!         than the target bmag.
! 3.2.1 Consider when the right grid point is equal to the target bmag.
               if ((abs(bmagin(i+1)-bmagi) < epsilon(0.) &
                    .or. abs(thetain(i+1)-thetai) < epsilon(0.)) &
                    .and. bprime(i+1) < -bpknob) then
                  call add_setset_root (thetain(i),thetain(i+1),bmagi,iset,321)
                  cycle grid_interval
               end if
               if (bmagin(i+1) > bmagi) then
! 3.2.2 Consider when the right grid point is greater than the target bmag.
! 3.2.2.1 If this interval bounds the starting theta point, that is what
!         the target point is, and it is already collected.
                  if (thetai >= thetain(i) .and. thetai <= thetain(i+1)) then
                     cycle grid_interval
                  end if
! 3.2.2.2 Otherwise, find and collect the target point.
                  call add_setset_root (thetain(i),thetain(i+1),bmagi,iset,3222)
                  cycle grid_interval
               end if
! 3.2.3 Check if the grid interval is concave up.
! 3.2.3.1 If not, skip to the next interval.
               if (bprime(i) <= 0.0 .or. bprime(i+1) >= 0.0) cycle grid_interval
! 3.2.3.2 Consider the case where the starting theta point is within
!         this interval.
               if (thetai > thetain(i) .and. thetai < thetain(i+1)) then
! 3.2.3.2.1 If the starting point is an extremum, skip to next interval.
                  if (essentialstart(iset)) cycle grid_interval
! 3.2.3.2.2 Otherwise, the other target point is right of the starting
!           point if the slope is positive, and left if negative.
                  if (bprime(i) > 0.0) then
                     call add_setset_root (thetain(i+1),thetai,bmagi,iset,323221)
                  else
                     call add_setset_root (thetain(i),thetai,bmagi,iset,323222)
                  end if
                  cycle grid_interval
               end if
! 3.2.3.3 If the maximum within this interval is greater than the target bmag,
!         then there are two target points in this interval, one on each side
!         of the maximum.
! 3.2.3.3.1 If this interval is on the edge, there will not be any maximum.
               if (i == 1 .or. i == nbmag-1) cycle grid_interval
! 3.2.3.3.2 Find the maximum in this interval.
               if (abs(bextr(i)) < epsilon(0.)) then
                  print *, "gridgen4.f90:gg4collect:3.2.3.3.2:"," missing extremum"
                  print *, "iset,i:",iset,i
                  print *, "bmagi:",bmagi
                  print *, "bprimel,bprimer:", bprime(i),bprime(i+1)
                  print *, "thetain(i):", thetain(i)
                  print *, "thetain(i+1):", thetain(i+1)
                  print *, "thetai:", thetai
                  stop
               end if
! 3.2.3.3.2.1 If the maximum is less than the target bmag, skip to
!             the next interval.
               if (bextr(i) <= bmagi) cycle grid_interval
! 3.2.3.3.2.2 Collect the point left of the maximum.
               call add_setset_root (thetain(i),thetaextr(i),bmagi,iset,323322)
! 3.2.3.3.2.3 Collect the point right of the maximum.
               call add_setset_root (thetain(i+1),thetaextr(i),bmagi,iset,323323)
               cycle grid_interval
            else if (abs(bmagin(i)-bmagi) < epsilon(0.) &
                 .or. abs(thetain(i)-thetai) < epsilon(0.)) then
! 3.3 Consider when then left grid point is equal to the target bmag.
! 3.3.1 Add the point if it is not the starting grid point.
               if (abs(thetai-thetain(i)) > epsilon(0.)) then
                  call add_setset (thetain(i), iset, 331)
               end if
! 3.3.2 Check if there is another matching target bmag in the interval.
               if (bprime(i) > 0.0 .and. bmagin(i+1) < bmagi &
                    .and. i /= 1 .and. i /= nbmag-1) then
                  call add_setset_root (thetain(i+1),thetain(i),bmagi,iset,3321)
                  cycle grid_interval
               else if (bprime(i) < 0.0 .and. bmagin(i+1) > bmagi &
                    .and. i /= 1 .and. i /= nbmag-1) then
                  call add_setset_root (thetain(i),thetain(i+1),bmagi,iset,3322)
                  cycle grid_interval
               end if
            end if

         end do grid_interval
      end do starting_points

! 4 Compatibility with the old gg4collect.
      allocate (bmagset(nstart))
      allocate (essentialset(nstart))

      nset = nstart
      bmagset = bmagstart(1:nstart)
      essentialset = essentialstart(1:nstart)

    end subroutine gg4collect

    subroutine old_gg4collect
      implicit none
      integer :: iset, i, ii, imatch
      real :: thetai, bmagi
      real :: bprimer, bprimel, bprime0
      
      allocate (bmagset(nbmag*(2+npadd)))
      allocate (essentialset(nbmag*(2+npadd)))
      allocate (thetasetset(nbmag*(2+npadd)*2))
      allocate (ibmagsetset(nbmag*(2+npadd)*2))
      allocate (icollsetset(nbmag*(2+npadd)*2))
      
      starting_points: do iset = 1, nstart
         thetai = thetastart(iset)
         bmagi = bmagint(thetai)
         
! 1 For extrema, check previous sets to attach to.
         if (essentialstart(iset)) then
            do i = 1, nset
! 1.1 Check if the extremum should be attached to each previous set
               if (abs(bmagi-bmagset(i)) < epsknob) then
! 1.1.1 If the extremum does belong in this set, eliminate points
!       near the extremum from this set.
                  where (ibmagsetset(1:nsetset) == i &
                       .and. abs(thetasetset(1:nsetset)-thetai) < epsknob)
                     ibmagsetset(1:nsetset) = 0
                  end where
! 1.1.2 Attach the extremum to this set.
                  call add_setset (thetai, i, 112)
                  cycle starting_points
               end if
            end do
         end if

! 2 Start a new set.
         call add_set (bmagi, essentialstart(iset))
         call add_setset (thetai, nset, 2)
         
! 3 Check each original grid interval for matching bmag.
         grid_interval: do i = 1, nbmag-1
! 3.0.1 Stoopid problems near -pi.
            if (nset == 1 .and. i == 1) cycle grid_interval
            if (bmagin(i) > bmagi .and. abs(thetain(i)-thetai) > epsilon(0.)) then
! 3.1 Consider when the left grid point is greater than the target bmag.
!    Then, there are three cases in which there are matching points.
!    (1) The right grid point is equal to the target bmag, and the slope
!        at the right point is positive.
!    (2) The right gridpoint is less than the target bmag.
!    (3) The right grid point is greater than the target bmag,
!        and the interval is concave down, and the minimum in
!        this interval, which is guaranteed to exist, is less than
!        the target bmag.
               bprimer = bmagp(thetain(i+1))
               if ((abs(bmagin(i+1)-bmagi) < epsilon(0.) &
                    .or. abs(thetain(i+1)-thetai) < epsilon(0.)) &
                    .and. bprimer > 0.0) &
                    then
! 3.1.1 Consider when the right grid point is equal to the target bmag.
                  call add_setset_root (thetain(i+1),thetain(i),bmagi,nset,311)
                  cycle grid_interval
               end if
               if (bmagin(i+1) < bmagi) then
! 3.1.2 Consider when the right grid point is less than the target bmag.
! 3.1.2.1 If this interval bounds the starting theta point, that is what
!         the target point is, and it is already collected.
                  if (thetai >= thetain(i) .and. thetai <= thetain(i+1)) then
                     cycle grid_interval
                  end if
! 3.1.2.2 Otherwise, find and collect the target point.
                  call add_setset_root (thetain(i+1),thetain(i),bmagi,nset,3122)
                  cycle grid_interval
               end if
! 3.1.3 Check if the grid interval is concave down.
               bprimel = bmagp (thetain(i))
               bprimer = bmagp (thetain(i+1))
! 3.1.3.1 If not, skip to the next interval.
               if (bprimel >= 0.0 .or. bprimer <= 0.0) cycle grid_interval
! 3.1.3.2 Consider the case where the starting theta point is within
!         this interval.
               if (thetai > thetain(i) .and. thetai < thetain(i+1)) then
! 3.1.3.2.1 If the starting point is an extremum, skip to next interval.
                  if (essentialstart(iset)) cycle grid_interval
! 3.1.3.2.2 Otherwise, the other target point is right of the starting
!           point is the slope is negative, and left if positive.
                  bprime0 = bmagp(thetai)
                  if (bprime0 < 0) then
                     call add_setset_root (thetai,thetain(i+1),bmagi,nset,313221)
                  else
                     call add_setset_root (thetai,thetain(i),bmagi,nset,313222)
                  end if
                  cycle grid_interval
               end if
! 3.1.3.3 If the minimum within this interval is less than the target bmag,
!         then there are two target points in this interval, one on each side
!         of the minimum.
! 3.1.3.3.1 If this interval is on the edge, there will not be any minimum.
               if (i == 1 .or. i == nbmag-1) cycle grid_interval
! 3.1.3.3.2 Find the minimum in this interval.
               imatch = 0
               do ii = 1, nstart
                ! All essential points are at the beginning.
                  if (.not. essentialstart(ii)) exit
                  if (thetain(i) < thetastart(ii) &
                       .and. thetastart(ii) < thetain(i+1)) &
                       then
                     if (imatch /= 0) then
                        print *, "gridgen4.f90:gg4collect:3.1.3.3.2:", &
                             " multiple extrema in interval"
                        print *, "iset,i,ii,imatch:",iset,i,ii,imatch
                        print *, "bmagi:",bmagi
                        print *, "bprimel,bprimer:", bprimel,bprimer
                        print *, "thetain(i):", thetain(i)
                        print *, "thetain(i+1):", thetain(i+1)
                        print *, "thetai:", thetai
                        stop
                     end if
                     imatch = ii
                  end if
               end do
               if (imatch == 0) then
                  print *, "gridgen4.f90:gg4collect:3.1.3.3.2:", &
                       " missing extremum"
                  print *, "iset,i,ii,imatch:",iset,i,ii,imatch
                  print *, "bmagi:",bmagi
                  print *, "bprimel,bprimer:", bprimel,bprimer
                  print *, "thetain(i):", thetain(i)
                  print *, "thetain(i+1):", thetain(i+1)
                  print *, "thetai:", thetai
                  stop
               end if
! 3.1.3.3.2.1 If the minimum is greater than the target bmag, skip to
!             the next interval.
               if (bmagint(thetastart(imatch)) > bmagi) cycle grid_interval
! 3.1.3.3.2.2 Collect the point left of the minimum.
               call add_setset_root (thetastart(imatch),thetain(i),&
                    bmagi,nset,313322)
! 3.1.3.3.2.3 Collect the point right of the minimum.
               call add_setset_root (thetastart(imatch),thetain(i+1),&
                    bmagi,nset,313323)
               cycle grid_interval
            else if (bmagin(i) < bmagi .and. abs(thetain(i)-thetai) > epsilon(0.)) then
! 3.2 Consider then the left grid point is less than the target bmag.
!     Then, there are three cases in which there are matching points.
!     (1) The right grid point is equal to the target bmag, and the
!         slope at the right point is negative.
!     (2) The right grid point is greater than the target bmag.
!     (3) The right grid point is less than the target bmag,
!         and the interval is concave up, and the maximum in
!         this interval, which is guaranteed to exist, is greater
!         than the target bmag.
               bprimer = bmagp(thetain(i+1))
! 3.2.1 Consider when the right grid point is equal to the target bmag.
               if ((abs(bmagin(i+1)-bmagi) < epsilon(0.) &
                    .or. abs(thetain(i+1)-thetai) < epsilon(0.)) &
                    .and. bprimer < 0.0) &
                    then
                  call add_setset_root (thetain(i),thetain(i+1),bmagi,nset,321)
                  cycle grid_interval
               end if
               if (bmagin(i+1) > bmagi) then
! 3.2.2 Consider when the right grid point is greater than the target bmag.
! 3.2.2.1 If this interval bounds the starting theta point, that is what
!         the target point is, and it is already collected.
                  if (thetai >= thetain(i) .and. thetai <= thetain(i+1)) then
                     cycle grid_interval
                  end if
! 3.2.2.2 Otherwise, find and collect the target point.
                  call add_setset_root (thetain(i),thetain(i+1),bmagi,nset,3222)
                  cycle grid_interval
               end if
! 3.2.3 Check if the grid interval is concave up.
               bprimel = bmagp(thetain(i))
               bprimer = bmagp(thetain(i+1))
! 3.2.3.1 If not, skip to the next interval.
               if (bprimel <= 0.0 .or. bprimer >= 0.0) cycle grid_interval
! 3.2.3.2 Consider the case where the starting theta point is within
!         this interval.
               if (thetai > thetain(i) .and. thetai < thetain(i+1)) then
! 3.2.3.2.1 If the starting point is an extremum, skip to next interval.
                  if (essentialstart(iset)) cycle grid_interval
! 3.2.3.2.2 Otherwise, the other target point is right of the starting
!           point if the slope is positive, and left if negative.
                  bprime0 = bmagp(thetai)
                  if (bprime0 > 0.0) then
                     call add_setset_root (thetain(i+1),thetai,bmagi,nset,323221)
                  else
                     call add_setset_root (thetain(i),thetai,bmagi,nset,323222)
                  end if
                  cycle grid_interval
               end if
! 3.2.3.3 If the maximum within this interval is greater than the target bmag,
!         then there are two target points in this interval, one on each side
!         of the maximum.
! 3.2.3.3.1 If this interval is on the edge, there will not be any maximum.
               if (i == 1 .or. i == nbmag-1) cycle grid_interval
! 3.2.3.3.2 Find the maximum in this interval.
               imatch = 0
               do ii = 1, nstart
                ! All essential points are at the beginning.
                  if (.not. essentialstart(ii)) exit
                  if (thetain(i) < thetastart(ii) &
                       .and. thetastart(ii) < thetain(i+1)) &
                       then
                     if (imatch /= 0) then
                        print *, "gridgen4.f90:gg4collect:3.2.3.3.2:", &
                             " multiple extrema in interval"
                        print *, "iset,i,ii,imatch:",iset,i,ii,imatch
                        print *, "bmagi:",bmagi
                        print *, "bprimel,bprimer:", bprimel,bprimer
                        print *, "thetain(i):", thetain(i)
                        print *, "thetain(i+1):", thetain(i+1)
                        print *, "thetai:", thetai
                        stop
                     end if
                     imatch = ii
                  end if
               end do
               if (imatch == 0) then
                  print *, "gridgen4.f90:gg4collect:3.2.3.3.2:", &
                       " missing extremum"
                  print *, "iset,i,ii,imatch:",iset,i,ii,imatch
                  print *, "bmagi:",bmagi
                  print *, "bprimel,bprimer:", bprimel,bprimer
                  print *, "thetain(i):", thetain(i)
                  print *, "thetain(i+1):", thetain(i+1)
                  print *, "thetai:", thetai
                  stop
               end if
! 3.2.3.3.2.1 If the maximum is less than the target bmag, skip to
!             the next interval.
               if (bmagint(thetastart(imatch)) <= bmagi) cycle grid_interval
! 3.2.3.3.2.2 Collect the point left of the maximum.
               call add_setset_root (thetain(i),thetastart(imatch),&
                    bmagi,nset,323322)
! 3.2.3.3.2.3 Collect the point right of the maximum.
               call add_setset_root (thetain(i+1),thetastart(imatch),&
                    bmagi,nset,323323)
               cycle grid_interval
            else if (abs(bmagin(i)-bmagi) < epsilon(0.) &
                 .or. abs(thetain(i)-thetai) < epsilon(0.)) then
! 3.3 Consider when then left grid point is equal to the target bmag.
! 3.3.1 Add the point if it is not the starting grid point.
               if (abs(thetai-thetain(i)) > epsilon(0.)) then
                  call add_setset (thetain(i), nset, 331)
! 3.3.2 Check if there is another matching target bmag in the interval.
                  bprime0 = bmagp(thetain(i))
                  if (bprime0 > 0.0 .and. bmagin(i+1) < bmagi &
                       .and. i /= 1 .and. i /= nbmag-1) &
                       then
                     call add_setset_root (thetain(i+1),thetain(i), &
                          bmagi,nset,3321)
                     cycle grid_interval
                  else if (bprime0 < 0.0 .and. bmagin(i+1) > bmagi &
                       .and. i /= 1 .and. i /= nbmag-1) &
                       then
                     call add_setset_root (thetain(i),thetain(i+1), &
                          bmagi,nset,3322)
                     cycle grid_interval
                  end if
               end if
            end if
         end do grid_interval
      end do starting_points
    end subroutine old_gg4collect

    subroutine add_set (bmag, essential)
      implicit none
      real, intent (in) :: bmag
      logical, intent (in) :: essential

      real, allocatable, dimension (:) :: bmagtmp
      logical, allocatable, dimension (:) :: essentialtmp
      
      if (nset >= size(bmagset)) then
         allocate (bmagtmp(nset))
         allocate (essentialtmp(nset))
         bmagtmp = bmagset(1:nset)
         essentialtmp = essentialset(1:nset)
         deallocate (bmagset, essentialset)
         allocate (bmagset(nset+nbmag))
         allocate (essentialset(nset+nbmag))
         bmagset(1:nset) = bmagtmp
         essentialset(1:nset) = essentialtmp
         deallocate (bmagtmp, essentialtmp)
      end if
      nset = nset + 1
      bmagset(nset) = bmag
      essentialset(nset) = essential
    end subroutine add_set

    subroutine add_setset (theta, ibmag, icoll)
      implicit none
      real, intent (in) :: theta
      integer, intent (in) :: ibmag, icoll

      real, allocatable, dimension (:) :: thetatmp
      integer, allocatable, dimension (:) :: ibmagtmp, icolltmp

      if (nsetset >= size(thetasetset)) then
         allocate (thetatmp(nsetset))
         allocate (ibmagtmp(nsetset), icolltmp(nsetset))
         thetatmp = thetasetset(1:nsetset)
         ibmagtmp = ibmagsetset(1:nsetset)
         icolltmp = icollsetset(1:nsetset)
         deallocate (thetasetset, ibmagsetset, icollsetset)
         allocate (thetasetset(nsetset+nbmag*(2+npadd)))
         allocate (ibmagsetset(nsetset+nbmag*(2+npadd)))
         allocate (icollsetset(nsetset+nbmag*(2+npadd)))
         thetasetset(1:nsetset) = thetatmp
         ibmagsetset(1:nsetset) = ibmagtmp
         icollsetset(1:nsetset) = icolltmp
         deallocate (thetatmp, ibmagtmp, icolltmp)
      end if
      nsetset = nsetset + 1
      thetasetset(nsetset) = theta
      ibmagsetset(nsetset) = ibmag
      icollsetset(nsetset) = icoll
      if (any(ibmagsetset(1:nsetset-1) == ibmag .and. &
           abs(thetasetset(1:nsetset-1)-theta) < epsknob)) &
           then
         ibmagsetset(nsetset) = 0
      end if
    end subroutine add_setset
    
    subroutine add_setset_root (thetal, thetag, bmagi, ibmag, icoll)
      implicit none
      real, intent (in) :: thetal, thetag, bmagi
      integer, intent (in) :: ibmag, icoll

      real :: thl, thg, theta0, bmag0
      integer :: i
      
      thl = thetal
      thg = thetag
      do i = 1, maxiter
         theta0 = 0.5*(thl+thg)
         bmag0 = bmagint(theta0)
         if (bmag0 > bmagi) then
            thg = theta0
         else if (bmag0 < bmagi) then
            thl = theta0
         else
            exit
         end if
      end do
      call add_setset (theta0, ibmag, icoll)
    end subroutine add_setset_root
    
    subroutine gg4sort
      implicit none
      allocate (ithetasort(nsetset), ibmagsort(nset))
      call heapsort (nsetset, thetasetset(1:nsetset), ithetasort)
      call heapsort (nset, bmagset(1:nset), ibmagsort)
    end subroutine gg4sort

    subroutine heapsort (n, v, index)
      implicit none
      integer, intent (in) :: n
      real, dimension (n), intent (in) :: v
      integer, dimension (n), intent (out) :: index
      integer :: i, j, l, ir, ira
      
      index = (/ (i, i=1,n) /)
      
      l = n/2+1
      ir = n
      
      do
         if (l > 1) then
            l = l - 1
            ira = index(l)
         else
            ira = index(ir)
            index(ir) = index(1)
            ir = ir - 1
            if (ir == 1) then
               index(1) = ira
               return
            end if
         end if
         i = l
         j = l + l
         do while (j <= ir)
            if (j < ir) then
               if (v(index(j)) < v(index(j+1))) j = j + 1
            end if
            if (v(ira) < v(index(j))) then
               index(i) = index(j)
               i = j
               j = j + j
            else
               j = ir + 1
            end if
         end do
         index(i) = ira
      end do
    end subroutine heapsort

    subroutine gg4metrics
      implicit none
      integer :: i
      
      allocate (thetares(nset), alambdares(nset))
      
      do i = 1, nset
         if (abs(bmagset(i)) > epsilon(0.)) then
            call get_thetares (i, thetares(i))
            call get_lambdares (i, alambdares(i))
         else
            thetares(i) = 0.0
            alambdares(i) = 0.0
         end if
      end do
    end subroutine gg4metrics
    
    subroutine get_thetares (iset, thetares)
      implicit none
      integer, intent (in) :: iset
      real, intent (out) :: thetares

      integer :: i, ileft, iright
      real :: dthetal, dthetar
      real :: res
      integer :: npts

! 1 Return large value for essential sets with odd number of points or
!   for set containing -pi.
      if (essentialset(iset)) then
! 1.1 Set containing -pi.
         if (iset == 1) then
            thetares = 2e20
            return
         end if
! 1.2 Sets with odd number of points.
! 1.2.1 Count points.
!       npts = count(ibmagsetset == iset)
         npts = count(ibmagsetset(1:nsetset) == iset)
         if (mod(npts,2) == 1) then
            thetares = 1e20
            return
         end if
         res = extrknob/real(npts*npts)
      else
         res = 0.0
      end if

! 2 Look through all points for points in the set.
      npts = 0
      points_in_set: do i = 1, nsetset
         if (ibmagsetset(ithetasort(i)) /= iset) cycle
! 2.1 Look for the point to the left.
         ileft = i
         do
            ileft = ileft - 1
            if (ileft < 1) cycle points_in_set
            if (ibmagsetset(ithetasort(ileft)) == 0) cycle
            if (abs(bmagset(ibmagsetset(ithetasort(ileft)))) > epsilon(0.)) exit
         end do
! 2.2 Look for the point to the right.
         iright = i
         do
            iright = iright + 1
            if (iright > nsetset) cycle points_in_set
            if (ibmagsetset(ithetasort(iright)) == 0) cycle
            if (abs(bmagset(ibmagsetset(ithetasort(iright)))) > epsilon(0.)) exit
         end do
! 2.3 Add contribution from this interval.
         dthetal=abs(thetasetset(ithetasort(i))-thetasetset(ithetasort(ileft)))
         dthetar=abs(thetasetset(ithetasort(i))-thetasetset(ithetasort(iright)))
         res = res + tfunc(thetasetset(ithetasort(i)))*dthetal*dthetar &
              /(dthetal+dthetar+1e-20)
         npts = npts + 1
      end do points_in_set
      if (npts > 0) then
         thetares = res/real(npts)
      else
         thetares = res
      end if
    end subroutine get_thetares
    
    subroutine get_lambdares (iset, alambdares)
      implicit none
      integer, intent (in) :: iset
      real, intent (out) :: alambdares

      integer :: i, iplus, iminus
      real :: al, alplus, alminus, dalplus, dalminus
      real :: res
      integer :: npts
      
      al = 1.0/bmagset(iset)
      do i = 1, nset
! 1 Look for target lambda
         if (ibmagsort(i) == iset) then
! 2 Look for bordering lambdas.
            alplus = 0.0
            do iplus = i-1, 1, -1
               if (ibmagsort(iplus) == 0) cycle
               if (abs(bmagset(ibmagsort(iplus))) > epsilon(0.)) then
                  alplus = 1.0/bmagset(ibmagsort(iplus))
                  exit
               end if
            end do
            alminus = 0.0
            do iminus = i+1, nset
               if (ibmagsort(iminus) == 0) cycle
               if (abs(bmagset(ibmagsort(iminus))) > epsilon(0.)) then
                  alminus = 1.0/bmagset(ibmagsort(iminus))
                  exit
               end if
            end do
            exit
         end if
      end do

! 3 Add up contributions to the result.
      res = 0.0
      npts = 0
      do i = 1, nbmag
         dalplus = abs(sqrt(max(1.0-alplus*bmagin(i),0.0)) &
              -sqrt(max(1.0-al*bmagin(i),0.0)))
         dalminus = abs(sqrt(max(1.0-alminus*bmagin(i),0.0)) &
              -sqrt(max(1.0-al*bmagin(i),0.0)))
         if (abs(dalplus+dalminus) > epsilon(0.)) then
            npts = npts + 1
            res = res + dalplus*dalminus/(dalplus+dalminus+1e-20)
         end if
      end do
      if (npts /= 0) then
         alambdares = res/real(npts)
      else
         alambdares = res
      end if
    end subroutine get_lambdares

    subroutine gg4remove
      implicit none
      integer :: idel, i
      integer :: ntheta_left, nlambda_left
      integer, dimension (nset) :: work
      
      nlambda_left = nset
      ntheta_left = nsetset
      do i = 1, nset
         if (.not. any(ibmagsetset(1:nsetset) == i)) then
            nlambda_left = nlambda_left - 1
         end if
      end do
! 1 Find the set with the minimum resolution metric.
      do while (ntheta_left > ntheta .or. nlambda_left > nbset)
         work(1:1) = minloc(thetares(1:nset)+alknob*alambdares(1:nset), &
              abs(bmagset(1:nset)) > epsilon(0.))
         idel = work(1)
         
! 2 Delete the set just found.
         if (idel == 0) then
            print *, "gridgen4.f:gg4remove:2: This cannot happen."
            stop
         end if
         call delete_set (idel, work, ntheta_left)
         nlambda_left = nlambda_left - 1
      end do
    end subroutine gg4remove
    
    subroutine delete_set (idel, work, ntheta_left)
      implicit none
      integer, intent (in) :: idel
      integer, intent (in out), dimension (nset) :: work
      integer, intent (out) :: ntheta_left

      integer :: i, j
      
      work = 0
! 1 Mark neighboring lambda sets to be recalculated.
      do i = 1, nset
         if (ibmagsort(i) == idel) then
            do j = i-1, 1, -1
               if (abs(bmagset(ibmagsort(j))) > epsilon(0.)) then
                  work(ibmagsort(j)) = ibmagsort(j)
                  exit
               end if
            end do
            do j = i+1, nset
               if (abs(bmagset(ibmagsort(j))) > epsilon(0.)) then
                  work(ibmagsort(j)) = ibmagsort(j)
                  exit
               end if
            end do
            exit
         end if
      end do

! 2 Mark lambda sets with neighboring theta points to have their resolution
!   metrics recalculated, counting the remaining theta points.
      ntheta_left = 0
      do i = 1, nsetset
         if (ibmagsetset(ithetasort(i)) == 0) cycle
         if (abs(bmagset(ibmagsetset(ithetasort(i)))) < epsilon(0.)) cycle
         if (ibmagsetset(ithetasort(i)) /= idel) then
            ntheta_left = ntheta_left + 1
            cycle
         end if
! 2.1 Found point to be deleted.
! 2.1.1 Mark set of neighboring point to the left to be recalculated.
         do j = i-1, 1, -1
            if (ibmagsetset(ithetasort(j)) == idel) cycle
            if (ibmagsetset(ithetasort(j)) == 0) cycle
            if (abs(bmagset(ibmagsetset(ithetasort(j)))) < epsilon(0.)) cycle
            work(ibmagsetset(ithetasort(j))) = ibmagsetset(ithetasort(j))
            exit
         end do
! 2.1.2 Mark set of neighboring point to the right to be recalculated.
         do j = i+1, nsetset
            if (ibmagsetset(ithetasort(j)) == idel) cycle
            if (ibmagsetset(ithetasort(j)) == 0) cycle
            if (abs(bmagset(ibmagsetset(ithetasort(j)))) < epsilon(0.)) cycle
            work(ibmagsetset(ithetasort(j))) = ibmagsetset(ithetasort(j))
            exit
         end do
      end do

! 3 Delete this set.
      bmagset(idel) = 0.0

! 4 Recalculate resolution metric for affected sets.
      do i = 1, nset
         if (work(i) /= 0) then
            call get_thetares (work(i), thetares(work(i)))
            call get_lambdares (work(i), alambdares(work(i)))
         end if
      end do
    end subroutine delete_set
    
    subroutine gg4results (iperiod)
      implicit none
      integer, intent (in) :: iperiod
      integer :: i, n
      
      n = 0
      do i = 1, nsetset
         if (ibmagsetset(ithetasort(i)) == 0) cycle
         if (abs(bmagset(ibmagsetset(ithetasort(i)))) < epsilon(0.)) cycle
         n = n + 1
         thetagrid(n) = thetasetset(ithetasort(i))
         bmaggrid(n) = bmagset(ibmagsetset(ithetasort(i)))
      end do
      
! Check point at +pi*iperiod
      if (abs(bmaggrid(n)-bmaggrid(1)) > epsilon(0.)) then
         n = n + 1
         thetagrid(n) = pi*iperiod
         bmaggrid(n) = bmaggrid(1)
      end if
      nthetaout = n

      n = 0
      do i = nset, 1, -1
         if (abs(bmagset(ibmagsort(i))) < epsilon(0.)) cycle
         n = n + 1
         bset(n) = bmagset(ibmagsort(i))
      end do
      nlambdaout = n
    end subroutine gg4results

    subroutine gg4debug (n, x, y, label)
      implicit none
      integer, intent (in) :: n
      real, dimension (:), intent (in) :: x, y
      character(*), intent (in) :: label

      integer :: i
      
      if (debug_output) then
         write (unit=debug_unit, fmt="('#',a)") label
         do i = 1, n
            write (unit=debug_unit, fmt="(i5,1x,2(g19.12,1x))") i, x(i), y(i)
         end do
      end if
    end subroutine gg4debug
    
    subroutine gg4debugi (n, x, i1, i2, label)
      implicit none
      integer, intent (in) :: n
      real, dimension (:), intent (in) :: x
      integer, dimension (:), intent (in) :: i1, i2
      character(*), intent (in) :: label

      integer :: i
      
      if (debug_output) then
         write (unit=debug_unit, fmt="('#',a)") label
         do i = 1, n
            write (unit=debug_unit, fmt="(i5,1x,g22.15,1x,i4,1x,i10)") &
                 i, x(i), i1(i), i2(i)
         end do
      end if
    end subroutine gg4debugi
    
    function bmagint (theta)
      use splines, only: fitp_curvp2
      implicit none
      real, intent (in) :: theta
      real :: bmagint
      bmagint = fitp_curvp2(theta,nbmag-1,thetain,bmagin,npi,bmagspl,tension)
    end function bmagint
    
    function bmagp (theta)
      implicit none
      real, intent (in) :: theta
      real :: bmagp
      real, parameter :: diff = 1e-5

      real :: th1, th2, bm1, bm2
      th1 = theta + 0.5*diff
      th2 = theta - 0.5*diff
      bm1 = bmagint (th1)
      bm2 = bmagint (th2)
      bmagp = (bm1-bm2)/diff
    end function bmagp
    
    real function tfunc (theta)
      implicit none
      real :: theta
      tfunc = 1.0 + deltaw*(1.0/(1.0+(theta-thetamax)**2/widthw**2) &
           +1.0/(1.0+(theta+thetamax)**2/widthw**2))
    end function tfunc
    
  end subroutine gridgen4_2

  subroutine dsmooth (n, xin, yin, yout)
    implicit none

    integer n
    real, dimension (:), intent (in) ::  xin, yin
    real, dimension (:), intent (out) :: yout

    real, dimension (:), allocatable :: yp, ypp, xypp
   
    real ypb, yppb, xb
    integer iter, i

    integer, parameter ::  maxiter = 20
    real, parameter ::  diff = 0.4

    allocate (yp(n), ypp(n), xypp(n))

    yout = yin

    do iter=1,maxiter
       do i=1,n-1
          yp(i) = (yout(i+1)-yout(i))/(xin(i+1)-xin(i))
       enddo
       do i=2,n-1
          xypp(i) = yout(i+1)-2.0*yout(i)+yout(i-1)
          ypp(i) = 2.0*(yp(i)-yp(i-1))/(xin(i+1)-xin(i-1))
       enddo
       ypb = 0.0
       yppb = 0.0
       xb = xin(n-2)-xin(2)
       do i=2,n-2
          ypb = ypb+abs(yp(i)/xb)
          yppb = yppb+abs(ypp(i)/xb)
       enddo
       if (ypb/yppb .gt. 2.0*xb/float(n)) goto 1
       do i=2,n-1
          yout(i) = yout(i) + diff*xypp(i)
       enddo
    enddo
1   continue

    deallocate (yp, ypp, xypp)

  end subroutine dsmooth

  subroutine d2smooth (n, xin, yin, yout)
    implicit none

    integer n
    real, dimension (:), intent (in) ::  xin, yin
    real, dimension (:), intent (out) :: yout

    real, dimension (:), allocatable :: yp, ypp, yppp, xypp
    real ypb, yppb, ypppb, xb
    integer iter,i
    integer, parameter ::  maxiter=20
    real, parameter :: diff = 0.4

    allocate (yp(n), ypp(n), yppp(n), xypp(n))

    yout = yin

    do iter=1,maxiter
       do i=1,n-1
          yp(i) = (yout(i+1)-yout(i))/(xin(i+1)-xin(i))
       enddo
       do i=2,n-1
          xypp(i) = yout(i+1)-2.0*yout(i)+yout(i-1)
          ypp(i) = 2.0*(yp(i)-yp(i-1))/(xin(i+1)-xin(i-1))
       enddo
       do i=2,n-2
          yppp(i) = (ypp(i+1)-ypp(i))/(xin(i+1)-xin(i))
       enddo
       ypb = 0.0
       yppb = 0.0
       ypppb = 0.0
       xb = xin(n-2)-xin(2)
       do i=2,n-2
          ypb = ypb+abs(yp(i)/xb)
          yppb = yppb+abs(ypp(i)/xb)
          ypppb = ypppb+abs(yppp(i)/xb)
       enddo
       if (ypb/yppb .gt. 4.0*xb/float(n) .and. yppb/ypppb .gt. 4.0*xb/float(n)) goto 1
       do i=2,n-1
          yout(i) = yout(i) + diff*xypp(i)
       enddo
    enddo
1   continue

    deallocate (yp, ypp, yppp, xypp)
    
  end subroutine d2smooth

  subroutine smooth (n, xin, yin, var, yout, ifail)
    implicit none

    integer n
    real, dimension (:), intent (in) :: xin, yin
    real, dimension (:), intent (out) :: yout
    real :: var

    integer ifail

! these next arrays should be double precision, which we usually get with 
! compiler options
    real, dimension (1) :: se
    real, dimension (:), allocatable :: x, f, y, df
    real, dimension (:,:), allocatable :: wk, c
    real :: dvar

    allocate (x(n), f(n), y(n), df(n), c(n, 3))
    allocate (wk(n+2,7))

    ifail = 0
    x = xin
    f = yin
    df = 1.
    dvar=var
    call cubgcv (x,f,df,n,y,c,n,dvar,0,se,wk,ifail)
    var=dvar
    yout = y

    deallocate (x, f, y, df, c, wk)

  end subroutine smooth

!     algorithm 642 collected algorithms from acm.
!     algorithm appeared in acm-trans. math. software, vol.12, no. 2,
!     jun., 1986, p. 150.
!   subroutine name     - cubgcv
!
!--------------------------------------------------------------------------
!
!   computer            - vax/double
!
!   author              - m.f.hutchinson
!                         csiro division of mathematics and statistics
!                         p.o. box 1965
!                         canberra, act 2601
!                         australia
!
!   latest revision     - 15 august 1985
!
!   purpose             - cubic spline data smoother
!
!   usage               - call cubgcv (x,f,df,n,y,c,ic,var,job,se,wk,ier)
!
!   arguments    x      - vector of length n containing the
!                           abscissae of the n data points
!                           (x(i),f(i)) i=1..n. (input) x
!                           must be ordered so that
!                           x(i) .lt. x(i+1).
!                f      - vector of length n containing the
!                           ordinates (or function values)
!                           of the n data points (input).
!                df     - vector of length n. (input/output)
!                           df(i) is the relative standard deviation
!                           of the error associated with data point i.
!                           each df(i) must be positive.  the values in
!                           df are scaled by the subroutine so that
!                           their mean square value is 1, and unscaled
!                           again on normal exit.
!                           the mean square value of the df(i) is returned
!                           in wk(7) on normal exit.
!                           if the absolute standard deviations are known,
!                           these should be provided in df and the error
!                           variance parameter var (see below) should then
!                           be set to 1.
!                           if the relative standard deviations are unknown,
!                           set each df(i)=1.
!                n      - number of data points (input).
!                           n must be .ge. 3.
!                y,c    - spline coefficients. (output) y
!                           is a vector of length n. c is
!                           an n-1 by 3 matrix. the value
!                           of the spline approximation at t is
!                           s(t)=((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
!                           where x(i).le.t.lt.x(i+1) and
!                           d = t-x(i).
!                ic     - row dimension of matrix c exactly
!                           as specified in the dimension
!                           statement in the calling program. (input)
!                var    - error variance. (input/output)
!                           if var is negative (i.e. unknown) then
!                           the smoothing parameter is determined
!                           by minimizing the generalized cross validation
!                           and an estimate of the error variance is
!                           returned in var.
!                           if var is non-negative (i.e. known) then the
!                           smoothing parameter is determined to minimize
!                           an estimate, which depends on var, of the true
!                           mean square error, and var is unchanged.
!                           in particular, if var is zero, then an
!                           interpolating natural cubic spline is calculated.
!                           var should be set to 1 if absolute standard
!                           deviations have been provided in df (see above).
!                job    - job selection parameter. (input)
!                         job = 0 should be selected if point standard error
!                           estimates are not required in se.
!                         job = 1 should be selected if point standard error
!                           estimates are required in se.
!                se     - vector of length n containing bayesian standard
!                           error estimates of the fitted spline values in y.
!                           se is not referenced if job=0. (output)
!                wk     - work vector of length 7*(n + 2). on normal exit the
!                           first 7 values of wk are assigned as follows:-
!
!                           wk(1) = smoothing parameter (= rho/(rho + 1))
!                           wk(2) = estimate of the number of degrees of
!                                   freedom of the residual sum of squares
!                           wk(3) = generalized cross validation
!                           wk(4) = mean square residual
!                           wk(5) = estimate of the true mean square error
!                                   at the data points
!                           wk(6) = estimate of the error variance
!                           wk(7) = mean square value of the df(i)
!
!                           if wk(1)=0 (rho=0) an interpolating natural cubic
!                           spline has been calculated.
!                           if wk(1)=1 (rho=infinite) a least squares
!                           regression line has been calculated.
!                           wk(2) is an estimate of the number of degrees of
!                           freedom of the residual which reduces to the
!                           usual value of n-2 when a least squares regression
!                           line is calculated.
!                           wk(3),wk(4),wk(5) are calculated with the df(i)
!                           scaled to have mean square value 1.  the
!                           unscaled values of wk(3),wk(4),wk(5) may be
!                           calculated by dividing by wk(7).
!                           wk(6) coincides with the output value of var if
!                           var is negative on input.  it is calculated with
!                           the unscaled values of the df(i) to facilitate
!                           comparisons with a priori variance estimates.
!
!                ier    - error parameter. (output)
!                         terminal error
!                           ier = 129, ic is less than n-1.
!                           ier = 130, n is less than 3.
!                           ier = 131, input abscissae are not
!                             ordered so that x(i).lt.x(i+1).
!                           ier = 132, df(i) is not positive for some i.
!                           ier = 133, job is not 0 or 1.
!
!   precision/hardware  - double
!
!   required routines   - spint1,spfit1,spcof1,sperr1
!
!   remarks      the number of arithmetic operations required by the
!                subroutine is proportional to n.  the subroutine
!                uses an algorithm developed by m.f. hutchinson and
!                f.r. de hoog, 'smoothing noisy data with spline
!                functions', numer. math. (in press)
!
!-----------------------------------------------------------------------
!
  subroutine cubgcv(x,f,df,n,y,c,ic,var,job,se,wk,ier)
!
!---specifications for arguments---
    integer n,ic,job,ier
!
! All real variables should be double precision: 
!
    real, dimension (:) :: x, f, df, y, se
    real, dimension (:,:) :: c
    real, dimension (0:n+1,7) :: wk
    real :: var
!
!---specifications for local variables---
    real :: delta,err,gf1,gf2,gf3,gf4,r1,r2,r3,r4,avh,avdf,avar,stat(6),p,q
    real :: ratio = 2.
    real :: tau = 1.618033989
    real :: zero = 0.
    real :: one = 1.

    integer :: i
!
!---initialize---
    ier = 133
    if (job.lt.0 .or. job.gt.1) go to 140
    call spint1(x,avh,f,df,avdf,n,y,c,ic,wk,wk(0,4),ier)
    if (ier.ne.0) go to 140
    avar = var
    if (var.gt.zero) avar = var*avdf*avdf
!
!---check for zero variance---
    if (abs(var)>epsilon(0.)) go to 10
    r1 = zero
    go to 90
!
!---find local minimum of gcv or the expected mean square error---
10  r1 = one
    r2 = ratio*r1
    call spfit1(x,avh,df,n,r2,p,q,gf2,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
20  call spfit1(x,avh,df,n,r1,p,q,gf1,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
    if (gf1.gt.gf2) go to 30
!
!---exit if p zero---
    if (p.le.zero) go to 100
    r2 = r1
    gf2 = gf1
    r1 = r1/ratio
    go to 20

30  r3 = ratio*r2
40  call spfit1(x,avh,df,n,r3,p,q,gf3,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
    if (gf3.gt.gf2) go to 50
!
!---exit if q zero---
    if (q.le.zero) go to 100
    r2 = r3
    gf2 = gf3
    r3 = ratio*r3
    go to 40

50  r2 = r3
    gf2 = gf3
    delta = (r2-r1)/tau
    r4 = r1 + delta
    r3 = r2 - delta
    call spfit1(x,avh,df,n,r3,p,q,gf3,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
    call spfit1(x,avh,df,n,r4,p,q,gf4,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
!
!---golden section search for local minimum---
60  if (gf3.gt.gf4) go to 70
    r2 = r4
    gf2 = gf4
    r4 = r3
    gf4 = gf3
    delta = delta/tau
    r3 = r2 - delta
    call spfit1(x,avh,df,n,r3,p,q,gf3,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
    go to 80
    
70  r1 = r3
    gf1 = gf3
    r3 = r4
    gf3 = gf4
    delta = delta/tau
    r4 = r1 + delta
    call spfit1(x,avh,df,n,r4,p,q,gf4,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
80  err = (r2-r1)/ (r1+r2)
    if (err*err+one.gt.one .and. err.gt.1.0e-6) go to 60
    r1 = (r1+r2)*0.5
!
!---calculate spline coefficients---
90  call spfit1(x,avh,df,n,r1,p,q,gf1,avar,stat,y,c,ic,wk(:,1:3),wk(:,4:5),wk(:,6),wk(:,7))
100 call spcof1(x,avh,f,df,n,p,q,y,c,ic,wk(:,6),wk(:,7))
!
!---optionally calculate standard error estimates---
    if (var.ge.zero) go to 110
    avar = stat(6)
    var = avar/ (avdf*avdf)
110 if (job.eq.1) call sperr1(x,avh,df,n,wk,p,avar,se)
!
!---unscale df---
    df = df * avdf
!
!--put statistics in wk---
    do i = 0,5
       wk(i,1) = stat(i+1)
    end do
    wk(5,1) = stat(6)/ (avdf*avdf)
    wk(6,1) = avdf*avdf
    go to 150
!
!---check for error condition---
140 continue
!     if (ier.ne.0) continue
150 return
  end subroutine cubgcv

  subroutine spint1(x,avh,y,dy,avdy,n,a,c,ic,r,t,ier)
!
! initializes the arrays c, r and t for one dimensional cubic
! smoothing spline fitting by subroutine spfit1.  the values
! df(i) are scaled so that the sum of their squares is n
! and the average of the differences x(i+1) - x(i) is calculated
! in avh in order to avoid underflow and overflow problems in
! spfit1.
!
! subroutine sets ier if elements of x are non-increasing,
! if n is less than 3, if ic is less than n-1 or if dy(i) is
! not positive for some i.
!
!---specifications for arguments---
    integer n,ic,ier
! All variables should be double precision
    real, dimension (:) :: x, y, dy, a
    real, dimension (:,:) :: c
    real, dimension (0:n+1,3) :: r
    real, dimension (0:n+1,2) :: t
    real :: avh, avdy

!
!---specifications for local variables---
    integer i
    real :: e,f,g,h
    real :: zero = 0.

!
!---initialization and input checking---
    ier = 0
    if (n.lt.3) go to 60
    if (ic.lt.n-1) go to 70
!
!---get average x spacing in avh---
    g = zero
    do i = 1,n - 1
       h = x(i+1) - x(i)
       if (h.le.zero) go to 80
       g = g + h
    end do
    avh = g/ (n-1)
!
!---scale relative weights---
    g = zero
    do i = 1,n
       if (dy(i).le.zero) go to 90
       g = g + dy(i)*dy(i)
    end do
    avdy = sqrt(g/n)
    
    dy = dy / avdy

!
!---initialize h,f---
    h = (x(2)-x(1))/avh
    f = (y(2)-y(1))/h
!
!---calculate a,t,r---
    do i = 2,n - 1
       g = h
       h = (x(i+1)-x(i))/avh
       e = f
       f = (y(i+1)-y(i))/h
       a(i) = f - e
       t(i,1) = 2.0d0* (g+h)/3.0d0
       t(i,2) = h/3.0d0
       r(i,3) = dy(i-1)/g
       r(i,1) = dy(i+1)/h
       r(i,2) = -dy(i)/g - dy(i)/h
    end do
!
!---calculate c = r'*r---
    r(n,2) = zero
    r(n,3) = zero
    r(n+1,3) = zero
    do i = 2,n - 1
       c(i,1) = r(i,1)*r(i,1) + r(i,2)*r(i,2) + r(i,3)*r(i,3)
       c(i,2) = r(i,1)*r(i+1,2) + r(i,2)*r(i+1,3)
       c(i,3) = r(i,1)*r(i+2,3)
    end do
    return
!
!---error conditions---
60  ier = 130
    return
    
70  ier = 129
    return

80  ier = 131
    return

   90 ier = 132
    return
  end subroutine spint1

  subroutine spfit1(x,avh,dy,n,rho,p,q,fun,var,stat,a,c,ic,r,t,u,v)
!
! fits a cubic smoothing spline to data with relative
! weighting dy for a given value of the smoothing parameter
! rho using an algorithm based on that of c.h. reinsch (1967),
! numer. math. 10, 177-183.
!
! the trace of the influence matrix is calculated using an
! algorithm developed by m.f.hutchinson and f.r.de hoog (numer.
! math., in press), enabling the generalized cross validation
! and related statistics to be calculated in order n operations.
!
! the arrays a, c, r and t are assumed to have been initialized
! by the subroutine spint1.  overflow and underflow problems are
! avoided by using p=rho/(1 + rho) and q=1/(1 + rho) instead of
! rho and by scaling the differences x(i+1) - x(i) by avh.
!
! the values in df are assumed to have been scaled so that the
! sum of their squared values is n.  the value in var, when it is
! non-negative, is assumed to have been scaled to compensate for
! the scaling of the values in df.
!
! the value returned in fun is an estimate of the true mean square
! when var is non-negative, and is the generalized cross validation
! when var is negative.
!
!---specifications for arguments---
    integer ic,n
! all variables should be double precision:
    real, dimension (:) :: x, dy, a
    real, dimension (6) :: stat
    real, dimension (ic,3) :: c
    real, dimension (0:n+1,3) :: r
    real, dimension (0:n+1,2) :: t
    real, dimension (0:n+1) :: u, v
    real :: fun, var, avh, p, q, rho
!
!---local variables---
    integer i
    real :: e,f,g,h,rho1
    real :: zero = 0.
    real :: one = 1.
    real :: two = 2.

!
!---use p and q instead of rho to prevent overflow or underflow---
    rho1 = one + rho
    p = rho/rho1
    q = one/rho1
    if (abs(rho1-one)<epsilon(0.)) p = zero
    if (abs(rho1-rho)<epsilon(0.)) q = zero
!
!---rational cholesky decomposition of p*c + q*t---
    f = zero
    g = zero
    h = zero
    do i = 0,1
       r(i,1) = zero
    end do
    do i = 2,n - 1
       r(i-2,3) = g*r(i-2,1)
       r(i-1,2) = f*r(i-1,1)
       r(i,1) = one/ (p*c(i,1)+q*t(i,1)-f*r(i-1,2)-g*r(i-2,3))
       f = p*c(i,2) + q*t(i,2) - h*r(i-1,2)
       g = h
       h = p*c(i,3)
    end do
!
!---solve for u---
    u(0) = zero
    u(1) = zero
    do i = 2,n - 1
       u(i) = a(i) - r(i-1,2)*u(i-1) - r(i-2,3)*u(i-2)
    end do
    u(n) = zero
    u(n+1) = zero
    do i = n - 1,2,-1
       u(i) = r(i,1)*u(i) - r(i,2)*u(i+1) - r(i,3)*u(i+2)
    end do
!
!---calculate residual vector v---
    e = zero
    h = zero
    do i = 1,n - 1
       g = h
       h = (u(i+1)-u(i))/ ((x(i+1)-x(i))/avh)
       v(i) = dy(i)* (h-g)
       e = e + v(i)*v(i)
    end do
    v(n) = dy(n)* (-h)
    e = e + v(n)*v(n)
!
!---calculate upper three bands of inverse matrix---
    r(n,1) = zero
    r(n,2) = zero
    r(n+1,1) = zero
    do i = n - 1,2,-1
       g = r(i,2)
       h = r(i,3)
       r(i,2) = -g*r(i+1,1) - h*r(i+1,2)
       r(i,3) = -g*r(i+1,2) - h*r(i+2,1)
       r(i,1) = r(i,1) - g*r(i,2) - h*r(i,3)
    end do
!
!---calculate trace---
    f = zero
    g = zero
    h = zero
    do i = 2,n - 1
       f = f + r(i,1)*c(i,1)
       g = g + r(i,2)*c(i,2)
       h = h + r(i,3)*c(i,3)
    end do
    f = f + two* (g+h)
!
!---calculate statistics---
    stat(1) = p
    stat(2) = f*p
    stat(3) = n*e/ (f*f)
    stat(4) = e*p*p/n
    stat(6) = e*p/f
    if (var.ge.zero) go to 80
    stat(5) = stat(6) - stat(4)
    fun = stat(3)
    go to 90

80  stat(5) = amax1(stat(4)-two*var*stat(2)/n+var,zero)
    fun = stat(5)
90  return
  end subroutine spfit1

  subroutine sperr1(x,avh,dy,n,r,p,var,se)
!
! calculates bayesian estimates of the standard errors of the fitted
! values of a cubic smoothing spline by calculating the diagonal elements
! of the influence matrix.
!
!---specifications for arguments---
    integer n
    real, dimension (:) :: x, dy, se
    real, dimension (0:n+1, 3) :: r
    real :: avh, p, var
!
!---specifications for local variables---
    integer i
    real :: f,g,h,f1,g1,h1
    real :: zero = 0.
    real :: one = 1.
!
!---initialize---
    h = avh/ (x(2)-x(1))
    se(1) = one - p*dy(1)*dy(1)*h*h*r(2,1)
    r(1,1) = zero
    r(1,2) = zero
    r(1,3) = zero
!
!---calculate diagonal elements---
    do i = 2,n - 1
       f = h
       h = avh/ (x(i+1)-x(i))
       g = -f - h
       f1 = f*r(i-1,1) + g*r(i-1,2) + h*r(i-1,3)
       g1 = f*r(i-1,2) + g*r(i,1) + h*r(i,2)
       h1 = f*r(i-1,3) + g*r(i,2) + h*r(i+1,1)
       se(i) = one - p*dy(i)*dy(i)* (f*f1+g*g1+h*h1)
    end do
    se(n) = one - p*dy(n)*dy(n)*h*h*r(n-1,1)
!
!---calculate standard error estimates---
    do i = 1,n
       se(i) = sqrt(amax1(se(i)*var,zero))*dy(i)
    end do
    
  end subroutine sperr1

  subroutine spcof1(x,avh,y,dy,n,p,q,a,c,ic,u,v)
!
! calculates coefficients of a cubic smoothing spline from
! parameters calculated by subroutine spfit1.
!
!---specifications for arguments---
    integer ic,n
    real, dimension (:) :: x, y, dy, a
    real, dimension (ic, 3) :: c
    real, dimension (0:n+1) :: u, v
    real :: p, q, avh
!
!---specifications for local variables---
    integer i
    real :: h,qh
!
!---calculate a---
    qh = q/ (avh*avh)
    do i = 1,n
       a(i) = y(i) - p*dy(i)*v(i)
       u(i) = qh*u(i)
    end do
!
!---calculate c---
    do i = 1,n - 1
       h = x(i+1) - x(i)
       c(i,3) = (u(i+1)-u(i))/ (3.0d0*h)
       c(i,1) = (a(i+1)-a(i))/h - (h*c(i,3)+u(i))*h
       c(i,2) = u(i)
    end do
  end subroutine spcof1

!---c cubgcv test driver
!---c ------------------
!---c
!---c author          - m.f.hutchinson
!---c                   csiro division of water and land resources
!---c                   gpo box 1666
!---c                   canberra act 2601
!---c                   australia
!---c
!---c latest revision - 7 august 1986
!---c
!---c computer        - vax/double
!---c
!---c usage           - main program
!---c
!---c required routines - cubgcv,spint1,spfit1,spcof1,sperr1,ggrand
!---c
!---c remarks   uses subroutine cubgcv to fit a cubic smoothing spline
!---c           to 50 data points which are generated by adding a random
!---c           variable with uniform density in the interval [-0.3,0.3]
!---c           to 50 points sampled from the curve  y=sin(3*pi*x/2).
!---c           random deviates in the interval [0,1] are generated by the
!---c           double precision function ggrand (similar to imsl function
!---c           ggubfs).  the abscissae are unequally spaced in [0,1].
!---c
!---c           point standard error estimates are returned in se by
!---c           setting job=1.  the error variance estimate is returned
!---c           in var.  it compares favourably with the true value of 0.03.
!---c           summary statistics from the array wk are written to
!---c           unit 6.  data values and fitted values with estimated
!---c           standard errors are also written to unit 6.
!---c
!---      parameter (n=50, ic=49)
!---c
!---      integer            job,ier
!---      double precision   x(n),f(n),y(n),df(n),c(ic,3),wk(7*(n+2)),
!---     *                   var,se(n)
!---      double precision   ggrand,dseed
!---c
!---c---initialize---
!---      dseed=1.2345d4
!---      job=1
!---      var=-1.0d0
!---c
!---c---calculate data points---
!---      do 10 i=1,n
!---      x(i)=(i - 0.5)/n + (2.0*ggrand(dseed) - 1.0)/(3.0*n)
!---      f(i)=dsin(4.71238*x(i)) + (2.0*ggrand(dseed) - 1.0)*0.3
!---      df(i)=1.0d0
!---  10  continue
!---c
!---c---fit cubic spline---
!---      call cubgcv(x,f,df,n,y,c,ic,var,job,se,wk,ier)
!---c
!---c---write out results---
!---      write(6,20)
!---  20  format(' cubgcv test driver results:')
!---      write(6,30) ier,var,wk(3),wk(4),wk(2)
!---  30  format(/' ier =',i4/' var =',f7.4/
!---     *        ' generalized cross validation =',f7.4/
!---     *        ' mean square residual         =',f7.4/
!---     *        ' residual degrees of freedom  =',f7.2)
!---      write(6,40)
!---  40  format(/' input data',17x,'output results'//
!---     *         '   i    x(i)    f(i)',6x,'    y(i)   se(i)',
!---     *          '      c(i,1)      c(i,2)      c(i,3)')
!---      do 60 i=1,n-1
!---      write(6,50) i,x(i),f(i),y(i),se(i),(c(i,j),j=1,3)
!---  50  format(i4,2f8.4,6x,2f8.4,3e12.4)
!---  60  continue
!---      write(6,50) n,x(n),f(n),y(n),se(n)
!---      stop
!---      end
!---      double precision function ggrand(dseed)
!---c
!---c double precision uniform random number generator
!---c
!---c constants: a = 7**5
!---c            b = 2**31 - 1
!---c            c = 2**31
!---c
!---c reference: imsl manual, chapter g - generation and testing of
!---c                                     random numbers
!---c
!---c---specifications for arguments---
!---      double precision dseed
!---c
!---c---specifications for local variables---
!---      double precision a,b,c,s
!---c
!---      data a,b,c/16807.0d0, 2147483647.0d0, 2147483648.0d0/
!---c
!---      s=dseed
!---      s=dmod(a*s, b)
!---      ggrand=s/c
!---      dseed=s
!---      return
!---      end
!---
!---cubgcv test driver results:
!---
!---ier =   0
!---var = 0.0279
!---generalized cross validation = 0.0318
!---mean square residual         = 0.0246
!---residual degrees of freedom  =  43.97
!---
!---input data                 output results
!---
!---  i    x(i)    f(i)          y(i)   se(i)      c(i,1)      c(i,2)      c(i,3)
!---  1  0.0046  0.2222        0.0342  0.1004  0.3630e+01  0.0000e+00  0.2542e+02
!---  2  0.0360 -0.1098        0.1488  0.0750  0.3705e+01  0.2391e+01 -0.9537e+01
!---  3  0.0435 -0.0658        0.1767  0.0707  0.3740e+01  0.2175e+01 -0.4233e+02
!---  4  0.0735  0.3906        0.2900  0.0594  0.3756e+01 -0.1642e+01 -0.2872e+02
!---  5  0.0955  0.6054        0.3714  0.0558  0.3642e+01 -0.3535e+01  0.2911e+01
!---  6  0.1078  0.3034        0.4155  0.0549  0.3557e+01 -0.3428e+01 -0.1225e+02
!---  7  0.1269  0.7386        0.4822  0.0544  0.3412e+01 -0.4131e+01  0.2242e+02
!---  8  0.1565  0.4616        0.5800  0.0543  0.3227e+01 -0.2143e+01  0.6415e+01
!---  9  0.1679  0.4315        0.6165  0.0543  0.3180e+01 -0.1923e+01 -0.1860e+02
!--- 10  0.1869  0.5716        0.6762  0.0544  0.3087e+01 -0.2985e+01 -0.3274e+02
!--- 11  0.2149  0.6736        0.7595  0.0542  0.2843e+01 -0.5733e+01 -0.4435e+02
!--- 12  0.2356  0.7388        0.8155  0.0539  0.2549e+01 -0.8486e+01 -0.5472e+02
!--- 13  0.2557  1.1953        0.8630  0.0537  0.2139e+01 -0.1180e+02 -0.9784e+01
!--- 14  0.2674  1.0299        0.8864  0.0536  0.1860e+01 -0.1214e+02  0.9619e+01
!--- 15  0.2902  0.7981        0.9225  0.0534  0.1322e+01 -0.1149e+02 -0.7202e+01
!--- 16  0.3155  0.8973        0.9485  0.0532  0.7269e+00 -0.1203e+02 -0.1412e+02
!--- 17  0.3364  1.2695        0.9583  0.0530  0.2040e+00 -0.1292e+02  0.2796e+02
!--- 18  0.3557  0.7253        0.9577  0.0527 -0.2638e+00 -0.1130e+02 -0.3453e+01
!--- 19  0.3756  1.2127        0.9479  0.0526 -0.7176e+00 -0.1151e+02  0.3235e+02
!--- 20  0.3881  0.7304        0.9373  0.0525 -0.9889e+00 -0.1030e+02  0.4381e+01
!--- 21  0.4126  0.9810        0.9069  0.0525 -0.1486e+01 -0.9977e+01  0.1440e+02
!--- 22  0.4266  0.7117        0.8842  0.0525 -0.1756e+01 -0.9373e+01 -0.8925e+01
!--- 23  0.4566  0.7203        0.8227  0.0524 -0.2344e+01 -0.1018e+02 -0.2278e+02
!--- 24  0.4704  0.9242        0.7884  0.0524 -0.2637e+01 -0.1112e+02 -0.4419e+01
!--- 25  0.4914  0.7345        0.7281  0.0523 -0.3110e+01 -0.1140e+02 -0.3562e+01
!--- 26  0.5084  0.7378        0.6720  0.0524 -0.3500e+01 -0.1158e+02  0.5336e+01
!--- 27  0.5277  0.7441        0.6002  0.0525 -0.3941e+01 -0.1127e+02  0.2479e+02
!--- 28  0.5450  0.5612        0.5286  0.0527 -0.4310e+01 -0.9980e+01  0.2920e+02
!--- 29  0.5641  0.5049        0.4429  0.0529 -0.4659e+01 -0.8309e+01  0.3758e+02
!--- 30  0.5857  0.4725        0.3390  0.0531 -0.4964e+01 -0.5878e+01  0.5563e+02
!--- 31  0.6159  0.1380        0.1850  0.0531 -0.5167e+01 -0.8307e+00  0.4928e+02
!--- 32  0.6317  0.1412        0.1036  0.0531 -0.5157e+01  0.1499e+01  0.5437e+02
!--- 33  0.6446 -0.1110        0.0371  0.0531 -0.5091e+01  0.3614e+01  0.3434e+02
!--- 34  0.6707 -0.2605       -0.0927  0.0532 -0.4832e+01  0.6302e+01  0.1164e+02
!--- 35  0.6853 -0.1284       -0.1619  0.0533 -0.4640e+01  0.6812e+01  0.1617e+02
!--- 36  0.7064 -0.3452       -0.2564  0.0536 -0.4332e+01  0.7834e+01  0.4164e+01
!--- 37  0.7310 -0.5527       -0.3582  0.0538 -0.3939e+01  0.8141e+01 -0.2214e+02
!--- 38  0.7531 -0.3459       -0.4415  0.0540 -0.3611e+01  0.6674e+01 -0.9205e+01
!--- 39  0.7686 -0.5902       -0.4961  0.0541 -0.3410e+01  0.6245e+01 -0.2193e+02
!--- 40  0.7952 -0.7644       -0.5828  0.0541 -0.3125e+01  0.4494e+01 -0.4649e+02
!--- 41  0.8087 -0.5392       -0.6242  0.0541 -0.3029e+01  0.2614e+01 -0.3499e+02
!--- 42  0.8352 -0.4247       -0.7031  0.0539 -0.2964e+01 -0.1603e+00  0.2646e+01
!--- 43  0.8501 -0.6327       -0.7476  0.0538 -0.2967e+01 -0.4132e-01  0.1817e+02
!--- 44  0.8726 -0.9983       -0.8139  0.0538 -0.2942e+01  0.1180e+01 -0.6774e+01
!--- 45  0.8874 -0.9082       -0.8574  0.0542 -0.2911e+01  0.8778e+00 -0.1364e+02
!--- 46  0.9139 -0.8930       -0.9340  0.0566 -0.2893e+01 -0.2044e+00 -0.8094e+01
!--- 47  0.9271 -1.0233       -0.9723  0.0593 -0.2903e+01 -0.5258e+00 -0.1498e+02
!--- 48  0.9473 -0.8839       -1.0313  0.0665 -0.2942e+01 -0.1433e+01  0.4945e+01
!--- 49  0.9652 -1.0172       -1.0843  0.0766 -0.2989e+01 -0.1168e+01  0.1401e+02
!--- 50  0.9930 -1.2715       -1.1679  0.0998
!---documentation:
!---c   computer            - vax/double
!---c
!---c   author              - m.f.hutchinson
!---c                         csiro division of mathematics and statistics
!---c                         p.o. box 1965
!---c                         canberra, act 2601
!---c                         australia
!---c
!---c   latest revision     - 15 august 1985
!---c
!---c   purpose             - cubic spline data smoother
!---c
!---c   usage               - call cubgcv (x,f,df,n,y,c,ic,var,job,se,wk,ier)
!---c
!---c   arguments    x      - vector of length n containing the
!---c                           abscissae of the n data points
!---c                           (x(i),f(i)) i=1..n. (input) x
!---c                           must be ordered so that
!---c                           x(i) .lt. x(i+1).
!---c                f      - vector of length n containing the
!---c                           ordinates (or function values)
!---c                           of the n data points (input).
!---c                df     - vector of length n. (input/output)
!---c                           df(i) is the relative standard deviation
!---c                           of the error associated with data point i.
!---c                           each df(i) must be positive.  the values in
!---c                           df are scaled by the subroutine so that
!---c                           their mean square value is 1, and unscaled
!---c                           again on normal exit.
!---c                           the mean square value of the df(i) is returned
!---c                           in wk(7) on normal exit.
!---c                           if the absolute standard deviations are known,
!---c                           these should be provided in df and the error
!---c                           variance parameter var (see below) should then
!---c                           be set to 1.
!---c                           if the relative standard deviations are unknown,
!---c                           set each df(i)=1.
!---c                n      - number of data points (input).
!---c                           n must be .ge. 3.
!---c                y,c    - spline coefficients. (output) y
!---c                           is a vector of length n. c is
!---c                           an n-1 by 3 matrix. the value
!---c                           of the spline approximation at t is
!---c                           s(t)=((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
!---c                           where x(i).le.t.lt.x(i+1) and
!---c                           d = t-x(i).
!---c                ic     - row dimension of matrix c exactly
!---c                           as specified in the dimension
!---c                           statement in the calling program. (input)
!---c                var    - error variance. (input/output)
!---c                           if var is negative (i.e. unknown) then
!---c                           the smoothing parameter is determined
!---c                           by minimizing the generalized cross validation
!---c                           and an estimate of the error variance is
!---c                           returned in var.
!---c                           if var is non-negative (i.e. known) then the
!---c                           smoothing parameter is determined to minimize
!---c                           an estimate, which depends on var, of the true
!---c                           mean square error, and var is unchanged.
!---c                           in particular, if var is zero, then an
!---c                           interpolating natural cubic spline is calculated.
!---c                           var should be set to 1 if absolute standard
!---c                           deviations have been provided in df (see above).
!---c                job    - job selection parameter. (input)
!---c                         job = 0 should be selected if point standard error
!---c                           estimates are not required in se.
!---c                         job = 1 should be selected if point standard error
!---c                           estimates are required in se.
!---c                se     - vector of length n containing bayesian standard
!---c                           error estimates of the fitted spline values in y.
!---c                           se is not referenced if job=0. (output)
!---c                wk     - work vector of length 7*(n + 2). on normal exit the
!---c                           first 7 values of wk are assigned as follows:-
!---c
!---c                           wk(1) = smoothing parameter (= rho/(rho + 1))
!---c                           wk(2) = estimate of the number of degrees of
!---c                                   freedom of the residual sum of squares
!---c                           wk(3) = generalized cross validation
!---c                           wk(4) = mean square residual
!---c                           wk(5) = estimate of the true mean square error
!---c                                   at the data points
!---c                           wk(6) = estimate of the error variance
!---c                           wk(7) = mean square value of the df(i)
!---c
!---c                           if wk(1)=0 (rho=0) an interpolating natural cubic
!---c                           spline has been calculated.
!---c                           if wk(1)=1 (rho=infinite) a least squares
!---c                           regression line has been calculated.
!---c                           wk(2) is an estimate of the number of degrees of
!---c                           freedom of the residual which reduces to the
!---c                           usual value of n-2 when a least squares regression
!---c                           line is calculated.
!---c                           wk(3),wk(4),wk(5) are calculated with the df(i)
!---c                           scaled to have mean square value 1.  the
!---c                           unscaled values of wk(3),wk(4),wk(5) may be
!---c                           calculated by dividing by wk(7).
!---c                           wk(6) coincides with the output value of var if
!---c                           var is negative on input.  it is calculated with
!---c                           the unscaled values of the df(i) to facilitate
!---c                           comparisons with a priori variance estimates.
!---c
!---c                ier    - error parameter. (output)
!---c                         terminal error
!---c                           ier = 129, ic is less than n-1.
!---c                           ier = 130, n is less than 3.
!---c                           ier = 131, input abscissae are not
!---c                             ordered so that x(i).lt.x(i+1).
!---c                           ier = 132, df(i) is not positive for some i.
!---c                           ier = 133, job is not 0 or 1.
!---c
!---c   precision/hardware  - double
!---c
!---c   required routines   - spint1,spfit1,spcof1,sperr1
!---c
!---c   remarks      the number of arithmetic operations required by the
!---c                subroutine is proportional to n.  the subroutine
!---c                uses an algorithm developed by m.f. hutchinson and
!---c                f.r. de hoog, 'smoothing noisy data with spline
!---c                functions', numer. math. (in press)
!---c
!---c-----------------------------------------------------------------------
!---c
!---
!---
!---algorithm
!---
!---cubgcv calculates a natural cubic spline curve which smoothes a given set
!---of data points, using statistical considerations to determine the amount
!---of smoothing required, as described in reference 2.  if the error variance
!---is known, it should be supplied to the routine in var.  the degree of
!---smoothing is then determined by minimizing an unbiased estimate of the true
!---mean square error.  on the other hand, if the error variance is not known, var
!---should be set to -1.0.  the routine then determines the degree of smoothing
!---by minimizing the generalized cross validation.  this is asymptotically the
!---same as minimizing the true mean square error (see reference 1).  in this
!---case, an estimate of the error variance is returned in var which may be
!---compared with any a priori approximate estimates.  in either case, an
!---estimate of the true mean square error is returned in wk(5).  this estimate,
!---however, depends on the error variance estimate, and should only be accepted
!---if the error variance estimate is reckoned to be correct.
!---
!---if job is set to 1, bayesian estimates of the standard error of each
!---smoothed data value are returned in the array se.  these also depend on
!---the error variance estimate and should only be accepted if the error
!---variance estimate is reckoned to be correct.  see reference 4.
!---
!---the number of arithmetic operations and the amount of storage required by
!---the routine are both proportional to n, so that very large data sets may be
!---analysed.  the data points do not have to be equally spaced or uniformly
!---weighted.  the residual and the spline coefficients are calculated in the
!---manner described in reference 3, while the trace and various statistics,
!---including the generalized cross validation, are calculated in the manner
!---described in reference 2.
!---
!---when var is known, any value of n greater than 2 is acceptable.  it is
!---advisable, however, for n to be greater than about 20 when var is unknown.
!---if the degree of smoothing done by cubgcv when var is unknown is not
!---satisfactory, the user should try specifying the degree of smoothing by
!---setting var to a reasonable value.
!---
!---references:
!---
!---1.  craven, peter and wahba, grace, "smoothing noisy data with spline
!---    functions", numer. math. 31, 377-403 (1979).
!---2.  hutchinson, m.f. and de hoog, f.r., "smoothing noisy data with spline
!---    functions", numer. math. (in press).
!---3.  reinsch, c.h., "smoothing by spline functions", numer. math. 10,
!---    177-183 (1967).
!---4.  wahba, grace, "bayesian 'confidence intervals' for the cross-validated
!---    smoothing spline", j.r.statist. soc. b 45, 133-150 (1983).
!---
!---
!---example
!---
!---a sequence of 50 data points are generated by adding a random variable with
!---uniform density in the interval [-0.3,0.3] to the curve y=sin(3*pi*x/2).
!---the abscissae are unequally spaced in [0,1].  point standard error estimates
!---are returned in se by setting job to 1.  the error variance estimate is
!---returned in var.  it compares favourably with the true value of 0.03.
!---the imsl function ggubfs is used to generate sample values of a uniform
!---variable on [0,1].
!---
!---
!---input:
!---
!---      integer          n,ic,job,ier
!---      double precision x(50),f(50),y(50),df(50),c(49,3),wk(400),
!---     *                 var,se(50)
!---      double precision ggubfs,dseed,dn
!---      data dseed/1.2345d4/
!---c
!---      n=50
!---      ic=49
!---      job=1
!---      var=-1.0d0
!---      dn=n
!---c
!---      do 10 i=1,n
!---      x(i)=(i - 0.5)/dn + (2.0*ggubfs(dseed) - 1.0)/(3.0*dn)
!---      f(i)=dsin(4.71238*x(i)) + (2.0*ggubfs(dseed) - 1.0)*0.3
!---      df(i)=1.0d0
!---  10  continue
!---      call cubgcv(x,f,df,n,y,c,ic,var,job,se,wk,ier)
!---       .
!---       .
!---       .
!---      end
!---
!---output:
!---
!---ier = 0
!---var = 0.0279
!---generalized cross validation = 0.0318
!---mean square residual = 0.0246
!---residual degrees of freedom = 43.97
!---for checking purposes the following output is given:
!---
!---x(1)  = 0.0046    f(1)  =  0.2222     y(1)  =  0.0342     se(1)  = 0.1004
!---x(21) = 0.4126    f(21) =  0.9810     y(21) =  0.9069     se(21) = 0.0525
!---x(41) = 0.8087    f(41) = -0.5392     y(41) = -0.6242     se(41) = 0.0541
!---

end module gridgen4mod
