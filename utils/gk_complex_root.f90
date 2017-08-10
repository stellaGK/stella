!*****************************************************************************!
!*****************************************************************************!
!**                        Complex Root Finding Routine                     **!
!*****************************************************************************!
!*****************************************************************************!
!
!------------------------------------------------------------------------------
!                             Copyright,  2005
!                                Greg Howes
!------------------------------------------------------------------------------
!  
!------------------------------------------------------------------------------
module gk_complex_root
  implicit none
  private

  public :: rtsec

 contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
!     NOTE: This routine was adapted from f77 routine by Eliot Quataert
   complex function rtsec(func,x1,x2,xacc,iflag)
     integer, parameter :: maxit=1000
     complex :: func, x1, xl, x2
     complex :: fl, f, swap, dx 
     real    :: xacc
     integer, intent(out), optional :: iflag
     integer :: j
     complex :: maxr,minr
     logical :: limits=.false.

     !Set limits for solution
     if (limits) then
        maxr=x2*10.
        minr=x1*0.1
     endif

     fl=func(x1)
     f=func(x2)
     if(abs(fl).lt.abs(f))then
        rtsec=x1
        xl=x2
        swap=fl
        fl=f
        f=swap
     else
        xl=x1
        rtsec=x2
     endif
     do  j=1,maxit
        iflag = j
!        write(*,'(a,i4)')'Iteration ',j
        if (abs(f-fl) .gt. 1.0E-37) then
           dx=(xl-rtsec)*f/(f-fl)
	else
           dx = (x2-x1)/25.0
	end if        
        xl=rtsec
        fl=f
!	if (Real(rtsec + dx) .gt. 0.075 .or. Real(rtsec+dx) .lt. 0) then
!		rtsec = (x1 + x2)/2.0
!	else
        !NOTE: Reduce jump by 0.5 to improve convergence GH
!        rtsec=rtsec+dx
        rtsec=rtsec+dx/2.
!	end if

        !Enforce limits
        if (limits) then
           !Real limit
           if (real(rtsec) .gt. real(maxr) .or. real(rtsec) .lt. real(minr)) &
              rtsec = cmplx(sqrt(real(minr)*real(maxr)), Aimag(rtsec))
           !NOTE: aimag(rtsec) should be negative!
           if (aimag(rtsec) .lt. aimag(maxr) .or. aimag(rtsec) .gt. aimag(minr)) &
              rtsec = cmplx(real(rtsec), -sqrt(aimag(minr)*aimag(maxr)))
        endif

        !LIMIT: Im(rtsec) < 0.
!	if (Aimag(rtsec) .gt. 0) then
!           rtsec = cmplx(Real(rtsec), -Aimag(rtsec))
!	end if
        !LIMIT: Re(rtsec) > 0.
!	if (Real(rtsec) .lt. 0) then
!           rtsec = cmplx(-Real(rtsec), Aimag(rtsec))
!	end if

        f=func(rtsec)
        if(abs(dx).lt.xacc.or. Abs(f) .eq. 0.) then
!        write(*,'(a,i4)')'Subroutine rtsec Converged on Iteration ',j
        return
        endif
     enddo
 !    stop
!        write(*,'(a,i4)')'Subroutine rtsec stopped at iteration ',j
     return
!     pause 'rtsec exceed maximum iterations'
   end function rtsec
!------------------------------------------------------------------------------
 end module gk_complex_root
  
