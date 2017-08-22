module finite_differences

  implicit none

  public :: third_order_upwind
  public :: third_order_upwind_zed

contains
  
  subroutine third_order_upwind (llim, f, del, sgn, df)
    
    implicit none
    
    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    complex, dimension (llim:), intent (out) :: df
    
    integer :: i

    ! zero BC, 1st order accurate upwind
    i = -llim*sgn
    df(i) = -f(i)*sgn/del
    ! zero BC, 3rd order accurate upwind
    i = -(llim+1)*sgn
    df(i) = -sgn*(2.*f(i-sgn)+3.*f(i)-6.*f(i+sgn))/(6.*del)
    ! 1st order accurate upwind
    i = llim*sgn
    df(i) = sgn*(f(i+sgn)-f(i))/del

    ! 3rd order accurate upwind
    do i = -(llim+2)*sgn, (llim+1)*sgn, -sgn
       df(i) = -sgn*(2.*f(i-sgn)+3*f(i)-6.*f(i+sgn)+f(i+2*sgn))/(6.*del)
    end do

  end subroutine third_order_upwind

  subroutine third_order_upwind_zed (llim, iseg, nseg, f, del, sgn, fl, fr, df)

    implicit none
    
    integer, intent (in) :: llim, iseg, nseg
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    complex, dimension (:), intent (in) :: fl, fr
    complex, dimension (llim:), intent (out) :: df

    integer :: i, istart, iend, ulim

    ulim = size(f)+llim-1    
    ! if sgn > 0, then stream speed is negative
    ! so sweep from more positive to more negative zed
    if (sgn > 0) then
       if (iseg == nseg) then
          i = ulim
          df(i) = -f(i)/del
          i = ulim-1
          df(i) = -(2.*f(i-1)+3.*f(i)-6.*f(i+1))/(6.*del)
       else
          i = ulim
          df(i) = -(2.*f(i-1)+3.*f(i)-6.*fr(1)+fr(2))/(6.*del)
          i = ulim-1
          df(i) = -(2.*f(i-1)+3.*f(i)-6.*f(i+1)+fr(1))/(6.*del)
       end if
       if (iseg == 1) then
          i = llim
          df(i) = (f(i+1)-f(i))/del
       else
          i = llim
          df(i) = -(2.*fl(2)+3*f(i)-6.*f(i+1)+f(i+2))/(6.*del)
       end if
       istart = ulim
       iend = llim
    else
       if (iseg == 1) then
          i = llim
          df(i) = f(i)/del
          i = llim+1
          df(i) = (2.*f(i+1)+3.*f(i)-6.*f(i-1))/(6.*del)
       else
          i = llim
          df(i) = (2.*f(i+1)+3*f(i)-6.*fl(2)+fl(1))/(6.*del)
          i = llim+1
          df(i) = (2.*f(i+1)+3*f(i)-6.*f(i-1)+fl(2))/(6.*del)
       end if
       if (iseg == nseg) then
          i = ulim
          df(i) = (f(i)-f(i-1))/del
       else
          i = ulim
          df(i) = (2.*fr(1)+3*f(i)-6.*f(i-1)+f(i-2))/(6.*del)
       end if
       istart = llim
       iend = ulim
    end if

    ! 3rd order accurate upwind
    do i = istart-2*sgn, iend+sgn, -sgn
       df(i) = -sgn*(2.*f(i-sgn)+3*f(i)-6.*f(i+sgn)+f(i+2*sgn))/(6.*del)
    end do

  end subroutine third_order_upwind_zed

end module finite_differences
