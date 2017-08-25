module finite_differences

  implicit none

  public :: third_order_upwind
  public :: third_order_upwind_zed
  public :: fd3pt
  public :: d2_3pt

  interface fd3pt
     module procedure fd3pt_real
     module procedure fd3pt_array
  end interface

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

  ! only good for equally-spaced grid-pts
  subroutine fd3pt_real (prof, profgrad, dr)
    
    implicit none
    
    real, dimension (:), intent (in) :: prof
    real, dimension (:), intent (out) :: profgrad
    real, intent (in) :: dr
    
    integer :: ix, npts
    real, dimension (:), allocatable :: aa, bb, cc
    
    npts = size(prof)
    allocate (aa(npts), bb(npts), cc(npts))
    
    aa = 1.0 ; bb = 4.0 ; cc = 1.0
    aa(1) = 0.0 ; bb(1) = 0.5 ; cc(1) = 0.5
    aa(npts) = 0.5 ; bb(npts) = 0.5 ; cc(npts) = 0.0
    
    do ix = 2, npts-1
       profgrad(ix) = 3.0 * (prof(ix+1) - prof(ix-1)) / dr
    end do
    profgrad(1) = (prof(2)-prof(1))/dr
    profgrad(npts) = (prof(npts)-prof(npts-1))/dr
    
    call tridag (aa, bb, cc, profgrad)
    
    deallocate (aa, bb, cc)
    
  end subroutine fd3pt_real
  
  subroutine fd3pt_array (prof, profgrad, dr)
    
    implicit none
    
    real, dimension (:), intent (in) :: prof, dr
    real, dimension (:), intent (out) :: profgrad
    
    integer :: ix, npts
    real, dimension (:), allocatable :: aa, bb, cc
    
    npts = size(prof)
    allocate (aa(npts), bb(npts), cc(npts))
    
    do ix = 2, npts-1
       profgrad(ix) = ((prof(ix)-prof(ix-1))*dr(ix)/dr(ix-1) &
            + (prof(ix+1)-prof(ix))*dr(ix-1)/dr(ix)) / (dr(ix-1)+dr(ix))
    end do
    profgrad(1) = (prof(2)-prof(1))/dr(1)
    profgrad(npts) = (prof(npts)-prof(npts-1))/dr(npts-1)
    
    deallocate (aa, bb, cc)
    
  end subroutine fd3pt_array

  ! second derivative using centered differences
  ! second order accurate
  subroutine d2_3pt (f, d2f, dr)

    implicit none

    real, dimension (:), intent (in) :: f
    real, dimension (:), intent (in) :: dr
    real, dimension (:), intent (out) :: d2f

    real :: a, b, c
    integer :: i, n

    n = size(f)

    do i = 2, n-1
       a = 2./(dr(i-1)*(dr(i)+dr(i-1)))
       b = -2./(dr(i-1)*dr(i))
       c = 2./(dr(i)*(dr(i)+dr(i-1)))
       d2f(i) = a*f(i-1)+b*f(i)+c*f(i+1)
    end do
    ! FLAG -- this is a hack
    ! do not anticipate needing 2nd derivatives
    ! at first and last grid points
    d2f(1) = d2f(2)
    d2f(n) = d2f(n-1)

  end subroutine d2_3pt

  subroutine tridag (aa, bb, cc, sol)
    
    implicit none
    
    real, dimension (:), intent (in) :: aa, bb, cc
    real, dimension (:), intent (in out) :: sol
    
    integer :: ix, npts
    real :: bet
    
    real, dimension (:), allocatable :: gam
    
    npts = size(aa)
    allocate (gam(npts))
    
    bet = bb(1)
    sol(1) = sol(1)/bet
    
    do ix = 2, npts
       gam(ix) = cc(ix-1)/bet
       bet = bb(ix) - aa(ix)*gam(ix)
       if (bet == 0.0) write (*,*) 'tridiagonal solve failed'
       sol(ix) = (sol(ix)-aa(ix)*sol(ix-1))/bet
    end do

    do ix = npts-1, 1, -1
       sol(ix) = sol(ix) - gam(ix+1)*sol(ix+1)
    end do

    deallocate (gam)

  end subroutine tridag
  
end module finite_differences
