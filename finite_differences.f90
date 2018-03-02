module finite_differences

  implicit none

  public :: first_order_upwind
  public :: third_order_upwind
  public :: third_order_upwind_zed
  public :: first_order_upwind_zed
  public :: second_order_centered_zed
  public :: fd3pt, fd5pt
  public :: d2_3pt
  public :: fd_variable_upwinding_vpa
  public :: fd_variable_upwinding_zed
  public :: fd_cell_centres_zed, cell_centres_zed

  interface fd3pt
     module procedure fd3pt_real
     module procedure fd3pt_array
  end interface

  interface fd5pt
     module procedure fd5pt_real
     module procedure fd5pt_array
  end interface

  interface first_order_upwind
     module procedure first_order_upwind_real
     module procedure first_order_upwind_complex
  end interface

  interface third_order_upwind
     module procedure third_order_upwind_complex
     module procedure third_order_upwind_real
  end interface

  interface tridag
     module procedure tridag_real
     module procedure tridag_complex
  end interface

  interface second_order_centered_zed
     module procedure second_order_centered_zed_real
     module procedure second_order_centered_zed_complex
  end interface

contains

  subroutine first_order_upwind_real (llim, f, del, sgn, df)
    
    implicit none
    
    integer, intent (in) :: llim
    real, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    real, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)

    if (sgn == -1) then
       istart = llim
       iend = llim+n-1
    else
       istart = llim+n-1
       iend = llim
    end if

    ! zero BC, 1st order accurate upwind
    df(istart) = -f(istart)*sgn/del
    do i = istart-sgn, iend, -sgn
       df(i) = sgn*(f(i+sgn)-f(i))/del
    end do

  end subroutine first_order_upwind_real
  
  subroutine first_order_upwind_complex (llim, f, del, sgn, df)
    
    implicit none
    
    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    complex, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)

    if (sgn == -1) then
       istart = llim
       iend = llim+n-1
    else
       istart = llim+n-1
       iend = llim
    end if

    ! zero BC, 1st order accurate upwind
    df(istart) = -f(istart)*sgn/del
    do i = istart-sgn, iend, -sgn
       df(i) = sgn*(f(i+sgn)-f(i))/del
    end do

  end subroutine first_order_upwind_complex
  
  subroutine third_order_upwind_complex (llim, f, del, sgn, df)
    
    implicit none
    
    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    complex, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)
    if (sgn == -1) then
       istart = llim
       iend = llim+n-1
    else
       istart = llim+n-1
       iend = llim
    end if

    ! zero BC, 1st order accurate upwind
    df(istart) = -f(istart)*sgn/del
    ! zero BC, 3rd order accurate upwind
    i = istart-sgn
    df(i) = -sgn*(2.*f(i-sgn)+3.*f(i)-6.*f(i+sgn))/(6.*del)
    ! 1st order accurate upwind
    df(iend) = sgn*(f(iend+sgn)-f(iend))/del

    ! 3rd order accurate upwind
    do i = istart-2*sgn, iend+sgn, -sgn
       df(i) = -sgn*(2.*f(i-sgn)+3*f(i)-6.*f(i+sgn)+f(i+2*sgn))/(6.*del)
    end do
!     ! zero BC, 1st order accurate upwind
!     i = -llim*sgn
!     df(i) = -f(i)*sgn/del
!     ! zero BC, 3rd order accurate upwind
!     i = -(llim+1)*sgn
!     df(i) = -sgn*(2.*f(i-sgn)+3.*f(i)-6.*f(i+sgn))/(6.*del)
!     ! 1st order accurate upwind
!     i = llim*sgn
!     df(i) = sgn*(f(i+sgn)-f(i))/del

!     ! 3rd order accurate upwind
!     do i = -(llim+2)*sgn, (llim+1)*sgn, -sgn
!        df(i) = -sgn*(2.*f(i-sgn)+3*f(i)-6.*f(i+sgn)+f(i+2*sgn))/(6.*del)
!     end do

  end subroutine third_order_upwind_complex

  subroutine third_order_upwind_real (llim, f, del, sgn, df)
    
    implicit none
    
    integer, intent (in) :: llim
    real, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    real, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)
    if (sgn == -1) then
       istart = llim
       iend = llim+n-1
    else
       istart = llim+n-1
       iend = llim
    end if

    ! zero BC, 1st order accurate upwind
    df(istart) = -f(istart)*sgn/del
    ! zero BC, 3rd order accurate upwind
    i = istart-sgn
    df(i) = -sgn*(2.*f(i-sgn)+3.*f(i)-6.*f(i+sgn))/(6.*del)
    ! 1st order accurate upwind
    df(iend) = sgn*(f(iend+sgn)-f(iend))/del

    ! 3rd order accurate upwind
    do i = istart-2*sgn, iend+sgn, -sgn
       df(i) = -sgn*(2.*f(i-sgn)+3*f(i)-6.*f(i+sgn)+f(i+2*sgn))/(6.*del)
    end do

!     ! zero BC, 1st order accurate upwind
!     i = -llim*sgn
!     df(i) = -f(i)*sgn/del
!     ! zero BC, 3rd order accurate upwind
!     i = -(llim+1)*sgn
!     df(i) = -sgn*(2.*f(i-sgn)+3.*f(i)-6.*f(i+sgn))/(6.*del)
!     ! 1st order accurate upwind
!     i = llim*sgn
!     df(i) = sgn*(f(i+sgn)-f(i))/del

!     ! 3rd order accurate upwind
!     do i = -(llim+2)*sgn, (llim+1)*sgn, -sgn
!        df(i) = -sgn*(2.*f(i-sgn)+3*f(i)-6.*f(i+sgn)+f(i+2*sgn))/(6.*del)
!     end do

  end subroutine third_order_upwind_real

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

  subroutine first_order_upwind_zed (llim, iseg, nseg, f, del, sgn, fl, fr, df)

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
          df(i) = (f(i+1)-f(i))/del
       else
          i = ulim
          df(i) = (fr(1)-f(i))/del
          i = ulim-1
          df(i) = (f(i+1)-f(i))/del
       end if
       i = llim
       df(i) = (f(i+1)-f(i))/del
       istart = ulim
       iend = llim
    else
       if (iseg == 1) then
          i = llim
          df(i) = f(i)/del
          i = llim+1
          df(i) = (f(i)-f(i-1))/del
       else
          i = llim
          df(i) = (f(i)-fl(2))/del
          i = llim+1
          df(i) = (f(i)-f(i-1))/del
       end if
       i = ulim
       df(i) = (f(i)-f(i-1))/del
       istart = llim
       iend = ulim
    end if

    ! 3rd order accurate upwind
    do i = istart-2*sgn, iend+sgn, -sgn
       df(i) = sgn*(f(i+sgn)-f(i))/del
    end do

  end subroutine first_order_upwind_zed

  subroutine second_order_centered_zed_real (llim, iseg, nseg, f, del, sgn, fl, fr, df)

    implicit none
    
    integer, intent (in) :: llim, iseg, nseg
    real, dimension (llim:), intent (in) :: f
    integer, intent (in) :: sgn
    real, intent (in) :: del
    real, dimension (:), intent (in) :: fl, fr
    real, dimension (llim:), intent (out) :: df

    integer :: i, ulim

    ulim = size(f)+llim-1    

    i = llim
    if (iseg == 1) then
       if (sgn>0) then
          ! sgn > 0 corresponds to negative advection speed
          ! upwind at boundary requires taking information from right
          df(i) = (f(i+1)-f(i))/del
       else
          ! sgn < 0 corresponds to positive advection speed
          ! upwind at boundary requires taking information from left
          ! with f=0 as incoming BC
          df(i) = f(i)/del
       end if
    else
       df(i) = 0.5*(f(i+1)-fl(2))/del
    end if

    i = ulim
    if (iseg == nseg) then
       if (sgn > 0) then
          ! sgn > 0 corresponds to negative advection speed
          ! upwind at boundary requires taking information from right
          ! and BC is f=0
          df(i) = -f(i)/del
       else
          ! sgn < 0 corresponds to positive advection speed
          ! upwind at boundary requires taking information from left
          df(i) = (f(i)-f(i-1))/del
       end if
    else
       df(i) = 0.5*(fr(1)-f(i-1))/del
    end if

    do i = llim+1, ulim-1
       df(i) = 0.5*(f(i+1)-f(i-1))/del
    end do

  end subroutine second_order_centered_zed_real

  subroutine second_order_centered_zed_complex (llim, iseg, nseg, f, del, sgn, fl, fr, df)

    implicit none
    
    integer, intent (in) :: llim, iseg, nseg
    complex, dimension (llim:), intent (in) :: f
    integer, intent (in) :: sgn
    real, intent (in) :: del
    complex, dimension (:), intent (in) :: fl, fr
    complex, dimension (llim:), intent (out) :: df

    integer :: i, ulim

    ulim = size(f)+llim-1    

    i = llim
    if (iseg == 1) then
       if (sgn>0) then
          ! sgn > 0 corresponds to negative advection speed
          ! upwind at boundary requires taking information from right
          df(i) = (f(i+1)-f(i))/del
       else
          ! sgn < 0 corresponds to positive advection speed
          ! upwind at boundary requires taking information from left
          ! with f=0 as incoming BC
          df(i) = f(i)/del
       end if
    else
       df(i) = 0.5*(f(i+1)-fl(2))/del
    end if

    i = ulim
    if (iseg == nseg) then
       if (sgn > 0) then
          ! sgn > 0 corresponds to negative advection speed
          ! upwind at boundary requires taking information from right
          ! and BC is f=0
          df(i) = -f(i)/del
       else
          ! sgn < 0 corresponds to positive advection speed
          ! upwind at boundary requires taking information from left
          df(i) = (f(i)-f(i-1))/del
       end if
    else
       df(i) = 0.5*(fr(1)-f(i-1))/del
    end if

    do i = llim+1, ulim-1
       df(i) = 0.5*(f(i+1)-f(i-1))/del
    end do

  end subroutine second_order_centered_zed_complex

  subroutine second_order_centered_vpa (llim, f, del, df)

    implicit none
    
    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    complex, dimension (llim:), intent (out) :: df

    integer :: i, ulim

    ulim = size(f)+llim-1    

    i = llim
    df(i) = 0.5*f(i+1)/del

    i = ulim
    df(i) = -0.5*f(i-1)/del

    do i = llim+1, ulim-1
       df(i) = 0.5*(f(i+1)-f(i-1))/del
    end do

  end subroutine second_order_centered_vpa

  subroutine fd_cell_centres_zed (llim, f, del, sgn, fl, fr, df)

    implicit none

    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    complex, intent (in) :: fl, fr
    complex, dimension (llim:), intent (out) :: df

    integer :: i, ulim

    ulim = size(f)+llim-1
       
    if (sgn > 0) then
       ! if sgn > 0, then stream speed is negative
       ! so sweep from more positive to more negative zed
       i = ulim
       df(i) = (fr-f(i))/del
       do i = ulim-1, llim, -1
          df(i) = (f(i+1)-f(i))/del
       end do
    else
       ! if sgn < 0, then stream speed is positive
       ! so sweep from more negative to more positive zed
       i = llim
       df(i) = (f(i)-fl)/del
       do i = llim+1, ulim
          df(i) = (f(i)-f(i-1))/del
       end do
    end if

  end subroutine fd_cell_centres_zed

  ! cell_centres_zed takes f at z grid locations
  ! and returns f at cell centres
  ! (with possible offset due to upwinding)
  subroutine cell_centres_zed (llim, f, upwnd, sgn, fl, fr, fc)

    implicit none

    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: upwnd
    integer, intent (in) :: sgn
    complex, intent (in) :: fl, fr
    complex, dimension (llim:), intent (out) :: fc

    integer :: i, ulim

    ulim = size(f)+llim-1

    if (sgn > 0) then
       ! if sgn > 0, then stream speed is negative
       ! so sweep from more positive to more negative zed
       i = ulim
       fc(i) = 0.5*((1.-upwnd)*fr + (1.+upwnd)*f(i))
       do i = ulim-1, llim, -1
          fc(i) = 0.5*((1.-upwnd)*f(i+1) + (1.+upwnd)*f(i))
       end do
    else
       ! if sgn < 0, then stream speed is positive
       ! so sweep from more negative to more positive zed
       i = llim
       fc(i) = 0.5*((1.+upwnd)*f(i)+(1.-upwnd)*fl)
       do i = llim+1, ulim
          fc(i) = 0.5*((1.+upwnd)*f(i)+(1.-upwnd)*f(i-1))
       end do
    end if

  end subroutine cell_centres_zed

  subroutine fd_variable_upwinding_zed (llim, iseg, nseg, f, del, sgn, upwnd, fl, fr, df)

    implicit none

    integer, intent (in) :: llim, iseg, nseg
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del, upwnd
    integer, intent (in) :: sgn
    complex, dimension (:), intent (in) :: fl, fr
    complex, dimension (llim:), intent (out) :: df

    integer :: i, istart, iend, ulim

    ! if upwnd is zero or if vpa=0, then use centered differences
    if (abs(upwnd) < epsilon(0.) .or. sgn == 0) then
       call second_order_centered_zed (llim, iseg, nseg, f, del, sgn, fl, fr, df)
    else
       ulim = size(f)+llim-1
       
       ! if sgn > 0, then stream speed is negative
       ! so sweep from more positive to more negative zed
       if (sgn > 0) then
          if (iseg == nseg) then
             i = ulim
             df(i) = (0.5*(upwnd-1.)*f(i-1)-upwnd*f(i))/del
          else
             i = ulim
             df(i) = (0.5*(upwnd-1.)*f(i-1)-upwnd*f(i)+0.5*(1.+upwnd)*fr(1))/del
          end if
          if (iseg == 1) then
             i = llim
             ! at left boundary, must upwind fully as no info for f(i-1)
             df(i) = (f(i+1)-f(i))/del
          else
             i = llim
             df(i) = (0.5*(1.+upwnd)*f(i+1)-upwnd*f(i)+0.5*(upwnd-1.)*fl(2))/del
          end if
          istart = ulim
          iend = llim
       else
          if (iseg == 1) then
             i = llim
             df(i) = (0.5*(1.-upwnd)*f(i+1)+upwnd*f(i))/del
          else
             i = llim
             df(i) = (0.5*(1.-upwnd)*f(i+1)+upwnd*f(i)-0.5*(1.+upwnd)*fl(2))/del
          end if
          if (iseg == nseg) then
             i = ulim
             ! if at rightmost zed, have no info for f(i+1) so must fully upwind
             df(i) = (f(i)-f(i-1))/del
          else
             i = ulim
             df(i) = (0.5*(1.-upwnd)*fr(1)+upwnd*f(i)-0.5*(1.+upwnd)*f(i-1))/del
          end if
          istart = llim
          iend = ulim
       end if

       ! mixed 2nd order centered and 1st order upwind scheme
       do i = istart-sgn, iend+sgn, -sgn
          df(i) = sgn*(0.5*(1.+upwnd)*f(i+sgn) - upwnd*f(i) + 0.5*(upwnd-1.)*f(i-sgn))/del
       end do

    end if

  end subroutine fd_variable_upwinding_zed

  subroutine fd_variable_upwinding_vpa (llim, f, del, sgn, upwnd, df)

    implicit none

    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del, upwnd
    integer, intent (in) :: sgn
    complex, dimension (llim:), intent (out) :: df

    integer :: i, istart, iend, ulim

    ! if upwnd is zero or if vpa=0, then use centered differences
    if (abs(upwnd) < epsilon(0.) .or. sgn == 0) then
       call second_order_centered_vpa (llim, f, del, df)
    else
       ulim = size(f)+llim-1
       
       ! if sgn > 0, then stream speed is negative
       ! so sweep from more positive to more negative zed
       if (sgn > 0) then
          i = ulim
          df(i) = (0.5*(upwnd-1.)*f(i-1)-upwnd*f(i))/del
          i = llim
          df(i) = (0.5*(1.+upwnd)*f(i+1)-upwnd*f(i))/del
          istart = ulim
          iend = llim
       else
          i = llim
          df(i) = (0.5*(1.-upwnd)*f(i+1)+upwnd*f(i))/del
          i = ulim
          df(i) = (upwnd*f(i)-0.5*(1.+upwnd)*f(i-1))/del
          istart = llim
          iend = ulim
       end if

       ! mixed centered and 1st order upwind scheme
       do i = istart-sgn, iend+sgn, -sgn
          df(i) = sgn*(0.5*(1.+upwnd)*f(i+sgn) - upwnd*f(i) + 0.5*(upwnd-1.)*f(i-sgn))/del
       end do

    end if

  end subroutine fd_variable_upwinding_vpa

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

  ! boundary points are 2nd-order accurate (2-pt compact difference)
  ! next to boundary points are 4th-order accurate (2-pt centered compact difference)
  ! interior points are 6th-order accurate (4-pt centered compact difference)
  subroutine fd5pt_real (prof, profgrad, dr)

    implicit none

    real, dimension (:), intent (in) :: prof
    real, dimension (:), intent (out) :: profgrad
    real, intent (in) :: dr

    integer :: ix, npts
    real, dimension (:), allocatable :: aa, bb, cc

    npts = size(prof)
    allocate (aa(npts), bb(npts), cc(npts))

    aa = 1.0 ; bb = 3.0 ; cc = 1.0
    aa(1) = 0.0 ; bb(1) = 0.5 ; cc(1) = 0.5
    aa(2) = 1.0 ; bb(2) = 4.0 ; cc(2) = 1.0
    aa(npts-1) = 1.0 ; bb(npts-1) = 4.0 ; cc(npts-1) = 1.0
    aa(npts) = 0.5 ; bb(npts) = 0.5 ; cc(npts) = 0.0

    do ix = 3, npts-2
       profgrad(ix) = (7.*(prof(ix+1) - prof(ix-1)) + 0.25*(prof(ix+2)-prof(ix-2))) / (3.*dr)
    end do
    profgrad(1) = (prof(2)-prof(1))/dr
    profgrad(2) = 3.0*(prof(3) - prof(1))/dr
    profgrad(npts-1) = 3.0*(prof(npts) - prof(npts-2))/dr
    profgrad(npts) = (prof(npts)-prof(npts-1))/dr

    call tridag (aa, bb, cc, profgrad)

    deallocate (aa, bb, cc)

  end subroutine fd5pt_real

  ! boundary points are 2nd-order accurate (2-pt compact difference)
  ! next to boundary points are 4th-order accurate (2-pt centered compact difference)
  ! interior points are 6th-order accurate (4-pt centered compact difference)
  subroutine fd5pt_array (prof, profgrad, dr)

    implicit none

    real, dimension (:), intent (in) :: prof, dr
    real, dimension (:), intent (out) :: profgrad

    integer :: ix, npts
    real, dimension (:), allocatable :: aa, bb, cc

    npts = size(prof)
    allocate (aa(npts), bb(npts), cc(npts))

    aa = 1.0 ; bb = 3.0 ; cc = 1.0
    aa(1) = 0.0 ; bb(1) = 0.5 ; cc(1) = 0.5
    aa(2) = 1.0 ; bb(2) = 4.0 ; cc(2) = 1.0
    aa(npts-1) = 1.0 ; bb(npts-1) = 4.0 ; cc(npts-1) = 1.0
    aa(npts) = 0.5 ; bb(npts) = 0.5 ; cc(npts) = 0.0

    do ix = 3, npts-2
       profgrad(ix) = (7.*(prof(ix+1) - prof(ix-1)) + 0.25*(prof(ix+2)-prof(ix-2))) / (3.*dr(ix))
    end do
    profgrad(1) = (prof(2)-prof(1))/dr(1)
    profgrad(2) = 3.0*(prof(3) - prof(1))/dr(2)
    profgrad(npts-1) = 3.0*(prof(npts) - prof(npts-2))/dr(npts-1)
    profgrad(npts) = (prof(npts)-prof(npts-1))/dr(npts)

    call tridag (aa, bb, cc, profgrad)

    deallocate (aa, bb, cc)

  end subroutine fd5pt_array


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

  subroutine tridag_real (aa, bb, cc, sol)
    
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

  end subroutine tridag_real

  subroutine tridag_complex (llim, aa, bb, cc, sol)
    
    implicit none
    
    integer, intent (in) :: llim
    real, dimension (llim:), intent (in) :: aa, bb, cc
    complex, dimension (llim:), intent (in out) :: sol
    
    integer :: ix, npts
    real :: bet
    
    real, dimension (:), allocatable :: gam
    
    npts = size(bb)
    allocate (gam(llim:llim+npts-1))
    
    bet = bb(llim)
    sol(llim) = sol(llim)/bet
    
    do ix = llim+1, llim+npts-1
       gam(ix) = cc(ix-1)/bet
       bet = bb(ix) - aa(ix)*gam(ix)
       if (bet == 0.0) write (*,*) 'tridiagonal solve failed'
       sol(ix) = (sol(ix)-aa(ix)*sol(ix-1))/bet
    end do

    do ix = llim+npts-2, llim, -1
       sol(ix) = sol(ix) - gam(ix+1)*sol(ix+1)
    end do

    deallocate (gam)

  end subroutine tridag_complex
  
end module finite_differences
