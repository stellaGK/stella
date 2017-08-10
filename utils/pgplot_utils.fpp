!
module pgplot_utils
  use constants, only: kind_rs, kind_rd
  implicit none
  private
  public :: palett, papergrid
  public :: pgkind

# ifdef PGDBLE
  integer, parameter :: pgkind=kind_rd
# else
  integer, parameter :: pgkind=kind_rs
# endif

contains
  
  subroutine palett(type, contra, bright)
    !-----------------------------------------------------------------------
    ! set a "palette" of colors in the range of color indices used by
    ! pgimag.
    !-----------------------------------------------------------------------
    implicit none

    integer :: i
    integer, intent(in):: type
    real (kind=pgkind) :: contra, bright

    real (kind=pgkind) :: &
         & gl(2)=(/0.0_pgkind,1.0_pgkind/), &
         & gr(2)=(/0.0_pgkind,1.0_pgkind/), &
         & gg(2)=(/0.0_pgkind,1.0_pgkind/), &
         & gb(2)=(/0.0_pgkind,1.0_pgkind/)

    real (kind=pgkind) :: &
         & rl(9)=(/ &
         & -0.5_pgkind, 0.0_pgkind,  0.17_pgkind, 0.33_pgkind, &
         & 0.50_pgkind, 0.67_pgkind, 0.83_pgkind, 1.0_pgkind,  &
         & 1.7_pgkind/), &
         & rr(9)=(/ &
         & 0.0_pgkind,  0.0_pgkind,  0.0_pgkind,  0.0_pgkind,  &
         & 0.6_pgkind,  1.0_pgkind,  1.0_pgkind,  1.0_pgkind,  &
         & 1.0_pgkind/), &
         & rg(9)=(/ &
         & 0.0_pgkind,  0.0_pgkind,  0.0_pgkind,  1.0_pgkind,  &
         & 1.0_pgkind,  1.0_pgkind,  0.6_pgkind,  0.0_pgkind,  &
         & 1.0_pgkind/), &
         & rb(9)=(/ &
         & 0.0_pgkind,  0.3_pgkind,  0.8_pgkind,  1.0_pgkind,  &
         & 0.3_pgkind,  0.0_pgkind,  0.0_pgkind,  0.0_pgkind,  &
         & 1.0_pgkind/)

    real (kind=pgkind) :: &
         & hl(5)=(/ &
         & 0.0_pgkind, 0.2_pgkind, 0.4_pgkind, 0.6_pgkind, 1.0_pgkind/), &
         & hr(5)=(/ &
         & 0.0_pgkind, 0.5_pgkind, 1.0_pgkind, 1.0_pgkind, 1.0_pgkind/), &
         & hg(5)=(/ &
         & 0.0_pgkind, 0.0_pgkind, 0.5_pgkind, 1.0_pgkind, 1.0_pgkind/), &
         & hb(5)=(/ &
         & 0.0_pgkind, 0.0_pgkind, 0.0_pgkind, 0.3_pgkind, 1.0_pgkind/)

    real (kind=pgkind) :: &
         & wl(10)=(/ &
         & 0.0_pgkind,  0.5_pgkind,  0.5_pgkind,  0.7_pgkind,  &
         & 0.7_pgkind,  0.85_pgkind, 0.85_pgkind, 0.95_pgkind, &
         & 0.95_pgkind, 1.0_pgkind/), &
         & wr(10)=(/ &
         & 0.0_pgkind,  1.0_pgkind,  0.0_pgkind,  0.0_pgkind,  &
         & 0.3_pgkind,  0.8_pgkind,  0.3_pgkind,  1.0_pgkind,  &
         & 1.0_pgkind,  1.0_pgkind/), &
         & wg(10)=(/ &
         & 0.0_pgkind,  0.5_pgkind,  0.4_pgkind,  1.0_pgkind,  &
         & 0.0_pgkind,  0.0_pgkind,  0.2_pgkind,  0.7_pgkind,  &
         & 1.0_pgkind,  1.0_pgkind/), &
         & wb(10)=(/ &
         & 0.0_pgkind,  0.0_pgkind,  0.0_pgkind,  0.0_pgkind,  &
         & 0.4_pgkind,  1.0_pgkind,  0.0_pgkind,  0.0_pgkind,  &
         & 0.95_pgkind, 1.0_pgkind/)

    real (kind=pgkind) :: &
         & al(20)=(/0.0_pgkind, (0.1_pgkind*i,0.1_pgkind*i,i=1,9), 1.0_pgkind/), &
         & ar(20)=(/(0.0_pgkind,i=1,2), (0.3_pgkind,i=1,2), (0.5_pgkind,i=1,2), &
         & (0.0_pgkind,i=1,8), (1.0_pgkind,i=1,6)/), &
         & ag(20)=(/(0.0_pgkind,i=1,2), (0.3_pgkind,i=1,2), (0.0_pgkind,i=1,4), &
         & (0.8_pgkind,i=1,2), (0.6_pgkind,i=1,2), (1.0_pgkind,i=1,4), &
         & (0.8_pgkind,i=1,2), (0.0_pgkind,i=1,2)/), &
         & ab(20)=(/(0.0_pgkind,i=1,2), (0.3_pgkind,i=1,2), (0.7_pgkind,i=1,4), &
         & (0.9_pgkind,i=1,2), (0.0_pgkind,i=1,10)/)

    integer, parameter :: nn=10
    real (kind=pgkind) :: &
         & wrl(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /), &
         & wrr(nn)=(/ (1._pgkind,i=1,nn) /), &
         & wrg(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /), &
         & wrb(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /)
    real (kind=pgkind) :: &
         & rbl(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /), &
         & rbr(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /), &
         & rbg(nn)=(/ (0._pgkind,i=1,nn) /), &
         & rbb(nn)=(/ (0._pgkind,i=1,nn) /)

    real (kind=pgkind) :: &
         & wbl(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /), &
         & wbr(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /), &
         & wbg(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /), &
         & wbb(nn)=(/ (1._pgkind,i=1,nn) /)
    real (kind=pgkind) :: &
         & bbl(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /), &
         & bbr(nn)=(/ (0._pgkind,i=1,nn) /), &
         & bbg(nn)=(/ (0._pgkind,i=1,nn) /), &
         & bbb(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /)

    real (kind=pgkind) :: &
         & wgl(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /), &
         & wgr(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /), &
         & wgg(nn)=(/ (1._pgkind,i=1,nn) /), &
         & wgb(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /)
    real (kind=pgkind) :: &
         & gbl(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /), &
         & gbr(nn)=(/ (0._pgkind,i=1,nn) /), &
         & gbg(nn)=(/ ((i-1)*1._pgkind/(nn-1),i=1,nn) /), &
         & gbb(nn)=(/ (0._pgkind,i=1,nn) /)

    real (kind=pgkind) :: &
         & brl(2*nn),brr(2*nn),brg(2*nn),brb(2*nn)

    brl(1:2*nn)=(/ ((i-1)*1._pgkind/(2*nn-1),i=1,2*nn) /)
    brr(1:nn)=wbr(1:nn); brr(nn+1:2*nn)=wrr(nn:1:-1)
    brg(1:nn)=wbg(1:nn); brg(nn+1:2*nn)=wrg(nn:1:-1)
    brb(1:nn)=wbb(1:nn); brb(nn+1:2*nn)=wrb(nn:1:-1)

    if (type.eq.1) then
       ! -- gray scale
       call pgctab(gl, gr, gg, gb, 2, contra, bright)
    else if (type.eq.2) then
       ! -- rainbow
       call pgctab(rl, rr, rg, rb, 9, contra, bright)
    else if (type.eq.3) then
       ! -- heat
       call pgctab(hl, hr, hg, hb, 5, contra, bright)
    else if (type.eq.4) then
       ! -- weird iraf
       call pgctab(wl, wr, wg, wb, 10, contra, bright)
    else if (type.eq.5) then
       ! -- aips
       call pgctab(al, ar, ag, ab, 20, contra, bright)
    else if (type.eq.6) then
       ! red -> white
       call pgctab(wrl, wrr, wrg, wrb, nn, contra, bright)
    else if (type.eq.7) then
       ! black -> red
       call pgctab(rbl, rbr, rbg, rbb, nn, contra, bright)
    else if (type.eq.8) then
       ! green -> white
       call pgctab(wgl, wgr, wgg, wgb, nn, contra, bright)
    else if (type.eq.9) then
       ! black -> green
       call pgctab(gbl, gbr, gbg, gbb, nn, contra, bright)
    else if (type.eq.10) then
       ! blue -> white
       call pgctab(wbl, wbr, wbg, wbb, nn, contra, bright)
    else if (type.eq.11) then
       ! black -> blue
       call pgctab(bbl, bbr, bbg, bbb, nn, contra, bright)
    else if (type.eq.12) then
       ! blue -> white -> red
       call pgctab(brl, brr, brg, brb, 2*nn, contra, bright)
    end if

  end subroutine palett

  subroutine pglineb(n,xpts,ypts,blnk)

    implicit none

    integer :: i
    integer, intent(in) :: n
    real (kind=pgkind), intent(in) :: xpts(*), ypts(*)
    integer :: n2
    real (kind=pgkind) :: xpts2(n), ypts2(n)
    real (kind=pgkind), intent(in) :: blnk
    real, parameter :: epsilon=1.e-20


    n2=0
    do i = 1, n
       if(abs(ypts(i)-blnk).lt.epsilon) cycle
       n2=n2+1
       xpts2(n2)=xpts(i)
       ypts2(n2)=ypts(i)
    end do

    call pgline(n2,xpts2,ypts2)

    return

  end subroutine pglineb

  subroutine papergrid(color)
  
    implicit none
  
    logical, intent(in) :: color
    real (kind=pgkind) :: x1,x2,y1,y2
    integer :: ls,lc

    call pgqls(ls)
    call pgqci(lc)
  
    call pgsls(4)
    if(color.eqv..true.) call pgsci(11)
    call pgsvp (0._pgkind,1._pgkind,0._pgkind,1._pgkind)
    call pgswin(0._pgkind,1._pgkind,0._pgkind,1._pgkind)
    call pgbox('G',.1_pgkind,0,'G',.1_pgkind,0)
    call pgsls(ls)
    if(color.eqv..true.) call pgsci(lc)

    return

  end subroutine papergrid

end module pgplot_utils
