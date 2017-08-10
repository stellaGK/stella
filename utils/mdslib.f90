! RN 2008/05/22: modularized
! this is only used in geo/eeq.f90
! I don't understand what this is doing
module mdslib
  implicit none
  private
  public :: mdsopen, mdsconnect, mdsclose

contains

  function mdsopen (c,int )
    integer :: mdsopen
    character*1 :: c
    integer, optional :: int

    mdsopen = 1

  end function mdsopen

  function mdsconnect (c, int)
    integer :: mdsconnect
    character*1 :: c
    integer, optional :: int

    mdsconnect = 1

  end function mdsconnect

!  function mdsclose (c, int)
  subroutine mdsclose (c, int)
!    integer :: mdsclose
    character*1 :: c
    integer :: int

!    mdsclose = 1

!  end function mdsclose
  end subroutine mdsclose

end module mdslib
