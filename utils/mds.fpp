# include "define.inc"

! RN 2008/05/30: Note for precision
! The real and complex variables which is explicitly declared to have
! a given kind should not be promoted to higher precision.
! Real (kind=kind_rs) should be always single precision, and so on.
! I guess MDSplus library does not handle quad precision real, so
! the promotion of precision may cause error. Currently, there's no
! problem because mds related routines are not used for double variables.
! (Note that some compilers (xl, for example) have an option to promote 
! precisions even when the kind is gievn.)
!-----------------------------------------------------------------------
!     file mdsplus_io.f.
!     performs mdsplus tree io
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!     0. mdsplus_io_mod.
!     1. getmdslogical.
!     2. getmdsinteger.
!     3. getmdsdouble.
!     4. getmdsdouble2darray.
!     5. getmdsdouble3darray.
!     6. getmdstext.
!     7. checkmds.
!     8. getmdserrortext.
!     9. getmdsreal
!-----------------------------------------------------------------------
!     subprogram 0. mdsplus_io_mod.
!     module declarations.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.
!-----------------------------------------------------------------------
module mdsio
  use constants, only: kind_rd, kind_rs
  implicit none
  private
# ifdef MDSPLUS
  external mdssetdefault
  logical  mdssetdefault
  external mdsvalue, mdsput
  logical  mdsvalue, mdsput

! these are stolen from  mdslib.inc
  integer, parameter :: IDTYPE_UNSIGNED_BYTE=2, IDTYPE_BYTE=6
  integer, parameter :: IDTYPE_UNSIGNED_SHORT=3, IDTYPE_SHORT=7
  integer, parameter :: IDTYPE_UNSIGNED_LONG=4, IDTYPE_LONG=8
  integer, parameter :: IDTYPE_UNSIGNED_LONGLONG=5, IDTYPE_LONGLONG=8
  integer, parameter :: IDTYPE_FLOAT=10, IDTYPE_DOUBLE=11
  integer, parameter :: IDTYPE_COMPLEX=12, IDTYPE_DOUBLE_COMPLEX=13
  integer, parameter :: IDTYPE_CSTRING=14
# endif
  public :: mds_read
  public :: mds_write
  public :: checkmds
  public :: getmdserrortext
  
  interface mds_read
     module procedure mds_r_log_0
     module procedure mds_r_char_1
     module procedure mds_r_int_0
     module procedure mds_r_int_1
     module procedure mds_r_int_2
     module procedure mds_r_int_3
     module procedure mds_r_int_4
     module procedure mds_r_int_5
     module procedure mds_r_int_6
     module procedure mds_r_int_7
     module procedure mds_r_real_0
     module procedure mds_r_real_1
     module procedure mds_r_real_2
     module procedure mds_r_real_3
     module procedure mds_r_real_4
     module procedure mds_r_real_5
     module procedure mds_r_real_6
     module procedure mds_r_real_7
     module procedure mds_r_doub_0
     module procedure mds_r_doub_1
     module procedure mds_r_doub_2
     module procedure mds_r_doub_3
     module procedure mds_r_doub_4
     module procedure mds_r_doub_5
     module procedure mds_r_doub_6
     module procedure mds_r_doub_7
     module procedure mds_r_cmplx_0
     module procedure mds_r_cmplx_1
     module procedure mds_r_cmplx_2
     module procedure mds_r_cmplx_3
     module procedure mds_r_cmplx_4
     module procedure mds_r_cmplx_5
     module procedure mds_r_cmplx_6
     module procedure mds_r_cmplx_7
     module procedure mds_r_dcmplx_0
     module procedure mds_r_dcmplx_1
     module procedure mds_r_dcmplx_2
     module procedure mds_r_dcmplx_3
     module procedure mds_r_dcmplx_4
     module procedure mds_r_dcmplx_5
     module procedure mds_r_dcmplx_6
     module procedure mds_r_dcmplx_7
  end interface

  interface mds_write
     module procedure mds_w_log_0
     module procedure mds_w_char_1
     module procedure mds_w_int_0
     module procedure mds_w_int_1
     module procedure mds_w_int_2
     module procedure mds_w_int_3
     module procedure mds_w_int_4
     module procedure mds_w_int_5
     module procedure mds_w_int_6
     module procedure mds_w_int_7
     module procedure mds_w_real_0
     module procedure mds_w_real_1
     module procedure mds_w_real_2
     module procedure mds_w_real_3
     module procedure mds_w_real_4
     module procedure mds_w_real_5
     module procedure mds_w_real_6
     module procedure mds_w_real_7
     module procedure mds_w_doub_0
     module procedure mds_w_doub_1
     module procedure mds_w_doub_2
     module procedure mds_w_doub_3
     module procedure mds_w_doub_4
     module procedure mds_w_doub_5
     module procedure mds_w_doub_6
     module procedure mds_w_doub_7
     module procedure mds_w_cmplx_0
     module procedure mds_w_cmplx_1
     module procedure mds_w_cmplx_2
     module procedure mds_w_cmplx_3
     module procedure mds_w_cmplx_4
     module procedure mds_w_cmplx_5
     module procedure mds_w_cmplx_6
     module procedure mds_w_cmplx_7
     module procedure mds_w_dcmplx_0
     module procedure mds_w_dcmplx_1
     module procedure mds_w_dcmplx_2
     module procedure mds_w_dcmplx_3
     module procedure mds_w_dcmplx_4
     module procedure mds_w_dcmplx_5
     module procedure mds_w_dcmplx_6
     module procedure mds_w_dcmplx_7
  end interface

contains
!---------------------------------------------
! mdsplus internal routines
!---------------------------------------------

  subroutine checkmds(status,msg)
    logical, intent(in) :: status
    character(*), intent(in) :: msg
# ifdef MDSPLUS
    character(512) text
    integer istat
    logical lstat
    equivalence(istat,lstat)      
    lstat = status
    if (iand(istat,1) .eq. 0) then
       call getmdserrortext(status,text)
       write(*,*) text
       write(*,*) msg
       !       call program_stop(msg//", "//trim(text))
    endif
# endif
  end subroutine checkmds

  subroutine getmdserrortext(status,text)
    logical, intent(in) :: status
    character(*) :: text ! intent?
# ifdef MDSPLUS
    logical loc_status
    integer istat
    equivalence (loc_status,istat)
    integer length
    integer :: descr
    loc_status = mdsvalue("getmsg($)",descr(IDTYPE_LONG,status,0),&
         &descr(IDTYPE_CSTRING,text,0,len(text),0),0,length)
    if (iand(istat,1) .eq. 0) then
       loc_status = status
       write(text,*) "error status = ",istat
    endif
    return
# endif
  end subroutine getmdserrortext


!--------------------------------------------
! mds_read subroutines
!--------------------------------------------

  subroutine mds_r_log_0(name,value)
    character(*), intent(in) :: name
    logical, intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    status = mdsvalue(name//char(0),descr(IDTYPE_LONG,value,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_log_0

  subroutine mds_r_char_1(name,value)
    character(*), intent(in) :: name
    character*(*) , intent(out):: value
# ifdef MDSPLUS
    logical :: status
    integer :: clen 
    integer :: descr
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_CSTRING,value,0,len(value)),0,clen)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_char_1

  subroutine mds_r_int_0(name,value)
    character(*), intent(in) :: name
    integer , intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    status = mdsvalue(name//char(0),descr(IDTYPE_LONG,value,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_int_0

  subroutine mds_r_int_1(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsvalue(name//char(0),descr(IDTYPE_LONG,value,dim1,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_int_1

  subroutine mds_r_int_2(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsvalue(name//char(0),descr(IDTYPE_LONG,value,dim1,dim2,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_int_2

  subroutine mds_r_int_3(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsvalue(name//char(0), &
         descr(IDTYPE_LONG,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_int_3

  subroutine mds_r_int_4(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_LONG,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_int_4

  subroutine mds_r_int_5(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_LONG,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_int_5

  subroutine mds_r_int_6(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_LONG,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_int_6

  subroutine mds_r_int_7(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_LONG,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_int_7

  subroutine mds_r_real_0(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    status = mdsvalue(name//char(0),descr(IDTYPE_FLOAT,value,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_real_0

  subroutine mds_r_real_1(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsvalue(name//char(0),descr(IDTYPE_FLOAT,value,dim1,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_real_1

  subroutine mds_r_real_2(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_FLOAT,value,dim1,dim2,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_real_2

  subroutine mds_r_real_3(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_FLOAT,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_real_3

  subroutine mds_r_real_4(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_FLOAT,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_real_4

  subroutine mds_r_real_5(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_FLOAT,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_real_5

  subroutine mds_r_real_6(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:,:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_FLOAT,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_real_6

  subroutine mds_r_real_7(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:,:,:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_FLOAT,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_real_7

  subroutine mds_r_doub_0(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd) , intent(out):: value
!    double precision, intent(out):: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    status = mdsvalue(name//char(0),descr(IDTYPE_DOUBLE,value,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_doub_0

  subroutine mds_r_doub_1(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsvalue(name//char(0),descr(IDTYPE_DOUBLE,value,dim1,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_doub_1

  subroutine mds_r_doub_2(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_DOUBLE,value,dim1,dim2,0,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_doub_2

  subroutine mds_r_doub_3(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_DOUBLE,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_doub_3

  subroutine mds_r_doub_4(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_DOUBLE,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_doub_4

  subroutine mds_r_doub_5(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_DOUBLE,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_doub_5

  subroutine mds_r_doub_6(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:,:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_DOUBLE,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_doub_6

  subroutine mds_r_doub_7(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:,:,:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_DOUBLE,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_doub_7

  subroutine mds_r_cmplx_0(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    status = mdsvalue(name//char(0),descr(IDTYPE_COMPLEX,value,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_cmplx_0

  subroutine mds_r_cmplx_1(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsvalue(name//char(0),descr(IDTYPE_COMPLEX,value,dim1,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_cmplx_1

  subroutine mds_r_cmplx_2(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_COMPLEX,value,dim1,dim2,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_cmplx_2

  subroutine mds_r_cmplx_3(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_COMPLEX,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_cmplx_3

  subroutine mds_r_cmplx_4(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_COMPLEX,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_cmplx_4

  subroutine mds_r_cmplx_5(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_COMPLEX,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_cmplx_5

  subroutine mds_r_cmplx_6(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:,:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_COMPLEX,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_cmplx_6

  subroutine mds_r_cmplx_7(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:,:,:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_COMPLEX,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_cmplx_7

  subroutine mds_r_dcmplx_0(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    status = mdsvalue(name//char(0),descr(IDTYPE_DOUBLE_COMPLEX,value,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_dcmplx_0

  subroutine mds_r_dcmplx_1(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsvalue(name//char(0),descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_dcmplx_1

  subroutine mds_r_dcmplx_2(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,dim2,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_dcmplx_2

  subroutine mds_r_dcmplx_3(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_dcmplx_3

  subroutine mds_r_dcmplx_4(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_dcmplx_4

  subroutine mds_r_dcmplx_5(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_dcmplx_5

  subroutine mds_r_dcmplx_6(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:,:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_dcmplx_6

  subroutine mds_r_dcmplx_7(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:,:,:,:,:,:,:), intent(out) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsvalue(name//char(0),&
         &descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error reading "//name)
    return
# endif
  end subroutine mds_r_dcmplx_7

!--------------------------------------------
! mds_write subroutines
!--------------------------------------------

  subroutine mds_w_log_0(name,value)
    character(*), intent(in) :: name
    logical :: value ! intent?
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",descr(IDTYPE_LONG,value,0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_log_0

  subroutine mds_w_char_1(name,value)
    character(*), intent(in) :: name
    character*(*), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",&
         &descr(IDTYPE_CSTRING,value,0,len(value)),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_char_1

  subroutine mds_w_int_0(name,value)
    character(*), intent(in) :: name
    integer , intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",descr(IDTYPE_LONG ,value, 0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_int_0

  subroutine mds_w_int_1(name,value)
    character(*), intent(in) :: name
    integer , dimension(:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsput(name//char(0),"$",descr(IDTYPE_LONG ,value, dim1, 0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_int_1

  subroutine mds_w_int_2(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_LONG,value,dim1,dim2,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_int_2

  subroutine mds_w_int_3(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_LONG,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_int_3

  subroutine mds_w_int_4(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_LONG,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_int_4

  subroutine mds_w_int_5(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_LONG,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_int_5

  subroutine mds_w_int_6(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_LONG,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_int_6

  subroutine mds_w_int_7(name,value)
    character(*), intent(in) :: name 
    integer , dimension(:,:,:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_LONG,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_int_7

  subroutine mds_w_real_0(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",descr(IDTYPE_FLOAT,value,0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_real_0

  subroutine mds_w_real_1(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs), dimension(:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsput(name//char(0),"$",descr(IDTYPE_FLOAT,value,dim1,0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_real_1

  subroutine mds_w_real_2(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_FLOAT,value,dim1,dim2,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_real_2

  subroutine mds_w_real_3(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_FLOAT,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_real_3

  subroutine mds_w_real_4(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_FLOAT,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_real_4

  subroutine mds_w_real_5(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_FLOAT,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_real_5

  subroutine mds_w_real_6(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:,:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_FLOAT,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_real_6

  subroutine mds_w_real_7(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rs) , dimension(:,:,:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_FLOAT,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_real_7

  subroutine mds_w_doub_0(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",descr(IDTYPE_DOUBLE,value,0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_doub_0

  subroutine mds_w_doub_1(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsput(name//char(0),"$",descr(IDTYPE_DOUBLE,value,dim1,0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_doub_1

  subroutine mds_w_doub_2(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsput(name//char(0),"$",&
         &descr(IDTYPE_DOUBLE,value,dim1,dim2,0,0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_doub_2

  subroutine mds_w_doub_3(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsput(name//char(0),"$",&
         &descr(IDTYPE_DOUBLE,value,dim1,dim2,dim3,0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_doub_3

  subroutine mds_w_doub_4(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_DOUBLE,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_doub_4

  subroutine mds_w_doub_5(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_DOUBLE,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_doub_5

  subroutine mds_w_doub_6(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:,:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_DOUBLE,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_doub_6

  subroutine mds_w_doub_7(name,value)
    character(*), intent(in) :: name
    real(kind=kind_rd), dimension(:,:,:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_DOUBLE,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_doub_7

  subroutine mds_w_cmplx_0(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",descr(IDTYPE_COMPLEX,value,0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_cmplx_0

  subroutine mds_w_cmplx_1(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsput(name//char(0),"$",descr(IDTYPE_COMPLEX,value,dim1,0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_cmplx_1

  subroutine mds_w_cmplx_2(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_COMPLEX,value,dim1,dim2,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_cmplx_2

  subroutine mds_w_cmplx_3(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_COMPLEX,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_cmplx_3

  subroutine mds_w_cmplx_4(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_COMPLEX,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_cmplx_4

  subroutine mds_w_cmplx_5(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_COMPLEX,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_cmplx_5

  subroutine mds_w_cmplx_6(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:,:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_COMPLEX,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_cmplx_6

  subroutine mds_w_cmplx_7(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rs), dimension(:,:,:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_COMPLEX,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_cmplx_7

  subroutine mds_w_dcmplx_0(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    status = mdsput(name//char(0),"$",descr(IDTYPE_DOUBLE_COMPLEX,value,0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_dcmplx_0

  subroutine mds_w_dcmplx_1(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: descr
    integer :: dim1
    dim1 = size(value,1)
    status = mdsput(name//char(0),"$",descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,0),0)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_dcmplx_1

  subroutine mds_w_dcmplx_2(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2
    dim1 = size(value,1)
    dim2 = size(value,2)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,dim2,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_dcmplx_2

  subroutine mds_w_dcmplx_3(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,dim2,dim3,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_dcmplx_3

  subroutine mds_w_dcmplx_4(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,dim2,dim3,dim4,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_dcmplx_4

  subroutine mds_w_dcmplx_5(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,dim2,dim3,dim4,dim5,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_dcmplx_5

  subroutine mds_w_dcmplx_6(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:,:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,dim2,dim3,dim4,dim5,dim6,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_dcmplx_6

  subroutine mds_w_dcmplx_7(name,value)
    character(*), intent(in) :: name
    complex(kind=kind_rd), dimension(:,:,:,:,:,:,:), intent(in) :: value
# ifdef MDSPLUS
    logical :: status
    integer :: len
    integer :: descr
    integer :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    dim1 = size(value,1)
    dim2 = size(value,2)
    dim3 = size(value,3)
    dim4 = size(value,4)
    dim5 = size(value,5)
    dim6 = size(value,6)
    dim7 = size(value,7)
    status = mdsput(name//char(0),&
         &"$",descr(IDTYPE_DOUBLE_COMPLEX,value,dim1,dim2,dim3,dim4,dim5,dim6,dim7,0),0,len)
    call checkmds(status,"error writing "//name)
    return
# endif
  end subroutine mds_w_dcmplx_7

end module mdsio

