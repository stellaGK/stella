# include "define.inc"

module mp
!
! <doc> Easier Fortran90 interface to the MPI Message Passing Library. </doc>
!
!     (c) Copyright 1991 to 1998 by Michael A. Beer, William D. Dorland, 
!     P. B. Snyder, Q. P. Liu, and Gregory W. Hammett. ALL RIGHTS RESERVED.
!     
! Note: mp_mpi_r8.f90 is a version of mp_mpi.f90 to use when compiling 
! with -r8 (where the default real type is taken to be 8 bytes).  Just 
! replaced all occurances of MPI_REAL with MPI_DOUBLE_PRECISION and 
! MPI_COMPLEX with MPI_DOUBLE_COMPLEX.
!
# ifdef MPI
# ifndef MPIINC
  use mpi
#endif
  ! TT: I experienced a problem on Ranger with mvapich-devel.
  ! TT: Compiler complained an inconsistent variable type for reorder below.
  ! TT: In this case, the problem was solved by commenting out the above line
  ! TT: and use mpif.h below.  (4/16/08)
# endif
  implicit none
  private

  public :: init_mp, finish_mp
  public :: broadcast, sum_reduce, sum_allreduce
  public :: max_reduce, max_allreduce
  public :: min_reduce, min_allreduce
  public :: nproc, iproc, proc0, job
  public :: send, ssend, receive
  public :: barrier
  public :: waitany
  public :: test_driver_flag
! JH> new abort method
  public :: mp_abort
! <JH
! MAB> needed by Trinity
  public :: scope, allprocs, subprocs
  public :: all_to_group, group_to_all
  public :: trin_flag
! <MAB
  public :: init_jobs
  public :: mp_comm
  public :: mp_info
  public :: mp_gather

# ifdef MPI
# ifdef MPIINC
! CMR: defined MPIINC for machines where need to include mpif.h
  include 'mpif.h'
#endif
  integer, pointer :: nproc
  integer, target :: ntot_proc, ngroup_proc

  integer, pointer :: iproc
  integer, target :: aproc, gproc

  logical, pointer :: proc0
  logical, target :: aproc0, gproc0

  integer, pointer :: mp_comm
  integer, target :: comm_all, comm_group

  integer, parameter :: mp_info = MPI_INFO_NULL

  integer :: job = 0
  integer (kind(MPI_REAL)) :: mpireal, mpicmplx
# else
  integer, parameter :: nproc = 1, iproc = 0
  logical, parameter :: proc0 = .true.

  integer, parameter :: mp_info = -1
  integer, parameter :: job = 0, mp_comm = -1
# endif
  integer, parameter :: allprocs = 0, subprocs = 1

! needed for Trinity -- MAB
  integer, dimension (:), allocatable :: grp0
  logical :: trin_flag = .false.

  logical :: test_driver_flag = .false.

  interface broadcast
     module procedure broadcast_integer 
     module procedure broadcast_integer_array 

     module procedure broadcast_real    
     module procedure broadcast_real_array   
     module procedure broadcast_real_2array 
     module procedure broadcast_real_3array 

     module procedure broadcast_complex 
     module procedure broadcast_complex_array
     module procedure broadcast_complex_2array
     module procedure broadcast_complex_3array

     module procedure broadcast_logical 
     module procedure broadcast_logical_array 

     module procedure bcastfrom_integer 
     module procedure bcastfrom_integer_array 

     module procedure bcastfrom_real    
     module procedure bcastfrom_real_array    

     module procedure bcastfrom_complex 
     module procedure bcastfrom_complex_array 

     module procedure bcastfrom_logical 
     module procedure bcastfrom_logical_array 

     module procedure broadcast_character
     module procedure bcastfrom_character
  end interface

  interface sum_reduce
     module procedure sum_reduce_integer
     module procedure sum_reduce_integer_array

     module procedure sum_reduce_real
     module procedure sum_reduce_real_array
     module procedure sum_reduce_real_2array
     module procedure sum_reduce_real_3array
     module procedure sum_reduce_real_4array
     module procedure sum_reduce_real_5array

     module procedure sum_reduce_complex
     module procedure sum_reduce_complex_array
     module procedure sum_reduce_complex_2array
     module procedure sum_reduce_complex_3array
     module procedure sum_reduce_complex_4array
     module procedure sum_reduce_complex_5array
  end interface
  
  !KDN 100526: Allows summing into alternate variable
  !rather than overwriting local data
!  interface sum_reduce_alt
!     module procedure sum_reduce_alt_complex_3array
!  end interface

  interface sum_allreduce
     module procedure sum_allreduce_integer
     module procedure sum_allreduce_integer_array

     module procedure sum_allreduce_real
     module procedure sum_allreduce_real_array
     module procedure sum_allreduce_real_2array
     module procedure sum_allreduce_real_3array
     module procedure sum_allreduce_real_4array
     module procedure sum_allreduce_real_5array

     module procedure sum_allreduce_complex
     module procedure sum_allreduce_complex_array
     module procedure sum_allreduce_complex_2array
     module procedure sum_allreduce_complex_3array
     module procedure sum_allreduce_complex_4array
     module procedure sum_allreduce_complex_5array
  end interface

  interface max_reduce
     module procedure max_reduce_integer
     module procedure max_reduce_integer_array

     module procedure max_reduce_real
     module procedure max_reduce_real_array
  end interface

  interface max_allreduce
     module procedure max_allreduce_integer
     module procedure max_allreduce_integer_array

     module procedure max_allreduce_real
     module procedure max_allreduce_real_array
  end interface

  interface min_reduce
     module procedure min_reduce_integer
     module procedure min_reduce_integer_array

     module procedure min_reduce_real
     module procedure min_reduce_real_array
  end interface

  interface min_allreduce
     module procedure min_allreduce_integer
     module procedure min_allreduce_integer_array

     module procedure min_allreduce_real
     module procedure min_allreduce_real_array
  end interface

  interface send
     module procedure send_integer
     module procedure send_integer_array

     module procedure send_real
     module procedure send_real_array

     module procedure send_complex
     module procedure send_complex_array
     module procedure nonblocking_send_complex_array

     module procedure send_logical
     module procedure send_logical_array
     
     module procedure send_character
  end interface

  interface receive
     module procedure receive_integer
     module procedure receive_integer_array

     module procedure receive_real
     module procedure receive_real_array

     module procedure receive_complex
     module procedure receive_complex_array
     module procedure receive_complex_2array
     module procedure nonblocking_receive_complex_array

     module procedure receive_logical
     module procedure receive_logical_array

     module procedure receive_character
  end interface

! MAB> needed for Trinity
! synchronous sends
  interface ssend
     module procedure ssend_integer
     module procedure ssend_integer_array

     module procedure ssend_real
     module procedure ssend_real_array

     module procedure ssend_complex
     module procedure ssend_complex_array
     module procedure ssend_complex_2array

     module procedure ssend_logical
     module procedure ssend_logical_array
  end interface

! send stuff from global proc0 to group proc0s
  interface all_to_group
     module procedure all_to_group_real
     module procedure all_to_group_real_array
  end interface

! send stuff from group proc0s to global proc0
  interface group_to_all
     module procedure group_to_all_real
     module procedure group_to_all_real_array
  end interface
! <MAB

contains

!MAB> changed to allow for calling within trinity (coupled flux tubes)
!  subroutine init_mp
  subroutine init_mp (comm_in)
# ifdef MPI
    use constants, only: pi, kind_rs, kind_rd
    use file_utils, only: error_unit
    implicit none
# endif
    integer, intent (in), optional :: comm_in ! MAB
# ifdef MPI
    integer :: ierror!, rank
    logical :: init ! MAB

    call mpi_initialized (init, ierror)  ! MAB
    if (.not. init) call mpi_init (ierror)  ! MAB
!    call mpi_init (ierror)
!MAB>
    if (present(comm_in)) then
       comm_all = comm_in
       if (.not. test_driver_flag) trin_flag = .true. 
       !EGH  we can pass in a communicator if we want to run tests
    else
       comm_all = mpi_comm_world
    end if
    call mpi_comm_size (comm_all, ntot_proc, ierror)
    call mpi_comm_rank (comm_all, aproc, ierror)
!    call mpi_comm_size (mpi_comm_world, ntot_proc, ierror)
!    call mpi_comm_rank (mpi_comm_world, aproc, ierror)
!    comm_all = mpi_comm_world
!<MAB
    aproc0 = aproc == 0

    call scope (allprocs)

    if ( (kind(pi)==kind_rs) .and. (kind_rs/=kind_rd) ) then
       mpireal = MPI_REAL
       mpicmplx = MPI_COMPLEX
    else if (kind(pi)==kind_rd) then
       mpireal = MPI_DOUBLE_PRECISION
       mpicmplx = MPI_DOUBLE_COMPLEX
    else
       write (error_unit(),*) 'ERROR: precision mismatch in mpi'
    end if
# endif

  end subroutine init_mp

  subroutine scope (focus)

    integer, intent (in) :: focus

# ifdef MPI
    if (focus == allprocs) then
       mp_comm => comm_all
       nproc => ntot_proc
       iproc => aproc
       proc0 => aproc0
    else
       mp_comm => comm_group
       nproc => ngroup_proc
       iproc => gproc
       proc0 => gproc0
    end if
# endif

  end subroutine scope

  subroutine finish_mp
# ifdef MPI
    implicit none
    integer :: ierror

    call mpi_finalize (ierror)
# endif
  end subroutine finish_mp

! ************** broadcasts *****************************

  subroutine broadcast_character (char)
    implicit none
    character(*), intent (in out) :: char
# ifdef MPI
    integer :: ierror
    call mpi_bcast (char, len(char), MPI_CHARACTER, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_character

  subroutine broadcast_integer (i)
    implicit none
    integer, intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_integer

  subroutine broadcast_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_integer_array

  subroutine broadcast_real (x)
    implicit none
    real, intent (in out) :: x
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, 1, mpireal, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_real

  subroutine broadcast_real_array (x)
    implicit none
    real, dimension (:), intent (in out) :: x
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, size(x), mpireal, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_real_array
  
  subroutine broadcast_real_2array(x)
  implicit none
  real, dimension(:,:), intent (in out) :: x
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, size(x), mpireal, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_real_2array

  subroutine broadcast_real_3array(x)
  implicit none
  real, dimension(:,:,:), intent (in out) :: x
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, size(x), mpireal, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_real_3array

  subroutine broadcast_complex (z)
    implicit none
    complex, intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, 1, mpicmplx, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_complex

  subroutine broadcast_complex_array (z)
    implicit none
    complex, dimension (:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_complex_array
  
  subroutine broadcast_complex_2array (z)
    implicit none
    complex, dimension (:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_complex_2array

  subroutine broadcast_complex_3array (z)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_complex_3array

  subroutine broadcast_logical (f)
    implicit none
    logical, intent (in out) :: f
# ifdef MPI
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_logical

  subroutine broadcast_logical_array (f)
    implicit none
    logical, dimension (:), intent (in out) :: f
# ifdef MPI
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, 0, mp_comm, ierror)
# endif
  end subroutine broadcast_logical_array

  subroutine bcastfrom_logical (f, src)
    implicit none
    logical, intent (in out) :: f
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, src, mp_comm, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_logical

  subroutine bcastfrom_logical_array (f, src)
    implicit none
    logical, dimension (:), intent (in out) :: f
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, src, mp_comm, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_logical_array

  subroutine bcastfrom_character (c, src)
    implicit none
    character(*), intent (in out) :: c
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (c, len(c), MPI_CHARACTER, src, mp_comm, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_character

  subroutine bcastfrom_integer (i, src)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, src, mp_comm, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_integer

  subroutine bcastfrom_integer_array (i, src)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, src, mp_comm, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_integer_array

  subroutine bcastfrom_real (x, src)
    implicit none
    real, intent (in out) :: x
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, 1, mpireal, src, mp_comm, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_real

  subroutine bcastfrom_real_array (x, src)
    implicit none
    real, dimension (:), intent (in out) :: x
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (x, size(x), mpireal, src, mp_comm, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_real_array

  subroutine bcastfrom_complex (z, src)
    implicit none
    complex, intent (in out) :: z
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, 1, mpicmplx, src, mp_comm, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_complex

  subroutine bcastfrom_complex_array (z, src)
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: src
# ifdef MPI
    integer :: ierror
    call mpi_bcast (z, size(z), mpicmplx, src, mp_comm, ierror)
# else
    if (src /= 0) call error ("broadcast from")
# endif
  end subroutine bcastfrom_complex_array

! ************** reductions ***********************

  subroutine sum_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (i, i, 1, MPI_INTEGER, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_integer

  subroutine sum_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (i, i, size(i), MPI_INTEGER, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_integer_array

!   subroutine sum_reduce_logical (a, dest)
!     implicit none
!     logical, intent (in out) :: a
!     integer, intent (in) :: dest
! # ifdef MPI
!     integer :: ierror
!     if(iproc.eq.dest)then
!        call mpi_reduce &
!             (MPI_IN_PLACE, a, 1, MPI_LOGICAL, MPI_LOR, dest, mp_comm, ierror)
!     else
!        call mpi_reduce &
!             (a, a, 1, MPI_LOGICAL, MPI_LOR, dest, mp_comm, ierror)
!     endif
! # else
!     if (dest /= 0) call error ("reduce to")
! # endif
!   end subroutine sum_reduce_logical

  subroutine sum_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, a, 1, mpireal, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (a, a, 1, mpireal, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real

  subroutine sum_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
         (a, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real_array

  subroutine sum_reduce_real_2array (a, dest)
    implicit none
    real, dimension (:,:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
         (a, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real_2array

  subroutine sum_reduce_real_3array (a, dest)
    implicit none
    real, dimension (:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
         (a, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real_3array

  subroutine sum_reduce_real_4array (a, dest)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
         (a, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real_4array

  subroutine sum_reduce_real_5array (a, dest)
    implicit none
    real, dimension (:,:,:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
         (a, a, size(a), mpireal, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_real_5array

  subroutine sum_reduce_complex (z, dest)
    implicit none
    complex, intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, 1, mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (z, z, 1, mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex

  subroutine sum_reduce_complex_array (z, dest)
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (z, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_array

  subroutine sum_reduce_complex_2array (z, dest)
    implicit none
    complex, dimension (:,:), intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (z, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_2array

  subroutine sum_reduce_complex_3array (z, dest)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (z, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_3array

  subroutine sum_reduce_complex_4array (z, dest)
    implicit none
    complex, dimension (:,:,:,:), intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (z, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_4array

  subroutine sum_reduce_complex_5array (z, dest)
    implicit none
    complex, dimension (:,:,:,:,:), intent (in out) :: z
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if (iproc.eq.dest) then
       call mpi_reduce &
            (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (z, z, size(z), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine sum_reduce_complex_5array

! Sum all z1 values into z2 at dest
!   subroutine sum_reduce_alt_complex_3array (z1,z2,dest)
!     implicit none
!     complex, dimension (:,:,:), intent (in out) :: z1
!     complex, dimension (:,:,:), intent (in out) :: z2
!     integer, intent (in) :: dest
! # ifdef MPI
!     integer :: ierror
!     call mpi_reduce &
!          (z1, z2, size(z1), mpicmplx, MPI_SUM, dest, mp_comm, ierror)
! # else
!     if (dest /= 0) call error ("reduce to")
! # endif
!   end subroutine sum_reduce_alt_complex_3array
  
  subroutine sum_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_integer

  subroutine sum_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_integer_array

  subroutine sum_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, 1, mpireal, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_real

  subroutine sum_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_real_array

  subroutine sum_allreduce_real_2array (a)
    implicit none
    real, dimension (:,:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_real_2array

  subroutine sum_allreduce_real_3array (a)
    implicit none
    real, dimension (:,:,:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_real_3array

  subroutine sum_allreduce_real_4array (a)
    implicit none
    real, dimension (:,:,:,:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_real_4array

  subroutine sum_allreduce_real_5array (a)
    implicit none
    real, dimension (:,:,:,:,:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_real_5array

  subroutine sum_allreduce_complex (z)
    implicit none
    complex, intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, 1, mpicmplx, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_complex

  subroutine sum_allreduce_complex_array (z)
    implicit none
    complex, dimension (:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_complex_array

  subroutine sum_allreduce_complex_2array (z)
    implicit none
    complex, dimension (:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_complex_2array

  subroutine sum_allreduce_complex_3array (z)
    implicit none
    complex, dimension (:,:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_complex_3array

  subroutine sum_allreduce_complex_4array (z)
    implicit none
    complex, dimension (:,:,:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_complex_4array

  subroutine sum_allreduce_complex_5array (z)
    implicit none
    complex, dimension (:,:,:,:,:), intent (in out) :: z
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, z, size(z), mpicmplx, MPI_SUM, mp_comm, ierror)
# endif
  end subroutine sum_allreduce_complex_5array

  subroutine max_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_MAX, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (i, i, 1, MPI_INTEGER, MPI_MAX, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine max_reduce_integer

  subroutine max_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_MAX, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (i, i, size(i), MPI_INTEGER, MPI_MAX, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine max_reduce_integer_array

  subroutine max_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, a, 1, mpireal, MPI_MAX, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (a, a, 1, mpireal, MPI_MAX, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine max_reduce_real

  subroutine max_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, a, size(a), mpireal, MPI_MAX, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (a, a, size(a), mpireal, MPI_MAX, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine max_reduce_real_array

  subroutine max_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_MAX, mp_comm, ierror)
# endif
  end subroutine max_allreduce_integer

  subroutine max_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_MAX, mp_comm, ierror)
# endif
  end subroutine max_allreduce_integer_array

  subroutine max_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, 1, mpireal, MPI_MAX, mp_comm, ierror)
# endif
  end subroutine max_allreduce_real

  subroutine max_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_MAX, mp_comm, ierror)
# endif
  end subroutine max_allreduce_real_array

  subroutine min_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_MIN, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (i, i, 1, MPI_INTEGER, MPI_MIN, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine min_reduce_integer

  subroutine min_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_MIN, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (i, i, size(i), MPI_INTEGER, MPI_MIN, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine min_reduce_integer_array

  subroutine min_reduce_real (a, dest)
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, a, 1, mpireal, MPI_MIN, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (a, a, 1, mpireal, MPI_MIN, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine min_reduce_real

  subroutine min_reduce_real_array (a, dest)
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
# ifdef MPI
    integer :: ierror
    if(iproc.eq.dest)then
       call mpi_reduce &
            (MPI_IN_PLACE, a, size(a), mpireal, MPI_MIN, dest, mp_comm, ierror)
    else
       call mpi_reduce &
            (a, a, size(a), mpireal, MPI_MIN, dest, mp_comm, ierror)
    endif
# else
    if (dest /= 0) call error ("reduce to")
# endif
  end subroutine min_reduce_real_array

  subroutine min_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_MIN, mp_comm, ierror)
# endif
  end subroutine min_allreduce_integer

  subroutine min_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_MIN, mp_comm, ierror)
# endif
  end subroutine min_allreduce_integer_array

  subroutine min_allreduce_real (a)
    implicit none
    real, intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, 1, mpireal, MPI_MIN, mp_comm, ierror)
# endif
  end subroutine min_allreduce_real

  subroutine min_allreduce_real_array (a)
    implicit none
    real, dimension (:), intent (in out) :: a
# ifdef MPI
    integer :: ierror
    call mpi_allreduce &
         (MPI_IN_PLACE, a, size(a), mpireal, MPI_MIN, mp_comm, ierror)
# endif
  end subroutine min_allreduce_real_array

! ********************* barrier **********************

  subroutine barrier
# ifdef MPI
    implicit none
    integer :: ierror
    call mpi_barrier (mp_comm, ierror)
# endif
  end subroutine barrier

! ********************* sends **********************

  subroutine send_integer (i, dest, tag)
    implicit none
    integer, intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, 1, MPI_INTEGER, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine send_integer

  subroutine send_integer_array (i, dest, tag)
    implicit none
    integer, dimension (:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, size(i), MPI_INTEGER, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine send_integer_array

  subroutine send_real (a, dest, tag)
    implicit none
    real, intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, 1, mpireal, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine send_real

  subroutine send_real_array (a, dest, tag)
    implicit none
    real, dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, size(a), mpireal, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine send_real_array

  subroutine send_complex (z, dest, tag)
    implicit none
    complex, intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, 1, mpicmplx, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine send_complex

  subroutine send_complex_array (z, dest, tag)
    implicit none
    complex, dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, size(z), mpicmplx, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine send_complex_array

  subroutine nonblocking_send_complex_array (z, dest, tag, request)
    implicit none
    complex, dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer, intent (out) :: request
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (z, size(z), mpicmplx, dest, tagp, mp_comm, request, ierror)
# else
    call error ("send")
# endif
  end subroutine nonblocking_send_complex_array

  subroutine send_logical (f, dest, tag)
    implicit none
    logical, intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, 1, MPI_LOGICAL, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine send_logical

  subroutine send_logical_array (f, dest, tag)
    implicit none
    logical, dimension (:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, size(f), MPI_LOGICAL, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine send_logical_array

  subroutine send_character (s, dest, tag)
    implicit none
    character(*), intent (in) :: s
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send &
         (s, len(s), MPI_CHARACTER, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine send_character

! MAB> needed for Trinity
! ********************* synchronous sends **********************

  subroutine ssend_integer (i, dest, tag)
    implicit none
    integer, intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (i, 1, MPI_INTEGER, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_integer

  subroutine ssend_integer_array (i, dest, tag)
    implicit none
    integer, dimension (:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (i, size(i), MPI_INTEGER, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_integer_array

  subroutine ssend_real (a, dest, tag)
    implicit none
    real, intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (a, 1, MPI_DOUBLE_PRECISION, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_real

  subroutine ssend_real_array (a, dest, tag)
    implicit none
    real, dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (a, size(a), MPI_DOUBLE_PRECISION, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_real_array

  subroutine ssend_complex (z, dest, tag)
    implicit none
    complex, intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (z, 1, MPI_DOUBLE_COMPLEX, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_complex

  subroutine ssend_complex_array (z, dest, tag)
    implicit none
    complex, dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_complex_array
  
  subroutine ssend_complex_2array (z, dest, tag)
    implicit none
    complex, dimension (:,:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_complex_2array
  
  subroutine ssend_logical (f, dest, tag)
    implicit none
    logical, intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (f, 1, MPI_LOGICAL, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_logical

  subroutine ssend_logical_array (f, dest, tag)
    implicit none
    logical, dimension (:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_ssend (f, size(f), MPI_LOGICAL, dest, tagp, mp_comm, ierror)
# else
    call error ("send")
# endif
  end subroutine ssend_logical_array

!   subroutine ssend_character (s, dest, tag)
!     implicit none
!     character(*), intent (in) :: s
!     integer, intent (in) :: dest
!     integer, intent (in), optional :: tag
! # ifdef MPI
!     integer :: ierror
!     integer :: tagp
!     tagp = 0
!     if (present(tag)) tagp = tag
!     call mpi_ssend &
!          (s, len(s), MPI_CHARACTER, dest, tagp, mp_comm, ierror)
! # else
!     call error ("send")
! # endif
!   end subroutine ssend_character
! <MAB

! ********************* receives  **********************

  subroutine receive_integer (i, src, tag)
    implicit none
# ifdef MPI
    integer, intent (out) :: i
# else
    integer :: i
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, 1, MPI_INTEGER, src, tagp, mp_comm, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_integer

  subroutine receive_integer_array (i, src, tag)
    implicit none
# ifdef MPI
    integer, dimension (:), intent (out) :: i
# else
    integer, dimension (:) :: i
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, size(i), MPI_INTEGER, src, tagp, mp_comm, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_integer_array

  subroutine receive_real (a, src, tag)
    implicit none
# ifdef MPI
    real, intent (out) :: a
# else
    real :: a
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, 1, mpireal, src, tagp, mp_comm, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_real

  subroutine receive_real_array (a, src, tag)
    implicit none
# ifdef MPI
    real, dimension (:), intent (out) :: a
# else
    real, dimension (:) :: a
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, size(a), mpireal, src, tagp, mp_comm, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_real_array

  subroutine receive_complex (z, src, tag)
    implicit none
# ifdef MPI
    complex, intent (out) :: z
# else
    complex :: z
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, 1, mpicmplx, src, tagp, mp_comm, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_complex

  subroutine receive_complex_array (z, src, tag)
    implicit none
# ifdef MPI
    complex, dimension (:), intent (out) :: z
# else
    complex, dimension (:) :: z
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, size(z), mpicmplx, src, tagp, mp_comm, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_complex_array

  subroutine receive_complex_2array (z, src, tag)
    implicit none
# ifdef MPI
    complex, dimension (:,:), intent (out) :: z
# else
    complex, dimension (:,:) :: z
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, size(z), mpicmplx, src, tagp, mp_comm, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_complex_2array
  
  subroutine nonblocking_receive_complex_array (z, src, tag, request)
    implicit none
# ifdef MPI
    complex, dimension (:), intent (inout) :: z
# else
    complex, dimension (:) :: z
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer, intent (out) :: request
# ifdef MPI
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (z, size(z), mpicmplx, src, tagp, mp_comm, &
        request, ierror)
# else
    call error ("receive")
# endif
  end subroutine nonblocking_receive_complex_array

  subroutine receive_logical (f, src, tag)
    implicit none
# ifdef MPI
    logical, intent (out) :: f
# else
    logical :: f
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, 1, MPI_LOGICAL, src, tagp, mp_comm, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_logical

  subroutine receive_logical_array (f, src, tag)
    implicit none
# ifdef MPI
    logical, dimension (:), intent (out) :: f
# else
    logical, dimension (:) :: f
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, size(f), MPI_LOGICAL, src, tagp, mp_comm, &
        status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_logical_array

  subroutine receive_character (s, src, tag)
    implicit none
# ifdef MPI
    character(*), intent (out) :: s
# else
    character(*) :: s
# endif
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
# ifdef MPI
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (s, len(s), MPI_CHARACTER, src, tagp, mp_comm, &
         status, ierror)
# else
    call error ("receive")
# endif
  end subroutine receive_character

  subroutine waitany (count, requests, requestindex, status)

    implicit none
    integer, intent(in) :: count
    integer, dimension(:), intent(inout) :: requests
    integer, intent(out) :: requestindex
# ifdef MPI
    integer, dimension(MPI_STATUS_SIZE), intent(out) :: status
# else
    integer, dimension(1), intent(out) :: status
# endif

# ifdef MPI
    integer :: ierror

    call mpi_waitany(count, requests, requestindex, status, ierror)
# else
    call error ("waitany")
# endif

  end subroutine waitany

  subroutine init_jobs (ncolumns, group0, ierr)

    implicit none  
# ifdef MPI
!    integer, parameter :: reorder=1
    ! TT: I changed variable definition by assuming integer 1 corresponds to
    ! TT: logical .true. but I am not sure if reorder is needed.
    ! TT: In any case this subroutine is only called when you use job fork.
    logical, parameter :: reorder=.true.
    integer :: ip, j, comm2d, id2d, ierr, nrows
# endif
    integer, intent(in) :: ncolumns
    integer, dimension(0:), intent (out) :: group0
# ifndef MPI
    integer :: ierr
    if (ncolumns /= 1) call error ("jobs")
# else
    integer, parameter :: ndim=2
    integer, dimension(ndim) :: dims
    integer, dimension(0:ndim-1) :: coords1d, coords2d
    logical, dimension(0:ndim-1) :: belongs
    logical, dimension(ndim) :: period
    
    logical :: isroot

    if (.not. allocated(grp0)) allocate (grp0(0:size(group0)-1))
    
! calculate dimensions  mpi processor grid will have and check that 
! ncolumns*nrows = number of processes

! nrows is # of processors per job (or group)    
    nrows = ntot_proc/ncolumns
    dims=(/ ncolumns, nrows /)     
    if(ntot_proc /= ncolumns*nrows) then
       ierr = 1
       if(aproc0) write(*,*) 'Number of processes must be divisible by number of groups'
       return
    endif
    ngroup_proc = nrows
    
    ! create 2d cartesian topology for processes
    
    period=(/ .false., .false. /)  !! no circular shift

!    call mpi_cart_create(mpi_comm_world, ndim, dims, period, reorder, comm2d, ierr)
    call mpi_cart_create(comm_all, ndim, dims, period, reorder, comm2d, ierr)
    call mpi_comm_rank(comm2d, id2d, ierr)
    call mpi_cart_coords(comm2d, id2d, ndim, coords2d, ierr)
    
! each processor knows which subgrid it is in from variable mpi_group
    job = coords2d(0)

! create 1d subgrids from 2d processor grid, variable belongs denotes
! whether processor grid is split by column or row

    belongs(1) = .true.    ! this dimension belongs to subgrid
    belongs(0) = .false.  

    call mpi_cart_sub(comm2d, belongs, comm_group, ierr)
    call mpi_comm_rank(comm_group, gproc, ierr)     
    call mpi_cart_coords(comm_group, gproc, 1, coords1d, ierr)
    gproc0 = (gproc == 0)
    
! find root process of each 1d subgrid and place in array group0 indexed 
! from 0 to subgrids-1
     
! MAB> following two lines were incorrect
!    j=1
!    group0(0) = 0
! replace with
    j = 0
    if (proc0 .and. gproc0) then
       group0(0) = 0
       j = 1
    end if
! <MAB
    do ip = 1, ntot_proc-1
       if (proc0) then
          call receive (isroot, ip)
          if (isroot) then
             group0(j) = ip
             j=j+1
          end if
       else if (ip == aproc) then
          call send (gproc0, 0)
       end if
       call barrier
    end do

!let all processors have the group0 array
    call broadcast (group0)

    grp0 = group0    

! TT> brought down here from init_job_name in file_utils.fpp
    call scope (subprocs)
! <TT

# endif
  end subroutine init_jobs

  subroutine all_to_group_real (all, group, njobs)
    
    implicit none

    real, dimension (:), intent (in) :: all
    real, intent (out) :: group
    integer, intent (in) :: njobs

    integer :: ik, tag, idx

    tag = 1000

# ifndef MPI
    call error("all_to_group")
# else
!    do ik = 0, njobs-1
!       if (proc0) then
!          if (iproc == grp0(ik)) then
!             group = all(ik+1)
!          else
!             call ssend (all(ik+1), grp0(ik), tag)
!          end if
!       else if (iproc == grp0(ik)) then
!          call receive (group, 0, tag)
!       end if
!!       call barrier
!    end do
    do ik = 0, njobs-1
       if (proc0) then
          idx = mod(ik,size(all))
          if (iproc == grp0(ik)) then
             group = all(idx+1)
          else
             call ssend (all(idx+1), grp0(ik), tag)
          end if
       else if (iproc == grp0(ik)) then
          call receive (group, 0, tag)
       end if
!       call barrier
    end do

# endif
  end subroutine all_to_group_real

  subroutine all_to_group_real_array (all, group, njobs)
    
    implicit none

    real, dimension (:,:), intent (in) :: all
    real, dimension (:), intent (out) :: group
    integer, intent (in) :: njobs

    integer :: ik, tag, idx

# ifndef MPI
    call error ("all_to_group")
# else
    tag = 1001

!    do ik = 0, njobs-1
!       if (proc0) then
!          if (iproc == grp0(ik)) then
!             group = all(ik+1,:)
!          else
!             call ssend (all(ik+1,:), grp0(ik), tag)
!          end if
!       else if (iproc == grp0(ik)) then
!          call receive (group, 0, tag)
!       end if
!!       call barrier
!    end do
    do ik = 0, njobs-1
       if (proc0) then
          idx = mod(ik,size(all,dim=1))
          if (iproc == grp0(ik)) then
             group = all(idx+1,:)
          else
             call ssend (all(idx+1,:), grp0(ik), tag)
          end if
       else if (iproc == grp0(ik)) then
          call receive (group, 0, tag)
       end if
!       call barrier
    end do

# endif
  end subroutine all_to_group_real_array

  subroutine group_to_all_real (group, all, njobs)
    
    implicit none

    real, intent (in) :: group
    real, dimension (:), intent (out) :: all
    integer, intent (in) :: njobs

    integer :: ik, tag, idx

    tag = 1002

# ifndef MPI
    call error("group_to_all")
# else
    do ik = 0, njobs-1
       if (iproc == grp0(ik)) then
          if (.not. proc0) then
             call ssend (group, 0, tag)
          else
             idx = mod(ik,size(all))
             all(idx+1) = group
          end if
       else if (proc0) then
          idx = mod(ik,size(all))
          call receive (all(idx+1), grp0(ik), tag)
       end if
!       call barrier
    end do

# endif    
  end subroutine group_to_all_real

  subroutine group_to_all_real_array (group, all, njobs)
    
    implicit none

    real, dimension (:), intent (in) :: group
    real, dimension (:,:), intent (out) :: all
    integer, intent (in) :: njobs

    integer :: ik, tag, idx

    tag = 1003

# ifndef MPI
    call error("group_to_all")
# else
    do ik = 0, njobs-1
       if (iproc == grp0(ik)) then
          if (.not. proc0) then
             call ssend (group, 0, tag)
          else
             idx = mod(ik,size(all))
             all(idx+1,:) = group
          end if
       else if (proc0) then
          idx = mod(ik,size(all))
          call receive (all(idx+1,:), grp0(ik), tag)
       end if
!       call barrier
    end do
    
# endif
  end subroutine group_to_all_real_array

  subroutine mp_abort (msg)
    use file_utils, only: error_unit, flush_output_file
    implicit none
    character(len=*), intent (in) :: msg
# ifdef MPI
    integer :: ierror
    integer, parameter :: error_code=MPI_ERR_UNKNOWN
# endif

    if (proc0) then
       write (error_unit(),*) "Error: "//msg
       call flush_output_file (error_unit())
    end if

# ifndef MPI
# ifndef NO_ABORT
    call abort
# endif
# else
    call mpi_abort(comm_all, error_code, ierror)
# endif
  end subroutine mp_abort


# ifndef MPI
  subroutine error (msg)
    use file_utils, only: error_unit
    implicit none
    character(len=*), intent (in) :: msg

    write (error_unit(),*) "mp error: "//msg
# ifndef NO_ABORT
    call abort
# endif
    stop
  end subroutine error
# endif

  ! this gathers a single integer from each processor into an array on proc0
  subroutine mp_gather (senddata, recvarray)

    implicit none

    integer, intent (in) :: senddata
    integer, dimension (:), intent (out) :: recvarray

    integer :: ierr

    call mpi_gather (senddata, 1, mpi_integer, recvarray, &
         1, mpi_integer, 0, mp_comm, ierr)

  end subroutine mp_gather

end module mp
