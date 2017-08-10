# include "define.inc"

module hdf_wrapper

  !
  ! hdf5 wrapper subroutines written by Tomo Tatsuno
  ! once ver. 1.8.0 becomes available, we can remove this module
  ! and use those High-level APIs.
  !
  ! I haven't written it carefully since it will go away soon.
  ! Honestly there must be more checks such as size mismatch
  ! between file and memory arrays, skipping after error etc...
  !

# ifdef HDF
# if FCOMPILER == _GFORTRAN_
  use hdf5
# else
  use hdf5, only: HID_T, HSIZE_T, H5T_NATIVE_INTEGER, H5S_SCALAR_F, SIZE_T
  use hdf5, only: h5screate_f, h5screate_simple_f, h5dcreate_f
  use hdf5, only: h5dwrite_f, h5dclose_f, h5sclose_f
  use hdf5, only: h5dopen_f, h5dread_f
# endif
# endif

  implicit none

  public :: hdf_init, hdf_finish
# ifdef HDF
  public :: hdf_write, hdf_read, hdf_error
  public :: hdf_file_real, hdf_mem_real

  private

  interface hdf_write
     module procedure hdf_write_string, hdf_write_logical_scalar
     module procedure hdf_write_integer_scalar, hdf_write_integer_rank1
     module procedure hdf_write_real_scalar, hdf_write_real_rank1
     module procedure hdf_write_real_rank2, hdf_write_real_rank3
     module procedure hdf_write_real_rank4, hdf_write_real_rank5
     module procedure hdf_write_complex_scalar, hdf_write_complex_rank1
     module procedure hdf_write_complex_rank2, hdf_write_complex_rank3
     module procedure hdf_write_complex_rank4
  end interface

  interface hdf_read
     module procedure hdf_read_string, hdf_read_logical_scalar
     module procedure hdf_read_integer_scalar
     module procedure hdf_read_real_scalar, hdf_read_real_rank1
     module procedure hdf_read_real_rank2, hdf_read_real_rank3
     module procedure hdf_read_real_rank4, hdf_read_real_rank5
     module procedure hdf_read_complex_scalar, hdf_read_complex_rank1
     module procedure hdf_read_complex_rank2, hdf_read_complex_rank3
     module procedure hdf_read_complex_rank4
  end interface

  logical, save :: stop_private
  integer (HID_T), save :: hdf_file_real, hdf_mem_real
  logical, save :: initialized = .false.
# endif
contains

  subroutine hdf_init (stop, dbl)

    use constants, only: kind_rs, kind_rd, pi
# if (defined HDF && FCOMPILER != _GFORTRAN_)
    use hdf5, only: H5T_NATIVE_REAL, H5T_NATIVE_DOUBLE, h5open_f
# endif
!    use mp, only: proc0

    logical, intent (in), optional :: stop
    logical, intent (in), optional :: dbl
# ifdef HDF
    integer :: ier
    call h5open_f (ier)

    stop_private = .false.
    if (present(stop)) stop_private = stop

    ! find and set code precision
    if ( (kind(pi)==kind_rs) .or. (kind_rs==kind_rd) ) then
       hdf_mem_real = H5T_NATIVE_REAL
!       if (proc0) print *, '# hdf_init: code is single precision'
    else if (kind(pi)==kind_rd) then
       hdf_mem_real = H5T_NATIVE_DOUBLE
!       if (proc0) print *, '# hdf_init: code is double precision'
    else
       write (*,*) 'ERROR: precision mismatch in hdf_init'
    end if

    ! set file precision
    hdf_file_real = H5T_NATIVE_REAL
    if (present(dbl)) then
       if (dbl .and. hdf_mem_real == H5T_NATIVE_DOUBLE) &
            hdf_file_real = H5T_NATIVE_DOUBLE
    end if

    initialized = .true.
# endif
  end subroutine hdf_init

  subroutine hdf_finish
# ifdef HDF
# if FCOMPILER != _GFORTRAN_
    use hdf5, only: h5close_f
# endif
    integer :: ier
    ! finalize if initialized
    if (.not. initialized) return
    call h5close_f (ier)
    initialized = .false.
# endif
  end subroutine hdf_finish

# ifdef HDF
  !
  ! write subroutines
  !
  subroutine hdf_write_string (loc, name, data, ier)

# if FCOMPILER != _GFORTRAN_
    use hdf5, only: H5T_NATIVE_CHARACTER, h5tcopy_f, h5tset_size_f, h5tclose_f
# endif
    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    character (*), intent (in) :: data
    integer, intent (inout) :: ier
    integer (HID_T) :: dtp, dsp, dst
    integer (SIZE_T) :: length
    integer (HSIZE_T), dimension (1) :: dim_scalar = (/1/)
    integer :: ierr

    dtp=0 ; dsp=0 ; dst=0
    call h5tcopy_f (H5T_NATIVE_CHARACTER, dtp, ier)
    if (ier < 0) call hdf_error (dset=name, message='copy dataset type')
    if (ier >= 0) then
       length = len(data)
       call h5tset_size_f (dtp, length, ier)
       if (ier < 0) call hdf_error (dset=name, message='set dataset size')
    end if
    if (ier >= 0) then
       call h5screate_f (H5S_SCALAR_F, dsp, ier)
       if (ier < 0) call hdf_error (dset=name, message='create dataspace')
    end if
    if (ier >= 0) then
       call h5dcreate_f (loc, name, dtp, dsp, dst, ier)
       if (ier < 0) call hdf_error (dset=name, message='create dataset')
    end if
    if (ier >= 0) then
       call h5dwrite_f (dst, dtp, data, dim_scalar, ier)
       if (ier < 0) call hdf_error (dset=name, message='write')
    end if
    if (dst /= 0) then
       call h5dclose_f (dst, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataset')
          if (ier >= 0) ier = ierr
       end if
    end if
    if (dsp /= 0) then
       call h5sclose_f (dsp, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataspace')
          if (ier >= 0) ier = ierr
       end if
    end if
    if (dtp /= 0) then
       call h5tclose_f (dtp, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataset type')
          if (ier >= 0) ier = ierr
       end if
    end if

  end subroutine hdf_write_string

  subroutine hdf_write_logical_scalar (loc, name, data, ier)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    logical, intent (in) :: data
    integer, intent (inout) :: ier
    character :: char

    write (char, '(l1)') data
    call hdf_write_string (loc, name, char, ier)

  end subroutine hdf_write_logical_scalar

  subroutine hdf_write_integer_scalar (loc, name, data, ier)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    integer, intent (in) :: data
    integer, intent (inout) :: ier
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (1) :: dim_scalar = (/1/)
    integer :: ierr

    dsp=0 ; dst=0
    call h5screate_f (H5S_SCALAR_F, dsp, ier)
    if (ier < 0) call hdf_error (dset=name, message='create dataspace')
    if (ier >= 0) then
       call h5dcreate_f (loc, name, H5T_NATIVE_INTEGER, dsp, dst, ier)
       if (ier < 0) call hdf_error (dset=name, message='create dataset')
    end if
    if (ier >= 0) then
       call h5dwrite_f (dst, H5T_NATIVE_INTEGER, data, dim_scalar, ier)
       if (ier < 0) call hdf_error (dset=name, message='write')
    end if
    if (dst /= 0) then
       call h5dclose_f (dst, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataset')
          if (ier >= 0) ier = ierr
       end if
    end if
    if (dsp /= 0) then
       call h5sclose_f (dsp, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataspace')
          if (ier >= 0) ier = ierr
       end if
    end if

  end subroutine hdf_write_integer_scalar

  subroutine hdf_write_integer_rank1 (loc, name, data, ier)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    integer, dimension (:), intent (in) :: data
    integer, intent (inout) :: ier
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (:), allocatable :: dims
    integer :: ierr, rank

    dsp=0 ; dst=0
    rank = size(shape(data))
    allocate (dims(rank))
    dims = shape(data)
    call h5screate_simple_f (rank, dims, dsp, ier)
    if (ier < 0) call hdf_error (dset=name, message='create dataspace')
    if (ier >= 0) then
       call h5dcreate_f (loc, name, H5T_NATIVE_INTEGER, dsp, dst, ier)
       if (ier < 0) call hdf_error (dset=name, message='create dataset')
    end if
    if (ier >= 0) then
       call h5dwrite_f (dst, H5T_NATIVE_INTEGER, data, dims, ier)
       if (ier < 0) call hdf_error (dset=name, message='write')
    end if
    if (dst /= 0) then
       call h5dclose_f (dst, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataset')
          if (ier >= 0) ier = ierr
       end if
    end if
    if (dsp /= 0) then
       call h5sclose_f (dsp, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataspace')
          if (ier >= 0) ier = ierr
       end if
    end if
    deallocate (dims)

  end subroutine hdf_write_integer_rank1

  subroutine hdf_write_real_scalar (loc, name, data, ier)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    real, intent (in) :: data
    integer, intent (inout) :: ier
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (1) :: dim_scalar = (/1/)
    integer :: ierr

    dsp=0 ; dst=0
    call h5screate_f (H5S_SCALAR_F, dsp, ier)
    if (ier < 0) call hdf_error (dset=name, message='create dataspace')
    if (ier >= 0) then
       call h5dcreate_f (loc, name, hdf_file_real, dsp, dst, ier)
       if (ier < 0) call hdf_error (dset=name, message='create dataset')
    end if
    if (ier >= 0) then
       call h5dwrite_f (dst, hdf_mem_real, data, dim_scalar, ier)
       if (ier < 0) call hdf_error (dset=name, message='write')
    end if
    if (dst /= 0) then
       call h5dclose_f (dst, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataset')
          if (ier >= 0) ier = ierr
       end if
    end if
    if (dsp /= 0) then
       call h5sclose_f (dsp, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataspace')
          if (ierr >= 0) ier = ierr
       end if
    end if

  end subroutine hdf_write_real_scalar

  subroutine hdf_write_real_rank1 (loc, name, data, ier)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    real, dimension (:), intent (in) :: data
    integer, intent (inout) :: ier
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (:), allocatable :: dims
    integer :: ierr, rank

    dsp=0 ; dst=0
    rank = size(shape(data))
    allocate (dims(rank))
    dims = shape(data)
    call h5screate_simple_f (rank, dims, dsp, ier)
    if (ier < 0) call hdf_error (dset=name, message='create dataspace')
    if (ier >= 0) then
       call h5dcreate_f (loc, name, hdf_file_real, dsp, dst, ier)
       if (ier < 0) call hdf_error (dset=name, message='create dataset')
    end if
    if (ier >= 0) then
       call h5dwrite_f (dst, hdf_mem_real, data, dims, ier)
       if (ier < 0) call hdf_error (dset=name, message='write')
    end if
    if (dst /= 0) then
       call h5dclose_f (dst, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataset')
          if (ier >= 0) ier = ierr
       end if
    end if
    if (dsp /= 0) then
       call h5sclose_f (dsp, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataspace')
          if (ierr >= 0) ier = ierr
       end if
    end if

  end subroutine hdf_write_real_rank1

  subroutine hdf_write_real_rank2 (loc, name, data, ier)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    real, dimension (:,:), intent (in) :: data
    integer, intent (inout) :: ier
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (:), allocatable :: dims
    integer :: ierr, rank

    dsp=0 ; dst=0
    rank = size(shape(data))
    allocate (dims(rank))
    dims = shape(data)
    call h5screate_simple_f (rank, dims, dsp, ier)
    if (ier < 0) call hdf_error (dset=name, message='create dataspace')
    if (ier >= 0) then
       call h5dcreate_f (loc, name, hdf_file_real, dsp, dst, ier)
       if (ier < 0) call hdf_error (dset=name, message='create dataset')
    end if
    if (ier >= 0) then
       call h5dwrite_f (dst, hdf_mem_real, data, dims, ier)
       if (ier < 0) call hdf_error (dset=name, message='write')
    end if
    if (dst /= 0) then
       call h5dclose_f (dst, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataset')
          if (ier >= 0) ier = ierr
       end if
    end if
    if (dsp /= 0) then
       call h5sclose_f (dsp, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataspace')
          if (ierr >= 0) ier = ierr
       end if
    end if

  end subroutine hdf_write_real_rank2

  subroutine hdf_write_real_rank3 (loc, name, data, ier)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    real, dimension (:,:,:), intent (in) :: data
    integer, intent (inout) :: ier
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (:), allocatable :: dims
    integer :: ierr, rank

    dsp=0 ; dst=0
    rank = size(shape(data))
    allocate (dims(rank))
    dims = shape(data)
    call h5screate_simple_f (rank, dims, dsp, ier)
    if (ier < 0) call hdf_error (dset=name, message='create dataspace')
    if (ier >= 0) then
       call h5dcreate_f (loc, name, hdf_file_real, dsp, dst, ier)
       if (ier < 0) call hdf_error (dset=name, message='create dataset')
    end if
    if (ier >= 0) then
       call h5dwrite_f (dst, hdf_mem_real, data, dims, ier)
       if (ier < 0) call hdf_error (dset=name, message='write')
    end if
    if (dst /= 0) then
       call h5dclose_f (dst, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataset')
          if (ier >= 0) ier = ierr
       end if
    end if
    if (dsp /= 0) then
       call h5sclose_f (dsp, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataspace')
          if (ierr >= 0) ier = ierr
       end if
    end if

  end subroutine hdf_write_real_rank3

  subroutine hdf_write_real_rank4 (loc, name, data, ier)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    real, dimension (:,:,:,:), intent (in) :: data
    integer, intent (inout) :: ier
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (:), allocatable :: dims
    integer :: ierr, rank

    dsp=0 ; dst=0
    rank = size(shape(data))
    allocate (dims(rank))
    dims = shape(data)
    call h5screate_simple_f (rank, dims, dsp, ier)
    if (ier < 0) call hdf_error (dset=name, message='create dataspace')
    if (ier >= 0) then
       call h5dcreate_f (loc, name, hdf_file_real, dsp, dst, ier)
       if (ier < 0) call hdf_error (dset=name, message='create dataset')
    end if
    if (ier >= 0) then
       call h5dwrite_f (dst, hdf_mem_real, data, dims, ier)
       if (ier < 0) call hdf_error (dset=name, message='write')
    end if
    if (dst /= 0) then
       call h5dclose_f (dst, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataset')
          if (ier >= 0) ier = ierr
       end if
    end if
    if (dsp /= 0) then
       call h5sclose_f (dsp, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataspace')
          if (ierr >= 0) ier = ierr
       end if
    end if

  end subroutine hdf_write_real_rank4

  subroutine hdf_write_real_rank5 (loc, name, data, ier)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    real, dimension (:,:,:,:,:), intent (in) :: data
    integer, intent (inout) :: ier
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (:), allocatable :: dims
    integer :: ierr, rank

    dsp=0 ; dst=0
    rank = size(shape(data))
    allocate (dims(rank))
    dims = shape(data)
    call h5screate_simple_f (rank, dims, dsp, ier)
    if (ier < 0) call hdf_error (dset=name, message='create dataspace')
    if (ier >= 0) then
       call h5dcreate_f (loc, name, hdf_file_real, dsp, dst, ier)
       if (ier < 0) call hdf_error (dset=name, message='create dataset')
    end if
    if (ier >= 0) then
       call h5dwrite_f (dst, hdf_mem_real, data, dims, ier)
       if (ier < 0) call hdf_error (dset=name, message='write')
    end if
    if (dst /= 0) then
       call h5dclose_f (dst, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataset')
          if (ier >= 0) ier = ierr
       end if
    end if
    if (dsp /= 0) then
       call h5sclose_f (dsp, ierr)
       if (ierr < 0) then
          call hdf_error (dset=name, message='close dataspace')
          if (ierr >= 0) ier = ierr
       end if
    end if

  end subroutine hdf_write_real_rank5

  subroutine hdf_write_complex_scalar (loc, name, data, ier)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    complex, intent (in) :: data
    integer, intent (inout) :: ier
    real, dimension (2) :: rtmp

    rtmp = (/ real(data), aimag(data) /)
    call hdf_write_real_rank1 (loc, name, rtmp, ier)

  end subroutine hdf_write_complex_scalar

  subroutine hdf_write_complex_rank1 (loc, name, data, ier)

    use convert, only: c2r
    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    complex, dimension (:), intent (in) :: data
    integer, intent (inout) :: ier
    real, dimension (:,:), allocatable :: rtmp

    allocate (rtmp(2,size(data)))
    call c2r (data, rtmp)
    call hdf_write_real_rank2 (loc, name, rtmp, ier)
    deallocate (rtmp)

  end subroutine hdf_write_complex_rank1

  subroutine hdf_write_complex_rank2 (loc, name, data, ier)

    use convert, only: c2r
    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    complex, dimension (:,:), intent (in) :: data
    integer, intent (inout) :: ier
    real, dimension (:,:,:), allocatable :: rtmp

    allocate (rtmp(2,size(data,1),size(data,2)))
    call c2r (data, rtmp)
    call hdf_write_real_rank3 (loc, name, rtmp, ier)
    deallocate (rtmp)

  end subroutine hdf_write_complex_rank2

  subroutine hdf_write_complex_rank3 (loc, name, data, ier)

    use convert, only: c2r
    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    complex, dimension (:,:,:), intent (in) :: data
    integer, intent (inout) :: ier
    real, dimension (:,:,:,:), allocatable :: rtmp

    allocate (rtmp(2,size(data,1),size(data,2),size(data,3)))
    call c2r (data, rtmp)
    call hdf_write_real_rank4 (loc, name, rtmp, ier)
    deallocate (rtmp)

  end subroutine hdf_write_complex_rank3

  subroutine hdf_write_complex_rank4 (loc, name, data, ier)

    use convert, only: c2r
    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    complex, dimension (:,:,:,:), intent (in) :: data
    integer, intent (inout) :: ier
    real, dimension (:,:,:,:,:), allocatable :: rtmp

    allocate (rtmp(2,size(data,1),size(data,2),size(data,3),size(data,4)))
    call c2r (data, rtmp)
    call hdf_write_real_rank5 (loc, name, rtmp, ier)
    deallocate (rtmp)

  end subroutine hdf_write_complex_rank4

  !
  ! read subroutines
  !
  subroutine hdf_read_string (loc, name, data)

# if FCOMPILER != _GFORTRAN_
    use hdf5, only: h5dget_type_f
# endif
    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    character (*), intent (out) :: data
    integer (HID_T) :: dsp, dst, dtp
    integer (HSIZE_T), dimension (1) :: dim_scalar = (/1/)
    integer :: ier

    data = ''
    call h5dopen_f (loc, name, dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dget_type_f (dst, dtp, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dread_f (dst, dtp, data, dim_scalar, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dclose_f (dst, ier)
    if (ier < 0) call hdf_error (dset=name)

  end subroutine hdf_read_string

  subroutine hdf_read_logical_scalar (loc, name, data)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    logical, intent (out) :: data
    character :: char

    call hdf_read_string (loc, name, char)
    read (char, '(l1)') data

  end subroutine hdf_read_logical_scalar

  subroutine hdf_read_integer_scalar (loc, name, data)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    integer, intent (out) :: data
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (1) :: dim_scalar = (/1/)
    integer :: ier

    call h5dopen_f (loc, name, dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dread_f (dst, H5T_NATIVE_INTEGER, data, dim_scalar, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dclose_f (dst, ier)
    if (ier < 0) call hdf_error (dset=name)

  end subroutine hdf_read_integer_scalar

  subroutine hdf_read_real_scalar (loc, name, data)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    real, intent (out) :: data
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (1) :: dim_scalar = (/1/)
    integer :: ier

    call h5dopen_f (loc, name, dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dread_f (dst, hdf_mem_real, data, dim_scalar, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dclose_f (dst, ier)
    if (ier < 0) call hdf_error (dset=name)

  end subroutine hdf_read_real_scalar

  subroutine hdf_read_real_rank1 (loc, name, data)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    real, dimension (:), intent (out) :: data
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (:), allocatable :: dims
    integer :: ier, rank

    rank = size(shape(data))
    allocate (dims(rank))
    dims = shape(data)
    call h5dopen_f (loc, name, dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dread_f (dst, hdf_mem_real, data, dims, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dclose_f (dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    deallocate (dims)

  end subroutine hdf_read_real_rank1

  subroutine hdf_read_real_rank2 (loc, name, data)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    real, dimension (:,:), intent (out) :: data
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (:), allocatable :: dims
    integer :: ier, rank

    rank = size(shape(data))
    allocate (dims(rank))
    dims = shape(data)
    call h5dopen_f (loc, name, dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dread_f (dst, hdf_mem_real, data, dims, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dclose_f (dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    deallocate (dims)

  end subroutine hdf_read_real_rank2

  subroutine hdf_read_real_rank3 (loc, name, data)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    real, dimension (:,:,:), intent (out) :: data
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (:), allocatable :: dims
    integer :: ier, rank

    rank = size(shape(data))
    allocate (dims(rank))
    dims = shape(data)
    call h5dopen_f (loc, name, dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dread_f (dst, hdf_mem_real, data, dims, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dclose_f (dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    deallocate (dims)

  end subroutine hdf_read_real_rank3

  subroutine hdf_read_real_rank4 (loc, name, data)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    real, dimension (:,:,:,:), intent (out) :: data
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (:), allocatable :: dims
    integer :: ier, rank

    rank = size(shape(data))
    allocate (dims(rank))
    dims = shape(data)
    call h5dopen_f (loc, name, dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dread_f (dst, hdf_mem_real, data, dims, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dclose_f (dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    deallocate (dims)

  end subroutine hdf_read_real_rank4

  subroutine hdf_read_real_rank5 (loc, name, data)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    real, dimension (:,:,:,:,:), intent (out) :: data
    integer (HID_T) :: dsp, dst
    integer (HSIZE_T), dimension (:), allocatable :: dims
    integer :: ier, rank

    rank = size(shape(data))
    allocate (dims(rank))
    dims = shape(data)
    call h5dopen_f (loc, name, dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dread_f (dst, hdf_mem_real, data, dims, ier)
    if (ier < 0) call hdf_error (dset=name)
    call h5dclose_f (dst, ier)
    if (ier < 0) call hdf_error (dset=name)
    deallocate (dims)

  end subroutine hdf_read_real_rank5

  subroutine hdf_read_complex_scalar (loc, name, data)

    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    complex, intent (out) :: data
    real, dimension (2) :: rtmp

    call hdf_read_real_rank1 (loc, name, rtmp)
    data = cmplx(rtmp(1),rtmp(2))

  end subroutine hdf_read_complex_scalar

  subroutine hdf_read_complex_rank1 (loc, name, data)

    use convert, only: r2c
    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    complex, dimension (:), intent (out) :: data
    real, dimension (:,:), allocatable :: rtmp

    allocate (rtmp(2,size(data)))
    call hdf_read_real_rank2 (loc, name, rtmp)
    call r2c (data, rtmp)

  end subroutine hdf_read_complex_rank1

  subroutine hdf_read_complex_rank2 (loc, name, data)

    use convert, only: r2c
    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    complex, dimension (:,:), intent (out) :: data
    real, dimension (:,:,:), allocatable :: rtmp

    allocate (rtmp(2,size(data,1),size(data,2)))
    call hdf_read_real_rank3 (loc, name, rtmp)
    call r2c (data, rtmp)

  end subroutine hdf_read_complex_rank2

  subroutine hdf_read_complex_rank3 (loc, name, data)

    use convert, only: r2c
    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    complex, dimension (:,:,:), intent (out) :: data
    real, dimension (:,:,:,:), allocatable :: rtmp

    allocate (rtmp(2,size(data,1),size(data,2),size(data,3)))
    call hdf_read_real_rank4 (loc, name, rtmp)
    call r2c (data, rtmp)

  end subroutine hdf_read_complex_rank3

  subroutine hdf_read_complex_rank4 (loc, name, data)

    use convert, only: r2c
    integer (HID_T), intent (in) :: loc
    character (*), intent (in) :: name
    complex, dimension (:,:,:,:), intent (out) :: data
    real, dimension (:,:,:,:,:), allocatable :: rtmp

    allocate (rtmp(2,size(data,1),size(data,2),size(data,3),size(data,4)))
    call hdf_read_real_rank5 (loc, name, rtmp)
    call r2c (data, rtmp)

  end subroutine hdf_read_complex_rank4

  !
  ! error output
  !
  subroutine hdf_error (file, grp, dset, message, iproc)

    use file_utils, only: error_unit
    character (*), intent (in), optional :: file
    character (*), intent (in), optional :: grp
    character (*), intent (in), optional :: dset
    character (*), intent (in), optional :: message
    integer, intent (in), optional :: iproc
    integer :: uer

    uer = error_unit()
!    call h5eprint_f (uer) ! file name has to be specified rather than unit no
    if (present(file)) write (uer,*) 'HDF error in file: ', trim (file)
    if (present(grp))  write (uer,*) 'HDF error in group: ', trim (grp)
    if (present(dset)) write (uer,*) 'HDF error in dataset: ', trim (dset)
    if (present(message)) write (uer,*) trim(message), ' failure'
    if (present(iproc)) write (uer,*) 'iproc: ', iproc

    if (stop_private) stop

  end subroutine hdf_error

# endif

end module hdf_wrapper
