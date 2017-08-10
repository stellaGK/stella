! Modifications for using FFTW version 3:
! (c) The Numerical Algorithms Group (NAG) Ltd, 2009 
!                                 on behalf of the HECToR project


# include "define.inc"

module fft_work

  use constants, only: kind_id

  implicit none

  public :: fft_type, delete_fft
  public :: init_ccfftw, init_crfftw, init_rcfftw, init_z
  public :: FFTW_FORWARD, FFTW_BACKWARD

  private

  interface init_crfftw
     module procedure init_crfftw_1d
     module procedure init_crfftw_2d
  end interface

  interface init_rcfftw
     module procedure init_rcfftw_1d
     module procedure init_rcfftw_2d
  end interface

  type :: fft_type
! TT>
!     integer :: n, plan, is, type
     integer :: n, is, type
     integer (kind_id) :: plan

! <TT
# if FFT == _FFTW3_
     integer :: howmany
     logical :: strided
# endif     
     real :: scale
  end type fft_type

# if FFT == _FFTW_

  ! parameters defined in fftw_f77.i
  integer, parameter :: fftw_estimate   =  0
  integer, parameter :: fftw_measure    =  1
  integer, parameter :: fftw_in_place   =  8
  integer, parameter :: fftw_use_wisdom = 16
  integer, parameter :: FFTW_FORWARD=-1, FFTW_BACKWARD=1
  integer, parameter :: FFTW_REAL_TO_COMPLEX=FFTW_FORWARD
  integer, parameter :: FFTW_COMPLEX_TO_REAL=FFTW_BACKWARD

# elif FFT == _FFTW3_

  ! read the parameters from the file of installation
#include "fftw3.f"

# else
  integer, parameter :: FFTW_FORWARD=-1, FFTW_BACKWARD=1
#endif

contains

# if FFT == _FFTW3_  
  subroutine init_z (fft, is, n, howmany)
# else
  subroutine init_z (fft, is, n)
# endif

    use mp, only: mp_abort

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
# if FFT == _FFTW3_
    integer, optional, intent (in) :: howmany
# endif
    integer :: j

# if FFT == _FFTW3_
    complex, dimension (:,:), allocatable :: dummy_in_data, dummy_out_data
    integer, dimension (1) :: array_n, embed
# endif    

    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 1
# if FFT == _FFTW3_
    if (present(howmany)) then
       fft%howmany = howmany
       fft%strided = .false.
    else
       call mp_abort("For FFTW3 howmany needs to be present in init_z")
    endif
# endif
    
# if FFT == _FFTW_
    j = fftw_measure + fftw_use_wisdom
    call fftw_f77_create_plan(fft%plan,n,is,j)
# elif FFT == _FFTW3_
    array_n = n
    embed   = n+1
    
    allocate (dummy_in_data(n+1, howmany), dummy_out_data(n+1, howmany))

    j = FFTW_PATIENT + FFTW_UNALIGNED
    call dfftw_plan_many_dft(fft%plan, 1, array_n, howmany, &
         dummy_in_data,  embed, 1, n+1, &
         dummy_out_data, embed, 1, n+1, &
         is, j)

    deallocate (dummy_in_data, dummy_out_data)
# endif

  end subroutine init_z


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# if FFT == _FFTW3_
  subroutine init_ccfftw (fft, is, n, howmany, data_array)
# else
  subroutine init_ccfftw (fft, is, n)
# endif
    use mp, only: mp_abort

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
# if FFT == _FFTW3_
    integer, optional, intent (in) :: howmany
    complex, optional, dimension(:,:), intent(inout) :: data_array 
    integer, dimension(1) :: array_n
# endif    

    integer :: j
    
    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 1
# if FFT == _FFTW3_
    if (present(howmany) .and. present(data_array)) then
       fft%howmany = howmany
       fft%strided = .false.
    else
       call mp_abort ("no howmany or data_array in init_ccfftw: FIX!")
    endif
# endif
    
    
# if FFT == _FFTW_
    j = fftw_in_place + fftw_measure + fftw_use_wisdom
    call fftw_f77_create_plan(fft%plan,n,is,j)
# elif FFT == _FFTW3_
    ! the planer expects this as an array of size 1
    array_n = n
    j = FFTW_PATIENT
    call dfftw_plan_many_dft(fft%plan, 1, array_n, howmany, &
         data_array, array_n, 1, n, &
         data_array, array_n, 1, n, &
         is, j)
# endif

  end subroutine init_ccfftw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# if FFT == _FFTW_

  subroutine init_rcfftw_1d (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    integer :: j

    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 0


    j = fftw_measure + fftw_use_wisdom
    call rfftwnd_f77_create_plan(fft%plan,1,N,is,j)

  end subroutine init_rcfftw_1d
  
# elif FFT == _FFTW3_

  subroutine init_rcfftw_1d (fft, is, n, howmany)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    integer, intent(in) :: howmany
    integer :: j


    ! a few things required for FFTW3 over FFTW2    

    integer, dimension(1) :: vector_sizes_real, vector_sizes_complex
    integer :: stride

    ! we need two dummy arrays for the planner.  Since at 
    ! present we are not using SSE instructions, these are save
    ! to be local to this routine

    real, dimension(:,:), allocatable :: dummy_real_data
    complex, dimension (:,:), allocatable :: dummy_complex_data



    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 0
    fft%howmany = howmany
    fft%strided = .false.

    allocate (dummy_real_data(N, max (1, howmany)))
    allocate (dummy_complex_data(N/2+1, max(1, howmany)))

    vector_sizes_real = N
    vector_sizes_complex = N/2+1
    j = FFTW_PATIENT ! + FFTW_UNALIGNED 
    call dfftw_plan_many_dft_r2c (fft%plan, 1, &
         vector_sizes_real, howmany, &
         dummy_real_data,    vector_sizes_real,    1, N,     &
         dummy_complex_data, vector_sizes_complex, 1, N/2+1, &
         j)

    deallocate (dummy_real_data, dummy_complex_data)


  end subroutine init_rcfftw_1d

!<GGH
#else
  subroutine init_rcfftw_1d (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    integer :: j

    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 0

  end subroutine init_rcfftw_1d

!>GGH

# endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! because of the interface we need different routines for FFTW2 and FFTW3

# if FFT == _FFTW_

  subroutine init_crfftw_1d (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    integer :: j


    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 0


    j = fftw_measure + fftw_use_wisdom
    call rfftwnd_f77_create_plan(fft%plan,1,N,is,j)


  end subroutine init_crfftw_1d


# elif FFT == _FFTW3_


  subroutine init_crfftw_1d (fft, is, n, howmany)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    integer, intent(in) :: howmany
    integer :: j

    ! a few things required for FFTW3 over FFTW2    

    integer, dimension(1) :: vector_sizes_real, vector_sizes_complex
    integer :: stride

    ! we need two dummy arrays for the planner.  Since at 
    ! present we are not using SSE instructions, these are save
    ! to be local to this routine

    real, dimension(:,:), allocatable :: dummy_real_data
    complex, dimension (:,:), allocatable :: dummy_complex_data


    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 0
    fft%howmany = howmany
    fft%strided = .false.

    allocate (dummy_real_data(N, max (1, howmany)))
    allocate (dummy_complex_data(N/2+1, max(1, howmany)))

    vector_sizes_real = N
    vector_sizes_complex = N/2+1
    j = FFTW_PATIENT ! + FFTW_UNALIGNED 
    call dfftw_plan_many_dft_c2r (fft%plan, 1, &
         vector_sizes_real, howmany, &
         dummy_complex_data, vector_sizes_complex, 1, N/2+1, &
         dummy_real_data,    vector_sizes_real,    1, N,     &
         j)

    deallocate (dummy_real_data, dummy_complex_data)

  end subroutine init_crfftw_1d

!<GGH
# else

  subroutine init_crfftw_1d (fft, is, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    integer :: j


    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 0

  end subroutine init_crfftw_1d
!>GGH

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# if FFT == _FFTW_

  subroutine init_rcfftw_2d (fft, is, m, n)
       
    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, m, n
    integer :: j

    fft%scale = 1./real(m*n)
    if (is > 0) fft%scale = 1.

    j = fftw_measure + fftw_use_wisdom
    call rfftw2d_f77_create_plan(fft%plan,m,n,is,j)

  end subroutine init_rcfftw_2d


# elif FFT == _FFTW3_


    subroutine init_rcfftw_2d (fft, is, m, n, howmany)
       
    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, m, n
    integer, intent (in) :: howmany
    integer :: j

! a few things required for FFTW3 over FFTW2    

    integer, parameter :: fft_dimension = 2
    integer  :: stride

    ! to pass to FFTW3, we need vector for actual sizes of the real 
    ! and complex data

    integer, dimension(fft_dimension) :: &
         vector_sizes_real, vector_sizes_complex 


    ! we need two dummy arrays for the planner.  Since at 
    ! present we are not using SSE instructions, these are save
    ! to be local to this routine

    real, dimension(:,:), allocatable :: dummy_real_data
    complex, dimension (:,:), allocatable :: dummy_complex_data



    fft%scale = 1./real(m*n)
    if (is > 0) fft%scale = 1.
    fft%howmany = howmany
    fft%strided = .true.

    stride  = howmany

    allocate (dummy_real_data(howmany, m*n))
    allocate (dummy_complex_data(howmany, (m/2+1)*n))

    vector_sizes_real(1) = m
    vector_sizes_real(2) = n

    vector_sizes_complex(1) = m/2+1
    vector_sizes_complex(2) = n

    j = FFTW_PATIENT + FFTW_UNALIGNED

    call dfftw_plan_many_dft_r2c(fft%plan, fft_dimension, &
         vector_sizes_real, howmany, &
         dummy_real_data, vector_sizes_real, stride, 1, &
         dummy_complex_data, vector_sizes_complex, stride, 1, &
         j)
    
    deallocate (dummy_real_data, dummy_complex_data)


  end subroutine init_rcfftw_2d
!<GGH
# else
  subroutine init_rcfftw_2d (fft, is, m, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, m, n
    integer :: j

    fft%scale = 1./real(m*n)
    if (is > 0) fft%scale = 1.

  end subroutine init_rcfftw_2d
!>GGH


# endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# if FFT == _FFTW_

    subroutine init_crfftw_2d (fft, is, m, n)
    
    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, m, n
    integer :: j

    fft%scale = 1./real(m*n)
    if (is > 0) fft%scale = 1.

    j = fftw_measure + fftw_use_wisdom
    call rfftw2d_f77_create_plan(fft%plan,m,n,is,j)

  end subroutine init_crfftw_2d


# elif FFT == _FFTW3_

  subroutine init_crfftw_2d (fft, is, m, n, howmany)

    
    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, m, n
    integer, intent (in) :: howmany
    integer :: j

! a few things required for FFTW3 over FFTW2    

    integer, parameter :: fft_dimension = 2
    integer  :: stride


    ! to pass to FFTW3, we need vector for actual sizes of the real 
    ! and complex data

    integer, dimension(fft_dimension) :: &
         vector_sizes_real, vector_sizes_complex 


    ! we need two dummy arrays for the planner.  Since at 
    ! present we are not using SSE instructions, these are save
    ! to be local to this routine

    real, dimension(:,:), allocatable :: dummy_real_data
    complex, dimension (:,:), allocatable :: dummy_complex_data



    fft%scale = 1./real(m*n)
    if (is > 0) fft%scale = 1.
    fft%howmany = howmany
    fft%strided = .true.

    stride  = howmany

    allocate (dummy_real_data(howmany, m*n))
    allocate (dummy_complex_data(howmany, (m/2+1)*n))

    vector_sizes_real(1) = m
    vector_sizes_real(2) = n

    vector_sizes_complex(1) = m/2+1
    vector_sizes_complex(2) = n

    j = FFTW_PATIENT + FFTW_UNALIGNED

    call dfftw_plan_many_dft_c2r(fft%plan, fft_dimension, &
         vector_sizes_real, howmany, &
         dummy_complex_data, vector_sizes_complex, stride, 1, &
         dummy_real_data, vector_sizes_real, stride, 1, &
         j)
    
    deallocate (dummy_real_data, dummy_complex_data)

  end subroutine init_crfftw_2d

!<GGH
# else

  subroutine init_crfftw_2d (fft, is, m, n)

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, m, n
    integer :: j

    fft%scale = 1./real(m*n)
    if (is > 0) fft%scale = 1.

  end subroutine init_crfftw_2d

!>GGH

# endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine delete_fft(fft)
    
    type (fft_type), intent (in out) :: fft

# if FFT == _FFTW_
    if (fft%type == 1) then
       call fftw_f77_destroy_plan(fft%plan)
    else
       call rfftw_f77_destroy_plan(fft%plan)
    end if
# elif FFT == _FFTW3_
       call dfftw_destroy_plan(fft%plan)
# endif

  end subroutine delete_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module fft_work
