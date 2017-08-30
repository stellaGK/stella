! Modifications for using FFTW version 3:
! (c) The Numerical Algorithms Group (NAG) Ltd, 2009 
!                                 on behalf of the HECToR project


# include "define.inc"

module fft_work

  use constants, only: kind_id

  implicit none

  public :: fft_type, delete_fft
  public :: init_ccfftw, init_crfftw, init_rcfftw
  public :: FFTW_FORWARD, FFTW_BACKWARD

  private

  type :: fft_type
     integer :: n, is, type
     integer (kind_id) :: plan
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# if FFT == _FFTW3_
  subroutine init_ccfftw (fft, is, n, howmany, data_in, data_out)
# else
  subroutine init_ccfftw (fft, is, n)
# endif
    use mp, only: mp_abort

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
# if FFT == _FFTW3_
    integer, optional, intent (in) :: howmany
    complex, dimension(:), intent(inout) :: data_in, data_out
    integer, dimension(1) :: array_n
# endif

    integer :: j
    
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
       call mp_abort ("no howmany in init_ccfftw: FIX!")
    endif
# endif
    
# if FFT == _FFTW_
    j = fftw_in_place + fftw_measure + fftw_use_wisdom
    call fftw_f77_create_plan(fft%plan,n,is,j)
# elif FFT == _FFTW3_
    ! the planner expects this as an array of size 1
    array_n = n
    j = FFTW_PATIENT
    call dfftw_plan_dft_1d(fft%plan, n, data_in, data_out, is, j)
!    call dfftw_plan_many_dft(fft%plan, 1, array_n, howmany, &
!         data_array, array_n, 1, n, &
!         data_array, array_n, 1, n, &
!         is, j)
# endif

  end subroutine init_ccfftw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# if FFT == _FFTW3_
  subroutine init_crfftw (fft, is, n, howmany, data_in, data_out)
# else
  subroutine init_crfftw (fft, is, n)
# endif
    use mp, only: mp_abort

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
# if FFT == _FFTW3_
    integer, optional, intent (in) :: howmany
    complex, dimension(:), intent(in out) :: data_in
    real, dimension (:), intent (in out) :: data_out
    integer, dimension(1) :: array_n
# endif

    integer :: j
    
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
       call mp_abort ("no howmany in init_ccfftw: FIX!")
    endif
# endif
    
# if FFT == _FFTW_
    j = fftw_in_place + fftw_measure + fftw_use_wisdom
    call fftw_f77_create_plan(fft%plan,n,is,j)
# elif FFT == _FFTW3_
    ! the planner expects this as an array of size 1
    array_n = n
    j = FFTW_PATIENT
    call dfftw_plan_dft_c2r_1d(fft%plan, n, data_in, data_out, j)
# endif

  end subroutine init_crfftw

# if FFT == _FFTW3_
  subroutine init_rcfftw (fft, is, n, howmany, data_in, data_out)
# else
  subroutine init_rcfftw (fft, is, n)
# endif
    use mp, only: mp_abort

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
# if FFT == _FFTW3_
    integer, optional, intent (in) :: howmany
    real, dimension(:), intent(in out) :: data_in
    complex, dimension (:), intent (in out) :: data_out
    integer, dimension(1) :: array_n
# endif

    integer :: j
    
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
       call mp_abort ("no howmany in init_ccfftw: FIX!")
    endif
# endif
    
# if FFT == _FFTW_
    j = fftw_in_place + fftw_measure + fftw_use_wisdom
    call fftw_f77_create_plan(fft%plan,n,is,j)
# elif FFT == _FFTW3_
    ! the planner expects this as an array of size 1
    array_n = n
    j = FFTW_PATIENT
    call dfftw_plan_dft_r2c_1d(fft%plan, n, data_in, data_out, j)
# endif

  end subroutine init_rcfftw

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
