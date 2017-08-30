# include "define.inc"

module stella_transforms

  use fft_work, only: fft_type

  implicit none

  public :: init_transforms, finish_transforms
  public :: transform_ky2y, transform_y2ky
  public :: transform_kx2x, transform_x2kx

  interface transform_ky2y
     module procedure transform_ky2y_4d
     module procedure transform_ky2y_2d
  end interface

  interface transform_y2ky
     module procedure transform_y2ky_4d
     module procedure transform_y2ky_2d
  end interface

  private

  type (fft_type) :: yf_fft, yb_fft
  type (fft_type) :: xf_fft, xb_fft

  logical :: transforms_initialized = .false.

  complex, dimension (:), allocatable :: fft_y_in, fft_y_out
  complex, dimension (:), allocatable :: fft_x_k
  real, dimension (:), allocatable :: fft_x_x

contains

  subroutine init_transforms

    use stella_layouts, only: init_stella_layouts

    implicit none

    if (transforms_initialized) return
    transforms_initialized = .true.

    call init_stella_layouts
    call init_y_fft
    call init_x_fft

  end subroutine init_transforms

  subroutine init_y_fft

    use stella_layouts, only: vmu_lo
    use fft_work, only: init_ccfftw

    implicit none

    logical :: initialized = .false.
# if FFT == _FFTW3_
    integer :: nb_ffts
# endif
    
    if (initialized) return
    initialized = .true.

    if (.not.allocated(fft_y_in)) allocate (fft_y_in(vmu_lo%ny))
    if (.not.allocated(fft_y_out)) allocate (fft_y_out(vmu_lo%ny))
# if FFT == _FFTW_
    call init_ccfftw (yf_fft,  1, vmu_lo%ny)
    call init_ccfftw (yb_fft, -1, vmu_lo%ny)
# elif FFT == _FFTW3_
    ! number of ffts to be calculated
    ! equal to number of elements off-processor * number of non-y elements on processor
    nb_ffts = (vmu_lo%ulim_alloc - vmu_lo%llim_proc + 1)*vmu_lo%nakx*vmu_lo%nzed

    call init_ccfftw (yf_fft,  1, vmu_lo%ny, nb_ffts, fft_y_in, fft_y_out)
    call init_ccfftw (yb_fft, -1, vmu_lo%ny, nb_ffts, fft_y_in, fft_y_out)
# endif

  end subroutine init_y_fft

  subroutine init_x_fft

    use stella_layouts, only: vmu_lo
    use fft_work, only: init_crfftw, init_rcfftw

    implicit none

    logical :: initialized = .false.
# if FFT == _FFTW3_
    integer :: nb_ffts
# endif
    
    if (initialized) return
    initialized = .true.

    if (.not.allocated(fft_x_k)) allocate (fft_x_k(vmu_lo%nx/2+1))
    if (.not.allocated(fft_x_x)) allocate (fft_x_x(vmu_lo%nx))
# if FFT == _FFTW_
    call init_crfftw (xf_fft,  1, vmu_lo%nx)
    call init_rcfftw (xb_fft, -1, vmu_lo%nx)
# elif FFT == _FFTW3_
    ! number of ffts to be calculated
    ! equal to number of elements off-processor * number of non-xy elements on processor
    nb_ffts = (vmu_lo%ulim_alloc - vmu_lo%llim_proc + 1)*vmu_lo%nzed

    call init_crfftw (xf_fft,  1, vmu_lo%nx, nb_ffts, fft_x_k, fft_x_x)
    call init_rcfftw (xb_fft, -1, vmu_lo%nx, nb_ffts, fft_x_x, fft_x_k)
# endif

  end subroutine init_x_fft

  subroutine transform_ky2y_4d (gky_unpad, gy)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:,-vmu_lo%nzgrid:,vmu_lo%llim_proc:), intent (in) :: gky_unpad
    complex, dimension (:,:,-vmu_lo%nzgrid:,vmu_lo%llim_proc:), intent (out) :: gy

    integer :: iky_max, ipad_up
    integer :: ikx, ig, ivmu

    ! first need to pad input array with zeros
    iky_max = vmu_lo%naky/2+1
    ipad_up = iky_max+vmu_lo%ny-vmu_lo%naky
    fft_y_in(iky_max+1:ipad_up) = 0.

    ! now fill in non-zero elements of array
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do ig = -vmu_lo%nzgrid, vmu_lo%nzgrid
          do ikx = 1, vmu_lo%nakx
             fft_y_in(:iky_max) = gky_unpad(:iky_max,ikx,ig,ivmu)
             fft_y_in(ipad_up+1:) = gky_unpad(iky_max+1:,ikx,ig,ivmu)
             call dfftw_execute_dft(yf_fft%plan, fft_y_in, fft_y_out)
             fft_y_out = fft_y_out*yf_fft%scale
             gy(:,ikx,ig,ivmu) = fft_y_out
          end do
       end do
    end do

  end subroutine transform_ky2y_4d

  subroutine transform_ky2y_2d (gky_unpad, gy)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:), intent (in) :: gky_unpad
    complex, dimension (:,:), intent (out) :: gy

    integer :: iky_max, ipad_up
    integer :: ikx

    ! first need to pad input array with zeros
    iky_max = vmu_lo%naky/2+1
    ipad_up = iky_max+vmu_lo%ny-vmu_lo%naky
    fft_y_in(iky_max+1:ipad_up) = 0.

    ! now fill in non-zero elements of array
    do ikx = 1, vmu_lo%nakx
       fft_y_in(:iky_max) = gky_unpad(:iky_max,ikx)
       fft_y_in(ipad_up+1:) = gky_unpad(iky_max+1:,ikx)
       call dfftw_execute_dft(yf_fft%plan, fft_y_in, fft_y_out)
       fft_y_out = fft_y_out*yf_fft%scale
       gy(:,ikx) = fft_y_out
    end do

  end subroutine transform_ky2y_2d

  subroutine transform_y2ky_4d (gy, gky)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:,-vmu_lo%nzgrid:,vmu_lo%llim_proc:), intent (in out) :: gy
    complex, dimension (:,:,-vmu_lo%nzgrid:,vmu_lo%llim_proc:), intent (out) :: gky

    integer :: iky_max, ipad_up
    integer :: ikx, ig, ivmu

    iky_max = vmu_lo%naky/2+1
    ipad_up = iky_max+vmu_lo%ny-vmu_lo%naky
    ! now fill in non-zero elements of array
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do ig = -vmu_lo%nzgrid, vmu_lo%nzgrid
          do ikx = 1, vmu_lo%nakx
             fft_y_in = gy(:,ikx,ig,ivmu)
             call dfftw_execute_dft(yb_fft%plan, fft_y_in, fft_y_out)
             fft_y_out = fft_y_out*yb_fft%scale
             gky(:iky_max,ikx,ig,ivmu) = fft_y_out(:iky_max)
             gky(iky_max+1:,ikx,ig,ivmu) = fft_y_out(ipad_up+1:)
          end do
       end do
    end do

  end subroutine transform_y2ky_4d

  subroutine transform_y2ky_2d (gy, gky)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:), intent (in out) :: gy
    complex, dimension (:,:), intent (out) :: gky

    integer :: iky_max, ipad_up
    integer :: ikx

    iky_max = vmu_lo%naky/2+1
    ipad_up = iky_max+vmu_lo%ny-vmu_lo%naky
    ! now fill in non-zero elements of array
    do ikx = 1, vmu_lo%nakx
       fft_y_in = gy(:,ikx)
       call dfftw_execute_dft(yb_fft%plan, fft_y_in, fft_y_out)
       fft_y_out = fft_y_out*yb_fft%scale
       gky(:iky_max,ikx) = fft_y_out(:iky_max)
       gky(iky_max+1:,ikx) = fft_y_out(ipad_up+1:)
    end do

  end subroutine transform_y2ky_2d

  subroutine transform_kx2x (gkx, gx)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:), intent (in) :: gkx
    real, dimension (:,:), intent (out) :: gx

    integer :: iy

    ! first need to pad input array with zeros
    fft_x_k(vmu_lo%nakx+1:) = 0.

    ! now fill in non-zero elements of array
    do iy = 1, vmu_lo%ny
       fft_x_k(:vmu_lo%nakx) = gkx(iy,:)
       call dfftw_execute_dft_c2r(xf_fft%plan, fft_x_k, fft_x_x)
       fft_x_x = fft_x_x*xf_fft%scale
       gx(iy,:) = fft_x_x
    end do

  end subroutine transform_kx2x

  subroutine transform_x2kx (gx, gkx)

    use stella_layouts, only: vmu_lo

    implicit none

    real, dimension (:,:), intent (in) :: gx
    complex, dimension (:,:), intent (out) :: gkx

    integer :: iy

    do iy = 1, vmu_lo%ny
       fft_x_x = gx(iy,:)
       call dfftw_execute_dft_r2c(xb_fft%plan, fft_x_x, fft_x_k)
       fft_x_k = fft_x_k*xb_fft%scale
       gkx(iy,:) = fft_x_k(:vmu_lo%nakx)
    end do

  end subroutine transform_x2kx

  subroutine finish_transforms

    implicit none

    call dfftw_destroy_plan (yf_fft%plan)
    call dfftw_destroy_plan (yb_fft%plan)
    call dfftw_destroy_plan (xf_fft%plan)
    call dfftw_destroy_plan (xb_fft%plan)
    if (allocated(fft_y_in)) deallocate (fft_y_in)
    if (allocated(fft_y_out)) deallocate (fft_y_out)
    if (allocated(fft_x_k)) deallocate (fft_x_k)
    if (allocated(fft_x_x)) deallocate (fft_x_x)
    transforms_initialized = .false.

  end subroutine finish_transforms

end module stella_transforms
