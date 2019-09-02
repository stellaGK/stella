module stella_transforms

  use fft_work, only: fft_type

  implicit none

  public :: init_transforms, finish_transforms
  public :: transform_ky2y, transform_y2ky
  public :: transform_kx2x, transform_x2kx
  public :: transform_kalpha2alpha, transform_alpha2kalpha

  interface transform_ky2y
     module procedure transform_ky2y_5d
     module procedure transform_ky2y_2d
  end interface

  interface transform_y2ky
     module procedure transform_y2ky_5d
     module procedure transform_y2ky_2d
  end interface

  private

  type (fft_type) :: yf_fft, yb_fft
  type (fft_type) :: xf_fft, xb_fft
  type (fft_type) :: alpha_f_fft, alpha_b_fft

  logical :: transforms_initialized = .false.

  complex, dimension (:), allocatable :: fft_y_in, fft_y_out
  complex, dimension (:), allocatable :: fft_x_k
  real, dimension (:), allocatable :: fft_x_x
  ! arrays for transforming from alpha-space to k-alpha space
  real, dimension (:), allocatable :: fft_alpha_alpha
  complex, dimension (:), allocatable :: fft_alpha_kalpha

contains

  subroutine init_transforms

    use physics_flags, only: full_flux_surface
    use stella_layouts, only: init_stella_layouts

    implicit none

    if (transforms_initialized) return
    transforms_initialized = .true.

    call init_stella_layouts
    call init_y_fft
    call init_x_fft
    if (full_flux_surface) call init_alpha_fft

  end subroutine init_transforms

  subroutine init_y_fft

    use stella_layouts, only: vmu_lo
    use fft_work, only: init_ccfftw

    implicit none

    logical :: initialized = .false.
    
    if (initialized) return
    initialized = .true.

    if (.not.allocated(fft_y_in)) allocate (fft_y_in(vmu_lo%ny))
    if (.not.allocated(fft_y_out)) allocate (fft_y_out(vmu_lo%ny))

    call init_ccfftw (yf_fft,  1, vmu_lo%ny, fft_y_in, fft_y_out)
    call init_ccfftw (yb_fft, -1, vmu_lo%ny, fft_y_in, fft_y_out)

  end subroutine init_y_fft

  subroutine init_x_fft

    use stella_layouts, only: vmu_lo
    use fft_work, only: init_crfftw, init_rcfftw

    implicit none

    logical :: initialized = .false.
    
    if (initialized) return
    initialized = .true.

    if (.not.allocated(fft_x_k)) allocate (fft_x_k(vmu_lo%nx/2+1))
    if (.not.allocated(fft_x_x)) allocate (fft_x_x(vmu_lo%nx))

    call init_crfftw (xf_fft,  1, vmu_lo%nx, fft_x_k, fft_x_x)
    call init_rcfftw (xb_fft, -1, vmu_lo%nx, fft_x_x, fft_x_k)

  end subroutine init_x_fft

  subroutine init_alpha_fft

    use fft_work, only: init_rcfftw, init_crfftw
    use stella_layouts, only: vmu_lo

    implicit none

    logical :: initialized = .false.
    
    if (initialized) return
    initialized = .true.

    if (.not.allocated(fft_alpha_kalpha)) allocate (fft_alpha_kalpha(vmu_lo%nalpha/2+1))
    if (.not.allocated(fft_alpha_alpha)) allocate (fft_alpha_alpha(vmu_lo%nalpha))

    call init_crfftw (alpha_f_fft,  1, vmu_lo%nalpha, fft_alpha_kalpha, fft_alpha_alpha)
    call init_rcfftw (alpha_b_fft, -1, vmu_lo%nalpha, fft_alpha_alpha, fft_alpha_kalpha)

  end subroutine init_alpha_fft

  subroutine transform_ky2y_5d (gky_unpad, gy)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:,-vmu_lo%nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: gky_unpad
    complex, dimension (:,:,-vmu_lo%nzgrid:,:,vmu_lo%llim_proc:), intent (out) :: gy

    integer :: iky_max, ipad_up
    integer :: ikx, iz, it, ivmu

    ! first need to pad input array with zeros
    iky_max = vmu_lo%naky
    ipad_up = iky_max+vmu_lo%ny-(2*vmu_lo%naky-1)

    ! now fill in non-zero elements of array
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do it = 1, vmu_lo%ntubes
          do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
             do ikx = 1, vmu_lo%nakx/2+1
                fft_y_in(iky_max+1:ipad_up) = 0.
                fft_y_in(:iky_max) = gky_unpad(:iky_max,ikx,iz,it,ivmu)
                fft_y_in(ipad_up+1:) = gky_unpad(iky_max+1:,ikx,iz,it,ivmu)
                call dfftw_execute_dft(yf_fft%plan, fft_y_in, fft_y_out)
                fft_y_out = fft_y_out*yf_fft%scale
                gy(:,ikx,iz,it,ivmu) = fft_y_out
             end do
          end do
       end do
    end do

  end subroutine transform_ky2y_5d

  subroutine transform_ky2y_2d (gky_unpad, gy)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:), intent (in) :: gky_unpad
    complex, dimension (:,:), intent (out) :: gy

    integer :: iky_max, ipad_up
    integer :: ikx

    ! first need to pad input array with zeros
    iky_max = vmu_lo%naky
    ipad_up = iky_max+vmu_lo%ny-(2*vmu_lo%naky-1)
!    fft_y_in(iky_max+1:ipad_up) = 0.

    ! now fill in non-zero elements of array
    do ikx = 1, vmu_lo%nakx/2+1
       fft_y_in(iky_max+1:ipad_up) = 0.
       fft_y_in(:iky_max) = gky_unpad(:iky_max,ikx)
       fft_y_in(ipad_up+1:) = gky_unpad(iky_max+1:,ikx)
       call dfftw_execute_dft(yf_fft%plan, fft_y_in, fft_y_out)
       fft_y_out = fft_y_out*yf_fft%scale
       gy(:,ikx) = fft_y_out
    end do

  end subroutine transform_ky2y_2d

  subroutine transform_y2ky_5d (gy, gky)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:,-vmu_lo%nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gy
    complex, dimension (:,:,-vmu_lo%nzgrid:,:,vmu_lo%llim_proc:), intent (out) :: gky

    integer :: iky_max, ipad_up
    integer :: ikx, iz, it, ivmu

    iky_max = vmu_lo%naky
    ipad_up = iky_max+vmu_lo%ny-(2*vmu_lo%naky-1)
    ! now fill in non-zero elements of array
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do it = 1, vmu_lo%ntubes
          do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
             do ikx = 1, vmu_lo%nakx/2+1
                fft_y_in = gy(:,ikx,iz,it,ivmu)
                call dfftw_execute_dft(yb_fft%plan, fft_y_in, fft_y_out)
                fft_y_out = fft_y_out*yb_fft%scale
                gky(:iky_max,ikx,iz,it,ivmu) = fft_y_out(:iky_max)
                gky(iky_max+1:,ikx,iz,it,ivmu) = fft_y_out(ipad_up+1:)
             end do
          end do
       end do
    end do

  end subroutine transform_y2ky_5d

  subroutine transform_y2ky_2d (gy, gky)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:), intent (in out) :: gy
    complex, dimension (:,:), intent (out) :: gky

    integer :: iky_max, ipad_up
    integer :: ikx

    iky_max = vmu_lo%naky
    ipad_up = iky_max+vmu_lo%ny-(2*vmu_lo%naky-1)
    ! now fill in non-zero elements of array
    do ikx = 1, vmu_lo%nakx/2+1
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

    ! now fill in non-zero elements of array
    do iy = 1, vmu_lo%ny
       ! first need to pad input array with zeros
       fft_x_k(vmu_lo%nakx/2+2:) = 0.
       fft_x_k(:vmu_lo%nakx/2+1) = gkx(iy,:)
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
       gkx(iy,:) = fft_x_k(:vmu_lo%nakx/2+1)
    end do

  end subroutine transform_x2kx

  subroutine transform_kalpha2alpha (gkalph, galph)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:), intent (in) :: gkalph
    real, dimension (:), intent (out) :: galph

    ! first need to pad input array with zeros
    fft_alpha_kalpha(vmu_lo%naky+1:) = 0.
    ! then fill in non-zero elements of array
    fft_alpha_kalpha(:vmu_lo%naky) = gkalph
    call dfftw_execute_dft_c2r(alpha_f_fft%plan, fft_alpha_kalpha, fft_alpha_alpha)
    fft_alpha_alpha = fft_alpha_alpha*alpha_f_fft%scale
    galph = fft_alpha_alpha

  end subroutine transform_kalpha2alpha

  subroutine transform_alpha2kalpha (galph, gkalph)

    use stella_layouts, only: vmu_lo

    implicit none

    real, dimension (:), intent (in) :: galph
    complex, dimension (:), intent (out) :: gkalph

    fft_alpha_alpha = galph
    call dfftw_execute_dft_r2c(alpha_b_fft%plan, fft_alpha_alpha, fft_alpha_kalpha)
    fft_alpha_kalpha = fft_alpha_kalpha*alpha_b_fft%scale
    ! filter out highest k-alpha modes to avoid aliasing
    gkalph = fft_alpha_kalpha(:vmu_lo%naky)

  end subroutine transform_alpha2kalpha

  subroutine finish_transforms

    use physics_flags, only: full_flux_surface

    implicit none

    call dfftw_destroy_plan (yf_fft%plan)
    call dfftw_destroy_plan (yb_fft%plan)
    call dfftw_destroy_plan (xf_fft%plan)
    call dfftw_destroy_plan (xb_fft%plan)
    if (full_flux_surface) then
       call dfftw_destroy_plan (alpha_f_fft%plan)
       call dfftw_destroy_plan (alpha_b_fft%plan)
    end if
    if (allocated(fft_y_in)) deallocate (fft_y_in)
    if (allocated(fft_y_out)) deallocate (fft_y_out)
    if (allocated(fft_x_k)) deallocate (fft_x_k)
    if (allocated(fft_x_x)) deallocate (fft_x_x)
    if (allocated(fft_alpha_alpha)) deallocate (fft_alpha_alpha)
    if (allocated(fft_alpha_kalpha)) deallocate (fft_alpha_kalpha)
    transforms_initialized = .false.

  end subroutine finish_transforms

end module stella_transforms
