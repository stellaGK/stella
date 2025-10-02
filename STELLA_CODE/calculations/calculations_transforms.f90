!###############################################################################
!                          FOURIER TRANSFORMATIONS
!###############################################################################
! 
! Apply discrete Fourier transformations of both real and complex data, with
! the help of the package FFTW, which stands for the "Fastest Fourier Transform
! in the West". The package takes advantage of using multiple processors.
! See module "utils/fft_work.f90" for details on the FFTW.
! 
! Stella considers (nx,ny) modes in xyz-space and (nakx,naky) modes in k-space.
! To make the nx/ny modes in a centered convolution free of aliasing errors
! we consider nakx/naky de-aliased modes, which is 2/3 of the total number of modes.
! The remainder of the vector is padded with zeros to get vectors of length nx/ny.
! Along naky we only consider positive modes due the reality condition
!      naky = (ny-1)/3 + 1          ordered like [0, dky, ..., kymax]
!      nakx = 2*((nx-1)/3) +  1     ordered like [0, dkx, ..., kxmax, -kxmax, ..., -dkx]
! 
! If the input data is real, the output satisfies the Hermitian redundancy
! and out[i] = conjugate(out[n-i]) so there are only n/2+1 redundant outputs
! 
! This module contains 6 Fourier transforms in sets of 2 where we first transform
! our quantity along x and then along y, or vice versa. Each of these transforms
! consists of a forward transform (k-space -> xyz-space) and a backward transform.
! 
!  (1) FFT along ny, followed by a FFT along nx, with zero-padding
!         I. Fourier transformation along ny while the x-axis is still in k-space
!               yf_fft      ky --C2C--> y      transform_ky2y_2d/5d (gky_unpad, gy)
!               yb_fft      y  --C2C--> ky     transform_y2ky_2d/5d (gy, gky)
!         II. Fourier transformation along nx while the y-axis is already in xyz-space
!               xf_fft      kx --C2R--> x      transform_kx2x (gkx, gx)
!               xb_fft      x  --RC2--> kx     transform_x2kx (gx, gkx)
! 
!  (2) FFT along nx, followed by a FFT along ny, with zero-padding
!         I. Fourier transformation along nx while the y-axis is still in k-space
!               xsf_fft    kx --C2C--> x       transform_kx2x_xfirst (gkx, gx)
!               xsb_fft    x  --C2C--> kx      transform_x2kx_xfirst (gx, gkx)
!         II. Fourier transformation along ny while the x-axis is already in xyz-space
!               ysf_fft    ky --C2R--> y       transform_ky2y_xfirst (gky, gy)
!               ysb_fft    y  --RC2--> ky      transform_y2ky_xfirst (gy, gky)
! 
!  (3) FFT along nakx, followed by a FFT along naky, without zero-padding
!         I. Fourier transformation along nakx while the y-axis is still in k-space
!               xfnp_fft   kx --C2C--> x       transform_kx2x_unpadded (gkx, gx)
!               xbnp_fft   x  --C2C--> kx      transform_x2kx_unpadded (gx, gkx)
!         II. Fourier transformation along naky while the x-axis is already in xyz-space
!               yfnp_fft   ky --C2R--> y       transform_ky2y_unpadded (gky, gy)
!               ybnp_fft   y  --RC2--> ky      transform_y2ky_unpadded (gy, gky)
! 
!  (4) FFT along nalpha = ny if full_flux_surface==True
!          alpha_f_fft   kalpha --C2R--> alpha    transform_kalpha2alpha (gkalph, galph)
!          alpha_b_fft   alpha  --RC2--> kalpha   transform_alpha2kalpha (galph, gkalph)
! 
!###############################################################################
module calculations_transforms

   use fft_work, only: fft_type

   implicit none

   ! Make the Fourier transformations available to other modules
   public :: init_transforms, finish_transforms
   public :: transform_ky2y, transform_y2ky
   public :: transform_kx2x, transform_x2kx
   public :: transform_ky2y_unpadded, transform_y2ky_unpadded
   public :: transform_kx2x_unpadded, transform_x2kx_unpadded
   public :: transform_kx2x_xfirst, transform_x2kx_xfirst
   public :: transform_ky2y_xfirst, transform_y2ky_xfirst
   public :: transform_kalpha2alpha, transform_alpha2kalpha

   private

   interface transform_ky2y
      module procedure transform_ky2y_5d
      module procedure transform_ky2y_2d
   end interface

   interface transform_y2ky
      module procedure transform_y2ky_5d
      module procedure transform_y2ky_2d
   end interface

   ! Stores the FFTW plans and the scaling factor
   type(fft_type) :: yf_fft, yb_fft
   type(fft_type) :: xf_fft, xb_fft
   type(fft_type) :: yfnp_fft, ybnp_fft
   type(fft_type) :: xfnp_fft, xbnp_fft
   type(fft_type) :: xsf_fft, xsb_fft
   type(fft_type) :: ysf_fft, ysb_fft
   type(fft_type) :: alpha_f_fft, alpha_b_fft

   ! FFT along ny while the x-axis is still in k-space
   complex, dimension(:), allocatable :: fft_y_in        !   dimension: ny
   complex, dimension(:), allocatable :: fft_y_out          ! dimension: ny
   
   ! FFT along nx while the y-axis is already in xyz-space
   complex, dimension(:), allocatable :: fft_x_k            ! dimension: nx/2+1
   real, dimension(:), allocatable :: fft_x_x               ! dimension: nx
   
   ! FFT along nx while the y-axis is still in k-space
   complex, dimension(:), allocatable :: fft_xs_k           ! dimension: nx
   complex, dimension(:), allocatable :: fft_xs_x           ! dimension: nx
   
   ! FFT along ny while the x-axis is already in xyz-space
   complex, dimension(:), allocatable :: fft_ys_k           ! dimension: ny/2+1
   real, dimension(:), allocatable :: fft_ys_y              ! dimension: ny
   
   ! FFT along nakx while the y-axis is still in k-space, without zero-padding
   complex, dimension(:), allocatable :: fftnp_x_k          ! dimension: nakx
   complex, dimension(:), allocatable :: fftnp_x_x          ! dimension: nakx
   
   ! FFT along naky while the x-axis is already in xyz-space, without zero-padding
   complex, dimension(:), allocatable :: fftnp_y_k          ! dimension: naky
   real, dimension(:), allocatable :: fftnp_y_y             ! dimension: 2*naky-1
      
   ! Arrays for transforming from alpha-space to k-alpha space
   real, dimension(:), allocatable :: fft_alpha_alpha       ! dimension: nalpha/2+1
   complex, dimension(:), allocatable :: fft_alpha_kalpha   ! dimension: nalpha

   ! Only initialise the FFTW plans once
   logical :: transforms_initialized = .false.

contains

!###############################################################################
!                              INITIALIZATIONS
!###############################################################################

  !=============================================================================
  !==================== INITIALIZE THE FOURIER TRANSFORMS !=====================
  !=============================================================================
  ! Allocate the input and output arrays of the Fourier transforms
  ! and create the FFT plans that are used to apply the transformations.
  !=============================================================================
   subroutine init_transforms

      use parameters_physics, only: full_flux_surface
      use parallelisation_layouts, only: read_parameters_parallelisation_layouts

      implicit none

      if (transforms_initialized) return
      transforms_initialized = .true.

      call read_parameters_parallelisation_layouts
      call init_y_fft
      call init_x_fft
      call init_x_xfirst_fft
      call init_y_xfirst_fft
      call init_unpadded_x_fft
      call init_unpadded_y_fft
      if (full_flux_surface) call init_alpha_fft

   end subroutine init_transforms

   !****************************************************************************
   !                 CREATE PLAN FOR A C2C FFT ALONG NY
   !****************************************************************************
   ! Fourier transformation along ny while the x-axis is still in k-space
   ! Create the FFT plans for a complex to complex fourier transform along ny:
   !    <yf_fft>: "forward" Fourier Transform for ky(complex) --> y(complex)
   !    <yb_fft>: "backward" Fourier Transform for y(complex) --> ky(complex)
   !****************************************************************************
   subroutine init_y_fft

      use parallelisation_layouts, only: vmu_lo
      use fft_work, only: init_ccfftw

      implicit none

      logical :: initialized = .false.

      if (initialized) return
      initialized = .true.

      ! Allocate the input and output arrays of the transformation
      ! Since the xyz-space is complex, the k-space contains all frequencies
      if (.not. allocated(fft_y_in)) allocate (fft_y_in(vmu_lo%ny))
      if (.not. allocated(fft_y_out)) allocate (fft_y_out(vmu_lo%ny))

      ! Create the plans for the Fourier transforms along ny
      call init_ccfftw(yf_fft,  1, vmu_lo%ny, fft_y_in, fft_y_out)
      call init_ccfftw(yb_fft, -1, vmu_lo%ny, fft_y_in, fft_y_out)

   end subroutine init_y_fft

   !****************************************************************************
   !                  CREATE PLAN FOR A R2C FFT ALONG NX 
   !****************************************************************************
   ! Fourier transformation along nx while the y-axis is already in xyz-space
   ! Create the FFT plans for a r2c or c2r fourier transform along nx:
   !    <xf_fft>: "forward" Fourier Transform for kx(complex) --> x(real)
   !    <xb_fft>: "backward" Fourier Transform for x(real) --> kx(complex)
   !****************************************************************************
   subroutine init_x_fft

      use parallelisation_layouts, only: vmu_lo
      use fft_work, only: init_crfftw, init_rcfftw

      implicit none

      logical :: initialized = .false.

      if (initialized) return
      initialized = .true.

      ! Allocate the input and output arrays of the transformation
      ! Since the xyz-space is real, the k-space only contains positive frequencies
      if (.not. allocated(fft_x_k)) allocate (fft_x_k(vmu_lo%nx / 2 + 1))  ! complex kx
      if (.not. allocated(fft_x_x)) allocate (fft_x_x(vmu_lo%nx))          ! real x

      ! Create the plans for the Fourier transforms along nx
      ! c2r DFTs are always FFTW_BACKWARD --> use the opposite naming convention
      ! r2c DFTs are always FFTW_FORWARD --> use the opposite naming convention
      call init_crfftw(xf_fft,  1, vmu_lo%nx, fft_x_k, fft_x_x)   ! FFTW_BACKWARD
      call init_rcfftw(xb_fft, -1, vmu_lo%nx, fft_x_x, fft_x_k)   ! FFTW_FORWARD

   end subroutine init_x_fft


   !****************************************************************************
   !                    CREATE PLAN FOR A C2C FFT ALONG NX 
   !****************************************************************************
   ! Fourier transformation along nx while the y-axis is still in k-space
   ! Create the FFT plans for a complex to complex fourier transform along nx:
   !    <xsf_fft>: "forward" Fourier Transform for kx(complex) --> x(complex)
   !    <xsb_fft>: "backward" Fourier Transform for x(complex) --> kx(complex)
   !****************************************************************************
   subroutine init_x_xfirst_fft

      use parallelisation_layouts, only: vmu_lo
      use fft_work, only: init_ccfftw

      implicit none

      logical :: initialized = .false.

      if (initialized) return
      initialized = .true.

      ! Allocate the input and output arrays of the transformation
      ! Since the xyz-space is complex, the k-space contains all frequencies
      if (.not. allocated(fft_xs_k)) allocate (fft_xs_k(vmu_lo%nx)) ! complex kx
      if (.not. allocated(fft_xs_x)) allocate (fft_xs_x(vmu_lo%nx)) ! complex x

      ! Create the plans for the Fourier transforms along nx
      call init_ccfftw(xsf_fft,  1, vmu_lo%nx, fft_xs_k, fft_xs_x)
      call init_ccfftw(xsb_fft, -1, vmu_lo%nx, fft_xs_x, fft_xs_k)

   end subroutine init_x_xfirst_fft

   !****************************************************************************
   !                    CREATE PLAN FOR A R2C FFT ALONG NY  
   !****************************************************************************
   ! Fourier transformation along ny while the x-axis is already in xyz-space
   ! Create the FFT plans for a r2c or c2r fourier transform along ny:
   !    <ysf_fft>: "forward" Fourier Transforms for ky(complex) --> y(real)
   !    <ysb_fft>: "backward" Fourier Transforms for y(real) --> ky(complex)
   !****************************************************************************
   subroutine init_y_xfirst_fft

      use parallelisation_layouts, only: vmu_lo
      use fft_work, only: init_crfftw, init_rcfftw

      implicit none

      logical :: initialized = .false.

      if (initialized) return
      initialized = .true.

      ! Allocate the input and output arrays of the transform
      ! Since the xyz-space is real, the k-space only contains positive frequencies
      if (.not. allocated(fft_ys_k)) allocate (fft_ys_k(vmu_lo%ny / 2 + 1))  ! complex ky
      if (.not. allocated(fft_ys_y)) allocate (fft_ys_y(vmu_lo%ny))          ! real y

      ! Create the plans for the Fourier transforms along ny
      ! c2r DFTs are always FFTW_BACKWARD --> use the opposite naming convention
      ! r2c DFTs are always FFTW_FORWARD --> use the opposite naming convention
      call init_crfftw(ysf_fft,  1, vmu_lo%ny, fft_ys_k, fft_ys_y)   ! FFTW_BACKWARD
      call init_rcfftw(ysb_fft, -1, vmu_lo%ny, fft_ys_y, fft_ys_k)   ! FFTW_FORWARD

   end subroutine init_y_xfirst_fft

   !****************************************************************************
   !                  CREATE PLAN FOR A C2C FFT ALONG NAKX 
   !****************************************************************************
   ! Fourier transformation along nakx while the y-axis is still in k-space
   ! Create the FFT plans for a c2c fourier transform along nakx:
   !    <xfnp_fft>: "forward" Fourier Transforms for kx(complex) --> x(complex)
   !    <xbnp_fft>: "backward" Fourier Transforms for x(complex) --> kx(complex)
   !****************************************************************************
   subroutine init_unpadded_x_fft

      use parallelisation_layouts, only: vmu_lo
      use fft_work, only: init_ccfftw, FFT_BACKWARD, FFT_FORWARD

      implicit none

      ! Allocate the input and output arrays of the transform without zero-padding
      ! Since the xyz-space is complex, the k-space contains all frequencies
      if (.not. allocated(fftnp_x_k)) allocate (fftnp_x_k(vmu_lo%nakx)) ! complex kx
      if (.not. allocated(fftnp_x_x)) allocate (fftnp_x_x(vmu_lo%nakx)) ! complex x

      ! Create the plans for the Fourier transforms along nakx
      ! forward Fourier uses FFTW_BACKWARD --> use the opposite naming convention
      ! backward Fourier uses FFT_FORWARD --> use the opposite naming convention
      call init_ccfftw(xfnp_fft, FFT_BACKWARD, vmu_lo%nakx, fftnp_x_k, fftnp_x_x)
      call init_ccfftw(xbnp_fft, FFT_FORWARD, vmu_lo%nakx, fftnp_x_x, fftnp_x_k)

   end subroutine init_unpadded_x_fft

   !****************************************************************************
   !                  CREATE PLAN FOR A C2R FFT ALONG NAKY  
   !****************************************************************************
   ! Fourier transformation along naky while the x-axis is already in xyz-space
   ! Create the FFT plans for a c2c fourier transform along naky:
   !    <yfnp_fft>: "forward" Fourier Transforms for ky(complex) --> y(real)
   !    <ybnp_fft>: "backward" Fourier Transforms for y(real) --> ky(complex)
   !****************************************************************************
   subroutine init_unpadded_y_fft

      use parallelisation_layouts, only: vmu_lo
      use fft_work, only: init_crfftw, init_rcfftw, FFT_BACKWARD, FFT_FORWARD

      implicit none

      ! Allocate the input and output arrays of the transform without zero-padding
      ! Since the xyz-space is real, the k-space only contains positive frequencies
      if (.not. allocated(fftnp_y_k)) allocate (fftnp_y_k(vmu_lo%naky))      ! complex ky
      if (.not. allocated(fftnp_y_y)) allocate (fftnp_y_y(2*vmu_lo%naky-1))  ! real y

    ! Create the plans for the Fourier transforms along nakx
    ! forward Fourier uses FFTW_BACKWARD --> use the opposite naming convention
    ! backward Fourier uses FFT_FORWARD --> use the opposite naming convention
      call init_crfftw(yfnp_fft, FFT_BACKWARD, 2*vmu_lo%naky-1, fftnp_y_k, fftnp_y_y)
      call init_rcfftw(ybnp_fft, FFT_FORWARD, 2*vmu_lo%naky-1, fftnp_y_y, fftnp_y_k)

   end subroutine init_unpadded_y_fft

   !****************************************************************************
   !                 CREATE PLAN FOR A C2R FFT ALONG NALPHA  
   !****************************************************************************
   ! Recall that if (full_flux_surface) then nalpha = ny
   ! Create the FFT plans for a c2r fourier transform along nalpha:
   !    <alpha_f_fft>: "forward" Fourier Transforms for kalpha(complex) --> alpha(real)
   !    <alpha_b_fft>: "backward" Fourier Transforms for alpha(real) --> kalpha(complex)
   !****************************************************************************
   subroutine init_alpha_fft

      use fft_work, only: init_rcfftw, init_crfftw
      use fft_work, only: fft_backward, fft_forward
      use parallelisation_layouts, only: vmu_lo

      implicit none

      logical :: initialized = .false.

      if (initialized) return
      initialized = .true.

      ! Allocate the input and output arrays of the transform
      ! Since the xyz-space is real, the k-space only contains positive frequencies
      if (.not. allocated(fft_alpha_kalpha)) allocate (fft_alpha_kalpha(vmu_lo%nalpha / 2 + 1))
      if (.not. allocated(fft_alpha_alpha)) allocate (fft_alpha_alpha(vmu_lo%nalpha))

      ! Create the plans for the Fourier transforms along nalpha
      call init_crfftw(alpha_f_fft, fft_backward, vmu_lo%nalpha, fft_alpha_kalpha, fft_alpha_alpha)
      call init_rcfftw(alpha_b_fft, fft_forward, vmu_lo%nalpha, fft_alpha_alpha, fft_alpha_kalpha)

   end subroutine init_alpha_fft

!###############################################################################
!                              TRANSFORMATIONS
!###############################################################################

   !****************************************************************************
   !                      FFT:   KY --C2C--> Y   ALONG NY  
   !****************************************************************************
   ! Fourier transformation along ny while the x-axis is still in k-space
   !     gky_unpad(2*naky-1, nakx/2+1)     with all ky and kx >= 0
   !     gy(ny, nakx/2+1)                  with all  y and kx >= 0
   !****************************************************************************
   subroutine transform_ky2y_5d(gky_unpad, gy)

      use parallelisation_layouts, only: vmu_lo

      implicit none

      complex, dimension(:, :, -vmu_lo%nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gky_unpad
      complex, dimension(:, :, -vmu_lo%nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: gy

      integer :: iky_max, ipad_up
      integer :: ikx, iz, it, ivmu

      ! Get the indices of the last positive ky <iky_max> and negative ky <ipad_up>
      iky_max = vmu_lo%naky
      ipad_up = iky_max + vmu_lo%ny - (2 * vmu_lo%naky - 1)

      ! Iterate over vmu, tubes, z and positive kx since we Fourier transform along ny
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, vmu_lo%ntubes
            do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
               do ikx = 1, vmu_lo%nakx / 2 + 1
                  ! Pad the middle of the array with zeros: zero padding in the wavenumber
                  ! domain results in an increased sampling rate in the real space domain
                  fft_y_in(iky_max + 1:ipad_up) = 0.
                  ! Fill in the non-zero elements with gky_unpad(2*naky-1, nakx/2+1)
                  ! Add the positive frequencies gky_unpad[:naky] to the start
                  ! And the negative frequencies gky_unpad[naky+1:] to the end
                  fft_y_in(:iky_max) = gky_unpad(:iky_max, ikx, iz, it, ivmu)
                  fft_y_in(ipad_up + 1:) = gky_unpad(iky_max + 1:, ikx, iz, it, ivmu)
                  ! Apply the "forward" FFT for ky(complex) --> y(complex) along ny
                  ! The result is still complex since the x-axis remains in reciprocal space
                  call dfftw_execute_dft(yf_fft%plan, fft_y_in, fft_y_out)
                  ! Do not rescale the Fourier components since yf_fft%scale=1
                  fft_y_out = fft_y_out * yf_fft%scale
                  ! Return gy(ny, nakx/2+1) with y in real space and x in reciprocal space
                  gy(:, ikx, iz, it, ivmu) = fft_y_out
               end do
            end do
         end do
      end do

   end subroutine transform_ky2y_5d

   !****************************************************************************
   !                  FFT:   KY --C2C--> Y   ALONG NY  
   !****************************************************************************
   ! Fourier transformation along ny while the x-axis is still in k-space
   !     gky_unpad(2*naky-1, nakx/2+1)     with all ky and kx >= 0
   !     gy(ny, nakx/2+1)                  with all  y and kx >= 0
   !****************************************************************************
   subroutine transform_ky2y_2d(gky_unpad, gy)

      use parallelisation_layouts, only: vmu_lo

      implicit none

      complex, dimension(:, :), intent(in) :: gky_unpad
      complex, dimension(:, :), intent(out) :: gy

      integer :: iky_max, ipad_up
      integer :: ikx

      ! Get the indices of the last positive ky <iky_max> and negative ky <ipad_up>
      iky_max = vmu_lo%naky
      ipad_up = iky_max + vmu_lo%ny - (2 * vmu_lo%naky - 1)

      ! Iterate over positive kx since we Fourier transform along ny
      do ikx = 1, vmu_lo%nakx / 2 + 1
         ! Pad the middle of the array with zeros: zero padding in the wavenumber
         ! domain results in an increased sampling rate in the real space domain
         fft_y_in(iky_max + 1:ipad_up) = 0.
         ! Fill in the non-zero elements with gky_unpad(2*naky-1, nakx/2+1)
         ! Add the positive frequencies gky_unpad[:naky] to the start
         ! And the negative frequencies gky_unpad[naky+1:] to the end
         fft_y_in(:iky_max) = gky_unpad(:iky_max, ikx)
         fft_y_in(ipad_up + 1:) = gky_unpad(iky_max + 1:, ikx)
         ! Apply the "forward" FFT for ky(complex) --> y(complex) along ny
         ! The result is still complex since the x-axis remains in reciprocal space
         call dfftw_execute_dft(yf_fft%plan, fft_y_in, fft_y_out)
         ! Do not rescale the Fourier components since yf_fft%scale=1
         fft_y_out = fft_y_out * yf_fft%scale
         ! Return gy(ny, nakx/2+1) with y in real space and x in reciprocal space
         gy(:, ikx) = fft_y_out
      end do

   end subroutine transform_ky2y_2d

   !****************************************************************************
   !                   FFT:   Y --C2C--> KY   ALONG NY  
   !****************************************************************************
   ! Fourier transformation along ny while the x-axis is already in k-space
   !     gy(ny, nakx/2+1, z, tubes, vmu)            with all  y and kx >= 0
   !     gky(2*naky-1, nakx/2+1, z, tubes, vmu)     with all ky and kx >= 0
   !****************************************************************************
   subroutine transform_y2ky_5d(gy, gky)

      use parallelisation_layouts, only: vmu_lo

      implicit none

      complex, dimension(:, :, -vmu_lo%nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gy
      complex, dimension(:, :, -vmu_lo%nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: gky

      integer :: iky_max, ipad_up
      integer :: ikx, iz, it, ivmu

      ! Get the indices of the last positive ky <iky_max> and negative ky <ipad_up>
      iky_max = vmu_lo%naky
      ipad_up = iky_max + vmu_lo%ny - (2 * vmu_lo%naky - 1)
      
    ! Iterate over vmu, tubes, z and positive kx since we Fourier transform along ny
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, vmu_lo%ntubes
            do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
               do ikx = 1, vmu_lo%nakx / 2 + 1
                  ! Fill in the non-zero elements with gy(ny, nakx/2+1)
                  fft_y_in = gy(:, ikx, iz, it, ivmu)
                  ! Apply the "backward" FFT for y(complex) --> ky(complex) along ny
                  call dfftw_execute_dft(yb_fft%plan, fft_y_in, fft_y_out)
                  ! Rescale the Fourier components with xf_fft%scale=1/ny
                  fft_y_out = fft_y_out * yb_fft%scale
                  ! Return gky(2*naky-1, nakx/2+1) with x and y in reciprocal space
                  ! The positive frequencies are located at the start
                  ! The negative frequencies are located at the end
                  gky(:iky_max, ikx, iz, it, ivmu) = fft_y_out(:iky_max)
                  gky(iky_max + 1:, ikx, iz, it, ivmu) = fft_y_out(ipad_up + 1:)
               end do
            end do
         end do
      end do

   end subroutine transform_y2ky_5d

   !****************************************************************************
   !                     FFT:   Y --C2C--> KY   ALONG NY  
   !****************************************************************************
   ! Fourier transformation along ny while the x-axis is already in k-space
   !     gy(ny, nakx/2+1)            with all  y and kx >= 0
   !     gky(2*naky-1, nakx/2+1)     with all ky and kx >= 0
   !****************************************************************************
   subroutine transform_y2ky_2d(gy, gky)

      use parallelisation_layouts, only: vmu_lo

      implicit none

      complex, dimension(:, :), intent(in out) :: gy    ! gy(ny, nakx/2+1)
      complex, dimension(:, :), intent(out) :: gky      ! gky(2*naky-1, nakx/2+1)

      integer :: iky_max, ipad_up
      integer :: ikx

      ! Get the indices of the last positive ky <iky_max> and negative ky <ipad_up>
      iky_max = vmu_lo%naky
      ipad_up = iky_max + vmu_lo%ny - (2 * vmu_lo%naky - 1)
      
      ! Iterate over positive kx since we Fourier transform along ny
      do ikx = 1, vmu_lo%nakx / 2 + 1
         ! Fill in the non-zero elements with gy(ny, nakx/2+1)
         fft_y_in = gy(:, ikx)
         ! Apply the "backward" FFT for y(complex) --> ky(complex) along ny
         call dfftw_execute_dft(yb_fft%plan, fft_y_in, fft_y_out)
         ! Rescale the Fourier components with xf_fft%scale=1/ny
         fft_y_out = fft_y_out * yb_fft%scale
         ! Return gky(2*naky-1, nakx/2+1) with x and y in reciprocal space
         ! The positive frequencies are located at the start
         ! The negative frequencies are located at the end
         gky(:iky_max, ikx) = fft_y_out(:iky_max)
         gky(iky_max + 1:, ikx) = fft_y_out(ipad_up + 1:)
      end do

   end subroutine transform_y2ky_2d
   
   !****************************************************************************
   !                       FFT:   KX --C2R--> X   ALONG NX  
   !****************************************************************************
   ! Fourier transformation along nx while the y-axis is already in xyz-space
   !     gkx(ny,nakx/2+1)      with all kx and all y
   !     gx(ny,nx)             with all  x and all y
   !****************************************************************************
   subroutine transform_kx2x(gkx, gx)

      use parallelisation_layouts, only: vmu_lo

      implicit none

      complex, dimension(:, :), intent(in) :: gkx     ! gkx(ny,nakx/2+1)
      real, dimension(:, :), intent(out) :: gx        ! gx(ny,nx)

      integer :: iy

      ! Iterate over ny since we Fourier transform along nx
      do iy = 1, vmu_lo%ny
         ! Pad the remainder of the array with zeros (since the real array is longer
         fft_x_k(vmu_lo%nakx / 2 + 2:) = 0.
         ! Fill in the non-zero elements with gkx(ny,nakx/2+1)
         fft_x_k(:vmu_lo%nakx / 2 + 1) = gkx(iy, :)
         ! Apply the "forward" FFT for kx(complex) --> x(real) along nx
         call dfftw_execute_dft_c2r(xf_fft%plan, fft_x_k, fft_x_x)
         ! Do not rescale the Fourier components since xf_fft%scale=1
         fft_x_x = fft_x_x * xf_fft%scale
         ! Return gx(ny,nx) in real space
         gx(iy, :) = fft_x_x
      end do

   end subroutine transform_kx2x

   !****************************************************************************
   !                    FFT:   X --R2C--> KX   ALONG NX  
   !****************************************************************************
   ! Fourier transformation along nx while the y-axis is already in xyz-space
   !     gx(ny,nx)             with all  x and all y
   !     gkx(ny,nakx/2+1)      with positive kx and all y
   !****************************************************************************
   subroutine transform_x2kx(gx, gkx)

      use parallelisation_layouts, only: vmu_lo

      implicit none

      real, dimension(:, :), intent(in) :: gx            ! gx(ny,nx)
      complex, dimension(:, :), intent(out) :: gkx       ! gkx(ny,nakx/2+1)

      integer :: iy

      ! Iterate over ny since we Fourier transform gx(ny,nx) along nx
      do iy = 1, vmu_lo%ny
         ! Fill in the array with gx(ny,nx)
         fft_x_x = gx(iy, :)
         ! Apply the "backward" FFT for x(real) -->  kx(complex) along nx
         call dfftw_execute_dft_r2c(xb_fft%plan, fft_x_x, fft_x_k)
         ! Rescale the Fourier components with xf_fft%scale=1/nx
         fft_x_k = fft_x_k * xb_fft%scale
         ! Return the positive wavenumbers gkx(ny,nakx/2+1) along x in reciprocal space
         gkx(iy, :) = fft_x_k(:vmu_lo%nakx / 2 + 1)
      end do

   end subroutine transform_x2kx

   !****************************************************************************
   !                    FFT:   KX --C2C--> X   ALONG NX  
   !****************************************************************************
   ! Fourier transformation along nx while the y-axis is still in k-space
   !     gkx(naky,nakx)      with all kx and ky > 0
   !     gx(naky,nx)         with all  x and ky > 0
   !****************************************************************************
   subroutine transform_kx2x_xfirst(gkx, gx)

      use parallelisation_layouts, only: vmu_lo

      implicit none

      complex, dimension(:, :), intent(in) :: gkx     ! gkx(naky,nakx)
      complex, dimension(:, :), intent(out) :: gx     ! gx(naky,nx)

      integer :: iky, ikx_max, ipad_up

      ! Get the indices of the last positive kx <ikx_max> and negative kx <ipad_up>
      ikx_max = vmu_lo%nakx / 2 + 1
      ipad_up = ikx_max + vmu_lo%nx - vmu_lo%nakx

      ! Iterate over naky since we Fourier transform gkx(nakx,nakx) along nx
      do iky = 1, vmu_lo%naky
          ! Pad the middle of the array with zeros: zero padding in the wavenumber
         ! domain results in an increased sampling rate in the real space domain
         fft_xs_k(ikx_max + 1:ipad_up) = 0.
         ! Fill in the non-zero elements with gkx(naky,nakx)
         ! Add the positive frequencies gkx[:vmu_lo%nakx/2+1] to the beginning
         ! And the negative frequencies gkx[vmu_lo%nakx/2+1:] to the end
         fft_xs_k(:ikx_max) = gkx(iky, :ikx_max)
         fft_xs_k(ipad_up + 1:) = gkx(iky, ikx_max + 1:)
         ! Apply the "forward" FFT for kx(complex) --> x(complex) along nx
         call dfftw_execute_dft(xsf_fft%plan, fft_xs_k, fft_xs_x)
         ! Do not rescale the Fourier components since xsf_fft%scale=1
         fft_xs_x = fft_xs_x * xsf_fft%scale
         ! Return gkx(naky,nx) with x in real space
         gx(iky, :) = fft_xs_x
      end do

   end subroutine transform_kx2x_xfirst

   !****************************************************************************
   !                   FFT:   X --C2C--> KX   ALONG NX 
   !****************************************************************************
   ! Fourier transformation along nx while the y-axis is still in k-space
   !     gx(naky,nx)         with all  x and ky > 0
   !     gkx(naky,nakx)      with all kx and ky > 0
   !****************************************************************************
   subroutine transform_x2kx_xfirst(gx, gkx)

      use parallelisation_layouts, only: vmu_lo

      implicit none

      complex, dimension(:, :), intent(in) :: gx            ! gx(ny,nx)
      complex, dimension(:, :), intent(out) :: gkx          ! gx(nakx,nakx)

      integer :: iky, ikx_max, ipad_up

      ! Get the indices of the last positive kx <ikx_max> and negative kx <ipad_up>
      ikx_max = vmu_lo%nakx / 2 + 1
      ipad_up = ikx_max + vmu_lo%nx - vmu_lo%nakx

      ! Iterate over naky since we Fourier transform gx(naky,nx) along nx 
      do iky = 1, vmu_lo%naky
         ! Fill in the array with gx(naky,nx)
         fft_xs_x = gx(iky, :)
         ! Apply the "backward" FFT for x(complex) --> kx(complex) along nx
         call dfftw_execute_dft(xsb_fft%plan, fft_xs_x, fft_xs_k)
         ! Rescale the Fourier components with xsb_fft%scale=1/nx
         fft_xs_k = fft_xs_k * xsb_fft%scale
         ! Return gkx(naky,nakx) with x and y in reciprocal space
         ! The positive frequencies are located at the start
         ! The negative frequencies are located at the end
         gkx(iky, :ikx_max) = fft_xs_k(:ikx_max)
         gkx(iky, ikx_max + 1:) = fft_xs_k(ipad_up + 1:)
      end do

   end subroutine transform_x2kx_xfirst

   !****************************************************************************
   !                   FFT:   KY --C2R--> Y   ALONG NY 
   !****************************************************************************
   ! Fourier transformation along ny while the x-axis is already in xyz-space
   !     gky(naky,nx)         with all x and ky > 0
   !     gy(ny,nx)            with all x and all y
   !****************************************************************************
   subroutine transform_ky2y_xfirst(gky, gy)

      use parallelisation_layouts, only: vmu_lo

      implicit none

      complex, dimension(:, :), intent(in) :: gky
      real, dimension(:, :), intent(out) :: gy

      integer :: ix

      ! Iterate over nx since we Fourier transform gky(naky,nx) along ny
      do ix = 1, vmu_lo%nx
         ! Pad the remainder of the array with zeros (since the real array is longer)
         fft_ys_k(vmu_lo%naky + 1:) = 0.
         ! Fill in the non-zero elements with gky(naky,nx)
         fft_ys_k(:vmu_lo%naky) = gky(:, ix)
         ! Apply the "forward" FFT for ky(complex) --> y(real) along ny
         call dfftw_execute_dft_c2r(ysf_fft%plan, fft_ys_k, fft_ys_y)
         ! Do not rescale the Fourier components since xf_fft%scale=1
         fft_ys_y = fft_ys_y * ysf_fft%scale
         ! Return gx(ny,nx) in real space
         gy(:, ix) = fft_ys_y
      end do

   !****************************************************************************
   !                   FFT:   KY --C2R--> Y   ALONG NY  
   !****************************************************************************
   ! Fourier transformation along ny while the x-axis is already in xyz-space
   !     gy(ny,nx)            with all x and all y
   !     gky(naky,nx)         with all x and ky > 0
   !****************************************************************************
   end subroutine transform_ky2y_xfirst

   subroutine transform_y2ky_xfirst(gy, gky)

      use parallelisation_layouts, only: vmu_lo

      implicit none

      real, dimension(:, :), intent(in) :: gy
      complex, dimension(:, :), intent(out) :: gky

      integer :: ix

      ! Iterate over nx since we Fourier transform gy(ny,nx) along ny
      do ix = 1, vmu_lo%nx
         ! Fill in the array with gy(ny,nx)
         fft_ys_y = gy(:, ix)
         ! Apply the "backward" FFT for y(complex) --> ky(real) along ny
         call dfftw_execute_dft_r2c(ysb_fft%plan, fft_ys_y, fft_ys_k)
         ! Rescale the Fourier components with xf_fft%scale=1/ny
         fft_ys_k = fft_ys_k * ysb_fft%scale
         ! Return the positive wavenumbers gky(naky,nx) along y in reciprocal space
         gky(:, ix) = fft_ys_k(:vmu_lo%naky)
      end do

   end subroutine transform_y2ky_xfirst
   
   !****************************************************************************
   !                   FFT:   KX --C2C--> X   ALONG NAKX 
   !****************************************************************************
   ! Fourier transformation along nakx while the y-axis is still in k-space
   !     gkx(naky,nakx)      with all kx and ky > 0
   !     gx(naky,nx)         with all  x and ky > 0
   !****************************************************************************
   subroutine transform_kx2x_unpadded(gkx, gx)

      implicit none

      complex, dimension(:, :), intent(in)  :: gkx
      complex, dimension(:, :), intent(out) :: gx

      integer :: iy

      ! Iterate over ky since we Fourier transform along nakx
      do iy = 1, size(gkx, 1)
         ! Fill in the array with gkx(naky,nakx)
         fftnp_x_k = gkx(iy, :)
         ! Apply the "forward" FT for kx(complex) --> x(complex) along nakx
         call dfftw_execute_dft(xfnp_fft%plan, fftnp_x_k, fftnp_x_x)
         ! Return gx(naky,nx) with x in real space and rescale with 1/nakx
         gx(iy, :) = fftnp_x_x * xfnp_fft%scale
      end do

   end subroutine transform_kx2x_unpadded

   !****************************************************************************
   !                 FFT:   X --C2C--> KX   ALONG NAKX  
   !****************************************************************************
   ! Fourier transformation along nakx while the y-axis is still in k-space
   !     gx(naky,nx)         with all  x and ky > 0
   !     gkx(naky,nakx)      with all kx and ky > 0
   !****************************************************************************
   subroutine transform_x2kx_unpadded(gx, gkx)

      implicit none

      complex, dimension(:, :), intent(in)  :: gx
      complex, dimension(:, :), intent(out) :: gkx

      integer :: iy

      ! Iterate over ky since we Fourier transform along nakx
      do iy = 1, size(gx, 1)
         ! Fill in the array with gx(naky,nx)
         fftnp_x_x = gx(iy, :)
         ! Apply the "backward" FT for x(complex) --> kx(complex) along nakx
         call dfftw_execute_dft(xbnp_fft%plan, fftnp_x_x, fftnp_x_k)
         ! Return gkx(naky,nakx) in reciprocal space and do not rescale
         gkx(iy, :) = fftnp_x_k * xbnp_fft%scale
      end do

   end subroutine transform_x2kx_unpadded

   !****************************************************************************
   !                  FFT:   KY --C2R--> Y   ALONG NAKY  
   !****************************************************************************
   ! Fourier transformation along naky while the x-axis is already in xyz-space
   !     gky(naky,nx)      with all kx and ky > 0
   !     gy(2*naky-1,nx)   with all  x and y
   !****************************************************************************
   subroutine transform_ky2y_unpadded(gky, gy)

      implicit none

      complex, dimension(:, :), intent(in) :: gky
      real, dimension(:, :), intent(out) :: gy

      integer :: ikx

      ! Iterate over x since we Fourier transform along naky
      do ikx = 1, size(gky, 1)
         ! Fill in the array with gky(naky,nx)
         fftnp_y_k = gky(:, ikx)
         ! Apply the "forward" FT for ky(complex) --> y(real) along naky
         call dfftw_execute_dft_c2r(yfnp_fft%plan, fftnp_y_k, fftnp_y_y)
         ! Return gkx(naky,nakx) in reciprocal space rescale with 1/naky
         gy(:, ikx) = fftnp_y_y * yfnp_fft%scale
      end do

   end subroutine transform_ky2y_unpadded

   !****************************************************************************
   !                    FFT:   Y --C2R--> KY   ALONG NAKY 
   !****************************************************************************
   ! Fourier transformation along naky while the x-axis is already in xyz-space
   !     gy(2*naky-1,nx)   with all  x and y
   !     gky(naky,nx)      with all kx and ky > 0
   !****************************************************************************
   subroutine transform_y2ky_unpadded(gy, gky)

      implicit none

      real, dimension(:, :), intent(in out) :: gy
      complex, dimension(:, :), intent(out) :: gky

      integer :: ikx

      ! Iterate over x since we Fourier transform along naky
      do ikx = 1, size(gy, 2)
         ! Fill in the array with gy(2*naky-1,nx)
         fftnp_y_k = gy(:, ikx)
         ! Apply the "backward" FT for y(complex) --> ky(real) along naky
         call dfftw_execute_dft_r2c(ybnp_fft%plan, fftnp_y_y, fftnp_y_k)
         ! Return gkx(naky,nakx) in reciprocal space do not rescale
         gky(:, ikx) = fftnp_y_y * ybnp_fft%scale
      end do

   end subroutine transform_y2ky_unpadded

   !****************************************************************************
   !****************************************************************************
   !****************************************************************************
   subroutine transform_kalpha2alpha(gkalph, galph)

      use parallelisation_layouts, only: vmu_lo

      implicit none

      complex, dimension(:), intent(in) :: gkalph
      real, dimension(:), intent(out) :: galph

      ! first need to pad input array with zeros
      fft_alpha_kalpha(vmu_lo%naky + 1:) = 0.
      ! then fill in non-zero elements of array
      fft_alpha_kalpha(:vmu_lo%naky) = gkalph
      call dfftw_execute_dft_c2r(alpha_f_fft%plan, fft_alpha_kalpha, fft_alpha_alpha)
      fft_alpha_alpha = fft_alpha_alpha * alpha_f_fft%scale
      galph = fft_alpha_alpha

   end subroutine transform_kalpha2alpha

   !****************************************************************************
   !****************************************************************************
   !****************************************************************************
   ! input galph array is real and contains values on the padded alpha grid
   ! gkalph is output array; it contains the Fourier coefficients of galph
   ! for positive ky values only (reality can be used to obtain the negative ky coefs)
   ! the highest 1/3 of the ky modes from the FFT have been discarded to avoid de-aliasing
   subroutine transform_alpha2kalpha(galph, gkalph)

      use parallelisation_layouts, only: vmu_lo

      implicit none

      real, dimension(:), intent(in) :: galph
      complex, dimension(:), intent(out) :: gkalph

      fft_alpha_alpha = galph
      call dfftw_execute_dft_r2c(alpha_b_fft%plan, fft_alpha_alpha, fft_alpha_kalpha)
      fft_alpha_kalpha = fft_alpha_kalpha * alpha_b_fft%scale
      ! filter out highest k-alpha modes to avoid aliasing
      gkalph = fft_alpha_kalpha(:vmu_lo%naky)

   end subroutine transform_alpha2kalpha

!###############################################################################
!                                 FINISH
!###############################################################################

   !****************************************************************************
   !                      FINISH THE FOURIER TRANSFORMS 
   !****************************************************************************
   ! Destroy the FFT plans and deallocate the arrays
   subroutine finish_transforms

      use parameters_physics, only: full_flux_surface

      implicit none

      call dfftw_destroy_plan(yf_fft%plan)
      call dfftw_destroy_plan(yb_fft%plan)
      call dfftw_destroy_plan(xf_fft%plan)
      call dfftw_destroy_plan(xb_fft%plan)
      call dfftw_destroy_plan(xsf_fft%plan)
      call dfftw_destroy_plan(xsb_fft%plan)
      call dfftw_destroy_plan(ysf_fft%plan)
      call dfftw_destroy_plan(ysb_fft%plan)
      call dfftw_destroy_plan(yfnp_fft%plan)
      call dfftw_destroy_plan(ybnp_fft%plan)
      call dfftw_destroy_plan(xfnp_fft%plan)
      call dfftw_destroy_plan(xbnp_fft%plan)
      if (full_flux_surface) then
         call dfftw_destroy_plan(alpha_f_fft%plan)
         call dfftw_destroy_plan(alpha_b_fft%plan)
      end if
      if (allocated(fft_y_in)) deallocate (fft_y_in)
      if (allocated(fft_y_out)) deallocate (fft_y_out)
      if (allocated(fft_x_k)) deallocate (fft_x_k)
      if (allocated(fft_x_x)) deallocate (fft_x_x)
      if (allocated(fft_xs_k)) deallocate (fft_xs_k)
      if (allocated(fft_xs_x)) deallocate (fft_xs_x)
      if (allocated(fft_ys_k)) deallocate (fft_ys_k)
      if (allocated(fft_ys_y)) deallocate (fft_ys_y)
      if (allocated(fft_alpha_alpha)) deallocate (fft_alpha_alpha)
      if (allocated(fft_alpha_kalpha)) deallocate (fft_alpha_kalpha)
      if (allocated(fftnp_y_k)) deallocate (fftnp_y_k)
      if (allocated(fftnp_y_y)) deallocate (fftnp_y_y)
      if (allocated(fftnp_x_k)) deallocate (fftnp_x_k)
      if (allocated(fftnp_x_x)) deallocate (fftnp_x_x)
      transforms_initialized = .false.

   end subroutine finish_transforms

end module calculations_transforms
