module multibox

  use fft_work, only: fft_type

  implicit none

  public :: init_multibox
  public :: finish_multibox
  public :: multibox_communicate
  public :: boundary_size
  public :: shear_rate

  private

  complex, dimension (:), allocatable :: g_buffer0
  complex, dimension (:), allocatable :: g_buffer1
  
  complex, dimension (:,:), allocatable :: fft_xky
  real, dimension (:,:), allocatable :: fft_xy

  ! for the unpadded FFTs
  type (fft_type) :: yf_fft, yb_fft
  type (fft_type) :: xf_fft, xb_fft

  complex, dimension (:), allocatable :: fft_x_k, fft_x_x
  complex, dimension (:), allocatable :: fft_y_k
  real, dimension (:), allocatable :: fft_y_y

  real :: shear_rate

  integer :: g_buff_size
  integer :: phi_buff_size
  logical :: mb_transforms_initialized = .false.
  logical :: smooth_ZFs


  integer :: temp_ind = 0

!>DSO 
! the multibox simulation has one parameter: boundary_size
! 
  integer :: boundary_size

contains

  subroutine init_multibox
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky
    use file_utils, only: input_unit_exist
    use file_utils, only: runtype_option_switch, runtype_multibox
    use mp, only: broadcast, proc0
    use fft_work, only: init_ccfftw, delete_fft,FFT_FORWARD, fft_type
    use linear_solve, only: lu_decomposition

    implicit none

#ifndef MPI
    return
#else

    integer :: in_file
    logical exist

    namelist /multibox_parameters/ boundary_size, shear_rate, smooth_ZFs

    if(runtype_option_switch /= runtype_multibox) return

    boundary_size = 4
    shear_rate = 0
    smooth_ZFs = .true.

    if (proc0) then
      in_file = input_unit_exist("multibox_parameters", exist)
      if (exist) read (in_file, nml=multibox_parameters)
    endif

    call broadcast(boundary_size)

    phi_buff_size = boundary_size*naky*ntubes*(2*nzgrid+1)
    g_buff_size   = phi_buff_size*(vmu_lo%ulim_alloc-vmu_lo%llim_proc+1)

    if (.not.allocated(g_buffer0)) allocate(g_buffer0(g_buff_size))
    if (.not.allocated(g_buffer1)) allocate(g_buffer1(g_buff_size))

    call init_mb_transforms

  
#endif
  end subroutine init_multibox

  subroutine finish_multibox

    implicit none

    if (allocated(g_buffer0))   deallocate (g_buffer0)
    if (allocated(g_buffer1))   deallocate (g_buffer1)

    call finish_mb_transforms

  end subroutine finish_multibox


  subroutine multibox_communicate

    use dist_fn_arrays, only: gnew
    use kt_grids, only: nakx,naky,naky_all
    use file_utils, only: runtype_option_switch, runtype_multibox
    use file_utils, only: get_unused_unit
    use fields_arrays, only: phi, apar
    use fields, only: advance_fields, fields_updated
    use job_manage, only: njobs
    use linear_solve, only: lu_back_substitution
    use stella_layouts, only: vmu_lo
    use mp, only: job, scope, mp_abort,  &
                  crossdomprocs, subprocs, &
                  send, receive, proc0

    implicit none

    integer :: num,ix,iy,iz,it,iv,offset
    integer :: ii,jj, temp_unit
    complex :: dzm,dzp
    character(len=1024) :: filename

#ifndef MPI
    return
#else
    if(runtype_option_switch /= runtype_multibox) return
    if(njobs /= 3) call mp_abort("Multibox only supports 3 domains at the moment.")


    if(mod(temp_ind,50)==0 .and. proc0) then
     ! call get_unused_unit(temp_unit)
      temp_unit=3023+job
      call transform_kx2x(phi(:,:,0,1),fft_xky)  
      call transform_ky2y(fft_xky,fft_xy)
      write (filename,"(A,I1,A,I0.6)") "phiout",job,"_",temp_ind
      open (unit=temp_unit, file=filename, status="replace",action="write")
      do ii=1,naky_all
        do jj=1,nakx
          write (temp_unit,*) jj,ii,fft_xy(ii,jj)
        enddo
        write (temp_unit,*) ""
      enddo
      close (unit=temp_unit)
    endif

! DSO - change communicator
    call scope(crossdomprocs)

    if(job==0 .or. job==(njobs-1)) then
      offset=0;
      if(job==0) offset=nakx-boundary_size
      num=1
      do iv = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, vmu_lo%ntubes
          do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
            call transform_kx2x(gnew(:,:,iz,it,iv),fft_xky)  
            do ix=1,boundary_size
              do iy=1,naky
                !if in the future the grids can have different naky, one will
                !have to divide by naky here, and multiply on the receiving end
                g_buffer0(num) = fft_xky(iy,ix + offset)
                num=num+1
              enddo
            enddo
          enddo
        enddo
      enddo
! DSO - send data
      call send(g_buffer0,1,43 + job)
    else
      offset = nakx - boundary_size
! DSO - receive the data

! left
      call receive(g_buffer0,0, 43)
! right
      call receive(g_buffer1,njobs-1, 43+njobs-1)

      num=1
      do iv = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, vmu_lo%ntubes
          do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
            call transform_kx2x(gnew(:,:,iz,it,iv),fft_xky)  
            do ix=1,boundary_size
              do iy=1,naky
                fft_xky(iy,ix)        = g_buffer0(num)
                fft_xky(iy,ix+offset) = g_buffer1(num)
                num=num+1
              enddo
            enddo
            if(smooth_ZFs) then
              dzm = fft_xky(1,boundary_size+1)      - fft_xky(1,boundary_size)
              dzp = fft_xky(1,nakx-boundary_size+1) - fft_xky(1,nakx-boundary_size)
              do ix=1,boundary_size
                fft_xky(1,ix)=fft_xky(1,ix) + dzm
                fft_xky(1,ix+offset)=fft_xky(1,ix+offset) - dzp
              enddo
            endif
            call transform_x2kx(fft_xky,gnew(:,:,iz,it,iv))  
          enddo
        enddo
      enddo
    endif

! DSO - change communicator
    call scope(subprocs)
    if(job==1) then
      fields_updated = .false.
      call advance_fields(gnew,phi,apar, dist='gbar')
    endif

    temp_ind=temp_ind+1

#endif
  end subroutine multibox_communicate


!!>DSO - The following subroutines are the _unpadded_ analogues of the ones found in
! stella_transforms.f90. These perhaps should also be moved to that file, but for now I
! want to keep all the multibox subroutines in one place.

  subroutine init_mb_transforms

    use stella_layouts, only: init_stella_layouts
    use kt_grids, only: nakx, naky, naky_all

    implicit none

    if (mb_transforms_initialized) return
    mb_transforms_initialized = .true.

    if (.not.allocated(fft_xky)) allocate (fft_xky(naky,nakx))
    if (.not.allocated(fft_xy)) allocate (fft_xy(naky_all,nakx))

    call init_x_fft
    call init_y_fft

  end subroutine init_mb_transforms

  subroutine init_x_fft

    use kt_grids, only: nakx
    use fft_work, only: init_ccfftw, FFT_BACKWARD, FFT_FORWARD

    implicit none

    if (.not.allocated(fft_x_k)) allocate (fft_x_k(nakx))
    if (.not.allocated(fft_x_x)) allocate (fft_x_x(nakx))

    call init_ccfftw (xf_fft, FFT_BACKWARD, nakx, fft_x_k, fft_x_x)
    call init_ccfftw (xb_fft, FFT_FORWARD , nakx, fft_x_x, fft_x_k)

  end subroutine init_x_fft

  subroutine init_y_fft

    use kt_grids, only: naky, naky_all
    use fft_work, only: init_crfftw, init_rcfftw, FFT_BACKWARD, FFT_FORWARD

    implicit none

    if (.not.allocated(fft_y_k)) allocate (fft_y_k(naky))
    if (.not.allocated(fft_y_y)) allocate (fft_y_y(naky_all))

    call init_crfftw (yf_fft, FFT_BACKWARD, naky_all, fft_y_k, fft_y_y)
    call init_rcfftw (yb_fft, FFT_FORWARD , naky_all, fft_y_y, fft_y_k)

  end subroutine init_y_fft


!  
!> transform routines start here
!

  subroutine transform_kx2x (gkx, gx)

    use kt_grids, only: naky

    implicit none

    complex, dimension (:,:), intent (in)  :: gkx
    complex, dimension (:,:), intent (out) :: gx

    integer :: iy

    do iy = 1, naky
       fft_x_k = gkx(iy,:)
       call dfftw_execute_dft(xf_fft%plan, fft_x_k, fft_x_x)
       gx(iy,:) = fft_x_x*xf_fft%scale
    end do

  end subroutine transform_kx2x

  subroutine transform_x2kx (gx, gkx)

    use kt_grids, only: naky

    implicit none

    complex, dimension (:,:), intent (in)  :: gx
    complex, dimension (:,:), intent (out) :: gkx

    integer :: iy

    do iy = 1, naky
       fft_x_x = gx(iy,:)
       call dfftw_execute_dft(xb_fft%plan, fft_x_x, fft_x_k)
       gkx(iy,:) = fft_x_k*xb_fft%scale
    end do

  end subroutine transform_x2kx

  subroutine transform_ky2y (gky, gy)

    use kt_grids, only: nakx

    implicit none

    complex, dimension (:,:), intent (in) :: gky
    real, dimension (:,:), intent (out) :: gy

    integer :: ikx

    do ikx = 1, nakx
       fft_y_k = gky(:,ikx)
       call dfftw_execute_dft_c2r(yf_fft%plan, fft_y_k, fft_y_y)
       gy(:,ikx) =fft_y_y*yf_fft%scale
    end do

  end subroutine transform_ky2y

  subroutine transform_y2ky (gy, gky)

    use kt_grids, only: nakx

    implicit none

    real, dimension (:,:), intent (in out) :: gy
    complex, dimension (:,:), intent (out) :: gky

    integer :: ikx

    do ikx = 1, nakx
       fft_y_k = gy(:,ikx)
       call dfftw_execute_dft_r2c(yb_fft%plan, fft_y_k, fft_y_y)
       gky(:,ikx) = fft_y_y*yb_fft%scale
    end do

  end subroutine transform_y2ky


  subroutine finish_mb_transforms

    implicit none

    call dfftw_destroy_plan (yf_fft%plan)
    call dfftw_destroy_plan (yb_fft%plan)
    call dfftw_destroy_plan (xf_fft%plan)
    call dfftw_destroy_plan (xb_fft%plan)
    if (allocated(fft_y_k)) deallocate (fft_y_k)
    if (allocated(fft_y_y)) deallocate (fft_y_y)
    if (allocated(fft_x_k)) deallocate (fft_x_k)
    if (allocated(fft_x_x)) deallocate (fft_x_x)
    if (allocated(fft_xky)) deallocate (fft_xky)
    if (allocated(fft_xy)) deallocate (fft_xy)

    mb_transforms_initialized = .false.

  end subroutine finish_mb_transforms

end module multibox
