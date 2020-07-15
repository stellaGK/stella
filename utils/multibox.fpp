module multibox

  use fft_work, only: fft_type
  use common_types, only: flux_surface_type

  implicit none

  public :: read_multibox_parameters
  public :: init_multibox
  public :: finish_multibox
  public :: multibox_communicate
  public :: add_multibox_krook
  public :: boundary_size
  public :: bs_fullgrid
  public :: g_exb, shear
  public :: xL, xR
  public :: kx0_L, kx0_R
  public :: RK_step
  public :: include_multibox_krook

  private

  complex, dimension (:), allocatable :: g_buffer0
  complex, dimension (:), allocatable :: g_buffer1
  complex, dimension (:), allocatable :: fsa_x
  real,    dimension (:), allocatable :: shear
  real,    dimension (:), allocatable :: copy_mask_left, copy_mask_right
  real,    dimension (:), allocatable :: krook_mask_left, krook_mask_right
  
  complex, dimension (:,:), allocatable :: fft_xky
  real, dimension (:,:), allocatable :: fft_xy

  ! for the unpadded FFTs
  type (fft_type) :: yf_fft, yb_fft
  type (fft_type) :: xf_fft, xb_fft

  complex, dimension (:), allocatable :: fft_x_k, fft_x_x
  complex, dimension (:), allocatable :: fft_y_k
  real, dimension (:), allocatable :: fft_y_y

  logical :: mb_transforms_initialized = .false.
  integer :: temp_ind = 0
  integer :: bs_fullgrid
  integer :: mb_debug_step
  

  real :: xL = 0., xR = 0.
  real :: kx0_L, kx0_R

!>DSO 
! the multibox simulation has one parameter: boundary_size
! 
  integer :: boundary_size, krook_size
  real :: g_exb, nu_krook_mb
  logical :: smooth_ZFs
  logical :: RK_step, include_multibox_krook
  integer :: krook_option_switch
  integer, parameter:: krook_option_default = 0, &
                       krook_option_linear  = 1, &
                       krook_option_expin   = 2, &
                       krook_option_expout  = 3 
  integer:: mb_zf_option_switch
  integer, parameter :: mb_zf_option_default = 0, &
                        mb_zf_option_no_ky0  = 1, &
                        mb_zf_option_no_fsa  = 2
  integer :: shear_option_switch
  integer, parameter:: shear_option_triangle  = 0, &
                       shear_option_sine      = 1, &
                       shear_option_flat      = 2 

contains

  subroutine read_multibox_parameters

    use file_utils, only: input_unit_exist, error_unit
    use file_utils, only: runtype_option_switch, runtype_multibox
    use text_options, only: text_option, get_option_value
    use mp, only: broadcast, proc0

    implicit none

    integer :: in_file, ierr
    logical exist
    

    type (text_option), dimension (4), parameter :: krook_opts = &
      (/ text_option('default', krook_option_default), &
         text_option('linear',  krook_option_linear) , &
         text_option('exp_in',  krook_option_expin) , &
         text_option('exp_out', krook_option_expout)/)
    type (text_option), dimension (3), parameter :: mb_zf_opts = &
      (/ text_option('default', mb_zf_option_default), &
         text_option('no_ky0',  mb_zf_option_no_ky0) , &
         text_option('no_fsa',  mb_zf_option_no_fsa)/)
    type (text_option), dimension (4), parameter :: shear_opts = &
      (/ text_option('default',  shear_option_triangle), &
         text_option('triangle', shear_option_triangle) , &
         text_option('sine',     shear_option_sine), &
         text_option('flat',     shear_option_flat)/)
    character(30) :: zf_option, krook_option, shear_option

    namelist /multibox_parameters/ boundary_size, krook_size, shear_option,& 
                                   g_exb, smooth_ZFs, zf_option, &
                                   krook_option, RK_step, nu_krook_mb, &
                                   mb_debug_step

    if(runtype_option_switch /= runtype_multibox) return

    boundary_size = 4
    krook_size = 0
    nu_krook_mb = 0.0
    g_exb = 0.
    mb_debug_step = 1000
    smooth_ZFs = .false.
    RK_step = .false.
    zf_option = 'default'
    shear_option = 'default'
    krook_option = 'default'
    
    if (proc0) then
      in_file = input_unit_exist("multibox_parameters", exist)
      if (exist) read (in_file, nml=multibox_parameters)

      ierr = error_unit()
      call get_option_value & 
        (krook_option, krook_opts, krook_option_switch, & 
         ierr, "krook_option in multibox_parameters")
      call get_option_value & 
        (zf_option, mb_zf_opts, mb_zf_option_switch, & 
         ierr, "zf_option in multibox_parameters")
      call get_option_value & 
        (shear_option, shear_opts, shear_option_switch, & 
         ierr, "shear_option in multibox_parameters")

       if(krook_size > boundary_size) krook_size = boundary_size 
    endif


    call broadcast(boundary_size)
    call broadcast(krook_size)
    call broadcast(nu_krook_mb)
    call broadcast(g_exb)
    call broadcast(smooth_ZFs)
    call broadcast(mb_zf_option_switch)
    call broadcast(krook_option_switch)
    call broadcast(shear_option_switch)
    call broadcast(RK_step)
    call broadcast(mb_debug_step)

    if(krook_option_switch.ne.krook_option_default) include_multibox_krook = .true.


  end subroutine read_multibox_parameters

  subroutine init_multibox 
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx,naky,akx, nx,x, x_clamped
    use file_utils, only: runtype_option_switch, runtype_multibox
    use job_manage, only: njobs
    use mp, only: scope, crossdomprocs, subprocs, &
                  send, receive, job

    implicit none

#ifndef MPI

    allocate (shear(1)); shear(1)=0.0
    return
#else

    integer :: g_buff_size
    integer :: phi_buff_size
    integer :: i

    real :: db

    if(runtype_option_switch /= runtype_multibox) return

    bs_fullgrid = nint((3.0*boundary_size)/2.0)

    phi_buff_size = boundary_size*naky*ntubes*(2*nzgrid+1)
    g_buff_size   = phi_buff_size*(vmu_lo%ulim_alloc-vmu_lo%llim_proc+1)

    if (.not.allocated(g_buffer0)) allocate(g_buffer0(g_buff_size))
    if (.not.allocated(g_buffer1)) allocate(g_buffer1(g_buff_size))
    if (.not.allocated(fsa_x) .and. (mb_zf_option_switch.eq.mb_zf_option_no_fsa)) then
      allocate(fsa_x(nakx))
    endif
    if (.not.allocated(copy_mask_left))  allocate(copy_mask_left(boundary_size));  copy_mask_left =1.0
    if (.not.allocated(copy_mask_right)) allocate(copy_mask_right(boundary_size)); copy_mask_right=1.0
    if (.not.allocated(krook_mask_left))  allocate(krook_mask_left(boundary_size));  krook_mask_left =0.0
    if (.not.allocated(krook_mask_right)) allocate(krook_mask_right(boundary_size)); krook_mask_right=0.0

    select case (krook_option_switch)
    case (krook_option_default)
      ! do nothing
    case (krook_option_linear)
      db = 1.0/krook_size
      do i = 1, krook_size
        krook_mask_right(i) = i*db
        copy_mask_right(i) = 0.0
      enddo
    case (krook_option_expin)
      db = 3.0/krook_size
      do i = 1, krook_size
        krook_mask_right(i) = 1.0-(1.0-exp(-(krook_size-i)*db))/(1.0-exp(-3.0))
        copy_mask_right(i) = 0.0
      enddo
    case (krook_option_expout)
      db = 3.0/krook_size
      do i = 1, krook_size
        krook_mask_right(i) = (1.0-exp(-i*db))/(1.0-exp(-3.0))
        copy_mask_right(i) = 0.0
      enddo
    end select

    do i = 1, boundary_size
      copy_mask_left(i)  = copy_mask_right(boundary_size - i + 1)
      krook_mask_left(i) = krook_mask_right(boundary_size - i + 1)
    enddo

    call scope(crossdomprocs)

    if(job==1) then
      xL = x(bs_fullgrid)
      xR = x(nx-bs_fullgrid+1)
      call send(xL,0)
      call send(xR,njobs-1)

      if(.not.allocated(x_clamped)) allocate(x_clamped(nx))
      x_clamped = x
      do i = 1, nx
        if(x_clamped(i) < xL) x_clamped(i) = xL
        if(x_clamped(i) > xR) x_clamped(i) = xR
      enddo

      call receive(kx0_L,0)
      call receive(kx0_R,njobs-1)
    else if(job==0) then
      call receive(xL,1)
      call send(akx(2),1)
    elseif(job==njobs-1) then
      call receive(xR,1)
      call send(akx(2),1)
    endif

    call scope(subprocs)

    call init_shear

    call init_mb_transforms
  
#endif
  end subroutine init_multibox

  subroutine finish_multibox

    use kt_grids, only: x_clamped

    implicit none

    if (allocated(g_buffer0))   deallocate (g_buffer0)
    if (allocated(g_buffer1))   deallocate (g_buffer1)
    if (allocated(fsa_x))       deallocate (fsa_x)
    if (allocated(x_clamped))   deallocate (x_clamped)
    if (allocated(shear))       deallocate (shear)

    if (allocated(copy_mask_left))   deallocate (copy_mask_left)
    if (allocated(copy_mask_right))  deallocate (copy_mask_right)
    if (allocated(krook_mask_left))  deallocate (krook_mask_left)
    if (allocated(krook_mask_right)) deallocate (krook_mask_right)

    call finish_mb_transforms

  end subroutine finish_multibox

  subroutine init_shear

    use kt_grids, only: akx, nx, x, dx
    use mp, only: job

    implicit none

    integer :: i, j, ccount

    if (.not.allocated(shear)) allocate(shear(nx)); shear = 0.0

    if(g_exb*g_exb > epsilon(0.0))  then
      select case (job)

      case (0)
        do i=1,nx
          select case (shear_option_switch)

          case (shear_option_triangle)
            j=i-bs_fullgrid-nx/2
            do while(j<1) 
              j = j + nx
            enddo
            shear(j) = g_exb*(xL + dx*(abs(i-nx/2) - abs(2*bs_fullgrid)))
          case (shear_option_sine)
            shear(i) = g_exb*(xL + sin(akx(2)*dx*(i-bs_fullgrid))/akx(2))
          case (shear_option_flat)
            shear(i) = g_exb*xL
          end select
        enddo
      case (1)
        select case (shear_option_switch)
          
        case (shear_option_triangle)
          do i=1,nx
            shear(i) = g_exb*x(i)
          enddo
        case (shear_option_sine)
          do i=1,nx
            shear(i) = g_exb*x(i)
          enddo
          do i = 1, bs_fullgrid
            shear(i)                = g_exb*(xL + sin(kx0_L*dx*(i-bs_fullgrid))/kx0_L)
            shear(nx-bs_fullgrid+i) = g_exb*(xR + sin(kx0_R*dx*(i-1))/kx0_R)
          enddo
        case (shear_option_flat)
          ccount = nx - 2*bs_fullgrid
          do i=1,ccount
            shear(i+bs_fullgrid) = g_exb*x(i+bs_fullgrid)
          enddo
          shear(1:bs_fullgrid)         = g_exb*xL
          shear((nx-bs_fullgrid+1):nx) = g_exb*xR
        end select
      case (2)
        do i=1,nx
          select case (shear_option_switch)
          
          case (shear_option_triangle)
            j=i+bs_fullgrid
            do while(j>nx) 
              j = j - nx
            enddo
            shear(j) = g_exb*(xR + dx*(abs(i-nx/2) - abs(nx/2 - 2*bs_fullgrid + 1)))
          case (shear_option_sine)
            shear(i) = g_exb*(xR + sin(akx(2)*dx*(i-1+bs_fullgrid))/akx(2))
          case (shear_option_flat)
            shear(i) = g_exb*xR
          end select
        enddo
      end select
    endif

    do i=1,nx
      write (*,*) i, shear(i), job
    enddo

  end subroutine init_shear

  subroutine multibox_communicate (gin)

    use kt_grids, only: nakx,naky,naky_all,nx,ny,dx,dy, zonal_mode
    use file_utils, only: runtype_option_switch, runtype_multibox
    use file_utils, only: get_unused_unit
    use fields_arrays, only: phi
    use fields, only: advance_fields, fields_updated
    use job_manage, only: njobs
    use stella_layouts, only: vmu_lo
    use stella_geometry, only: dl_over_b
    use zgrid, only: nzgrid
    use mp, only: job, scope, mp_abort,  &
                  crossdomprocs, subprocs, allprocs, &
                  send, receive, proc0

    implicit none

    integer :: num,ia, ix,iky,iz,it,iv,offset
    integer :: ii,jj, temp_unit
    real :: afacx, afacy
    complex :: dzm,dzp
    character(len=512) :: filename

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (inout) :: gin

#ifndef MPI
    return
#else
    if(runtype_option_switch /= runtype_multibox) return
    if(njobs /= 3) call mp_abort("Multibox only supports 3 domains at the moment.")


    if(mod(temp_ind,mb_debug_step)==0 .and. proc0) then
     ! call get_unused_unit(temp_unit)
      temp_unit=3023+job
      afacx = real(nx)/real(nakx)
      afacy = real(ny)/real(2*naky-1)
      call transform_kx2x(phi(:,:,0,1),fft_xky)  
      call transform_ky2y(fft_xky,fft_xy)
      write (filename,"(A,I1,A,I0.6)") "phiout",job,"_",temp_ind
      open (unit=temp_unit, file=filename, status="replace",& 
            action="write",form="unformatted",access="stream")
      write (temp_unit) real(nakx,4)
      do ii=1,nakx
        write(temp_unit) real(afacx*dx*(ii-1),4)
      enddo
      do ii=1,naky_all
        write (temp_unit) real(afacy*dy*(ii-1),4)
        do jj=1,nakx
          write (temp_unit) real(fft_xy(ii,jj),4)
        enddo
      enddo
      close (unit=temp_unit)
    endif

! DSO - change communicator
    call scope(crossdomprocs)

    ia=1

    if(job==0 .or. job==(njobs-1)) then
      offset=0;
      ! DSO the next line might seem backwards, but this makes it easier to stitch together imaages
      ! FLAG DSO - might do something weird with magnetic shear
      if(job==njobs-1) offset=nakx-boundary_size
      num=1
      do iv = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, vmu_lo%ntubes

          !this is where the FSA goes
          if(mb_zf_option_switch .eq. mb_zf_option_no_fsa) then
            do ix= 1,nakx
              fft_x_k(ix) = sum(dl_over_b(ia,:)*gin(1,ix,:,it,iv))
            enddo
            call dfftw_execute_dft(xf_fft%plan, fft_x_k, fft_x_x)
            fsa_x = fft_x_x*xf_fft%scale
          endif

          do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
            call transform_kx2x(gin(:,:,iz,it,iv),fft_xky)  
            do ix=1,boundary_size
              do iky=1,naky
                !DSO if in the future the grids can have different naky, one will
                !have to divide by naky here, and multiply on the receiving end
                g_buffer0(num) = fft_xky(iky,ix + offset)
                if((iky.eq. 1) .and. (zonal_mode(iky))) then
                  if(    mb_zf_option_switch .eq. mb_zf_option_no_ky0) then
                    g_buffer0(num) = 0
                  elseif(mb_zf_option_switch .eq. mb_zf_option_no_fsa) then
                    g_buffer0(num) = fft_xky(iky,ix + offset) - fsa_x(ix + offset)
                  endif
                endif
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
            call transform_kx2x(gin(:,:,iz,it,iv),fft_xky)  
            do ix=1,boundary_size
              do iky=1,naky
                fft_xky(iky,ix)        = fft_xky(iky,ix)*(1-copy_mask_left(ix)) &
                                       + g_buffer0(num)*copy_mask_left(ix)
                                        
                fft_xky(iky,ix+offset) = fft_xky(iky,ix+offset)*(1-copy_mask_right(ix)) &
                                       + g_buffer1(num)*copy_mask_right(ix)
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
            call transform_x2kx(fft_xky,gin(:,:,iz,it,iv))  
          enddo
        enddo
      enddo
    endif

! DSO - change communicator
    call scope(subprocs)

    if(job==1) fields_updated = .false.
    !if(job==1) then
     ! fields_updated = .false.
      !call advance_fields(gnew,phi,apar, dist='gbar')
    !endif

    temp_ind=temp_ind+1

#endif
  end subroutine multibox_communicate

  subroutine add_multibox_krook (g, rhs)

    use stella_time, only: code_dt
    use stella_layouts, only: vmu_lo
    use kt_grids, only:  nakx, naky
    use zgrid, only: nzgrid, ntubes
    use mp, only: job


    implicit none

    integer ::  iky, ix, iz, it, ivmu, num, offset

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: rhs

    complex, allocatable, dimension (:,:) :: g0x, g0k
    if(job /= 1) return

    allocate (g0k(naky,nakx))
    allocate (g0x(naky,nakx))
    
    offset = nakx - boundary_size

    num=1
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      do it = 1, ntubes
        do iz = -nzgrid, nzgrid
          g0x = 0.0
          call transform_kx2x(g(:,:,iz,it,ivmu),fft_xky)  
          do ix=1,boundary_size
            do iky=1,naky
                g0x(iky,ix)        = (fft_xky(iky,ix)        - g_buffer0(num))*krook_mask_left(ix)
                g0x(iky,ix+offset) = (fft_xky(iky,ix+offset) - g_buffer1(num))*krook_mask_right(ix)
                num=num+1
            enddo
          enddo
          call transform_x2kx(g0x,g0k)  
          rhs(:,:,iz,it,ivmu) = rhs(:,:,iz,it,ivmu) - code_dt*nu_krook_mb*g0k
        enddo
      enddo
    enddo

    deallocate(g0k,g0x)

  end subroutine add_multibox_krook



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

!   subroutine transform_y2ky (gy, gky)
!
!    use kt_grids, only: nakx
!
!    implicit none
!
!    real, dimension (:,:), intent (in out) :: gy
!    complex, dimension (:,:), intent (out) :: gky
!
!    integer :: ikx
!
!    do ikx = 1, nakx
!       fft_y_k = gy(:,ikx)
!       call dfftw_execute_dft_r2c(yb_fft%plan, fft_y_y, fft_y_k)
!       gky(:,ikx) = fft_y_y*yb_fft%scale
!    end do
!
!  end subroutine transform_y2ky


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
