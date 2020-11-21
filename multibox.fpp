module multibox

  use fft_work, only: fft_type

  implicit none

  public :: read_multibox_parameters
  public :: init_multibox
  public :: finish_multibox
  public :: multibox_communicate
  public :: communicate_multibox_parameters
  public :: add_multibox_krook
  public :: boundary_size
  public :: bs_fullgrid
  public :: xL, xR, xL_d, xR_d
  public :: rhoL, rhoR
  public :: kx0_L, kx0_R
  public :: RK_step, comm_at_init
  public :: include_multibox_krook


  private

  complex, dimension (:), allocatable :: g_buffer0
  complex, dimension (:), allocatable :: g_buffer1
  complex, dimension (:), allocatable :: fsa_x
  real,    dimension (:), allocatable :: copy_mask_left, copy_mask_right
  real,    dimension (:), allocatable :: krook_mask_left, krook_mask_right
  real,    dimension (:), allocatable :: krook_fac
  real, dimension(:), allocatable :: x_mb

  real :: dx_mb

  complex, dimension (:,:), allocatable :: fft_kxky, fft_xky
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
  integer :: x_fft_size
  

  real :: xL = 0., xR = 0.
  real :: rhoL = 0., rhoR = 0.
  real :: xL_d = 0., xR_d = 0.
  real :: kx0_L, kx0_R

  integer :: boundary_size, krook_size
  real :: nu_krook_mb, krook_exponent
  logical :: smooth_ZFs
  logical :: RK_step, include_multibox_krook, comm_at_init
  integer :: krook_option_switch
  integer, parameter:: krook_option_default = 0, &
                       krook_option_linear  = 0, &
                       krook_option_exp   = 1, &
                       krook_option_exp_rev  = 2 
  integer:: mb_zf_option_switch
  integer, parameter :: mb_zf_option_default = 0, &
                        mb_zf_option_no_ky0  = 1, &
                        mb_zf_option_no_fsa  = 2
  integer :: LR_debug_switch
  integer, parameter:: LR_debug_option_default  = 0, &
                       LR_debug_option_L        = 1, &
                       LR_debug_option_R        = 2 

contains

  subroutine read_multibox_parameters

    use file_utils, only: input_unit_exist, error_unit
    use file_utils, only: runtype_option_switch, runtype_multibox
    use text_options, only: text_option, get_option_value
    use mp, only: broadcast, proc0
    use kt_grids, only: nx, nakx
    use job_manage, only: njobs
    use mp, only: scope, crossdomprocs, subprocs, &
                  send, receive, job

    implicit none

    integer :: in_file, ierr
    integer :: nakxl, nxl, nakxr, nxr, fac
    logical exist
    

    type (text_option), dimension (4), parameter :: krook_opts = &
      (/ text_option('default', krook_option_default), &
         text_option('linear',  krook_option_linear) , &
         text_option('exp',         krook_option_exp) , &
         text_option('exp_reverse', krook_option_exp_rev)/)
    type (text_option), dimension (3), parameter :: mb_zf_opts = &
      (/ text_option('default', mb_zf_option_default), &
         text_option('no_ky0',  mb_zf_option_no_ky0) , &
         text_option('no_fsa',  mb_zf_option_no_fsa)/)
    type (text_option), dimension (3), parameter :: LR_db_opts = &
      (/ text_option('default', LR_debug_option_default), &
         text_option('L',       LR_debug_option_L) , &
         text_option('R',       LR_debug_option_R)/)
    character(30) :: zf_option, krook_option, LR_debug_option

    namelist /multibox_parameters/ boundary_size, krook_size, & 
                                   smooth_ZFs, zf_option, LR_debug_option, &
                                   krook_option, RK_step, nu_krook_mb, &
                                   mb_debug_step, krook_exponent, comm_at_init

    if(runtype_option_switch /= runtype_multibox) return

    boundary_size = 4
    krook_size = 0
    krook_exponent = 0.0
    nu_krook_mb = 0.0
    mb_debug_step = 1000
    smooth_ZFs = .false.
    comm_at_init = .false.
    RK_step = .false.
    zf_option = 'default'
    krook_option = 'default'
    LR_debug_option = 'default'
    
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
        (LR_debug_option, LR_db_opts, LR_debug_switch, & 
         ierr, "LR_debug_option in multibox_parameters")

       if(krook_size > boundary_size) krook_size = boundary_size 
    endif

    call broadcast(boundary_size)
    call broadcast(krook_size)
    call broadcast(nu_krook_mb)
    call broadcast(smooth_ZFs)
    call broadcast(mb_zf_option_switch)
    call broadcast(krook_option_switch)
    call broadcast(krook_exponent)
    call broadcast(LR_debug_switch)
    call broadcast(RK_step)
    call broadcast(mb_debug_step)
    call broadcast(comm_at_init)

    call scope(crossdomprocs)

    if(job==1) then
      call receive(nakxl,0)
      call receive(nxl  ,0)
      call receive(nakxr,njobs-1)
      call receive(nxr  ,njobs-1)

      ! the following assumes nx in the center domain is some
      ! integer multiple of nx in the left or right domain.
      ! Also assumes dx is the same in every domain, which should
      ! be the case
      fac=nx/nxl
      x_fft_size=nakxl*fac
    else
      call send(nakx,1)
      call send(nx,1)
      x_fft_size = nakx
    endif

    call scope(subprocs)

    if(abs(nu_krook_mb) > epsilon(0.0)) then
      include_multibox_krook = .true.
    endif

  end subroutine read_multibox_parameters

  subroutine init_multibox 
    use constants, only: pi
    use stella_layouts, only: vmu_lo
    use stella_geometry, only: get_x_to_rho
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx,naky, akx, aky, nx,x, x_d, rho, x0, rho_clamped,rho_d_clamped
    use file_utils, only: runtype_option_switch, runtype_multibox
    use job_manage, only: njobs
    use mp, only: scope, crossdomprocs, subprocs, &
                  send, receive, job

    implicit none


    integer :: g_buff_size
    integer :: phi_buff_size
    integer :: i

    real, dimension (:), allocatable :: x_clamped, x_d_clamped

    real :: db


    if(runtype_option_switch /= runtype_multibox) return

    bs_fullgrid = nint((3.0*boundary_size)/2.0)

    phi_buff_size = boundary_size*naky*ntubes*(2*nzgrid+1)
    g_buff_size   = phi_buff_size*(vmu_lo%ulim_alloc-vmu_lo%llim_proc+1)

    if (.not.allocated(g_buffer0)) allocate(g_buffer0(g_buff_size))
    if (.not.allocated(g_buffer1)) allocate(g_buffer1(g_buff_size))
    if (.not.allocated(fsa_x) .and. (mb_zf_option_switch.eq.mb_zf_option_no_fsa)) then
      allocate(fsa_x(nakx)); fsa_x=0.0
    endif
    if (.not.allocated(copy_mask_left))  allocate(copy_mask_left(boundary_size));  copy_mask_left =1.0
    if (.not.allocated(copy_mask_right)) allocate(copy_mask_right(boundary_size)); copy_mask_right=1.0
    if (.not.allocated(krook_mask_left))  allocate(krook_mask_left(boundary_size));  krook_mask_left =0.0
    if (.not.allocated(krook_mask_right)) allocate(krook_mask_right(boundary_size)); krook_mask_right=0.0

    select case (krook_option_switch)
    case (krook_option_linear)
      db = 1.0/krook_size
      do i = 1, krook_size
        krook_mask_right(i) = i*db
        copy_mask_right(i) = 0.0
      enddo
    case (krook_option_exp)
      db = 3.0/krook_size
      do i = 1, krook_size
        krook_mask_right(i) = 1.0-(1.0-exp(-(krook_size-i)*db))/(1.0-exp(-3.0))
        copy_mask_right(i) = 0.0
      enddo
    case (krook_option_exp_rev)
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

    if(.not.allocated(krook_fac)) allocate (krook_fac(naky))

    if(naky>1) krook_fac = (aky/aky(2))**krook_exponent

#ifdef MPI
    call scope(crossdomprocs)

    if(.not.allocated(x_mb)) allocate(x_mb(x_fft_size))
    dx_mb = (2*pi*x0)/x_fft_size

    if(job==1) then

      xL = x(bs_fullgrid)
      xR = x(nx-bs_fullgrid+1)
      rhoL = rho(bs_fullgrid)
      rhoR = rho(nx-bs_fullgrid+1)
      xL_d = x_d(boundary_size)
      xR_d = x_d(x_fft_size-boundary_size+1)

      allocate(x_clamped(nx))
      allocate(x_d_clamped(nakx))

      if(LR_debug_switch == LR_debug_option_L) then
        x_clamped   = xL
        x_d_clamped = xL_d
      else if(LR_debug_switch == LR_debug_option_R) then
        x_clamped   = xR
        x_d_clamped = xR_d
      else
        x_clamped   = x
        x_d_clamped = x_d
        do i = 1, nx
          if(x_clamped(i) < xL) x_clamped(i) = xL
          if(x_clamped(i) > xR) x_clamped(i) = xR
        enddo
        do i = 1, nakx
          if(x_d_clamped(i) < xL_d) x_d_clamped(i) = xL_d
          if(x_d_clamped(i) > xR_d) x_d_clamped(i) = xR_d
        enddo
      endif

      call get_x_to_rho(1, x_clamped, rho_clamped)
      call get_x_to_rho(1, x_d_clamped, rho_d_clamped)

      deallocate (x_clamped, x_d_clamped)

    elseif(job==0) then
      do i = 1, x_fft_size
        x_mb(i) = (i-1)*dx_mb
      enddo

      call receive(xL,1)
      call receive(xL_d,1)
      call send(akx(2),1)
    elseif(job==njobs-1) then
      do i = 1, x_fft_size
        x_mb(i) = (i-1)*dx_mb
      enddo

      call receive(xR,1)
      call receive(xR_d,1)
      call send(akx(2),1)
    endif

    call scope(subprocs)
#endif

    call init_mb_transforms

  end subroutine init_multibox

  subroutine communicate_multibox_parameters
    use job_manage, only: njobs
    use mp, only: scope, crossdomprocs, subprocs, &
                  send, receive, job

    implicit none

#ifdef MPI
    if(job==1) then
      call scope(crossdomprocs)

      call send(xL  ,0)
      call send(xL_d,0)
      call send(xR  ,njobs-1)
      call send(xR_d,njobs-1)

      call receive(kx0_L,0)
      call receive(kx0_R,njobs-1)

      call scope(subprocs)

    endif
#endif
  end subroutine communicate_multibox_parameters

  subroutine finish_multibox

    implicit none

    if (allocated(g_buffer0))   deallocate (g_buffer0)
    if (allocated(g_buffer1))   deallocate (g_buffer1)
    if (allocated(fsa_x))       deallocate (fsa_x)

    if (allocated(copy_mask_left))   deallocate (copy_mask_left)
    if (allocated(copy_mask_right))  deallocate (copy_mask_right)
    if (allocated(krook_mask_left))  deallocate (krook_mask_left)
    if (allocated(krook_mask_right)) deallocate (krook_mask_right)

    if (allocated(fft_kxky)) deallocate (fft_kxky)
    if (allocated(fft_xky))  deallocate (fft_xky)
    if (allocated(fft_xy))   deallocate (fft_xy)

    call finish_mb_transforms

  end subroutine finish_multibox

  subroutine multibox_communicate (gin)

    use constants, only: zi
    use kt_grids, only: nakx,naky,naky_all, aky, nx,ny,dx,dy, zonal_mode
    use file_utils, only: runtype_option_switch, runtype_multibox
    use file_utils, only: get_unused_unit
    use fields_arrays, only: phi, phi_corr_QN
    use fields, only: advance_fields, fields_updated
    use job_manage, only: njobs
    use physics_flags, only: radial_variation,prp_shear_enabled, hammett_flow_shear
    use physics_parameters, only: g_exb, g_exbfac
    use stella_layouts, only: vmu_lo
    use stella_geometry, only: dl_over_b
    use zgrid, only: nzgrid
    use flow_shear, only:  shift_state
    use mp, only: job, scope, mp_abort,  &
                  crossdomprocs, subprocs, allprocs, &
                  send, receive, proc0

    implicit none

    integer :: num,ia, ix,iky,iz,it,iv,offset
    integer :: ii,jj, temp_unit
    real :: afacx, afacy
    complex :: dzm,dzp
    character(len=512) :: filename

    complex, dimension (:,:), allocatable :: prefac
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (inout) :: gin

#ifndef MPI
    return
#else
    if(runtype_option_switch /= runtype_multibox) return
    if(LR_debug_switch /= LR_debug_option_default) return
    if(njobs /= 3) call mp_abort("Multibox only supports 3 domains at the moment.")

    allocate (prefac(naky,x_fft_size)); prefac = 1.0

    if(prp_shear_enabled.and.hammett_flow_shear) then
      prefac = exp(-zi*g_exb*g_exbfac*spread(x_mb,1,naky)*spread(aky*shift_state,2,x_fft_size))
    endif

    if(mod(temp_ind,mb_debug_step)==0 .and. proc0) then
     ! call get_unused_unit(temp_unit)
      temp_unit=3023+job
      afacx = real(nx)/real(x_fft_size)
      afacy = real(ny)/real(2*naky-1)
      fft_kxky=phi(:,:,0,1)
      if(radial_variation) then
        fft_kxky = fft_kxky + phi_corr_QN(:,:,0,1)
      endif 
      call transform_kx2x(fft_kxky,fft_xky)  
      fft_xky=fft_xky*prefac
      call transform_ky2y(fft_xky,fft_xy)
      write (filename,"(A,I1,A,I0.6)") "phiout",job,"_",temp_ind
      open (unit=temp_unit, file=filename, status="replace",& 
            action="write",form="unformatted",access="stream")
      write (temp_unit) real(x_fft_size,4)
      do ii=1,x_fft_size
        write(temp_unit) real(afacx*dx*(ii-1),4)
      enddo
      do ii=1,naky_all
        write (temp_unit) real(afacy*dy*(ii-1),4)
        do jj=1,x_fft_size
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
          if(zonal_mode(1) .and. mb_zf_option_switch .eq. mb_zf_option_no_fsa) then
            do ix= 1,nakx
              fsa_x(ix) = sum(dl_over_b(ia,:)*gin(1,ix,:,it,iv))
            enddo
          endif

          do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
            fft_kxky = gin(:,:,iz,it,iv)

            if(zonal_mode(1)) then
              if(    mb_zf_option_switch .eq. mb_zf_option_no_ky0) then
                fft_kxky(1,:) = 0.0
              elseif(mb_zf_option_switch .eq. mb_zf_option_no_fsa) then
                fft_kxky(1,:) = fft_kxky(1,:) - fsa_x
              endif
            endif

            call transform_kx2x(fft_kxky,fft_xky)
            fft_xky = fft_xky*prefac
            do ix=1,boundary_size
              do iky=1,naky
                !DSO if in the future the grids can have different naky, one will
                !have to divide by naky here, and multiply on the receiving end
                g_buffer0(num) = fft_xky(iky,ix + offset)
                num=num+1
              enddo
            enddo
          enddo
        enddo
      enddo
! DSO - send data
      call send(g_buffer0,1,43 + job)
    else
      offset = x_fft_size - boundary_size
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
              dzm = fft_xky(1,boundary_size+1)            - fft_xky(1,boundary_size)
              dzp = fft_xky(1,x_fft_size-boundary_size+1) - fft_xky(1,x_fft_size-boundary_size)
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

    deallocate (prefac)

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
    allocate (g0x(naky,x_fft_size))
    
    offset = x_fft_size - boundary_size

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
          rhs(:,:,iz,it,ivmu) = rhs(:,:,iz,it,ivmu) &
                              - code_dt*nu_krook_mb*spread(krook_fac,2,nakx)*g0k
        enddo
      enddo
    enddo

    deallocate(g0k,g0x)

  end subroutine add_multibox_krook

!!>DSO - The following subroutines are the _unpadded_ analogues of the ones found in
! stella_transforms.f90. 
! The ones uses here ensure that the grid spacing in real space is consistent between 
! domains (since we do not keep the checkboard mode)

  subroutine init_mb_transforms

    use stella_layouts, only: init_stella_layouts
    use kt_grids, only: nakx, naky, naky_all

    implicit none

    if (mb_transforms_initialized) return
    mb_transforms_initialized = .true.

    if (.not.allocated(fft_kxky)) allocate (fft_kxky(naky  ,nakx))
    if (.not.allocated(fft_xky))  allocate (fft_xky(naky   ,x_fft_size))
    if (.not.allocated(fft_xy))   allocate (fft_xy(naky_all,x_fft_size))

    call init_x_fft
    call init_y_fft

  end subroutine init_mb_transforms

  subroutine init_x_fft

    use fft_work, only: init_ccfftw, FFT_BACKWARD, FFT_FORWARD

    implicit none

    if (.not.allocated(fft_x_k)) allocate (fft_x_k(x_fft_size))
    if (.not.allocated(fft_x_x)) allocate (fft_x_x(x_fft_size))

    call init_ccfftw (xf_fft, FFT_BACKWARD, x_fft_size, fft_x_k, fft_x_x)
    call init_ccfftw (xb_fft, FFT_FORWARD , x_fft_size, fft_x_x, fft_x_k)

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

    use kt_grids, only: naky, ikx_max

    implicit none

    complex, dimension (:,:), intent (in)  :: gkx
    complex, dimension (:,:), intent (out) :: gx

    integer :: iy

    do iy = 1, naky
       fft_x_k = 0.
       fft_x_k(:ikx_max) = gkx(iy,:ikx_max)
       fft_x_k((x_fft_size-ikx_max+2):) = gkx(iy,ikx_max+1:)
       call dfftw_execute_dft(xf_fft%plan, fft_x_k, fft_x_x)
       gx(iy,:) = fft_x_x*xf_fft%scale
    end do

  end subroutine transform_kx2x

  subroutine transform_x2kx (gx, gkx)

    use kt_grids, only: naky, ikx_max

    implicit none

    complex, dimension (:,:), intent (in)  :: gx
    complex, dimension (:,:), intent (out) :: gkx

    integer :: iy

    do iy = 1, naky
       fft_x_x = gx(iy,:)
       call dfftw_execute_dft(xb_fft%plan, fft_x_x, fft_x_k)
       gkx(iy,:ikx_max)   = fft_x_k(:ikx_max)*xb_fft%scale
       gkx(iy,ikx_max+1:) = fft_x_k((x_fft_size-ikx_max+2):)*xb_fft%scale
    end do

  end subroutine transform_x2kx

  subroutine transform_ky2y (gky, gy)

    implicit none

    complex, dimension (:,:), intent (in) :: gky
    real, dimension (:,:), intent (out) :: gy

    integer :: ikx

    do ikx = 1, x_fft_size
       fft_y_k = gky(:,ikx)
       call dfftw_execute_dft_c2r(yf_fft%plan, fft_y_k, fft_y_y)
       gy(:,ikx) =fft_y_y*yf_fft%scale
    end do

  end subroutine transform_ky2y

!   subroutine transform_y2ky (gy, gky)
!
!    implicit none
!
!    real, dimension (:,:), intent (in out) :: gy
!    complex, dimension (:,:), intent (out) :: gky
!
!    integer :: ikx
!
!    do ikx = 1, x_fft_size
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
