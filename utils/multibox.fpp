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
  public :: kx0_L, kx0_R
  public :: RK_step
  public :: include_multibox_krook


  private

  complex, dimension (:), allocatable :: g_buffer0
  complex, dimension (:), allocatable :: g_buffer1
  complex, dimension (:), allocatable :: fsa_x
  real,    dimension (:), allocatable :: copy_mask_left, copy_mask_right
  real,    dimension (:), allocatable :: krook_mask_left, krook_mask_right
  real,    dimension (:), allocatable :: krook_fac
  
  complex, dimension (:,:), allocatable :: fft_xky
  real, dimension (:,:), allocatable :: fft_xy

  integer :: temp_ind = 0
  integer :: bs_fullgrid
  integer :: mb_debug_step
  

  real :: xL = 0., xR = 0.
  real :: xL_d = 0., xR_d = 0.
  real :: kx0_L, kx0_R

!>DSO 
! the multibox simulation has one parameter: boundary_size
! 
  integer :: boundary_size, krook_size
  real :: nu_krook_mb, krook_exponent
  logical :: smooth_ZFs
  logical :: RK_step, include_multibox_krook
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

    implicit none

    integer :: in_file, ierr
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
                                   mb_debug_step, krook_exponent

    if(runtype_option_switch /= runtype_multibox) return

    boundary_size = 4
    krook_size = 0
    krook_exponent = 0.0
    nu_krook_mb = 0.0
    mb_debug_step = 1000
    smooth_ZFs = .false.
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

    if(abs(nu_krook_mb) > epsilon(0.0)) then
      include_multibox_krook = .true.
    endif

  end subroutine read_multibox_parameters

  subroutine init_multibox 
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx,naky,naky_all, akx, aky, nx,x, x_clamped, x_d
    use file_utils, only: runtype_option_switch, runtype_multibox
    use job_manage, only: njobs
    use mp, only: scope, crossdomprocs, subprocs, &
                  send, receive, job

    implicit none


    integer :: g_buff_size
    integer :: phi_buff_size
    integer :: i

    real :: db

    if (.not.allocated(fft_xky)) allocate (fft_xky(naky,nakx))
    if (.not.allocated(fft_xy)) allocate (fft_xy(naky_all,nakx))

    if(.not.allocated(x_clamped)) allocate(x_clamped(nx)); x_clamped = 0.

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

      do i = 1, krook_size
        write (*,*) krook_mask_right(i), job
      enddo

    do i = 1, boundary_size
      copy_mask_left(i)  = copy_mask_right(boundary_size - i + 1)
      krook_mask_left(i) = krook_mask_right(boundary_size - i + 1)
    enddo

    if(.not.allocated(krook_fac)) allocate (krook_fac(naky))

    krook_fac = (aky/aky(2))**krook_exponent

#ifdef MPI
    call scope(crossdomprocs)

    if(job==1) then
      xL = x(bs_fullgrid)
      xR = x(nx-bs_fullgrid+1)
      xL_d = x_d(boundary_size)
      xR_d = x_d(nakx-boundary_size+1)

      if(LR_debug_switch == LR_debug_option_L) then
        x_clamped = xL
      else if(LR_debug_switch == LR_debug_option_R) then
        x_clamped = xR
      else
        x_clamped = x
        do i = 1, nx
          if(x_clamped(i) < xL) x_clamped(i) = xL
          if(x_clamped(i) > xR) x_clamped(i) = xR
        enddo
      endif
    elseif(job==0) then
      call receive(xL,1)
      call receive(xL_d,1)
      call send(akx(2),1)
    elseif(job==njobs-1) then
      call receive(xR,1)
      call receive(xR_d,1)
      call send(akx(2),1)
    endif

    call scope(subprocs)
#endif
  end subroutine init_multibox

  subroutine communicate_multibox_parameters
    use job_manage, only: njobs
    use mp, only: scope, crossdomprocs, subprocs, &
                  send, receive, job

    implicit none
#ifndef MPI
    return
#else
    if(job==1) then
      call scope(crossdomprocs)

      call send(xL,0)
      call send(xR,njobs-1)
      call send(xL_d,0)
      call send(xR_d,njobs-1)

      call receive(kx0_L,0)
      call receive(kx0_R,njobs-1)

      call scope(subprocs)
    endif
#endif
  end subroutine communicate_multibox_parameters

  subroutine finish_multibox

    use kt_grids, only: x_clamped

    implicit none

    if (allocated(g_buffer0))   deallocate (g_buffer0)
    if (allocated(g_buffer1))   deallocate (g_buffer1)
    if (allocated(fsa_x))       deallocate (fsa_x)
    if (allocated(x_clamped))   deallocate (x_clamped)

    if (allocated(copy_mask_left))   deallocate (copy_mask_left)
    if (allocated(copy_mask_right))  deallocate (copy_mask_right)
    if (allocated(krook_mask_left))  deallocate (krook_mask_left)
    if (allocated(krook_mask_right)) deallocate (krook_mask_right)

    if (allocated(fft_xky)) deallocate (fft_xky)
    if (allocated(fft_xy)) deallocate (fft_xy)

  end subroutine finish_multibox

  subroutine multibox_communicate (gin)

    use kt_grids, only: nakx,naky,naky_all,nx,ny,dx,dy, zonal_mode
    use file_utils, only: runtype_option_switch, runtype_multibox
    use file_utils, only: get_unused_unit
    use fields_arrays, only: phi
    use fields, only: advance_fields, fields_updated
    use job_manage, only: njobs
    use stella_layouts, only: vmu_lo
!   use stella_geometry, only: dl_over_b
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
    use stella_transforms, only: transform_ky2y_unpadded
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
    if(LR_debug_switch /= LR_debug_option_default) return
    if(njobs /= 3) call mp_abort("Multibox only supports 3 domains at the moment.")


    if(mod(temp_ind,mb_debug_step)==0 .and. proc0) then
     ! call get_unused_unit(temp_unit)
      temp_unit=3023+job
      afacx = real(nx)/real(nakx)
      afacy = real(ny)/real(2*naky-1)
      call transform_kx2x_unpadded(phi(:,:,0,1),fft_xky)  
      call transform_ky2y_unpadded(fft_xky,fft_xy)
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
!           do ix= 1,nakx
!             fft_x_k(ix) = sum(dl_over_b(ia,:)*gin(1,ix,:,it,iv))
!           enddo
!           call dfftw_execute_dft(xf_fft%plan, fft_x_k, fft_x_x)
!           fsa_x = fft_x_x*xf_fft%scale
          endif

          do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
            call transform_kx2x_unpadded(gin(:,:,iz,it,iv),fft_xky)  
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
            call transform_kx2x_unpadded(gin(:,:,iz,it,iv),fft_xky)  
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
            call transform_x2kx_unpadded(fft_xky,gin(:,:,iz,it,iv))  
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
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
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
          call transform_kx2x_unpadded(g(:,:,iz,it,ivmu),fft_xky)  
          do ix=1,boundary_size
            do iky=1,naky
                g0x(iky,ix)        = (fft_xky(iky,ix)        - g_buffer0(num))*krook_mask_left(ix)
                g0x(iky,ix+offset) = (fft_xky(iky,ix+offset) - g_buffer1(num))*krook_mask_right(ix)
                num=num+1
            enddo
          enddo
          call transform_x2kx_unpadded(g0x,g0k)  
          rhs(:,:,iz,it,ivmu) = rhs(:,:,iz,it,ivmu) &
                              - code_dt*nu_krook_mb*spread(krook_fac,2,nakx)*g0k
        enddo
      enddo
    enddo

    deallocate(g0k,g0x)

  end subroutine add_multibox_krook

end module multibox
