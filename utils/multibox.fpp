module multibox

  use fft_work, only: fft_type

  implicit none

  public :: init_multibox
  public :: finish_multibox
  public :: multibox_communicate
  public :: boundary_size
  public :: shear_rate

  private

  real, dimension (:), allocatable :: g_buffer0
  real, dimension (:), allocatable :: g_buffer1
  real, dimension (:), allocatable :: phi_buffer0
  real, dimension (:), allocatable :: phi_buffer1
  
  ! for swapped kxky space; could move to kt_grids.f90
  real, dimension (:), allocatable :: akx_2

  real, dimension (:,:), allocatable :: kernels
  real, dimension (:,:), allocatable :: conv_solv
  real, dimension (:), allocatable :: rhs
  real, dimension (:), allocatable :: solution

  integer, dimension (:), allocatable ::  bndry_ind
  integer, dimension (:), allocatable :: nbndry_ind
  integer, dimension (:), allocatable :: conv_idx

  complex, dimension (:,:), allocatable :: fft_swap
  complex, dimension (:,:), allocatable :: fft_kxy
  real, dimension (:,:), allocatable :: fft_xy

  ! for the unpadded FFTs
  type (fft_type) :: yf_fft, yb_fft
  type (fft_type) :: xf_fft, xb_fft

  complex, dimension (:), allocatable :: fft_y_in, fft_y_out
  complex, dimension (:), allocatable :: fft_x_k
  real, dimension (:), allocatable :: fft_x_x

  real :: conv_d
  real :: shear_rate

  integer :: g_buff_size
  integer :: phi_buff_size
  integer :: naky_f
  integer :: iL,iR
  logical :: mb_transforms_initialized = .false.
  logical :: mb_straight_copy = .true.

  integer :: temp_ind = 0

!>DSO 
! the multibox simulation has two parameters: bbits and  b_offset
! 
  integer :: bbits, b_offset, boundary_size

contains

  subroutine init_multibox
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx,naky,x0,akx,ikx_max
    use file_utils, only: input_unit_exist
    use file_utils, only: runtype_option_switch, runtype_multibox
    use mp, only: broadcast, proc0 
    use constants, only: pi, zi
    use fft_work, only: init_ccfftw, delete_fft,FFT_FORWARD, fft_type
    use linear_solve, only: lu_decomposition

    implicit none

#ifndef MPI
    return
#else

    integer :: i,j,k, in_file
    real dx, xm, xp, conv_fac
    real, dimension (nakx) :: x
    logical exist
    complex, dimension (nakx) :: fft_in
    complex, dimension (nakx) :: fft_out
    complex, dimension (nakx) :: expm, expp
    type (fft_type) :: conv_fft


    namelist /multibox_parameters/ bbits, b_offset, boundary_size, shear_rate

    if(runtype_option_switch /= runtype_multibox) return

    bbits = 1 ! just the 0th derivative as default
    b_offset = 0 
    boundary_size = 4
    shear_rate = 0


    if (proc0) then
      in_file = input_unit_exist("multibox_parameters", exist)
      if (exist) read (in_file, nml=multibox_parameters)
    endif

    call broadcast(bbits)
    call broadcast(b_offset)
    call broadcast(boundary_size)

    !calculate how many derivates need to be sent based on the bitmask
    if(.not.mb_straight_copy) then
      j=1
      boundary_size=0
      do while (j.le.bbits)
        if(IAND(bbits,j)>0) boundary_size=boundary_size+1
        j=j*2
      enddo
    endif

    naky_f = 2*naky-1

    phi_buff_size = boundary_size*naky_f*ntubes*(2*nzgrid+1)
    g_buff_size   = phi_buff_size*(vmu_lo%ulim_alloc-vmu_lo%llim_proc+1)

    if (.not.allocated(g_buffer0)) allocate(g_buffer0(g_buff_size))
    if (.not.allocated(g_buffer1)) allocate(g_buffer1(g_buff_size))
    if (.not.allocated(phi_buffer0)) allocate(phi_buffer0(phi_buff_size))
    if (.not.allocated(phi_buffer1)) allocate(phi_buffer1(phi_buff_size))

    if(.not. mb_straight_copy) then
      if (.not.allocated(kernels)) allocate(kernels(2*(boundary_size-1),nakx))
      if (.not.allocated(conv_solv)) allocate(conv_solv(2*(boundary_size-1),2*(boundary_size-1)))
      if (.not.allocated(conv_idx)) allocate(conv_idx(2*(boundary_size-1)))
      if (.not.allocated(rhs))  allocate(rhs(2*(boundary_size-1)))
      if (.not.allocated(solution))  allocate(solution(2*(boundary_size-1)))
      if (.not.allocated(bndry_ind))  allocate( bndry_ind(2*(boundary_size-1)))
      if (.not.allocated(nbndry_ind)) allocate(nbndry_ind(nakx-2*(boundary_size-1)))
      if (.not.allocated(akx_2)) allocate(akx_2(ikx_max))

      akx_2 = akx(1:ikx_max)
      iL = boundary_size - b_offset
      iR = nakx - boundary_size + 1 + b_offset

      !precompute convolution arrays and perform LU decomposition of our matrix

      call init_ccfftw(conv_fft, FFT_FORWARD, nakx, fft_in,fft_out)

      dx=2*pi*x0/nakx
      j=1
      k=1
      do i = 1, nakx
        x(i) = (i-1)*dx
        if((i==iL)            .or. & ! 0th derivative in left boundary
          ((i>boundary_size).and.(i<nakx-boundary_size+1)) .or. & ! physical domain
           (i==iR)) then    ! 0th derivative in right boundary
          nbndry_ind(j) = i
          j=j+1
        else
          bndry_ind(k) = i
          k=k+1
        endif
      enddo

      xm = x(iL)
      xp = x(iR)

      do i = 1, nakx
        expm(i) = exp(zi*akx(i)*xm)
        expp(i) = exp(zi*akx(i)*xp)
      enddo

      conv_fac = real(nakx)

      i=1
      j=2
      k=1
      do while (j.le.bbits)
        if(IAND(bbits,j)>0) then
          !left boundary convolution
          fft_in = (zi*akx)**i*expm
          call dfftw_execute_dft(conv_fft%plan,fft_in,fft_out)
          kernels(k,:) = real(fft_out)/conv_fac
          conv_solv(k,:) =  kernels(k,bndry_ind)
          k=k+1

          ! right boundary convolution
          fft_in = (zi*akx)**i*expp
          call dfftw_execute_dft(conv_fft%plan,fft_in,fft_out)
          kernels(k,:) = real(fft_out)/conv_fac
          conv_solv(k,:) =  kernels(k,bndry_ind)
          k=k+1
        endif
      
        i=i+1
        j=j*2
      end do

      call lu_decomposition(conv_solv,conv_idx,conv_d)
      call delete_fft(conv_fft)

    endif

    call init_mb_transforms

  
#endif
  end subroutine init_multibox

  subroutine finish_multibox

    implicit none

    if (allocated(g_buffer0))   deallocate (g_buffer0)
    if (allocated(g_buffer1))   deallocate (g_buffer1)
    if (allocated(phi_buffer0)) deallocate (phi_buffer0)
    if (allocated(phi_buffer1)) deallocate (phi_buffer1)
    if (allocated(kernels)) deallocate (kernels)
    if (allocated(conv_solv)) deallocate (conv_solv)
    if (allocated(conv_idx)) deallocate(conv_idx)
    if (allocated(rhs))  deallocate(rhs)
    if (allocated(solution))  deallocate(solution)
    if (allocated(bndry_ind))  deallocate( bndry_ind)
    if (allocated(nbndry_ind)) deallocate(nbndry_ind)
    if (allocated(akx_2)) deallocate(akx_2)

    call finish_mb_transforms

  end subroutine finish_multibox


  subroutine multibox_communicate

    use dist_fn_arrays, only: gnew
    use kt_grids, only: swap_kxky, swap_kxky_back,nakx, zonal_mode
    use file_utils, only: runtype_option_switch, runtype_multibox
    use fields_arrays, only: phi, apar
    use fields, only: advance_fields, fields_updated
    use job_manage, only: njobs
    use constants, only: zi
    use linear_solve, only: lu_back_substitution
    use stella_layouts, only: vmu_lo
    use mp, only: job, scope, mp_abort,  &
                  crossdomprocs, subprocs, &
                  send, receive

    implicit none

    integer :: deriv,ibit,num,ix,iy,iz,it,iv,offset
    integer :: ii,jj,fout
    character(len=1024) :: filename

#ifndef MPI
    return
#else
    if(runtype_option_switch /= runtype_multibox) return
    if(njobs /= 3) call mp_abort("Multibox only supports 3 domains at the moment.")

    if(mod(temp_ind,50)==0) then
      call swap_kxky(phi(:,:,0,1),fft_swap)
      call transform_ky2y_2d(fft_swap,fft_kxy)
      call transform_kx2x(fft_kxy,fft_xy)  
      write (filename,"(A,I1,A,I0.6)") "phiout",job,"_",temp_ind
      open (unit=fout, file=filename, status="replace",action="write")
      do ii=1,naky_f
        do jj=1,nakx
          write (fout,*) jj,ii,fft_xy(ii,jj)
        enddo
          write (fout,*) ""
      enddo
      close (unit=fout)
    endif

! DSO - change communicator
      call scope(crossdomprocs)

    if(mb_straight_copy) then
      if(job==0 .or. job==(njobs-1)) then
        offset=0;
        if(job==0) offset=nakx-boundary_size
        num=1
        do iv = vmu_lo%llim_proc, vmu_lo%ulim_proc
          do it = 1, vmu_lo%ntubes
            do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
              call swap_kxky(gnew(:,:,iz,it,iv),fft_swap)
              call transform_ky2y_2d(fft_swap,fft_kxy)
              call transform_kx2x(fft_kxy,fft_xy)  
              do ix=1,boundary_size
                do iy=1,naky_f
                  g_buffer0(num) = fft_xy(iy,ix + offset)
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
              call swap_kxky(gnew(:,:,iz,it,iv),fft_swap)
              call transform_ky2y_2d(fft_swap,fft_kxy)
              call transform_kx2x(fft_kxy,fft_xy)  
              do ix=1,boundary_size
                do iy=1,naky_f
                  fft_xy(iy,ix)        = g_buffer0(num)
                  fft_xy(iy,ix+offset) = g_buffer1(num)
                  num=num+1
                enddo
              enddo
              call transform_x2kx(fft_xy,fft_kxy)  
              call transform_y2ky_2d(fft_kxy, fft_swap)
              call swap_kxky_back(fft_swap,gnew(:,:,iz,it,iv))
            enddo
          enddo
        enddo
      endif
    else



! DSO - pack up data to be sent
! first do phi

      if(job==0 .or. job==(njobs-1)) then
        !zi*spread(akx,1,naky)
        num=1
        do it = 1, vmu_lo%ntubes
          do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid

            call swap_kxky(phi(:,:,iz,it),fft_swap)
            call transform_ky2y_2d(fft_swap,fft_kxy)
            !should perhaps enforce reality here
            do iy=1,naky_f
              deriv=0
              ibit=1
              do while (ibit.le.bbits)
                if(IAND(bbits,ibit)>0) then
                  fft_x_k = (zi*akx_2)**deriv*fft_kxy(iy,:)
                  call dfftw_execute_dft_c2r(xf_fft%plan, fft_x_k, fft_x_x)
                  fft_x_x = fft_x_x*xf_fft%scale
                  phi_buffer0(num) = fft_x_x(1) 
                  num=num+1
                endif
                deriv=deriv+1
                ibit=2*ibit
              enddo
            enddo

          !if(iz==3) then
          !    fft_buffer = spread((zi*akx_2)**i,1,naky_f)*fft_kxy
          !    call transform_kx2x(fft_buffer,fft_xy)
          !    do num=1,naky_f
          !      phi_buffer0(k) = fft_xy(num,1) 
          !      k=k+1
          !  fft_x_k = (zi*akx_2)**1*fft_kxy(5,:)
          !  call dfftw_execute_dft_c2r(xf_fft%plan, fft_x_k, fft_x_x)
          !  fft_x_x = fft_x_x*xf_fft%scale
          !  write(*,*) job, fft_x_x(1)
          !endif

          enddo
        enddo
! DSO - send data
        call send(phi_buffer0,1,43 + job)
      else
! DSO - receive the data

! left
        call receive(phi_buffer0,0, 43)
! right
        call receive(phi_buffer1,njobs-1, 43+njobs-1)

        num=1
        do it = 1, vmu_lo%ntubes
          do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
            call swap_kxky(phi(:,:,iz,it),fft_swap)
            call transform_ky2y_2d(fft_swap,fft_kxy)
            call transform_kx2x(fft_kxy,fft_xy)  
            do iy=1,naky_f
              fft_xy(iy,iL) = phi_buffer0(num)
              fft_xy(iy,iR) = phi_buffer1(num)
              num=num+1
              do ibit = 1, (boundary_size-1)
                rhs(2*(ibit-1) + 1) = phi_buffer0(num) - &
                        sum(fft_xy(iy,nbndry_ind)*kernels(2*(ibit-1)+1,nbndry_ind))
                rhs(2*(ibit-1) + 2) = phi_buffer1(num) - &
                        sum(fft_xy(iy,nbndry_ind)*kernels(2*(ibit-1)+2,nbndry_ind))
                num=num+1
              enddo
              if(IAND(bbits,1) < bbits) then ! derivatives used for boundary condition 
                call lu_back_substitution(conv_solv,conv_idx,rhs)
                fft_xy(iy,bndry_ind) = rhs
              endif
            enddo

            call transform_x2kx(fft_xy,fft_kxy)  
            call transform_y2ky_2d(fft_kxy, fft_swap)
            call swap_kxky_back(fft_swap,phi(:,:,iz,it))

          !if(iz==3) then
          !  call swap_kxky(phi(:,:,iz,it),fft_swap)
          !  call transform_ky2y_2d(fft_swap,fft_kxy)
          !  call transform_kx2x(spread((zi*akx_2)**1,1,naky_f)*fft_kxy,fft_xy)  
          !  write(*,*) job, fft_xy(5,iR)
          !endif
          enddo
        enddo
      endif

! now do g
      if(job==0 .or. job==(njobs-1)) then
        num=1
        do iv = vmu_lo%llim_proc, vmu_lo%ulim_proc
          do it = 1, vmu_lo%ntubes
            do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid

              call swap_kxky(gnew(:,:,iz,it,iv),fft_swap)
              call transform_ky2y_2d(fft_swap,fft_kxy)
          !should perhaps enforce reality here

              do iy=1,naky_f
                deriv=0
                ibit=1
                do while (ibit.le.bbits)
                  if(IAND(bbits,ibit)>0) then
                    fft_x_k = (zi*akx_2)**deriv*fft_kxy(iy,:)
                    call dfftw_execute_dft_c2r(xf_fft%plan, fft_x_k, fft_x_x)
                    fft_x_x = fft_x_x*xf_fft%scale
                    g_buffer0(num) = fft_x_x(1) 
                    num=num+1
                  endif
                  deriv=deriv+1
                  ibit=2*ibit
                enddo
              enddo
            enddo
          enddo
        enddo
! DSO - send data
        call send(g_buffer0,1,43 + job)
      else
! DSO - receive the data

! left
        call receive(g_buffer0,0, 43)
! right
        call receive(g_buffer1,njobs-1, 43+njobs-1)

        num=1
        do iv = vmu_lo%llim_proc, vmu_lo%ulim_proc
          do it = 1, vmu_lo%ntubes
            do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
              call swap_kxky(gnew(:,:,iz,it,iv),fft_swap)
              call transform_ky2y_2d(fft_swap,fft_kxy)
              call transform_kx2x(fft_kxy,fft_xy)  
              do iy=1,naky_f
                fft_xy(iy,iL) = g_buffer0(num)
                fft_xy(iy,iR) = g_buffer1(num)
                num=num+1
                do ibit = 1, (boundary_size-1)
                  rhs(2*(ibit-1) + 1) = g_buffer0(num) - &
                          sum(fft_xy(iy,nbndry_ind)*kernels(2*(ibit-1)+1,nbndry_ind))
                  rhs(2*(ibit-1) + 2) = g_buffer1(num) - &
                          sum(fft_xy(iy,nbndry_ind)*kernels(2*(ibit-1)+2,nbndry_ind))
                    num=num+1
                enddo
                if(IAND(bbits,1) < bbits) then ! derivatives used for boundary condition 
                  call lu_back_substitution(conv_solv,conv_idx,rhs)
                  fft_xy(iy,bndry_ind) = rhs
                endif
              enddo
              call transform_x2kx(fft_xy,fft_kxy)  
              call transform_y2ky_2d(fft_kxy, fft_swap)
              call swap_kxky_back(fft_swap,gnew(:,:,iz,it,iv))
            enddo
          enddo
        enddo
      endif

      if(zonal_mode(1)) then
        phi(1,:,:,:) = 0.0
        gnew(1,:,:,:,:) = 0.0
      endif
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

    use physics_flags, only: full_flux_surface
    use stella_layouts, only: init_stella_layouts
    use kt_grids, only: nakx, ikx_max

    implicit none

    if (mb_transforms_initialized) return
    mb_transforms_initialized = .true.

    if (.not.allocated(fft_swap)) allocate (fft_swap(naky_f,ikx_max))
    if (.not.allocated(fft_kxy)) allocate (fft_kxy(naky_f,ikx_max))
    if (.not.allocated(fft_xy)) allocate (fft_xy(naky_f,nakx))

    call init_y_fft
    call init_x_fft

  end subroutine init_mb_transforms

  subroutine init_y_fft

    use fft_work, only: init_ccfftw, FFT_BACKWARD, FFT_FORWARD

    implicit none

    if (.not.allocated(fft_y_in)) allocate (fft_y_in(naky_f))
    if (.not.allocated(fft_y_out)) allocate (fft_y_out(naky_f))

    call init_ccfftw (yf_fft, FFT_BACKWARD, naky_f, fft_y_in, fft_y_out)
    call init_ccfftw (yb_fft, FFT_FORWARD , naky_f, fft_y_in, fft_y_out)

  end subroutine init_y_fft

  subroutine init_x_fft

    use kt_grids, only: nakx, ikx_max
    use fft_work, only: init_crfftw, init_rcfftw, FFT_BACKWARD, FFT_FORWARD

    implicit none

    if (.not.allocated(fft_x_k)) allocate (fft_x_k(ikx_max))
    if (.not.allocated(fft_x_x)) allocate (fft_x_x(nakx))

    call init_crfftw (xf_fft, FFT_BACKWARD, nakx, fft_x_k, fft_x_x)
    call init_rcfftw (xb_fft, FFT_FORWARD , nakx, fft_x_x, fft_x_k)

  end subroutine init_x_fft

!  
!> transform routines start here
!

  subroutine transform_ky2y_2d (gky, gy)

    use kt_grids, only: ikx_max

    implicit none

    complex, dimension (:,:), intent (in) :: gky
    complex, dimension (:,:), intent (out) :: gy

    integer :: ikx

    do ikx = 1, ikx_max
       fft_y_in = gky(:,ikx)
       call dfftw_execute_dft(yf_fft%plan, fft_y_in, fft_y_out)
       gy(:,ikx) =fft_y_out*yf_fft%scale
    end do

  end subroutine transform_ky2y_2d

  subroutine transform_y2ky_2d (gy, gky)

    use kt_grids, only: ikx_max

    implicit none

    complex, dimension (:,:), intent (in out) :: gy
    complex, dimension (:,:), intent (out) :: gky

    integer :: ikx

    do ikx = 1, ikx_max
       fft_y_in = gy(:,ikx)
       call dfftw_execute_dft(yb_fft%plan, fft_y_in, fft_y_out)
       gky(:,ikx) = fft_y_out*yb_fft%scale
    end do

  end subroutine transform_y2ky_2d

  subroutine transform_kx2x (gkx, gx)

    implicit none

    complex, dimension (:,:), intent (in) :: gkx
    real, dimension (:,:), intent (out) :: gx

    integer :: iy

    do iy = 1, naky_f
       fft_x_k = gkx(iy,:)
       call dfftw_execute_dft_c2r(xf_fft%plan, fft_x_k, fft_x_x)
       gx(iy,:) = fft_x_x*xf_fft%scale
    end do

  end subroutine transform_kx2x

  subroutine transform_x2kx (gx, gkx)

    implicit none

    real, dimension (:,:), intent (in) :: gx
    complex, dimension (:,:), intent (out) :: gkx

    integer :: iy

    do iy = 1, naky_f
       fft_x_x = gx(iy,:)
       call dfftw_execute_dft_r2c(xb_fft%plan, fft_x_x, fft_x_k)
       gkx(iy,:) = fft_x_k*xb_fft%scale
    end do

  end subroutine transform_x2kx

  subroutine finish_mb_transforms

    use physics_flags, only: full_flux_surface

    implicit none

    call dfftw_destroy_plan (yf_fft%plan)
    call dfftw_destroy_plan (yb_fft%plan)
    call dfftw_destroy_plan (xf_fft%plan)
    call dfftw_destroy_plan (xb_fft%plan)
    if (allocated(fft_y_in)) deallocate (fft_y_in)
    if (allocated(fft_y_out)) deallocate (fft_y_out)
    if (allocated(fft_x_k)) deallocate (fft_x_k)
    if (allocated(fft_x_x)) deallocate (fft_x_x)
    if (allocated(fft_swap)) deallocate (fft_swap)
    if (allocated(fft_kxy)) deallocate (fft_kxy)
    if (allocated(fft_xy)) deallocate (fft_xy)

    mb_transforms_initialized = .false.

  end subroutine finish_mb_transforms

end module multibox

          !i=0 j=1 k=1
          !do while (j.le.bbits)
          !  if(IAND(bbits,j)>0) then
          !    fft_buffer = spread((zi*akx_2)**i,1,naky_f)*fft_kxy
          !    call transform_kx2x(fft_buffer,fft_xy)
          !    do num=1,naky_f
          !      phi_buffer0(k) = fft_xy(num,1) 
          !      k=k+1
          !    enddo
          !  endif
          !  i=i+1
          !  j=2*j
          !enddo
