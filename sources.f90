module sources

  implicit none

  public :: init_sources, finish_sources
  public :: include_krook_operator, update_tcorr_krook
  public :: remove_zero_projection, project_out_zero
  public :: add_krook_operator
  public :: tcorr_source
  public :: int_krook, int_proj, int_QN

  private

  logical :: include_krook_operator, remove_zero_projection
  logical :: krook_odd, exclude_boundary_regions
  real :: nu_krook, tcorr_source, int_krook, int_proj, int_QN
  integer:: ikxmax_source


contains

  subroutine init_sources

    use kt_grids, only: nakx
    use zgrid, only: nzgrid, ntubes
    use stella_layouts, only: vmu_lo
    use dist_fn_arrays, only: g_krook, g_proj
    use fields_arrays, only : phi_proj

    implicit none

    call read_parameters

    if (include_krook_operator.and..not.allocated(g_krook)) then
      allocate (g_krook(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g_krook = 0.
    endif

    if (remove_zero_projection.and..not.allocated(g_proj)) then
      allocate (g_proj(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g_proj = 0.
    endif

    if (.not.allocated(phi_proj)) then
      allocate (phi_proj(nakx,1,ntubes)); phi_proj = 0.
    endif

    int_krook = 0.
    int_proj  = 0.

  end subroutine init_sources

  subroutine read_parameters

    use file_utils, only: input_unit_exist
    use physics_flags, only: full_flux_surface, radial_variation
    use mp, only: proc0, broadcast
    use kt_grids, only: ikx_max, periodic_variation

    implicit none

    namelist /sources/ &
         include_krook_operator, nu_krook, tcorr_source, remove_zero_projection, &
         ikxmax_source, krook_odd, exclude_boundary_regions

    integer :: in_file
    logical :: dexist

    if (proc0) then
       include_krook_operator = .false.
       exclude_boundary_regions = radial_variation.and..not.periodic_variation
       remove_zero_projection = .false.
       nu_krook = 0.05
       tcorr_source =0.02
       ikxmax_source = 1 ! kx=0
       if(periodic_variation) ikxmax_source = 2 ! kx=0 and kx=1
       krook_odd = .true. ! damp only the odd mode that can affect profiles

       in_file = input_unit_exist("sources", dexist)
       if (dexist) read (unit=in_file, nml=sources)
    end if

    ikxmax_source = min(ikxmax_source,ikx_max)

    call broadcast (include_krook_operator)
    call broadcast (exclude_boundary_regions)
    call broadcast (nu_krook)
    call broadcast (tcorr_source)
    call broadcast (ikxmax_source)
    call broadcast (remove_zero_projection)
    call broadcast (krook_odd)

  end subroutine read_parameters

  subroutine finish_sources

    use dist_fn_arrays, only: g_krook, g_proj

    implicit none

    if(allocated(g_krook)) deallocate(g_krook)
    if(allocated(g_proj))  deallocate(g_proj)

  end subroutine finish_sources

  subroutine add_krook_operator (g, gke_rhs)

    use zgrid, only: nzgrid, ntubes
    use constants, only: zi
    use kt_grids, only: akx, nakx, zonal_mode
    use stella_layouts, only: vmu_lo
    use stella_time, only: code_dt
    use dist_fn_arrays, only: g_krook
    use multibox, only: boundary_size
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

    implicit none

    real :: exp_fac
    complex :: tmp
    integer :: ikx, iz, it, ia, ivmu

    !complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), optional, intent (in) :: f0
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gke_rhs

    complex, dimension (:,:), allocatable :: g0k, g0x

    ia = 1
    if(.not.zonal_mode(1)) return

    exp_fac = exp(-code_dt/tcorr_source)

    !TODO: add number and momentum conservation
    if (exclude_boundary_regions) then
      allocate(g0k(1,nakx))
      allocate(g0x(1,nakx))
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            g0k(1,:) = g(1,:,iz,it,ivmu)
            call transform_kx2x_unpadded(g0k,g0x)
            tmp = sum(g0x(1,(boundary_size+1):(nakx-boundary_size)))/real(nakx-2*boundary_size)
            if(tcorr_source.le.epsilon(0.0)) then
              g0x = tmp
            else
              g0x = (code_dt*tmp + exp_fac*int_krook*g_krook(1,iz,it,ivmu)) &
                  / (code_dt     + exp_fac*int_krook)
            endif
            g0x(1,1:boundary_size) = 0.0
            g0x(1,(nakx-boundary_size+1):nakx) = 0.0
            call transform_x2kx_unpadded (g0x,g0k)
            gke_rhs(1,:,iz,it,ivmu) = gke_rhs(1,:,iz,it,ivmu) - code_dt*nu_krook*g0k(1,:)
          enddo
        enddo
      enddo
      deallocate(g0k, g0x)
    else
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
              if(abs(akx(ikx)).gt.akx(ikxmax_source)) cycle
              tmp = g(1,ikx,iz,it,ivmu)
              if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) tmp = zi*aimag(tmp)
              if(tcorr_source.le.epsilon(0.0)) then
                gke_rhs(1,ikx,iz,it,ivmu) = gke_rhs(1,ikx,iz,it,ivmu) - code_dt*nu_krook*tmp
              else
                gke_rhs(1,ikx,iz,it,ivmu) = gke_rhs(1,ikx,iz,it,ivmu) - code_dt*nu_krook &
                                          * (code_dt*tmp + exp_fac*int_krook*g_krook(ikx,iz,it,ivmu)) &
                                          / (code_dt     + exp_fac*int_krook)
              endif
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine add_krook_operator

  subroutine update_tcorr_krook (g)

    use constants, only: zi
    use dist_fn_arrays, only: g_krook
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: akx, nakx, zonal_mode
    use stella_layouts, only: vmu_lo
    use stella_time, only: code_dt
    use multibox, only: boundary_size
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (in) :: g
    complex, dimension (:,:), allocatable :: g0k, g0x

    integer :: ivmu, iz, it, ikx, ia
    real :: int_krook_old, exp_fac
    complex :: tmp

    if(.not.zonal_mode(1)) return

    exp_fac = exp(-code_dt/tcorr_source)

    ia = 1

    int_krook_old = int_krook
    int_krook =  code_dt + exp_fac*int_krook_old

    if (exclude_boundary_regions) then
      allocate(g0k(1,nakx))
      allocate(g0x(1,nakx))
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            g0k(1,:) = g(1,:,iz,it,ivmu)
            call transform_kx2x_unpadded(g0k,g0x)
            tmp = sum(g0x(1,(boundary_size+1):(nakx-boundary_size)))/real(nakx-2*boundary_size)
            g_krook(:,iz,it,ivmu) = (code_dt*tmp + exp_fac*int_krook_old*g_krook(:,iz,it,ivmu))/int_krook
          enddo
        enddo
      enddo
      deallocate(g0k, g0x)
    else
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
              tmp = g(1,ikx,iz,it,ivmu)
              if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) tmp = zi*aimag(tmp)
              g_krook(ikx,iz,it,ivmu) = (code_dt*tmp + exp_fac*int_krook_old*g_krook(ikx,iz,it,ivmu))/int_krook
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine update_tcorr_krook

  subroutine project_out_zero (g)

    use zgrid, only: nzgrid, ntubes
    use constants, only: zi
    use kt_grids, only: zonal_mode, akx, nakx
    use stella_layouts, only: vmu_lo
    use stella_time, only: code_dt
    use dist_fn_arrays, only: g_proj
    use multibox, only: boundary_size
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

    implicit none

    real :: exp_fac
    complex :: tmp
    integer :: ikx, iz, it, ia, ivmu

    complex, dimension (:,:), allocatable :: g0k, g0x
    complex, dimension (:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (inout) :: g

    ia = 1
    if(.not.zonal_mode(1)) return

    exp_fac = exp(-code_dt/tcorr_source)

    if (exclude_boundary_regions) then !here we do not require ikxmax_source
      allocate (g0k(1,nakx))
      allocate (g0x(1,nakx))
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            g0k(1,:) = g(:,iz,it,ivmu)
            call transform_kx2x_unpadded (g0k,g0x)
            tmp = sum(g0x(1,(boundary_size+1):(nakx-boundary_size)))/real(nakx-2*boundary_size)
            if(tcorr_source.le.epsilon(0.)) then
              g0x = tmp
            else
              g0x = (code_dt*tmp + exp_fac*int_proj*g_proj(1,iz,it,ivmu)) &
                  / (code_dt     + exp_fac*int_proj)
            endif
            g_proj(1,iz,it,ivmu) = g0x(1,1)
            g0x(1,1:boundary_size) = 0.0
            g0x(1,(nakx-boundary_size+1):nakx) = 0.0
            call transform_x2kx_unpadded (g0x,g0k)
            g(:,iz,it,ivmu) = g0k(1,:)
          enddo
        enddo
      enddo
      deallocate (g0k, g0x)
    else
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
              if(abs(akx(ikx)).gt.akx(ikxmax_source)) then
                g(ikx,iz,it,ivmu) = 0.0
              else
                tmp = g(ikx,iz,it,ivmu)
                if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) tmp = zi*aimag(tmp)
                if(tcorr_source.le.epsilon(0.)) then
                  g(ikx,iz,it,ivmu) = tmp
                else
                  g(ikx,iz,it,ivmu) = (code_dt*tmp + exp_fac*int_proj*g_proj(ikx,iz,it,ivmu)) &
                                    / (code_dt     + exp_fac*int_proj)
                endif
              endif
              if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) then
                g_proj(ikx,iz,it,ivmu) = zi*aimag(g(ikx,iz,it,ivmu))
              else
                g_proj(ikx,iz,it,ivmu) = g(ikx,iz,it,ivmu)
              endif
            enddo
          enddo
        enddo
      enddo
    endif

    int_proj = code_dt + exp_fac*int_proj

  end subroutine project_out_zero

end module sources
