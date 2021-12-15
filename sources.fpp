module sources

#if defined MPI && defined ISO_C_BINDING
  use mpi
#endif

  implicit none

  public :: init_sources, finish_sources
  public :: init_quasineutrality_source
  public :: init_source_timeaverage
  public :: update_quasineutrality_source
  public :: include_qn_source
  public :: include_krook_operator, update_tcorr_krook
  public :: remove_zero_projection, project_out_zero
  public :: add_krook_operator
  public :: tcorr_source, exclude_boundary_regions, exp_fac
  public :: int_krook, int_proj
  public :: qn_source_initialized
  public :: time_sources

  private

  logical :: include_krook_operator, remove_zero_projection
  logical :: krook_odd, exclude_boundary_regions
  logical :: from_zero
  logical :: conserve_momentum, conserve_density
  integer:: ikxmax_source
  real :: nu_krook, tcorr_source, int_krook, int_proj
  real :: exp_fac

  logical :: qn_source_initialized, include_qn_source

  logical :: debug = .false.

#if defined MPI && defined ISO_C_BINDING
  logical :: qn_window_initialized = .false.
#endif

  real, dimension (2,2) :: time_sources = 0.

contains

  subroutine init_sources

    use mp, only: job, proc0
    use run_parameters, only: fphi
    use run_parameters, only: ky_solve_radial, ky_solve_real
    use kt_grids, only: naky, nakx, zonal_mode
    use zgrid, only: nzgrid, ntubes
    use stella_layouts, only: vmu_lo
    use dist_fn_arrays, only: g_krook, g_proj, g_symm
    use fields_arrays, only : phi_proj, phi_proj_stage
    use physics_flags, only: radial_variation
    use species, only: spec, has_electron_species
    use physics_flags, only: adiabatic_option_switch
    use physics_flags, only: adiabatic_option_fieldlineavg
    use file_utils, only: runtype_option_switch, runtype_multibox

    implicit none

    logical :: has_elec, adia_elec
    real :: fac

    debug = debug .and. proc0

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
      allocate (phi_proj(nakx,-nzgrid:nzgrid,ntubes)); phi_proj = 0.
    endif
    if (.not.allocated(phi_proj_stage)) then
      allocate (phi_proj_stage(nakx,-nzgrid:nzgrid,ntubes)); phi_proj_stage = 0.
    endif

    if ((conserve_momentum.or.conserve_density).and..not.allocated(g_symm)) then
      allocate (g_symm(naky, nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    endif

    fac = 1.
    if (from_zero) fac = 0.

    if (int_krook.lt.0.) int_krook = fac*tcorr_source
    if (int_proj.lt.0.)  int_proj  = fac*tcorr_source

    include_qn_source = .false.
    if (fphi > epsilon(0.0).and.radial_variation.and.ky_solve_radial.gt.0) then
      has_elec  = has_electron_species(spec)
      adia_elec = .not.has_elec.and.zonal_mode(1) &
                  .and.adiabatic_option_switch == adiabatic_option_fieldlineavg
      if(adia_elec) then
        if(runtype_option_switch.ne.runtype_multibox.or.(job.eq.1.and..not.ky_solve_real)) then
          include_qn_source = .true.
        endif
      endif
    endif

  end subroutine init_sources

  subroutine read_parameters

    use file_utils, only: input_unit_exist
    use physics_flags, only: full_flux_surface, radial_variation
    use mp, only: proc0, broadcast
    use kt_grids, only: ikx_max, periodic_variation
    use fields_arrays, only: tcorr_source_qn, exclude_boundary_regions_qn

    implicit none

    namelist /sources/ &
         include_krook_operator, nu_krook, tcorr_source, remove_zero_projection, &
         ikxmax_source, krook_odd, exclude_boundary_regions, &
         tcorr_source_qn, exclude_boundary_regions_qn, from_zero, &
         conserve_momentum, conserve_density

    integer :: in_file
    logical :: dexist

    if (proc0) then
       include_krook_operator = .false.
       exclude_boundary_regions = radial_variation.and..not.periodic_variation
       exclude_boundary_regions_qn = exclude_boundary_regions
       remove_zero_projection = .false.
       nu_krook = 0.05
       tcorr_source = 0.02
       tcorr_source_qn = -1.0
       ikxmax_source = 1 ! kx=0
       if(periodic_variation) ikxmax_source = 2 ! kx=0 and kx=1
       krook_odd = .true. ! damp only the odd mode that can affect profiles
       from_zero = .true.

       conserve_momentum = .false.
       conserve_density = .false.

       in_file = input_unit_exist("sources", dexist)
       if (dexist) read (unit=in_file, nml=sources)

       if (tcorr_source_qn.lt.0) tcorr_source_qn = tcorr_source
    end if

    ikxmax_source = min(ikxmax_source,ikx_max)

    int_proj = -1.
    int_krook= -1.

    call broadcast (include_krook_operator)
    call broadcast (exclude_boundary_regions)
    call broadcast (exclude_boundary_regions_qn)
    call broadcast (nu_krook)
    call broadcast (tcorr_source)
    call broadcast (tcorr_source_qn)
    call broadcast (ikxmax_source)
    call broadcast (remove_zero_projection)
    call broadcast (krook_odd)
    call broadcast (from_zero)
    call broadcast (conserve_momentum)
    call broadcast (conserve_density)

  end subroutine read_parameters

  subroutine init_source_timeaverage

    use stella_time, only: code_dt
    use fields_arrays, only: tcorr_source_qn, exp_fac_qn

    implicit none

    exp_fac    = exp(-code_dt/tcorr_source)
    exp_fac_qn = exp(-code_dt/tcorr_source_qn)

  end subroutine init_source_timeaverage

  subroutine finish_sources

    use dist_fn_arrays, only: g_krook, g_proj, g_symm
    use fields_arrays, only : phi_proj, phi_proj_stage
#if defined MPI && defined ISO_C_BINDING
    use fields_arrays, only : qn_zf_window
#else
    use fields_arrays, only : phizf_solve, phi_ext
#endif

    implicit none

    integer :: ierr

    if (allocated(g_krook))  deallocate (g_krook)
    if (allocated(g_proj))   deallocate (g_proj)
    if (allocated(g_symm))   deallocate (g_symm)
    if (allocated(phi_proj)) deallocate (phi_proj)
    if (allocated(phi_proj_stage)) deallocate (phi_proj_stage)

#if defined MPI && defined ISO_C_BINDING
    if (qn_zf_window.ne.MPI_WIN_NULL) then 
      call mpi_win_free(qn_zf_window, ierr)
      qn_window_initialized = .false.
    endif
#else    
    if (associated(phizf_solve%zloc)) deallocate (phizf_solve%zloc)
    if (associated(phizf_solve%idx)) deallocate (phizf_solve%idx)
    if (associated(phi_ext)) deallocate (phi_ext)
#endif

  end subroutine finish_sources

  subroutine add_krook_operator (g, gke_rhs)

    use mp, only: proc0
    use job_manage, only: time_message
    use zgrid, only: nzgrid, ntubes
    use constants, only: pi, zi
    use kt_grids, only: naky, akx, nakx, zonal_mode, boundary_size
    use stella_layouts, only: vmu_lo
    use stella_time, only: code_dt
    use dist_fn_arrays, only: g_krook, g_symm
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
    use physics_flags, only: radial_variation

    implicit none

    complex :: tmp
    integer :: ikx, jkx, iz, it, ia, ivmu, npts

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (in), target :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gke_rhs
    
    complex, dimension (:,:,:,:,:), pointer :: g_work

    complex, dimension (:,:), allocatable :: g0k, g0x, g1x
    real, dimension (:), allocatable :: basis_func

    ia = 1
    if (.not.zonal_mode(1)) return

    if (proc0) call time_message(.false.,time_sources(:,1), ' sources')

    g_work => g
    if (conserve_momentum.or.conserve_density) then
      g_work => g_symm
      g_work = g
    endif

    if (conserve_momentum) call enforce_momentum_conservation (g_work)
    if (conserve_density)  call enforce_density_conservation (g_work(1,:,:,:,:))
    
    if (exclude_boundary_regions) then
      npts = nakx - 2*boundary_size
      allocate (g0k(1,nakx))
      allocate (g0x(1,nakx))
      allocate (g1x(1,nakx))
      allocate (basis_func(npts))
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            g0k(1,:) = g_work(1,:,iz,it,ivmu)
            g1x = 0.
            call transform_kx2x_unpadded(g0k,g0x)
            do ikx = 1, ikxmax_source
              if (ikx.eq.1) then
                basis_func = 1.0
                tmp = sum(g0x(1,(boundary_size+1):(nakx-boundary_size)))/real(npts)
              else
                do jkx = 1, npts
                  basis_func(jkx) = sin(2.0*pi*(ikx-1)*jkx/real(npts+1))
                enddo
                tmp = 2.0*sum(basis_func*g0x(1,(boundary_size+1):(nakx-boundary_size)))/real(npts+1)
              endif
              if(tcorr_source.gt.epsilon(0.0)) then
                tmp = (code_dt*tmp + exp_fac*int_krook*g_krook(ikx,iz,it,ivmu)) &
                    / (code_dt     + exp_fac*int_krook)
              endif
              do jkx = 1, npts
                g1x(1,boundary_size+jkx) = g1x(1,boundary_size+jkx) + tmp*basis_func(jkx) 
              enddo
            enddo
            call transform_x2kx_unpadded (g1x,g0k)
            gke_rhs(1,:,iz,it,ivmu) = gke_rhs(1,:,iz,it,ivmu) - code_dt*nu_krook*g0k(1,:)
          enddo
        enddo
      enddo
      deallocate (g0k, g0x, g1x, basis_func)
    else
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
              if(abs(akx(ikx)).gt.akx(ikxmax_source)) cycle
              tmp = g_work(1,ikx,iz,it,ivmu)
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

    if (proc0) call time_message(.false.,time_sources(:,1), ' sources')

  end subroutine add_krook_operator

  subroutine update_tcorr_krook (g)

    use mp, only: proc0
    use job_manage, only: time_message
    use constants, only: pi, zi
    use dist_fn_arrays, only: g_krook, g_symm
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: akx, nakx, zonal_mode, boundary_size
    use stella_layouts, only: vmu_lo
    use stella_time, only: code_dt
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), target, intent (in) :: g
    complex, dimension (:,:), allocatable :: g0k, g0x

    complex, dimension (:,:,:,:,:), pointer :: g_work

    integer :: ivmu, iz, it, ikx, jkx, ia, npts
    real :: int_krook_old
    complex :: tmp

    if(.not.zonal_mode(1)) return

    if (proc0) call time_message(.false.,time_sources(:,1), ' sources')

    ia = 1

    g_work => g
    if (conserve_momentum.or.conserve_density) then
      g_work => g_symm
      g_work = g
    endif

    if (conserve_momentum) call enforce_momentum_conservation (g_work)
    if (conserve_density)  call enforce_density_conservation (g_work(1,:,:,:,:))

    int_krook_old = int_krook
    int_krook =  code_dt + exp_fac*int_krook_old

    if (exclude_boundary_regions) then
      npts = nakx - 2*boundary_size
      allocate (g0k(1,nakx))
      allocate (g0x(1,nakx))
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            g0k(1,:) = g_work(1,:,iz,it,ivmu)
            call transform_kx2x_unpadded(g0k,g0x)
            do ikx = 1, ikxmax_source
              if (ikx.eq.1) then
                tmp = sum(g0x(1,(boundary_size+1):(nakx-boundary_size)))/real(npts)
              else
                tmp = 0.
                do jkx = 1, npts
                  tmp = tmp + sin(2.0*pi*(ikx-1)*jkx/real(npts+1))*g0x(1,boundary_size+jkx)
                enddo
                tmp = 2.0*tmp/real(npts+1)
              endif
              g_krook(ikx,iz,it,ivmu) = (code_dt*tmp + exp_fac*int_krook_old*g_krook(ikx,iz,it,ivmu)) &
                                        /int_krook
            enddo
          enddo
        enddo
      enddo
      deallocate (g0k, g0x)
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

    if (proc0) call time_message(.false.,time_sources(:,1), ' sources')

  end subroutine update_tcorr_krook


  subroutine enforce_momentum_conservation (g_work)

    use mp, only: proc0
    use job_manage, only: time_message
    use redistribute, only: scatter, gather
    use stella_layouts, only: vmu_lo, kxkyz_lo
    use stella_layouts, only: imu_idx, is_idx, iv_idx
    use vpamu_grids , only: nvgrid, nvpa, nmu
    use dist_redistribute, only: kxkyz2vmu
    use dist_fn_arrays, only: gvmu
    use zgrid, only: nzgrid

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (inout) :: g_work

    integer :: ikxkyz, imu, iv, iv2
    complex :: tmp


    if (proc0) call time_message(.false.,time_sources(:,2), ' source_redist')
    call scatter (kxkyz2vmu, g_work, gvmu)
    if (proc0) call time_message(.false.,time_sources(:,2), ' source_redist')

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
      do imu = 1, nmu
        do iv = 1, nvgrid
          iv2 = nvpa-iv+1
          tmp = 0.5*(gvmu(iv,imu,ikxkyz) + gvmu(iv2,imu,ikxkyz))
          gvmu(iv ,imu,ikxkyz) = tmp
          gvmu(iv2,imu,ikxkyz) = tmp
        end do
      end do
    end do

    if (proc0) call time_message(.false.,time_sources(:,2), ' source_redist')
    call gather (kxkyz2vmu, gvmu, g_work)
    if (proc0) call time_message(.false.,time_sources(:,2), ' source_redist')

  end subroutine enforce_momentum_conservation

  subroutine enforce_density_conservation (g_work)

    use mp, only: sum_allreduce
    use species, only: spec, nspec
    use physics_flags, only: radial_variation
    use vpamu_grids , only: integrate_species, mu, vpa, vperp2
    use vpamu_grids , only: maxwell_vpa, maxwell_mu, maxwell_fac
    use kt_grids, only: nakx, rho_d_clamped
    use stella_layouts, only: vmu_lo, imu_idx, is_idx, iv_idx
    use stella_geometry, only: bmag, dBdrho, dl_over_b, d_dl_over_b_drho
    use gyro_averages, only: gyro_average, aj0x, aj1x
    use dist_fn_arrays, only: kperp2, dkperp2dr
    use zgrid, only: nzgrid, ntubes
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

    implicit none

    complex, dimension (:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (inout) :: g_work
    integer :: ia, ikx, it, iz, imu, iv, ivmu, is

    complex, dimension (:,:), allocatable :: gyro_g, g0k, g0x
    complex, dimension (:), allocatable :: g_fsa
    real :: energy

    ia = 1

    allocate (gyro_g(nakx,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate (g0k(1,nakx))
    allocate (g0x(1,nakx))
    allocate (g_fsa(nakx)); g_fsa = 0.

    do it = 1, ntubes
      do iz = -nzgrid, nzgrid
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          iv = iv_idx(vmu_lo,ivmu)
          gyro_g(:,ivmu) = g_work(:,iz,it,ivmu)*aj0x(1,:,iz,ivmu)
          g0k = 0.0
          if(radial_variation) then
            g0k(1,:) = gyro_g(:,ivmu) &
                     * (-0.5*aj1x(1,:,iz,ivmu)/aj0x(1,:,iz,ivmu)*(spec(is)%smz)**2 &
                     * (kperp2(1,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                     * (dkperp2dr(1,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                     + dBdrho(iz)/bmag(ia,iz))

            call transform_kx2x_unpadded (g0k,g0x)
            g0x(1,:) = rho_d_clamped*g0x(1,:)
            call transform_x2kx_unpadded (g0x,g0k)
          endif
          gyro_g(:,ivmu) = gyro_g(:,ivmu) + g0k(1,:)
        enddo

        do ikx = 1, nakx
          call integrate_species (gyro_g(ikx,:), iz, spec%dens_psi0,g0k(1,ikx),reduce_in=.false.)
        enddo
        call sum_allreduce (g0k)

        !we now have delta n. Flux surface average
        call transform_kx2x_unpadded (g0k,g0x)
        g_fsa = g_fsa + dl_over_b(ia,iz)*g0x(1,:)
        if(radial_variation) g_fsa = g_fsa + d_dl_over_b_drho(ia,iz)*rho_d_clamped*g0x(1,:)
      enddo
    enddo
    g_fsa = g_fsa/ntubes

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      is = is_idx(vmu_lo,ivmu)
      imu = imu_idx(vmu_lo,ivmu)
      iv = iv_idx(vmu_lo,ivmu)
      do it = 1, ntubes
        do iz = -nzgrid, nzgrid
          g0x(1,:) = g_fsa
          !multiply by f0
          g0x(1,:) = g0x(1,:)*maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is)
          if (radial_variation) then !variation in the density cancels
            energy = (vpa(iv)**2 + vperp2(ia,iz,imu))*(spec(is)%temp_psi0/spec(is)%temp)
            g0x(1,:) = g0x(1,:)*(1.0 - rho_d_clamped*(spec(is)%tprim*(energy-1.5) + 2.*mu(imu)*dBdrho(iz)))
          endif

          call transform_x2kx_unpadded(g0x,g0k)
          g_work(:,iz,it,ivmu) = g_work(:,iz,it,ivmu) - g0k(1,:)
        enddo
      enddo
    enddo
    deallocate (gyro_g,g0k,g0x,g_fsa)

  end subroutine enforce_density_conservation
  
  subroutine project_out_zero (gold, gnew)

    use mp, only: proc0
    use job_manage, only: time_message
    use zgrid, only: nzgrid, ntubes
    use constants, only: pi, zi
    use kt_grids, only: zonal_mode, akx, nakx, boundary_size
    use stella_layouts, only: vmu_lo
    use stella_time, only: code_dt
    use dist_fn_arrays, only: g_proj, g_symm
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded


    implicit none

    complex :: tmp
    integer :: ikx, jkx, iz, it, ia, ivmu, npts

    complex, dimension (:,:), allocatable :: g0k, g0x, g1x
    real, dimension (:), allocatable :: basis_func
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (in) :: gold
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (inout) :: gnew
    complex, allocatable, dimension (:,:,:,:) :: g

    ia = 1
    if(.not.zonal_mode(1)) return

    if (proc0) call time_message(.false.,time_sources(:,1), ' sources')

    allocate (g(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    !divide by code_dt to ensure time averaging is performed correctly
    if (conserve_momentum) then
      g_symm = (gnew - gold)/code_dt
      call enforce_momentum_conservation (g_symm)
      g = g_symm(1,:,:,:,:)

    else
      g = (gnew(1,:,:,:,:) - gold(1,:,:,:,:))/code_dt
    endif

    if (conserve_density) call enforce_density_conservation (g)


    if (exclude_boundary_regions) then
      npts = nakx - 2*boundary_size
      allocate (g0k(1,nakx))
      allocate (g0x(1,nakx))
      allocate (g1x(1,nakx))
      allocate (basis_func(npts))
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            g0k(1,:) = g(:,iz,it,ivmu)
            g1x = 0.
            call transform_kx2x_unpadded (g0k,g0x)
            do ikx = 1, ikxmax_source
              !physical region should have an odd number of collocation points
              if (ikx.eq.1) then
                basis_func = 1.0
                tmp = sum(g0x(1,(boundary_size+1):(nakx-boundary_size)))/real(npts)
              else
                ! here we use a Fourier basis due to periodicity, 
                ! though we could use Legendre polynomials
                ! NB: Only a constant or linear function (or nearly linear, i.e. first
                ! sine harmonic) make physical sense as sources, so ikxmax_source <= 2
                do jkx = 1, npts
                  basis_func(jkx) = sin(2.0*pi*(ikx-1)*jkx/real(npts+1))
                enddo
                tmp = 2.0*sum(basis_func*g0x(1,(boundary_size+1):(nakx-boundary_size)))/real(npts+1)
              endif
              if(tcorr_source.gt.epsilon(0.)) then
                tmp = (code_dt*tmp + exp_fac*int_proj*g_proj(ikx,iz,it,ivmu)) &
                    / (code_dt     + exp_fac*int_proj)
                g_proj(ikx,iz,it,ivmu) = tmp
              endif
              do jkx = 1, npts
                g1x(1,boundary_size+jkx) = g1x(1,boundary_size+jkx) + tmp*basis_func(jkx) 
              enddo
            enddo
            call transform_x2kx_unpadded (g1x,g0k)
            g(:,iz,it,ivmu) = g0k(1,:)
          enddo
        enddo
      enddo
      deallocate (g0k, g0x, g1x, basis_func)
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

    gnew(1,:,:,:,:) = gnew(1,:,:,:,:) - code_dt*g

    deallocate (g)

    if (proc0) call time_message(.false.,time_sources(:,1), ' sources')

  end subroutine project_out_zero

  subroutine init_quasineutrality_source

#if defined MPI && defined ISO_C_BINDING
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use fields_arrays, only: qn_zf_window
    use mp, only: sgproc0, curr_focus, mp_comm, sharedsubprocs, comm_sgroup
    use mp, only: scope, real_size, nbytes_real
    use mp_lu_decomposition, only: lu_decomposition_local, lu_inverse_local
    use mpi
#endif
    use physics_flags, only: radial_variation
    use stella_geometry, only: dl_over_b, d_dl_over_b_drho
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
    use zgrid, only: nzgrid, nztot
    use kt_grids, only: naky, nakx, rho_d_clamped, boundary_size
    use linear_solve, only: lu_decomposition
    use fields_arrays, only: phizf_solve, c_mat, theta, phi_ext
    use fields_arrays, only: tcorr_source_qn, exclude_boundary_regions_qn, exp_fac_qn
    

    implicit none

    integer :: iz, ikx, ia, jkx, jz
    integer :: inmat, jnmat, nmat_zf
    real :: dum
#if defined MPI && ISO_C_BINDING
    integer :: prior_focus, ierr, temp_window
    integer :: disp_unit = 1
    integer*8 :: cur_pos
    integer (kind=MPI_ADDRESS_KIND) :: win_size
    type(c_ptr) :: cptr
    complex, dimension (:,:), pointer :: temp_mat
#endif
    complex, dimension (:,:), allocatable :: g0k, g0x, g1k

    ia = 1

    if (qn_source_initialized) return
    qn_source_initialized = .true.

    if (include_qn_source) then
      nmat_zf = nakx*(nztot-1)
#if defined MPI && ISO_C_BINDING
      if((.not.qn_window_initialized).or.(qn_zf_window.eq.MPI_WIN_NULL)) then
        prior_focus = curr_focus
        call scope (sharedsubprocs)
        win_size = 0
        if(sgproc0) then
          win_size = int(nmat_zf,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND & 
                   + int(nmat_zf*(nmat_zf+1),MPI_ADDRESS_KIND)*2*real_size !complex size
        endif

        if (debug) write (6,*) 'sources::init_quasineutrality_source::win_allocate'
        call mpi_win_allocate_shared(win_size,disp_unit,MPI_INFO_NULL, &
                                     mp_comm,cptr,qn_zf_window,ierr)

        if (.not.sgproc0) then
          !make sure all the procs have the right memory address
          call mpi_win_shared_query(qn_zf_window,0,win_size,disp_unit,cptr,ierr)
        end if
        call mpi_win_fence(0,qn_zf_window,ierr)
        cur_pos = transfer(cptr,cur_pos)

        !allocate the memory
        if (.not.associated(phizf_solve%zloc))  then
          cptr = transfer(cur_pos,cptr)
          call c_f_pointer (cptr,phizf_solve%zloc,(/nmat_zf,nmat_zf/))
        endif
        cur_pos = cur_pos + nmat_zf**2*2*nbytes_real

        if (.not.associated(phi_ext))  then
          cptr = transfer(cur_pos,cptr)
          call c_f_pointer (cptr,phi_ext,(/nmat_zf/))
        endif
        cur_pos = cur_pos + nmat_zf*2*nbytes_real

        if (.not.associated(phizf_solve%idx))  then
          cptr = transfer(cur_pos,cptr)
          call c_f_pointer (cptr,phizf_solve%idx,(/nmat_zf/))
        endif

        call mpi_win_fence(0,qn_zf_window,ierr)

        call scope (prior_focus)

        qn_window_initialized = .true.
      endif
#else
      if (debug) write (6,*) 'sources::init_quasineutrality_source::allocate_phizf'
      if (.not.associated(phizf_solve%zloc)) allocate (phizf_solve%zloc(nmat_zf,nmat_zf))
      if (.not.associated(phizf_solve%idx))  allocate (phizf_solve%idx(nmat_zf))
      if (.not.associated(phi_ext)) allocate (phi_ext(nmat_zf))
#endif

#if defined MPI && ISO_C_BINDING
      if(sgproc0) then
#endif
        allocate (g0k(1,nakx))
        allocate (g0x(1,nakx))
        allocate (g1k(1,nakx))

        phizf_solve%zloc = 0.

        !get the big matrix
        do jz = -nzgrid, nzgrid-1
          do jkx = 1, nakx
            jnmat = jkx + nakx*(jz+nzgrid)
            ! C.phi
            do ikx = 1, nakx
              inmat = ikx + nakx*(jz+nzgrid)
              phizf_solve%zloc(inmat,jnmat) = phizf_solve%zloc(inmat,jnmat) + c_mat(ikx,jkx)
            enddo

            ! -C.<phi>_\psi
            g0k = 0.0; g0k(1,jkx) = 1.0
            call transform_kx2x_unpadded (g0k,g0x)
            g0x(1,:) = (dl_over_b(ia,jz) + d_dl_over_b_drho(ia,jz)*rho_d_clamped)*g0x(1,:)
            call transform_x2kx_unpadded(g0x,g0k)

            !set the gauge potential
            if(jkx.eq.1) g0k(1,1) = 0. 

            do ikx = 1, nakx
              g1k(1,ikx) = sum(c_mat(ikx,:)*g0k(1,:))
            enddo

            do iz = -nzgrid, nzgrid-1
              do ikx = 1, nakx
                inmat = ikx + nakx*(iz+nzgrid)
                phizf_solve%zloc(inmat,jnmat) = phizf_solve%zloc(inmat,jnmat) - g1k(1,ikx)
              enddo
            enddo

            ! get theta.phi
            g1k(1,:) = theta(:,jkx,jz)

            ! +theta.phi
            do ikx = 1, nakx
              inmat = ikx + nakx*(jz+nzgrid)
              phizf_solve%zloc(inmat,jnmat) = phizf_solve%zloc(inmat,jnmat) + g1k(1,ikx)
            enddo

            ! -<<theta.phi>_psi>_T
            call transform_kx2x_unpadded (g1k, g0x)
            g0x(1,:) = (dl_over_b(ia,jz) + d_dl_over_b_drho(ia,jz)*rho_d_clamped)*g0x(1,:)

            if (exclude_boundary_regions_qn) then
              g0x(1,:) = sum(g0x(1,(boundary_size+1):(nakx-boundary_size))) &
                       / (nakx - 2*boundary_size)
              g0x(1,1:boundary_size) = 0.0
              g0x(1,(nakx-boundary_size+1):) = 0.0
            else
              g0x(1,:) = sum(g0x(1,:))/nakx
            endif

            call transform_x2kx_unpadded(g0x,g0k)

            if (tcorr_source_qn.gt.epsilon(0.)) then
              g0k = (1. - exp_fac_qn)*g0k
            endif

            do iz = -nzgrid, nzgrid-1
              do ikx = 1, nakx
                inmat = ikx + nakx*(iz+nzgrid)
                phizf_solve%zloc(inmat,jnmat) = phizf_solve%zloc(inmat,jnmat) &
                                              - g0k(1,ikx)
              enddo
            enddo
          enddo
        enddo
        deallocate (g0k,g1k,g0x)
#if defined MPI && ISO_C_BINDING
      endif

      call mpi_win_fence(0,qn_zf_window,ierr)
      if (debug) write (6,*) 'sources::init_quasineutrality_source::lu_decomposition'
      call lu_decomposition_local(comm_sgroup, 0, qn_zf_window, &
                                  phizf_solve%zloc, phizf_solve%idx, dum)

      call mpi_win_fence(0,qn_zf_window,ierr)

      if (debug) write (6,*) 'sources::init_quasineutrality_source::temp_mat'
      call scope (sharedsubprocs)
      win_size = 0
      if (sgproc0) then
        win_size = int(nmat_zf**2,MPI_ADDRESS_KIND)*2*real_size !complex size
      endif

      if (debug) write (6,*) 'sources::init_quasineutrality_source::win_allocate'
      call mpi_win_allocate_shared(win_size,disp_unit,MPI_INFO_NULL, &
                                   mp_comm,cptr,temp_window,ierr)

      if (.not.sgproc0) then
        !make sure all the procs have the right memory address
        call mpi_win_shared_query(temp_window,0,win_size,disp_unit,cptr,ierr)
      end if

      call mpi_win_fence(0,temp_window,ierr) !DSO - not sure if this fence is necessary
      call c_f_pointer (cptr,temp_mat,(/nmat_zf,nmat_zf/))

      if (sgproc0) temp_mat = phizf_solve%zloc

      call mpi_win_fence(0,temp_window,ierr)
      call scope (prior_focus)

      ! inverse is calculated since it is more straightforward to parallelize
      ! inverse calculation/matrix multiplication than the lu back substitution
      if (debug) write (6,*) 'sources::init_quasineutrality_source::lu_inverse'
      call lu_inverse_local(comm_sgroup, qn_zf_window, &
                            temp_mat, phizf_solve%idx, phizf_solve%zloc)
      call mpi_win_free (temp_window, ierr)
#else
      if (debug) write (6,*) 'sources::init_quasineutrality_source::lu_decomposition'
      call lu_decomposition(phizf_solve%zloc, phizf_solve%idx, dum)
#endif
    endif
  end subroutine init_quasineutrality_source

  subroutine update_quasineutrality_source

    use fields_arrays, only: phi_proj, phi_proj_stage
    use fields_arrays, only: tcorr_source_qn, exp_fac_qn

    implicit none

    if (tcorr_source_qn.lt.epsilon(0.)) then
      phi_proj = phi_proj_stage
    else
      phi_proj = exp_fac_qn*phi_proj + (1.-exp_fac_qn)*phi_proj_stage
    endif

  end subroutine update_quasineutrality_source

end module sources
