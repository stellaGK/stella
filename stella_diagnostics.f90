module stella_diagnostics

  implicit none

  public :: init_stella_diagnostics, finish_stella_diagnostics
  public :: diagnose_stella

  private

  integer :: ntg_out
  integer :: nwrite, nsave, nmovie
  logical :: save_for_restart
  logical :: write_phi_vs_time
  logical :: write_gvmus
  logical :: write_gzvs

  integer :: stdout_unit

  ! arrays needed for averaging in x,y,z
  real, dimension (:), allocatable :: fac
  real, dimension (:), allocatable :: dl_over_b

  real, dimension (:,:,:), allocatable :: pflux, vflux, qflux, exchange
  real, dimension (:), allocatable :: pflux_avg, vflux_avg, qflux_avg, heat_avg

  integer :: nout = 1
  logical :: diagnostics_initialized = .false.

  logical :: debug = .false.

contains

  subroutine init_stella_diagnostics (nstep)

    use zgrid, only: init_zgrid
    use kt_grids, only: init_kt_grids
    use run_parameters, only: init_run_parameters
    use species, only: init_species
    use dist_fn, only: init_dist_fn
    use init_g, only: init_init_g
    use stella_io, only: init_stella_io
    use mp, only: broadcast, proc0

    implicit none

    integer, intent (in) :: nstep

    integer :: nmovie_tot

    if (diagnostics_initialized) return
    diagnostics_initialized = .true.
    
    debug = debug .and. proc0
    
    call init_zgrid
    call init_kt_grids
    call init_run_parameters
    call init_species
    call init_init_g
    call init_dist_fn
    
    call read_parameters
    call allocate_arrays
    
    call broadcast (nwrite)
    call broadcast (nmovie)
    call broadcast (nsave)
    call broadcast (save_for_restart)
    call broadcast (write_phi_vs_time)
    call broadcast (write_gvmus)
    call broadcast (write_gzvs)
    
    nmovie_tot = nstep/nmovie
    
    ! note that init_averages needs ntg_out, defined in read_parameters
    call init_averages
    call init_stella_io (write_phi_vs_time, write_gvmus, write_gzvs)
    call open_loop_ascii_files
    
  end subroutine init_stella_diagnostics
  
  subroutine read_parameters

    use mp, only: proc0
    use file_utils, only: input_unit_exist
    use zgrid, only: nperiod, ntheta

    implicit none

    logical :: exist
    integer :: in_file

    namelist /stella_diagnostics_knobs/ nwrite, nmovie, nsave, &
         save_for_restart, write_phi_vs_time, write_gvmus, write_gzvs

    if (proc0) then
       nwrite = 100
       nmovie = 100000
       nsave = -1
       save_for_restart = .false.
       write_phi_vs_time = .false.
       write_gvmus = .false.
       write_gzvs = .false.

       in_file = input_unit_exist ("stella_diagnostics_knobs", exist)
       if (exist) read (unit=in_file, nml=stella_diagnostics_knobs)

       if (.not. save_for_restart) nsave = -1
    end if
    ntg_out = ntheta/2 + (nperiod-1)*ntheta

  end subroutine read_parameters

  subroutine allocate_arrays

    use species, only: nspec
    use kt_grids, only: nakx, naky

    implicit none

    if (.not.allocated(pflux)) allocate(pflux (nakx,naky,nspec)) ; pflux = 0.
    if (.not.allocated(qflux)) allocate(qflux (nakx,naky,nspec)) ; qflux = 0.
    if (.not.allocated(vflux)) allocate(vflux (nakx,naky,nspec)) ; vflux = 0.
    if (.not.allocated(exchange)) allocate(exchange (nakx,naky,nspec)) ; exchange = 0.
    if (.not.allocated(pflux_avg)) allocate(pflux_avg(nspec)) ; pflux_avg = 0.
    if (.not.allocated(qflux_avg)) allocate(qflux_avg(nspec)) ; qflux_avg = 0.
    if (.not.allocated(vflux_avg)) allocate(vflux_avg(nspec)) ; vflux_avg = 0.
    if (.not.allocated(heat_avg)) allocate(heat_avg(nspec)) ; heat_avg = 0.

  end subroutine allocate_arrays

  subroutine init_averages

    use zgrid, only: delthet
    use geometry, only: jacob
    use kt_grids, only: akx, nakx

    implicit none

    if (.not.allocated(dl_over_b)) then
       allocate (dl_over_b(-ntg_out:ntg_out))
       dl_over_b = delthet(-ntg_out:ntg_out)*jacob(-ntg_out:ntg_out)
       dl_over_b = dl_over_b / sum(dl_over_b)
    end if

    if (.not.allocated(fac)) then
       allocate (fac(nakx)) ; fac = 2.0
       if (akx(1)<epsilon(0.)) fac(1) = 1.0
    end if

  end subroutine init_averages

  subroutine open_loop_ascii_files

    use file_utils, only: open_output_file

    implicit none

    call open_output_file (stdout_unit,'.out')

  end subroutine open_loop_ascii_files

  subroutine close_loop_ascii_files
    
    use file_utils, only: close_output_file
    
    implicit none
    
    call close_output_file (stdout_unit)

  end subroutine close_loop_ascii_files

  subroutine diagnose_stella (istep)

    use mp, only: proc0
    use fields_arrays, only: phi, apar
    use dist_fn_arrays, only: gvmu, gnew
    use stella_io, only: write_time_nc
    use stella_io, only: write_phi_nc
    use stella_io, only: write_gvmus_nc
    use stella_io, only: write_gzvs_nc
    use stella_time, only: code_time
    use zgrid, only: ntgrid
    use vpamu_grids, only: nvgrid, nmu
    use species, only: nspec

    implicit none

    integer, intent (in) :: istep
    
    real :: phi2, apar2
    real, dimension (:,:,:), allocatable :: gvmus
    real, dimension (:,:,:), allocatable :: gzvs

    ! only write data to file every nwrite time steps
    if (mod(istep,nwrite) /= 0) return

    if (proc0) then
       call volume_average (phi, phi2)
       call volume_average (apar, apar2)
       write (*,'(a7,i7,a6,e12.4,a10,e12.4,a11,e12.4)') 'istep=', istep, &
            'time=', code_time, '|phi|^2=', phi2, '|apar|^2= ', apar2
       call write_loop_ascii_files (istep, phi2, apar2)
    end if

    if (proc0) then
       if (debug) write (*,*) 'stella_diagnostics::write_time_nc'
       call write_time_nc (nout, code_time)
       if (write_phi_vs_time) then
          if (debug) write (*,*) 'stella_diagnostics::diagnose_stella::write_phi_nc'
          call write_phi_nc (nout, phi)
       end if
    end if
    if (write_gvmus) then
       allocate (gvmus(2*nvgrid+1,nmu,nspec))
       if (debug) write (*,*) 'stella_diagnostics::diagnose_stella::get_gvmus'
       ! note that gvmus is h at this point
       call get_gvmus (gvmu, gvmus)
       if (debug) write (*,*) 'stella_diagnostics::diagnose_stella::write_gvmus_nc'
       if (proc0) call write_gvmus_nc (nout, gvmus)
       deallocate (gvmus)
    end if
    if (write_gzvs) then
       allocate (gzvs(2*ntgrid+1,2*nvgrid+1,nspec))
       if (debug) write (*,*) 'stella_diagnostics::diagnose_stella::get_gzvs'
       call get_gzvs (gnew, gzvs)
       if (debug) write (*,*) 'stella_diagnostics::diagnose_stella::write_gzvs_nc'
       if (proc0) call write_gzvs_nc (nout, gzvs)
       deallocate (gzvs)
    end if

    nout = nout + 1

  end subroutine diagnose_stella

  subroutine volume_average (unavg, avg)

    use zgrid, only: ntgrid
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-ntgrid:), intent (in) :: unavg
    real, intent (out) :: avg

    integer :: iky, ikx, ig

    avg = 0.
    do ig = -ntgrid, ntgrid
       do ikx = 1, nakx
          do iky = 1, naky
             avg = avg + real(unavg(iky,ikx,ig)*conjg(unavg(iky,ikx,ig)))*fac(ikx)*dl_over_b(ig)
          end do
       end do
    end do

  end subroutine volume_average

  ! get_gvmus takes g(kx,ky,theta) and returns average over theta of int dxdy g(x,y,theta)^2
  ! SHOULD MODIFY TO TAKE ADVANTAGE OF FACT THAT G(KY,KX,Z) LOCAL IS AVAILABLE
  subroutine get_gvmus (g, gv)

    use mp, only: nproc, sum_reduce
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: is_idx, ikx_idx, ig_idx
    use vpamu_grids, only: nvgrid, nmu

    implicit none

    complex, dimension (-nvgrid:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    real, dimension (:,:,:), intent (out) :: gv

    integer :: ikxkyz, iv, is, imu, ig, ikx, ivp

    ! when doing volume averages, note the following:
    ! int dxdy g(x,y)^2 = sum_kx |g(kx=0,ky)|^2 + 2 * sum_{kx,ky} |g(kx>0,ky)|^2
    ! factor of 2 accounted for in fac

    gv = 0.
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       is = is_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       ig = ig_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          do iv = -nvgrid, nvgrid
             ivp = iv+nvgrid+1
             gv(ivp,imu,is) = gv(ivp,imu,is) + real(g(iv,imu,ikxkyz)*conjg(g(iv,imu,ikxkyz)))*fac(ikx)*dl_over_b(ig)
          end do
       end do
    end do

    if (nproc > 1) call sum_reduce (gv,0)

  end subroutine get_gvmus

  ! get_gzvs takes g(kx,ky,theta,vpa,mu,s) and returns int dmudxdy g(x,y,theta,vpa,mu,s)^2
  subroutine get_gzvs (g, gz)

    use stella_layouts, only: vmu_lo
    use zgrid, only: ntgrid
    use vpamu_grids, only: nvgrid
    use vpamu_grids, only: integrate_mu
    use kt_grids, only: nakx, naky

    implicit none

    complex, dimension (:,:,-ntgrid:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (:,:,:), intent (out) :: gz

    integer :: ivmu, ig, ikx, iky, igp

    real, dimension (:,:), allocatable :: gtmp

    allocate (gtmp(-ntgrid:ntgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    ! when doing volume averages, note the following:
    ! int dxdy g(x,y)^2 = sum_kx |g(kx,ky=0)|^2 + 2 * sum_{kx,ky} |g(kx,ky>0)|^2
    ! factor of 2 accounted for in fac

    gtmp = 0.
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do ikx = 1, nakx
          do iky = 1, naky
             gtmp(:,ivmu) = gtmp(:,ivmu) + real(g(iky,ikx,:,ivmu)*conjg(g(iky,ikx,:,ivmu)))*fac(ikx)
          end do
       end do
    end do

    do ig = -ntgrid, ntgrid
       igp = ig+ntgrid+1
       call integrate_mu (ig, gtmp(ig,:), gz(igp,:,:))
    end do

    deallocate (gtmp)

  end subroutine get_gzvs

  subroutine finish_stella_diagnostics

    use mp, only: proc0
    use stella_io, only: finish_stella_io

    implicit none

    if (proc0) then
       call write_final_ascii_files
       call close_loop_ascii_files
    end if
    call finish_stella_io
    call finish_averages
    call deallocate_arrays

    nout = 1
    diagnostics_initialized = .false.

  end subroutine finish_stella_diagnostics

  subroutine write_loop_ascii_files (istep, phi2, apar2)

    use stella_time, only: code_time

    implicit none
    
    integer, intent (in) :: istep
    real, intent (in) :: phi2, apar2

    write (stdout_unit,'(a7,i7,a6,e12.4,a10,e12.4,a11,e12.4)'), 'istep=', istep, &
         'time=', code_time, '|phi|^2=', phi2, '|apar|^2= ', apar2

  end subroutine write_loop_ascii_files

  subroutine write_final_ascii_files

    use file_utils, only: open_output_file, close_output_file
    use fields_arrays, only: phi, apar
    use zgrid, only: ntgrid
    use zgrid, only: theta
    use kt_grids, only: naky, nakx
    use kt_grids, only: aky, akx, theta0

    implicit none

    integer :: tmpunit
    integer :: iky, ikx, ig

    call open_output_file (tmpunit,'.final_fields')
    write (tmpunit,'(8a12)') '# theta', 'thet-thet0', 'aky', 'akx', &
         'real(phi)', 'imag(phi)', 'real(apar)', 'imag(apar)'
    do iky = 1, naky
       do ikx = 1, nakx
          do ig = -ntgrid, ntgrid
             write (tmpunit,'(8e12.4)') theta(ig), theta(ig)-theta0(iky,ikx), aky(iky), akx(ikx), &
                  real(phi(iky,ikx,ig)), aimag(phi(iky,ikx,ig)), &
                  real(apar(iky,ikx,ig)), aimag(apar(iky,ikx,ig))
          end do
          write (tmpunit,*)
       end do
    end do
    call close_output_file (tmpunit)
    
  end subroutine write_final_ascii_files

  subroutine finish_averages

    implicit none

    if (allocated(dl_over_b)) deallocate (dl_over_b)
    if (allocated(fac)) deallocate (fac)

  end subroutine finish_averages

  subroutine deallocate_arrays

    implicit none

    if (allocated(pflux)) deallocate (pflux)
    if (allocated(qflux)) deallocate (qflux)
    if (allocated(vflux)) deallocate (vflux)
    if (allocated(exchange)) deallocate (exchange)
    if (allocated(pflux_avg)) deallocate (pflux_avg)
    if (allocated(qflux_avg)) deallocate (qflux_avg)
    if (allocated(vflux_avg)) deallocate (vflux_avg)
    if (allocated(heat_avg)) deallocate (heat_avg)

  end subroutine deallocate_arrays

end module stella_diagnostics
