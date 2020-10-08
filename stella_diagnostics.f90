module stella_diagnostics

  implicit none

  public :: init_stella_diagnostics, finish_stella_diagnostics
  public :: diagnose_stella
  public :: nsave

  private

  interface fieldline_average
     module procedure fieldline_average_real
     module procedure fieldline_average_complex
  end interface

  integer :: ntg_out
  integer :: nwrite, nsave, nmovie, navg
  logical :: save_for_restart
  logical :: write_omega
  logical :: write_moments
  logical :: write_phi_vs_time
  logical :: write_gvmus
  logical :: write_gzvs
  logical :: write_kspectra
!  logical :: write_symmetry

  integer :: stdout_unit, fluxes_unit, omega_unit

  ! arrays needed for averaging in x,y,z
  real, dimension (:), allocatable :: fac

  real, dimension (:,:,:), allocatable :: pflux, vflux, qflux, exchange
  real, dimension (:), allocatable :: pflux_avg, vflux_avg, qflux_avg, heat_avg

  ! needed for calculating growth rates and frequencies
  complex, dimension (:,:,:), allocatable :: omega_vs_time

  integer :: nout = 1
  logical :: diagnostics_initialized = .false.

  logical :: debug = .false.

contains

  subroutine init_stella_diagnostics (restart, nstep, tstart)

    use zgrid, only: init_zgrid
    use kt_grids, only: init_kt_grids
    use physics_parameters, only: init_physics_parameters
    use stella_geometry, only: geo_surf, twist_and_shift_geo_fac, q_as_x
    use run_parameters, only: init_run_parameters
    use species, only: init_species
    use dist_fn, only: init_dist_fn
    use init_g, only: init_init_g
    use stella_io, only: init_stella_io, get_nout
    use mp, only: broadcast, proc0

    implicit none

    integer, intent (in) :: nstep
    logical, intent (in) :: restart
    real, intent (in) :: tstart

    integer :: nmovie_tot

    if (diagnostics_initialized) return
    diagnostics_initialized = .true.
    
    debug = debug .and. proc0
    
    call init_zgrid
    call init_physics_parameters
    call init_kt_grids (geo_surf, twist_and_shift_geo_fac, q_as_x)
    call init_run_parameters
    call init_species
    call init_init_g
    call init_dist_fn
    
    call read_parameters
    call allocate_arrays
    
    call broadcast (nwrite)
    call broadcast (navg)
    call broadcast (nmovie)
    call broadcast (nsave)
    call broadcast (save_for_restart)
    call broadcast (write_omega)
    call broadcast (write_kspectra)
    call broadcast (write_moments)
    call broadcast (write_phi_vs_time)
    call broadcast (write_gvmus)
    call broadcast (write_gzvs)
!    call broadcast (write_symmetry)
    
    nmovie_tot = nstep/nmovie
    
    call init_averages
    call init_stella_io (restart, write_phi_vs_time, write_kspectra, &
!         write_gvmus, write_gzvs, write_symmetry, write_moments)
         write_gvmus, write_gzvs, write_moments)
    call open_loop_ascii_files(restart)

    if(proc0) call get_nout(tstart,nout)
    call broadcast (nout)
  end subroutine init_stella_diagnostics
  
  subroutine read_parameters

    use mp, only: proc0
    use file_utils, only: input_unit_exist
    use zgrid, only: nperiod, nzed

    implicit none

    logical :: exist
    integer :: in_file

    namelist /stella_diagnostics_knobs/ nwrite, navg, nmovie, nsave, &
         save_for_restart, write_phi_vs_time, write_gvmus, write_gzvs, &
!         write_omega, write_kspectra, write_symmetry, write_moments
         write_omega, write_kspectra, write_moments

    if (proc0) then
       nwrite = 50
       navg = 50
       nmovie = 10000
       nsave = -1
       save_for_restart = .false.
       write_omega = .false.
       write_phi_vs_time = .false.
       write_gvmus = .false.
       write_gzvs = .false.
       write_kspectra = .false.
       write_moments = .false.
!       write_symmetry = .false.

       in_file = input_unit_exist ("stella_diagnostics_knobs", exist)
       if (exist) read (unit=in_file, nml=stella_diagnostics_knobs)

       if (.not. save_for_restart) nsave = -1
    end if
    ntg_out = nzed/2 + (nperiod-1)*nzed

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
    if (.not.allocated(omega_vs_time)) then
       if (write_omega) then
          allocate (omega_vs_time(navg,naky,nakx))
          omega_vs_time = 0.
       else
          allocate (omega_vs_time(1,1,1))
       end if
    end if

  end subroutine allocate_arrays

  subroutine init_averages

    use kt_grids, only: aky, naky

    implicit none

    if (.not.allocated(fac)) then
       allocate (fac(naky)) ; fac = 2.0
       if (aky(1)<epsilon(0.)) fac(1) = 1.0
    end if

  end subroutine init_averages

  subroutine open_loop_ascii_files(restart)

    use file_utils, only: open_output_file
    use species, only: nspec

    implicit none

    logical, intent (in) :: restart
    character (3) :: nspec_str
    character (100) :: str

    logical :: overwrite

    overwrite = .not.restart

    call open_output_file (stdout_unit,'.out',overwrite)
    call open_output_file (fluxes_unit,'.fluxes',overwrite)
    write (nspec_str,'(i3)') nspec*12
    str = trim('(2a12,2a'//trim(nspec_str)//')')
    write (fluxes_unit,str) '#time', 'pflx', 'vflx', 'qflx'
    if (write_omega) then
       call open_output_file (omega_unit,'.omega',overwrite)
       write (omega_unit,'(7a12)') '#time', 'ky', 'kx', &
            'Re[om]', 'Im[om]', 'Re[omavg]', 'Im[omavg]'
    end if

  end subroutine open_loop_ascii_files

  subroutine close_loop_ascii_files
    
    use file_utils, only: close_output_file
    
    implicit none
    
    call close_output_file (stdout_unit)
    call close_output_file (fluxes_unit)
    if (write_omega) call close_output_file (omega_unit)

  end subroutine close_loop_ascii_files

  subroutine diagnose_stella (istep)

    use mp, only: proc0,job
    use constants, only: zi
    use redistribute, only: scatter
    use fields_arrays, only: phi, apar
    use fields_arrays, only: phi_old, phi_corr_QN
    use dist_fn_arrays, only: gvmu, gnew
!    use g_tofrom_h, only: g_to_h
    use stella_io, only: write_time_nc
    use stella_io, only: write_phi2_nc
    use stella_io, only: write_phi_nc
    use stella_io, only: write_gvmus_nc
    use stella_io, only: write_gzvs_nc
    use stella_io, only: write_kspectra_nc
    use stella_io, only: write_moments_nc
    use stella_io, only: sync_nc
!    use stella_io, only: write_symmetry_nc
    use stella_time, only: code_time, code_dt
    use run_parameters, only: fphi
    use zgrid, only: nztot, nzgrid, ntubes
    use vpamu_grids, only: nmu, nvpa
    use species, only: nspec
    use kt_grids, only: naky, nakx, nx
    use dist_redistribute, only: kxkyz2vmu
    use physics_flags, only: radial_variation

    implicit none

    integer, intent (in) :: istep
    
    real :: phi2, apar2
    real :: zero
    real, dimension (:,:,:), allocatable :: gvmus
    real, dimension (:,:,:,:), allocatable :: gzvs
!    real, dimension (:,:,:), allocatable :: pflx_zvpa, vflx_zvpa, qflx_zvpa
    real, dimension (:), allocatable :: part_flux, mom_flux, heat_flux
    real, dimension (:,:), allocatable :: phi2_vs_kxky
    complex, dimension (:,:,:,:,:), allocatable :: density, upar, temperature

    complex, dimension (:,:), allocatable :: omega_avg
    complex, dimension (:,:), allocatable :: phiavg, phioldavg
    complex, dimension (:,:,:,:), allocatable :: phi_out

    ! calculation of omega requires computation of omega more
    ! frequently than every nwrite time steps
    if (write_omega .and. proc0) then
       zero = 100.*epsilon(0.)
       if (istep > 0) then
          allocate (phiavg(naky,nakx))
          allocate (phioldavg(naky,nakx))
          call fieldline_average (phi, phiavg)
          call fieldline_average (phi_old, phioldavg)
          where (abs(phiavg) < zero .or. abs(phioldavg) < zero)
             omega_vs_time(mod(istep,navg)+1,:,:) = 0.0
          elsewhere
             omega_vs_time(mod(istep,navg)+1,:,:) = log(phiavg/phioldavg)*zi/code_dt
          end where
          deallocate (phiavg, phioldavg)
       end if
    end if

    ! only write data to file every nwrite time steps
    if (mod(istep,nwrite) /= 0) return

    allocate(phi_out(naky,nakx,-nzgrid:nzgrid,ntubes))
    phi_out = phi
    if(radial_variation) then
      phi_out = phi_out + phi_corr_QN
    endif

    allocate (part_flux(nspec))
    allocate (mom_flux(nspec))
    allocate (heat_flux(nspec))

    ! obtain turbulent fluxes
    if(radial_variation) then
      call get_fluxes_vmulo (gnew, phi_out, part_flux, mom_flux, heat_flux)
    else
      call scatter (kxkyz2vmu, gnew, gvmu)
!     call g_to_h (gvmu, phi, fphi)
      call get_fluxes (gvmu, part_flux, mom_flux, heat_flux)
!     call g_to_h (gvmu, phi, -fphi)
    endif

    if (proc0) then
       if(write_omega) then
         allocate (omega_avg(naky,nakx))
         omega_avg = sum(omega_vs_time,dim=1)/real(navg)
       else
         allocate (omega_avg(1,1))
       endif
       call volume_average (phi_out, phi2)
       call volume_average (apar, apar2)
       write (*,'(a7,i7,a6,e12.4,a4,e12.4,a10,e12.4,a11,e12.4,a6,i3)') 'istep=', istep, &
            'time=', code_time, 'dt=', code_dt, '|phi|^2=', phi2, '|apar|^2= ', apar2, "job=", job
       call write_loop_ascii_files (istep, phi2, apar2, part_flux, mom_flux, heat_flux, &
            omega_vs_time(mod(istep,navg)+1,:,:), omega_avg)

       ! do not need omega_avg again this time step
       deallocate (omega_avg)
    end if


    if (proc0) then
       if (debug) write (*,*) 'stella_diagnostics::write_time_nc'
       call write_time_nc (nout, code_time)
       call write_phi2_nc (nout, phi2)
       if (write_phi_vs_time) then
          if (debug) write (*,*) 'stella_diagnostics::diagnose_stella::write_phi_nc'
          call write_phi_nc (nout, phi_out)
       end if
       if (write_kspectra) then
          if (debug) write (*,*) 'stella_diagnostics::diagnose_stella::write_kspectra'
          allocate (phi2_vs_kxky(naky,nakx))
          call fieldline_average (real(phi_out*conjg(phi_out)),phi2_vs_kxky)
          call write_kspectra_nc (nout, phi2_vs_kxky)
          deallocate (phi2_vs_kxky)
       end if
    end if
    if (write_moments) then
       if (debug) write (*,*) 'stella_diagnostics::diagnose_stella::write_moments'
       allocate (density(naky,nakx,nztot,ntubes,nspec))
       allocate (upar(naky,nakx,nztot,ntubes,nspec))
       allocate (temperature(naky,nakx,nztot,ntubes,nspec))
       call get_moments (gnew, density, upar, temperature)
       if (proc0) call write_moments_nc (nout, density, upar, temperature)
       deallocate (density, upar, temperature)
    end if
    if (write_gvmus) then
       allocate (gvmus(nvpa,nmu,nspec))
       if (debug) write (*,*) 'stella_diagnostics::diagnose_stella::get_gvmus'
       ! note that gvmus is h at this point
       call get_gvmus (gvmu, gvmus)
       if (debug) write (*,*) 'stella_diagnostics::diagnose_stella::write_gvmus_nc'
       if (proc0) call write_gvmus_nc (nout, gvmus)
       deallocate (gvmus)
    end if
    if (write_gzvs) then
       allocate (gzvs(ntubes,nztot,nvpa,nspec))
       if (debug) write (*,*) 'stella_diagnostics::diagnose_stella::get_gzvs'
       call get_gzvs (gnew, gzvs)
       if (debug) write (*,*) 'stella_diagnostics::diagnose_stella::write_gzvs_nc'
       if (proc0) call write_gzvs_nc (nout, gzvs)
       deallocate (gzvs)
    end if
!     if (write_symmetry) then
!        allocate (pflx_zvpa(nztot,nvpa,nspec))
!        allocate (vflx_zvpa(nztot,nvpa,nspec))
!        allocate (qflx_zvpa(nztot,nvpa,nspec))
!        call get_fluxes_vs_zvpa (gnew, pflx_zvpa, vflx_zvpa, qflx_zvpa)
!        deallocate (pflx_zvpa, vflx_zvpa, qflx_zvpa)
!     end if

    if (proc0) call sync_nc

    deallocate (part_flux, mom_flux, heat_flux)
    deallocate(phi_out)

    nout = nout + 1

  end subroutine diagnose_stella

  subroutine fieldline_average_real (unavg, avg)

    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky
    use stella_geometry, only: dl_over_b

    implicit none

    real, dimension (:,:,-nzgrid:,:), intent (in) :: unavg
    real, dimension (:,:), intent (out) :: avg

    integer :: it, ia

    ia = 1

    avg = 0.0
    do it = 1, ntubes
       avg = avg + sum(spread(spread(dl_over_b(ia,:),1,naky),2,nakx)*unavg(:,:,:,it),dim=3)
    end do
    avg = avg/real(ntubes)
    
  end subroutine fieldline_average_real

  subroutine fieldline_average_complex (unavg, avg)

    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky
    use stella_geometry, only: dl_over_b

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: unavg
    complex, dimension (:,:), intent (out) :: avg

    integer :: it, ia

    ia = 1

    avg = 0.0
    do it = 1, ntubes
       avg = avg + sum(spread(spread(dl_over_b(ia,:),1,naky),2,nakx)*unavg(:,:,:,it),dim=3)
    end do
    avg = avg/real(ntubes)

  end subroutine fieldline_average_complex

  subroutine volume_average (unavg, avg)

    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx
    use stella_geometry, only: dl_over_b

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: unavg
    real, intent (out) :: avg

    integer :: iky, ikx, iz, it, ia

    ia = 1

    avg = 0.
    do it = 1, ntubes
       do iz = -nzgrid, nzgrid
          do ikx = 1, nakx
             do iky = 1, naky
                avg = avg + real(unavg(iky,ikx,iz,it)*conjg(unavg(iky,ikx,iz,it)))*fac(iky)*dl_over_b(ia,iz)
             end do
          end do
       end do
    end do
    avg = avg/real(ntubes)

  end subroutine volume_average

  ! assumes that the non-Boltzmann part of df is passed in (aka h)
  subroutine get_fluxes (g, pflx, vflx, qflx)

    use mp, only: sum_reduce
    use constants, only: zi
    use fields_arrays, only: phi, apar
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use species, only: spec
    use stella_geometry, only: jacob, grho, bmag, btor
    use stella_geometry, only: gds21, gds22
    use stella_geometry, only: geo_surf
    use zgrid, only: delzed, nzgrid, ntubes
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: vperp2, vpa
    use run_parameters, only: fphi, fapar
    use kt_grids, only: aky, theta0
    use gyro_averages, only: gyro_average, gyro_average_j1

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (out) :: pflx, vflx, qflx

    integer :: ikxkyz, iky, ikx, iz, it, is, ia
    real, dimension (:), allocatable :: flx_norm
    complex, dimension (:,:), allocatable :: gtmp1, gtmp2, gtmp3

    allocate (flx_norm(-nzgrid:nzgrid))
    allocate (gtmp1(nvpa,nmu), gtmp2(nvpa,nmu), gtmp3(nvpa,nmu))

    pflx = 0. ; vflx = 0. ; qflx = 0.

    flx_norm = jacob(1,:)*delzed
    flx_norm = flx_norm/sum(flx_norm*grho(1,:))

    ia = 1
    ! get electrostatic contributions to fluxes
    if (fphi > epsilon(0.0)) then
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          
          ! get particle flux
          call gyro_average (g(:,:,ikxkyz), ikxkyz, gtmp1)
          call get_one_flux (iky, iz, flx_norm(iz), gtmp1, phi(iky,ikx,iz,it), pflx(is))

          ! get heat flux
          ! NEEDS TO BE MODIFIED TO TREAT ENERGY = ENERGY(ALPHA)
          gtmp1 = gtmp1*(spread(vpa**2,2,nmu)+spread(vperp2(1,iz,:),1,nvpa))
          call get_one_flux (iky, iz, flx_norm(iz), gtmp1, phi(iky,ikx,iz,it), qflx(is))

          ! get momentum flux
          ! parallel component
          gtmp1 = g(:,:,ikxkyz)*spread(vpa,2,nmu)*geo_surf%rmaj*btor(iz)/bmag(ia,iz)
          call gyro_average (gtmp1, ikxkyz, gtmp2)
          gtmp1 = -g(:,:,ikxkyz)*zi*aky(iky)*spread(vperp2(ia,iz,:),1,nvpa)*geo_surf%rhoc &
               * (gds21(ia,iz)+theta0(iky,ikx)*gds22(ia,iz))*spec(is)%smz &
               / (geo_surf%qinp*geo_surf%shat*bmag(ia,iz)**2)
          call gyro_average_j1 (gtmp1, ikxkyz, gtmp3)
          gtmp1 = gtmp2 + gtmp3

          call get_one_flux (iky, iz, flx_norm(iz), gtmp1, phi(iky,ikx,iz,it), vflx(is))
       end do
    end if

    if (fapar > epsilon(0.0)) then
       ! particle flux
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          
          ! Apar contribution to particle flux
          gtmp1 = -g(:,:,ikxkyz)*spec(is)%stm*spread(vpa,2,nmu)
          call gyro_average (gtmp1, ikxkyz, gtmp2)
          call get_one_flux (iky, iz, flx_norm(iz), gtmp2, apar(iky,ikx,iz,it), pflx(is))
          
          ! Apar contribution to heat flux
          gtmp2 = gtmp2*(spread(vpa**2,2,nmu)+spread(vperp2(ia,iz,:),1,nvpa))
          call get_one_flux (iky, iz, flx_norm(iz), gtmp2, apar(iky,ikx,iz,it), qflx(is))
          
          ! Apar contribution to momentum flux
          ! parallel component
          gtmp1 = -spread(vpa**2,2,nmu)*spec(is)%stm*g(:,:,ikxkyz) &
               * geo_surf%rmaj*btor(iz)/bmag(1,iz)
          call gyro_average (gtmp1, ikxkyz, gtmp2)
          ! perp component
          gtmp1 = spread(vpa,2,nmu)*spec(is)%stm*g(:,:,ikxkyz) &
               * zi*aky(iky)*spread(vperp2(ia,iz,:),1,nvpa)*geo_surf%rhoc &
               * (gds21(ia,iz)+theta0(iky,ikx)*gds22(ia,iz))*spec(is)%smz &
               / (geo_surf%qinp*geo_surf%shat*bmag(ia,iz)**2)
          call gyro_average_j1 (gtmp1, ikxkyz, gtmp3)
          ! FLAG -- NEED TO ADD IN CONTRIBUTION FROM BOLTZMANN PIECE !!

          gtmp1 = gtmp2 + gtmp3

          call get_one_flux (iky, iz, flx_norm(iz), gtmp1, apar(iky,ikx,iz,it), vflx(is))
       end do
    end if

    call sum_reduce (pflx, 0) ; pflx = pflx*spec%dens_psi0
    call sum_reduce (qflx, 0) ; qflx = qflx*spec%dens_psi0*spec%temp_psi0
    call sum_reduce (vflx, 0) ; vflx = vflx*spec%dens_psi0*sqrt(spec%mass*spec%temp_psi0)

    ! normalise to account for contributions from multiple flux tubes
    ! in flux tube train
    pflx = pflx/real(ntubes)
    qflx = qflx/real(ntubes)
    vflx = vflx/real(ntubes)

    deallocate (gtmp1, gtmp2, gtmp3)
    deallocate (flx_norm)

  end subroutine get_fluxes

  subroutine get_fluxes_vmulo (g, phi, pflx, vflx, qflx)

    use mp, only: sum_reduce
    use constants, only: zi
    use dist_fn_arrays, only: g1, g2, kperp2, dkperp2dr
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use species, only: spec
    use stella_geometry, only: jacob, grho, bmag, btor
    use stella_geometry, only: drhodpsi
    use stella_geometry, only: gds21, gds22
    use stella_geometry, only: dgds21dr, dgds22dr
    use stella_geometry, only: geo_surf
    use stella_geometry, only: dBdrho, dIdrho, rho_to_x
    use stella_geometry, only: dl_over_b, d_dl_over_b_drho
    use zgrid, only: delzed, nzgrid, ntubes
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: vperp2, vpa, mu
    use run_parameters, only: fphi, fapar
    use kt_grids, only: aky, theta0, naky, nakx, nx, x_clamped
    use physics_flags, only: radial_variation
    use gyro_averages, only: gyro_average, gyro_average_j1, aj0x, aj1x
    use stella_transforms, only: transform_x2kx_xfirst, transform_kx2x_xfirst

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phi
    real, dimension (:), intent (out) :: pflx, vflx, qflx

    integer :: ivmu, imu, iv, iz, it, is, ia
    real, dimension (:), allocatable :: flx_norm
    complex, dimension (:,:), allocatable :: gtmp1, gtmp2, gtmp3
    complex, dimension (:,:), allocatable :: g0k, g0x

    allocate (flx_norm(-nzgrid:nzgrid))
    allocate (gtmp1(nvpa,nmu), gtmp2(nvpa,nmu), gtmp3(nvpa,nmu))

    pflx = 0. ; vflx = 0. ; qflx = 0.

    flx_norm = jacob(1,:)*delzed
    flx_norm = flx_norm/sum(flx_norm*grho(1,:))

    if(radial_variation) then
      allocate (g0k(naky,nakx))
      allocate (g0x(naky,nx))
    endif

    ia = 1
    ! FLAG - electrostatic for now
    ! get electrostatic contributions to fluxes

    if (fphi > epsilon(0.0)) then
       ia = 1

       !get particle flux
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          iv = iv_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          is = is_idx(vmu_lo,ivmu)

          call gyro_average (g(:,:,:,:,ivmu), ivmu, g1(:,:,:,:,ivmu))

          if(radial_variation) then
            do it = 1, ntubes
              do iz= -nzgrid, nzgrid
                g0k = g1(:,:,iz,it,ivmu) &
                  * (-0.5*aj1x(:,:,iz,ivmu)/aj0x(:,:,iz,ivmu)*(spec(is)%smz)**2 & 
                  * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                  * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                  + dBdrho(iz)/bmag(ia,iz) + d_dl_over_b_drho(ia,iz)/dl_over_b(ia,iz))
               
                call transform_kx2x_xfirst (g0k,g0x)
                g0x = rho_to_x*spread(x_clamped,1,naky)*g0x
                call transform_x2kx_xfirst (g0x,g0k)

                g1(:,:,iz,it,ivmu) = g1(:,:,iz,it,ivmu) + g0k
              enddo
            enddo
          endif
       enddo
       call get_one_flux_vmulo (flx_norm, spec%dens_psi0, g1, phi, pflx)

       !get heat flux
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          iv = iv_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          is = is_idx(vmu_lo,ivmu)

          call gyro_average (g(:,:,:,:,ivmu), ivmu, g1(:,:,:,:,ivmu))
          do it = 1, ntubes
            do iz= -nzgrid, nzgrid

              g1(:,:,iz,it,ivmu) = g1(:,:,iz,it,ivmu)*(vpa(iv)**2+vperp2(ia,iz,imu))

              if(radial_variation) then
                g0k = g1(:,:,iz,it,ivmu) & 
                  * (-0.5*aj1x(:,:,iz,ivmu)/aj0x(:,:,iz,ivmu)*(spec(is)%smz)**2 & 
                     * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                     * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                     + dBdrho(iz)/bmag(ia,iz) &
                     + 2.0*mu(imu)*dBdrho(iz)/(vpa(iv)**2+vperp2(ia,iz,imu)) &
                     + d_dl_over_b_drho(ia,iz)/dl_over_b(ia,iz))

                call transform_kx2x_xfirst (g0k,g0x)
                g0x = rho_to_x*spread(x_clamped,1,naky)*g0x
                call transform_x2kx_xfirst (g0x,g0k)

                g1(:,:,iz,it,ivmu) = g1(:,:,iz,it,ivmu) + g0k
              endif
            enddo
          enddo
       enddo
       call get_one_flux_vmulo (flx_norm,spec%dens_psi0*spec%temp_psi0, g1, phi, qflx)

       ! get momentum flux
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          iv = iv_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          is = is_idx(vmu_lo,ivmu)
          do it = 1, ntubes
            do iz= -nzgrid, nzgrid
            ! parallel component
              g0k = g(:,:,iz,it,ivmu)*vpa(iv)*geo_surf%rmaj*btor(iz)/bmag(ia,iz)
              call gyro_average (g0k, iz, ivmu, g1(:,:,iz,it,ivmu))

              if(radial_variation) then
                g0k = g1(:,:,iz,it,ivmu) &
                  * (-0.5*aj1x(:,:,iz,ivmu)/aj0x(:,:,iz,ivmu)*(spec(is)%smz)**2 & 
                  * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                  * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                  + dIdrho/(geo_surf%rmaj*btor(iz)) & 
                  + d_dl_over_b_drho(ia,iz)/dl_over_b(ia,iz))

                call transform_kx2x_xfirst (g0k,g0x)
                g0x = rho_to_x*spread(x_clamped,1,naky)*g0x
                call transform_x2kx_xfirst (g0x,g0k)

                g1(:,:,iz,it,ivmu) = g1(:,:,iz,it,ivmu) + g0k

              endif

              ! perpendicular component
              g0k = -g(:,:,iz,it,ivmu)*zi*spread(aky,2,nakx)*vperp2(ia,iz,imu)*geo_surf%rhoc &
                * (gds21(ia,iz)+theta0*gds22(ia,iz))*spec(is)%smz &
                / (geo_surf%qinp*geo_surf%shat*bmag(ia,iz)**2)

              call gyro_average_j1 (g0k, iz, ivmu, g2(:,:,iz,it,ivmu))
              if(radial_variation) then
                g0k = g2(:,:,iz,it,ivmu) &
                    * ( (dgds21dr(ia,iz)+theta0*dgds22dr(ia,iz))/(gds21(ia,iz)+theta0*gds22(ia,iz)) &
                       - geo_surf%d2qdr2*geo_surf%rhoc/(geo_surf%shat*geo_surf%qinp) & 
                       - geo_surf%d2psidr2*drhodpsi &
                       + (0.5*aj0x(:,:,iz,ivmu)/aj1x(:,:,iz,ivmu) - 1.0) & 
                       * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                       + d_dl_over_b_drho(ia,iz)/dl_over_b(ia,iz))

                call transform_kx2x_xfirst (g0k,g0x)
                g0x = rho_to_x*spread(x_clamped,1,naky)*g0x
                call transform_x2kx_xfirst (g0x,g0k)

                g2(:,:,iz,it,ivmu) = g2(:,:,iz,it,ivmu) + g0k
              endif
          enddo
        enddo
      enddo

      g1 = g1 + g2
      call get_one_flux_vmulo (flx_norm,spec%dens_psi0*sqrt(spec%mass*spec%temp_psi0), g1, phi, vflx)

    end if

    ! normalise to account for contributions from multiple flux tubes
    ! in flux tube train
    pflx = pflx/real(ntubes)
    qflx = qflx/real(ntubes)
    vflx = vflx/real(ntubes)

    deallocate (flx_norm)

  end subroutine get_fluxes_vmulo

  subroutine get_one_flux (iky, iz, norm, gin, fld, flxout)

    use vpamu_grids, only: integrate_vmu
    use kt_grids, only: aky

    implicit none

    integer, intent (in) :: iky, iz
    real, intent (in) :: norm
    complex, dimension (:,:), intent (in) :: gin
    complex, intent (in) :: fld
    real, intent (in out) :: flxout
    
    complex :: flx

    call integrate_vmu (gin,iz,flx)
    flxout = flxout &
         + 0.5*fac(iky)*aky(iky)*aimag(flx*conjg(fld))*norm

  end subroutine get_one_flux

  subroutine get_one_flux_vmulo (norm, weights, gin, fld, flxout)

    use vpamu_grids, only: integrate_vmu
    use stella_layouts, only: vmu_lo
    use kt_grids, only: aky, nakx, naky
    use zgrid, only: nzgrid, ntubes
    use species, only: nspec

    implicit none

    real, dimension (-nzgrid:), intent (in) :: norm
    real, dimension (:), intent (in) :: weights
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: gin
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: fld
    real, dimension (:), intent (in out) :: flxout
    
    complex, dimension (:,:,:,:,:), allocatable :: totals

    integer :: is, it, iz, ikx

    allocate (totals(naky,nakx,-nzgrid:nzgrid,ntubes,nspec))

    call integrate_vmu (gin,weights,totals)
    do is = 1, nspec
      do it = 1, ntubes
        do iz= -nzgrid, nzgrid
          do ikx = 1, nakx
            flxout(is) = flxout(is) &
              + sum(0.5*fac*aky*aimag(totals(:,ikx,iz,it,is)*conjg(fld(:,ikx,iz,it)))*norm(iz))
          enddo
        enddo
      enddo
    enddo

    deallocate (totals)

  end subroutine get_one_flux_vmulo

!   subroutine get_fluxes_vs_zvpa (g, pflx, vflx, qflx)

!     use zgrid, only: nzgrid, delzed
!     use stella_layouts, only: vmu_lo
!     use stella_geometry, only: jacob, grho

!     implicit none

!     complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
!     real, dimension (-nzgrid:,:,:), intent (out) :: pflx, vflx, qflx

!     real, dimension (:), allocatable :: flx_norm
! !    real, dimension (:,:), allocatable :: gtmp

! !    allocate (gtmp(-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

!     allocate (flx_norm(-nzgrid:nzgrid))

! !    gtmp = 0.

!     pflx = 0. ; vflx = 0. ; qflx = 0.

!     flx_norm = jacob(1,:)*delzed
!     flx_norm = flx_norm/sum(flx_norm*grho(1,:))

! !     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
! !        do ikx = 1, nakx
! !           do iky = 1, naky
! !        = g(:,:,:,ivmu)*aj0x(:,:,:,ivmu)
! !        integrate_mu (,pflx)
! !        pflx(iky,ikx,iz) = pflx + 0.5*fac(iky)*aky(iky)*aimag(pflx*conjg(phi))*flx_norm(iz)
! !     end do

! !    deallocate (gtmp)
!     deallocate (flx_norm)

!   end subroutine get_fluxes_vs_zvpa

  subroutine get_moments (g, dens, upar, temp)
    
    use zgrid, only: nzgrid, ntubes
    use species, only: spec
    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: vpa, vperp2, mu
    use vpamu_grids, only: maxwell_mu, ztmax, maxwell_fac
    use kt_grids, only: naky, nakx, nx, x_clamped
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use dist_fn_arrays, only: g1, g2, kperp2, dkperp2dr
    use stella_geometry, only: bmag, dBdrho, rho_to_x
    use gyro_averages, only: aj0x, aj1x, gyro_average
    use fields_arrays, only: phi, phi_corr_QN
    use physics_flags, only: radial_variation
    use stella_transforms, only: transform_x2kx_xfirst, transform_kx2x_xfirst

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,:,:,:), intent (out) :: dens, upar, temp

    complex, dimension (:,:), allocatable :: g0k, g0x

    integer :: ivmu, iv, imu, is, ia
    integer :: iz, it

    if(radial_variation) then
      allocate (g0k(naky,nakx))
      allocate (g0x(naky,nx))
    endif

    ! f = h - Ze*phi/T * F0
    ! g = h - Ze*<phi>/T * F0
    ! f = g + Ze*(<phi>-phi)/T * F0
    ia = 1
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       call gyro_average (g(:,:,:,:,ivmu), ivmu, g1(:,:,:,:,ivmu))
       ! FLAG -- AJ0X NEEDS DEALING WITH BELOW
       g2(:,:,:,:,ivmu) = g1(:,:,:,:,ivmu) + ztmax(iv,is) &
            * spread(spread(spread(maxwell_mu(ia,:,imu,is),1,naky),2,nakx) &
            * maxwell_fac(is)*(aj0x(:,:,:,ivmu)**2-1.0),4,ntubes)*phi

       if(radial_variation) then
         do it = 1, ntubes
           do iz= -nzgrid, nzgrid
             !phi
             g0k = ztmax(iv,is)*maxwell_mu(ia,iz,imu,is) &
               * maxwell_fac(is)*(aj0x(:,:,iz,ivmu)**2-1.0)*phi(:,:,iz,it) &
               *(-spec(is)%tprim*(vpa(iv)**2+vperp2(ia,iz,imu)-2.5) &
                 -spec(is)%fprim+(dBdrho(iz)/bmag(ia,iz))*(1.0 - 2.0*mu(imu)*bmag(ia,iz)) &
                 -aj1x(:,:,iz,ivmu)*aj0x(:,:,iz,ivmu)*(spec(is)%smz)**2 & 
                 * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                 * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                    /(aj0x(:,:,iz,ivmu)**2 - 1.0))

             !g
             g0k = g0k + g1(:,:,iz,it,ivmu) &
               * (-0.5*aj1x(:,:,iz,ivmu)/aj0x(:,:,iz,ivmu)*(spec(is)%smz)**2 & 
               * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
               * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
               + dBdrho(iz)/bmag(ia,iz))

             call transform_kx2x_xfirst (g0k,g0x)
             g0x = rho_to_x*spread(x_clamped,1,naky)*g0x
             call transform_x2kx_xfirst (g0x,g0k)

             !phi QN
             g0k = g0k + ztmax(iv,is)*maxwell_mu(ia,iz,imu,is) &
               * maxwell_fac(is)*(aj0x(:,:,iz,ivmu)**2-1.0)*phi_corr_QN(:,:,iz,it)

             g2(:,:,iz,it,ivmu) = g2(:,:,iz,it,ivmu) + g0k
           enddo
         enddo
       endif
    end do
    call integrate_vmu (g2, spec%dens, dens)

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       g2(:,:,:,:,ivmu) = (g1(:,:,:,:,ivmu) + ztmax(iv,is) &
            * spread(spread(spread(maxwell_mu(ia,:,imu,is),1,naky),2,nakx) &
            * maxwell_fac(is)*(aj0x(:,:,:,ivmu)**2-1.0),4,ntubes)*phi) &
            *(vpa(iv)**2+spread(spread(spread(vperp2(1,:,imu),1,naky),2,nakx),4,ntubes))/1.5
       if(radial_variation) then
         do it = 1, ntubes
           do iz= -nzgrid, nzgrid
             !phi
             g0k = ztmax(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is) &
               *(vpa(iv)**2 + vperp2(ia,iz,imu))/1.5 & 
               *(aj0x(:,:,iz,ivmu)**2-1.0)*phi(:,:,iz,it) &
               *(-spec(is)%tprim*(vpa(iv)**2+vperp2(ia,iz,imu)-2.5) &
                 -spec(is)%fprim+(dBdrho(iz)/bmag(ia,iz))*(1.0 - 2.0*mu(imu)*bmag(ia,iz)) &
                 -aj1x(:,:,iz,ivmu)*aj0x(:,:,iz,ivmu)*(spec(is)%smz)**2 & 
                 * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                 * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                    /(aj0x(:,:,iz,ivmu)**2 - 1.0) &
                 + 2.0*mu(imu)*dBdrho(iz)/(vpa(iv)**2+vperp2(ia,iz,imu)))


             !g
             g0k = g0k + g1(:,:,iz,it,ivmu)*(vpa(iv)**2+vperp2(ia,iz,imu))/1.5 & 
               * (-0.5*aj1x(:,:,iz,ivmu)/aj0x(:,:,iz,ivmu)*(spec(is)%smz)**2 & 
               * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
               * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                 + dBdrho(iz)/bmag(ia,iz) &
                 + 2.0*mu(imu)*dBdrho(iz)/(vpa(iv)**2+vperp2(ia,iz,imu)))

             call transform_kx2x_xfirst (g0k,g0x)
             g0x = rho_to_x*spread(x_clamped,1,naky)*g0x
             call transform_x2kx_xfirst (g0x,g0k)

             !phi QN
             g0k = g0k + ztmax(iv,is)*maxwell_mu(ia,iz,imu,is) &
               * (vpa(iv)**2 + vperp2(ia,iz,imu))/1.5 & 
               * maxwell_fac(is)*(aj0x(:,:,iz,ivmu)**2-1.0)*phi_corr_QN(:,:,iz,it)

             g2(:,:,iz,it,ivmu) = g2(:,:,iz,it,ivmu) + g0k
           enddo
         enddo
       endif
    end do
    ! integrate to get dTs/Tr
!    call integrate_vmu (g2, spec%temp, temp)
    call integrate_vmu (g2, spec%temp_psi0*spec%dens, temp)

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       g2(:,:,:,:,ivmu) = vpa(iv)*g1(:,:,:,:,ivmu)
       if(radial_variation) then
         do it = 1, ntubes
           do iz= -nzgrid, nzgrid
             !g
             g0k = vpa(iv)*g1(:,:,iz,it,ivmu) &
               * (-0.5*aj1x(:,:,iz,ivmu)/aj0x(:,:,iz,ivmu)*(spec(is)%smz)**2 & 
               * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
               * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
               + dBdrho(iz)/bmag(ia,iz))

             call transform_kx2x_xfirst (g0k,g0x)
             g0x = rho_to_x*spread(x_clamped,1,naky)*g0x
             call transform_x2kx_xfirst (g0x,g0k)

             g2(:,:,iz,it,ivmu) = g2(:,:,iz,it,ivmu) + g0k
           enddo
         enddo
       endif
    end do
    call integrate_vmu (g2, spec%stm_psi0, upar)

    if(allocated(g0k)) deallocate(g0k)
    if(allocated(g0x)) deallocate(g0x)


  end subroutine get_moments

  ! get_gvmus takes g(kx,ky,z) and returns average over z of int dxdy g(x,y,z)^2
  ! SHOULD MODIFY TO TAKE ADVANTAGE OF FACT THAT G(KY,KX,Z) LOCAL IS AVAILABLE
  subroutine get_gvmus (g, gv)

    use mp, only: nproc, sum_reduce
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: is_idx, iky_idx, iz_idx
    use zgrid, only: ntubes
    use vpamu_grids, only: nvpa, nmu
    use stella_geometry, only: dl_over_b

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    real, dimension (:,:,:), intent (out) :: gv

    integer :: ikxkyz, iv, is, imu, iz, iky, ia!, ivp

    ia = 1

    ! when doing volume averages, note the following:
    ! int dxdy g(x,y)^2 = sum_ky |g(ky=0,kx)|^2 + 2 * sum_{kx,ky} |g(ky>0,kx)|^2
    ! factor of 2 accounted for in fac

    gv = 0.
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       is = is_idx(kxkyz_lo,ikxkyz)
       iky = iky_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          do iv = 1, nvpa
             gv(iv,imu,is) = gv(iv,imu,is) + real(g(iv,imu,ikxkyz)*conjg(g(iv,imu,ikxkyz)))*fac(iky)*dl_over_b(ia,iz)
          end do
       end do
    end do
    gv = gv/real(ntubes)

    if (nproc > 1) call sum_reduce (gv,0)

  end subroutine get_gvmus

  ! get_gzvs takes g(kx,ky,z,vpa,mu,s) and returns int dmudxdy g(x,y,z,vpa,mu,s)^2
  subroutine get_gzvs (g, gz)

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: integrate_mu
    use kt_grids, only: nakx, naky

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (:,:,:,:), intent (out) :: gz

    integer :: ivmu, iz, it, ikx, iky, izp

    real, dimension (:,:,:), allocatable :: gtmp

    allocate (gtmp(-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    ! when doing volume averages, note the following:
    ! int dxdy g(x,y)^2 = sum_kx |g(kx,ky=0)|^2 + 2 * sum_{kx,ky} |g(kx,ky>0)|^2
    ! factor of 2 accounted for in fac

    gtmp = 0.
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do ikx = 1, nakx
          do iky = 1, naky
             gtmp(:,:,ivmu) = gtmp(:,:,ivmu) + real(g(iky,ikx,:,:,ivmu)*conjg(g(iky,ikx,:,:,ivmu)))*fac(iky)
          end do
       end do
    end do
    
    do it = 1, ntubes
       do iz = -nzgrid, nzgrid
          izp = iz+nzgrid+1
          call integrate_mu (iz, gtmp(iz,it,:), gz(it,izp,:,:))
       end do
    end do

    deallocate (gtmp)

  end subroutine get_gzvs

  subroutine finish_stella_diagnostics(istep)

    use mp, only: proc0
    use redistribute, only: scatter
    use stella_io, only: finish_stella_io
    use run_parameters, only: fphi, fapar
    use stella_time, only: code_dt, code_time
    use stella_save, only: stella_save_for_restart
    use dist_redistribute, only: kxkyz2vmu
    use dist_fn_arrays, only: gnew, gvmu

    implicit none

    integer :: istatus
    integer, intent (in) :: istep

    if (proc0) then
       call write_final_ascii_files
       call close_loop_ascii_files
    end if
    if (save_for_restart) then
        call scatter (kxkyz2vmu, gnew, gvmu)
        call stella_save_for_restart (gvmu, istep, code_time, code_dt, istatus, fphi, fapar, .true.)
    end if
    call finish_stella_io
    call finish_averages
    call deallocate_arrays

    nout = 1
    diagnostics_initialized = .false.

  end subroutine finish_stella_diagnostics

  subroutine write_loop_ascii_files (istep, phi2, apar2, pflx, vflx, qflx, om, om_avg)

    use stella_time, only: code_time
    use species, only: nspec
    use kt_grids, only: naky, nakx
    use kt_grids, only: aky, akx

    implicit none
    
    integer, intent (in) :: istep
    real, intent (in) :: phi2, apar2
    real, dimension (:), intent (in) :: pflx, vflx, qflx
    complex, dimension (:,:), intent (in) :: om, om_avg

    character (3) :: nspec_str
    character (100) :: str
    integer :: ikx, iky

    write (stdout_unit,'(a7,i7,a6,e12.4,a10,e12.4,a11,e12.4)') 'istep=', istep, &
         'time=', code_time, '|phi|^2=', phi2, '|apar|^2= ', apar2

    call flush(stdout_unit)

    write (nspec_str,'(i3)') 3*nspec+1
    str = trim('('//trim(nspec_str)//'e12.4)')
    write (fluxes_unit,str) code_time, pflx, vflx, qflx

    call flush(stdout_unit)
    call flush(fluxes_unit)

    if (write_omega .and. istep > 0) then
       do iky = 1, naky
          do ikx = 1, nakx
             write (omega_unit,'(7e16.8)') code_time, aky(iky), akx(ikx),&
                  real(om(iky,ikx)), aimag(om(iky,ikx)), &
                  real(om_avg(iky,ikx)), aimag(om_avg(iky,ikx))
          end do
          if (nakx > 1) write (omega_unit,*)
       end do
       if (naky > 1) write (omega_unit,*)
       call flush(omega_unit)
    end if

  end subroutine write_loop_ascii_files

  subroutine write_final_ascii_files

    use file_utils, only: open_output_file, close_output_file
    use fields_arrays, only: phi, apar
    use zgrid, only: nzgrid, ntubes
    use zgrid, only: zed
    use kt_grids, only: naky, nakx
    use kt_grids, only: aky, akx, theta0
    use stella_geometry, only: zed_eqarc

    implicit none

    integer :: tmpunit
    integer :: iky, ikx, iz, it

    call open_output_file (tmpunit,'.final_fields')
    write (tmpunit,'(9a14)') '# z', 'z-thet0', 'aky', 'akx', &
         'real(phi)', 'imag(phi)', 'real(apar)', 'imag(apar)', &
         'z_eqarc-thet0'
    do iky = 1, naky
       do ikx = 1, nakx
          do it = 1, ntubes
             do iz = -nzgrid, nzgrid
                write (tmpunit,'(9es15.4e3,i3)') zed(iz), zed(iz)-theta0(iky,ikx), aky(iky), akx(ikx), &
                  real(phi(iky,ikx,iz,it)), aimag(phi(iky,ikx,iz,it)), &
                  real(apar(iky,ikx,iz,it)), aimag(apar(iky,ikx,iz,it)), zed_eqarc(iz)-theta0(iky,ikx), it
             end do
             write (tmpunit,*)
          end do
       end do
    end do
    call close_output_file (tmpunit)
    
  end subroutine write_final_ascii_files

  subroutine finish_averages

    implicit none

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
    if (allocated(omega_vs_time)) deallocate (omega_vs_time)

  end subroutine deallocate_arrays

end module stella_diagnostics
