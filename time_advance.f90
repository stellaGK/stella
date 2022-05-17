module time_advance

  public :: init_time_advance, finish_time_advance
  public :: advance_stella
  public :: time_gke, time_parallel_nl
  public :: checksum

  private

  interface get_dgdy
     module procedure get_dgdy_2d
     module procedure get_dgdy_3d
     module procedure get_dgdy_4d
  end interface

  interface get_dgdx
     module procedure get_dgdx_2d
     module procedure get_dgdx_3d
     module procedure get_dgdx_4d
  end interface

  interface checksum
     module procedure checksum_field
     module procedure checksum_dist
  end interface

  logical :: time_advance_initialized = .false.

  logical :: wdriftinit   = .false.
  logical :: wstarinit    = .false.
  logical :: parnlinit    = .false.
  logical :: readinit     = .false.
  logical :: radialinit   = .false.
  logical :: driftimpinit = .false.

  ! if .true., dist fn is represented on alpha grid
  ! if .false., dist fn is given on k-alpha grid
  ! default is .false.; will only ever be set to
  ! .true. during full_flux_surface simulations
!  logical :: alpha_space = .false.

  integer :: explicit_option_switch
  integer, parameter :: explicit_option_rk3 = 1, &
       explicit_option_rk2 = 2, &
       explicit_option_rk4 = 3

  real :: xdriftknob, ydriftknob, wstarknob
  logical :: flip_flop
  logical :: leapfrog_this_timestep, runge_kutta_terms_this_timestep

  complex, dimension (:,:,:), allocatable :: gamtot_drifts!, apar_denom_drifts
  complex, dimension (:,:), allocatable :: gamtot3_drifts

  ! factor multiplying parallel nonlinearity
  real, dimension (:,:), allocatable :: par_nl_fac, d_par_nl_fac_dr
  ! factor multiplying higher order linear term in parallel acceleration
  real, dimension (:,:), allocatable :: par_nl_curv, d_par_nl_curv_dr
  real, dimension (:), allocatable :: par_nl_driftx, par_nl_drifty
  real, dimension (:), allocatable :: d_par_nl_driftx_dr, d_par_nl_drifty_dr

  ! needed for timing various pieces of gke solve
  real, dimension (2,10) :: time_gke = 0.
  real, dimension (2,2) :: time_parallel_nl = 0.

  logical :: debug = .false.

contains

  subroutine init_time_advance

    use mp, only: proc0
    use stella_transforms, only: init_transforms
    use run_parameters, only: drifts_implicit
    use physics_flags, only: radial_variation
    use physics_flags, only: include_parallel_nonlinearity
    use neoclassical_terms, only: init_neoclassical_terms
    use dissipation, only: init_collisions, include_collisions
    use parallel_streaming, only: init_parallel_streaming
    use mirror_terms, only: init_mirror
    use flow_shear, only: init_flow_shear

    implicit none

    if (time_advance_initialized) return
    time_advance_initialized = .true.

    debug = debug .and. proc0

    if (debug) write (6,*) 'time_advance::init_time_advance::read_parameters'
    call read_parameters
    if (debug) write (6,*) 'time_advance::init_time_advance::allocate_arrays'
    call allocate_arrays
    if (debug) write (6,*) 'time_advance::init_time_advance::init_neoclassical_terms'
    call init_neoclassical_terms
    if (debug) write (6,*) 'time_advance::init_time_advance::init_mirror'
    call init_mirror
    if (debug) write (6,*) 'time_advance::init_time_advance::init_parstream'
    call init_parallel_streaming
    if (debug) write (6,*) 'time_advance::init_time_advance::init_wdrift'
    call init_wdrift
    if (debug) write (6,*) 'time_advance::init_time_advance::init_wstar'
    call init_wstar
    if (debug) write (6,*) 'time_advance::init_time_advance::init_flow_shear'
    call init_flow_shear
    if (debug) write (6,*) 'time_advance::init_time_advance::init_parallel_nonlinearity'
    if (include_parallel_nonlinearity) call init_parallel_nonlinearity
    if (debug) write (6,*) 'time_advance::init_time_advance::init_radial_variation'
    if (radial_variation) call init_radial_variation
    if (debug) write (6,*) 'time_advance::init_time_advance::init_drifts_implicit'
    if (drifts_implicit) call init_drifts_implicit
    if (include_collisions) then
      if (debug) write (6,*) 'time_advance::init_time_advance::init_collisions'
      call init_collisions
    endif
    if (debug) write (6,*) 'time_advance::init_time_advance::init_cfl'
    call init_cfl

    !call write_drifts

  end subroutine init_time_advance

  subroutine read_parameters

    use file_utils, only: error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    use run_parameters, only: none_implicit

    implicit none

    logical :: taexist

    type (text_option), dimension (4), parameter :: explicitopts = &
         (/ text_option('default', explicit_option_rk3), &
         text_option('rk3', explicit_option_rk3), &
         text_option('rk2', explicit_option_rk2), &
         text_option('rk4', explicit_option_rk4) /)
    character(10) :: explicit_option

    namelist /time_advance_knobs/ xdriftknob, ydriftknob, wstarknob, explicit_option, flip_flop

    integer :: ierr, in_file

    if (readinit) return
    readinit = .true.

    if (proc0) then
       explicit_option = 'default'
       xdriftknob = 1.0
       ydriftknob = 1.0
       wstarknob = 1.0
       flip_flop = .false.

       in_file = input_unit_exist("time_advance_knobs", taexist)
       if (taexist) read (unit=in_file, nml=time_advance_knobs)

       ierr = error_unit()
       call get_option_value &
            (explicit_option, explicitopts, explicit_option_switch, &
            ierr, "explicit_option in time_advance_knobs")
    end if

    call broadcast (explicit_option_switch)
    call broadcast (xdriftknob)
    call broadcast (ydriftknob)
    call broadcast (wstarknob)
    call broadcast (flip_flop)

    if (none_implicit) flip_flop = .false.

  end subroutine read_parameters

! subroutine write_drifts
!   use dist_fn_arrays, only: wdriftx_g, wdrifty_g
!   use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
!   use dist_fn_arrays, only: wstar, wstarp
!   use dist_fn_arrays, only: wdriftpx_g, wdriftpy_g
!   use dist_fn_arrays, only: wdriftpx_phi, wdriftpy_phi
!   use zgrid, only: nzgrid
!   use stella_layouts, only: vmu_lo

!   use file_utils, only: run_name
!
!   implicit none

!   integer ia, iz, ivmu
!   character(len=512) :: filename

!   ia=1

!   filename=trim(run_name)//".drifts"
!   open(3345,file=trim(filename),status='unknown')
!   do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!     do iz= -nzgrid, nzgrid
!       write(3345,'(10e25.8)') wstar(ia,iz,ivmu),wstarp(ia,iz,ivmu), &
!                              wdriftx_g(ia,iz,ivmu), wdriftpx_g(ia,iz,ivmu), &
!                              wdrifty_g(ia,iz,ivmu), wdriftpy_g(ia,iz,ivmu), &
!                              wdriftx_phi(ia,iz,ivmu), wdriftpx_phi(ia,iz,ivmu), &
!                              wdrifty_phi(ia,iz,ivmu), wdriftpy_phi(ia,iz,ivmu)
!     enddo
!   enddo
!   close (3345)

! end subroutine write_drifts

  subroutine init_wdrift

    use dist_fn_arrays, only: wdriftx_g, wdrifty_g
    use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use species, only: spec
    use zgrid, only: nzgrid
    use kt_grids, only: nalpha
    use stella_geometry, only: cvdrift, gbdrift
    use stella_geometry, only: cvdrift0, gbdrift0
    use stella_geometry, only: gds23, gds24
    use stella_geometry, only: geo_surf, q_as_x
    use stella_geometry, only: dxdXcoord, drhodpsi, dydalpha
    use vpamu_grids, only: vpa, vperp2
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use neoclassical_terms, only: include_neoclassical_terms
    use neoclassical_terms, only: dphineo_dzed, dphineo_drho, dphineo_dalpha
    use neoclassical_terms, only: dfneo_dvpa, dfneo_dzed, dfneo_dalpha

    implicit none

    integer :: ivmu, iv, imu, is, ia
    real :: fac
    real, dimension (:,:), allocatable :: wcvdrifty, wgbdrifty
    real, dimension (:,:), allocatable :: wcvdriftx, wgbdriftx

    if (wdriftinit) return
    wdriftinit = .true.

    if (.not.allocated(wdriftx_phi)) then
      allocate (wdriftx_phi(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      wdriftx_phi = 0.0
    endif
    if (.not.allocated(wdrifty_phi)) then
      allocate (wdrifty_phi(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      wdrifty_phi = 0.0
    endif
    if (.not.allocated(wdriftx_g)) then
      allocate (wdriftx_g(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      wdriftx_g = 0.0
    endif
    if (.not.allocated(wdrifty_g)) then
      allocate (wdrifty_g(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      wdrifty_g = 0.0
    endif

    allocate (wcvdrifty(nalpha,-nzgrid:nzgrid))
    allocate (wgbdrifty(nalpha,-nzgrid:nzgrid))
    allocate (wcvdriftx(nalpha,-nzgrid:nzgrid))
    allocate (wgbdriftx(nalpha,-nzgrid:nzgrid))

    ia = 1
    ! FLAG -- need to deal with shat=0 case.  ideally move away from q as x-coordinate
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)

       fac = -ydriftknob*0.5*code_dt*spec(is)%tz_psi0
       ! this is the curvature drift piece of wdrifty with missing factor of vpa
       ! vpa factor is missing to avoid singularity when including
       ! non-Maxwellian corrections to equilibrium
       wcvdrifty = fac*cvdrift*vpa(iv)
       ! this is the grad-B drift piece of wdrifty
       wgbdrifty = fac*gbdrift*0.5*vperp2(:,:,imu)
       wdrifty_g(:,:,ivmu) = wcvdrifty*vpa(iv) + wgbdrifty
       ! if including neoclassical correction to equilibrium Maxwellian,
       ! then add in v_E^{nc} . grad y dg/dy coefficient here
       if (include_neoclassical_terms) then
          wdrifty_g(:,:,ivmu) = wdrifty_g(:,:,ivmu)+code_dt*0.5*(gds23*dphineo_dzed &
               + drhodpsi*dydalpha*dphineo_drho)
       end if

       wdrifty_phi(:,:,ivmu) = spec(is)%zt*(wgbdrifty + wcvdrifty*vpa(iv)) &
            * maxwell_vpa(iv,is)*maxwell_mu(:,:,imu,is)*maxwell_fac(is)
       ! if including neoclassical corrections to equilibrium,
       ! add in -(Ze/m) * v_curv/vpa . grad y d<phi>/dy * dF^{nc}/dvpa term
       ! and v_E . grad z dF^{nc}/dz (here get the dphi/dy part of v_E)
       if (include_neoclassical_terms) then
          wdrifty_phi(:,:,ivmu) = wdrifty_phi(:,:,ivmu) &
               - 0.5*spec(is)%zt*dfneo_dvpa(:,:,ivmu)*wcvdrifty &
               - code_dt*0.5*dfneo_dzed(:,:,ivmu)*gds23
       end if

       if(q_as_x) then
         fac = -xdriftknob*0.5*code_dt*spec(is)%tz_psi0
       else
         fac = -xdriftknob*0.5*code_dt*spec(is)%tz_psi0/geo_surf%shat
       endif
       ! this is the curvature drift piece of wdriftx with missing factor of vpa
       ! vpa factor is missing to avoid singularity when including
       ! non-Maxwellian corrections to equilibrium
       wcvdriftx = fac*cvdrift0*vpa(iv)
       ! this is the grad-B drift piece of wdriftx
       wgbdriftx = fac*gbdrift0*0.5*vperp2(:,:,imu)
       wdriftx_g(:,:,ivmu) = wcvdriftx*vpa(iv) + wgbdriftx
       ! if including neoclassical correction to equilibrium Maxwellian,
       ! then add in v_E^{nc} . grad x dg/dx coefficient here
       if (include_neoclassical_terms) then
          wdriftx_g(:,:,ivmu) = wdriftx_g(:,:,ivmu)+code_dt*0.5*(gds24*dphineo_dzed &
               - dxdXcoord*dphineo_dalpha)
       end if
       wdriftx_phi(:,:,ivmu) = spec(is)%zt*(wgbdriftx + wcvdriftx*vpa(iv)) &
            * maxwell_vpa(iv,is)*maxwell_mu(:,:,imu,is)*maxwell_fac(is)
       ! if including neoclassical corrections to equilibrium,
       ! add in (Ze/m) * v_curv/vpa . grad x d<phi>/dx * dF^{nc}/dvpa term
       ! and v_E . grad z dF^{nc}/dz (here get the dphi/dx part of v_E)
       ! and v_E . grad alpha dF^{nc}/dalpha (dphi/dx part of v_E)
       if (include_neoclassical_terms) then
          wdriftx_phi(:,:,ivmu) = wdriftx_phi(:,:,ivmu) &
               - 0.5*spec(is)%zt*dfneo_dvpa(:,:,ivmu)*wcvdriftx &
               + code_dt*0.5*(dfneo_dalpha(:,:,ivmu)*dxdXcoord-dfneo_dzed(:,:,ivmu)*gds24)
       end if

    end do

    deallocate (wcvdriftx, wgbdriftx, wcvdrifty, wgbdrifty)

  end subroutine init_wdrift

  subroutine init_wstar

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use species, only: spec
    use zgrid, only: nzgrid
    use kt_grids, only: nalpha
    use stella_geometry, only: dydalpha, drhodpsi
    use vpamu_grids, only: vperp2, vpa
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use dist_fn_arrays, only: wstar
    use neoclassical_terms, only: include_neoclassical_terms
    use neoclassical_terms, only: dfneo_drho

    implicit none

    integer :: is, imu, iv, ivmu
    real, dimension (:,:), allocatable :: energy

    if (wstarinit) return
    wstarinit = .true.

    if (.not.allocated(wstar)) &
         allocate (wstar(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc)) ; wstar = 0.0

    allocate (energy(nalpha,-nzgrid:nzgrid))

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       iv = iv_idx(vmu_lo,ivmu)
       energy = (vpa(iv)**2 + vperp2(:,:,imu))*(spec(is)%temp_psi0/spec(is)%temp)
       if (include_neoclassical_terms) then
          wstar(:,:,ivmu) = dydalpha*drhodpsi*wstarknob*0.5*code_dt &
               * (maxwell_vpa(iv,is)*maxwell_mu(:,:,imu,is)*maxwell_fac(is) &
               * (spec(is)%fprim+spec(is)%tprim*(energy-1.5)) &
               - dfneo_drho(:,:,ivmu))
       else
          wstar(:,:,ivmu) = dydalpha*drhodpsi*wstarknob*0.5*code_dt &
               * maxwell_vpa(iv,is)*maxwell_mu(:,:,imu,is)*maxwell_fac(is) &
               * (spec(is)%fprim+spec(is)%tprim*(energy-1.5))
       end if
    end do

    deallocate (energy)

  end subroutine init_wstar

  subroutine init_drifts_implicit

    use constants, only: zi
    use mp, only: sum_allreduce, mp_abort
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use gyro_averages, only: aj0x
    use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
    use dist_fn_arrays, only: wdriftx_g, wdrifty_g
    use dist_fn_arrays, only: wstar
    use fields_arrays, only: gamtot
    use fields, only: efac
    use run_parameters, only: fphi, fapar, time_upwind
    use species, only: spec, has_electron_species
    use stella_geometry, only: dl_over_b
    use zgrid, only: nzgrid
    use vpamu_grids, only: integrate_species
    use species, only: spec
    use kt_grids, only: naky, nakx, aky, akx, zonal_mode
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg

    implicit none

      integer :: ivmu, iz, ikx, is, ia, iv, imu
      complex :: tmps
      complex, dimension (:,:,:), allocatable :: g0
      complex, dimension (:,:), allocatable :: wd_g, wd_phi, wstr, tmp

      if (driftimpinit) return
      driftimpinit = .true.

      ia = 1

      allocate (wd_g(naky,nakx))
      allocate (wd_phi(naky,nakx))
      allocate (wstr(naky,nakx))
      allocate (tmp(naky,nakx))

      if (.not.allocated(gamtot_drifts)) &
           allocate (gamtot_drifts(naky,nakx,-nzgrid:nzgrid))
      gamtot_drifts = 0.
      if (.not.allocated(gamtot3_drifts)) &
           allocate (gamtot3_drifts(nakx,-nzgrid:nzgrid))
      gamtot3_drifts = 0.
!     if (.not.allocated(apar_denom_drifts)) &
!          allocate (apar_denom_wstar(naky,nakx,-nzgrid:nzgrid))
!     apar_denom_wstar = 0.

      if (fphi > epsilon(0.0)) then
        allocate (g0(naky,nakx,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        do iz = -nzgrid, nzgrid
          do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv  = iv_idx(vmu_lo,ivmu)
            imu = imu_idx(vmu_lo,ivmu)
            is  = is_idx(vmu_lo,ivmu)

            !there terms already contain a factor of code_dt as well as
            !a negative sign to account for RHS
            wd_g   = -zi*(spread(akx,1,naky)*wdriftx_g(ia,iz,ivmu) &
                        + spread(aky,2,nakx)*wdrifty_g(ia,iz,ivmu))

            wd_phi = -zi*(spread(akx,1,naky)*wdriftx_phi(ia,iz,ivmu) &
                        + spread(aky,2,nakx)*wdrifty_phi(ia,iz,ivmu))

            wstr   = -zi*spread(aky,2,nakx)*wstar(ia,iz,ivmu)

            g0(:,:,ivmu) = 0.5*(1.0+time_upwind)*aj0x(:,:,iz,ivmu)**2 &
                           *(wd_phi + wstr)/(1.0 + 0.5*(1.0+time_upwind)*wd_g)
          enddo
          call integrate_species (g0, iz, spec%z*spec%dens_psi0, gamtot_drifts(:,:,iz))
        enddo

        gamtot_drifts = gamtot_drifts + gamtot

        deallocate (g0)

        if (.not.has_electron_species(spec)) then
          ! no need to do anything extra for ky /= 0 because
          ! already accounted for in gamtot_h
          if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            if (zonal_mode(1)) then
              do ikx = 1, nakx
                ! do not need kx=ky=0 mode
                if (abs(akx(ikx)) < epsilon(0.)) cycle
                tmps = 1.0/efac - sum(dl_over_b(ia,:)/gamtot_drifts(1,ikx,:))
                gamtot3_drifts(ikx,:) = 1./(gamtot_drifts(1,ikx,:)*tmps)
              end do
              if (akx(1) < epsilon(0.)) gamtot3_drifts(1,:) = 0.0
            end if
          end if
        end if
      end if

      deallocate (wd_g, wd_phi, wstr, tmp)


      ! FLAG -- NEED TO SORT OUT FINITE FAPAR FOR GSTAR
       if (fapar > epsilon(0.)) then
          write (*,*) 'APAR NOT SETUP FOR GSTAR YET. aborting'
          call mp_abort ('APAR NOT SETUP FOR GSTAR YET. aborting')
       end if

  !        allocate (g0(-nvgrid:nvgrid,nmu))
  !        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
  !           iky = iky_idx(kxkyz_lo,ikxkyz)
  !           ikx = ikx_idx(kxkyz_lo,ikxkyz)
  !           ig = iz_idx(kxkyz_lo,ikxkyz)
  !           is = is_idx(kxkyz_lo,ikxkyz)
  !           g0 = spread(vpa**2,2,nmu)*spread(aj0v(:,ikxkyz)**2,1,nvpa)*anon(ig,:,:)
  !           wgt = 2.0*beta*spec(is)%z*spec(is)%z*spec(is)%dens/spec(is)%mass
  !           call integrate_vmu (g0, ig, tmp)
  !           apar_denom(iky,ikx,ig) = apar_denom(iky,ikx,ig) + tmp*wgt
  !        end do
  !        call sum_allreduce (apar_denom)
  !        apar_denom = apar_denom + kperp2

  !        deallocate (g0)
  !     end if

  end subroutine init_drifts_implicit

  subroutine init_parallel_nonlinearity

    use physics_parameters, only: rhostar
    use species, only: spec, nspec
    use zgrid, only: nztot, nzgrid
    use stella_geometry, only: geo_surf, drhodpsi, q_as_x
    use stella_geometry, only: gradpar, dbdzed, bmag
    use stella_geometry, only: cvdrift, cvdrift0
    use stella_geometry, only: dIdrho, dgradpardrho, dBdrho, d2Bdrdth
    use stella_geometry, only: dcvdriftdrho, dcvdrift0drho
    use physics_flags, only: radial_variation

    implicit none

    if (.not. allocated(par_nl_fac)) allocate (par_nl_fac(-nzgrid:nzgrid,nspec))
    ! this is the factor multiplying -dphi/dz * dg/dvpa in the parallel nonlinearity
    par_nl_fac = 0.5*rhostar*spread(spec%stm_psi0*spec%zt_psi0,1,nztot)*spread(gradpar,2,nspec)

    if (.not. allocated(par_nl_curv)) allocate (par_nl_curv(-nzgrid:nzgrid,nspec))
    ! ydriftknob is here because this term comes from bhat x curvature . grad B
    par_nl_curv = -ydriftknob*rhostar*geo_surf%rgeo*geo_surf%betaprim*drhodpsi &
         *spread(dbdzed(1,:)*gradpar/bmag(1,:),2,nspec)/spread(spec%zt_psi0,1,nztot)

    if (.not. allocated(par_nl_drifty)) allocate (par_nl_drifty(-nzgrid:nzgrid))
    par_nl_drifty = 0.25*rhostar*cvdrift(1,:)
    if (.not. allocated(par_nl_driftx)) allocate (par_nl_driftx(-nzgrid:nzgrid))
    if(q_as_x) then
      par_nl_driftx = 0.25*rhostar*cvdrift0(1,:)
    else
      par_nl_driftx = 0.25*rhostar*cvdrift0(1,:)/geo_surf%shat
    endif

    if (radial_variation) then
      if (.not. allocated(d_par_nl_fac_dr)) allocate (d_par_nl_fac_dr(-nzgrid:nzgrid,nspec))
      ! this is the factor multiplying -dphi/dz * dg/dvpa in the parallel nonlinearity
      d_par_nl_fac_dr = 0.5*rhostar*spread(spec%stm_psi0*spec%zt_psi0,1,nztot)*spread(dgradpardrho,2,nspec)

      if (.not. allocated(d_par_nl_curv_dr)) allocate (d_par_nl_curv_dr(-nzgrid:nzgrid,nspec))
      ! ydriftknob is here because this term comes from bhat x curvature . grad B
      ! handle terms with no zeroes
      d_par_nl_curv_dr = par_nl_curv*(dIdrho/geo_surf%rgeo - drhodpsi*geo_surf%d2psidr2 &
                         - spread(dBdrho/bmag(1,:) + dgradpardrho/gradpar,2,nspec))
      ! handle terms with possible zeroes
      d_par_nl_curv_dr = d_par_nl_curv_dr &
                       - ((ydriftknob*rhostar*geo_surf%rgeo*drhodpsi*spread(gradpar/bmag(1,:),2,nspec)) &
                          /spread(spec%zt_psi0,1,nztot)) &
                         *(  geo_surf%betadbprim*spread(dbdzed(1,:),2,nspec) &
                           + geo_surf%betaprim*spread(d2Bdrdth,2,nspec))

      if (.not. allocated(d_par_nl_drifty_dr)) allocate (d_par_nl_drifty_dr(-nzgrid:nzgrid))
      d_par_nl_drifty_dr = 0.25*rhostar*dcvdriftdrho(1,:)
      if (.not. allocated(d_par_nl_drifty_dr)) allocate (d_par_nl_driftx_dr(-nzgrid:nzgrid))
      if(q_as_x) then
        d_par_nl_driftx_dr = 0.25*rhostar*dcvdrift0drho(1,:)
      else
        d_par_nl_driftx_dr = 0.25*rhostar*dcvdrift0drho(1,:)/geo_surf%shat
      endif
    endif

    parnlinit = .true.

  end subroutine init_parallel_nonlinearity

  subroutine init_radial_variation
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use species, only: spec, pfac
    use zgrid, only: nzgrid
    use kt_grids, only: nalpha
    use stella_geometry, only: drhodpsi, dydalpha, gfac
    use stella_geometry, only: dBdrho, geo_surf, q_as_x
    use stella_geometry, only: dcvdriftdrho, dcvdrift0drho
    use stella_geometry, only: dgbdriftdrho, dgbdrift0drho
    use vpamu_grids, only: vperp2, vpa, mu
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use dist_fn_arrays, only: wstarp
    use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
    use dist_fn_arrays, only: wdriftpx_g, wdriftpy_g
    use dist_fn_arrays, only: wdriftpx_phi, wdriftpy_phi!, adiabatic_phi
!   use neoclassical_terms, only: include_neoclassical_terms

    implicit none

    integer :: is, imu, iv, ivmu
    real :: fac
    real, dimension (:,:), allocatable :: energy

    real, dimension (:,:), allocatable :: wcvdrifty, wgbdrifty
    real, dimension (:,:), allocatable :: wcvdriftx, wgbdriftx

!wstar

    if (radialinit) return
    radialinit = .true.

    allocate (wcvdrifty(nalpha,-nzgrid:nzgrid))
    allocate (wgbdrifty(nalpha,-nzgrid:nzgrid))
    allocate (wcvdriftx(nalpha,-nzgrid:nzgrid))
    allocate (wgbdriftx(nalpha,-nzgrid:nzgrid))

    if (.not.allocated(wstarp)) &
      allocate (wstarp(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc)) ; wstarp = 0.0
    if (.not.allocated(wdriftpx_phi)) &
       allocate (wdriftpx_phi(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
!   if (.not.allocated(adiabatic_phi)) &
!      allocate (adiabatic_phi(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    if (.not.allocated(wdriftpy_phi)) &
      allocate (wdriftpy_phi(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    if (.not.allocated(wdriftpx_g)) &
      allocate (wdriftpx_g(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    if (.not.allocated(wdriftpy_g)) &
      allocate (wdriftpy_g(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate (energy(nalpha,-nzgrid:nzgrid))

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       iv = iv_idx(vmu_lo,ivmu)
       energy = (vpa(iv)**2 + vperp2(:,:,imu))*(spec(is)%temp_psi0/spec(is)%temp)
       !FLAG DSO - THIS NEEDS TO BE ADDED SOMEDAY!
       !if (include_neoclassical_terms) then
       !   wstarp(:,:,ivmu) = dydalpha*drhodpsi*wstarknob*0.5*code_dt &
       !        * (maxwell_vpa(iv)*maxwell_mu(:,:,imu) &
       !        * (spec(is)%fprim+spec(is)%tprim*(energy-1.5)) &
       !        - dfneo_drho(:,:,ivmu))
       !else
       !recall that fprim = -dn/dr and trpim = -dt/dr

       wstarp(:,:,ivmu) = -wstarknob*0.5*code_dt &
          * dydalpha*drhodpsi*maxwell_vpa(iv,is)*maxwell_mu(:,:,imu,is)*maxwell_fac(is) &
          * (pfac*(spec(is)%d2ndr2-(spec(is)%fprim)**2-(spec(is)%tprim)**2*energy) &
           + pfac*(spec(is)%d2Tdr2-(spec(is)%tprim)**2)*(energy-1.5) &
           - gfac*2*spec(is)%tprim*mu(imu)*spread(dBdrho,1,nalpha) &
           + (spec(is)%fprim+spec(is)%tprim*(energy-1.5)) &
           * (  pfac*(spec(is)%fprim+spec(is)%tprim*(energy-1.5))  &
              + gfac*2*mu(imu)*spread(dBdrho,1,nalpha) &
              + gfac*drhodpsi*geo_surf%d2psidr2))

       !end if

       !wdrift
       fac = -ydriftknob*0.5*code_dt*spec(is)%tz_psi0
       ! this is the curvature drift piece of wdrifty with missing factor of vpa
       ! vpa factor is missing to avoid singularity when including
       ! non-Maxwellian corrections to equilibrium
       wcvdrifty = gfac*fac*dcvdriftdrho*vpa(iv)
       ! this is the grad-B drift piece of wdrifty
       wgbdrifty = gfac*fac*dgbdriftdrho*0.5*vperp2(:,:,imu)
       wdriftpy_g(:,:,ivmu) = wcvdrifty*vpa(iv) + wgbdrifty

       wdriftpy_phi(:,:,ivmu) = spec(is)%zt*(wgbdrifty + wcvdrifty*vpa(iv)) &
            * maxwell_vpa(iv,is)*maxwell_mu(:,:,imu,is)*maxwell_fac(is) &
            - wdrifty_phi(:,:,ivmu)*(pfac*(spec(is)%fprim + spec(is)%tprim*(energy-2.5)) &
              + gfac*2.*mu(imu)*spread(dBdrho,1,nalpha))


       if(q_as_x) then
         fac = -xdriftknob*0.5*code_dt*spec(is)%tz_psi0
       else
         fac = -xdriftknob*0.5*code_dt*spec(is)%tz_psi0/geo_surf%shat
       endif
       ! this is the curvature drift piece of wdriftx with missing factor of vpa
       ! vpa factor is missing to avoid singularity when including
       ! non-Maxwellian corrections to equilibrium
       wcvdriftx = gfac*fac*dcvdrift0drho*vpa(iv)
       ! this is the grad-B drift piece of wdriftx
       wgbdriftx = gfac*fac*dgbdrift0drho*0.5*vperp2(:,:,imu)
       wdriftpx_g(:,:,ivmu) = wgbdriftx + wcvdriftx*vpa(iv)

       wdriftpx_phi(:,:,ivmu) = spec(is)%zt*(wgbdriftx + wcvdriftx*vpa(iv)) &
            * maxwell_vpa(iv,is)*maxwell_mu(:,:,imu,is)*maxwell_fac(is) &
            - wdriftx_phi(:,:,ivmu)*(pfac*(spec(is)%fprim + spec(is)%tprim*(energy-2.5)) &
              + gfac*2.*mu(imu)*spread(dBdrho,1,nalpha))


!      !the next piece is everything under the x derivative, as this needs to be
!      !transformed separately
!      wdriftpx_phi(:,:,ivmu) = spec(is)%zt*(wgbdriftx + wcvdriftx*vpa(iv))  &
!           * maxwell_vpa(iv,is)*maxwell_mu(:,:,imu,is)*maxwell_fac(is)

!      !this is variation in the Maxwellian part of the adiabatic response of phi,
!      !which needs to be transformed separately before differentiation wrt x
!      !the gyroaveraging and quasineutrality is already done in fields
!      adiabatic_phi(:,:,ivmu) = -(pfac*(spec(is)%fprim+spec(is)%tprim*(energy-2.5)) &
!                                 +gfac*2.*mu(imu)*spread(dBdrho,1,nalpha))

    end do

    deallocate (energy, wcvdriftx, wgbdriftx, wcvdrifty, wgbdrifty)

  end subroutine init_radial_variation


  subroutine allocate_arrays

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx
    use dist_fn_arrays, only: g0, g1, g2, g3

    implicit none

    if (.not.allocated(g0)) &
         allocate (g0(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    g0 = 0.
    if (.not.allocated(g1)) &
         allocate (g1(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    g1 = 0.
    if (.not.allocated(g2)) &
         allocate (g2(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    g2 = 0.
    if (.not.allocated(g3) .and. explicit_option_switch==explicit_option_rk4) then
       allocate (g3(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       g3 = 0.
    else
       allocate (g3(1,1,1,1,1))
    end if

  end subroutine allocate_arrays

  subroutine init_cfl

    use mp, only: proc0, nproc, max_allreduce, min_allreduce
    use mp, only: scope, allprocs, subprocs
    use dist_fn_arrays, only: wdriftx_g, wdrifty_g
    use stella_time, only: cfl_dt, code_dt, write_dt
    use run_parameters, only: cfl_cushion
    use physics_flags, only: radial_variation, prp_shear_enabled
    use zgrid, only: delzed
    use vpamu_grids, only: dvpa
    use kt_grids, only: akx, aky, nx, rho
    use run_parameters, only: stream_implicit, mirror_implicit, drifts_implicit
    use parallel_streaming, only: stream
    use parallel_streaming, only: stream_rad_var1, stream_rad_var2
    use mirror_terms, only: mirror
    use flow_shear, only: prl_shear, shift_times
    use file_utils, only: runtype_option_switch, runtype_multibox
    use dissipation, only: include_collisions, collisions_implicit
    use dissipation, only: vpa_operator, mu_operator
    use dissipation, only: cfl_dt_vpadiff, cfl_dt_mudiff

    implicit none

    real :: cfl_dt_mirror, cfl_dt_stream, cfl_dt_shear
    real :: cfl_dt_wdriftx, cfl_dt_wdrifty
    real :: zero
    real :: wdriftx_max, wdrifty_max

    ! avoid divide by zero in cfl_dt terms below
    zero = 100.*epsilon(0.)

    ! FLAG -- assuming equal spacing in zed!

    ! Bob: Not clear why cfl_dt is set this way. If no terms are explicit, then
    ! cfl_dt is arbitrarily large (there's no CFL constraint).
    ! if (cfl_dt.lt.0) cfl_dt = code_dt/cfl_cushion
    ! Replace with very large number
    if (cfl_dt.lt.0) cfl_dt = 99999E99

    if(.not.drifts_implicit) then
      ! get the local max value of wdriftx on each processor
      wdriftx_max = maxval(abs(wdriftx_g))
      ! compare these max values across processors to get global max
      if (nproc > 1) then
        call max_allreduce (wdriftx_max)
      endif
      ! NB: wdriftx_g has code_dt built-in, which accounts for code_dt factor here
      cfl_dt_wdriftx = abs(code_dt)/max(maxval(abs(akx))*wdriftx_max,zero)
      cfl_dt = cfl_dt_wdriftx
    endif

    cfl_dt_shear = abs(code_dt)/max(maxval(abs(aky))*maxval(abs(prl_shear)),zero)
    cfl_dt = min(cfl_dt,cfl_dt_shear)

    if(prp_shear_enabled) then
      cfl_dt_shear = minval(shift_times)
      cfl_dt = min(cfl_dt,cfl_dt_shear)
    endif


    if (.not.stream_implicit) then
       ! NB: stream has code_dt built-in, which accounts for code_dt factor here
       cfl_dt_stream = abs(code_dt)*delzed(0)/max(maxval(abs(stream)),zero)
       cfl_dt = min(cfl_dt,cfl_dt_stream)
    end if

    if (.not.mirror_implicit) then
       ! NB: mirror has code_dt built-in, which accounts for code_dt factor here
       cfl_dt_mirror = abs(code_dt)*dvpa/max(maxval(abs(mirror)),zero)
       cfl_dt = min(cfl_dt,cfl_dt_mirror)
    end if

    if (radial_variation) then
      !while other quantities should go here, parallel streaming with electrons
      !is what will limit us
      cfl_dt_stream = abs(code_dt)*delzed(0)/max(maxval(abs(stream_rad_var1)),zero)
      cfl_dt_stream = cfl_dt_stream/abs(rho(nx)+zero)
      cfl_dt = min(cfl_dt,cfl_dt_stream)

      cfl_dt_stream = abs(code_dt)*delzed(0)/max(maxval(abs(stream_rad_var2)),zero)
      cfl_dt_stream = cfl_dt_stream/abs(rho(nx)+zero)
      cfl_dt = min(cfl_dt,cfl_dt_stream)

    end if

    if (include_collisions.and..not.collisions_implicit) then
      if (vpa_operator) cfl_dt = min(cfl_dt,cfl_dt_vpadiff)
      if (mu_operator)  cfl_dt = min(cfl_dt,cfl_dt_mudiff)
    endif

    if(.not.drifts_implicit) then
      ! get the local max value of wdrifty on each processor
      wdrifty_max = maxval(abs(wdrifty_g))
      ! compare these max values across processors to get global max
      if (nproc > 1) then
        call max_allreduce (wdrifty_max)
      endif
      ! NB: wdrifty_g has code_dt built-in, which accounts for code_dt factor here
      cfl_dt_wdrifty = abs(code_dt)/max(maxval(abs(aky))*wdrifty_max,zero)
      cfl_dt = min(cfl_dt,cfl_dt_wdrifty)
    endif

    if(runtype_option_switch == runtype_multibox) call scope(allprocs)
    call min_allreduce (cfl_dt)
    if(runtype_option_switch == runtype_multibox) call scope(subprocs)

    if (proc0) then
       write (*,'(A)') "############################################################"
       write (*,'(A)') "                        CFL CONDITION"
       write (*,'(A)') "############################################################"
       write (*,'(A16)') 'LINEAR CFL_DT: '
       if (.not.drifts_implicit) write (*,'(A12,ES12.4)') '   wdriftx: ', cfl_dt_wdriftx
       if (.not.drifts_implicit) write (*,'(A12,ES12.4)') '   wdrifty: ', cfl_dt_wdrifty
       if (.not.stream_implicit) write (*,'(A12,ES12.4)') '   stream: ', cfl_dt_stream
       if (.not.mirror_implicit) write (*,'(A12,ES12.4)') '   mirror: ', cfl_dt_mirror
       write (*,*)
    end if

    if (abs(code_dt) > cfl_dt*cfl_cushion) then
       if (proc0) then
          write (*,*) 'CHANGING TIME STEP:'
          write (*,'(A16, ES10.2E2)') "   code_dt:"//REPEAT(' ',50),code_dt
          write (*,'(A16, ES10.2E2)') "   cfl_dt:"//REPEAT(' ',50),cfl_dt
          write (*,'(A16, ES10.2E2)') "   cfl_cushion:"//REPEAT(' ',50),cfl_cushion
          write (*,'(A65)') '     ==> User-specified delt is larger than cfl_dt*cfl_cushion.'//REPEAT(' ',50)
          write (*,'(A49,ES12.4)') '     ==> Changing code_dt to cfl_dt*cfl_cushion ='//REPEAT(' ',50), cfl_dt*cfl_cushion
       end if
       code_dt = sign(1.0,code_dt)*cfl_dt*cfl_cushion
       call reset_dt
    else if (proc0) then
       call write_dt
       write (*,*)
    end if

  end subroutine init_cfl

  subroutine reset_dt

    use parallel_streaming, only: parallel_streaming_initialized
    use parallel_streaming, only: init_parallel_streaming
    use dissipation, only: init_collisions, collisions_initialized, include_collisions
    use run_parameters, only: stream_implicit, driftkinetic_implicit, drifts_implicit
    use response_matrix, only: response_matrix_initialized
    use response_matrix, only: init_response_matrix
    use mirror_terms, only: mirror_initialized
    use mirror_terms, only: init_mirror
    use flow_shear, only: flow_shear_initialized
    use flow_shear, only: init_flow_shear
    use physics_flags, only: radial_variation

    implicit none

    ! need to recompute mirror and streaming terms
    ! to account for updated code_dt
    wdriftinit = .false.
    wstarinit = .false.
    radialinit = .false.
    driftimpinit = .false.
    flow_shear_initialized = .false.
    mirror_initialized = .false.
    parallel_streaming_initialized = .false.

    call init_wstar
    call init_wdrift
    call init_mirror
    call init_parallel_streaming
    call init_flow_shear
    if (radial_variation) call init_radial_variation
    if (drifts_implicit) call init_drifts_implicit
    if (include_collisions) then
      collisions_initialized = .false.
      call init_collisions
    endif
    ! do not try to re-init response matrix
    ! before it has been initialized the first time
    if ((stream_implicit.or.driftkinetic_implicit) .and. response_matrix_initialized) then
       response_matrix_initialized = .false.
       call init_response_matrix
    end if

  end subroutine reset_dt

  subroutine advance_stella (istep)

    use dist_fn_arrays, only: golder, gold, gnew
    use fields_arrays, only: phi, apar
    use fields_arrays, only: phi_old
    use fields, only: advance_fields, fields_updated
    use run_parameters, only: none_implicit, use_leapfrog_splitting, nisl_nonlinear
    use run_parameters, only:  leapfrog_drifts, leapfrog_nonlinear
    use run_parameters, only: stream_implicit, mirror_implicit, drifts_implicit
    use physics_flags, only: nonlinear
    use physics_flags, only: include_parallel_nonlinearity
    use physics_flags, only: include_mirror
    use physics_flags, only: include_drifts
    use physics_flags, only: include_parallel_streaming
    use physics_flags, only: full_flux_surface, radial_variation
    use physics_parameters, only: g_exb
    use multibox, only: RK_step
    use dissipation, only: include_krook_operator, update_delay_krook
    use dissipation, only: remove_zero_projection, project_out_zero
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx
    use stella_layouts, only: vmu_lo
    use dissipation, only: include_collisions, collisions_implicit
    use dissipation, only: include_krook_operator
    use file_utils, only: runtype_option_switch, runtype_multibox
    use multibox, only: include_multibox_krook

    implicit none

    integer, intent (in) :: istep
    complex, allocatable, dimension (:,:,:,:) :: g1
    logical  :: restart_time_step, time_advance_successful, reverse_implicit_order

    if(.not.RK_step) then
      if (debug) write (*,*) 'time_advance::multibox'
      call mb_communicate(gnew)
    endif

    ! save value of phi
    ! for use in diagnostics (to obtain frequency)
    phi_old = phi
    ! We can use the leapfrog step provided:
    ! 1) istep > 1 (need gold and golder)
    ! 2) The timestep hasn't just been reset (because then golder, gold, gnew
    ! aren't equally spaced in time)
    ! 3) The simulation hasn't just been restarted - currently, save_for_restart
    ! only saves gold, so we don't have golder.
    leapfrog_this_timestep = .false.

    ! Flag which is set to true once we've taken a step without needing to
    ! reset dt.
    time_advance_successful = .false.
    ! if CFL condition is violated by nonlinear term
    ! then must modify time step size and restart time step
    ! assume false and test.
    ! We don't want to do any operations after dt is reset (we discard any
    ! updates to g and the timestep again), so we'll be frequently checking
    ! this flag.
    restart_time_step = .false.
    !! Original
    !if (use_leapfrog_splitting .and. istep > 1) then
    !!! Hack - in a sim we're doing, we restart from istep=10000. Want the first step
    ! (istep=10001) to be RK
    if (use_leapfrog_splitting .and. istep > 1) then
      leapfrog_this_timestep = .true.
    end if

    ! Work out whether there's any terms which are to be evolved by a
    ! Runge-Kutta scheme this timestep.
    ! TODO: Tidy up this logic, maybe put in its own subroutine.
    if ((nonlinear) .and. ( (.not. leapfrog_this_timestep))) then
      runge_kutta_terms_this_timestep = .true.
    else if (include_parallel_nonlinearity) then
      runge_kutta_terms_this_timestep = .true.
    else if ((g_exb**2).gt.epsilon(0.0)) then
      runge_kutta_terms_this_timestep = .true.
    else if (include_mirror.and..not.mirror_implicit) then
      runge_kutta_terms_this_timestep = .true.
    else if ( (include_drifts .and. .not.drifts_implicit) .and. (.not. leapfrog_drifts .or. .not. leapfrog_this_timestep)) then
      runge_kutta_terms_this_timestep = .true.
    else if (include_collisions.and..not.collisions_implicit) then
      runge_kutta_terms_this_timestep = .true.
    else if (include_parallel_streaming.and.(.not.stream_implicit)) then
      runge_kutta_terms_this_timestep = .true.
    else if (radial_variation) then
      runge_kutta_terms_this_timestep = .true.
    else if (include_krook_operator) then
      runge_kutta_terms_this_timestep = .true.
    else if (runtype_option_switch == runtype_multibox .and. include_multibox_krook) then
      runge_kutta_terms_this_timestep = .true.
    else
      runge_kutta_terms_this_timestep = .false.
    end if

    ! For reversing the order in which the implicit operators are applied
    ! (either ABC...  or ...CBA )
    ! The order is never reversed for Lie splitting (always ABC...), but
    ! is reversed for flip-flop or Leapfrog splitting.
    reverse_implicit_order = .false.

    if (leapfrog_this_timestep) then
      ! Try taking a leapfrog step, but if dt is reset, we can't take a Leapfrog
      ! step (need to do an "ordinary" single-step approach).
      ! Need to check if dt is reset after the explicit step, and after the
      ! Leapfrog step - either could contain the nonlinearity, so either could
      ! cause dt to be reset.

      ! We advance golder by a step; we need to get the fields at golder, so
      ! set the flag to ensure fields get re-calculated. NB we could skip the
      ! field solve because these fields have been calculated previously, but at the
      ! expense of memory (storing another set of fields) and code refactoring
      fields_updated = .false.
      if (runge_kutta_terms_this_timestep) call advance_explicit (golder, restart_time_step)

      if (.not. restart_time_step) then
        ! NB reverse_implicit_order is .false. at this point
        if (.not.none_implicit) call advance_implicit (phi, apar, reverse_implicit_order, golder)
        ! advance_leapfrog uses the updated golder (g^{n-1} advanced by single step
        ! operators) and gnew = g^{n}. It returns gnew = golder + 2*rhs(gnew)
        ! To get rhs(gnew), we need the fields at gnew; set the flag to ensure
        ! fields get re-calculated.
        ! NB we could skip the
        ! field solve because these fields have been calculated previously, but at the
        ! expense of memory (storing another set of fields) and code refactoring
        fields_updated = .false.
        call advance_leapfrog(gold, golder, restart_time_step)

        if (.not. restart_time_step) then
          reverse_implicit_order = .true.  ! Swap the order in which we apply the implicit operators
          if (.not.none_implicit) call advance_implicit (phi, apar, reverse_implicit_order, golder)
          if (runge_kutta_terms_this_timestep) call advance_explicit (golder, restart_time_step)
        end if
      end if

      if (restart_time_step) then
        ! Need to re-do the step with Lie splitting
        ! This flag tells us we need to treat whatever terms we were going to
        ! treat with the leapfrog approach with a non-leapfrog scheme.
        leapfrog_this_timestep = .false.
        runge_kutta_terms_this_timestep = .true.
        ! We're discarding changes to gnew and starting over, so fields will
        ! need to be updated.
        fields_updated = .false.
      else
        time_advance_successful = .true.
        ! Update gnew
        gnew = golder
      end if
    end if

    if ((nisl_nonlinear) .and. (istep .eq. 1) )then
      ! Store the value of g(t=0) in golder
      golder = gold
      call advance_initial_nisl_step()

    !!! Bob: Should be in a separate branch?
    else
      ! Perform the Lie or flip-flopping step until we've done it without the
      ! timestep changing.
      ! If we're using NISL for the ExB nonlinearity, we want to take many steps

      ! with RK3, such that in each RK3 step the CFL constraint is not violated,
      ! but the RK3 steps all add up to a step of code_dt
      do while (.not. time_advance_successful)
        ! If we've already attempted a time advance then we've updated gnew,
        ! so reset it now.
        gnew = gold
        ! Ensure fields are consistent with gnew.
        call advance_fields (gnew, phi, apar, dist='gbar')


        ! Reset the flag that checks if we've reset dt
        restart_time_step = .false. ! Becomes true if we reset dt

        ! reverse the order of operations every time step
        ! as part of alternating direction operator splitting
        ! this is needed to ensure 2nd order accuracy in time
        if (mod(istep,2)==1 .or. .not.flip_flop) then
          reverse_implicit_order = .false.
          ! advance the explicit parts of the GKE
          call advance_explicit (gnew, restart_time_step)

          ! enforce periodicity for zonal mode
          !    if (zonal_mode(1)) gnew(1,:,-nzgrid,:) = gnew(1,:,nzgrid,:)

          ! use operator splitting to separately evolve
          ! all terms treated implicitly
          if (.not. restart_time_step .and. .not. none_implicit) call advance_implicit (phi, apar, reverse_implicit_order, gnew)

        else
          reverse_implicit_order = .true.
          if (.not.none_implicit) call advance_implicit (phi, apar, reverse_implicit_order, gnew)
          call advance_explicit (gnew, restart_time_step)
        end if

        if (.not. restart_time_step) then
          time_advance_successful = .true.
        else
          ! We're discarding changes to gnew and starting over, so fields will
          ! need to be re-calculated
          fields_updated = .false.
        end if
      end do

      if(remove_zero_projection) then
        allocate (g1(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        g1 = gnew(1,:,:,:,:) - gold(1,:,:,:,:)
        call project_out_zero(g1)
        gnew(1,:,:,:,:) = gnew(1,:,:,:,:) - g1
        deallocate (g1)
      end if

      !update the delay parameters for the Krook operator
      if(include_krook_operator) call update_delay_krook(gnew)

      golder = gold
      gold = gnew

      ! Ensure fields are updated so that omega calculation is correct.
      call advance_fields (gnew, phi, apar, dist='gbar')
    end if

  end subroutine advance_stella

  subroutine advance_initial_nisl_step()

    ! TODO: Get rid of unneeded variables
    use dist_fn_arrays, only: gold, gnew
    use fields_arrays, only: phi, apar
    use fields_arrays, only: phi_old
    use fields, only: advance_fields, fields_updated
    use run_parameters, only: none_implicit, use_leapfrog_splitting, nisl_nonlinear
    use run_parameters, only:  leapfrog_drifts, exact_exb_nonlinear_solution_first_step
    use run_parameters, only: stream_implicit, mirror_implicit, drifts_implicit
    use physics_flags, only: nonlinear
    use physics_flags, only: include_parallel_nonlinearity
    use physics_flags, only: include_mirror
    use physics_flags, only: include_drifts
    use physics_flags, only: include_parallel_streaming
    use physics_flags, only: full_flux_surface, radial_variation
    use physics_parameters, only: g_exb
    use stella_time, only: code_dt
    use multibox, only: RK_step
    use dissipation, only: include_krook_operator, update_delay_krook
    use dissipation, only: remove_zero_projection, project_out_zero
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx
    use stella_layouts, only: vmu_lo
    use dissipation, only: include_collisions, collisions_implicit
    use dissipation, only: include_krook_operator
    use file_utils, only: runtype_option_switch, runtype_multibox
    use multibox, only: include_multibox_krook

    implicit none

    complex, allocatable, dimension (:,:,:,:) :: g1
    logical  :: restart_time_step, time_advance_successful, reverse_implicit_order
    real :: total_dt, nisl_step_dt, tolerance

    !!! Take a series of steps using an explicit calculation of the nonlinear term.
    !!! We don't know how many steps this will be, because the CFL constraint
    !!! depends upon how big vexb is.

    if (exact_exb_nonlinear_solution_first_step) then
      call advance_ExB_nonlinearity_exact (gold, gnew)
      gold = gnew
      fields_updated = .false.
    else

      total_dt = 0


      ! We may end up changing dt, but at the end we need to change delt back to the
      ! original value, so store original value
      nisl_step_dt = code_dt

      ! Advance until total_dt = nisl_step_dt
      ! Include a tolerance in case the series of smaller dts are
      ! ~ machine precision less than code_dt
      tolerance = 1E-9
      do while (total_dt .lt. (nisl_step_dt-tolerance))
        ! Take a single time step. As before, keep going until
        ! restart_time_step = .false.
        time_advance_successful = .false.
        do while (.not. time_advance_successful)
          ! We're definitely not flip-flopping here, so no need to
          ! reverse order of operations
          reverse_implicit_order = .false.
          ! advance the explicit parts of the GKE
          !
          call advance_explicit (gnew, restart_time_step)
          ! enforce periodicity for zonal mode
          !    if (zonal_mode(1)) gnew(1,:,-nzgrid,:) = gnew(1,:,nzgrid,:)

          ! use operator splitting to separately evolve
          ! all terms treated implicitly
          if (.not. restart_time_step .and. .not. none_implicit) call advance_implicit (phi, apar, reverse_implicit_order, gnew)
          if (.not. restart_time_step) then
            time_advance_successful = .true.
          else
            ! Because we're unsuccessful, need to revert gnew to gold
            gnew = gold
            ! We're discarding changes to gnew and starting over, so fields will
            ! need to be re-calculated
            fields_updated = .false.
          end if
        end do

        if(remove_zero_projection) then
          allocate (g1(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
          g1 = gnew(1,:,:,:,:) - gold(1,:,:,:,:)
          call project_out_zero(g1)
          gnew(1,:,:,:,:) = gnew(1,:,:,:,:) - g1
          deallocate (g1)
        end if

        !update the delay parameters for the Krook operator
        ! Bob: What does this do and does it need to be done?
        if(include_krook_operator) call update_delay_krook(gnew)

        ! Update gold. We don't need to update golder because it's not used here,
        ! and we don't want to update it because we want to use g(t=0) as golder
        ! in subsequent steps
        gold = gnew

        ! Ensure fields are updated so that omega calculation is correct.
        call advance_fields (gnew, phi, apar, dist='gbar')
        ! Update the total time advanced
        total_dt = total_dt + code_dt
      end do

      ! We've now taken a series of steps and have gold, corresponding to
      ! g(t=nisl_step_dt). If the timestep has been changed, change it back
      ! to nisl_step_dt
      if (code_dt .lt. nisl_step_dt) then
        code_dt = nisl_step_dt
        call reset_dt
      end if
    end if
  end subroutine advance_initial_nisl_step

!  subroutine advance_explicit (phi, apar, g)
  subroutine advance_explicit (g, restart_time_step)

    use mp, only: proc0
    use job_manage, only: time_message
    use zgrid, only: nzgrid
    use extended_zgrid, only: periodic
    use kt_grids, only: naky
    use stella_layouts, only: vmu_lo, iv_idx
    use parallel_streaming, only: stream_sign
    implicit none

    !    complex, dimension (:,:,-nzgrid:), intent (in out) :: phi, apar
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
    logical, intent (in out) :: restart_time_step

    integer :: ivmu, iv, sgn, iky
    ! start the timer for the explicit part of the solve
    if (proc0) call time_message(.false.,time_gke(:,8),' explicit')

    select case (explicit_option_switch)
    case (explicit_option_rk2)
       ! SSP RK2
       call advance_explicit_rk2 (g, restart_time_step)
    case (explicit_option_rk3)
       ! default is SSP RK3
       call advance_explicit_rk3 (g, restart_time_step)
    case (explicit_option_rk4)
       ! RK4
       call advance_explicit_rk4 (g, restart_time_step)
    end select

    ! enforce periodicity for periodic (including zonal) modes
    do iky = 1, naky
       if (periodic(iky)) then
          do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
             iv = iv_idx(vmu_lo,ivmu)
             ! stream_sign > 0 corresponds to dz/dt < 0
             sgn = stream_sign(iv)
             g(iky,:,sgn*nzgrid,:,ivmu) = g(iky,:,-sgn*nzgrid,:,ivmu)
          end do
       end if
    end do

    ! stop the timer for the explicit part of the solve
    if (proc0) call time_message(.false.,time_gke(:,8),' explicit')

  end subroutine advance_explicit

  !> Advance g by the Leapfrog technique:
  !> g(i+1) = g(i-1) + (dg/dt)_leapfrog(i) * 2dt
  !> where (dg/dt)_leapfrog(i) are terms in GKE evaluated at ti (i.e.
  !> pertaining to gold) - the user can select these terms to be any of the
  !> following:
  !> nonlinear (ExB advection) term
  !> drifts (magnetic drifts + wstar drifts)
  !> NB more terms could be implemented in this way if desired
  !> This subroutine takes g(i-1) and g(i) as inputs; g(i-1) is overwritten
  !> by g(i+1)
  subroutine advance_leapfrog (gold, golder, restart_time_step)

    !use mp, only: proc0, min_allreduce
    !use job_manage, only: time_message
    !use stella_time, only: cfl_dt, code_dt, code_dt_max
    use run_parameters, only: nisl_nonlinear, leapfrog_nonlinear, leapfrog_drifts, exact_exb_nonlinear_solution
    !use physics_parameters, only: g_exb, g_exbfac
    use zgrid, only: nzgrid, ntubes
    use stella_layouts, only: vmu_lo
    use kt_grids, only: naky, nakx
    use fields_arrays, only: phi, apar
    use fields, only: fields_updated, advance_fields
    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: gold
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: golder
    logical, intent (in out)  :: restart_time_step

    complex, dimension (:,:,:,:,:), allocatable :: rhs

    ! if CFL condition is violated by nonlinear term
    ! then must modify time step size and restart time step
    ! assume false and test
    restart_time_step = .false.
    ! Get the fields corresponding to g(i)
    call advance_fields (gold, phi, apar, dist='gbar')
    ! There are 2 different "flavours" of leapfrog option we can use:
    ! (1) Use NISL for the leapfrog step - has to be done as a single
    ! operation (SL scheme)
    !
    ! (2) Use "normal" (explicit) method, whereby (dg/dt) is explicitly evaluated
    if (nisl_nonlinear) then
      if (exact_exb_nonlinear_solution) then
        call advance_ExB_nonlinearity_exact (gold, golder)
      else
        call advance_ExB_nonlinearity_nisl (gold, golder)
      end if
      fields_updated = .false.
    else
      ! Calculate the rhs from the nonlinear and/or drift terms.
      ! There are already methods which calculate this, assuming that the
      ! step to be taking is dt - calculate the RHS using these and then
      ! double to get a step of 2dt.
      allocate (rhs(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        ! Set RHS to zero, then add terms relating to nonlinear and/or wdrift
        ! piece.
        rhs = 0
        ! Calculate the nonlinear term and add to RHS
        if (leapfrog_nonlinear) call advance_ExB_nonlinearity (gold, rhs, restart_time_step)
        if (.not. restart_time_step) then
          if (leapfrog_drifts) then
            ! calculate and add alpha-component of magnetic drift term to RHS of GK eqn
            call advance_wdrifty_explicit (gold, phi, rhs)

            ! calculate and add psi-component of magnetic drift term to RHS of GK eqn
            call advance_wdriftx_explicit (gold, phi, rhs)

            ! calculate and add omega_* term to RHS of GK eqn
            call advance_wstar_explicit (rhs)
          end if
        end if


      ! g(i+1) = g(i-1) + (dg/dt)_leapfrog(i) * 2dt
      golder = golder+2*rhs
      fields_updated = .false.

      deallocate(rhs)
    end if

  end subroutine advance_leapfrog

  ! strong stability-preserving RK2
  subroutine advance_explicit_rk2 (g, restart_time_step)

    use dist_fn_arrays, only: g0, g1
    use zgrid, only: nzgrid
    use stella_layouts, only: vmu_lo
    use multibox, only: RK_step

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
    logical, intent (in out) :: restart_time_step

    if(RK_step) call mb_communicate (g)

    g0 = g

    ! SSP rk3 algorithm to advance explicit part of code
    ! if GK equation written as dg/dt = rhs - vpar . grad h,
    ! solve_gke returns rhs*dt

    ! First step of RK2
    call solve_gke (g0, g1, restart_time_step)
    ! Don't bother to try the subsequent steps if we've reset dt, because
    ! we'll be re-doing the timestep anyway.
    if (.not. restart_time_step) then
      ! Second step of RK2
      g1 = g0 + g1
      if(RK_step) call mb_communicate (g1)
      call solve_gke (g1, g, restart_time_step)

      ! this is gbar at intermediate time level
      g = 0.5*g0 + 0.5*(g1 + g)
    end if

  end subroutine advance_explicit_rk2

  ! strong stability-preserving RK3
  subroutine advance_explicit_rk3 (g, restart_time_step)

    use dist_fn_arrays, only: g0, g1, g2
    use zgrid, only: nzgrid
    use stella_layouts, only: vmu_lo
    use multibox, only: RK_step

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
    logical, intent (in out) :: restart_time_step

    if(RK_step) call mb_communicate (g)

    g0 = g

    ! SSP rk3 algorithm to advance explicit part of code
    ! if GK equation written as dg/dt = rhs - vpar . grad h,
    ! solve_gke returns rhs*dt

    ! First step of RK3
    call solve_gke (g0, g1, restart_time_step)
    if (.not. restart_time_step) then
       ! Second step of RK3
       g1 = g0 + g1
       if(RK_step) call mb_communicate (g1)
       call solve_gke (g1, g2, restart_time_step)

       if (.not. restart_time_step) then
          ! Third step of RK3
          g2 = g1 + g2
          if(RK_step) call mb_communicate (g2)
          call solve_gke (g2, g, restart_time_step)

          ! this is gbar at intermediate time level
          g = g0/3. + 0.5*g1 + (g2 + g)/6.

       end if
    end if

  end subroutine advance_explicit_rk3

  subroutine advance_explicit_rk4 (g, restart_time_step)

    use dist_fn_arrays, only: g0, g1, g2, g3
    use zgrid, only: nzgrid
    use stella_layouts, only: vmu_lo
    use multibox, only: RK_step

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
    logical, intent (in out) :: restart_time_step

    if(RK_step) call mb_communicate(g)

    g0 = g

    ! SSP rk3 algorithm to advance explicit part of code
    ! if GK equation written as dg/dt = rhs - vpar . grad h,
    ! solve_gke returns rhs*dt

    ! First step of RK4
    call solve_gke (g0, g1, restart_time_step)

    if (.not. restart_time_step) then
       ! Step 2 of RK4
       ! g1 is h*k1
       g3 = g0 + 0.5*g1
       if(RK_step) call mb_communicate(g3)
       call solve_gke (g3, g2, restart_time_step)
       g1 = g1 + 2.*g2

       if (.not. restart_time_step) then
          ! Step 3 of RK4
          ! g2 is h*k2
          g2 = g0+0.5*g2
          if(RK_step) call mb_communicate(g2)
          call solve_gke (g2, g3, restart_time_step)
          g1 = g1 + 2.*g3

          if (.not. restart_time_step) then
             ! Step 4 of RK4
             ! g3 is h*k3
             g3 = g0+g3
             if(RK_step) call mb_communicate(g3)
             call solve_gke (g3, g, restart_time_step)
             g1 = g1 + g

            ! this is gbar at intermediate time level
            g = g0 + g1/6.
          end if
       end if
    end if

  end subroutine advance_explicit_rk4

  subroutine solve_gke (gin, rhs_ky, restart_time_step)

    use job_manage, only: time_message
    use fields_arrays, only: phi, apar
    use stella_layouts, only: vmu_lo
    use stella_transforms, only: transform_y2ky
    use redistribute, only: gather, scatter
    use physics_flags, only: include_parallel_nonlinearity
    use physics_flags, only: include_parallel_streaming
    use physics_flags, only: include_mirror
    use physics_flags, only: include_drifts
    use physics_flags, only: nonlinear
    use physics_flags, only: full_flux_surface, radial_variation
    use run_parameters, only: nisl_nonlinear, leapfrog_nonlinear, leapfrog_drifts
    use physics_parameters, only: g_exb
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, ny
    use run_parameters, only: stream_implicit, mirror_implicit, drifts_implicit
    use dissipation, only: include_collisions, advance_collisions_explicit, collisions_implicit
    use dissipation, only: include_krook_operator, add_krook_operator
    use parallel_streaming, only: advance_parallel_streaming_explicit
    use fields, only: advance_fields, fields_updated, get_radial_correction
    use mirror_terms, only: advance_mirror_explicit
    use flow_shear, only: advance_parallel_flow_shear
    use file_utils, only: runtype_option_switch, runtype_multibox
    use multibox, only: include_multibox_krook, add_multibox_krook

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: gin
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (out), target :: rhs_ky
    logical, intent (in out) :: restart_time_step

    complex, dimension (:,:,:,:,:), allocatable, target :: rhs_y
    complex, dimension (:,:,:,:,:), pointer :: rhs

    rhs_ky = 0.

    if (full_flux_surface) then
       allocate (rhs_y(ny,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       rhs_y = 0.
       rhs => rhs_y
    else
       rhs => rhs_ky
    end if
    ! start with gbar in k-space and (ky,kx,z) local
    ! obtain fields corresponding to gbar
    call advance_fields (gin, phi, apar, dist='gbar')

    if(radial_variation) call get_radial_correction(gin, phi, dist='gbar')

    ! default is to continue with same time step size
    ! if estimated CFL condition for nonlinear terms is violated
    ! then restart_time_step will be set to .true.
    restart_time_step = .false.
    ! calculate and add ExB nonlinearity to RHS of GK eqn
    ! do this first, as the CFL condition may require a change in time step
    ! and thus recomputation of mirror, wdrift, wstar, and parstream
    ! We want to calculate the nonlinear term here if the following is true:
    !  (1) nonlinear = .true.
    !  (2) Either (leapfrog_nonlinear = .false. and nisl_nonlinear = .false.) or (leapfrog_this_timestep = .false.)
    if (nonlinear) then
      if ( ((.not. leapfrog_nonlinear) .and. (.not. nisl_nonlinear)) &
        .or. ( .not. leapfrog_this_timestep)) call advance_ExB_nonlinearity (gin, rhs, restart_time_step)
    end if

    if (include_parallel_nonlinearity .and. .not.restart_time_step) &
         call advance_parallel_nonlinearity (gin, rhs, restart_time_step)

    if (.not. restart_time_step) then

       if ((g_exb**2).gt.epsilon(0.0)) call advance_parallel_flow_shear (rhs)

       ! calculate and add mirror term to RHS of GK eqn
       if (include_mirror.and..not.mirror_implicit) then
          call advance_mirror_explicit (gin, rhs)
       end if

       if ((include_drifts .and. .not. drifts_implicit) .and. (.not. leapfrog_drifts .or. .not. leapfrog_this_timestep)) then
         ! calculate and add alpha-component of magnetic drift term to RHS of GK eqn
         call advance_wdrifty_explicit (gin, phi, rhs)

         ! calculate and add psi-component of magnetic drift term to RHS of GK eqn
         call advance_wdriftx_explicit (gin, phi, rhs)

         ! calculate and add omega_* term to RHS of GK eqn
         call advance_wstar_explicit (rhs)
       end if

       if (include_collisions.and..not.collisions_implicit) call advance_collisions_explicit (gin, phi, rhs)

       if (full_flux_surface) then
          call transform_y2ky (rhs_y, rhs_ky)
          deallocate (rhs_y)
       end if

       ! calculate and add parallel streaming term to RHS of GK eqn
       if (include_parallel_streaming.and.(.not.stream_implicit)) &
            call advance_parallel_streaming_explicit (gin, rhs_ky)

       if (radial_variation) call advance_radial_variation(gin,rhs)

       if (include_krook_operator) call add_krook_operator(gin,rhs)

       if (runtype_option_switch == runtype_multibox .and. include_multibox_krook) &
         call add_multibox_krook(gin,rhs)

       ! calculate and add omega_* term to RHS of GK eqn
!       if (wstar_explicit) call advance_wstar_explicit (rhs_ky)
!       call advance_wstar_explicit (rhs_ky)
       ! calculate and add collision term to RHS of GK eqn
       !    call advance_collisions
    end if

    fields_updated = .false.

    nullify (rhs)

  end subroutine solve_gke

  subroutine advance_wstar_explicit (gout)

    use mp, only: proc0, mp_abort
    use job_manage, only: time_message
    use fields, only: get_dchidy
    use fields_arrays, only: phi, apar
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx
    use physics_flags, only: full_flux_surface

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout

    complex, dimension (:,:,:,:,:), allocatable :: g0

    if (proc0) call time_message(.false.,time_gke(:,6),' wstar advance')

    allocate (g0(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    ! omega_* stays in ky,kx,z space with ky,kx,z local
    ! get d<chi>/dy
    if (debug) write (*,*) 'time_advance::solve_gke::get_dchidy'
    if (full_flux_surface) then
       ! FLAG -- ADD SOMETHING HERE
       call mp_abort ('wstar term not yet setup for full_flux_surface = .true. aborting.')
    else
       call get_dchidy (phi, apar, g0)
       ! multiply with omega_* coefficient and add to source (RHS of GK eqn)
       if (debug) write (*,*) 'time_advance::solve_gke::add_wstar_term'
       call add_wstar_term (g0, gout)
    end if
    deallocate (g0)

    if (proc0) call time_message(.false.,time_gke(:,6),' wstar advance')

  end subroutine advance_wstar_explicit

  subroutine advance_wdrifty_explicit (g, phi, gout)

    use mp, only: proc0
    use stella_layouts, only: vmu_lo
    use job_manage, only: time_message
    use stella_transforms, only: transform_ky2y
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky, ny
    use physics_flags, only: full_flux_surface
    use gyro_averages, only: gyro_average
    use dist_fn_arrays, only: wdrifty_g, wdrifty_phi

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phi
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout

    integer :: ivmu
    complex, dimension (:,:,:,:), allocatable :: dphidy
    complex, dimension (:,:,:,:,:), allocatable :: g0k, g0y

    ! alpha-component of magnetic drift (requires ky -> y)
    if (proc0) call time_message(.false.,time_gke(:,4),' dgdy advance')

    allocate (dphidy(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (g0k(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    if (debug) write (*,*) 'time_advance::solve_gke::get_dgdy'
    call get_dgdy (g, g0k)
    call get_dgdy (phi, dphidy)

    if (full_flux_surface) then
       allocate (g0y(ny,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       ! transform dg/dy from k-space to y-space
       call transform_ky2y (g0k, g0y)
       ! add vM . grad y dg/dy term to equation
       call add_dg_term_global (g0y, wdrifty_g, gout)
       ! get <dphi/dy> in k-space
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          call gyro_average (dphidy, ivmu, g0k(:,:,:,:,ivmu))
       end do
       ! transform d<phi>/dy from k-space to y-space
       call transform_ky2y (g0k, g0y)
       ! add vM . grad y d<phi>/dy term to equation
       call add_dphi_term_global (g0y, wdrifty_phi, gout)
       deallocate (g0y)
    else
       if (debug) write (*,*) 'time_advance::solve_gke::add_dgdy_term'
       ! add vM . grad y dg/dy term to equation
       call add_dg_term (g0k, wdrifty_g(1,:,:), gout)
       ! get <dphi/dy> in k-space
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          call gyro_average (dphidy, ivmu, g0k(:,:,:,:,ivmu))
       end do
       ! add vM . grad y d<phi>/dy term to equation
       call add_dphi_term (g0k, wdrifty_phi(1,:,:), gout)
    end if
    deallocate (g0k, dphidy)

    if (proc0) call time_message(.false.,time_gke(:,4),' dgdy advance')

  end subroutine advance_wdrifty_explicit

  subroutine advance_wdriftx_explicit (g, phi, gout)

    use mp, only: proc0
    use stella_layouts, only: vmu_lo
    use job_manage, only: time_message
    use stella_transforms, only: transform_ky2y
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky, ny, akx
    use physics_flags, only: full_flux_surface
    use gyro_averages, only: gyro_average
    use dist_fn_arrays, only: wdriftx_g, wdriftx_phi

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phi
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout

    integer :: ivmu
    complex, dimension (:,:,:,:), allocatable :: dphidx
    complex, dimension (:,:,:,:,:), allocatable :: g0k, g0y

    ! psi-component of magnetic drift (requires ky -> y)
    if (proc0) call time_message(.false.,time_gke(:,5),' dgdx advance')

    ! do not calculate if wdriftx terms are all zero
    if (maxval(abs(akx))<epsilon(0.)) then
       if (proc0) call time_message(.false.,time_gke(:,5),' dgdx advance')
       return
    end if

    allocate (dphidx(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (g0k(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    if (debug) write (*,*) 'time_advance::solve_gke::get_dgdx'
    call get_dgdx (g, g0k)
    call get_dgdx (phi, dphidx)

    if (full_flux_surface) then
       allocate (g0y(ny,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       ! transform dg/dx from k-space to y-space
       call transform_ky2y (g0k, g0y)
       ! add vM . grad x dg/dx term to equation
       call add_dg_term_global (g0y, wdriftx_g, gout)
       ! get <dphi/dx> in k-space
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          call gyro_average (dphidx, ivmu, g0k(:,:,:,:,ivmu))
       end do
       ! transform d<phi>/dx from k-space to y-space
       call transform_ky2y (g0k, g0y)
       ! add vM . grad x d<phi>/dx term to equation
       call add_dphi_term_global (g0y, wdriftx_phi, gout)
       deallocate (g0y)
    else
       if (debug) write (*,*) 'time_advance::solve_gke::add_dgdx_term'
       ! add vM . grad x dg/dx term to equation
       call add_dg_term (g0k, wdriftx_g(1,:,:), gout)
       ! get <dphi/dx> in k-space
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          call gyro_average (dphidx, ivmu, g0k(:,:,:,:,ivmu))
       end do
       ! add vM . grad x d<phi>/dx term to equation
       call add_dphi_term (g0k, wdriftx_phi(1,:,:), gout)
    end if
    deallocate (g0k, dphidx)

    if (proc0) call time_message(.false.,time_gke(:,5),' dgdx advance')

  end subroutine advance_wdriftx_explicit

  subroutine advance_ExB_nonlinearity (g, gout, restart_time_step)

    use mp, only: proc0, min_allreduce
    use mp, only: scope, allprocs, subprocs
    use stella_layouts, only: vmu_lo, imu_idx, is_idx
    use job_manage, only: time_message
    use gyro_averages, only: gyro_average
    use fields, only: get_dchidx, get_dchidy
    use fields_arrays, only: phi, apar, shift_state
    use fields_arrays, only: phi_corr_QN,   phi_corr_GA
    !   use fields_arrays, only: apar_corr_QN, apar_corr_GA
    use stella_transforms, only: transform_y2ky,transform_x2kx
    use stella_transforms, only: transform_y2ky_xfirst, transform_x2kx_xfirst
    use stella_time, only: cfl_dt, code_dt, code_dt_max
    use run_parameters, only: cfl_cushion, delt_adjust, fphi, nisl_nonlinear
    use physics_parameters, only: g_exb, g_exbfac
    use zgrid, only: nzgrid, ntubes
    use stella_geometry, only: exb_nonlin_fac, exb_nonlin_fac_p, gfac
    use kt_grids, only: nakx, naky, nx, ny, ikx_max
    use kt_grids, only: akx, aky, rho_clamped
    use physics_flags, only: full_flux_surface, radial_variation
    use physics_flags, only: prp_shear_enabled, hammett_flow_shear
    use physics_flags, only: override_vexb, vexb_x, vexb_y
    use kt_grids, only: x, y, swap_kxky, swap_kxky_back
    use constants, only: pi, zi
    use file_utils, only: runtype_option_switch, runtype_multibox

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout
    logical, intent (out) :: restart_time_step


    complex, dimension (:,:), allocatable :: g0k, g0a, g0k_swap
    complex, dimension (:,:), allocatable :: g0kxy, g0xky, prefac
    real, dimension (:,:), allocatable :: g0xy, g1xy, bracket

    integer :: ivmu, iz, it, ia, imu, is
    logical :: yfirst

    ! alpha-component of magnetic drift (requires ky -> y)
    if (proc0) call time_message(.false.,time_gke(:,7),' ExB nonlinear advance')
    if (debug) write (*,*) 'time_advance::solve_gke::advance_ExB_nonlinearity::get_dgdy'

    restart_time_step = .false.
    yfirst = .not.prp_shear_enabled
    allocate (g0k(naky,nakx))
    allocate (g0a(naky,nakx))
    allocate (g0xy(ny,nx))
    allocate (g1xy(ny,nx))
    allocate (bracket(ny,nx))
    allocate (prefac(naky,nx))

    if(yfirst) then
      allocate (g0k_swap(2*naky-1,ikx_max))
      allocate (g0kxy(ny,ikx_max))
    else
      allocate (g0xky(naky,nx))
    endif

    prefac = 1.0
    if(prp_shear_enabled.and.hammett_flow_shear) then
      prefac = exp(-zi*g_exb*g_exbfac*spread(x,1,naky)*spread(aky*shift_state,2,nx))
    endif

    ia=1
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do it = 1, ntubes
          do iz = -nzgrid, nzgrid
             call get_dgdy (g(:,:,iz,it,ivmu), g0k)
             call forward_transform(g0k,g0xy)

             call get_dchidx (iz, ivmu, phi(:,:,iz,it), apar(:,:,iz,it), g0k)
             if(prp_shear_enabled.and.hammett_flow_shear) then
               call get_dchidy (iz, ivmu, phi(:,:,iz,it), apar(:,:,iz,it), g0a)
               g0k = g0k - g_exb*g_exbfac*spread(shift_state,2,nakx)*g0a
             endif
             call forward_transform(g0k,g1xy)
             if (override_vexb) then
               g1xy = vexb_y
               bracket = - g0xy*g1xy ! -ve sign because we're on the RHS
             else
               g1xy = g1xy*exb_nonlin_fac
               bracket = g0xy*g1xy
             end if

             cfl_dt = min(cfl_dt,2.*pi/(maxval(abs(g1xy))*aky(naky)))

             if(radial_variation) then
               bracket = bracket + gfac*g0xy*g1xy*exb_nonlin_fac_p*spread(rho_clamped,1,ny)
               call gyro_average (phi_corr_QN(:,:,iz,it),iz,ivmu,g0a)
               g0a = fphi*(g0a + phi_corr_GA(:,:,iz,it,ivmu))
               call get_dgdx(g0a,g0k)
               call forward_transform(g0k,g1xy)
               g1xy = g1xy*exb_nonlin_fac
               bracket = bracket + g0xy*g1xy
             endif

             cfl_dt = min(cfl_dt,2.*pi/(maxval(abs(g1xy))*aky(naky)))

             call get_dgdx (g(:,:,iz,it,ivmu), g0k)
             if(prp_shear_enabled.and.hammett_flow_shear) then
               call get_dgdy (g(:,:,iz,it,ivmu), g0a)
               g0k = g0k - g_exb*g_exbfac*spread(shift_state,2,nakx)*g0a
             endif
             call forward_transform(g0k,g0xy)
             call get_dchidy (iz, ivmu, phi(:,:,iz,it), apar(:,:,iz,it), g0k)
             call forward_transform(g0k,g1xy)
             if (override_vexb) then
               g1xy = vexb_x
               bracket = bracket - g0xy*g1xy ! -ve sign because we're on the RHS
             else
               g1xy = g1xy*exb_nonlin_fac
               bracket = bracket - g0xy*g1xy
             end if

             cfl_dt = min(cfl_dt,2.*pi/(maxval(abs(g1xy))*akx(ikx_max)))

             if(radial_variation) then
               bracket = bracket - gfac*g0xy*g1xy*exb_nonlin_fac_p*spread(rho_clamped,1,ny)
               call gyro_average (phi_corr_QN(:,:,iz,it),iz,ivmu,g0a)
               g0a = fphi*(g0a + phi_corr_GA(:,:,iz,it,ivmu))
               call get_dgdy(g0a,g0k)
               call forward_transform(g0k,g1xy)
               g1xy = g1xy*exb_nonlin_fac
               bracket = bracket - g0xy*g1xy
             endif

             cfl_dt = min(cfl_dt,2.*pi/(maxval(abs(g1xy))*akx(ikx_max)))

             if(yfirst) then
                call transform_x2kx (bracket, g0kxy)
                if (full_flux_surface) then
                  gout(:,:,iz,it,ivmu) = g0kxy
                else
                  call transform_y2ky (g0kxy, g0k_swap)
                  call swap_kxky_back (g0k_swap, gout(:,:,iz,it,ivmu))
                endif
             else
                call transform_y2ky_xfirst (bracket, g0xky)
                g0xky = g0xky/prefac
                call transform_x2kx_xfirst (g0xky, gout(:,:,iz,it,ivmu))
             end if
          end do
       end do
       ! enforce periodicity for zonal mode
       ! FLAG -- THIS IS PROBABLY NOT NECESSARY (DONE AT THE END OF EXPLICIT ADVANCE)
       ! AND MAY INDEED BE THE WRONG THING TO DO
       gout(1,:,-nzgrid,:,ivmu) = 0.5*(gout(1,:,nzgrid,:,ivmu)+gout(1,:,-nzgrid,:,ivmu))
       gout(1,:,nzgrid,:,ivmu) = gout(1,:,-nzgrid,:,ivmu)
    end do

    deallocate (g0k, g0a, g0xy, g1xy, bracket)
    if (allocated(g0k_swap)) deallocate(g0k_swap)
    if (allocated(g0xky)) deallocate(g0xky)
    if (allocated(g0kxy)) deallocate(g0kxy)

    if(runtype_option_switch == runtype_multibox) call scope(allprocs)

    call min_allreduce (cfl_dt)
    if(runtype_option_switch == runtype_multibox) call scope(subprocs)


    if (code_dt > cfl_dt*cfl_cushion) then
      if (.not. nisl_nonlinear) then
        if (proc0) then
          write (*,*) ' '
          write (*,*) 'CHANGING TIME STEP:'
          write (*,'(A16, ES10.2E2)') "   code_dt:"//REPEAT(' ',50),code_dt
          write (*,'(A16, ES10.2E2)') "   cfl_dt:"//REPEAT(' ',50),cfl_dt
          write (*,'(A16, ES10.2E2)') "   cfl_cushion:"//REPEAT(' ',50),cfl_cushion
          write (*,'(A16, ES10.2E2)') "   delt_adjust:"//REPEAT(' ',50),delt_adjust
          write (*,'(A65)') '     ==> User-specified delt is larger than cfl_dt*cfl_cushion.'//REPEAT(' ',50)
          write (*,'(A61,ES12.4)') '     ==> Changing code_dt to cfl_dt*cfl_cushion/delt_adjust ='//REPEAT(' ',50), cfl_dt*cfl_cushion/delt_adjust
          write(*,*) ' '
        end if
        code_dt = cfl_dt*cfl_cushion/delt_adjust
      else
        !!! If we're using NISL, then we hit this branch only when we start (when we're
        !!! taking the first timesteps explictly to start the Leapfrog-like NISL.)
        !!! Rather than setting code_dt = cfl_dt*cfl_cushion/delt_adjust, we want to
        !!! set it to so that old_dt = N*new_dt for integer N.
        ! Want new_dt <= cfl_dt*cfl_cushion and old_dt = N*new_dt for integer N
        !max_dt = cfl_dt*cfl_cushion
        code_dt = code_dt / ceiling(code_dt/(cfl_dt*cfl_cushion))
      end if

      call reset_dt
      restart_time_step = .true.

    else if (code_dt < min(cfl_dt*cfl_cushion/delt_adjust,code_dt_max)) then
      ! If NISL, we don't want this to adjust
      if (.not. nisl_nonlinear) then
        if (proc0) then
          write (*,*) ' '
          write (*,*) 'CHANGING TIME STEP:'
          write (*,'(A16, ES10.2E2)') "   code_dt:"//REPEAT(' ',50),code_dt
          write (*,'(A16, ES10.2E2)') "   cfl_dt:"//REPEAT(' ',50),cfl_dt
          write (*,'(A16, ES10.2E2)') "   cfl_cushion:"//REPEAT(' ',50),cfl_cushion
          write (*,'(A16, ES10.2E2)') "   delt_adjust:"//REPEAT(' ',50),delt_adjust
          write (*,'(A65)') '     ==> User-specified delt is larger than cfl_dt*cfl_cushion.'//REPEAT(' ',50)
          write (*,'(A61,ES12.4)') '     ==> Changing code_dt to cfl_dt*cfl_cushion/delt_adjust ='//REPEAT(' ',50), cfl_dt*cfl_cushion/delt_adjust
          write(*,*) ' '
        end if
         code_dt = min(cfl_dt*cfl_cushion/delt_adjust,code_dt_max)
         call reset_dt
         ! FLAG -- NOT SURE THIS IS CORRECT
         gout = code_dt*gout
      end if
    else
       gout = code_dt*gout
    end if

    if (proc0) call time_message(.false.,time_gke(:,7),' ExB nonlinear advance')

    contains

    subroutine forward_transform (gk,gx)

      use stella_transforms, only: transform_ky2y, transform_kx2x
      use stella_transforms, only: transform_ky2y_xfirst, transform_kx2x_xfirst

      implicit none

      complex, dimension(:,:), intent (in) :: gk
      real, dimension(:,:), intent (out) :: gx

      if(yfirst) then
        ! we have i*ky*g(kx,ky) for ky >= 0 and all kx
        ! want to do 1D complex to complex transform in y
        ! which requires i*ky*g(kx,ky) for all ky and kx >= 0
        ! use g(kx,-ky) = conjg(g(-kx,ky))
        ! so i*(-ky)*g(kx,-ky) = -i*ky*conjg(g(-kx,ky)) = conjg(i*ky*g(-kx,ky))
        ! and i*kx*g(kx,-ky) = i*kx*conjg(g(-kx,ky)) = conjg(i*(-kx)*g(-kx,ky))
        ! and i*(-ky)*J0(kx,-ky)*phi(kx,-ky) = conjg(i*ky*J0(-kx,ky)*phi(-kx,ky))
        ! and i*kx*J0(kx,-ky)*phi(kx,-ky) = conjg(i*(-kx)*J0(-kx,ky)*phi(-kx,ky))
        ! i.e., can calculate dg/dx, dg/dy, d<phi>/dx and d<phi>/dy
        ! on stella (kx,ky) grid, then conjugate and flip sign of (kx,ky)
        ! NB: J0(kx,ky) = J0(-kx,-ky)
        ! TODO DSO: coordinate change for shearing
        call swap_kxky (gk, g0k_swap)
        call transform_ky2y (g0k_swap, g0kxy)
        call transform_kx2x (g0kxy, gx)
      else
        call transform_kx2x_xfirst(gk, g0xky)
        g0xky = g0xky*prefac
        call transform_ky2y_xfirst(g0xky, gx)
      endif

    end subroutine forward_transform

  end subroutine advance_ExB_nonlinearity



  !> A time-advanced golder is passed in, along with gold (which hasn't been
  !> evolved at all). We apply a NISL-Leapfrog step to replace golder with
  !> nonlinearly-advanced golder, that is we update golder as follows:
  !>
  !>  golder(x,y) = iFFt(golder(kx,ky))
  !>  golder(x,y) = golder(x-p*dx,y-q*dy) + 2*rhs{gold(x-p*dx/2,y-q*dy/2)}
  !>  golder(kx,ky) = FFT(golder(x,y))
  !>
  !> NB it's likely that sequentially advecting in x,y leads to O(dt) errors,
  !> so advect in both directions at the same time.
  subroutine advance_ExB_nonlinearity_exact (gold, golder)

    use mp, only: proc0, min_allreduce
    use mp, only: scope, allprocs, subprocs
    use stella_layouts, only: vmu_lo, imu_idx, is_idx
    use job_manage, only: time_message
    use gyro_averages, only: gyro_average
    use fields, only: get_dchidx, get_dchidy
    use fields_arrays, only: phi, apar, shift_state
    use stella_transforms, only: transform_y2ky,transform_x2kx
    use stella_transforms, only: transform_y2ky_xfirst, transform_x2kx_xfirst
    use stella_time, only: code_dt
    use run_parameters, only: fphi
    use physics_flags, only: override_vexb, vexb_x, vexb_y
    use physics_parameters, only: g_exb, g_exbfac
    use run_parameters, only: add_nl_source_in_real_space, no_extra_padding
    use zgrid, only: nzgrid, ntubes
    use stella_geometry, only: exb_nonlin_fac, gfac, dxdXcoord, dydalpha
    use kt_grids, only: nakx, naky, nx, ny, ikx_max
    use kt_grids, only: akx, aky, rho_clamped
    use kt_grids, only: x, y, swap_kxky, swap_kxky_back
    use kt_grids, only: dx, dy, x0, y0
    use constants, only: pi, zi
    use file_utils, only: runtype_option_switch, runtype_multibox

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: gold
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: golder

    integer :: ivmu, iz, it, ia, imu, is, ix, iy, p, q, iky, ikx
    ia=1
    if (override_vexb) then
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        imu = imu_idx(vmu_lo,ivmu)
        is = is_idx(vmu_lo,ivmu)
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            ! Hack - get analytic solution
            ! gnew = gold*exp(-i(k_x U_x + k_y U_y)t)
            do iky = 1, naky
              do ikx = 1, nakx
                golder(iky, ikx, iz, it, ivmu) = gold(iky, ikx, iz, it, ivmu)*exp(cmplx(0,-1)*(akx(ikx)*vexb_x+aky(iky)*vexb_y)*code_dt)
              end do
            end do
          end do
        end do
      end do
    else
      stop "Analytic solution impossible for non-constant vexb"
    end if

  end subroutine advance_ExB_nonlinearity_exact


  !> A time-advanced golder is passed in, along with gold (which hasn't been
  !> evolved at all). We apply a NISL-Leapfrog step to replace golder with
  !> nonlinearly-advanced golder, that is we update golder as follows:
  !>
  !>  golder(x,y) = iFFt(golder(kx,ky))
  !>  golder(x,y) = golder(x-p*dx,y-q*dy) + 2*rhs{gold(x-p*dx/2,y-q*dy/2)}
  !>  golder(kx,ky) = FFT(golder(x,y))
  !>
  !> NB it's likely that sequentially advecting in x,y leads to O(dt) errors,
  !> so advect in both directions at the same time.
  subroutine advance_ExB_nonlinearity_nisl (gold, golder, single_step)

    use mp, only: proc0, min_allreduce
    use mp, only: scope, allprocs, subprocs
    use stella_layouts, only: vmu_lo, imu_idx, is_idx
    use job_manage, only: time_message
    use gyro_averages, only: gyro_average
    use fields, only: get_dchidx, get_dchidy
    use fields_arrays, only: phi, apar, shift_state
    use stella_transforms, only: transform_y2ky,transform_x2kx
    use stella_transforms, only: transform_y2ky_xfirst, transform_x2kx_xfirst
    use stella_time, only: code_dt
    use run_parameters, only: fphi
    use physics_flags, only: override_vexb, vexb_x, vexb_y
    use physics_parameters, only: g_exb, g_exbfac
    use run_parameters, only: add_nl_source_in_real_space, no_extra_padding
    use zgrid, only: nzgrid, ntubes
    use stella_geometry, only: exb_nonlin_fac, gfac, dxdXcoord, dydalpha
    use kt_grids, only: nakx, naky, nx, ny, ikx_max
    use kt_grids, only: akx, aky, rho_clamped
    use kt_grids, only: x, y, swap_kxky, swap_kxky_back
    use kt_grids, only: dx, dy, x0, y0
    use constants, only: pi, zi
    use file_utils, only: runtype_option_switch, runtype_multibox

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: gold
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: golder
    logical, optional, intent(in) :: single_step  ! First timestep or from restart needs the single-step version

    !complex, dimension (:,:,:,:,:), allocatable :: gout, g
    complex, dimension (:,:), allocatable :: g0k, g0k_swap, rhs_array_fourier !, g0k_swap_extra_padding
    complex, dimension (:,:), allocatable :: g0kxy, g0kxy_extra_padding
    real, dimension (:,:), allocatable :: vchiold_x, vchiold_y, dgold_dy, dgold_dx, golderxy, gnewxy, dgold_dx_normal, dgold_dy_normal, rhs_array!, bracket

    integer :: ivmu, iz, it, ia, imu, is, ix, iy, p, q, xidx_departure, yidx_departure, yidx_for_upsampled_array, xidx_for_upsampled_array
    integer :: upsampled_xidx, upsampled_yidx
    logical :: yfirst
    logical :: single_step_local
    real :: y_departure, x_departure, yval, xval, rhs, velocity_x, velocity_y
    real :: max_velocity_x, max_velocity_y
    ! alpha-component of magnetic drift (requires ky -> y)
    if (proc0) call time_message(.false.,time_gke(:,7),' ExB nonlinear advance, nisl')

    if (debug) write (*,*) 'time_advance::solve_gke::advance_ExB_nonlinearity::get_dgdy'

    ! Assume no perp shear; yfirst = true

    single_step_local = .false.
    if(present(single_step)) single_step_local = single_step

    allocate (g0k(naky,nakx))
    allocate (rhs_array_fourier(naky,nakx))
    allocate (dgold_dx(2*ny,2*nx))
    if (no_extra_padding) then
      allocate (dgold_dx_normal(ny,nx))
      allocate (dgold_dy_normal(ny,nx))
    end if
    allocate (dgold_dy(2*ny,2*nx))
    allocate (vchiold_x(ny,nx))
    allocate (vchiold_y(ny,nx))
    allocate (golderxy(ny,nx))
    allocate (gnewxy(ny,nx))
    allocate (rhs_array(ny,nx))
    !allocate (bracket(ny,nx))

    allocate (g0k_swap(2*naky-1,ikx_max))
    allocate (g0kxy(ny,ikx_max))
    allocate (g0kxy_extra_padding(2*ny,ikx_max))

    ia=1

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      imu = imu_idx(vmu_lo,ivmu)
      is = is_idx(vmu_lo,ivmu)
      do it = 1, ntubes
        do iz = -nzgrid, nzgrid

          ! Need to do the following:
          ! 1) Get dgold/dx(x,y), dgold/dy(x,y), vchiold_x(x,y), vchiold_y(x,y)
          ! 2) Get golder(x,y)
          ! 3) Calculate departure points, and hence velocities and p, q
          ! 4) golder(x,y) = golder(x-p*dx,y-q*dy) + 2*rhs{gold(x-p*dx/2,y-q*dy/2)}
          ! 5) F.T. golder to get golder(kx,ky)

          ! Step 1) Get dgold/dx(x,y), dgold/dy(x,y), vchiold_x(x,y), vchiold_y(x,y)
          call get_dgdy (gold(:,:,iz,it,ivmu), g0k)
          if (no_extra_padding) then
            call forward_transform(g0k,dgold_dy_normal)
          else
            call forward_transform_extra_padding(g0k,dgold_dy)
          end if
          ! For testing - get dg/dy in the non-upsampled and upsampled forms and print
          ! call forward_transform(g0k,golderxy)

          call get_dchidx (iz, ivmu, phi(:,:,iz,it), apar(:,:,iz,it), g0k)
          call forward_transform(g0k,vchiold_y)

          ! NB we haven't flipped the signs (in contrast to advance_exb_nonlinearity)
          ! because this is the actual velocity - we haven't moved to the RHS yet.
          vchiold_y = -vchiold_y*exb_nonlin_fac

          call get_dgdx (gold(:,:,iz,it,ivmu), g0k)
          if (no_extra_padding) then
            call forward_transform(g0k,dgold_dx_normal)
          else
            call forward_transform_extra_padding(g0k,dgold_dx)
          end if
          call get_dchidy (iz, ivmu, phi(:,:,iz,it), apar(:,:,iz,it), g0k)
          call forward_transform(g0k,vchiold_x)

          vchiold_x = vchiold_x*exb_nonlin_fac

          ! Step 2) Get golder(x,y)
          call forward_transform(golder(:,:,iz,it,ivmu),golderxy)

          if (single_step_local) then
            max_velocity_x = dx/(2*code_dt)
            max_velocity_y = dy/(2*code_dt)
          else
            max_velocity_x = dx/(4*code_dt)
            max_velocity_y = dy/(4*code_dt)
          end if

          if (override_vexb) then
            vchiold_y = vexb_y
            vchiold_x = vexb_x
          end if
          ! Step 3) Calculate departure points, and hence velocities and p, q
          do iy = 1, ny
            do ix = 1, nx
              if (override_vexb) then
                velocity_x = vexb_x
                velocity_y = vexb_y
                x_departure = x(ix) - 2*velocity_x*code_dt
                y_departure = y(iy) - 2*velocity_y*code_dt
              else
                call get_approx_departure_point(vchiold_y, vchiold_x, iy, ix, y_departure, x_departure, single_step_local)
                velocity_x = (x(ix) - x_departure)/(2*code_dt)
                velocity_y = (y(iy) - y_departure)/(2*code_dt)
              end if

              ! Check velocities are sensible
              ! Is this tolerance too high?
              if (override_vexb) then
                if (abs(velocity_x - vexb_x) > 1E-10) then
                  write(*,*) "velocity_x, vexb_x = ", velocity_x, vexb_x
                  stop "Aborting"
                end if
                if (abs(velocity_y - vexb_y) > 1E-10) then
                  write(*,*) "velocity_y, vexb_y = ", velocity_y, vexb_y
                  stop "Aborting"
                end if
              end if
              p = nint((x(ix) - x_departure)/dx)
              q = nint((y(iy) - y_departure)/dy)
              ! if ((p /= 0) .or. (q /=0)) then
              !   stop "Not equal to zero!"
              !
              ! end if
              if (override_vexb) then
                if (p == 0) then
                  velocity_x = vexb_x
                end if

                if (q == 0) then
                  velocity_y = vexb_y
                end if
              end if

              ! Find idxs of gridpoint closest to departure point.
              ! Normal idxs must be in the range(1,2,3,...,nx(or ny))
              ! Upsampled idxs must be in the range(1,2,3,...,2*nx(or 2*ny))
              xidx_departure = modulo((ix - p - 1),nx) + 1
              yidx_departure = modulo((iy - q - 1),ny) + 1
              if (single_step_local) then
                ! These idxs are at the time location told, which is the same
                ! as the time location of the g which we're evolving (although it's
                ! called golder). So e.g. for ix=10, p=1, we'd want xidx_departure=9 and
                ! xidx_for_upsampled_array = 18.
                ! So the formula is xidx_for_upsampled_array=(2*ix - 2*p), but correctly moduloed such that
                ! it lies in the range (1,2,3,...,2*nx)
                xidx_for_upsampled_array = modulo((2*ix - 2*p - 2),(2*nx)) + 1
                yidx_for_upsampled_array = modulo((2*iy - 2*q - 2),(2*ny)) + 1
              else
                ! These idxs are halfway between the 2 time locations (tolder, tnew)
                ! So e.g. for ix=10, p=1, we'd want xidx_departure=9 and
                ! xidx_for_upsampled_array = 19, because the upsampled idxs
                ! for ix=10 and 9 are 20 and 18 respectively.
                ! So the formula is xidx_for_upsampled_array=(2*ix - p), but correctly moduloed such that
                ! it lies in the range (1,2,3,...,2*nx)
                xidx_for_upsampled_array = modulo((2*ix - p - 2),(2*nx)) + 1
                yidx_for_upsampled_array = modulo((2*iy - q - 2),(2*ny)) + 1
              end if

              ! Residual velocities
              if (single_step_local) then
                velocity_x = velocity_x - p*dx/(code_dt)
                velocity_y = velocity_y - q*dy/(code_dt)
              else
                velocity_x = velocity_x - p*dx/(2*code_dt)
                velocity_y = velocity_y - q*dy/(2*code_dt)
              end if

              ! Can we check velocity_x, velocity_y to see if they're breaking a CFL condition?
              ! Are the residual velocities larger than expected?
              if (abs(velocity_x) > max_velocity_x) then
                write(*,*) "velocity_x, max_velocity_x =  ", velocity_x, max_velocity_x
                write(*,*) "velocity_y, max_velocity_y =  ", velocity_y, max_velocity_y
                stop "Stopping early"
              end if
              if (abs(velocity_y) > max_velocity_y) then
                write(*,*) "velocity_x, max_velocity_x =  ", velocity_x, max_velocity_x
                write(*,*) "velocity_y, max_velocity_y =  ", velocity_y, max_velocity_y
                stop "Stopping early"
              end if
              !! Update golder(x,y)
              if (single_step_local) then
                if (no_extra_padding) then
                  write(*,*) "Options no_extra_padding,  not yet implemented for the single-step version"
                  stop "Stopping early"
                end if
                !  The single-step version: We only have one timestep to work with.
                rhs = - code_dt * (velocity_x * dgold_dx(yidx_for_upsampled_array, xidx_for_upsampled_array) &
                        + velocity_y * dgold_dy(yidx_for_upsampled_array, xidx_for_upsampled_array) )
              else
                ! The 3-step version

                !!! Original implemenation  - calculate RHS, then update g on an idx-by-idx basis

                if (.not. no_extra_padding) then
                  !! "Correct" RHS implementation; use upsampled dg/dx, dg/dy
                  ! -ve sign because it's on the RHS
                  rhs = - 2*code_dt*(velocity_x * dgold_dx(yidx_for_upsampled_array, xidx_for_upsampled_array) &
                          + velocity_y * dgold_dy(yidx_for_upsampled_array, xidx_for_upsampled_array) )
                else
                  !! Implementation which should be equivalent to the "Leapfrog" implementation
                  !! Only valid if no SL
                  if ((p .ne. 0) .or. (q .ne. 0)) then
                    write(*,*) "p, q = ", p, q
                    stop "stopping"
                  end if
                  if ((yidx_departure .ne. iy) .or. (xidx_departure .ne. ix) ) then
                    write(*,*) "yidx_departure, iy, xidx_departure, ix = ", yidx_departure, iy, xidx_departure, ix
                    stop "stopping"
                  end if
                  rhs = -2*code_dt*(velocity_x * dgold_dx_normal(iy, ix) &
                          + velocity_y * dgold_dy_normal(iy, ix) )
                end if
              end if
              if (add_nl_source_in_real_space) then
                gnewxy(iy, ix) = golderxy(yidx_departure, xidx_departure) + rhs
              else
                rhs_array(iy,ix) = rhs
                gnewxy(iy, ix) = golderxy(yidx_departure, xidx_departure)
              end if

            end do
          end do

          !! Invert to get golder(kx,ky)
          call transform_x2kx (gnewxy, g0kxy)
          call transform_y2ky (g0kxy, g0k_swap)
          call swap_kxky_back (g0k_swap, golder(:,:,iz,it,ivmu))

          !! Invert rhs_array to get rhs_array in Fourier space
          if (.not. add_nl_source_in_real_space) then
            call transform_x2kx (rhs_array, g0kxy)
            call transform_y2ky (g0kxy, g0k_swap)
            call swap_kxky_back (g0k_swap, rhs_array_fourier)
            golder(:,:,iz,it,ivmu) = golder(:,:,iz,it,ivmu) + rhs_array_fourier
          end if
        end do
      end do
      ! ! enforce periodicity for zonal mode
      ! ! FLAG -- THIS IS PROBABLY NOT NECESSARY (DONE AT THE END OF EXPLICIT ADVANCE)
      ! ! AND MAY INDEED BE THE WRONG THING TO DO
      golder(1,:,-nzgrid,:,ivmu) = 0.5*(golder(1,:,nzgrid,:,ivmu)+golder(1,:,-nzgrid,:,ivmu))
      golder(1,:,nzgrid,:,ivmu) = golder(1,:,-nzgrid,:,ivmu)
    end do

    deallocate (g0k, dgold_dx, dgold_dy, vchiold_x, vchiold_y, golderxy)
    if (allocated(g0k_swap)) deallocate(g0k_swap)
    if (allocated(g0kxy)) deallocate(g0kxy)
    if (allocated(rhs_array_fourier)) deallocate(rhs_array_fourier)
    if (allocated(rhs_array)) deallocate(rhs_array)
    if (allocated(dgold_dx_normal)) deallocate(dgold_dx_normal)
    if (allocated(dgold_dy_normal)) deallocate(dgold_dy_normal)

    if(runtype_option_switch == runtype_multibox) call scope(allprocs)

    if(runtype_option_switch == runtype_multibox) call scope(subprocs)

    if (proc0) call time_message(.false.,time_gke(:,7),' ExB nonlinear advance, nisl')

    contains

    subroutine forward_transform (gk,gx)

      use stella_transforms, only: transform_ky2y, transform_kx2x
      use stella_transforms, only: transform_ky2y_xfirst, transform_kx2x_xfirst

      implicit none

      complex, dimension(:,:), intent (in) :: gk
      real, dimension(:,:), intent (out) :: gx

      ! we have i*ky*g(kx,ky) for ky >= 0 and all kx
      ! want to do 1D complex to complex transform in y
      ! which requires i*ky*g(kx,ky) for all ky and kx >= 0
      ! use g(kx,-ky) = conjg(g(-kx,ky))
      ! so i*(-ky)*g(kx,-ky) = -i*ky*conjg(g(-kx,ky)) = conjg(i*ky*g(-kx,ky))
      ! and i*kx*g(kx,-ky) = i*kx*conjg(g(-kx,ky)) = conjg(i*(-kx)*g(-kx,ky))
      ! and i*(-ky)*J0(kx,-ky)*phi(kx,-ky) = conjg(i*ky*J0(-kx,ky)*phi(-kx,ky))
      ! and i*kx*J0(kx,-ky)*phi(kx,-ky) = conjg(i*(-kx)*J0(-kx,ky)*phi(-kx,ky))
      ! i.e., can calculate dg/dx, dg/dy, d<phi>/dx and d<phi>/dy
      ! on stella (kx,ky) grid, then conjugate and flip sign of (kx,ky)
      ! NB: J0(kx,ky) = J0(-kx,-ky)
      ! TODO DSO: coordinate change for shearing
      call swap_kxky (gk, g0k_swap)
      call transform_ky2y (g0k_swap, g0kxy)
      call transform_kx2x (g0kxy, gx)

    end subroutine forward_transform

    !> Perform the Fourier transform, but with extra padding with zeroes so that
    !> the size of the untransformed and transformed arrays are twice as large.
    !> This gets us the values in real space on a grid which is twice as fine,
    !> with spectral accuracy.
    subroutine forward_transform_extra_padding (gk,gx_extra_padding)

      use stella_transforms, only: transform_ky2y_extra_padding, transform_kx2x_extra_padding
      !use stella_transforms, only: transform_ky2y_xfirst, transform_kx2x_xfirst

      implicit none

      complex, dimension(:,:), intent (in) :: gk
      real, dimension(:,:), intent (out) :: gx_extra_padding

      ! we have i*ky*g(kx,ky) for ky >= 0 and all kx
      ! want to do 1D complex to complex transform in y
      ! which requires i*ky*g(kx,ky) for all ky and kx >= 0
      ! use g(kx,-ky) = conjg(g(-kx,ky))
      ! so i*(-ky)*g(kx,-ky) = -i*ky*conjg(g(-kx,ky)) = conjg(i*ky*g(-kx,ky))
      ! and i*kx*g(kx,-ky) = i*kx*conjg(g(-kx,ky)) = conjg(i*(-kx)*g(-kx,ky))
      ! and i*(-ky)*J0(kx,-ky)*phi(kx,-ky) = conjg(i*ky*J0(-kx,ky)*phi(-kx,ky))
      ! and i*kx*J0(kx,-ky)*phi(kx,-ky) = conjg(i*(-kx)*J0(-kx,ky)*phi(-kx,ky))
      ! i.e., can calculate dg/dx, dg/dy, d<phi>/dx and d<phi>/dy
      ! on stella (kx,ky) grid, then conjugate and flip sign of (kx,ky)
      ! NB: J0(kx,ky) = J0(-kx,-ky)
      ! TODO DSO: coordinate change for shearing
      call swap_kxky (gk, g0k_swap)
      call transform_ky2y_extra_padding (g0k_swap, g0kxy_extra_padding)
      call transform_kx2x_extra_padding (g0kxy_extra_padding, gx_extra_padding)

    end subroutine forward_transform_extra_padding

    subroutine get_approx_departure_point (vy, vx, yidx_in, xidx_in, y_departure, x_departure, single_step)

      use kt_grids, only: x, y
      use kt_grids, only: dx, dy
      use stella_time, only: code_dt

      implicit none

      real, dimension(:,:), intent (in) :: vy, vx
      integer, intent (in) :: yidx_in, xidx_in
      logical, intent(in) :: single_step
      !real, intent (in) :: y_in, x_in
      real, intent (out) :: y_departure, x_departure

      integer :: max_iterations, counter, xidx, yidx
      real :: u_x, u_y, u_y_new, u_x_new, time_remaining, xboundary, yboundary, xnorm, ynorm
      real :: min_vy_magnitude, ydist_dt, min_vx_magnitude, xdist_dt, very_small_dt

      y_departure = y(yidx_in)
      x_departure = x(xidx_in)

      max_iterations = 100
      very_small_dt = 1E-12

      xidx = xidx_in
      yidx = yidx_in
      u_x = vx(yidx,xidx)
      u_y = vy(yidx,xidx)
      if (single_step) then
        time_remaining = code_dt
      else
        time_remaining = 2*code_dt
      end if

      counter=0
      do while ((time_remaining > 0) .and. (counter < max_iterations))
        counter = counter + 1
        !!! Take a step in (x,y) based on the velocity at [x, y, t^n].

        !!! Find the boundaries in x and y (at which the cell velocity changes)
        !! Old - velocities upsampled, which means the cells are dx/2 * dy/2.
        !! Cell boundaries at -dx/4, dx/4, 3dx/4, . . .
        ! xnorm = (2*x_departure/dx - 0.5); ynorm = (2*y_departure/dy - 0.5)

        !! New - velocities not upsampled, so cells are dx*dy. Cell boundaries
        !! at -dx/2, dx/2, 3dx/2, . . .
        ! "Normalised" x value & y value. This is an integer when we're
        ! sitting on a cell boundary.
        xnorm = (x_departure/dx - 0.5); ynorm = (y_departure/dy - 0.5)

        ! Our boundaries are nonperiodic.
        if (u_x > 0) then
          ! floor/ceiling gets us the normalised x value of the nearest
          ! cell boundary. We then convert back to "ordinary" x value.
          xboundary = ((floor(xnorm) + 0.5) * dx)
        else
          xboundary = ((ceiling(xnorm)  + 0.5) * dx)
        end if

        if (u_y > 0) then
          yboundary = ((floor(ynorm) + 0.5) * dy)
        else
          yboundary = ((ceiling(ynorm)  + 0.5) * dy)
        end if

        ! Calculate the time required to reach the nearest boundaries in y and x
        ! Careful though - if our velocities are too small we'll get <>dist_dt=inf
        min_vx_magnitude = (dx/code_dt)*1e-10   ! We know it's very unlikely that the trajectory with this velocity will make it into the next cell
        min_vy_magnitude = (dy/code_dt)*1e-10   ! We know it's very unlikely that the trajectory with this velocity will make it into the next cell

        if (abs(u_y) < min_vy_magnitude) then
          ydist_dt = 10*code_dt    ! This ensures ydist_dt > time_remaining
        else
          ydist_dt = -(yboundary - y_departure)/u_y
        end if

        if (abs(u_x) < min_vx_magnitude) then
          xdist_dt = 10*code_dt ! This ensures xdist_dt > time_remaining
        else
          xdist_dt = -(xboundary - x_departure)/u_x
        end if

        if ((ydist_dt > time_remaining) .and. (xdist_dt > time_remaining)) then
          ! Stop before we hit the next boundary
          ! Update location
          y_departure = (y_departure - u_y * time_remaining)
          x_departure = (x_departure - u_x * time_remaining)
          time_remaining = 0
        else
          ! Hit the next boundary before we run out of time.
          if (ydist_dt < xdist_dt) then
            ! We've hit the boundary in y - update the non-periodic
            ! position of the particle and the periodic position of the
            ! particle.
            y_departure = (yboundary - u_y*very_small_dt )   ! slightly overstep, so we're definitely in the next cell
            x_departure = (x_departure - u_x * (ydist_dt + very_small_dt))

            ! Update the values of u_x, u_y
            if (u_y > 0) then
              ! +ve u_y, so going back in time we're going in the -ve y direction; our new cell is
              ! below the old cell
              yidx = yidx - 1
            else
              ! -ve u_y, so going back in time we're going in the +ve y direction; our new cell is
              ! below the old cell
              yidx = yidx + 1
            end if

            ! The yidx range is (1, 2, 3, . . . ny), so anything
            ! beyond this range must be mapped into this range using periodicity
            ! modulo operation maps into range(0, ny-1), so need to correct for this
            yidx = modulo((yidx - 1),(ny))+1

            ! Update the velocities
            u_x = vx(yidx, xidx)
            ! Update u_y, but if the sign is different, we're going to
            ! bounce back and forth (this indicates the velocity falling
            ! to zero somewhere between the 2 cell centres). To avoid the "bouncing",
            ! set velocity to zero.
            u_y_new = vy(yidx, xidx)
            if ((u_y * u_y_new) < 0) then
              ! Opposite signs, so set u_y=0
              u_y=0
            else
              u_y = u_y_new
            end if

            ! Update time_remaining. Include the "very small dt" contribution
            time_remaining = time_remaining - (ydist_dt + very_small_dt)
          else
            ! Hit the boundary in x
            x_departure = (xboundary - u_x*very_small_dt)    ! slightly overstep, so we're definitely in the next cell
            y_departure = (y_departure - u_y*(xdist_dt + very_small_dt))

            ! Update the values of u_x, u_y
            if (u_x > 0) then
              ! +ve u_x, so going back in time we're going in the -ve x direction; our new cell is
              ! to the left of the old cell
              xidx = xidx - 1
            else
              ! -ve u_y, so going back in time we're going in the +ve y direction; our new cell is
              ! below the old cell
              xidx = xidx + 1
            end if

            ! The xidx range is (1, 2, 3, . . . nx), so anything
            ! beyond this range must be mapped into this range using periodicity
            ! modulo operation maps into range(0, nx-1), so need to correct for this
            xidx = modulo((xidx - 1),(nx))+1

            ! Update velocities
            u_y = vy(yidx, xidx)
            ! Update u_x, but if the sign is different, we're going to
            ! bounce back and forth (this indicates the velocity falling
            ! to zero somewhere between the 2 cell centres). To avoid the "bouncing",
            ! set velocity to zero.
            u_x_new = vx(yidx, xidx)
            if ((u_x * u_x_new) < 0) then
              ! Opposite signs, so set u_x=0 (ux passes through zero)
              u_x=0
            else
                u_x = u_x_new
            end if
            ! Update time_remaining. Include the "very small dt" contribution
            time_remaining = time_remaining - (xdist_dt + very_small_dt)
          end if
        end if
      end do

    end subroutine get_approx_departure_point

  end subroutine advance_ExB_nonlinearity_nisl

  subroutine advance_parallel_nonlinearity (g, gout, restart_time_step)

    use constants, only: zi
    use mp, only: proc0, min_allreduce, mp_abort
    use mp, only: scope, allprocs,subprocs
    use stella_layouts, only: vmu_lo, xyz_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use job_manage, only: time_message
    use finite_differences, only: second_order_centered_zed
    use finite_differences, only: third_order_upwind
    use redistribute, only: gather, scatter
    use fields_arrays, only: phi, phi_corr_QN,  phi_corr_GA
    use stella_transforms, only: transform_ky2y, transform_y2ky
    use stella_transforms, only: transform_kx2x, transform_x2kx
    use stella_time, only: cfl_dt, code_dt, code_dt_max
    use run_parameters, only: cfl_cushion, delt_adjust
    use zgrid, only: nzgrid, delzed, ntubes
    use extended_zgrid, only: neigen, nsegments, ikxmod
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: periodic
    use physics_flags, only: full_flux_surface, radial_variation
    use kt_grids, only: akx, aky, nakx, naky, nx, ny, ikx_max
    use kt_grids, only: swap_kxky, swap_kxky_back, rho_clamped
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: dvpa, vpa, mu
    use gyro_averages, only: gyro_average
    use parallel_streaming, only: stream_sign
    use dist_redistribute, only: xyz2vmu
    use file_utils, only: runtype_option_switch, runtype_multibox
    use extended_zgrid, only: fill_zed_ghost_zones

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout
    logical, intent (out) :: restart_time_step

    integer :: ivmu, ixyz
    integer :: iz, it, iv, imu, is
    integer :: iky, ie, iseg
    integer :: advect_sign
    real, dimension (:), allocatable :: dgdv
    real, dimension (:,:,:,:,:), allocatable :: g0xy
    real, dimension (:,:,:), allocatable :: gxy_vmulocal
    real, dimension (:,:), allocatable :: g1xy, advect_speed
    complex, dimension (2) :: gleft, gright
    complex, dimension (:,:,:,:), allocatable :: phi_gyro, dphidz
    complex, dimension (:,:), allocatable :: g0k, g0kxy, g0k_swap
    complex, dimension (:,:), allocatable :: tmp

    ! alpha-component of magnetic drift (requires ky -> y)
    if (proc0) call time_message(.false.,time_parallel_nl(:,1),' parallel nonlinearity advance')

    restart_time_step = .false.

    ! overview:
    ! need g and d<phi>/dz in (x,y) space in
    ! order to upwind dg/dvpa
    ! 1) transform d<phi>/dz from (kx,ky) to (x,y). layout: vmu_lo
    ! 2) need sign of parnl advection in xyz_lo (since dg/dvpa
    !    requires vpa local), so d<phi>/dz(vmu_lo) --> d<phi>/dz(xyz_lo)
    ! 3) transform g from (kx,ky) to (x,y). layout: vmu_lo
    ! 4) dg/dvpa requires vpa local, so g(vmu_lo) --> g(xyz_lo)
    ! 5) calculate dg/dvpa
    ! 6) multiply dg/dvpa with d<phi>/dz
    ! 7) product(xyz_lo) --> product(vmu_lo)
    ! 8) inverse transform product(vmu_lo)

    allocate (g0k(naky,nakx))
    allocate (g0xy(ny,nx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate (g0kxy(ny,ikx_max))
    if(radial_variation) allocate (g1xy(ny,nx))
    allocate (phi_gyro(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (dphidz(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (g0k_swap(2*naky-1,ikx_max))
    allocate (tmp(size(gout,1),size(gout,2)))

    ! get d<phi>/dz in vmu_lo
    ! we will need to transform it to real-space
    ! as its sign is needed for upwinding of dg/dvpa
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)

       ! construct <phi>
       dphidz = phi
       if(radial_variation) dphidz = dphidz + phi_corr_QN
       call gyro_average (dphidz, ivmu, phi_gyro)
       if(radial_variation) phi_gyro = phi_gyro + phi_corr_GA(:,:,:,:,ivmu)

       do iky = 1, naky
          do it = 1, ntubes
             do ie = 1, neigen(iky)
                do iseg = 1, nsegments(ie,iky)
                   ! first fill in ghost zones at boundaries in g(z)
                   call fill_zed_ghost_zones (it, iseg, ie, iky, phi_gyro, gleft, gright)
                   ! now get d<phi>/dz
                   call second_order_centered_zed (iz_low(iseg), iseg, nsegments(ie,iky), &
                        phi_gyro(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it), &
                        delzed(0), stream_sign(iv), gleft, gright, periodic(iky), &
                        dphidz(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it))
                end do
             end do
          end do
       end do

       if(radial_variation) then
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               ! use reality to swap from ky >= 0, all kx to kx >= 0 , all ky
               call swap_kxky (dphidz(:,:,iz,it), g0k_swap)
               ! transform in y
               call transform_ky2y (g0k_swap, g0kxy)
               ! transform in x
               call transform_kx2x (g0kxy, g1xy)
               g0xy(:,:,iz,it,ivmu) = g1xy*(par_nl_fac(iz,is) + d_par_nl_fac_dr(iz,is)*spread(rho_clamped,1,ny))

               g0k = zi*spread(aky,2,nakx)*phi_gyro(:,:,iz,it)
               call swap_kxky(g0k,g0k_swap)
               call transform_ky2y (g0k_swap, g0kxy)
               call transform_kx2x (g0kxy, g1xy)
               g0xy(:,:,iz,it,ivmu) = g0xy(:,:,iz,it,ivmu)  &
                                    + vpa(iv)*g1xy*(par_nl_drifty(iz) + d_par_nl_drifty_dr(iz)*spread(rho_clamped,1,ny))

               g0k = zi*spread(akx,1,naky)*phi_gyro(:,:,iz,it)
               call swap_kxky(g0k,g0k_swap)
               call transform_ky2y (g0k_swap, g0kxy)
               call transform_kx2x (g0kxy, g1xy)
               g0xy(:,:,iz,it,ivmu) = g0xy(:,:,iz,it,ivmu)  &
                                    + vpa(iv)*g1xy*(par_nl_driftx(iz) + d_par_nl_driftx_dr(iz)*spread(rho_clamped,1,ny))

               g0xy(:,:,iz,it,ivmu) = g0xy(:,:,iz,it,ivmu)  &
                                    + vpa(iv)*mu(imu)*(par_nl_curv(iz,is) + d_par_nl_curv_dr(iz,is)*spread(rho_clamped,1,ny))

            end do
         end do
       else
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0k = dphidz(:,:,iz,it)*par_nl_fac(iz,is) + vpa(iv)*mu(imu)*par_nl_curv(iz,is) &
                   + zi*vpa(iv)*phi_gyro(:,:,iz,it)*(  spread(akx,1,naky)*par_nl_driftx(iz)   &
                                                     + spread(aky,2,nakx)*par_nl_drifty(iz))
               ! use reality to swap from ky >= 0, all kx to kx >= 0 , all ky
               call swap_kxky (g0k, g0k_swap)
               ! transform in y
               call transform_ky2y (g0k_swap, g0kxy)
               ! transform in x
               call transform_kx2x (g0kxy, g0xy(:,:,iz,it,ivmu))
            end do
         end do
       endif
    end do

    ! do not need phi_gyro or dphidz  again so deallocate
    deallocate (phi_gyro, dphidz)
    deallocate (g0k)
    if (allocated(g1xy)) deallocate (g1xy)

    allocate (gxy_vmulocal(nvpa,nmu,xyz_lo%llim_proc:xyz_lo%ulim_alloc))
    allocate (advect_speed(nmu,xyz_lo%llim_proc:xyz_lo%ulim_alloc))

    ! we now have the advection velocity in vpa in (x,y) space
    ! next redistribute it so that (vpa,mu) are local
    if (proc0) call time_message(.false.,time_parallel_nl(:,2),' parallel nonlinearity redist')
    call scatter (xyz2vmu, g0xy, gxy_vmulocal)
    if (proc0) call time_message(.false.,time_parallel_nl(:,2),' parallel nonlinearity redist')
    ! advect_speed does not depend on vpa
    advect_speed = gxy_vmulocal(1,:,:)


    ! transform g from (kx,ky) to (x,y)
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do it = 1, ntubes
          do iz = -nzgrid, nzgrid
             call swap_kxky (g(:,:,iz,it,ivmu), g0k_swap)
             ! transform in y
             call transform_ky2y (g0k_swap, g0kxy)
             ! transform in x
             call transform_kx2x (g0kxy, g0xy(:,:,iz,it,ivmu))
          end do
       end do
    end do

    ! redistribute so that (vpa,mu) local
    if (proc0) call time_message(.false.,time_parallel_nl(:,2),' parallel nonlinearity redist')
    call scatter (xyz2vmu, g0xy, gxy_vmulocal)
    if (proc0) call time_message(.false.,time_parallel_nl(:,2),' parallel nonlinearity redist')

    allocate (dgdv(nvpa))

    ! we now need to form dg/dvpa and obtain product of dg/dvpa with advection speed
    do ixyz = xyz_lo%llim_proc, xyz_lo%ulim_proc
       do imu = 1, nmu
          ! advect_sign set to +/- 1 depending on sign of the parallel nonlinearity
          ! advection velocity
          ! NB: advect_sign = -1 corresponds to positive advection velocity
          advect_sign = int(sign(1.0,advect_speed(imu,ixyz)))
          call third_order_upwind (1,gxy_vmulocal(:,imu,ixyz),dvpa,advect_sign,dgdv)
          gxy_vmulocal(:,imu,ixyz) = dgdv*advect_speed(imu,ixyz)
          cfl_dt = min(cfl_dt,dvpa/abs(advect_speed(imu,ixyz)))
       end do
    end do

    ! finished with dgdv and advect_speed
    deallocate (dgdv, advect_speed)

    ! now that we have the full parallel nonlinearity in (x,y)-space
    ! need to redistribute so that (x,y) local for transforms
    if (proc0) call time_message(.false.,time_parallel_nl(:,2),' parallel nonlinearity redist')
    call gather (xyz2vmu, gxy_vmulocal, g0xy)
    if (proc0) call time_message(.false.,time_parallel_nl(:,2),' parallel nonlinearity redist')

    ! finished with gxy_vmulocal
    deallocate (gxy_vmulocal)

    ! g0xy is parallel nonlinearity term with (x,y) on processor
    ! need to inverse Fourier transform
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do it = 1, ntubes
          do iz = -nzgrid, nzgrid
             call transform_x2kx (g0xy(:,:,iz,it,ivmu), g0kxy)
             if (full_flux_surface) then
                gout(:,:,iz,it,ivmu) = gout(:,:,iz,it,ivmu) + code_dt*g0kxy
             else
                call transform_y2ky (g0kxy, g0k_swap)
                call swap_kxky_back (g0k_swap, tmp)
                gout(:,:,iz,it,ivmu) = gout(:,:,iz,it,ivmu) + code_dt*tmp
             end if
          end do
       end do
    end do
    deallocate (g0k_swap, g0kxy, g0xy)

    if(runtype_option_switch == runtype_multibox) call scope(allprocs)

    call min_allreduce (cfl_dt)

    if(runtype_option_switch == runtype_multibox) call scope(subprocs)

    if (code_dt > cfl_dt*cfl_cushion) then
       if (proc0) then
          write (*,*) ' '
          write (*,*) 'CHANGING TIME STEP:'
          write (*,'(A16, ES10.2E2)') "   code_dt:"//REPEAT(' ',50),code_dt
          write (*,'(A16, ES10.2E2)') "   cfl_dt:"//REPEAT(' ',50),cfl_dt
          write (*,'(A16, ES10.2E2)') "   cfl_cushion:"//REPEAT(' ',50),cfl_cushion
          write (*,'(A16, ES10.2E2)') "   delt_adjust:"//REPEAT(' ',50),delt_adjust
          write (*,'(A65)') '     ==> User-specified delt is larger than cfl_dt*cfl_cushion.'//REPEAT(' ',50)
          write (*,'(A61,ES12.4)') '     ==> Changing code_dt to cfl_dt*cfl_cushion/delt_adjust ='//REPEAT(' ',50), cfl_dt*cfl_cushion/delt_adjust
          write(*,*) ' '
       end if
       code_dt = cfl_dt*cfl_cushion/delt_adjust
       call reset_dt
       restart_time_step = .true.
    else if (code_dt < min(cfl_dt*cfl_cushion/delt_adjust,code_dt_max)) then
       if (proc0) then
          write (*,*) ' '
          write (*,*) 'CHANGING TIME STEP:'
          write (*,'(A16, ES10.2E2)') "   code_dt:"//REPEAT(' ',50),code_dt
          write (*,'(A16, ES10.2E2)') "   cfl_dt:"//REPEAT(' ',50),cfl_dt
          write (*,'(A16, ES10.2E2)') "   cfl_cushion:"//REPEAT(' ',50),cfl_cushion
          write (*,'(A16, ES10.2E2)') "   delt_adjust:"//REPEAT(' ',50),delt_adjust
          write (*,'(A65)') '     ==> User-specified delt is larger than cfl_dt*cfl_cushion.'//REPEAT(' ',50)
          write (*,'(A61,ES12.4)') '     ==> Changing code_dt to cfl_dt*cfl_cushion/delt_adjust ='//REPEAT(' ',50), cfl_dt*cfl_cushion/delt_adjust
          write(*,*) ' '
       end if
       code_dt = min(cfl_dt*cfl_cushion/delt_adjust,code_dt_max)
       call reset_dt
       !    else
       !       gout = code_dt*gout
    end if

    if (proc0) call time_message(.false.,time_parallel_nl(:,1),' parallel nonlinearity advance')

  end subroutine advance_parallel_nonlinearity

  subroutine advance_radial_variation (g, gout)

    use mp, only: mp_abort, proc0
    use job_manage, only: time_message
    use fields, only: get_dchidy
    use fields_arrays, only: phi, apar
    use fields_arrays, only: phi_corr_QN,  phi_corr_GA
    !   use fields_arrays, only: apar_corr_QN, apar_corr_GA
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky, multiply_by_rho
    use gyro_averages, only: gyro_average, gyro_average_j1
    use run_parameters, only: fphi
    use physics_flags, only: full_flux_surface
    use physics_flags, only: include_parallel_streaming, include_mirror
    use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
    use dist_fn_arrays, only: wdriftpx_g, wdriftpy_g
    use dist_fn_arrays, only: wdriftpx_phi, wdriftpy_phi !, adiabatic_phi
    use dist_fn_arrays, only: wstar, wstarp
    use mirror_terms, only: add_mirror_radial_variation
    use flow_shear, only: prl_shear, prl_shear_p
    use parallel_streaming, only: add_parallel_streaming_radial_variation

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout

    integer :: ia, ivmu, iv, imu, is, iz, it

    complex, dimension (:,:), allocatable :: g0k, g1k, g0a
    complex, dimension (:,:,:,:,:), allocatable :: g_corr

    allocate (g0k(naky,nakx))

    allocate (g1k(naky,nakx))
    allocate (g0a(naky,nakx))


    if (debug) write (*,*) 'time_advance::solve_gke::advance_radial_variation'

    if (proc0) call time_message(.false.,time_gke(:,10),' radial variation advance')

    if(include_mirror .or. include_parallel_streaming) then
      allocate (g_corr(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g_corr = 0.
    endif

    !grab the mirror and parallel streaming corrections here to save on FFTs
    if (include_mirror) then
      call add_mirror_radial_variation(g,g_corr)
    endif
    if (include_parallel_streaming) then
      call add_parallel_streaming_radial_variation(g,g_corr,gout)
    endif


    if (full_flux_surface) then
       ! FLAG -- ADD SOMETHING HERE
       call mp_abort ('wstarp term not yet setup for full_flux_surface = .true. aborting.')
    endif

    ia = 1
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            g0k = 0.

            !wstar variation
            call get_dchidy (iz, ivmu, phi(:,:,iz,it), apar(:,:,iz,it), g0a)
            g0k = g0k + g0a*wstarp(ia,iz,ivmu)

            !radial variation in ExB nonlinearity is handled in advance_ExB_nonlinearity

            !wdrift(x/y) - g

            call get_dgdx(g(:,:,iz,it,ivmu),g0a)
            g0k = g0k + g0a*wdriftpx_g(ia,iz,ivmu)

            call get_dgdy(g(:,:,iz,it,ivmu),g0a)
            g0k = g0k + g0a*wdriftpy_g(ia,iz,ivmu)

            !wdrift - phi
            call get_dgdx(phi(:,:,iz,it),g1k)
            !wdriftx variation
            call gyro_average (g1k,iz,ivmu,g0a )
            g0k = g0k + g0a*wdriftpx_phi(ia,iz,ivmu)

            call get_dgdy(phi(:,:,iz,it),g1k)
            !wdrifty variation
            call gyro_average (g1k,iz,ivmu,g0a)
            g0k = g0k + g0a*wdriftpy_phi(ia,iz,ivmu)

            !prl_shear variation
            g0k = g0k + g0a*prl_shear_p(ia,iz,ivmu)

            !mirror term and/or parallel streaming
            if(include_mirror.or.include_parallel_streaming) then
              g0k = g0k + g_corr(:,:,iz,it,ivmu)
            endif

            !inverse and forward transforms
            call multiply_by_rho(g0k)


            !
            !quasineutrality/gyroaveraging
            !
            call gyro_average (phi_corr_QN(:,:,iz,it),iz,ivmu,g0a)
            g0a = fphi*(g0a + phi_corr_GA(:,:,iz,it,ivmu))

            !wstar - gyroaverage/quasineutrality variation
            call get_dgdy(g0a,g1k)
            g0k = g0k + g1k*wstar(ia,iz,ivmu)

            !wdrifty gyroaverage/quasineutrality variation
            g0k = g0k + g1k*wdrifty_phi(ia,iz,ivmu)

            !prl_shear gyroaverage/quasineutrality variation
            g0k = g0k + g1k*prl_shear(ia,iz,ivmu)

            !wdriftx gyroaverage/quasineutrality variation
            call get_dgdx(g0a,g1k)
            g0k = g0k + g1k*wdriftx_phi(ia,iz,ivmu)

            !           !wdriftx F_M/T_s variation
!           call gyro_average (phi(:,:,iz,it),iz,ivmu,g0a)
!           g0a = adiabatic_phi(ia,iz,ivmu)*g0a
!           call multiply_by_rho(g0a)
!           call get_dgdx(g0a,g1k)
!           g0k = g0k + g1k*wdriftx_phi(ia,iz,ivmu)

            gout(:,:,iz,it,ivmu) = gout(:,:,iz,it,ivmu) + g0k
          end do
       end do
    end do

    deallocate (g0k, g1k, g0a)
    if(allocated(g_corr)) deallocate(g_corr)

    if (proc0) call time_message(.false.,time_gke(:,10),' radial variation advance')

  end subroutine advance_radial_variation


  subroutine get_dgdy_2d (g, dgdy)

    use constants, only: zi
    use kt_grids, only: nakx, aky

    implicit none

    complex, dimension (:,:), intent (in) :: g
    complex, dimension (:,:), intent (out) :: dgdy

    dgdy = zi*spread(aky,2,nakx)*g

  end subroutine get_dgdy_2d

  subroutine get_dgdy_3d (g, dgdy)

    use constants, only: zi
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, aky

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: dgdy

    dgdy = zi*spread(spread(spread(aky,2,nakx),3,2*nzgrid+1),4,ntubes)*g

  end subroutine get_dgdy_3d

  subroutine get_dgdy_4d (g, dgdy)

    use constants, only: zi
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, aky

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (out) :: dgdy

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       dgdy(:,:,:,:,ivmu) = zi*spread(spread(spread(aky,2,nakx),3,2*nzgrid+1),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine get_dgdy_4d

  subroutine add_dg_term (g, wdrift_in, src)

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (-nzgrid:,vmu_lo%llim_proc:), intent (in) :: wdrift_in
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,:,ivmu) = src(:,:,:,:,ivmu) &
            + spread(spread(spread(wdrift_in(:,ivmu),1,naky),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_dg_term

  subroutine add_dg_term_global (g, wdrift_in, src)

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: wdrift_in
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,:,ivmu) = src(:,:,:,:,ivmu) - spread(spread(wdrift_in(:,:,ivmu),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_dg_term_global

  subroutine get_dgdx_2d (g, dgdx)

    use constants, only: zi
    use kt_grids, only: naky, akx

    implicit none

    complex, dimension (:,:), intent (in) :: g
    complex, dimension (:,:), intent (out) :: dgdx

    dgdx = zi*spread(akx,1,naky)*g

  end subroutine get_dgdx_2d

  subroutine get_dgdx_3d (g, dgdx)

    use constants, only: zi
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, akx

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: dgdx

    dgdx = zi*spread(spread(spread(akx,1,naky),3,2*nzgrid+1),4,ntubes)*g

  end subroutine get_dgdx_3d

  subroutine get_dgdx_4d (g, dgdx)

    use constants, only: zi
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, akx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (out) :: dgdx

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       dgdx(:,:,:,:,ivmu) = zi*spread(spread(spread(akx,1,naky),3,2*nzgrid+1),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine get_dgdx_4d

  subroutine add_dphi_term (g, wdrift, src)

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (-nzgrid:,vmu_lo%llim_proc:), intent (in) :: wdrift
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,:,ivmu) = src(:,:,:,:,ivmu) &
            + spread(spread(spread(wdrift(:,ivmu),1,naky),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_dphi_term

  subroutine add_dphi_term_global (g, wdrift, src)

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: wdrift
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,:,ivmu) = src(:,:,:,:,ivmu) &
            + spread(spread(wdrift(:,:,ivmu),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_dphi_term_global

  subroutine add_wstar_term (g, src)

    use dist_fn_arrays, only: wstar
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,:,ivmu) = src(:,:,:,:,ivmu) &
            + spread(spread(spread(wstar(1,:,ivmu),1,naky),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_wstar_term

  subroutine advance_implicit (phi, apar, reverse_implicit_order, g)
!  subroutine advance_implicit (phi, apar, g)

    use mp, only: proc0
    use job_manage, only: time_message
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use dissipation, only: hyper_dissipation, advance_hyper_dissipation
    use physics_flags, only: include_parallel_streaming, include_drifts
    use physics_flags, only: radial_variation, full_flux_surface
    use physics_flags, only: include_mirror, prp_shear_enabled
    use run_parameters, only: stream_implicit, mirror_implicit, drifts_implicit
    use parallel_streaming, only: advance_parallel_streaming_implicit
    use fields, only: advance_fields, fields_updated
    use mirror_terms, only: advance_mirror_implicit
    use dissipation, only: collisions_implicit, include_collisions
    use dissipation, only: advance_collisions_implicit
    use run_parameters, only: driftkinetic_implicit
    use flow_shear, only: advance_perp_flow_shear
    use multibox, only: RK_step

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar
    logical, intent(in) :: reverse_implicit_order
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
!    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out), target :: g

!    complex, dimension (:,:,:,:,:), pointer :: gk, gy
!    complex, dimension (:,:,:,:,:), allocatable, target :: g_dual

!    ! the 'g' that enters this subroutine may be in alpha-space or kalpha-space
!    ! figure out which it is
!    if (size(g,1) == naky) then
!       alpha_space = .false.
!       gk => g
!       if (full_flux_surface) then
!          allocate (g_dual(nalpha,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
!          gy => g_dual
!       end if
!    else
!       alpha_space = .true.
!       allocate (g_dual(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
!       gy => g
!       gk => g_dual
!    end if
    ! start the timer for the implicit part of the solve
    if (proc0) call time_message(.false.,time_gke(:,9),' implicit')

    ! reverse the order of operations every time step
    ! as part of alternating direction operator splitting
    ! this is needed to ensure 2nd order accuracy in time
!    if (mod(istep,2)==0) then
       ! g^{*} (coming from explicit solve) is input
       ! get g^{**}, with g^{**}-g^{*} due to mirror term

    if(RK_step) call mb_communicate (g)

    ! Apply the operators sequentially, order depending on the reverse_implicit_order
    ! logical
    if (.not. reverse_implicit_order) then

       if (prp_shear_enabled) then
          call advance_perp_flow_shear(g)
          fields_updated = .false.
       end if
       if (hyper_dissipation) then
!          ! for hyper-dissipation, need to be in k-alpha space
!          if (alpha_space) call transform_y2ky (gy, gk)
          call advance_hyper_dissipation (g)
          fields_updated = .false.
       end if
       if (collisions_implicit .and. include_collisions) then
          call advance_fields (g, phi, apar, dist='gbar')
          call advance_collisions_implicit (mirror_implicit, phi, apar, g)
          fields_updated = .false.
       end if
       if (mirror_implicit .and. include_mirror) then
!          if (full_flux_surface) then
!             allocate (gy(ny,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
!             if (.not.alpha_space) call transform_ky2y (g, gy)
!          else
!             g_mirror => g
!          end if
          call advance_mirror_implicit (collisions_implicit, g)
          fields_updated = .false.
       end if
       ! get updated fields corresponding to advanced g
       ! note that hyper-dissipation and mirror advances
       ! depended only on g and so did not need field update
       call advance_fields (g, phi, apar, dist='gbar')

       ! g^{**} is input
       ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
       if ((stream_implicit.or.driftkinetic_implicit) .and. include_parallel_streaming) then
            call advance_parallel_streaming_implicit (g, phi, apar)
            if(radial_variation.or.full_flux_surface) fields_updated = .false.
       endif
       call advance_fields (g, phi, apar, dist='gbar')
       if (drifts_implicit .and. include_drifts) call advance_drifts_implicit (g, phi, apar)
    else
       ! get updated fields corresponding to advanced g
       ! note that hyper-dissipation and mirror advances
       ! depended only on g and so did not need field update
       call advance_fields (g, phi, apar, dist='gbar')
       if (drifts_implicit) call advance_drifts_implicit (g, phi, apar)

       ! g^{**} is input
       ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
       if ((stream_implicit.or.driftkinetic_implicit) .and. include_parallel_streaming) then
          call advance_parallel_streaming_implicit (g, phi, apar)
          if(radial_variation.or.full_flux_surface) fields_updated = .false.
       endif

       if (mirror_implicit .and. include_mirror) then
          call advance_mirror_implicit (collisions_implicit, g)
          fields_updated = .false.
       end if

       if (collisions_implicit .and. include_collisions) then
          call advance_fields (g, phi, apar, dist='gbar')
          call advance_collisions_implicit (mirror_implicit, phi, apar, g)
          fields_updated = .false.
       end if

       if (hyper_dissipation) then
          call advance_hyper_dissipation (g)
          fields_updated = .false.
       end if

       if (prp_shear_enabled) then
          call advance_perp_flow_shear(g)
          fields_updated = .false.
       end if

    end if

!       call checksum (phi, phitot)
!       call checksum (g, gtot)
!       if (proc0) write (*,*) 'stream', phitot, gtot

    ! stop the timer for the implict part of the solve
    if (proc0) call time_message(.false.,time_gke(:,9),' implicit')
  end subroutine advance_implicit

  subroutine advance_drifts_implicit (g, phi, apar)

    use constants, only: zi
    use stella_layouts, only: vmu_lo
    use stella_geometry, only: dl_over_b
    use run_parameters, only: fphi, fapar, time_upwind
    use dist_fn_arrays, only: g1
    use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
    use dist_fn_arrays, only: wdriftx_g, wdrifty_g
    use dist_fn_arrays, only: wstar
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg
    use gyro_averages, only: aj0x, gyro_average
    use kt_grids, only: akx, aky, nakx, naky, zonal_mode
    use zgrid, only: nzgrid, ntubes
    use species, only: spec, has_electron_species
    use fields, only: advance_fields
    use vpamu_grids, only: integrate_species

    implicit none

    integer :: ivmu, iz, it, ia, ikx
    complex :: tmp

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar

    complex, dimension (:,:), allocatable :: wd_g, wd_phi, wstr
    complex, dimension (:,:,:), allocatable :: gyro_g

    ia = 1
    allocate (wd_g(naky,nakx))
    allocate (wd_phi(naky,nakx))
    allocate (wstr(naky,nakx))

    ! given g^{*}, obtain phi^{*} and apar^{*}
    call advance_fields (g, phi, apar, dist='gbar')

    ! solve for g^inh
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      do it = 1, ntubes
        do iz = -nzgrid, nzgrid
          wd_g   = -zi*(spread(akx,1,naky)*wdriftx_g(ia,iz,ivmu) &
                      + spread(aky,2,nakx)*wdrifty_g(ia,iz,ivmu))

          wd_phi = -zi*(spread(akx,1,naky)*wdriftx_phi(ia,iz,ivmu) &
                      + spread(aky,2,nakx)*wdrifty_phi(ia,iz,ivmu))

          wstr   = -zi*spread(aky,2,nakx)*wstar(ia,iz,ivmu)

          g1(:,:,iz,it,ivmu) = (g(:,:,iz,it,ivmu)*(1.0 - 0.5*(1.0 - time_upwind)*wd_g) &
                              - 0.5*(1.0 - time_upwind)*(wd_phi + wstr) &
                                *aj0x(:,:,iz,ivmu)*fphi*phi(:,:,iz,ia)) &
                                /(1.0 + 0.5*(1.0 + time_upwind)*wd_g)
        enddo
      enddo
    enddo

    !we have g_inh, now get phi
    if (fphi > epsilon(0.0)) then
      allocate (gyro_g(naky,nakx,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      do it = 1, ntubes
        do iz = -nzgrid, nzgrid
          do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            call gyro_average (g1(:,:,iz,it,ivmu), iz, ivmu, gyro_g(:,:,ivmu))
          enddo
          call integrate_species (gyro_g, iz, spec%z*spec%dens_psi0, phi(:,:,iz,it))
        enddo
        phi(:,:,:,it) = phi(:,:,:,it) / gamtot_drifts
        if(any(real(gamtot_drifts(1,1,:)).lt.epsilon(0.))) phi(1,1,:,it) = 0.0

        if (.not.has_electron_species(spec)) then
          ! no need to do anything extra for ky /= 0 because
          ! already accounted for in gamtot_h
          if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            if (zonal_mode(1)) then
              do ikx = 1, nakx
                tmp = sum(dl_over_b(ia,:)*phi(1,ikx,:,it))
                phi(1,ikx,:,it) = phi(1,ikx,:,it) + tmp*gamtot3_drifts(ikx,:)
              enddo
            endif
          endif
        endif
      enddo
      deallocate (gyro_g)
    endif

    !finally, get g
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      do it = 1, ntubes
        do iz = -nzgrid, nzgrid
          !these terms already contain a factor of code_dt and a
          ! negative sign
          wd_g   = -zi*(spread(akx,1,naky)*wdriftx_g(ia,iz,ivmu) &
                      + spread(aky,2,nakx)*wdrifty_g(ia,iz,ivmu))

          wd_phi = -zi*(spread(akx,1,naky)*wdriftx_phi(ia,iz,ivmu) &
                      + spread(aky,2,nakx)*wdrifty_phi(ia,iz,ivmu))

          wstr   = -zi*spread(aky,2,nakx)*wstar(ia,iz,ivmu)

          g(:,:,iz,it,ivmu) = g1(:,:,iz,it,ivmu) &
                              - 0.5*(1.0 + time_upwind)*(wd_phi + wstr) &
                              *aj0x(:,:,iz,ivmu)*fphi*phi(:,:,iz,it) &
                             /(1.0 + 0.5*(1.0 + time_upwind)*wd_g)
        enddo
      enddo
    enddo

    deallocate (wd_g, wd_phi, wstr)

  end subroutine advance_drifts_implicit

  subroutine mb_communicate (g_in)

    use mp, only: job
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use multibox, only: multibox_communicate
    use fields, only: fields_updated,advance_fields
    use fields_arrays, only: phi, apar
    use file_utils, only: runtype_option_switch, runtype_multibox

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g_in

    if(runtype_option_switch.ne.runtype_multibox) return

    if(job.ne.1) then
      call advance_fields(g_in,phi,apar,dist='gbar')
    endif

    call multibox_communicate(g_in)

    if(job.eq.1) then
      fields_updated = .false.
      call advance_fields(g_in,phi,apar,dist='gbar')
    endif

  end subroutine mb_communicate

  subroutine checksum_field (field, total)

    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky
    use extended_zgrid, only: neigen, nsegments, ikxmod
    use extended_zgrid, only: iz_low, iz_up

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: field
    real, intent (out) :: total

    integer :: it, iky, ie, iseg
    integer :: ikx

    total = 0.

    do iky = 1, naky
       do it = 1, ntubes
          do ie = 1, neigen(iky)
             iseg = 1
             ikx = ikxmod(iseg,ie,iky)
             total = total + sum(cabs(field(iky,ikx,iz_low(iseg):iz_up(iseg),it)))
             if (nsegments(ie,iky) > 1) then
                do iseg = 2, nsegments(ie,iky)
                   ikx = ikxmod(iseg,ie,iky)
                   total = total + sum(cabs(field(iky,ikx,iz_low(iseg)+1:iz_up(iseg),it)))
                end do
             end if
          end do
       end do
    end do

  end subroutine checksum_field

  subroutine checksum_dist (dist, total)

    use mp, only: sum_allreduce
    use zgrid, only: nzgrid
    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: dist
    real, intent (out) :: total

    integer :: ivmu
    real :: subtotal

    total = 0.

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       call checksum (dist(:,:,:,:,ivmu), subtotal)
       total = total + subtotal
    end do

    call sum_allreduce (total)

  end subroutine checksum_dist

  subroutine finish_time_advance

    use stella_transforms, only: finish_transforms
    use physics_flags, only: full_flux_surface
    use extended_zgrid, only: finish_extended_zgrid
    use parallel_streaming, only: finish_parallel_streaming
    use mirror_terms, only: finish_mirror
    use flow_shear, only : finish_flow_shear
    use neoclassical_terms, only: finish_neoclassical_terms
    use dissipation, only: finish_dissipation

    implicit none

    if (full_flux_surface) call finish_transforms
    call finish_dissipation
    call finish_parallel_nonlinearity
    call finish_wstar
    call finish_wdrift
    call finish_drifts_implicit
    call finish_parallel_streaming
    call finish_flow_shear
    call finish_mirror
    call finish_neoclassical_terms
    call deallocate_arrays

    time_advance_initialized = .false.
    readinit = .false.

  end subroutine finish_time_advance

  subroutine finish_parallel_nonlinearity

    implicit none

    if (allocated(par_nl_fac)) deallocate (par_nl_fac)
    if (allocated(par_nl_curv)) deallocate (par_nl_curv)
    if (allocated(par_nl_driftx)) deallocate (par_nl_driftx)
    if (allocated(par_nl_drifty)) deallocate (par_nl_drifty)

    parnlinit = .false.

  end subroutine finish_parallel_nonlinearity

  subroutine finish_wdrift

    use dist_fn_arrays, only: wdriftx_g, wdrifty_g
    use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
    use dist_fn_arrays, only: wdriftpx_g, wdriftpy_g
    use dist_fn_arrays, only: wdriftpx_phi, wdriftpy_phi
!   use dist_fn_arrays, only: adiabatic_phi

    implicit none

    if (allocated(wdriftx_g)) deallocate (wdriftx_g)
    if (allocated(wdrifty_g)) deallocate (wdrifty_g)
    if (allocated(wdriftx_phi)) deallocate (wdriftx_phi)
    if (allocated(wdrifty_phi)) deallocate (wdrifty_phi)
    if (allocated(wdriftpx_g)) deallocate (wdriftpx_g)
    if (allocated(wdriftpy_g)) deallocate (wdriftpy_g)
    if (allocated(wdriftpx_phi)) deallocate (wdriftpx_phi)
    if (allocated(wdriftpy_phi)) deallocate (wdriftpy_phi)
!   if (allocated(adiabatic_phi)) deallocate (adiabatic_phi)

    wdriftinit = .false.

  end subroutine finish_wdrift

  subroutine finish_wstar

    use dist_fn_arrays, only: wstar, wstarp

    implicit none

    if (allocated(wstar)) deallocate (wstar)
    if (allocated(wstarp)) deallocate (wstarp)

    wstarinit = .false.

  end subroutine finish_wstar

  subroutine finish_drifts_implicit

    implicit none

    if (allocated(gamtot_drifts)) deallocate (gamtot_drifts)
    if (allocated(gamtot3_drifts)) deallocate (gamtot3_drifts)

    driftimpinit = .false.

  end subroutine finish_drifts_implicit

  subroutine deallocate_arrays

    use dist_fn_arrays, only: g0, g1, g2, g3

    implicit none

    if (allocated(g0)) deallocate (g0)
    if (allocated(g1)) deallocate (g1)
    if (allocated(g2)) deallocate (g2)
    if (allocated(g3)) deallocate (g3)

  end subroutine deallocate_arrays

end module time_advance
