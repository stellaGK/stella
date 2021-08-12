module parallel_streaming

  implicit none

  public :: init_z_equation, finish_z_equation
  public :: stream, stream_c, stream_sign
  public :: advance_parallel_streaming_explicit
  public :: z_tridiagonal_solve
  public :: advance_z_implicit
  public :: parallel_streaming_initialized
  public :: add_parallel_streaming_radial_variation
  public :: stream_rad_var1
  public :: stream_rad_var2
  public :: time_parallel_streaming

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!! From original parallel_streaming
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! public :: init_parallel_streaming, finish_parallel_streaming
  ! public :: advance_parallel_streaming_explicit
  ! public :: advance_parallel_streaming_implicit
  ! public :: stream_tridiagonal_solve
  ! public :: parallel_streaming_initialized
  ! public :: stream, stream_c, stream_sign
  ! public :: stream_rad_var1
  ! public :: stream_rad_var2

  private

  interface center_zed
     module procedure center_zed_segment_real
     module procedure center_zed_extended
  end interface

  logical :: parallel_streaming_initialized = .false.

  integer, dimension (:), allocatable :: stream_sign
  real, dimension (:,:,:), allocatable :: stream, stream_c
  real, dimension (:,:), allocatable :: dgdt_tri_a, stream_tri_a, wdrift_tri_a
  real, dimension (:,:), allocatable :: dgdt_tri_b, stream_tri_b, wdrift_tri_b
  real, dimension (:,:), allocatable :: dgdt_tri_c, stream_tri_c, wdrift_tri_c
  real, dimension (:,:), allocatable :: gradpar_c
  real, dimension (:,:,:), allocatable :: wdriftx_g_c, wdrifty_g_c, wdriftx_phi_c, wdrifty_phi_c
  real, dimension (:,:,:), allocatable :: wstar_c
  real, dimension (:,:,:), allocatable :: stream_rad_var1
  real, dimension (:,:,:), allocatable :: stream_rad_var2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!! From original parallel_streaming
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! integer, dimension (:), allocatable :: stream_sign
  ! real, dimension (:,:,:), allocatable :: stream, stream_c
  ! real, dimension (:,:), allocatable :: stream_tri_a1, stream_tri_a2
  ! real, dimension (:,:), allocatable :: stream_tri_b1, stream_tri_b2
  ! real, dimension (:,:), allocatable :: stream_tri_c1, stream_tri_c2
  ! real, dimension (:,:), allocatable :: gradpar_c

  real, dimension (2) :: time_parallel_streaming

contains

  !> Initialisation of PDE-in-z time advance requires the following:
  !> 1) Calculate the entries of the bidiagonal matrix on the LHS of the GKE
  !> 2) Initialise stream_sign
  !> 3) Center the terms which appear in the PDE (accounts for z-upwinding)
  subroutine init_z_equation

    use mp, only: mp_abort
    use stella_time, only: code_dt
    use species, only: spec, nspec
    use vpamu_grids, only: nvpa
    use vpamu_grids, only: vpa
    use stella_layouts, only: vmu_lo, iv_idx
    use zgrid, only: nzgrid, nztot
    use kt_grids, only: nalpha
    use stella_geometry, only: gradpar
    use physics_flags, only: include_parallel_streaming, radial_variation
    use run_parameters, only: stream_implicit, driftkinetic_implicit
    use run_parameters, only: drifts_implicit_in_z
    use dist_fn_arrays, only: wdrifty_phi, wdriftx_phi, wdrifty_g, wdriftx_g, wstar

    implicit none

    integer :: iv, is, ivmu, ia

    ! RJD: taken from init_parallel_streaming.

    if (parallel_streaming_initialized) return
    parallel_streaming_initialized = .true.

    if (.not.allocated(stream)) allocate (stream(-nzgrid:nzgrid,nvpa,nspec)) ; stream = 0.
    if (.not.allocated(stream_sign)) allocate (stream_sign(nvpa)) ; stream_sign = 0

    ! sign of stream corresponds to appearing on RHS of GK equation
    ! i.e., this is the factor multiplying dg/dz on RHS of equation
    if (include_parallel_streaming) then
       stream = -code_dt*spread(spread(spec%stm_psi0,1,nztot),2,nvpa) &
            * spread(spread(vpa,1,nztot)*spread(gradpar,2,nvpa),3,nspec)
    else
       stream = 0.0
    end if

    ! stream_sign set to +/- 1 depending on the sign of the parallel streaming term.
    ! NB: stream_sign = -1 corresponds to positive advection velocity
    do iv = 1, nvpa
       stream_sign(iv) = int(sign(1.0,stream(0,iv,1)))
    end do

    ! If we're treating streaming implicitly, we need to center in z i.e.
    ! evaluate stream at the cell center rather than the grifpoint.
    ! If explicit, it seems we don't want to center streaming term (TODO: find
    ! out why this is and document it.)
    if (stream_implicit .or. driftkinetic_implicit) then
       ! call init_invert_stream_operator
       if (.not.allocated(stream_c)) allocate (stream_c(-nzgrid:nzgrid,nvpa,nspec))
       stream_c = stream
       do is = 1, nspec
          do iv = 1, nvpa
             call center_zed (iv, stream_c(:,iv,is))
          end do
       end do
       if (.not.allocated(gradpar_c)) allocate (gradpar_c(-nzgrid:nzgrid,-1:1))
       gradpar_c = spread(gradpar,2,3)
       ! get gradpar centred in zed for negative vpa (affects upwinding)
       call center_zed(1,gradpar_c(:,-stream_sign(1)))
       ! get gradpar centred in zed for positive vpa (affects upwinding)
       call center_zed(nvpa,gradpar_c(:,-stream_sign(nvpa)))
       stream = stream_c
    end if

    if (drifts_implicit_in_z) then
      ! Center wdriftx,y_g,phi - need the centered values to get the RHS and LHS
      ! of the GKE
      if (.not.allocated(wdriftx_g_c)) allocate(wdriftx_g_c(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not.allocated(wdrifty_g_c)) allocate(wdrifty_g_c(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not.allocated(wdriftx_phi_c)) allocate(wdriftx_phi_c(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not.allocated(wdrifty_phi_c)) allocate(wdrifty_phi_c(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      wdriftx_g_c = wdriftx_g
      wdrifty_g_c = wdrifty_g
      wdriftx_phi_c = wdriftx_phi
      wdrifty_phi_c = wdrifty_phi

      ia = 1

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        iv = iv_idx(vmu_lo,ivmu)
        call center_zed(iv, wdriftx_g_c(ia,:,ivmu))
        call center_zed(iv, wdrifty_g_c(ia,:,ivmu))
        call center_zed(iv, wdriftx_phi_c(ia,:,ivmu))
        call center_zed(iv, wdrifty_phi_c(ia,:,ivmu))
      end do
    end if

    ! Could include in the earlier if block, but we might want to separate wdrift
    ! and wstar in the future
    if (drifts_implicit_in_z) then
      ! Center wstar
      if (.not.allocated(wstar_c)) allocate(wstar_c(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      wstar_c = wstar
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        iv = iv_idx(vmu_lo,ivmu)
        call center_zed(iv, wstar_c(ia,:,ivmu))
      end do
    end if

    if(radial_variation) then
      call mp_abort("radial variation not currently supported electromagnetically, aborting")
      ! Taken from original parallel_streaming
      ! allocate (energy(-nzgrid:nzgrid))
      !
      ! if(.not.allocated(stream_rad_var1)) then
      !   allocate(stream_rad_var1(-nzgrid:nzgrid,nvpa,nspec))
      ! endif
      ! if(.not.allocated(stream_rad_var2)) then
      !   allocate(stream_rad_var2(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      !   stream_rad_var2 = 0.0
      ! endif
      ! ia=1
      ! stream_rad_var1 = -code_dt*spread(spread(spec%stm_psi0,1,nztot),2,nvpa) &
      !       * gfac*spread(spread(vpa,1,nztot)*spread(dgradpardrho,2,nvpa),3,nspec)
      ! do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      !   is  = is_idx(vmu_lo,ivmu)
      !   imu = imu_idx(vmu_lo,ivmu)
      !   iv  = iv_idx(vmu_lo,ivmu)
      !   energy = (vpa(iv)**2 + vperp2(ia,:,imu))*(spec(is)%temp_psi0/spec(is)%temp)
      !   stream_rad_var2(ia,:,ivmu) = &
      !           +code_dt*spec(is)%stm_psi0*vpa(iv)*gradpar &
      !           *spec(is)%zt*maxwell_vpa(iv,is)*maxwell_mu(ia,:,imu,is)*maxwell_fac(is) &
      !           *(  pfac*(spec(is)%fprim + spec(is)%tprim*(energy-2.5)) &
      !             + gfac*2*mu(imu)*dBdrho)
      ! enddo
      ! deallocate (energy)
    endif

    ! We always initialise the operator inversion - is this necessary?
    ! We could only do so if stream_implicit or driftkinetic_implicit is true.
    ! Currently, the explicit calculation doesn't do a sweep/inversion - is
    ! this consistent?
    call init_invert_operator

  end subroutine init_z_equation

  !> Initialise the coefficients for the response matrix
  !> solve. Whenever we calculate g, we're solving a matrix equation with a
  !> bidiagonal matrix. However, depending on the stream sign, this is either
  !> upper bidiagonal or lower bidiagnonal.
  !>
  !> The general equation we're solving is:
  !>   / b c 0 0 . . .  . . 0 \     / g(1) \       / RHS(1)  \
  !>  /  a b c 0 . . .      .  \   |  g(2)  |     |  RHS(2)  |
  !> /   0 a b c . . .      .   \  |  g(3)  |     |  RHS(3)  |
  !>|        . . .               | |    .   |   = |     .    |
  !> \   .      . . . 0 a b c 0 /  |    .   |     |     .    |
  !>  \  .       . . . 0 a b c /   |    .   |     |     .    |
  !>   \ 0 . .  . . . 0 0 a b /     \   .  /       \    .   /
  !>
  !> a, b, c correspond to dg/dt, streaming and wdrift terms on the LHS
  !> of the GKE. Each of these physics terms has a physics piece and a
  !> piece corresponding to the disretisation in z.
  !> Here we calculate the coefficients for both signs.
  !> For +ve stream sign (-ve advection velocity), a=0 (upper bidiagonal matrix)
  !> For -ve stream sign (+ve advection velocity), c=0 (lower bidiagonal matrix)
  subroutine init_invert_operator

    use zgrid, only: delzed
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: nsegments
    use run_parameters, only: zed_upwind, time_upwind

    implicit none

    integer :: nz, nseg_max

    nz = maxval(iz_up-iz_low)
    nseg_max = maxval(nsegments)

    if (.not.allocated(dgdt_tri_a)) then
       allocate (dgdt_tri_a(nz*nseg_max+1,-1:1)) ; dgdt_tri_a = 0.
       allocate (stream_tri_a(nz*nseg_max+1,-1:1)) ; stream_tri_a = 0.
       allocate (wdrift_tri_a(nz*nseg_max+1,-1:1)) ; wdrift_tri_a = 0.
       allocate (dgdt_tri_b(nz*nseg_max+1,-1:1)) ; dgdt_tri_b = 1.
       allocate (stream_tri_b(nz*nseg_max+1,-1:1)) ; stream_tri_b = 0.
       allocate (wdrift_tri_b(nz*nseg_max+1,-1:1)) ; wdrift_tri_b = 0.
       allocate (dgdt_tri_c(nz*nseg_max+1,-1:1)) ; dgdt_tri_c = 0.
       allocate (stream_tri_c(nz*nseg_max+1,-1:1)) ; stream_tri_c = 0.
       allocate (wdrift_tri_c(nz*nseg_max+1,-1:1)) ; wdrift_tri_c = 0.
    end if

    ! corresponds to sign of stream term positive on RHS of equation
    ! i.e., negative parallel advection speed
    ! NB: assumes equal spacing in zed
    dgdt_tri_b(:,1) = 0.5*(1.0+zed_upwind)
    stream_tri_b(:,1) = -1.0/delzed(0)
    dgdt_tri_c(:nz*nseg_max,1) = 0.5*(1.0-zed_upwind)
    stream_tri_c(:nz*nseg_max,1) = 1.0/delzed(0)
    ! corresponds to sign of stream term negative on RHS of equation
    ! NB: assumes equal spacing in zed
    dgdt_tri_b(:,-1) = 0.5*(1.0+zed_upwind)
    stream_tri_b(:,-1) = 1.0/delzed(0)
    dgdt_tri_a(2:,-1) = 0.5*(1.0-zed_upwind)
    stream_tri_a(2:,-1) = -1.0/delzed(0)

    stream_tri_a = 0.5*(1.0+time_upwind)*stream_tri_a
    stream_tri_b = 0.5*(1.0+time_upwind)*stream_tri_b
    stream_tri_c = 0.5*(1.0+time_upwind)*stream_tri_c
    ! The wdrift term contains both a z-upwinding term and a
    ! time-upwinding term.
    wdrift_tri_a = 0.5*(1.0+time_upwind)*dgdt_tri_a
    wdrift_tri_b = 0.5*(1.0+time_upwind)*dgdt_tri_b
    wdrift_tri_c = 0.5*(1.0+time_upwind)*dgdt_tri_c

  end subroutine init_invert_operator

  ! solve (I + (1+alph)/2*dt*vpa . grad)g^{n+1} = RHS
  ! g = RHS is input and overwritten by g = g^{n+1}
  subroutine invert_gke (ivmu, g)

    use mp, only : mp_abort
    use zgrid, only: nzgrid, ntubes
    use extended_zgrid, only: neigen
    use extended_zgrid, only: nsegments
    use extended_zgrid, only: nzed_segment
    use extended_zgrid, only: map_to_extended_zgrid
    use extended_zgrid, only: map_from_extended_zgrid
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, is_idx
    use kt_grids, only: naky
    use kt_grids, only: zonal_mode

    implicit none

    integer, intent (in) :: ivmu
    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: g

    integer :: iv, is
    integer :: iky, ie, it
    integer :: ulim, sgn
    complex, dimension (:), allocatable :: gext

    iv = iv_idx(vmu_lo,ivmu)
    is = is_idx(vmu_lo,ivmu)
    sgn = stream_sign(iv)

    do iky = 1, naky
       if (zonal_mode(iky)) then
          call mp_abort("sweep_zed_zonal not currently set up for electromagnetic, aborting")
          call sweep_zed_zonal (iv, is, sgn, g(iky,:,:,:))
       else
          do it = 1, ntubes
             do ie = 1, neigen(iky)
                allocate (gext(nsegments(ie,iky)*nzed_segment+1))
                ! get g on extended domain in zed
                call map_to_extended_zgrid (it, ie, iky, g(iky,:,:,:), gext, ulim)
                ! solve (I + (1+alph)/2*dt*vpa . grad)g_{inh}^{n+1} = RHS
                call z_tridiagonal_solve (iky, ie, ivmu, gext(:ulim))
                ! extract g from extended domain in zed
                call map_from_extended_zgrid (it, ie, iky, gext, g(iky,:,:,:))
                deallocate (gext)
             end do
          end do
       end if
    end do

  end subroutine invert_gke

  subroutine advance_z_implicit (g, phi, apar, bpar)

    use mp, only: proc0, mp_abort
    use job_manage, only: time_message
    use stella_layouts, only: vmu_lo, iv_idx
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx
    use dist_fn_arrays, only: g1
    use run_parameters, only: stream_matrix_inversion
    use fields, only: advance_fields, fields_updated

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar, bpar

    integer :: ivmu, iv
    complex, dimension (:,:,:,:), allocatable :: phi1, apar1, bpar1

    if (proc0) call time_message(.false.,time_parallel_streaming,' Stream advance')

    allocate (phi1(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (apar1(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (bpar1(naky,nakx,-nzgrid:nzgrid,ntubes))
    ! Bob: To make electromagnetic, need to replace phi with chi. However, because
    ! of how the gyroaverages are implemented, need to perform the replacement later.
    ! save the incoming g and phi, as they will be needed later
    g1 = g
    phi1 = phi
    apar1 = apar
    bpar1 = bpar

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)

       ! obtain RHS of inhomogeneous GK eqn;
       ! i.e., (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g_{inh}^{n+1}
       ! = (1-(1-alph)/2*dt*vpa*gradpar*d/dz)g^{n}
       ! + (1-alph)/2*dt*Ze*dlnF0/dE*exp(-vpa^2)*vpa*b.gradz*d<phi^{n}>/dz
       call get_gke_rhs (ivmu, g1(:,:,:,:,ivmu), phi1, phi, apar1, apar, bpar1, bpar, g(:,:,:,:,ivmu), eqn='inhomogeneous')

       if (stream_matrix_inversion) then
          ! solve (I + (1+alph)/2*dt*vpa . grad)g_{inh}^{n+1} = RHS
          ! g = RHS is input and overwritten by g = g_{inh}^{n+1}
          call invert_gke (ivmu, g(:,:,:,:,ivmu))
       else
          call mp_abort("sweep_g_zed currently not supported")
          ! call sweep_g_zed (ivmu, g(:,:,:,:,ivmu))
       end if
    end do

    fields_updated = .false.

    ! we now have g_{inh}^{n+1}
    ! calculate associated fields (phi_{inh}^{n+1})
    call advance_fields (g, phi, apar, bpar, dist='gbar')

    ! solve response_matrix*phi^{n+1} = phi_{inh}^{n+1}
    ! phi = phi_{inh}^{n+1} is input and overwritten by phi = phi^{n+1}
    call invert_implicitz_response (phi, apar, bpar)

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)

       ! now have phi^{n+1} for non-negative kx
       ! obtain RHS of GK eqn;
       ! i.e., (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g^{n+1}
       ! = (1-(1-alph)/2*dt*vpa*gradpar*d/dz)g^{n}
       ! + dt*Ze*dlnF0/dE*exp(-vpa^2)*vpa*b.gradz*d/dz((1+alph)/2*<phi^{n+1}>+(1-alph)/2*<phi^{n}>)
       call get_gke_rhs (ivmu, g1(:,:,:,:,ivmu), phi1, phi, apar1, apar, bpar1, bpar, g(:,:,:,:,ivmu), eqn='full')

       if (stream_matrix_inversion) then
          ! solve (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g^{n+1} = RHS
          ! g = RHS is input and overwritten by g = g^{n+1}
          call invert_gke (ivmu, g(:,:,:,:,ivmu))
       else
          call mp_abort("sweep_g_zed currently not supported")
          ! call sweep_g_zed (ivmu, g(:,:,:,:,ivmu))
       end if
    end do

    deallocate (phi1)

    if (proc0) call time_message(.false.,time_parallel_streaming,' Stream advance')

  end subroutine advance_z_implicit

  subroutine advance_parallel_streaming_explicit (g, gout)

    use mp, only: proc0, mp_abort
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use job_manage, only: time_message
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use species, only: spec
    use gyro_averages, only: gyro_average
    use fields, only: get_gyroaverage_chi
    use fields_arrays, only: phi, apar, bpar
    use run_parameters, only: driftkinetic_implicit, fapar, fbpar

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout

    integer :: ivmu, iv, imu, is, ia
    complex, dimension (:,:,:,:), allocatable :: g0, g1

    allocate (g0(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (g1(naky,nakx,-nzgrid:nzgrid,ntubes))
    ! parallel streaming stays in ky,kx,z space with ky,kx,z local
    if (proc0) call time_message(.false.,time_parallel_streaming,' Stream advance')

    ia = 1
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
    ! obtain <chi> (or <phi>-phi if driftkinetic_implicit=T)
       if (driftkinetic_implicit) then
         if ((fapar > epsilon(0.)) .or. (fbpar > epsilon(0.))) then
            call mp_abort ('driftkinetic_implicit not set up for apar, bpar. aborting')
         end if
         call gyro_average (phi, ivmu, g0(:,:,:,:))
         g0(:,:,:,:) = g0(:,:,:,:) - phi
       else
         call get_gyroaverage_chi(ivmu, phi, apar, bpar, g0)
       end if

    ! get d<phi>/dz, with z the parallel coordinate and store in g1
       call get_dgdz (g0, ivmu, g1)
    !    call get_dgdz_centered (g0, ivmu, g1)
    ! only want to treat vpar . grad (<phi>-phi)*F0 term explicitly
       if (driftkinetic_implicit) then
         g0 = 0.
       else
         call get_dgdz (g(:,:,:,:,ivmu), ivmu, g0)
         !call get_dgdz_centered (g(:,:,:,:,ivmu), ivmu, g0)
       end if

       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       g0(:,:,:,:) = g0(:,:,:,:) + g1(:,:,:,:)*spec(is)%zt &
            *maxwell_vpa(iv,is)*spread(spread(spread(maxwell_mu(ia,:,imu,is),1,naky),2,nakx),4,ntubes) &
            *maxwell_fac(is)

       ! multiply dg/dz with vpa*(b . grad z) and add to source (RHS of GK equation)
       call add_stream_term (g0, ivmu, gout(:,:,:,:,ivmu))
    end do
    if (proc0) call time_message(.false.,time_parallel_streaming,' Stream advance')
    deallocate (g0, g1)

  end subroutine advance_parallel_streaming_explicit

  !> Get the RHS of the GKE in z
  !> This looks like:
  !> RHS = g^n_(i*) - dt*(1-tupwind)/2*vpa*stm*gradpar_(i*)(dg^n/dz)_(i*)
  !>       + (1-tupwind)/2*zi*(ky*(wdrifty_g)_(i*) + kx*(wdriftx_g)_(i*))*g^n_(i*)
  !>       - dt*vpa*stm*gradpar_(i*)*Z/T*e**(-vpa**2)*e**(-vperp**2)_(i*)(d<chi^(n*)>/dz)_(i*)
  !>       + zi*(ky*(wdrifty_phi)_(i*) + kx*(wdriftx_phi)_(i*))<chi^(n*)>_(i*)
  !>       + zi*ky*(wstar)_(i*)*<chi^(n*)>_(i*)
  !> Currently this is designed to be flexible such that the streaming, wdrift and wstar
  !> terms can individually be included or not included. Also the "eqn" variable
  !> determines whether we are solving the inhomogeneous equation (setting the
  !> chi^(n+1) terms to zero) or the full GKE (including contributions from
  !> t^n and t^(n+1)).
  !> TODO: also include an option to solve the homogeneous equation (i.e. only
  !> including the chi^(n+1) terms) - this can then be called by get_dgdfield_matrix_column
  !> in the response matrix calculation, and remove code duplication in response_matrix.fpp
  !>
  !> May be able to reduce memory usage by more clever use of dummy variables
  subroutine get_gke_rhs (ivmu, gold, phiold, phi, aparold, apar, bparold, bpar, rhs, eqn)

    use mp, only : mp_abort
    use constants, only: zi
    use stella_time, only: code_dt
    use zgrid, only: nzgrid, ntubes
    use species, only: spec
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use kt_grids, only: naky, nakx, aky, akx
    ! use gyro_averages, only: gyro_average
    use vpamu_grids, only: vpa, mu
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    ! use stella_geometry, only: dbdzed
    use neoclassical_terms, only: include_neoclassical_terms
    use neoclassical_terms, only: dfneo_dvpa
    use run_parameters, only: time_upwind
    use run_parameters, only: driftkinetic_implicit, drifts_implicit_in_z, stream_implicit
    use run_parameters, only: maxwellian_inside_zed_derivative
    use fields, only: get_chi, get_gyroaverage_chi
    !use dist_fn_arrays, only: wdrifty_phi, wdriftx_phi, wdrifty_g, wdriftx_g, wstar
    !use stella_geometry, only: gradpar

    implicit none

    integer, intent (in) :: ivmu
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: gold
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phiold, phi, aparold, apar, bparold, bpar
    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: rhs
    character (*), intent (in) :: eqn

    integer :: iv, imu, is, iz, ia
    real :: tupwnd1, tupwnd2, stream_fac
    real, dimension (:), allocatable :: vpadf0dE_fac
    real, dimension (:), allocatable :: dummy1 !, dummy2
    complex, dimension (:,:,:,:), allocatable :: g, chi, dgdz, dchidz
    complex, dimension (:,:,:,:), allocatable :: field1, field2, field3

    allocate (vpadf0dE_fac(-nzgrid:nzgrid))
    allocate (dummy1(-nzgrid:nzgrid))
    ! allocate (dummy2(-nzgrid:nzgrid))
    allocate (g(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (dgdz(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (chi(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (dchidz(naky,nakx,-nzgrid:nzgrid,ntubes))

    ia = 1

    tupwnd1 = 0.5*(1.0-time_upwind)
    if (eqn=='full') then
       tupwnd2 = 0.5*(1.0+time_upwind)
    else
       ! Getting the inhomogeneous RHS - no contribution from t^(n+1)
       tupwnd2 = 0.0
    end if

    iv = iv_idx(vmu_lo,ivmu)
    imu = imu_idx(vmu_lo,ivmu)
    is = is_idx(vmu_lo,ivmu)

    !--------------------------------------------------------------------!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
    ! !!!First calculate g^n_(i*), <chi^(n*)>_(i*), (dg^n/dz)_(i*), !!!! !
    ! !!!  (d<chi^(n*)>/dz)_(i*)                                    !!!! !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
    !--------------------------------------------------------------------!

    ! obtain g^n_(i*) and store in g
    g = gold
    call center_zed (iv,g)

    ! obtain (dg^{n}/dz)_(i*) and store in dgdz
    ! NB: could eliminate this calculation at the expense of memory
    ! as this was calculated previously
    call get_dzed (iv,gold,dgdz)

    ! obtain <chi^(n*)>_(i) and store in chi
    ! First, get phi, apar, bpar time-centered:
    ! phi = (1+tupwind)/2*phi^{n+1} + (1-tupwind)/2*phi^{n}
    allocate (field1(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (field2(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (field3(naky,nakx,-nzgrid:nzgrid,ntubes))
    field1 = tupwnd1*phiold+tupwnd2*phi
    field2 = tupwnd1*aparold+tupwnd2*apar
    field3 = tupwnd1*bparold+tupwnd2*bpar

    ! set chi to be chi or <chi> depending on whether parallel streaming is
    ! implicit or only implicit in the kperp = 0 (drift kinetic) piece
    ! TODO: Check that the drift kinetic piece is correct.
    if (driftkinetic_implicit) then
       call mp_abort("Need to check that the electromagnetic drift kinetic implementation is correct.")
       !! Old:
       !g = field
       !! New:
       call get_chi(field1, field2, field3, ivmu, chi)
    else
      !! Old:
      ! call gyro_average (field, ivmu, g)
      !! New: Get gy_chi as a function of (kx, ky, z, tube), for a particular ivmu == a particular (iv, imu, is)
      call get_gyroaverage_chi(ivmu, field1, field2, field3, chi)

    end if
    deallocate (field1)
    deallocate (field2)
    deallocate (field3)

    ! obtain (d<chi^(n*)>/dz)_(i*) and store in dchidz
    call get_dzed (iv,chi,dchidz)

    ! obtain <chi^(n*)>_(i*) and store in chi
    call center_zed (iv,chi)

    !!! Now calculate the RHS term by term and store in rhs
    ! Start with g^n_(i*)
    rhs = g

    ! Add the streaming term
    ! streaming_term = - dt*(1-tupwind)/2*vpa*stm*gradpar_(i*)(dg^n/dz)_(i*)
    !                  - dt*vpa*stm*gradpar_(i*)*Z/T*e**(-vpa**2)*e**(-vperp**2)_(i*)(d<chi^(n*)>/dz)_(i*)
    if (stream_implicit) then

      ! stream_fac appears in front of both the g and chi terms.
      stream_fac = code_dt*spec(is)%stm_psi0

      ! NB: could do this once at beginning of simulation to speed things up
      ! this is vpa*Z/T*exp(-vpa^2)
      vpadf0dE_fac = vpa(iv)*spec(is)%zt*maxwell_vpa(iv,is)
      ! if including neoclassical correction to equilibrium distribution function
      ! then must also account for -vpa*dF_neo/dvpa*Z/T
      ! CHECK TO ENSURE THAT DFNEO_DVPA EXCLUDES EXP(-MU*B/T) FACTOR !!
      if (include_neoclassical_terms) then
         call mp_abort("Check that neoclassical terms are correct in electromagnetic implementation")
         do iz = -nzgrid, nzgrid
            vpadf0dE_fac(iz) = vpadf0dE_fac(iz)-0.5*dfneo_dvpa(1,iz,ivmu)*spec(is)%zt
         end do
         call center_zed (iv,vpadf0dE_fac)
      end if

      if (maxwellian_inside_zed_derivative) then
        call mp_abort("maxwellian_inside_zed_derivative not currently supported electromagnetically")
      else
        ! center Maxwellian factor in mu
        ! and store in dummy variable gp
        ! Could do this in initialisation?
        dummy1 = maxwell_mu(ia,:,imu,is)*maxwell_fac(is)
        call center_zed (iv,dummy1)
        ! multiply dchidz by Maxwellian factor and store in dchidz
        dchidz = dchidz*spread(spread(spread(dummy1,1,naky),2,nakx),4,ntubes)
      end if

      !! We've done this in initialisation, so don't need to do it here.
      ! obtain gradpar_(i*) and store in dummy variable dummy1
      ! dummy1 = gradpar
      ! call center_zed(iv, dummy1)

      do iz = -nzgrid, nzgrid
         rhs(:,:,iz,:) = rhs(:,:,iz,:) - stream_fac*gradpar_c(iz, stream_sign(iv)) &
              * (tupwnd1*vpa(iv)*dgdz(:,:,iz,:) + vpadf0dE_fac(iz)*dchidz(:,:,iz,:))
      end do

    end if

    ! Add the wdrift term
    ! wdrift_term = + (1-tupwind)/2*zi*(ky*(wdrifty_g)_(i*) + kx*(wdriftx_g)_(i*))*g^n_(i*)
    !               + zi*(ky*(wdrifty_phi)_(i*) + kx*(wdriftx_phi)_(i*))<chi^(n*)>_(i*)
    if (drifts_implicit_in_z) then
      !! We've done this in initialisation, so don't need to do it here.
      ! Center wdrifty_g and wdriftx_g. Store in dummy variables dummy1 and dummy2
      ! dummy1 = wdrifty_g(ia,:,ivmu)
      ! call center_zed(iv, dummy1)
      ! dummy2 = wdriftx_g(ia,:,ivmu)
      ! call center_zed(iv, dummy2)

      ! Add to RHS
      do iz = -nzgrid, nzgrid
        rhs(:,:,iz,:) = rhs(:,:,iz,:) + tupwnd1*zi*(spread(spread(aky,2,nakx),3,ntubes)*wdrifty_g_c(ia,iz,ivmu) &
                          + spread(spread(akx,1,naky),3,ntubes)*wdriftx_g_c(ia,iz,ivmu))*g(:,:,iz,:)
      end do
      !! We've done this in initialisation, so don't need to do it here.
      ! Center wdrifty_phi and wdriftx_phi. Store in dummy variables dummy1 and dummy2
      ! dummy1 = wdrifty_phi(ia,:,ivmu)
      ! call center_zed(iv, dummy1)
      ! dummy2 = wdriftx_phi(ia,:,ivmu)
      ! call center_zed(iv, dummy2)

      ! Add to RHS
      do iz = -nzgrid, nzgrid
        rhs(:,:,iz,:) = rhs(:,:,iz,:) + tupwnd1*zi*(spread(spread(aky,2,nakx),3,ntubes)*wdrifty_phi_c(ia,iz,ivmu) &
                          + spread(spread(akx,1,naky),3,ntubes)*wdriftx_phi_c(ia,iz,ivmu))*chi(:,:,iz,:)
      end do

    end if

    ! Add the wstar term
    ! wstar_term = + zi*ky*(wstar)_(i*)*<chi^(n*)>_(i*)
    ! Could include in the earlier if block, but we might want to separate wdrift
    ! and wstar in the future
    if (drifts_implicit_in_z) then
      !! We've done this in initialisation
      ! Center wstar, store in dummy variable dummy1
      !dummy1 = wstar(ia,:,ivmu)
      !call center_zed(iv, dummy1)

      ! Add to RHS
      do iz = -nzgrid, nzgrid
        rhs(:,:,iz,:) = rhs(:,:,iz,:) + spread(spread(aky,2,nakx),3,ntubes)*wstar_c(ia,iz,ivmu)*chi(:,:,iz,:)
      end do

    end if

    deallocate (vpadf0dE_fac)
    deallocate (dummy1)
    !deallocate (dummy2)
    deallocate (g)
    deallocate (dgdz)
    deallocate (chi)
    deallocate (dchidz)

  end subroutine get_gke_rhs

  !> Based on stream_tridiagonal_solve, this solves the bidiagonal matrix equation
  !> but with a combination of streaming and drifts.
  !> The equation we seek to solve is
  !> The general equation we're solving is:
  !>   MATRIX . g = RHS
  !> where MATRIX is an lower(upper) bidiagonal matrix for positive (negative)
  !> streaming velocity, corresponding to negative (postive) streaming sign.
  !> The RHS comes in as the variable g, and is overwritten by g. Both
  !> are defined on the extended z grid rather than the regular z grid.
  subroutine z_tridiagonal_solve (iky, ie, ivmu, g)

    use finite_differences, only: tridag
    use constants, only: zi
    use zgrid, only: nzgrid
    use extended_zgrid, only: iz_low, iz_up, ikxmod
    use extended_zgrid, only: nsegments
    use extended_zgrid, only: nzed_segment
    use kt_grids, only: aky, akx
    use mp, only : mp_abort
    use run_parameters, only: drifts_implicit_in_z, stream_implicit
    use dist_fn_arrays, only: wdrifty_phi, wdriftx_phi, wdrifty_g, wdriftx_g
    use stella_layouts, only: iv_idx, is_idx, vmu_lo

    implicit none

    integer, intent (in) :: iky, ie, ivmu
    complex, dimension (:), intent (in out) :: g

    integer :: iseg, llim, ulim, n, iv, is, ia
    integer :: nz, nseg_max, sgn, n_ext
    complex, dimension (:), allocatable :: a, b, c
    !real, dimension (:), allocatable :: dummy1, dummy2
    complex, dimension(:), allocatable :: wdrift_term

    ia = 1
    iv = iv_idx(vmu_lo,ivmu)
    is = is_idx(vmu_lo,ivmu)

    ! avoid double-counting at boundaries between 2pi segments
    nz = nzed_segment
    nseg_max = nsegments(ie,iky)
    sgn = stream_sign(iv)

    n_ext = nseg_max*nz+1
    allocate (a(n_ext))
    allocate (b(n_ext))
    allocate (c(n_ext))

    allocate(wdrift_term(-nzgrid:nzgrid))
    ! allocate(dummy1(-nzgrid:nzgrid))
    ! allocate(dummy2(-nzgrid:nzgrid))

    !!! Go through each segment, for each one calculating the contributions
    !!! to the LHS from the dg/dt term, the streaming term and the wdrift
    !!! term.
    iseg = 1
    llim = 1 ; ulim = nz+1
    ! Begin with the dg/dt terms
    a(llim:ulim) = dgdt_tri_a(llim:ulim,sgn)
    b(llim:ulim) = dgdt_tri_b(llim:ulim,sgn)
    c(llim:ulim) = dgdt_tri_c(llim:ulim,sgn)

    ! Add the streaming terms
    if (stream_implicit) then
      a(llim:ulim) = a(llim:ulim) - stream(iz_low(iseg):iz_up(iseg),iv,is)*stream_tri_a(llim:ulim,sgn)
      b(llim:ulim) = b(llim:ulim) - stream(iz_low(iseg):iz_up(iseg),iv,is)*stream_tri_b(llim:ulim,sgn)
      c(llim:ulim) = c(llim:ulim) - stream(iz_low(iseg):iz_up(iseg),iv,is)*stream_tri_c(llim:ulim,sgn)
    end if

    ! Add the wdrift terms
    if (drifts_implicit_in_z) then
      !! We've done this in the intialisation, so don't need to do here.
      ! Center wdrifty_g and wdriftx_g. Store in dummy variables dummy1 and dummy2
      ! dummy1 = wdrifty_g(ia,iz_low(iseg):iz_up(iseg),ivmu)
      ! call center_zed(iv, dummy1)
      ! dummy2 = wdriftx_g(ia,iz_low(iseg):iz_up(iseg),ivmu)
      ! call center_zed(iv, dummy2)
      ! Since we're on the extended zgrid, kx may vary from segment to segment
      ! (if BCs are linked). Use ikxmod to get the value of kx for this
      ! particular segment.
      wdrift_term = zi*(aky(iky)*wdrifty_g_c(ia,iz_low(iseg):iz_up(iseg),ivmu) &
                    + akx(ikxmod(iseg,ie,iky))*wdriftx_g_c(ia,iz_low(iseg):iz_up(iseg),ivmu))
      a(llim:ulim) = a(llim:ulim) - wdrift_tri_a(llim:ulim,sgn)*wdrift_term
      b(llim:ulim) = b(llim:ulim) - wdrift_tri_b(llim:ulim,sgn)*wdrift_term
      c(llim:ulim) = c(llim:ulim) - wdrift_tri_c(llim:ulim,sgn)*wdrift_term
    end if

    if (nsegments(ie,iky) > 1) then
       do iseg = 2, nsegments(ie,iky)
          llim = ulim+1
          ulim = llim+nz-1
          a(llim:ulim) = dgdt_tri_a(llim:ulim,sgn)
          b(llim:ulim) = dgdt_tri_b(llim:ulim,sgn)
          c(llim:ulim) = dgdt_tri_c(llim:ulim,sgn)

          if (stream_implicit) then
            a(llim:ulim) = a(llim:ulim) - stream(iz_low(iseg)+1:iz_up(iseg),iv,is)*stream_tri_a(llim:ulim,sgn)
            b(llim:ulim) = b(llim:ulim) - stream(iz_low(iseg)+1:iz_up(iseg),iv,is)*stream_tri_b(llim:ulim,sgn)
            c(llim:ulim) = c(llim:ulim) - stream(iz_low(iseg)+1:iz_up(iseg),iv,is)*stream_tri_c(llim:ulim,sgn)
          end if

          if (drifts_implicit_in_z) then
            !! We've done this in the intialisation, so don't need to do here.
            ! dummy1 = wdrifty_g(ia,iz_low(iseg)+1:iz_up(iseg),ivmu)
            ! call center_zed(iv, dummy1)
            ! dummy2 = wdriftx_g(ia,iz_low(iseg)+1:iz_up(iseg),ivmu)
            ! call center_zed(iv, dummy2)
            wdrift_term = zi*(aky(iky)*wdrifty_g_c(ia,iz_low(iseg)+1:iz_up(iseg),ivmu) &
                          + akx(ikxmod(iseg,ie,iky))*wdriftx_g_c(ia,iz_low(iseg)+1:iz_up(iseg),ivmu))
            a(llim:ulim) = a(llim:ulim) - wdrift_tri_a(llim:ulim,sgn)*wdrift_term
            b(llim:ulim) = b(llim:ulim) - wdrift_tri_b(llim:ulim,sgn)*wdrift_term
            c(llim:ulim) = c(llim:ulim) - wdrift_tri_c(llim:ulim,sgn)*wdrift_term
          end if
       end do
    end if

    ! Unclear why these are needed - seem to reproduce the calculations above
    ! n = size(stream_tri_a1,1)
    ! a(ulim) = stream_tri_a1(n,sgn)-stream(iz_up(nsegments(ie,iky)),iv,is)*stream_tri_a2(n,sgn)
    ! b(ulim) = stream_tri_b1(n,sgn)-stream(iz_up(nsegments(ie,iky)),iv,is)*stream_tri_b2(n,sgn)
    ! c(ulim) = 0. ! this line should not be necessary, as c(ulim) should not be accessed by tridag
    call tridag (1, a(:ulim), b(:ulim), c(:ulim), g)

    deallocate (a, b, c)
    deallocate (wdrift_term)
    !deallocate (dummy1, dummy2)

  end subroutine z_tridiagonal_solve

  !> Given the "inhomogeneous fields" phi, apar, bpar at t^(n+1) (i.e. those
  !> fields corresponding to g_inh^(n+1) ), find the full fields at t^(n+1)
  !> by solving the response matrix equation.
  !> This subroutine takes the inhomogeneous fields phi, apar, bpar, and
  !> replaces them with the full fields.
  subroutine invert_implicitz_response (phi, apar, bpar)

    use linear_solve, only: lu_back_substitution
    use zgrid, only: nzgrid, ntubes, nztot
    use extended_zgrid, only: neigen
    use extended_zgrid, only: nsegments
    use extended_zgrid, only: nzed_segment
    use extended_zgrid, only: map_to_extended_zgrid
    use extended_zgrid, only: map_from_extended_zgrid
    use extended_zgrid, only: ikxmod
    use kt_grids, only: naky, zonal_mode
    use fields_arrays, only: response_matrix

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar, bpar

    integer :: iky, ie, it, ulim, nz_segment
    integer :: ikx
    complex, dimension (:), allocatable :: fields

    ! need to put the fields into extended zed grid
    do iky = 1, naky
       ! avoid double counting of periodic endpoints for zonal modes
       if (zonal_mode(iky)) then
          allocate(fields(3*(nztot-1)))
          do it = 1, ntubes
             do ie = 1, neigen(iky)
                ikx = ikxmod(1,ie,iky)
                ! Because we want to avoid double-counting the endpoints,
                ! populate fields array ignoring the endpoint in
                ! phi, apar, bpar.
                fields(:nztot-1) = phi(iky,ikx,:nzgrid-1,it)
                fields(nztot:2*nztot-2) = apar(iky,ikx,:nzgrid-1,it)
                fields(2*nztot-1:3*nztot-3) = bpar(iky,ikx,:nzgrid-1,it)
                call lu_back_substitution (response_matrix(iky)%eigen(ie)%zloc, &
                     response_matrix(iky)%eigen(ie)%idx, fields)
                phi(iky,ikx,:nzgrid-1,it) = fields(:nztot-1)
                phi(iky,ikx,nzgrid,it) = phi(iky,ikx,-nzgrid,it)
                apar(iky,ikx,:nzgrid-1,it) = fields(nztot:2*nztot-2)
                apar(iky,ikx,nzgrid,it) = apar(iky,ikx,-nzgrid,it)
                bpar(iky,ikx,:nzgrid-1,it) = fields(2*nztot-1:3*nztot-3)
                bpar(iky,ikx,nzgrid,it) = bpar(iky,ikx,-nzgrid,it)
             end do
          end do
          deallocate(fields)
       else
          do it = 1, ntubes
             do ie = 1, neigen(iky)
                ! solve response_matrix*phi^{n+1} = phi_{inh}^{n+1}
                nz_segment = (nsegments(ie,iky)*nzed_segment+1)
                allocate (fields(3*nz_segment))
                call map_to_extended_zgrid (it, ie, iky, phi(iky,:,:,:), fields(:nz_segment), ulim)
                call map_to_extended_zgrid (it, ie, iky, apar(iky,:,:,:), fields(nz_segment+1:2*nz_segment), ulim)
                call map_to_extended_zgrid (it, ie, iky, bpar(iky,:,:,:), fields(2*nz_segment+1:3*nz_segment), ulim)
                call lu_back_substitution (response_matrix(iky)%eigen(ie)%zloc, &
                     response_matrix(iky)%eigen(ie)%idx, fields)
                call map_from_extended_zgrid (it, ie, iky, fields(:nz_segment), phi(iky,:,:,:))
                call map_from_extended_zgrid (it, ie, iky, fields(nz_segment+1:2*nz_segment), apar(iky,:,:,:))
                call map_from_extended_zgrid (it, ie, iky, fields(2*nz_segment+1:3*nz_segment), bpar(iky,:,:,:))
                deallocate (fields)
             end do
          end do
       end if
    end do

  end subroutine invert_implicitz_response

  subroutine add_stream_term (g, ivmu, src)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, is_idx
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx!, zonal_mode

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: src
    integer, intent (in) :: ivmu

    integer :: iv, is

    iv = iv_idx(vmu_lo,ivmu)
    is = is_idx(vmu_lo,ivmu)
    src(:,:,:,:) = src(:,:,:,:) + spread(spread(spread(stream(:,iv,is),1,naky),2,nakx),4,ntubes)*g(:,:,:,:)
  !  if (zonal_mode(1)) src(1,:,-nzgrid,:) = src(1,:,nzgrid,:)

  end subroutine add_stream_term

  ! Taken from parallel_streaming
  ! TODO: Make electromagnetic
  subroutine sweep_zed_zonal (iv, is, sgn, g)

    use zgrid, only: nzgrid, delzed, nztot, ntubes
    use kt_grids, only: nakx
    use run_parameters, only: zed_upwind, time_upwind

    implicit none

    integer, intent (in) :: iv, is, sgn
    complex, dimension (:,-nzgrid:,:), intent (in out) :: g

    integer :: iz, iz1, iz2
    real :: fac1, fac2
    complex, dimension (:), allocatable :: gcf
    complex, dimension (:,:,:), allocatable :: gpi

    allocate (gpi(nakx,-nzgrid:nzgrid,ntubes))
    allocate (gcf(-nzgrid:nzgrid))
    ! ky=0 is 2pi periodic (no extended zgrid)
    ! decompose into complementary function + particular integral
    ! zero BC for particular integral
    ! unit BC for complementary function (no source)
    if (sgn < 0) then
       iz1 = -nzgrid ; iz2 = nzgrid
    else
       iz1 = nzgrid ; iz2 = -nzgrid
    end if
    gpi(:,iz1,:) = 0. ; gcf(iz1) = 1.
    do iz = iz1-sgn, iz2, -sgn
       fac1 = 1.0+zed_upwind+sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
       fac2 = 1.0-zed_upwind-sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
       gpi(:,iz,:) = (-gpi(:,iz+sgn,:)*fac2 + 2.0*g(:,iz,:))/fac1
       gcf(iz) = -gcf(iz+sgn)*fac2/fac1
    end do
    ! g = g_PI + (g_PI(pi)/(1-g_CF(pi))) * g_CF
    g = gpi + (spread(gpi(:,iz2,:),2,nztot)/(1.-gcf(iz2)))*spread(spread(gcf,1,nakx),3,ntubes)
    deallocate (gpi, gcf)

  end subroutine sweep_zed_zonal

  subroutine finish_z_equation

    implicit none

    if (allocated(dgdt_tri_a)) then
       deallocate (dgdt_tri_a)
       deallocate (dgdt_tri_b)
       deallocate (dgdt_tri_c)
    end if

    if (allocated(stream_tri_a)) then
      deallocate (stream_tri_a)
      deallocate (stream_tri_b)
      deallocate (stream_tri_c)
    end if

    if (allocated(wdrift_tri_a)) then
      deallocate (wdrift_tri_a)
      deallocate (wdrift_tri_b)
      deallocate (wdrift_tri_c)
    end if

    if (allocated(stream_c)) then
      deallocate(stream_c)
    end if

    if (allocated(wdriftx_g_c)) then
      deallocate(wdriftx_g_c)
      deallocate(wdrifty_g_c)
      deallocate(wdriftx_phi_c)
      deallocate(wdrifty_phi_c)
    end if

    if (allocated(gradpar_c)) then
      deallocate(gradpar_c)
    end if

  end subroutine finish_z_equation

  !> TODO: Make electromagnetic.
  subroutine add_parallel_streaming_radial_variation (g, gout, rhs)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use job_manage, only: time_message
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use species, only: spec
    use gyro_averages, only: gyro_average, gyro_average_j1
    use fields_arrays, only: phi, phi_corr_QN, phi_corr_GA

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout

    !the next input/output is for quasineutrality and gyroaveraging corrections
    !that go directly in the RHS (since they don't require further FFTs)
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: rhs

    integer :: ivmu, iv, imu, is, it, ia, iz

    complex, dimension (:,:,:,:), allocatable :: g0, g1, g2, g3
    complex, dimension (:,:), allocatable :: g0k

    allocate (g0(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (g1(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (g2(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (g3(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (g0k(naky,nakx)); g0k =0

    ! parallel streaming stays in ky,kx,z space with ky,kx,z local

    ia = 1
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
    ! obtain <phi>
    ! get d<phi>/dz, with z the parallel coordinate and store in g1
       call gyro_average (phi, ivmu, g0)
       call get_dgdz_variable (g0, ivmu, g1)

    ! get variation in gyroaveraging and store in g2
       call get_dgdz_variable (phi_corr_GA(:,:,:,:,ivmu), ivmu, g2)

    ! get variation in quasineutrality and store in g3
       call gyro_average (phi_corr_QN, ivmu, g0)
       call get_dgdz_variable (g0,  ivmu, g3)

       call get_dgdz_variable (g(:,:,:,:,ivmu),  ivmu, g0)

       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do it = 1, ntubes
         do iz = -nzgrid,nzgrid

    !!#1 - variation in gradpar
           g0k = g0(:,:,iz,it) &
               + g1(:,:,iz,it)*spec(is)%zt*maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is)

           g0k = g0k*stream_rad_var1(iz,iv,is)

    !!#2 - variation in F_s/T_s
           g0k = g0k + g1(:,:,iz,it)*stream_rad_var2(ia,iz,ivmu)

           gout(:,:,iz,it,ivmu) = gout(:,:,iz,it,ivmu) + g0k

    !!#3 - variation in the gyroaveraging and quasineutrality of phi
    !!     These variations already have the linear part calculated, so
    !!     ad it into the rhs directly
           g0k = spec(is)%zt*stream(iz,iv,is)*maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is) &
                 *(g2(:,:,iz,it) + g3(:,:,iz,it))

           rhs(:,:,iz,it,ivmu) = rhs(:,:,iz,it,ivmu) + g0k

         enddo
       enddo
    end do
    deallocate (g0, g1, g2, g3, g0k)

  end subroutine add_parallel_streaming_radial_variation

  subroutine center_zed_extended (iv, g)

    use finite_differences, only: cell_centres_zed
    use kt_grids, only: naky, nakx
    use zgrid, only: nzgrid, ntubes
    use extended_zgrid, only: neigen, nsegments
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: fill_zed_ghost_zones
    use run_parameters, only: zed_upwind

    implicit none

    integer, intent (in) :: iv
    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: g

    integer :: iky, ie, iseg, it
    complex, dimension (2) :: gleft, gright
    complex, dimension (:,:,:), allocatable :: gc

    allocate (gc(nakx,-nzgrid:nzgrid,ntubes))

    do iky = 1, naky
       do it = 1, ntubes
          do ie = 1, neigen(iky)
             do iseg = 1, nsegments(ie,iky)
                ! first fill in ghost zones at boundaries in g(z)
                call fill_zed_ghost_zones (it, iseg, ie, iky, g, gleft, gright)
                ! get cell centres values
                call cell_centres_zed (iz_low(iseg), &
                     g(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it), &
                     zed_upwind, stream_sign(iv), gleft(2), gright(1), &
                     gc(ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it))
             end do
          end do
       end do
       g(iky,:,:,:) = gc
    end do

    deallocate (gc)

  end subroutine center_zed_extended

  subroutine center_zed_segment_real (iv, g)

    use zgrid, only: nzgrid
    use run_parameters, only: zed_upwind

    integer, intent (in) :: iv
    real, dimension (-nzgrid:), intent (in out) :: g

    if (stream_sign(iv) > 0) then
       g(:nzgrid-1) = 0.5*((1.+zed_upwind)*g(:nzgrid-1) + (1.-zed_upwind)*g(-nzgrid+1:))
       g(nzgrid) = g(-nzgrid)
    else
       g(-nzgrid+1:) = 0.5*((1.-zed_upwind)*g(:nzgrid-1) + (1.+zed_upwind)*g(-nzgrid+1:))
       g(-nzgrid) = g(nzgrid)
    end if

  end subroutine center_zed_segment_real

  subroutine get_dzed (iv, g, dgdz)

    use finite_differences, only: fd_cell_centres_zed
    use kt_grids, only: naky
    use zgrid, only: nzgrid, delzed, ntubes
    use extended_zgrid, only: neigen, nsegments
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: fill_zed_ghost_zones

    implicit none

    integer, intent (in) :: iv
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: dgdz

    integer :: iky, ie, iseg, it
    complex, dimension (2) :: gleft, gright

    do it = 1, ntubes
       do iky = 1, naky
          do ie = 1, neigen(iky)
             do iseg = 1, nsegments(ie,iky)
                ! first fill in ghost zones at boundaries in g(z)
                call fill_zed_ghost_zones (it, iseg, ie, iky, g, gleft, gright)
                ! get finite difference approximation for dg/dz at cell centres
                ! iv > nvgrid corresponds to positive vpa, iv <= nvgrid to negative vpa
                call fd_cell_centres_zed (iz_low(iseg), &
                     g(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it), &
                     delzed(0), stream_sign(iv), gleft(2), gright(1), &
                     dgdz(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it))
             end do
          end do
       end do
    end do

  end subroutine get_dzed

  subroutine get_dgdz (g, ivmu, dgdz)

    use finite_differences, only: third_order_upwind_zed
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx
    use zgrid, only: nzgrid, delzed, ntubes
    use extended_zgrid, only: neigen, nsegments
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: fill_zed_ghost_zones
    use extended_zgrid, only: periodic
    use kt_grids, only: naky

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: dgdz
    integer, intent (in) :: ivmu

    integer :: iseg, ie, it, iky, iv
    complex, dimension (2) :: gleft, gright

    ! FLAG -- assuming delta zed is equally spaced below!
    iv = iv_idx(vmu_lo,ivmu)
    do iky = 1, naky
      do it = 1, ntubes
        do ie = 1, neigen(iky)
          do iseg = 1, nsegments(ie,iky)
            ! first fill in ghost zones at boundaries in g(z)
            call fill_zed_ghost_zones (it, iseg, ie, iky, g(:,:,:,:), gleft, gright)
            ! now get dg/dz
            call third_order_upwind_zed (iz_low(iseg), iseg, nsegments(ie,iky), &
                 g(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it), &
                 delzed(0), stream_sign(iv), gleft, gright, periodic(iky), &
                 dgdz(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it))
            end do
          end do
      end do
    end do

  end subroutine get_dgdz

  subroutine get_dgdz_variable (g, ivmu, dgdz)

     use finite_differences, only: fd_variable_upwinding_zed
     use stella_layouts, only: vmu_lo
     use stella_layouts, only: iv_idx
     use zgrid, only: nzgrid, delzed, ntubes
     use extended_zgrid, only: neigen, nsegments
     use extended_zgrid, only: iz_low, iz_up
     use extended_zgrid, only: ikxmod
     use extended_zgrid, only: fill_zed_ghost_zones
     use extended_zgrid, only: periodic
     use run_parameters, only: zed_upwind
     use kt_grids, only: naky

     implicit none

     complex, dimension (:,:,-nzgrid:,:), intent (in) :: g
     complex, dimension (:,:,-nzgrid:,:), intent (out) :: dgdz
     integer, intent (in) :: ivmu

     integer :: iseg, ie, iky, iv, it
     complex, dimension (2) :: gleft, gright
     ! FLAG -- assuming delta zed is equally spaced below!
      iv = iv_idx(vmu_lo,ivmu)
      do iky = 1, naky
        do it = 1, ntubes
          do ie = 1, neigen(iky)
            do iseg = 1, nsegments(ie,iky)
               ! first fill in ghost zones at boundaries in g(z)
               call fill_zed_ghost_zones (it, iseg, ie, iky, g(:,:,:,:), gleft, gright)
               ! now get dg/dz
               call fd_variable_upwinding_zed (iz_low(iseg), iseg, nsegments(ie,iky), &
                    g(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it), &
                    delzed(0), stream_sign(iv), zed_upwind,gleft, gright, periodic(iky), &
                    dgdz(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it))
            end do
          end do
        enddo
      end do
  end subroutine get_dgdz_variable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The following are from the original parallel_streaming,
!!!! should now be redundant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Initialise the coefficients for the response matrix
!   !> solve. Whenever we calculate g, we're solving a matrix equation with a
!   !> bidiagonal matrix. However, depending on the stream sign, this is either
!   !> upper bidiagonal or lower bidiagnonal.
!   !>
!   !> The general equation we're solving is:
!   !>   / b c 0 0 . . .  . . 0 \     / g(1) \       / RHS(1)  \
!   !>  /  a b c 0 . . .      .  \   |  g(2)  |     |  RHS(2)  |
!   !> /   0 a b c . . .      .   \  |  g(3)  |     |  RHS(3)  |
!   !>|        . . .               | |    .   |   = |     .    |
!   !> \   .      . . . 0 a b c 0 /  |    .   |     |     .    |
!   !>  \  .       . . . 0 a b c /   |    .   |     |     .    |
!   !>   \ 0 . .  . . . 0 0 a b /     \   .  /       \    .   /
!   !>
!   !> For +ve stream sign (-ve advection velocity), a=0 (upper bidiagonal matrix)
!   !> For -ve stream sign (+ve advection velocity), c=0 (lower bidiagonal matrix)
!   !> Here we calculate the coefficients for both signs.
!   !> The "1" and "2" correspond to the discretisation of the d/dt and the d/dz terms
!   !> respectively
!   subroutine init_invert_stream_operator
!
!     use zgrid, only: delzed
!     use extended_zgrid, only: iz_low, iz_up
!     use extended_zgrid, only: nsegments
!     use run_parameters, only: zed_upwind, time_upwind
!
!     implicit none
!
!     integer :: nz, nseg_max
!
!     nz = maxval(iz_up-iz_low)
!     nseg_max = maxval(nsegments)
!
!     if (.not.allocated(stream_tri_a1)) then
!        allocate (stream_tri_a1(nz*nseg_max+1,-1:1)) ; stream_tri_a1 = 0.
!        allocate (stream_tri_a2(nz*nseg_max+1,-1:1)) ; stream_tri_a2 = 0.
!        allocate (stream_tri_b1(nz*nseg_max+1,-1:1)) ; stream_tri_b1 = 1.
!        allocate (stream_tri_b2(nz*nseg_max+1,-1:1)) ; stream_tri_b2 = 0.
!        allocate (stream_tri_c1(nz*nseg_max+1,-1:1)) ; stream_tri_c1 = 0.
!        allocate (stream_tri_c2(nz*nseg_max+1,-1:1)) ; stream_tri_c2 = 0.
!     end if
!
!     ! corresponds to sign of stream term positive on RHS of equation
!     ! i.e., negative parallel advection speed
!     ! NB: assumes equal spacing in zed
!     stream_tri_b1(:,1) = 0.5*(1.0+zed_upwind)
!     stream_tri_b2(:,1) = -1.0/delzed(0)
!     stream_tri_c1(:nz*nseg_max,1) = 0.5*(1.0-zed_upwind)
!     stream_tri_c2(:nz*nseg_max,1) = 1.0/delzed(0)
!     ! corresponds to sign of stream term negative on RHS of equation
!     ! NB: assumes equal spacing in zed
!     stream_tri_b1(:,-1) = 0.5*(1.0+zed_upwind)
!     stream_tri_b2(:,-1) = 1.0/delzed(0)
!     stream_tri_a1(2:,-1) = 0.5*(1.0-zed_upwind)
!     stream_tri_a2(2:,-1) = -1.0/delzed(0)
!
!     stream_tri_a2 = 0.5*(1.0+time_upwind)*stream_tri_a2
!     stream_tri_b2 = 0.5*(1.0+time_upwind)*stream_tri_b2
!     stream_tri_c2 = 0.5*(1.0+time_upwind)*stream_tri_c2
!
!   end subroutine init_invert_stream_operator
!
!   subroutine advance_parallel_streaming_implicit (g, phi, apar, bpar)
!
!     use mp, only: proc0
!     use job_manage, only: time_message
!     use stella_layouts, only: vmu_lo, iv_idx
!     use zgrid, only: nzgrid, ntubes
!     use kt_grids, only: naky, nakx
!     use dist_fn_arrays, only: g1
!     use run_parameters, only: stream_matrix_inversion
!     use fields, only: advance_fields, fields_updated
!
!     implicit none
!
!     complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
!     complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar, bpar
!
!     integer :: ivmu, iv
!     complex, dimension (:,:,:,:), allocatable :: phi1, apar1, bpar1
!
!     if (proc0) call time_message(.false.,time_parallel_streaming,' Stream advance')
!
!     allocate (phi1(naky,nakx,-nzgrid:nzgrid,ntubes))
!     allocate (apar1(naky,nakx,-nzgrid:nzgrid,ntubes))
!     allocate (bpar1(naky,nakx,-nzgrid:nzgrid,ntubes))
!     ! Bob: To make electromagnetic, need to replace phi with chi. However, because
!     ! of how the gyroaverages are implemented, need to perform the replacement later.
!     ! save the incoming g and phi, as they will be needed later
!     g1 = g
!     phi1 = phi
!     apar1 = apar
!     bpar1 = bpar
!
!     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!        iv = iv_idx(vmu_lo,ivmu)
!
!        ! obtain RHS of inhomogeneous GK eqn;
!        ! i.e., (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g_{inh}^{n+1}
!        ! = (1-(1-alph)/2*dt*vpa*gradpar*d/dz)g^{n}
!        ! + (1-alph)/2*dt*Ze*dlnF0/dE*exp(-vpa^2)*vpa*b.gradz*d<phi^{n}>/dz
!        call get_gke_rhs (ivmu, g1(:,:,:,:,ivmu), phi1, phi, apar1, apar, bpar1, bpar, g(:,:,:,:,ivmu), eqn='inhomogeneous')
!
!        if (stream_matrix_inversion) then
!           ! solve (I + (1+alph)/2*dt*vpa . grad)g_{inh}^{n+1} = RHS
!           ! g = RHS is input and overwritten by g = g_{inh}^{n+1}
!           call invert_parstream (ivmu, g(:,:,:,:,ivmu))
!        else
!           call sweep_g_zed (ivmu, g(:,:,:,:,ivmu))
!        end if
!     end do
!
!     fields_updated = .false.
!
!     ! we now have g_{inh}^{n+1}
!     ! calculate associated fields (phi_{inh}^{n+1})
!     call advance_fields (g, phi, apar, bpar, dist='gbar')
!
!     ! solve response_matrix*phi^{n+1} = phi_{inh}^{n+1}
!     ! phi = phi_{inh}^{n+1} is input and overwritten by phi = phi^{n+1}
!     call invert_parstream_response (phi, apar, bpar)
!
!     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!        iv = iv_idx(vmu_lo,ivmu)
!
!        ! now have phi^{n+1} for non-negative kx
!        ! obtain RHS of GK eqn;
!        ! i.e., (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g^{n+1}
!        ! = (1-(1-alph)/2*dt*vpa*gradpar*d/dz)g^{n}
!        ! + dt*Ze*dlnF0/dE*exp(-vpa^2)*vpa*b.gradz*d/dz((1+alph)/2*<phi^{n+1}>+(1-alph)/2*<phi^{n}>)
!        call get_gke_rhs (ivmu, g1(:,:,:,:,ivmu), phi1, phi, apar1, apar, bpar1, bpar, g(:,:,:,:,ivmu), eqn='full')
!
!        if (stream_matrix_inversion) then
!           ! solve (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g^{n+1} = RHS
!           ! g = RHS is input and overwritten by g = g^{n+1}
!           call invert_parstream (ivmu, g(:,:,:,:,ivmu))
!        else
!           call sweep_g_zed (ivmu, g(:,:,:,:,ivmu))
!        end if
!     end do
!
!     deallocate (phi1)
!
!     if (proc0) call time_message(.false.,time_parallel_streaming,' Stream advance')
!
!   end subroutine advance_parallel_streaming_implicit
!
!   subroutine get_gke_rhs (ivmu, gold, phiold, phi, aparold, apar, bparold, bpar, g, eqn)
!
!     use stella_time, only: code_dt
!     use zgrid, only: nzgrid, ntubes
!     use species, only: spec
!     use stella_layouts, only: vmu_lo
!     use stella_layouts, only: iv_idx, imu_idx, is_idx
!     use kt_grids, only: naky, nakx
!     use gyro_averages, only: gyro_average
!     use vpamu_grids, only: vpa, mu
!     use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
!     use stella_geometry, only: dbdzed
!     use neoclassical_terms, only: include_neoclassical_terms
!     use neoclassical_terms, only: dfneo_dvpa
!     use run_parameters, only: time_upwind
!     use run_parameters, only: driftkinetic_implicit
!     use run_parameters, only: maxwellian_inside_zed_derivative
!     !! Bob: Need to compute chi & its gyroaverage
!     use fields, only: get_chi, get_gyroaverage_chi
!
!     implicit none
!
!     integer, intent (in) :: ivmu
!     complex, dimension (:,:,-nzgrid:,:), intent (in) :: gold
!     complex, dimension (:,:,-nzgrid:,:), intent (in) :: phiold, phi, aparold, apar, bparold, bpar
!     complex, dimension (:,:,-nzgrid:,:), intent (in out) :: g
!     character (*), intent (in) :: eqn
!
!     integer :: iv, imu, is, iz, ia
!     real :: tupwnd1, tupwnd2, fac
!     real, dimension (:), allocatable :: vpadf0dE_fac
!     real, dimension (:), allocatable :: gp
!     complex, dimension (:,:,:,:), allocatable :: dgdz, dphidz
!     complex, dimension (:,:,:,:), allocatable :: field1, field2, field3
!
!     allocate (vpadf0dE_fac(-nzgrid:nzgrid))
!     allocate (gp(-nzgrid:nzgrid))
!     allocate (dgdz(naky,nakx,-nzgrid:nzgrid,ntubes))
!     allocate (dphidz(naky,nakx,-nzgrid:nzgrid,ntubes))
!
!     ia = 1
!
!     tupwnd1 = 0.5*(1.0-time_upwind)
!     if (eqn=='full') then
!        tupwnd2 = 0.5*(1.0+time_upwind)
!     else
!        tupwnd2 = 0.0
!     end if
!
!     ! now have phi^{n+1} for non-negative kx
!     ! obtain RHS of GK eqn;
!     ! i.e., (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g^{n+1}
!     ! = (1-(1-alph)/2*dt*vpa*gradpar*d/dz)g^{n}
!     ! + dt*Ze*dlnF0/dE*exp(-vpa^2)*vpa*b.gradz*d/dz((1+alph)/2*<phi^{n+1}>+(1-alph)/2*<phi^{n}>
!     iv = iv_idx(vmu_lo,ivmu)
!     imu = imu_idx(vmu_lo,ivmu)
!     is = is_idx(vmu_lo,ivmu)
!
!     ! obtain dg^{n}/dz and store in dgdz
!     ! NB: could eliminate this calculation at the expense of memory
!     ! as this was calculated previously
!     call get_dzed (iv,gold,dgdz)
!
!     allocate (field1(naky,nakx,-nzgrid:nzgrid,ntubes))
!     allocate (field2(naky,nakx,-nzgrid:nzgrid,ntubes))
!     allocate (field3(naky,nakx,-nzgrid:nzgrid,ntubes))
!     ! old: get <phi> = (1+alph)/2*<phi^{n+1}> + (1-alph)/2*<phi^{n}>
!     field1 = tupwnd1*phiold+tupwnd2*phi
!     ! Bob: do same for apar, bpar.
!     field2 = tupwnd1*aparold+tupwnd2*apar
!     field3 = tupwnd1*bparold+tupwnd2*bpar
!
!     ! set g to be phi or <phi> depending on whether parallel streaming is
!     ! implicit or only implicit in the kperp = 0 (drift kinetic) piece
!     ! Bob: we want g to be chi or <chi>
!     if (driftkinetic_implicit) then
!        !! Old:
!        !g = field
!        !! New:
!        call get_chi(field1, field2, field3, ivmu, g)
!     else
!       !! Old:
!       ! call gyro_average (field, ivmu, g)
!       !! New: Get gy_chi as a function of (kx, ky, z, tube), for a particular ivmu == a particular (iv, imu, is)
!       call get_gyroaverage_chi(ivmu, field1, field2, field3, g)
!
!     end if
!     deallocate (field1)
!     deallocate (field2)
!     deallocate (field3)
!
!     if (maxwellian_inside_zed_derivative) then
!        ! obtain d(exp(-mu*B/T)*<phi>)/dz and store in dphidz
!        g = g*spread(spread(spread(maxwell_mu(ia,:,imu,is)*maxwell_fac(is),1,naky),2,nakx),4,ntubes)
!        ! Bob: this is actually dchidz
!        call get_dzed (iv,g,dphidz)
!        ! get <phi>*exp(-mu*B/T)*dB/dz at cell centres
!        g = g*spread(spread(spread(dbdzed(ia,:),1,naky),2,nakx),4,ntubes)
!        call center_zed (iv,g)
!        ! construct d(<phi>*exp(-mu*B/T))/dz + 2*mu*<phi>*exp(-mu*B/T)*dB/dz
!        ! = d<phi>/dz * exp(-mu*B/T)
!        dphidz = dphidz + 2.0*mu(imu)*g
!     else
!        ! obtain d<phi>/dz and store in dphidz
!        call get_dzed (iv,g,dphidz)
!        ! center Maxwellian factor in mu
!        ! and store in dummy variable gp
!        gp = maxwell_mu(ia,:,imu,is)*maxwell_fac(is)
!        call center_zed (iv,gp)
!        ! multiply by Maxwellian factor
!        dphidz = dphidz*spread(spread(spread(gp,1,naky),2,nakx),4,ntubes)
!     end if
!
!     ! NB: could do this once at beginning of simulation to speed things up
!     ! this is vpa*Z/T*exp(-vpa^2)
!     vpadf0dE_fac = vpa(iv)*spec(is)%zt*maxwell_vpa(iv,is)
!     ! if including neoclassical correction to equilibrium distribution function
!     ! then must also account for -vpa*dF_neo/dvpa*Z/T
!     ! CHECK TO ENSURE THAT DFNEO_DVPA EXCLUDES EXP(-MU*B/T) FACTOR !!
!     if (include_neoclassical_terms) then
!        do iz = -nzgrid, nzgrid
!           vpadf0dE_fac(iz) = vpadf0dE_fac(iz)-0.5*dfneo_dvpa(1,iz,ivmu)*spec(is)%zt
!        end do
!        call center_zed (iv,vpadf0dE_fac)
!     end if
!
!     g = gold
!     call center_zed (iv,g)
!
!     if (stream_sign(iv) > 0) then
!        gp = gradpar_c(:,-1)
!     else
!        gp = gradpar_c(:,1)
!     end if
!
!     ! construct RHS of GK eqn
!     fac = code_dt*spec(is)%stm_psi0
!     do iz = -nzgrid, nzgrid
!        g(:,:,iz,:) = g(:,:,iz,:) - fac*gp(iz) &
!             * (tupwnd1*vpa(iv)*dgdz(:,:,iz,:) + vpadf0dE_fac(iz)*dphidz(:,:,iz,:))
!     end do
!
!     deallocate (vpadf0dE_fac, gp)
!     deallocate (dgdz, dphidz)
!
!   end subroutine get_gke_rhs
!
!   ! solve (I + (1+alph)/2*dt*vpa . grad)g^{n+1} = RHS
!   ! g = RHS is input and overwritten by g = g^{n+1}
!   subroutine invert_parstream (ivmu, g)
!
!     use zgrid, only: nzgrid, ntubes
!     use extended_zgrid, only: neigen
!     use extended_zgrid, only: nsegments
!     use extended_zgrid, only: nzed_segment
!     use extended_zgrid, only: map_to_extended_zgrid
!     use extended_zgrid, only: map_from_extended_zgrid
!     use stella_layouts, only: vmu_lo
!     use stella_layouts, only: iv_idx, is_idx
!     use kt_grids, only: naky
!     use kt_grids, only: zonal_mode
!
!     implicit none
!
!     integer, intent (in) :: ivmu
!     complex, dimension (:,:,-nzgrid:,:), intent (in out) :: g
!
!     integer :: iv, is
!     integer :: iky, ie, it
!     integer :: ulim, sgn
!     complex, dimension (:), allocatable :: gext
!
!     iv = iv_idx(vmu_lo,ivmu)
!     is = is_idx(vmu_lo,ivmu)
!     sgn = stream_sign(iv)
!
!     do iky = 1, naky
!        if (zonal_mode(iky)) then
!           call sweep_zed_zonal (iv, is, sgn, g(iky,:,:,:))
!        else
!           do it = 1, ntubes
!              do ie = 1, neigen(iky)
!                 allocate (gext(nsegments(ie,iky)*nzed_segment+1))
!                 ! get g on extended domain in zed
!                 call map_to_extended_zgrid (it, ie, iky, g(iky,:,:,:), gext, ulim)
!                 ! solve (I + (1+alph)/2*dt*vpa . grad)g_{inh}^{n+1} = RHS
!                 call stream_tridiagonal_solve (iky, ie, iv, is, gext(:ulim))
!                 ! extract g from extended domain in zed
!                 call map_from_extended_zgrid (it, ie, iky, gext, g(iky,:,:,:))
!                 deallocate (gext)
!              end do
!           end do
!        end if
!     end do
!
!   end subroutine invert_parstream
!
!   ! Bob: stream_tridiagonal_solve redundant (replaced by z_tridiagonal_solve)
!   subroutine stream_tridiagonal_solve (iky, ie, iv, is, g)
!
!     use finite_differences, only: tridag
!     use extended_zgrid, only: iz_low, iz_up
!     use extended_zgrid, only: nsegments
!     use extended_zgrid, only: nzed_segment
!
!     implicit none
!
!     integer, intent (in) :: iky, ie, iv, is
!     complex, dimension (:), intent (in out) :: g
!
!     integer :: iseg, llim, ulim, n
!     integer :: nz, nseg_max, sgn, n_ext
!     real, dimension (:), allocatable :: a, b, c
!
!     ! avoid double-counting at boundaries between 2pi segments
!     nz = nzed_segment
!     nseg_max = nsegments(ie,iky)
!     sgn = stream_sign(iv)
!
!     n_ext = nseg_max*nz+1
!     allocate (a(n_ext))
!     allocate (b(n_ext))
!     allocate (c(n_ext))
!
!     iseg = 1
!     llim = 1 ; ulim = nz+1
!     a(llim:ulim) = stream_tri_a1(llim:ulim,sgn) &
!          - stream(iz_low(iseg):iz_up(iseg),iv,is)*stream_tri_a2(llim:ulim,sgn)
!     b(llim:ulim) = stream_tri_b1(llim:ulim,sgn) &
!          -stream(iz_low(iseg):iz_up(iseg),iv,is)*stream_tri_b2(llim:ulim,sgn)
!     c(llim:ulim) = stream_tri_c1(llim:ulim,sgn) &
!          -stream(iz_low(iseg):iz_up(iseg),iv,is)*stream_tri_c2(llim:ulim,sgn)
!
!     if (nsegments(ie,iky) > 1) then
!        do iseg = 2, nsegments(ie,iky)
!           llim = ulim+1
!           ulim = llim+nz-1
!           a(llim:ulim) = stream_tri_a1(llim:ulim,sgn) &
!                - stream(iz_low(iseg)+1:iz_up(iseg),iv,is)*stream_tri_a2(llim:ulim,sgn)
!           b(llim:ulim) = stream_tri_b1(llim:ulim,sgn) &
!                - stream(iz_low(iseg)+1:iz_up(iseg),iv,is)*stream_tri_b2(llim:ulim,sgn)
!           c(llim:ulim) = stream_tri_c1(llim:ulim,sgn) &
!                - stream(iz_low(iseg)+1:iz_up(iseg),iv,is)*stream_tri_c2(llim:ulim,sgn)
!        end do
!     end if
!     n = size(stream_tri_a1,1)
!     a(ulim) = stream_tri_a1(n,sgn)-stream(iz_up(nsegments(ie,iky)),iv,is)*stream_tri_a2(n,sgn)
!     b(ulim) = stream_tri_b1(n,sgn)-stream(iz_up(nsegments(ie,iky)),iv,is)*stream_tri_b2(n,sgn)
!     c(ulim) = 0. ! this line should not be necessary, as c(ulim) should not be accessed by tridag
!     call tridag (1, a(:ulim), b(:ulim), c(:ulim), g)
!
!     deallocate (a, b, c)
!
!   end subroutine stream_tridiagonal_solve
!
!   ! g= RHS of gke is input
!   ! g = g^{n+1} is output
!   subroutine sweep_g_zed (ivmu, g)
!
!     use zgrid, only: nzgrid, delzed, ntubes
!     use extended_zgrid, only: neigen, nsegments, nzed_segment
!     use extended_zgrid, only: map_to_extended_zgrid
!     use extended_zgrid, only: map_from_extended_zgrid
!     use kt_grids, only: naky
!     use kt_grids, only: zonal_mode
!     use stella_layouts, only: vmu_lo
!     use stella_layouts, only: iv_idx, is_idx
!     use run_parameters, only: zed_upwind, time_upwind
!
!     implicit none
!
!     integer, intent (in) :: ivmu
!     complex, dimension (:,:,-nzgrid:,:), intent (in out) :: g
!
!     integer :: iv, is
!     integer :: iky, ie, it
!     integer :: ulim, sgn
!     integer :: iz, izext, iz1, iz2
!     real :: fac1, fac2
!     complex, dimension (:), allocatable :: gext
!
!     iv = iv_idx (vmu_lo,ivmu)
!     is = is_idx (vmu_lo,ivmu)
!     sgn = stream_sign(iv)
!     ! will sweep to right (positive vpa) or left (negative vpa)
!     ! and solve for g on the extended z-grid
!     do iky = 1, naky
!        if (zonal_mode(iky)) then
!           call sweep_zed_zonal (iv, is, sgn, g(iky,:,:,:))
!        else
!           do it = 1, ntubes
!              do ie = 1, neigen(iky)
!                 allocate (gext(nsegments(ie,iky)*nzed_segment+1))
!                 ! get g on extended domain in zed
!                 call map_to_extended_zgrid (it, ie, iky, g(iky,:,:,:), gext, ulim)
!                 if (sgn < 0) then
!                    iz1 = 1 ; iz2 = ulim
!                 else
!                    iz1 = ulim ; iz2 = 1
!                 end if
!                 izext = iz1 ; iz = sgn*nzgrid
!                 fac1 = 1.0+zed_upwind+sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
!                 gext(izext) = gext(izext)*2.0/fac1
!                 do izext = iz1-sgn, iz2, -sgn
!                    if (iz == -sgn*nzgrid) then
!                       iz = sgn*nzgrid-sgn
!                    else
!                       iz = iz - sgn
!                    end if
!                    fac1 = 1.0+zed_upwind+sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
!                    fac2 = 1.0-zed_upwind-sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
!                    gext(izext) = (-gext(izext+sgn)*fac2 + 2.0*gext(izext))/fac1
!                 end do
!                 ! extract g from extended domain in zed
!                 call map_from_extended_zgrid (it, ie, iky, gext, g(iky,:,:,:))
!                 deallocate (gext)
!              end do
!           end do
!        end if
!     end do
!
!   end subroutine sweep_g_zed
!
!   subroutine invert_parstream_response (phi, apar, bpar)
!
!     use linear_solve, only: lu_back_substitution
!     use zgrid, only: nzgrid, ntubes, nztot
!     use extended_zgrid, only: neigen
!     use extended_zgrid, only: nsegments
!     use extended_zgrid, only: nzed_segment
!     use extended_zgrid, only: map_to_extended_zgrid
!     use extended_zgrid, only: map_from_extended_zgrid
!     use extended_zgrid, only: ikxmod
!     use kt_grids, only: naky, zonal_mode
!     use fields_arrays, only: response_matrix
!
!     implicit none
!
!     complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar, bpar
!
!     integer :: iky, ie, it, ulim, nz_segment
!     integer :: ikx
!     complex, dimension (:), allocatable :: fields
!
!     ! need to put the fields into extended zed grid
!     do iky = 1, naky
!        ! avoid double counting of periodic endpoints for zonal modes
!        if (zonal_mode(iky)) then
!           allocate(fields(3*(nztot-1)))
!           do it = 1, ntubes
!              do ie = 1, neigen(iky)
!                 ikx = ikxmod(1,ie,iky)
!                 ! Because we want to avoid double-counting the endpoints,
!                 ! populate fields array ignoring the endpoint in
!                 ! phi, apar, bpar.
!                 fields(:nztot-1) = phi(iky,ikx,:nzgrid-1,it)
!                 fields(nztot:2*nztot-2) = apar(iky,ikx,:nzgrid-1,it)
!                 fields(2*nztot-1:3*nztot-3) = bpar(iky,ikx,:nzgrid-1,it)
!                 call lu_back_substitution (response_matrix(iky)%eigen(ie)%zloc, &
!                      response_matrix(iky)%eigen(ie)%idx, fields)
!                 phi(iky,ikx,:nzgrid-1,it) = fields(:nztot-1)
!                 phi(iky,ikx,nzgrid,it) = phi(iky,ikx,-nzgrid,it)
!                 apar(iky,ikx,:nzgrid-1,it) = fields(nztot:2*nztot-2)
!                 apar(iky,ikx,nzgrid,it) = apar(iky,ikx,-nzgrid,it)
!                 bpar(iky,ikx,:nzgrid-1,it) = fields(2*nztot-1:3*nztot-3)
!                 bpar(iky,ikx,nzgrid,it) = bpar(iky,ikx,-nzgrid,it)
!              end do
!           end do
!           deallocate(fields)
!        else
!           do it = 1, ntubes
!              do ie = 1, neigen(iky)
!                 ! solve response_matrix*phi^{n+1} = phi_{inh}^{n+1}
!                 nz_segment = (nsegments(ie,iky)*nzed_segment+1)
!                 allocate (fields(3*nz_segment))
!                 call map_to_extended_zgrid (it, ie, iky, phi(iky,:,:,:), fields(:nz_segment), ulim)
!                 call map_to_extended_zgrid (it, ie, iky, apar(iky,:,:,:), fields(nz_segment+1:2*nz_segment), ulim)
!                 call map_to_extended_zgrid (it, ie, iky, bpar(iky,:,:,:), fields(2*nz_segment+1:3*nz_segment), ulim)
!                 call lu_back_substitution (response_matrix(iky)%eigen(ie)%zloc, &
!                      response_matrix(iky)%eigen(ie)%idx, fields)
!                 call map_from_extended_zgrid (it, ie, iky, fields(:nz_segment), phi(iky,:,:,:))
!                 call map_from_extended_zgrid (it, ie, iky, fields(nz_segment+1:2*nz_segment), apar(iky,:,:,:))
!                 call map_from_extended_zgrid (it, ie, iky, fields(2*nz_segment+1:3*nz_segment), bpar(iky,:,:,:))
!                 deallocate (fields)
!              end do
!           end do
!        end if
!     end do
!
!   end subroutine invert_parstream_response
!
!   subroutine finish_parallel_streaming
!
!     use run_parameters, only: stream_implicit, driftkinetic_implicit
!
!     implicit none
!
!     if (allocated(stream)) deallocate (stream)
!     if (allocated(stream_c)) deallocate (stream_c)
!     if (allocated(stream_sign)) deallocate (stream_sign)
!     if (allocated(gradpar_c)) deallocate (gradpar_c)
!     if (allocated(stream_rad_var1)) deallocate (stream_rad_var1)
!     if (allocated(stream_rad_var2)) deallocate (stream_rad_var2)
!
!     if (stream_implicit .or. driftkinetic_implicit) call finish_invert_stream_operator
!
!     parallel_streaming_initialized = .false.
!
!   end subroutine finish_parallel_streaming
!
!   subroutine finish_invert_stream_operator
!
!     implicit none
!
!     if (allocated(stream_tri_a1)) then
!        deallocate (stream_tri_a1)
!        deallocate (stream_tri_a2)
!        deallocate (stream_tri_b1)
!        deallocate (stream_tri_b2)
!        deallocate (stream_tri_c1)
!        deallocate (stream_tri_c2)
!     end if
!
!   end subroutine finish_invert_stream_operator

  !   subroutine init_parallel_streaming
  !
  !     use finite_differences, only: fd3pt
  !     use stella_time, only: code_dt
  !     use stella_layouts, only: vmu_lo
  !     use stella_layouts, only: iv_idx, imu_idx, is_idx
  !     use species, only: spec, nspec, pfac
  !     use vpamu_grids, only: nvpa
  !     use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
  !     use vpamu_grids, only: vperp2, vpa, mu
  !     use kt_grids, only: nalpha
  !     use zgrid, only: nzgrid, nztot
  !     use stella_geometry, only: gradpar, dgradpardrho, dBdrho, gfac
  !     use run_parameters, only: stream_implicit, driftkinetic_implicit
  !     use physics_flags, only: include_parallel_streaming, radial_variation
  !
  !     implicit none
  !
  !     integer :: iv, imu, is, ivmu, ia
  !
  !     real, dimension (:), allocatable :: energy
  !
  !     if (parallel_streaming_initialized) return
  !     parallel_streaming_initialized = .true.
  !
  !     if (.not.allocated(stream)) allocate (stream(-nzgrid:nzgrid,nvpa,nspec)) ; stream = 0.
  !     if (.not.allocated(stream_sign)) allocate (stream_sign(nvpa)) ; stream_sign = 0
  !
  !     ! sign of stream corresponds to appearing on RHS of GK equation
  !     ! i.e., this is the factor multiplying dg/dz on RHS of equation
  !     if (include_parallel_streaming) then
  !        stream = -code_dt*spread(spread(spec%stm_psi0,1,nztot),2,nvpa) &
  !             * spread(spread(vpa,1,nztot)*spread(gradpar,2,nvpa),3,nspec)
  !     else
  !        stream = 0.0
  !     end if
  !
  !
  !     if(radial_variation) then
  !       allocate (energy(-nzgrid:nzgrid))
  !
  !       if(.not.allocated(stream_rad_var1)) then
  !         allocate(stream_rad_var1(-nzgrid:nzgrid,nvpa,nspec))
  !       endif
  !       if(.not.allocated(stream_rad_var2)) then
  !         allocate(stream_rad_var2(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
  !         stream_rad_var2 = 0.0
  !       endif
  !       ia=1
  !       stream_rad_var1 = -code_dt*spread(spread(spec%stm_psi0,1,nztot),2,nvpa) &
  !             * gfac*spread(spread(vpa,1,nztot)*spread(dgradpardrho,2,nvpa),3,nspec)
  !       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
  !         is  = is_idx(vmu_lo,ivmu)
  !         imu = imu_idx(vmu_lo,ivmu)
  !         iv  = iv_idx(vmu_lo,ivmu)
  !         energy = (vpa(iv)**2 + vperp2(ia,:,imu))*(spec(is)%temp_psi0/spec(is)%temp)
  !         stream_rad_var2(ia,:,ivmu) = &
  !                 +code_dt*spec(is)%stm_psi0*vpa(iv)*gradpar &
  !                 *spec(is)%zt*maxwell_vpa(iv,is)*maxwell_mu(ia,:,imu,is)*maxwell_fac(is) &
  !                 *(  pfac*(spec(is)%fprim + spec(is)%tprim*(energy-2.5)) &
  !                   + gfac*2*mu(imu)*dBdrho)
  !       enddo
  !       deallocate (energy)
  !     endif
  !
  !     ! stream_sign set to +/- 1 depending on the sign of the parallel streaming term.
  !     ! NB: stream_sign = -1 corresponds to positive advection velocity
  !     do iv = 1, nvpa
  !        stream_sign(iv) = int(sign(1.0,stream(0,iv,1)))
  !     end do
  !     ! vpa = 0 is special case
  ! !    stream_sign(0) = 0
  !
  !     ! RJD: We no longer need to worry about this because it's handled by
  !     ! implicit_z.
  !     if (stream_implicit .or. driftkinetic_implicit) then
  !        call init_invert_stream_operator
  !        if (.not.allocated(stream_c)) allocate (stream_c(-nzgrid:nzgrid,nvpa,nspec))
  !        stream_c = stream
  !        do is = 1, nspec
  !           do iv = 1, nvpa
  !              call center_zed (iv, stream_c(:,iv,is))
  !           end do
  !        end do
  !        if (.not.allocated(gradpar_c)) allocate (gradpar_c(-nzgrid:nzgrid,-1:1))
  !        gradpar_c = spread(gradpar,2,3)
  !        ! get gradpar centred in zed for negative vpa (affects upwinding)
  !        call center_zed(1,gradpar_c(:,-stream_sign(1)))
  !        ! get gradpar centred in zed for positive vpa (affects upwinding)
  !        call center_zed(nvpa,gradpar_c(:,-stream_sign(nvpa)))
  !        stream = stream_c
  !     end if
  !
  !   end subroutine init_parallel_streaming

end module parallel_streaming
