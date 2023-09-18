module parallel_streaming

   implicit none

   public :: init_parallel_streaming, finish_parallel_streaming
   public :: advance_parallel_streaming_explicit
   public :: add_parallel_streaming_radial_variation
   public :: stream_tridiagonal_solve
   public :: parallel_streaming_initialized
   public :: stream, stream_c, stream_sign, gradpar_c
   public :: time_parallel_streaming
   public :: stream_rad_var1
   public :: stream_rad_var2
   public :: center_zed, get_dzed
   public :: get_zed_derivative_extended_domain

   private

   interface center_zed
      module procedure center_zed_segment_real
      module procedure center_zed_segment_complex
      module procedure center_zed_extended
   end interface center_zed

   logical :: parallel_streaming_initialized = .false.

   integer, dimension(:), allocatable :: stream_sign
   real, dimension(:, :, :, :), allocatable :: stream
   real, dimension(:, :, :), allocatable :: stream_c
   real, dimension(:, :, :), allocatable :: stream_rad_var1
   real, dimension(:, :, :), allocatable :: stream_rad_var2
   real, dimension(:, :), allocatable :: stream_tri_a1, stream_tri_a2
   real, dimension(:, :), allocatable :: stream_tri_b1, stream_tri_b2
   real, dimension(:, :), allocatable :: stream_tri_c1, stream_tri_c2
   real, dimension(:, :), allocatable :: gradpar_c

   real, dimension(2, 3) :: time_parallel_streaming = 0.

contains

   subroutine init_parallel_streaming

      use finite_differences, only: fd3pt
      use stella_time, only: code_dt
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use species, only: spec, nspec, pfac
      use vpamu_grids, only: nvpa, nvpa
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: vperp2, vpa, mu
      use kt_grids, only: nalpha
      use zgrid, only: nzgrid, nztot
      use stella_geometry, only: gradpar, dgradpardrho, dBdrho, gfac, b_dot_grad_z
      use run_parameters, only: stream_implicit, driftkinetic_implicit
      use physics_flags, only: include_parallel_streaming, radial_variation

      implicit none

      integer :: iv, imu, is, ivmu
      integer :: ia, iz

      real, dimension(:), allocatable :: energy

      if (parallel_streaming_initialized) return
      parallel_streaming_initialized = .true.

      if (.not. allocated(stream)) allocate (stream(nalpha, -nzgrid:nzgrid, nvpa, nspec)); stream = 0.
      if (.not. allocated(stream_sign)) allocate (stream_sign(nvpa)); stream_sign = 0

      ! sign of stream corresponds to appearing on RHS of GK equation
      ! i.e., this is the factor multiplying dg/dz on RHS of equation
      if (include_parallel_streaming) then
         do iv = 1, nvpa
            do iz = -nzgrid, nzgrid
               do ia = 1, nalpha
                  stream(ia, iz, iv, :) = -code_dt * b_dot_grad_z(ia, iz) * vpa(iv) * spec%stm_psi0
               end do
            end do
         end do
      else
         stream = 0.0
      end if

      if (radial_variation) then
         allocate (energy(-nzgrid:nzgrid))

         if (.not. allocated(stream_rad_var1)) then
            allocate (stream_rad_var1(-nzgrid:nzgrid, nvpa, nspec))
         end if
         if (.not. allocated(stream_rad_var2)) then
            allocate (stream_rad_var2(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            stream_rad_var2 = 0.0
         end if
         ia = 1
         stream_rad_var1 = -code_dt * spread(spread(spec%stm_psi0, 1, nztot), 2, nvpa) &
                           * gfac * spread(spread(vpa, 1, nztot) * spread(dgradpardrho, 2, nvpa), 3, nspec)
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            energy = (vpa(iv)**2 + vperp2(ia, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
            stream_rad_var2(ia, :, ivmu) = &
               +code_dt * spec(is)%stm_psi0 * vpa(iv) * gradpar &
               * spec(is)%zt * maxwell_vpa(iv, is) * maxwell_mu(ia, :, imu, is) * maxwell_fac(is) &
               * (pfac * (spec(is)%fprim + spec(is)%tprim * (energy - 2.5)) &
                  + gfac * 2 * mu(imu) * dBdrho)
         end do
         deallocate (energy)
      end if

      !> stream_sign set to +/- 1 depending on the sign of the parallel streaming term.
      !> NB: stream_sign = -1 corresponds to positive advection velocity
      !> only need to consider ia=1, iz=0 and is=1 because alpha, z and species dependences
      !> do not lead to change in sign of the streaming pre-factor
      do iv = 1, nvpa
         stream_sign(iv) = int(sign(1.0, stream(1, 0, iv, 1)))
      end do

      if (stream_implicit .or. driftkinetic_implicit) then
         call init_invert_stream_operator
         if (.not. allocated(stream_c)) allocate (stream_c(-nzgrid:nzgrid, nvpa, nspec))
         stream_c = stream(1, :, :, :)
         do is = 1, nspec
            do iv = 1, nvpa
               call center_zed(iv, stream_c(:, iv, is), -nzgrid)
            end do
         end do
         if (.not. allocated(gradpar_c)) allocate (gradpar_c(-nzgrid:nzgrid, -1:1))
         gradpar_c = spread(gradpar, 2, 3)
         !> get gradpar centred in zed for negative vpa (affects upwinding)
         call center_zed(1, gradpar_c(:, -stream_sign(1)), -nzgrid)
         !> get gradpar centred in zed for positive vpa (affects upwinding)
         call center_zed(nvpa, gradpar_c(:, -stream_sign(nvpa)), -nzgrid)
         stream = spread(stream_c, 1, nalpha)
      end if

   end subroutine init_parallel_streaming

   subroutine init_invert_stream_operator

      use zgrid, only: delzed
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: nsegments
      use run_parameters, only: zed_upwind_plus, zed_upwind_minus, time_upwind_plus

      implicit none

      integer :: nz, nseg_max

      nz = maxval(iz_up - iz_low)
      nseg_max = maxval(nsegments)

      if (.not. allocated(stream_tri_a1)) then
         allocate (stream_tri_a1(nz * nseg_max + 1, -1:1)); stream_tri_a1 = 0.
         allocate (stream_tri_a2(nz * nseg_max + 1, -1:1)); stream_tri_a2 = 0.
         allocate (stream_tri_b1(nz * nseg_max + 1, -1:1)); stream_tri_b1 = 1.
         allocate (stream_tri_b2(nz * nseg_max + 1, -1:1)); stream_tri_b2 = 0.
         allocate (stream_tri_c1(nz * nseg_max + 1, -1:1)); stream_tri_c1 = 0.
         allocate (stream_tri_c2(nz * nseg_max + 1, -1:1)); stream_tri_c2 = 0.
      end if

      ! corresponds to sign of stream term positive on RHS of equation
      ! i.e., negative parallel advection speed
      ! NB: assumes equal spacing in zed
      stream_tri_b1(:, 1) = zed_upwind_plus
      stream_tri_b2(:, 1) = -1.0 / delzed(0)
      stream_tri_c1(:nz * nseg_max, 1) = zed_upwind_minus
      stream_tri_c2(:nz * nseg_max, 1) = 1.0 / delzed(0)
      ! corresponds to sign of stream term negative on RHS of equation
      ! NB: assumes equal spacing in zed
      stream_tri_b1(:, -1) = zed_upwind_plus
      stream_tri_b2(:, -1) = 1.0 / delzed(0)
      stream_tri_a1(2:, -1) = zed_upwind_minus
      stream_tri_a2(2:, -1) = -1.0 / delzed(0)

      stream_tri_a2 = time_upwind_plus * stream_tri_a2
      stream_tri_b2 = time_upwind_plus * stream_tri_b2
      stream_tri_c2 = time_upwind_plus * stream_tri_c2

   end subroutine init_invert_stream_operator

   subroutine advance_parallel_streaming_explicit(g, phi, gout)

      use mp, only: proc0
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use job_manage, only: time_message
      use stella_transforms, only: transform_ky2y
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, naky_all, nakx, ikx_max, ny
      use kt_grids, only: swap_kxky
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use species, only: spec
      use physics_flags, only: full_flux_surface
      use gyro_averages, only: gyro_average
      use run_parameters, only: driftkinetic_implicit

      use fields, only: advance_fields, fields_updated
      use fields_arrays, only: apar
      use gyro_averages, only: j0_ffs

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      integer :: ivmu, iv, imu, is, ia, iz, it
      complex, dimension(:, :, :, :), allocatable :: g0, dgphi_dz
      complex, dimension(:, :, :, :), allocatable :: g0y, g1y
      complex, dimension(:, :), allocatable :: g0_swap

      !> if flux tube simulation parallel streaming stays in ky,kx,z space with ky,kx,z local
      !> if full flux surface (flux annulus), will need to calculate in y space

      !> start the timer for the parallel streaming part of the time advance
      if (proc0) call time_message(.false., time_parallel_streaming(:, 1), ' Stream advance')

      !> allocate arrays needed for intermmediate calculations
      allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dgphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
      !> if simulating a full flux surface, will also need version of the above arrays
      !> that is Fourier transformed to y-space
      if (full_flux_surface) then
         allocate (g0_swap(naky_all, ikx_max))
         allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes))
         allocate (g1y(ny, ikx_max, -nzgrid:nzgrid, ntubes))
      end if

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         !> get (iv,imu,is) indices corresponding to ivmu super-index
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         !> obtain <phi> (or <phi>-phi if driftkinetic_implicit=T)
         if (full_flux_surface) then
            call gyro_average(phi, g0(:, :, :, :), j0_ffs(:, :, :, ivmu))
         else
            call gyro_average(phi, ivmu, g0(:, :, :, :))
         end if

         if (driftkinetic_implicit) g0(:, :, :, :) = g0(:, :, :, :) - phi

         !> get d<phi>/dz, with z the parallel coordinate and store in dgphi_dz
         !> note that this should be a centered difference to avoid numerical
         !> unpleasantness to do with inexact cancellations in later velocity integration
         !> see appendix of the stella JCP 2019 for details
         call get_dgdz_centered(g0, ivmu, dgphi_dz)

         !> if driftkinetic_implicit=T, then only want to treat vpar . grad (<phi>-phi)*F0 term explicitly;
         !> in this case, zero out dg/dz term (or d(g/F)/dz for full-flux-surface)
         if (driftkinetic_implicit) then
            g0 = 0.
         else
            !> compute dg/dz in k-space and store in g0
            call get_dgdz(g(:, :, :, :, ivmu), ivmu, g0)
            !> if simulating a full flux surface, need to obtain the contribution from parallel streaming
            !> in y-space, so FFT d(g/F)/dz from ky to y
            if (full_flux_surface) then
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
                     call swap_kxky(g0(:, :, iz, it), g0_swap)
                     call transform_ky2y(g0_swap, g0y(:, :, iz, it))
                  end do
               end do
            end if
            ! ! if simulating a full flux surface, must calculate F * d/dz (g/F) rather than dg/dz
            ! ! since F=F(y) in this case, avoid multiple Fourier transforms by applying chain rule
            ! ! to z derivative: F * d/dz (g/F) = dg/dz - g * d ln F / dz = dg/dz + g * mu/T * dB/dz
            ! if (full_flux_surface) then
            !    ! transform g and dg/dz from ky to y space and store in g0y and g1y, respectively
            !    g1y = g(:,:,:,:,ivmu)
            !    do it = 1, ntubes
            !       do iz = -nzgrid, nzgrid
            !          call transform_ky2y (g1y(:,:,iz,it), g0y(:,:,iz,it))
            !          ! no longer need g1y so re-use as FFT of dg/dz (g0)
            !          call transform_ky2y (g0(:,:,iz,it), g1y(:,:,iz,it))
            !       end do
            !    end do
            !    ! overwrite g0y with dg/dz + g * mu/T * dB/dz
            !    g0y = g1y + 2.0*mu(imu)*spread(spread(dBdzed,2,nakx),4,ntubes) * g0y
            !    ! g1y no longer needed so can over-write with d<phi>/dz below
            ! end if
         end if

         if (full_flux_surface) then
            !> transform d<phi>/dz (fully explicit) or d(<phi>-phi)/dz (if driftkinetic_implicit)
            !> from kalpha (ky) to alpha (y) space and store in g1y
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call swap_kxky(dgphi_dz(:, :, iz, it), g0_swap)
                  call transform_ky2y(g0_swap, g1y(:, :, iz, it))
               end do
            end do
            ! ! over-write g0y with F * d/dz (g/F) + ZeF/T * d<phi>/dz (or <phi>-phi for driftkinetic_implicit).
            g0y(:,:,:,:) = g0y(:,:,:,:) + g1y(:,:,:,:)*spec(is)%zt*maxwell_fac(is) &
                 * maxwell_vpa(iv,is)*spread(spread(maxwell_mu(:,:,imu,is),2,nakx),4,ntubes)*maxwell_fac(is)

            !> multiply d(g/F)/dz and d<phi>/dz terms with vpa*(b . grad z) and add to source (RHS of GK equation)
            call add_stream_term_ffs(g0y, ivmu, gout(:, :, :, :, ivmu))
         else
            ia = 1
            g0(:, :, :, :) = g0(:, :, :, :) + dgphi_dz(:, :, :, :) * spec(is)%zt * maxwell_fac(is) &
                                * maxwell_vpa(iv, is) * spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx), 4, ntubes)

            ! multiply dg/dz with vpa*(b . grad z) and add to source (RHS of GK equation)
            call add_stream_term(g0, ivmu, gout(:, :, :, :, ivmu))
         end if

      end do

      !> deallocate intermediate arrays used in this subroutine
      deallocate (g0, dgphi_dz)
      if (full_flux_surface) deallocate (g0y, g1y, g0_swap)

      !> finish timing the subroutine
      if (proc0) call time_message(.false., time_parallel_streaming(:, 1), ' Stream advance')

   end subroutine advance_parallel_streaming_explicit

   subroutine add_parallel_streaming_radial_variation(g, gout, rhs)

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

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      !the next input/output is for quasineutrality and gyroaveraging corrections
      !that go directly in the RHS (since they don't require further FFTs)
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: rhs

      integer :: ivmu, iv, imu, is, it, ia, iz

      complex, dimension(:, :, :, :), allocatable :: g0, g1, g2, g3
      complex, dimension(:, :), allocatable :: g0k

      allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g1(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g2(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g3(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g0k(naky, nakx)); g0k = 0

      ! parallel streaming stays in ky,kx,z space with ky,kx,z local

      ia = 1
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! obtain <phi>
         ! get d<phi>/dz, with z the parallel coordinate and store in g1
         call gyro_average(phi, ivmu, g0)
         call get_dgdz_centered(g0, ivmu, g1)

         ! get variation in gyroaveraging and store in g2
         call get_dgdz_centered(phi_corr_GA(:, :, :, :, ivmu), ivmu, g2)

         ! get variation in quasineutrality and store in g3
         call gyro_average(phi_corr_QN, ivmu, g0)
         call get_dgdz_centered(g0, ivmu, g3)

         call get_dgdz(g(:, :, :, :, ivmu), ivmu, g0)

         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid

    !!#1 - variation in gradpar
               g0k = g0(:, :, iz, it) &
                     + g1(:, :, iz, it) * spec(is)%zt * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)

               g0k = g0k * stream_rad_var1(iz, iv, is)

    !!#2 - variation in F_s/T_s
               g0k = g0k + g1(:, :, iz, it) * stream_rad_var2(ia, iz, ivmu)

               gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + g0k

    !!#3 - variation in the gyroaveraging and quasineutrality of phi
    !!     These variations already have the linear part calculated, so
    !!     ad it into the rhs directly
               g0k = spec(is)%zt * stream(1, iz, iv, is) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is) &
                     * (g2(:, :, iz, it) + g3(:, :, iz, it))

               rhs(:, :, iz, it, ivmu) = rhs(:, :, iz, it, ivmu) + g0k

            end do
         end do
      end do
      deallocate (g0, g1, g2, g3, g0k)

   end subroutine add_parallel_streaming_radial_variation

   subroutine get_dgdz(g, ivmu, dgdz)

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

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdz
      integer, intent(in) :: ivmu

      integer :: iseg, ie, it, iky, iv
      complex, dimension(2) :: gleft, gright

      ! FLAG -- assuming delta zed is equally spaced below!
      iv = iv_idx(vmu_lo, ivmu)
      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! first fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g(:, :, :, :), gleft, gright)
                  ! now get dg/dz
                  call third_order_upwind_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
                                              g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                                              delzed(0), stream_sign(iv), gleft, gright, periodic(iky), &
                                              dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
      end do

   end subroutine get_dgdz

   subroutine get_dgdz_centered(g, ivmu, dgdz)

      use finite_differences, only: second_order_centered_zed
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

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdz
      integer, intent(in) :: ivmu

      integer :: iseg, ie, iky, iv, it
      complex, dimension(2) :: gleft, gright
      ! FLAG -- assuming delta zed is equally spaced below!
      iv = iv_idx(vmu_lo, ivmu)
      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! first fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g(:, :, :, :), gleft, gright)
                  ! now get dg/dz
                  call second_order_centered_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
                                                 g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                                                 delzed(0), stream_sign(iv), gleft, gright, periodic(iky), &
                                                 dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
      end do
   end subroutine get_dgdz_centered

! subroutine get_dgdz_variable (g, ivmu, dgdz)

!    use finite_differences, only: fd_variable_upwinding_zed
!    use stella_layouts, only: vmu_lo
!    use stella_layouts, only: iv_idx
!    use zgrid, only: nzgrid, delzed, ntubes
!    use extended_zgrid, only: neigen, nsegments
!    use extended_zgrid, only: iz_low, iz_up
!    use extended_zgrid, only: ikxmod
!    use extended_zgrid, only: fill_zed_ghost_zones
!    use extended_zgrid, only: periodic
!    use run_parameters, only: zed_upwind
!    use kt_grids, only: naky

!    implicit none

!    complex, dimension (:,:,-nzgrid:,:), intent (in) :: g
!    complex, dimension (:,:,-nzgrid:,:), intent (out) :: dgdz
!    integer, intent (in) :: ivmu

!    integer :: iseg, ie, iky, iv, it
!    complex, dimension (2) :: gleft, gright
!    ! FLAG -- assuming delta zed is equally spaced below!
!     iv = iv_idx(vmu_lo,ivmu)
!     do iky = 1, naky
!       do it = 1, ntubes
!         do ie = 1, neigen(iky)
!           do iseg = 1, nsegments(ie,iky)
!              ! first fill in ghost zones at boundaries in g(z)
!              call fill_zed_ghost_zones (it, iseg, ie, iky, g(:,:,:,:), gleft, gright)
   ! now get dg/dz
!              call fd_variable_upwinding_zed (iz_low(iseg), iseg, nsegments(ie,iky), &
!                   g(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it), &
!                   delzed(0), stream_sign(iv), zed_upwind,gleft, gright, periodic(iky), &
!                   dgdz(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it))
!           end do
!         end do
!       enddo
!     end do
! end subroutine get_dgdz_variable

   subroutine add_stream_term(g, ivmu, src)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, is_idx
      use zgrid, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: src
      integer, intent(in) :: ivmu

      integer :: iz, iv, is

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      do iz = -nzgrid, nzgrid
         src(:, :, iz, :) = src(:, :, iz, :) + stream(1, iz, iv, is) * g(:, :, iz, :)
      end do

   end subroutine add_stream_term

   subroutine add_stream_term_ffs(g, ivmu, src)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, is_idx
      use zgrid, only: nzgrid
      use kt_grids, only: ny

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: src
      integer, intent(in) :: ivmu

      integer :: iz, iy, iv, is

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      do iz = -nzgrid, nzgrid
         do iy = 1, ny
            src(iy, :, iz, :) = src(iy, :, iz, :) + stream(iy, iz, iv, is) * g(iy, :, iz, :)
         end do
      end do

   end subroutine add_stream_term_ffs

   subroutine stream_tridiagonal_solve(iky, ie, iv, is, g)

      use finite_differences, only: tridag
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: nsegments
      use extended_zgrid, only: nzed_segment

      implicit none

      integer, intent(in) :: iky, ie, iv, is
      complex, dimension(:), intent(in out) :: g

      integer :: iseg, llim, ulim, n
      integer :: nz, nseg_max, sgn, n_ext
      integer :: ia
      real, dimension(:), allocatable :: a, b, c

      ia = 1

      ! avoid double-counting at boundaries between 2pi segments
      nz = nzed_segment
      nseg_max = nsegments(ie, iky)
      sgn = stream_sign(iv)

      n_ext = nseg_max * nz + 1
      allocate (a(n_ext))
      allocate (b(n_ext))
      allocate (c(n_ext))

      iseg = 1
      llim = 1; ulim = nz + 1
      a(llim:ulim) = stream_tri_a1(llim:ulim, sgn) &
                     - stream(ia, iz_low(iseg):iz_up(iseg), iv, is) * stream_tri_a2(llim:ulim, sgn)
      b(llim:ulim) = stream_tri_b1(llim:ulim, sgn) &
                     - stream(ia, iz_low(iseg):iz_up(iseg), iv, is) * stream_tri_b2(llim:ulim, sgn)
      c(llim:ulim) = stream_tri_c1(llim:ulim, sgn) &
                     - stream(ia, iz_low(iseg):iz_up(iseg), iv, is) * stream_tri_c2(llim:ulim, sgn)

      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            llim = ulim + 1
            ulim = llim + nz - 1
            a(llim:ulim) = stream_tri_a1(llim:ulim, sgn) &
                           - stream(ia, iz_low(iseg) + 1:iz_up(iseg), iv, is) * stream_tri_a2(llim:ulim, sgn)
            b(llim:ulim) = stream_tri_b1(llim:ulim, sgn) &
                           - stream(ia, iz_low(iseg) + 1:iz_up(iseg), iv, is) * stream_tri_b2(llim:ulim, sgn)
            c(llim:ulim) = stream_tri_c1(llim:ulim, sgn) &
                           - stream(ia, iz_low(iseg) + 1:iz_up(iseg), iv, is) * stream_tri_c2(llim:ulim, sgn)
         end do
      end if
      n = size(stream_tri_a1, 1)
      a(ulim) = stream_tri_a1(n, sgn) - stream(ia, iz_up(nsegments(ie, iky)), iv, is) * stream_tri_a2(n, sgn)
      b(ulim) = stream_tri_b1(n, sgn) - stream(ia, iz_up(nsegments(ie, iky)), iv, is) * stream_tri_b2(n, sgn)
      c(ulim) = 0. ! this line should not be necessary, as c(ulim) should not be accessed by tridag
      call tridag(1, a(:ulim), b(:ulim), c(:ulim), g)

      deallocate (a, b, c)

   end subroutine stream_tridiagonal_solve

   subroutine get_dzed(iv, g, dgdz)

      use finite_differences, only: fd_cell_centres_zed
      use kt_grids, only: naky
      use zgrid, only: nzgrid, delzed, ntubes
      use extended_zgrid, only: neigen, nsegments
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: fill_zed_ghost_zones

      implicit none

      integer, intent(in) :: iv
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdz

      integer :: iky, ie, iseg, it
      complex, dimension(2) :: gleft, gright

      do it = 1, ntubes
         do iky = 1, naky
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! first fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g, gleft, gright)
                  ! get finite difference approximation for dg/dz at cell centres
                  ! iv > nvgrid corresponds to positive vpa, iv <= nvgrid to negative vpa
                  call fd_cell_centres_zed(iz_low(iseg), &
                                           g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                                           delzed(0), stream_sign(iv), gleft(2), gright(1), &
                                           dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
      end do

   end subroutine get_dzed

   subroutine get_zed_derivative_extended_domain(iv, f, f_left, f_right, df_dz)

      use zgrid, only: delzed
      use finite_differences, only: fd_cell_centres_zed

      implicit none

      integer, intent(in) :: iv
      complex, dimension(:), intent(in) :: f
      complex, intent(in) :: f_left, f_right
      complex, dimension(:), intent(out) :: df_dz

      call fd_cell_centres_zed(1, f, delzed(0), stream_sign(iv), f_left, f_right, df_dz)

   end subroutine get_zed_derivative_extended_domain

   subroutine center_zed_extended(iv, g)

      use finite_differences, only: cell_centres_zed
      use kt_grids, only: naky, nakx
      use zgrid, only: nzgrid, ntubes
      use extended_zgrid, only: neigen, nsegments
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: fill_zed_ghost_zones
      use run_parameters, only: zed_upwind

      implicit none

      integer, intent(in) :: iv
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: g

      integer :: iky, ie, iseg, it
      complex, dimension(2) :: gleft, gright
      complex, dimension(:, :, :), allocatable :: gc

      allocate (gc(nakx, -nzgrid:nzgrid, ntubes))

      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! first fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g, gleft, gright)
                  ! get cell centres values
                  call cell_centres_zed(iz_low(iseg), &
                                        g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                                        zed_upwind, stream_sign(iv), gleft(2), gright(1), &
                                        gc(ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
         g(iky, :, :, :) = gc
      end do

      deallocate (gc)

   end subroutine center_zed_extended

   !> center_zed_segment_real takes as arguments the vpa index (iv)
   !> the z-depenendent real function f, and the starting iz index for the array f (llim),
   !> and overwrites f with the cell-centered version
   ! it is assumed that any function passed to this subroutine is periodic
   subroutine center_zed_segment_real(iv, f, llim)

      use zgrid, only: nzgrid
      use run_parameters, only: zed_upwind_plus, zed_upwind_minus

      integer, intent(in) :: iv, llim
      real, dimension(llim:), intent(in out) :: f

      integer :: ulim

      ulim = llim + size(f) - 1

      if (stream_sign(iv) > 0) then
         f(:ulim - 1) = zed_upwind_plus * f(:ulim - 1) + zed_upwind_minus * f(llim + 1:)
         f(ulim) = f(llim)
      else
         f(llim + 1:) = zed_upwind_minus * f(:ulim - 1) + zed_upwind_plus * f(llim + 1:)
         f(llim) = f(ulim)
      end if

   end subroutine center_zed_segment_real

   !> center_zed_segment_complex takes as arguments the vpa index (iv)
   !> the z-depenendent conplex function f, and the starting iz index for the array f (llim),
   !> and overwrites f with the cell-centered version;
   subroutine center_zed_segment_complex(iv, f, llim, periodic)

      use zgrid, only: nzgrid
      use run_parameters, only: zupwnd_p => zed_upwind_plus
      use run_parameters, only: zupwnd_m => zed_upwind_minus

      integer, intent(in) :: iv, llim
      complex, dimension(llim:), intent(in out) :: f
      logical, intent(in) :: periodic

      integer :: ulim

      ulim = llim + size(f) - 1

      ! stream_sign > 0 means negative advection velocity
      if (stream_sign(iv) > 0) then
         f(:ulim - 1) = zupwnd_p * f(:ulim - 1) + zupwnd_m * f(llim + 1:)
         if (periodic) then
            f(ulim) = f(llim)
         else
            f(ulim) = zupwnd_p * f(ulim)
         end if
      else
         f(llim + 1:) = zupwnd_m * f(:ulim - 1) + zupwnd_p * f(llim + 1:)
         if (periodic) then
            f(llim) = f(ulim)
         else
            f(llim) = zupwnd_p * f(llim)
         end if
      end if

   end subroutine center_zed_segment_complex

   subroutine center_zed_midpoint(iv, g)

      use zgrid, only: nzgrid

      integer, intent(in) :: iv
      real, dimension(-nzgrid:), intent(in out) :: g

      if (stream_sign(iv) > 0) then
         g(:nzgrid - 1) = 0.5 * (g(:nzgrid - 1) + g(-nzgrid + 1:))
         g(nzgrid) = g(-nzgrid)
      else
         g(-nzgrid + 1:) = 0.5 * (g(:nzgrid - 1) + g(-nzgrid + 1:))
         g(-nzgrid) = g(nzgrid)
      end if

   end subroutine center_zed_midpoint

   subroutine finish_parallel_streaming

      use run_parameters, only: stream_implicit, driftkinetic_implicit

      implicit none

      if (allocated(stream)) deallocate (stream)
      if (allocated(stream_c)) deallocate (stream_c)
      if (allocated(stream_sign)) deallocate (stream_sign)
      if (allocated(gradpar_c)) deallocate (gradpar_c)
      if (allocated(stream_rad_var1)) deallocate (stream_rad_var1)
      if (allocated(stream_rad_var2)) deallocate (stream_rad_var2)

      if (stream_implicit .or. driftkinetic_implicit) call finish_invert_stream_operator

      parallel_streaming_initialized = .false.

   end subroutine finish_parallel_streaming

   subroutine finish_invert_stream_operator

      implicit none

      if (allocated(stream_tri_a1)) then
         deallocate (stream_tri_a1)
         deallocate (stream_tri_a2)
         deallocate (stream_tri_b1)
         deallocate (stream_tri_b2)
         deallocate (stream_tri_c1)
         deallocate (stream_tri_c2)
      end if

   end subroutine finish_invert_stream_operator

end module parallel_streaming
