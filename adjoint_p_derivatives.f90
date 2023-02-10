module adjoint_p_derivatives

   implicit none

   public :: perturb_p
   public :: allocate_unpert
   public :: deallocate_p_derivatives
   public :: get_denominator
   public :: lagrangian_integrals

   public :: integrate_unpert

   private

contains

   subroutine allocate_unpert

      use stella_layouts, only: vmu_lo
      use kt_grids, only: nakx, naky
      use zgrid, only: nzgrid, ntubes

      use adjoint_distfn_arrays, only: g_unpert
      use adjoint_field_arrays, only: q_unpert
      use adjoint_field_arrays, only: derivative, denominator

      implicit none

      if (.not. allocated(g_unpert)) then
         allocate (g_unpert(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         g_unpert = 0.0
      end if

      if (.not. allocated(q_unpert)) then
         allocate (q_unpert(naky, nakx, -nzgrid:nzgrid, ntubes))
         q_unpert = 0.0
      end if

      if (.not. allocated(denominator)) allocate (denominator(naky, nakx)); denominator = 0.0
      if (.not. allocated(derivative)) allocate (derivative(naky, nakx)); derivative = 0.0

   end subroutine allocate_unpert

   subroutine perturb_p(g_term, q_term, unpert)

      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx

      use physics_flags, only: include_parallel_streaming
      use physics_flags, only: include_mirror

      use millerlocal, only: del

      use adjoint_distfn_arrays, only: g_unpert
      use adjoint_field_arrays, only: q_unpert
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: g_term
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: q_term
      logical, optional, intent(in) :: unpert

    !! Adjoint - g_term is the gyrokinetic equation contribution to Langrangian integrand
    !! Adjoint - q_term is the quasineutrality contribution to Lagrangian integrand

      g_term = 0.0
      q_term = 0.0

      if (include_mirror) call get_mirror_term(g_term)
      if (include_parallel_streaming) call get_prll_strm_term(g_term)
      call get_wstar_term(g_term)
      call get_wdrift_term(g_term)
      call get_omega_term(g_term)

      call get_quasin_term(q_term)

   end subroutine perturb_p

  !! Mirror Term
   subroutine get_mirror_term(gout)

      use vpamu_grids, only: nmu
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx, nalpha

      use species, only: nspec
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx

      use physics_flags, only: full_flux_surface
      use mirror_terms, only: get_dgdvpa_explicit
      use dist_fn_arrays, only: gvmu
      use adjoint_distfn_arrays, only: gsave

      use redistribute, only: gather, scatter
      use dist_redistribute, only: kxkyz2vmu

      use mirror_terms, only: mirror
      use stella_time, only: code_dt

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
      real, dimension(:, :, :, :), allocatable :: mirror_ad
      complex, dimension(:, :, :, :, :), allocatable :: dgdvpa

      allocate (mirror_ad(nalpha, -nzgrid:nzgrid, nmu, nspec))
      allocate (dgdvpa(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    !! Adjoint - get dg/dvp
      call scatter(kxkyz2vmu, gsave, gvmu)
      call get_dgdvpa_explicit(gvmu)
      call gather(kxkyz2vmu, gvmu, dgdvpa)

    !! Adjoint - get mirror term coefficent; v_th*mu* b_hat . grad B
    !!           with correct sign
      mirror_ad = -mirror / code_dt

    !! Adjoint - add mirror term to lagrangian integrand
      call add_mirror_term(dgdvpa, mirror_ad, gout)

      deallocate (mirror_ad)
      deallocate (dgdvpa)

   end subroutine get_mirror_term

   subroutine add_mirror_term(g, mirror_ad, mirr_out)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: imu_idx, is_idx
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(:, -nzgrid:, :, :), intent(in) :: mirror_ad
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: mirr_out

      integer :: imu, is, ivmu

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         mirr_out(:, :, :, :, ivmu) = mirr_out(:, :, :, :, ivmu) + &
                                      spread(spread(spread(mirror_ad(1, :, imu, is), 1, naky), 2, nakx), 4, ntubes) * g(:, :, :, :, ivmu)
      end do

   end subroutine add_mirror_term

  !! Parallel Streaming Term
   subroutine get_prll_strm_term(gout)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use vpamu_grids, only: nvpa
      use kt_grids, only: nakx, naky, nalpha
      use zgrid, only: nzgrid, ntubes
      use species, only: spec, nspec

      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac

      use gyro_averages, only: gyro_average
      use parallel_streaming, only: get_dgdz_centered, stream
      use stella_time, only: code_dt

      use adjoint_distfn_arrays, only: gsave
      use adjoint_field_arrays, only: phi_save

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      real, dimension(:, :, :, :), allocatable :: stream_ad
      complex, dimension(:, :, :, :), allocatable :: par
      complex, dimension(:, :, :, :), allocatable :: gyrophi
      complex, dimension(:, :, :, :), allocatable :: dphidz, dgdz

      integer :: is, imu, iv, ia, ivmu

      allocate (stream_ad(nalpha, -nzgrid:nzgrid, nvpa, nspec)); stream_ad = 0.0
      allocate (gyrophi(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dphidz(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dgdz(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (par(naky, nakx, -nzgrid:nzgrid, ntubes))

    !! Adjoint - get stream co-efficient; v_th* vpa* b_hat . grad z
    !!           with correct sign
      stream_ad = -stream / code_dt

      ia = 1
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         gyrophi = 0.

         call gyro_average(phi_save, ivmu, gyrophi)

         dphidz = 0.0
         dgdz = 0.0

       !! Adjoint - get d/dzed of g and J0*phi
         call get_dgdz_centered(gyrophi, ivmu, dphidz)
         call get_dgdz_centered(gsave(:, :, :, :, ivmu), ivmu, dgdz)

       !! Adjoint - get part of parallel streaming term without co-efficient
       !!         - dg/dz + Z/T * F_0 * d(phi *J_0)/dz
         par = 0.

         par(:, :, :, :) = dgdz(:, :, :, :) + dphidz(:, :, :, :) * spec(is)%zt &
                           * maxwell_vpa(iv, is) * spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx), 4, ntubes) &
                           * maxwell_fac(is)

       !! Adjoint - Multiply by co-efficient and add to Lagrangian integrand
         call add_stream_term(par, stream_ad(1, :, :, :), ivmu, gout(:, :, :, :, ivmu))
      end do

      deallocate (stream_ad)
      deallocate (gyrophi)
      deallocate (dphidz)
      deallocate (dgdz)
      deallocate (par)

   end subroutine get_prll_strm_term

   subroutine add_stream_term(g, stream, ivmu, stream_out)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, is_idx
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      real, dimension(:, :, :), intent(in) :: stream
      integer, intent(in) :: ivmu

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: stream_out

      integer :: iv, is

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      stream_out(:, :, :, :) = stream_out(:, :, :, :) + &
                               spread(spread(spread(stream(:, iv, is), 1, naky), 2, nakx), 4, ntubes) * g(:, :, :, :)

   end subroutine add_stream_term

  !! Omega* Term
   subroutine get_wstar_term(gout)

      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx, nalpha
      use stella_layouts, only: vmu_lo

      use dist_fn_arrays, only: wstar
      use fields_arrays, only: apar

      use adjoint_field_arrays, only: phi_save
      use stella_time, only: code_dt
      use fields, only: get_dchidy

      use stella_layouts, only: iv_idx, imu_idx, is_idx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
      complex, dimension(:, :, :, :, :), allocatable :: dchidy
      real, dimension(:, :, :), allocatable :: wstar_ad

      allocate (dchidy(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      allocate (wstar_ad(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      dchidy = 0.
    !! Adjoint - get omega* frequency with correct sign
      wstar_ad = -wstar / code_dt

    !! Adjoint - i*iky*J0*phi
      apar = 0.0
      call get_dchidy(phi_save, apar, dchidy)
      apar = 0.0

    !! Adjoint - add omega* term to Lagrangian integrand
      call add_wstar_term(dchidy, wstar_ad, gout)

      deallocate (wstar_ad)
      deallocate (dchidy)

   end subroutine get_wstar_term

   subroutine add_wstar_term(g, wstar, wstar_out)

      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(in) :: wstar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: wstar_out
      integer :: ivmu

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         wstar_out(:, :, :, :, ivmu) = wstar_out(:, :, :, :, ivmu) + &
                                       spread(spread(spread(wstar(1, :, ivmu), 1, naky), 2, nakx), 4, ntubes) * g(:, :, :, :, ivmu)
      end do

   end subroutine add_wstar_term

  !! Omega drift terms
   subroutine get_wdrift_term(gout)

      use dist_fn_arrays, only: wdriftx_g, wdrifty_g
      use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi

      use adjoint_field_arrays, only: phi_save
      use adjoint_distfn_arrays, only: gsave

      use gyro_averages, only: gyro_average
      use stella_time, only: code_dt

      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx, nalpha
      use constants, only: zi

      use time_advance, only: get_dgdy, get_dgdx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
      complex, dimension(:, :, :, :, :), allocatable :: gyro
      complex, dimension(:, :, :, :, :), allocatable :: gy, gx
      real, dimension(:, :, :), allocatable :: drifty_g, driftx_g
      real, dimension(:, :, :), allocatable :: drifty_phi, driftx_phi

      integer :: ivmu

      allocate (gyro(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); gyro = 0.0
      allocate (gy(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); gy = 0.0
      allocate (gx(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); gx = 0.0

      allocate (drifty_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); drifty_g = 0.0
      allocate (driftx_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); driftx_g = 0.0
      allocate (drifty_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); drifty_phi = 0.0
      allocate (driftx_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); driftx_phi = 0.0

      drifty_g = -wdrifty_g / code_dt
      driftx_g = -wdriftx_g / code_dt
      drifty_phi = -wdrifty_phi / code_dt
      driftx_phi = -wdriftx_phi / code_dt

      call get_dgdy(gsave, gy)
      call get_dgdx(gsave, gx)

    !! Adjoint - add g drift terms to lagrangian integrand
      call add_drift(gy, drifty_g(1, :, :), gout)
      call add_drift(gx, driftx_g(1, :, :), gout)

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         call gyro_average(phi_save, ivmu, gyro(:, :, :, :, ivmu))
      end do

      call get_dgdy(gyro, gy)
      call get_dgdx(gyro, gx)

    !! Adjoint - add phi drift terms to lagrangian integrand
      call add_drift(gy, drifty_phi(1, :, :), gout)
      call add_drift(gx, driftx_phi(1, :, :), gout)

      deallocate (gy, gx)
      deallocate (drifty_g, driftx_g)
      deallocate (drifty_phi, driftx_phi)
      deallocate (gyro)

   end subroutine get_wdrift_term

   subroutine add_drift(g, wdrift, src)

      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(-nzgrid:, vmu_lo%llim_proc:), intent(in) :: wdrift
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: src

      integer :: ivmu

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         src(:, :, :, :, ivmu) = src(:, :, :, :, ivmu) + &
                                 spread(spread(spread(wdrift(:, ivmu), 1, naky), 2, nakx), 4, ntubes) * g(:, :, :, :, ivmu)
      end do

   end subroutine add_drift

  !! Omega Term - gamma_0 * g
   subroutine get_omega_term(gout)

      use adjoint_distfn_arrays, only: gsave
      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes, nztot
      use adjoint_field_arrays, only: omega_g

      use kt_grids, only: naky, nakx
      use stella_layouts, only: iv_idx, imu_idx, is_idx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
      integer :: ivmu, iv, imu, is

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         gout(:, :, :, :, ivmu) = gout(:, :, :, :, ivmu) + &
                                  spread(spread(omega_g(:, :), 3, nztot), 4, ntubes) * gsave(:, :, :, :, ivmu)

      end do

   end subroutine get_omega_term

  !! Quasineutrality Term
   subroutine get_quasin_term(qout)

      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx

      use adjoint_distfn_arrays, only: gsave
      use adjoint_field_arrays, only: phi_save
      use fields_arrays, only: gamtot, apar

      use fields, only: advance_fields, fields_updated

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(out) :: qout
      complex, dimension(:, :, :, :), allocatable :: phi_from_g

      allocate (phi_from_g(naky, nakx, -nzgrid:nzgrid, ntubes))
      phi_from_g = 0.0

    !! Adjoint - this is int(J_0 * gsave)/gamtot
      fields_updated = .false.
      call advance_fields(gsave, phi_from_g, apar, dist='gbar')

    !! Adjoint - Note that phi_save has a minus sign as gamtot is defined with
    !!           a minus sign as it appears in Quasineutrality
      qout = (phi_from_g - phi_save) * spread(gamtot, 4, ntubes)

      deallocate (phi_from_g)

   end subroutine get_quasin_term

  !! Get int(lam*g)
   subroutine get_denominator

      use adjoint_distfn_arrays, only: lam_save, gsave
      use adjoint_field_arrays, only: denominator

      use vpamu_grids, only: integrate_species
      use species, only: nspec
      use kt_grids, only: naky, nakx
      use zgrid, only: nzgrid, ntubes

      use mp, only: sum_allreduce
      use volume_averages, only: fieldline_average

      use stella_geometry, only: gradpar
      implicit none

      complex, dimension(:, :, :, :), allocatable :: int_tubes
      real, dimension(:), allocatable :: wgts
      integer :: iz, it

      real :: unp

      allocate (int_tubes(naky, nakx, -nzgrid:nzgrid, ntubes)); int_tubes = 0.0
      allocate (wgts(nspec))
      wgts = 1.
      denominator = 0.

      unp = sum(1 / (gradpar(:)), dim=1)

    !! Adjoint - Velocity integral
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
!          call integrate_species(gsave(:,:,iz,it,:)*lam_save(:,:,iz,it,:),&
            call integrate_species(conjg(lam_save(:, :, iz, it, :)) * gsave(:, :, iz, it, :), &
                                   iz, wgts, int_tubes(:, :, iz, it), it, reduce_in=.false.)
         end do
      end do
      call sum_allreduce(int_tubes)

    !! Adjoint - zed integral
      call fieldline_average(int_tubes, denominator)

      denominator = denominator / unp

      deallocate (int_tubes)
      deallocate (wgts)

   end subroutine get_denominator

  !! Integrals in numerator to get Lagrangian
   subroutine lagrangian_integrals(gin, qin, gout)
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use species, only: nspec

      use adjoint_distfn_arrays, only: lam_save
      use adjoint_field_arrays, only: chi_save

      use stella_layouts, only: vmu_lo

      use volume_averages, only: fieldline_average
      use vpamu_grids, only: integrate_species
      use mp, only: sum_allreduce

      use stella_geometry, only: gradpar
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: qin
      complex, dimension(:, :), allocatable, intent(out) :: gout

      complex, dimension(:, :), allocatable :: qout
      complex, dimension(:, :, :, :), allocatable :: int_tubes
      real, dimension(:), allocatable :: wgts

      real :: unp
      integer :: it, iz

      allocate (gout(naky, nakx)); gout = 0.0
      allocate (int_tubes(naky, nakx, -nzgrid:nzgrid, ntubes)); int_tubes = 0.
      allocate (qout(naky, nakx)); qout = 0.0

      allocate (wgts(nspec)); wgts = 1.0

      unp = sum(1 / (gradpar(:)), dim=1)

      int_tubes = 0.0
    !! Adjoint - Velocity Integral
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
!          call integrate_species(gin(:,:,iz,it,:)*lam_save(:,:,iz,it,:),&
            call integrate_species(conjg(lam_save(:, :, iz, it, :)) * gin(:, :, iz, it, :), &
                                   iz, wgts, int_tubes(:, :, iz, it), it, reduce_in=.False.)
         end do
      end do
      call sum_allreduce(int_tubes)

      call fieldline_average(int_tubes, gout)
      call fieldline_average(conjg(chi_save) * qin, qout)
      !    call fieldline_average (qin*chi_save, qout)

      gout = (gout + qout) / unp

      deallocate (int_tubes)
      deallocate (wgts)
      deallocate (qout)
   end subroutine lagrangian_integrals

   subroutine integrate_unpert(check_out)

      use stella_geometry, only: gradpar

      use adjoint_distfn_arrays, only: gsave
!    use adjoint_field_arrays, only: q_unpert
      use adjoint_distfn_arrays, only: lam_save
      use adjoint_field_arrays, only: chi_save

      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use species, only: nspec

      use volume_averages, only: fieldline_average
      use vpamu_grids, only: integrate_species
      use mp, only: sum_allreduce

      implicit none

      complex, intent(out) :: check_out
      complex, dimension(:, :), allocatable :: check

      complex, dimension(:, :), allocatable :: qout
      real, dimension(:), allocatable :: wgts
      complex, dimension(:, :, :, :), allocatable :: dum

      !real :: unp
      integer :: iz, it

      allocate (check(naky, nakx)); check = 0.0
      allocate (dum(naky, nakx, -nzgrid:nzgrid, ntubes)); dum = 0.
      allocate (wgts(nspec)); wgts = 1.0
      allocate (qout(naky, nakx)); qout = 0.0

      !unp = sum(1/(gradpar(:)), dim = 1)

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            call integrate_species(conjg(lam_save(:, :, iz, it, :)) * gsave(:, :, iz, it, :), &
                                   iz, wgts, dum(:, :, iz, it), it, reduce_in=.False.)
         end do
      end do
      call sum_allreduce(dum)
      call fieldline_average(dum, check)

      check_out = check(1, 1)!/unp

      deallocate (qout)
      deallocate (dum)
      deallocate (wgts)

   end subroutine integrate_unpert

   ! !! Integrals in numerator to get Lagrangian
  !! Deallocate Arrays
   subroutine deallocate_p_derivatives

      use adjoint_field_arrays, only: derivative, denominator
      use adjoint_distfn_arrays, only: g_unpert
      use adjoint_field_arrays, only: q_unpert

      implicit none

      deallocate (g_unpert, q_unpert)
      deallocate (denominator)
      deallocate (derivative)

   end subroutine deallocate_p_derivatives

end module adjoint_p_derivatives
