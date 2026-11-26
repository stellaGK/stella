!###############################################################################
!            Evaluate GVEC geometry directly, via the GVEC Fortran library
!###############################################################################
!
! This is the direct alternative to geometry_gvec_read_netCDF_file: instead of
! reading a geometry netCDF written by gvec_to_stella.py, it links against GVEC
! and evaluates the equilibrium in-place, so a GVEC state file can be given to
! stella exactly as a VMEC wout file is.
!
! It produces the SAME outputs as read_gvec_geometry, so everything downstream
! (the arc-length grid, the (alpha,psi_t) -> (x,y) transform, the boundary
! conditions) is shared between the two paths, and the two can be compared
! directly on the same equilibrium.
!
! Requires stella to be built with GVEC support (STELLA_ENABLE_GVEC=ON); without
! it the routine below aborts with an explanatory message and stella still
! builds and runs for every other geometry option.
!
! Notes on the GVEC interface
! ---------------------------
! * The interface used is MODgvec_gvec_to_gene, GVEC's existing library API for
!   gyrokinetic codes.
! * <FillStrings> must be called on the GVEC parameter (.ini) file BEFORE the
!   state is read whenever the equilibrium uses a general frame mapping, since
!   those h-maps read extra parameters of their own (the axisNB frame reads
!   hmap_ncfile). Omitting it is why the interface is documented as supporting
!   only the cylindrical frame; with it, the general frames work.
! * gvec_to_gene_scalars returns the TOTAL toroidal flux (it multiplies by 2*pi
!   internally), unlike the pygvec profile, which is per radian.
!
!###############################################################################
module gvec_evaluate_geometry

   implicit none

   public :: evaluate_gvec_geometry

   private

   logical, parameter :: debug = .false.

contains

   !============================================================================
   ! Evaluate all geometry arrays stella needs from a GVEC state file.
   ! Output arguments match read_gvec_geometry exactly.
   !============================================================================
   subroutine evaluate_gvec_geometry( &
      ! Input
      state_file, parameter_file, sfl_coord, sfl_factor, &
      nalpha, nzgrid_file, verbose, &
      s_in, alpha0_in, zeta_center_in, nfield_periods_in, &
      ! Scalar outputs
      s, safety_factor_q, shat, L_reference, B_reference, nfp_out, &
      sign_toroidal_flux, nfp_sim_out, dzdzeta, alpha0_out, zeta_center_out, &
      ! 1-D outputs
      alpha, zeta, &
      ! 2-D outputs (nalpha, -nzgrid_file:nzgrid_file)
      theta, grho, bmag, b_dot_gradzeta, &
      B_sub_theta, B_sub_zeta, &
      grad_alpha_grad_alpha, grad_alpha_grad_psit, grad_psit_grad_psit, &
      gds23_psitalpha, gds24_psitalpha, gds25_psitalpha, gds26_psitalpha, &
      gbdrift_alpha, B_times_gradB_dot_gradpsit, &
      gbdrift_kappa_alpha, B_times_kappa_dot_gradpsit, &
      psit_displacement_fac, &
      gradzeta_gradalpha_R2overB2, gradzeta_gradpsit_R2overB2, &
      b_dot_gradzeta_RR, &
      ierr)

#ifdef GVEC
      use MODgvec_ReadInTools,  only: FillStrings
      use MODgvec_gvec_to_gene, only: init_gvec_to_gene, gvec_to_gene_scalars
      use MODgvec_gvec_to_gene, only: gvec_to_gene_profile, gvec_to_gene_coords
      use MODgvec_gvec_to_gene, only: gvec_to_gene_metrics, finalize_gvec_to_gene
#endif
      use mp,        only: mp_abort
      use constants, only: pi
      use debug_flags, only: print_extra_info_to_terminal

      implicit none

      ! Input
      character(*), intent(in) :: state_file, parameter_file
      integer,      intent(in) :: sfl_coord, sfl_factor
      integer,      intent(in) :: nalpha, nzgrid_file
      logical,      intent(in) :: verbose
      real,         intent(in) :: s_in, alpha0_in, zeta_center_in, nfield_periods_in

      ! Scalar outputs
      integer, intent(out) :: sign_toroidal_flux, nfp_out, ierr, nfp_sim_out
      real,    intent(out) :: s, safety_factor_q, shat
      real,    intent(out) :: L_reference, B_reference, dzdzeta
      real,    intent(out) :: alpha0_out, zeta_center_out

      ! 1-D outputs
      real, dimension(:),             intent(out) :: alpha
      real, dimension(-nzgrid_file:), intent(out) :: zeta

      ! 2-D outputs
      real, dimension(:, -nzgrid_file:), intent(out) :: theta, grho, bmag
      real, dimension(:, -nzgrid_file:), intent(out) :: b_dot_gradzeta
      real, dimension(:, -nzgrid_file:), intent(out) :: B_sub_theta, B_sub_zeta
      real, dimension(:, -nzgrid_file:), intent(out) :: grad_alpha_grad_alpha
      real, dimension(:, -nzgrid_file:), intent(out) :: grad_alpha_grad_psit
      real, dimension(:, -nzgrid_file:), intent(out) :: grad_psit_grad_psit
      real, dimension(:, -nzgrid_file:), intent(out) :: gds23_psitalpha
      real, dimension(:, -nzgrid_file:), intent(out) :: gds24_psitalpha
      real, dimension(:, -nzgrid_file:), intent(out) :: gds25_psitalpha
      real, dimension(:, -nzgrid_file:), intent(out) :: gds26_psitalpha
      real, dimension(:, -nzgrid_file:), intent(out) :: gbdrift_alpha
      real, dimension(:, -nzgrid_file:), intent(out) :: B_times_gradB_dot_gradpsit
      real, dimension(:, -nzgrid_file:), intent(out) :: gbdrift_kappa_alpha
      real, dimension(:, -nzgrid_file:), intent(out) :: B_times_kappa_dot_gradpsit
      real, dimension(:, -nzgrid_file:), intent(out) :: psit_displacement_fac
      real, dimension(:, -nzgrid_file:), intent(out) :: gradzeta_gradalpha_R2overB2
      real, dimension(:, -nzgrid_file:), intent(out) :: gradzeta_gradpsit_R2overB2
      real, dimension(:, -nzgrid_file:), intent(out) :: b_dot_gradzeta_RR

#ifdef GVEC
      ! Local variables
      integer, parameter :: dp = kind(1.0d0)
      real(dp), parameter :: mu_0 = 4.0d-7 * 3.14159265358979323846d0

      integer  :: nz_file, ia, iz, nfp_gvec
      real(dp) :: Fa, minor_r, PhiPrime_edge, q_edge, major_R
      real(dp) :: rho0, q0, q_prime, pres, dp_drho
      real(dp) :: iota0, diota_drho, dp_ds
      real(dp) :: phi_edge_over_2pi, a_ref, B_ref
      real(dp) :: nfp_sim_real
      logical  :: flip_tor, flip_pol
      real(dp) :: alpha0_use, zeta_center_use

      real(dp), allocatable :: thet2d(:, :), zeta2d(:, :), theta_out2d(:, :)
      real(dp), allocatable :: cart(:, :, :)
      real(dp), allocatable :: grad_rho(:, :, :), grad_ts(:, :, :), grad_zt(:, :, :)
      real(dp), allocatable :: Bvec(:, :, :), grad_absB(:, :, :)
      real(dp), allocatable :: zeta1d(:)

      real(dp) :: gr(3), gts(3), gz(3), bb(3), gB(3)
      real(dp) :: grad_s_v(3), grad_psit(3), grad_alpha(3), BcGB(3), BcGs(3)
      real(dp) :: modB, Rcyl, zmc, bm, bm2
      real(dp) :: gag, gap, gpp, gza, gzp, gta, gtp
      real(dp) :: BcGBdGA, BcGBdGP, BcGsdGA, BcGsdGP
      real(dp) :: Btg_dGA, Btg_dGP
#endif

      !-------------------------------------------------------------------------

      ierr = 0

#ifndef GVEC
      write (*, '(A)') "  ERROR: stella was built without GVEC support, so a GVEC"
      write (*, '(A)') "         state file cannot be read directly."
      write (*, '(A)') "         Either rebuild with -DSTELLA_ENABLE_GVEC=ON, or use"
      write (*, '(A)') "         gvec_to_stella.py and set <gvec_filename> instead."
      ierr = 1
      call mp_abort('evaluate_gvec_geometry: stella was compiled without GVEC support')
      ! Silence unused-argument warnings in the no-GVEC build.
      alpha = 0.0; zeta = 0.0; theta = 0.0; grho = 0.0; bmag = 0.0
      b_dot_gradzeta = 0.0; B_sub_theta = 0.0; B_sub_zeta = 0.0
      grad_alpha_grad_alpha = 0.0; grad_alpha_grad_psit = 0.0
      grad_psit_grad_psit = 0.0
      gds23_psitalpha = 0.0; gds24_psitalpha = 0.0
      gds25_psitalpha = 0.0; gds26_psitalpha = 0.0
      gbdrift_alpha = 0.0; B_times_gradB_dot_gradpsit = 0.0
      gbdrift_kappa_alpha = 0.0; B_times_kappa_dot_gradpsit = 0.0
      psit_displacement_fac = 0.0
      gradzeta_gradalpha_R2overB2 = 0.0; gradzeta_gradpsit_R2overB2 = 0.0
      b_dot_gradzeta_RR = 0.0
      s = 0.0; safety_factor_q = 0.0; shat = 0.0
      L_reference = 0.0; B_reference = 0.0; dzdzeta = 1.0
      nfp_out = 1; nfp_sim_out = 1; sign_toroidal_flux = 1
      alpha0_out = 0.0; zeta_center_out = 0.0
#else

      nz_file = 2 * nzgrid_file + 1

      ! <alpha0> and <zeta_center> carry the sentinel -999 meaning "take it from
      ! the geometry file". There is no file on this path, so the sentinels must
      ! be resolved to their actual defaults; using them literally would place
      ! the field line and the toroidal domain at alpha = zeta = -999.
      alpha0_use      = alpha0_in
      zeta_center_use = zeta_center_in
      if (alpha0_use      < -900.0) alpha0_use      = 0.0d0
      if (zeta_center_use < -900.0) zeta_center_use = 0.0d0

      if (verbose .and. print_extra_info_to_terminal) then
         write (*, '(A)') " "
         write (*, '(A)') "######################################################"
         write (*, '(A)') "         GVEC GEOMETRY (direct, via libGVEC)"
         write (*, '(A)') "######################################################"
         write (*, '(A)') " "
         write (*, '(2A)') "  GVEC state file:     ", trim(state_file)
         write (*, '(2A)') "  GVEC parameter file: ", trim(parameter_file)
      end if

      ! Load the GVEC parameters, then the state
      ! A general frame mapping reads extra parameters of its own (the axisNB
      ! frame reads hmap_ncfile), so the parameter file must be processed first.
      ! Paths inside it are resolved relative to the working directory.
      if (len_trim(parameter_file) > 0) then
         call FillStrings(trim(parameter_file))
      end if

      call init_gvec_to_gene(trim(state_file), SFLcoord_in=sfl_coord, &
                             factorSFL_in=sfl_factor)

      ! Scalars
      ! NOTE: Fa is the TOTAL toroidal flux; gvec_to_gene_scalars multiplies the
      ! per-radian profile by 2*pi internally.
      call gvec_to_gene_scalars(Fa, minor_r, PhiPrime_edge, q_edge, nfp_gvec, major_R)

      phi_edge_over_2pi = Fa / (2.0d0 * pi)
      a_ref = minor_r
      B_ref = 2.0d0 * abs(phi_edge_over_2pi) / a_ref**2

      ! Surface quantities
      s    = s_in
      rho0 = sqrt(s_in)
      call gvec_to_gene_profile(rho0, q=q0, q_prime=q_prime, p=pres, p_prime=dp_drho)

      iota0 = 1.0d0 / q0
      ! q_prime is dq/drho, so diota/drho = -(dq/drho)/q^2
      diota_drho = -q_prime / (q0 * q0)
      ! stella works in s_vmec = rho^2, so d/ds = (1/2rho) d/drho
      dp_ds = dp_drho / (2.0d0 * rho0)

      ! Handedness
      ! The formulas of stella's VMEC geometry hold in the frame of the VMEC
      ! angles, whose Jacobian is negative, with the toroidal flux of the
      ! right-handed frame as psi_t. The logical coordinates of GVEC are
      ! right-handed, so the toroidal angle is reversed, as in gvec_to_stella.py.
      ! This reverses the sign of iota, but not the toroidal flux
      ! Phi' = B^zeta * Jacobian, so phi_edge_over_2pi keeps its sign.
      flip_tor = .true.
      flip_pol = .false.
      iota0 = -iota0
      diota_drho = -diota_drho

      safety_factor_q = 1.0d0 / iota0
      ! shat = -rho * (diota/drho) / iota
      shat = -rho0 * diota_drho / iota0

      L_reference = a_ref
      B_reference = B_ref
      sign_toroidal_flux = int(sign(1.0d0, phi_edge_over_2pi))
      nfp_out = nfp_gvec

      ! Field periods simulated
      if (nfield_periods_in > 0.0) then
         nfp_sim_real = nfield_periods_in
      else
         nfp_sim_real = real(nfp_gvec, dp)
      end if
      nfp_sim_out = nint(nfp_sim_real)
      dzdzeta = real(nfp_gvec, dp) / nfp_sim_real

      alpha0_out      = real(alpha0_use)
      zeta_center_out = real(zeta_center_use)

      if (verbose .and. print_extra_info_to_terminal) then
         write (*, '(A,F8.5)') "    s             = ", s
         write (*, '(A,F8.5)') "    q             = ", safety_factor_q
         write (*, '(A,F8.5)') "    shat          = ", shat
         write (*, '(A,F8.5)') "    L_reference   = ", L_reference
         write (*, '(A,F8.5)') "    B_reference   = ", B_reference
         write (*, '(A,I3)')   "    nfp           = ", nfp_out
         write (*, '(A,I3)')   "    sign_torflux  = ", sign_toroidal_flux
         write (*, '(A,2L2)')  "    flip tor/pol  = ", flip_tor, flip_pol
         write (*, '(A)') " "
      end if

      ! Field-line grid
      allocate (zeta1d(nz_file))
      do iz = 1, nz_file
         zeta1d(iz) = zeta_center_use + pi * real(iz - 1 - nzgrid_file, dp) &
                      * nfp_sim_real / (real(nfp_gvec, dp) * real(nzgrid_file, dp))
      end do
      do ia = 1, nalpha
         alpha(ia) = alpha0_use + real(ia - 1, dp) * 2.0d0 * pi / real(nalpha, dp)
      end do
      zeta(-nzgrid_file:nzgrid_file) = real(zeta1d(:))

      allocate (thet2d(nalpha, nz_file), zeta2d(nalpha, nz_file))
      allocate (theta_out2d(nalpha, nz_file))
      do ia = 1, nalpha
         do iz = 1, nz_file
            zeta2d(ia, iz) = zeta1d(iz)
            thet2d(ia, iz) = alpha(ia) + iota0 * zeta1d(iz)
         end do
      end do

      ! The grid above is in the FLIPPED (stella-facing) frame; GVEC must be
      ! evaluated in its own, so undo the flip on the sample points.
      if (flip_tor) zeta2d = -zeta2d
      if (flip_pol) thet2d = -thet2d

      ! Evaluate GVEC
      allocate (cart(3, nalpha, nz_file))
      allocate (grad_rho(3, nalpha, nz_file), grad_ts(3, nalpha, nz_file))
      allocate (grad_zt(3, nalpha, nz_file), Bvec(3, nalpha, nz_file))
      allocate (grad_absB(3, nalpha, nz_file))

      if (debug) write (*, *) 'gvec_evaluate_geometry::gvec_to_gene_coords'
      call gvec_to_gene_coords(nalpha, nz_file, rho0, thet2d, zeta2d, &
                               theta_out2d, cart)
      if (debug) write (*, *) 'gvec_evaluate_geometry::gvec_to_gene_metrics'
      call gvec_to_gene_metrics(nalpha, nz_file, rho0, thet2d, zeta2d, &
                                grad_rho, grad_ts, grad_zt, Bvec, grad_absB)

      ! Express the evaluated fields in the flipped frame. Quantities carrying a
      ! free theta or zeta index change sign; true vectors (position, grad rho,
      ! B, grad|B|) do not.
      if (flip_tor) grad_zt = -grad_zt
      if (flip_pol) then
         grad_ts     = -grad_ts
         theta_out2d = -theta_out2d
      end if

      ! Assemble stella's coefficients
      do ia = 1, nalpha
         do iz = 1, nz_file

            gr  = grad_rho(:, ia, iz)
            gts = grad_ts(:, ia, iz)
            gz  = grad_zt(:, ia, iz)
            bb  = Bvec(:, ia, iz)
            gB  = grad_absB(:, ia, iz)

            modB = sqrt(dot_product(bb, bb))
            Rcyl = sqrt(cart(1, ia, iz)**2 + cart(2, ia, iz)**2)
            zmc  = zeta1d(iz) - zeta_center_use

            ! grad s_vmec = 2 rho grad rho ;  grad psi_t = (Phi_edge/2pi) grad s
            grad_s_v  = 2.0d0 * rho0 * gr
            grad_psit = phi_edge_over_2pi * grad_s_v

            ! grad alpha = grad theta* - iota grad zeta - (zeta-zeta_c) (diota/drho) grad rho
            grad_alpha = gts - iota0 * gz - zmc * diota_drho * gr

            bm  = modB / B_ref
            bm2 = bm * bm

            bmag(ia, iz - nzgrid_file - 1) = real(bm)
            theta(ia, iz - nzgrid_file - 1) = real(theta_out2d(ia, iz))
            grho(ia, iz - nzgrid_file - 1) = real(a_ref * sqrt(dot_product(gr, gr)))

            ! b . grad zeta = (a/|B|) B . grad zeta
            b_dot_gradzeta(ia, iz - nzgrid_file - 1) = &
               real(a_ref * dot_product(bb, gz) / modB)

            gag = a_ref**2 * dot_product(grad_alpha, grad_alpha)
            gap = dot_product(grad_alpha, grad_psit) / B_ref
            gpp = dot_product(grad_psit, grad_psit) / (a_ref**2 * B_ref**2)

            gza = a_ref * dot_product(gz, grad_alpha)
            gzp = dot_product(gz, grad_psit) / (a_ref * B_ref)
            gta = a_ref * dot_product(gts, grad_alpha)
            gtp = dot_product(gts, grad_psit) / (a_ref * B_ref)

            grad_alpha_grad_alpha(ia, iz - nzgrid_file - 1) = real(gag)
            grad_alpha_grad_psit(ia, iz - nzgrid_file - 1)  = real(gap)
            grad_psit_grad_psit(ia, iz - nzgrid_file - 1)   = real(gpp)

            gds23_psitalpha(ia, iz - nzgrid_file - 1) = &
               real(-(gza * gap - gzp * gag) / bm2)
            gds24_psitalpha(ia, iz - nzgrid_file - 1) = &
               real(-(gza * gpp - gzp * gap) / bm2)
            gds25_psitalpha(ia, iz - nzgrid_file - 1) = &
               real(-real(sign_toroidal_flux, dp) * (gta * gap - gtp * gag) / bm2)
            gds26_psitalpha(ia, iz - nzgrid_file - 1) = &
               real(-0.5d0 * real(sign_toroidal_flux, dp) * (gta * gpp - gtp * gap) / bm2)

            ! Magnetic drifts
            BcGB(1) = bb(2) * gB(3) - bb(3) * gB(2)
            BcGB(2) = bb(3) * gB(1) - bb(1) * gB(3)
            BcGB(3) = bb(1) * gB(2) - bb(2) * gB(1)
            BcGs(1) = bb(2) * grad_s_v(3) - bb(3) * grad_s_v(2)
            BcGs(2) = bb(3) * grad_s_v(1) - bb(1) * grad_s_v(3)
            BcGs(3) = bb(1) * grad_s_v(2) - bb(2) * grad_s_v(1)

            BcGBdGA = dot_product(BcGB, grad_alpha)
            BcGBdGP = dot_product(BcGB, grad_psit)
            BcGsdGA = dot_product(BcGs, grad_alpha)
            BcGsdGP = dot_product(BcGs, grad_psit)

            Btg_dGA = a_ref**2 * B_ref / modB**3 * BcGBdGA
            Btg_dGP = BcGBdGP / modB**3

            gbdrift_alpha(ia, iz - nzgrid_file - 1) = real(Btg_dGA)
            B_times_gradB_dot_gradpsit(ia, iz - nzgrid_file - 1) = real(Btg_dGP)

            ! Curvature drift adds the pressure-gradient term
            gbdrift_kappa_alpha(ia, iz - nzgrid_file - 1) = &
               real(Btg_dGA + a_ref**2 * B_ref * mu_0 / modB**4 * dp_ds * BcGsdGA)
            B_times_kappa_dot_gradpsit(ia, iz - nzgrid_file - 1) = &
               real(Btg_dGP + mu_0 / modB**4 * dp_ds * BcGsdGP)

            ! R^2/B^2 terms
            gradzeta_gradalpha_R2overB2(ia, iz - nzgrid_file - 1) = &
               real(gza * Rcyl**2 / bm2)
            gradzeta_gradpsit_R2overB2(ia, iz - nzgrid_file - 1) = &
               real(gzp * Rcyl**2 / bm2)
            b_dot_gradzeta_RR(ia, iz - nzgrid_file - 1) = &
               real(a_ref * dot_product(bb, gz) / modB * Rcyl**2)

         end do
      end do

      ! Not available through the gvec_to_gene interface:
      !   B_sub_theta / B_sub_zeta   need the covariant basis vectors e_theta,
      !                              e_zeta, which the interface does not return.
      !   psit_displacement_fac      needs dLambda/drho.
      ! None of these is used by the GVEC branch of stella (B_sub_* are allocated
      ! and deallocated unused, and x_displacement_fac appears only in
      ! commented-out code in calculations_gyro_averages), so they are zeroed
      ! rather than partially filled, which would be silently wrong.
      B_sub_theta = 0.0
      B_sub_zeta = 0.0
      psit_displacement_fac = 0.0

      call finalize_gvec_to_gene()

      deallocate (thet2d, zeta2d, theta_out2d, cart, zeta1d)
      deallocate (grad_rho, grad_ts, grad_zt, Bvec, grad_absB)
#endif

   end subroutine evaluate_gvec_geometry

end module gvec_evaluate_geometry
