!###############################################################################
!################################ GVEC GEOMETRY ################################
!###############################################################################
!
! Geometry module for reading equilibria from GVEC via a pre-computed netCDF
! file produced by the Python script gvec_to_stella.py.
!
! Workflow:
!   1. Run GVEC to obtain an equilibrium state file.
!   2. Run:  python gvec_to_stella.py STATE_FILE --s 0.5 --nalpha 1 \
!                --nzgrid_file <nzgrid*factor>  -o my_gvec_geometry.nc
!   3. In the stella input file set:
!        &geometry_options
!          geometry_option = 'gvec'
!        /
!        &geometry_gvec
!          gvec_filename = 'my_gvec_geometry.nc'
!          z_grid_refinement_factor = <factor>   ! default 1
!        /
!
! When z_grid_refinement_factor > 1, the file must have been produced with
! --nzgrid_file = nzgrid * z_grid_refinement_factor so that arc-length
! coordinate (zed_equal_arc = .true.) can be computed accurately.
!
! Normalisation: all arrays follow stella's VMEC conventions.  See
! geometry_vmec_read_netCDF_file.f90 for full definitions.
!
! Written 2025 - mirrors geometry_vmec.f90.
!###############################################################################
module gvec_geometry

   use debug_flags, only: debug => geometry_debug
   use namelist_geometry, only: radial_coordinate_sgnpsitpsit
   use namelist_geometry, only: radial_coordinate_minuspsit
   use namelist_geometry, only: radial_coordinate_r

   implicit none

   public :: read_gvec_parameters
   public :: get_gvec_geometry

   public :: radial_coordinate_switch
   public :: radial_coordinate_sgnpsitpsit
   public :: radial_coordinate_minuspsit
   public :: radial_coordinate_r

   private

   ! Module-level parameters (read from namelist)
   integer  :: radial_coordinate_switch
   integer  :: z_grid_refinement_factor
   logical  :: verbose
   character(2000) :: gvec_filename
   character(2000) :: gvec_state_file, gvec_parameter_file
   integer :: sfl_coord, sfl_factor
   real     :: torflux_namelist, alpha0_namelist, zeta_center_namelist, nfield_periods_namelist

contains

   !============================================================================
   !===================== READ GVEC PARAMETERS =================================
   !============================================================================
   subroutine read_gvec_parameters

      use mp,     only: mp_abort
      use grids_z, only: zed_equal_arc
      use namelist_geometry, only: read_namelist_geometry_gvec

      implicit none

      call read_namelist_geometry_gvec( &
         z_grid_refinement_factor, radial_coordinate_switch, verbose, gvec_filename, &
         torflux_namelist, alpha0_namelist, zeta_center_namelist, nfield_periods_namelist, &
         gvec_state_file, gvec_parameter_file, sfl_coord, sfl_factor)

      ! Two ways to supply GVEC geometry:
      !   <gvec_state_file>  a GVEC state file, evaluated directly through the
      !                      GVEC library (needs STELLA_ENABLE_GVEC)
      !   <gvec_filename>    a geometry netCDF written by gvec_to_stella.py
      ! Using the direct path requires <torflux>, since the flux surface is no
      ! longer recorded in a pre-computed file.
      if (len_trim(gvec_state_file) > 0) then
         if (torflux_namelist <= 0.0) then
            call mp_abort('&geometry_gvec: <torflux> must be set (0 < torflux <= 1) '// &
                          'when reading a GVEC state file directly. Aborting.')
         end if
      end if

      if (.not. zed_equal_arc) then
         if (z_grid_refinement_factor > 1) then
            write (*, *) 'There is no reason to use z_grid_refinement_factor > 1 unless zed_equal_arc=T'
            write (*, *) 'Setting z_grid_refinement_factor = 1'
            z_grid_refinement_factor = 1
         end if
      end if

   end subroutine read_gvec_parameters


   !============================================================================
   !================== GET GEOMETRY ARRAYS FROM GVEC ==========================
   !============================================================================
   subroutine get_gvec_geometry( &
      nzgrid, nalpha, naky, surf, grho, bmag, &
      b_dot_gradz_avg, b_dot_gradz, &
      grad_alpha_grad_alpha, grad_alpha_grad_psit, grad_psit_grad_psit, &
      gds23_psitalpha, gds24_psitalpha, gds25_psitalpha, gds26_psitalpha, &
      gbdrift_alpha, B_times_gradB_dot_gradx_psi, &
      B_times_kappa_dot_grady_alpha, B_times_kappa_dot_gradx_psi, &
      gradzeta_gradpsit_R2overB2, gradzeta_gradalpha_R2overB2, b_dot_gradzeta_RR, &
      sign_torflux, theta_gvec, dzdzeta, L_reference, B_reference, &
      alpha, zeta, field_period_ratio, psit_displacement_fac)

      use mp,        only: mp_abort
      use constants, only: pi
      use splines,   only: geo_spline

      use grids_z,   only: zed_equal_arc, zed
      use grids_z,   only: get_total_arc_length, get_arc_length_grid

      use common_types,       only: flux_surface_type
      use file_utils,         only: open_output_file
      use debug_flags,        only: print_extra_info_to_terminal

      use geometry_gvec_read_netCDF_file, only: read_gvec_geometry
      use gvec_evaluate_geometry, only: evaluate_gvec_geometry

      implicit none

      ! Arguments
      integer, intent(in)  :: nzgrid, nalpha, naky
      integer, intent(out) :: sign_torflux
      type(flux_surface_type), intent(out) :: surf
      real,    intent(out) :: dzdzeta, L_reference, B_reference, field_period_ratio
      real, dimension(-nzgrid:),    intent(out) :: b_dot_gradz_avg
      real, dimension(:, -nzgrid:), intent(out) :: grho, bmag, b_dot_gradz
      real, dimension(:, -nzgrid:), intent(out) :: grad_alpha_grad_alpha
      real, dimension(:, -nzgrid:), intent(out) :: grad_alpha_grad_psit
      real, dimension(:, -nzgrid:), intent(out) :: grad_psit_grad_psit
      real, dimension(:, -nzgrid:), intent(out) :: gds23_psitalpha, gds24_psitalpha
      real, dimension(:, -nzgrid:), intent(out) :: gds25_psitalpha, gds26_psitalpha
      real, dimension(:, -nzgrid:), intent(out) :: gbdrift_alpha
      real, dimension(:, -nzgrid:), intent(out) :: B_times_gradB_dot_gradx_psi
      real, dimension(:, -nzgrid:), intent(out) :: B_times_kappa_dot_grady_alpha
      real, dimension(:, -nzgrid:), intent(out) :: B_times_kappa_dot_gradx_psi
      real, dimension(:, -nzgrid:), intent(out) :: gradzeta_gradpsit_R2overB2
      real, dimension(:, -nzgrid:), intent(out) :: gradzeta_gradalpha_R2overB2
      real, dimension(:, -nzgrid:), intent(out) :: b_dot_gradzeta_RR
      real, dimension(:, -nzgrid:), intent(out) :: theta_gvec, zeta, psit_displacement_fac
      real, dimension(:), intent(out) :: alpha

      ! Local variables
      integer :: ierr, ia, iz, nzgrid_file, nzgrid_file_expected, tmpunit
      integer :: zetamax_idx
      real    :: dzeta_file, zmin, nfp_real, nfp_sim_real
      integer :: nfp_out, nfp_sim

      ! Fine-grid versions of all geometry arrays (nzgrid_file = nzgrid * factor)
      real, allocatable :: zeta_file(:)
      real, allocatable :: alpha_file(:)
      real, allocatable :: arc_length(:, :)
      real, allocatable :: zed_domain_size(:)
      real, allocatable :: b_dot_gradzeta_file(:, :)
      real, allocatable :: b_dot_grad_arclength(:, :)
      real, allocatable :: b_dot_gradzeta_avg(:)
      real, allocatable :: b_dot_grad_arclength_avg(:)
      real, allocatable :: theta_file(:, :)
      real, allocatable :: grho_file(:, :)
      real, allocatable :: bmag_file(:, :)
      real, allocatable :: gag_file(:, :), gap_file(:, :), gpp_file(:, :)
      real, allocatable :: gds23_file(:, :), gds24_file(:, :)
      real, allocatable :: gds25_file(:, :), gds26_file(:, :)
      real, allocatable :: gbdrift_file(:, :)
      real, allocatable :: Btg_dGP_file(:, :)
      real, allocatable :: Btk_dGA_file(:, :), Btk_dGP_file(:, :)
      real, allocatable :: psit_disp_file(:, :)
      real, allocatable :: gza_R2oB2_file(:, :), gzp_R2oB2_file(:, :)
      real, allocatable :: bdgz_RR_file(:, :)
      real, allocatable :: B_sub_theta_file(:, :), B_sub_zeta_file(:, :)
      real    :: dzdzeta_file, s_file, alpha0_file, zeta_center_file
      integer :: nfp_sim_file

      if (debug) write (*, *) 'gvec_geometry::get_gvec_geometry'

      ! Grid sizes
      nzgrid_file = nzgrid * z_grid_refinement_factor

      ! Allocate fine-grid arrays
      allocate (zeta_file(-nzgrid_file:nzgrid_file));            zeta_file = 0.0
      allocate (alpha_file(nalpha));                             alpha_file = 0.0
      allocate (arc_length(nalpha, -nzgrid_file:nzgrid_file));  arc_length = 0.0
      allocate (zed_domain_size(nalpha));                        zed_domain_size = 0.0
      allocate (b_dot_gradzeta_file(nalpha, -nzgrid_file:nzgrid_file)); b_dot_gradzeta_file = 0.0
      allocate (b_dot_grad_arclength(nalpha, -nzgrid:nzgrid));  b_dot_grad_arclength = 0.0
      allocate (b_dot_gradzeta_avg(-nzgrid_file:nzgrid_file));  b_dot_gradzeta_avg = 0.0
      allocate (b_dot_grad_arclength_avg(-nzgrid:nzgrid));      b_dot_grad_arclength_avg = 0.0
      allocate (theta_file(nalpha, -nzgrid_file:nzgrid_file));  theta_file = 0.0
      allocate (grho_file(nalpha, -nzgrid_file:nzgrid_file));   grho_file = 0.0
      allocate (bmag_file(nalpha, -nzgrid_file:nzgrid_file));   bmag_file = 0.0
      allocate (gag_file(nalpha, -nzgrid_file:nzgrid_file));    gag_file = 0.0
      allocate (gap_file(nalpha, -nzgrid_file:nzgrid_file));    gap_file = 0.0
      allocate (gpp_file(nalpha, -nzgrid_file:nzgrid_file));    gpp_file = 0.0
      allocate (gds23_file(nalpha, -nzgrid_file:nzgrid_file));  gds23_file = 0.0
      allocate (gds24_file(nalpha, -nzgrid_file:nzgrid_file));  gds24_file = 0.0
      allocate (gds25_file(nalpha, -nzgrid_file:nzgrid_file));  gds25_file = 0.0
      allocate (gds26_file(nalpha, -nzgrid_file:nzgrid_file));  gds26_file = 0.0
      allocate (gbdrift_file(nalpha, -nzgrid_file:nzgrid_file)); gbdrift_file = 0.0
      allocate (Btg_dGP_file(nalpha, -nzgrid_file:nzgrid_file)); Btg_dGP_file = 0.0
      allocate (Btk_dGA_file(nalpha, -nzgrid_file:nzgrid_file)); Btk_dGA_file = 0.0
      allocate (Btk_dGP_file(nalpha, -nzgrid_file:nzgrid_file)); Btk_dGP_file = 0.0
      allocate (psit_disp_file(nalpha, -nzgrid_file:nzgrid_file)); psit_disp_file = 0.0
      allocate (gza_R2oB2_file(nalpha, -nzgrid_file:nzgrid_file)); gza_R2oB2_file = 0.0
      allocate (gzp_R2oB2_file(nalpha, -nzgrid_file:nzgrid_file)); gzp_R2oB2_file = 0.0
      allocate (bdgz_RR_file(nalpha, -nzgrid_file:nzgrid_file));   bdgz_RR_file = 0.0
      allocate (B_sub_theta_file(nalpha, -nzgrid_file:nzgrid_file)); B_sub_theta_file = 0.0
      allocate (B_sub_zeta_file(nalpha, -nzgrid_file:nzgrid_file));  B_sub_zeta_file = 0.0

      ! Obtain the geometry, either directly from GVEC or from a netCDF file
      ierr = 0
      if (len_trim(gvec_state_file) > 0) then

         if (debug) write (*, *) 'gvec_geometry::call evaluate_gvec_geometry'
         call evaluate_gvec_geometry( &
            gvec_state_file, gvec_parameter_file, sfl_coord, sfl_factor, &
            nalpha, nzgrid_file, verbose, &
            torflux_namelist, alpha0_namelist, zeta_center_namelist, &
            nfield_periods_namelist, &
            s_file, surf%qinp, surf%shat, L_reference, B_reference, nfp_out, &
            sign_torflux, nfp_sim_file, dzdzeta_file, alpha0_file, zeta_center_file, &
            alpha_file, zeta_file, &
            theta_file, grho_file, bmag_file, b_dot_gradzeta_file, &
            B_sub_theta_file, B_sub_zeta_file, &
            gag_file, gap_file, gpp_file, &
            gds23_file, gds24_file, gds25_file, gds26_file, &
            gbdrift_file, Btg_dGP_file, Btk_dGA_file, Btk_dGP_file, &
            psit_disp_file, gza_R2oB2_file, gzp_R2oB2_file, bdgz_RR_file, &
            ierr)

         if (ierr /= 0) call mp_abort('evaluate_gvec_geometry returned an error')

      else

      if (debug) write (*, *) 'gvec_geometry::call read_gvec_geometry'
      call read_gvec_geometry( &
         gvec_filename, nalpha, nzgrid_file, verbose, &
         s_file, surf%qinp, surf%shat, L_reference, B_reference, nfp_out, &
         sign_torflux, nfp_sim_file, dzdzeta_file, alpha0_file, zeta_center_file, &
         alpha_file, zeta_file, &
         theta_file, grho_file, bmag_file, b_dot_gradzeta_file, &
         B_sub_theta_file, B_sub_zeta_file, &
         gag_file, gap_file, gpp_file, &
         gds23_file, gds24_file, gds25_file, gds26_file, &
         gbdrift_file, Btg_dGP_file, Btk_dGA_file, Btk_dGP_file, &
         psit_disp_file, gza_R2oB2_file, gzp_R2oB2_file, bdgz_RR_file, &
         ierr)

      if (ierr /= 0) call mp_abort('read_gvec_geometry returned an error')

      end if

      ! Validate namelist parameters against the file
      ! If the user specified torflux / alpha0 / zeta_center / nfield_periods in
      ! &geometry_gvec, check that they match what was used to generate the file.
      ! Sentinel values (-1 / -999) mean "not specified - use whatever is in the file".
      call validate_gvec_namelist_vs_file( &
         gvec_filename, s_file, alpha0_file, zeta_center_file, nfp_sim_file, nfp_out, &
         torflux_namelist, alpha0_namelist, zeta_center_namelist, nfield_periods_namelist)

      surf%rhoc = sqrt(s_file)    ! rhoc = rho = sqrt(normalised toroidal flux)
      surf%rhotor = surf%rhoc
      nfp_real = real(nfp_out)

      ! nfp_sim is inferred from dzdzeta: dzdzeta = nfp / nfp_sim
      nfp_sim_real = nfp_real / dzdzeta_file
      nfp_sim = nint(nfp_sim_real)
      field_period_ratio = nfp_sim_real / nfp_real

      ! zeta step on the fine grid
      dzeta_file = (zeta_file(nzgrid_file) - zeta_file(-nzgrid_file)) &
                   / real(2 * nzgrid_file)

      ! Arc-length coordinate
      if (zed_equal_arc) then

         if (debug) write (*, *) 'gvec_geometry::zed_equal_arc'
         if (z_grid_refinement_factor <= 1) then
            write (*, *) 'WARNING: zed_equal_arc=T but z_grid_refinement_factor=1.'
            write (*, *) 'For accurate arc-length, re-run gvec_to_stella.py with'
            write (*, *) '--nzgrid_file = <nzgrid * factor>  (factor >= 4 recommended)'
         end if

         zetamax_idx = nzgrid_file

         do ia = 1, nalpha
            ! Total (un-normalised) arc length: int 1/(b.grad zeta) dzeta
            call get_total_arc_length( &
               zetamax_idx, b_dot_gradzeta_file(ia, -zetamax_idx:zetamax_idx), &
               dzeta_file, zed_domain_size(ia))

            ! The arc-length grid runs over [-L/2, +L/2], so the lower boundary
            ! is -zed_domain_size/2 -- NOT the lower zeta of the domain. Both
            ! index arguments must be the FINE extent, since <b_dot_gradzeta_file>
            ! and <arc_length> are dimensioned (-nzgrid_file:nzgrid_file);
            ! passing the coarse <nzgrid> mis-declares the dummy array bounds.
            zmin = -zed_domain_size(ia) * 0.5
            call get_arc_length_grid( &
               zetamax_idx, nzgrid_file, zmin, b_dot_gradzeta_file(ia, :), &
               dzeta_file, arc_length(ia, :))
         end do

         ! Compress/expand the arc length onto [-pi, pi], which is the range of
         ! stella's <zed>. Without this the spline abscissa is in physical arc
         ! length while the target grid is in [-pi, pi], and every interpolation
         ! extrapolates off the end of the data. (Mirrors get_vmec_geometry.)
         arc_length = 2.0 * pi * arc_length &
                      / spread(zed_domain_size, 2, 2 * nzgrid_file + 1)

         ! Interpolate all arrays from (zeta_file) to arc-length grid (zed)
         ! alpha and zeta are taken from the file (they carry the field-line/domain info)
         alpha(:) = alpha_file(:)
         ! In arc-length mode the zed grid is uniform in the normalised arc length, so
         ! the toroidal angle zeta at the zed points is interpolated from the file, as for VMEC.
         do ia = 1, nalpha
            call geo_spline(arc_length(ia, :), zeta_file, zed, zeta(ia, :))
         end do

         do ia = 1, nalpha
            b_dot_grad_arclength(ia, :) = 2.0 * pi / zed_domain_size(ia)

            call geo_spline(arc_length(ia, :), bmag_file(ia, :),  zed, bmag(ia, :))
            call geo_spline(arc_length(ia, :), grho_file(ia, :),  zed, grho(ia, :))
            call geo_spline(arc_length(ia, :), theta_file(ia, :), zed, theta_gvec(ia, :))
            call geo_spline(arc_length(ia, :), gag_file(ia, :),   zed, grad_alpha_grad_alpha(ia, :))
            call geo_spline(arc_length(ia, :), gap_file(ia, :),   zed, grad_alpha_grad_psit(ia, :))
            call geo_spline(arc_length(ia, :), gpp_file(ia, :),   zed, grad_psit_grad_psit(ia, :))
            call geo_spline(arc_length(ia, :), gds23_file(ia, :), zed, gds23_psitalpha(ia, :))
            call geo_spline(arc_length(ia, :), gds24_file(ia, :), zed, gds24_psitalpha(ia, :))
            call geo_spline(arc_length(ia, :), gds25_file(ia, :), zed, gds25_psitalpha(ia, :))
            call geo_spline(arc_length(ia, :), gds26_file(ia, :), zed, gds26_psitalpha(ia, :))
            call geo_spline(arc_length(ia, :), gbdrift_file(ia, :),   zed, gbdrift_alpha(ia, :))
            call geo_spline(arc_length(ia, :), Btg_dGP_file(ia, :),  zed, B_times_gradB_dot_gradx_psi(ia, :))
            call geo_spline(arc_length(ia, :), Btk_dGA_file(ia, :),  zed, B_times_kappa_dot_grady_alpha(ia, :))
            call geo_spline(arc_length(ia, :), Btk_dGP_file(ia, :),  zed, B_times_kappa_dot_gradx_psi(ia, :))
            call geo_spline(arc_length(ia, :), psit_disp_file(ia, :),  zed, psit_displacement_fac(ia, :))
            call geo_spline(arc_length(ia, :), gza_R2oB2_file(ia, :),  zed, gradzeta_gradalpha_R2overB2(ia, :))
            call geo_spline(arc_length(ia, :), gzp_R2oB2_file(ia, :),  zed, gradzeta_gradpsit_R2overB2(ia, :))
            call geo_spline(arc_length(ia, :), bdgz_RR_file(ia, :),    zed, b_dot_gradzeta_RR(ia, :))

            ! Interpolate b_dot_gradzeta onto the arc-length grid (stored in b_dot_gradz)
            ! so we can use it to rescale gds23 and gds24 in the same loop pass.
            call geo_spline(arc_length(ia, :), b_dot_gradzeta_file(ia, :), zed, b_dot_gradz(ia, :))

            ! Apply the dl/dζ → dl/dl scaling for gds23/gds24:
            !     dz/dζ  = b_dot_grad_arclength / b_dot_gradzeta_on_arc_grid
            do iz = -nzgrid, nzgrid
               if (abs(b_dot_gradz(ia, iz)) > 0.0) then
                  gds23_psitalpha(ia, iz) = gds23_psitalpha(ia, iz) &
                     * b_dot_grad_arclength(ia, iz) / b_dot_gradz(ia, iz)
                  gds24_psitalpha(ia, iz) = gds24_psitalpha(ia, iz) &
                     * b_dot_grad_arclength(ia, iz) / b_dot_gradz(ia, iz)
               end if
            end do
         end do

         ! b_dot_gradz = b · ∇l in arc-length coordinate
         b_dot_grad_arclength_avg = sum(b_dot_grad_arclength, 1) / real(nalpha)

         ! dzdzeta does not apply in arc-length (zed = l/total_length * 2pi)
         dzdzeta = 1.0

         b_dot_gradz_avg = b_dot_grad_arclength_avg
         b_dot_gradz     = b_dot_grad_arclength

      else

         ! Non arc-length: copy the fine-grid arrays directly
         if (debug) write (*, *) 'gvec_geometry::not_zed_equal_arc'

         zeta(:, :) = spread(zeta_file(-nzgrid:nzgrid), 1, nalpha)
         alpha(:)   = alpha_file(:)
         theta_gvec  = theta_file(:, -nzgrid:nzgrid)
         grho        = grho_file(:, -nzgrid:nzgrid)
         bmag        = bmag_file(:, -nzgrid:nzgrid)
         grad_alpha_grad_alpha = gag_file(:, -nzgrid:nzgrid)
         grad_alpha_grad_psit  = gap_file(:, -nzgrid:nzgrid)
         grad_psit_grad_psit   = gpp_file(:, -nzgrid:nzgrid)
         gds23_psitalpha       = gds23_file(:, -nzgrid:nzgrid)
         gds24_psitalpha       = gds24_file(:, -nzgrid:nzgrid)
         gds25_psitalpha       = gds25_file(:, -nzgrid:nzgrid)
         gds26_psitalpha       = gds26_file(:, -nzgrid:nzgrid)
         gbdrift_alpha             = gbdrift_file(:, -nzgrid:nzgrid)
         B_times_gradB_dot_gradx_psi = Btg_dGP_file(:, -nzgrid:nzgrid)
         B_times_kappa_dot_grady_alpha = Btk_dGA_file(:, -nzgrid:nzgrid)
         B_times_kappa_dot_gradx_psi = Btk_dGP_file(:, -nzgrid:nzgrid)
         psit_displacement_fac = psit_disp_file(:, -nzgrid:nzgrid)
         gradzeta_gradalpha_R2overB2 = gza_R2oB2_file(:, -nzgrid:nzgrid)
         gradzeta_gradpsit_R2overB2  = gzp_R2oB2_file(:, -nzgrid:nzgrid)
         b_dot_gradzeta_RR   = bdgz_RR_file(:, -nzgrid:nzgrid)

         b_dot_gradzeta_file(:, :) = b_dot_gradzeta_file(:, :) * dzdzeta_file
         b_dot_gradzeta_avg = sum(b_dot_gradzeta_file(:, -nzgrid:nzgrid), 1) / real(nalpha)
         b_dot_gradz_avg = b_dot_gradzeta_avg(-nzgrid:nzgrid)
         b_dot_gradz     = b_dot_gradzeta_file(:, -nzgrid:nzgrid)

         ! dzdzeta: the scaling factor between z (stella's parallel coordinate,
         ! restricted to [-pi,pi]) and zeta (physical toroidal angle).
         ! z = zeta * (nfp/nfp_sim), so dz/dzeta = nfp/nfp_sim.
         dzdzeta = dzdzeta_file

         ! Apply the dzdzeta scaling to gds23, gds24
         gds23_psitalpha = gds23_psitalpha * dzdzeta
         gds24_psitalpha = gds24_psitalpha * dzdzeta

      end if

      ! Populate flux-surface struct
      ! surf%rhoc was set above; compute remaining scalar quantities

      ! grho = sqrt(grad_psit_grad_psit) / rhotor  (a |∇ρ|)
      ! Already stored in the file; nothing to recompute.

      surf%psitor_lcfs = 0.5 * sign_torflux
      surf%drhotordrho = 1.0
      surf%rhoc_psi0   = surf%rhoc
      surf%qinp_psi0   = surf%qinp
      surf%shat_psi0   = surf%shat

      ! zed0_fac: needed for the ballooning-mode kx-ky mapping
      ! = -zed(nzgrid) / zeta(1,nzgrid) * q
      if (abs(zeta(1, nzgrid)) > 0.0) then
         surf%zed0_fac = -zed(nzgrid) / zeta(1, nzgrid) * surf%qinp
      else
         surf%zed0_fac = 0.0
      end if

      ! theta: the straight-field-line angle, i.e.
      !     theta_pest(alpha, zeta) = alpha + iota * zeta
      ! Overwrite theta_gvec (the GVEC logical angle) with theta_pest to match
      ! the convention used by stella's VMEC branch.
      ! theta_gvec already holds the GVEC logical angle theta; we convert it to
      ! the PEST straight-field-line angle here.  Alternatively, one could keep
      ! the logical angle for diagnostic purposes, but the VMEC branch stores
      ! theta_pest (= alpha + zeta/q) as the "theta_vmec" array.
      theta_gvec = spread(alpha, 2, 2 * nzgrid + 1) + zeta / surf%qinp

      ! Write geometry summary to output file
      if (verbose .and. print_extra_info_to_terminal) then
         call open_output_file(tmpunit, '.gvec_run_data')
         write (tmpunit, '(6a12)') '#rhotor', 'qinp', 'shat', 'aref', 'Bref', 'dzdzeta'
         write (tmpunit, '(6e12.4)') surf%rhoc, surf%qinp, surf%shat, &
                                     L_reference, B_reference, dzdzeta
         close (tmpunit)
      end if

      ! Deallocate
      deallocate (zeta_file, alpha_file, arc_length, zed_domain_size)
      deallocate (b_dot_gradzeta_file, b_dot_grad_arclength)
      deallocate (b_dot_gradzeta_avg, b_dot_grad_arclength_avg)
      deallocate (theta_file, grho_file, bmag_file)
      deallocate (gag_file, gap_file, gpp_file)
      deallocate (gds23_file, gds24_file, gds25_file, gds26_file)
      deallocate (gbdrift_file, Btg_dGP_file, Btk_dGA_file, Btk_dGP_file)
      deallocate (psit_disp_file, gza_R2oB2_file, gzp_R2oB2_file, bdgz_RR_file)
      deallocate (B_sub_theta_file, B_sub_zeta_file)

   end subroutine get_gvec_geometry


   !============================================================================
   !================ VALIDATE NAMELIST PARAMETERS VS FILE =====================
   !============================================================================
   ! Called after reading the netCDF file.  If the user specified any of
   ! torflux / alpha0 / zeta_center / nfield_periods in &geometry_gvec, this
   ! routine checks them against the values in the file and aborts with a
   ! helpful message if they don't match.
   !
   ! Sentinel values (torflux=-1, alpha0=-999, zeta_center=-999, nfield_periods=-1)
   ! mean "not set by user - skip check for this parameter".
   !============================================================================
   subroutine validate_gvec_namelist_vs_file( &
      gvec_filename, s_file, alpha0_file, zeta_center_file, nfp_sim_file, nfp_file, &
      torflux_nl, alpha0_nl, zeta_center_nl, nfield_periods_nl)

      use mp, only: mp_abort
      implicit none

      character(*), intent(in) :: gvec_filename
      real,         intent(in) :: s_file, alpha0_file, zeta_center_file
      integer,      intent(in) :: nfp_sim_file, nfp_file
      real,         intent(in) :: torflux_nl, alpha0_nl, zeta_center_nl, nfield_periods_nl

      real, parameter :: tol = 1.0e-4   ! tolerance for floating-point comparisons
      logical :: ok
      character(len=2000) :: cmd

      ok = .true.

      ! Build the prototype gvec_to_stella.py command to show in error messages
      write (cmd, '(A)') "  python gvec_to_stella.py STATE_FILE"

      if (torflux_nl > 0.0 .and. abs(torflux_nl - s_file) > tol) then
         write (*, '(A)') " "
         write (*, '(A)') "  ERROR (geometry_gvec): torflux mismatch"
         write (*, '(A,F8.5,A,F8.5)') "    namelist: torflux = ", torflux_nl, &
                                       "   file: s = ", s_file
         write (*, '(A)') "  Re-run gvec_to_stella.py with:"
         write (*, '(A,F6.4)') "    --s ", torflux_nl
         ok = .false.
      end if

      if (alpha0_nl > -998.0 .and. abs(alpha0_nl - alpha0_file) > tol) then
         write (*, '(A)') " "
         write (*, '(A)') "  ERROR (geometry_gvec): alpha0 mismatch"
         write (*, '(A,F8.5,A,F8.5)') "    namelist: alpha0 = ", alpha0_nl, &
                                       "   file: alpha(1) = ", alpha0_file
         write (*, '(A)') "  Re-run gvec_to_stella.py with:"
         write (*, '(A,F8.5)') "    --alpha0 ", alpha0_nl
         ok = .false.
      end if

      if (zeta_center_nl > -998.0 .and. abs(zeta_center_nl - zeta_center_file) > tol) then
         write (*, '(A)') " "
         write (*, '(A)') "  ERROR (geometry_gvec): zeta_center mismatch"
         write (*, '(A,F8.5,A,F8.5)') "    namelist: zeta_center = ", zeta_center_nl, &
                                       "   file: zeta_center = ", zeta_center_file
         write (*, '(A)') "  Re-run gvec_to_stella.py with:"
         write (*, '(A,F8.5)') "    --zeta_center ", zeta_center_nl
         ok = .false.
      end if

      if (nfield_periods_nl > 0.0) then
         if (abs(nfield_periods_nl - real(nfp_sim_file)) > tol) then
            write (*, '(A)') " "
            write (*, '(A)') "  ERROR (geometry_gvec): nfield_periods mismatch"
            write (*, '(A,F6.1,A,I4)') "    namelist: nfield_periods = ", nfield_periods_nl, &
                                        "   file: nfp_sim = ", nfp_sim_file
            write (*, '(A)') "  Re-run gvec_to_stella.py with:"
            write (*, '(A,I4)') "    --nfield_periods ", nint(nfield_periods_nl)
            ok = .false.
         end if
      end if

      if (.not. ok) then
         write (*, '(A)') " "
         write (*, '(A)') "  The parameters above must match the values used when generating"
         write (*, '(2A)') "  the geometry file: ", trim(gvec_filename)
         write (*, '(A)') "  Full example gvec_to_stella.py command (adjust STATE_FILE):"
         write (*, '(A)') "    python gvec_to_stella.py STATE_FILE.dat \\"
         if (torflux_nl   > 0.0)    write (*, '(A,F6.4,A)')  "        --s ", torflux_nl, " \\"
         if (alpha0_nl    > -998.0) write (*, '(A,F8.5,A)')  "        --alpha0 ", alpha0_nl, " \\"
         if (zeta_center_nl > -998.0) write (*, '(A,F8.5,A)')"        --zeta_center ", zeta_center_nl, " \\"
         if (nfield_periods_nl > 0.0) write (*, '(A,I4,A)') "        --nfield_periods ", nint(nfield_periods_nl), " \\"
         write (*, '(A)')  "        --nzgrid <nzgrid> --nzgrid_file <nzgrid_file> -o new_geometry.nc"
         write (*, '(A)') " "
         call mp_abort('geometry_gvec mismatch: please re-run gvec_to_stella.py')
      end if

   end subroutine validate_gvec_namelist_vs_file


end module gvec_geometry
