!###############################################################################
!                    Define a sheared slab (and curvy slab) equilibrium
!###############################################################################
!
! Implements a sheared slab magnetic geometry for stella.
!
! The magnetic field is uniform with B = B_ref * z_hat, and the magnetic
! shear enters through the metric tensor. The parallel coordinate z ranges
! from -nperiod*pi to +nperiod*pi and is normalised by the perpendicular
! scale length a_ref (so z here is the field-line arc-length divided by a_ref).
!
! In the sheared slab the (x, y, z) coordinates are related to the
! Cartesian coordinates (X, Y, Z) by:
!   x = X / a_ref        (radial direction)
!   y = Y / a_ref        (bi-normal direction)
!   z = Z / a_ref        (parallel direction, proportional to arc length)
!
! The metric coefficients for this geometry are:
!   |grad y|^2 = 1 + (shat * z)^2
!   grad x . grad y = -shat * z
!   |grad x|^2 = 1
!   b . grad z = 1 / (q * R_maj)
!
! For the "curvy slab" option, cylindrical curvature is added:
!   B x grad B . grad y = 1 / R_maj
!   B x kappa . grad y = 1 / R_maj + betaprim
!
! All other drift terms are zero for flat slab (curvy_slab = .false.),
! except B x kappa . grad y = betaprim (pressure-driven drift).
!
! References:
!   - Cowley, Kulsrud & Sudan, Phys. Fluids B 3, 2767 (1991) [slab GK]
!   - Hammett et al., Plasma Phys. Control. Fusion 35, 973 (1993) [GS2 slab]
!   - Kotschenreuther et al., Comp. Phys. Comm. 88, 128 (1995) [GS2]
!   - Parra & Catto, Plasma Phys. Control. Fusion 50, 065014 (2008) [slab metric]
!
!###############################################################################
module slab_geometry

   implicit none

   public :: get_slab_geometry_coefficients

   private

   real :: shat_slab, qinp_slab, rmaj_slab, betaprim_slab
   logical :: curvy_slab_flag

contains

   !****************************************************************************
   !               Calculate geometric quantities for a slab equilibrium
   !****************************************************************************
   subroutine get_slab_geometry_coefficients(nzgrid, zed, bmag, b_dot_gradz_avg, grad_rho, surf, &
      grad_y_dot_grad_y, grad_x_dot_grad_y, grad_x_dot_grad_x, &
      B_times_gradB_dot_gradx, B_times_gradB_dot_grady, &
      B_times_kappa_dot_gradx, B_times_kappa_dot_grady, btor, rmajor)

      use common_types, only: flux_surface_type
      use namelist_geometry, only: read_namelist_geometry_slab

      implicit none

      ! Arguments
      integer, intent(in) :: nzgrid
      real, dimension(-nzgrid:), intent(in) :: zed
      real, dimension(-nzgrid:), intent(out) :: bmag, b_dot_gradz_avg, grad_rho
      real, dimension(-nzgrid:), intent(out) :: grad_y_dot_grad_y, grad_x_dot_grad_y
      real, dimension(-nzgrid:), intent(out) :: grad_x_dot_grad_x
      real, dimension(-nzgrid:), intent(out) :: B_times_gradB_dot_gradx, B_times_gradB_dot_grady
      real, dimension(-nzgrid:), intent(out) :: B_times_kappa_dot_gradx, B_times_kappa_dot_grady
      real, dimension(-nzgrid:), intent(out) :: btor, rmajor
      type(flux_surface_type), intent(out) :: surf

      integer :: iz

      !-------------------------------------------------------------------------

      ! Read the slab namelist from the input file
      call read_namelist_geometry_slab(shat_slab, qinp_slab, rmaj_slab, betaprim_slab, curvy_slab_flag)

      !-------------------------------------------------------------------------
      ! Magnetic field magnitude: B = B_ref everywhere in a slab
      !-------------------------------------------------------------------------
      bmag = 1.0

      !-------------------------------------------------------------------------
      ! Parallel streaming: b . grad z = 1 / (q * R_maj)
      ! In stella's normalisation, with z ~ theta and safety factor q,
      ! the field-line length per radian is q*R, giving b.grad(z) = 1/(q*R).
      !-------------------------------------------------------------------------
      b_dot_gradz_avg = 1.0 / (qinp_slab * rmaj_slab)

      !-------------------------------------------------------------------------
      ! Metric tensor elements for a sheared slab
      ! With magnetic shear shat, the perpendicular metric is:
      !   |grad y|^2  = 1 + (shat * z)^2
      !   grad x.grad y = -shat * z
      !   |grad x|^2  = 1
      !-------------------------------------------------------------------------
      do iz = -nzgrid, nzgrid
         grad_y_dot_grad_y(iz) = 1.0 + (shat_slab * zed(iz))**2
         grad_x_dot_grad_y(iz) = -shat_slab * zed(iz)
         grad_x_dot_grad_x(iz) = 1.0
      end do

      !-------------------------------------------------------------------------
      ! grad rho = |grad x| = 1 in normalised slab coordinates
      !-------------------------------------------------------------------------
      grad_rho = 1.0

      !-------------------------------------------------------------------------
      ! Drift terms
      ! For a flat slab: grad B = 0, curvature = 0 --> all drift terms zero
      !   except for the pressure-driven term betaprim in B_times_kappa_dot_grady
      ! For a curvy slab (cylindrical curvature with major radius rmaj):
      !   B x grad B . grad y = 1 / R_maj  (grad-B drift)
      !   B x kappa . grad y  = 1 / R_maj + betaprim  (curvature + pressure drift)
      ! The x-components (gradx) are zero in both cases (no toroidal variation).
      !-------------------------------------------------------------------------
      B_times_gradB_dot_gradx = 0.0
      B_times_kappa_dot_gradx = 0.0

      if (curvy_slab_flag) then
         ! Add cylindrical curvature: 1/R_maj contribution
         B_times_gradB_dot_grady = 1.0 / rmaj_slab
         B_times_kappa_dot_grady = 1.0 / rmaj_slab + betaprim_slab
      else
         ! Pure flat slab: no curvature or grad-B drift
         B_times_gradB_dot_grady = 0.0
         B_times_kappa_dot_grady = betaprim_slab
      end if

      !-------------------------------------------------------------------------
      ! Toroidal field and major radius
      ! btor = 0 means all flow shear is treated as perpendicular
      ! rmajor = rmaj_slab for curvy slab, 1.0 for flat slab (avoids divide by zero)
      !-------------------------------------------------------------------------
      btor = 0.0
      if (curvy_slab_flag) then
         rmajor = rmaj_slab
      else
         rmajor = 1.0
      end if

      !-------------------------------------------------------------------------
      ! Flux surface parameters
      ! These are set to physically meaningful values for slab geometry.
      ! rhoc = 1.0 (we are at the reference surface)
      ! kappa = 1.0 (circular cross-section, though not really meaningful for slab)
      !-------------------------------------------------------------------------
      surf%shat = shat_slab
      surf%qinp = qinp_slab
      surf%rmaj = rmaj_slab
      surf%rgeo = rmaj_slab
      surf%rhoc = 1.0
      surf%kappa = 1.0
      surf%kapprim = 0.0
      surf%tri = 0.0
      surf%triprim = 0.0
      surf%shift = 0.0
      surf%betaprim = betaprim_slab
      surf%betadbprim = 0.0
      surf%d2qdr2 = 0.0
      surf%d2psidr2 = 0.0
      surf%dpsitordrho = 1.0
      surf%rhotor = 1.0
      surf%drhotordrho = 1.0
      surf%psitor_lcfs = 1.0

      ! psi0 variants (reference surface values -- same as current surface for slab)
      surf%rhoc_psi0 = surf%rhoc
      surf%qinp_psi0 = surf%qinp
      surf%shat_psi0 = surf%shat

   end subroutine get_slab_geometry_coefficients

end module slab_geometry
