!###############################################################################
!############################### STELLA GEOMETRY ###############################
!###############################################################################
! 
! Routines for calculating the geometric quantities needed by stella.
! 
! This routine will call the <vmec_geometry> or <geometry_miller> modules.
! Which each uses specific (alpha, psi) and (x,y) coordinates. Nonetheless,
! stella is completely general, it just needs geometric variables as a function 
! of <psi> as well as the factor <dxdpsi> and <dpsitdpsi>.
! 
!---------------------------- Geometric quantities -----------------------------
! 
! The following geometric quantities are defined in this module:
!    <bmag>(ia,iz) = B / B_ref
!    <gradx_dot_gradx>(ia,iz) = |∇x|²
!    <grady_dot_grady>(ia,iz) = |∇y|²
!    <gradx_dot_grady>(ia,iz) = ∇x . ∇y
!    <B_times_gradB_dot_gradx>(ia,iz) = B × ∇B · ∇x (a*B_ref/B^3)
!    <B_times_gradB_dot_grady>(ia,iz) = B × ∇B · ∇y (a*B_ref/B^3)
!    <B_times_kappa_dot_gradx>(ia,iz) = B × κ · ∇x (a*B_ref/B^2)
!    <B_times_kappa_dot_grady>(ia,iz) = B × κ · ∇y (a*B_ref/B^2)
!    <b_dot_gradz>(ia,iz) = b · ∇z
!    <b_dot_gradz_avg>(iz) = sum_alpha b · ∇z J dalpha / sum_alpha J dalpha
! 
! Geometric quantities needed for radial variation:
!    <d_gradxdotgradx_drho> = d(∇x . ∇y)/drho
!    <d_gradydotgrady_drho> = d(|∇x|²)/drho
!    <d_gradxdotgrady_drho> = d(|∇y|²)/drho
! 
! The normalized derivatives are defined as,
!    <dxdpsi> = (rho_r/a) (d tilde{x} / d tilde{psi})
!    <dydalpha> = (rho_r/a) (d tilde{y} / d tilde{alpha})
! 
!--------------------------- Backwards Compatibility ---------------------------
! 
! An overview of the name changes implemented in September 2025 are given here,
!    - gds22         -->   gradx_dot_gradx * shat * shat
!    - gds2          -->   grady_dot_grady
!    - gds21         -->   gradx_dot_grady * shat
!    - gbdrift0      -->   B_times_gradB_dot_gradx * 2. * shat
!    - gbdrift       -->   B_times_gradB_dot_grady * 2.
!    - cvdrift0      -->   B_times_kappa_dot_gradx * 2. * shat
!    - cvdrift       -->   B_times_kappa_dot_grady * 2.
!    - gradpar       -->   b_dot_gradz_avg
!    - dgradpardrho  -->   d_bdotgradz_drho
!    - gradpar_c     -->   b_dot_gradz_centerdinz (gk_parallel_streaming.f90)
!    - d_gradydotgrady_drho       -->   d_gradydotgrady_drho
!    - d_gradxdotgradx_drho      -->   d_gradxdotgradx_drho
!    - d_gradxdotgrady_drho      -->   d_gradxdotgrady_drho
! 
!###############################################################################
module geometry

   ! Load debug flags
   use debug_flags, only: debug => geometry_debug
   use common_types, only: flux_surface_type
   
   ! Read the parameters for <geo_option_switch> from namelist_geometry.f90
   use namelist_geometry, only: geo_option_local
   use namelist_geometry, only: geo_option_inputprof
   use namelist_geometry, only: geo_option_vmec
   use namelist_geometry, only: geo_option_multibox
   use namelist_geometry, only: geo_option_zpinch

   implicit none

   ! Make routines available to other modules
   public :: init_geometry, finish_init_geometry, finish_geometry
   public :: communicate_geo_multibox, x_displacement_fac
   
   ! Geometric quantities for the gyrokinetic equations
   public :: bmag, dbdzed, btor, bmag_psi0, grho, grho_norm, grad_x
   public :: dcvdriftdrho, dcvdrift0drho, dgbdriftdrho, dgbdrift0drho
   public :: grady_dot_grady, gradx_dot_grady, gradx_dot_gradx, gds23, gds24, gds25, gds26
   public :: B_times_kappa_dot_grady, B_times_kappa_dot_gradx
   public :: B_times_gradB_dot_grady, B_times_gradB_dot_gradx
   public :: d_gradydotgrady_drho, d_gradxdotgrady_drho, d_gradxdotgradx_drho
   public :: exb_nonlin_fac, exb_nonlin_fac_p, flux_fac
   public :: jacob, djacdrho, drhodpsi, drhodpsip, drhodpsip_psi0
   public :: dl_over_b, d_dl_over_b_drho
   public :: dBdrho, d2Bdrdth, d_bdotgradz_drho, dIdrho
   public :: geo_surf, Rmajor, dzetadz
   public :: theta_vmec, zeta, alpha 
   public :: dxdpsi, dydalpha, clebsch_factor
   public :: aref, bref, twist_and_shift_geo_fac
   public :: q_as_x, get_x_to_rho, gfac
   public :: dVolume, grad_x_grad_y_end
   
   ! Flux tube only needs b_dot_gradz_avg(z)
   public :: b_dot_gradz_avg
   
   ! Full flux surface needs b_dot_gradz(alpha, z)
   public :: b_dot_gradz
   
   ! Extended z-grid for final fields diagnostics
   public :: b_dot_gradz_avg_eqarc, zed_eqarc
   
   ! Geometric quantities for momentum flux
   public :: gradzeta_gradx_R2overB2
   public :: gradzeta_grady_R2overB2
   public :: b_dot_gradzeta_RR
   
   ! Used in kt_grids.f90
   public :: geo_option_switch, geo_option_vmec

   private
   
   !----------------------------------------------------------------------------
   
   ! Define the geometric quantities of the chosen magnetic flux surface in <geo_surf>:
   ! rmaj, rgeo, kappa, kapprim, tri, triprim, rhoc, dr, shift, qinp
   ! shat, betaprim, betadbprim, d2qdr2, d2psidr2, dpsitordrho, d2psitordrho2
   ! rhotor, drhotordrho, psitor_lcfs, zed0_fac, rhoc_psi0, qinp_psi0, shat_psi0
   type(flux_surface_type) :: geo_surf

   ! Geometric quantities
   real :: grad_x_grad_y_end, clebsch_factor
   real :: aref, bref, dxdpsi, dydalpha
   real :: dqdrho, dIdrho, grho_norm
   real :: drhodpsi, drhodpsip, drhodpsip_psi0, shat, qinp
   real :: exb_nonlin_fac, exb_nonlin_fac_p, flux_fac
   real :: b_dot_gradz_avg_eqarc, dzetadz
   real :: twist_and_shift_geo_fac, gfac
   integer :: sign_torflux
   integer :: geo_option_switch

   ! Geometric quantities for the gyrokinetic equations
   real, dimension(:), allocatable :: zed_eqarc, alpha
   real, dimension(:), allocatable :: b_dot_gradz_avg
   real, dimension(:), allocatable :: dBdrho, d2Bdrdth, d_bdotgradz_drho, btor, Rmajor 
   real, dimension(:, :), allocatable :: bmag, bmag_psi0, dbdzed 
   real, dimension(:, :), allocatable :: B_times_kappa_dot_grady, B_times_kappa_dot_gradx, B_times_gradB_dot_grady, B_times_gradB_dot_gradx
   real, dimension(:, :), allocatable :: dcvdriftdrho, dcvdrift0drho, dgbdriftdrho, dgbdrift0drho
   real, dimension(:, :), allocatable :: grady_dot_grady, gradx_dot_grady, gradx_dot_gradx, gds23, gds24, gds25, gds26
   real, dimension(:, :), allocatable :: d_gradydotgrady_drho, d_gradxdotgrady_drho, d_gradxdotgradx_drho, x_displacement_fac
   real, dimension(:, :), allocatable :: jacob, djacdrho, grho, grad_x
   real, dimension(:, :), allocatable :: dl_over_b, d_dl_over_b_drho 
   real, dimension(:, :), allocatable :: theta_vmec, zeta
   real, dimension(:, :, :), allocatable :: dVolume
   
   ! Geometric quantities for full flux surface
   real, dimension(:, :), allocatable :: b_dot_gradz
   
   ! Geometric quantities for the momentum flux
   real, dimension(:, :), allocatable :: gradzeta_gradx_R2overB2
   real, dimension(:, :), allocatable :: gradzeta_grady_R2overB2
   real, dimension(:, :), allocatable :: b_dot_gradzeta_RR

 
   ! The geometric quantities can be read from an old geometry file
   logical :: overwrite_bmag, overwrite_b_dot_gradzeta, overwrite_geometry
   logical :: overwrite_grady_dot_grady, overwrite_gradx_dot_grady, overwrite_gradx_dot_gradx
   logical :: overwrite_gds23, overwrite_gds24, overwrite_B_times_gradB_dot_grady
   logical :: overwrite_B_times_kappa_dot_grady, overwrite_B_times_gradB_dot_gradx, q_as_x
   character(100) :: geometry_file
  
   ! Only initialise once
   logical :: initialised_geometry = .false.   

contains

   !============================================================================
   !========================= INITIALIZE THE GEOMETRY ==========================
   !============================================================================
   subroutine init_geometry(nalpha, naky)

      ! Grids
      use grids_z, only: nzgrid, zed, delzed, shat_zero, grad_x_grad_y_zero
      use grids_z, only: boundary_option_switch, boundary_option_self_periodic
      use grids_z, only: boundary_option_linked, boundary_option_linked_stellarator

      ! VMEC equilibria
      use vmec_geometry, only: read_vmec_parameters, get_vmec_geometry 

      ! Flags
      use parameters_multibox, only: include_geometric_variation

      ! Routines
      use namelist_geometry, only: read_namelist_geometry_from_txt, read_namelist_geometry_options
      use mp, only: proc0

      implicit none

      integer, intent(in) :: nalpha, naky
      real :: bmag_z0
      integer :: iy, ia, iz

      !---------------------------------------------------------------------- 

      ! Only initialize once
      if (initialised_geometry) return
      initialised_geometry = .true.

      ! Default is no re-scaling of zed
      dzetadz = 1.0
      
      ! Track the code 
      if (debug) write (*, *) 'geometry::init_geometry'

      ! Only calculate the geometry on proc0
      if (proc0) then
         
         ! Read the <geo_knobs> namelist in the input file
         call read_namelist_geometry_options (geo_option_switch, q_as_x)

         ! Only read in the geometry file if the geometry_option is set to be an input profile
         if (geo_option_switch==geo_option_inputprof) then
            call read_namelist_geometry_from_txt(geometry_file, overwrite_bmag, overwrite_b_dot_gradzeta, &
               overwrite_grady_dot_grady, overwrite_gradx_dot_grady, overwrite_gradx_dot_gradx, overwrite_gds23, overwrite_gds24, &
               overwrite_B_times_gradB_dot_grady, overwrite_B_times_kappa_dot_grady, &
               overwrite_B_times_gradB_dot_gradx, overwrite_geometry)
         end if

         ! Use Miller parameters or VMEC to get the geometry needed for stella
         if (geo_option_switch==geo_option_local)     call get_geometry_arrays_from_Miller(nalpha)
         if (geo_option_switch==geo_option_inputprof) call get_geometry_arrays_from_Miller(nalpha)
         if (geo_option_switch==geo_option_multibox)  call get_geometry_arrays_from_Miller(nalpha)
         if (geo_option_switch==geo_option_vmec)      call get_geometry_arrays_from_VMEC(nalpha, naky) 
         if (geo_option_switch==geo_option_zpinch)    call get_geometry_arrays_from_zpinch(nalpha)
         
         ! Overwrite the selected geometric coefficients
         if (overwrite_geometry) call overwrite_selected_geometric_coefficients(nalpha)

         ! <exb_nonlin_fac> = -(0.5/C)*Bref*(dx/dpsi)(dy/dalpha) = - 0.5 * (1/<clebsch_factor>) * <dxdpsi> * <dydalpha>
         exb_nonlin_fac = -0.5 / clebsch_factor * dxdpsi * dydalpha
         exb_nonlin_fac_p = 0.0

         ! <flux_fac> = -(0.5/C)*(drho/dpsi)*(dy/dalpha)
         flux_fac = -0.5 / clebsch_factor * drhodpsi * dydalpha

         ! If <radial_variation> = True than <q_as_x> = True
         ! The following will get multiplied by <exb_nonlin_fac> in <advance_exb_nonlinearity>
         if (q_as_x) exb_nonlin_fac_p = geo_surf%d2qdr2 / dqdrho - geo_surf%d2psidr2 * drhodpsip
         
         ! In the diagnostics of final fields we want out extented z-grid in arc length units
         ! So here we calculate <zed_eqarc> which is z = arc length
         call get_b_dot_gradz_avg_eqarc(b_dot_gradz_avg, zed, delzed, b_dot_gradz_avg_eqarc)
         call get_zed_eqarc(b_dot_gradz_avg, delzed, zed, b_dot_gradz_avg_eqarc, zed_eqarc)

      end if
 
      !=========================================================================
      !=================== CALCULATIONS ON ALL PROCESSORS  =====================
      !=========================================================================
      
      ! Track the code 
      if (debug) write (*, *) 'geometry::init_geometry::calculate_on_all_processors'

      ! We calculated the geometric quantities on proc0, now allocate the arrays on 
      ! the other processors, and broadcast from proc0 to the other processors
      if (.not. proc0) call allocate_arrays(nalpha, nzgrid) 
      call broadcast_arrays

      ! <gfac> will allow us to include the geometric variation
      ! By default <include_geometric_variation> = True
      if (include_geometric_variation) gfac = 1.0
      if (.not. include_geometric_variation) gfac = 0.0
 
      ! <dqdrho> = dq/drho = hat{s} * q/rho  with  hat{s} = r/q (dq/dr) = rho/q (dq/drho)
      ! FLAG DSO - the following assumes a linear relation from q to rho, but this will not be correct if d2qdrho != 0
      dqdrho = geo_surf%shat * geo_surf%qinp / geo_surf%rhoc

      ! <jacob> is the Jacobian from Cartesian coordinates to (y,x,z) coordinates with B = C ∇ψ x ∇α
      ! jacob^{-1} = (∇y x ∇x) . ∇z = (dy/dalpha) (dx/dpsi) (∇α x ∇ψ) . ∇z)
      !            = - (1/C) (dy/dalpha) (dx/dpsi) (B . ∇z)
      !            = - (1/C) (d(y/a)/dalpha) (d(x/a)/d(psi/a^2Bref)) ((B/Bref) . ∇z)
      !            = - (a^2Bref/C) (d(y/a)/dalpha) (d(x/a)/d(psi)) ((bmag*b/Bref) . ∇z)  
      ! Warning <jacob> was not correct for <radial_variation>, nonetheless, it was only ever used 
      ! in both the numerator and denominator of averages -- and so any constant factors cancelled out 
      jacob = -clebsch_factor / (dydalpha * dxdpsi * b_dot_gradz * bmag)
   
      ! <dl_over_b> = dl/J are the integration weights along the field line 
      ! For flux tube simulations with psi = psit and psi = psip it reduces to <dl_over_b> = dl/B (?)
      dl_over_b = spread(delzed, 1, nalpha) * jacob

      ! Avoid double counting at the end points for ky = 0 modes (which leads to destabilization of the zonal modes)
      ! FLAG DSO - while this is correct for ky = 0 modes and sufficient for output, if dl_over_b is applied to
      ! non-zero ky modes, a more sophisticated approach will be required that takes into account the sign of v_parallel
      dl_over_b(:, nzgrid) = 0.

      ! Correction to flux-surface-averaging for adiabatic electrons
      ! FLAG COOKIE TO GEORGIA - <djacdrho> is not defined for VMEC?
      d_dl_over_b_drho = spread(delzed, 1, nalpha) * djacdrho
      d_dl_over_b_drho(:, nzgrid) = 0
      d_dl_over_b_drho = d_dl_over_b_drho - dl_over_b &
            * spread(sum(d_dl_over_b_drho, dim=2) / sum(dl_over_b, dim=2), 2, 2 * nzgrid + 1)
      d_dl_over_b_drho = gfac * d_dl_over_b_drho / spread(sum(dl_over_b, dim=2), 2, 2 * nzgrid + 1)

      ! Normalize dl/B by int dl/B
      dl_over_b = dl_over_b / spread(sum(dl_over_b, dim=2), 2, 2 * nzgrid + 1)

      ! We normalize the fluxes with sum( dl/J * |nabla rho| )
      grho_norm = sum(dl_over_b(1, :) * grho(1, :))

      ! FLAG - would probably be better to compute this in the various geometry
      ! subroutines (Miller, VMEC, etc.), as there B is likely calculated on a finer z-grid
      do iy = 1, nalpha
         call get_dzed(nzgrid, delzed, bmag(iy, :), dbdzed(iy, :))
      end do

      ! Change the boundary conditions if the shear is too low or if |∇x . ∇y| is too low at the ends of the field line 
      select case (boundary_option_switch)

      ! If the magnetic shear is almost zero, override the parallel
      ! boundary condition so that it is periodic if using the standard
      ! twist and shift bc, in which kx_shift is proportional to shat
      case (boundary_option_linked)
         if (abs(geo_surf%shat) <= shat_zero) then
            write (*, *) 'Using periodic boundary conditions as shat < shat_zero' 
            boundary_option_switch = boundary_option_self_periodic
         end if

      ! If the magnetic |∇x . ∇y| is almost zero, override parallel boundary 
      ! condition so that it is periodic if using the stellarator symmetric
      !twist and shift bc, in which kx_shift is proportional to |∇x . ∇y|
      ! TODO: this will fail for Miller + boundary_option_linked_stellarator since grad_x_grad_y_end isn't set
      case (boundary_option_linked_stellarator)
         if (abs(grad_x_grad_y_end) <= grad_x_grad_y_zero) then
            write (*, *) 'Using periodic boundary conditions as grad_x_grad_y_end < grad_x_grad_y_zero'
            boundary_option_switch = boundary_option_self_periodic 
         end if

      end select

      if (proc0) call write_geometric_coefficients(nalpha)
      
      ! Deallocate the local arrays within the Miller module
      call finish_init_geometry

   end subroutine init_geometry

   !======================================================================
   !====================== READ GEOMETRY FROM VMEC =======================
   !======================================================================
   ! Regardless of the choice of the coordinates we have
   !     hat{s} = (r/q) (dq/dr) = (rho/q) (dq/drho) = -(rho/iota)(diota/drho)
   !     r = a * sqrt(psi_t/psi_{t,LCFS}) = a * rho = a * sqrt(s)
   !     drho/dψt = sgn(ψt) / (a^2*Bref*rho)
   !     alpha = theta - iota*zeta
   !     iota = 1/q
   ! 
   ! The original stella implementation used the following coordinates,
   ! However, for clockwise VMECs the coordinate x would point inwards
   !     ψ = -ψt = -psi_toroidal = - enclosed toroidal flux divided by 2pi
   !     x = - sgn(ψt)/(r*Bref) (ψ - ψ0)
   !     y = r (α - α0)
   !     dx/dψ = -sgn(ψt)/(r*Bref)           (ρref/a) dx̃/dψ̃ = -sgn(ψt)/rho
   !     dy/dα = r                           (ρref/a) dỹ/dα̃ = r/a = rho
   !     dψ/dψt = -1                         dψ̃/dψ̃t = -1
   !     drho/dψ = -sgn(ψt)/(a^2*Bref*rho)   drho/dψ̃ = -sgn(ψt)/rho
   !     B = - ∇ψ x ∇α                       B/Bref = - (a∇)ψ̃ x (a∇)α̃
   ! 
   ! We implemented a better option, so that the coordinate x always points radially outwards
   !     ψ = sgn(ψt) * ψt
   !     x = 1/(r*Bref) (ψ - ψ0)
   !     y = r (α - α0)
   !     dx/dψ = 1/(r*Bref)                  (ρref/a) dx̃/dψ̃ = a/r = 1/rho
   !     dy/dα = r                           (ρref/a) dỹ/dα̃ = r/a = rho
   !     dψ/dψt = sgn(ψt)                    dψ̃/dψ̃t = sgn(ψt)
   !     drho/dψ = 1/(a^2*Bref*rho)          drho/dψ̃ = 1/rho
   !     B = sgn(ψt) ∇ψ x ∇α                 B/Bref = sgn(ψt) (a∇)ψ̃ x (a∇)α̃
   ! 
   ! The most intuitive option is to choose psi = r
   !     ψ = r
   !     x = (ψ - ψ0)
   !     y = r (α - α0)
   !     dx/dψ = 1                           (ρref/a) dx̃/dψ̃ = 1
   !     dy/dα = r                           (ρref/a) dỹ/dα̃ = r/a = rho
   !     dψ/dψt = sgn(ψt)/(a*Bref*rho)       dψ̃/dψ̃t = sgn(ψt)/rho
   !     drho/dψ = 1/a                       drho/dψ̃ = 1
   !     B = sgn(ψt)*r*Bref ∇ψ x ∇α          B/Bref = sgn(ψt)*(r/a) (a∇)(r/a) x (a∇)α̃ 
   ! 
   ! We need to define the normalized dx/dψ, dy/dalpha, dpsi_t/dpsi, drho/dψ...
   !     <dxdpsi> = (ρref/a)(dx̃/dψ̃) = (ρref/a)(d(x/ρref)/d(ψ/(a^2 Bref)) = a Bref (dx/dψ)     (if psi is a flux)
   !     <dxdpsi> = (ρref/a)(dx̃/dψ̃) = (ρref/a)(d(x/ρref)/d(ψ/a) = (dx/dψ)                     (if psi is a length)
   !     <dydalpha> = (ρref/a)(dỹ/dα̃) = (ρref/a)(d(y/ρref)/dα) = (1/a) (dy/dα)
   ! 
   !======================================================================
   subroutine get_geometry_arrays_from_VMEC(nalpha, naky)

      use mp, only: mp_abort
      use vmec_geometry, only: read_vmec_parameters, get_vmec_geometry
      use vmec_geometry, only: radial_coordinate_switch, radial_coordinate_sgnpsitpsit
      use vmec_geometry, only: radial_coordinate_minuspsit, radial_coordinate_r
      use debug_flags, only: const_alpha_geo
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      integer, intent(in) :: nalpha, naky

      ! Local geometric arrays to construct <grady_dot_grady>, <gradx_dot_grady>, <gradx_dot_gradx>, ...
      real, dimension(:, :), allocatable :: psit_displacement_fac, grad_alpha_grad_alpha
      real, dimension(:, :), allocatable :: grad_alpha_grad_psit, grad_alpha_grad_psi 
      real, dimension(:, :), allocatable :: grad_psit_grad_psit, grad_psi_grad_psi
      real, dimension(:, :), allocatable :: B_times_gradB_dot_gradalpha, B_times_kappa_dot_gradalpha
      real, dimension(:, :), allocatable :: B_times_gradB_dot_gradpsit, B_times_kappa_dot_gradpsit
      real, dimension(:, :), allocatable :: gds23_alphapsit, gds24_alphapsit
      real, dimension(:, :), allocatable :: gds25_alphapsit, gds26_alphapsit
      real, dimension(:, :), allocatable :: grad_x_grad_x, grad_y_grad_y, grad_y_grad_x
      real, dimension(:, :), allocatable :: gradzeta_gradpsit_R2overB2, gradzeta_gradalpha_R2overB2

      ! Local variables
      real :: rho, shat, iota, field_period_ratio
      real :: dpsidpsit, dpsidx, dpsitdx, dxdpsit 

      !----------------------------------------------------------------------
      
      ! Pretty sure q_as_x is not implemented for VMEC
      if (q_as_x) then
         call mp_abort('q_as_x = True is not implemented for VMEC. Aborting.')
      end if 
      
      ! Read the <vmec_parameters> namelist in the input file
      call read_vmec_parameters()

      ! Allocate geometry arrays 
      call allocate_arrays(nalpha, nzgrid)
      call allocate_temporary_arrays(nalpha, nzgrid)

      ! Call the <vmec_geometry> module to calculate the geometric coefficients
      ! needed by stella, based on the VMEC equilibrium file
      call get_vmec_geometry(nzgrid, nalpha, naky, geo_surf, grho, bmag, &
               b_dot_gradz_avg, b_dot_gradz, & 
               grad_alpha_grad_alpha, grad_alpha_grad_psit, grad_psit_grad_psit, &
               gds23_alphapsit, gds24_alphapsit, gds25_alphapsit, gds26_alphapsit, & 
               B_times_gradB_dot_gradalpha, B_times_gradB_dot_gradpsit, B_times_kappa_dot_gradalpha, B_times_kappa_dot_gradpsit, &
               gradzeta_gradpsit_R2overB2, gradzeta_gradalpha_R2overB2, b_dot_gradzeta_RR, &
               sign_torflux, theta_vmec, dzetadz, aref, bref, alpha, zeta, &
               field_period_ratio, psit_displacement_fac)

      ! Flux surface quantities that we need
      rho = geo_surf%rhotor 
      shat = geo_surf%shat
      iota = 1/geo_surf%qinp

      ! Define the (psi,alpha) and (x,y) coordinates
      if (radial_coordinate_switch==radial_coordinate_sgnpsitpsit) then
         dxdpsi = 1. / rho
         dydalpha = rho
         dpsidpsit = sign_torflux
         drhodpsi = 1. / rho
         clebsch_factor = sign_torflux
      else if (radial_coordinate_switch==radial_coordinate_minuspsit) then
         dxdpsi = -sign_torflux / rho
         dydalpha = rho
         dpsidpsit = -1
         drhodpsi = -sign_torflux / rho
         clebsch_factor = -1
      else if (radial_coordinate_switch==radial_coordinate_r) then
         dxdpsi = 1  
         dydalpha = rho
         dpsidpsit = sign_torflux / rho
         drhodpsi = 1
         clebsch_factor = sign_torflux * rho
      end if

      ! Define the inverse or derived relations
      dxdpsit = dxdpsi * dpsidpsit
      dpsidx = 1. / dxdpsi
      dpsitdx = 1. / dxdpsit

      ! Now that we have defined our coordinates, transform psi_t to psi or x
      grad_alpha_grad_psi = dpsidpsit * grad_alpha_grad_psit
      grad_psi_grad_psi = dpsidpsit * dpsidpsit * grad_psit_grad_psit
      grad_y_grad_x = dxdpsi * dydalpha * grad_alpha_grad_psi
      grad_x_grad_x = grad_psi_grad_psi * dxdpsi**2
      grad_y_grad_y = grad_alpha_grad_alpha * dydalpha**2
      gds23 = gds23_alphapsit * dydalpha * dydalpha * dxdpsit
      gds24 = gds24_alphapsit * dydalpha * dxdpsit * dxdpsit
      gds25 = gds25_alphapsit * dydalpha * dydalpha * dxdpsit
      gds26 = gds26_alphapsit * dydalpha * dxdpsit * dxdpsit
      B_times_gradB_dot_grady = B_times_gradB_dot_gradalpha * dydalpha
      B_times_gradB_dot_gradx = B_times_gradB_dot_gradpsit * dxdpsit
      B_times_kappa_dot_grady = B_times_kappa_dot_gradalpha * dydalpha
      B_times_kappa_dot_gradx = B_times_kappa_dot_gradpsit * dxdpsit

      ! <grad_x> = |∇x|
      ! <gradx_dot_grady> = ∇x . ∇y
      ! <grady_dot_grady> = |∇y|^2
      ! <gradx_dot_gradx> = |∇x|^2
      grad_x = sqrt(abs(grad_x_grad_x))
      grady_dot_grady = grad_y_grad_y
      gradx_dot_grady = grad_y_grad_x
      gradx_dot_gradx = (grad_x)**2
      
      ! For the momentum flux we need (R^2/B^2) ∇ζ . ∇y and (R^2/B^2) ∇ζ . ∇x
      gradzeta_gradx_R2overB2 = gradzeta_gradpsit_R2overB2 * dxdpsit
      gradzeta_grady_R2overB2 = gradzeta_gradalpha_R2overB2 * dydalpha

      ! We want |ds/dx|*sqrt((dR/ds)^2+(dZ/ds)^2) so do dψ̃t/dx̃ * |ds̃/dψ̃t|*sqrt((dR/ds)^2+(dZ/ds)^2)  
      x_displacement_fac = dpsitdx*psit_displacement_fac

      ! Calculate the <twist_and_shift_geo_fac> needed for the boundary conditions
      call calculate_twist_and_shift_geo_fac()

      ! Arrays as a function of <psi_t> are not needed in stella since we 
      ! only need arrays as a function of <x> and <y>
      call deallocate_temporary_arrays

      ! Test FFS implementation by setting all geometric coefficients
      ! to their values at a given alpha; i.e., make the system axisymmetric
      if (const_alpha_geo) call set_ffs_geo_coefs_constant(nalpha)

      ! Variables for extra stella options (e.g., radial variation, ...)
      ! TODO: WARNING: FLAG: <btor> is not defined for VMEC
      ! It is used to calculate the momentum flux
      btor = -1000.
      drhodpsip = -1000.
      drhodpsip_psi0 = -1000.
      bmag_psi0 = bmag 
      
   contains

      !**********************************************************************
      !                     BOUNDARY CONDITIONS IN VMEC                     !
      !**********************************************************************
      ! Calculate <twist_and_shift_geo_fac> = dkx/dky * jtwist 
      ! Minus its sign gives the direction of the shift in kx for the twist-and-shift BC
      ! Note that the default BC are the unconnected BC in the <zgrid> module 
      ! 
      ! Regardless of the choice of coordinates we have 
      !     kψ = (dx/dψ) kx
      !     kα = (dy/dα) ky
      !     hat{s} = -(rho/iota)(diota/drho) 
      ! 
      ! For VMEC we use
      !     alpha = theta - iota*zeta
      !     z = zeta = 2*pi*P = 2*pi*field_period_ratio = 2*pi*nfield_periods/nfp
      ! 
      ! Standard twist-and-shift boundary conditions  
      !     dkpsi/dkalpha * jtwist = 2*pi*P*(diota/dψ) = 2*pi*P * (drho/dψ) (diota/drho) 
      !     dkx/dky * jtwist = -2*pi*P * (dψ/dx) (dy/dα) * (drho/dψ) * hat{s} (iota/rho)
      ! 
      ! Stellarator symmetric twist-and-shift boundary conditions  
      !     dkpsi/dkalpha * jtwist = -2 (∇ψ . ∇α) / |∇ψ|^2   (at the last z-point)
      !     dkx/dky * jtwist = -2  * (dψ/dx) (dy/dα) * (∇ψ . ∇α) / |∇ψ|^2
      !**********************************************************************
      subroutine calculate_twist_and_shift_geo_fac()
  
         use grids_z, only: boundary_option_switch, boundary_option_linked_stellarator 
         use constants, only: pi
         use debug_flags, only: print_extra_info_to_terminal
         implicit none  

         !---------------------------------------------------------------------- 

         ! Print the boundary condition information to the output file
         if (print_extra_info_to_terminal) then 
            write (*, '(A)') "############################################################"
            write (*, '(A)') "                     BOUNDARY CONDITIONS"
            write (*, '(A)') "############################################################"
         end if
         select case (boundary_option_switch)

         ! Stellarator symmetric twist-and-shift BC (see Martin et al, 2018)
         case (boundary_option_linked_stellarator)

            ! dkx/dky * jtwist = -2  * (dψ/dx) (dy/dα) * (∇ψ . ∇α) / |∇ψ|^2
            if (print_extra_info_to_terminal) then; write (*, *) ' '; write (*, *) 'Stellarator symmetric twist and shift BC selected'; end if 
            twist_and_shift_geo_fac = -2. * dpsidx * dydalpha * (grad_alpha_grad_psi(1, nzgrid)) / (grad_psi_grad_psi(1, nzgrid))

         ! Standard twist-and-shift boundary conditions or unconnected BC (then twist_and_shift_geo_fac doesn't matter)
         case default

            ! dkx/dky * jtwist = -2*pi*P * (dψ/dx) (dy/dα) * (drho/dψ) * hat{s} (iota/rho)
            if (print_extra_info_to_terminal) then; write (*, *) ' '; write (*, *) 'Standard twist and shift BC selected'; end if 
            twist_and_shift_geo_fac = -2. * pi * field_period_ratio * dpsidx * dydalpha * drhodpsi * shat * iota / rho  

         end select

         ! If ∇x . ∇y is very small at the ends of the flux tube, use periodic boundary conditions
         grad_x_grad_y_end = grad_y_grad_x(1, nzgrid)

         ! Print this factor to the output file
         if (print_extra_info_to_terminal) write (*, *) 'twist_and_shift_geo_fac: ', twist_and_shift_geo_fac; write (*, *) ' '
   
      end subroutine calculate_twist_and_shift_geo_fac
      

      !**********************************************************************
      !                      ALLOCATE TEMPORARY ARRAYS                      !
      !**********************************************************************
      subroutine allocate_temporary_arrays(nalpha, nzgrid)

         implicit none

         integer, intent(in) :: nalpha, nzgrid

         ! We need the following arrays to construct <grady_dot_grady>, <gradx_dot_grady>, <gradx_dot_gradx>, ...
         allocate (psit_displacement_fac(nalpha, -nzgrid:nzgrid))
         allocate (grad_alpha_grad_alpha(nalpha, -nzgrid:nzgrid))
         allocate (grad_alpha_grad_psi(nalpha, -nzgrid:nzgrid))
         allocate (grad_alpha_grad_psit(nalpha, -nzgrid:nzgrid))
         allocate (grad_psit_grad_psit(nalpha, -nzgrid:nzgrid))
         allocate (grad_psi_grad_psi(nalpha, -nzgrid:nzgrid))
         allocate (B_times_gradB_dot_gradalpha(nalpha, -nzgrid:nzgrid))
         allocate (B_times_kappa_dot_gradalpha(nalpha, -nzgrid:nzgrid))
         allocate (B_times_gradB_dot_gradpsit(nalpha, -nzgrid:nzgrid))
         allocate (B_times_kappa_dot_gradpsit(nalpha, -nzgrid:nzgrid))
         allocate (gds23_alphapsit(nalpha, -nzgrid:nzgrid))
         allocate (gds24_alphapsit(nalpha, -nzgrid:nzgrid))
         allocate (gds25_alphapsit(nalpha, -nzgrid:nzgrid))
         allocate (gds26_alphapsit(nalpha, -nzgrid:nzgrid))
         allocate (grad_x_grad_x(nalpha, -nzgrid:nzgrid))
         allocate (grad_y_grad_y(nalpha, -nzgrid:nzgrid))
         allocate (grad_y_grad_x(nalpha, -nzgrid:nzgrid))
         
         ! Needed for the momentum flux diagnostic for non-axisymmetric devices
         allocate (gradzeta_gradpsit_R2overB2(nalpha, -nzgrid:nzgrid)); gradzeta_gradpsit_R2overB2 = 0.0 ! Temp array
         allocate (gradzeta_gradalpha_R2overB2(nalpha, -nzgrid:nzgrid)); gradzeta_gradalpha_R2overB2 = 0.0 ! Temp array

      end subroutine allocate_temporary_arrays

      !**********************************************************************
      !                     DEALLOCATE TEMPORARY ARRAYS                     !
      !**********************************************************************
      subroutine deallocate_temporary_arrays

         implicit none

         deallocate (psit_displacement_fac, grad_alpha_grad_alpha)
         deallocate (grad_alpha_grad_psi, grad_psi_grad_psi)
         deallocate (grad_alpha_grad_psit, grad_psit_grad_psit)
         deallocate (B_times_gradB_dot_gradalpha, B_times_kappa_dot_gradalpha)
         deallocate (B_times_gradB_dot_gradpsit, B_times_kappa_dot_gradpsit)
         deallocate (gds23_alphapsit, gds24_alphapsit)
         deallocate (gds25_alphapsit, gds26_alphapsit)
         deallocate (grad_y_grad_y, grad_x_grad_x, grad_y_grad_x)
         deallocate (gradzeta_gradpsit_R2overB2, gradzeta_gradalpha_R2overB2)

      end subroutine deallocate_temporary_arrays

   end subroutine get_geometry_arrays_from_VMEC

   !=========================================================================
   !========================== MILLER EQUILIBRIUM  ==========================
   !=========================================================================
   ! The default option uses the following coordinates:
   !     psi = psi_poloidal = enclosed poloidal flux divided by 2pi 
   !     x = q/(r*B_zeta) (psip - psip0)
   !     y = (1/Bref)(dpsip/dr) (alpha - alpha0)
   !     r/q = (1/Bref)(dpsip/dr)
   !     alpha = zeta - q * theta
   !     B = - nabla psi x nabla alpha
   ! 
   ! For <radial_variation> = True, we use <q_as_x> = True
   ! See "A novel approach to radially global gyrokinetic simulation 
   ! using the flux-tube code stella" by D. A. St-Onge
   !     psi = q = psit/psip
   !     x = (1/Bref)(dpsip/dr) (q - q0)
   !     y = (1/Bref)(dpsip/dr) (alpha - alpha0)
   !     alpha =  zeta - q * theta
   !     B = - nabla psi x nabla alpha
   !======================================================================
   subroutine get_geometry_arrays_from_Miller(nalpha)

      ! Grids
      use grids_z, only: zed, nzed, nzgrid
      use grids_z, only: zed_equal_arc
      use constants, only: pi
      
      ! Miller geometry
      use geometry_miller, only: read_miller_parameters
      use geometry_miller, only: get_miller_geometry
      
      ! Radial variation
      use geometry_miller, only: communicate_parameters_multibox
      use geometry_inputprofiles_interface, only: read_inputprof_geo

      implicit none

      integer, intent(in) :: nalpha 
      real :: dpsipdrho, dpsipdrho_psi0

      !---------------------------------------------------------------------- 
     
      ! Read the <millergeo_parameters> namelist in the input file  
      ! <nzed> and <nzgrid> are inputs, and <geo_surf> is returned 
      ! which is a dictionary that contains all the geometry information 
      if (debug) write (*, *) 'geometry::Miller::read_miller_parameters'
      call read_miller_parameters(nzed, nzgrid, geo_surf)

      ! Allocate geometry arrays for stella
      if (debug) write (*, *) 'geometry::Miller::allocate_arrays'
      call allocate_arrays(nalpha, nzgrid)

      ! Overwrite parameters from the input file with those from the input.profiles file
      ! We use <rhoc> from the input file to select the surface
      if (geo_option_switch==geo_option_inputprof) then
         call read_inputprof_geo(geo_surf)
      end if

      ! Multi box simulations that include radial variation
      ! WARNING: I took this line out of a seperate loop, perhaps 
      ! <geo_option_multibox> is broken now ... check with an older version
      if (geo_option_switch==geo_option_multibox) then
         call communicate_parameters_multibox(surf=geo_surf)
      end if

      ! Call the <geometry_miller.f90> module to calculate the geometric coefficients 
      ! needed by stella, based on the local Miller parameters. For Miller geometries
      ! each field line <alpha> has the same geometry, hence we will only pass on the 
      ! ialpha=1 arrays to the get_miller_geometry() routine. Note that in Miller z = theta.
      if (debug) write (*, *) 'geometry::Miller::get_miller_geometry'
      call get_miller_geometry(nzgrid, zed, zed_equal_arc, &
                dpsipdrho, dpsipdrho_psi0, dIdrho, grho(1, :), bmag(1, :), bmag_psi0(1, :), &
                grady_dot_grady(1, :), gradx_dot_grady(1, :), gradx_dot_gradx(1, :), gds23(1, :), gds24(1, :), b_dot_gradz(1, :), &
                B_times_gradB_dot_gradx(1, :), B_times_gradB_dot_grady(1, :), &
                B_times_kappa_dot_gradx(1, :), B_times_kappa_dot_grady(1, :), &
                dBdrho, d2Bdrdth, d_bdotgradz_drho, btor, rmajor, &
                dcvdrift0drho(1, :), dcvdriftdrho(1, :), dgbdrift0drho(1, :), dgbdriftdrho(1, :), &
                d_gradydotgrady_drho(1, :), d_gradxdotgrady_drho(1, :), d_gradxdotgradx_drho(1, :), djacdrho(1, :))
      if (debug) write (*, *) 'geometry::Miller::get_miller_geometry_finished'

      ! <drhodpsip> = drho/d(psip/a^2*Bref) = (a^2*Bref) * drho/dpsip = (a*Bref) * dr/dpsip
      drhodpsip = 1./dpsipdrho 
      drhodpsip_psi0 = 1./dpsipdrho_psi0
      drhodpsi = drhodpsip

      ! Miller parameters use psi = psip; dx/dpsi = dr/dpsip = q/(r*Br) and rho = r/a 
      ! <dxdpsi> = (rhor/a)(dx̃/dψ̃) = (rhor/a)(d(x/rhor)/d(psi/a^2Br) = a*Br*(dx/dpsi) = a*q/r = q/rho 
      dxdpsi = geo_surf%qinp_psi0 / geo_surf%rhoc_psi0 

      ! <dydalpha> = (rhor/a)(d(y/rhor)/dalpha) = (1/a)(dy/dalpha) = 1/(a*Bref) * (dpsi/dr) 
      ! <dpsipdrho> = dψ̃/drho = d(psip/(a^2*Br))/d(r/a) = 1/(a*Bref) * dpsip/dr 
      dydalpha = dpsipdrho

      ! For psi = psi_p we have B = - ∇ψ x ∇α 
      clebsch_factor = -1.

      ! Since the profile gradients are given with respect to r we also need (dr/dx) evaluated at the field line
      ! We defined x = q0/(r0*Bref) (psip - psip0) and r^2 = 2*q0*psi_p/Bref
      ! dr/dpsip = sqrt(2*q0/Bref) dsqrt(psip)/dpsip = sqrt(2*q0/Bref) (1/2) (1/sqrt(psip0)) = sqrt(q0/2*Bref*psip0)
      ! dx/dpsip = q0/(r0*Bref) dpsip/dpsip = q0/(r0*Bref) = q0/Bref sqrt(Bref/2*q0*psip0)
      ! dr/dx = dr/dpsip * dpsip/dx = sqrt(q0/2*Bref*psip0) * Bref/q0 * sqrt(2*q0*psip0/Bref) = 1
      ! TODO: drdx = 1

      ! <grad_x> = |grad x| = dx/drho * |grad rho| = dx/dpsi * dpsi/drho * |grad rho| 
      !          = a*Br*(dx/dpsi) * 1/(a^2*Bref) * dpsip/drho * a |grad rho| = <dxdpsi> * <dpsipdrho> * <grho>
      ! <dpsipdrho> = 1/(a*Bref) * dpsip/dr = 1/(a^2*Bref) * dpsip/drho   
      ! <dxdpsi> = a*Br*(dx/dpsi)
      ! <grho> = a |grad rho| 
      grad_x = dxdpsi * dpsipdrho_psi0 * grho

      ! abs(<twist_and_shift_geo_fac>) is dkx/dky * jtwist
      ! minus its sign gives the direction of the shift in kx for the twist-and-shift BC
      ! For x = q/(r*B_zeta) (psip - psip0) --> dkx/dky * jtwist = 2*pi*shat
      twist_and_shift_geo_fac = 2.0 * pi * geo_surf%shat_psi0 

      ! <aref> and <bref> should not be needed, so set to 1
      aref = 1.0; bref = 1.0

      ! We use zed = theta, so zeta = q*theta = q*zed
      zeta(1, :) = zed * geo_surf%qinp_psi0

      ! If <radial_variation> = True than <q_as_x> = True
      ! x = (1/Bref)(dpsip/dr) (q - q0) and psi = q
      ! hat{s} = r/q (dq/dr) = rho/q (dq/drho)
      ! B = - (dpsi_p/dq) ∇ψ x ∇α 
      ! <dqdrho> = dq/drho = hat{s} * q/rho 
      ! <dpsipdrho> = dψ̃/drho = d(psip/(a^2*Br))/d(r/a) = 1/(a*Bref) * dpsip/dr 
      ! <dxdpsi> = (rhor/a)(dx̃/dq) = (rhor/a)(d(x/rhor)/dq = (1/a)*(dx/dq) = 1/(a*Bref) * dpsip/dr  
      ! <twist_and_shift_geo_fac> = dkx/dky * jtwist = 2*pi 
      ! <clebsch_factor> = - d(psip/(a^2*Br))/dq = - 1/(a^2*Bref) (dpsi_p/q) = - 1/(a^2*Bref) (dpsi_p/drho) (drho/dq) 
      !                 = - 1/(a*Bref) (dpsi_p/dr) (drho/dq) = - <dpsipdrho> * 1/<dqdrho>
      ! <grad_x> = |grad x| = dx/drho * |grad rho| = dx/dq * dq/drho * |grad rho| 
      !          = (1/a) dx/dq * dq/drho * a |grad rho| = <dxdpsi> * <dqdrho> * <grho>
      ! WARNING: <grad_x> was wrong in pervious version for <q_as_x> = True, but it was not used
      if (q_as_x) then
         dqdrho = geo_surf%shat_psi0 * geo_surf%qinp_psi0 / geo_surf%rhoc_psi0 
         dxdpsi = dpsipdrho_psi0
         dydalpha = dpsipdrho_psi0
         clebsch_factor = - dpsipdrho_psi0 * (1/dqdrho)
         twist_and_shift_geo_fac = 2.0 * pi
         grad_x = dxdpsi * dqdrho * grho
         drhodpsi = 1/dqdrho
      end if

      ! For Miller each field line <alpha> has the same geometry so
      ! ensure that all arrays are filled with the ialpha = 1 information
      if (debug) write (*, *) 'geometry::Miller::spread_geometry'
      bmag = spread(bmag(1, :), 1, nalpha)
      grho = spread(grho(1, :), 1, nalpha) 
      bmag_psi0 = spread(bmag_psi0(1, :), 1, nalpha)
      grady_dot_grady = spread(grady_dot_grady(1, :), 1, nalpha)
      gradx_dot_grady = spread(gradx_dot_grady(1, :), 1, nalpha)
      gradx_dot_gradx = spread(gradx_dot_gradx(1, :), 1, nalpha)
      gds23 = spread(gds23(1, :), 1, nalpha)
      gds24 = spread(gds24(1, :), 1, nalpha)
      B_times_gradB_dot_gradx = spread(B_times_gradB_dot_gradx(1, :), 1, nalpha)
      B_times_gradB_dot_grady = spread(B_times_gradB_dot_grady(1, :), 1, nalpha)
      B_times_kappa_dot_gradx = spread(B_times_kappa_dot_gradx(1, :), 1, nalpha)
      B_times_kappa_dot_grady = spread(B_times_kappa_dot_grady(1, :), 1, nalpha)
      dcvdrift0drho = spread(dcvdrift0drho(1, :), 1, nalpha)
      dcvdriftdrho = spread(dcvdriftdrho(1, :), 1, nalpha)
      dgbdrift0drho = spread(dgbdrift0drho(1, :), 1, nalpha)
      dgbdriftdrho = spread(dgbdriftdrho(1, :), 1, nalpha)
      d_gradydotgrady_drho = spread(d_gradydotgrady_drho(1, :), 1, nalpha)
      d_gradxdotgrady_drho = spread(d_gradxdotgrady_drho(1, :), 1, nalpha)
      d_gradxdotgradx_drho = spread(d_gradxdotgradx_drho(1, :), 1, nalpha)
      djacdrho = spread(djacdrho(1, :), 1, nalpha)
      b_dot_gradz = spread(b_dot_gradz(1, :), 1, nalpha)
      zeta = spread(zeta(1, :), 1, nalpha)
      b_dot_gradz_avg = b_dot_gradz(1, :) ! For Miller <b_dot_gradz> = <b_dot_gradz_avg> for ialpha=1

      ! For the momentum flux we need (R^2/B^2) ∇ζ . ∇y and (R^2/B^2) ∇ζ . ∇x
      ! For Miller or axi-symmetric devices we have: ∇ζ . ∇ψp = 0 and ∇ζ . ∇α = ∇ζ . ∇ζ = (1/R^2) 
      !      (R^2/B^2) * ∇ζ . ∇α * (dy/dα) = (R^2/B^2) (1/R^2) (dy/dα) = (1/B^2) (dy/dα) = geo_surf%rhoc / (geo_surf%qinp * bmag**2)
      gradzeta_grady_R2overB2 = geo_surf%rhoc / (geo_surf%qinp * bmag**2)
      gradzeta_gradx_R2overB2 = 0.0

      ! For the momentum flux we need R^2 * b . ∇ζ
      ! Note that in Miller z = theta, so b_dot_gradz = b_dot_grad_theta
      !      R^2 * b . ∇ζ = R^2 * (1/B) (∇ζ x ∇ψ + I ∇ζ) . ∇ζ =  R^2/B * I * ∇ζ . ∇ζ
      !                  = R^2/B * I * (1/R^2) = I/B = (R/B) (I/R) = (R/B) * Btor
      b_dot_gradzeta_RR = geo_surf%rmaj * spread(btor, 1, nalpha) / bmag 
      if (debug) write (*, *) 'geometry::Miller::get_geometry_arrays_from_Miller_finished'

   end subroutine get_geometry_arrays_from_Miller

   !=========================================================================
   !========================== ZPINCH EQUILIBRIUM  ==========================
   !=========================================================================
   subroutine get_geometry_arrays_from_zpinch (nalpha)

      use zpinch, only: get_zpinch_geometry_coefficients
      use grids_z, only: nzgrid
      
      implicit none

      integer, intent (in) :: nalpha
      real :: dpsipdrho, dpsipdrho_psi0

      !-------------------------------------------------------------------------

      ! Allocate the arrays
      call allocate_arrays(nalpha, nzgrid)

      ! Calculate the geometric coefficients for a z-pinch magnetic equilibrium
      call get_zpinch_geometry_coefficients(nzgrid, bmag(1, :), b_dot_gradz_avg, grho(1, :), geo_surf, &
          grady_dot_grady(1, :), gradx_dot_grady(1, :), gradx_dot_gradx(1, :), B_times_gradB_dot_gradx(1, :), &
          B_times_gradB_dot_grady(1, :), B_times_kappa_dot_gradx(1, :), B_times_kappa_dot_grady(1, :), btor, rmajor)

      ! Note that <b_dot_gradz> is the alpha-dependent b . grad z,
      ! and <b_dot_gradz_avg> is the constant-in-alpha part of it.
      ! For a z-pinch, <b_dot_gradz> is independent of alpha.
      b_dot_gradz(1, :) = b_dot_gradz_avg

      ! Effectively choose psi = x * B * a_ref = x * B * L_B
      dpsipdrho = 1.0
      dpsipdrho_psi0 = 1.0
      bmag_psi0 = bmag
      
      ! Note that psi here is meaningless
      drhodpsi = 1./dpsipdrho
      drhodpsip_psi0 = 1./dpsipdrho_psi0
      dxdpsi = 1.0
      dydalpha = 1.0
      sign_torflux = -1
      clebsch_factor = sign_torflux
      grad_x = sqrt(gradx_dot_gradx)
      
      ! There is no magnetic shear in the z-pinch and thus no need for twist-and-shift
      twist_and_shift_geo_fac = 1.0
      
      ! The reference quantities aref and bref should not be needed, so set to 1
      aref = 1.0; bref = 1.0
      
      ! The array zeta should not be needed
      zeta(1, :) = 0.0
     
   end subroutine get_geometry_arrays_from_zpinch
   
   !=========================================================================
   !================== OVERWRITE GEOMETRIC COEFFICIENTS  ====================
   !=========================================================================
   subroutine overwrite_selected_geometric_coefficients(nalpha)

      use file_utils, only: get_unused_unit
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      integer, intent(in) :: nalpha
      
      ! Local variables
      integer :: geofile_unit
      character(100) :: dum_char
      real :: dum_real
      integer :: ia, iz
      real :: bmag_file, b_dot_gradzeta_file
      real :: grady_dot_grady_file, gradx_dot_grady_file, gradx_dot_gradx_file, gds23_file, gds24_file
      real :: B_times_gradB_dot_grady_file, B_times_kappa_dot_grady_file, B_times_gradB_dot_gradx_file

      !-------------------------------------------------------------------------

      ! Open the geometry file
      call get_unused_unit(geofile_unit)
      open (geofile_unit, file=trim(geometry_file), status='old', action='read')

      ! Deal with the first lines
      read (geofile_unit, fmt=*) dum_char
      read (geofile_unit, fmt=*) dum_char
      read (geofile_unit, fmt=*) dum_char

      ! Overwrite the following geometric quantities:
      !   - bmag, b_dot_gradzeta, grady_dot_grady, gradx_dot_grady, gradx_dot_gradx, gds23, gds24,
      !   - B_times_gradB_dot_grady, B_times_kappa_dot_grady, B_times_gradB_dot_gradx, and B_times_kappa_dot_gradx
      do ia = 1, nalpha
         do iz = -nzgrid, nzgrid
         
            ! Read the (alpha,z) point in the geometry file
            read (geofile_unit, fmt='(13e12.4)') dum_real, dum_real, dum_real, bmag_file, b_dot_gradzeta_file, &
               grady_dot_grady_file, gradx_dot_grady_file, gradx_dot_gradx_file, gds23_file, &
               gds24_file, B_times_gradB_dot_grady_file, B_times_kappa_dot_grady_file, B_times_gradB_dot_gradx_file
               
            ! If an overwrite flag is turned on, overwrite the geometric quantity
            if (overwrite_bmag) bmag(ia, iz) = bmag_file
            if (overwrite_grady_dot_grady) grady_dot_grady(ia, iz) = grady_dot_grady_file
            if (overwrite_gradx_dot_grady) gradx_dot_grady(ia, iz) = gradx_dot_grady_file
            if (overwrite_gradx_dot_gradx) gradx_dot_gradx(ia, iz) = gradx_dot_gradx_file
            if (overwrite_gds23) gds23(ia, iz) = gds23_file
            if (overwrite_gds24) gds24(ia, iz) = gds24_file
            if (overwrite_B_times_gradB_dot_grady) B_times_gradB_dot_grady(ia, iz) = B_times_gradB_dot_grady_file
            if (overwrite_B_times_kappa_dot_grady) B_times_kappa_dot_grady(ia, iz) = B_times_kappa_dot_grady_file
            if (overwrite_B_times_gradB_dot_gradx) B_times_gradB_dot_gradx(ia, iz) = B_times_gradB_dot_gradx_file
            
            ! Assume we are only reading in for a single alpha. 
            ! Usually, b_dot_gradzeta is the average of all b_dot_gradz values.
            if (overwrite_b_dot_gradzeta) b_dot_gradz(1, iz) = b_dot_gradzeta_file
            
         end do
      end do
      
      ! For any static ideal MHD equilibrium B × κ · ∇ψ = b × ∇B · ∇ψ
      B_times_kappa_dot_gradx = B_times_gradB_dot_gradx

      ! Close the geometry file
      close (geofile_unit)

   end subroutine overwrite_selected_geometric_coefficients

   !=========================================================================
   !============= MAKE GEOMETRIC COEFFICIENTS CONSTANT IN ALPHA =============
   !=========================================================================
   subroutine set_ffs_geo_coefs_constant(nalpha)

      implicit none

      integer, intent(in) :: nalpha

      !-------------------------------------------------------------------------

      call set_coef_constant(B_times_gradB_dot_gradx, nalpha)
      call set_coef_constant(B_times_kappa_dot_gradx, nalpha)
      call set_coef_constant(B_times_gradB_dot_grady, nalpha)
      call set_coef_constant(B_times_kappa_dot_grady, nalpha)
      call set_coef_constant(grad_x, nalpha)
      call set_coef_constant(grho, nalpha)
      call set_coef_constant(bmag, nalpha)
      call set_coef_constant(bmag_psi0, nalpha)
      call set_coef_constant(grady_dot_grady, nalpha)
      call set_coef_constant(gradx_dot_grady, nalpha)
      call set_coef_constant(gradx_dot_gradx, nalpha)
      call set_coef_constant(gds23, nalpha)
      call set_coef_constant(gds24, nalpha)
      call set_coef_constant(gds25, nalpha)
      call set_coef_constant(gds26, nalpha)
      call set_coef_constant(theta_vmec, nalpha)
      call set_coef_constant(x_displacement_fac, nalpha)
      call set_coef_constant(zeta, nalpha)
      call set_coef_constant(b_dot_gradz, nalpha)
      call set_coef_constant(gradzeta_gradx_R2overB2, nalpha)
      call set_coef_constant(gradzeta_grady_R2overB2, nalpha)
      call set_coef_constant(b_dot_gradzeta_RR, nalpha)

   end subroutine set_ffs_geo_coefs_constant

   subroutine set_coef_constant(coef, nalpha)

      use grids_z, only: nzgrid

      implicit none

      real, dimension(:, -nzgrid:), intent(in out) :: coef
      integer, intent(in) :: nalpha

      coef = spread(coef(1, :), 1, nalpha)

   end subroutine set_coef_constant

   !============================================================================ 
   !============================ ALLOCATE ARRAYS  ==============================
   !============================================================================
   subroutine allocate_arrays(nalpha, nzgrid)

      implicit none

      integer, intent(in) :: nalpha, nzgrid

      !-------------------------------------------------------------------------

      ! It is possible to call <init_geometry> multiple times so we do not want to
      ! create more copies of the geometric arrays if they already exist, hence
      ! we need to use if (.not. allocated(...)) statements
      if (.not. allocated(bmag)) allocate (bmag(nalpha, -nzgrid:nzgrid)); bmag = 0.0 
      if (.not. allocated(bmag_psi0)) allocate (bmag_psi0(nalpha, -nzgrid:nzgrid)); bmag_psi0 = 0.0
      if (.not. allocated(grady_dot_grady)) allocate (grady_dot_grady(nalpha, -nzgrid:nzgrid)); grady_dot_grady = 0.0
      if (.not. allocated(gradx_dot_grady)) allocate (gradx_dot_grady(nalpha, -nzgrid:nzgrid)); gradx_dot_grady = 0.0
      if (.not. allocated(gradx_dot_gradx)) allocate (gradx_dot_gradx(nalpha, -nzgrid:nzgrid)); gradx_dot_gradx = 0.0
      if (.not. allocated(gds23)) allocate (gds23(nalpha, -nzgrid:nzgrid)); gds23 = 0.0
      if (.not. allocated(gds24)) allocate (gds24(nalpha, -nzgrid:nzgrid)); gds24 = 0.0
      if (.not. allocated(gds25)) allocate (gds25(nalpha, -nzgrid:nzgrid)); gds25 = 0.0
      if (.not. allocated(gds26)) allocate (gds26(nalpha, -nzgrid:nzgrid)); gds26 = 0.0
      if (.not. allocated(d_gradydotgrady_drho)) allocate (d_gradydotgrady_drho(nalpha, -nzgrid:nzgrid)); d_gradydotgrady_drho = 0.0
      if (.not. allocated(d_gradxdotgrady_drho)) allocate (d_gradxdotgrady_drho(nalpha, -nzgrid:nzgrid)); d_gradxdotgrady_drho = 0.0
      if (.not. allocated(d_gradxdotgradx_drho)) allocate (d_gradxdotgradx_drho(nalpha, -nzgrid:nzgrid)); d_gradxdotgradx_drho = 0.0
      if (.not. allocated(B_times_gradB_dot_grady)) allocate (B_times_gradB_dot_grady(nalpha, -nzgrid:nzgrid)); B_times_gradB_dot_grady = 0.0
      if (.not. allocated(B_times_gradB_dot_gradx)) allocate (B_times_gradB_dot_gradx(nalpha, -nzgrid:nzgrid)); B_times_gradB_dot_gradx = 0.0
      if (.not. allocated(B_times_kappa_dot_grady)) allocate (B_times_kappa_dot_grady(nalpha, -nzgrid:nzgrid)); B_times_kappa_dot_grady = 0.0
      if (.not. allocated(B_times_kappa_dot_gradx)) allocate (B_times_kappa_dot_gradx(nalpha, -nzgrid:nzgrid)); B_times_kappa_dot_gradx = 0.0
      if (.not. allocated(dgbdriftdrho)) allocate (dgbdriftdrho(nalpha, -nzgrid:nzgrid)); dgbdriftdrho = 0.0
      if (.not. allocated(dcvdriftdrho)) allocate (dcvdriftdrho(nalpha, -nzgrid:nzgrid)); dcvdriftdrho = 0.0
      if (.not. allocated(dgbdrift0drho)) allocate (dgbdrift0drho(nalpha, -nzgrid:nzgrid)); dgbdrift0drho = 0.0
      if (.not. allocated(dcvdrift0drho)) allocate (dcvdrift0drho(nalpha, -nzgrid:nzgrid)); dcvdrift0drho = 0.0
      if (.not. allocated(dbdzed)) allocate (dbdzed(nalpha, -nzgrid:nzgrid)); dbdzed = 0.0
      if (.not. allocated(theta_vmec)) allocate (theta_vmec(nalpha, -nzgrid:nzgrid)); theta_vmec = 0.0
      if (.not. allocated(jacob)) allocate (jacob(nalpha, -nzgrid:nzgrid)); jacob = 0.0
      if (.not. allocated(djacdrho)) allocate (djacdrho(nalpha, -nzgrid:nzgrid)); djacdrho = 0.0
      if (.not. allocated(grho)) allocate (grho(nalpha, -nzgrid:nzgrid)); grho = 0.0
      if (.not. allocated(grad_x)) allocate (grad_x(nalpha, -nzgrid:nzgrid)); grad_x = 0.0
      if (.not. allocated(dl_over_b)) allocate (dl_over_b(nalpha, -nzgrid:nzgrid)); dl_over_b = 0.0
      if (.not. allocated(d_dl_over_b_drho)) allocate (d_dl_over_b_drho(nalpha, -nzgrid:nzgrid)); d_dl_over_b_drho = 0.0
      if (.not. allocated(b_dot_gradz)) allocate (b_dot_gradz(nalpha, -nzgrid:nzgrid)); b_dot_gradz = 0.0
      if (.not. allocated(b_dot_gradz_avg)) allocate (b_dot_gradz_avg(-nzgrid:nzgrid)); b_dot_gradz_avg = 0.0
      if (.not. allocated(zed_eqarc)) allocate (zed_eqarc(-nzgrid:nzgrid)); zed_eqarc = 0.0
      if (.not. allocated(btor)) allocate (btor(-nzgrid:nzgrid)); btor = 0.0
      if (.not. allocated(rmajor)) allocate (rmajor(-nzgrid:nzgrid)); rmajor = 0.0
      if (.not. allocated(dBdrho)) allocate (dBdrho(-nzgrid:nzgrid)); dBdrho = 0.0
      if (.not. allocated(d2Bdrdth)) allocate (d2Bdrdth(-nzgrid:nzgrid)); d2Bdrdth = 0.0
      if (.not. allocated(d_bdotgradz_drho)) allocate (d_bdotgradz_drho(-nzgrid:nzgrid)); d_bdotgradz_drho = 0.0
      if (.not. allocated(alpha)) allocate (alpha(nalpha)); alpha = 0.0
      if (.not. allocated(zeta)) allocate (zeta(nalpha, -nzgrid:nzgrid)); zeta = 0.0
      if (.not. allocated(x_displacement_fac)) allocate (x_displacement_fac(nalpha, -nzgrid:nzgrid)); x_displacement_fac = 0.0
      
      ! Needed for the momentum flux diagnostic for non-axisymmetric devices
      if (.not. allocated(gradzeta_gradx_R2overB2)) allocate (gradzeta_gradx_R2overB2(nalpha, -nzgrid:nzgrid)); gradzeta_gradx_R2overB2 = 0.0
      if (.not. allocated(gradzeta_grady_R2overB2)) allocate (gradzeta_grady_R2overB2(nalpha, -nzgrid:nzgrid)); gradzeta_grady_R2overB2 = 0.0
      if (.not. allocated(b_dot_gradzeta_RR)) allocate (b_dot_gradzeta_RR(nalpha, -nzgrid:nzgrid)); b_dot_gradzeta_RR = 0.0 

   end subroutine allocate_arrays

   !============================================================================
   !============================= BROADCAST ARRAYS =============================
   !============================================================================
   subroutine broadcast_arrays

      use mp, only: broadcast

      implicit none

      ! Flags
      call broadcast(q_as_x)

      ! Switch between coordinates
      call broadcast(clebsch_factor)
      call broadcast(drhodpsi)
      call broadcast(drhodpsip)
      call broadcast(drhodpsip_psi0)
      call broadcast(dxdpsi)
      call broadcast(dydalpha)

      ! Z-grid
      call broadcast(zeta)
      call broadcast(alpha)
      call broadcast(dzetadz)
      call broadcast(twist_and_shift_geo_fac)
      call broadcast(grad_x_grad_y_end)

      ! Geometric variables
      call broadcast(rmajor)
   
      ! Factors needed in the gyrokinetic equation
      call broadcast(exb_nonlin_fac)
      call broadcast(exb_nonlin_fac_p)

      ! Factor needed in the fluxes
      call broadcast(flux_fac)

      ! Geometric arrays
      call broadcast(dIdrho)
      call broadcast(grho)
      call broadcast(grad_x)
      call broadcast(bmag)
      call broadcast(bmag_psi0)
      call broadcast(btor)
      call broadcast(b_dot_gradz)
      call broadcast(b_dot_gradz_avg)
      call broadcast(grady_dot_grady)
      call broadcast(gradx_dot_grady)
      call broadcast(gradx_dot_gradx)
      call broadcast(gds23)
      call broadcast(gds24)
      call broadcast(gds25)
      call broadcast(gds26)
      call broadcast(d_gradydotgrady_drho)
      call broadcast(d_gradxdotgrady_drho)
      call broadcast(d_gradxdotgradx_drho)
      call broadcast(B_times_gradB_dot_gradx)
      call broadcast(B_times_gradB_dot_grady)
      call broadcast(B_times_kappa_dot_gradx)
      call broadcast(B_times_kappa_dot_grady)
      call broadcast(dgbdrift0drho)
      call broadcast(dgbdriftdrho)
      call broadcast(dcvdrift0drho)
      call broadcast(dcvdriftdrho)
      call broadcast(dBdrho)
      call broadcast(d2Bdrdth)
      call broadcast(d_bdotgradz_drho)
      call broadcast(djacdrho)
      
      ! Arrays for the momentum flux
      call broadcast(gradzeta_gradx_R2overB2)
      call broadcast(gradzeta_grady_R2overB2)
      call broadcast(b_dot_gradzeta_RR)
   
      ! Geometric variables at the chosen radial location
      call broadcast(qinp)
      call broadcast(shat)
      call broadcast(geo_surf%rmaj)
      call broadcast(geo_surf%rgeo)
      call broadcast(geo_surf%kappa)
      call broadcast(geo_surf%kapprim)
      call broadcast(geo_surf%tri)
      call broadcast(geo_surf%triprim)
      call broadcast(geo_surf%rhoc)
      call broadcast(geo_surf%rhoc_psi0)
      call broadcast(geo_surf%dr)
      call broadcast(geo_surf%shift)
      call broadcast(geo_surf%qinp)
      call broadcast(geo_surf%qinp_psi0)
      call broadcast(geo_surf%shat)
      call broadcast(geo_surf%shat_psi0)
      call broadcast(geo_surf%betaprim)
      call broadcast(geo_surf%betadbprim)
      call broadcast(geo_surf%d2qdr2)
      call broadcast(geo_surf%d2psidr2)
      call broadcast(geo_surf%dpsitordrho)
      call broadcast(geo_surf%rhotor)
      call broadcast(geo_surf%psitor_lcfs)
      call broadcast(geo_surf%drhotordrho)

      ! Reference quantities
      call broadcast(aref)
      call broadcast(bref)

   end subroutine broadcast_arrays

   !============================================================================ 
   !============== COMMUNICATE GEOMETRY FOR A MULTIBOX SIMULATION ==============
   !============================================================================
   subroutine communicate_geo_multibox(l_edge, r_edge)

      use mp, only: proc0
      use geometry_miller, only: communicate_parameters_multibox

      implicit none

      real, intent(in) :: l_edge, r_edge

      !-------------------------------------------------------------------------

      if (proc0) then
         call communicate_parameters_multibox(geo_surf, gfac * l_edge, gfac * r_edge)
      end if

   end subroutine communicate_geo_multibox

   !============================================================================ 
   !============================== CALCULATE DZED ==============================
   !============================================================================
   ! given function f(z:-pi->pi), calculate z derivative
   ! second order accurate, with equal grid spacing assumed
   ! assumes periodic in z -- may need to change this in future
   !============================================================================
   subroutine get_dzed(nz, dz, f, df)

      implicit none

      integer, intent(in) :: nz
      real, dimension(-nz:), intent(in) :: dz, f
      real, dimension(-nz:), intent(out) :: df

      !-------------------------------------------------------------------------

      df(-nz + 1:nz - 1) = (f(-nz + 2:) - f(:nz - 2)) / (dz(:nz - 2) + dz(-nz + 1:nz - 1))

      ! TODO-GA hack to avoid non-periodicity in full-flux-surface case
      ! if (full_flux_annulus .and. .not. const_alpha_geo) then
      !   df(-nz) = (f(-nz + 1) - f(-nz)) / dz(-nz)
      !  df(nz) = (f(nz) - f(nz - 1)) / dz(nz - 1)
      !else
      ! assume periodicity in the B-field
      df(-nz) = (f(-nz + 1) - f(nz - 1)) / (dz(-nz) + dz(nz - 1))
      df(nz) = df(-nz)
      !end if

   end subroutine get_dzed

   !============================================================================ 
   !========================= CALCULATE b_dot_gradzeta EQARC ==========================
   !============================================================================
   subroutine get_b_dot_gradz_avg_eqarc(gp, z, dz, gp_eqarc)

      use constants, only: pi
      use grids_z, only: nzgrid

      implicit none

      real, dimension(-nzgrid:), intent(in) :: gp, z, dz
      real, intent(out) :: gp_eqarc

      !-------------------------------------------------------------------------

      ! First get int dz b . grad z
      call integrate_zed(dz, 1./gp, gp_eqarc)
      
      ! Then take (zmax-zmin)/int (dz b . gradz) to get b . grad z'
      gp_eqarc = (z(nzgrid) - z(-nzgrid)) / gp_eqarc

   end subroutine get_b_dot_gradz_avg_eqarc

   !============================================================================ 
   !=========================== CALCULATE ZED EQARC ============================
   !============================================================================
   subroutine get_zed_eqarc(gp, dz, z, gp_eqarc, z_eqarc)

      use grids_z, only: nzgrid

      implicit none

      real, dimension(-nzgrid:), intent(in) :: gp, dz, z
      real, intent(in) :: gp_eqarc
      real, dimension(-nzgrid:), intent(out) :: z_eqarc

      integer :: iz

      !-------------------------------------------------------------------------

      z_eqarc(-nzgrid) = z(-nzgrid)
      do iz = -nzgrid + 1, nzgrid
         call integrate_zed(dz(:iz), 1./gp(:iz), z_eqarc(iz))
      end do
      z_eqarc(-nzgrid + 1:) = z(-nzgrid) + z_eqarc(-nzgrid + 1:) * gp_eqarc

   end subroutine get_zed_eqarc

   !============================================================================ 
   !============================== INTEGRATE ZED ===============================
   !============================================================================
   ! Trapezoidal rule to integrate in zed.
   !============================================================================
   subroutine integrate_zed(dz, f, intf)

      use grids_z, only: nzgrid

      implicit none

      real, dimension(-nzgrid:), intent(in) :: dz
      real, dimension(-nzgrid:), intent(in) :: f
      real, intent(out) :: intf

      integer :: iz, iz_max

      !-------------------------------------------------------------------------

      iz_max = -nzgrid + size(dz) - 1

      intf = 0.
      do iz = -nzgrid + 1, iz_max
         intf = intf + dz(iz) * (f(iz - 1) + f(iz))
      end do
      intf = 0.5 * intf

   end subroutine integrate_zed

   !============================================================================ 
   !================================= X TO RHO =================================
   !============================================================================
   subroutine get_x_to_rho(llim, x_in, rho_out)

      use parameters_physics, only: rhostar

      implicit none

      integer, intent(in) :: llim
      real, dimension(:), intent(in) :: x_in
      real, dimension(:), intent(out) :: rho_out
      
      integer :: ix, ulim
      real :: a, b, c

      !-------------------------------------------------------------------------

      ulim = size(x_in) + llim - 1

      if (q_as_x) then
         a = 0.5 * geo_surf%d2qdr2 / dqdrho
         b = 1.0
         c = -rhostar / (dqdrho * dxdpsi)
      else
         a = 0.5 * geo_surf%d2psidr2 * drhodpsip
         b = 1.0
         c = -rhostar * drhodpsip / dxdpsi
      end if

      do ix = llim, ulim
         if (abs(4.0 * a * c * x_in(ix)) < 1.e-6) then
            rho_out(ix) = -(c * x_in(ix)) / b - a * (c * x_in(ix))**2 / b**3
         else
            rho_out(ix) = (-b + sqrt(b**2 - 4.*a * c * x_in(ix))) / (2.*a)
         end if
      end do

   end subroutine get_x_to_rho

   !============================================================================ 
   !============================== WRITE GEOMETRY ==============================
   !============================================================================
   subroutine write_geometric_coefficients(nalpha)

      use file_utils, only: open_output_file, close_output_file
      use grids_z, only: nzgrid, zed

      implicit none

      integer, intent(in) :: nalpha
      integer :: geometry_unit
      integer :: ia, iz

      !-------------------------------------------------------------------------

      ! Open a *.geometry text file
      call open_output_file(geometry_unit, '.geometry')

      ! Write the most important geometric variables to a text file
      write (geometry_unit, '(a1,11a14)') '#', 'rhoc', 'qinp', 'shat', 'rhotor', &
         'aref', 'bref', 'dxdpsi', 'dydalpha', 'exb_nonlin', 'flux_fac'
      write (geometry_unit, '(a1,11e14.4)') '#', geo_surf%rhoc, geo_surf%qinp, &
         geo_surf%shat, geo_surf%rhotor, aref, bref, dxdpsi, dydalpha, exb_nonlin_fac, flux_fac
      write (geometry_unit, *)

      ! Write the most important geometric arrays to a text file
      write (geometry_unit, '(15a12)') '# alpha', 'zed', 'zeta', 'bmag', 'bdot_grad_z', 'grady_dot_grady', 'gradx_dot_grady', 'gradx_dot_gradx', &
         'gds23', 'gds24', 'B_times_gradB_dot_grady', 'B_times_kappa_dot_grady', 'B_times_gradB_dot_gradx', 'bmag_psi0', 'btor'
      do ia = 1, nalpha
         do iz = -nzgrid, nzgrid
            write (geometry_unit, '(15e12.4)') alpha(ia), zed(iz), zeta(ia, iz), bmag(ia, iz), b_dot_gradz(ia, iz), &
               grady_dot_grady(ia, iz), gradx_dot_grady(ia, iz), gradx_dot_gradx(ia, iz), gds23(ia, iz), &
               gds24(ia, iz), B_times_gradB_dot_grady(ia, iz), B_times_kappa_dot_grady(ia, iz), B_times_gradB_dot_gradx(ia, iz), &
               bmag_psi0(ia, iz), btor(iz)
         end do
         write (geometry_unit, *)
      end do

      ! Close the *.geometry text file
      call close_output_file(geometry_unit)

   end subroutine write_geometric_coefficients

   !============================================================================ 
   !====================== FINISH INITIALIZING THE GEOMETRY ====================
   !============================================================================
   ! TODO - Is this really necessary, can Miller not deallocate everything immediatly?
   !============================================================================
   subroutine finish_init_geometry

      use mp, only: proc0
      use geometry_miller, only: finish_miller_geometry

      implicit none

      if (proc0) then
         select case (geo_option_switch)
         case (geo_option_local)
            call finish_miller_geometry
         case (geo_option_multibox)
            call finish_miller_geometry
         case (geo_option_inputprof)
            call finish_miller_geometry
         case (geo_option_vmec)
         end select
      end if

   end subroutine finish_init_geometry

   !============================================================================ 
   !============================ FINISH THE GEOMETRY ===========================
   !============================================================================
   subroutine finish_geometry

      implicit none

      if (allocated(zed_eqarc)) deallocate (zed_eqarc)
      if (allocated(grho)) deallocate (grho)
      if (allocated(grad_x)) deallocate (grad_x)
      if (allocated(bmag)) deallocate (bmag)
      if (allocated(bmag_psi0)) deallocate (bmag_psi0)
      if (allocated(btor)) deallocate (btor)
      if (allocated(rmajor)) deallocate (rmajor)
      if (allocated(dbdzed)) deallocate (dbdzed)
      if (allocated(jacob)) deallocate (jacob)
      if (allocated(djacdrho)) deallocate (djacdrho)
      if (allocated(b_dot_gradz)) deallocate (b_dot_gradz)
      if (allocated(b_dot_gradz_avg)) deallocate (b_dot_gradz_avg) 
      if (allocated(dl_over_b)) deallocate (dl_over_b)
      if (allocated(d_dl_over_b_drho)) deallocate (d_dl_over_b_drho)
      if (allocated(grady_dot_grady)) deallocate (grady_dot_grady)
      if (allocated(gradx_dot_grady)) deallocate (gradx_dot_grady)
      if (allocated(gradx_dot_gradx)) deallocate (gradx_dot_gradx)
      if (allocated(gds23)) deallocate (gds23)
      if (allocated(gds24)) deallocate (gds24)
      if (allocated(gds25)) deallocate (gds25)
      if (allocated(gds26)) deallocate (gds26)
      if (allocated(d_gradydotgrady_drho)) deallocate (d_gradydotgrady_drho)
      if (allocated(d_gradxdotgrady_drho)) deallocate (d_gradxdotgrady_drho)
      if (allocated(d_gradxdotgradx_drho)) deallocate (d_gradxdotgradx_drho)
      if (allocated(B_times_gradB_dot_grady)) deallocate (B_times_gradB_dot_grady)
      if (allocated(B_times_gradB_dot_gradx)) deallocate (B_times_gradB_dot_gradx)
      if (allocated(B_times_kappa_dot_grady)) deallocate (B_times_kappa_dot_grady)
      if (allocated(B_times_kappa_dot_gradx)) deallocate (B_times_kappa_dot_gradx)
      if (allocated(dgbdriftdrho)) deallocate (dgbdriftdrho)
      if (allocated(dcvdriftdrho)) deallocate (dcvdriftdrho)
      if (allocated(dgbdrift0drho)) deallocate (dgbdrift0drho)
      if (allocated(dcvdrift0drho)) deallocate (dcvdrift0drho)
      if (allocated(dBdrho)) deallocate (dBdrho)
      if (allocated(d2Bdrdth)) deallocate (d2Bdrdth)
      if (allocated(d_bdotgradz_drho)) deallocate (d_bdotgradz_drho)
      if (allocated(theta_vmec)) deallocate (theta_vmec)
      if (allocated(alpha)) deallocate (alpha)
      if (allocated(zeta)) deallocate (zeta)
      if (allocated(x_displacement_fac)) deallocate (x_displacement_fac)
      
      ! Arrays for the momentum flux 
      if (allocated(gradzeta_gradx_R2overB2)) deallocate (gradzeta_gradx_R2overB2)
      if (allocated(gradzeta_grady_R2overB2)) deallocate (gradzeta_grady_R2overB2)
      if (allocated(b_dot_gradzeta_RR)) deallocate (b_dot_gradzeta_RR)

      ! Only initialise once
      initialised_geometry = .false.

   end subroutine finish_geometry

end module geometry
