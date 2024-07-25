
!###############################################################################
!############################### STELLA GEOMETRY ###############################
!###############################################################################
! 
! Routines for calculating the geometry needed by stella.
! 
! This routine will call the <vmec_geometry> or <miller_geometry> modules.
! Which each uses specific (x, psi) coordinates. Nonetheless, stella is 
! completely general, it just needs geometric variables as a function of
! <psi> as well as the factor <dxdpsi> and <dpsitdpsi>.
! 
!###############################################################################


module stella_geometry

   use common_types, only: flux_surface_type

   implicit none

   public :: init_geometry, finish_init_geometry, finish_geometry
   public :: communicate_geo_multibox, x_displacement_fac 
   public :: bmag, dbdzed, btor, bmag_psi0, grho, grho_norm, grad_x
   public :: dcvdriftdrho, dcvdrift0drho, dgbdriftdrho, dgbdrift0drho
   public :: gds2, gds21, gds22, gds23, gds24, gds25, gds26, gradpar
   public :: b_dot_grad_zeta
   public :: cvdrift, cvdrift0, gbdrift, gbdrift0
   public :: dgds2dr, dgds21dr, dgds22dr
   public :: exb_nonlin_fac, exb_nonlin_fac_p, flux_fac
   public :: jacob, djacdrho, drhodpsi, drhodpsip, drhodpsip_psi0
   public :: dl_over_b, d_dl_over_b_drho
   public :: dBdrho, d2Bdrdth, dgradpardrho, dIdrho
   public :: geo_surf, Rmajor, dzetadz
   public :: theta_vmec, zeta, alpha 
   public :: dxdpsi, dydalpha, clebsch_factor
   public :: aref, bref, twist_and_shift_geo_fac
   public :: q_as_x, get_x_to_rho, gfac
   public :: dVolume, grad_x_grad_y_end  
   public :: b_dot_grad_z_averaged_eqarc, zed_eqarc ! Extended z-grid for final fields diagnostics
   public :: b_dot_grad_z ! Full flux surface
   public :: gradzeta_gradx_RRoverBB ! Momentum flux
   public :: gradzeta_grady_RRoverBB ! Momentum flux
   public :: b_dot_grad_zeta_RR ! Momentum flux
   public :: geo_option_switch, geo_option_vmec ! Used in kt_grids.f90

   private
   
   ! Debugging
   logical :: debug = .true.

   type(flux_surface_type) :: geo_surf

   real :: grad_x_grad_y_end, clebsch_factor
   real :: aref, bref, dxdpsi, dydalpha
   real :: dqdrho, dIdrho, grho_norm
   real :: drhodpsi, drhodpsip, drhodpsip_psi0, shat, qinp
   real :: exb_nonlin_fac, exb_nonlin_fac_p, flux_fac
   real :: b_dot_grad_z_averaged_eqarc, dzetadz
   real :: twist_and_shift_geo_fac, gfac

   ! Geometric quantities for the gyrokinetic equations (electrostatic, fluxtube)
   real, dimension(:), allocatable :: zed_eqarc, b_dot_grad_zeta, alpha
   real, dimension(:), allocatable :: gradpar, b_dot_grad_z_averaged
   real, dimension(:), allocatable :: dBdrho, d2Bdrdth, dgradpardrho, btor, Rmajor 
   real, dimension(:, :), allocatable :: bmag, bmag_psi0, dbdzed 
   real, dimension(:, :), allocatable :: cvdrift, cvdrift0, gbdrift, gbdrift0
   real, dimension(:, :), allocatable :: dcvdriftdrho, dcvdrift0drho, dgbdriftdrho, dgbdrift0drho
   real, dimension(:, :), allocatable :: gds2, gds21, gds22, gds23, gds24, gds25, gds26
   real, dimension(:, :), allocatable :: dgds2dr, dgds21dr, dgds22dr, x_displacement_fac
   real, dimension(:, :), allocatable :: jacob, djacdrho, grho, grad_x
   real, dimension(:, :), allocatable :: dl_over_b, d_dl_over_b_drho 
   real, dimension(:, :), allocatable :: theta_vmec, zeta
   real, dimension(:, :, :), allocatable :: dVolume
   
   ! Geometric quantities for full flux surface 
   real, dimension(:, :), allocatable :: b_dot_grad_z
   
   ! Geometric quantities for the momentum flux
   real, dimension(:, :), allocatable :: gradzeta_gradx_RRoverBB
   real, dimension(:, :), allocatable :: gradzeta_grady_RRoverBB
   real, dimension(:, :), allocatable :: b_dot_grad_zeta_RR

   integer :: sign_torflux
   integer :: geo_option_switch
   integer, parameter :: geo_option_local = 1
   integer, parameter :: geo_option_inputprof = 2
   integer, parameter :: geo_option_vmec = 3
   integer, parameter :: geo_option_multibox = 4
 
   logical :: overwrite_bmag, overwrite_b_dot_grad_zeta, overwrite_geometry
   logical :: overwrite_gds2, overwrite_gds21, overwrite_gds22
   logical :: overwrite_gds23, overwrite_gds24, overwrite_gbdrift
   logical :: overwrite_cvdrift, overwrite_gbdrift0, q_as_x
   character(100) :: geo_file
  
   logical :: geoinit = .false.
   logical :: set_bmag_const
   

contains

   !============================================================================
   !========================= INITIALIZE THE GEOMETRY ==========================
   !============================================================================
   subroutine init_geometry(nalpha, naky)

      ! Zgrid
      use zgrid, only: nzgrid, zed, delzed, shat_zero, grad_x_grad_y_zero
      use zgrid, only: boundary_option_switch, boundary_option_self_periodic
      use zgrid, only: boundary_option_linked, boundary_option_linked_stellarator 

      ! VMEC equilibria
      use vmec_geometry, only: read_vmec_parameters, get_vmec_geometry 

      ! Flags
      use physics_flags, only: include_geometric_variation

      ! Routines
      use file_utils, only: get_unused_unit
      use mp, only: proc0

      implicit none

      integer, intent(in) :: nalpha, naky
      real :: bmag_z0
      integer :: iy, ia, iz 

      !---------------------------------------------------------------------- 

      ! Only initialize once
      if (geoinit) return
      geoinit = .true.

      ! Default is no re-scaling of zed
      dzetadz = 1.0
      
      ! Track the code 
      if (debug) write (*, *) 'stella_geometry::init_geometry'

      ! Only calculate the geometry on proc0
      if (proc0) then
         
         ! Read the <geo_knobs> namelist in the input file
         call read_parameters 
 
         ! Use Miller parameters or VMEC to get the geometry needed for stella 
         if (geo_option_switch==geo_option_local)     call get_geometry_arrays_from_Miller(nalpha)
         if (geo_option_switch==geo_option_inputprof) call get_geometry_arrays_from_Miller(nalpha)
         if (geo_option_switch==geo_option_multibox)  call get_geometry_arrays_from_Miller(nalpha)
         if (geo_option_switch==geo_option_vmec)      call get_geometry_arrays_from_VMEC(nalpha, naky) 

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
         call get_b_dot_grad_z_averaged_eqarc(b_dot_grad_z_averaged, zed, delzed, b_dot_grad_z_averaged_eqarc)
         call get_zed_eqarc(b_dot_grad_zeta, delzed, zed, b_dot_grad_z_averaged_eqarc, zed_eqarc)
		   
		   ! A lot of modules use <gradpar> even though <b_dot_grad_z_averaged> is a better name  
         gradpar = b_dot_grad_z_averaged

      end if
 
      !=========================================================================
      !=================== CALCULATIONS ON ALL PROCESSORS  =====================
      !=========================================================================
      
      ! Track the code 
      if (debug) write (*, *) 'stella_geometry::init_geometry::calculate_on_all_processors'

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
      jacob = -clebsch_factor / (dydalpha * dxdpsi * b_dot_grad_z * bmag)
	
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

      ! AVB: temporary, set bmag = constant in z for Spitzer problem
      if (set_bmag_const) then
         bmag_z0 = bmag(1, 0)
         print *, ''
         print *, '! SETTING BMAG = CONSTANT IN Z'
         print *, ''
         do ia = 1, nalpha
            do iz = -nzgrid, nzgrid
               bmag(ia, iz) = bmag_z0
            end do
         end do
      end if
   

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

      use vmec_geometry, only: read_vmec_parameters, get_vmec_geometry 
      use vmec_geometry, only: radial_coordinate_option, radial_coordinate_sgnpsitpsit
      use vmec_geometry, only: radial_coordinate_minuspsit, radial_coordinate_r 
      use physics_flags, only: const_alpha_geo
      use zgrid, only: nzgrid 

      implicit none

      integer, intent(in) :: nalpha, naky

      real, dimension(:, :), allocatable :: psit_displacement_fac, grad_alpha_grad_alpha
      real, dimension(:, :), allocatable :: grad_alpha_grad_psit, grad_alpha_grad_psi 
      real, dimension(:, :), allocatable :: grad_psit_grad_psit, grad_psi_grad_psi
      real, dimension(:, :), allocatable :: gbdrift_alpha, cvdrift_alpha
      real, dimension(:, :), allocatable :: gbdrift0_psit, cvdrift0_psit
      real, dimension(:, :), allocatable :: gds23_alphapsit, gds24_alphapsit
      real, dimension(:, :), allocatable :: gds25_alphapsit, gds26_alphapsit
      real, dimension(:, :), allocatable :: grad_x_grad_x, grad_y_grad_y, grad_y_grad_x
      real, dimension(:, :), allocatable :: gradzeta_gradpsit_RRoverBB, gradzeta_gradalpha_RRoverBB, b_dot_grad_zeta_RR

      real :: rho, shat, iota, field_period_ratio
      real :: dpsidpsit, dpsidx, dpsitdx, dxdpsit 

      !---------------------------------------------------------------------- 

      ! Read the <vmec_parameters> namelist in the input file    
      call read_vmec_parameters()

      ! Allocate geometry arrays 
      call allocate_arrays(nalpha, nzgrid)
      call allocate_temporary_arrays(nalpha, nzgrid)

      ! Call the <vmec_geometry> module to calculate the geometric coefficients 
      ! needed by stella, based on the VMEC equilibrium file
      call get_vmec_geometry(nzgrid, nalpha, naky, geo_surf, grho, bmag, &
               b_dot_grad_z_averaged, b_dot_grad_z, & 
               grad_alpha_grad_alpha, grad_alpha_grad_psit, grad_psit_grad_psit, &
               gds23_alphapsit, gds24_alphapsit, gds25_alphapsit, gds26_alphapsit, & 
               gbdrift_alpha, gbdrift0_psit, cvdrift_alpha, cvdrift0_psit, &
               gradzeta_gradpsit_RRoverBB, gradzeta_gradalpha_RRoverBB, b_dot_grad_zeta_RR, &
               sign_torflux, theta_vmec, dzetadz, aref, bref, alpha, zeta, &
               field_period_ratio, psit_displacement_fac)

      ! Flux surface quantities that we need
      rho = geo_surf%rhotor 
      shat = geo_surf%shat
      iota = geo_surf%qinp

      ! Define the (psi,alpha) and (x,y) coordinates
      if (radial_coordinate_option==radial_coordinate_sgnpsitpsit) then
         dxdpsi = 1. / rho    
         dydalpha = rho
         dpsidpsit = sign_torflux
         drhodpsi = 1. / rho
         clebsch_factor = sign_torflux
      else if (radial_coordinate_option==radial_coordinate_minuspsit) then
         dxdpsi = -sign_torflux / rho    
         dydalpha = rho
         dpsidpsit = -1
         drhodpsi = -sign_torflux / rho
         clebsch_factor = -1
      else if (radial_coordinate_option==radial_coordinate_r) then
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
      gbdrift = gbdrift_alpha * dydalpha; gbdrift0 = gbdrift0_psit * dxdpsit 
      cvdrift = cvdrift_alpha * dydalpha; cvdrift0 = cvdrift0_psit * dxdpsit

      ! <grad_x> = |∇x|    <gds21> = hat{s} ∇x . ∇y  
      ! <gds2> = |∇y|^2    <gds22> = shat^2 * |∇x|^2 
      grad_x = sqrt(abs(grad_x_grad_x))
      gds2 = grad_y_grad_y
      gds21 = geo_surf%shat * grad_y_grad_x
      gds22 = (geo_surf%shat * grad_x)**2
      
      ! For the momentum flux we need (R^2/B^2) ∇ζ . ∇y and (R^2/B^2) ∇ζ . ∇x
      gradzeta_gradx_RRoverBB = gradzeta_gradpsit_RRoverBB * dxdpsit
      gradzeta_grady_RRoverBB = gradzeta_gradalpha_RRoverBB * dydalpha

      ! We want |ds/dx|*sqrt((dR/ds)^2+(dZ/ds)^2) so do dψ̃t/dx̃ * |ds̃/dψ̃t|*sqrt((dR/ds)^2+(dZ/ds)^2)  
      x_displacement_fac = dpsitdx*psit_displacement_fac

      ! Calculate the <twist_and_shift_geo_fac> needed for the boundary conditions
      call calculate_twist_and_shift_geo_fac()

      ! Deallocate arrays
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
 
         use zgrid, only: boundary_option_switch, boundary_option_linked_stellarator 
         use zgrid, only: shat_zero, grad_x_grad_y_zero
         use constants, only: pi
  
         implicit none  

         !---------------------------------------------------------------------- 

         ! Print the boundary condition information to the output file
         write (*, '(A)') "############################################################"
         write (*, '(A)') "                     BOUNDARY CONDITIONS"
         write (*, '(A)') "############################################################"
         select case (boundary_option_switch)

         ! Stellarator symmetric twist-and-shift BC (see Martin et al, 2018)
         case (boundary_option_linked_stellarator)

            ! dkx/dky * jtwist = -2  * (dψ/dx) (dy/dα) * (∇ψ . ∇α) / |∇ψ|^2
            write (*, *) ' '; write (*, *) 'Stellarator symmetric twist and shift BC selected' 
            twist_and_shift_geo_fac = -2. * dpsidx * dydalpha * (grad_alpha_grad_psi(1, nzgrid)) / (grad_psi_grad_psi(1, nzgrid))

         ! Standard twist-and-shift boundary conditions or unconnected BC (then twist_and_shift_geo_fac doesn't matter)
         case default

            ! dkx/dky * jtwist = -2*pi*P * (dψ/dx) (dy/dα) * (drho/dψ) * hat{s} (iota/rho)
            write (*, *) ' '; write (*, *) 'Standard twist and shift BC selected'
            twist_and_shift_geo_fac = -2. * pi * field_period_ratio * dpsidx * dydalpha * drhodpsi * shat * iota / rho

         end select

         ! If ∇x . ∇y is very small at the ends of the flux tube, use periodic boundary conditions
         grad_x_grad_y_end = grad_y_grad_x(1, nzgrid)

         ! Print this factor to the output file
         write (*, *) 'twist_and_shift_geo_fac: ', twist_and_shift_geo_fac; write (*, *) ' '
   
      end subroutine calculate_twist_and_shift_geo_fac
      

      !**********************************************************************
      !                      ALLOCATE TEMPORARY ARRAYS                      !
      !**********************************************************************
      subroutine allocate_temporary_arrays(nalpha, nzgrid)

         implicit none

         integer, intent(in) :: nalpha, nzgrid

         allocate (psit_displacement_fac(nalpha, -nzgrid:nzgrid)) 
         allocate (grad_alpha_grad_alpha(nalpha, -nzgrid:nzgrid))
         allocate (grad_alpha_grad_psi(nalpha, -nzgrid:nzgrid))
         allocate (grad_alpha_grad_psit(nalpha, -nzgrid:nzgrid))
         allocate (grad_psit_grad_psit(nalpha, -nzgrid:nzgrid))
         allocate (grad_psi_grad_psi(nalpha, -nzgrid:nzgrid))
         allocate (gbdrift_alpha(nalpha, -nzgrid:nzgrid))
         allocate (cvdrift_alpha(nalpha, -nzgrid:nzgrid))
         allocate (gbdrift0_psit(nalpha, -nzgrid:nzgrid))
         allocate (cvdrift0_psit(nalpha, -nzgrid:nzgrid))
         allocate (gds23_alphapsit(nalpha, -nzgrid:nzgrid))
         allocate (gds24_alphapsit(nalpha, -nzgrid:nzgrid))
         allocate (gds25_alphapsit(nalpha, -nzgrid:nzgrid))
         allocate (gds26_alphapsit(nalpha, -nzgrid:nzgrid)) 
         allocate (grad_x_grad_x(nalpha, -nzgrid:nzgrid))
         allocate (grad_y_grad_y(nalpha, -nzgrid:nzgrid))
         allocate (grad_y_grad_x(nalpha, -nzgrid:nzgrid))

      end subroutine allocate_temporary_arrays

      !**********************************************************************
      !                     DEALLOCATE TEMPORARY ARRAYS                     !
      !**********************************************************************
      subroutine deallocate_temporary_arrays

         implicit none

         deallocate (grad_y_grad_y, grad_x_grad_x, grad_y_grad_x)
         deallocate (psit_displacement_fac, grad_alpha_grad_alpha)
         deallocate (grad_alpha_grad_psi, grad_psit_grad_psit)
         deallocate (gbdrift_alpha, cvdrift_alpha)
         deallocate (gbdrift0_psit, cvdrift0_psit)
         deallocate (gds23_alphapsit, gds24_alphapsit)
         deallocate (gds25_alphapsit, gds26_alphapsit)

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

      use miller_geometry, only: read_local_parameters, get_local_geo
      use miller_geometry, only: communicate_parameters_multibox
      use inputprofiles_interface, only: read_inputprof_geo
      use zgrid, only: zed, nzed, nzgrid, zed_equal_arc
      use constants, only: pi

      implicit none   

      integer, intent(in) :: nalpha 
      real :: dpsipdrho, dpsipdrho_psi0

      !---------------------------------------------------------------------- 

      ! Read the <millergeo_parameters> namelist in the input file  
      ! <nzed> and <nzgrid> are inputs, and <geo_surf> is returned 
      ! which is a dictionary that contains all the geometry information 
      call read_local_parameters(nzed, nzgrid, geo_surf)

      ! Allocate geometry arrays for stella
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

      ! Call the <miller_geometry.f90> module to calculate the geometric coefficients 
      ! needed by stella, based on the local Miller parameters. For Miller geometries
      ! each field line <alpha> has the same geometry, hence we will only pass on the 
      ! ialpha=1 arrays to the get_local_geo() routine. Note that in Miller z = theta.
      call get_local_geo(nzed, nzgrid, zed, zed_equal_arc, &
                dpsipdrho, dpsipdrho_psi0, dIdrho, grho(1, :), bmag(1, :), bmag_psi0(1, :), &
                gds2(1, :), gds21(1, :), gds22(1, :), gds23(1, :), gds24(1, :), b_dot_grad_z(1, :), &
                gbdrift0(1, :), gbdrift(1, :), cvdrift0(1, :), cvdrift(1, :), &
                dBdrho, d2Bdrdth, dgradpardrho, btor, rmajor, &
                dcvdrift0drho(1, :), dcvdriftdrho(1, :), dgbdrift0drho(1, :), dgbdriftdrho(1, :), &
                dgds2dr(1, :), dgds21dr(1, :), dgds22dr(1, :), djacdrho(1, :))

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
      bmag = spread(bmag(1, :), 1, nalpha)
      bmag_psi0 = spread(bmag_psi0(1, :), 1, nalpha)
      gds2 = spread(gds2(1, :), 1, nalpha)
      gds21 = spread(gds21(1, :), 1, nalpha)
      gds22 = spread(gds22(1, :), 1, nalpha)
      gds23 = spread(gds23(1, :), 1, nalpha)
      gds24 = spread(gds24(1, :), 1, nalpha)
      gbdrift0 = spread(gbdrift0(1, :), 1, nalpha)
      gbdrift = spread(gbdrift(1, :), 1, nalpha)
      cvdrift0 = spread(cvdrift0(1, :), 1, nalpha)
      cvdrift = spread(cvdrift(1, :), 1, nalpha)
      dcvdrift0drho = spread(dcvdrift0drho(1, :), 1, nalpha)
      dcvdriftdrho = spread(dcvdriftdrho(1, :), 1, nalpha)
      dgbdrift0drho = spread(dgbdrift0drho(1, :), 1, nalpha)
      dgbdriftdrho = spread(dgbdriftdrho(1, :), 1, nalpha)
      dgds2dr = spread(dgds2dr(1, :), 1, nalpha)
      dgds21dr = spread(dgds21dr(1, :), 1, nalpha)
      dgds22dr = spread(dgds22dr(1, :), 1, nalpha)
      djacdrho = spread(djacdrho(1, :), 1, nalpha)
      b_dot_grad_z = spread(b_dot_grad_z(1, :), 1, nalpha)
      zeta = spread(zeta(1, :), 1, nalpha)

      ! For the momentum flux we need (R^2/B^2) ∇ζ . ∇y and (R^2/B^2) ∇ζ . ∇x
      ! For Miller or axi-symmetric devices we have: ∇ζ . ∇ψp = 0 and ∇ζ . ∇α = ∇ζ . ∇ζ = (1/R^2) 
      !		(R^2/B^2) * ∇ζ . ∇α * (dy/dα) = (R^2/B^2) (1/R^2) (dy/dα) = (1/B^2) (dy/dα) = geo_surf%rhoc / (geo_surf%qinp * bmag**2)
      gradzeta_grady_RRoverBB = geo_surf%rhoc / (geo_surf%qinp * bmag**2)
      gradzeta_gradx_RRoverBB = 0.0

      ! For the momentum flux we need R^2 * b . ∇ζ
      ! Note that in Miller z = theta, so b_dot_grad_z = b_dot_grad_theta
      !		R^2 * b . ∇ζ = R^2 * (1/B) (∇ζ x ∇ψ + I ∇ζ) . ∇ζ =  R^2/B * I * ∇ζ . ∇ζ
      !                  = R^2/B * I * (1/R^2) = I/B = (R/B) (I/R) = (R/B) * Btor
      b_dot_grad_zeta_RR = geo_surf%rmaj * spread(btor, 1, nalpha) / bmag 

   end subroutine get_geometry_arrays_from_Miller

   !=========================================================================
   !================== OVERWRITE GEOMETRIC COEFFICIENTS  ====================
   !=========================================================================
   subroutine overwrite_selected_geometric_coefficients(nalpha)

      use file_utils, only: get_unused_unit
      use zgrid, only: nzgrid

      implicit none

      integer, intent(in) :: nalpha
      integer :: geofile_unit
      character(100) :: dum_char
      real :: dum_real

      integer :: ia, iz
      real :: bmag_file, b_dot_grad_zeta_file
      real :: gds2_file, gds21_file, gds22_file, gds23_file, gds24_file
      real :: gbdrift_file, cvdrift_file, gbdrift0_file

      call get_unused_unit(geofile_unit)
      open (geofile_unit, file=trim(geo_file), status='old', action='read')

      read (geofile_unit, fmt=*) dum_char
      read (geofile_unit, fmt=*) dum_char
      read (geofile_unit, fmt=*) dum_char

      ! overwrite bmag, b_dot_grad_zeta, gds2, gds21, gds22, gds23, gds24, gbdrift, cvdrift, gbdrift0, and cvdrift0
      ! with values from file
      do ia = 1, nalpha
         do iz = -nzgrid, nzgrid
            read (geofile_unit, fmt='(13e12.4)') dum_real, dum_real, dum_real, bmag_file, b_dot_grad_zeta_file, &
               gds2_file, gds21_file, gds22_file, gds23_file, &
               gds24_file, gbdrift_file, cvdrift_file, gbdrift0_file
            if (overwrite_bmag) bmag(ia, iz) = bmag_file
            if (overwrite_b_dot_grad_zeta) b_dot_grad_zeta(iz) = b_dot_grad_zeta_file
            ! assuming we are only reading in for a single alpha. Usually, b_dot_grad_zeta is the average of all b_dot_grad_z values.
            if (overwrite_b_dot_grad_zeta) b_dot_grad_z(1, iz) = b_dot_grad_zeta_file 
            if (overwrite_gds2) gds2(ia, iz) = gds2_file
            if (overwrite_gds21) gds21(ia, iz) = gds21_file
            if (overwrite_gds22) gds22(ia, iz) = gds22_file
            if (overwrite_gds23) gds23(ia, iz) = gds23_file
            if (overwrite_gds24) gds24(ia, iz) = gds24_file
            if (overwrite_gbdrift) gbdrift(ia, iz) = gbdrift_file
            if (overwrite_cvdrift) cvdrift(ia, iz) = cvdrift_file
            if (overwrite_gbdrift0) gbdrift0(ia, iz) = gbdrift0_file
         end do
      end do
      cvdrift0 = gbdrift0

      close (geofile_unit)

   end subroutine overwrite_selected_geometric_coefficients

   !=========================================================================
   !============= MAKE GEOMETRIC COEFFICIENTS CONSTANT IN ALPHA =============
   !=========================================================================
   subroutine set_ffs_geo_coefs_constant(nalpha)

      implicit none

      integer, intent(in) :: nalpha

      call set_coef_constant(gbdrift0, nalpha)
      call set_coef_constant(cvdrift0, nalpha)
      call set_coef_constant(gbdrift, nalpha)
      call set_coef_constant(cvdrift, nalpha)
      call set_coef_constant(grad_x, nalpha)
      call set_coef_constant(grho, nalpha)
      call set_coef_constant(bmag, nalpha)
      call set_coef_constant(bmag_psi0, nalpha)
      call set_coef_constant(gds2, nalpha)
      call set_coef_constant(gds21, nalpha)
      call set_coef_constant(gds22, nalpha)
      call set_coef_constant(gds23, nalpha)
      call set_coef_constant(gds24, nalpha)
      call set_coef_constant(gds25, nalpha)
      call set_coef_constant(gds26, nalpha)
      call set_coef_constant(theta_vmec, nalpha)
      call set_coef_constant(x_displacement_fac, nalpha)
      call set_coef_constant(zeta, nalpha)
      call set_coef_constant(b_dot_grad_z, nalpha)
      call set_coef_constant(gradzeta_gradx_RRoverBB, nalpha)
      call set_coef_constant(gradzeta_grady_RRoverBB, nalpha)
      call set_coef_constant(b_dot_grad_zeta_RR, nalpha)

   end subroutine set_ffs_geo_coefs_constant

   subroutine set_coef_constant(coef, nalpha)

      use zgrid, only: nzgrid

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

      if (.not. allocated(bmag)) allocate (bmag(nalpha, -nzgrid:nzgrid)); bmag = 0.0 
      if (.not. allocated(bmag_psi0)) allocate (bmag_psi0(nalpha, -nzgrid:nzgrid)); bmag_psi0 = 0.0
      if (.not. allocated(gds2)) allocate (gds2(nalpha, -nzgrid:nzgrid)); gds2 = 0.0
      if (.not. allocated(gds21)) allocate (gds21(nalpha, -nzgrid:nzgrid)); gds21 = 0.0
      if (.not. allocated(gds22)) allocate (gds22(nalpha, -nzgrid:nzgrid)); gds22 = 0.0
      if (.not. allocated(gds23)) allocate (gds23(nalpha, -nzgrid:nzgrid)); gds23 = 0.0
      if (.not. allocated(gds24)) allocate (gds24(nalpha, -nzgrid:nzgrid)); gds24 = 0.0
      if (.not. allocated(gds25)) allocate (gds25(nalpha, -nzgrid:nzgrid)); gds25 = 0.0
      if (.not. allocated(gds26)) allocate (gds26(nalpha, -nzgrid:nzgrid)); gds26 = 0.0
      if (.not. allocated(dgds2dr)) allocate (dgds2dr(nalpha, -nzgrid:nzgrid)); dgds2dr = 0.0
      if (.not. allocated(dgds21dr)) allocate (dgds21dr(nalpha, -nzgrid:nzgrid)); dgds21dr = 0.0
      if (.not. allocated(dgds22dr)) allocate (dgds22dr(nalpha, -nzgrid:nzgrid)); dgds22dr = 0.0
      if (.not. allocated(gbdrift)) allocate (gbdrift(nalpha, -nzgrid:nzgrid)); gbdrift = 0.0
      if (.not. allocated(gbdrift0)) allocate (gbdrift0(nalpha, -nzgrid:nzgrid)); gbdrift0 = 0.0
      if (.not. allocated(cvdrift)) allocate (cvdrift(nalpha, -nzgrid:nzgrid)); cvdrift = 0.0
      if (.not. allocated(cvdrift0)) allocate (cvdrift0(nalpha, -nzgrid:nzgrid)); cvdrift0 = 0.0
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
      if (.not. allocated(b_dot_grad_z)) allocate (b_dot_grad_z(nalpha, -nzgrid:nzgrid)); b_dot_grad_z = 0.0
      if (.not. allocated(gradpar)) allocate (gradpar(-nzgrid:nzgrid)); gradpar = 0.0
      if (.not. allocated(zed_eqarc)) allocate (zed_eqarc(-nzgrid:nzgrid)); zed_eqarc = 0.0
      if (.not. allocated(btor)) allocate (btor(-nzgrid:nzgrid)); btor = 0.0
      if (.not. allocated(rmajor)) allocate (rmajor(-nzgrid:nzgrid)); rmajor = 0.0
      if (.not. allocated(dBdrho)) allocate (dBdrho(-nzgrid:nzgrid)); dBdrho = 0.0
      if (.not. allocated(d2Bdrdth)) allocate (d2Bdrdth(-nzgrid:nzgrid)); d2Bdrdth = 0.0
      if (.not. allocated(dgradpardrho)) allocate (dgradpardrho(-nzgrid:nzgrid)); dgradpardrho = 0.0
      if (.not. allocated(alpha)) allocate (alpha(nalpha)); alpha = 0.0
      if (.not. allocated(zeta)) allocate (zeta(nalpha, -nzgrid:nzgrid)); zeta = 0.0
      if (.not. allocated(x_displacement_fac)) allocate (x_displacement_fac(nalpha, -nzgrid:nzgrid)); x_displacement_fac = 0.0

      ! Needed for the momentum flux diagnostic for non-axisymmetric devices
      if (.not. allocated(gradzeta_gradx_RRoverBB)) allocate (gradzeta_gradx_RRoverBB(nalpha, -nzgrid:nzgrid)); gradzeta_gradx_RRoverBB = 0.0
      if (.not. allocated(gradzeta_grady_RRoverBB)) allocate (gradzeta_grady_RRoverBB(nalpha, -nzgrid:nzgrid)); gradzeta_grady_RRoverBB = 0.0
      if (.not. allocated(b_dot_grad_zeta_RR)) allocate (b_dot_grad_zeta_RR(nalpha, -nzgrid:nzgrid)); b_dot_grad_zeta_RR = 0.0

   end subroutine allocate_arrays

   !============================================================================
   !============= READ THE <GEO_KNOBS> NAMELIST FROM THE INPUT FILE ============
   !============================================================================
   subroutine read_parameters

      ! Flags
      use physics_flags, only: radial_variation

      ! Multibox runs
      use file_utils, only: runtype_option_Switch, runtype_multibox

      ! Routines
      use text_options, only: text_option, get_option_value
      use file_utils, only: error_unit, input_unit_exist
      use mp, only: job

      implicit none

      character(20) :: geo_option
      integer :: in_file, ierr
      logical :: exist

      ! Text options for <geo_option> in the <geo_knobs> namrlist
      type(text_option), dimension(5), parameter :: geoopts = (/ &
             text_option('default', geo_option_local), &
             text_option('miller', geo_option_local), &
             text_option('local', geo_option_local), &
             text_option('input.profiles', geo_option_inputprof), &
             text_option('vmec', geo_option_vmec)/)

      !---------------------------------------------------------------------- 

      ! Define the variables in the namelist
      namelist /geo_knobs/ geo_option, geo_file, overwrite_bmag, overwrite_b_dot_grad_zeta, &
         overwrite_gds2, overwrite_gds21, overwrite_gds22, overwrite_gds23, overwrite_gds24, &
         overwrite_gbdrift, overwrite_cvdrift, overwrite_gbdrift0, q_as_x, set_bmag_const

      ! Assign default variables
      geo_option = 'local'
      overwrite_bmag = .false.
      overwrite_b_dot_grad_zeta = .false.
      overwrite_gds2 = .false.
      overwrite_gds21 = .false.
      overwrite_gds22 = .false.
      overwrite_gds23 = .false.
      overwrite_gds24 = .false.
      overwrite_gbdrift = .false.
      overwrite_cvdrift = .false.
      overwrite_gbdrift0 = .false.
      set_bmag_const = .false.
      geo_file = 'input.geometry'

      ! The following is True by default in radial variation runs
      q_as_x = radial_variation 

      ! Read the <geo_knobs> namelist in the input file
      in_file = input_unit_exist("geo_knobs", exist)
      if (exist) read (unit=in_file, nml=geo_knobs)

      ! Read <geo_option> in the input file
      ierr = error_unit()
      call get_option_value (geo_option, geoopts, geo_option_switch, ierr, "geo_option in geo_knobs")

      ! Multibox run
      if (radial_variation .and. runtype_option_switch == runtype_multibox .and. job /= 1) then
         geo_option_switch = geo_option_multibox
      end if

      ! If any geometric array needs to be overwritten, set <overwrite_geometry> = True
      overwrite_geometry = overwrite_bmag .or. overwrite_b_dot_grad_zeta &
            .or. overwrite_gds2 .or. overwrite_gds21 .or. overwrite_gds22 &
            .or. overwrite_gds23 .or. overwrite_gds24 &
            .or. overwrite_cvdrift .or. overwrite_gbdrift .or. overwrite_gbdrift0

   end subroutine read_parameters

   !============================================================================
   !============================= BROADCAST ARRAYS =============================
   !============================================================================
   subroutine broadcast_arrays

      use mp, only: broadcast

      implicit none

      ! Flags 
      call broadcast(q_as_x)
      call broadcast(set_bmag_const)

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
      call broadcast(gradpar)
      call broadcast(b_dot_grad_zeta)
      call broadcast(b_dot_grad_z)
      call broadcast(b_dot_grad_z_averaged) 
      call broadcast(gds2)
      call broadcast(gds21)
      call broadcast(gds22)
      call broadcast(gds23)
      call broadcast(gds24)
      call broadcast(gds25)
      call broadcast(gds26)
      call broadcast(dgds2dr)
      call broadcast(dgds21dr)
      call broadcast(dgds22dr)
      call broadcast(gbdrift0)
      call broadcast(gbdrift)
      call broadcast(cvdrift0)
      call broadcast(cvdrift)
      call broadcast(dgbdrift0drho)
      call broadcast(dgbdriftdrho)
      call broadcast(dcvdrift0drho)
      call broadcast(dcvdriftdrho)
      call broadcast(dBdrho)
      call broadcast(d2Bdrdth)
      call broadcast(dgradpardrho)
      call broadcast(djacdrho)
      
      ! Arrays for the momentum flux
      call broadcast(gradzeta_gradx_RRoverBB)
      call broadcast(gradzeta_grady_RRoverBB)
      call broadcast(b_dot_grad_zeta_RR)
   
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

      use miller_geometry, only: communicate_parameters_multibox
      use mp, only: proc0

      implicit none

      real, intent(in) :: l_edge, r_edge

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
   subroutine get_dzed(nz, dz, f, df)

      implicit none

      integer, intent(in) :: nz
      real, dimension(-nz:), intent(in) :: dz, f
      real, dimension(-nz:), intent(out) :: df

      df(-nz + 1:nz - 1) = (f(-nz + 2:) - f(:nz - 2)) / (dz(:nz - 2) + dz(-nz + 1:nz - 1))

      ! TODO-GA hack to avoid non-periodicity in full-flux-surface case
		! if (full_flux_surface .and. .not. const_alpha_geo) then
      !   df(-nz) = (f(-nz + 1) - f(-nz)) / dz(-nz)
      !  df(nz) = (f(nz) - f(nz - 1)) / dz(nz - 1)
      !else
      ! assume periodicity in the B-field
      df(-nz) = (f(-nz + 1) - f(nz - 1)) / (dz(-nz) + dz(nz - 1))
      df(nz) = df(-nz)
   	!end if

   end subroutine get_dzed

   !============================================================================ 
   !========================= CALCULATE b_dot_grad_zeta EQARC ==========================
   !============================================================================
   subroutine get_b_dot_grad_z_averaged_eqarc(gp, z, dz, gp_eqarc)

      use constants, only: pi
      use zgrid, only: nzgrid

      implicit none

      real, dimension(-nzgrid:), intent(in) :: gp, z, dz
      real, intent(out) :: gp_eqarc

      ! first get int dz b . grad z
      call integrate_zed(dz, 1./gp, gp_eqarc)
      ! then take (zmax-zmin)/int (dz b . gradz)
      ! to get b . grad z'
      gp_eqarc = (z(nzgrid) - z(-nzgrid)) / gp_eqarc

   end subroutine get_b_dot_grad_z_averaged_eqarc

   !============================================================================ 
   !=========================== CALCULATE ZED EQARC ============================
   !============================================================================
   subroutine get_zed_eqarc(gp, dz, z, gp_eqarc, z_eqarc)

      use zgrid, only: nzgrid

      implicit none

      real, dimension(-nzgrid:), intent(in) :: gp, dz, z
      real, intent(in) :: gp_eqarc
      real, dimension(-nzgrid:), intent(out) :: z_eqarc

      integer :: iz

      z_eqarc(-nzgrid) = z(-nzgrid)
      do iz = -nzgrid + 1, nzgrid
         call integrate_zed(dz(:iz), 1./gp(:iz), z_eqarc(iz))
      end do
      z_eqarc(-nzgrid + 1:) = z(-nzgrid) + z_eqarc(-nzgrid + 1:) * gp_eqarc

   end subroutine get_zed_eqarc

   !============================================================================ 
   !============================== INTEGRATE ZED ===============================
   !============================================================================
   ! trapezoidal rule to integrate in zed
   subroutine integrate_zed(dz, f, intf)

      use zgrid, only: nzgrid

      implicit none

      real, dimension(-nzgrid:), intent(in) :: dz
      real, dimension(-nzgrid:), intent(in) :: f
      real, intent(out) :: intf

      integer :: iz, iz_max

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

      use physics_parameters, only: rhostar

      implicit none

      integer, intent(in) :: llim

      real, dimension(:), intent(in) :: x_in
      real, dimension(:), intent(out) :: rho_out

      integer :: ix, ulim

      real :: a, b, c

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
      use zgrid, only: nzgrid, zed

      implicit none

      integer, intent(in) :: nalpha
      integer :: geometry_unit
      integer :: ia, iz

      call open_output_file(geometry_unit, '.geometry')

      write (geometry_unit, '(a1,11a14)') '#', 'rhoc', 'qinp', 'shat', 'rhotor', 'aref', 'bref', 'dxdpsi', 'dydalpha', &
         'exb_nonlin', 'flux_fac'
      write (geometry_unit, '(a1,11e14.4)') '#', geo_surf%rhoc, geo_surf%qinp, geo_surf%shat, geo_surf%rhotor, aref, bref, &
         dxdpsi, dydalpha, exb_nonlin_fac, flux_fac
      write (geometry_unit, *)

      write (geometry_unit, '(15a12)') '# alpha', 'zed', 'zeta', 'bmag', 'bdot_grad_z', 'gds2', 'gds21', 'gds22', &
         'gds23', 'gds24', 'gbdrift', 'cvdrift', 'gbdrift0', 'bmag_psi0', 'btor'
      do ia = 1, nalpha
         do iz = -nzgrid, nzgrid
            !          write (geometry_unit,'(15e12.4)') alpha(ia), zed(iz), zeta(ia,iz), bmag(ia,iz), b_dot_grad_zeta(iz), &
            write (geometry_unit, '(15e12.4)') alpha(ia), zed(iz), zeta(ia, iz), bmag(ia, iz), b_dot_grad_z(ia, iz), &
               gds2(ia, iz), gds21(ia, iz), gds22(ia, iz), gds23(ia, iz), &
               gds24(ia, iz), gbdrift(ia, iz), cvdrift(ia, iz), gbdrift0(ia, iz), &
               bmag_psi0(ia, iz), btor(iz)
         end do
         write (geometry_unit, *)
      end do

      call close_output_file(geometry_unit)

   end subroutine write_geometric_coefficients

   !============================================================================ 
   !====================== FINISH INITIALIZING THE GEOMETRY ====================
   !============================================================================
   subroutine finish_init_geometry

      use mp, only: proc0
      use miller_geometry, only: finish_local_geo

      implicit none

      if (proc0) then
         select case (geo_option_switch)
         case (geo_option_local)
            call finish_local_geo
         case (geo_option_multibox)
            call finish_local_geo
         case (geo_option_inputprof)
            call finish_local_geo
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
      if (allocated(gradpar)) deallocate (gradpar) 
      if (allocated(b_dot_grad_z)) deallocate (b_dot_grad_z)
      if (allocated(dl_over_b)) deallocate (dl_over_b)
      if (allocated(d_dl_over_b_drho)) deallocate (d_dl_over_b_drho)
      if (allocated(gds2)) deallocate (gds2)
      if (allocated(gds21)) deallocate (gds21)
      if (allocated(gds22)) deallocate (gds22)
      if (allocated(gds23)) deallocate (gds23)
      if (allocated(gds24)) deallocate (gds24)
      if (allocated(gds25)) deallocate (gds25)
      if (allocated(gds26)) deallocate (gds26)
      if (allocated(dgds2dr)) deallocate (dgds2dr)
      if (allocated(dgds21dr)) deallocate (dgds21dr)
      if (allocated(dgds22dr)) deallocate (dgds22dr)
      if (allocated(gbdrift)) deallocate (gbdrift)
      if (allocated(gbdrift0)) deallocate (gbdrift0)
      if (allocated(cvdrift)) deallocate (cvdrift)
      if (allocated(cvdrift0)) deallocate (cvdrift0)
      if (allocated(dgbdriftdrho)) deallocate (dgbdriftdrho)
      if (allocated(dcvdriftdrho)) deallocate (dcvdriftdrho)
      if (allocated(dgbdrift0drho)) deallocate (dgbdrift0drho)
      if (allocated(dcvdrift0drho)) deallocate (dcvdrift0drho)
      if (allocated(dBdrho)) deallocate (dBdrho)
      if (allocated(d2Bdrdth)) deallocate (d2Bdrdth)
      if (allocated(dgradpardrho)) deallocate (dgradpardrho)
      if (allocated(theta_vmec)) deallocate (theta_vmec)

      if (allocated(alpha)) deallocate (alpha)
      if (allocated(zeta)) deallocate (zeta)

      if (allocated(x_displacement_fac)) deallocate (x_displacement_fac)
      
      ! Arrays for the momentum flux 
      if (allocated(gradzeta_gradx_RRoverBB)) deallocate (gradzeta_gradx_RRoverBB)
      if (allocated(gradzeta_grady_RRoverBB)) deallocate (gradzeta_grady_RRoverBB)
      if (allocated(b_dot_grad_zeta_RR)) deallocate (b_dot_grad_zeta_RR)

      geoinit = .false.

   end subroutine finish_geometry

end module stella_geometry
