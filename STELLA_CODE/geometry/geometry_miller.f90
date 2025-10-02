!###############################################################################
!      Construct a magnetic equilibrium based on a set of Miller parameters     
!###############################################################################
! 
! This module constructs a magnetic equilibrium based on a set of Miller parameters.
! 
!--------------------------------- Mathematics ---------------------------------
! 
! Miller solves the Grad Shafranov equation locally: 
! 
!                       R^2 ∇ · (∇ψ/R^2) = −4π R^2 p′− II′
! 
! From this we can get all of the local geometric variables.
! 
! The Miller equilibrium is a local a formalism capable of describing the local 
! magnetic geometry of a flux surface within axisymmetric systems. This approach
! ensures that the Grad-Shafranov equation is locally satisfied in ψ.
! 
! The definition of the Miller equilibrium is given in cyclindrical coordinates:
! 
!               R(r,θ) = R0(r) + r cos[ θ + sin θ arcsin δ(r) ]
!               Z(r,θ) = κ(r) r sin θ
! 
! Note that stella defines <tri> = arcsin δ(r).
! 
! 
! The GS equation can be expressed using this R and Z, giving a long expression, 
! which will ~eventually~ be in the "stella_manual". For this we require terms like 
! ∂Z/∂θ, ∂R/∂θ, ∂Z/∂r, ∂R/∂r etc. (plus second derivatives). These will depend on
! quantities like R0, r, δ and κ. e.g.
! 
!              ∂R/∂r = R0' + cos[θ+ sin(θδ)]− rsin{θsin[θ+ sin(θδ)]}δ′
!              ∂Z/∂r = (κ′r+ κ) sinθ
!              ∂R/∂θ = −sin[θ+ sin(θδ)][1 + cos(θδ)]
!              ∂Z/∂θ = κ * r * cosθ
! 
! Once we have these derivatives we can calclulate the Jacobian for the system and 
! and  other geometric quantities, such as,
! 
!                     |∇r|^2 = R^2/Jr^2 [(∂Z/∂θ)^2 + (∂R/∂θ)^2]
! 
! etc. These all depend on the form of R and Z, and require the Jacobians.
! 
!------------------------------- Input Variables -------------------------------
! 
!&geometry_miller
!   rhoc = 0.5
!   rmaj = 3.0
!   shift = 0.0
!   qinp = 1.4
!   shat = 0.8
!   kappa = 0.0
!   kapprim = 0.0
!   tri = 0.0
!   triprim = 0.0
!   rgeo = 3.0
!   betaprim = 0.0
!   betadbprim = 0.0
!   d2qdr2 = 0.0
!   d2psidr2 = 0.0
!   nzed_local = 128.0
!   read_profile_variation = .false.
!   write_profile_variation = .false.
!/
! 
! Definitions of input Miller quantities:
!     - rhoc         -->   Minor radius, r
!     - rmaj         -->   Major radius, R0
!     - shift        -->   Shafranov shift, R0'
!     - kappa        -->   Elongation, κ
!     - kappprim     -->   Radial derivative of elongation, κ'
!     - tri          -->   Triangularity, δ
!     - triprim      -->   Radial derivative of triangularity, δ'
!     - qinp         -->   Safety factor, q
!     - shat         -->   Magnetic shear, s
!     - rgeo         -->   Sets the refernce magnetic field, R_geo
!     - betaprim     -->   Radial derivative of plasma beta, β'
!                          β' = 4*π*p/B0^2 * (- ∂lnp/∂r )
!                          where p is the total pressure and B0 is the 
!                          magnetic field strength
! 
! Other reference variables that are needed (but inferred):
!     - bref         -->   Magnetic field strength reference, B_ref
!     - aref         -->   Minor radius reference length, a_ref
! 
! 
!------------------------------- Code Specifics --------------------------------
! 
! This module first reads in the miller parameters in the input file,
! and stores them as local%(name of variable) or surf%(name of variable). 
! Once the input variables have been read, this module:
! 
!     - Defines dq/∂r = s * q / r
!     - Gets the forms of R and Z (using functions called Rpos and Zpos, and are
!       stored as Rr and Zr)
!     - Then all the necessary derivatives, and second derivatives are computed, 
!       along with the Jacobians
!     - With these derivatives, geometric quantities can be calculated, e.g. 
!                 <bmag> = B/B0 = sqrt(I^2 + |∇ψ|^2)/R
!                 <b_dot_gradtheta> = b . ∇θ = dψ/∂r / (B/B0 * Jacobian_r)
!     - These are the quantities used in geometry.f90 to compute the geometric
!       coefficients that appear in the code. 
!
!---------------------------- Geometric quantities -----------------------------
! 
!    <b_dot_gradtheta>(ia,iz) = b · ∇θ
!    <b_dot_gradB> = b . ∇B
!    <d_bdotgradtheta_drho> = d (b . ∇θ) / drho
!    <d_bdotgradB_drho> = d (b . ∇B) / drho
! 
!--------------------------- Backwards Compatibility ---------------------------
! 
! An overview of the name changes implemented in September 2025 are given here,
!    - gradpar           -->   b_dot_gradtheta
!    - gradparB          -->   b_dot_gradB
!    - dgradpardrho      -->   d_bdotgradtheta_drho
!    - dgradparBdrho     -->   d_bdotgradB_drho
! 
!###############################################################################
module geometry_miller

   ! Load debug flags
   use debug_flags, only: debug => geometry_debug
   use common_types, only: flux_surface_type

   implicit none

   ! Make routines available to other modules
   public :: read_miller_parameters
   public :: get_miller_geometry
   public :: communicate_parameters_multibox
   public :: finish_miller_geometry
   public :: local

   private

   !----------------------------- Input Variables ------------------------------
   
   ! Variables read from the "geometry_miller" namelist in the input file
   real :: rhoc, rgeo, rmaj, shift
   real :: kappa, kapprim
   real :: tri, triprim
   real :: betaprim, betadbprim
   real :: qinp, shat, d2qdr2
   real :: dpsip_drho, d2psidr2, dpsip_drho_psi0
   real :: psitor_lcfs, rhotor
   real :: drhotordrho, dI_drho, rhoc0
   real :: dI = 0.0
   logical :: write_profile_variation
   logical :: read_profile_variation
   logical :: load_psi0_variables
   
   !----------------------------- Local Variables ------------------------------

   ! Local variables for the z-gid
   integer :: nzed_local
   integer :: nz, nzgrid_2pi
   integer :: nperiod

   ! Geometric quantities
   real, dimension(:), allocatable :: gradrho, bmag, gradrho_psi0, bmag_psi0
   real, dimension(:), allocatable :: b_dot_gradtheta, b_dot_gradtheta_arc, arc
   real, dimension(:), allocatable :: grady_dot_grady, gradx_dot_grady, gradx_dot_gradx
   real, dimension(:), allocatable :: B_times_gradB_dot_gradx, B_times_gradB_dot_grady
   real, dimension(:), allocatable :: B_times_kappa_dot_gradx, B_times_kappa_dot_grady
   real, dimension(:), allocatable :: gds23, gds24
   real, dimension(:), allocatable :: d2R_dtheta2, d2Z_dtheta2, d2R_dr_dtheta, d2Z_dr_dtheta
   real, dimension(:), allocatable :: gradpsi, dB_drho, d2B_dR_dtheta
   real, dimension(:), allocatable :: d_bdotgradtheta_drho, d_bdotgradB_drho, dB_dtheta, b_dot_gradB
   real, dimension(:), allocatable :: d_Btimeskappadotgradx_drho, d_BtimesgradBdotgradx_drho, theta
   real, dimension(:), allocatable :: vartheta, dvartheta_dr, gradrho_gradtheta
   real, dimension(:), allocatable :: gradalpha_times_B_dot_gradtheta, d2vartheta_dr2
   real, dimension(:), allocatable :: gradtheta2, gradalpha_gradtheta, gradrho_gradalpha, gradalpha2
   real, dimension(:), allocatable :: d2B_dr2, d2R_dr2, d2Z_dr2, drz, drz_dtheta
   real, dimension(:), allocatable :: d3R_dr2_dtheta, d3Z_dr2_dtheta
   real, dimension(:), allocatable :: d2gradpsi_dr2, d_gradalphatimesBdotgradtheta_drho
   real, dimension(:), allocatable :: d_Btimeskappadotgrady_drho, d_BtimesgradBdotgrady_drho
   real, dimension(:), allocatable :: d_gradydotgrady_drho, d_gradxdotgrady_drho, d_gradxdotgradx_drho
   real, dimension(:), allocatable :: d_gradr2_drho, d_gradpsi2_drho
   real, dimension(:), allocatable :: d_nablarnablatheta_drho, d_nablatheta2_drho
   real, dimension(:), allocatable :: d_nablaalphanablarho_drho, d_nablaalphanablatheta_drho, d_nablaalpha2_drho
   real, dimension(:), allocatable :: jacrho, delta_theta, djac_drho, djacr_drho
   real, dimension(:), allocatable :: d2jac_dr2, dR_drho, dZ_drho, dR_dtheta, dZ_dtheta
   real, dimension(:), allocatable :: d2R, d2Z
   real, dimension(:, :), allocatable :: Rr, Zr
   real :: rgeo_temp, dqdr, d2I_dr2
   
   ! Store the characterics of the flux surface as local%(name of variable)
   type(flux_surface_type) :: local

   ! Only initialise once
   logical :: initialised_miller = .false.

contains

!===============================================================================
!================================= INITIALISE ==================================
!===============================================================================

   !****************************************************************************
   !                            READ LOCAL PARAMETERS                           
   !****************************************************************************
   subroutine read_miller_parameters(nzed, nzgrid, local_out)

      use file_utils, only: input_unit_exist
      use common_types, only: flux_surface_type
      use namelist_geometry, only: read_namelist_geometry_miller

      implicit none

      ! Arguments
      type(flux_surface_type), intent(out) :: local_out
      integer, intent(in) :: nzed, nzgrid

      ! Local variables
      real :: dum
      integer :: j

      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_miller) return
      initialised_miller = .true.
         
      ! Read the "geometry_miller" namelist in the input file
      call read_namelist_geometry_miller(rhoc, rmaj, shift, qinp, shat, &
         kappa, kapprim, tri, triprim, rgeo, betaprim, &
         betadbprim, d2qdr2, d2psidr2, nzed_local, &
         read_profile_variation, write_profile_variation)

      ! The following variables are only needed for sfincs, when we are not 
      ! reading the geometry info from a file.
      rhotor = rhoc
      psitor_lcfs = 1.0
      drhotordrho = 1.0

      ! Save the input variables 
      local%rhoc = rhoc
      local%rmaj = rmaj
      local%rgeo = rgeo
      local%shift = shift
      local%kappa = kappa
      local%kapprim = kapprim
      local%qinp = qinp
      local%shat = shat
      local%tri = tri
      local%triprim = triprim
      local%betaprim = betaprim
      local%betadbprim = betadbprim
      local%d2qdr2 = d2qdr2
      local%d2psidr2 = d2psidr2
      local%zed0_fac = 1.0
      local%dr = 1.e-3 * (rhoc / rmaj)
      local%rhotor = rhotor
      local%psitor_lcfs = psitor_lcfs
      local%drhotordrho = drhotordrho
      local%dpsitordrho = 0.0
      local%d2psitordrho2 = 0.0

      ! Variables for multibox simulations with radial variation
      local%rhoc_psi0 = rhoc
      local%qinp_psi0 = qinp
      local%shat_psi0 = shat
      
      ! Note that <nzed> is an input parameter, denoting the number of z-points in a single period or 2*pi segment.
      ! It is used to calculate the total number of positive z-points as <nzgrid> = <nzed> / 2 + (<nperiod> - 1) * <nzed>.
      ! and the total number of z-points is given by <nztot> = 2 * <nzgrid> + 1
      ! Instead of importing <np> = <nperiod> we calculate it based on <nzed> and <nzgrid>.
      nzgrid_2pi = nzed / 2
      nperiod = (nzgrid - nzgrid_2pi) / nzed + 1

      ! Now we switch to using (possible higher resolution) local grid.
      ! The number of positive z-points in a single period or 2*pi segment is given by <nzgrid_2pi>,
      ! and the total number of positive z-points is given by <nz>, equivalent to <nzgrid> on the local grid.
      nzgrid_2pi = nzed_local / 2
      nz = nzgrid_2pi + nzed_local * (nperiod - 1)
     
      ! Set variables for radial variation simulations
      call read_miller_parameters_radial_variation

      ! Give the flux surface variables to the geometry.f90 module as <local_out> = <geo_surf>
      local_out = local

   contains
   
      ! ------------------------------------------------------------------------
      !                              Radial Variation
      ! ------------------------------------------------------------------------
      subroutine read_miller_parameters_radial_variation

         implicit none
         
         !----------------------------------------------------------------------
            
         ! The following variables are only needed for radial variation simulations.
         ! If <read_profile_variation>=True, it reads <rhoc0>, and sets <load_psi0_variables>=False. 
         rhoc0 = 0.5
         load_psi0_variables = .true.
      
         ! Initialise arrays to zero. These will be overwritten if reading in from file.
         ! They are only relevant for radial_variation tests.
         allocate (d2R(-nz:nz)); d2R = 0.
         allocate (d2Z(-nz:nz)); d2Z = 0.
         allocate (bmag_psi0(-nz:nz)); bmag_psi0 = 0.
         allocate (gradrho_psi0(-nz:nz)); gradrho_psi0 = 0.

         ! Read a radial profile variation file
         if (read_profile_variation) then
         
            ! Read the "RZ.in" file
            open (1002, file='RZ.in', status='old')
            read (1002, '(12e13.5)') rhoc0, dI, qinp, shat, d2qdr2, kappa, kapprim, &
               tri, triprim, betaprim, betadbprim, dpsip_drho_psi0
            do j = -nz, nz
               read (1002, '(5e13.5)') dum, d2R(j), d2Z(j), bmag_psi0(j), gradrho_psi0(j)
            end do
            close (1002)
            
            ! Calculate some geometric quantities
            local%qinp = qinp + shat * qinp / rhoc0 * (local%rhoc - rhoc0) + 0.5 * (local%rhoc - rhoc0)**2 * d2qdr2
            local%shat = (local%rhoc / local%qinp) * (shat * qinp / rhoc0 + (local%rhoc - rhoc0) * d2qdr2)
            local%kappa = kappa + kapprim * (local%rhoc - rhoc0)
            local%tri = tri + triprim * (local%rhoc - rhoc0)
            local%betaprim = betaprim + betadbprim * (local%rhoc - rhoc0)

            ! Set some geometric quantities at the psi0 flux surface
            local%rhoc_psi0 = rhoc0
            local%qinp_psi0 = qinp
            local%shat_psi0 = shat

            ! We do not need to load the variables at the psi0 flux surface later
            load_psi0_variables = .false.
            
         end if
      
      end subroutine read_miller_parameters_radial_variation

   end subroutine read_miller_parameters
   
   
!===============================================================================
!=============================== MILLER GEOMETRY ===============================
!===============================================================================

   !****************************************************************************
   !                          CALCULATE MILLER GEOMETRY                         
   !****************************************************************************
   ! This is the main routine that is called from geometry.f90
   !****************************************************************************
   subroutine get_miller_geometry(nzgrid, zed_in, zed_equal_arc, &
      dpsip_drho_out, dpsip_drho_psi0_out, dI_drho_out, gradrho_out, &
      bmag_out, bmag_psi0_out, &
      grady_dot_grady_out, gradx_dot_grady_out, gradx_dot_gradx_out, &
      gds23_out, gds24_out, b_dot_gradtheta_out, &
      B_times_gradB_dot_gradx_out, B_times_gradB_dot_grady_out, &
      B_times_kappa_dot_gradx_out, B_times_kappa_dot_grady_out, &
      dB_drho_out, d2B_dR_dtheta_out, d_bdotgradtheta_drho_out, &
      btor_out, rmajor_out, &
      d_Btimeskappadotgradx_drho_out, d_Btimeskappadotgrady_drho_out, &
      d_BtimesgradBdotgradx_drho_out, d_BtimesgradBdotgrady_drho_out, &
      d_gradydotgrady_drho_out, d_gradxdotgrady_drho_out, &
      d_gradxdotgradx_drho_out, djac_drho_out)

      use constants, only: pi
      use splines, only: geo_spline
      use file_utils, only: run_name

      implicit none

      ! Inputs
      integer, intent(in) :: nzgrid
      real, dimension(-nzgrid:), intent(in) :: zed_in
      logical, intent(in) :: zed_equal_arc
      
      ! Outputs 
      real, intent(out) :: dpsip_drho_out, dpsip_drho_psi0_out, dI_drho_out
      real, dimension(-nzgrid:), intent(out) :: gradrho_out
      real, dimension(-nzgrid:), intent(out) :: bmag_out, bmag_psi0_out
      real, dimension(-nzgrid:), intent(out) :: grady_dot_grady_out, gradx_dot_grady_out
      real, dimension(-nzgrid:), intent(out) :: gradx_dot_gradx_out, gds23_out, gds24_out
      real, dimension(-nzgrid:), intent(out) :: b_dot_gradtheta_out, B_times_gradB_dot_gradx_out
      real, dimension(-nzgrid:), intent(out) :: B_times_gradB_dot_grady_out, B_times_kappa_dot_gradx_out
      real, dimension(-nzgrid:), intent(out) :: B_times_kappa_dot_grady_out
      real, dimension(-nzgrid:), intent(out) :: dB_drho_out, d2B_dR_dtheta_out, d_bdotgradtheta_drho_out
      real, dimension(-nzgrid:), intent(out) :: btor_out, rmajor_out
      real, dimension(-nzgrid:), intent(out) :: d_Btimeskappadotgradx_drho_out, d_Btimeskappadotgrady_drho_out
      real, dimension(-nzgrid:), intent(out) :: d_BtimesgradBdotgradx_drho_out, d_BtimesgradBdotgrady_drho_out
      real, dimension(-nzgrid:), intent(out) :: d_gradydotgrady_drho_out, d_gradxdotgrady_drho_out
      real, dimension(-nzgrid:), intent(out) :: d_gradxdotgradx_drho_out, djac_drho_out

      ! Local variables - used for calculations
      integer :: nr, n
      integer :: i, j
      real :: rmin, dum
      real :: integral_jacrho_over_R2
      real, dimension(3) :: dr
      real, allocatable, dimension(:) :: zed_arc
      character(len=512) :: filename

      !-------------------------------------------------------------------------
      
      ! Debug message
      if (debug) write (*, *) 'geometry_miller::get_miller_geometry'
      
      ! Number of grid points used for radial derivatives (-dr, 0, +dr)
      nr = 3

      ! Allocate arrays
      call allocate_arrays(nr, nz)

      ! Get dq/∂r = q/r * s 
      dqdr = local%shat * local%qinp / local%rhoc

      ! Store (-dr, 0, +dr) as these are needed for the radial derivatives
      dr(1) = -local%dr
      dr(2) = 0.
      dr(3) = local%dr

      ! ------------------------------------------------------------------------
      !                          Function Calculations                          
      ! ------------------------------------------------------------------------
      ! Compute the functions R(r,θ) and Z(r,θ), including their radial derivatives.
      !    R(r,θ) = R0(r) + r cos[ θ + sin θ arcsin δ(r) ]
      !    Z(r,θ) = κ(r) r sin θ
      ! ------------------------------------------------------------------------
      do j = -nz, nz
         theta(j) = j * (2 * nperiod - 1) * pi / real(nz)
         do i = 1, 3
            rmin = local%rhoc + dr(i)
            Rr(i, j) = Rpos(rmin, theta(j), j)
            Zr(i, j) = Zpos(rmin, theta(j), j)
         end do
      end do

      ! ------------------------------------------------------------------------
      !                          Calculate Derivatives                          
      ! ------------------------------------------------------------------------
      ! The following sets of routines compute derivatives of R and Z with 
      ! respect to r and θ (get first and second derivatives).
      ! ------------------------------------------------------------------------
      
      ! Debug message
      if (debug) write (*, *) 'geometry_miller::radial_derivatives'

      ! Allocate arrays
      if (.not. allocated(delta_theta)) allocate (delta_theta(-nz:nz - 1))
      
      ! Get dθ as a function of theta - needed for theta derivatives
      delta_theta = theta(-nz + 1:) - theta(:nz - 1)

      ! Get ∂R/∂r and ∂Z/∂r
      call get_drho(Rr, dR_drho)
      call get_drho(Zr, dZ_drho)

      ! Get ∂R/∂θ and ∂Z/∂θ
      call get_dtheta(Rr(2, :), dR_dtheta)
      call get_dtheta(Zr(2, :), dZ_dtheta)

      ! Get second derivatives of R and Z with respect to theta, ∂²R/∂θ² and ∂²Z/∂θ²
      call get_d2dtheta2(Rr(2, :), d2R_dtheta2)
      call get_d2dtheta2(Zr(2, :), d2Z_dtheta2)
      
      ! Get mixed theta and rho derivatives of R and Z, ∂²R/∂θ∂r and ∂²Z/∂θ∂r
      call get_dtheta(dR_drho, d2R_dr_dtheta)
      call get_dtheta(dZ_drho, d2Z_dr_dtheta)

      ! ------------------------------------------------------------------------
      !                         Calculate the Jacobian                          
      ! ------------------------------------------------------------------------
      ! Get the Jacobian <jacrho> of the transformation from (r, θ, ζ) -> (R, Z, ζ)
      ! As opposed to the Jacobian <jac> for the transformation from (ψ, θ, ζ) -> (R, Z, ζ)
      !   Jρ = R * (∂R/∂r * ∂Z/∂θ - ∂R/∂θ * ∂Z/dr) 
      !   <jacrho> = Jρ/a² = (R/a) * (∂(R/a)/∂ρ * ∂(Z/a)/∂θ - ∂(R/a)/∂θ * ∂(Z/a)/dρ)
      ! ------------------------------------------------------------------------

      ! Calculate the Jacobian from (r, θ, ζ) -> (R, Z, ζ)
      jacrho = Rr(2, :) * (dR_drho * dZ_dtheta - dR_dtheta * dZ_drho)

      ! ------------------------------------------------------------------------
      !                          Calculate dψ_pol/dr                            
      ! ------------------------------------------------------------------------
      ! For an axisymmetric magnetic field we have
      !   B = ∇φ × ∇ψ + R * Btor * ∇φ
      ! 
      ! The safety factor q is defined as,
      !   q = ∂ψ_tor/∂ψ_pol = (∂ψ_tor/∂r) / (∂ψ_pol/∂r)
      !   <qinp> = <dpsitordrho> / <dpsipoldrho>
      ! 
      ! The Jacobian from (r, θ, ζ) -> (R, Z, ζ) is defined as,
      !   Jρ = R * (∂R/∂r * ∂Z/∂θ - ∂R/∂θ * ∂Z/dr) 
      ! 
      ! The derivate dψ_pol/dρ can be calculated as,
      !   dψ_{pol,SI}/dr = Rgeo * Btor / (2*pi*q) int ( Jρ / R² ) dtheta
      !   dψ_pol/dρ = (Rgeo/a) / (2*pi*q) int ( (Jρ/a²) / (R/a)² ) dtheta
      !   <dpsipoldrho> = <rgeo> * (1/<qinp) * 1/(2*pi) * int ( <jacrho> / <Rr>² ) dtheta
      !                 = <rgeo> * (1/<qinp) * <integral_jacrho_over_R2>
      ! 
      ! Here we use the following normalisations,
      !   <psitor> = ψ_pol = ψ_{pol,SI} / (a² Bref) = ψ_{pol,SI} / (a² Btor)
      !   <rho> = ρ = r/a
      !   <jacrho> = Jρ/a² = (R/a) * (∂(R/a)/∂ρ * ∂(Z/a)/∂θ - ∂(R/a)/∂θ * ∂(Z/a)/dρ)
      ! 
      ! And we've defined,
      !    <rgeo> = <bi> = (Rgeo/a) = I(ψ)/(a*Bref)
      !    <integral_jacrho_over_R2> = 1 / (2*pi) int ( <jacrho> / (R/a)² ) dtheta
      ! 
      ! If we are given <dpsitordrho>, we can use it to calculate  <rgeo>
      !   <rgeo> = <dpsipoldrho> * <qinp> / <integral_jacrho_over_R2>
      !          = <dpsitordrho> / <integral_jacrho_over_R2>
      ! 
      ! 
      ! ------------------------------------------------------------------------
      
      ! Calculate 1 / (2*pi) int ( <jacrho> / (R/a)² ) dtheta
      ! Note that theta_integrate returns the theta integral from 0 -> 2*pi
      call theta_integrate(jacrho(-nzgrid_2pi:nzgrid_2pi) / Rr(2, -nzgrid_2pi:nzgrid_2pi)**2, integral_jacrho_over_R2)
      integral_jacrho_over_R2 = integral_jacrho_over_R2 / (2.*pi)

      ! If using input.profiles, we are given <dpsitordrho> and use it to compute <rgeo>
      !   (Rgeo/a) = (∂ψ_tor/∂r) / <integral_jacrho_over_R2>
      !   ∂ψ_pol/∂ρ = (∂ψ_tor/∂ρ) / q 
      !   <rgeo> = <bi> = (Rgeo/a) = I(ψ)/(a*Bref)
      if (abs(local%dpsitordrho) > epsilon(0.)) then
         local%rgeo = local%dpsitordrho / integral_jacrho_over_R2
         dpsip_drho = local%dpsitordrho / local%qinp
         local%d2psidr2 = (local%d2psitordrho2 - local%dpsitordrho * local%shat / local%rhoc) / local%qinp
         rgeo_temp = local%rgeo
         
      ! Otherwise, we are given <rgeo> and must use it to compute <dpsip_drho>
      ! Note that without radial variation, dI = 0, and we have <rgeo> = (Rgeo/a)
      ! Calculate dψ_pol/dρ = (Rgeo/a) * (1/q) * <integral_jacrho_over_R2>
      else
         rgeo_temp = local%rgeo + dI * (rhoc - rhoc0)
         dpsip_drho = integral_jacrho_over_R2 * rgeo_temp / local%qinp
      end if
      
      ! ------------------------------------------------------------------------
      !                       Calculate Geometric Factors                       
      ! ------------------------------------------------------------------------
      ! Calculate the following geometric quantities:
      !   ∇r = R * sqrt[ (∂R/∂θ)^2 + (∂Z/∂r)^2 ] / Jρ
      !   ∇ψ_{pol,SI} = ∇r * ∂ψ_{pol,SI}/∂r
      !   r_z = 1/Jρ * (∂R/∂r ∂R/∂θ + ∂Z/∂r ∂Z/∂θ
      ! 
      ! Taking into account the normalisations, this gives,
      !   <gradrho> = ∇ρ = (R/a) * sqrt[ (∂(R/a)/∂θ)^2 + (∂(Z/a)/∂θ)^2 ] / (Jρ/a²)
      !   <gradpsi> = ∇ψ = ∇ρ * ∂ψ/∂ρ
      !   <drz> = r_z/a = 1/(Jρ/a²) * (∂(R/a)/∂ρ ∂(R/a)/∂θ + ∂(Z/a)/∂ρ ∂(Z/a)/∂θ)
      ! 
      ! Here we use the following normalisations,
      !   <psitor> = ψ_pol = ψ_{pol,SI} / (a² Bref) = ψ_{pol,SI} / (a² Btor)
      !   <rho> = ρ = r/a
      !   <jacrho> = Jρ/a² = (R/a) * (∂(R/a)/∂ρ * ∂(Z/a)/∂θ - ∂(R/a)/∂θ * ∂(Z/a)/dρ)
      ! ------------------------------------------------------------------------
      
      ! Calculate ∇ρ = (R/a) * sqrt[ (∂(R/a)/∂θ)^2 + (∂(Z/a)/∂θ)^2 ] / (Jρ/a²)
      gradrho = Rr(2, :) * sqrt(dR_dtheta**2 + dZ_dtheta**2) / jacrho
      
      ! Calculate ∇ψ = ∇ρ * ∂ψ/∂ρ
      gradpsi = gradrho * dpsip_drho

      ! Calculate <drz> = r_z/a = 1/(Jρ/a²) * (∂(R/a)/∂ρ ∂(R/a)/∂θ + ∂(Z/a)/∂ρ ∂(Z/a)/∂θ)
      ! Quantity needed in calculation of ∂I/∂r and djacrho/∂r
      drz = (dR_drho * dR_dtheta + dZ_drho * dZ_dtheta) / jacrho
      
      ! Calculate <drz_dtheta> = ∂(r_z/a)/∂θ
      call get_dtheta(drz, drz_dtheta)

      ! Get ∂I/∂r
      call get_dI_drho(dpsip_drho, gradrho, dI_drho)
      dI_drho_out = dI_drho

      ! Get ∂<jacrho>/∂ρ * ∂ψ/∂ρ and ∂<jacrho>/∂ρ
      call get_djac_drho(dpsip_drho, dI_drho, gradrho)

      ! Get ∂^2R/∂r^2 and ∂^2Z/∂r^2
      call get_d2RZdr2

      ! For radial variation we use <d2R> instead of <d2R_dr2>
      d2R = d2R_dr2
      d2Z = d2Z_dr2

      ! Get theta derivative of ∂^2R/∂r^2 and ∂^2Z/∂r^2
      call get_dtheta(d2R_dr2, d3R_dr2_dtheta)
      call get_dtheta(d2Z_dr2, d3Z_dr2_dtheta)

      ! Calculate the magnitude of B (normalized by B(psi,theta corresponding to Rgeo))
      ! <bmag> = B/B0 = sqrt(I^2 + |∇ψ|^2)/R
      bmag = sqrt(rgeo_temp**2 + gradpsi**2) / Rr(2, :)

      ! ------------------------------------------------------------------------
      !                              Radial Variation
      ! ------------------------------------------------------------------------
      
      ! The next lines are for multibox runs
      if (load_psi0_variables) then
         dpsip_drho_psi0 = dpsip_drho
         bmag_psi0 = bmag
         gradrho_psi0 = gradrho
      end if

      if (write_profile_variation) then
         open (1002, file='RZ.out', status='unknown')
         write (1002, '(12e13.5)') local%rhoc, dI_drho, local%qinp, local%shat, local%d2qdr2, &
            local%kappa, local%kapprim, &
            local%tri, local%triprim, &
            local%betaprim, local%betadbprim, dpsip_drho
         do j = -nz, nz
            write (1002, '(5e13.5)') theta(j), d2R_dr2(j), d2Z_dr2(j), bmag(j), gradrho(j)
         end do
         close (1002)
      end if
     
      ! ------------------------------------------------------------------------
      !                     Geometric Quantities + Their Derivatives
      ! ------------------------------------------------------------------------

      ! Get ∂B/∂θ
      call get_dtheta(bmag, dB_dtheta)

      ! Calculate <b_dot_gradtheta> = b . ∇θ (formerly <gradpar>)
      b_dot_gradtheta = dpsip_drho / (bmag * jacrho)
      
      ! Calculate <b_dot_gradB> = b . ∇B (formerly <gradparB>)
      b_dot_gradB = b_dot_gradtheta * dB_dtheta

      ! Get ∂|∇r|^2/∂r and ∂|∇ψ|^2/∂r
      call get_d_gradr2_drho(dpsip_drho, gradrho)

      ! Get ∂B/∂r and ∂^2 B/∂r^2
      call get_dB_drho(bmag, dI_drho)

      ! Calculate <d_bdotgradtheta_drho> = ∂(b . ∇θ)/∂r (formerly <dgradpardrho>)
      d_bdotgradtheta_drho = -b_dot_gradtheta * (dB_drho / bmag + djac_drho / jacrho)

      ! Get ∂(∂B/∂r)/∂θ
      call get_dtheta(dB_drho, d2B_dR_dtheta)

      ! Calculate <d_bdotgradB_drho> = ∂(b . ∇B)/∂r (formerly <dgradparBdrho>)
      d_bdotgradB_drho = d_bdotgradtheta_drho * dB_dtheta + b_dot_gradtheta * d2B_dR_dtheta

      ! Obtain vartheta: ϑ = I/(q * (∂ψ/∂r))  int_0^θ dθ' Jacrho / R^2
      call get_vartheta(dpsip_drho)

      ! Obtain ∂ϑ/∂r
      call get_dvartheta_dr(dpsip_drho, dI_drho)

      ! Get |∇θ|^2 , ∇r.∇θ , ∇α.∇θ , ∇r.∇α , |∇α|^2
      call get_graddotgrad(dpsip_drho, gradrho)
      
      ! The following routine calculates: 
      !   <grady_dot_grady> = ∇y.∇y
      !   <gradx_dot_grady> = ∇x.∇y 
      !   <gradx_dot_gradx> = ∇x.∇x 
      !   <gds23> = ∇θ . [∇α x (∇r x ∇α)] * (∂ψ_N/∂r)^2 / B^2
      !   <gds24> = ∇θ . [∇r x (∇r x ∇α)] * (∂ψ_N/∂r)^2 * (q/r) / B^2
      call get_gds(grady_dot_grady, gradx_dot_grady, gradx_dot_gradx, gds23, gds24)

      ! Calculate <gradalpha_times_B_dot_gradtheta> = <cross> = (∇α x B) . ∇θ
      gradalpha_times_B_dot_gradtheta = dpsip_drho * (gradrho_gradalpha * gradalpha_gradtheta - gradalpha2 * gradrho_gradtheta)

      ! Note that the definitions of <B_times_gradB_dot_grady>, <B_times_gradB_dot_gradx>, <dgbdriftdr> 
      ! and <dgbdrift0dr> are such that it gets multiplied by vperp2, not mu.  
      ! This is in contrast to Michael's GS3 notes

      ! This is bhat/B x (∇B/B) . ∇α * 2 * ∂ψ_N/∂r
      ! We redefined B_times_gradB_dot_grady = gbdrift / 2
      B_times_gradB_dot_grady = (-dB_drho + gradalpha_times_B_dot_gradtheta * dB_dtheta * dpsip_drho / bmag**2) / bmag

      ! This is bhat/B x (bhat . ∇ bhat) . ∇α * 2 * ∂ψ_N/∂r
      ! This is assuming:  β' = 4*π*p/B0^2 * (- ∂lnp/∂r ),  where p is the total pressure and B0 is the 
      ! magnetic field strength.
      ! We redefined B_times_kappa_dot_grady = cvdrift / 2
      B_times_kappa_dot_grady = (B_times_gradB_dot_grady + local%betaprim / bmag**2)

      ! This is 2 *(bhat/B x ∇B / B) . (∇q) * ∂ψ_N/∂r / (bhat . ∇B)
      ! same as usual GS2 definition once bhat . ∇B is added in below
      ! We redefined B_times_kappa_dot_gradx = cvdrift0 / 2 / shat
      B_times_kappa_dot_gradx = -2.*rgeo_temp * dqdr / bmag**2 / 2. / local%shat

      ! This is 2*∂ψ_N/∂r times the rho derivative (bhat/B x ∇B / B) . (∇q)
      d_Btimeskappadotgradx_drho = B_times_kappa_dot_gradx * 2. * local%shat * (d_bdotgradB_drho + b_dot_gradB * &
         (dI_drho / rgeo_temp - 2.*dB_drho / bmag - local%d2psidr2 / dpsip_drho)) &
         - 2.*rgeo_temp * b_dot_gradB * local%d2qdr2 / bmag**2

      ! This is 2*∂ψ_N/∂r/B times the rho derivative of (bhat x ∇B/B) . (∇q)
      ! Note that there's an extra factor of 1/B that's not expanded due to v_perp -> mu
      d_BtimesgradBdotgradx_drho = B_times_kappa_dot_gradx * 2. * local%shat * &
         (d_bdotgradB_drho + b_dot_gradB * (dI_drho / rgeo_temp - dB_drho / bmag - local%d2psidr2 / dpsip_drho)) &
         - 2.*rgeo_temp * b_dot_gradB * local%d2qdr2 / bmag**2

      B_times_kappa_dot_gradx = B_times_kappa_dot_gradx * b_dot_gradB

      ! This is 2 * ∂ψ_N/∂r * (bhat/B x ∇B/B) . (∇q)
      B_times_gradB_dot_gradx = B_times_kappa_dot_gradx

      ! Get ∂^2I/∂r^2 and ∂^2 Jacrho / dr^2
      if (debug) write (*, *) 'geometry_miller::get_d2I_dr2_d2jac_dr2'
      call get_d2I_dr2_d2jac_dr2(gradrho, dI_drho)

      ! Get ∂^2ϑ/∂r^2
      call get_d2vartheta_dr2(dpsip_drho, dI_drho)

      ! Get ∂^2B/∂r^2
      call get_d2B_dr2(bmag, dI_drho)

      ! Get ∂[(∇α x B) . ∇θ]/∂r (and others - defined in routine)
      call get_d_gradalphatimesBdotgradtheta_drho(dpsip_drho, dI_drho, gradrho)

      ! <d_BtimesgradBdotgrady_drho> = ∂/∂r [(bhat/B x (∇B) . ∇α) * 2 * ∂ψ_N/∂r] / B
      ! Note that there's an extra factor of 1/B that's not expanded due to v_perp -> mu
      d_BtimesgradBdotgrady_drho = 2.0 * (local%d2psidr2 * dB_drho / dpsip_drho - d2B_dr2 &
          + dpsip_drho * (d_gradalphatimesBdotgradtheta_drho * dB_dtheta + gradalpha_times_B_dot_gradtheta * (d2B_dR_dtheta &
          - 2.*dB_dtheta * dB_drho / bmag)) / bmag**2) / bmag
          
      ! <d_Btimeskappadotgrady_drho> = ∂/∂r (bhat/B x [bhat . ∇ bhat] . ∇α) * 2 * ∂ψ_N/∂r
      d_Btimeskappadotgrady_drho = d_BtimesgradBdotgrady_drho - B_times_gradB_dot_grady * 2. * dB_drho / bmag &
         + 2.0 * local%betadbprim / bmag**2 - 4.0 * local%betaprim * dB_drho / bmag**3 &
         - 2.0 * local%betaprim * local%d2psidr2 / dpsip_drho

      ! The next two sets of lines are corrections needed for the side boxes in a multibox simulation
      ! gbdrift  = gbdrift *(dpsip_drho_psi0/dpsip_drho)*(bmag/bmag_psi0)
      ! B_times_∇B_dot_gradx = B_times_∇B_dot_gradx*(dpsip_drho_psi0/dpsip_drho)*(bmag/bmag_psi0)
      B_times_gradB_dot_grady = B_times_gradB_dot_grady * (dpsip_drho_psi0 / dpsip_drho)
      B_times_gradB_dot_gradx = B_times_gradB_dot_gradx * (dpsip_drho_psi0 / dpsip_drho)
      B_times_kappa_dot_grady = B_times_kappa_dot_grady * (dpsip_drho_psi0 / dpsip_drho)
      B_times_kappa_dot_gradx = B_times_kappa_dot_gradx * (dpsip_drho_psi0 / dpsip_drho)

      ! <d_BtimesgradBdotgrady_drho>  = d_BtimesgradBdotgrady_drho *(dpsip_drho_psi0/dpsip_drho)*(bmag/bmag_psi0)
      ! <d_BtimesgradBdotgradx_drho> = d_BtimesgradBdotgradx_drho*(dpsip_drho_psi0/dpsip_drho)*(bmag/bmag_psi0)
      d_BtimesgradBdotgrady_drho = d_BtimesgradBdotgrady_drho * (dpsip_drho_psi0 / dpsip_drho)
      d_BtimesgradBdotgradx_drho = d_BtimesgradBdotgradx_drho * (dpsip_drho_psi0 / dpsip_drho)
      d_Btimeskappadotgrady_drho = d_Btimeskappadotgrady_drho * (dpsip_drho_psi0 / dpsip_drho)
      d_Btimeskappadotgradx_drho = d_Btimeskappadotgradx_drho * (dpsip_drho_psi0 / dpsip_drho)

      ! Interpolate here - put onto simulation zed grid
      call interpolate_functions

      ! Get the toroidal component of the magnetic field
      ! <btor> = B_toroidal/Bref = I/R Bref = <rgeo> * a/R
      btor_out = rgeo_temp / rmajor_out

      dpsip_drho_out = dpsip_drho
      dpsip_drho_psi0_out = dpsip_drho_psi0
      
      ! The <d_BtimesgradBdotgradx_drho> and <d_nablarnablatheta_drho> variables differ on macos-14 (CMake) 
      ! with respect to the other operating systems
      do i = -nz, nz
         d_BtimesgradBdotgradx_drho(i) = round(d_BtimesgradBdotgradx_drho(i), 9) 
         d_nablarnablatheta_drho(i) = round(d_nablarnablatheta_drho(i), 9)
      end do

      ! Write geometry txt files
      call write_geometry_miller_txt_files

   contains

      ! ------------------------------------------------------------------------
      !                          Interpolate functions                          
      ! ------------------------------------------------------------------------
      subroutine interpolate_functions
         
         implicit none
         
         !----------------------------------------------------------------------

         if (zed_equal_arc) then
            if (debug) write (*, *) 'geometry_miller::zed_equal_arc=.true.'
            call theta_integrate(1./b_dot_gradtheta, dum)
            b_dot_gradtheta_arc = (theta(nz) - theta(-nz)) / ((2 * nperiod - 1) * dum)
            call theta_integrate_indef(b_dot_gradtheta_arc / b_dot_gradtheta, arc)

            allocate (zed_arc(-nzgrid:nzgrid))

            call geo_spline(arc, theta, zed_in, zed_arc)
            call geo_spline(theta, gradrho_psi0, zed_arc, gradrho_out) ! gradrho is used to normalize fluxes
            call geo_spline(theta, bmag, zed_arc, bmag_out)
            call geo_spline(theta, bmag_psi0, zed_arc, bmag_psi0_out)
            call geo_spline(theta, grady_dot_grady, zed_arc, grady_dot_grady_out)
            call geo_spline(theta, gradx_dot_grady, zed_arc, gradx_dot_grady_out)
            call geo_spline(theta, gradx_dot_gradx, zed_arc, gradx_dot_gradx_out)
            call geo_spline(theta, gds23, zed_arc, gds23_out)
            call geo_spline(theta, gds24, zed_arc, gds24_out)
            call geo_spline(theta, b_dot_gradtheta_arc, zed_arc, b_dot_gradtheta_out)
            call geo_spline(theta, B_times_gradB_dot_grady, zed_arc, B_times_gradB_dot_grady_out)
            call geo_spline(theta, B_times_gradB_dot_gradx, zed_arc, B_times_gradB_dot_gradx_out)
            call geo_spline(theta, B_times_kappa_dot_grady, zed_arc, B_times_kappa_dot_grady_out)
            call geo_spline(theta, B_times_kappa_dot_gradx, zed_arc, B_times_kappa_dot_gradx_out)
            call geo_spline(theta, dB_drho, zed_arc, dB_drho_out)
            call geo_spline(theta, d2B_dR_dtheta, zed_arc, d2B_dR_dtheta_out)
            call geo_spline(theta, d_bdotgradtheta_drho, zed_arc, d_bdotgradtheta_drho_out)
            call geo_spline(theta, Rr(2, :), zed_arc, rmajor_out)
            call geo_spline(theta, d_Btimeskappadotgrady_drho, zed_arc, d_Btimeskappadotgrady_drho_out)
            call geo_spline(theta, d_BtimesgradBdotgrady_drho, zed_arc, d_BtimesgradBdotgrady_drho_out)
            call geo_spline(theta, d_Btimeskappadotgradx_drho, zed_arc, d_Btimeskappadotgradx_drho_out)
            call geo_spline(theta, d_BtimesgradBdotgradx_drho, zed_arc, d_BtimesgradBdotgradx_drho_out)
            call geo_spline(theta, d_gradydotgrady_drho, zed_arc, d_gradydotgrady_drho_out)
            call geo_spline(theta, d_gradxdotgrady_drho, zed_arc, d_gradxdotgrady_drho_out)
            call geo_spline(theta, d_gradxdotgradx_drho, zed_arc, d_gradxdotgradx_drho_out)
            call geo_spline(theta, djac_drho / dpsip_drho, zed_arc, djac_drho_out)

            deallocate (zed_arc)
         else
            if (debug) write (*, *) 'geometry_miller::zed_equal_arc=.false.'
            call geo_spline(theta, gradrho_psi0, zed_in, gradrho_out) !gradrho is used to normalize fluxes
            call geo_spline(theta, bmag, zed_in, bmag_out)
            call geo_spline(theta, bmag_psi0, zed_in, bmag_psi0_out)
            call geo_spline(theta, grady_dot_grady, zed_in, grady_dot_grady_out)
            call geo_spline(theta, gradx_dot_grady, zed_in, gradx_dot_grady_out)
            call geo_spline(theta, gradx_dot_gradx, zed_in, gradx_dot_gradx_out)
            call geo_spline(theta, gds23, zed_in, gds23_out)
            call geo_spline(theta, gds24, zed_in, gds24_out)
            call geo_spline(theta, b_dot_gradtheta, zed_in, b_dot_gradtheta_out)
            call geo_spline(theta, B_times_gradB_dot_grady, zed_in, B_times_gradB_dot_grady_out)
            call geo_spline(theta, B_times_gradB_dot_gradx, zed_in, B_times_gradB_dot_gradx_out)
            call geo_spline(theta, B_times_kappa_dot_grady, zed_in, B_times_kappa_dot_grady_out)
            call geo_spline(theta, B_times_kappa_dot_gradx, zed_in, B_times_kappa_dot_gradx_out)
            call geo_spline(theta, dB_drho, zed_in, dB_drho_out)
            call geo_spline(theta, d2B_dR_dtheta, zed_in, d2B_dR_dtheta_out)
            call geo_spline(theta, d_bdotgradtheta_drho, zed_in, d_bdotgradtheta_drho_out)
            call geo_spline(theta, Rr(2, :), zed_in, rmajor_out)
            call geo_spline(theta, d_Btimeskappadotgrady_drho, zed_in, d_Btimeskappadotgrady_drho_out)
            call geo_spline(theta, d_BtimesgradBdotgrady_drho, zed_in, d_BtimesgradBdotgrady_drho_out)
            call geo_spline(theta, d_Btimeskappadotgradx_drho, zed_in, d_Btimeskappadotgradx_drho_out)
            call geo_spline(theta, d_BtimesgradBdotgradx_drho, zed_in, d_BtimesgradBdotgradx_drho_out)
            call geo_spline(theta, d_gradydotgrady_drho, zed_in, d_gradydotgrady_drho_out)
            call geo_spline(theta, d_gradxdotgrady_drho, zed_in, d_gradxdotgrady_drho_out)
            call geo_spline(theta, d_gradxdotgradx_drho, zed_in, d_gradxdotgradx_drho_out)
            call geo_spline(theta, djac_drho / dpsip_drho, zed_in, djac_drho_out)
         end if
      end subroutine interpolate_functions

      ! ------------------------------------------------------------------------
      !                           Write to text files                           
      ! ------------------------------------------------------------------------
      subroutine write_geometry_miller_txt_files

         implicit none

         if (debug) write (*, *) 'geometry_miller::write_geometry_miller_txt_files'
         filename = "geometry_miller."//trim(run_name)//".input"
         open (1002, file=trim(filename), status='unknown')
         write (1002, '(5a16)') '#1.rhoc', '2.rmaj', '3.rgeo', '4.shift', '5.qinp'
         write (1002, '(5e16.8)') local%rhoc, local%rmaj, local%rgeo, local%shift, local%qinp
         write (1002, *)
         write (1002, '(5a16)') '#6.shat', '7.kappa', '8.kapprim', '9.tri', '10.triprim'
         write (1002, '(5e16.8)') local%shat, local%kappa, local%kapprim, local%tri, local%triprim
         write (1002, *)
         write (1002, '(5a16)') '11.betaprim', '12.dpsitordrho', '13.rhotor', '14.drhotordrho', '15.d2qdr2'
         write (1002, '(5e16.8)') local%betaprim, local%dpsitordrho, local%rhotor, local%drhotordrho, local%d2qdr2
         write (1002, *)
         write (1002, '(3a16)') '16.d2psidr2', '17.betadbprim', '18.psitor_lcfs'
         write (1002, '(3e16.8)') local%d2psidr2, local%betadbprim, local%psitor_lcfs
         close (1002)
         filename = "geometry_miller."//trim(run_name)//".output"
         open (1001, file=trim(filename), status='unknown')
         write (1001, '(a9,e18.9,a11,e18.9,a11,e18.9)') '#dI/dr: ', dI_drho, 'd2I/dr2: ', d2I_dr2, 'dpsi/dr: ', dpsip_drho
         write (1001, '(58a15)') '#1.theta', '2.R', '3.dR/dr', '4.d2R_dr2', '5.dR/dth', &
            '6.d2R_dr_dtheta', '7.dZ/dr', '8.d2Z_dr2', '9.dZ/dth', '10.d2Z_dr_dtheta', &
            '11.bmag', '12.dBdr', '13.d2B_dr2', '14.dB/dth', '15.d2B_dR_dtheta', &
            '16.vartheta', '17.dvartheta_dr', '18.d2vartheta_dr2', '19.jacr', '20.djacrdr', &
            '21.djac_drho', '22.d2jac_dr2', '23.gradrho2', '24.d_gradr2_drho', '25.gthet2', &
            '26.d_nablatheta2_drho', '27.grgthet', '28.d_nablarnablatheta_drho', &
            '29.galphgth', '30.d_nablaalphanablatheta_drho', &
            '31.grgalph', '32.d_nablaalphanablarho_drho', '33.galph2', &
            '34.d_nablaalpha2_drho', '35.gradalpha_times_B_dot_gradtheta', &
            '36.d_gradalphatimesBdotgradtheta_drho', '37.B_times_gradB_dot_gradx', &
            '38.dgbdrift0', '39.B_times_kappa_dot_gradx', '40.dcvdrift0', &
            '41.B_times_gradB_dot_grady', '42.dgbdrift', '43.B_times_kappa_dot_grady', '44.dcvdrift', '45.drz_dtheta', &
            '46.b_dot_gradtheta', '47.dgpardr', '48.b_dot_gradB', '49.dgparBdr', '50.grady_dot_grady', &
            '51.d_gradydotgrady_drho', '52.gradx_dot_grady', '53.d_gradxdotgrady_drho', &
            '54.gradx_dot_gradx', '55.d_gradxdotgradx_drho', &
            '56.gds23', '57.gds24', '58.Zr'

         if (debug) write (*, *) 'geometry_miller::write_geometry_miller_txt_files::start_loop'
         if (debug) then 
            i = nz
            write(*,*) 'a', theta(i), Rr(2, i), dR_drho(i), d2R_dr2(i), dR_dtheta(i)
            write(*,*) 'b', d2R_dr_dtheta(i), dZ_drho(i), d2Z_dr2(i), dZ_dtheta(i), d2Z_dr_dtheta(i)
            write(*,*) 'c', bmag(i), dB_drho(i), d2B_dr2(i), dB_dtheta(i), d2B_dR_dtheta(i)
            write(*,*) 'd', vartheta(i), dvartheta_dr(i), d2vartheta_dr2(i), jacrho(i), djacr_drho(i)
            write(*,*) 'e', djac_drho(i), d2jac_dr2(i), gradrho(i)**2, d_gradr2_drho(i), gradtheta2(i)
            write(*,*) 'f', d_nablatheta2_drho(i), gradrho_gradtheta(i), &
               d_nablarnablatheta_drho(i), gradalpha_gradtheta(i), d_nablaalphanablatheta_drho(i)
            write(*,*) 'g', gradrho_gradalpha(i), d_nablaalphanablarho_drho(i), &
               gradalpha2(i), d_nablaalpha2_drho(i), gradalpha_times_B_dot_gradtheta(i)
            write(*,*) 'h', d_gradalphatimesBdotgradtheta_drho(i), B_times_gradB_dot_gradx(i), &
               d_BtimesgradBdotgradx_drho(i), B_times_kappa_dot_gradx(i), d_Btimeskappadotgradx_drho(i)
            write(*,*) 'i', B_times_gradB_dot_grady(i), d_BtimesgradBdotgrady_drho(i), &
               B_times_kappa_dot_grady(i), d_Btimeskappadotgrady_drho(i), drz_dtheta(i)
            write(*,*) 'j', b_dot_gradtheta(i), d_bdotgradtheta_drho(i), b_dot_gradB(i), d_bdotgradB_drho(i), grady_dot_grady(i)
            write(*,*) 'k', d_gradydotgrady_drho(i), gradx_dot_grady(i), d_gradxdotgrady_drho(i), &
               gradx_dot_gradx(i), d_gradxdotgradx_drho(i), gds23(i), gds24(i)
            write(*,*) 'l', Zr(2, i)
         end if
         do i = -nz, nz
            write (1001, '(59e18.9)') theta(i), Rr(2, i), dR_drho(i), d2R_dr2(i), dR_dtheta(i), &
               d2R_dr_dtheta(i), dZ_drho(i), d2Z_dr2(i), dZ_dtheta(i), d2Z_dr_dtheta(i), &
               bmag(i), dB_drho(i), d2B_dr2(i), dB_dtheta(i), d2B_dR_dtheta(i), &
               vartheta(i), dvartheta_dr(i), d2vartheta_dr2(i), jacrho(i), djacr_drho(i), &
               djac_drho(i), d2jac_dr2(i), gradrho(i)**2, d_gradr2_drho(i), gradtheta2(i), &
               d_nablatheta2_drho(i), gradrho_gradtheta(i), d_nablarnablatheta_drho(i), gradalpha_gradtheta(i), d_nablaalphanablatheta_drho(i), &
               gradrho_gradalpha(i), d_nablaalphanablarho_drho(i), gradalpha2(i), d_nablaalpha2_drho(i), gradalpha_times_B_dot_gradtheta(i), &
               d_gradalphatimesBdotgradtheta_drho(i), B_times_gradB_dot_gradx(i), &
               d_BtimesgradBdotgradx_drho(i), B_times_kappa_dot_gradx(i), d_Btimeskappadotgradx_drho(i), &
               B_times_gradB_dot_grady(i), d_BtimesgradBdotgrady_drho(i), B_times_kappa_dot_grady(i), d_Btimeskappadotgrady_drho(i), drz_dtheta(i), &
               b_dot_gradtheta(i), d_bdotgradtheta_drho(i), b_dot_gradB(i), d_bdotgradB_drho(i), grady_dot_grady(i), &
               d_gradydotgrady_drho(i), gradx_dot_grady(i), d_gradxdotgrady_drho(i), gradx_dot_gradx(i), d_gradxdotgradx_drho(i), gds23(i), gds24(i), &
               Zr(2, i)
         end do
         close (1001)
         if (debug) write (*, *) 'geometry_miller::write_geometry_miller_txt_files_finished'

      end subroutine write_geometry_miller_txt_files

   end subroutine get_miller_geometry

   !****************************************************************************
   !                              Allocate arrays                               
   !****************************************************************************
   ! Periodic quantities can be computed on 2*pi grid and replicated.
   !****************************************************************************
   subroutine allocate_arrays(nr, nz)

      implicit none

      integer, intent(in) :: nr, nz

      !-------------------------------------------------------------------------

      ! Geometric quantities needed for the gyrokinetic equation (appearing in 
      ! the gyrokinetic equation, in |k_perp| and the drift frequency).
      allocate (bmag(-nz:nz)); bmag = 0.0
      allocate (b_dot_gradB(-nz:nz)); b_dot_gradB = 0.0
      allocate (b_dot_gradtheta(-nz:nz)); b_dot_gradtheta = 0.0
      allocate (b_dot_gradtheta_arc(-nz:nz)); b_dot_gradtheta_arc = 0.0
      allocate (grady_dot_grady(-nz:nz)); grady_dot_grady = 0.0
      allocate (gradx_dot_grady(-nz:nz)); gradx_dot_grady = 0.0
      allocate (gradx_dot_gradx(-nz:nz)); gradx_dot_gradx = 0.0
      allocate (B_times_gradB_dot_gradx(-nz:nz)); B_times_gradB_dot_gradx = 0.0
      allocate (B_times_gradB_dot_grady(-nz:nz)); B_times_gradB_dot_grady = 0.0
      allocate (B_times_kappa_dot_gradx(-nz:nz)); B_times_kappa_dot_gradx = 0.0
      allocate (B_times_kappa_dot_grady(-nz:nz)); B_times_kappa_dot_grady = 0.0
      
      ! For the neoclassical terms
      allocate (gds23(-nz:nz)); gds23 = 0.0
      allocate (gds24(-nz:nz)); gds24 = 0.0
      
      ! The R(r, θ) and Z(r, θ) positions of a Miller equilibrium
      ! Here we have nr = 3 for the radial derivatives
      allocate (Rr(nr, -nz:nz)); Rr = 0.0
      allocate (Zr(nr, -nz:nz)); Zr = 0.0
      
      ! Derivatives R(r, θ) and Z(r, θ) with respect to r and θ
      allocate (dR_drho(-nz:nz)); dR_drho = 0.0
      allocate (dZ_drho(-nz:nz)); dZ_drho = 0.0
      allocate (dR_dtheta(-nz:nz)); dR_dtheta  = 0.0
      allocate (dZ_dtheta(-nz:nz)); dZ_dtheta= 0.0
      allocate (d2R_dr_dtheta(-nz:nz)); d2R_dr_dtheta = 0.0
      allocate (d2Z_dr_dtheta(-nz:nz)); d2Z_dr_dtheta = 0.0
      allocate (d2R_dtheta2(-nz:nz)); d2R_dtheta2 = 0.0
      allocate (d2Z_dtheta2(-nz:nz)); d2Z_dtheta2 = 0.0
      allocate (d2R_dr2(-nz:nz)); d2R_dr2 = 0.0
      allocate (d2Z_dr2(-nz:nz)); d2Z_dr2 = 0.0
      allocate (d3R_dr2_dtheta(-nz:nz)); d3R_dr2_dtheta = 0.0
      allocate (d3Z_dr2_dtheta(-nz:nz)); d3Z_dr2_dtheta = 0.0
      
      ! The r_z variable
      allocate (drz(-nz:nz)); drz = 0.0
      allocate (drz_dtheta(-nz:nz)); drz_dtheta = 0.0
      
      ! Derivatives of the magnetic field
      allocate (dB_drho(-nz:nz)); dB_drho = 0.0
      allocate (dB_dtheta(-nz:nz)); dB_dtheta = 0.0
      allocate (d2B_dR_dtheta(-nz:nz)); d2B_dR_dtheta = 0.0
      allocate (d2B_dr2(-nz:nz)); d2B_dr2 = 0.0
      
      ! The theta and arc_length dimensions
      allocate (theta(-nz:nz)); theta = 0.0
      allocate (arc(-nz:nz)); arc = 0.0
      allocate (vartheta(-nz:nz)); vartheta = 0.0
      allocate (dvartheta_dr(-nz:nz)); dvartheta_dr = 0.0
      allocate (d2vartheta_dr2(-nz:nz)); d2vartheta_dr2 = 0.0
      
      ! Jacobians from  (r, θ, ζ) -> (R, Z, ζ)
      allocate (jacrho(-nz:nz)); jacrho = 0.0
      allocate (djac_drho(-nz:nz)); djac_drho = 0.0
      allocate (djacr_drho(-nz:nz)); djacr_drho = 0.0
      allocate (d2jac_dr2(-nz:nz)); d2jac_dr2 = 0.0
      
      ! Gradients of ∇α and ∇r
      allocate (gradrho(-nz:nz)); gradrho = 0.0
      allocate (gradpsi(-nz:nz)); gradpsi = 0.0
      allocate (gradalpha_gradtheta(-nz:nz)); gradalpha_gradtheta = 0.0
      allocate (gradrho_gradalpha(-nz:nz)); gradrho_gradalpha = 0.0
      allocate (gradrho_gradtheta(-nz:nz)); gradrho_gradtheta = 0.0
      allocate (gradtheta2(-nz:nz)); gradtheta2 = 0.0
      allocate (gradalpha2(-nz:nz)); gradalpha2 = 0.0
      allocate (gradalpha_times_B_dot_gradtheta(-nz:nz)); gradalpha_times_B_dot_gradtheta = 0.0
      allocate (d2gradpsi_dr2(-nz:nz)); d2gradpsi_dr2 = 0.0
      
      ! Radial variation
      allocate (d_bdotgradB_drho(-nz:nz)); d_bdotgradB_drho = 0.0
      allocate (d_gradydotgrady_drho(-nz:nz)); d_gradydotgrady_drho = 0.0
      allocate (d_gradxdotgrady_drho(-nz:nz)); d_gradxdotgrady_drho = 0.0
      allocate (d_gradxdotgradx_drho(-nz:nz)); d_gradxdotgradx_drho = 0.0
      allocate (d_BtimesgradBdotgradx_drho(-nz:nz)); d_BtimesgradBdotgradx_drho = 0.0
      allocate (d_BtimesgradBdotgrady_drho(-nz:nz)); d_BtimesgradBdotgrady_drho = 0.0
      allocate (d_Btimeskappadotgradx_drho(-nz:nz)); d_Btimeskappadotgradx_drho = 0.0
      allocate (d_Btimeskappadotgrady_drho(-nz:nz)); d_Btimeskappadotgrady_drho = 0.0
      allocate (d_bdotgradtheta_drho(-nz:nz)); d_bdotgradtheta_drho = 0.0
      allocate (d_gradalphatimesBdotgradtheta_drho(-nz:nz)); d_gradalphatimesBdotgradtheta_drho = 0.0
      allocate (d_gradpsi2_drho(-nz:nz)); d_gradpsi2_drho = 0.0
      allocate (d_gradr2_drho(-nz:nz)); d_gradr2_drho = 0.0
      allocate (d_nablarnablatheta_drho(-nz:nz)); d_nablarnablatheta_drho = 0.0
      allocate (d_nablatheta2_drho(-nz:nz)); d_nablatheta2_drho = 0.0
      allocate (d_nablaalphanablarho_drho(-nz:nz)); d_nablaalphanablarho_drho = 0.0
      allocate (d_nablaalphanablatheta_drho(-nz:nz)); d_nablaalphanablatheta_drho = 0.0
      allocate (d_nablaalpha2_drho(-nz:nz)); d_nablaalpha2_drho = 0.0

   end subroutine allocate_arrays
   
   
!===============================================================================
!=================================== FINISH ====================================
!===============================================================================

   !****************************************************************************
   !                         Finish the Miller geometry                         
   !****************************************************************************
   subroutine finish_miller_geometry

      implicit none

      call deallocate_arrays

   contains
   
      !---------------------- Deallocate temporary arrays ----------------------
      subroutine deallocate_arrays
      
         implicit none
         
         ! Geometric quantities needed for the gyrokinetic equation (appearing in 
         ! the gyrokinetic equation, in |k_perp| and the drift frequency).
         deallocate (bmag)
         deallocate (b_dot_gradB)
         deallocate (b_dot_gradtheta)
         deallocate (b_dot_gradtheta_arc)
         deallocate (grady_dot_grady)
         deallocate (gradx_dot_grady)
         deallocate (gradx_dot_gradx)
         deallocate (B_times_gradB_dot_gradx)
         deallocate (B_times_gradB_dot_grady)
         deallocate (B_times_kappa_dot_gradx)
         deallocate (B_times_kappa_dot_grady)
         
         ! For the neoclassical terms
         deallocate (gds23)
         deallocate (gds24)
         
         ! The R(r, θ) and Z(r, θ) positions of a Miller equilibrium
         deallocate (Rr)
         deallocate (Zr)
         
         ! Derivatives R(r, θ) and Z(r, θ) with respect to r and θ
         deallocate (dR_drho)
         deallocate (dZ_drho)
         deallocate (dR_dtheta)
         deallocate (dZ_dtheta)
         deallocate (d2R_dr_dtheta)
         deallocate (d2Z_dr_dtheta)
         deallocate (d2R_dr2)
         deallocate (d2Z_dr2)
         deallocate (d2R_dtheta2)
         deallocate (d2Z_dtheta2)
         deallocate (d3R_dr2_dtheta)
         deallocate (d3Z_dr2_dtheta)
         deallocate (d2R, d2Z)
         
         ! The r_z variable
         deallocate (drz)
         deallocate (drz_dtheta)
      
         ! Derivatives of the magnetic field
         deallocate (dB_drho)
         deallocate (dB_dtheta)
         deallocate (d2B_dR_dtheta)
         deallocate (d2B_dr2)
         
         ! The theta and arc_length dimensions
         deallocate (theta)
         deallocate (arc)
         deallocate (vartheta)
         deallocate (dvartheta_dr)
         deallocate (d2vartheta_dr2)
         
         ! Jacobians from  (r, θ, ζ) -> (R, Z, ζ)
         deallocate (jacrho)
         deallocate (djac_drho)
         deallocate (djacr_drho)
         deallocate (d2jac_dr2)
      
         ! Gradients of ∇α and ∇r
         deallocate (gradrho)
         deallocate (gradpsi)
         deallocate (gradrho_gradtheta)
         deallocate (gradalpha_gradtheta)
         deallocate (gradrho_gradalpha)
         deallocate (gradtheta2)
         deallocate (gradalpha2)
         deallocate (gradalpha_times_B_dot_gradtheta)
         deallocate (d2gradpsi_dr2)
      
         ! Radial variation
         deallocate (d_bdotgradB_drho)
         deallocate (d_bdotgradtheta_drho)
         deallocate (d_gradxdotgradx_drho)
         deallocate (d_gradydotgrady_drho)
         deallocate (d_gradxdotgrady_drho)
         deallocate (d_Btimeskappadotgradx_drho)
         deallocate (d_Btimeskappadotgrady_drho)
         deallocate (d_BtimesgradBdotgradx_drho)
         deallocate (d_BtimesgradBdotgrady_drho)
         deallocate (d_nablarnablatheta_drho)
         deallocate (d_nablatheta2_drho)
         deallocate (d_nablaalphanablarho_drho)
         deallocate (d_nablaalphanablatheta_drho)
         deallocate (d_nablaalpha2_drho)
         deallocate (d_gradalphatimesBdotgradtheta_drho)
         deallocate (d_gradpsi2_drho)
         deallocate (d_gradr2_drho)
         
         ! Radial variation
         if (allocated(delta_theta)) deallocate (delta_theta)
         if (allocated(bmag_psi0)) deallocate (bmag_psi0)
         if (allocated(gradrho_psi0)) deallocate (gradrho_psi0)
         
      end subroutine deallocate_arrays

   end subroutine finish_miller_geometry
   
!===============================================================================
!================================ CALCULATIONS =================================
!===============================================================================

   !****************************************************************************
   !                      Subroutine for computing ∂ . /∂r
   !****************************************************************************
   ! Takes in f(r), with r given at three radial locations
   ! and returns df = df/∂r at the middle radius
   !****************************************************************************
   subroutine get_drho(f, df)

      implicit none

      real, dimension(:, -nz:), intent(in) :: f
      real, dimension(-nz:), intent(out) :: df

      !-------------------------------------------------------------------------
      
      df = 0.5 * (f(3, :) - f(1, :)) / local%dr
      
   end subroutine get_drho

   !****************************************************************************
   !                       Subroutine for computing ∂ . /∂θ
   !****************************************************************************
   ! Given a function f(theta:-pi->pi), calculate theta derivative
   ! second order accurate, with equal grid spacing assumed
   ! assumes periodic in theta -- may need to change this in future
   !****************************************************************************
   subroutine get_dtheta(f, df)

      implicit none

      real, dimension(-nz:), intent(in) :: f
      real, dimension(-nz:), intent(out) :: df

      !-------------------------------------------------------------------------
      
      ! Assuming equal grid spacing in theta here
      df(-nz + 1:nz - 1) = (f(-nz + 2:) - f(:nz - 2)) / (delta_theta(:nz - 2) + delta_theta(-nz + 1:))

      ! Use periodicity at boundary
      df(-nz) = (f(-nz + 1) - f(nz - 1)) / (delta_theta(-nz) + delta_theta(nz - 1))
      df(nz) = df(-nz)

   end subroutine get_dtheta

   !****************************************************************************
   !                       Subroutine for computing ∂^2 . /∂θ^2
   !****************************************************************************
   ! Given a function f(theta), calculate second derivative of f with respect
   ! to theta. This is second order accurate, with equal grid spacing assumed
   !****************************************************************************
   subroutine get_d2dtheta2(f, d2f)

      implicit none

      real, dimension(-nz:), intent(in) :: f
      real, dimension(-nz:), intent(out) :: d2f

      !-------------------------------------------------------------------------
      
      ! Assuming equal grid spacing in theta here
      d2f(-nz + 1:nz - 1) = (f(:nz - 2) - 2.*f(-nz + 1:nz - 1) + f(-nz + 2:)) / delta_theta(-nz + 1:nz - 1)**2

      ! Use periodicity at boundary
      d2f(-nz) = (f(nz - 1) - 2.*f(-nz) + f(-nz + 1)) / delta_theta(-nz + 1)**2
      d2f(nz) = d2f(-nz)
      
   end subroutine get_d2dtheta2

   !****************************************************************************
   !                                Compute ∂I/∂r
   !****************************************************************************
   subroutine get_dI_drho(dpsip_drho, gradrho, dI_drho)

      use constants, only: pi

      implicit none

      real, intent(in) :: dpsip_drho
      real, dimension(-nz:), intent(in) :: gradrho
      real, intent(out) :: dI_drho

      real :: num1, num2, denom
      real, dimension(:), allocatable :: dum

      !-------------------------------------------------------------------------

      allocate (dum(-nz:nz)); dum = 0.

      dum = jacrho * (1.0 + (rgeo_temp / gradpsi)**2) / Rr(2, :)**2
      call theta_integrate(dum(-nzgrid_2pi:nzgrid_2pi), denom)

      dum = jacrho * (2.*dR_drho / Rr(2, :) + dqdr / local%qinp) / Rr(2, :)**2
      call theta_integrate(dum(-nzgrid_2pi:nzgrid_2pi), num1)

      ! Here, <betaprim> = β' = 4*π*p/B0^2 * (- ∂lnp/∂r )
      ! where p is the total pressure and B0 is the magnetic field strength
      dum = (-2.*(dR_dtheta * d2R_dr_dtheta + dZ_dtheta * d2Z_dr_dtheta) / jacrho &
         + drz_dtheta + local%betaprim * jacrho / dpsip_drho**2) / gradrho**2
      call theta_integrate(dum(-nzgrid_2pi:nzgrid_2pi), num2)

      dI_drho = rgeo_temp * (num1 + num2) / denom

      deallocate (dum)

   end subroutine get_dI_drho

   !****************************************************************************
   !                              Compute ∂J/∂r
   !****************************************************************************
   subroutine get_djac_drho(dpsip_drho, dI_drho, gradrho)

      implicit none

      real, intent(in) :: dpsip_drho, dI_drho
      real, dimension(-nz:), intent(in) :: gradrho

      !-------------------------------------------------------------------------

      ! This is ∂ψ/∂r * ∂(Jacrho)/∂r
      ! Here, <betaprim> = β' = 4*π*p/B0^2 * (- ∂lnp/∂r )
      ! where p is the total pressure and B0 is the magnetic field strength
      djac_drho = (Rr(2, :) / gradrho)**2 * (2.*(dR_dtheta * d2R_dr_dtheta + dZ_dtheta * d2Z_dr_dtheta) / jacrho &
         - drz_dtheta + jacrho * (rgeo_temp * dI_drho / Rr(2, :)**2 - local%betaprim) / dpsip_drho**2)

      ! This is ∂(Jacrho)/∂r
      djacr_drho = djac_drho + jacrho * local%d2psidr2 / dpsip_drho

   end subroutine get_djac_drho

   !****************************************************************************
   !                        Compute ∂^2 R/∂r^2 and ∂^2 Z/∂r^2
   !****************************************************************************
   subroutine get_d2RZdr2

      implicit none

      !-------------------------------------------------------------------------

      ! Get factor common to both d^2 R/∂r^2 and d^2 Z/∂r^2
      d2R_dr2 = ((djacr_drho - jacrho * dR_drho / Rr(2, :)) / Rr(2, :) &
          - dR_drho * d2Z_dr_dtheta + dZ_drho * d2R_dr_dtheta) / (dR_dtheta**2 + dZ_dtheta**2)

      d2Z_dr2 = -d2R_dr2 * dR_dtheta
      d2R_dr2 = d2R_dr2 * dZ_dtheta

   end subroutine get_d2RZdr2

   !****************************************************************************
   !                     Compute ∂(|∇r|^2)/∂r and ∂(|∇ψ|^2)/∂r
   !****************************************************************************
   subroutine get_d_gradr2_drho(dpsip_drho, gradrho)

      implicit none

      real, intent(in) :: dpsip_drho
      real, dimension(-nz:), intent(in) :: gradrho

      !-------------------------------------------------------------------------
      ! ∂(|∇r|^2)/∂r
      d_gradr2_drho = 2.*(gradrho**2 * (dR_drho / Rr(2, :) - djacr_drho / jacrho) &
         + (Rr(2, :) / jacrho)**2 * (dR_dtheta * d2R_dr_dtheta + d2Z_dr_dtheta * dZ_dtheta))
         
      ! ∂(|∇ψ|^2)/∂r
      d_gradpsi2_drho = 2.*(gradpsi**2 * (dR_drho / Rr(2, :) - djac_drho / jacrho) &
         + (Rr(2, :) / jacrho)**2 * (dR_dtheta * d2R_dr_dtheta + d2Z_dr_dtheta * dZ_dtheta) * dpsip_drho**2)

   end subroutine get_d_gradr2_drho

   !****************************************************************************
   !                   Compute: ∇θ.∇θ , ∇r.∇θ , ∇α.∇θ , ∇r.∇α , ∇α.∇α
   !****************************************************************************
   ! Define: 
   ! <gradtheta2> = ∇θ.∇θ = |∇θ|^2
   ! <gradrho_gradtheta> = ∇r.∇θ
   ! <gradalpha_gradtheta> = ∇α.∇θ
   ! <gradrho_gradalpha> = ∇r.∇α
   ! <gradalpha2> = ∇α.∇α = |∇α|^2
   !****************************************************************************
   subroutine get_graddotgrad(dpsip_drho, gradrho)

      implicit none

      real, intent(in) :: dpsip_drho
      real, dimension(-nz:), intent(in) :: gradrho

      !-------------------------------------------------------------------------

      ! ∇θ.∇θ
      gradtheta2 = (Rr(2, :) / jacrho)**2 * (dR_drho**2 + dZ_drho**2)

      ! ∇r.∇θ
      gradrho_gradtheta = -(Rr(2, :) / jacrho)**2 * (dR_drho * dR_dtheta + dZ_drho * dZ_dtheta)

      ! ∇α.∇θ
      gradalpha_gradtheta = -(vartheta * dqdr + local%qinp * dvartheta_dr) * gradrho_gradtheta &
         - rgeo_temp * jacrho / (dpsip_drho * Rr(2, :)**2) * gradtheta2
      
      ! ∇r.∇α
      gradrho_gradalpha = -(vartheta * dqdr + local%qinp * dvartheta_dr) * gradrho**2 &
         - rgeo_temp * jacrho / (dpsip_drho * Rr(2, :)**2) * gradrho_gradtheta
      
      ! ∇α.∇α
      gradalpha2 = (1./Rr(2, :)**2) + ((vartheta * dqdr + local%qinp * dvartheta_dr) * gradrho)**2 &
         + 2.*rgeo_temp * jacrho * (vartheta * dqdr + local%qinp * dvartheta_dr) * gradrho_gradtheta / (dpsip_drho * Rr(2, :)**2) &
         + (rgeo_temp * jacrho / (dpsip_drho * Rr(2, :)**2))**2 * gradtheta2

   end subroutine get_graddotgrad

   !****************************************************************************
   !                           Compute Gradient Factors
   !****************************************************************************
   ! Define:
   ! <grady_dot_grady> = ∇y.∇y = |∇α|^2 * (∂ψ_N/∂r)^2
   ! <gradx_dot_grady> = ∇x.∇y = (∇q . ∇α) * (∂ψ_N/∂r)^2
   ! <gradx_dot_gradx> = ∇x.∇x = |∇q|^2 * (∂ψ_N/∂r)^2
   ! <gds23> = ∇θ . [∇α x (∇r x ∇α)] * (∂ψ_N/∂r)^2 / B^2
   !         = (∇r . ∇θ * |∇α|^2 - ∇α . ∇θ * ∇r . ∇α) * (∂ψ_N/∂r)^2 / B^2
   ! <gds24> = ∇θ . [∇r x (∇r x ∇α)] * (∂ψ_N/∂r)^2 * (q/r) / B^2
   !         = (∇r . ∇θ * ∇r . ∇α - ∇α . ∇θ * |∇r|^2) * (∂ψ_N/∂r)^2 / B^2 * q/rho
   !****************************************************************************
   subroutine get_gds(grady_dot_grady, gradx_dot_grady, gradx_dot_gradx, gds23, gds24)

      implicit none

      real, dimension(-nz:), intent(out) :: grady_dot_grady, gradx_dot_grady, gradx_dot_gradx, gds23, gds24

      !-------------------------------------------------------------------------

      ! Debug
      if (debug) write (*, *) 'geometry_miller::get_gds'
      
      ! <grady_dot_grady> = ∇y.∇y = ∇|α|^2 * (∂ψ_N/∂r)^2
      ! Note: the (∂ψ_N/∂r) factor accounts for ky normalization
      grady_dot_grady = gradalpha2 * dpsip_drho_psi0**2

      ! <gradx_dot_grady> = ∇x.∇y = (∇q . ∇α) * (∂ψ_N/∂r)^2
      gradx_dot_grady = gradrho_gradalpha * dqdr * dpsip_drho_psi0**2 / shat

      ! <gradx_dot_gradx> = ∇x.∇x = |∇q|^2 * (∂ψ_N/∂r)^2
      gradx_dot_gradx = (gradrho * dpsip_drho_psi0 * dqdr)**2 / (shat)**2

      ! <gds23> = (∇r . ∇θ * |∇α|^2 - ∇α . ∇θ * ∇r . ∇α) * (∂ψ_N/∂r)^2 / B^2
      gds23 = (gradrho_gradtheta * gradalpha2 - gradalpha_gradtheta * gradrho_gradalpha) * (dpsip_drho_psi0 / bmag)**2

      ! <gds24> = (∇r . ∇θ * ∇r . ∇α - ∇α . ∇θ * |∇r|^2) * (∂ψ_N/∂r)^2 / B^2 * q/rho
      gds24 = (gradrho_gradtheta * gradrho_gradalpha - gradalpha_gradtheta * gradrho**2) &
         * (dpsip_drho_psi0 / bmag)**2 * (local%qinp_psi0 / local%rhoc_psi0)

      ! Note that kperp2 = (n0/a)^2 * (∂ψ_N/∂r)^2 * [ ∇y.∇y + 2*theta0*(∇y.∇x)*shat + theta0^2*(∇x.∇x)*shat^2 ]
      ! theta0 = kx/(ky*shat)

   end subroutine get_gds

   !****************************************************************************
   !                             Compute ∂B/∂r
   !****************************************************************************
   subroutine get_dB_drho(bmag, dI_drho)

      implicit none

      real, dimension(-nz:), intent(in) :: bmag
      real, intent(in) :: dI_drho

      !-------------------------------------------------------------------------

      ! ∂B/∂r
      dB_drho = (rgeo_temp * dI_drho + 0.5 * d_gradpsi2_drho) / (bmag * Rr(2, :)**2) &
               - bmag * dR_drho / Rr(2, :)

   end subroutine get_dB_drho

   !****************************************************************************
   !                                Compute ϑ (varthetaa)
   !****************************************************************************
   ! varthetaa is defined through α = φ − qϑ, which is used to get: 
   !                          ϑ = I/q int_0^θ dθ' Jacrho / R^2
   ! This is an integral over θ from 0 to θ
   !****************************************************************************
   subroutine get_vartheta(dpsip_drho)

      implicit none

      real, intent(in) :: dpsip_drho

      !-------------------------------------------------------------------------

      call theta_integrate_indef(jacrho / Rr(2, :)**2, vartheta)
      vartheta = rgeo_temp * vartheta / (dpsip_drho * local%qinp)

   end subroutine get_vartheta

   !****************************************************************************
   !                                   Compute ∂ϑ/∂r
   !****************************************************************************
   ! ∂ϑ/∂r= (I'/I − q'/q) * ϑ + I/q int_0^θ dθ' ∂(Jacrho / R^2)/∂r
   !****************************************************************************
   subroutine get_dvartheta_dr(dpsip_drho, dI_drho)

      implicit none

      real, intent(in) :: dpsip_drho, dI_drho

      real, dimension(-nz:nz) :: dum

      !-------------------------------------------------------------------------

      dum = rgeo_temp * jacrho * (dI_drho / rgeo_temp - dqdr / local%qinp + djac_drho / jacrho &
         - 2.*dR_drho / Rr(2, :)) / Rr(2, :)**2
      call theta_integrate_indef(dum, dvartheta_dr)

      ! ∂ϑ/∂r
      dvartheta_dr = dvartheta_dr / (dpsip_drho * local%qinp)

   end subroutine get_dvartheta_dr

   !****************************************************************************
   !                           Compute ∂^2I/∂r^2, ∂^2Jacrho/∂r^2
   !****************************************************************************
   subroutine get_d2I_dr2_d2jac_dr2(gradrho, dI_drho)

      use constants, only: pi

      implicit none

      real, dimension(-nz:), intent(in) :: gradrho
      real, intent(in) :: dI_drho

      real :: denom, num1, num2, num3, num4
      real, dimension(-nz:nz) :: tmp, tmp2

      !-------------------------------------------------------------------------

      ! Denom is the denominator in the expression for ∂^2 I /∂r^2
      tmp = jacrho / Rr(2, :)**2 * (1.0 + (rgeo_temp / gradpsi)**2)
      call theta_integrate(tmp(-nzgrid_2pi:nzgrid_2pi), denom)
      denom = denom / rgeo_temp

      d2jac_dr2 = dI_drho * rgeo_temp * jacrho / gradpsi**2 &
        * (dI_drho / rgeo_temp + djacr_drho / jacrho - d_gradpsi2_drho / gradpsi**2 - 2.*dR_drho / Rr(2, :))

      tmp = -d2jac_dr2 / Rr(2, :)**2 - dI_drho * jacrho / (rgeo_temp * Rr(2, :)**2) &
         * (djacr_drho / jacrho - dI_drho / rgeo_temp - 2.*dR_drho / Rr(2, :))
         
      call theta_integrate(tmp(-nzgrid_2pi:nzgrid_2pi), num1)

      tmp = (d2R_dr2 * dR_dtheta + dR_drho * d2R_dr_dtheta + d2Z_dr2 * dZ_dtheta + dZ_drho * d2Z_dr_dtheta) / jacrho &
            - djacr_drho * (dR_drho * dR_dtheta + dZ_drho * dZ_dtheta) / jacrho**2
            
      call get_dtheta(tmp, tmp2)
      
      tmp = (tmp2 - 2./jacrho * (-djacr_drho / jacrho * (dR_dtheta * d2R_dr_dtheta + dZ_dtheta * d2Z_dr_dtheta) &
          + d2R_dr_dtheta**2 + dR_dtheta * d3R_dr2_dtheta + d2Z_dr_dtheta**2 + dZ_dtheta * d3Z_dr2_dtheta)) / gradrho**2 &
          - d_gradr2_drho * (drz_dtheta - 2./jacrho * (dR_dtheta * d2R_dr_dtheta + dZ_dtheta * d2Z_dr_dtheta)) / gradrho**4
          
      call theta_integrate(tmp(-nzgrid_2pi:nzgrid_2pi), num2)
      
      d2jac_dr2 = d2jac_dr2 - tmp * Rr(2, :)**2

      tmp = jacrho * (local%betadbprim + local%betaprim * (djacr_drho / jacrho - d_gradpsi2_drho / gradpsi**2)) / gradpsi**2
      
      call theta_integrate(tmp(-nzgrid_2pi:nzgrid_2pi), num3)
      
      !FLAG - next negative sign?
      d2jac_dr2 = d2jac_dr2 - tmp * Rr(2, :)**2

      tmp = jacrho / Rr(2, :)**2 * (2.*d2R_dr2 / Rr(2, :) - 2.*(dR_drho / Rr(2, :))**2 &
         + local%d2qdr2 / local%qinp - (dqdr / local%qinp)**2 + (2 * dR_drho / Rr(2, :) + dqdr / local%qinp) &
         * (djacr_drho / jacrho - 2.*dR_drho / Rr(2, :)))
         
      call theta_integrate(tmp(-nzgrid_2pi:nzgrid_2pi), num4)

      d2I_dr2 = (num1 + num2 + num3 + num4) / denom

      d2jac_dr2 = d2jac_dr2 + rgeo_temp * jacrho / gradpsi**2 * d2I_dr2 + 2.*djac_drho * dR_drho / Rr(2, :)

   end subroutine get_d2I_dr2_d2jac_dr2

   !****************************************************************************
   !                              Compute ∂^2ϑ/∂r^2
   !****************************************************************************
   subroutine get_d2vartheta_dr2(dpsip_drho, dI_drho)

      implicit none

      real, intent(in) :: dpsip_drho, dI_drho

      real, dimension(-nz:nz) :: dum

      !-------------------------------------------------------------------------

      dum = rgeo_temp * jacrho / (local%qinp * dpsip_drho * Rr(2, :)**2) * ((dI_drho / rgeo_temp - dqdr / local%qinp &
         + djac_drho / jacrho - 2.*dR_drho / Rr(2, :))**2 &
         + d2I_dr2 / rgeo_temp - (dI_drho / rgeo_temp)**2 - local%d2qdr2 / local%qinp &
         + (dqdr / local%qinp)**2 + d2jac_dr2 / jacrho - (djac_drho / jacrho)**2 &
         - djac_drho * local%d2psidr2 / (dpsip_drho * jacrho) &
         - 2.*d2R_dr2 / Rr(2, :) + 2.*(dR_drho / Rr(2, :))**2)

      call theta_integrate_indef(dum, d2vartheta_dr2)

   end subroutine get_d2vartheta_dr2

   !****************************************************************************
   !                             Compute ∂^2B/∂r^2
   !****************************************************************************
   subroutine get_d2B_dr2(bmag, dI_drho)

      implicit none

      real, dimension(-nz:), intent(in) :: bmag
      real, intent(in) :: dI_drho

      !-------------------------------------------------------------------------

      d2gradpsi_dr2 = 2.*(dR_drho / Rr(2, :) - djac_drho / jacrho) * d_gradpsi2_drho &
         + 2.*gradpsi**2 * (d2R_dr2 / Rr(2, :) - (dR_drho / Rr(2, :))**2 - d2jac_dr2 / jacrho + djac_drho * djacr_drho / jacrho**2) &
         + 2.*(Rr(2, :) * gradpsi / jacrho)**2 * (d2R_dr_dtheta**2 + dR_dtheta * d3R_dr2_dtheta + d2Z_dr_dtheta**2 + dZ_dtheta * d3Z_dr2_dtheta &
         + 2.*(dR_dtheta * d2R_dr_dtheta + dZ_dtheta * d2Z_dr_dtheta) * (dR_drho / Rr(2, :) - djac_drho / jacrho))

      ! Get ∂/∂r (∂B/∂r)
      d2B_dr2 = -dB_drho * dR_drho / Rr(2, :) + bmag * (dR_drho / Rr(2, :))**2 &
         - bmag * d2R_dr2 / Rr(2, :) + 0.5 * (2.*(dI_drho**2 + rgeo_temp * d2I_dr2) + d2gradpsi_dr2) / (bmag * Rr(2, :)**2) &
         - (dB_drho + bmag * dR_drho / Rr(2, :)) * (2.*dR_drho / Rr(2, :) + dB_drho / bmag)

   end subroutine get_d2B_dr2

   !****************************************************************************
   !                  Get More Derivatives with respect to rho                  
   !****************************************************************************
   ! Things computed in this routine are: 
   ! <d_nablarnablatheta_drho> = ∂(∇r . ∇θ)/∂r
   ! <d_nablatheta2_drho> = ∂(|∇θ|^2)/∂r
   ! <d_nablaalpha2_drho> = ∂(|∇α|^2)/∂r
   ! <d_nablaalphanablarho_drho> = ∂(∇α . ∇r)/∂r
   ! <d_nablaalphanablatheta_drho> = ∂(∇α . ∇θ)/∂r
   ! <d_gradalphatimesBdotgradtheta_drho> = ∂[(∇α x B) . ∇θ)]/∂r
   ! <d_gradydotgrady_drho> = (∂ψ/∂r)^2 * ∂(|∇α|^2)/∂r
   ! <d_gradxdotgrady_drho> = (∂ψ/∂r)^2 * ∂(∇α . ∇q)/∂r
   ! <d_gradxdotgradx_drho> = (∂ψ/∂r)^2 * ∂(|∇q|^2)/∂r
   !****************************************************************************
   subroutine get_d_gradalphatimesBdotgradtheta_drho(dpsip_drho, dI_drho, gradrho)

      implicit none

      real, intent(in) :: dpsip_drho, dI_drho
      real, dimension(-nz:), intent(in) :: gradrho

      !-------------------------------------------------------------------------

      ! dgr2 = ∂/∂r (|∇r|^2)
      ! dgr2 = 2.*(Rr(2,:)/jacrho)**2*((dR_drho/Rr(2,:)-djac_drho/jacrho)*(dR_dtheta**2+dZ_dtheta**2) &
      !      + dR_dtheta*d2R_dr_dtheta + dZ_dtheta*d2Z_dr_dtheta)

      ! <d_nablarnablatheta_drho> = ∂(∇r . ∇θ)/∂r
      d_nablarnablatheta_drho = 2.*gradrho_gradtheta * (dR_drho / Rr(2, :) - djacr_drho / jacrho) &
           - (Rr(2, :) / jacrho)**2 * (d2R_dr2 * dR_dtheta + dR_drho * d2R_dr_dtheta + d2Z_dr2 * dZ_dtheta + dZ_drho * d2Z_dr_dtheta)

      ! <d_nablatheta2_drho> = ∂(|∇θ|^2)/∂r
      d_nablatheta2_drho = 2.*(Rr(2, :) / jacrho)**2 * ((dR_drho / Rr(2, :) - djacr_drho / jacrho) * (dR_drho**2 + dZ_drho**2) &
         + dR_drho * d2R_dr2 + dZ_drho * d2Z_dr2)
         
      ! <d_nablaalpha2_drho> = ∂(|∇α|^2)/∂r
      ! will later multiply it by:   0.5 * (∂ψ_N/∂r)^2
      d_nablaalpha2_drho = -2 * dR_drho / Rr(2, :)**3 + d_gradr2_drho * (vartheta * dqdr + local%qinp * dvartheta_dr)**2 &
          + (2.0 * gradrho**2 * (vartheta * dqdr + local%qinp * dvartheta_dr) &
          + 2.*rgeo_temp * jacrho * gradrho_gradtheta / (dpsip_drho * Rr(2, :)**2)) &
          * (local%d2qdr2 * vartheta + 2.*dqdr * dvartheta_dr + local%qinp * d2vartheta_dr2) &
          + 2.*(vartheta * dqdr + local%qinp * dvartheta_dr) * rgeo_temp * jacrho / (dpsip_drho * Rr(2, :)**2) &
          * (d_nablarnablatheta_drho + gradrho_gradtheta * (dI_drho / rgeo_temp + djac_drho / jacrho - 2.*dR_drho / Rr(2, :))) &
          + (rgeo_temp * jacrho / (dpsip_drho * Rr(2, :)**2))**2 * &
          (d_nablatheta2_drho + 2.*gradtheta2 * (dI_drho / rgeo_temp + djac_drho / jacrho - 2.*dR_drho / Rr(2, :)))

      ! <d_nablaalphanablarho_drho> = ∂(∇α . ∇r)/∂r
      d_nablaalphanablarho_drho = -gradrho**2 * (2.*dvartheta_dr * dqdr + vartheta * local%d2qdr2 + local%qinp * d2vartheta_dr2) &
        - d_gradr2_drho * (vartheta * dqdr + local%qinp * dvartheta_dr) - rgeo_temp * jacrho / (dpsip_drho * Rr(2, :)**2) &
        * (d_nablarnablatheta_drho + gradrho_gradtheta * (dI_drho / rgeo_temp + djac_drho / jacrho - 2.*dR_drho / Rr(2, :)))

      ! <d_nablaalphanablatheta_drho> = ∂(∇α . ∇θ)/∂r
      d_nablaalphanablatheta_drho = -gradrho_gradtheta * (2.*dvartheta_dr * dqdr + vartheta * local%d2qdr2 + local%qinp * d2vartheta_dr2) &
        - d_nablarnablatheta_drho * (vartheta * dqdr + local%qinp * dvartheta_dr) - rgeo_temp * jacrho / (dpsip_drho * Rr(2, :)**2) &
        * (d_nablatheta2_drho + gradtheta2 * (dI_drho / rgeo_temp + djac_drho / jacrho - 2.*dR_drho / Rr(2, :)))

      ! <d_gradalphatimesBdotgradtheta_drho> = ∂[(∇α x B) . ∇θ)]/∂r
      d_gradalphatimesBdotgradtheta_drho = dpsip_drho * (d_nablaalphanablarho_drho * gradalpha_gradtheta &
        + gradrho_gradalpha * d_nablaalphanablatheta_drho &
        - d_nablaalpha2_drho * gradrho_gradtheta - gradalpha2 * d_nablarnablatheta_drho) &
        + local%d2psidr2 * gradalpha_times_B_dot_gradtheta / dpsip_drho

      ! <d_gradydotgrady_drho> = (∂ψ/∂r)^2 * ∂(|∇α|^2)/∂r
      d_gradydotgrady_drho = d_nablaalpha2_drho * dpsip_drho_psi0**2

      ! <d_gradxdotgrady_drho> = (∂ψ/∂r)^2 * ∂(∇α . ∇q)/∂r
      ! Note that there will be multiplication by 2 in dist_fn.fpp
      d_gradxdotgrady_drho = (d_nablaalphanablarho_drho * dqdr + local%d2qdr2 * gradrho_gradalpha) * dpsip_drho_psi0**2

      ! <d_gradxdotgradx_drho> = (∂ψ/∂r)^2 * ∂(|∇q|^2)/∂r
      d_gradxdotgradx_drho = (dqdr**2 * d_gradr2_drho + 2.*gradrho**2 * dqdr * local%d2qdr2) * dpsip_drho_psi0**2

      ! note that dkperp2/∂r = (n0/a)^2*(∂ψ_N/∂r)^2*(d_gradydotgrady_drho + 2*theta0*d_gradxdotgrady_drho + theta0^2*d_gradxdotgradx_drho)

   end subroutine get_d_gradalphatimesBdotgradtheta_drho

   !****************************************************************************
   !                       Definite Integral in θ - from 0 to 2*π
   !****************************************************************************
   subroutine theta_integrate(integrand, integral)

      implicit none

      real, dimension(-nzgrid_2pi:), intent(in) :: integrand
      real, intent(out) :: integral

      !-------------------------------------------------------------------------

      ! Use trapezoidal rule to integrate in theta
      integral = 0.5 * sum(delta_theta(-nzgrid_2pi:nzgrid_2pi - 1) * (integrand(-nzgrid_2pi:nzgrid_2pi - 1) &
         + integrand(-nzgrid_2pi + 1:nzgrid_2pi)))

   end subroutine theta_integrate

   !****************************************************************************
   !                            Indefinite Integral in θ 
   !****************************************************************************
   subroutine theta_integrate_indef(integrand, integral)

      implicit none

      real, dimension(-nz:), intent(in) :: integrand
      real, dimension(-nz:), intent(out) :: integral

      integer :: i

      !-------------------------------------------------------------------------
      
      ! Use trapezoidal rule to integrate in theta
      integral(0) = 0.0
      do i = 1, nz
         integral(i) = integral(i - 1) + 0.5 * delta_theta(i - 1) * (integrand(i - 1) + integrand(i))
      end do
      do i = -1, -nz, -1
         integral(i) = integral(i + 1) - 0.5 * delta_theta(i) * (integrand(i + 1) + integrand(i))
      end do
      
   end subroutine theta_integrate_indef

   !****************************************************************************
   !                            COMMUNICATE PARAMETERS 
   !****************************************************************************
   ! Only needed for radial_variation simulations
   !****************************************************************************
   subroutine communicate_parameters_multibox(surf, drl, drr)
   
      use mp, only: job, scope, mp_abort
      use mp, only: crossdomprocs, subprocs
      use mp, only: send, receive
      use job_manage, only: njobs
      use common_types, only: flux_surface_type

      implicit none

      ! Arguments
      real, optional, intent(in) :: drl, drr
      type(flux_surface_type), intent(inout) :: surf

      ! Local variables
      real :: lrhoc, lqinp, lshat, lkappa, ltri, lbetaprim
      real :: rrhoc, rqinp, rshat, rkappa, rtri, rbetaprim
      real :: dqdr, rhoc_psi0, qinp_psi0, shat_psi0

      !-------------------------------------------------------------------------

      !FLAG DSO -  I think d2psidrho2 needs to be communicated, but
      !            I'm unsure what quantity needs to be updated

      if (debug) write (*, *) 'geometry_miller::communicate_parameters_multibox'
      if (job == 1) then
         dqdr = local%shat * local%qinp / local%rhoc

         lrhoc = local%rhoc + drl
         lqinp = local%qinp + drl * dqdr + 0.5 * drl**2 * local%d2qdr2
         lshat = (lrhoc / lqinp) * (dqdr + drl * local%d2qdr2)
         lkappa = kappa + drl * kapprim
         ltri = tri + drl * triprim
         lbetaprim = betaprim + drl * betadbprim

         rrhoc = local%rhoc + drr
         rqinp = local%qinp + drr * dqdr + 0.5 * drr**2 * local%d2qdr2
         rshat = (rrhoc / rqinp) * (dqdr + drr * local%d2qdr2)
         rkappa = kappa + drr * kapprim
         rtri = tri + drr * triprim
         rbetaprim = betaprim + drr * betadbprim
      end if

      call scope(crossdomprocs)

      if (job == 1) then
         call send(lrhoc, 0, 120)
         call send(lqinp, 0, 121)
         call send(lshat, 0, 122)
         call send(lkappa, 0, 123)
         call send(ltri, 0, 124)
         call send(lbetaprim, 0, 125)
         call send(local%rhoc, 0, 126)
         call send(d2R, 0, 127)
         call send(d2Z, 0, 128)
         call send(dI_drho, 0, 129)
         call send(rhoc, 0, 130)
         call send(qinp, 0, 131)
         call send(shat, 0, 132)
         call send(dpsip_drho, 0, 133)
         call send(bmag, 0, 134)
         call send(gradrho, 0, 135)

         call send(rrhoc, njobs - 1, 220)
         call send(rqinp, njobs - 1, 221)
         call send(rshat, njobs - 1, 222)
         call send(rkappa, njobs - 1, 223)
         call send(rtri, njobs - 1, 224)
         call send(rbetaprim, njobs - 1, 225)
         call send(local%rhoc, njobs - 1, 226)
         call send(d2R, njobs - 1, 227)
         call send(d2Z, njobs - 1, 228)
         call send(dI_drho, njobs - 1, 229)
         call send(rhoc, njobs - 1, 230)
         call send(qinp, njobs - 1, 231)
         call send(shat, njobs - 1, 232)
         call send(dpsip_drho, njobs - 1, 233)
         call send(bmag, njobs - 1, 234)
         call send(gradrho, njobs - 1, 235)
         rhoc_psi0 = rhoc
         qinp_psi0 = qinp
         shat_psi0 = shat
         local%rhoc_psi0 = rhoc_psi0
         local%qinp_psi0 = qinp_psi0
         local%shat_psi0 = shat_psi0
      elseif (job == 0) then
         call receive(rhoc, 1, 120)
         call receive(qinp, 1, 121)
         call receive(shat, 1, 122)
         call receive(kappa, 1, 123)
         call receive(tri, 1, 124)
         call receive(betaprim, 1, 125)
         call receive(rhoc0, 1, 126)
         call receive(d2R, 1, 127)
         call receive(d2Z, 1, 128)
         call receive(dI, 1, 129)
         call receive(rhoc_psi0, 1, 130)
         call receive(qinp_psi0, 1, 131)
         call receive(shat_psi0, 1, 132)
         call receive(dpsip_drho_psi0, 1, 133)
         call receive(bmag_psi0, 1, 134)
         call receive(gradrho_psi0, 1, 135)
         local%rhoc = rhoc
         local%qinp = qinp
         local%shat = shat
         local%kappa = kappa
         local%tri = tri
         local%betaprim = betaprim
         local%rhoc_psi0 = rhoc_psi0
         local%qinp_psi0 = qinp_psi0
         local%shat_psi0 = shat_psi0

         load_psi0_variables = .false.
      elseif (job == njobs - 1) then
         call receive(rhoc, 1, 220)
         call receive(qinp, 1, 221)
         call receive(shat, 1, 222)
         call receive(kappa, 1, 223)
         call receive(tri, 1, 224)
         call receive(betaprim, 1, 225)
         call receive(rhoc0, 1, 226)
         call receive(d2R, 1, 227)
         call receive(d2Z, 1, 228)
         call receive(dI, 1, 229)
         call receive(rhoc_psi0, 1, 230)
         call receive(qinp_psi0, 1, 231)
         call receive(shat_psi0, 1, 232)
         call receive(dpsip_drho_psi0, 1, 233)
         call receive(bmag_psi0, 1, 234)
         call receive(gradrho_psi0, 1, 235)
         local%rhoc = rhoc
         local%qinp = qinp
         local%shat = shat
         local%kappa = kappa
         local%tri = tri
         local%betaprim = betaprim
         local%rhoc_psi0 = rhoc_psi0
         local%qinp_psi0 = qinp_psi0
         local%shat_psi0 = shat_psi0

         load_psi0_variables = .false.
      end if

      surf%rhoc = local%rhoc
      surf%qinp = local%qinp
      surf%shat = local%shat
      surf%kappa = local%kappa
      surf%tri = local%tri
      surf%betaprim = local%betaprim
      surf%rhoc_psi0 = rhoc_psi0
      surf%qinp_psi0 = qinp_psi0
      surf%shat_psi0 = shat_psi0

      call scope(subprocs)

   end subroutine communicate_parameters_multibox

   !****************************************************************************
   !                               Function for R                               
   !****************************************************************************
   ! The position R(r,θ) is defined as,
   !    R(r,θ) = R0(r) + r cos( θ + sin θ arcsin δ(r) )
   ! 
   ! In order to calculate radial derivatives we also need,
   !    ∂R/∂r = R'0 + cos ( θ + sin θ arcsin δ(r) ) - r * sin θ * sin ( θ + sin θ * arcsin δ(r) ) * δ'
   ! 
   ! For short notation, let's define
   !    <cos_theta> = cos ( θ + sin θ arcsin δ(r) )
   !    <sin_theta> = sin θ * sin ( θ + sin θ * arcsin δ(r) ) * δ'
   ! 
   ! Therefore, the position R(r,θ) is calculated as,
   !    R = R0 + r * <cos_theta> + R'0 * dr + <cos_theta> * dr - r * <sin_theta> * dr
   ! 
   ! Recall the definitions of the Miller parameters,
   !    <rhoc> = r
   !    <tri> = arcsin δ(r)
   !    <triprim> = δ'
   !    <rmaj> = R0(r)
   !    <shift> = R'0(r) = Horizontal Shafranov shift
   !****************************************************************************
   function Rpos(r, theta, iz)

      use constants, only: pi

      ! Arguments
      integer, intent(in) :: iz
      real, intent(in) :: r, theta
      
      ! Local variables
      real :: cos_theta
      real :: sin_theta
      real :: dr
      integer :: i
      
      ! Result of the function
      real :: Rpos

      !-------------------------------------------------------------------------

      ! Allow for strange specification of R_ψ
      if (iz == nz + 1) then
         i = -nz
      else
         i = iz
      end if
      
      ! Within the Miller routine we have defined
      !   <dr> = 1.e-3 * (rhoc / rmaj)
      ! The radial location which is passed in here, is,
      !   <r> = local%rhoc + dr(i)
      ! Therefore, we can obtain <dr> again through,
      dr = r - local%rhoc
      
      ! Calculate cos ( θ + sin θ arcsin δ(r) )
      cos_theta = cos(theta +  sin(theta) * local%tri)
      
      ! Calculate - sin θ * sin ( θ + sin θ * arcsin δ(r) ) * δ'
      sin_theta = - sin(theta) * sin(theta + sin(theta) * local%tri) * local%triprim

      ! Calculate R = R0 + r * <cos_theta> + R'0 * dr + <cos_theta> * dr - r * <sin_theta> * dr
      ! We also add (1/2)*(r-r0)^2 * ∂^2 R/ ∂r|_r0. Note that d2R=0 unless read_profile_variation = T in input file
      Rpos = local%rmaj + local%shift * dr + cos_theta * local%rhoc &
         + (cos_theta + local%rhoc * sin_theta) * dr + 0.5 * (r - rhoc0)**2 * d2R(i)

   end function Rpos

   !****************************************************************************
   !                               Function for Z                               
   !****************************************************************************
   ! The position Z(r,θ) is defined as,
   !      Z(r,θ) = κ(r) * r * sin θ
   ! 
   ! In order to calculate radial derivatives we also need,
   !    ∂Z/∂r = (κ' * r + κ) * sin θ
   ! 
   ! Therefore, the position Z(r,θ) is calculated as,
   !    Z = κ(r) * r * sin θ + (κ' * r + κ) * sin θ * dr
   ! 
   ! Recall the definitions of the Miller parameters,
   !    <rhoc> = r
   !    <kappa> = κ(r)
   !    <kapprim> = κ'(r)
   !****************************************************************************
   function Zpos(r, theta, iz)

      ! Arguments
      integer, intent(in) :: iz
      real, intent(in) :: r, theta
      
      ! Local variables
      real :: dr
      integer :: i
      
      ! Result of the function
      real :: Zpos

      !-------------------------------------------------------------------------

      ! Allow for strange specification of Z_psi
      if (iz == nz + 1) then
         i = -nz
      else
         i = iz
      end if
      
      ! Within the Miller routine we have defined
      !   <dr> = 1.e-3 * (rhoc / rmaj)
      ! The radial location which is passed in here, is,
      !   <r> = local%rhoc + dr(i)
      ! Therefore, we can obtain <dr> again through,
      dr = r - local%rhoc
      
      ! Calculate Z = κ(r) * r * sin θ + (κ' * r + κ) * sin θ * dr
      ! Note that d2Z=0 unless read_profile_variation = T in input file
      Zpos = local%kappa * local%rhoc * sin(theta) + &
         (local%kapprim * local%rhoc + local%kappa) * sin(theta) * dr + 0.5 * (r - rhoc0)**2 * d2Z(i)

   end function Zpos

   !****************************************************************************
   !                             Function to modify θ
   !****************************************************************************
   function mod2pi(theta)

      real, intent(in) :: theta
      real :: pi, th, mod2pi
      real, parameter :: theta_tol = 1.e-6
      logical :: out

      !-------------------------------------------------------------------------

      pi = 2.*acos(0.)

      if (theta <= pi .and. theta >= -pi) then
         mod2pi = theta
         return
      end if

      if (theta - theta_tol <= pi .and. theta >= -pi) then
         mod2pi = pi
         return
      end if

      if (theta <= pi .and. theta + theta_tol >= -pi) then
         mod2pi = -pi
         return
      end if

      th = theta
      out = .true.
      do while (out)
         if (th > pi) th = th - 2.*pi
         if (th < -pi) th = th + 2.*pi
         if (th <= pi .and. th >= -pi) out = .false.
      end do
      mod2pi = th

   end function mod2pi

   !****************************************************************************
   !                                Round Variables
   !****************************************************************************
   ! This is a routine that is required for the automatic tests to work as some
   ! compilers use different rounding rules
   !****************************************************************************s
   function round(val, n)
   
      implicit none

      real :: val, round
      real :: scaled, remainder
      integer :: n
      integer :: sgn

      !-------------------------------------------------------------------------

      scaled = val*10.0**n
      sgn = sign(1.0, scaled)
      remainder = modulo(abs(scaled), 10.0)
      round = (scaled - sgn * remainder )/ 10.0**n

   end function round
   
end module geometry_miller
