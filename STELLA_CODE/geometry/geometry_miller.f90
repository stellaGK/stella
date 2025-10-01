!###############################################################################
!      Construct a magnetic equilibrium based on a set of Miller parameters     
!###############################################################################
! 
! This module constructs a magnetic equilibrium based on a set of Miller parameters.
! 
!------------------------------- Input Variables -------------------------------
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
!--------------------------------- Mathematics ---------------------------------
! Miller solves the Grad Shafranov equation locally: 
!                       R^2 ∇ · (∇ψ/R^2) = −4π R^2 p′− II′
!
! From this we can get all of the local geometry variables. 
! 
! The Miller equilibrium is a local a formalism capable of describing the local 
! magnetic geometry of a flux surface within axisymmetric systems. This approach
! ensures that the Grad-Shafranov equation is locally satisfied in ψ.
!
! The definition of the Miller equilibrium is given in cyclindrical coordinates:
!              R(r,θ) = R0(r) + rcos[θ+ sin (θarcsin \bar{δ}(r)) ]
!              Z(r,θ) = κ(r)rsin θ
! stella defines δ = arcsin \bar{δ}
!
! 
! The GS equation can be expressed using this R and Z, giving a long expression, 
! which will ~eventually~ be in the "stella_manual". For this we require terms like 
! ∂Z/∂θ, ∂R/∂θ, ∂Z/∂r, ∂R/∂r etc. (plus second derivatives). These will depend on
! quantities like R0, r, δ and κ. e.g.
!              ∂R/∂r = R0' + cos[θ+ sin(θδ)]− rsin{θsin[θ+ sin(θδ)]}δ′
!              ∂Z/∂r = (κ′r+ κ) sinθ
!              ∂R/∂θ = −sin[θ+ sin(θδ)][1 + cos(θδ)]
!              ∂Z/∂θ = κ * r * cosθ
!
! Once we have these derivatives we can also get the Jacobian for the system and 
! other geometric quantities, such as
!                           |∇r|^2 = R^2/Jr^2 [(∂Z/∂θ)^2 + (∂R/∂θ)^2]
! etc. These all depend on the form of R and Z, and require the Jacobians. 
!
!------------------------------- Code Specifics --------------------------------
! Code first reads in the input miller parameters and stores them as 
! local%(name of variable).
! 
!     - Defines dq/∂r = s * q / r
!     - Gets the forms of R and Z (using functions called Rpos and Zpos, and are
!                                  stored as Rr and Zr)
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
   public :: read_local_parameters
   public :: communicate_parameters_multibox
   public :: get_local_geo
   public :: finish_local_geo
   public :: local

   private

   real :: rhoc, rmaj, shift
   real :: kappa, kapprim
   real :: tri, triprim
   real :: betaprim, betadbprim
   real :: qinp, shat, d2qdr2
   real :: rgeo
   real :: dpsipdrho, d2psidr2, dpsipdrho_psi0
   real :: psitor_lcfs
   real :: rhotor, drhotordrho, dIdrho, dI
   real :: rhoc0
   logical :: write_profile_variation, read_profile_variation
   logical :: load_psi0_variables

   integer :: nzed_local
   integer :: nz, nz2pi

   real :: bi, dqdr, d2Idr2
   real, dimension(:), allocatable :: grho, bmag, grho_psi0, bmag_psi0
   real, dimension(:), allocatable :: b_dot_gradtheta, b_dot_gradtheta_arc, arc
   real, dimension(:), allocatable :: grady_dot_grady, gradx_dot_grady, gradx_dot_gradx
   real, dimension(:), allocatable :: gds23, gds24
   real, dimension(:), allocatable :: B_times_gradB_dot_gradx, B_times_gradB_dot_grady
   real, dimension(:), allocatable :: B_times_kappa_dot_gradx, B_times_kappa_dot_grady
   real, dimension(:), allocatable :: d2Rdth2, d2Zdth2, d2Rdrdth, d2Zdrdth
   real, dimension(:), allocatable :: gpsi, dBdrho, d2Bdrdth
   real, dimension(:), allocatable :: d_bdotgradtheta_drho, d_bdotgradB_drho, dBdth, b_dot_gradB
   real, dimension(:), allocatable :: dcvdrift0drho, dgbdrift0drho, theta
   real, dimension(:), allocatable :: varthet, dvarthdr, gradrho_gradthet, cross, d2varthdr2
   real, dimension(:), allocatable :: gradthet2, gradalph_gradthet, gradrho_gradalph, gradalph2
   real, dimension(:), allocatable :: d2Bdr2, d2Rdr2, d2Zdr2, drz, drzdth
   real, dimension(:), allocatable :: d2Rdr2dth, d2Zdr2dth, d2gpsidr2, dcrossdr
   real, dimension(:), allocatable :: dcvdriftdrho, dgbdriftdrho
   real, dimension(:), allocatable :: dgds2dr, dgds21dr, dgds22dr
   real, dimension(:), allocatable :: dgr2dr, dgpsi2dr
   real, dimension(:), allocatable :: dgrgt, dgt2, dgagr, dgagt, dga2
   real, dimension(:, :), allocatable :: Rr, Zr

   real, dimension(:), allocatable :: jacrho, delthet, djacdrho, djacrdrho
   real, dimension(:), allocatable :: d2jacdr2, dRdrho, dZdrho, dRdth, dZdth

   real, dimension(:), allocatable :: d2R, d2Z

   type(flux_surface_type) :: local

   ! Only initialise once
   logical :: initialised_miller = .false.

contains

   !============================================================================
   !=========================== READ LOCAL PARAMETERS ==========================
   !============================================================================ 
   subroutine read_local_parameters(nzed, nzgrid, local_out)

      use file_utils, only: input_unit_exist
      use common_types, only: flux_surface_type
      use namelist_geometry, only: read_namelist_geometry_miller

      implicit none

      type(flux_surface_type), intent(out) :: local_out
      integer, intent(in) :: nzed, nzgrid

      real :: dum
      integer :: np, j

      !-------------------------------------------------------------------------

      call read_namelist_geometry_miller(rhoc, rmaj, shift, qinp, shat, &
            kappa, kapprim, tri, triprim, rgeo, betaprim, &
            betadbprim, d2qdr2, d2psidr2, &
            nzed_local, read_profile_variation, write_profile_variation)

      call init_local_defaults

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

      ! Following two variables are not inputs from Miller input
      local%dr = 1.e-3 * (rhoc / rmaj)
      local%rhotor = rhotor
      local%psitor_lcfs = psitor_lcfs
      local%drhotordrho = drhotordrho
      local%dpsitordrho = 0.0
      local%d2psitordrho2 = 0.0

      ! The next three variablaes are for multibox simulations
      ! with radial variation
      local%rhoc_psi0 = rhoc
      local%qinp_psi0 = qinp
      local%shat_psi0 = shat

      ! First get nperiod corresponding to input number of grid points
      nz2pi = nzed / 2
      np = (nzgrid - nz2pi) / nzed + 1

      ! Now switch to using (possible higher resolution) local grid
      nz2pi = nzed_local / 2
      ! This is the equivalent of nzgrid on the local grid
      nz = nz2pi + nzed_local * (np - 1)

      ! ------------------------------------------------------------------------
      !                              Radial Variation
      ! ------------------------------------------------------------------------
      ! Initialise to zero. These will be overwritten if reading in from file.
      ! They are only relevant for radial_variation tests.
      allocate (d2R(-nz:nz))
      allocate (d2Z(-nz:nz))
      allocate (bmag_psi0(-nz:nz))
      allocate (grho_psi0(-nz:nz))
      d2R = 0.; d2Z = 0.; dI = 0.

      if (read_profile_variation) then
         open (1002, file='RZ.in', status='old')
         read (1002, '(12e13.5)') rhoc0, dI, qinp, shat, d2qdr2, kappa, kapprim, tri, triprim, &
            betaprim, betadbprim, dpsipdrho_psi0
         do j = -nz, nz
            read (1002, '(5e13.5)') dum, d2R(j), d2Z(j), bmag_psi0(j), grho_psi0(j)
         end do
         close (1002)
         local%qinp = qinp + shat * qinp / rhoc0 * (local%rhoc - rhoc0) &
                      + 0.5 * (local%rhoc - rhoc0)**2 * d2qdr2
         local%shat = (local%rhoc / local%qinp) &
                      * (shat * qinp / rhoc0 + (local%rhoc - rhoc0) * d2qdr2)
         local%kappa = kappa + kapprim * (local%rhoc - rhoc0)
         local%tri = tri + triprim * (local%rhoc - rhoc0)
         local%betaprim = betaprim + betadbprim * (local%rhoc - rhoc0)

         local%rhoc_psi0 = rhoc0
         local%qinp_psi0 = qinp
         local%shat_psi0 = shat

         load_psi0_variables = .false.
      end if

      local_out = local

   contains

      subroutine init_local_defaults

         implicit none

         if (initialised_miller) return
         initialised_miller = .true.
         
         rhoc0 = 0.5

         load_psi0_variables = .true.

         ! Only needed for sfincs when not using
         ! geo info from file
         rhotor = rhoc
         psitor_lcfs = 1.0
         drhotordrho = 1.0

      end subroutine init_local_defaults

   end subroutine read_local_parameters

   !============================================================================
   !========================= CALCULATE MILLER GEOMETRY ========================
   !============================================================================ 
   ! This is the main routine that is called from geometry.f90
   subroutine get_local_geo(nzed, nzgrid, zed_in, zed_equal_arc, &
      dpsipdrho_out, dpsipdrho_psi0_out, dIdrho_out, grho_out, &
      bmag_out, bmag_psi0_out, &
      grady_dot_grady_out, gradx_dot_grady_out, gradx_dot_gradx_out, &
      gds23_out, gds24_out, b_dot_gradtheta_out, &
      B_times_gradB_dot_gradx_out, B_times_gradB_dot_grady_out, &
      B_times_kappa_dot_gradx_out, B_times_kappa_dot_grady_out, &
      dBdrho_out, d2Bdrdth_out, d_bdotgradtheta_drho_out, &
      btor_out, rmajor_out, &
      dcvdrift0drho_out, dcvdriftdrho_out, &
      dgbdrift0drho_out, dgbdriftdrho_out, &
      dgds2dr_out, dgds21dr_out, &
      dgds22dr_out, djacdrho_out)

      use constants, only: pi
      use splines, only: geo_spline
      use file_utils, only: run_name

      implicit none

      ! Inputs
      integer, intent(in) :: nzed, nzgrid
      real, dimension(-nzgrid:), intent(in) :: zed_in
      logical, intent(in) :: zed_equal_arc
      ! Outputs 
      real, intent(out) :: dpsipdrho_out, dpsipdrho_psi0_out, dIdrho_out
      real, dimension(-nzgrid:), intent(out) :: grho_out
      real, dimension(-nzgrid:), intent(out) :: bmag_out, bmag_psi0_out
      real, dimension(-nzgrid:), intent(out) :: grady_dot_grady_out, gradx_dot_grady_out
      real, dimension(-nzgrid:), intent(out) :: gradx_dot_gradx_out, gds23_out, gds24_out
      real, dimension(-nzgrid:), intent(out) :: b_dot_gradtheta_out, B_times_gradB_dot_gradx_out
      real, dimension(-nzgrid:), intent(out) :: B_times_gradB_dot_grady_out, B_times_kappa_dot_gradx_out
      real, dimension(-nzgrid:), intent(out) :: B_times_kappa_dot_grady_out
      real, dimension(-nzgrid:), intent(out) :: dBdrho_out, d2Bdrdth_out, d_bdotgradtheta_drho_out
      real, dimension(-nzgrid:), intent(out) :: btor_out, rmajor_out
      real, dimension(-nzgrid:), intent(out) :: dcvdrift0drho_out, dcvdriftdrho_out
      real, dimension(-nzgrid:), intent(out) :: dgbdrift0drho_out, dgbdriftdrho_out
      real, dimension(-nzgrid:), intent(out) :: dgds2dr_out, dgds21dr_out
      real, dimension(-nzgrid:), intent(out) :: dgds22dr_out, djacdrho_out

      ! Local variables - used for calculations
      integer :: nr, np
      integer :: i, j
      real :: rmin, dum
      real, dimension(3) :: dr
      real, allocatable, dimension(:) :: zed_arc
      character(len=512) :: filename
      integer :: n

      !-------------------------------------------------------------------------
      
      ! Debug message
      if (debug) write (*, *) 'geometry_miller::get_local_geo'
      
      ! Number of grid points used for radial derivatives (-dr, 0, +dr)
      nr = 3

      ! First get nperiod corresponding to input number of grid points
      nz2pi = nzed / 2
      np = (nzgrid - nz2pi) / nzed + 1

      ! Now switch to using (possible higher resolution) local grid to compute
      ! the geometry variables on. This is to ensure that the geometry is 
      ! smooth and accurately captured
      nz2pi = nzed_local / 2
      ! This is the equivalent of nzgrid on the local grid
      nz = nz2pi + nzed_local * (np - 1)

      ! Allocate arrays
      if (debug) write (*, *) 'geometry_miller::get_local_geo'
      call allocate_arrays(nr, nz)

      ! Get dq/∂r = q/r * s 
      dqdr = local%shat * local%qinp / local%rhoc

      ! Store (-dr, 0, +dr) as these are needed for the radial derivatives
      dr(1) = -local%dr
      dr(2) = 0.
      dr(3) = local%dr

      ! ------------------------------------------------------------------------
      !                           Function Calculations
      ! ------------------------------------------------------------------------
      ! Compute the functions R(r,θ) and Z(r,θ)
      !              R(r,θ) = R0(r) + rcos[θ+ sin (θarcsin \bar{δ}(r)) ]
      !              Z(r,θ) = κ(r)rsin θ
      do j = -nz, nz
         theta(j) = j * (2 * np - 1) * pi / real(nz)
         do i = 1, 3
            rmin = local%rhoc + dr(i)
            Rr(i, j) = Rpos(rmin, theta(j), j)
            Zr(i, j) = Zpos(rmin, theta(j), j)
         end do
      end do

      ! ------------------------------------------------------------------------
      !                          Derivative Calculations
      ! ------------------------------------------------------------------------
      ! The following sets of routines compute derivatives of R and Z with 
      ! respect to r and θ (get first and second derivatives).
      ! ------------------------------------------------------------------------

      if (.not. allocated(delthet)) allocate (delthet(-nz:nz - 1))
      ! Get dθ theta as a function of theta - needed for theta derivatives
      delthet = theta(-nz + 1:) - theta(:nz - 1)
      
      ! Debug message
      if (debug) write (*, *) 'geometry_miller::radial_derivatives'

      ! Get ∂R/∂r and ∂Z/∂r
      call get_drho(Rr, dRdrho)
      call get_drho(Zr, dZdrho)

      ! Get ∂R/∂θ and ∂Z/∂θ
      call get_dthet(Rr(2, :), dRdth)
      call get_dthet(Zr(2, :), dZdth)

      ! Get second derivatives of R and Z with respect to theta
      call get_d2dthet2(Rr(2, :), d2Rdth2)
      call get_d2dthet2(Zr(2, :), d2Zdth2)
      ! Get mixed theta and rho derivatives of R and Z
      call get_dthet(dRdrho, d2Rdrdth)
      call get_dthet(dZdrho, d2Zdrdth)

      ! ------------------------------------------------------------------------
      !                      Calculation Geometric Factors 
      ! ------------------------------------------------------------------------
      ! These routines compute the geometric factors that are common in a lot of
      ! computations for the final geometric quantites - e.g. the Jacobian.
      ! ------------------------------------------------------------------------

      ! Get the Jacobian of the transformation from (r, θ, ζ) -> (R, Z, ζ)
      ! This is called jacr or jacrho in following comments as opposed
      ! to jacobian, which is for tranformation from (ψ, θ, ζ) -> (R, Z, ζ)
      ! <Jacrho> = R * (∂R/∂r * ∂Z/∂θ - ∂R/∂θ * ∂Z/dr)
      if (debug) write (*, *) 'geometry_miller::get_jacrho'
      call get_jacrho

      ! theta_integrate returns integral from 0 -> 2*pi
      ! Note that <dpsipdrho> here is an intermediary
      ! that requires manipulation to get final <dpsipdrho>
      call theta_integrate(jacrho(-nz2pi:nz2pi) / Rr(2, -nz2pi:nz2pi)**2, dpsipdrho)
      dpsipdrho = dpsipdrho / (2.*pi)

      ! Define: <bi> = I/(Btor(psi,theta of Rgeo)*a) 
      !              = Rgeo/a                ->  Where I = Btor * R is a flux function
      !         <dpsipdrho> = ∂ψ_N/∂r        -> The radial derivative of the poloidal flux function
      !                     = ∂ψ_tor/∂r / q  -> Where ψ_tor is the toroidal flux function
      ! 
      if (abs(local%dpsitordrho) > epsilon(0.)) then
         ! If using input.profiles, we are given dpsitordrho 
         ! and must use it to compute Rgeo.
         local%rgeo = local%dpsitordrho / dpsipdrho
         dpsipdrho = local%dpsitordrho / local%qinp
         local%d2psidr2 = (local%d2psitordrho2 - local%dpsitordrho * local%shat / local%rhoc) &
                          / local%qinp
         bi = local%rgeo
      else
         ! Otherwise, we are given rgeo and must use it to compute dpsipdrho
         bi = local%rgeo + dI * (rhoc - rhoc0)
         dpsipdrho = dpsipdrho * bi / local%qinp
      end if

      ! Get |∇r| and |∇ψ|
      call get_gradrho(dpsipdrho, grho)

      ! Quantity needed in calculation of ∂I/∂r and djacrho/∂r
      drz = (dRdrho * dRdth + dZdrho * dZdth) / jacrho
      call get_dthet(drz, drzdth)

      ! Get ∂I/∂r
      call get_dIdrho(dpsipdrho, grho, dIdrho)
      dIdrho_out = dIdrho

      ! Get ∂Jacrho/∂r * ∂ψ/∂r and ∂Jacrho/∂r
      call get_djacdrho(dpsipdrho, dIdrho, grho)

      ! Get ∂^2R/∂r^2 and ∂^2Z/∂r^2
      call get_d2RZdr2

      d2R = d2Rdr2
      d2Z = d2Zdr2

      ! Get theta derivative of ∂^2R/∂r^2 and ∂^2Z/∂r^2
      call get_dthet(d2Rdr2, d2Rdr2dth)
      call get_dthet(d2Zdr2, d2Zdr2dth)

      ! Calculate the magnitude of B (normalized by B(psi,theta corresponding to Rgeo))
      ! <bmag> = B/B0 = sqrt(I^2 + |∇ψ|^2)/R
      bmag = sqrt(bi**2 + gpsi**2) / Rr(2, :)

      ! ------------------------------------------------------------------------
      !                              Radial Variation
      ! ------------------------------------------------------------------------
      ! The next lines are for multibox runs
      if (load_psi0_variables) then
         dpsipdrho_psi0 = dpsipdrho
         bmag_psi0 = bmag
         grho_psi0 = grho
      end if

      if (write_profile_variation) then
         open (1002, file='RZ.out', status='unknown')
         write (1002, '(12e13.5)') local%rhoc, dIdrho, local%qinp, local%shat, local%d2qdr2, &
            local%kappa, local%kapprim, &
            local%tri, local%triprim, &
            local%betaprim, local%betadbprim, dpsipdrho
         do j = -nz, nz
            write (1002, '(5e13.5)') theta(j), d2Rdr2(j), d2Zdr2(j), bmag(j), grho(j)
         end do
         close (1002)
      end if
     
      ! ------------------------------------------------------------------------
      !                     Geometric Quantities + Their Derivatives
      ! ------------------------------------------------------------------------

      ! Get ∂B/∂θ
      call get_dthet(bmag, dbdth)

      ! Calculate <b_dot_gradtheta> = b . ∇θ (formerly <gradpar>)
      b_dot_gradtheta = dpsipdrho / (bmag * jacrho)
      
      ! Calculate <b_dot_gradB> = b . ∇B (formerly <gradparB>)
      b_dot_gradB = b_dot_gradtheta * dBdth

      ! Get ∂|∇r|^2/∂r and ∂|∇ψ|^2/∂r
      call get_dgr2dr(dpsipdrho, grho)

      ! Get ∂B/∂r and ∂^2 B/∂r^2
      call get_dBdrho(bmag, dIdrho)

      ! Calculate <d_bdotgradtheta_drho> = ∂(b . ∇θ)/∂r (formerly <dgradpardrho>)
      d_bdotgradtheta_drho = -b_dot_gradtheta * (dBdrho / bmag + djacdrho / jacrho)

      ! Get ∂(∂B/∂r)/∂θ
      call get_dthet(dBdrho, d2Bdrdth)

      ! Calculate <d_bdotgradB_drho> = ∂(b . ∇B)/∂r (formerly <dgradparBdrho>)
      d_bdotgradB_drho = d_bdotgradtheta_drho * dBdth + b_dot_gradtheta * d2Bdrdth

      ! Obtain varthet: ϑ = I/(q * (∂ψ/∂r))  int_0^θ dθ' Jacrho / R^2
      call get_varthet(dpsipdrho)

      ! Obtain ∂ϑ/∂r
      call get_dvarthdr(dpsipdrho, dIdrho)

      ! Get |∇θ|^2 , ∇r.∇θ , ∇α.∇θ , ∇r.∇α , |∇α|^2
      call get_graddotgrad(dpsipdrho, grho)
      
      if (debug) write (*, *) 'geometry_miller::get_gds'
      ! The following routine gets: 
      !           <grady_dot_grady> = ∇y.∇y
      !           <gradx_dot_grady> = ∇x.∇y 
      !           <gradx_dot_gradx> = ∇x.∇x 
      !           <gds23> = ∇θ . [∇α x (∇r x ∇α)] * (∂ψ_N/∂r)^2 / B^2
      !           <gds24> = ∇θ . [∇r x (∇r x ∇α)] * (∂ψ_N/∂r)^2 * (q/r) / B^2
      call get_gds(grady_dot_grady, gradx_dot_grady, gradx_dot_gradx, gds23, gds24)

      ! This is (∇α x B) . ∇θ
      cross = dpsipdrho * (gradrho_gradalph * gradalph_gradthet - gradalph2 * gradrho_gradthet)

      ! Note that the definitions of <B_times_gradB_dot_grady>, <B_times_gradB_dot_gradx>, <dgbdriftdr> 
      ! and <dgbdrift0dr> are such that it gets multiplied by vperp2, not mu.  
      ! This is in contrast to Michael's GS3 notes

      ! This is bhat/B x (∇B/B) . ∇α * 2 * ∂ψ_N/∂r
      ! We redefined B_times_gradB_dot_grady = gbdrift / 2
      B_times_gradB_dot_grady = (-dBdrho + cross * dBdth * dpsipdrho / bmag**2) / bmag

      ! This is bhat/B x (bhat . ∇ bhat) . ∇α * 2 * ∂ψ_N/∂r
      ! This is assuming:  β' = 4*π*p/B0^2 * (- ∂lnp/∂r ),  where p is the total pressure and B0 is the 
      ! magnetic field strength.
      ! We redefined B_times_kappa_dot_grady = cvdrift / 2
      B_times_kappa_dot_grady = (B_times_gradB_dot_grady + local%betaprim / bmag**2)

      ! This is 2 *(bhat/B x ∇B / B) . (∇q) * ∂ψ_N/∂r / (bhat . ∇B)
      ! same as usual GS2 definition once bhat . ∇B is added in below
      ! We redefined B_times_kappa_dot_gradx = cvdrift0 / 2 / shat
      B_times_kappa_dot_gradx = -2.*bi * dqdr / bmag**2 / 2. / local%shat

      ! This is 2*∂ψ_N/∂r times the rho derivative (bhat/B x ∇B / B) . (∇q)
      dcvdrift0drho = B_times_kappa_dot_gradx * 2. * local%shat * (d_bdotgradB_drho + b_dot_gradB * &
            (dIdrho / bi - 2.*dBdrho / bmag - local%d2psidr2 / dpsipdrho)) &
            - 2.*bi * b_dot_gradB * local%d2qdr2 / bmag**2

      ! This is 2*∂ψ_N/∂r/B times the rho derivative of (bhat x ∇B/B) . (∇q)
      ! Note that there's an extra factor of 1/B that's not expanded due to v_perp -> mu
      dgbdrift0drho = B_times_kappa_dot_gradx * 2. * local%shat * &
            (d_bdotgradB_drho + b_dot_gradB * (dIdrho / bi - dBdrho / bmag - local%d2psidr2 / dpsipdrho)) &
            - 2.*bi * b_dot_gradB * local%d2qdr2 / bmag**2

      B_times_kappa_dot_gradx = B_times_kappa_dot_gradx * b_dot_gradB

      ! This is 2 * ∂ψ_N/∂r * (bhat/B x ∇B/B) . (∇q)
      B_times_gradB_dot_gradx = B_times_kappa_dot_gradx

      ! Get ∂^2I/∂r^2 and ∂^2 Jacrho / dr^2
      if (debug) write (*, *) 'geometry_miller::get_d2Idr2_d2jacdr2'
      call get_d2Idr2_d2jacdr2(grho, dIdrho)

      ! Get ∂^2ϑ/∂r^2
      call get_d2varthdr2(dpsipdrho, dIdrho)

      ! Get ∂^2B/∂r^2
      call get_d2Bdr2(bmag, dIdrho)

      ! Get ∂[(∇α x B) . ∇θ]/∂r (and others - defined in routine)
      call get_dcrossdr(dpsipdrho, dIdrho, grho)

      ! <dgbdriftdrho> = ∂/∂r [(bhat/B x (∇B) . ∇α) * 2 * ∂ψ_N/∂r] / B
      ! Note that there's an extra factor of 1/B that's not expanded due to v_perp -> mu
      dgbdriftdrho = 2.0 * (local%d2psidr2 * dBdrho / dpsipdrho - d2Bdr2 &
                            + dpsipdrho * (dcrossdr * dBdth + cross * (d2Bdrdth - 2.*dBdth * dBdrho / bmag)) / bmag**2) / bmag
      ! <dcvdriftdrho> = ∂/∂r (bhat/B x [bhat . ∇ bhat] . ∇α) * 2 * ∂ψ_N/∂r
      dcvdriftdrho = dgbdriftdrho - B_times_gradB_dot_grady * 2. * dBdrho / bmag &
                     + 2.0 * local%betadbprim / bmag**2 - 4.0 * local%betaprim * dBdrho / bmag**3 &
                     - 2.0 * local%betaprim * local%d2psidr2 / dpsipdrho

      ! The next two sets of lines are corrections needed for the side boxes in a multibox simulation
      ! gbdrift  = gbdrift *(dpsipdrho_psi0/dpsipdrho)*(bmag/bmag_psi0)
      ! B_times_∇B_dot_gradx = B_times_∇B_dot_gradx*(dpsipdrho_psi0/dpsipdrho)*(bmag/bmag_psi0)
      B_times_gradB_dot_grady = B_times_gradB_dot_grady * (dpsipdrho_psi0 / dpsipdrho)
      B_times_gradB_dot_gradx = B_times_gradB_dot_gradx * (dpsipdrho_psi0 / dpsipdrho)
      B_times_kappa_dot_grady = B_times_kappa_dot_grady * (dpsipdrho_psi0 / dpsipdrho)
      B_times_kappa_dot_gradx = B_times_kappa_dot_gradx * (dpsipdrho_psi0 / dpsipdrho)

      ! <dgbdriftdrho>  = dgbdriftdrho *(dpsipdrho_psi0/dpsipdrho)*(bmag/bmag_psi0)
      ! <dgbdrift0drho> = dgbdrift0drho*(dpsipdrho_psi0/dpsipdrho)*(bmag/bmag_psi0)
      dgbdriftdrho = dgbdriftdrho * (dpsipdrho_psi0 / dpsipdrho)
      dgbdrift0drho = dgbdrift0drho * (dpsipdrho_psi0 / dpsipdrho)
      dcvdriftdrho = dcvdriftdrho * (dpsipdrho_psi0 / dpsipdrho)
      dcvdrift0drho = dcvdrift0drho * (dpsipdrho_psi0 / dpsipdrho)

      ! Interpolate here - put onto simulation zed grid
      call interpolate_functions

      ! Get the toroidal component of the magnetic field
      ! <btor> = B_toroidal/Bref = I/R Bref = <rgeo> * a/R
      btor_out = bi / rmajor_out

      dpsipdrho_out = dpsipdrho
      dpsipdrho_psi0_out = dpsipdrho_psi0

      ! Round functions - This is needed for the automatic tests
      call round_functions
      
      ! The <dgbdrift0drho> and <dgrgt> variables differ on macos-14 (CMake) 
      ! with respect to the other operating systems
      do i = -nz, nz
         dgbdrift0drho(i) = round(dgbdrift0drho(i), 9) 
         dgrgt(i) = round(dgrgt(i), 9)
      end do

      ! Write geometry txt files
      call write_geometry_miller_txt_files
      
      initialised_miller = .false.

   contains

      subroutine interpolate_functions
         
         implicit none

         if (zed_equal_arc) then
            if (debug) write (*, *) 'geometry_miller::zed_equal_arc=.true.'
            call theta_integrate(1./b_dot_gradtheta, dum)
            b_dot_gradtheta_arc = (theta(nz) - theta(-nz)) / ((2 * np - 1) * dum)
            call theta_integrate_indef(b_dot_gradtheta_arc / b_dot_gradtheta, arc)

            allocate (zed_arc(-nzgrid:nzgrid))

            call geo_spline(arc, theta, zed_in, zed_arc)
            call geo_spline(theta, grho_psi0, zed_arc, grho_out) ! grho is used to normalize fluxes
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
            call geo_spline(theta, dBdrho, zed_arc, dBdrho_out)
            call geo_spline(theta, d2Bdrdth, zed_arc, d2Bdrdth_out)
            call geo_spline(theta, d_bdotgradtheta_drho, zed_arc, d_bdotgradtheta_drho_out)
            call geo_spline(theta, Rr(2, :), zed_arc, rmajor_out)
            call geo_spline(theta, dcvdriftdrho, zed_arc, dcvdriftdrho_out)
            call geo_spline(theta, dgbdriftdrho, zed_arc, dgbdriftdrho_out)
            call geo_spline(theta, dcvdrift0drho, zed_arc, dcvdrift0drho_out)
            call geo_spline(theta, dgbdrift0drho, zed_arc, dgbdrift0drho_out)
            call geo_spline(theta, dgds2dr, zed_arc, dgds2dr_out)
            call geo_spline(theta, dgds21dr, zed_arc, dgds21dr_out)
            call geo_spline(theta, dgds22dr, zed_arc, dgds22dr_out)
            call geo_spline(theta, djacdrho / dpsipdrho, zed_arc, djacdrho_out)

            deallocate (zed_arc)
         else
            if (debug) write (*, *) 'geometry_miller::zed_equal_arc=.false.'
            call geo_spline(theta, grho_psi0, zed_in, grho_out) !grho is used to normalize fluxes
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
            call geo_spline(theta, dBdrho, zed_in, dBdrho_out)
            call geo_spline(theta, d2Bdrdth, zed_in, d2Bdrdth_out)
            call geo_spline(theta, d_bdotgradtheta_drho, zed_in, d_bdotgradtheta_drho_out)
            call geo_spline(theta, Rr(2, :), zed_in, rmajor_out)
            call geo_spline(theta, dcvdriftdrho, zed_in, dcvdriftdrho_out)
            call geo_spline(theta, dgbdriftdrho, zed_in, dgbdriftdrho_out)
            call geo_spline(theta, dcvdrift0drho, zed_in, dcvdrift0drho_out)
            call geo_spline(theta, dgbdrift0drho, zed_in, dgbdrift0drho_out)
            call geo_spline(theta, dgds2dr, zed_in, dgds2dr_out)
            call geo_spline(theta, dgds21dr, zed_in, dgds21dr_out)
            call geo_spline(theta, dgds22dr, zed_in, dgds22dr_out)
            call geo_spline(theta, djacdrho / dpsipdrho, zed_in, djacdrho_out)
         end if
      end subroutine interpolate_functions

      subroutine round_functions

         implicit none 

         ! We test stella on macos-12, macos-13, macos-14, ubuntu-20.04, ubuntu-22.04 and ubuntu-24.04
         ! and the different operating systems give differences in the last digits, so round the values.
         n = 12
         do i = -nz, nz
            theta(i) = round(theta(i), n)
            Rr(2, i) = round(Rr(2, i), n)
            dRdrho(i) = round(dRdrho(i), n)
            d2Rdr2(i) = round(d2Rdr2(i), n) 
            dRdth(i) = round(dRdth(i), n) 
            d2Rdrdth(i) = round(d2Rdrdth(i), n) 
            dZdrho(i) = round(dZdrho(i), n)
            d2Zdr2(i) = round(d2Zdr2(i), n)
            dZdth(i) = round(dZdth(i), n)
            d2Zdrdth(i) = round(d2Zdrdth(i), n) 
            bmag(i) = round(bmag(i), n) 
            dBdrho(i) = round(dBdrho(i), n)
            d2Bdr2(i) = round(d2Bdr2(i), n)
            dBdth(i) = round(dBdth(i), n) 
            d2Bdrdth(i) = round(d2Bdrdth(i), n)
            varthet(i) = round(varthet(i), n) 
            dvarthdr(i) = round(dvarthdr(i), n)
            d2varthdr2(i) = round(d2varthdr2(i), n)
            jacrho(i) = round(jacrho(i), n)
            djacrdrho(i) = round(djacrdrho(i), n) 
            djacdrho(i) = round(djacdrho(i), n) 
            d2jacdr2(i) = round(d2jacdr2(i), n)
            grho(i) = round(grho(i), n)
            dgr2dr(i) = round(dgr2dr(i), n)
            gradthet2(i) = round(gradthet2(i), n)
            dgt2(i) = round(dgt2(i), n)
            gradrho_gradthet(i) = round(gradrho_gradthet(i), n)
            gradalph_gradthet(i) = round(gradalph_gradthet(i), n)
            dgagt(i) = round(dgagt(i), n) 
            gradrho_gradalph(i) = round(gradrho_gradalph(i), n) 
            dgagr(i) = round(dgagr(i), n) 
            gradalph2(i) = round(gradalph2(i), n) 
            dga2(i) = round(dga2(i), n)
            cross(i) = round(cross(i), n) 
            dcrossdr(i) = round(dcrossdr(i), n) 
            B_times_gradB_dot_gradx(i) = round(B_times_gradB_dot_gradx(i), n) 
            B_times_kappa_dot_gradx(i) = round(B_times_kappa_dot_gradx(i), n) 
            dcvdrift0drho(i) = round(dcvdrift0drho(i), n)
            B_times_gradB_dot_grady(i) = round(B_times_gradB_dot_grady(i), n)
            dgbdriftdrho(i) = round(dgbdriftdrho(i), n)
            B_times_kappa_dot_grady(i) = round(B_times_kappa_dot_grady(i), n)
            dcvdriftdrho(i) = round(dcvdriftdrho(i), n)
            drzdth(i) = round(drzdth(i), n)
            b_dot_gradtheta(i) = round(b_dot_gradtheta(i), n)
            d_bdotgradtheta_drho(i) = round(d_bdotgradtheta_drho(i), n)
            b_dot_gradB(i) = round(b_dot_gradB(i), n)
            d_bdotgradB_drho(i) = round(d_bdotgradB_drho(i), n)
            grady_dot_grady(i) = round(grady_dot_grady(i), n)
            dgds2dr(i) = round(dgds2dr(i), n)
            gradx_dot_grady(i) = round(gradx_dot_grady(i), n)
            dgds21dr(i) = round(dgds21dr(i), n)
            gradx_dot_gradx(i) = round(gradx_dot_gradx(i), n)
            dgds22dr(i) = round(dgds22dr(i), n)
            gds23(i) = round(gds23(i), n)
            gds24(i) = round(gds24(i), n) 
            Zr(2, i) = round(Zr(2,i), n)
         end do
      end subroutine round_functions

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
         write (1002, '(5a16)') '11.betaprim', '12.dpsitordrho', '13.rhotor', &
            '14.drhotordrho', '15.d2qdr2'
         write (1002, '(5e16.8)') local%betaprim, local%dpsitordrho, local%rhotor, &
            local%drhotordrho, local%d2qdr2
         write (1002, *)
         write (1002, '(3a16)') '16.d2psidr2', '17.betadbprim', '18.psitor_lcfs'
         write (1002, '(3e16.8)') local%d2psidr2, local%betadbprim, local%psitor_lcfs
         close (1002)
         filename = "geometry_miller."//trim(run_name)//".output"
         open (1001, file=trim(filename), status='unknown')
         write (1001, '(a9,e18.9,a11,e18.9,a11,e18.9)') '#dI/dr: ', dIdrho, 'd2I/dr2: ', d2Idr2, 'dpsi/dr: ', dpsipdrho
         write (1001, '(58a15)') '#1.theta', '2.R', '3.dR/dr', '4.d2Rdr2', '5.dR/dth', &
            '6.d2Rdrdth', '7.dZ/dr', '8.d2Zdr2', '9.dZ/dth', '10.d2Zdrdth', &
            '11.bmag', '12.dBdr', '13.d2Bdr2', '14.dB/dth', '15.d2Bdrdth', &
            '16.varthet', '17.dvarthdr', '18.d2varthdr2', '19.jacr', '20.djacrdr', &
            '21.djacdrho', '22.d2jacdr2', '23.grho2', '24.dgr2dr', '25.gthet2', &
            '26.dgt2', '27.grgthet', '28.dgrgt', '29.galphgth', '30.dgagt', &
            '31.grgalph', '32.dgagr', '33.galph2', '34.dga2', '35.cross', &
            '36.dcrossdr', '37.B_times_gradB_dot_gradx', '38.dgbdrift0', '39.B_times_kappa_dot_gradx', '40.dcvdrift0', &
            '41.B_times_gradB_dot_grady', '42.dgbdrift', '43.B_times_kappa_dot_grady', '44.dcvdrift', '45.drzdth', &
            '46.b_dot_gradtheta', '47.dgpardr', '48.b_dot_gradB', '49.dgparBdr', '50.grady_dot_grady', &
            '51.dgds2dr', '52.gradx_dot_grady', '53.dgds21dr', '54.gradx_dot_gradx', '55.dgds22dr', &
            '56.gds23', '57.gds24', '58.Zr'

         if (debug) write (*, *) 'geometry_miller::write_geometry_miller_txt_files::start_loop'
         if (debug) then 
            i = nz
            write(*,*) 'a', theta(i), Rr(2, i), dRdrho(i), d2Rdr2(i), dRdth(i)
            write(*,*) 'b', d2Rdrdth(i), dZdrho(i), d2Zdr2(i), dZdth(i), d2Zdrdth(i)
            write(*,*) 'c', bmag(i), dBdrho(i), d2Bdr2(i), dBdth(i), d2Bdrdth(i)
            write(*,*) 'd', varthet(i), dvarthdr(i), d2varthdr2(i), jacrho(i), djacrdrho(i)
            write(*,*) 'e', djacdrho(i), d2jacdr2(i), grho(i)**2, dgr2dr(i), gradthet2(i)
            write(*,*) 'f', dgt2(i), gradrho_gradthet(i), dgrgt(i), gradalph_gradthet(i), dgagt(i)
            write(*,*) 'g', gradrho_gradalph(i), dgagr(i), gradalph2(i), dga2(i), cross(i)
            write(*,*) 'h', dcrossdr(i), B_times_gradB_dot_gradx(i), dgbdrift0drho(i), B_times_kappa_dot_gradx(i), dcvdrift0drho(i)
            write(*,*) 'i', B_times_gradB_dot_grady(i), dgbdriftdrho(i), B_times_kappa_dot_grady(i), dcvdriftdrho(i), drzdth(i)
            write(*,*) 'j', b_dot_gradtheta(i), d_bdotgradtheta_drho(i), b_dot_gradB(i), d_bdotgradB_drho(i), grady_dot_grady(i)
            write(*,*) 'k', dgds2dr(i), gradx_dot_grady(i), dgds21dr(i), gradx_dot_gradx(i), dgds22dr(i), gds23(i), gds24(i)
            write(*,*) 'l', Zr(2, i)
         end if
         do i = -nz, nz
            write (1001, '(59e18.9)') theta(i), Rr(2, i), dRdrho(i), d2Rdr2(i), dRdth(i), &
               d2Rdrdth(i), dZdrho(i), d2Zdr2(i), dZdth(i), d2Zdrdth(i), &
               bmag(i), dBdrho(i), d2Bdr2(i), dBdth(i), d2Bdrdth(i), &
               varthet(i), dvarthdr(i), d2varthdr2(i), jacrho(i), djacrdrho(i), &
               djacdrho(i), d2jacdr2(i), grho(i)**2, dgr2dr(i), gradthet2(i), &
               dgt2(i), gradrho_gradthet(i), dgrgt(i), gradalph_gradthet(i), dgagt(i), &
               gradrho_gradalph(i), dgagr(i), gradalph2(i), dga2(i), cross(i), &
               dcrossdr(i), B_times_gradB_dot_gradx(i), dgbdrift0drho(i), B_times_kappa_dot_gradx(i), dcvdrift0drho(i), &
               B_times_gradB_dot_grady(i), dgbdriftdrho(i), B_times_kappa_dot_grady(i), dcvdriftdrho(i), drzdth(i), &
               b_dot_gradtheta(i), d_bdotgradtheta_drho(i), b_dot_gradB(i), d_bdotgradB_drho(i), grady_dot_grady(i), &
               dgds2dr(i), gradx_dot_grady(i), dgds21dr(i), gradx_dot_gradx(i), dgds22dr(i), gds23(i), gds24(i), &
               Zr(2, i)
         end do
         close (1001)
         if (debug) write (*, *) 'geometry_miller::write_geometry_miller_txt_files_finished'

      end subroutine write_geometry_miller_txt_files

   end subroutine get_local_geo

   !============================================================================
   !============================== ALLOCATE ARRAYS =============================
   !============================================================================ 
   subroutine allocate_arrays(nr, nz)

      implicit none

      integer, intent(in) :: nr, nz

      !-------------------------------------------------------------------------

      ! Periodic quantities can be computed on 2*pi grid and replicated
      allocate (grho(-nz:nz), bmag(-nz:nz), b_dot_gradtheta(-nz:nz)); grho = 0.0; bmag = 0.0; b_dot_gradtheta = 0.0
      allocate (grady_dot_grady(-nz:nz), gradx_dot_grady(-nz:nz), gradx_dot_gradx(-nz:nz), gds23(-nz:nz), gds24(-nz:nz))
      grady_dot_grady = 0.0; gradx_dot_grady = 0.0; gradx_dot_gradx = 0.0; gds23 = 0.0; gds24 = 0.0
      allocate (B_times_gradB_dot_gradx(-nz:nz), B_times_gradB_dot_grady(-nz:nz)); B_times_gradB_dot_gradx = 0.0; B_times_gradB_dot_grady = 0.0
      allocate (B_times_kappa_dot_gradx(-nz:nz), B_times_kappa_dot_grady(-nz:nz)); B_times_kappa_dot_gradx = 0.0; B_times_kappa_dot_grady = 0.0
      allocate (Rr(nr, -nz:nz), Zr(nr, -nz:nz)); Rr = 0.0; Zr = 0.0
      allocate (jacrho(-nz:nz), djacdrho(-nz:nz), djacrdrho(-nz:nz), d2jacdr2(-nz:nz))
      jacrho = 0.0; djacdrho = 0.0; djacrdrho = 0.0; d2jacdr2 = 0.0
      allocate (d2Rdrdth(-nz:nz), d2Zdrdth(-nz:nz), gpsi(-nz:nz)); d2Rdrdth = 0.0; d2Zdrdth = 0.0; gpsi = 0.0
      allocate (dBdrho(-nz:nz), d_bdotgradtheta_drho(-nz:nz)); dBdrho = 0.0; d_bdotgradtheta_drho = 0.0
      allocate (d2Bdrdth(-nz:nz), d_bdotgradB_drho(-nz:nz), dBdth(-nz:nz), b_dot_gradB(-nz:nz))
      d2Bdrdth = 0.0; d_bdotgradB_drho = 0.0; dBdth = 0.0; b_dot_gradB = 0.0
      allocate (dcvdrift0drho(-nz:nz), dgbdrift0drho(-nz:nz)); dcvdrift0drho = 0.0; dgbdrift0drho = 0.0
      allocate (theta(-nz:nz)); theta = 0.0
      allocate (b_dot_gradtheta_arc(-nz:nz)); b_dot_gradtheta_arc = 0.0
      allocate (arc(-nz:nz)); arc = 0.0
      allocate (dRdrho(-nz:nz), dZdrho(-nz:nz), dRdth(-nz:nz), dZdth(-nz:nz))
      dRdrho = 0.0; dZdrho = 0.0; dRdth  = 0.0; dZdth= 0.0
      allocate (gradrho_gradthet(-nz:nz), gradthet2(-nz:nz), dgr2dr(-nz:nz), dgpsi2dr(-nz:nz))
      gradrho_gradthet = 0.0; gradthet2 = 0.0; dgr2dr = 0.0; dgpsi2dr = 0.0
      allocate (dgrgt(-nz:nz), dgt2(-nz:nz), dgagr(-nz:nz), dgagt(-nz:nz), dga2(-nz:nz))
      dgrgt = 0.0; dgt2 = 0.0; dgagr = 0.0; dgagt = 0.0; dga2 = 0.0
      allocate (d2Rdr2(-nz:nz), d2Zdr2(-nz:nz), d2Bdr2(-nz:nz)); d2Rdr2 = 0.0; d2Zdr2 = 0.0; d2Bdr2 = 0.0
      allocate (drz(-nz:nz), drzdth(-nz:nz), d2Rdr2dth(-nz:nz), d2Zdr2dth(-nz:nz))
      drz = 0.0; drzdth = 0.0; d2Rdr2dth = 0.0; d2Zdr2dth = 0.0
      allocate (d2Rdth2(-nz:nz), d2Zdth2(-nz:nz)); d2Rdth2 = 0.0; d2Zdth2 = 0.0
      allocate (d2gpsidr2(-nz:nz)); d2gpsidr2 = 0.0
      allocate (gradalph_gradthet(-nz:nz), gradalph2(-nz:nz), gradrho_gradalph(-nz:nz))
      gradalph_gradthet = 0.0; gradalph2 = 0.0; gradrho_gradalph = 0.0
      allocate (dgds2dr(-nz:nz), dgds21dr(-nz:nz)); dgds2dr = 0.0; dgds21dr = 0.0
      allocate (dgds22dr(-nz:nz)); dgds22dr = 0.0
      allocate (dcvdriftdrho(-nz:nz), dgbdriftdrho(-nz:nz)); dcvdriftdrho = 0.0; dgbdriftdrho = 0.0
      allocate (varthet(-nz:nz), dvarthdr(-nz:nz), d2varthdr2(-nz:nz)); varthet = 0.0; dvarthdr = 0.0; d2varthdr2 = 0.0
      allocate (cross(-nz:nz)); cross = 0.0
      allocate (dcrossdr(-nz:nz)); dcrossdr = 0.0

   end subroutine allocate_arrays

   !===============================================================================
   !==================================== FINISH ===================================
   !===============================================================================
   subroutine finish_local_geo

      implicit none

      call deallocate_arrays

   contains

      deallocate (grho)
      deallocate (bmag)
      deallocate (b_dot_gradtheta)
      deallocate (grady_dot_grady)
      deallocate (gradx_dot_grady)
      deallocate (gradx_dot_gradx)
      deallocate (gds23)
      deallocate (gds24)
      deallocate (B_times_gradB_dot_gradx)
      deallocate (B_times_gradB_dot_grady)
      deallocate (B_times_kappa_dot_gradx)
      deallocate (B_times_kappa_dot_grady)
      deallocate (Rr, Zr)
      deallocate (jacrho, djacdrho, djacrdrho, d2jacdr2)
      deallocate (d2Rdrdth, d2Zdrdth, gpsi)
      deallocate (dBdrho, d_bdotgradtheta_drho)
      deallocate (d2Bdrdth, d_bdotgradB_drho, dBdth, b_dot_gradB)
      deallocate (dcvdrift0drho, dgbdrift0drho)
      deallocate (theta)
      deallocate (b_dot_gradtheta_arc)
      deallocate (arc)
      deallocate (dRdrho, dZdrho, dRdth, dZdth)
      deallocate (gradrho_gradthet, gradthet2, dgr2dr, dgpsi2dr)
      deallocate (dgrgt, dgt2, dgagr, dgagt, dga2)
      deallocate (d2Rdr2, d2Zdr2, d2Bdr2)
      deallocate (drz, drzdth, d2Rdr2dth, d2Zdr2dth)
      deallocate (d2Rdth2, d2Zdth2)
      deallocate (d2gpsidr2)
      deallocate (gradalph_gradthet, gradalph2, gradrho_gradalph)
      deallocate (dgds2dr, dgds21dr)
      deallocate (dgds22dr)
      deallocate (dcvdriftdrho, dgbdriftdrho)
      deallocate (varthet, dvarthdr, d2varthdr2)
      deallocate (cross)
      deallocate (dcrossdr)
      deallocate (d2R, d2Z)
      if (allocated(delthet)) deallocate (delthet)
      if (allocated(bmag_psi0)) deallocate (bmag_psi0)
      if (allocated(grho_psi0)) deallocate (grho_psi0)

   end subroutine finish_local_geo
   
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
   subroutine get_dthet(f, df)

      implicit none

      real, dimension(-nz:), intent(in) :: f
      real, dimension(-nz:), intent(out) :: df

      !-------------------------------------------------------------------------
      
      ! Assuming equal grid spacing in theta here
      df(-nz + 1:nz - 1) = (f(-nz + 2:) - f(:nz - 2)) / (delthet(:nz - 2) + delthet(-nz + 1:))

      ! Use periodicity at boundary
      df(-nz) = (f(-nz + 1) - f(nz - 1)) / (delthet(-nz) + delthet(nz - 1))
      df(nz) = df(-nz)

   end subroutine get_dthet

   !****************************************************************************
   !                       Subroutine for computing ∂^2 . /∂θ^2
   !****************************************************************************
   ! Given a function f(theta), calculate second derivative of f with respect
   ! to theta. This is second order accurate, with equal grid spacing assumed
   !****************************************************************************
   subroutine get_d2dthet2(f, d2f)

      implicit none

      real, dimension(-nz:), intent(in) :: f
      real, dimension(-nz:), intent(out) :: d2f

      !-------------------------------------------------------------------------
      
      ! Assuming equal grid spacing in theta here
      d2f(-nz + 1:nz - 1) = (f(:nz - 2) - 2.*f(-nz + 1:nz - 1) + f(-nz + 2:)) / delthet(-nz + 1:nz - 1)**2

      ! Use periodicity at boundary
      d2f(-nz) = (f(nz - 1) - 2.*f(-nz) + f(-nz + 1)) / delthet(-nz + 1)**2
      d2f(nz) = d2f(-nz)
      
   end subroutine get_d2dthet2

   !****************************************************************************
   !                               Compute Jacrho 
   !****************************************************************************
   ! Get the Jacobian of the transformation from (r, θ, ζ) -> (R, Z, ζ)
   ! This is called jacr or Jacrho in the comments as opposed
   ! to jacobian, which is for tranformation from (ψ, θ, ζ) -> (R, Z, ζ)
   !****************************************************************************
   subroutine get_jacrho

      implicit none

      !-------------------------------------------------------------------------

      ! <Jacrho> = R * (∂R/∂r * ∂Z/∂θ - ∂R/∂θ * ∂Z/dr)
      jacrho = Rr(2, :) * (dRdrho * dZdth - dRdth * dZdrho)

   end subroutine get_jacrho

   !****************************************************************************
   !                             Compute ∇r and ∇ψ
   !****************************************************************************
   ! Define: 
   ! <grho> = ∇r
   ! <gpsi> = ∇ψ
   !****************************************************************************
   subroutine get_gradrho(dpsipdrho, grho)

      implicit none

      real, intent(in) :: dpsipdrho
      real, dimension(-nz:), intent(out) :: grho

      !-------------------------------------------------------------------------

      ! ∇r = R * sqrt[ (∂Z/∂θ)^2 + (∂R/∂r)^2 ] / Jacrho
      grho = Rr(2, :) * sqrt(dRdth**2 + dZdth**2) / jacrho
      ! ∇ψ = ∇r * ∂ψ/∂r
      gpsi = grho * dpsipdrho

   end subroutine get_gradrho

   !****************************************************************************
   !                                Compute ∂I/∂r
   !****************************************************************************
   subroutine get_dIdrho(dpsipdrho, grho, dIdrho)

      use constants, only: pi

      implicit none

      real, intent(in) :: dpsipdrho
      real, dimension(-nz:), intent(in) :: grho
      real, intent(out) :: dIdrho

      real :: num1, num2, denom
      real, dimension(:), allocatable :: dum

      !-------------------------------------------------------------------------

      allocate (dum(-nz:nz)); dum = 0.

      dum = jacrho * (1.0 + (bi / gpsi)**2) / Rr(2, :)**2
      call theta_integrate(dum(-nz2pi:nz2pi), denom)

      dum = jacrho * (2.*dRdrho / Rr(2, :) + dqdr / local%qinp) / Rr(2, :)**2
      call theta_integrate(dum(-nz2pi:nz2pi), num1)

      ! Here, <betaprim> = β' = 4*π*p/B0^2 * (- ∂lnp/∂r )
      ! where p is the total pressure and B0 is the magnetic field strength
      dum = (-2.*(dRdth * d2Rdrdth + dZdth * d2Zdrdth) / jacrho &
             + drzdth + local%betaprim * jacrho / dpsipdrho**2) / grho**2
      call theta_integrate(dum(-nz2pi:nz2pi), num2)

      dIdrho = bi * (num1 + num2) / denom

      deallocate (dum)

   end subroutine get_dIdrho

   !****************************************************************************
   !                              Compute ∂J/∂r
   !****************************************************************************
   subroutine get_djacdrho(dpsipdrho, dIdrho, grho)

      implicit none

      real, intent(in) :: dpsipdrho, dIdrho
      real, dimension(-nz:), intent(in) :: grho

      !-------------------------------------------------------------------------

      ! This is ∂ψ/∂r * ∂(Jacrho)/∂r
      ! Here, <betaprim> = β' = 4*π*p/B0^2 * (- ∂lnp/∂r )
      ! where p is the total pressure and B0 is the magnetic field strength
      djacdrho = (Rr(2, :) / grho)**2 * (2.*(dRdth * d2Rdrdth + dZdth * d2Zdrdth) / jacrho &
                                         - drzdth + jacrho * (bi * dIdrho / Rr(2, :)**2 - local%betaprim) / dpsipdrho**2)

      ! This is ∂(Jacrho)/∂r
      djacrdrho = djacdrho + jacrho * local%d2psidr2 / dpsipdrho

   end subroutine get_djacdrho

   !****************************************************************************
   !                        Compute ∂^2 R/∂r^2 and ∂^2 Z/∂r^2
   !****************************************************************************
   subroutine get_d2RZdr2

      implicit none

      !-------------------------------------------------------------------------

      ! Get factor common to both d^2 R/∂r^2 and d^2 Z/∂r^2
      d2Rdr2 = ((djacrdrho - jacrho * dRdrho / Rr(2, :)) / Rr(2, :) &
                - dRdrho * d2Zdrdth + dZdrho * d2Rdrdth) / (dRdth**2 + dZdth**2)

      d2Zdr2 = -d2Rdr2 * dRdth
      d2Rdr2 = d2Rdr2 * dZdth

   end subroutine get_d2RZdr2

   !****************************************************************************
   !                     Compute ∂(|∇r|^2)/∂r and ∂(|∇ψ|^2)/∂r
   !****************************************************************************
   subroutine get_dgr2dr(dpsipdrho, grho)

      implicit none

      real, intent(in) :: dpsipdrho
      real, dimension(-nz:), intent(in) :: grho

      !-------------------------------------------------------------------------
      ! ∂(|∇r|^2)/∂r
      dgr2dr = 2.*(grho**2 * (dRdrho / Rr(2, :) - djacrdrho / jacrho) &
                   + (Rr(2, :) / jacrho)**2 * (dRdth * d2Rdrdth + d2Zdrdth * dZdth))
      ! ∂(|∇ψ|^2)/∂r
      dgpsi2dr = 2.*(gpsi**2 * (dRdrho / Rr(2, :) - djacdrho / jacrho) &
                     + (Rr(2, :) / jacrho)**2 * (dRdth * d2Rdrdth + d2Zdrdth * dZdth) * dpsipdrho**2)

   end subroutine get_dgr2dr

   !****************************************************************************
   !                   Compute: ∇θ.∇θ , ∇r.∇θ , ∇α.∇θ , ∇r.∇α , ∇α.∇α
   !****************************************************************************
   ! Define: 
   ! <gradthet2> = ∇θ.∇θ = |∇θ|^2
   ! <gradrho_gradthet> = ∇r.∇θ
   ! <gradalph_gradthet> = ∇α.∇θ
   ! <gradrho_gradalph> = ∇r.∇α
   ! <gradalph2> = ∇α.∇α = |∇α|^2
   !****************************************************************************
   subroutine get_graddotgrad(dpsipdrho, grho)

      implicit none

      real, intent(in) :: dpsipdrho
      real, dimension(-nz:), intent(in) :: grho

      !-------------------------------------------------------------------------

      ! ∇θ.∇θ
      gradthet2 = (Rr(2, :) / jacrho)**2 * (dRdrho**2 + dZdrho**2)

      ! ∇r.∇θ
      gradrho_gradthet = -(Rr(2, :) / jacrho)**2 * (dRdrho * dRdth + dZdrho * dZdth)

      ! ∇α.∇θ
      gradalph_gradthet = -(varthet * dqdr + local%qinp * dvarthdr) * gradrho_gradthet &
                          - bi * jacrho / (dpsipdrho * Rr(2, :)**2) * gradthet2
      
      ! ∇r.∇α
      gradrho_gradalph = -(varthet * dqdr + local%qinp * dvarthdr) * grho**2 &
                         - bi * jacrho / (dpsipdrho * Rr(2, :)**2) * gradrho_gradthet
      
      ! ∇α.∇α
      gradalph2 = (1./Rr(2, :)**2) + ((varthet * dqdr + local%qinp * dvarthdr) * grho)**2 &
                  + 2.*bi * jacrho * (varthet * dqdr + local%qinp * dvarthdr) * gradrho_gradthet / (dpsipdrho * Rr(2, :)**2) &
                  + (bi * jacrho / (dpsipdrho * Rr(2, :)**2))**2 * gradthet2

   end subroutine get_graddotgrad

   !****************************************************************************
   !                           Compute Gradient Factors
   !****************************************************************************
   ! Define:
   ! <grady_dot_grady> = ∇y.∇y = ∇α|^2 * (∂ψ_N/∂r)^2
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

      ! <grady_dot_grady> = ∇y.∇y = ∇α|^2 * (∂ψ_N/∂r)^2
      ! Note: the (∂ψ_N/∂r) factor accounts for ky normalization
      grady_dot_grady = gradalph2 * dpsipdrho_psi0**2

      ! <gradx_dot_grady> = ∇x.∇y = (∇q . ∇α) * (∂ψ_N/∂r)^2
      gradx_dot_grady = gradrho_gradalph * dqdr * dpsipdrho_psi0**2 / shat

      ! <gradx_dot_gradx> = ∇x.∇x = |∇q|^2 * (∂ψ_N/∂r)^2
      gradx_dot_gradx = (grho * dpsipdrho_psi0 * dqdr)**2 / (shat)**2

      ! <gds23> = (∇r . ∇θ * |∇α|^2 - ∇α . ∇θ * ∇r . ∇α) * (∂ψ_N/∂r)^2 / B^2
      gds23 = (gradrho_gradthet * gradalph2 - gradalph_gradthet * gradrho_gradalph) * (dpsipdrho_psi0 / bmag)**2

      ! <gds24> = (∇r . ∇θ * ∇r . ∇α - ∇α . ∇θ * |∇r|^2) * (∂ψ_N/∂r)^2 / B^2 * q/rho
      gds24 = (gradrho_gradthet * gradrho_gradalph - gradalph_gradthet * grho**2) &
              * (dpsipdrho_psi0 / bmag)**2 * (local%qinp_psi0 / local%rhoc_psi0)

      ! Note that kperp2 = (n0/a)^2 * (∂ψ_N/∂r)^2 * [ ∇y.∇y + 2*theta0*(∇y.∇x)*shat + theta0^2*(∇x.∇x)*shat^2 ]
      ! theta0 = kx/(ky*shat)

   end subroutine get_gds

   !****************************************************************************
   !                             Compute ∂B/∂r
   !****************************************************************************
   subroutine get_dBdrho(bmag, dIdrho)

      implicit none

      real, dimension(-nz:), intent(in) :: bmag
      real, intent(in) :: dIdrho

      !-------------------------------------------------------------------------

      ! ∂B/∂r
      dBdrho = (bi * dIdrho + 0.5 * dgpsi2dr) / (bmag * Rr(2, :)**2) &
               - bmag * dRdrho / Rr(2, :)

   end subroutine get_dBdrho

   !****************************************************************************
   !                                Compute ϑ (vartheta)
   !****************************************************************************
   ! Vartheta is defined through α = φ − qϑ, which is used to get: 
   !                          ϑ = I/q int_0^θ dθ' Jacrho / R^2
   ! This is an integral over θ from 0 to θ
   !****************************************************************************
   subroutine get_varthet(dpsipdrho)

      implicit none

      real, intent(in) :: dpsipdrho

      !-------------------------------------------------------------------------

      call theta_integrate_indef(jacrho / Rr(2, :)**2, varthet)
      varthet = bi * varthet / (dpsipdrho * local%qinp)

   end subroutine get_varthet

   !****************************************************************************
   !                                   Compute ∂ϑ/∂r
   !****************************************************************************
   ! ∂ϑ/∂r= (I'/I − q'/q) * ϑ + I/q int_0^θ dθ' ∂(Jacrho / R^2)/∂r
   !****************************************************************************
   subroutine get_dvarthdr(dpsipdrho, dIdrho)

      implicit none

      real, intent(in) :: dpsipdrho, dIdrho

      real, dimension(-nz:nz) :: dum

      !-------------------------------------------------------------------------

      dum = bi * jacrho * (dIdrho / bi - dqdr / local%qinp + djacdrho / jacrho &
                           - 2.*dRdrho / Rr(2, :)) / Rr(2, :)**2
      call theta_integrate_indef(dum, dvarthdr)

      ! ∂ϑ/∂r
      dvarthdr = dvarthdr / (dpsipdrho * local%qinp)

   end subroutine get_dvarthdr

   !****************************************************************************
   !                           Compute ∂^2I/∂r^2, ∂^2Jacrho/∂r^2
   !****************************************************************************
   subroutine get_d2Idr2_d2jacdr2(grho, dIdrho)

      use constants, only: pi

      implicit none

      real, dimension(-nz:), intent(in) :: grho
      real, intent(in) :: dIdrho

      real :: denom, num1, num2, num3, num4
      real, dimension(-nz:nz) :: tmp, tmp2

      !-------------------------------------------------------------------------

      ! Denom is the denominator in the expression for ∂^2 I /∂r^2
      tmp = jacrho / Rr(2, :)**2 * (1.0 + (bi / gpsi)**2)
      call theta_integrate(tmp(-nz2pi:nz2pi), denom)
      denom = denom / bi

      d2jacdr2 = dIdrho * bi * jacrho / gpsi**2 &
                 * (dIdrho / bi + djacrdrho / jacrho - dgpsi2dr / gpsi**2 &
                    - 2.*dRdrho / Rr(2, :))

      tmp = -d2jacdr2 / Rr(2, :)**2 - dIdrho * jacrho / (bi * Rr(2, :)**2) &
            * (djacrdrho / jacrho - dIdrho / bi - 2.*dRdrho / Rr(2, :))
      call theta_integrate(tmp(-nz2pi:nz2pi), num1)

      tmp = (d2Rdr2 * dRdth + dRdrho * d2Rdrdth + d2Zdr2 * dZdth + dZdrho * d2Zdrdth) / jacrho &
            - djacrdrho * (dRdrho * dRdth + dZdrho * dZdth) / jacrho**2
      call get_dthet(tmp, tmp2)
      tmp = (tmp2 - 2./jacrho * (-djacrdrho / jacrho * (dRdth * d2Rdrdth + dZdth * d2Zdrdth) &
                                 + d2Rdrdth**2 + dRdth * d2Rdr2dth + d2Zdrdth**2 + dZdth * d2Zdr2dth)) / grho**2 &
            - dgr2dr * (drzdth - 2./jacrho * (dRdth * d2Rdrdth + dZdth * d2Zdrdth)) / grho**4
      call theta_integrate(tmp(-nz2pi:nz2pi), num2)
      d2jacdr2 = d2jacdr2 - tmp * Rr(2, :)**2

      tmp = jacrho * (local%betadbprim + local%betaprim * (djacrdrho / jacrho - dgpsi2dr / gpsi**2)) / gpsi**2
      call theta_integrate(tmp(-nz2pi:nz2pi), num3)
      
      !FLAG - next negative sign?
      d2jacdr2 = d2jacdr2 - tmp * Rr(2, :)**2

      tmp = jacrho / Rr(2, :)**2 * (2.*d2Rdr2 / Rr(2, :) - 2.*(dRdrho / Rr(2, :))**2 &
                                    + local%d2qdr2 / local%qinp - (dqdr / local%qinp)**2 + (2 * dRdrho / Rr(2, :) + dqdr / local%qinp) &
                                    * (djacrdrho / jacrho - 2.*dRdrho / Rr(2, :)))
      call theta_integrate(tmp(-nz2pi:nz2pi), num4)

      d2Idr2 = (num1 + num2 + num3 + num4) / denom

      d2jacdr2 = d2jacdr2 + bi * jacrho / gpsi**2 * d2Idr2 + 2.*djacdrho * dRdrho / Rr(2, :)

   end subroutine get_d2Idr2_d2jacdr2

   !****************************************************************************
   !                              Compute ∂^2ϑ/∂r^2
   !****************************************************************************
   subroutine get_d2varthdr2(dpsipdrho, dIdrho)

      implicit none

      real, intent(in) :: dpsipdrho, dIdrho

      real, dimension(-nz:nz) :: dum

      !-------------------------------------------------------------------------

      dum = bi * jacrho / (local%qinp * dpsipdrho * Rr(2, :)**2) * ((dIdrho / bi - dqdr / local%qinp &
            + djacdrho / jacrho - 2.*dRdrho / Rr(2, :))**2 &
            + d2Idr2 / bi - (dIdrho / bi)**2 - local%d2qdr2 / local%qinp &
            + (dqdr / local%qinp)**2 + d2jacdr2 / jacrho - (djacdrho / jacrho)**2 &
            - djacdrho * local%d2psidr2 / (dpsipdrho * jacrho) &
            - 2.*d2Rdr2 / Rr(2, :) + 2.*(dRdrho / Rr(2, :))**2)

      call theta_integrate_indef(dum, d2varthdr2)

   end subroutine get_d2varthdr2

   !****************************************************************************
   !                             Compute ∂^2B/∂r^2
   !****************************************************************************
   subroutine get_d2Bdr2(bmag, dIdrho)

      implicit none

      real, dimension(-nz:), intent(in) :: bmag
      real, intent(in) :: dIdrho

      !-------------------------------------------------------------------------

      d2gpsidr2 = 2.*(dRdrho / Rr(2, :) - djacdrho / jacrho) * dgpsi2dr &
                  + 2.*gpsi**2 * (d2Rdr2 / Rr(2, :) - (dRdrho / Rr(2, :))**2 - d2jacdr2 / jacrho + djacdrho * djacrdrho / jacrho**2) &
                  + 2.*(Rr(2, :) * gpsi / jacrho)**2 * (d2Rdrdth**2 + dRdth * d2Rdr2dth + d2Zdrdth**2 + dZdth * d2Zdr2dth &
                  + 2.*(dRdth * d2Rdrdth + dZdth * d2Zdrdth) * (dRdrho / Rr(2, :) - djacdrho / jacrho))

      ! Get ∂/∂r (∂B/∂r)
      d2Bdr2 = -dBdrho * dRdrho / Rr(2, :) + bmag * (dRdrho / Rr(2, :))**2 &
               - bmag * d2Rdr2 / Rr(2, :) + 0.5 * (2.*(dIdrho**2 + bi * d2Idr2) &
                                                   + d2gpsidr2) / (bmag * Rr(2, :)**2) &
               - (dBdrho + bmag * dRdrho / Rr(2, :)) * (2.*dRdrho / Rr(2, :) + dBdrho / bmag)

   end subroutine get_d2Bdr2

   !****************************************************************************
   !                     Get More Derivatives with respect to r
   !****************************************************************************
   ! Things computed in this routine are: 
   ! <dgrgt> = ∂(∇r . ∇θ)/∂r
   ! <dgt2> = ∂(|∇θ|^2)/∂r
   ! <dga2> = ∂(|∇α|^2)/∂r
   ! <dgagr> = ∂(∇α . ∇r)/∂r
   ! <dgagt> = ∂(∇α . ∇θ)/∂r
   ! <dcrossdr> = ∂[(∇α x B) . ∇θ)]/∂r
   ! <dgds2dr> = (∂ψ/∂r)^2 * ∂(|∇α|^2)/∂r
   ! <dgds21dr> = (∂ψ/∂r)^2 * ∂(∇α . ∇q)/∂r
   ! <dgds22dr> = (∂ψ/∂r)^2 * ∂(|∇q|^2)/∂r
   !****************************************************************************
   subroutine get_dcrossdr(dpsipdrho, dIdrho, grho)

      implicit none

      real, intent(in) :: dpsipdrho, dIdrho
      real, dimension(-nz:), intent(in) :: grho

      !-------------------------------------------------------------------------

      ! dgr2 = ∂/∂r (|∇r|^2)
      ! dgr2 = 2.*(Rr(2,:)/jacrho)**2*((dRdrho/Rr(2,:)-djacdrho/jacrho)*(dRdth**2+dZdth**2) &
      !      + dRdth*d2Rdrdth + dZdth*d2Zdrdth)

      ! <dgrgt> = ∂(∇r . ∇θ)/∂r
      dgrgt = 2.*gradrho_gradthet * (dRdrho / Rr(2, :) - djacrdrho / jacrho) &
              - (Rr(2, :) / jacrho)**2 * (d2Rdr2 * dRdth + dRdrho * d2Rdrdth + d2Zdr2 * dZdth + dZdrho * d2Zdrdth)

      ! <dgt2> = ∂(|∇θ|^2)/∂r
      dgt2 = 2.*(Rr(2, :) / jacrho)**2 * ((dRdrho / Rr(2, :) - djacrdrho / jacrho) * (dRdrho**2 + dZdrho**2) &
                                          + dRdrho * d2Rdr2 + dZdrho * d2Zdr2)
      ! <dga2> = ∂(|∇α|^2)/∂r
      ! will later multiply it by:   0.5 * (∂ψ_N/∂r)^2
      dga2 = -2 * dRdrho / Rr(2, :)**3 + dgr2dr * (varthet * dqdr + local%qinp * dvarthdr)**2 &
             + (2.0 * grho**2 * (varthet * dqdr + local%qinp * dvarthdr) &
                + 2.*bi * jacrho * gradrho_gradthet / (dpsipdrho * Rr(2, :)**2)) &
             * (local%d2qdr2 * varthet + 2.*dqdr * dvarthdr + local%qinp * d2varthdr2) &
             + 2.*(varthet * dqdr + local%qinp * dvarthdr) * bi * jacrho / (dpsipdrho * Rr(2, :)**2) &
             * (dgrgt + gradrho_gradthet * (dIdrho / bi + djacdrho / jacrho - 2.*dRdrho / Rr(2, :))) &
             + (bi * jacrho / (dpsipdrho * Rr(2, :)**2))**2 * (dgt2 + 2.*gradthet2 * (dIdrho / bi + djacdrho / jacrho &
                                                                                     - 2.*dRdrho / Rr(2, :)))

      ! <dgagr> = ∂(∇α . ∇r)/∂r
      dgagr = -grho**2 * (2.*dvarthdr * dqdr + varthet * local%d2qdr2 + local%qinp * d2varthdr2) &
              - dgr2dr * (varthet * dqdr + local%qinp * dvarthdr) - bi * jacrho / (dpsipdrho * Rr(2, :)**2) &
              * (dgrgt + gradrho_gradthet * (dIdrho / bi + djacdrho / jacrho - 2.*dRdrho / Rr(2, :)))

      ! <dgagt> = ∂(∇α . ∇θ)/∂r
      dgagt = -gradrho_gradthet * (2.*dvarthdr * dqdr + varthet * local%d2qdr2 + local%qinp * d2varthdr2) &
              - dgrgt * (varthet * dqdr + local%qinp * dvarthdr) - bi * jacrho / (dpsipdrho * Rr(2, :)**2) &
              * (dgt2 + gradthet2 * (dIdrho / bi + djacdrho / jacrho - 2.*dRdrho / Rr(2, :)))

      ! <dcrossdr> = ∂[(∇α x B) . ∇θ)]/∂r
      dcrossdr = dpsipdrho * (dgagr * gradalph_gradthet + gradrho_gradalph * dgagt &
                             - dga2 * gradrho_gradthet - gradalph2 * dgrgt) + local%d2psidr2 * cross / dpsipdrho

      ! <dgds2dr> = (∂ψ/∂r)^2 * ∂(|∇α|^2)/∂r
      dgds2dr = dga2 * dpsipdrho_psi0**2

      ! <dgds21dr> = (∂ψ/∂r)^2 * ∂(∇α . ∇q)/∂r
      ! Note that there will be multiplication by 2 in dist_fn.fpp
      dgds21dr = (dgagr * dqdr + local%d2qdr2 * gradrho_gradalph) * dpsipdrho_psi0**2

      ! <dgds22dr> = (∂ψ/∂r)^2 * ∂(|∇q|^2)/∂r
      dgds22dr = (dqdr**2 * dgr2dr + 2.*grho**2 * dqdr * local%d2qdr2) * dpsipdrho_psi0**2

      ! note that dkperp2/∂r = (n0/a)^2*(∂ψ_N/∂r)^2*(dgds2dr + 2*theta0*dgds21dr + theta0^2*dgds22dr)

   end subroutine get_dcrossdr

   !****************************************************************************
   !                       Definite Integral in θ - from 0 to 2*π
   !****************************************************************************
   subroutine theta_integrate(integrand, integral)

      implicit none

      real, dimension(-nz2pi:), intent(in) :: integrand
      real, intent(out) :: integral

      !-------------------------------------------------------------------------

      ! Use trapezoidal rule to integrate in theta
      integral = 0.5 * sum(delthet(-nz2pi:nz2pi - 1) * (integrand(-nz2pi:nz2pi - 1) + integrand(-nz2pi + 1:nz2pi)))

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
         integral(i) = integral(i - 1) + 0.5 * delthet(i - 1) * (integrand(i - 1) + integrand(i))
      end do
      do i = -1, -nz, -1
         integral(i) = integral(i + 1) - 0.5 * delthet(i) * (integrand(i + 1) + integrand(i))
      end do
      
   end subroutine theta_integrate_indef

   !****************************************************************************
   !                            COMMUNICATE PARAMETERS 
   !****************************************************************************
   ! Only needed for radial_variation simulations
   !****************************************************************************
   subroutine communicate_parameters_multibox(surf, drl, drr)
      use mp, only: job, scope, mp_abort, &
                    crossdomprocs, subprocs, &
                    send, receive
      use job_manage, only: njobs
      use common_types, only: flux_surface_type

      implicit none

      real, optional, intent(in) :: drl, drr
      type(flux_surface_type), intent(inout) :: surf

      real :: lrhoc, lqinp, lshat, lkappa, ltri, lbetaprim
      real :: rrhoc, rqinp, rshat, rkappa, rtri, rbetaprim
      real :: dqdr
      real :: rhoc_psi0, qinp_psi0, shat_psi0

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
         call send(dIdrho, 0, 129)
         call send(rhoc, 0, 130)
         call send(qinp, 0, 131)
         call send(shat, 0, 132)
         call send(dpsipdrho, 0, 133)
         call send(bmag, 0, 134)
         call send(grho, 0, 135)

         call send(rrhoc, njobs - 1, 220)
         call send(rqinp, njobs - 1, 221)
         call send(rshat, njobs - 1, 222)
         call send(rkappa, njobs - 1, 223)
         call send(rtri, njobs - 1, 224)
         call send(rbetaprim, njobs - 1, 225)
         call send(local%rhoc, njobs - 1, 226)
         call send(d2R, njobs - 1, 227)
         call send(d2Z, njobs - 1, 228)
         call send(dIdrho, njobs - 1, 229)
         call send(rhoc, njobs - 1, 230)
         call send(qinp, njobs - 1, 231)
         call send(shat, njobs - 1, 232)
         call send(dpsipdrho, njobs - 1, 233)
         call send(bmag, njobs - 1, 234)
         call send(grho, njobs - 1, 235)
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
         call receive(dpsipdrho_psi0, 1, 133)
         call receive(bmag_psi0, 1, 134)
         call receive(grho_psi0, 1, 135)
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
         call receive(dpsipdrho_psi0, 1, 233)
         call receive(bmag_psi0, 1, 234)
         call receive(grho_psi0, 1, 235)
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
   !              R(r,θ) = R0(r) + rcos[θ+ sin (θarcsin \bar{δ}(r)) ]
   !****************************************************************************
   function Rpos(r, theta, j)

      use constants, only: pi

      integer, intent(in) :: j
      real, intent(in) :: r, theta
      real :: Rpos
      real :: g, gp, dr
      integer :: i

      !-------------------------------------------------------------------------

      dr = r - local%rhoc

      ! For Y Xiao:
      !    g = local%delp/local%rhoc + local%d * sin(theta)**2
      !    Rpos = local%rmaj*(1.+r*(cos(theta)-g)-g*dr)

      g = cos(theta + local%tri * sin(theta))
      gp = -sin(theta + local%tri * sin(theta)) &
           * local%triprim * sin(theta)

      ! Allow for strange specification of R_ψ
      if (j == nz + 1) then
         i = -nz
      else
         i = j
      end if

      ! Second line here is (1/2)*(r-r0)^2 * ∂^2 R/ ∂r|_r0
      ! Note that d2R=0 unless read_profile_variation = T in input file
      Rpos = local%rmaj + local%shift * dr + g * local%rhoc + (g + local%rhoc * gp) * dr &
             + 0.5 * (r - rhoc0)**2 * d2R(i)

   end function Rpos

   !****************************************************************************
   !                                Function for Z 
   !****************************************************************************
   !                              Z(r,θ) = κ(r)rsin θ
   !****************************************************************************
   function Zpos(r, theta, j)

      integer, intent(in) :: j
      real, intent(in) :: r, theta
      real :: Zpos, dr
      integer :: i

      !-------------------------------------------------------------------------

      ! Allow for strange specification of Z_psi
      if (j == nz + 1) then
         i = -nz
      else
         i = j
      end if

      dr = r - local%rhoc
      ! Note that d2Z=0 unless read_profile_variation = T in input file
      Zpos = local%kappa * sin(theta) * local%rhoc + (local%rhoc * local%kapprim + local%kappa) * sin(theta) * dr &
             + 0.5 * (r - rhoc0)**2 * d2Z(i)

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
