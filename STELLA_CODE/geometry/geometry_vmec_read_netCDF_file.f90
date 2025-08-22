!###############################################################################
!################################ VMEC GEOMETRY ################################
!###############################################################################
! 
! Routines for calculating the geometry needed by stella, from a VMEC file.
! Inside the <geometry> module we call:
! 
! call get_vmec_geometry(&
!            nzgrid, nalpha, naky, geo_surf, grho, bmag, gradpar, &
!            b_dot_grad_z, grad_alpha_grad_alpha, &
!            grad_alpha_grad_psit, grad_psit_grad_psitt, &
!            gds23, gds24, gds25, gds26, gbdrift_alpha, gbdrift0_psit, &
!            cvdrift_alpha, cvdrift0_psit, sign_torflux, &
!            theta, dzetadz, aref, bref, alpha, zeta, &
!            field_period_ratio, psit_displacement_fac)
! 
! The VMEC module will calculate the geometric arrays with psi_t as the
! the radial coordinate, and zeta as the parallel coordinate, on a z-grid 
! with <z_grid_refinement_factor> more z-points than the real stella grid. 
! This module can change the parallel coordinate to the normalized arc-length,
! and it interpolates the VMEC z-grid to the stella z-grid.
! 
! Note that I removed <b_dot_grad_zeta_prefac> and <z_grid_scalefac>.
!  
! Initial VMEC geometry code was written by Matt Landreman, University of Maryland in August 2017. 
! Modified in 2018-2019 by Michael Barnes, and cleaned in 2024 by Hanne Thienpondt.
! 
!###############################################################################

module geometry_vmec_read_netCDF_file

   implicit none

   public :: calculate_vmec_geometry

   private

   logical, parameter :: debug = .false.

   logical :: lasym
   integer :: nfp, isigng
   integer :: ns, mnmax, mnmax_nyq
   integer :: mpol, ntor
   real :: Aminor
   real, dimension(:), allocatable :: xm, xn
   real, dimension(:), allocatable :: xm_nyq, xn_nyq
   real, dimension(:, :), allocatable :: rmnc, rmns
   real, dimension(:, :), allocatable :: lmnc, lmns
   real, dimension(:, :), allocatable :: zmnc, zmns
   real, dimension(:, :), allocatable :: bmnc, bmns
   real, dimension(:, :), allocatable :: gmnc, gmns
   real, dimension(:, :), allocatable :: bsupumnc, bsupumns
   real, dimension(:, :), allocatable :: bsupvmnc, bsupvmns
   real, dimension(:, :), allocatable :: bsubumnc, bsubumns
   real, dimension(:, :), allocatable :: bsubvmnc, bsubvmns
   real, dimension(:, :), allocatable :: bsubsmnc, bsubsmns
   real, dimension(:), allocatable :: phi, phip, iotas, iotaf, presf

   ! Radial integration weights and surfaces
   real, dimension(2) :: weight_full, weight_half
   integer, dimension(2) :: index_full, index_half

contains

   !============================================================================
   !=========================== READ VMEC EQUILIBRIUM ==========================
   !============================================================================ 
   ! Read in everything from the vmec wout file using libstell.
   ! <vmec_filename> is the vmec wout_* file that will be read.
   !============================================================================ 
   subroutine read_vmec_equilibrium(vmec_filename, verbose, ierr)

      ! Rename all variables from the VMEC file as '*_vmec' 
      use read_wout_mod, only: read_wout_file, read_wout_deallocate
      use read_wout_mod, only: nfp_vmec => nfp
      use read_wout_mod, only: lasym_vmec => lasym
      use read_wout_mod, only: isigng_vmec => isigng
      use read_wout_mod, only: Aminor_vmec => Aminor
      use read_wout_mod, only: ns_vmec => ns
      use read_wout_mod, only: mnmax_nyq_vmec => mnmax_nyq
      use read_wout_mod, only: mnmax_vmec => mnmax
      use read_wout_mod, only: mpol_vmec => mpol
      use read_wout_mod, only: ntor_vmec => ntor
      use read_wout_mod, only: xm_vmec => xm
      use read_wout_mod, only: xn_vmec => xn
      use read_wout_mod, only: xm_nyq_vmec => xm_nyq
      use read_wout_mod, only: xn_nyq_vmec => xn_nyq
      use read_wout_mod, only: phi_vmec => phi
      use read_wout_mod, only: phip_vmec => phip
      use read_wout_mod, only: lmnc_vmec => lmnc
      use read_wout_mod, only: lmns_vmec => lmns
      use read_wout_mod, only: rmnc_vmec => rmnc
      use read_wout_mod, only: rmns_vmec => rmns
      use read_wout_mod, only: zmnc_vmec => zmnc
      use read_wout_mod, only: zmns_vmec => zmns
      use read_wout_mod, only: bmnc_vmec => bmnc
      use read_wout_mod, only: bmns_vmec => bmns
      use read_wout_mod, only: gmnc_vmec => gmnc
      use read_wout_mod, only: gmns_vmec => gmns
      use read_wout_mod, only: bsupumnc_vmec => bsupumnc
      use read_wout_mod, only: bsupvmnc_vmec => bsupvmnc
      use read_wout_mod, only: bsupumns_vmec => bsupumns
      use read_wout_mod, only: bsupvmns_vmec => bsupvmns
      use read_wout_mod, only: bsubumnc_vmec => bsubumnc
      use read_wout_mod, only: bsubvmnc_vmec => bsubvmnc
      use read_wout_mod, only: bsubumns_vmec => bsubumns
      use read_wout_mod, only: bsubvmns_vmec => bsubvmns
      use read_wout_mod, only: bsubsmnc_vmec => bsubsmnc
      use read_wout_mod, only: bsubsmns_vmec => bsubsmns
      use read_wout_mod, only: iotas_vmec => iotas
      use read_wout_mod, only: iotaf_vmec => iotaf
      use read_wout_mod, only: presf_vmec => presf
      use mp, only: mp_abort

      implicit none
      
      ! TODO-GA: import this from run_parameters as "use run_parameters, only: print_extra_info_to_terminal"
      logical :: print_extra_info_to_terminal = .false.

      ! vmec_filename is the vmec wout_* file that will be read.
      logical, intent(in) :: verbose
      character(*), intent(in) :: vmec_filename

      integer, intent(out) :: ierr
      integer :: iopen

      !---------------------------------------------------------------------- 

      ! Track the code 
      if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::read_vmec_equilibrium'
 
      ! Print some key information about the VMEC file to the command prompt
      if (verbose .and. print_extra_info_to_terminal) then
         write (*, '(A)') "############################################################"
         write (*, '(A)') "                       MAGNETIC FIELD"
         write (*, '(A)') "############################################################"
         write (*, *) ' '; write (*, *) "The VMEC data is read from: '", trim(vmec_filename), "'."
      end if

      ! Read the VMEC file
      call read_wout_file(vmec_filename, ierr, iopen)

      ! Return errors if the reading of the VMEC fails
      if (iopen /= 0) then
         write(*,*) ' '; write(*,*) 'In the input file, you specified <vmec_filename> = ', trim(vmec_filename), '.'
         write(*,*) 'This file is missing in this folder. Aborting.'; write(*,*) ' ' 
         ierr = iopen
         call mp_abort('read_vmec_equilibrium returned error.') 
      end if

      if (ierr /= 0) then
         write(*,*) 'Error reading wout file'
         call mp_abort('read_vmec_equilibrium returned error.') 
      end if

      ! Make the following VMEC variables available to the entire module
      nfp = nfp_vmec
      lasym = lasym_vmec
      isigng = isigng_vmec
      aminor = aminor_vmec
      ns = ns_vmec
      mnmax = mnmax_vmec
      mnmax_nyq = mnmax_nyq_vmec
      mpol = mpol_vmec
      ntor = ntor_vmec

      ! Print some important characteristics of the VMEC to the command prompt
      if (verbose .and. print_extra_info_to_terminal) then
         write (*, *) " "
         write (*, *) "  Characteristics of the magnetic field:"
         write (*, '(A52, I2)') "      Number of field periods of the machine (nfp): ", nfp
         write (*, '(A40, L1)') "      Stellarator-asymmetric? (lasym):  ", lasym
      end if

      ! Allocate and assign the VMEC array, so they are available to the entire module
      if (.not. allocated(rmnc)) then

         ! Basic VMEC arrays and symmetric data
         allocate (xm(mnmax)); xm = xm_vmec
         allocate (xn(mnmax)); xn = xn_vmec
         allocate (xm_nyq(mnmax_nyq)); xm_nyq = xm_nyq_vmec
         allocate (xn_nyq(mnmax_nyq)); xn_nyq = xn_nyq_vmec
         allocate (rmnc(mnmax, ns)); rmnc = rmnc_vmec
         allocate (lmns(mnmax, ns)); lmns = lmns_vmec
         allocate (zmns(mnmax, ns)); zmns = zmns_vmec
         allocate (bmnc(mnmax_nyq, ns)); bmnc = bmnc_vmec
         allocate (gmnc(mnmax_nyq, ns)); gmnc = gmnc_vmec
         allocate (bsupumnc(mnmax_nyq, ns)); bsupumnc = bsupumnc_vmec
         allocate (bsupvmnc(mnmax_nyq, ns)); bsupvmnc = bsupvmnc_vmec
         allocate (bsubumnc(mnmax_nyq, ns)); bsubumnc = bsubumnc_vmec
         allocate (bsubvmnc(mnmax_nyq, ns)); bsubvmnc = bsubvmnc_vmec
         allocate (bsubsmns(mnmax_nyq, ns)); bsubsmns = bsubsmns_vmec
         allocate (phi(ns)); phi = phi_vmec
         allocate (phip(ns)); phip = phip_vmec
         allocate (iotas(ns)); iotas = iotas_vmec
         allocate (iotaf(ns)); iotaf = iotaf_vmec
         allocate (presf(ns)); presf = presf_vmec

         ! If the VMEC is not symmetric we also need the following terms
         if (lasym) then
            allocate (rmns(mnmax, ns)); rmns = rmns_vmec
            allocate (lmnc(mnmax, ns)); lmnc = lmnc_vmec
            allocate (zmnc(mnmax, ns)); zmnc = zmnc_vmec
            allocate (bmns(mnmax_nyq, ns)); bmns = bmns_vmec
            allocate (gmns(mnmax_nyq, ns)); gmns = gmns_vmec
            allocate (bsupumns(mnmax_nyq, ns)); bsupumns = bsupumns_vmec
            allocate (bsupvmns(mnmax_nyq, ns)); bsupvmns = bsupvmns_vmec
            allocate (bsubumns(mnmax_nyq, ns)); bsubumns = bsubumns_vmec
            allocate (bsubvmns(mnmax_nyq, ns)); bsubvmns = bsubvmns_vmec
            allocate (bsubsmnc(mnmax_nyq, ns)); bsubsmnc = bsubsmnc_vmec
         end if

      end if

      ! Deallocate all arrays opened externally in <read_wout_mod>
      ! since they have now been allocated locally inside this module
      call read_wout_deallocate(ierr)
      if (ierr /= 0) then
         write(*,*) "Warning: error returned when deallocating wout arrays. ", ierr
      end if
      ierr = 0

   end subroutine read_vmec_equilibrium

 
   !============================================================================
   !=============== CALCULATE GEOMETRIC ARRAYS NEEDED FOR STELLA ===============
   !============================================================================ 
   subroutine calculate_vmec_geometry(&
                  ! Input parameters
                  vmec_filename, nalpha, alpha0, nzgrid, zeta_center, rectangular_cross_section, & 
                  number_of_field_periods_inputfile, s_inputfile, vmec_surface_option, verbose, & 
                  ! Output parameters
                  s, safety_factor_q, shat, L_reference, B_reference, nfp_out, & 
                  sign_toroidal_flux, alpha, zeta, bmag, b_dot_grad_zeta, grad_alpha_grad_alpha, &
                  grad_alpha_grad_psit, grad_psit_grad_psit, gds23_psitalpha, gds24_psitalpha, & 
                  gds25_psitalpha, gds26_psitalpha, gbdrift_alpha, gbdrift0_psit, cvdrift_alpha, & 
                  cvdrift0_psit,theta, B_sub_zeta, B_sub_theta, psit_displacement_fac, &
                  gradzeta_gradpsit_R2overB2, gradzeta_gradalpha_R2overB2, &
                  b_dot_grad_zeta_RR, ierr)

      use mp, only: mp_abort

      implicit none
 
      !*************************************************************************
      !                            Input parameters                            !
      !*************************************************************************
      ! <nalpha> is the number of grid points in the alpha coordinate (= number of field lines) 
      ! <alpha0> is the first alpha value to include in the alpha grid
      ! The zeta grid has <nzgrid>*2+1 points, including the "repeated" point at index -nzgrid and +nzgrid.
      ! The zeta domain is centered at <zeta_center>. Setting <zeta_center> = 2*pi*N/<nfp> for any integer N should
      ! yield identical results to setting <zeta_center> = 0, where <nfp> is the number of field periods (as in VMEC).
      ! <s_inputfile> determines which flux surface from the VMEC file will be used
      ! for the computation. This parameter should lie in the interval [0,1].
      ! If <verbose> is .true. in the vmec_parameters namelist, lots of diagnostic information is printed.
      ! 
      ! If <number_of_field_periods_inputfile> > 0, then this parameter does what you think:
      ! the extent of the toroidal in zeta will be 2*pi*<number_of_field_periods_inputfile>/<nfp>.
      ! If <number_of_field_periods_inputfile> <= 0, the entire 2*pi toroidal domain will be included.
      ! 
      ! If <vmec_surface_option> = 0, the magnetic surface specified by <s_inputfile> will 
      ! be used, by interpolating between the surfaces available in the vmec file.
      ! If <vmec_surface_option> = 1, the magnetic surface on vmec's HALF radial mesh will be used that is closest to <s_inputfile>.
      ! If <vmec_surface_option> = 2, the magnetic surface on vmec's FULL radial mesh will be used that is closest to <s_inputfile>.
      ! Other values of <vmec_surface_option> will cause the program to abort with an error.

      logical, intent(in) :: rectangular_cross_section, verbose 
      integer, intent(in) :: nalpha, nzgrid, vmec_surface_option
      real, intent(in) :: alpha0, zeta_center, s_inputfile 
      real, intent(in) :: number_of_field_periods_inputfile 
      character(*), intent(in) :: vmec_filename

      !*************************************************************************
      !                            Output quantities                           !
      !************************************************************************* 
      ! On exit, <s> = psi_toroidal / psi_{toroidal,edge} holds the 
      ! flux surface that was actually used for the geometry, 
      !  
      ! The reference length (in meters) and reference magnetic field (in Tesla) are 
      !     <L_reference> = a = Aminor = Aminor_p =minor radius calculated by VMEC
      !     <B_reference> = Bref = 2 * abs(edge_toroidal_flux_over_2pi) / (L_reference * L_reference)
      ! 
      ! VMEC uses the radial coordinate <s> which is related to the coordinates <rho>, <r> and <psi_t>, 
      !     s = psi_toroidal / psi_{toroidal,edge}
      !     rho = sqrt(psi_toroidal / psi_{toroidal,edge})  
      !     r = a * sqrt(s) = a * rho 
      ! 
      ! Other VMEC quantities of interest
      !     <nfp> is the number of field periods of the device (e.g. 5 in W7-X)
      !     <iota> = rotational transform
      !     <safety_factor_q> = 1/iota
      !     <shat> = magnetic shear = (r/q)(dq/dr) = -(r/iota)(diota/dr) = -2(s/iota)(diota/ds)
      ! 
      ! Coordinate grids
      !     <theta_p> = PEST toroidal angle
      !     <zeta> = ζ = toroidal angle zeta
      !     <alpha> = theta_p - iota * zeta
      !  
      ! Geometric arrays
      !     <b_dot_grad_zeta> = b . ∇ζ  (increases in the counter-clockwise direction)
      !     <grad_alpha_grad_alpha> = a^2 ∇α . ∇α  
      !     <grad_alpha_grad_psit> = (1/Bref) ∇α . ∇ψt  
      !     <grad_psit_grad_psit> = 1/(a^2 Bref^2) ∇ψt . ∇ψt  
      !     <gbdrift_alpha> = 2*a^2*Bref/B^3 * B x ∇B . ∇α 
      !     <cvdrift_alpha> = 2*a*Bref/B^2 * B x kappa . ∇α 
      !     <kappa> = (bhat . ∇bhat) 

      integer, intent(out) :: sign_toroidal_flux, ierr  
      real, intent(out) :: L_reference, B_reference, nfp_out
      real, intent(out) :: s, safety_factor_q, shat
      real, dimension(:, -nzgrid:), intent(out) :: theta, bmag, b_dot_grad_zeta, psit_displacement_fac
      real, dimension(:, -nzgrid:), intent(out) :: grad_alpha_grad_alpha, grad_alpha_grad_psit, grad_psit_grad_psit
      real, dimension(:, -nzgrid:), intent(out) :: gds23_psitalpha, gds24_psitalpha, gds25_psitalpha, gds26_psitalpha
      real, dimension(:, -nzgrid:), intent(out) :: gbdrift_alpha, cvdrift_alpha, gbdrift0_psit, cvdrift0_psit, B_sub_theta, B_sub_zeta 
      real, dimension(:, -nzgrid:), intent(out) :: gradzeta_gradpsit_R2overB2, gradzeta_gradalpha_R2overB2 
      real, dimension(:, -nzgrid:), intent(out) :: b_dot_grad_zeta_RR
      real, dimension(-nzgrid:), intent(out) :: zeta
      real, dimension(:), intent(out) :: alpha 

      !*************************************************************************
      !                           Internal variables                           !
      !*************************************************************************
      !-------------------------------------------------------------------------
      !*********************************************************************
      ! VMEC variables of interest:
      ! ns = number of flux surfaces used by VMEC
      ! nfp = number of field periods, e.g. 5 for W7-X, 4 for HSX
      ! iotas = rotational transform (1/q) on the half grid.
      ! iotaf = rotational transform on the full grid.
      ! presf = pressure on the full grid.
      ! 
      ! All VMEC quantities (B, pressure, etc) are in SI units.
      ! 
      ! In VMEC, quantities on the half grid have the same number of array elements (ns) as quantities on the full grid,
      ! but the first array element is 0.
      ! 
      !*********************************************************************  

      real, parameter :: pi = 3.1415926535897932d+0
      real, parameter :: zero = 0.0d+0
      real, parameter :: one = 1.0d+0
      real, parameter :: mu_0 = 4 * pi * (1.0d-7)

      integer :: j, ialpha
      real :: temp, edge_toroidal_flux_over_2pi
      real :: number_of_field_periods, sqrt_s, a, Bref
      real :: iota, d_pressure_d_s, d_iota_d_s  
      real, dimension(:, :), allocatable :: R, Z, B, sqrt_g 
      real, dimension(:, :), allocatable :: d_B_d_theta, d_B_d_zeta, d_B_d_s
      real, dimension(:, :), allocatable :: d_R_d_theta, d_R_d_zeta, d_R_d_s
      real, dimension(:, :), allocatable :: d_Z_d_theta, d_Z_d_zeta, d_Z_d_s
      real, dimension(:, :), allocatable :: d_X_d_theta, d_X_d_zeta, d_X_d_s
      real, dimension(:, :), allocatable :: d_Y_d_theta, d_Y_d_zeta, d_Y_d_s
      real, dimension(:, :), allocatable :: d_Lambda_d_theta, d_Lambda_d_zeta, d_Lambda_d_s
      real, dimension(:, :), allocatable :: B_sub_s, B_sup_theta, B_sup_zeta
      real, dimension(:, :), allocatable :: grad_s_X, grad_s_Y, grad_s_Z
      real, dimension(:, :), allocatable :: grad_theta_X, grad_theta_Y, grad_theta_Z
      real, dimension(:, :), allocatable :: grad_theta_pest_X, grad_theta_pest_Y, grad_theta_pest_Z
      real, dimension(:, :), allocatable :: grad_zeta_X, grad_zeta_Y, grad_zeta_Z
      real, dimension(:, :), allocatable :: grad_psit_X, grad_psit_Y, grad_psit_Z
      real, dimension(:, :), allocatable :: grad_alpha_X, grad_alpha_Y, grad_alpha_Z
      real, dimension(:, :), allocatable :: B_cross_grad_s_dot_grad_alpha, B_cross_grad_B_dot_grad_psit
      real, dimension(:, :), allocatable :: B_cross_grad_B_dot_grad_alpha
      real, dimension(:, :), allocatable :: grad_B_X, grad_B_Y, grad_B_Z
      real, dimension(:, :), allocatable :: B_X, B_Y, B_Z
      real, dimension(:, :), allocatable :: gradzeta_gradalpha, gradzeta_gradpsit
      real, dimension(:, :), allocatable :: gradthetap_gradalpha, gradthetap_gradpsit
      real :: ds 

      !-------------------------------------------------------------------------

      ! Read in equilibrium information from VMEC file, this is stored as a set 
      ! of global variables in <read_wout_mod> in mini_libstell, which will be accessible
      call read_vmec_equilibrium(vmec_filename, verbose, ierr)

      ! Check the input variables
      call sanity_checks_inputs(nalpha, nzgrid, s_inputfile)

      ! Do some sanity checking to ensure the VMEC arrays have some expected properties
      call sanity_checks_vmec()

      ! Length of the simulated field line 
      number_of_field_periods = number_of_field_periods_inputfile
      if (number_of_field_periods_inputfile <= 0) number_of_field_periods = nfp 

      ! Silly step because <nfp_out> is a real, and <nfp> is an integer
      nfp_out = nfp 

      ! Define the toroidal flux at the last closed flux surface 
      edge_toroidal_flux_over_2pi = phi(ns) / (2 * pi) * isigng 
      sign_toroidal_flux = int(sign(1.1, edge_toroidal_flux_over_2pi)) 

      ! Set reference length and magnetic field for stella's normalization
      ! Using the choices made by Pavlos Xanthopoulos in GIST
      B_reference = 2 * abs(edge_toroidal_flux_over_2pi) / (Aminor * Aminor); 
      L_reference = Aminor; a = L_reference; bref = B_reference

      ! The radial coordinate used in VMEC is the normalized toroidal flux
      !     <s> = psi_toroidal / psi_{toroidal,edge} 
      ! Determine which flux surface to use, based on <s_inputfile> and <vmec_surface_option>. 
      ! In general, we get quantities for stella by linear interpolation, taking a weighted average of the quantity from
      ! 2 surfaces in the VMEC file. For any VMEC quantity Q on the full/half grid, the value used in stella will be
      !     Q_stella = Q(index_full(1))*weight_full(1) + Q(index_full(2))*weight_full(2) 
      !     Q_stella = Q(index_half(1))*weight_half(1) + Q(index_half(2))*weight_half(2)
      call get_chosen_flux_surface(s_inputfile, s, sqrt_s, vmec_surface_option, ierr)
 
      ! Evaluate several radial-profile functions at the flux surface we ended up choosing
      call calculate_quantities_on_fluxsurface(s, iota, safety_factor_q, shat) 

      ! Set up the <alpha> and <zeta> coordinate grids, where <zeta> covers <number_of_field_periods>
      alpha = [(alpha0 + ((j - 1) * 2 * pi) / nalpha, j=1, nalpha)]
      zeta = [(zeta_center + (pi * j * number_of_field_periods) / (nfp * nzgrid), j=-nzgrid, nzgrid)]

      ! We need to know the <theta> coordinates that lie on these <alpha> field lines 
      !     <theta_pest> = alpha + iota * zeta = straight-field-line angle
      !     <theta> = theta_pest - Lambda = cylindrical angle theta 

      call calculate_theta(nzgrid, zeta, nalpha, alpha, iota, theta, ierr)

      ! Allocate the geometry arrays, and initialize them to zero
      call allocate_geometry_arrays() 

      ! Perform the sine/consine fourier transformations along the field lines, 
      ! e.g., we construct R(ialpha, izeta), using theta(ialpha, izeta).
      ! We obtain the following quantities as a fuction of (ialpha, izeta), 
      !     R, Z, B, sqrt_g, B_sup_theta, B_sup_zeta, B_sub_theta, B_sub_zeta, B_sub_s
      ! As well as the derivatives with respect to (s, zeta, theta),
      !     d_R_d_s, d_Z_d_s, d_Lambda_d_s, d_B_d_s
      !     d_R_d_zeta, d_Z_d_zeta, d_Lambda_d_zeta, d_B_d_zeta
      !     d_R_d_theta, d_Z_d_theta, d_Lambda_d_theta, d_B_d_theta
      call perform_sinecosine_fouriertransforms() 

      ! Use R(ialpha,izeta) and Z(ialpha,izeta), to compute 
      !     X = R * cos(zeta); Y = R * sin(zeta)
      !     grad_s_X, grad_s_Y, grad_s_Z
      !     grad_zeta_X, grad_zeta_Y, grad_zeta_Z
      !     grad_theta_X, grad_theta_Y, grad_theta_Z
      !     grad_theta_pest_X, grad_theta_pest_Y, grad_theta_pest_Z
      !     grad_psit_X, grad_psit_Y, grad_psit_Z
      !     grad_alpha_X, grad_alpha_Y, grad_alpha_Z
      !     grad_B_X, grad_B_Y, grad_B_Z; B_x, B_y, B_z
      call calculate_cartesian_gradient_vectors()

      ! Calculate the triple products B x ∇s . ∇α and B x ∇B . ∇α
      call calculate_triple_products()

      ! Track the code  
      if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::calculate_geometry_for_stella' 

      ! <bmag> = B̃ = B/Bref 
      bmag = B / Bref

      ! <b_dot_grad_zeta> = (a/B) B . ∇ζ = (a/B) B_zeta 
      b_dot_grad_zeta = (a/B) * B_sup_zeta 

      ! <grad_alpha_grad_alpha> = (a^2) ∇α . ∇α
      ! <grad_alpha_grad_psit> = (1/Bref) ∇α . ∇ψt 
      ! <grad_psit_grad_psit> = (1/(a^2*Bref^2)) ∇ψt . ∇ψt 
      grad_alpha_grad_alpha = (a*a) * (grad_alpha_X * grad_alpha_X + grad_alpha_Y * grad_alpha_Y + grad_alpha_Z * grad_alpha_Z) 
      grad_alpha_grad_psit = (1./Bref) * (grad_alpha_X * grad_psit_X + grad_alpha_Y * grad_psit_Y + grad_alpha_Z * grad_psit_Z)   
      grad_psit_grad_psit = 1/(a*a*Bref*Bref) * (grad_psit_X * grad_psit_X + grad_psit_Y * grad_psit_Y + grad_psit_Z * grad_psit_Z)  

      ! Note that we don't normalize ∇ζ, only ∇ψt and ∇α (not even the ∇ in ∇ζ)
      ! <gradzeta_gradpsit> = 1/(a*Bref) ∇ζ . ∇ψt = ∇ζ . tilde{∇} tilde{ψt} 
      ! <gradzeta_gradalpha> = a ∇ζ . ∇α
      ! <gradthetap_gradpsit> = 1/(a*Bref) ∇θp . ∇ψt
      ! <gradthetap_gradalpha> = a ∇θp . ∇α
      gradzeta_gradpsit = 1./(a*Bref) * (grad_zeta_X * grad_psit_X + grad_zeta_Y * grad_psit_Y + grad_zeta_Z * grad_psit_Z)
      gradzeta_gradalpha = a * (grad_zeta_X * grad_alpha_X + grad_zeta_Y * grad_alpha_Y + grad_zeta_Z * grad_alpha_Z)
      gradthetap_gradpsit = 1./(a*Bref) * (grad_theta_pest_X * grad_psit_X + grad_theta_pest_Y * grad_psit_Y + grad_theta_pest_Z * grad_psit_Z)
      gradthetap_gradalpha = a * (grad_theta_pest_X * grad_alpha_X + grad_theta_pest_Y * grad_alpha_Y + grad_theta_pest_Z * grad_alpha_Z)
      
      ! For the momentum flux we need (R^2/B^2) ∇ζ . ∇α and (R^2/B^2) ∇ζ . ∇ψt
      gradzeta_gradpsit_R2overB2 = gradzeta_gradpsit * R**2 / (bmag*bmag)
      gradzeta_gradalpha_R2overB2 = gradzeta_gradalpha * R**2 / (bmag*bmag)
      b_dot_grad_zeta_RR = b_dot_grad_zeta * R**2

      ! Define <gds23_psitalpha> = -1/B̃^2 * [(∇α . ∇ζ) * (∇ψt . ∇α)  - (∇ψt . ∇ζ) * |∇α|^2]
      ! Define <gds24_psitalpha> = -1/B̃^2 * [(∇α . ∇ζ) * |∇ψt|^2     - (∇ψt . ∇ζ) * (∇ψt . ∇α)]
      ! Define <gds25_psitalpha> = -sgn(ψt)/B̃^2 * [(∇α . ∇θp) * (∇ψt . ∇α) - (∇ψt . ∇θp) * |∇α|^2]
      ! Define <gds26_psitalpha> = -0.5*sgn(ψt)/B̃^2 * [(∇α . ∇θp) * |∇ψt|^2    - (∇ψt . ∇θp) * (∇ψt . ∇α)]
      gds23_psitalpha = -1./(bmag*bmag) * (gradzeta_gradalpha * grad_alpha_grad_psit - gradzeta_gradpsit * grad_alpha_grad_alpha)
      gds24_psitalpha = -1./(bmag*bmag) * (gradzeta_gradalpha * grad_psit_grad_psit - gradzeta_gradpsit * grad_alpha_grad_psit)
      gds25_psitalpha = -sign_toroidal_flux/(bmag*bmag) * (gradthetap_gradalpha * grad_alpha_grad_psit - gradthetap_gradpsit * grad_alpha_grad_alpha)
      gds26_psitalpha = -sign_toroidal_flux/(2.*bmag*bmag) * (gradthetap_gradalpha * grad_psit_grad_psit - gradthetap_gradpsit * grad_alpha_grad_psit)

      ! <gbdrift_alpha> = 2*a^2*Bref/B^3 (B x ∇B . ∇α)
      ! <gbdrift0_psit> = 2*hat{s}/B^3 (B x ∇B . ∇ψt) 
      ! <cvdrift_alpha> = 2*a^2*Bref/B^2 (B x kappa . ∇α) = gbdrift_alpha + 2*a^2*Bref*mu0/B^4 (dp/ds) B x ∇s . ∇α
      ! <cvdrift0_psit> = 2*hat{s}/B^2 (B x kappa . ∇ψt) = 2*hat{s}/B^3 (B x ∇B . ∇ψt) = gbdrift0_psit
      gbdrift_alpha = 2.*a*a*Bref/(B*B*B) * B_cross_grad_B_dot_grad_alpha
      gbdrift0_psit = 2.*shat/(B*B*B) * B_cross_grad_B_dot_grad_psit 
      cvdrift_alpha = gbdrift_alpha + 2.*a*a*Bref*mu_0/(B*B*B*B) * d_pressure_d_s * B_cross_grad_s_dot_grad_alpha 
      cvdrift0_psit = gbdrift0_psit

      ! Ratio of the physical displacement due to movement in the stella 
      ! x-coordinate to the x-coordinate itself: |ds/dx|*sqrt((dR/ds)^2+(dZ/ds)^2)
      ! We do not yet define x, so calculate the psit_displacement: |d(s/a)/d(psit/a^2Bref)|*sqrt((dR/ds)^2+(dZ/ds)^2)
      ! And we have d(s/a)/d(psit/a^2*Bref) = a*Bref ds/dpsit = a*Bref (1/psi_{t,LCFS})
      psit_displacement_fac = (a*Bref/edge_toroidal_flux_over_2pi) * sqrt(d_R_d_s**2 + d_Z_d_s**2)  

      ! Deallocate the local geometry arrays
      call deallocate_geometry_arrays()

      ! Print some information to the output file if <verbose> = True
      call print_variables_to_outputfile()

  contains

!###############################################################################
!############################# CALCULATE GEOMETRY ##############################
!############################################################################### 

      !*************************************************************************
      !                 CALCULATE TRIPLE PRODUCTS OF B, ∇s, ∇α                !
      !*************************************************************************
      subroutine calculate_triple_products()
  
         implicit none 

         integer :: izeta
         real, dimension(:, :), allocatable :: B_cross_grad_s_dot_grad_alpha_v2
         real, dimension(:, :), allocatable :: B_cross_grad_s_dot_grad_alpha_v3
         real, dimension(:, :), allocatable :: B_cross_grad_B_dot_grad_alpha_v2

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::calculate_triple_products' 

         ! Allocate local arrays 
         allocate (B_cross_grad_B_dot_grad_alpha_v2(nalpha, -nzgrid:nzgrid))
         allocate (B_cross_grad_s_dot_grad_alpha_v2(nalpha, -nzgrid:nzgrid))
         allocate (B_cross_grad_s_dot_grad_alpha_v3(nalpha, -nzgrid:nzgrid))

         ! Calculate B x ∇B . ∇ψt, using psi_t = s*psi_{t,LCFS} and 1/sqrt(g) = ∇s . ∇θ x ∇ζ
         ! B x ∇B . ∇ψt = psi_{t,LCFS} B x (dB/ds ∇s + dB/dθ ∇θ + dB/dζ ∇ζ) . ∇s 
         !               = psi_{t,LCFS} (dB/dθ B x ∇θ . ∇s + dB/dζ B x ∇ζ . ∇s)  
         !               = psi_{t,LCFS} (dB/dθ B_zeta ∇ζ x ∇θ . ∇s + dB/dζ B_theta ∇θ x ∇ζ . ∇s)  
         !               = psi_{t,LCFS}/sqrt(g) (- dB/dθ B_zeta + dB/dζ B_theta)  
         B_cross_grad_B_dot_grad_psit = edge_toroidal_flux_over_2pi / sqrt_g * (B_sub_theta * d_B_d_zeta - B_sub_zeta * d_B_d_theta) 

         ! Calculate B x ∇s . ∇α, using 1/sqrt(g) = ∇s . ∇θ x ∇ζ and B = B_s ∇s + B_theta ∇θ + B_zeta ∇ζ
         ! B x ∇s . ∇α = (B_s ∇s + B_theta ∇θ + B_zeta ∇ζ) x ∇s . ∇α 
         !              = (B_theta ∇θ + B_zeta ∇ζ) x ∇s . ∇(θ + λ - iota*ζ)
         !              = B_theta ∇θ x ∇s . ∇(λ - iota*ζ) + B_zeta ∇ζ x ∇s . ∇(θ + λ) 
         !              = B_theta (dλ/dζ - iota) ∇θ x ∇s . ∇ζ + B_zeta (1 + dλ/dθ) ∇ζ x ∇s . ∇θ
         !              = 1/sqrt(g) ( - B_theta (dλ/dζ - iota) + B_zeta (1 + dλ/dθ) )
         B_cross_grad_s_dot_grad_alpha = (B_sub_zeta * (1 + d_Lambda_d_theta) - B_sub_theta * (d_Lambda_d_zeta - iota)) / sqrt_g

         ! We can also directly calculate B x ∇s . ∇α  
         B_cross_grad_s_dot_grad_alpha_v2 = 0 &
               + B_X * grad_s_Y * grad_alpha_Z + B_Y * grad_s_Z * grad_alpha_X + B_Z * grad_s_X * grad_alpha_Y &
               - B_Z * grad_s_Y * grad_alpha_X - B_X * grad_s_Z * grad_alpha_Y - B_Y * grad_s_X * grad_alpha_Z

         ! Calculate B x ∇s . ∇α, using psi_t = s*psi_{t,LCFS} and B = ∇psi_t x ∇α = psi_{t,LCFS} ∇s x ∇α
         !  B x ∇s . ∇α = ∇s x ∇α . B = (1/psi_{t,LCFS}) B^2
         B_cross_grad_s_dot_grad_alpha_v3 = (1/edge_toroidal_flux_over_2pi) * (B_X * B_X + B_Y * B_Y + B_Z * B_Z)

         ! Check that the three methods of calculating B x ∇s . ∇α agree
         call check_that_arrays_match(B_cross_grad_s_dot_grad_alpha, B_cross_grad_s_dot_grad_alpha_v2, 1.0e-2, 'B_cross_grad_s_dot_grad_alpha')
         call check_that_arrays_match(B_cross_grad_s_dot_grad_alpha, B_cross_grad_s_dot_grad_alpha_v3, 1.0e-2, 'B_cross_grad_s_dot_grad_alpha')

         ! Calculate B x ∇B . ∇α, assuming alpha = theta_pest - iota * zeta  
         if (.not. rectangular_cross_section) then
            do izeta = -nzgrid, nzgrid
               B_cross_grad_B_dot_grad_alpha(:, izeta) = 0 &
                     + (B_sub_s(:, izeta) * d_B_d_theta(:, izeta) * (d_Lambda_d_zeta(:, izeta) - iota) &
                     + B_sub_theta(:, izeta) * d_B_d_zeta(:, izeta) * (d_Lambda_d_s(:, izeta) - zeta(izeta) * d_iota_d_s) &
                     + B_sub_zeta(:, izeta) * d_B_d_s(:, izeta) * (1 + d_Lambda_d_theta(:, izeta)) &
                     - B_sub_zeta(:, izeta) * d_B_d_theta(:, izeta) * (d_Lambda_d_s(:, izeta) - zeta(izeta) * d_iota_d_s) &
                     - B_sub_theta(:, izeta) * d_B_d_s(:, izeta) * (d_Lambda_d_zeta(:, izeta) - iota) &
                     - B_sub_s(:, izeta) * d_B_d_zeta(:, izeta) * (1 + d_Lambda_d_theta(:, izeta))) / sqrt_g(:, izeta)
            end do
         end if 

         ! Calculate B x ∇B . ∇α, assuming alpha = theta_pest - iota * (zeta - zeta_center)
         if (rectangular_cross_section) then
            do izeta = -nzgrid, nzgrid
               B_cross_grad_B_dot_grad_alpha(:, izeta) = 0 &
                     + (B_sub_s(:, izeta) * d_B_d_theta(:, izeta) * (d_Lambda_d_zeta(:, izeta) - iota) &
                     + B_sub_theta(:, izeta) * d_B_d_zeta(:, izeta) * (d_Lambda_d_s(:, izeta) - (zeta(izeta) - zeta_center) * d_iota_d_s) &
                     + B_sub_zeta(:, izeta) * d_B_d_s(:, izeta) * (1 + d_Lambda_d_theta(:, izeta)) &
                     - B_sub_zeta(:, izeta) * d_B_d_theta(:, izeta) * (d_Lambda_d_s(:, izeta) - (zeta(izeta) - zeta_center) * d_iota_d_s) &
                     - B_sub_theta(:, izeta) * d_B_d_s(:, izeta) * (d_Lambda_d_zeta(:, izeta) - iota) &
                     - B_sub_s(:, izeta) * d_B_d_zeta(:, izeta) * (1 + d_Lambda_d_theta(:, izeta))) / sqrt_g(:, izeta)
            end do
         end if 

         ! Calculate B x ∇B . ∇α directly
         B_cross_grad_B_dot_grad_alpha_v2 = 0 &
                     + B_X * grad_B_Y * grad_alpha_Z + B_Y * grad_B_Z * grad_alpha_X + B_Z * grad_B_X * grad_alpha_Y &
                     - B_Z * grad_B_Y * grad_alpha_X - B_X * grad_B_Z * grad_alpha_Y - B_Y * grad_B_X * grad_alpha_Z

         ! Check that the two methods of calculating B x ∇B . ∇α  agree
         call check_that_arrays_match(B_cross_grad_B_dot_grad_alpha, B_cross_grad_B_dot_grad_alpha_v2, 1.0e-2, 'B_cross_grad_B_dot_grad_alpha')

         ! Deallocate local arrays 
         deallocate (B_cross_grad_B_dot_grad_alpha_v2)
         deallocate (B_cross_grad_s_dot_grad_alpha_v2)
         deallocate (B_cross_grad_s_dot_grad_alpha_v3)

      end subroutine calculate_triple_products

      !*************************************************************************
      !       CALCULATE THE CARTESIAN COMPONENTS OF THE GRADIENT VECTORS       !
      !*************************************************************************
      ! Use R(ialpha,izeta) and Z(ialpha,izeta), to compute 
      !     X = R * cos(zeta)
      !     Y = R * sin(zeta)
      ! 
      ! We can write the magnetic field in covariant representation
      ! (the components transform like the transformation of the reference axes)
      !     B = B_s ∇s + B_theta ∇θ + B_zeta ∇ζ
      ! Where be assume that B_s = 0 for a vacuum magnetic field
      !     B = B_theta ∇θ + B_zeta ∇ζ
      ! 
      ! We can also write the magnetic field in contravariant representation
      ! (the components transform inversely to the transformation of the reference axes)
      !     B = (1/sqrt(g)) (dpsi_t/ds) (iota * e_theta + e_zeta)
      ! So we have a tangent basis (∇s, ∇θ, ∇ζ) and the dual basis (e_s, e_theta, e_zeta)
      !     ∇s . e_s = 1         ∇θ . e_theta = 1     ∇ζ . e_zeta = 1
      !     ∇s . e_theta = 0     ∇s . e_zeta = 0      ∇θ . e_s = 0      ...
      !     e_s = (dX/ds) e_X + (dY/ds) e_Y + (dZ/ds) e_Z
      !     e_ζ = (dX/dζ) e_X + (dY/dζ) e_Y + (dZ/dζ) e_Z
      !     e_θ = (dX/dθ) e_X + (dY/dθ) e_Y + (dZ/dθ) e_Z
      ! 
      ! The two bases are related through the dual relations
      !     sqrt(g) = 1 / (∇s . ∇θ x ∇ζ)  
      !     e_i = sqrt(g) (e^j x e^k)
      !     e_s = sqrt(g) (∇θ x ∇ζ)
      !     e_ζ = sqrt(g) (∇s x ∇θ)
      !     e_θ = sqrt(g) (∇ζ x ∇s)
      !     ∇s = 1/sqrt(g) (e_θ x e_ζ)
      !     ∇ζ = 1/sqrt(g) (e_s x e_θ)
      !     ∇θ = 1/sqrt(g) (e_ζ x e_s)
      ! 
      ! If we multiply the covariant representation with the contravariant representation we get
      !     B^2 = (1/sqrt(g)) (dpsi_t/ds) (iota*B_theta + B_zeta)
      ! The Jacobian sqrt(g) of the magnetic coordinate system can thus be calculated as
      !     sqrt(g) = (dpsi_t/ds) (iota*B_theta + B_zeta)/B^2
      ! 
      ! Use the dual relations to get the Cartesian components of ∇s, ∇θ, and ∇ζ
      !     ∇s = 1/sqrt(g) (e_θ x e_ζ) = 1/sqrt(g) ((dX/dθ)e_X + (dY/dθ)e_Y + (dZ/dθ)e_Z) x ((dX/dζ)e_X + (dY/dζ)e_Y + (dZ/dζ)e_Z)
      !     ∇s . e_X = 1/sqrt(g) ( dY/dθ * dZ/dζ - dZ/dθ * dY/dζ ) 
      !     ∇s . e_Y = 1/sqrt(g) ( dZ/dθ * dX/dζ - dX/dθ * dZ/dζ )
      !     ∇s . e_Z = 1/sqrt(g) ( dX/dθ * dY/dζ - dY/dθ * dX/dζ ) 
      !*************************************************************************

      subroutine calculate_cartesian_gradient_vectors()
  
         implicit none 

         integer :: izeta
         real :: cos_angle, sin_angle

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::calculate_cartesian_gradient_vectors' 
   
         ! Iterate over the <zeta> grid
         do izeta = -nzgrid, nzgrid

            ! We need cos(zeta) and sin(zeta) to transform (R) to (X,Y)
            cos_angle = cos(zeta(izeta))
            sin_angle = sin(zeta(izeta))

            ! X = R*cos(ζ); dX/dζ = (dR/dζ)*cos(ζ) - R*sin(ζ); dX/dθ = dR/dθ * cos(ζ)
            d_X_d_s(:, izeta) = d_R_d_s(:, izeta) * cos_angle
            d_X_d_zeta(:, izeta) = d_R_d_zeta(:, izeta) * cos_angle - R(:, izeta) * sin_angle
            d_X_d_theta(:, izeta) = d_R_d_theta(:, izeta) * cos_angle

            ! Y = R*sin(ζ); dY/dζ = (dR/dζ)*sin(ζ) + R*cos(ζ); dY/dθ = dR/dθ * sin(ζ)
            d_Y_d_s(:, izeta) = d_R_d_s(:, izeta) * sin_angle
            d_Y_d_zeta(:, izeta) = d_R_d_zeta(:, izeta) * sin_angle + R(:, izeta) * cos_angle
            d_Y_d_theta(:, izeta) = d_R_d_theta(:, izeta) * sin_angle

         end do

         ! Use the dual relations to get the Cartesian components of ∇s, ∇θ, and ∇ζ
         !     ∇s = 1/sqrt(g) (e_θ x e_ζ) = 1/sqrt(g) ((dX/dθ)e_X + (dY/dθ)e_Y + (dZ/dθ)e_Z) x ((dX/dζ)e_X + (dY/dζ)e_Y + (dZ/dζ)e_Z)
         !     ∇s . e_X = 1/sqrt(g) ( dY/dθ * dZ/dζ - dZ/dθ * dY/dζ ) 
         !     ∇s . e_Y = 1/sqrt(g) ( dZ/dθ * dX/dζ - dX/dθ * dZ/dζ )
         !     ∇s . e_Z = 1/sqrt(g) ( dX/dθ * dY/dζ - dY/dθ * dX/dζ ) 
         grad_s_X = (d_Y_d_theta * d_Z_d_zeta - d_Z_d_theta * d_Y_d_zeta) / sqrt_g
         grad_s_Y = (d_Z_d_theta * d_X_d_zeta - d_X_d_theta * d_Z_d_zeta) / sqrt_g
         grad_s_Z = (d_X_d_theta * d_Y_d_zeta - d_Y_d_theta * d_X_d_zeta) / sqrt_g
         grad_zeta_X = (d_Y_d_s * d_Z_d_theta - d_Z_d_s * d_Y_d_theta) / sqrt_g
         grad_zeta_Y = (d_Z_d_s * d_X_d_theta - d_X_d_s * d_Z_d_theta) / sqrt_g
         grad_zeta_Z = (d_X_d_s * d_Y_d_theta - d_Y_d_s * d_X_d_theta) / sqrt_g
         grad_theta_X = (d_Y_d_zeta * d_Z_d_s - d_Z_d_zeta * d_Y_d_s) / sqrt_g
         grad_theta_Y = (d_Z_d_zeta * d_X_d_s - d_X_d_zeta * d_Z_d_s) / sqrt_g
         grad_theta_Z = (d_X_d_zeta * d_Y_d_s - d_Y_d_zeta * d_X_d_s) / sqrt_g

         ! Calculate <grad_theta_pest> = ∇(θ + λ) with λ = λ(s, θ, ζ)
         !     Recall that ∇f = (df/ds) ∇s + (df/dθ) ∇θ + (df/dζ) ∇ζ 
         !     ∇theta_pest = (dλ/ds) ∇s + (1 + dλ/dθ) ∇s + 0*∇ζ
         grad_theta_pest_X = d_Lambda_d_s * grad_s_X + (1.0 + d_Lambda_d_theta) * grad_theta_X + d_Lambda_d_zeta * grad_zeta_X  
         grad_theta_pest_Y = d_Lambda_d_s * grad_s_Y + (1.0 + d_Lambda_d_theta) * grad_theta_Y + d_Lambda_d_zeta * grad_zeta_Y 
         grad_theta_pest_Z = d_Lambda_d_s * grad_s_Z + (1.0 + d_Lambda_d_theta) * grad_theta_Z + d_Lambda_d_zeta * grad_zeta_Z 

         ! Sanity checks: dζ/dX = -sin(ζ)/R; dζ/dY = cos(ζ)/R and dζ/dZ = 0
         call sanity_check_grad_zeta()

         ! In VMEC s = psi_t/psi_{t,LCFS} so psi_t = s*psi_{t,LCFS}
         ! ∇psi_t = (dpsi_t/ds) * ∇s + 0*∇θ + 0*∇ζ = d(psi_{t,LCFS}*s)/ds * ∇s = psi_{t,LCFS} * ∇s 
         grad_psit_X = grad_s_X * edge_toroidal_flux_over_2pi
         grad_psit_Y = grad_s_Y * edge_toroidal_flux_over_2pi
         grad_psit_Z = grad_s_Z * edge_toroidal_flux_over_2pi

         ! ∇alpha = ∇(theta_pest - iota * zeta) = ∇(theta + Lambda - iota * zeta)
         !        = (dλ/ds - ζ (diota/ds) ∇s + (1 + dλ/dθ) ∇θ + (dλ/dζ - iota) ∇ζ
         ! First add the part which depends on zeta, assuming alpha = theta_pest - iota * zeta
         if (.not. rectangular_cross_section) then
            do izeta = -nzgrid, nzgrid
               grad_alpha_X(:, izeta) = (d_Lambda_d_s(:, izeta) - (zeta(izeta)) * d_iota_d_s) * grad_s_X(:, izeta)
               grad_alpha_Y(:, izeta) = (d_Lambda_d_s(:, izeta) - (zeta(izeta)) * d_iota_d_s) * grad_s_Y(:, izeta)
               grad_alpha_Z(:, izeta) = (d_Lambda_d_s(:, izeta) - (zeta(izeta)) * d_iota_d_s) * grad_s_Z(:, izeta)
            end do
         end if 

         ! If we want to ensure that the perpendicular cross-section of the flux tube is rectangular, 
         ! then we need to set alpha = theta_pest - iota * (zeta - zeta_center) instead 
         if (rectangular_cross_section) then
            do izeta = -nzgrid, nzgrid
               grad_alpha_X(:, izeta) = (d_Lambda_d_s(:, izeta) - (zeta(izeta) - zeta_center) * d_iota_d_s) * grad_s_X(:, izeta)
               grad_alpha_Y(:, izeta) = (d_Lambda_d_s(:, izeta) - (zeta(izeta) - zeta_center) * d_iota_d_s) * grad_s_Y(:, izeta)
               grad_alpha_Z(:, izeta) = (d_Lambda_d_s(:, izeta) - (zeta(izeta) - zeta_center) * d_iota_d_s) * grad_s_Z(:, izeta)
            end do
         end if 

         ! Now add the final two terms in ∇alpha = (dλ/ds - ζ (diota/ds) ∇s + (1 + dλ/dθ) ∇θ + (dλ/dζ - iota) ∇ζ
         grad_alpha_X = grad_alpha_X + (1 + d_Lambda_d_theta) * grad_theta_X + (-iota + d_Lambda_d_zeta) * grad_zeta_X
         grad_alpha_Y = grad_alpha_Y + (1 + d_Lambda_d_theta) * grad_theta_Y + (-iota + d_Lambda_d_zeta) * grad_zeta_Y
         grad_alpha_Z = grad_alpha_Z + (1 + d_Lambda_d_theta) * grad_theta_Z + (-iota + d_Lambda_d_zeta) * grad_zeta_Z

         ! ∇B = (dB/ds) ∇s + (dB/dθ) ∇θ + (dB/dζ) ∇ζ
         grad_B_X = d_B_d_s * grad_s_X + d_B_d_theta * grad_theta_X + d_B_d_zeta * grad_zeta_X
         grad_B_Y = d_B_d_s * grad_s_Y + d_B_d_theta * grad_theta_Y + d_B_d_zeta * grad_zeta_Y
         grad_B_Z = d_B_d_s * grad_s_Z + d_B_d_theta * grad_theta_Z + d_B_d_zeta * grad_zeta_Z

         ! Recall that psi_t = s*psi_{t,LCFS}, and alpha = theta + Lambda - iota * zeta and
         !     ∇psi_t and ∇s are parallel vectors, thus ∇psi_t x ∇s = 0
         !     sqrt(g) = 1 / (∇s . ∇θ x ∇ζ) and ∇s = 1/sqrt(g) (e_θ x e_ζ)
         !     e_ζ = sqrt(g) (∇s x ∇θ) = (dX/dζ) e_X + (dY/dζ) e_Y + (dZ/dζ) e_Z
         !     e_θ = sqrt(g) (∇ζ x ∇s) = (dX/dθ) e_X + (dY/dθ) e_Y + (dZ/dθ) e_Z 
         ! Using the formulas above we can calculate the cartesian components of the magnetic field as
         !       B = ∇psi_t x ∇alpha = ∇psi_t x ( (dλ/ds - ζ (diota/ds) ∇s + (1 + dλ/dθ) ∇θ + (dλ/dζ - iota) ∇ζ )
         !         = psi_{t,LCFS} ( (1 + dλ/dθ) ∇s x ∇θ + (dλ/dζ - iota) ∇s x ∇ζ )
         !         = psi_{t,LCFS}/sqrt(g) ( (1 + dλ/dθ) e_ζ - (dλ/dζ - iota) e_θ )
         ! B . e_x = psi_{t,LCFS}/sqrt(g) ( (1 + dλ/dθ) (dX/dζ) - (dλ/dζ - iota) (dX/dθ) )
         B_X = edge_toroidal_flux_over_2pi / sqrt_g * ((1 + d_Lambda_d_theta) * d_X_d_zeta - (d_Lambda_d_zeta - iota) * d_X_d_theta) 
         B_Y = edge_toroidal_flux_over_2pi / sqrt_g * ((1 + d_Lambda_d_theta) * d_Y_d_zeta - (d_Lambda_d_zeta - iota) * d_Y_d_theta)  
         B_Z = edge_toroidal_flux_over_2pi / sqrt_g * ((1 + d_Lambda_d_theta) * d_Z_d_zeta - (d_Lambda_d_zeta - iota) * d_Z_d_theta) 

         ! Check that sqrt(g) = 1 / (∇s . ∇θ x ∇ζ) and 1/sqrt(g) = (∇s . ∇θ x ∇ζ) 
         call sanity_check_jacobian()

         ! Check that B_sub_s = B . e_s; B_sub_zeta = B . e_ζ and B_sub_theta = B . e_θ
         ! Check that B_sup_s = B . ∇s = 0; B_sup_zeta = B . ∇ζ and B_sup_theta = B . ∇θ
         call sanity_check_Bcomponents()

      end subroutine calculate_cartesian_gradient_vectors

      !*************************************************************************
      !            PERFORM THE SINE/COSINE FOURIER TRANSFORMATIONS             !
      !*************************************************************************
      ! We constructed the arrays <zeta> and <theta> along the chosen
      ! field lines, so now we can perform the sine/consine fourier 
      ! transformations along these field lines.
      ! 
      ! Note that R, Z, and Lambda use the (xm,xn) mode numbers, while all the
      ! other quantities use the (xm_nyq,xn_nyq) mode numbers. 
      ! B and Lambda are on the half mesh, so their radial derivatives are on the full mesh.
      ! R and Z are on the full mesh, so their radial derivatives are on the half mesh.
      ! 
      ! For the derivatives with respect to <theta> and <zeta> we have, e.g.,
      !     R = sum_{imn} R_{imn}*cos(angle)
      !     dR/dtheta = sum_{imn} d(R_{imn}*cos(angle))/dtheta = - m * sum_{imn} R_{imn}*sin(angle)
      !     dR/dzeta = sum_{imn} d(R_{imn}*cos(angle))/dzeta = n * sum_{imn} R_{imn}*sin(angle)
      !************************************************************************* 

      subroutine perform_sinecosine_fouriertransforms()
  
         implicit none 
 
         real, dimension(:), allocatable :: d_Lambda_d_s_mnc, d_Lambda_d_s_mns
         real, dimension(:), allocatable :: d_B_d_s_mnc, d_B_d_s_mns
         real, dimension(:), allocatable :: d_R_d_s_mnc, d_R_d_s_mns
         real, dimension(:), allocatable :: d_Z_d_s_mnc, d_Z_d_s_mns
      
         real :: angle, cos_angle, sin_angle
         integer :: imn, ialpha, izeta
         integer :: m, n

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::perform_sinecosine_fouriertransforms' 

         ! Allocate the local arrays, needed for the radial derivatives
         allocate (d_Lambda_d_s_mnc(ns)); allocate (d_Lambda_d_s_mns(ns))
         allocate (d_B_d_s_mnc(ns)); allocate (d_B_d_s_mns(ns))
         allocate (d_R_d_s_mnc(ns)); allocate (d_R_d_s_mns(ns))
         allocate (d_Z_d_s_mnc(ns)); allocate (d_Z_d_s_mns(ns))

         !---------------------------------------------------------------------- 
         ! Perform the sine/consine fourier transformations along <mnmax>
         !---------------------------------------------------------------------- 

         ! Iterate over the <mnmax> mode numbers to construct R, Z and Lambda
         do imn = 1, mnmax

            ! Get the (m, n) mode numbers
            m = int(xm(imn))
            n = int(xn(imn))

            !-----------------------------------------------------
            ! First, consider just the stellarator-symmetric terms
            !-----------------------------------------------------

            ! For this <imn> mode we construct the radial derivatives that we will need.
            ! Note that B and Lambda are on the half mesh, so their radial derivatives are on the full mesh.
            ! R and Z are on the full mesh, so their radial derivatives are on the half mesh.
            d_Lambda_d_s_mns(2:ns - 1) = (lmns(imn, 3:ns) - lmns(imn, 2:ns - 1)) / ds
            d_R_d_s_mnc(2:ns) = (rmnc(imn, 2:ns) - rmnc(imn, 1:ns - 1)) / ds
            d_Z_d_s_mns(2:ns) = (zmns(imn, 2:ns) - zmns(imn, 1:ns - 1)) / ds

            ! We use a simplistic extrapolation at the endpoints.
            d_Lambda_d_s_mns(1) = d_Lambda_d_s_mns(2)
            d_Lambda_d_s_mns(ns) = d_Lambda_d_s_mns(ns - 1)
            d_R_d_s_mnc(1) = 0; d_Z_d_s_mns(1) = 0

            ! For each alpha-fieldline we know the (zeta, theta)-coordinates and we can calculate <angle> 
            do ialpha = 1, nalpha
               do izeta = -nzgrid, nzgrid

                  ! Calculate the angle used in the sine/consine fourier transformations
                  angle = m * theta(ialpha, izeta) - n * zeta(izeta)
                  cos_angle = cos(angle); sin_angle = sin(angle)

                  ! Handle R, which is on the full mesh
                  call radial_interpolation(rmnc(imn,:), temp, 'full') 
                  R(ialpha, izeta) = R(ialpha, izeta) + temp * cos_angle
                  d_R_d_theta(ialpha, izeta) = d_R_d_theta(ialpha, izeta) - temp * m * sin_angle
                  d_R_d_zeta(ialpha, izeta) = d_R_d_zeta(ialpha, izeta) + temp * n * sin_angle

                  ! Handle Z, which is on the full mesh
                  call radial_interpolation(zmns(imn,:), temp, 'full') 
                  Z(ialpha, izeta) = Z(ialpha, izeta) + temp * sin_angle
                  d_Z_d_theta(ialpha, izeta) = d_Z_d_theta(ialpha, izeta) + temp * m * cos_angle
                  d_Z_d_zeta(ialpha, izeta) = d_Z_d_zeta(ialpha, izeta) - temp * n * cos_angle

                  ! Handle Lambda, which is on the half mesh
                  call radial_interpolation(lmns(imn,:), temp, 'half') 
                  d_Lambda_d_theta(ialpha, izeta) = d_Lambda_d_theta(ialpha, izeta) + m * temp * cos_angle
                  d_Lambda_d_zeta(ialpha, izeta) = d_Lambda_d_zeta(ialpha, izeta) - n * temp * cos_angle

                  ! Handle dR/ds, since R is on the full mesh, its radial derivative is on the half mesh
                  call radial_interpolation(d_R_d_s_mnc(:), temp, 'half')  
                  d_R_d_s(ialpha, izeta) = d_R_d_s(ialpha, izeta) + temp * cos_angle

                  ! Handle dZ/ds, since Z is on the full mesh, its radial derivative is on the half mesh
                  call radial_interpolation(d_Z_d_s_mns(:), temp, 'half')   
                  d_Z_d_s(ialpha, izeta) = d_Z_d_s(ialpha, izeta) + temp * sin_angle

                  ! Handle dLambda/ds, since Lambda is on the half mesh, its radial derivative is on the full mesh
                  call radial_interpolation(d_Lambda_d_s_mns(:), temp, 'full')    
                  d_Lambda_d_s(ialpha, izeta) = d_Lambda_d_s(ialpha, izeta) + temp * sin_angle

               end do
            end do

            !-----------------------------------------------------
            ! Now consider the stellarator-asymmetric terms.
            !-----------------------------------------------------

            if (lasym) then

               ! For this <imn> mode we construct the radial derivatives that we will need.
               ! Note that B and Lambda are on the half mesh, so their radial derivatives are on the full mesh.
               ! R and Z are on the full mesh, so their radial derivatives are on the half mesh.
               d_Lambda_d_s_mnc(2:ns - 1) = (lmnc(imn, 3:ns) - lmnc(imn, 2:ns - 1)) / ds
               d_R_d_s_mns(2:ns) = (rmns(imn, 2:ns) - rmns(imn, 1:ns - 1)) / ds
               d_Z_d_s_mnc(2:ns) = (zmnc(imn, 2:ns) - zmnc(imn, 1:ns - 1)) / ds

               ! We use a simplistic extrapolation at the endpoints.
               d_Lambda_d_s_mnc(1) = d_Lambda_d_s_mnc(2)
               d_Lambda_d_s_mnc(ns) = d_Lambda_d_s_mnc(ns - 1)
               d_R_d_s_mns(1) = 0; d_Z_d_s_mnc(1) = 0 

               ! For each alpha-fieldline we know the (zeta, theta)-coordinates and we can calculate <angle> 
               do ialpha = 1, nalpha
                  do izeta = -nzgrid, nzgrid

                     ! Calculate the angle used in the sine/consine fourier transformations
                     angle = m * theta(ialpha, izeta) - n * zeta(izeta)
                     cos_angle = cos(angle); sin_angle = sin(angle) 

                     ! Handle R, which is on the full mesh
                     call radial_interpolation(rmns(imn,:), temp, 'full')  
                     R(ialpha, izeta) = R(ialpha, izeta) + temp * sin_angle
                     d_R_d_theta(ialpha, izeta) = d_R_d_theta(ialpha, izeta) + temp * m * cos_angle
                     d_R_d_zeta(ialpha, izeta) = d_R_d_zeta(ialpha, izeta) - temp * n * cos_angle

                     ! Handle Z, which is on the full mesh
                     call radial_interpolation(zmnc(imn,:), temp, 'full')   
                     Z(ialpha, izeta) = Z(ialpha, izeta) + temp * cos_angle  
                     d_Z_d_theta(ialpha, izeta) = d_Z_d_theta(ialpha, izeta) - temp * m * sin_angle
                     d_Z_d_zeta(ialpha, izeta) = d_Z_d_zeta(ialpha, izeta) + temp * n * sin_angle

                     ! Handle Lambda, which is on the half mesh
                     call radial_interpolation(lmnc(imn,:), temp, 'half')   
                     d_Lambda_d_theta(ialpha, izeta) = d_Lambda_d_theta(ialpha, izeta) - temp * m * sin_angle
                     d_Lambda_d_zeta(ialpha, izeta) = d_Lambda_d_zeta(ialpha, izeta) + temp * n * sin_angle

                     ! Handle dR/ds, since R is on the full mesh, its radial derivative is on the half mesh
                     call radial_interpolation(d_R_d_s_mns(:), temp, 'half')   
                     d_R_d_s(ialpha, izeta) = d_R_d_s(ialpha, izeta) + temp * sin_angle

                     ! Handle dZ/ds, since Z is on the full mesh, its radial derivative is on the half mesh
                     call radial_interpolation(d_Z_d_s_mnc(:), temp, 'half')   
                     d_Z_d_s(ialpha, izeta) = d_Z_d_s(ialpha, izeta) + temp * cos_angle

                     ! Handle dLambda/ds, since Lambda is on the half mesh, its radial derivative is on the full mesh
                     call radial_interpolation(d_Lambda_d_s_mnc(:), temp, 'full')    
                     d_Lambda_d_s(ialpha, izeta) = d_Lambda_d_s(ialpha, izeta) + temp * cos_angle

                  end do
               end do
            end if
         end do

         !---------------------------------------------------------------------- 
         ! Perform the sine/consine fourier transformations along <mnmax_nyq>
         !---------------------------------------------------------------------- 

         ! Iterate over the <mnmax> mode numbers to construct B, ... 
         do imn = 1, mnmax_nyq  

            ! Get the (m, n) mode numbers
            m = int(xm_nyq(imn))
            n = int(xn_nyq(imn))

            !-----------------------------------------------------
            ! First, consider just the stellarator-symmetric terms
            !-----------------------------------------------------

            ! For this <imn> mode we construct the radial derivatives that we will need.
            ! Note that B is on the half mesh, so their radial derivatives are on the full mesh.
            d_B_d_s_mnc(2:ns - 1) = (bmnc(imn, 3:ns) - bmnc(imn, 2:ns - 1)) / ds
            d_B_d_s_mnc(1) = d_B_d_s_mnc(2)
            d_B_d_s_mnc(ns) = d_B_d_s_mnc(ns - 1)

            ! For each alpha-fieldline we know the (zeta, theta)-coordinates and we can calculate <angle> 
            do ialpha = 1, nalpha
               do izeta = -nzgrid, nzgrid

                  ! Calculate the angle used in the sine/consine fourier transformations
                  angle = m * theta(ialpha, izeta) - n * zeta(izeta)
                  cos_angle = cos(angle); sin_angle = sin(angle)

                  ! Handle |B|, which is on the half mesh
                  call radial_interpolation(bmnc(imn,:), temp, 'half')  
                  B(ialpha, izeta) = B(ialpha, izeta) + temp * cos_angle
                  d_B_d_theta(ialpha, izeta) = d_B_d_theta(ialpha, izeta) - m * temp * sin_angle
                  d_B_d_zeta(ialpha, izeta) = d_B_d_zeta(ialpha, izeta) + n * temp * sin_angle

                  ! Handle Jacobian, which is on the half mesh
                  call radial_interpolation(gmnc(imn,:), temp, 'half')   
                  sqrt_g(ialpha, izeta) = sqrt_g(ialpha, izeta) + temp * cos_angle

                  ! Handle B sup theta, which is on the half mesh
                  call radial_interpolation(bsupumnc(imn,:), temp, 'half')    
                  B_sup_theta(ialpha, izeta) = B_sup_theta(ialpha, izeta) + temp * cos_angle

                  ! Handle B sup zeta, which is on the half mesh
                  call radial_interpolation(bsupvmnc(imn,:), temp, 'half')     
                  B_sup_zeta(ialpha, izeta) = B_sup_zeta(ialpha, izeta) + temp * cos_angle

                  ! Handle B sub theta, which is on the half mesh
                  call radial_interpolation(bsubumnc(imn,:), temp, 'half')  
                  B_sub_theta(ialpha, izeta) = B_sub_theta(ialpha, izeta) + temp * cos_angle

                  ! Handle B sub zeta, which is on the half mesh
                  call radial_interpolation(bsubvmnc(imn,:), temp, 'half') 
                  B_sub_zeta(ialpha, izeta) = B_sub_zeta(ialpha, izeta) + temp * cos_angle

                  ! Handle B sub psi, unlike the other components of B, this one is on the full mesh
                  call radial_interpolation(bsubsmns(imn,:), temp, 'full')  
                  B_sub_s(ialpha, izeta) = B_sub_s(ialpha, izeta) + temp * sin_angle

                  ! Handle dB/ds, since bmnc is on the half mesh, its radial derivative is on the full mesh
                  call radial_interpolation(d_B_d_s_mnc(:), temp, 'full')  
                  d_B_d_s(ialpha, izeta) = d_B_d_s(ialpha, izeta) + temp * cos_angle

               end do
            end do

            !-----------------------------------------------------
            ! Now consider the stellarator-asymmetric terms.
            !-----------------------------------------------------

            if (lasym) then

               ! For this <imn> mode we construct the radial derivatives that we will need.
               ! Note that B is on the half mesh, so their radial derivatives are on the full mesh.  
               d_B_d_s_mns(2:ns - 1) = (bmns(imn, 3:ns) - bmns(imn, 2:ns - 1)) / ds 
               d_B_d_s_mns(1) = d_B_d_s_mns(2)
               d_B_d_s_mns(ns) = d_B_d_s_mns(ns - 1)

               ! For each alpha-fieldline we know the (zeta, theta)-coordinates and we can calculate <angle> 
               do ialpha = 1, nalpha
                  do izeta = -nzgrid, nzgrid

                     ! Calculate the angle used in the sine/consine fourier transformations
                     angle = m * theta(ialpha, izeta) - n * zeta(izeta)
                     cos_angle = cos(angle)
                     sin_angle = sin(angle) 

                     ! Handle |B|, which is on the half mesh
                     call radial_interpolation(bmns(imn,:), temp, 'half')  
                     B(ialpha, izeta) = B(ialpha, izeta) + temp * sin_angle
                     d_B_d_theta(ialpha, izeta) = d_B_d_theta(ialpha, izeta) + m * temp * cos_angle
                     d_B_d_zeta(ialpha, izeta) = d_B_d_zeta(ialpha, izeta) - n * temp * cos_angle

                     ! Handle Jacobian, which is on the half mesh
                     call radial_interpolation(gmns(imn,:), temp, 'half')   
                     sqrt_g(ialpha, izeta) = sqrt_g(ialpha, izeta) + temp * sin_angle

                     ! Handle B sup theta, which is on the half mesh
                     call radial_interpolation(bsupumns(imn,:), temp, 'half')
                     B_sup_theta(ialpha, izeta) = B_sup_theta(ialpha, izeta) + temp * sin_angle

                     ! Handle B sup zeta, which is on the half mesh
                     call radial_interpolation(bsupvmns(imn,:), temp, 'half') 
                     B_sup_zeta(ialpha, izeta) = B_sup_zeta(ialpha, izeta) + temp * sin_angle

                     ! Handle B sub theta, which is on the half mesh
                     call radial_interpolation(bsubumns(imn,:), temp, 'half')  
                     B_sub_theta(ialpha, izeta) = B_sub_theta(ialpha, izeta) + temp * sin_angle

                     ! Handle B sub zeta, which is on the half mesh
                     call radial_interpolation(bsubvmns(imn,:), temp, 'half')   
                     B_sub_zeta(ialpha, izeta) = B_sub_zeta(ialpha, izeta) + temp * sin_angle

                     ! Handle B sub psi, unlike the other components of B, this one is on the full mesh
                     call radial_interpolation(bsubsmnc(imn,:), temp, 'full')    
                     B_sub_s(ialpha, izeta) = B_sub_s(ialpha, izeta) + temp * cos_angle

                     ! Handle dB/ds, since bmns is on the half mesh, its radial derivative is on the full mesh
                     call radial_interpolation(d_B_d_s_mns(:), temp, 'full')    
                     d_B_d_s(ialpha, izeta) = d_B_d_s(ialpha, izeta) + temp * sin_angle

                  end do
               end do
            end if
         end do

         ! Deallocate the local arrays 
         deallocate (d_Lambda_d_s_mnc); deallocate (d_Lambda_d_s_mns)
         deallocate (d_B_d_s_mnc); deallocate (d_B_d_s_mns)
         deallocate (d_R_d_s_mnc); deallocate (d_R_d_s_mns)
         deallocate (d_Z_d_s_mnc); deallocate (d_Z_d_s_mns)

         ! Sanity check: iota = (B dot grad theta_pest) / (B dot grad zeta)  
         call sanity_check_iota()

      end subroutine perform_sinecosine_fouriertransforms

!###############################################################################
!############################# LITTLE CALCULATIONS #############################
!###############################################################################


      !*************************************************************************
      !                     QUANTITIES ON THE FLUX SURFACE                     !
      !*************************************************************************
      ! Evaluate several radial-profile functions at the flux surface we ended up choosing
      ! in the <get_chosen_flux_surface> routine. Here <ns>, <iotas>, <iotaf> and 
      ! <presf> and is a module variables
      !*************************************************************************

      subroutine calculate_quantities_on_fluxsurface(s, iota, safety_factor_q, shat)

         use debug_flags, only: print_extra_info_to_terminal

         implicit none 
  
         real, intent(in) :: s
         real, intent(out) :: iota, safety_factor_q, shat
 
         real, dimension(:), allocatable :: d_pressure_d_s_on_half_grid, d_iota_d_s_on_half_grid

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::calculate_quantities_on_fluxsurface' 

         ! Allocate the local arrays
         allocate (d_iota_d_s_on_half_grid(ns)); d_iota_d_s_on_half_grid = 0
         allocate (d_pressure_d_s_on_half_grid(ns)); d_pressure_d_s_on_half_grid = 0

         ! Construct the (diota/ds) and (dpressure/ds) arrays on the half grid 
         d_iota_d_s_on_half_grid(2:ns) = (iotaf(2:ns) - iotaf(1:ns - 1)) / ds 
         d_pressure_d_s_on_half_grid(2:ns) = (presf(2:ns) - presf(1:ns - 1)) / ds

         ! Calculate <iota>, <d_iota_d_s> and <d_pressure_d_s> on the chosen flux surface <s>
         call radial_interpolation(iotas, iota, 'half') 
         call radial_interpolation(d_iota_d_s_on_half_grid, d_iota_d_s, 'half')  
         call radial_interpolation(d_pressure_d_s_on_half_grid, d_pressure_d_s, 'half')   
         
         ! The radial coordinate is r = a sqrt(s) = a sqrt(psi_toroidal / psi_{toroidal,edge})  
         !     <shat> = (r/q)(dq/dr) = - (r/iota) (diota/dr) = -2 (s/iota) (diota/ds)
         !     <safety_factor_q> = q = 1/iota
         shat = (-2 * s / iota) * d_iota_d_s
         safety_factor_q = 1 / iota

         ! Deallocate the local arrays
         deallocate (d_iota_d_s_on_half_grid)
         deallocate (d_pressure_d_s_on_half_grid)

         ! Print these quantities to the command prompt
         if (verbose .and. print_extra_info_to_terminal) then
            write (*, *) " "
            write (*, *) "  Radial-profile functions at the chosen flux surface:"
            write (*, '(A21, F15.12)') "      iota:"//repeat(' ', 50), iota
            write (*, '(A21, ES20.12E3)') "      ds:"//repeat(' ', 50), ds
            write (*, '(A21, ES20.12E3)') "      diota/ds:"//repeat(' ', 50), d_iota_d_s
            write (*, '(A21, ES20.12E3)') "      dpressure/ds:"//repeat(' ', 50), d_pressure_d_s
            write (*, *) " "
         end if

      end subroutine calculate_quantities_on_fluxsurface

      !*************************************************************************
      !                        Choose the flux surface                         !
      !*************************************************************************
      ! Determine which flux surface to use, based on <s_inputfile> 
      ! and <vmec_surface_option>. Possible values of vmec_surface_option:
      !     0 = Use the exact radius requested.
      !     1 = Use the nearest value of the VMEC half grid.
      !     2 = Use the nearest value of the VMEC full grid.
      ! Here <ns> is a module variable, and <zero>, <weight_full>, 
      ! <weight_half>, <index_full>, <index_half>
      ! are subroutine variables.
      !*************************************************************************
      subroutine get_chosen_flux_surface(s_inputfile, s, sqrt_s, vmec_surface_option, ierr)

         use mp, only: mp_abort

         implicit none 
 
         integer, intent(in) :: vmec_surface_option
         real, intent(in) :: s_inputfile 
         real, intent(out) :: s, sqrt_s 
         integer, intent(inout) :: ierr

         real, dimension(:), allocatable :: dr2, s_full_grid, s_half_grid 
         integer :: j, index_of_minimum_error
         real :: min_dr2

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::get_chosen_flux_surface' 

         ! Select the flux surface                    
         ! The radial coordinate used in VMEC is the normalized toroidal flux
         !     <s> = psi_toroidal / psi_{toroidal,edge} 
         ! There are <ns> flux surfaces, <s> ranges from 0 to 1, and <ds> is the step size
         allocate (s_full_grid(ns))
         s_full_grid = [(real(j - 1) / (ns - 1), j=1, ns)] 
         ds = s_full_grid(2) - s_full_grid(1)

         ! Build an array of the half grid points
         allocate (s_half_grid(ns - 1))
         do j = 1, ns - 1
            s_half_grid(j) = (s_full_grid(j) + s_full_grid(j + 1)) * (0.5d+0)
         end do

         ! We select the flux surface based on the <vmec_surface_option>
         select case (vmec_surface_option)

            ! Use exact radius requested.
            case (0) 
               s = s_inputfile

            ! Use nearest value of the VMEC half grid
            ! Compute differences <dr2> and find the index of minimum error
            case (1)
               allocate (dr2(ns - 1))
               dr2 = (s_half_grid - s_inputfile)**2
               index_of_minimum_error = 1
               min_dr2 = dr2(1)
               do j = 2, ns - 1
                  if (dr2(j) < min_dr2) then
                     index_of_minimum_error = j
                     min_dr2 = dr2(j)
                  end if
               end do
               s = s_half_grid(index_of_minimum_error)
               deallocate (dr2)

            ! Use nearest value of the VMEC full grid
            ! Compute differences <dr2> and find the index of minimum error
            case (2)
               allocate (dr2(ns))
               dr2 = (s_full_grid - s_inputfile)**2
               index_of_minimum_error = 1
               min_dr2 = dr2(1) 
               do j = 2, ns
                  if (dr2(j) < min_dr2) then
                     index_of_minimum_error = j
                     min_dr2 = dr2(j)
                  end if
               end do
               s = s_full_grid(index_of_minimum_error)
               deallocate (dr2)

            ! Exit stella if another <vmec_surface_option> was given
            case default
               write(*,*) "Error! vmec_surface_option must be 0, 1, or 2. It is instead ", vmec_surface_option
               ierr = 114; call mp_abort("Error! vmec_surface_option must be 0, 1, or 2.") 

         end select

         ! Calculate rho = sqrt(s) = sqrt(psi_toroidal / psi_{toroidal,edge})
         sqrt_s = sqrt(s)

         ! Get radial integration weights on the full grid    
         ! In general, we get quantities for stella by linear interpolation, taking a weighted average of the quantity from
         ! 2 surfaces in the VMEC file. Sometimes the weights are 0 and 1, i.e., no interpolation is needed.
         ! For any VMEC quantity Q on the full grid, the value used in stella will be
         !     Q_stella = Q(index_full(1))*weight_full(1) + Q(index_full(2))*weight_full(2)
         ! For any VMEC quantity Q on the half grid, the value used in stella will be
         !     Q_stella = Q(index_half(1))*weight_half(1) + Q(index_half(2))*weight_half(2)
         ! Exit stella if <s> falls outside of [0,1]
         if (s > 1) then
            write(*,*) "Error! s = psi_toroidal / psi_{toroidal,edge} cannot be >1"
            ierr = 115; call mp_abort("Error! s = psi_toroidal / psi_{toroidal,edge} cannot be >1") 
         elseif (s < 0) then
            write(*,*) "Error! s = psi_toroidal / psi_{toroidal,edge} cannot be <0"
            ierr = 116; call mp_abort("Error! s = psi_toroidal / psi_{toroidal,edge} cannot be >1") 

         ! Special case if we select the last closed flux surface
         elseif (s == 1) then
            index_full(1) = ns - 1
            index_full(2) = ns
            weight_full(1) = zero

         ! Generally, <s> is >= 0 and <1 
         else
            index_full(1) = floor(s * (ns - 1)) + 1
            index_full(2) = index_full(1) + 1
            weight_full(1) = index_full(1) - s * (ns - one)
         end if

         ! Get radial integration weights on the half grid         
         ! Special case if <s> is smaller than our first half grid point
         ! We start at element 2 since element 1 is always 0 for quantities on the half grid.
         if (s < s_half_grid(1)) then
            write(*,*) "Warning: extrapolating beyond the end of VMEC's half grid."
            write(*,*) "(Extrapolating towards the magnetic axis.) Results are likely to be inaccurate." 
            index_half(1) = 2; index_half(2) = 3
            weight_half(1) = (s_half_grid(2) - s) / (s_half_grid(2) - s_half_grid(1))

         ! Special case if <s> is larger than our last half grid point
         elseif (s > s_half_grid(ns - 1)) then
            write(*,*) "Warning: extrapolating beyond the end of VMEC's half grid."
            write(*,*) "(Extrapolating towards the last closed flux surface.) Results may be inaccurate."
            index_half(1) = ns - 1; index_half(2) = ns
            weight_half(1) = (s_half_grid(ns - 1) - s) / (s_half_grid(ns - 1) - s_half_grid(ns - 2))

         ! Special case if <s> is the last half grid point
         elseif (s == s_half_grid(ns - 1)) then 
            index_half(1) = ns - 1; index_half(2) = ns
            weight_half(1) = zero

         ! Generally, <s> is inside the half grid.
         else 
            index_half(1) = floor(s * (ns - 1) + 0.5d+0) + 1
            if (index_half(1) < 2) index_half(1) = 2 
            index_half(2) = index_half(1) + 1
            weight_half(1) = index_half(1) - s * (ns - one) - (0.5d+0)
         end if

         ! We integrate between two flux surfaces, so second weights = 1 - first weight
         weight_full(2) = one - weight_full(1)
         weight_half(2) = one - weight_half(1)

         ! Deallocate the local arrays
         deallocate (s_full_grid)
         deallocate (s_half_grid)

      end subroutine get_chosen_flux_surface

!###############################################################################
!################################ SANITY CHECKS ################################
!###############################################################################

      !*************************************************************************
      !                  Print information of the output file                  !
      !*************************************************************************

      subroutine print_variables_to_outputfile()

         use debug_flags, only: print_extra_info_to_terminal
         
         implicit none

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::print_variables_to_outputfile'

         ! Only continue if <verbose> is True and <print_extra_info_to_terminal> is True
         if (.not. (verbose .and. print_extra_info_to_terminal)) return  

         ! Length of the simulated field line 
         if (number_of_field_periods_inputfile <= 0)  then
            write(*,*) "   Since <nfield_periods> in the <vmec_parameters> knob was <= 0, it is being set to nfp =", nfp
         end if 

         ! Set reference length and magnetic field for stella's normalization
         ! Using the choices made by Pavlos Xanthopoulos in GIST  
         ! Sign of the toroidal flux at the last closed flux surface  
         ! So we know whether the VMEC is clockwise of counter-clockwise 
         write (*, *) "  Reference values for the stella normalization:"
         write (*, '(A42, F15.12, A7)') "      Reference length (minor radius a):"//repeat(' ', 50), L_reference, " meters"
         write (*, '(A42, F15.12, A6)') "      Reference magnetic field strength:"//repeat(' ', 50), B_reference, " Tesla"
         write (*, '(A43, I2)') "      Sign of the toroidal flux:"//repeat(' ', 50), sign_toroidal_flux
         write (*, *) "  "

      end subroutine print_variables_to_outputfile


      !*************************************************************************
      !                        Test the input variables                        !
      !*************************************************************************
      subroutine sanity_checks_inputs(nalpha, nzgrid, s_inputfile)

         use mp, only: mp_abort

         implicit none
 
         integer, intent(in) :: nalpha, nzgrid
         real, intent(in) :: s_inputfile  

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::sanity_checks_inputs'

         ! Check the variables in the input file
         if (nalpha < 1) then
            write(*,*) "Error! nalpha must be >= 1. Instead it is", nalpha
            ierr = 100; call mp_abort("Error! nalpha must be >= 1") 
         end if

         if (nzgrid < 1) then
            write(*,*) "Error! nzgrid must be >= 1. Instead it is", nzgrid
            ierr = 101; call mp_abort("Error! nzgrid must be >= 1") 
         end if

         if (s_inputfile <= 0) then
            write(*,*) "Error! s_inputfile must be > 0. Instead it is", s_inputfile
            ierr = 102; call mp_abort("Error! s_inputfile must be > 0") 
         end if

         if (s_inputfile > 1) then
            write(*,*) "Error! s_inputfile must be <= 1. Instead it is", s_inputfile
            ierr = 103; call mp_abort("Error! s_inputfile must be <= 1") 
         end if

      end subroutine sanity_checks_inputs

      !*************************************************************************
      !                          Test the VMEC arrays                          !
      !*************************************************************************
      ! Do some sanity checking to ensure the VMEC arrays have some expected properties.
      ! Here <xm>, <xn>, <xm_nyq>, <xn_nyq>, <phi>, <ns> are module variables
      subroutine sanity_checks_vmec()

         use mp, only: mp_abort

         implicit none

         real :: dphi

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::sanity_checks_vmec'
      
         ! There is a bug in libstell read_wout_file for ASCII-format wout files, in which the xm_nyq and xn_nyq 
         ! arrays are sometimes not populated. The next few lines here provide a workaround. 
         if (maxval(abs(xm_nyq)) < 1 .and. maxval(abs(xn_nyq)) < 1) then
            if (mnmax_nyq == mnmax) then
               if (verbose) write(*,*) "xm_nyq and xn_nyq arrays are not populated in the wout file. Using xm and xn instead."
               xm_nyq = xm; xn_nyq = xn
            else
               write(*,*) "Error! xm_nyq and xn_nyq arrays are not populated in the wout file, and mnmax_nyq != mnmax."
               ierr = 104; call mp_abort("Error! xm_nyq and xn_nyq arrays are not populated in the wout file, and mnmax_nyq != mnmax.") 
            end if
         end if

         ! <phi> is vmec's array of the toroidal flux (not divided by 2*pi!) on vmec's radial grid.
         ! <phi> needs to start with a zero, and the grid needs to be uniformly spaced
         if (abs(phi(1)) > 1d-14) then
            write(*,*) "Error! VMEC phi array does not begin with 0."
            write(*,*) "phi:", phi
            ierr = 105; call mp_abort("Error! VMEC phi array does not begin with 0.") 
         end if 
         dphi = phi(2) - phi(1)
         do j = 3, ns
            if (abs(phi(j) - phi(j - 1) - dphi) > 1d-11) then
               write(*,*) "Error! VMEC phi array is not uniformly spaced."
               write(*,*) "phi:", phi
               ierr = 106; call mp_abort("Error! VMEC phi array is not uniformly spaced.") 
            end if
         end do

         ! <phips> needs to be constant and equal to -phi(ns)/(2*pi)
         ! Here <ns> is the number of flux surfaces included in the VMEC 
         ! <phips> is defined on the half-mesh, so skip first point.
         do j = 2, ns
            if (abs(phip(j) + phi(ns) / (2 * pi)) > 1d-11) then
               write(*,*) "Error! VMEC phips array is not constant and equal to -phi(ns)/(2*pi)."
               write(*,*) "phip(s):", phip
               ierr = 107; call mp_abort("Error! VMEC phips array is not constant and equal to -phi(ns)/(2*pi).") 
            end if
         end do

         ! The first mode in the <xm> and <xn> arrays should be m=n=0:
         if (xm(1) /= 0) then
            write(*,*) "First element of xm in the wout file should be 0."
            ierr = 108; call mp_abort("First element of xm in the wout file should be 0.") 
         end if
         if (xn(1) /= 0) then
            write(*,*) "First element of xn in the wout file should be 0."
            ierr = 109; call mp_abort("First element of xn in the wout file should be 0.") 
         end if
         if (xm_nyq(1) /= 0) then
            write(*,*) "First element of xm_nyq in the wout file should be 0."
            ierr = 110; call mp_abort("First element of xm_nyq in the wout file should be 0.") 
         end if
         if (xn_nyq(1) /= 0) then
            write(*,*) "First element of xn_nyq in the wout file should be 0."
            ierr = 111; call mp_abort("First element of xn_nyq in the wout file should be 0.") 
         end if

         ! Lambda should be on the half mesh, so its value at radial index 1 should be 0 for all (m,n)
         if (maxval(abs(lmns(:, 1))) > 0) then
            write(*,*) "Error! Expected lmns to be on the half mesh, but its value at radial index 1 is nonzero."
            write(*,*) "Here comes lmns(:,1):", lmns(:, 1)
            ierr = 112; call mp_abort("Error! Expected lmns to be on the half mesh, but its value at radial index 1 is nonzero.") 
         end if
         if (lasym) then
            if (maxval(abs(lmnc(:, 1))) > 0) then
               write(*,*) "Error! Expected lmnc to be on the half mesh, but its value at radial index 1 is nonzero."
               write(*,*) "Here comes lmnc(:,1):", lmnc(:, 1)
               ierr = 113; call mp_abort("Error! Expected lmnc to be on the half mesh, but its value at radial index 1 is nonzero.") 
            end if
         end if

      end subroutine sanity_checks_vmec

      !**********************************************************************
      !                              Test iota                              !
      !**********************************************************************
      ! If the conversion to <theta_pest> has been done correctly, we should find that
      !     iota = (B dot grad theta_pest) / (B dot grad zeta)    
      !**********************************************************************

      subroutine sanity_check_iota()

         implicit none

         real, dimension(:, :), allocatable :: temp1_vs_alphazeta, temp2_vs_alphazeta

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::sanity_check_iota'
         
         ! Allocate temporary arrays
         allocate (temp1_vs_alphazeta(nalpha, -nzgrid:nzgrid)) 
         allocate (temp2_vs_alphazeta(nalpha, -nzgrid:nzgrid)) 
   
         ! Calculate (B dot grad theta_pest) / (B dot grad zeta)  
         temp1_vs_alphazeta = (B_sup_theta * (1 + d_Lambda_d_theta) + B_sup_zeta * d_Lambda_d_zeta) / B_sup_zeta

         ! Turn iota into an array vs (alpha, zeta) 
         temp2_vs_alphazeta = iota
   
         ! Compare iota with (B dot grad theta_pest) / (B dot grad zeta)    
         call check_that_arrays_match(temp1_vs_alphazeta, temp2_vs_alphazeta, 0.01, 'iota')

         ! Deallocate temporary arrays
         deallocate (temp1_vs_alphazeta)
         deallocate (temp2_vs_alphazeta)

      end subroutine sanity_check_iota

      !**********************************************************************
      !                               Test ∇ζ                               !
      !**********************************************************************
      ! We have the following relations
      !     X = R*cos(ζ)      ζ = arccos(X/R)
      !     Y = R*sin(ζ)      ζ = arcsin(Y/R)
      !     R^2 = X^2 + Y^2
      ! Therefore, we should have
      !     dζ/dX = d(arccos(X/R))/dX = -Y/R^2 = -sin(ζ)/R  
      !     dζ/dY = d(arcsin(Y/R))/dX = X/R^2 = cos(ζ)/R  
      !     dζ/dZ = 0
      !**********************************************************************

      subroutine sanity_check_grad_zeta()

         implicit none

         real, dimension(:, :), allocatable :: minus_sinzeta_over_R, plus_coszeta_over_R
         integer :: izeta

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::sanity_check_grad_zeta'
         
         ! Allocate temporary arrays 
         allocate (minus_sinzeta_over_R(nalpha, -nzgrid:nzgrid)) 
         allocate (plus_coszeta_over_R(nalpha, -nzgrid:nzgrid)) 

         ! Calculate -sin(ζ)/R and cos(ζ)/R  
         do izeta = -nzgrid, nzgrid
            minus_sinzeta_over_R(:, izeta) = -sin(zeta(izeta)) / R(:, izeta)
            plus_coszeta_over_R(:, izeta) = cos(zeta(izeta)) / R(:, izeta)
         end do

         ! Compare dζ/dX with -sin(ζ)/R; compare dζ/dY with cos(ζ)/R; dζ/dZ should be 0 
         call check_that_arrays_match(grad_zeta_X, minus_sinzeta_over_R, 1.0e-2, 'grad_zeta_X')
         call check_that_arrays_match(grad_zeta_Y, plus_coszeta_over_R, 1.0e-2, 'grad_zeta_Y')
         call check_that_array_is_zero(grad_zeta_Z, 1.0e-14, 'grad_zeta_Z')


         ! We might as well use the exact values  
         grad_zeta_X = minus_sinzeta_over_R 
         grad_zeta_Y = plus_coszeta_over_R  
         grad_zeta_Z = 0

         ! Deallocate temporary arrays
         deallocate (minus_sinzeta_over_R)
         deallocate (plus_coszeta_over_R)

      end subroutine sanity_check_grad_zeta


      !**********************************************************************
      !                            Test Jacobian                            !
      !**********************************************************************
      ! sqrt(g) = 1 / (∇s . ∇θ x ∇ζ) 
      !         = dX/ds dY/dθ dZ/dζ + dY/ds dZ/dθ dX/dζ + dZ/ds dX/dθ dY/dζ
      !         - dX/ds dZ/dθ dY/dζ - dY/ds dX/dθ dZ/dζ - dZ/ds dY/dθ dX/dζ
      ! 1/sqrt(g) = (∇s . ∇θ x ∇ζ) 
      !           = (∇s.e_X) (∇θ.e_Y) (∇ζ.e_Z) + (∇s.e_Y) (∇θ.e_Z) (∇ζ.e_X) + (∇s.e_Z) (∇θ.e_X) (∇ζ.e_Y) 
      !           - (∇s.e_X) (∇θ.e_Z) (∇ζ.e_Y) + (∇s.e_Y) (∇θ.e_X) (∇ζ.e_Z) + (∇s.e_Z) (∇θ.e_Y) (∇ζ.e_X) 
      !**********************************************************************
      subroutine sanity_check_jacobian()

         implicit none

         real, dimension(:, :), allocatable :: temp_vs_alphazeta 

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::sanity_check_jacobian'

         ! Allocate temporary arrays 
         allocate (temp_vs_alphazeta(nalpha, -nzgrid:nzgrid))

         ! Compare sqrt(g) and 1 / (∇s . ∇θ x ∇ζ) 
         temp_vs_alphazeta = d_X_d_s * d_Y_d_theta * d_Z_d_zeta + d_Y_d_s * d_Z_d_theta * d_X_d_zeta &
                            + d_Z_d_s * d_X_d_theta * d_Y_d_zeta - d_Z_d_s * d_Y_d_theta * d_X_d_zeta &
                            - d_X_d_s * d_Z_d_theta * d_Y_d_zeta - d_Y_d_s * d_X_d_theta * d_Z_d_zeta
         call check_that_arrays_match(sqrt_g, temp_vs_alphazeta, 3.0e-3, 'sqrt_g')

         ! Compare 1/sqrt(g) and (∇s . ∇θ x ∇ζ) 
         temp_vs_alphazeta = grad_s_X * grad_theta_Y * grad_zeta_Z + grad_s_Y * grad_theta_Z * grad_zeta_X &
                            + grad_s_Z * grad_theta_X * grad_zeta_Y - grad_s_Z * grad_theta_Y * grad_zeta_X &
                            - grad_s_X * grad_theta_Z * grad_zeta_Y - grad_s_Y * grad_theta_X * grad_zeta_Z
         call check_that_arrays_match(1/sqrt_g, temp_vs_alphazeta, 1.0e-2, '1/sqrt_g')

         ! Deallocate temporary arrays
         deallocate (temp_vs_alphazeta)

      end subroutine sanity_check_jacobian

      !**********************************************************************
      !           Test covariant and contravariant componets of B           !
      !**********************************************************************
      ! We have the following dual basis
      !     e_s = (dX/ds) e_X + (dY/ds) e_Y + (dZ/ds) e_Z
      !     e_ζ = (dX/dζ) e_X + (dY/dζ) e_Y + (dZ/dζ) e_Z
      !     e_θ = (dX/dθ) e_X + (dY/dθ) e_Y + (dZ/dθ) e_Z
      ! Therefore, the covariant components of the magnetic field are 
      !     B_sub_s = B . e_s = B_X * dX/ds + B_Y * dY/ds + B_Z * dZ/ds
      !     B_sub_zeta = B . e_ζ = B_X * dX/dzeta + B_Y * dY/dzeta + B_Z * dZ/dzeta
      !     B_sub_theta = B . e_θ = B_X * dX/dtheta + B_Y * dY/dtheta + B_Z * dZ/dtheta
      ! And the contravariant components are
      !     B_sup_s = B . ∇s = (B.e_X)(∇s.e_X) + (B.e_Y)(∇s.e_Y) + (B.e_Z)(∇s.e_Z) 
      !     B_sup_zeta = B . ∇ζ = (B.e_X)(∇ζ.e_X) + (B.e_Y)(∇ζ.e_Y) + (B.e_Z)(∇ζ.e_Z) 
      !     B_sup_theta = B . ∇θ = (B.e_X)(∇θ.e_X) + (B.e_Y)(∇θ.e_Y) + (B.e_Z)(∇θ.e_Z) 
      !**********************************************************************
      subroutine sanity_check_Bcomponents()

         implicit none

         real, dimension(:, :), allocatable :: temp_vs_alphazeta 

         !----------------------------------------------------------------------  

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::sanity_check_Bcomponents'

         ! Allocate temporary arrays 
         allocate (temp_vs_alphazeta(nalpha, -nzgrid:nzgrid))

         ! B_sub_s = B . e_s = B_X * dX/ds B_Y * dY/ds B_Z * dZ/ds
         temp_vs_alphazeta = B_X * d_X_d_s + B_Y * d_Y_d_s + B_Z * d_Z_d_s
         call check_that_arrays_match(temp_vs_alphazeta, B_sub_s, 1.0e-2, 'B_sub_s')

         ! B_sub_zeta = B . e_ζ = B_X * dX/dzeta + B_Y * dY/dzeta + B_Z * dZ/dzeta
         temp_vs_alphazeta = B_X * d_X_d_zeta + B_Y * d_Y_d_zeta + B_Z * d_Z_d_zeta
         call check_that_arrays_match(temp_vs_alphazeta, B_sub_zeta, 1.0e-2, 'B_sub_zeta')

         ! B_sub_theta = B . e_θ = B_X * dX/dtheta + B_Y * dY/dtheta + B_Z * dZ/dtheta
         temp_vs_alphazeta = B_X * d_X_d_theta + B_Y * d_Y_d_theta + B_Z * d_Z_d_theta
         call check_that_arrays_match(temp_vs_alphazeta, B_sub_theta, 1.0e-2, 'B_sub_theta')

         ! B_sup_s = B . ∇s = (B.e_X)(∇s.e_X) + (B.e_Y)(∇s.e_Y) + (B.e_Z)(∇s.e_Z) = 0
         temp_vs_alphazeta = B_X * grad_s_X + B_Y * grad_s_Y + B_Z * grad_s_Z
         call check_that_array_is_zero(temp_vs_alphazeta, 1.0e-2, 'B_sup_s')
         
         ! B_sup_zeta = B . ∇ζ = (B.e_X)(∇ζ.e_X) + (B.e_Y)(∇ζ.e_Y) + (B.e_Z)(∇ζ.e_Z) 
         temp_vs_alphazeta = B_X * grad_zeta_X + B_Y * grad_zeta_Y + B_Z * grad_zeta_Z
         call check_that_arrays_match(temp_vs_alphazeta, B_sup_zeta, 1.0e-2, 'B_sup_zeta')
      
         ! B_sup_theta = B . ∇θ = (B.e_X)(∇θ.e_X) + (B.e_Y)(∇θ.e_Y) + (B.e_Z)(∇θ.e_Z) 
         temp_vs_alphazeta = B_X * grad_theta_X + B_Y * grad_theta_Y + B_Z * grad_theta_Z
         call check_that_arrays_match(temp_vs_alphazeta, B_sup_theta, 1.0e-2, 'B_sup_theta')

         ! Deallocate temporary arrays
         deallocate (temp_vs_alphazeta)

      end subroutine sanity_check_Bcomponents

      !*************************************************************************
      !                      Check that two arrays match                       !
      !*************************************************************************
      ! Verify that array1 = array2 to within a relative tolerance.
      ! Verify that |array1| = 0 to within a relative tolerance.
      !*************************************************************************

      subroutine check_that_arrays_match(array1, array2, tolerance, name)

         implicit none

         real, intent(in) :: tolerance
         character(len=*), intent(in) :: name
         real, dimension(nalpha, -nzgrid:nzgrid), intent(in) :: array1, array2

         integer :: ierr_local = 0
         real :: max_difference

         !---------------------------------------------------------------------- 
 
         ! Check the difference between the two arrays
         max_difference = maxval(abs(array1 - array2)) / maxval(abs(array1) + abs(array2))

         ! If the relative difference is bigger than the <tolarance>, give an error
         if (max_difference > tolerance) then
            ierr_local = 1
            write(*,*) "Error! Two methods for computing ", trim(name), " disagree."
            write(*,*) "Here comes method 1:"
            do ialpha = 1, nalpha
               write(*,*)array1(ialpha, :)
            end do
            write(*,*) "Here comes method 2:"
            do ialpha = 1, nalpha
               write(*,*)array2(ialpha, :)
            end do
            write(*,*) "Here comes the difference:"
            do ialpha = 1, nalpha
               write(*,*) array1(ialpha, :) - array2(ialpha, :)
            end do
         end if 

         ! Increase the total error if this check failed
         if (ierr_local /= 0) ierr = ierr + 1 

      end subroutine check_that_arrays_match

      subroutine check_that_array_is_zero(array1, tolerance, name)

         implicit none

         real, intent(in) :: tolerance
         character(len=*), intent(in) :: name
         real, dimension(nalpha, -nzgrid:nzgrid), intent(in) :: array1

         integer :: ierr_local = 0
         real :: max_value

         !---------------------------------------------------------------------- 
   
         ! Get the maximum value inside |array1|
         max_value = maxval(abs(array1))

         ! Check that the array is zero
         if (max_value > tolerance) then
            ierr_local = 1
            write(*,*) "Error! ", trim(name), " should be 0, but instead it is:"
            do ialpha = 1, nalpha
               write(*,*)array1(ialpha, :)
            end do
         end if

         ! Increase the total error if this check failed
         if (ierr_local /= 0) ierr = ierr + 1 

      end subroutine check_that_array_is_zero


      !*************************************************************************
      !                       Initialize geometry arrays                       !
      !*************************************************************************
      subroutine allocate_geometry_arrays()

         implicit none 

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::allocate_geometry_arrays'

         ! Set the already allocated arrays to zero (these enter the <calculate_vmec_geometry> routine) 
         gds23_psitalpha = 0.0; gds24_psitalpha = 0.0; gds25_psitalpha = 0.0; gds26_psitalpha = 0.0
         bmag = 0; b_dot_grad_zeta = 0.0; B_sub_theta = 0; B_sub_zeta = 0.0
         gbdrift_alpha = 0.0; gbdrift0_psit = 0.0; cvdrift_alpha = 0; cvdrift0_psit = 0.0
         grad_alpha_grad_alpha = 0.0; grad_alpha_grad_psit = 0.0; grad_psit_grad_psit = 0.0
         gradzeta_gradpsit_R2overB2 = 0.0; gradzeta_gradpsit_R2overB2 = 0.0
         
         ! Allocate the arrays versus (nalpha, nz) along the chosen field line
         allocate (B(nalpha, -nzgrid:nzgrid)); B = 0.0
         allocate (sqrt_g(nalpha, -nzgrid:nzgrid)); sqrt_g = 0.0
         allocate (R(nalpha, -nzgrid:nzgrid)); R = 0.0
         allocate (Z(nalpha, -nzgrid:nzgrid)); Z = 0.0
         allocate (d_B_d_theta(nalpha, -nzgrid:nzgrid)); d_B_d_theta = 0.0
         allocate (d_B_d_zeta(nalpha, -nzgrid:nzgrid)); d_B_d_zeta = 0.0
         allocate (d_B_d_s(nalpha, -nzgrid:nzgrid)); d_B_d_s = 0.0
         allocate (d_R_d_theta(nalpha, -nzgrid:nzgrid)); d_R_d_theta = 0.0
         allocate (d_R_d_zeta(nalpha, -nzgrid:nzgrid)); d_R_d_zeta = 0.0
         allocate (d_R_d_s(nalpha, -nzgrid:nzgrid)); d_R_d_s = 0.0
         allocate (d_Z_d_theta(nalpha, -nzgrid:nzgrid)); d_Z_d_theta = 0.0
         allocate (d_Z_d_zeta(nalpha, -nzgrid:nzgrid)); d_Z_d_zeta = 0.0
         allocate (d_Z_d_s(nalpha, -nzgrid:nzgrid)); d_Z_d_s = 0.0
         allocate (d_Lambda_d_theta(nalpha, -nzgrid:nzgrid)); d_Lambda_d_theta = 0.0
         allocate (d_Lambda_d_zeta(nalpha, -nzgrid:nzgrid)); d_Lambda_d_zeta = 0.0
         allocate (d_Lambda_d_s(nalpha, -nzgrid:nzgrid)); d_Lambda_d_s = 0.0
         allocate (B_sub_s(nalpha, -nzgrid:nzgrid)); B_sub_s = 0.0
         allocate (B_sup_theta(nalpha, -nzgrid:nzgrid)); B_sup_theta = 0.0
         allocate (B_sup_zeta(nalpha, -nzgrid:nzgrid)); B_sup_zeta = 0.0         
         allocate (B_cross_grad_B_dot_grad_alpha(nalpha, -nzgrid:nzgrid)); B_cross_grad_B_dot_grad_alpha = 0.0        
         allocate (B_cross_grad_s_dot_grad_alpha(nalpha, -nzgrid:nzgrid)); B_cross_grad_s_dot_grad_alpha = 0.0        
         allocate (B_cross_grad_B_dot_grad_psit(nalpha, -nzgrid:nzgrid)) ; B_cross_grad_B_dot_grad_psit = 0.0        
         allocate (gradzeta_gradalpha(nalpha, -nzgrid:nzgrid)); gradzeta_gradalpha = 0.0        
         allocate (gradzeta_gradpsit(nalpha, -nzgrid:nzgrid)); gradzeta_gradpsit = 0.0        
         allocate (gradthetap_gradalpha(nalpha, -nzgrid:nzgrid)); gradthetap_gradalpha = 0.0        
         allocate (gradthetap_gradpsit(nalpha, -nzgrid:nzgrid)); gradthetap_gradpsit = 0.0                   
          
         ! Allocate the arrays versus (nalpha, nz) in Cartesian coordinates
         allocate (d_X_d_s(nalpha, -nzgrid:nzgrid)); d_X_d_s = 0.0
         allocate (d_X_d_theta(nalpha, -nzgrid:nzgrid)); d_X_d_theta = 0.0
         allocate (d_X_d_zeta(nalpha, -nzgrid:nzgrid)); d_X_d_zeta = 0.0
         allocate (d_Y_d_s(nalpha, -nzgrid:nzgrid)); d_Y_d_s = 0.0
         allocate (d_Y_d_theta(nalpha, -nzgrid:nzgrid)); d_Y_d_theta = 0.0
         allocate (d_Y_d_zeta(nalpha, -nzgrid:nzgrid)); d_Y_d_zeta = 0.0 
         allocate (grad_s_X(nalpha, -nzgrid:nzgrid)); grad_s_X = 0.0
         allocate (grad_s_Y(nalpha, -nzgrid:nzgrid)); grad_s_Y = 0.0
         allocate (grad_s_Z(nalpha, -nzgrid:nzgrid)); grad_s_Z = 0.0
         allocate (grad_theta_X(nalpha, -nzgrid:nzgrid)); grad_theta_X = 0.0
         allocate (grad_theta_Y(nalpha, -nzgrid:nzgrid)); grad_theta_Y = 0.0
         allocate (grad_theta_Z(nalpha, -nzgrid:nzgrid)); grad_theta_Z = 0.0
         allocate (grad_theta_pest_X(nalpha, -nzgrid:nzgrid)); grad_theta_pest_X = 0.0
         allocate (grad_theta_pest_Y(nalpha, -nzgrid:nzgrid)); grad_theta_pest_Y = 0.0
         allocate (grad_theta_pest_Z(nalpha, -nzgrid:nzgrid)); grad_theta_pest_Z = 0.0
         allocate (grad_zeta_X(nalpha, -nzgrid:nzgrid)); grad_zeta_X = 0.0
         allocate (grad_zeta_Y(nalpha, -nzgrid:nzgrid)); grad_zeta_Y = 0.0
         allocate (grad_zeta_Z(nalpha, -nzgrid:nzgrid)); grad_zeta_Z = 0.0
         allocate (grad_psit_X(nalpha, -nzgrid:nzgrid)); grad_psit_X = 0.0
         allocate (grad_psit_Y(nalpha, -nzgrid:nzgrid)); grad_psit_Y = 0.0
         allocate (grad_psit_Z(nalpha, -nzgrid:nzgrid)); grad_psit_Z = 0.0
         allocate (grad_alpha_X(nalpha, -nzgrid:nzgrid)); grad_alpha_X = 0.0
         allocate (grad_alpha_Y(nalpha, -nzgrid:nzgrid)); grad_alpha_Y = 0.0
         allocate (grad_alpha_Z(nalpha, -nzgrid:nzgrid)); grad_alpha_Z = 0.0
         allocate (B_X(nalpha, -nzgrid:nzgrid)); B_X = 0.0
         allocate (B_Y(nalpha, -nzgrid:nzgrid)); B_Y = 0.0
         allocate (B_Z(nalpha, -nzgrid:nzgrid)); B_Z = 0.0
         allocate (grad_B_X(nalpha, -nzgrid:nzgrid)); grad_B_X = 0.0
         allocate (grad_B_Y(nalpha, -nzgrid:nzgrid)); grad_B_Y = 0.0
         allocate (grad_B_Z(nalpha, -nzgrid:nzgrid)); grad_B_Z = 0.0

      end subroutine 

      !*************************************************************************
      !                       Deallocate geometry arrays                       !
      !*************************************************************************
      subroutine deallocate_geometry_arrays()

         implicit none 

         !---------------------------------------------------------------------- 

         ! Track the code  
         if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::deallocate_geometry_arrays'

         deallocate (R, Z, B, sqrt_g) 
         deallocate (d_X_d_s, d_X_d_zeta, d_X_d_theta) 
         deallocate (d_Y_d_s, d_Y_d_zeta, d_Y_d_theta) 
         deallocate (d_R_d_s, d_R_d_zeta, d_R_d_theta) 
         deallocate (d_Z_d_s, d_Z_d_zeta, d_Z_d_theta)
         deallocate (d_B_d_s, d_B_d_zeta, d_B_d_theta)  
         deallocate (d_Lambda_d_s, d_Lambda_d_zeta, d_Lambda_d_theta) 
         deallocate (B_sub_s, B_sup_zeta, B_sup_theta)  
         deallocate (grad_s_X, grad_s_Y, grad_s_Z)
         deallocate (grad_theta_X, grad_theta_Y, grad_theta_Z)
         deallocate (grad_theta_pest_X, grad_theta_pest_Y, grad_theta_pest_Z)
         deallocate (grad_zeta_X, grad_zeta_Y, grad_zeta_Z)
         deallocate (grad_psit_X, grad_psit_Y, grad_psit_Z)
         deallocate (grad_alpha_X, grad_alpha_Y, grad_alpha_Z)
         deallocate (B_X, B_Y, B_Z)
         deallocate (grad_B_X, grad_B_Y, grad_B_Z)
         deallocate (B_cross_grad_B_dot_grad_alpha)
         deallocate (B_cross_grad_s_dot_grad_alpha)
         deallocate (B_cross_grad_B_dot_grad_psit)
         deallocate (gradzeta_gradalpha, gradzeta_gradpsit)  
         deallocate (gradthetap_gradalpha, gradthetap_gradpsit)

         if (allocated(rmnc)) then
            deallocate (xm, xn, xm_nyq, xn_nyq)
            deallocate (rmnc, lmns, zmns, bmnc, gmnc)
            deallocate (bsupumnc, bsupvmnc, bsubumnc, bsubvmnc, bsubsmns)
            deallocate (phi, phip, iotas, iotaf, presf)
         end if
         if (allocated(rmns)) then
            deallocate (rmns, lmnc, zmnc, bmns, gmns)
            deallocate (bsupumns, bsupvmns, bsubumns, bsubvmns, bsubsmnc)
         end if

      end subroutine deallocate_geometry_arrays

   end subroutine calculate_vmec_geometry

!###############################################################################
!################################ CALCULATIONS #################################
!###############################################################################

   !*************************************************************************
   !                INTERPOLATE QUANTITIES ON THE RADIAL GRID               !
   !*************************************************************************

   subroutine radial_interpolation(quantity, interpolated_quantity, grid)

      implicit none 
   
      character(4), intent(in) :: grid
      real, dimension(:), intent(in) :: quantity 
      real, intent(out) :: interpolated_quantity 

      !---------------------------------------------------------------------- 

      if (grid=='half') then
         interpolated_quantity = quantity(index_half(1)) * weight_half(1) + quantity(index_half(2)) * weight_half(2)
      else if (grid=='full') then 
         interpolated_quantity = quantity(index_full(1)) * weight_full(1) + quantity(index_full(2)) * weight_full(2)
      else 
         write(*,*) 'Warning: "grid" needs to be specified!'
      end if

   end subroutine radial_interpolation

!###############################################################################
!########################## CALCULATE THE THETA GRID ###########################
!###############################################################################
! We constructed the <alpha> grid based on the input variable <alpha0>, which 
! tells us which field lines we have selected. And we contructed the <zeta> grid
! so that the cylindrical toroidal angle <zeta> covers <number_of_field_periods>,
! centered around <zeta_center>. From the definition of <alpha> we know that 
!     <theta_pest> = alpha + iota * zeta = straight-field-line angle
! However, we need to know the cylindrical poloidal angle <theta> which 
! is used in VMEC, and from which we constructed <theta_pest> through, 
!     <theta_pest> = theta + Lambda 
! Here <Lambda> is given by the Sine/Cosine Fourier transformation of <lmns>,  
! where we have <imn> number of <xm> and <xn> Fourier coefficients, 
!     angle(imn) = xm(imn)*theta + xn(imn)*zeta
!     lambda(s, theta, zeta) = sum_{imn} lmns(imn, s) * sin(angle(imn))
! From the constructed <alpha> and <zeta> grids we thus now the <theta_pest> grid, 
! and we use a root solver, to try different values of <theta>, and see
! whether the corresponding lambda(s, theta, zeta) is such that 
!     <theta_pest> = theta + Lambda 
! Note that these routines have access to the module variables, but not
! to the <calculate_vmec_geometry()> subroutine variables.
!     weight_full, weight_half, index_full, index_half, lasym, nfp, 
!     isigng, ns, mnmax, mnmax_nyq, mpol, ntor, Aminor,xm, xn, xm_nyq, xn_nyq, 
!     rmnc, rmns, lmnc, lmns, zmnc, zmns, bmnc, bmns, gmnc, gmns,
!     bsupumnc, bsupumns, bsupvmnc, bsupvmns, bsubumnc, bsubumns, bsubvmnc, 
!     bsubvmns, bsubsmnc, bsubsmns, phi, phip, iotas, iotaf, presf
!###############################################################################

   !============================================================================
   !=================== FOR EACH (ALPHA,ZETA) FIND theta ==================
   !============================================================================
   subroutine calculate_theta(nzgrid, zeta, nalpha, alpha, iota, theta, ierr)

      use mp, only: mp_abort

      implicit none 

      real, intent(in) :: iota
      integer, intent(in) :: nzgrid, nalpha
      real, dimension(:), intent(in) :: alpha                     ! zeta(nalpha)
      real, dimension(-nzgrid:), intent(in) :: zeta               ! zeta(nzgrid)
      real, dimension(:, -nzgrid:), intent(out) :: theta          ! theta(nalpha, nzgrid)
      integer, intent(inout) :: ierr

      integer :: izeta, ialpha
      logical :: theta_converged
      real :: theta_pest_target
      real :: theta_min
      real :: theta_max
      real :: zeta0

      !------------------------------------------------------------------------- 

      ! Track the code  
      if (debug) write (*, *) 'geometry_vmec_read_netCDF_file::calculate_theta' 

      ! For each (alpha, zeta) we know <theta_pest> = alpha + iota * zeta
      ! and we use a root solver to find <theta> = theta_pest - Lambda
      do izeta = -nzgrid, nzgrid 

         ! Define the current zeta explicitly, since <fzero_residual()> will use it
         zeta0 = zeta(izeta)
         do ialpha = 1, nalpha
            
            ! For this (alpha, zeta) we have the following <theta_pest>
            theta_pest_target = alpha(ialpha) + iota * zeta0

            ! Guess that <theta> will be within 0.3 radians of <theta_pest>
            theta_min = theta_pest_target - 0.3
            theta_max = theta_pest_target + 0.3

            ! Use a root solver to find the <theta> for which 
            !     <theta_pest> = <theta> + Lambda(s, theta, zeta)
            call get_root(theta_min, theta_max, theta(ialpha, izeta), theta_converged)

            ! Sometimes the root solved can not find a corresponding <theta>
            if (.not. theta_converged) then
               write (*, *) "ERROR: could not find root needed to compute theta. Aborting."
               ierr = 117; call mp_abort("ERROR: could not find root needed to compute theta. Aborting.") 
            end if

         end do
      end do 

   contains

      !============================================================================
      !=================== RESIDUAL FUNCTION FOR THE ROOT SOLVER ==================
      !============================================================================ 
      ! The root solver will try to minimize the error function <fzero_residual>
      ! Since we want to find a <theta> for which,
      !     <theta_pest> = <theta> + Lambda(s, <theta>, zeta)
      ! The residual function that we want to minimize is,
      !     <fzero_residual> = <theta> + Lambda(s, <theta>, zeta) - <theta_pest> 
      ! Note that <lmns> and <lmnc> use the non-Nyquist <xm>, <xn>, and <mnmax>.
      ! and that <lmns> and <lmnc> are on the radial half grid.  
      !============================================================================ 
      function fzero_residual(theta_try) 

         implicit none

         real, intent(in) :: theta_try

         real :: lambda 
         real :: fzero_residual
         real :: angle, sinangle, cosangle
         real :: lmns_temp, lmnc_temp
         integer :: imn

         ! Initialize
         lambda = 0

         ! Calculate Lambda(s, <theta>, zeta) through the Sine/Cosine Fourier transformation of <lmns>
         do imn = 1, mnmax

            ! Calculate angle(imn) = xm(imn)*theta + xn(imn)*zeta
            angle = xm(imn) * theta_try - xn(imn) * zeta0
            sinangle = sin(angle); cosangle = cos(angle)

            ! Get <lmns> and <lmnc> on the chosen flux surface
            call radial_interpolation(lmns(imn,:), lmns_temp, 'half') 
            if (lasym) call radial_interpolation(lmnc(imn,:), lmnc_temp, 'half') 

            ! Add the <imn> component of Lambda to <fzero_residual>
            lambda = lambda + lmns_temp*sinangle
            if (lasym) lambda = lambda + lmnc_temp*cosangle
            
         end do

         ! <fzero_residual> = <theta> + Lambda(s, <theta>, zeta) - <theta_pest> 
         fzero_residual = theta_try + lambda - theta_pest_target

      end function fzero_residual

      !============================================================================
      !============ USE A ROOT SOLVER TO FIND THE theta THAT WORKS ===========
      !============================================================================ 

      subroutine get_root(a0, b0, root, converged)

         implicit none

         real, intent(in) :: a0, b0
         real, intent(out) :: root
         logical, intent(out) :: converged

         integer, parameter :: itmax_bracket = 10
         integer, parameter :: itmax_root = 10
         real, parameter :: tol = 1.0e-10
         integer :: it
         real :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm, eps

         a = a0
         b = b0
         fa = fzero_residual(a)
         fb = fzero_residual(b)
         do it = 1, itmax_bracket
            eps = epsilon(a)
            if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
               write (*, *)
               write (*, *) 'in calculate_vmec_geometry, theta_min=', a, ' and theta_max=', b, ' do not bracket root.'
               write (*, *) 'f(theta_min)=', fa, 'and f(theta_max)=', fb, '.'
               a = a - 0.3
               b = b + 0.3
               write (*, *) 'Trying again with values ', a, ' and ', b, ' .'
               fa = fzero_residual(a)
               fb = fzero_residual(b)
            else
               exit
            end if
         end do

         c = b
         fc = fb
         do it = 1, itmax_root
            if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
               c = a
               fc = fa
               d = b - a
               e = d
            end if
            if (abs(fc) < abs(fb)) then
               a = b
               b = c
               c = a
               fa = fb
               fb = fc
               fc = fa
            end if
            tol1 = 2.0 * eps * abs(b) + 0.5 * tol
            xm = 0.5 * (c - b)
            if (abs(xm) <= tol1 .or. fb == 0.0) then
               root = b
               converged = .true.
               exit
            end if
            if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
               s = fb / fa
               if (a == c) then
                  p = 2.0 * xm * s
                  q = 1.0 - s
               else
                  q = fa / fc
                  r = fb / fc
                  p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                  q = (q - 1.0) * (r - 1.0) * (s - 1.0)
               end if
               if (p > 0.0) q = -q
               p = abs(p)
               if (2.0 * p < min(3.0 * xm * q - abs(tol1 * q), abs(e * q))) then
                  e = d
                  d = p / q
               else
                  d = xm
                  e = d
               end if
            else
               d = xm
               e = d
            end if
            a = b
            fa = fb
            b = b + merge(d, sign(tol1, xm), abs(d) > tol1)
            fb = fzero_residual(b)
         end do

      end subroutine get_root

   end subroutine calculate_theta

end module geometry_vmec_read_netCDF_file
