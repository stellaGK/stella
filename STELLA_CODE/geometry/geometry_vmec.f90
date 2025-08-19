
!###############################################################################
!################################ VMEC GEOMETRY ################################
!###############################################################################
! 
! Routines for calculating the geometry needed by stella, from a VMEC file.
! Inside the <geometry> module we call:
! 
! call get_vmec_geometry(&
!            nzgrid, nalpha, naky, geo_surf, grho, bmag, b_dot_grad_z, &
!            b_dot_grad_z_averaged, grad_alpha_grad_alpha, &
!            grad_alpha_grad_psit, grad_psit_grad_psit, &
!            gds23_psitalpha , gds24_psitalpha , gds25_psitalpha , gds26_psitalpha , gbdrift_alpha, gbdrift0_psi, &
!            cvdrift_alpha, cvdrift0_psi, sign_torflux, &
!            theta_vmec, dzetadz, aref, bref, alpha, zeta, &
!            field_period_ratio, psit_displacement_fac)
! 
! The VMEC module will calculate the geometric arrays with psi_t as the
! the radial coordinate, and zeta as the parallel coordinate, on a z-grid 
! with <zgrid_refinement_factor> more z-points than the real stella grid. 
! This module can change the parallel coordinate to the normalized arc-length,
! and it interpolates the VMEC z-grid to the stella z-grid.
! 
! Initial VMEC geometry code was written by Matt Landreman, University of Maryland in August 2017. 
! Modified in 2018-2019 by Michael Barnes, and cleaned in 2024 by Hanne Thienpondt.
! 
! Changes
! -------
! 07/2024 Removed <b_dot_grad_zeta_prefac> and <zgrid_scalefac> 
! 07/2024 Changed <gradpar> to <b_dot_grad_z_averaged> and other sensible name changed
! 
!###############################################################################

module vmec_geometry

   implicit none

	! Routines
   public :: read_vmec_parameters
   public :: get_vmec_geometry

   ! Radial coordinate options
   public :: radial_coordinate_option
   public :: radial_coordinate_sgnpsitpsit
   public :: radial_coordinate_minuspsit
   public :: radial_coordinate_r

   private
   integer :: radial_coordinate_option
   integer, parameter :: radial_coordinate_sgnpsitpsit = 1
   integer, parameter :: radial_coordinate_minuspsit = 2
   integer, parameter :: radial_coordinate_r = 3
 
   real :: alpha0, nfield_periods
   real :: zeta_center, torflux

   integer :: n_tolerated_test_arrays_inconsistencies
   integer :: zgrid_refinement_factor
   integer :: surface_option

   logical :: verbose, rectangular_cross_section
   character(2000) :: vmec_filename
   character(20) :: radial_coordinate 

contains

   !============================================================================
   !========================== READ VMEC PARAMETERS ============================
   !============================================================================  
   subroutine read_vmec_parameters

      use text_options, only: text_option, get_option_value
      use file_utils, only: input_unit_exist, error_unit
      use zgrid, only: zed_equal_arc
      use mp, only: mp_abort

      implicit none

      integer :: in_file, ierr
      logical :: exist
      
      ! Backwards comptability, these parameters have been removed
      real :: zgrid_scalefac

      ! Allow text options for <radial_coordinate> to choose sgn(psi_t)*psi_t; -psi_t or r
      type(text_option), dimension(5), parameter :: radial_coordinate_options = & 
                (/text_option('default', radial_coordinate_sgnpsitpsit), &
                  text_option('sgn(psi_t) psi_t', radial_coordinate_sgnpsitpsit), &
                  text_option('original', radial_coordinate_minuspsit), &
                  text_option('- psi_t', radial_coordinate_minuspsit), &
                  text_option('r', radial_coordinate_r)/)

      !---------------------------------------------------------------------- 

      ! Define the variables in the namelist
      namelist /vmec_parameters/ alpha0, zeta_center, rectangular_cross_section, nfield_periods, &
         torflux, zgrid_refinement_factor, surface_option, radial_coordinate, &
         verbose, vmec_filename, n_tolerated_test_arrays_inconsistencies, &
         ! Backwards compatibility for old stella code
         zgrid_scalefac

      ! Assign default variables
      call init_vmec_defaults

      ! Read the <vmec_parameters> namelist in the input file 
      in_file = input_unit_exist("vmec_parameters", exist)
      if (exist) read (unit=in_file, nml=vmec_parameters)

      ! Read the text option for <radial_coordinate> 
      ierr = error_unit()
      call get_option_value(radial_coordinate, radial_coordinate_options, &
                  radial_coordinate_option, ierr, "radial_coordinate in vmec_parameters")

      ! If we set <zed_equal_arc> = True, we also define <zgrid_refinement_factor>
      if (.not. zed_equal_arc) then
         if (zgrid_refinement_factor > 1) then
            write (*, *) 'There is no reason to use zgrid_refinement_factor > 1 unless zed_equal_arc=T'
            write (*, *) 'Setting zgrid_refinement_factor = 1'
            zgrid_refinement_factor = 1
         end if
      end if

      ! Allow <n_tolerated_test_arrays_inconsistencies> errors when reading the VMEC
      if (n_tolerated_test_arrays_inconsistencies < 0) then
         write (*, *) 'n_tolerated_test_arrays_inconsistencies = ', n_tolerated_test_arrays_inconsistencies
         call mp_abort('n_tolerated_test_arrays_inconsistencies should always be >= 0.  aborting')
      end if
      
      ! Backwards compatibility
      if (zgrid_scalefac>0) write(*,*) 'WARNING: The parameter <zgrid_scalefac> in the <zgrid_parameters> knob has been removed.'
   
   contains

       !=========================================================================
       !========================= DEFAULT VMEC PARAMETERS =======================
       !=========================================================================  
       subroutine init_vmec_defaults

          use zgrid, only: zed_equal_arc

          implicit none

          !----------------------------------------------------------------------- 

          ! Default parameters for VMEC
          vmec_filename = 'equilibria/wout_w7x_standardConfig.nc'
          alpha0 = 0.0
          zeta_center = 0.0
          nfield_periods = -1.0
          torflux = 0.6354167d+0
          surface_option = 0
          verbose = .true.
          n_tolerated_test_arrays_inconsistencies = 0
          zgrid_refinement_factor = 1
          radial_coordinate = 'sgn(psi_t) psi_t'

          ! If we use the normalized arc-length as the parallel coordinate, 
          ! then use <zgrid_refinement_factor> more z-points to calculate
          ! the geometry arrays in VMEC, by using a smaller step dzeta.
          if (zed_equal_arc) then
             zgrid_refinement_factor = 4
          else
             zgrid_refinement_factor = 1
          end if

          ! For alpha=0 the perpendicular cross-section is rectangular at zeta_center=0
          ! For alpha!=0 it is not in the original stella, since alpha = theta_p - iota*zeta
          ! To make the cross-section rectangular we need to define alpa = theta_p - iota(zeta-zeta_center) 
          ! This is done by toggling <rectangular_cross_section> to .true.
          rectangular_cross_section = .false.
          
          ! Backwards compatibility
          zgrid_scalefac = -1.0

      end subroutine init_vmec_defaults

   end subroutine read_vmec_parameters

   !============================================================================
   !==================== GET THE GEOMETRY VECTORS FROM VMEC ====================
   !============================================================================   
   subroutine get_vmec_geometry(nzgrid, nalpha, naky, surf, grho, bmag, &
                     b_dot_grad_z_averaged, b_dot_grad_z, &
                     grad_alpha_grad_alpha, grad_alpha_grad_psit, grad_psit_grad_psit, &
                     gds23_psitalpha , gds24_psitalpha, gds25_psitalpha, gds26_psitalpha, &
                     gbdrift_alpha, gbdrift0_psi, cvdrift_alpha, cvdrift0_psi, &
                     gradzeta_gradpsit_R2overB2, gradzeta_gradalpha_R2overB2, b_dot_grad_zeta_RR, &
                     sign_torflux, theta_vmec, dzetadz, L_reference, B_reference, alpha, zeta, &
                     field_period_ratio, psit_displacement_fac)

      use constants, only: pi
      use common_types, only: flux_surface_type
      use splines, only: geo_spline
      use physics_parameters, only: full_flux_surface
      use debug_flags, only: const_alpha_geo 
      use zgrid, only: zed_equal_arc, get_total_arc_length, get_arc_length_grid
      use zgrid, only: zed 
      use geometry_vmec_read_netCDF_file, only: calculate_vmec_geometry
      use file_utils, only: open_output_file
      use mp, only: mp_abort

      implicit none
      
      ! TODO-HT TODO-GA Circular dependency, change to use run_parameters, only: print_extra_info_to_terminal
      logical :: print_extra_info_to_terminal = .false. 

      integer, intent(in) :: nzgrid, nalpha, naky 
      integer, intent(out) :: sign_torflux
      type(flux_surface_type), intent(out) :: surf      
      real, dimension(:), intent(out) :: alpha
      real, intent(out) :: dzetadz, L_reference, B_reference, field_period_ratio
      real, dimension(-nzgrid:), intent(out) :: b_dot_grad_z_averaged
      real, dimension(:, -nzgrid:), intent(out) :: grho, bmag, b_dot_grad_z, &
               grad_alpha_grad_alpha, grad_alpha_grad_psit, grad_psit_grad_psit, &
               gds23_psitalpha , gds24_psitalpha , gds25_psitalpha , gds26_psitalpha , gbdrift_alpha, gbdrift0_psi, &
               cvdrift_alpha, cvdrift0_psi, theta_vmec, zeta, psit_displacement_fac, &
               gradzeta_gradpsit_R2overB2, gradzeta_gradalpha_R2overB2, b_dot_grad_zeta_RR


      ! These routines are always called on the first processor only
      logical, parameter :: debug = .false.

      integer :: ierr
      integer :: tmpunit
      integer :: i, j, ia, iz
      integer :: nzgrid_vmec
      integer :: zetamax_idx 

      real :: dzeta_vmec, zmin, nfp
      real, dimension(nalpha, -nzgrid:nzgrid) :: theta
      real, dimension(nalpha, -nzgrid:nzgrid) :: B_sub_theta_vmec, B_sub_zeta 
      real, dimension(:), allocatable :: zeta_vmec, zed_domain_size
      real, dimension(:, :), allocatable :: thetamod_vmec, arc_length 
      real, dimension(:, :), allocatable :: B_sub_theta_vmec_mod, B_sub_zeta_mod  
      real, dimension(:, :), allocatable :: bmag_vmec, b_dot_grad_zeta_vmec
      real, dimension(:, :), allocatable :: grad_alpha_grad_alpha_vmec, grad_alpha_grad_psit_vmec, grad_psit_grad_psit_vmec
      real, dimension(:, :), allocatable :: gds23_psitalpha_vmec, gds24_psitalpha_vmec, gds25_psitalpha_vmec, gds26_psitalpha_vmec
      real, dimension(:, :), allocatable :: gbdrift_alpha_vmec, gbdrift0_psi_vmec, cvdrift_alpha_vmec, cvdrift0_psi_vmec 
      real, dimension(:, :), allocatable :: psit_displacement_fac_vmec       
      real, dimension(:, :), allocatable :: gradzeta_gradpsit_R2overB2_vmec
      real, dimension(:, :), allocatable :: gradzeta_gradalpha_R2overB2_vmec
      real, dimension(:, :), allocatable :: b_dot_grad_zeta_RR_vmec
      real, dimension(:, :), allocatable :: b_dot_grad_zeta, b_dot_grad_arclength
      real, dimension(:), allocatable :: b_dot_grad_zeta_averaged, b_dot_grad_arclength_averaged 


      !---------------------------------------------------------------------- 

      ! If desired, increase the number of sampled zeta grid points in VMEC data to increase
      ! the accuracy of later integration in zeta and interpolation onto stella zed grid 
      nzgrid_vmec = nzgrid * zgrid_refinement_factor 

      ! Allocate VMEC geometry arrays of size 2*<nzgrid_vmec>+1
      ! The '_vmec' indicated that the arrays are defined on the extended z-grid
      if (debug) write (*, *) 'get_vmec_geometry::allocate_arrays'
      allocate (zeta_vmec(-nzgrid_vmec:nzgrid_vmec)); zeta_vmec = 0.0
      allocate (thetamod_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); thetamod_vmec = 0.0
      allocate (B_sub_zeta_mod(nalpha, -nzgrid_vmec:nzgrid_vmec)); B_sub_zeta_mod = 0.0  
      allocate (B_sub_theta_vmec_mod(nalpha, -nzgrid_vmec:nzgrid_vmec)); B_sub_theta_vmec_mod = 0.0
      allocate (bmag_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); bmag_vmec = 0.0
      allocate (b_dot_grad_zeta_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); b_dot_grad_zeta_vmec = 0.0
      allocate (grad_alpha_grad_alpha_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); grad_alpha_grad_alpha_vmec = 0.0
      allocate (grad_alpha_grad_psit_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); grad_alpha_grad_psit_vmec = 0.0
      allocate (grad_psit_grad_psit_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); grad_psit_grad_psit_vmec = 0.0
      allocate (gds23_psitalpha_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); gds23_psitalpha_vmec = 0.0
      allocate (gds24_psitalpha_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); gds24_psitalpha_vmec = 0.0
      allocate (gds25_psitalpha_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); gds25_psitalpha_vmec = 0.0
      allocate (gds26_psitalpha_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); gds26_psitalpha_vmec = 0.0
      allocate (gbdrift_alpha_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); gbdrift_alpha_vmec = 0.0
      allocate (gbdrift0_psi_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); gbdrift0_psi_vmec = 0.0
      allocate (cvdrift_alpha_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); cvdrift_alpha_vmec = 0.0
      allocate (cvdrift0_psi_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); cvdrift0_psi_vmec = 0.0
      allocate (psit_displacement_fac_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); psit_displacement_fac_vmec = 0.0
      allocate (gradzeta_gradpsit_R2overB2_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); gradzeta_gradpsit_R2overB2_vmec = 0.0
      allocate (gradzeta_gradalpha_R2overB2_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); gradzeta_gradalpha_R2overB2_vmec = 0.0
      allocate (b_dot_grad_zeta_RR_vmec(nalpha, -nzgrid_vmec:nzgrid_vmec)); b_dot_grad_zeta_RR_vmec = 0.0
      allocate (arc_length(nalpha, -nzgrid_vmec:nzgrid_vmec)); arc_length = 0.0
      allocate (zed_domain_size(nalpha)); zed_domain_size = 0.0 
      
      ! We will define b_dot_grad_* arrays on the (-nzgrid:nzgrid) to define <b_dot_grad_z> and <b_dot_grad_z_averaged>
      allocate (b_dot_grad_arclength_averaged(-nzgrid:nzgrid)); b_dot_grad_arclength_averaged = 0.0
      allocate (b_dot_grad_arclength(nalpha, -nzgrid:nzgrid)); b_dot_grad_arclength = 0.0
      allocate (b_dot_grad_zeta_averaged(-nzgrid:nzgrid)); b_dot_grad_zeta_averaged = 0.0
      allocate (b_dot_grad_zeta(nalpha, -nzgrid:nzgrid)); b_dot_grad_zeta = 0.0

      ! Calculate the geometry arrays from the VMEC file using geometry_vmec_read_netCDF_file.f90
      ! Some quantities will be assigned to the module variables such as <nfp>
      if (debug) write (*, *) 'get_vmec_geometry::calculate_vmec_geometry'
      call calculate_vmec_geometry(&
               ! Input parameters
               vmec_filename, nalpha, alpha0, nzgrid_vmec, zeta_center, &
               rectangular_cross_section, nfield_periods, torflux, surface_option, verbose, &
               ! Output parameters
               surf%rhoc, surf%qinp, surf%shat, L_reference, B_reference, nfp, &
               sign_torflux, alpha, zeta_vmec, bmag_vmec, b_dot_grad_zeta_vmec, grad_alpha_grad_alpha_vmec, &
               grad_alpha_grad_psit_vmec, grad_psit_grad_psit_vmec, gds23_psitalpha_vmec, gds24_psitalpha_vmec, &
               gds25_psitalpha_vmec, gds26_psitalpha_vmec, gbdrift_alpha_vmec, gbdrift0_psi_vmec, cvdrift_alpha_vmec, &
               cvdrift0_psi_vmec, thetamod_vmec, B_sub_zeta_mod, B_sub_theta_vmec_mod, psit_displacement_fac_vmec, &
               gradzeta_gradpsit_R2overB2_vmec, gradzeta_gradalpha_R2overB2_vmec, &
               b_dot_grad_zeta_RR_vmec, ierr)

      ! Stop stella if we had too many errors when calculating the geometry arrays (TODO: I think this is broken now)
      if (ierr /= 0) then
         if (ierr > n_tolerated_test_arrays_inconsistencies .or. ierr < 0) then
            call mp_abort('calculate_vmec_geometry returned error.')
         end if
      end if

      ! Get ratio of number of simulated field periods to the number of field periods of the device
      field_period_ratio = nfield_periods / real(nfp) 

      ! Step size in the <zeta_vmec> grid
      dzeta_vmec = zeta_vmec(1) - zeta_vmec(0)

      ! Interpolate geometric quantities from (zeta,alpha) grid to (zed,alpha) grid, 
      ! with zed = z = l the normalised arc-length. First we need to get zed(zeta,alpha)
      ! We rewrite b . ∇ζ as b . ∇l with l the normalised arc-length
      !     <b_dot_grad_zeta_vmec> = b . ∇ζ = b . ∇l (dζ/dl)
      !     int 1 / (b . ∇ζ) dζ = int 1 / (b . ∇l) dl
      if (zed_equal_arc) then
      
			! If we choose z = arc length then the geometric coefficients from VMEC need to be interpolated to the arc length grid
         if (debug) write (*, *) 'get_vmec_geometry::zed_equal_arc'

         ! Index for the max zeta of the nominal zeta grid ranging from [-zetamax_idx, zetamax_idx]
         !     <zetamax_idx> = <nzgrid> * <zgrid_refinement_factor> = <nzgrid_vmec>
         zetamax_idx = nzgrid_vmec

         ! For each field line (or each alpha), we calculate the total arc length and <b_dot_grad_z>
         do ia = 1, nalpha

            ! We can calculate the total arc-length (not the normalized one!) as,
            !     <zed_domain_size> = total arc-length = int 1 / (b . ∇ζ) dζ
            call get_total_arc_length(zetamax_idx, b_dot_grad_zeta_vmec(ia, -zetamax_idx:zetamax_idx), dzeta_vmec, zed_domain_size(ia))

            ! We create the arc-length vector which ranges from
            !     <arc_length> = [-total arc-length/2, +total arc-length/2] = [-<zed_domain_size>/2, <zed_domain_size>/2]
            zmin = -zed_domain_size(ia) * 0.5
            call get_arc_length_grid(zetamax_idx, nzgrid_vmec, zmin, b_dot_grad_zeta_vmec(ia, :), dzeta_vmec, arc_length(ia, :))

            ! We know that b . ∇l is a constant along l, so write A = b . ∇l 
            !     b . ∇ζ = b . ∇l (dζ/dl) = A (dζ/dl)
            !     int 1 / (b . ∇ζ) dζ = int 1 / (b . ∇l) dl = (1/A) int dl
            !     <zed_domain_size> = total arc-length = int 1 / (b . ∇ζ) dζ
            !     A = (1/<zed_domain_size>) * int dl 
            ! At this point we assume that l is the the normalised arc-length, so that z = l lies in the range [-pi:pi]
            !     int dl = 2 * pi 
            ! Therefore, we can calculate b . ∇l as, 
            !     <b_dot_grad_arclength> = b . ∇l = A = int dl / <zed_domain_size> = 2 * pi / <zed_domain_size>
            ! Note that <b_dot_grad_arclength> is calculated on the stella grid, which has <zgrid_refinement_factor> less points than the VMEC grid
            b_dot_grad_arclength(ia, :) = 2.0 * pi / zed_domain_size(ia)

         end do

         ! Rescale z = arc_length so that it is arc-length compressed (or expanded) to the range [-pi:pi] 
         ! Note that <arc_length> has dimensions (nalpha, -nzgrid:nzgrid) and not (nalpha, -nzgrid_vmec:nzgrid_vmec))
         arc_length = 2.0 * pi * arc_length / spread(zed_domain_size, 2, 2 * nzgrid_vmec + 1)

         ! Now that we have normalized arc-length(alpha,zeta), interpolate from the VMEC zeta grid (which is 
         ! irregular in l) to the normalized arc-length grid (which is irregular in zeta)
			! In this case we have nzgrid_vmec > nzgrid so the resulting arrays will have less z-points
         if (debug) write (*, *) 'get_vmec_geometry::geo_spline'
         do ia = 1, nalpha
            call geo_spline(arc_length(ia, :), zeta_vmec, zed, zeta(ia, :))
            call geo_spline(arc_length(ia, :), b_dot_grad_zeta_vmec(ia, :), zed, b_dot_grad_zeta(ia, :))
            call geo_spline(arc_length(ia, :), bmag_vmec(ia, :), zed, bmag(ia, :))
            call geo_spline(arc_length(ia, :), grad_alpha_grad_alpha_vmec(ia, :), zed, grad_alpha_grad_alpha(ia, :))
            call geo_spline(arc_length(ia, :), grad_alpha_grad_psit_vmec(ia, :), zed, grad_alpha_grad_psit(ia, :))
            call geo_spline(arc_length(ia, :), grad_psit_grad_psit_vmec(ia, :), zed, grad_psit_grad_psit(ia, :))
            call geo_spline(arc_length(ia, :), gds23_psitalpha_vmec(ia, :), zed, gds23_psitalpha (ia, :))
            call geo_spline(arc_length(ia, :), gds24_psitalpha_vmec(ia, :), zed, gds24_psitalpha (ia, :))
            call geo_spline(arc_length(ia, :), gds25_psitalpha_vmec(ia, :), zed, gds25_psitalpha (ia, :))
            call geo_spline(arc_length(ia, :), gds26_psitalpha_vmec(ia, :), zed, gds26_psitalpha (ia, :))
            call geo_spline(arc_length(ia, :), gbdrift_alpha_vmec(ia, :), zed, gbdrift_alpha(ia, :))
            call geo_spline(arc_length(ia, :), gbdrift0_psi_vmec(ia, :), zed, gbdrift0_psi(ia, :))
            call geo_spline(arc_length(ia, :), cvdrift_alpha_vmec(ia, :), zed, cvdrift_alpha(ia, :))
            call geo_spline(arc_length(ia, :), cvdrift0_psi_vmec(ia, :), zed, cvdrift0_psi(ia, :))
            call geo_spline(arc_length(ia, :), thetamod_vmec(ia, :), zed, theta_vmec(ia, :))
            call geo_spline(arc_length(ia, :), B_sub_zeta_mod(ia, :), zed, B_sub_zeta(ia, :)) 
            call geo_spline(arc_length(ia, :), B_sub_theta_vmec_mod(ia, :), zed, B_sub_theta_vmec(ia, :)) 
            call geo_spline(arc_length(ia, :), psit_displacement_fac_vmec(ia, :), zed, psit_displacement_fac(ia, :))            
            call geo_spline(arc_length(ia, :), gradzeta_gradpsit_R2overB2_vmec(ia, :), zed, gradzeta_gradpsit_R2overB2(ia, :))
            call geo_spline(arc_length(ia, :), gradzeta_gradalpha_R2overB2_vmec(ia, :), zed, gradzeta_gradalpha_R2overB2(ia, :))
            call geo_spline(arc_length(ia, :), b_dot_grad_zeta_RR_vmec(ia, :), zed, b_dot_grad_zeta_RR(ia, :))

            ! Here we still have that <b_dot_grad_zeta> = b . ∇ζ but we want it to be b . ∇l = b . ∇ζ * dl/dζ.
            ! we have constructed <b_dot_grad_z> = b . ∇l to be a function purely of alpha,
            ! so dl/dζ = b_dot_grad_arclength(alpha) / b_dot_grad_zeta(alpha,zeta)
            ! and ∇ζ * dl/dζ = ∇l, so multiply <gds23_psitalpha > and <gds24_psitalpha > with dl/dζ
            gds23_psitalpha (ia, :) = gds23_psitalpha (ia, :) * b_dot_grad_arclength(ia, :) / b_dot_grad_zeta(ia, :)
            gds24_psitalpha (ia, :) = gds24_psitalpha (ia, :) * b_dot_grad_arclength(ia, :) / b_dot_grad_zeta(ia, :)
            
         end do

         ! Define <b_dot_grad_arclength_averaged> to be the average value of <b_dot_grad_arclength> in alpha
         ! Note that <b_dot_grad_arclength> and <b_dot_grad_arclength_averaged> have dimensions (nalpha, -nzgrid:nzgrid) 
         b_dot_grad_arclength_averaged = sum(b_dot_grad_arclength, 1) / size(b_dot_grad_arclength, 1)

         ! We now have geometric coefficients on an alpha-grid. As we will be multiplying this 
         ! with functions of g and phi, we must take care to avoid aliasing. This is 
         ! accomplished by filtering out the highest third of the wavenumber spectra.
         if (debug) write (*, *) 'get_vmec_geometry::filter_geo_coef'
         if (full_flux_surface .and. .not. const_alpha_geo) then
            do iz = -nzgrid, nzgrid
               call filter_geo_coef(naky, bmag(:, iz))
               call filter_geo_coef(naky, grad_alpha_grad_alpha(:, iz))
               call filter_geo_coef(naky, grad_alpha_grad_psit(:, iz))
               call filter_geo_coef(naky, grad_psit_grad_psit(:, iz))
               call filter_geo_coef(naky, gds23_psitalpha (:, iz))
               call filter_geo_coef(naky, gds24_psitalpha (:, iz))
               call filter_geo_coef(naky, gds25_psitalpha (:, iz))
               call filter_geo_coef(naky, gds26_psitalpha (:, iz))
               call filter_geo_coef(naky, gbdrift_alpha(:, iz))
               call filter_geo_coef(naky, gbdrift0_psi(:, iz))
               call filter_geo_coef(naky, cvdrift_alpha(:, iz))
               call filter_geo_coef(naky, cvdrift0_psi(:, iz))
               call filter_geo_coef(naky, b_dot_grad_arclength(:, iz))
            end do
         end if
      end if 

      ! If <zed_equal_arc> = .false., then the zed coordinate is the same as VMEC's zeta coordinate,
      ! so no need to interpolate onto the stella grid, with <zgrid_refinement_factor> less points than the VMEC grid
      if (.not. zed_equal_arc) then

			! If we choose z = zeta then the geometric coefficients from VMEC are already defined on the correct z-grid
			! In this case we also have nzgrid_vmec = nzgrid so we have the correct number of z-points
         if (debug) write (*, *) 'get_vmec_geometry::not_zed_equal_arc'
         zeta = spread(zeta_vmec, 1, nalpha)
         bmag = bmag_vmec
         b_dot_grad_zeta_averaged = b_dot_grad_zeta_vmec(1, :)
         b_dot_grad_zeta = b_dot_grad_zeta_vmec
         grad_alpha_grad_alpha = grad_alpha_grad_alpha_vmec
         grad_alpha_grad_psit = grad_alpha_grad_psit_vmec
         grad_psit_grad_psit = grad_psit_grad_psit_vmec
         gds23_psitalpha = gds23_psitalpha_vmec
         gds24_psitalpha = gds24_psitalpha_vmec
         gds25_psitalpha = gds25_psitalpha_vmec
         gds26_psitalpha = gds26_psitalpha_vmec
         gbdrift_alpha = gbdrift_alpha_vmec
         gbdrift0_psi = gbdrift0_psi_vmec
         cvdrift_alpha = cvdrift_alpha_vmec
         cvdrift0_psi = cvdrift0_psi_vmec
         gradzeta_gradpsit_R2overB2 = gradzeta_gradpsit_R2overB2_vmec
         gradzeta_gradalpha_R2overB2 = gradzeta_gradalpha_R2overB2_vmec
         b_dot_grad_zeta_RR = b_dot_grad_zeta_RR_vmec
         
         ! TODO-HT Not sure what these are for
         psit_displacement_fac = psit_displacement_fac_vmec
         theta_vmec = thetamod_vmec
         
         ! TODO-HT these are not used and can be removed I think
         B_sub_theta_vmec = B_sub_theta_vmec_mod ! JFP
         B_sub_zeta = B_sub_zeta_mod ! JFP
 
         ! The parallel coordinate is zeta, but we want to use z = zeta/P so that z 
         ! is compressed (or expanded) to the range [-pi,pi], hence dzeta/dz = P = nfp/nfield_periods
         dzetadz = real(nfp) / nfield_periods  
         b_dot_grad_zeta_averaged = b_dot_grad_zeta_averaged * dzetadz
         b_dot_grad_zeta = b_dot_grad_zeta * dzetadz
         gds23_psitalpha  = gds23_psitalpha  * dzetadz
         gds24_psitalpha  = gds24_psitalpha  * dzetadz 

      end if
      
      ! Choose z-coordinate
      if (zed_equal_arc) then
         b_dot_grad_z_averaged = b_dot_grad_arclength_averaged
         b_dot_grad_z = b_dot_grad_arclength
      else if (.not.zed_equal_arc) then 
         b_dot_grad_z_averaged = b_dot_grad_zeta_averaged
         b_dot_grad_z = b_dot_grad_zeta
      end if
     
      ! The arrays over the extended zeta-grid are no longer needed, so deallocate
      deallocate (grad_alpha_grad_alpha_vmec, grad_alpha_grad_psit_vmec, grad_psit_grad_psit_vmec)
      deallocate (B_sub_theta_vmec_mod, B_sub_zeta_mod, bmag_vmec, b_dot_grad_zeta_vmec)
      deallocate (gds23_psitalpha_vmec, gds24_psitalpha_vmec, gds25_psitalpha_vmec, gds26_psitalpha_vmec)
      deallocate (b_dot_grad_arclength_averaged, b_dot_grad_arclength, b_dot_grad_zeta_averaged, b_dot_grad_zeta) 
      deallocate (zed_domain_size, zeta_vmec, thetamod_vmec)
      deallocate (gbdrift_alpha_vmec, gbdrift0_psi_vmec)
      deallocate (cvdrift_alpha_vmec, cvdrift0_psi_vmec)
      deallocate (psit_displacement_fac_vmec, arc_length) 
      deallocate (gradzeta_gradpsit_R2overB2_vmec)
      deallocate (gradzeta_gradalpha_R2overB2_vmec) 
      deallocate (b_dot_grad_zeta_RR_vmec) 

      !> calculate_vmec_geometry returns psitor/psitor_lcfs as rhoc
      !> stella uses rhoc = rho = sqrt(psitor/psitor_lcfs) = rhotor
      surf%rhoc = sqrt(surf%rhoc)
      surf%rhotor = surf%rhoc 

      ! Use rho = sqrt(psi_t / psi_{t,LCFS}) and Bref = 2 |psi_LCFS|/a^2
      ! drho/dpsi_t = 1/(2*sqrt(psi_t*psi_{t,LCFS})) * 2 |psi_LCFS|/ (a^2*Bref) = sgn(psi_t)/(rho*a^2*Bref)
      ! <grho> = a * |grad rho| = a * |drho/dpsi_t| * |grad psi_t|
      !     = (a^2*Bref) |drho/dpsi_t| * 1/(a*Bref) |grad psi_t|
      !     = (a^2*Bref) 1/(rho*a^2*Bref) * sqrt(grad_psit_grad_psit)
      !     = 1/rho * sqrt(grad_psit_grad_psit) 
      grho = sqrt(grad_psit_grad_psit) / surf%rhotor

      ! Something for the <sfincs_interface> module?
      surf%psitor_lcfs = 0.5 * sign_torflux
      surf%drhotordrho = 1.0
      
      ! Something for the <radial_variation> module?
      surf%rhoc_psi0 = surf%rhoc
      surf%qinp_psi0 = surf%qinp
      surf%shat_psi0 = surf%shat

      ! When we plot the modes on the extended z-grid from the ballooning transformation
      ! We need to move the z-domains by a factor zed0, i.e., we plot |phi|^2 versus zed(iz) - zed0(iky, ikx)
      ! with zed0(ky,kx) = theta0 * geo_surf%zed0_fac; theta0(ky,kx) = kx/(ky*shat) and theta = q*zeta
      ! The (zed/theta)-value at the end of the z-domain is
      surf%zed0_fac = -zed(nzgrid) / zeta(1, nzgrid) * surf%qinp

      ! VMEC theta is the cylindrical theta-angle (not straight-field-line coordinate) scaled between [-pi:pi]
      theta_vmec = theta_vmec / nfp

      ! The straight-field-line angle and field line label are
      !     theta_pest = theta_vmec + Lambda(psi,alpha,theta_vmec)
      !     alpha = theta_pest - iota*zeta
      ! so theta is theta_pest up to constant (alpha)  
      theta = spread(alpha, 2, 2 * nzgrid + 1) + zeta / surf%qinp

      ! Write the VMEC geometry arrays to a text file
      call open_output_file(tmpunit, '.vmec.geo')
      write (tmpunit, '(6a12)') '#rhotor', 'qinp', 'shat', 'aref', 'Bref', 'dzetadz'
      write (tmpunit, '(6e12.4)') surf%rhoc, surf%qinp, surf%shat, L_reference, B_reference, dzetadz
      write (tmpunit, '(17a12)') '#    alpha', 'zeta', 'bmag', 'b_dot_grad_z_avg', 'bdot_grad_z', 'grad_alpha2', &
         'gd_alph_psi', 'grad_psi2', 'gds23_psitalpha ', 'gds24_psitalpha ', 'gbdriftalph', 'gbdrift0psi', 'cvdriftalph', &
         'cvdrift0psi', 'theta_vmec', 'B_sub_theta', 'B_sub_zeta' 
      do j = -nzgrid, nzgrid
         do i = 1, nalpha
            write (tmpunit, '(17e12.4)') alpha(i), zeta(i, j), bmag(i, j), b_dot_grad_z_averaged(j), b_dot_grad_z(i, j), &
               grad_alpha_grad_alpha(i, j), grad_alpha_grad_psit(i, j), grad_psit_grad_psit(i, j), &
               gds23_psitalpha (i, j), gds24_psitalpha (i, j),  gbdrift_alpha(i, j), gbdrift0_psi(i, j), & 
               cvdrift_alpha(i, j), cvdrift0_psi(i, j), theta_vmec(i, j), B_sub_theta_vmec(i, j), B_sub_zeta(i, j)  
         end do
         write (tmpunit, *)
      end do
      close (tmpunit)

      ! Write some information to the command prompt
      if (verbose .and. print_extra_info_to_terminal) then 
         if (radial_coordinate_option==radial_coordinate_r) write(*,*) '  The radial coordinate psi is chosen to be psi = r = a*sqrt(psi_t/psi_{t,LCFS})'
         if (radial_coordinate_option==radial_coordinate_minuspsit) write(*,*) '  The radial coordinate psi is chosen to be psi = -psi_t'
         if (radial_coordinate_option==radial_coordinate_sgnpsitpsit) write(*,*) '  The radial coordinate psi is chosen to be psi = sgn(psi_t) psi_t'
         write(*,*) ' '
      end if

   end subroutine get_vmec_geometry

   !============================================================================
   !===================== FILTER THE GEOMETRY COEFFICIENTS =====================
   !============================================================================ 
   subroutine filter_geo_coef(naky, geocoef)

      use stella_transforms, only: transform_alpha2kalpha, transform_kalpha2alpha

      implicit none

      integer, intent(in) :: naky
      real, dimension(:), intent(in out) :: geocoef

      complex, dimension(:), allocatable :: fourier

      ! Filtering and padding are built-in to the Fourier transform routines below
      allocate (fourier(naky))
      call transform_alpha2kalpha(geocoef, fourier)
      call transform_kalpha2alpha(fourier, geocoef)
      deallocate (fourier)

   end subroutine filter_geo_coef

end module vmec_geometry
