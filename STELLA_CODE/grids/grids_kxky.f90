!###############################################################################
!                             (KX,KY) GRIDS MODULE                              
!###############################################################################
! This module initialises the kx and ky grids for stella.
! It reads the relevant namelists within namelist_kxky_grid.f90,
! and sets up the kx and ky arrays accordingly.
! It also sets up the x and y arrays if running in box mode.
! 
! The range mode is typically used for linear simulations, while the box
! mode is typically ised for nonlinear simulations.
! 
!                                 RANGE
! 
! Range mode is used when one wants to specify a range of kx and ky values
! and the number of modes in each direction.
! 
! In this mode, ky is assumed to be positive and the user specifies
! naky, aky_min and aky_max. The ky values are then spaced linearly or exponentially
! between aky_min and aky_max depending on the value of kyspacing_option_switch.
! 
! The user also specifies nakx and either akx_min and akx_max or
! theta0_min and theta0_max. If the latter are specified, then
! akx_min and akx_max are determined using the relation
! theta0 = kx/(ky*shat) (or theta0 = kx/ky if q_as_x = .true.).
! 
!                                  BOX
! 
! In Box mode, the code sets up a two-dimensional grid in real space (x, y) and
! the corresponding Fourier space (kx, ky). This mode is typically used for
! nonlinear simulations. The user specifies the number of grid points (nx, ny),
! the box sizes (x0, y0), and other parameters such as jtwist and phase_shift_angle.
! The code then calculates the kx and ky arrays, as well as the real-space x and y
! arrays, ensuring consistency with the chosen boundary conditions (e.g., twist-and-shift,
! periodic, or linked). Additional arrays for radial coordinates (rho, rho_d) and
! their clamped versions are also set up if radial variation is enabled.
! 
!----------------------------------- Issues ------------------------------------
! 
! TODO - For the "range" mode we assume only a single ky mode has been specified
! if a range of kx modes are specified, because the ballooning angle is defined
! as theta0 = kx / (hat{s} * ky_min) therefore, it is badly defined for ky values
! that are not ky = ky_min. A range of ky-modes only works for kx = 0.
! 
!###############################################################################
module grids_kxky
   
   ! Read the parameters for <grid_option_switch> from namelist_kxky_grid.f90
   use namelist_kxky_grid, only: grid_option_range
   use namelist_kxky_grid, only: grid_option_box
   
   ! Read the parameters for <kyspacing_option_switch> from namelist_kxky_grid.f90
   use namelist_kxky_grid, only: kyspacing_linear
   use namelist_kxky_grid, only: kyspacing_exponential

   implicit none
   
   ! Although the parameters are available through namelist_kxky_grid
   ! make them available through grids_kxky as well
   public :: grid_option_switch
   public :: grid_option_range
   public :: grid_option_box
   public :: kyspacing_option_switch
   public :: kyspacing_linear
   public :: kyspacing_exponential

   ! Make the routines available to other modules
   public :: init_grids_kxky
   public :: finish_grids_kxky
   public :: read_parameters_kxky_grids
   
   ! Make the parameters available to other modules
   public :: aky, akx
   public :: aky_all, aky_all_ordered
   public :: theta0, zed0
   public :: zonal_mode
   public :: x, x_d, y
   public :: dy, dx
   public :: rho, rho_d, rho_clamped, rho_d_clamped
   public :: g0x
   public :: box

   ! For Box/Range
   public :: naky, nakx
   public :: nx, ny
   public :: nalpha, naky_all, ikx_max
   public :: reality
   public :: phase_shift_angle
   public :: jtwist, jtwistfac, ikx_twist_shift
   public :: centered_in_rho
   public :: periodic_variation, randomize_phase_shift
   public :: aky_min, aky_max
   public :: akx_min, akx_max
   public :: theta0_min, theta0_max
   public :: x0, y0

   ! Check initialisation
   public :: initialised_grids_kxky
   
   private 

   ! For the (kx,ky) grids
   real, dimension(:), allocatable :: aky, akx
   real, dimension(:), allocatable :: aky_all, aky_all_ordered
   real, dimension(:, :), allocatable :: theta0, zed0
   real, dimension(:), allocatable :: x, x_d, y
   logical, dimension(:), allocatable :: zonal_mode

   ! For radial variation
   real, dimension(:), allocatable :: rho, rho_d, rho_clamped, rho_d_clamped
   complex, dimension(:, :), allocatable :: g0x

   ! Required for flux calculations
   real :: dx, dy

   ! Internal calculations
   real :: dkx, dky, dx_d
   logical :: box
   
   ! Parameters
   integer :: grid_option_switch
   integer :: naky, nakx, nx, ny
   integer :: nalpha, naky_all, ikx_max
   logical :: reality = .false.
   real :: phase_shift_angle
   integer :: jtwist
   real :: jtwistfac
   real :: ikx_twist_shift
   logical :: centered_in_rho, periodic_variation, randomize_phase_shift

   ! For Range
   real :: aky_min, aky_max
   real :: akx_min, akx_max
   real :: theta0_min, theta0_max
   integer :: kyspacing_option_switch

   ! For Box
   real :: x0, y0
   
   ! Only initialise once
   logical :: initialised_grids_kxky
   logical :: initialised_read_parameters_kxky_grids

contains

!###############################################################################
!################################ READ NAMELIST ################################
!###############################################################################

   subroutine read_parameters_kxky_grids

      use mp, only: proc0, mp_abort
      use namelist_kxky_grid, only: read_namelist_kxky_grid_option
      use namelist_kxky_grid, only: read_namelist_kxky_grid_box
      use namelist_kxky_grid, only:read_namelist_kxky_grid_range
      
      implicit none
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_read_parameters_kxky_grids) return
      initialised_read_parameters_kxky_grids = .true.

      ! Read the "kxky_grid_option" namelist in the input file
      call read_namelist_kxky_grid_option (grid_option_switch)
      
      ! Read the "kxky_grid_range" or "kxky_grid_box" namelist in the input file
      if (proc0) then
         select case (grid_option_switch)
         case (grid_option_range)
            call read_namelist_kxky_grid_range (nalpha, naky, nakx, aky_min, aky_max, & 
               akx_min, akx_max, theta0_min, theta0_max, &
               kyspacing_option_switch, phase_shift_angle, ikx_max, naky_all)
         case (grid_option_box)
            call read_namelist_kxky_grid_box (nx, ny, ikx_max, naky_all, naky, nakx, nalpha, &
               x0, y0, jtwist, jtwistfac, phase_shift_angle, &
               centered_in_rho, randomize_phase_shift, periodic_variation, reality)
         end select
      end if
      
      ! Broadcast the parameters to all processors
      call broadcast_parameters

   contains
    
      ! Broadcast parameters to all processes
      subroutine broadcast_parameters
         
         use mp, only: broadcast
         
         implicit none

         call broadcast(grid_option_switch)
         call broadcast(naky)
         call broadcast(nakx)
         call broadcast(ny)
         call broadcast(nx)
         call broadcast(nalpha)
         call broadcast(naky_all)
         call broadcast(ikx_max)
         call broadcast(reality)
         call broadcast(phase_shift_angle)
         call broadcast(jtwist)
         call broadcast(jtwistfac)
         call broadcast(ikx_twist_shift)
         call broadcast(centered_in_rho)
         call broadcast(periodic_variation)
         call broadcast(randomize_phase_shift)
         call broadcast(aky_min)
         call broadcast(aky_max)
         call broadcast(akx_min)
         call broadcast(akx_max)
         call broadcast(theta0_min)
         call broadcast(theta0_max)
         call broadcast(kyspacing_option_switch)
         call broadcast(x0)
         call broadcast(y0)

      end subroutine broadcast_parameters
      
   end subroutine read_parameters_kxky_grids


!###############################################################################
!########################### INITIALISE (KX,KY) GRIDS ##########################
!###############################################################################

   !****************************************************************************
   !                        Initialise the (kx,ky) grids                       !
   !****************************************************************************
   subroutine init_grids_kxky

      use mp, only: mp_abort
      use common_types, only: flux_surface_type
      use grids_z, only: read_parameters_z_grid
      
      implicit none
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_grids_kxky) return
      initialised_grids_kxky = .true.

      ! Make sure the parameters in the z-grid namelists have been read
      call read_parameters_z_grid
      
      ! Construct the (kx,ky) grid using the "range" or "box" mode
      select case (grid_option_switch)
      case (grid_option_range)
         call init_grids_kxky_range
      case (grid_option_box)
         call init_grids_kxky_box
      end select

      ! Determine if <iky> corresponds to a zonal mode, i.e., whether ky=0
      if (.not. allocated(zonal_mode)) allocate (zonal_mode(naky))
      zonal_mode = .false.
      if (abs(aky(1)) < epsilon(0.)) zonal_mode(1) = .true.

   contains

      !*************************************************************************
      !                Initialise the (kx,ky) grids: Range mode                !
      !*************************************************************************
      subroutine init_grids_kxky_range

         use common_types, only: flux_surface_type
         use geometry, only: geo_surf, q_as_x
         use grids_z, only: shat_zero
           
         implicit none

         ! Local variables
         integer :: i, j
         real :: tfac
         real :: zero
         real :: dkx = 0.0
         real :: dky = 0.0
         real :: dtheta0 = 0.0
      
         !----------------------------------------------------------------------

         ! Recall that we are in "range" mode and not in "box" mode
         box = .false.
         
         ! Allocate the arrays (akx, aky, aky_all, aky_all_ordered, theta0, zed0)
         call allocate_arrays
         
         ! Create a range of ky-modes
         if (naky > 1) then
            select case (kyspacing_option_switch)
            case (kyspacing_linear)
               dky = (aky_max - aky_min) / real(naky - 1)
               aky = (/(aky_min + dky * real(i), i=0, naky - 1)/)
            case (kyspacing_exponential)
               dky = (log(aky_max) - log(aky_min)) / real(naky - 1)
               aky = (/(exp(log(aky_min) + dky * real(i)), i=0, naky - 1)/)
            end select
            
         ! Simulate a single ky-mode
         else
            aky = (/(aky_min, i=0, naky - 1)/)
         end if

         ! Initialise the values for akx and theta0
         akx = 0.0
         theta0 = 0.0

         ! Define the ballooning angle as theta0 = kx / (ky*shat)
         if (.not. q_as_x) tfac = geo_surf%shat
         
         ! If x = q we define the ballooning angle as theta0 = kx / ky
         if (q_as_x) tfac = 1.0
   
         ! Define a numerical zero
         zero = 100.*epsilon(0.)

         ! If <theta0_min> and <theta0_max> have been specified,
         ! use them to determine <akx_min> and <akx_max>
         if (theta0_max > theta0_min - zero) then
         
            ! Positive global magnetic shear
            if (geo_surf%shat > epsilon(0.)) then
               akx_min = theta0_min * tfac * aky(1)
               akx_max = theta0_max * tfac * aky(1)
               
            ! Negative global magnetic shear (link theta0_max to akx_min)
            else
               akx_min = theta0_max * tfac * aky(1)
               akx_max = theta0_min * tfac * aky(1)
            end if
            
         end if
   
         ! If <shat> is bigger than <shat_zero>, we assume linked boundary conditions are used
         if (abs(geo_surf%shat) > shat_zero) then
         
            ! If <theta0_min> and <theta0_max> had been specified in the input file
            ! then we have already calculated the corresponding <akx_min> and <akx_max>.
            ! If instead, <akx_min> and <akx_max> are specified in the input file instead of
            ! <theta0_min> and <theta0_max>, use them to calculate <theta0_min> and <theta0_max>.
            ! Next, use these values to construct the <akx> and <theta0> arrays.
            if (theta0_min > theta0_max + zero .and. abs(aky(1)) > zero) then
            
               ! Ballooning angle is theta0 = kx / (ky*shat); or theta0 = kx / ky if x = q
               ! The minimum and maximum balooning angle occur for the smallest value of ky
               theta0_min = akx_min / (tfac * aky(1))
               theta0_max = akx_max / (tfac * aky(1))
               
               ! Determine the step size in theta0; if nakx==1 then dtheta0 = 0
               if (nakx > 1) dtheta0 = (theta0_max - theta0_min) / real(nakx - 1)

               ! We assume that only a single ky value has been specified
               ! For this ky-value we construct a single chain of (kx,ky) modes
               ! TODO - Then why do we sum over ky?
               do j = 1, naky
                  theta0(j, :) = (/(theta0_min + dtheta0 * real(i), i=0, nakx - 1)/)
               end do
               
               ! Ballooning angle is theta0 = kx / (ky*shat); or theta0 = kx / ky if x = q
               ! Therefore, kx = ky * shat * theta0; or kx = theta0 * ky if x = q
               akx = theta0(1, :) * tfac * aky(1)
               
            ! If <theta0_min> and <theta0_max> had been specified in the input file
            ! then we have already calculated the corresponding <akx_min> and <akx_max>.
            ! Use these values to construct the <akx> and <theta0> arrays.
            else if (akx_max > akx_min - zero .or. nakx == 1) then

               ! Get the step size in kx from <akx_min> and <akx_max>; if nakx==1 then dkx = 0
               if (nakx > 1) dkx = (akx_max - akx_min) / real(nakx - 1)
               
               ! Construct the range of kx-modes
               akx = (/(akx_min + dkx * real(i), i=0, nakx - 1)/)
               
               ! Get the step size in theta0 from <theta0_max> and <theta0_min>; if nakx==1 then dtheta0 = 0
               if (nakx > 1) dtheta0 = (theta0_max - theta0_min) / real(nakx - 1)
               
               ! Construct the theta0 array for positive global magnetic shear
               if (geo_surf%shat > epsilon(0.)) then
                  do j = 1, naky
                     theta0(j, :) = (/(theta0_min + dtheta0 * real(i), i=0, nakx - 1)/)
                  end do
                  
               ! Construct the theta0 array for negative global magnetic shear (link theta0_max to akx_min)
               else
                  do j = 1, naky
                     theta0(j, :) = (/(theta0_min + dtheta0 * real(i), i=nakx - 1, 0, -1)/)
                  end do
               end if
               
            else
               call mp_abort('The choice ky=0 is inconsistent with akx_min different from akx_max. Aborting.')
            end if

         ! If <shat> is smaller than <shat_zero>, periodic boundary conditions are enforced
         ! These boundary conditions are used for periodic finite kx ballooning space runs with shat=0
         ! In this case the theta0 array contains only zeros since theta0 = kx / (hat{s} * ky) = 0 if shat=0
         else
            if (nakx > 1) dkx = (akx_max - akx_min) / real(nakx - 1)
            akx = (/(akx_min + dkx * real(i), i=0, nakx - 1)/)
         end if

         ! When we plot the modes on the extended z-grid from the ballooning transformation,
         ! we need to move the z-domains by a factor zed0, i.e., we plot |phi|^2 versus (zed(iz) - zed0(iky, ikx))
         ! with zed0(ky,kx) = theta0 * <zed0_fac> and <zed0_fac> = max(zed) / max(theta)
         ! For Miller geometries, the parallel coordinate is z = theta, thus we simply have <zed0_fac> = 1.
         ! For VMEC geometries, z = zeta, thus <zed0_fac> = max(zed) / max(theta) = max(zed) / max(zeta) * q
         ! where we used the safety factor q = dzeta/dtheta to convert theta to zeta.
         zed0 = theta0 * geo_surf%zed0_fac

         ! The following variables are more important for "box" simulations
         ikx_max = nakx
         naky_all = naky
         
      end subroutine init_grids_kxky_range

    
      !*************************************************************************
      !                 Initialise the (kx,ky) grids: Box mode                 !
      !*************************************************************************
      subroutine init_grids_kxky_box
           
         ! Parallelisation
         use mp, only: mp_abort
         use mp, only: proc0
         use mp, only: broadcast
         
         ! Calculations
         use constants, only: pi, zi
         use ran, only: ranf
         
         ! Flags
         use parameters_physics, only: full_flux_surface
         
         ! Parallel boundary condition
         use grids_z, only: nperiod
         use grids_z, only: boundary_option_switch
         use grids_z, only: boundary_option_linked
         use grids_z, only: boundary_option_linked_stellarator
         
         ! Geometry
         use parameters_physics, only: rhostar
         use common_types, only: flux_surface_type
         use geometry, only: twist_and_shift_geo_fac
         use geometry, only: geo_surf
         use geometry, only: dydalpha
         use geometry, only: q_as_x
         use geometry, only: geo_option_switch
         use geometry, only: geo_option_vmec
         
         implicit none
   
         ! Local variables
         integer :: ikx, iky
         integer :: ikyneg
         real :: norm
      
         !----------------------------------------------------------------------
   
         ! Recall that we are in "box" mode and not in "range" mode
         box = .true.
         
         ! Allocate the arrays (akx, aky, aky_all, aky_all_ordered, theta0, zed0)
         ! as well as (x_d, rho, rho_d, rho_clamped, rho_d_clamped, x, y)
         call allocate_arrays
         
         ! Set jtwist if it has not been specified in the input file
         ! The variable <jtwist> determines the number of eigenmodes at ky_min
         if (jtwist < 1) then
            jtwist = max(1, int(abs(twist_and_shift_geo_fac) + 0.5))
            jtwist = max(1, int(jtwistfac * jtwist + 0.5))
         end if
         
         ! Define a signed version of jtwist, with sign determined by, e.g., magnetic shear
         ikx_twist_shift = -jtwist * int(sign(1.0, twist_and_shift_geo_fac))
   
         ! Set y0 for a full-flux-surface simulation based on <rhostar>
         if (y0 < 0.) then
         
            ! When simulating a flux annulus, y0 is determined by the physical 
            ! extent of the device, i.e., y0 = geo_surf%rhotor/rhostar
            if (full_flux_surface) then
               if (rhostar > 0.) then
                  y0 = geo_surf%rhotor / rhostar
               else
                  call mp_abort('Must set rhostar if simulating a full flux surface. Aborting.')
               end if
               
            ! When simulating a flux tube, it makes no sense to have y0 < 0.0, abort.
            else
               call mp_abort('For a flux tube simulation, it is mandotary to set y0. Aborting.')
            end if
            
         end if

         ! The grid spacing in ky is determined by the width of the flux tube in real space
         ! with y0 = Ly / (2*pi*rho_ref) = 1 / (ky_SI_min * rho_ref) = 1 / <dky>
         dky = 1./y0
         
         ! Use the twist-and-shift parallel boundary condition to calculate the grid 
         ! spacing in kx, i.e., kx = ky * twist_shift_geo_fac / jtwist.
         ! For periodic boundary conditions (default case) we simply set dkx = dky if x0 is not specified.
         select case (boundary_option_switch)
         case (boundary_option_linked)
            dkx = (2 * nperiod - 1) * dky * abs(twist_and_shift_geo_fac) / real(jtwist)
         case (boundary_option_linked_stellarator)
            dkx = dky * abs(twist_and_shift_geo_fac) / real(jtwist)
         case default
            if (x0 < epsilon(0.0)) then
               dkx = dky
            else
               dkx = 1./x0
            end if
         end select
         
         ! Now that we have calculated <dkx> we can set the width of the flux tube box in real space
         ! with x0 = Lx / (2*pi*rho_ref) = 1 / (kx_SI_min * rho_ref) = 1 / <dkx>
         x0 = 1./dkx

         ! Contruct the ky-array, which goes from zero to aky_max. Negative ky-values are
         ! omitted thanks to the reality condition: phi(kx,ky,x) = conj(phi(-kx,-ky,z)
         do iky = 1, naky
            aky(iky) = real(iky - 1) * dky
         end do
         
         ! Construct the full ky-array, running from -aky_max to aky_max as <aky_all_ordered>
         ! Note that <naky_all> = 2 * naky - 1 has been calculated inside namelist_kxky_grid.f90
         ! Also contruct <aky_all>, stored in the same order as akx (0 -> aky_max, -aky_max -> -dky)
         
         ! First set the arrays equal to the aky array for ky >= 0
         aky_all_ordered(naky:naky_all) = aky
         aky_all(:naky) = aky
         
         ! Next fill in ky < 0, with <ikyneg> the ky index corresponding to +ky in original array
         do iky = naky + 1, naky_all 
            ikyneg = naky_all - iky + 2
            aky_all(iky) = -aky(ikyneg)
         end do
         aky_all_ordered(:naky - 1) = aky_all(naky + 1:)
   
         ! Construct the kx array, stored as (0 -> akx_max, -akx_max -> -dkx)
         do ikx = 1, ikx_max
            akx(ikx) = real(ikx - 1) * dkx
         end do
         do ikx = ikx_max + 1, nakx
            akx(ikx) = real(ikx - nakx - 1) * dkx
         end do
   
         ! Define the ballooning angle as theta0 = kx / (ky*shat)
         ! If x = q we define the ballooning angle as theta0 = kx / ky
         if (.not. q_as_x) then
            do ikx = 1, nakx
               theta0(2:, ikx) = akx(ikx) / (aky(2:) * geo_surf%shat)
            end do
         else if (q_as_x) then
            do ikx = 1, nakx
               theta0(2:, ikx) = akx(ikx) / aky(2:)
            end do
         end if
         
         ! For zonal modes, i.e., for ky = 0, set theta0 = 0.
         theta0(1, :) = 0.0

         ! We omitted negative ky-values using the reality condition: phi(kx,ky,x) = conj(phi(-kx,-ky,z)
         ! When summing over all modes, we sum ky=0 once and ky>0 twice to take this into account
         norm = 1.
         if (naky > 1) norm = aky(2)
         
         ! TODO - write comments on phase_shift_angle
         if (rhostar > 0.) then
            if (geo_option_switch == geo_option_vmec) then
               phase_shift_angle = -2.*pi * (2 * nperiod - 1) * dydalpha / (rhostar * geo_surf%qinp_psi0)
            else
               phase_shift_angle = -2.*pi * (2 * nperiod - 1) * geo_surf%qinp_psi0 * dydalpha / rhostar
            end if
         else if (randomize_phase_shift) then
            if (proc0) phase_shift_angle = 2.*pi * ranf() / norm
            call broadcast(phase_shift_angle)
         else
            phase_shift_angle = 2.*pi * phase_shift_angle / norm
         end if
   
         ! Calculate extra (kx,ky) related grids for radial variation runs
         call init_grids_kxky_box_radial_variation()
         
      end subroutine init_grids_kxky_box
      
      !*************************************************************************
      !      Initialise the (kx,ky) grids: Box mode for radial variation       !
      !*************************************************************************
      subroutine init_grids_kxky_box_radial_variation
      
         ! Parallelisation
         use mp, only: broadcast
         
         ! Calculations
         use constants, only: pi
         
         ! Geometry
         use parameters_physics, only: rhostar
         use geometry, only: geo_surf
         use geometry, only: q_as_x
         use geometry, only: drhodpsi
         use geometry, only: dxdpsi
         use geometry, only: get_x_to_rho

         ! Radial variation
         use file_utils, only: runtype_option_switch
         use file_utils, only: runtype_multibox
         use write_radial_grid, only: dump_radial_grid
         use parameters_physics, only: radial_variation
         
         implicit none
         
         ! Local variables
         integer :: ikx, iky
         real :: x_shift
         real :: dqdrho
         real :: pfac
      
         !----------------------------------------------------------------------
      
         ! Calculate the step size in real space
         dx = (2 * pi * x0) / nx
         dy = (2 * pi * y0) / ny
   
         x_shift = pi * x0
         pfac = 1.0
         if (periodic_variation) pfac = 0.5
         if (centered_in_rho) then
            if (q_as_x) then
               dqdrho = geo_surf%shat * geo_surf%qinp / geo_surf%rhoc
               x_shift = pi * x0 * (1.0 - 0.5 * pfac * rhostar * pi * x0 * geo_surf%d2qdr2 / (dqdrho**2 * dxdpsi))
            else
               x_shift = pi * x0 * (1.0 - 0.5 * pfac * rhostar * pi * x0 * geo_surf%d2psidr2 * drhodpsi**2 / dxdpsi)
            end if
         end if
   
         do ikx = 1, nx
            if (radial_variation .or. runtype_option_switch == runtype_multibox) then
               if (periodic_variation) then
                  if (ikx <= nx / 2) then
                     x(ikx) = (ikx - 1) * dx - 0.5 * x_shift
                  else
                     x(ikx) = x(nx - ikx + 1)
                  end if
               else
                  x(ikx) = (ikx - 0.5) * dx - x_shift
               end if
            else
               x(ikx) = (ikx - 1) * dx
            end if
         end do
   
         dx_d = (2 * pi * x0) / nakx
         do ikx = 1, nakx
            if (radial_variation .or. runtype_option_switch == runtype_multibox) then
               if (periodic_variation) then
                  if (ikx <= (nakx + 1) / 2) then
                     x_d(ikx) = (ikx - 1) * dx_d - 0.5 * x_shift
                  else
                     x_d(ikx) = x_d(nakx - ikx + 1)
                  end if
               else
                  x_d(ikx) = (ikx - 0.5) * dx_d - x_shift
               end if
            else
               x_d(ikx) = (ikx - 1) * dx_d
            end if
         end do
   
         call get_x_to_rho(1, x, rho)
         call get_x_to_rho(1, x_d, rho_d)

         rho_clamped = rho
         rho_d_clamped = rho_d
         
         zed0 = theta0 * geo_surf%zed0_fac
   
         if (radial_variation) call dump_radial_grid (x, rho, nx)
   
         if (radial_variation .and. (any((rho + geo_surf%rhoc) < 0.0) .or. any((rho + geo_surf%rhoc) > 1.0))) then
            call mp_abort('rho(x) is beyond range [0,1]. Try changing rhostar or q/psi profiles')
         end if
   
         do iky = 1, ny
            y(iky) = (iky - 1) * dy
         end do
         
         ! Broadcast the radial variation parameters to all processors
         call broadcast(x_d)
         call broadcast(rho)
         call broadcast(rho_d)
         call broadcast(rho_clamped)
         call broadcast(rho_d_clamped)
         
      end subroutine init_grids_kxky_box_radial_variation
       
      !*************************************************************************
      !              Allocate arrays needed for the (kx,ky) grids              !
      !*************************************************************************
      subroutine allocate_arrays
         
         implicit none
      
         !----------------------------------------------------------------------
         
         ! Arrays needed for both "range" and "box" mode
         if(.not. allocated(akx)) allocate (akx(nakx))
         if(.not. allocated(aky)) allocate (aky(naky))
         if(.not. allocated(aky_all)) allocate(aky_all(naky_all))
         if(.not. allocated(aky_all_ordered)) allocate (aky_all_ordered(naky_all))
         if(.not. allocated(theta0)) allocate (theta0(naky, nakx))
         if(.not. allocated(zed0)) allocate (zed0(naky, nakx))
         
         ! Arrays needed only for "box" mode
         if (box) then
            if (.not. allocated(x_d)) allocate (x_d(nakx))
            if (.not. allocated(rho)) allocate (rho(nx))
            if (.not. allocated(rho_d)) allocate (rho_d(nakx))
            if (.not. allocated(rho_clamped)) allocate (rho_clamped(nx))
            if (.not. allocated(rho_d_clamped)) allocate (rho_d_clamped(nakx))
            if (.not. allocated(x)) allocate (x(nx))
            if (.not. allocated(y)) allocate (y(ny))
         end if
      
      end subroutine allocate_arrays
      
   end subroutine init_grids_kxky
   
!###############################################################################
!############################ FINISH (KX,KY) GRIDS #############################
!###############################################################################

   subroutine finish_grids_kxky

      implicit none

      if (allocated(aky)) deallocate (aky)
      if (allocated(aky_all)) deallocate (aky_all)
      if (allocated(aky_all_ordered)) deallocate (aky_all_ordered)
      if (allocated(akx)) deallocate (akx)
      if (allocated(theta0)) deallocate (theta0)
      if (allocated(zed0)) deallocate (zed0)
      if (allocated(x)) deallocate (x)
      if (allocated(y)) deallocate (y)
      if (allocated(x_d)) deallocate (x_d)
      if (allocated(rho)) deallocate (rho)
      if (allocated(rho_d)) deallocate (rho_d)
      if (allocated(rho_clamped)) deallocate (rho_clamped)
      if (allocated(rho_d_clamped)) deallocate (rho_d_clamped)
      if (allocated(g0x)) deallocate (g0x)

      reality = .false.
      initialised_grids_kxky = .false.
      
   end subroutine finish_grids_kxky

end module grids_kxky
