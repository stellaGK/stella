!###############################################################################
!                            VELOCITY (VPA, MU) GRIDS                           
!###############################################################################
! This module initiates the grids for the parallel and perpendicular velocity
! and the associated integration weights to perform integrations in v-space. The
! integrations in v-space (along mu and vpa) and over the species are performed
! as weighted sums with the weights <wgts_mu>; <wgts_vpa> and <wgts_s>.
! 
!                               VPA AND MU GRIDS
! 
! Read the parameters in the namelist <velocity_grids> and initialize the
! vpa and mu grids and the integration weights. Calculate the Maxwell factors:
!       <maxwell_fac> = (n/n_psi0)*(T_psi0/T)**3/2
!       <maxwell_vpa> = exp(-vpa**2*(T_psi0/T))
!       <maxwell_mu> = exp(-2*mu*B*(T_psi0/T))
!
! The quantity vpa is the parallel velocity at grid points and wgts_vpa are the
! integration weights assigned to the parallel velocity grid points.
! The velocity grid goes from -vpa_max to vpa_max without a point at vpa = 0.
! The integration weights corresponding to the vpa grid points are obtained through
! the Simpson's 3/8 rule at lower/upper boundary and composite Simpson elsewhere.
! 
! The quantity mu is related to the perpendicular velocity at grid points and
! wgts_mu are the weights assigned to the perpendicular velocity grid points.
! Either we create an equally spaced grid with mu_max = vperp_max**2/(2*max(B)),
! or we create a dynamically spaced grid based on the Gauss-Laguerre quadrature.
! 
!                                INTEGRATIONS
! 
! Perform integrations over v-space. Integrating over the perpendicular velocity
! is equivalent to taking weighted sums over mu = [1 --> nmu] of the quantity <g>
! with weights wgts_mu(ia,iz,imu). Integrating over the parallel velocity is
! equivalent to taking weighted sums over v = [-nvgrid --> nvgrid] of the quantity
! <g> with weights wgts_vpa(iv). Integrating over the species is equivalent to
! summing the contributions of each species with the weights ws(is).
! 
! The interface "integrate_mu" integrates over the perpendicular velocity.
! The interface "integrate_vmu" integrates over both velocities.
! The interface "integrate_species" integrates over both velocities and the species.
! 
! Each interface has a subroutine to deal with complex and real distribution
! functions <g>. Moreover the integration can be performed locally by summing
! over all (imu,ivmu,is) or non-locally by summing over the ivmu points and then
! summing the contributions of all processors. Finally there are "single" and
! "block" routines where <g> in the block routines has indices (kx,ky) as well.
! 
!                             INTEGRATION WEIGHTS
! 
! Integrations over velocity space always have the following form:
!     int dmu int dvpa 2 bmag(z) / sqrt(pi)
! 
! We will absorb the factor [ 2 bmag(z) / sqrt(pi) ] into the integration weights:
!     <wgts_vpa>[ivpa] = dvpa[ivpa] / sqrt(pi)
!     <wgts_mu>[ialpha, iz, imu] = dmu[imu] * 2 * bmag[iz]
! 
!###############################################################################
module grids_velocity

   implicit none

   ! Make routines accesible to other modules
   public :: init_velocity_grids
   public :: finish_velocity_grids
   public :: read_parameters_velocity_grids
   
   ! Grid points
   public :: vpa, nvgrid, nvpa
   public :: mu, nmu
   public :: vperp2
   
   ! Integration weights
   public :: wgts_vpa, dvpa
   public :: wgts_mu, wgts_mu_bare, dmu
   
   ! The factor exp(v²) = exp(v²_parallel) * exp(v²_perp)
   public :: maxwell_vpa, maxwell_mu
   public :: maxwell_fac, ztmax
   public :: maxwell_mu_avg
   
   ! The following factors are used in the collision operators
   public :: int_unit, int_vpa2, int_vperp2, int_vfrth
   public :: dmu_cell, mu_cell
   public :: set_vpa_weights

   private

   ! Velocity grids
   integer :: nvgrid, nvpa, nmu
   real :: vpa_max, vperp_max, dvpa
   real, dimension(:), allocatable :: dmu
   real, dimension(:), allocatable :: vpa
   real, dimension(:), allocatable :: mu
   
   ! Integration weights
   real, dimension(:), allocatable :: wgts_vpa, wgts_vpa_default
   real, dimension(:), allocatable :: wgts_mu_bare
   real, dimension(:, :, :), allocatable :: wgts_mu
   
   ! The factor exp(v²) = exp(v²_parallel) * exp(v²_perp)
   real, dimension(:), allocatable :: maxwell_fac
   real, dimension(:, :), allocatable :: maxwell_vpa
   real, dimension(:, :, :, :), allocatable :: maxwell_mu
   real, dimension(:, :, :, :), allocatable :: maxwell_mu_avg
   
   ! Other arrays used in the collision operators
   real, dimension(:, :, :), allocatable :: int_unit, int_vpa2, int_vperp2, int_vfrth
   real, dimension(:, :), allocatable :: ztmax
   real, dimension(:), allocatable :: dmu_cell, mu_cell
   
   ! Flags
   logical :: equally_spaced_mu_grid
   logical :: conservative_wgts_vpa

   ! Arrays related to the (vpa,mu) grid that are declared here
   ! but allocated and filled elsewhere because they depend on z, etc.
   real, dimension(:, :, :), allocatable :: vperp2

   ! Only initialise once
   logical :: initialised_velocity_grids = .false.
   logical :: initialised_read_velocity_grids = .false.

contains

!###############################################################################
!################################ READ NAMELIST ################################
!###############################################################################

   subroutine read_parameters_velocity_grids

      use namelist_velocity_grids, only: read_namelist_velocity_grids

      implicit none
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_read_velocity_grids) return
      initialised_read_velocity_grids = .true.
      
      ! Read the "velocity_grids" namelist in the input file
      call read_namelist_velocity_grids(nvgrid, nmu, vpa_max, vperp_max, &
         equally_spaced_mu_grid, conservative_wgts_vpa)
         
      ! Broadcast the parameters to all processors
      call broadcast_velocity_grids

      ! Calculate <nvpa> based on <nvgrid>
      nvpa = 2 * nvgrid

   contains
   
      !-------------------------- Broadcast parameters--------------------------
      subroutine broadcast_velocity_grids
      
         use mp, only: broadcast

         implicit none

         call broadcast(nvgrid)
         call broadcast(nmu)
         call broadcast(vpa_max)
         call broadcast(vperp_max)
         call broadcast(equally_spaced_mu_grid)
         call broadcast(conservative_wgts_vpa)

      end subroutine broadcast_velocity_grids

   end subroutine read_parameters_velocity_grids

!###############################################################################
!########################### INITIALISE VELOCITY GRIDS #########################
!###############################################################################

   subroutine init_velocity_grids

      use grids_species, only: read_parameters_species
      use parameters_numerical, only: read_parameters_numerical

      implicit none
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_velocity_grids) return
      initialised_velocity_grids = .true.
      
      ! Make sure the dependencies of the velocity grids are initialised
      call read_parameters_numerical
      call read_parameters_species

      ! Set up the vpa and mu grid points and integration weights
      call init_vpa_grid
      call init_mu_grid
      call init_maxwellian_factors
      call init_velocity_parameters_for_collision_operators
      
   end subroutine init_velocity_grids


   !****************************************************************************
   !                     Initialise parallel velocity grid                     !
   !****************************************************************************
   ! The parallel velocity grid goes from -<vpa_max> to <vpa_max>, with no 
   ! point at <vpa=0>. The lack of a point at <vpa=0> avoids treating the
   ! <vpa=z=0> phase space location, which is isolated from all other phase 
   ! space points in the absence of collisions.
   !****************************************************************************
   subroutine init_vpa_grid

      implicit none
      
      call allocate_arrays
      call initialise_vpa_grid
      call initialise_vpa_integration_weights
      
   contains
   
      !---------------------------- Allocate arrays ----------------------------
      subroutine allocate_arrays
      
         use grids_species, only: nspec

         implicit none

         if (.not. allocated(vpa)) then
            allocate (vpa(nvpa)); vpa = 0.0
            allocate (wgts_vpa(nvpa)); wgts_vpa = 0.0
            allocate (wgts_vpa_default(nvpa)); wgts_vpa_default = 0.0
            allocate (maxwell_vpa(nvpa, nspec)); maxwell_vpa = 0.0
            allocate (ztmax(nvpa, nspec)); ztmax = 0.0
         end if
         
      end subroutine allocate_arrays
      
      !-------------------------- Initialise vpa grid --------------------------
      subroutine initialise_vpa_grid

         implicit none
         
         integer :: iv
      
         !----------------------------------------------------------------------

         ! Equal grid spacing in parallel velocity
         dvpa = 2.*vpa_max / (nvpa - 1)

         ! Obtain vpa grid for vpa > 0
         do iv = nvgrid + 1, nvpa
            vpa(iv) = real(iv - nvgrid - 0.5) * dvpa
         end do
         
         ! Fill in vpa grid for vpa < 0
         vpa(:nvgrid) = -vpa(nvpa:nvgrid + 1:-1)

      end subroutine initialise_vpa_grid
      
      !-------------------------------------------------------------------------
      !------------------- Initialise vpa integration weights ------------------
      !-------------------------------------------------------------------------
      ! Integrations over velocity space always have the following form:
      !     int dmu int dvpa 2 bmag(z) / sqrt(pi)
      ! 
      ! We will absorb the factor [ 1 / sqrt(pi) ] into the integration weight:
      !     <wgts_vpa>[ivpa] = <Simpsons_weights>[ivpa] / sqrt(pi)
      !-------------------------------------------------------------------------
      subroutine initialise_vpa_integration_weights
      
         use constants, only: pi

         implicit none
      
         !---------------------------------------------------------------------- 
         
         ! Conservative weights conserve the density in the collision operator
         if (conservative_wgts_vpa) then
            wgts_vpa = dvpa / sqrt(pi)
            
         ! By default, we use the Simpson's rule for the integration weights
         else 
            call initialise_vpa_integration_weights_Simpsons
         end if
   
         ! The collision operator will switch between the chosen integration 
         ! weights, which can be conservative or Simpson's weights, and 
         ! conservative weights, so we need to store the original weights
         wgts_vpa_default = wgts_vpa
         
      end subroutine initialise_vpa_integration_weights
      
      !----------------------------------------------------------------------
      !----------- Calculate parallel velocity integration weights ----------
      !----------------------------------------------------------------------
      ! Get integration weights corresponding to vpa grid points. For now use 
      ! Use Simpson's rule: i.e. subdivide grid into 3-point segments, with each 
      ! segment spanning vpa_low to vpa_up, then the contribution of each 
      ! segment to the integral is (vpa_up - vpa_low) * (f1 + 4*f2 + f3) / 6
      ! Inner boundary points are used in two segments, so they get double the weight.
      !----------------------------------------------------------------------
      subroutine initialise_vpa_integration_weights_Simpsons
      
         use mp, only: mp_abort
         use constants, only: pi
         use parameters_numerical, only: maxwellian_normalization

         implicit none
          
         integer :: idx, iseg, nvpa_seg
         real :: del
      
         !----------------------------------------------------------------------

         ! The Simpson's rule only works for a minimum of 6 grid points.
         if (nvpa < 6) call mp_abort('stella does not currently support nvgrid < 3. Aborting.')

         ! Use simpson 3/8 rule at lower boundary and composite Simpson elsewhere
         del = 0.375 * dvpa
         wgts_vpa(1) = del
         wgts_vpa(2:3) = 3.*del
         wgts_vpa(4) = del
         
         ! Composite simpson
         nvpa_seg = (nvpa - 4) / 2
         del = dvpa / 3.
         do iseg = 1, nvpa_seg
            idx = 2 * (iseg - 1) + 4
            wgts_vpa(idx) = wgts_vpa(idx) + del
            wgts_vpa(idx + 1) = wgts_vpa(idx + 1) + 4.*del
            wgts_vpa(idx + 2) = wgts_vpa(idx + 2) + del
         end do

         ! For the sake of symmetry, do the same thing with 3/8 rule at
         ! upper boundary, and composite elsewhere.
         del = 0.375 * dvpa
         wgts_vpa(nvpa - 3) = wgts_vpa(nvpa - 3) + del
         wgts_vpa(nvpa - 2:nvpa - 1) = wgts_vpa(nvpa - 2:nvpa - 1) + 3.*del
         wgts_vpa(nvpa) = wgts_vpa(nvpa) + del
         nvpa_seg = (nvpa - 4) / 2
         del = dvpa / 3.
         do iseg = 1, nvpa_seg
            idx = 2 * (iseg - 1) + 1
            wgts_vpa(idx) = wgts_vpa(idx) + del
            wgts_vpa(idx + 1) = wgts_vpa(idx + 1) + 4.*del
            wgts_vpa(idx + 2) = wgts_vpa(idx + 2) + del
         end do

         ! Divide by 2 to account for double-counting
         wgts_vpa = 0.5 * wgts_vpa
         
         ! Add the 1/sqrt(pi) factor of velocity integrations to <wgts_vpa>
         wgts_vpa = wgts_vpa / sqrt(pi)

         ! If maxwellian_normalization = .true., then the evolved pdf is normalized
         ! by a Maxwellian; this normalisation must be accounted for in the velocity
         ! space integrals, so include exp(-vpa^2) factor in the vpa weights.
         ! NB: the species index of maxwell_vpa is not needed for the radially local version 
         ! of the code and would otherwise add a species index to wgts_vpa, so currently 
         ! maxwellian_normalization is not supported for the radially global version of the code.
         if (maxwellian_normalization) wgts_vpa = wgts_vpa * maxwell_vpa(:, 1)
         
      end subroutine initialise_vpa_integration_weights_Simpsons

   end subroutine init_vpa_grid


   !****************************************************************************
   !                Initialise perpendicular (mu) velocity grid                !
   !****************************************************************************
   ! Integrations over velocity space always have the following form:
   !     int dmu int dvpa 2 bmag(z) / sqrt(pi)
   ! 
   ! We will absorb the factor [ 2 bmag(z) ] into the mu integration weights:
   !     <wgts_mu>[ialpha, iz, imu] = <wgts_mu_bare>[imu] * 2 * bmag[iz]
   !****************************************************************************
   subroutine init_mu_grid
   
      ! Flags
      use parameters_numerical, only: maxwellian_normalization
      
      ! Geometry
      use geometry, only: bmag
      
      ! Grids
      use grids_z, only: nztot
      use grids_kxky, only: nalpha
      
      implicit none
      
      !-------------------------------------------------------------------------
      
      ! Allocate arrays related to the mu-grid
      call allocate_arrays
      
      ! Initialise the <mu>-grid and the integration weights called <wgts_mu_bare>
      ! Either the mu-grid is equally spaced or a Gauss-Laguerre quadrature is used
      if (equally_spaced_mu_grid) then
         call initialise_mu_and_integration_weights_equally_spaced
      else
         call initialise_mu_and_integration_weights_gauss_laguerre
      end if

      ! Add the constant factor [ 2 bmag(z) ] of the velocity integration to the mu integration weights
      wgts_mu = 2.*spread(spread(wgts_mu_bare, 1, nalpha), 2, nztot) * spread(bmag, 3, nmu)

      ! If <maxwellian_normalization>, the evolved pdf is normalized by a Maxwwellian;
      ! in this case, the velocity integration must account for the Maxwellian.
      ! NB: the species index on maxwell_mu is only needed for radially global simulations,
      ! which are not currently supported for maxwellian_normalization = .true.
      if (maxwellian_normalization) wgts_mu = wgts_mu * maxwell_mu(:, :, :, 1)
      
   contains
   
      !---------------------------- Allocate arrays ----------------------------
      subroutine allocate_arrays
      
         use grids_species, only: nspec
         use grids_z, only: nzgrid
      
         implicit none

         if (.not. allocated(mu)) then
            allocate (mu(nmu)); mu = 0.0
            allocate (wgts_mu(nalpha, -nzgrid:nzgrid, nmu)); wgts_mu = 0.0
            allocate (wgts_mu_bare(nmu)); wgts_mu_bare = 0.0
            allocate (maxwell_mu(nalpha, -nzgrid:nzgrid, nmu, nspec)); maxwell_mu = 0.0
            allocate (maxwell_mu_avg(nalpha, -nzgrid:nzgrid, nmu, nspec)); maxwell_mu_avg = 0.0
            allocate (dmu(nmu - 1))
            allocate (mu_cell(nmu))
            allocate (dmu_cell(nmu))
         end if
         
      end subroutine allocate_arrays
       
      !-------------------------------------------------------------------------
      !------------ Calculate mu integration weights: equally spaced -----------
      !-------------------------------------------------------------------------
      ! Get an equally space grid in mu with max value
      !     <mu_max> = vperp_max**2/(2*max(bmag))
      ! Put first grid point at dmu/2 to avoid mu=0 special point
      !     dmu/2 + (nmu-1)*dmu = mu_max
      !     so dmu = mu_max/(nmu-1/2)
      !-------------------------------------------------------------------------
      subroutine initialise_mu_and_integration_weights_equally_spaced
         
         use geometry, only: bmag_psi0
         
         implicit none
         
         integer :: imu
         real :: mu_max
      
         !----------------------------------------------------------------------
         
         mu_max = vperp_max**2 / (2.*maxval(bmag_psi0))
         dmu = mu_max / (nmu - 0.5)
         mu(1) = 0.5 * dmu(1)
         do imu = 2, nmu
            mu(imu) = mu(1) + (imu - 1) * dmu(1)
         end do
         wgts_mu_bare = dmu(1)
      
      end subroutine initialise_mu_and_integration_weights_equally_spaced
       
      !-------------------------------------------------------------------------
      !------------ Calculate mu integration weights: Gauss-Laguerre -----------
      !-------------------------------------------------------------------------
      subroutine initialise_mu_and_integration_weights_gauss_laguerre

         use gauss_quad, only: get_laguerre_grids
         use geometry, only: bmag_psi0
         
         implicit none
      
         !----------------------------------------------------------------------
         
         real, dimension(:), allocatable :: wgts_GL
         real, dimension(:), allocatable :: x_GL
         real :: mu_max, x_max
         
         ! Allocate temporary arrays
         allocate (x_GL(nmu)); x_GL = 0.0
         allocate (wgts_GL(nmu)); wgts_GL = 0.0
         
         ! Calculate a Gauss-Laguerre quadrature with nmu points
         call get_laguerre_grids(x_GL, wgts_GL)
         
         ! Get the maximum value of the x-points of the Gauss-Laguerre quadrature
         x_max = x_GL(nmu)
         
         ! Make sure <vperp_max> is defined
         if (vperp_max < 0) vperp_max = sqrt(x_max)
         
         ! The mu grid ranges from mu_min to mu_max with
         mu_max = vperp_max**2 / ( 2.*minval(bmag_psi0) )
      
         ! Rescale the Gauss-Laguerre weights to <mu_max>
         ! Moreover, the quadrature has a factor exp(-x_GL) which we need to remove
         wgts_mu_bare = wgts_GL * exp(x_GL) / x_max * mu_max

         ! Rescale the x-grid from the Gauss-Laguerre quadrature to a mu-grid
         mu = x_GL / x_max * mu_max

         ! Calculate the step size between mu points
         dmu(:nmu - 1) = mu(2:) - mu(:nmu - 1)
         
         ! Note that we leave dmu(nmu) uninitialized. It should never be used,
         ! so we want valgrind or similar to return error if it is used.
      
         ! Deallocate temporary arrays
         deallocate (x_GL)
         deallocate (wgts_GL)
      
      end subroutine initialise_mu_and_integration_weights_gauss_laguerre

   end subroutine init_mu_grid
   
   
   !****************************************************************************
   !                       Initialise Maxwellian factors                       !
   !****************************************************************************
   subroutine init_maxwellian_factors
      
      use parameters_physics, only: full_flux_surface
      use grids_species, only: spec, nspec
      use grids_species, only: species_option_switch
      use grids_species, only: species_option_multibox
      use grids_z, only: nztot
      use grids_kxky, only: nalpha
      use geometry, only: bmag
   
      implicit none
      
      !----------------------------------------------------------------------

      ! Calculate the vpa part of the v-space Maxwellian
      maxwell_vpa = exp(-spread(vpa * vpa, 2, nspec))
      
      ! Allow for radial variation in the species parameters
      if (species_option_switch == species_option_multibox) then
         maxwell_vpa = exp(-spread(vpa * vpa, 2, nspec) * spread(spec%temp_psi0 / spec%temp, 1, nvpa))
      end if
      
      ! <ztmax> is the Maxwellian in <vpa>, multiplied by charge number over normalized temperature
      ztmax = spread(spec%zt, 1, nvpa) * maxwell_vpa
   
      ! Calculate the mu part of the v-space Maxwellian
      maxwell_mu = exp(-2.*spread(spread(spread(mu, 1, nalpha), 2, nztot) * spread(bmag, 3, nmu), 4, nspec) &
         * spread(spread(spread(spec%temp_psi0 / spec%temp, 1, nalpha), 2, nztot), 3, nmu))
      if (full_flux_surface) maxwell_mu_avg = spread(sum(maxwell_mu, dim = 1), 1, nalpha) / nalpha

      ! Calculate <maxwell_fac>. Note that <maxwell_fac> = 1 unless radially global
      if (.not. allocated(maxwell_fac)) allocate (maxwell_fac(nspec)); maxwell_fac = 1.0
      maxwell_fac = spec%dens / spec%dens_psi0 * (spec%temp_psi0 / spec%temp)**1.5
   
   end subroutine init_maxwellian_factors


   !****************************************************************************
   !           Initialise velocity parameters for collision operators          !
   !****************************************************************************
   subroutine init_velocity_parameters_for_collision_operators
   
      implicit none
      
      real, dimension(:), allocatable :: dmu_ghost
      
      !-------------------------------------------------------------------------
      
      ! Allocate temporary array
      allocate (dmu_ghost(nmu))

      ! For the collision operators we need <mu_cell> and <dmu_cell>
      ! Add ghost cell at mu=0 and beyond mu_max for purposes of differentiation
      ! note assuming here that grid spacing for ghost cell is equal to
      ! grid spacing for last non-ghost cell
      dmu_ghost(:nmu - 1) = dmu; dmu_ghost(nmu) = dmu(nmu - 1)
      
      ! This is mu at cell centres (including to left and right of mu grid boundary points)
      mu_cell(:nmu - 1) = 0.5 * (mu(:nmu - 1) + mu(2:))
      mu_cell(nmu) = mu(nmu) + 0.5 * dmu(nmu - 1)
      
      ! This is mu_{j+1/2} - mu_{j-1/2}
      dmu_cell(1) = mu_cell(1)
      dmu_cell(2:) = mu_cell(2:) - mu_cell(:nmu - 1)
      
      ! Deallocate temporary array
      deallocate (dmu_ghost)
   
   end subroutine init_velocity_parameters_for_collision_operators


   !****************************************************************************
   !               Redefine parallel velocity integration weights              !
   !****************************************************************************
   ! The collision module is able to change the velocity weights.
   ! Added option for density conserving form of collision operator.d
   !****************************************************************************
   subroutine set_vpa_weights(conservative)

      use constants, only: pi

      implicit none

      logical, intent(in) :: conservative
      
      !-------------------------------------------------------------------------

      if (conservative) then
         wgts_vpa = dvpa / sqrt(pi)
      else if (conservative_wgts_vpa) then 
         wgts_vpa = dvpa / sqrt(pi)
      else if ((.not. conservative_wgts_vpa) .and. (.not. conservative)) then
         wgts_vpa = wgts_vpa_default
      end if

   end subroutine set_vpa_weights
   
   !******************************************************************************
   !                           FINALISE VELOCITY GRIDS
   !******************************************************************************
   subroutine finish_velocity_grids

      implicit none

      call finish_vpa_grid
      call finish_mu_grid
      call finish_maxwellians
      call finish_arrays_for_collision_operator

      initialised_velocity_grids = .false.

   end subroutine finish_velocity_grids

   subroutine finish_vpa_grid

      implicit none

      if (allocated(vpa)) deallocate (vpa)
      if (allocated(wgts_vpa)) deallocate (wgts_vpa)
      if (allocated(wgts_vpa_default)) deallocate (wgts_vpa_default)
      if (allocated(ztmax)) deallocate (ztmax)

   end subroutine finish_vpa_grid

   subroutine finish_mu_grid

      implicit none

      if (allocated(mu)) deallocate (mu)
      if (allocated(wgts_mu)) deallocate (wgts_mu)
      if (allocated(wgts_mu_bare)) deallocate (wgts_mu_bare)
      if (allocated(dmu)) deallocate (dmu)

   end subroutine finish_mu_grid

   subroutine finish_maxwellians

      implicit none
      
      if (allocated(maxwell_mu)) deallocate (maxwell_mu)
      if (allocated(maxwell_mu_avg)) deallocate (maxwell_mu_avg)
      if (allocated(maxwell_vpa)) deallocate (maxwell_vpa)
      if (allocated(maxwell_fac)) deallocate (maxwell_fac)

   end subroutine finish_maxwellians

   subroutine finish_arrays_for_collision_operator

      implicit none

      if (allocated(int_unit)) deallocate (int_unit)
      if (allocated(int_vpa2)) deallocate (int_vpa2)
      if (allocated(int_vperp2)) deallocate (int_vperp2)
      if (allocated(int_vfrth)) deallocate (int_vfrth)
      if (allocated(mu_cell)) deallocate (mu_cell)
      if (allocated(dmu_cell)) deallocate (dmu_cell)

   end subroutine finish_arrays_for_collision_operator

end module grids_velocity
