!###############################################################################
!                                 EXTENDED Z-GRID                               
!###############################################################################
! 
! In the (kx,ky) directions, periodic boundary conditions are assumed.
! This automatically statisfied thanks to the Fourier decomposition along (x,y).
!     phi(Lx,y,z) = phi(-Lx,y,z)
!     phi(x,Ly,z) = phi(x,-Ly,z)
! 
!                      TWIST AND SHIFT BOUNDARY CONDITIONS
! 
! In the parallel direction, along the field lines, or in the z-direction, the 
! (kx,ky) grid at one end of the flux tube needs to be connected to the (kx,ky)
! grid on the other end of the flux tube. Since the (kx,ky) directions are not
! always perpendicular, i.e., their angle changes along the field line due to 
! the local magnetic shear. It is possible that the (kx,ky) grid do not match
! at the ends of the flux tube. To deal with the effect of shear we can apply:
!     - Standard twist-and-shift boundary conditions 
!     - Stellarator-symmetric twist-and-shift boundary conditions based on the local magnetic shear
! 
! The standard twist-and-shift boundary conditions calculated how the (kx,ky) grid
! is connected at the flux tube ends based on he global magnetic shear <shat>. In
! contrast the stellarator-symmetric twist-and-shift boundary conditions calculate
! conntections based on the local magnetic shear. The standard twist-and-shift 
! boundary conditions are only valid in axysimmetric devices such as tokamaks.
! It is also commonly applied to stellarators, even they they are formally
! incorrect in stellarators. In stellarators once has to use the Stellarator-symmetric
! twist-and-shift boundary conditions. Moreover, in low-shear devices they can be 
! a lot more efficient than the Standard twist-and-shift boundary conditions.
! 
! The twist and shift boundary conditions create chains of connected (kx,ky) modes.
! Since these modes are connected along the z-direction, they form a chain of connected
! 2*pi segments on a so-called extended z-grid. Each chain is a so-called eigenmode
! of the system, and the number of eigenmodes is stored in <neigen>. These chains
! conntect different (kx,ky) modes with the same ky-value but a different kx-value.
! Therefore, for each ky-value we will have <neigen(ky)> eigenmodes. The number of 
! connections goes down with ky. So the larger ky, the less connections are made and
! the more eigenmodes are present. At the lowest ky-value, typically all kx-modes
! are connected to each other and we have <neigen> = 1 (if <jtwist> = 1).
! 
! The zonal modes (modes with ky=0) are always connected onto themselves. So
! <neigen> = <nakx> for ky=0. Also for large ky values each (kx,ky) mode
! will be connected onto themselves, since they can not make connections.
! 
!                         PERIODIC BOUNDARY CONDITIONS
! 
! If the global magnetic shear is very low, or the (kx,ky) directions are the flux
! tube end are approximatly perpendicular, one can apply:
!     - Periodic boundary conditions
! which connects each (kx,ky) mode onto itself.
! 
!###############################################################################
module grids_extended_zgrid

   ! Load debug flags
   use debug_flags, only: debug => extended_grid_debug

   implicit none

   ! Make routines and parameters available to other modules
   public :: nsegments
   public :: neigen, neigen_max
   public :: ikxmod
   public :: iz_low, iz_mid, iz_up
   public :: it_right
   public :: periodic
   public :: phase_shift
   public :: nzed_segment
   public :: fill_zed_ghost_zones, fill_zext_ghost_zones
   public :: init_extended_zgrid, finish_extended_zgrid
   public :: map_to_extended_zgrid
   public :: map_from_extended_zgrid
   public :: map_to_iz_ikx_from_izext
   public :: nzext, nzext_max
   public :: izext_iz, izext_ikx, izext_dtube
   public :: ndup, dup_iz, dup_ikx, dup_dtube, dup_src
   public :: tube_index
   public :: get_chain_with_ghosts, enforce_chain_copies
   public :: zed_operator_along_chains
   public :: zop_upwind3, zop_centered2, zop_fd_cell_centers, zop_cell_centers, zop_fourth

   private

   ! Chains of modes are contructed, where modes are connected along z that have 
   ! the same ky-value but a different kx-value. The following arrays keep track
   !  of  connections between different 2pi segments of the z-grid.
   integer :: nzed_segment, neigen_max
   integer, dimension(:), allocatable :: neigen
   integer, dimension(:), allocatable :: iz_low, iz_mid, iz_up
   integer, dimension(:, :), allocatable :: nsegments
   integer, dimension(:, :, :), allocatable :: ikxmod
   
   ! Arrays to indicate which flux tube index to connect to on the left and on the right,
   ! as a function of current flux tube index. Pre-compute to avoid conditionals in loops
   integer, dimension(:), allocatable :: it_left, it_right
   complex, dimension(:), allocatable :: phase_shift

   ! Keep trach whether a ky-mode is periodic or not
   logical, dimension(:), allocatable :: periodic

   ! Each chain (ie, iky) is also described point by point on the extended z-grid, which is used by
   ! the implicit solver and the response matrix. This description is general: the (iz, ikx) location
   ! may change at any point of the chain, not only between 2pi segments, as needed for the
   ! non-twisting flux tube. For the twist-and-shift chains it is built from the segments above.
   !     <nzext>(ie, iky)                          number of points of the chain on the extended z-grid
   !     izext_iz, izext_ikx(izext, ie, iky)       (iz, ikx) location of each point of the chain
   !     izext_dtube(izext, ie, iky)               number of flux tubes to the right (flux tube train)
   !     phase factors                             gext = g * phase_to and g = gext * phase_from
   ! Grid points that are copies of a point of the chain (the shared end points of 2pi segments)
   ! are listed separately, they are filled by <map_from_extended_zgrid>:
   !     <ndup>(ie, iky)                           number of copies
   !     dup_iz, dup_ikx, dup_dtube(idup, ie, iky) location of each copy
   !     dup_src(idup, ie, iky)                    point of the chain (izext) that is copied
   ! The points of all chains are stored one after the other in 1D arrays (a chain of a low ky can
   ! be much longer than one of a high ky), starting after <point_offset>(ie, iky) and <dup_offset>(ie, iky).
   integer :: nzext_max
   integer, dimension(:, :), allocatable :: nzext, ndup, point_offset, dup_offset
   integer, dimension(:), allocatable :: point_iz, point_ikx, point_dtube
   integer, dimension(:), allocatable :: copy_iz, copy_ikx, copy_dtube, copy_src
   complex, dimension(:), allocatable :: point_phase_to, point_phase_from, copy_phase_from

   ! Only initialise once
   logical :: initialised_extended_zgrid = .false.

   ! Finite-difference operators along z that <zed_operator_along_chains> can apply
   integer, parameter :: zop_upwind3 = 1, zop_centered2 = 2, zop_fd_cell_centers = 3
   integer, parameter :: zop_cell_centers = 4, zop_fourth = 5

contains

!###############################################################################
!######################### INITIALISE EXTENDED Z-GRID ##########################
!###############################################################################

   subroutine init_extended_zgrid

      use grids_z, only: boundary_option_switch
      use grids_z, only: boundary_option_self_periodic
      use grids_z, only: boundary_option_zero
      use grids_z, only: boundary_option_linked
      use grids_z, only: boundary_option_linked_stellarator
      use grids_z, only: boundary_option_non_twisting
      use grids_z, only: nperiod, nzgrid, nzed, ntubes
      use grids_kxky, only: nakx, naky, ikx_max
      use grids_kxky, only: jtwist, ikx_twist_shift, phase_shift_angle
      use grids_kxky, only: aky
      use constants, only: zi

      implicit none

      integer :: iseg, iky, ie, ikx, it
      integer :: nseg_max
      integer, dimension(:), allocatable :: ikx_shift_end
      integer, dimension(:, :), allocatable :: ikx_shift

      if (initialised_extended_zgrid) return
      initialised_extended_zgrid = .true.
 
      if (debug) write (*, *) 'extended_zgrid::allocate_arrays'  
      if (.not. allocated(neigen)) allocate (neigen(naky))
      if (.not. allocated(periodic)) allocate (periodic(naky)); periodic = .false.
      if (.not. allocated(phase_shift)) allocate (phase_shift(naky))

      ! Periodic boundary condition along z
      if (boundary_option_switch == boundary_option_self_periodic) then
         periodic = .true.
         
      ! Zonal modes (modes with ky=0) are also periodic, i.e., they connect onto themselves
      else
         where (abs(aky) < epsilon(0.0)) periodic = .true.
      end if

      ! phase shift due to the twist-and-shift boundary condition
      ! Usually set to zero for standard local simulation, but can
      ! have an effect for global simulations and simulations with low
      ! magnetic shear that use periodic boundary conditions everywhere
      phase_shift = exp(zi * aky * phase_shift_angle)

      ! For unconnected boundary conditions ('default', 'zero' or 'unconnected') the
      ! 2pi segments tile one contiguous field line at fixed kx, so the segment joins
      ! are interior points (z = +-pi, +-3pi, ...) and no phase shift belongs there.
      ! Applying it would corrupt the parallel derivative at these joins whenever
      ! rhostar > 0, since <phase_shift_angle> is then set from rhostar in grids_kxky.
      if (boundary_option_switch == boundary_option_zero) phase_shift = (1.0, 0.0)

      ! The non-twisting flux tube builds its chains point by point
      if (boundary_option_switch == boundary_option_non_twisting) then
         ! The phase shift of global simulations is not implemented for the non-twisting flux tube
         phase_shift = (1.0, 0.0)
         call init_extended_zgrid_non_twisting
         call init_tube_connections
         return
      end if

      if (debug) write (*, *) 'extended_zgrid::boundary_option_switch'  
      if (boundary_option_switch == boundary_option_linked .or. boundary_option_switch == boundary_option_linked_stellarator) then

         ! All periodic modes (e.g., the zonal mode) have no connections
         do iky = 1, naky
            if (periodic(iky)) then
               neigen(iky) = nakx
            else
               ! TODO: check whether range grids with linked boundary conditions are correct.
               ! For the range grid jtwist = 0 (see read_parameters_kxky_grids in grids_kxky.f90),
               ! hence neigen = 0 for all ky > 0, i.e., no mode is connected to any other.
               neigen(iky) = min((iky - 1) * jtwist, nakx)
            end if
         end do

         neigen_max = maxval(neigen)

         if (.not. allocated(ikx_shift_end)) then
            allocate (ikx_shift_end(neigen_max)); ikx_shift_end = 0
            allocate (ikx_shift(nakx, naky)); ikx_shift = 0
         end if

         ! phi(kx-kx_shift,-nzgrid) = phi(kx,nzgrid) from twist-and-shift BC
         ! for positive (negative) magnetic shear, kx_shift is positive (negative),
         ! so start at most positive (negative) kx and
         ! progress to smaller (larger) kx values as connections occur

         ! figure out how much to shift ikx by to get to the end of the kx chain
         ! for positive (negative) magnetic shear, this is the left-most (right-most) theta-theta0
         ! in each set of connected 2pi segments
         ! note that theta0 goes from 0 to theta0_max and then from theta0_min back
         ! to -dtheta0

         do ikx = 1, neigen_max
            ! first ikx_max=nakx/2+1 theta0s are 0 and all positive theta0 values
            ! remainder are negative theta0s
            ! theta_0 = kx / ky / shat
            ! if ky > 0, then most positive theta_0 corresponds to most positive kx

            ! first consider case where shift in kx is negative (corresponds to positive magnetic shear)
            if (ikx_twist_shift < 0) then
               if (ikx <= ikx_max) then
                  ikx_shift_end(ikx) = ikx_max - 2 * ikx + 1
               else
                  ikx_shift_end(ikx) = 3 * ikx_max - 2 * ikx
               end if
               ! then consider case where shift in kx is positive
            else if (ikx_twist_shift > 0) then
               if (ikx < ikx_max) then
                  if (ikx + ikx_max <= nakx) then
                     ikx_shift_end(ikx) = ikx_max
                  else
                     ikx_shift_end(ikx) = ikx - nakx
                  end if
               else
                  ikx_shift_end(ikx) = 1 - ikx_max
               end if
            end if
            ! note that zero shift case is taken care of by initialization of ikx_shift_end
         end do

         do iky = 1, naky
            ! ikx_shift is how much to shift each ikx by to connect
            ! to the next theta0 (from most positive to most negative for positive magnetic shear
            ! and vice versa for negative magnetic shear)

            ! first consider shift in index for case where shift is negative
            ! (corresponds to positive magnetic shear)
            if (ikx_twist_shift < 0) then
               ! if ky > 0, then going to more negative theta0
               ! corresponds to going to more negative kx
               do ikx = 1, ikx_max
                  ! if theta0 is sufficiently positive, shifting to more
                  ! negative theta0 corresponds to decreasing ikx
                  if (ikx - neigen(iky) > 0) then
                     ikx_shift(ikx, iky) = -neigen(iky)
                     ! if a positive theta0 connects to a negative theta0
                     ! must do more complicated mapping of ikx
                  else if (ikx - neigen(iky) + nakx >= ikx_max + 1) then
                     ikx_shift(ikx, iky) = nakx - neigen(iky)
                  end if
               end do
               ! if theta0 is negative, then shifting to more negative
               ! theta0 corresponds to decreasing ikx
               do ikx = ikx_max + 1, nakx
                  ! if theta0 is sufficiently negative, it has no
                  ! more negative theta0 with which it can connect
                  if (ikx - neigen(iky) > ikx_max) then
                     ikx_shift(ikx, iky) = -neigen(iky)
                  end if
                  ! theta0 is positive
               end do
            else if (ikx_twist_shift > 0) then
               ! if ky > 0, then going to more positive theta0
               ! corresponds to going to more positive kx
               do ikx = 1, ikx_max
                  ! if shift in kx, kx_shift, is less than kx-kx_max,
                  ! then shift by the appropriate amount
                  if (ikx + neigen(iky) <= ikx_max) then
                     ikx_shift(ikx, iky) = neigen(iky)
                  end if
                  ! otherwise, no kx on grid to connect with
               end do
               do ikx = ikx_max + 1, nakx
                  ! if kx+kx_shift < 0, then simple shift by neigen
                  if (ikx + neigen(iky) <= nakx) then
                     ikx_shift(ikx, iky) = neigen(iky)
                     ! if 0 < kx+kx_shift <= kx_max, then more complicated shift
                     ! to positive set of kx values
                  else if (ikx - ikx_max + neigen(iky) <= nakx) then
                     ikx_shift(ikx, iky) = neigen(iky) - nakx
                  end if
                  ! otherwise, no kx on grid with which to connect
               end do
            end if
         end do

         if (.not. allocated(nsegments)) allocate (nsegments(neigen_max, naky))

         do iky = 1, naky
            if (neigen(iky) == 0) then
               nsegments(:, iky) = 1
            else
               nsegments(:, iky) = (nakx - 1) / neigen(iky)

               do ie = 1, mod(nakx - 1, neigen(iky)) + 1
                  nsegments(ie, iky) = nsegments(ie, iky) + 1
               end do
            end if
         end do

         nseg_max = maxval(nsegments) 
         if ((.not. allocated(iz_low)) .and. (neigen_max>0)) then
            allocate (iz_low(nseg_max)); iz_low = -nzgrid
            allocate (iz_mid(nseg_max)); iz_mid = 0
            allocate (iz_up(nseg_max)); iz_up = nzgrid
         end if

      else

         neigen = nakx; neigen_max = nakx

         if (.not. allocated(ikx_shift_end)) then
            allocate (ikx_shift_end(neigen_max))
            allocate (ikx_shift(nakx, naky))
         end if
         ikx_shift = 0; ikx_shift_end = 0

         if (.not. allocated(nsegments)) then
            allocate (nsegments(neigen_max, naky))
         end if

         ! this is the number of 2pi poloidal segments in the extended theta domain,
         ! which is needed in initializing the reponse matrix and doing the implicit sweep
         nsegments = 2 * (nperiod - 1) + 1

         nseg_max = maxval(nsegments)

         if (.not. allocated(iz_low)) then
            allocate (iz_low(nseg_max))
            allocate (iz_mid(nseg_max))
            allocate (iz_up(nseg_max))
         end if

         ! iz_low(j) is the ig index corresponding to the inboard midplane from below (theta=-pi) within the jth segment
         ! iz_mid(j) is the ig index corresponding to the outboard midplane (theta=0) within the jth segment
         do iseg = 1, nseg_max
            iz_low(iseg) = -nzgrid + (iseg - 1) * nzed
            iz_mid(iseg) = iz_low(iseg) + nzed / 2
            iz_up(iseg) = iz_low(iseg) + nzed
         end do

      end if

      if (debug) write (*, *) 'extended_zgrid::ikxmod_1'  
      if (.not. allocated(ikxmod)) then
         allocate (ikxmod(nseg_max, neigen_max, naky))
         ! initialize ikxmod to nakx
         ! should not be necessary but just in case one tries to access
         ! a value beyond nsegments(ie,iky)
         ikxmod = nakx
      end if

      if (debug) write (*, *) 'extended_zgrid::ikxmod_2'  
      do iky = 1, naky
         ! only do the following once for each independent set of theta0s
         ! the assumption here is that all kx are on processor and sequential
         do ie = 1, neigen(iky)
            ! remap to start at theta0 = theta0_max (theta0_min) for negative (positive) kx shift
            ! for this set of connected theta0s
            iseg = 1
            ikxmod(iseg, ie, iky) = ie + ikx_shift_end(ie)
            if (nsegments(ie, iky) > 1) then
               do iseg = 2, nsegments(ie, iky)
                  ikxmod(iseg, ie, iky) = ikxmod(iseg - 1, ie, iky) + ikx_shift(ikxmod(iseg - 1, ie, iky), iky)
               end do
            end if
         end do
      end do

      if (debug) write (*, *) 'extended_zgrid::deallocate'  
      if (allocated(ikx_shift_end)) deallocate (ikx_shift_end)
      if (allocated(ikx_shift)) deallocate (ikx_shift)

      call init_tube_connections

      ! this is the number of unique zed values in all segments but the first
      ! the first has one extra unique zed value (all others have one grid common
      ! with the previous segment due to periodicity)
      if (neigen_max>0) then 
         nzed_segment = iz_up(1) - iz_low(1)
      else 
         nzed_segment = 0
      end if 

      ! Describe each chain point by point on the extended z-grid
      call build_chain_points_from_segments
      
   end subroutine init_extended_zgrid

   !****************************************************************************
   !         Flux tube connections (flux tube train with ntubes > 1)
   !****************************************************************************
   subroutine init_tube_connections

      use grids_z, only: ntubes

      implicit none

      integer :: it

      if (.not. allocated(it_left)) allocate (it_left(ntubes))
      if (.not. allocated(it_right)) allocate (it_right(ntubes))

      it_right(ntubes) = 1
      if (ntubes > 1) then
         do it = 1, ntubes - 1
            it_right(it) = it + 1
         end do
      end if

      it_left(1) = ntubes
      if (ntubes > 1) then
         do it = 2, ntubes
            it_left(it) = it - 1
         end do
      end if

   end subroutine init_tube_connections

   !****************************************************************************
   !                 Chains of the non-twisting flux tube
   !****************************************************************************
   ! In the non-twisting flux tube (Ball and Brunner 2021) the kx index of a field
   ! line changes along z: the field line keeps the ballooning wavenumber
   ! akx + m0_ntft(ky, z) dkx, i.e. the radial mode number m + m0 (eq. B.12), so from
   ! z(iz) to z(iz+1) its grid mode number changes from m to m - [m0(iz+1) - m0(iz)].
   ! It leaves the grid (zero beyond, as in the conventional flux tube) where |m| > m_max,
   ! which can happen anywhere along z. At the end of the tube z(nzgrid) is the same
   ! point as z(-nzgrid) with the same grid index (eq. B.15), so a chain continues there
   ! in the next flux tube of the train. The unique points of a ky are (iz, m) with
   ! iz = -nzgrid, ..., nzgrid-1; the grid points at iz = nzgrid are copies.
   !
   ! A chain starts at a point without a predecessor on the grid and follows its
   ! successors until it leaves the grid. On the extended z-grid it is stored as the
   ! sequence of its points; when it passes the end of the tube the point is stored
   ! at iz = nzgrid (so that the extended grid is equally spaced) and its copy at
   ! iz = -nzgrid of the next tube. Without twist (n_twist_ntft = 0, e.g. ky = 0)
   ! every field line closes on itself after one pass: such a ky is periodic, and
   ! its chains are stored like the periodic chains of the conventional flux tube.
   !****************************************************************************
   subroutine init_extended_zgrid_non_twisting

      use mp, only: mp_abort, proc0
      use grids_z, only: nzgrid
      use grids_kxky, only: naky, nakx, ikx_max
      use grids_kxky, only: m0_ntft, n_twist_ntft

      implicit none

      integer :: iky, ie, mmax, nchains, max_chains
      integer, dimension(:), allocatable :: dm
      logical, dimension(:, :), allocatable :: visited
      integer, dimension(:, :, :), allocatable :: starts
      ! State of the chain that is being walked: number of points and copies, and whether to store them
      integer :: np, nc
      logical :: store

      !-------------------------------------------------------------------------

      ! Radial mode numbers on the grid: m = -mmax, ..., mmax
      mmax = ikx_max - 1
      if (nakx /= 2 * mmax + 1) call mp_abort('The non-twisting flux tube needs a symmetric kx grid. Aborting.')

      allocate (dm(-nzgrid:nzgrid - 1))
      allocate (visited(-nzgrid:nzgrid - 1, -mmax:mmax))
      allocate (starts(2, nakx * 2 * nzgrid, naky))
      if (.not. allocated(neigen)) allocate (neigen(naky))

      ! First pass: find the chains of every ky (their first point) and whether the ky is periodic
      do iky = 1, naky
         call set_dm(iky)
         periodic(iky) = n_twist_ntft(iky) == 0
         visited = .false.; nchains = 0
         call find_chains(iky)
         neigen(iky) = nchains
      end do
      neigen_max = maxval(neigen)
      max_chains = max(neigen_max, 1)

      ! The segment description is not used by the non-twisting flux tube, but keep it allocated:
      ! one 'segment' per chain, whose kx is the kx of the first point of the chain
      if (.not. allocated(nsegments)) allocate (nsegments(max_chains, naky)); nsegments = 1
      if (.not. allocated(ikxmod)) allocate (ikxmod(1, max_chains, naky)); ikxmod = 1
      if (.not. allocated(iz_low)) allocate (iz_low(1), iz_mid(1), iz_up(1))
      iz_low = -nzgrid; iz_mid = 0; iz_up = nzgrid
      nzed_segment = 2 * nzgrid

      ! Second pass: count the points of every chain, allocate, and store them
      call allocate_chain_offsets
      store = .false.
      do iky = 1, naky
         call set_dm(iky)
         do ie = 1, neigen(iky)
            call walk_chain(iky, ie, starts(1, ie, iky), starts(2, ie, iky))
            nzext(ie, iky) = np
            ndup(ie, iky) = nc
         end do
      end do
      call allocate_chain_points
      store = .true.
      do iky = 1, naky
         call set_dm(iky)
         do ie = 1, neigen(iky)
            call walk_chain(iky, ie, starts(1, ie, iky), starts(2, ie, iky))
            ikxmod(1, ie, iky) = izext_ikx(1, ie, iky)
         end do
      end do

      deallocate (dm, visited, starts)

      ! Summary of the chains
      if (proc0) then
         write (*, '(A)') ' Non-twisting flux tube (Ball and Brunner 2021): chains along z for every ky'
         write (*, '(A)') '    iky   n_twist   chains   periodic   min/max points per chain   max |m0|'
         do iky = 1, naky
            write (*, '(I7, I10, I9, L11, I14, A, I7, I13)') iky, n_twist_ntft(iky), neigen(iky), periodic(iky), &
               minval(nzext(:neigen(iky), iky)), ' /', maxval(nzext(:neigen(iky), iky)), maxval(abs(m0_ntft(iky, :)))
         end do
      end if

   contains

      ! Change of the grid mode number of a field line from z(iz) to z(iz+1): m -> m - dm(iz)
      subroutine set_dm(iky)
         integer, intent(in) :: iky
         dm = m0_ntft(iky, -nzgrid + 1:nzgrid) - m0_ntft(iky, -nzgrid:nzgrid - 1)
      end subroutine set_dm

      ! kx index of the radial mode number m
      integer function ikx_of(m)
         integer, intent(in) :: m
         if (m >= 0) then
            ikx_of = m + 1
         else
            ikx_of = nakx + m + 1
         end if
      end function ikx_of

      ! Find the first point of every chain of this ky
      subroutine find_chains(iky)

         integer, intent(in) :: iky
         integer :: iz, m, izp, mp

         ! Open chains start at a point whose predecessor is not on the grid
         do m = -mmax, mmax
            do iz = -nzgrid, nzgrid - 1
               if (iz > -nzgrid) then
                  izp = iz - 1
               else
                  izp = nzgrid - 1
               end if
               mp = m + dm(izp)
               if (abs(mp) > mmax) then
                  if (periodic(iky)) call mp_abort('Non-twisting flux tube: a field line leaves the kx grid &
                     &although the tube has no twist; increase nx. Aborting.')
                  nchains = nchains + 1
                  starts(:, nchains, iky) = (/iz, m/)
                  call mark_chain(iky, iz, m)
               end if
            end do
         end do
         ! Without twist, the remaining points lie on closed field lines; start them at iz = -nzgrid
         do m = -mmax, mmax
            if (visited(-nzgrid, m)) cycle
            if (.not. periodic(iky)) call mp_abort('Non-twisting flux tube: closed field line in a twisted tube. Aborting.')
            nchains = nchains + 1
            starts(:, nchains, iky) = (/-nzgrid, m/)
            call mark_chain(iky, -nzgrid, m)
         end do
         ! Every point must belong to one chain
         if (.not. all(visited)) call mp_abort('Non-twisting flux tube: not every point is on a chain. Aborting.')

      end subroutine find_chains

      ! Mark the points of the chain that starts at (iz_start, m_start)
      subroutine mark_chain(iky, iz_start, m_start)

         integer, intent(in) :: iky, iz_start, m_start
         integer :: iz, m, dtube

         iz = iz_start; m = m_start; dtube = 0
         do
            if (visited(iz, m)) then
               if (periodic(iky) .and. iz == -nzgrid .and. m == m_start .and. dtube > 0) exit
               call mp_abort('Non-twisting flux tube: a point lies on two chains. Aborting.')
            end if
            visited(iz, m) = .true.
            m = m - dm(iz)
            if (abs(m) > mmax) exit
            iz = iz + 1
            if (iz == nzgrid) then
               iz = -nzgrid; dtube = dtube + 1
            end if
         end do

      end subroutine mark_chain

      ! Walk along the chain (ie, iky) from its first point; count its points (np) and copies (nc),
      ! and store them if <store>
      subroutine walk_chain(iky, ie, iz_start, m_start)

         integer, intent(in) :: iky, ie, iz_start, m_start
         integer :: iz, m, dtube

         np = 0; nc = 0; dtube = 0
         iz = iz_start; m = m_start
         call add_point(iky, ie, iz, m, dtube)
         ! An open chain that starts at the end of the tube: its copy at z(nzgrid) of the previous tube
         if (iz == -nzgrid .and. .not. periodic(iky)) call add_copy(iky, ie, nzgrid, m, -1)
         do
            m = m - dm(iz)
            if (abs(m) > mmax) exit
            iz = iz + 1
            if (iz == nzgrid) then
               ! Passing the end of the tube: store the point at z(nzgrid), as for a periodic chain
               call add_point(iky, ie, nzgrid, m, dtube)
               if (periodic(iky)) exit
               ! and its copy at z(-nzgrid) of the next tube, from where the chain continues
               call add_copy(iky, ie, -nzgrid, m, dtube + 1)
               iz = -nzgrid; dtube = dtube + 1
            else
               call add_point(iky, ie, iz, m, dtube)
            end if
         end do

      end subroutine walk_chain

      subroutine add_point(iky, ie, iz, m, dtube)
         integer, intent(in) :: iky, ie, iz, m, dtube
         np = np + 1
         if (.not. store) return
         point_iz(point_offset(ie, iky) + np) = iz
         point_ikx(point_offset(ie, iky) + np) = ikx_of(m)
         point_dtube(point_offset(ie, iky) + np) = dtube
      end subroutine add_point

      ! The copy at z(nzgrid) of the previous tube is a copy of the first point; otherwise of the last point
      subroutine add_copy(iky, ie, iz, m, dtube)
         integer, intent(in) :: iky, ie, iz, m, dtube
         nc = nc + 1
         if (.not. store) return
         copy_iz(dup_offset(ie, iky) + nc) = iz
         copy_ikx(dup_offset(ie, iky) + nc) = ikx_of(m)
         copy_dtube(dup_offset(ie, iky) + nc) = dtube
         if (iz == nzgrid) then
            copy_src(dup_offset(ie, iky) + nc) = 1
         else
            copy_src(dup_offset(ie, iky) + nc) = np
         end if
      end subroutine add_copy

   end subroutine init_extended_zgrid_non_twisting

   !****************************************************************************
   !            Describe the twist-and-shift chains point by point
   !****************************************************************************
   ! The chain (ie, iky) consists of <nsegments> 2pi segments: the first one covers
   ! iz_low(1):iz_up(1) and every next segment iz_low(iseg)+1:iz_up(iseg), since its
   ! first point is shared with the last point of the previous segment. That shared
   ! point is stored on the grid twice, as the last point of segment iseg-1 and as
   ! the first point of segment iseg, which is listed as a copy (dup). The phase
   ! factors are accumulated in the same way as in the segment loops, so that the
   ! mapping to and from the extended z-grid is unchanged.
   !****************************************************************************
   subroutine build_chain_points_from_segments

      use grids_kxky, only: naky

      implicit none

      integer :: iky, ie, iseg, iz, ip, ic
      complex :: shift_to, shift_from

      !-------------------------------------------------------------------------

      ! Count the points and copies of every chain
      call allocate_chain_offsets
      do iky = 1, naky
         do ie = 1, neigen(iky)
            nzext(ie, iky) = nsegments(ie, iky) * nzed_segment + 1
            ndup(ie, iky) = nsegments(ie, iky) - 1
         end do
      end do
      call allocate_chain_points

      do iky = 1, naky
         do ie = 1, neigen(iky)
            ip = point_offset(ie, iky); ic = dup_offset(ie, iky)
            shift_to = 1.; shift_from = 1.
            do iseg = 1, nsegments(ie, iky)
               if (iseg > 1) then
                  shift_to = shift_to / phase_shift(iky)
                  shift_from = shift_from * phase_shift(iky)
                  ! The first point of this segment is a copy of the last point of the previous one
                  ic = ic + 1
                  copy_iz(ic) = iz_low(iseg)
                  copy_ikx(ic) = ikxmod(iseg, ie, iky)
                  copy_dtube(ic) = iseg - 1
                  copy_src(ic) = ip - point_offset(ie, iky)
                  copy_phase_from(ic) = shift_from
               end if
               do iz = iz_low(iseg) + merge(0, 1, iseg == 1), iz_up(iseg)
                  ip = ip + 1
                  point_iz(ip) = iz
                  point_ikx(ip) = ikxmod(iseg, ie, iky)
                  point_dtube(ip) = iseg - 1
                  point_phase_to(ip) = shift_to
                  point_phase_from(ip) = shift_from
               end do
            end do
         end do
      end do

   end subroutine build_chain_points_from_segments

   !****************************************************************************
   !        Allocate the number of points per chain (filled by the caller)
   !****************************************************************************
   subroutine allocate_chain_offsets

      use grids_kxky, only: naky

      implicit none

      call deallocate_chain_points
      allocate (nzext(max(neigen_max, 1), naky)); nzext = 0
      allocate (ndup(max(neigen_max, 1), naky)); ndup = 0

   end subroutine allocate_chain_offsets

   !****************************************************************************
   !  Allocate the 1D arrays of points and copies, once <nzext> and <ndup> are known
   !****************************************************************************
   subroutine allocate_chain_points

      use grids_kxky, only: naky

      implicit none

      integer :: iky, ie, np, nc

      !-------------------------------------------------------------------------

      allocate (point_offset(max(neigen_max, 1), naky)); point_offset = 0
      allocate (dup_offset(max(neigen_max, 1), naky)); dup_offset = 0
      np = 0; nc = 0
      do iky = 1, naky
         do ie = 1, size(nzext, 1)
            point_offset(ie, iky) = np; dup_offset(ie, iky) = nc
            np = np + nzext(ie, iky); nc = nc + ndup(ie, iky)
         end do
      end do
      nzext_max = max(maxval(nzext), 1)
      allocate (point_iz(max(np, 1))); point_iz = 0
      allocate (point_ikx(max(np, 1))); point_ikx = 1
      allocate (point_dtube(max(np, 1))); point_dtube = 0
      allocate (point_phase_to(max(np, 1))); point_phase_to = 1.
      allocate (point_phase_from(max(np, 1))); point_phase_from = 1.
      allocate (copy_iz(max(nc, 1))); copy_iz = 0
      allocate (copy_ikx(max(nc, 1))); copy_ikx = 1
      allocate (copy_dtube(max(nc, 1))); copy_dtube = 0
      allocate (copy_src(max(nc, 1))); copy_src = 1
      allocate (copy_phase_from(max(nc, 1))); copy_phase_from = 1.

   end subroutine allocate_chain_points

   subroutine deallocate_chain_points

      implicit none

      if (allocated(nzext)) deallocate (nzext, ndup)
      if (allocated(point_offset)) deallocate (point_offset, dup_offset)
      if (allocated(point_iz)) deallocate (point_iz, point_ikx, point_dtube, point_phase_to, point_phase_from)
      if (allocated(copy_iz)) deallocate (copy_iz, copy_ikx, copy_dtube, copy_src, copy_phase_from)

   end subroutine deallocate_chain_points

   !****************************************************************************
   !              Location of a point of a chain, or of a copy
   !****************************************************************************
   pure integer function izext_iz(izext, ie, iky)
      integer, intent(in) :: izext, ie, iky
      izext_iz = point_iz(point_offset(ie, iky) + izext)
   end function izext_iz

   pure integer function izext_ikx(izext, ie, iky)
      integer, intent(in) :: izext, ie, iky
      izext_ikx = point_ikx(point_offset(ie, iky) + izext)
   end function izext_ikx

   pure integer function izext_dtube(izext, ie, iky)
      integer, intent(in) :: izext, ie, iky
      izext_dtube = point_dtube(point_offset(ie, iky) + izext)
   end function izext_dtube

   pure integer function dup_iz(idup, ie, iky)
      integer, intent(in) :: idup, ie, iky
      dup_iz = copy_iz(dup_offset(ie, iky) + idup)
   end function dup_iz

   pure integer function dup_ikx(idup, ie, iky)
      integer, intent(in) :: idup, ie, iky
      dup_ikx = copy_ikx(dup_offset(ie, iky) + idup)
   end function dup_ikx

   pure integer function dup_dtube(idup, ie, iky)
      integer, intent(in) :: idup, ie, iky
      dup_dtube = copy_dtube(dup_offset(ie, iky) + idup)
   end function dup_dtube

   pure integer function dup_src(idup, ie, iky)
      integer, intent(in) :: idup, ie, iky
      dup_src = copy_src(dup_offset(ie, iky) + idup)
   end function dup_src

   !****************************************************************************
   !            Flux tube index <dtube> tubes to the right of tube <it>
   !****************************************************************************
   pure integer function tube_index(it, dtube)

      use grids_z, only: ntubes

      implicit none

      integer, intent(in) :: it, dtube

      tube_index = modulo(it - 1 + dtube, ntubes) + 1

   end function tube_index

!###############################################################################
!################################## FUNCTIONS ##################################
!###############################################################################

   !****************************************************************************
   !                                Fill ghost zones
   !****************************************************************************
   subroutine fill_zed_ghost_zones(it, iseg, ie, iky, g, gleft, gright)

      use grids_z, only: nzgrid

      implicit none

      integer, intent(in) :: it, iseg, ie, iky
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:), intent(out) :: gleft, gright

      integer :: nseg

      ! stream_sign > 0 --> stream speed < 0

      nseg = nsegments(ie, iky)

      if (iseg == 1) then
         if (periodic(iky)) then
            gleft = phase_shift(iky) * g(iky, ikxmod(iseg, ie, iky), iz_up(nseg) - 2:iz_up(nseg) - 1, it)
         else
            gleft = 0.0
         end if
      else
         gleft = phase_shift(iky) * g(iky, ikxmod(iseg - 1, ie, iky), iz_up(iseg - 1) - 2:iz_up(iseg - 1) - 1, it_left(it))
      end if

      if (nseg > iseg) then
         ! connect to segment with larger theta-theta0 (on right)
         gright = g(iky, ikxmod(iseg + 1, ie, iky), iz_low(iseg + 1) + 1:iz_low(iseg + 1) + 2, it_right(it)) / phase_shift(iky)
      else
         ! apply periodic BC where necessary and zero BC otherwise
         if (periodic(iky)) then
            gright = g(iky, ikxmod(iseg, ie, iky), iz_low(1) + 1:iz_low(1) + 2, it) / phase_shift(iky)
         else
            gright = 0.0
         end if
      end if

   end subroutine fill_zed_ghost_zones

   !****************************************************************************
   !                                Fill ghost zones
   !****************************************************************************
   subroutine fill_zext_ghost_zones(iky, pdf_ext, pdf_left, pdf_right)

      implicit none

      integer, intent(in) :: iky
      complex, dimension(:), intent(in) :: pdf_ext
      complex, intent(out) :: pdf_left, pdf_right

      integer :: nz_ext

      ! n_zext is the number of grid points in this extended zed domain
      nz_ext = size(pdf_ext)

      ! if periodic BCs are applied in zed, then ghost zones at ends of extended domain
      ! should be filled using periodicity (with any appropriate phase shift)
      ! otherwise, zero incoming BC is used
      if (periodic(iky)) then
         pdf_left = pdf_ext(nz_ext - 1) * phase_shift(iky)
         pdf_right = pdf_ext(2) / phase_shift(iky)
      else
         pdf_left = 0.0
         pdf_right = 0.0
      end if

   end subroutine fill_zext_ghost_zones

   !****************************************************************************
   !                           Map to extended z-grid
   !****************************************************************************
   subroutine map_to_extended_zgrid(it, ie, iky, g, gext, ulim)

      use grids_z, only: nzgrid

      implicit none

      integer, intent(in) :: it, ie, iky
      complex, dimension(:, -nzgrid:, :), intent(in) :: g
      complex, dimension(:), intent(out) :: gext
      integer, intent(out) :: ulim

      integer :: izext, ip

      ulim = nzext(ie, iky)
      do izext = 1, ulim
         ip = point_offset(ie, iky) + izext
         gext(izext) = g(point_ikx(ip), point_iz(ip), tube_index(it, point_dtube(ip))) * point_phase_to(ip)
      end do

   end subroutine map_to_extended_zgrid

   !****************************************************************************
   !                           Map from extended z-grid
   !****************************************************************************
   subroutine map_from_extended_zgrid(it, ie, iky, gext, g)

      use grids_z, only: nzgrid

      implicit none

      integer, intent(in) :: it, ie, iky
      complex, dimension(:), intent(in) :: gext
      complex, dimension(:, -nzgrid:, :), intent(in out) :: g

      integer :: izext, idup, ip, ic

      ! The copies (e.g. the first point of every 2pi segment but the first) and the points of the chain
      ! are different grid locations, so the order in which they are written does not matter
      do idup = 1, ndup(ie, iky)
         ic = dup_offset(ie, iky) + idup
         g(copy_ikx(ic), copy_iz(ic), tube_index(it, copy_dtube(ic))) = gext(copy_src(ic)) * copy_phase_from(ic)
      end do
      do izext = 1, nzext(ie, iky)
         ip = point_offset(ie, iky) + izext
         g(point_ikx(ip), point_iz(ip), tube_index(it, point_dtube(ip))) = gext(izext) * point_phase_from(ip)
      end do

   end subroutine map_from_extended_zgrid

   !****************************************************************************
   !             A chain on the extended z-grid, with ghost points
   !****************************************************************************
   ! Map the chain (ie, iky) of g(ky, kx, z, tube) to the extended z-grid, and fill two
   ! ghost points beyond each end: by periodicity for periodic chains (whose last point
   ! is the same as the first), and zero otherwise. Used to take z-derivatives along
   ! the chains of the non-twisting flux tube with the same finite-difference routines.
   !****************************************************************************
   subroutine get_chain_with_ghosts(it, ie, iky, g, gext, gleft, gright)

      use grids_z, only: nzgrid

      implicit none

      integer, intent(in) :: it, ie, iky
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:), intent(out) :: gext
      complex, dimension(2), intent(out) :: gleft, gright

      integer :: n

      call map_to_extended_zgrid(it, ie, iky, g(iky, :, :, :), gext, n)
      if (periodic(iky) .and. n > 3) then
         gleft = gext(n - 2:n - 1) * phase_shift(iky)
         gright = gext(2:3) / phase_shift(iky)
      else
         gleft = 0.0
         gright = 0.0
      end if

   end subroutine get_chain_with_ghosts

   !****************************************************************************
   !          Apply a finite-difference operator in z along every chain
   !****************************************************************************
   ! Used by the non-twisting flux tube, whose chains change kx at any z point, instead
   ! of the loops over 2pi segments. Every chain is mapped to the extended z-grid with
   ! two ghost points at each end (see <get_chain_with_ghosts>), the operator is applied
   ! as on a single segment (iseg = nseg = 1), and the result is mapped back, including
   ! the copies of the points. Chains shorter than five points (possible for large ky)
   ! are padded with zeros, which is the zero boundary condition beyond their ends, since
   ! the one-sided end formulas need more points.
   !     <scheme>   zop_upwind3: third order upwind dg/dz; zop_centered2: second order
   !                centered dg/dz; zop_fd_cell_centers: dg/dz at cell centres;
   !                zop_cell_centers: g at cell centres; zop_fourth: d^4 g/dz^4
   !     <sgn>      sign of the streaming speed (upwinding), <upwnd> = zed_upwind
   !****************************************************************************
   subroutine zed_operator_along_chains(g, dg, scheme, sgn, upwnd)

      use grids_z, only: nzgrid, ntubes, delzed
      use grids_kxky, only: naky
      use calculations_finite_differences, only: third_order_upwind_zed, second_order_centered_zed
      use calculations_finite_differences, only: fd_cell_centers_zed, cell_centers_zed
      use calculations_finite_differences, only: fourth_derivative_second_centered_zed

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: dg
      integer, intent(in) :: scheme, sgn
      real, intent(in) :: upwnd

      integer :: iky, ie, it, n, np, pad
      logical :: per
      complex, dimension(2) :: gl, gr
      complex, dimension(:), allocatable :: gext, dext, gpad, dpad

      allocate (gext(nzext_max), dext(nzext_max), gpad(nzext_max + 4), dpad(nzext_max + 4))
      do iky = 1, naky
         per = periodic(iky)
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               n = nzext(ie, iky)
               call get_chain_with_ghosts(it, ie, iky, g, gext(:n), gl, gr)
               pad = 0
               if (n < 5 .and. .not. per) pad = 2
               np = n + 2 * pad
               gpad(:np) = 0.0
               gpad(pad + 1:pad + n) = gext(:n)
               select case (scheme)
               case (zop_upwind3)
                  call third_order_upwind_zed(1, 1, 1, gpad(:np), delzed(0), sgn, gl, gr, per, dpad(:np))
               case (zop_centered2)
                  call second_order_centered_zed(1, 1, 1, gpad(:np), delzed(0), sgn, gl, gr, per, dpad(:np))
               case (zop_fd_cell_centers)
                  call fd_cell_centers_zed(1, gpad(:np), delzed(0), sgn, gl(2), gr(1), dpad(:np))
               case (zop_cell_centers)
                  call cell_centers_zed(1, gpad(:np), upwnd, sgn, gl(2), gr(1), dpad(:np))
               case (zop_fourth)
                  call fourth_derivative_second_centered_zed(1, 1, 1, gpad(:np), delzed(0), gl, gr, per, dpad(:np))
               end select
               dext(:n) = dpad(pad + 1:pad + n)
               call map_from_extended_zgrid(it, ie, iky, dext(:n), dg(iky, :, :, :))
            end do
         end do
      end do
      deallocate (gext, dext, gpad, dpad)

   end subroutine zed_operator_along_chains

   !****************************************************************************
   !        Make the copies of the points of every chain equal to the points
   !****************************************************************************
   ! The points of a chain that are stored twice on the (kx, z) grid (the shared end points
   ! of the 2pi segments, or of the non-twisting flux tube) and the last point of a
   ! periodic chain are set from the point of the chain that they copy.
   !****************************************************************************
   subroutine enforce_chain_copies(fld)

      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: fld

      integer :: iky, ie, it, n
      complex, dimension(:), allocatable :: fext

      allocate (fext(nzext_max))
      do iky = 1, naky
         do ie = 1, neigen(iky)
            n = nzext(ie, iky)
            do it = 1, ntubes
               call map_to_extended_zgrid(it, ie, iky, fld(iky, :, :, :), fext, n)
               if (periodic(iky)) fext(n) = fext(1) / phase_shift(iky)
               call map_from_extended_zgrid(it, ie, iky, fext(:n), fld(iky, :, :, :))
            end do
         end do
      end do
      deallocate (fext)

   end subroutine enforce_chain_copies

   !****************************************************************************
   !                                      Map 
   !****************************************************************************
   subroutine map_to_iz_ikx_from_izext(iky, ie, iz_from_izext, ikx_from_izext)

      implicit none

      integer, intent(in) :: iky, ie
      integer, dimension(:), intent(out) :: iz_from_izext, ikx_from_izext

      integer :: n

      n = nzext(ie, iky)
      iz_from_izext(:n) = point_iz(point_offset(ie, iky) + 1:point_offset(ie, iky) + n)
      ikx_from_izext(:n) = point_ikx(point_offset(ie, iky) + 1:point_offset(ie, iky) + n)

   end subroutine map_to_iz_ikx_from_izext

!###############################################################################
!########################### FINISH EXTENDED Z-GRID ############################
!###############################################################################

   subroutine finish_extended_zgrid

      implicit none

      if (allocated(neigen)) deallocate (neigen)
      if (allocated(periodic)) deallocate (periodic)
      if (allocated(nsegments)) deallocate (nsegments)
      if (allocated(iz_low)) deallocate (iz_low, iz_mid, iz_up)
      if (allocated(ikxmod)) deallocate (ikxmod)
      if (allocated(it_right)) deallocate (it_right)
      if (allocated(it_left)) deallocate (it_left)
      if (allocated(phase_shift)) deallocate (phase_shift)
      call deallocate_chain_points

      initialised_extended_zgrid = .false.

   end subroutine finish_extended_zgrid

end module grids_extended_zgrid
