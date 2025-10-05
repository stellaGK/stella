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

   ! Only initialise once
   logical :: initialised_extended_zgrid = .false.

contains

!###############################################################################
!######################### INITIALISE EXTENDED Z-GRID ##########################
!###############################################################################

   subroutine init_extended_zgrid

      use mp, only: broadcast
      use grids_z, only: boundary_option_switch
      use grids_z, only: boundary_option_self_periodic
      use grids_z, only: boundary_option_linked
      use grids_z, only: boundary_option_linked_stellarator
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

      if (debug) write (*, *) 'extended_zgrid::boundary_option_switch'  
      if (boundary_option_switch == boundary_option_linked .or. boundary_option_switch == boundary_option_linked_stellarator) then

         ! All periodic modes (e.g., the zonal mode) have no connections
         do iky = 1, naky
            if (periodic(iky)) then
               neigen(iky) = nakx
            else
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

      ! this is the number of unique zed values in all segments but the first
      ! the first has one extra unique zed value (all others have one grid common
      ! with the previous segment due to periodicity)
      if (neigen_max>0) then 
         nzed_segment = iz_up(1) - iz_low(1)
      else 
         nzed_segment = 0
      end if 

      call broadcast (phase_shift)
      call broadcast (periodic)
      
   end subroutine init_extended_zgrid

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

      integer :: iseg, ikx, itmod
      integer :: llim
      complex :: curr_shift

      complex :: tmp

      ! avoid double-counting at boundaries between 2pi segments
      iseg = 1
      curr_shift = 1.
      ikx = ikxmod(iseg, ie, iky)
      llim = 1; ulim = nzed_segment + 1

      if (periodic(iky)) then
         gext(llim:ulim - 1) = g(ikx, -nzgrid:nzgrid - 1, it)
         gext(ulim) = g(ikx, -nzgrid, it) / phase_shift(iky)
      else
         gext(llim:ulim) = g(ikx, iz_low(iseg):iz_up(iseg), it) * curr_shift
         if (nsegments(ie, iky) > 1) then
            itmod = it
            do iseg = 2, nsegments(ie, iky)
               curr_shift = curr_shift / phase_shift(iky)
               ikx = ikxmod(iseg, ie, iky)
               itmod = it_right(itmod)

               tmp = gext(ulim)

               llim = ulim
               ulim = llim + nzed_segment
               gext(llim:ulim) = g(ikx, iz_low(iseg):iz_up(iseg), itmod) * curr_shift

            end do
         end if
      end if 

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

      integer :: iseg, ikx, itmod
      integer :: llim, ulim
      complex :: curr_shift

      complex :: tmp 

      iseg = 1
      curr_shift = 1.
      ikx = ikxmod(iseg, ie, iky)
      llim = 1; ulim = nzed_segment + 1

      if (periodic(iky)) then
         g(ikx, iz_low(iseg):iz_up(iseg) - 1, it) = gext(llim:ulim - 1)
         g(ikx, iz_up(iseg), it) = g(ikx, iz_low(iseg), it) / phase_shift(iky)
      else
         g(ikxmod(1, ie, iky), iz_low(1), :) = 0.0
         g(ikxmod(nsegments(ie, iky), ie, iky), iz_up(nsegments(ie, iky)), : ) = 0.0

         g(ikx, iz_low(iseg):iz_up(iseg), it) = gext(llim:ulim)
         if (nsegments(ie, iky) > 1) then
            itmod = it
            do iseg = 2, nsegments(ie, iky)
               curr_shift = curr_shift * phase_shift(iky)
               ikx = ikxmod(iseg, ie, iky)
               itmod = it_right(itmod)

               llim = ulim
               ulim = llim + nzed_segment
               g(ikx, iz_low(iseg):iz_up(iseg), itmod) = gext(llim:ulim) * curr_shift
               g(ikxmod(iseg - 1, ie, iky), iz_up(iseg - 1), itmod) = g(ikx, iz_low(iseg), itmod) * 0.5 * (phase_shift(iky) + 1)
!               g(ikxmod(iseg - 1, ie, iky), iz_up(iseg - 1), itmod) = 0.5 * ( g(ikxmod(iseg - 1, ie, iky), iz_up(iseg - 1), itmod) + g(ikx, iz_low(iseg), itmod) * phase_shift(iky))
               
            end do
         end if
      end if

   end subroutine map_from_extended_zgrid

   !****************************************************************************
   !                                      Map 
   !****************************************************************************
   subroutine map_to_iz_ikx_from_izext(iky, ie, iz_from_izext, ikx_from_izext)

      implicit none

      integer, intent(in) :: iky, ie
      integer, dimension(:), intent(out) :: iz_from_izext, ikx_from_izext

      integer :: iseg
      integer :: llim, ulim
      integer :: izext

      iseg = 1
      llim = 1; ulim = nzed_segment + 1
      ikx_from_izext(llim:ulim) = ikxmod(iseg, ie, iky)
      do izext = llim, ulim
         iz_from_izext(izext) = izext - llim + iz_low(iseg)
      end do
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            llim = ulim + 1
            ulim = llim + nzed_segment - 1
            ikx_from_izext(llim:ulim) = ikxmod(iseg, ie, iky)
            do izext = llim, ulim
               iz_from_izext(izext) = izext - llim + iz_low(iseg) + 1
            end do
         end do
      end if

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

      initialised_extended_zgrid = .false.

   end subroutine finish_extended_zgrid

end module grids_extended_zgrid
