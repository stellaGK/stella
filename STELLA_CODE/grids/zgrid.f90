module zgrid

   implicit none

   public :: init_zgrid, finish_zgrid
   public :: nzgrid
   public :: nztot, nz2pi
   public :: zed
   public :: delzed
   public :: get_total_arc_length
   public :: get_arc_length_grid
   public :: boundary_option_switch
   public :: boundary_option_zero
   public :: boundary_option_self_periodic
   public :: boundary_option_linked
   public :: boundary_option_linked_stellarator
   public :: nzed, nperiod, ntubes, shat_zero, dkx_over_dky
   public :: boundary_option, zed_equal_arc, grad_x_grad_y_zero

   private

   integer :: nzed, nzgrid, nztot, nz2pi
   integer :: nperiod, ntubes
   logical :: zed_equal_arc
   real :: shat_zero, grad_x_grad_y_zero, dkx_over_dky
   real, dimension(:), allocatable :: zed, delzed

   integer :: boundary_option_switch
   integer, parameter :: boundary_option_zero = 1, &
                         boundary_option_self_periodic = 2, &
                         boundary_option_linked = 3, &
                         boundary_option_linked_stellarator = 4

   logical :: initialised = .false.
   character(20) :: boundary_option

contains

   !============================================================================
   !============================= INITIALISE Z GRID ============================
   !============================================================================ 

   subroutine init_zgrid

      use mp, only: proc0
      use constants, only: pi

      use input_file_z_grid, only: read_namelist_z_grid, read_namelist_z_boundary_condition

      implicit none

      integer :: i

      if (initialised) return
      initialised = .true.

      if (proc0) then
         call read_namelist_z_grid(nzed, nperiod, ntubes, zed_equal_arc)
         call read_namelist_z_boundary_condition (boundary_option_switch, shat_zero, & 
                                                   grad_x_grad_y_zero, dkx_over_dky)
      end if
      
      ! <nzed> specifies the grid points left and right of z=0 in a single segment
      ! whereas <nzgrid> is the total number of grid points
      nzgrid = nzed / 2 + (nperiod - 1) * nzed
      call broadcast_parameters
      call compute_useful_quantities
      
   contains

      subroutine compute_useful_quantities

         implicit none
         

         if (.not. allocated(zed)) allocate (zed(-nzgrid:nzgrid))
         if (.not. allocated(delzed)) allocate (delzed(-nzgrid:nzgrid))

         zed = (/(i * pi / real(nzed / 2), i=-nzgrid, nzgrid)/)
         delzed(:nzgrid - 1) = zed(-nzgrid + 1:) - zed(:nzgrid - 1)
         delzed(nzgrid) = delzed(-nzgrid)

         nztot = 2 * nzgrid + 1
         ! number of zed in a 2*pi segment, including points at +/- pi
         nz2pi = 2 * (nzed / 2) + 1

      end subroutine compute_useful_quantities

      subroutine broadcast_parameters

         use mp, only: broadcast

         implicit none

         call broadcast(nzed)
         call broadcast(nzgrid)
         call broadcast(nperiod)
         call broadcast(ntubes)
         call broadcast(zed_equal_arc)
         call broadcast(shat_zero)
         call broadcast(boundary_option_switch)
         call broadcast(grad_x_grad_y_zero)
         call broadcast(dkx_over_dky)

      end subroutine broadcast_parameters

   end subroutine init_zgrid

   !============================================================================
   !============================= FINIALISE Z GRID ============================
   !============================================================================ 

   subroutine finish_zgrid

      implicit none

      if (allocated(zed)) deallocate (zed)
      if (allocated(delzed)) deallocate (delzed)

      initialised = .false.

   end subroutine finish_zgrid

   !============================================================================
   !======================== CALCULATE TOTAL ARC LENGTH ========================
   !============================================================================  
   subroutine get_total_arc_length(nz, gp, dz, length)

      implicit none

      integer, intent(in) :: nz
      real, dimension(-nz:), intent(in) :: gp
      real, intent(in) :: dz
      real, intent(out) :: length

      !----------------------------------------------------------------------

      ! <nz> is the number of positive or negative z-points, i.e., f(-nz:nz)
      ! <dz> is the step size along z  
      ! 1/<gp> is the function to integration along z
      call integrate_zed(nz, dz, 1./gp, length)

   end subroutine get_total_arc_length

   !============================================================================
   !======================== CALCULATE ARC LENGTH GRID =========================
   !============================================================================  
   subroutine get_arc_length_grid(nz_max, nzext_max, zboundary, gp, dz, zarc)

      implicit none

      integer, intent(in) :: nz_max, nzext_max
      real, intent(in) :: zboundary, dz
      real, dimension(-nzext_max:), intent(in) :: gp 
      real, dimension(-nzext_max:), intent(out) :: zarc

      integer :: iz

      !---------------------------------------------------------------------- 

      zarc(-nz_max) = zboundary
      if (nz_max /= nzext_max) then
         do iz = -nzext_max, -nz_max - 1
            call integrate_zed(nzext_max, dz, 1./gp(iz:-nz_max), zarc(iz))
            zarc(iz) = zarc(-nz_max) - zarc(iz)
         end do
      end if
      ! TODO: this seems very inefficient -- could just add incremental change at each zed,
      ! rather than recomputing from the boundary each time
      do iz = -nz_max + 1, nzext_max
         call integrate_zed(nz_max, dz, 1./gp(-nz_max:iz), zarc(iz))
         zarc(iz) = zarc(-nz_max) + zarc(iz)
      end do

   end subroutine get_arc_length_grid

   !============================================================================
   !============================ INTEGRATE ALONG Z =============================
   !============================================================================  
   ! Use the trapezoidal rule to integrate in zed
   subroutine integrate_zed(nz, dz, f, intf)

      implicit none

      integer, intent(in) :: nz
      real, intent(in) :: dz
      real, dimension(-nz:), intent(in) :: f
      real, intent(out) :: intf

      integer :: iz, iz_max

      !---------------------------------------------------------------------- 

      ! <nz> is the number of positive or negative z-points, i.e., f(-nz:nz)
      iz_max = -nz + size(f) - 1

      ! Initialize the intgration
      intf = 0.

      ! For each z-point add ( f(iz) + f(iz-1) ) * dz / 2
      do iz = -nz + 1, iz_max
         intf = intf + dz * (f(iz - 1) + f(iz))
      end do
      intf = 0.5 * intf

   end subroutine integrate_zed

end module zgrid
