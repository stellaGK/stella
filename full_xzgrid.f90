module full_xzgrid

   implicit none

   public :: init_full_xzgrid, finish_full_xzgrid
   public :: map_to_full_xzgrid
   public :: map_from_full_xzgrid
   public :: xz_idx

   interface map_to_full_xzgrid
      module procedure map_to_full_xzgrid
      module procedure map_from_ezgrid_to_full_xzgrid
   end interface map_to_full_xzgrid

   interface map_from_full_xzgrid
      module procedure map_from_full_xzgrid
      module procedure map_to_ezgrid_from_full_xzgrid
   end interface map_from_full_xzgrid

   integer, dimension(:), allocatable :: nelements

   logical :: full_xzgrid_initialized = .false.

contains

   subroutine init_full_xzgrid

      use kt_grids, only: nakx, naky
      use zgrid, only: nztot, ntubes
      use extended_zgrid, only: periodic

      implicit none

      integer :: iky, pm

      if (full_xzgrid_initialized) return
      full_xzgrid_initialized = .true.

      if (.not. allocated(nelements)) allocate (nelements(naky))

      do iky = 1, naky
         pm = 0
         if (periodic(iky)) pm = 1
         nelements(iky) = nakx * ((nztot - 1) * ntubes + 1 - pm)
      end do

   end subroutine init_full_xzgrid

   subroutine map_to_full_xzgrid(iky, g, g_full)

      use kt_grids, only: nakx
      use zgrid, only: nzgrid, ntubes
      use extended_zgrid, only: periodic

      implicit none

      integer, intent(in) :: iky
      complex, dimension(:, -nzgrid:, :), intent(in) :: g
      complex, dimension(:), intent(out) :: g_full

      integer :: ikx, iz, it, pm

      pm = 0
      if (periodic(iky)) pm = 1

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid - pm
            do ikx = 1, nakx
               g_full(xz_idx(ikx, iz, it)) = g(ikx, iz, it)
            end do
         end do
      enddo

   end subroutine map_to_full_xzgrid

   subroutine map_from_ezgrid_to_full_xzgrid(it, ie, iky, g_ext, g_full)

      use extended_zgrid, only: ikxmod, nzed_segment, nsegments
      use extended_zgrid, only: iz_low, it_right

      implicit none

      integer, intent(in) :: it, ie, iky
      complex, dimension(:), intent(in) :: g_ext
      complex, dimension(:), intent(out) :: g_full

      integer :: idx, iseg, ikx, itmod
      integer :: llim, ulim

      ! avoid double-counting at boundaries between 2pi segments
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      llim = 1; ulim = nzed_segment + 1
      do idx = llim, ulim
         g_full(xz_idx(ikx, iz_low(iseg) + idx - llim, it)) = g_ext(idx)
      end do
      if (nsegments(ie, iky) > 1) then
         itmod = it
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            itmod = it_right(itmod)
            llim = ulim + 1
            ulim = llim + nzed_segment - 1
            g_full(xz_idx(ikx, iz_low(iseg), itmod)) = g_ext(llim - 1)
            do idx = llim, ulim
               g_full(xz_idx(ikx, iz_low(iseg) + idx - llim + 1, itmod)) = g_ext(idx)
            end do
         end do
      end if

   end subroutine map_from_ezgrid_to_full_xzgrid

   subroutine map_from_full_xzgrid(iky, g_full, g)

      use kt_grids, only: nakx
      use zgrid, only: nzgrid, ntubes
      use extended_zgrid, only: periodic

      implicit none

      integer, intent(in) :: iky
      complex, dimension(:), intent(in) :: g_full
      complex, dimension(:, -nzgrid:, :), intent(in out) :: g

      integer :: iz, ikx, it, pm

      pm = 0
      if (periodic(iky)) pm = 1
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid - pm
            do ikx = 1, nakx
               g(ikx, iz, it) = g_full(xz_idx(ikx, iz, it))
            end do
         end do
      end do
      if (periodic(iky)) g(:, nzgrid, :) = g(:, -nzgrid, :)

   end subroutine map_from_full_xzgrid

   subroutine map_to_ezgrid_from_full_xzgrid(it, ie, iky, g_full, g_ext)

      use extended_zgrid, only: ikxmod, nzed_segment, nsegments
      use extended_zgrid, only: iz_low, it_right

      implicit none

      integer, intent(in) :: it, ie, iky
      complex, dimension(:), intent(out) :: g_ext
      complex, dimension(:), intent(in) :: g_full

      integer :: idx, iseg, ikx, itmod
      integer :: llim, ulim

      ! avoid double-counting at boundaries between 2pi segments
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      llim = 1; ulim = nzed_segment + 1
      do idx = llim, ulim
         g_ext(idx) = g_full(xz_idx(ikx, iz_low(iseg) + idx - llim, it))
      end do
      if (nsegments(ie, iky) > 1) then
         itmod = it
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            itmod = it_right(itmod)
            llim = ulim + 1
            ulim = llim + nzed_segment - 1
            do idx = llim, ulim
               g_ext(idx) = g_full(xz_idx(ikx, iz_low(iseg) + idx - llim + 1, itmod))
            end do
         end do
      end if

   end subroutine map_to_ezgrid_from_full_xzgrid

   subroutine finish_full_xzgrid

      implicit none

      if (allocated(nelements)) deallocate (nelements)

      full_xzgrid_initialized = .false.

   end subroutine finish_full_xzgrid

   elemental function xz_idx(ikx, iz, it)

      use kt_grids, only: nakx
      use zgrid, only: nzgrid, ntubes, nztot

      implicit none

      integer xz_idx
      integer, intent(in) :: ikx, iz, it

      xz_idx = ikx + nakx * (iz + nzgrid + (nztot - 1) * (it - 1))

   end function xz_idx

end module full_xzgrid
