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

      use kt_grids, only: nakx, naky, zonal_mode
      use zgrid, only: nztot, ntubes
      use extended_zgrid, only: periodic

      implicit none

      integer :: iky, pm, zm

      if (full_xzgrid_initialized) return
      full_xzgrid_initialized = .true.

      if (.not. allocated(nelements)) allocate (nelements(naky))

      do iky = 1, naky
         pm = 0; zm = 0
         if (periodic(iky)) pm = 1
         if (zonal_mode(iky)) zm = 1
         nelements(iky) = (nakx - zm) * ((nztot - 1) * ntubes + 1 - pm)
      end do

   end subroutine init_full_xzgrid

   subroutine map_to_full_xzgrid(iky, g, g_full)

      use kt_grids, only: nakx, zonal_mode
      use zgrid, only: nzgrid, ntubes
      use extended_zgrid, only: periodic

      implicit none

      integer, intent(in) :: iky
      complex, dimension(:, -nzgrid:, :), intent(in) :: g
      complex, dimension(:), intent(out) :: g_full

      integer :: ikx, iz, it, pm, zm

      pm = 0; zm = 0
      if (periodic(iky)) pm = 1
      if (zonal_mode(iky)) zm = 1

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid - pm
            do ikx = 1 + zm, nakx
               g_full(xz_idx(ikx, iz, it, zm)) = g(ikx, iz, it)
            end do
         end do
      end do

   end subroutine map_to_full_xzgrid

   subroutine map_from_ezgrid_to_full_xzgrid(it, ie, iky, g_ext, g_full)

      use kt_grids, only: zonal_mode
      use extended_zgrid, only: ikxmod, nzed_segment, nsegments
      use extended_zgrid, only: iz_low, it_right

      implicit none

      integer, intent(in) :: it, ie, iky
      complex, dimension(:), intent(in) :: g_ext
      complex, dimension(:), intent(out) :: g_full

      integer :: idx, iseg, ikx, itmod
      integer :: llim, ulim, zm

      ! avoid double-counting at boundaries between 2pi segments
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      if (zonal_mode(iky) .and. ikx == 1) return
      zm = 0
      if (zonal_mode(iky)) zm = 1
      llim = 1; ulim = nzed_segment + 1
      do idx = llim, ulim
         g_full(xz_idx(ikx, iz_low(iseg) + idx - llim, it, zm)) = g_ext(idx)
      end do
      if (nsegments(ie, iky) > 1) then
         itmod = it
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            itmod = it_right(itmod)
            llim = ulim + 1
            ulim = llim + nzed_segment - 1
            g_full(xz_idx(ikx, iz_low(iseg), itmod, zm)) = g_ext(llim - 1)
            do idx = llim, ulim
               g_full(xz_idx(ikx, iz_low(iseg) + idx - llim + 1, itmod, zm)) = g_ext(idx)
            end do
         end do
      end if

   end subroutine map_from_ezgrid_to_full_xzgrid

   subroutine map_from_full_xzgrid(iky, g_full, g)

      use kt_grids, only: nakx, zonal_mode
      use zgrid, only: nzgrid, ntubes
      use extended_zgrid, only: periodic

      implicit none

      integer, intent(in) :: iky
      complex, dimension(:), intent(in) :: g_full
      complex, dimension(:, -nzgrid:, :), intent(in out) :: g

      integer :: iz, ikx, it, pm, zm

      pm = 0; zm = 0
      if (periodic(iky)) pm = 1
      if (zonal_mode(iky)) zm = 1
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid - pm
            do ikx = 1 + zm, nakx
               g(ikx, iz, it) = g_full(xz_idx(ikx, iz, it, zm))
            end do
         end do
      end do
      if (periodic(iky)) g(:, nzgrid, :) = g(:, -nzgrid, :)
      if (zm == 1) g(zm, :, :) = 0.0

   end subroutine map_from_full_xzgrid

   subroutine map_to_ezgrid_from_full_xzgrid(it, ie, iky, g_full, g_ext)

      use kt_grids, only: zonal_mode
      use extended_zgrid, only: ikxmod, nzed_segment, nsegments
      use extended_zgrid, only: iz_low, it_right

      implicit none

      integer, intent(in) :: it, ie, iky
      complex, dimension(:), intent(out) :: g_ext
      complex, dimension(:), intent(in) :: g_full

      integer :: idx, iseg, ikx, itmod
      integer :: llim, ulim, zm

      ! avoid double-counting at boundaries between 2pi segments
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      if (zonal_mode(iky) .and. ikx == 1) return
      zm = 0
      if (zonal_mode(iky)) zm = 1
      llim = 1; ulim = nzed_segment + 1
      do idx = llim, ulim
         g_ext(idx) = g_full(xz_idx(ikx, iz_low(iseg) + idx - llim, it, zm))
      end do
      if (nsegments(ie, iky) > 1) then
         itmod = it
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            itmod = it_right(itmod)
            llim = ulim + 1
            ulim = llim + nzed_segment - 1
            do idx = llim, ulim
               g_ext(idx) = g_full(xz_idx(ikx, iz_low(iseg) + idx - llim + 1, itmod, zm))
            end do
         end do
      end if

   end subroutine map_to_ezgrid_from_full_xzgrid

   subroutine finish_full_xzgrid

      implicit none

      if (allocated(nelements)) deallocate (nelements)

      full_xzgrid_initialized = .false.

   end subroutine finish_full_xzgrid

   elemental function xz_idx(ikx, iz, it, zm)

      use kt_grids, only: nakx
      use zgrid, only: nzgrid, ntubes, nztot

      implicit none

      integer xz_idx
      integer, intent(in) :: ikx, iz, it, zm

      xz_idx = ikx - zm + (nakx - zm) * (iz + nzgrid + (nztot - 1) * (it - 1))

   end function xz_idx

end module full_xzgrid
