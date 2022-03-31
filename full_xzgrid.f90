module full_xzgrid

   implicit none

   public :: init_full_xzgrid, finish_full_xzgrid
   public :: map_to_full_xzgrid
   public :: map_from_full_xzgrid

   interface map_to_full_xzgrid
      module procedure map_to_full_xzgrid
      module procedure map_from_ezgrid_to_full_xzgrid
   end interface map_to_full_xzgrid


   interface map_from_full_xzgrid
      module procedure map_from_full_xzgrid
      module procedure map_to_ezgrid_from_full_xzgrid
   end interface map_from_full_xzgrid

   integer, dimension (:), allocatable :: nelements

   logical :: full_xzgrid_initialized = .false.

contains

   subroutine init_full_xzgrid

      use zgrid, only: nztot
      use kt_grids, only: nakx, naky

      implicit none

      if (full_xzgrid_initialized) return
      full_xzgrid_initialized = .true.

      if (.not. allocated(nelements)) allocate (nelements(naky))

      nelements = nztot*nakx

   end subroutine init_full_xzgrid

   subroutine map_to_full_xzgrid (it, iky, g, g_full)

      use zgrid, only: nzgrid
      use kt_grids, only: nakx

      implicit none

      integer, intent (in) :: it, iky
      complex, dimension (:,-nzgrid:,:), intent (in) :: g
      complex, dimension (:), intent (out) :: g_full

      integer :: ikx, iz

      !integer :: zm = 0
      !if(aky(iky).lt.epsilon(0.0)) zm = 1

      do iz = -nzgrid, nzgrid
         do ikx = 1, nakx
            g_full(xz_idx(ikx, iz)) = g(ikx, iz, it)
         enddo
      enddo

   end subroutine map_to_full_xzgrid

   subroutine map_from_ezgrid_to_full_xzgrid (it, ie, iky, g_ext, g_full)

      use zgrid, only: nzgrid
      use extended_zgrid, only: ikxmod, nzed_segment, nsegments
      use extended_zgrid, only: iz_low, it_right

      implicit none

      integer, intent (in) :: it, ie, iky
      complex, dimension (:), intent (in) :: g_ext
      complex, dimension (:), intent (out) :: g_full

      integer :: idx, iseg, ikx, itmod
      integer :: llim, ulim

      ! avoid double-counting at boundaries between 2pi segments
      iseg = 1
      ikx = ikxmod(iseg,ie,iky)
      llim = 1 ; ulim = nzed_segment+1
      do idx = llim, ulim
         g_full(xz_idx(ikx, iz_low(iseg) + idx - llim)) = g_ext(idx)
      enddo
      if (nsegments(ie,iky) > 1) then
         itmod = it
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg,ie,iky)
            itmod = it_right(itmod)
            llim = ulim+1
            ulim = llim+nzed_segment-1
            g_full(xz_idx(ikx, iz_low(iseg))) = g_ext(llim - 1)
            do idx = llim, ulim
               g_full(xz_idx(ikx, iz_low(iseg) + idx - llim + 1)) = g_ext(idx)
            enddo
         end do
      end if

   end subroutine map_from_ezgrid_to_full_xzgrid

   subroutine map_from_full_xzgrid (it, iky, g_full, g)

      use zgrid, only: nzgrid
      use kt_grids, only: nakx

      implicit none

      integer, intent (in) :: it, iky
      complex, dimension (:), intent (in) :: g_full
      complex, dimension (:,-nzgrid:,:), intent (in out) :: g

      integer :: iz, ikx

      do iz = -nzgrid, nzgrid
         do ikx = 1, nakx
            g(ikx, iz, it) = g_full(xz_idx(ikx, iz))
         enddo
      enddo

  end subroutine map_from_full_xzgrid

  subroutine map_to_ezgrid_from_full_xzgrid (it, ie, iky, g_full, g_ext)
      use zgrid, only: nzgrid
      use extended_zgrid, only: ikxmod, nzed_segment, nsegments
      use extended_zgrid, only: iz_low, it_right

      implicit none

      integer, intent (in) :: it, ie, iky
      complex, dimension (:), intent (out) :: g_ext
      complex, dimension (:), intent (in) :: g_full

      integer :: idx, iseg, ikx, itmod
      integer :: llim, ulim

      ! avoid double-counting at boundaries between 2pi segments
      iseg = 1
      ikx = ikxmod(iseg,ie,iky)
      llim = 1 ; ulim = nzed_segment+1
      do idx = llim, ulim
         g_ext(idx) = g_full(xz_idx(ikx, iz_low(iseg) + idx - llim))
      enddo
      if (nsegments(ie,iky) > 1) then
         itmod = it
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg,ie,iky)
            itmod = it_right(itmod)
            llim = ulim+1
            ulim = llim+nzed_segment-1
            do idx = llim, ulim
               g_ext(idx) = g_full(xz_idx(ikx, iz_low(iseg) + idx - llim + 1))
            enddo
         end do
      end if

  end subroutine map_to_ezgrid_from_full_xzgrid

   subroutine finish_full_xzgrid

      implicit none

      if (allocated(nelements)) deallocate (nelements)

      full_xzgrid_initialized = .false.

   end subroutine finish_full_xzgrid

   elemental function xz_idx(ikx, iz)

      use zgrid, only: nzgrid
      use kt_grids, only: nakx

      implicit none

      integer xz_idx
      integer, intent(in) :: ikx, iz

      xz_idx = ikx + (iz + nzgrid) * nakx

   end function xz_idx

end module full_xzgrid
