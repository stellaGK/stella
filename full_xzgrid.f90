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
    module procedure map_to_full_xzgrid
    module procedure map_from_ezgrid_to_full_xzgrid
  end interface map_from_full_xzgrid

contains

  subroutine init_full_xzgrid

    use zgrid, only: boundary_option_switch
    use zgrid, only: boundary_option_self_periodic
    use zgrid, only: boundary_option_linked
    use zgrid, only: nperiod, nzgrid, nzed, ntubes
    use kt_grids, only: nakx, naky
    use kt_grids, only: jtwist, ikx_twist_shift
    use kt_grids, only: aky, ikx_max
    use species, only: nspec

    implicit none

  end subroutine init_extended_zgrid

  subroutine fill_zed_ghost_zones (it, iseg, ie, iky, g, gleft, gright)

    use zgrid, only: nzgrid
    use kt_grids, only: zonal_mode

    implicit none

    integer, intent (in) :: it, iseg, ie, iky
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: g
    complex, dimension (:), intent (out) :: gleft, gright

    ! stream_sign > 0 --> stream speed < 0

    if (iseg == 1) then
       ! if zonal mode, then periodic BC instead of zero BC
       if (zonal_mode(iky)) then
          gleft = g(iky,ikxmod(iseg,ie,iky),iz_up(iseg)-2:iz_up(iseg)-1,it)
       else
          gleft = 0.0
       end if
    else
       gleft = g(iky,ikxmod(iseg-1,ie,iky),iz_up(iseg-1)-2:iz_up(iseg-1)-1,it_left(it))
    end if
    
    if (nsegments(ie,iky) > iseg) then
       ! connect to segment with larger theta-theta0 (on right)
       gright = g(iky,ikxmod(iseg+1,ie,iky),iz_low(iseg+1)+1:iz_low(iseg+1)+2,it_right(it))
    else
       ! apply periodic BC to zonal mode and zero BC otherwise
       if (zonal_mode(iky)) then
          gright = g(iky,ikxmod(iseg,ie,iky),iz_low(iseg)+1:iz_low(iseg)+2,it)
       else
          gright = 0.0
       end if
    end if
    
  end subroutine fill_zed_ghost_zones

  subroutine map_to_xzgrid (it, ie, iky, g, g_full, ulim)

    use zgrid, only: nzgrid

    implicit none

    integer, intent (in) :: it, ie, iky
    complex, dimension (:,-nzgrid:,:), intent (in) :: g
    complex, dimension (:), intent (out) :: gext
    integer, intent (out) :: ulim

    integer :: iseg, ikx, itmod
    integer :: llim, zm = 0

    if(ky.eq.0) zm = 1

    ! avoid double-counting at boundaries between 2pi segments
    do iz = -nzgrid, nzgrid
      do ikx = 1+zm, nakx
          gext(llim:ulim) = g(ikx,iz_low(iseg)+1:iz_up(iseg),itmod)
       enddo
     enddo
    end if

  end subroutine map_to_xzgrid

  subroutine map_from_ezgrid_to_xzgrid (it, ie, iky, g_ext, g_full, ulim)

    use zgrid, only: nzgrid
    use extended_zgrid, only: ikxmod

    implicit none

    integer, intent (in) :: it, ie, iky
    complex, dimension (:), intent (out) :: gext
    complex, dimension (:,-nzgrid:,:), intent (in) :: gext
    integer, intent (out) :: ulim

    integer :: iseg, ikx, itmod
    integer :: llim

    ! avoid double-counting at boundaries between 2pi segments
    iseg = 1
    ikx = ikxmod(iseg,ie,iky)
    llim = 1 ; ulim = nzed_segment+1
    gext(llim:ulim) = g(ikx,iz_low(iseg):iz_up(iseg),it)
    if (nsegments(ie,iky) > 1) then
       itmod = it
       do iseg = 2, nsegments(ie,iky)
          ikx = ikxmod(iseg,ie,iky)
          itmod = it_right(itmod)
          llim = ulim+1
          ulim = llim+nzed_segment-1
          gext(llim:ulim) = g(ikx,iz_low(iseg)+1:iz_up(iseg),itmod)
       end do
    end if

  end subroutine map_from_ezgrid_to_xzgrid

  subroutine map_from_xzgrid (it, ie, iky, gext, g)

    use zgrid, only: nzgrid

    implicit none

    integer, intent (in) :: it, ie, iky
    complex, dimension (:), intent (in) :: gext
    complex, dimension (:,-nzgrid:,:), intent (in out) :: g

    integer :: iseg, ikx, itmod
    integer :: llim, ulim

    iseg = 1
    ikx = ikxmod(iseg,ie,iky)
    llim = 1 ; ulim = nzed_segment+1
    g(ikx,iz_low(iseg):iz_up(iseg),it) = gext(llim:ulim)
    if (nsegments(ie,iky) > 1) then
       itmod = it
       do iseg = 2, nsegments(ie,iky)
          llim = ulim+1
          ulim = llim+nzed_segment-1
          ikx = ikxmod(iseg,ie,iky)
          itmod = it_right(itmod)
          g(ikx,iz_low(iseg),itmod) = gext(llim-1)
          g(ikx,iz_low(iseg)+1:iz_up(iseg),itmod) = gext(llim:ulim)
       end do
    end if

  end subroutine map_from_xzgrid

  subroutine finish_extended_zgrid

    implicit none

    if (allocated(neigen)) deallocate (neigen)
    if (allocated(periodic)) deallocate (periodic)
    if (allocated(nsegments)) deallocate (nsegments)
    if (allocated(iz_low)) deallocate (iz_low, iz_mid, iz_up)
    if (allocated(ikxmod)) deallocate (ikxmod)
    if (allocated(it_right)) deallocate (it_right)
    if (allocated(it_left)) deallocate (it_left)

    extended_zgrid_initialized = .false.

  end subroutine finish_extended_zgrid

end module extended_zgrid
