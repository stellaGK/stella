module extended_zgrid

  implicit none

  public :: nsegments
  public :: neigen
  public :: ikxmod
  public :: iz_low, iz_mid, iz_up
  public :: periodic
  public :: nzed_segment
  public :: fill_zed_ghost_zones
  public :: init_extended_zgrid, finish_extended_zgrid
  public :: map_to_extended_zgrid
  public :: map_from_extended_zgrid

  ! these arrays needed to keep track of connections between different
  ! 2pi segments
  integer :: nzed_segment
  integer, dimension (:), allocatable :: neigen
  integer, dimension (:), allocatable :: iz_low, iz_mid, iz_up
  integer, dimension (:,:), allocatable :: nsegments
  integer, dimension (:,:,:), allocatable :: ikxmod

  ! FLAG -- NEED TO IMPLEMENT PERIODIC FOR ZONAL FLOW
  logical, dimension (:), allocatable :: periodic

  logical :: extended_zgrid_initialized = .false.

contains

  subroutine init_extended_zgrid

    use zgrid, only: boundary_option_switch
    use zgrid, only: boundary_option_self_periodic
    use zgrid, only: boundary_option_linked
    use zgrid, only: nperiod, nzgrid, nzed
    use kt_grids, only: nakx, naky
    use kt_grids, only: jtwist_out, aky, ikx_max
    use species, only: nspec

    implicit none

    integer :: iseg, iky, ie, ntg, ikx
    integer :: nseg_max, neigen_max
    integer, dimension (:), allocatable :: ikx_shift_left
    integer, dimension (:,:), allocatable :: ikx_shift

    if (extended_zgrid_initialized) return
    extended_zgrid_initialized = .true.
    
    ntg = nzed/2

    if (.not. allocated(neigen)) allocate (neigen(naky))
    if (.not. allocated(periodic)) allocate (periodic(naky)) ; periodic = .false.

    if (boundary_option_switch==boundary_option_self_periodic) then
       periodic = .true.
    else
       where (abs(aky) < epsilon(0.0)) periodic = .true.
    end if

    select case (boundary_option_switch)
    case (boundary_option_linked)

       ! if linked BC, then iky=1 corresponds to ky=0 which has no connections
       neigen(1) = nakx
       if (naky > 1) then
          do iky = 2, naky
             ! must link different kx values at theta = +/- pi
             ! neigen is the number of independent eigenfunctions along the field line
             neigen(iky) = min((iky-1)*jtwist_out,nakx)
          end do
       end if

       neigen_max = maxval(neigen)

       if (.not. allocated(ikx_shift_left)) then
          allocate (ikx_shift_left(neigen_max)) ; ikx_shift_left = 0
          allocate (ikx_shift(nakx,naky)) ; ikx_shift = 0
       end if

       ! figure out how much to shift ikx by to get to
       ! the left-most (theta-theta0) in each set of connected 2pi segments
       ! note that theta0 goes from 0 to theta0_max and then from theta0_min back
       ! to -dtheta0
       do ikx = 1, neigen_max
          ! first ikx_max=nakx/2+1 theta0s are 0 and all positive theta0 values
          ! remainder are negative theta0s
          ! theta_0 = kx / ky / shat
          ! if ky > 0, then most positive theta_0 corresponds to most positive kx
          if (ikx <= ikx_max) then
             ikx_shift_left(ikx) = ikx_max-2*ikx+1
          else
             ikx_shift_left(ikx) = 3*ikx_max-2*ikx
          end if
       end do

       do iky = 1, naky
          ! ikx_shift is how much to shift each ikx by to connect
          ! to the next theta0 (from most positive to most negative)

          ! if ky > 0, then going to more negative theta0
          ! corresponds to going to more negative kx
          do ikx = 1, ikx_max
             ! if theta0 is sufficiently positive, shifting to more
             ! negative theta0 corresponds to decreasing ikx
             if (ikx-neigen(iky) > 0) then
                ikx_shift(ikx,iky) = -neigen(iky)
                ! if a positive theta0 connects to a negative theta0
                ! must do more complicated mapping of ikx
             else if (ikx-neigen(iky)+nakx >= ikx_max+1) then
                ikx_shift(ikx,iky) = nakx - neigen(iky)
             end if
          end do
          ! if theta0 is negative, then shifting to more negative
          ! theta0 corresponds to decreasing ikx
          do ikx = ikx_max+1, nakx
             ! if theta0 is sufficiently negative, it has no
             ! more negative theta0 with which it can connect
             if (ikx-neigen(iky) >= ikx_max) then
                ikx_shift(ikx,iky) = -neigen(iky)
             end if
             ! theta0 is positive
          end  do
       end do

       if (.not. allocated(nsegments)) allocate (nsegments(neigen_max,naky))

       do iky = 1, naky
          if (neigen(iky) == 0) then
             nsegments(:,iky) = 1
          else
             nsegments(:,iky) = (nakx-1)/neigen(iky)

             do ie = 1, mod(nakx-1,neigen(iky))+1
                nsegments(ie,iky) = nsegments(ie,iky) + 1
             end do
          end if
       end do

       nseg_max = maxval(nsegments)

       if (.not. allocated(iz_low)) then
          allocate (iz_low(nseg_max)) ; iz_low = -nzgrid
          allocate (iz_mid(nseg_max)) ; iz_mid = 0
          allocate (iz_up(nseg_max)) ; iz_up = nzgrid
       end if
       
    case default
       
       neigen = nakx ; neigen_max = nakx
       
       if (.not. allocated(ikx_shift_left)) then
          allocate (ikx_shift_left(neigen_max))
          allocate (ikx_shift(nakx,naky))
       end if
       ikx_shift = 0 ; ikx_shift_left = 0
       
       if (.not. allocated(nsegments)) then
          allocate (nsegments(neigen_max,naky))
       end if
       
       ! this is the number of 2pi poloidal segments in the extended theta domain,
       ! which is needed in initializing the reponse matrix and doing the implicit sweep
       nsegments = 2*(nperiod-1) + 1

       nseg_max = maxval(nsegments)
       
       if (.not. allocated(iz_low)) then
          allocate (iz_low(nseg_max))
          allocate (iz_mid(nseg_max))
          allocate (iz_up(nseg_max))
       end if

       ! iz_low(j) is the ig index corresponding to the inboard midplane from below (theta=-pi) within the jth segment
       ! iz_mid(j) is the ig index corresponding to the outboard midplane (theta=0) within the jth segment
       do iseg = 1, nseg_max
          iz_low(iseg) = -nzgrid + (iseg-1)*nzed
          iz_mid(iseg) = iz_low(iseg) + nzed/2
          iz_up(iseg) = iz_low(iseg) + nzed
       end do

    end select

    if (.not. allocated(ikxmod)) then
       allocate (ikxmod(nseg_max,neigen_max,naky))
       ! initialize ikxmod to nakx
       ! should not be necessary but just in case one tries to access
       ! a value beyond nsegments(ie,iky)
       ikxmod = nakx
    end if
    do iky = 1, naky
       ! only do the following once for each independent set of theta0s
       ! the assumption here is that all kx are on processor and sequential
       do ie = 1, neigen(iky)
          ! remap to start at theta0 = theta0_max
          ! (so that theta-theta0 is most negative)
          ! for this set of connected theta0s
          iseg = 1
          ikxmod(iseg,ie,iky) = ie + ikx_shift_left(ie)
          if (nsegments(ie,iky) > 1) then
             do iseg = 2, nsegments(ie,iky)
                ikxmod(iseg,ie,iky) = ikxmod(iseg-1,ie,iky) + ikx_shift(ikxmod(iseg-1,ie,iky),iky)
             end do
          end if
       end do
    end do

    if (allocated(ikx_shift_left)) deallocate (ikx_shift_left)
    if (allocated(ikx_shift)) deallocate (ikx_shift)

    ! this is the number of unique zed values in all segments but the first
    ! the first has one extra unique zed value (all others have one grid common
    ! with the previous segment due to periodicity)
    nzed_segment = iz_up(1)-iz_low(1)

  end subroutine init_extended_zgrid

  subroutine fill_zed_ghost_zones (iseg, ie, iky, g, gleft, gright)

    use zgrid, only: nzgrid
    use kt_grids, only: zonal_mode

    implicit none

    integer, intent (in) :: iseg, ie, iky
    complex, dimension (:,:,-nzgrid:), intent (in) :: g
    complex, dimension (:), intent (out) :: gleft, gright

    ! stream_sign > 0 --> stream speed < 0

    if (iseg == 1) then
       ! if zonal mode, then periodic BC instead of zero BC
       if (zonal_mode(iky)) then
          gleft = g(iky,ikxmod(iseg,ie,iky),iz_up(iseg)-2:iz_up(iseg)-1)
       else
          gleft = 0.0
       end if
    else
       gleft = g(iky,ikxmod(iseg-1,ie,iky),iz_up(iseg-1)-2:iz_up(iseg-1)-1)
    end if
    
    if (nsegments(ie,iky) > iseg) then
       ! connect to segment with larger theta-theta0 (on right)
       gright = g(iky,ikxmod(iseg+1,ie,iky),iz_low(iseg+1)+1:iz_low(iseg+1)+2)
    else
       ! apply periodic BC to zonal mode and zero BC otherwise
       if (zonal_mode(iky)) then
          gright = g(iky,ikxmod(iseg,ie,iky),iz_low(iseg)+1:iz_low(iseg)+2)
       else
          gright = 0.0
       end if
    end if
    
  end subroutine fill_zed_ghost_zones

  subroutine map_to_extended_zgrid (ie, iky, g, gext, ulim)

    use zgrid, only: nzgrid

    implicit none

    integer, intent (in) :: ie, iky
    complex, dimension (:,-nzgrid:), intent (in) :: g
    complex, dimension (:), intent (out) :: gext
    integer, intent (out) :: ulim

    integer :: iseg, ikx
    integer :: llim

    ! avoid double-counting at boundaries between 2pi segments
    iseg = 1
    ikx = ikxmod(iseg,ie,iky)
    llim = 1 ; ulim = nzed_segment+1
    gext(llim:ulim) = g(ikx,iz_low(iseg):iz_up(iseg))
    if (nsegments(ie,iky) > 1) then
       do iseg = 2, nsegments(ie,iky)
          ikx = ikxmod(iseg,ie,iky)
          llim = ulim+1
          ulim = llim+nzed_segment-1
          gext(llim:ulim) = g(ikx,iz_low(iseg)+1:iz_up(iseg))
       end do
    end if

  end subroutine map_to_extended_zgrid

  subroutine map_from_extended_zgrid (ie, iky, gext, g)

    use zgrid, only: nzgrid

    implicit none

    integer, intent (in) :: ie, iky
    complex, dimension (:), intent (in) :: gext
    complex, dimension (:,-nzgrid:), intent (in out) :: g

    integer :: iseg, ikx
    integer :: llim, ulim

    iseg = 1
    ikx = ikxmod(iseg,ie,iky)
    llim = 1 ; ulim = nzed_segment+1
    g(ikx,iz_low(iseg):iz_up(iseg)) = gext(llim:ulim)
    if (nsegments(ie,iky) > 1) then
       do iseg = 2, nsegments(ie,iky)
          llim = ulim+1
          ulim = llim+nzed_segment-1
          ikx = ikxmod(iseg,ie,iky)
          g(ikx,iz_low(iseg)) = gext(llim-1)
          g(ikx,iz_low(iseg)+1:iz_up(iseg)) = gext(llim:ulim)
       end do
    end if

  end subroutine map_from_extended_zgrid

  subroutine finish_extended_zgrid

    implicit none

    if (allocated(neigen)) deallocate (neigen)
    if (allocated(periodic)) deallocate (periodic)
    if (allocated(nsegments)) deallocate (nsegments)
    if (allocated(iz_low)) deallocate (iz_low, iz_mid, iz_up)
    if (allocated(ikxmod)) deallocate (ikxmod)

    extended_zgrid_initialized = .false.

  end subroutine finish_extended_zgrid

end module extended_zgrid
