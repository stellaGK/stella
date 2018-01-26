module extended_zgrid

  implicit none

  public :: nsegments, nsegments_poskx
  public :: neigen
  public :: ikxmod
  public :: iz_low, iz_mid, iz_up
  public :: periodic
  public :: nzed_segment

  ! these arrays needed to keep track of connections between different
  ! 2pi segments
  integer :: nzed_segment
  integer, dimension (:), allocatable :: neigen
  integer, dimension (:), allocatable :: iz_low, iz_mid, iz_up
  integer, dimension (:,:), allocatable :: nsegments, nsegments_poskx
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
    use kt_grids, only: ntheta0, nakx, naky
    use kt_grids, only: jtwist_out, aky
    use species, only: nspec
    use vpamu_grids, only: nmu, nvgrid

    implicit none

    integer :: iseg, iky, ie, ntg, ikx, ikxshiftend
    integer :: nseg_max, neigen_max
    integer :: iky_max
    integer, dimension (:), allocatable :: ikx_shift_left_kypos, ikx_shift_left_kyneg
    integer, dimension (:,:), allocatable :: ikx_shift

    if (extended_zgrid_initialized) return
    extended_zgrid_initialized = .true.
    
    ntg = nzed/2
    ! iky_max is the index of the most positive ky
    iky_max = naky/2+1

    if (.not. allocated(neigen)) allocate (neigen(naky))
    if (.not. allocated(periodic)) allocate (periodic(naky)) ; periodic = .false.

    if (boundary_option_switch==boundary_option_self_periodic) then
       periodic = .true.
    else
       where (abs(aky) < epsilon(0.0)) periodic = .true.
    end if

    select case (boundary_option_switch)
    case (boundary_option_linked)

       ! ntheta0 = 2*(nakx-1) + 1
       ! nakx includes kx >= 0
       ! ntheta0 also includes kx < 0
       neigen(1) = ntheta0
       if (naky > 1) then
          do iky = 2, iky_max
             ! must link different kx values at theta = +/- pi
             ! neigen is the number of independent eigenfunctions along the field line
             neigen(iky) = min((iky-1)*jtwist_out,ntheta0)
          end do
          ! number of eigenfunctions for -ky is same as for +ky
          neigen(iky_max+1:) = neigen(iky_max:2:-1)
       end if

       neigen_max = maxval(neigen)

       if (.not. allocated(ikx_shift_left_kypos)) then
          allocate (ikx_shift_left_kypos(neigen_max)) ; ikx_shift_left_kypos = 0
          allocate (ikx_shift_left_kyneg(neigen_max)) ; ikx_shift_left_kyneg = 0
          allocate (ikx_shift(ntheta0,naky)) ; ikx_shift = 0
       end if

       ! figure out how much to shift ikx by to get to
       ! the left-most (theta-theta0) in each set of connected 2pi segments
       ! note that theta0 goes from 0 to theta0_max and then from theta0_min back
       ! to -dtheta0
       do ikx = 1, neigen_max
          ! first ntheta0/2+1 theta0s are 0 and all positive theta0 values
          ! remainder are negative theta0s
          ! note that ntheta0 is always positive for box
          ! theta_0 = kx / ky / shat
          ! if ky > 0, then most positive theta_0 corresponds to most positive kx
          if (ikx <= nakx) then
             ikx_shift_left_kypos(ikx) = nakx-2*ikx+1
          else
             ikx_shift_left_kypos(ikx) = 3*nakx-2*ikx
          end if
          ! if ky < 0, most positive theta_0 corresponds to most negative kx
          if (ikx < nakx) then
             ikx_shift_left_kyneg(ikx) = nakx
          else
             ikx_shift_left_kyneg(ikx) = 1-nakx
          end if
       end do

       do iky = 1, naky
          ! ikx_shift is how much to shift each ikx by to connect
          ! to the next theta0 (from most positive to most negative)

          ! if ky < 0, then going to more negative theta0
          ! corresponds to going to more positive kx
          if (aky(iky) < 0.0) then
             ! first treat kx positive
             ! connect to more positive kx neigen away
             ! but only if not trying to connect to kx
             ! so positive that ikx is not on grid
             if (nakx - neigen(iky) > 0) ikx_shift(:nakx-neigen(iky),iky) = neigen(iky)
             ! next treat kx negative
             ! if kx sufficiently negative, then 
             ! shifting by neigen keeps kx negative
             do ikx = nakx+1, ntheta0
                if (ikx+neigen(iky) <= ntheta0) then
                   ikx_shift(ikx,iky) = neigen(iky)
                   ! if theta0 not sufficiently negative,
                   ! then must shift to postive theta0
                else if (ikx+neigen(iky) <= ntheta0+nakx) then
                   ikx_shift(ikx,iky) = neigen(iky) - ntheta0
                end if
             end do
          else
             ! if ky > 0, then going to more negative theta0
             ! corresponds to going to more negative kx
             do ikx = 1, nakx
                ! if theta0 is sufficiently positive, shifting to more
                ! negative theta0 corresponds to decreasing ikx
                if (ikx-neigen(iky) > 0) then
                   ikx_shift(ikx,iky) = -neigen(iky)
                   ! if a positive theta0 connects to a negative theta0
                   ! must do more complicated mapping of ikx
                else if (ikx-neigen(iky)+ntheta0 >= nakx+1) then
                   ikx_shift(ikx,iky) = ntheta0 - neigen(iky)
                end if
             end do
             ! if theta0 is negative, then shifting to more negative
             ! theta0 corresponds to decreasing ikx
             do ikx = nakx+1, ntheta0
                ! if theta0 is sufficiently negative, it has no
                ! more negative theta0 with which it can connect
                if (ikx-neigen(iky) >= nakx) then
                   ikx_shift(ikx,iky) = -neigen(iky)
                end if
                ! theta0 is positive
             end  do
          end if
       end do

       if (.not. allocated(nsegments)) then
          allocate (nsegments(neigen_max,naky))
          allocate (nsegments_poskx(neigen_max,naky))
       end if

       do iky = 1, naky
          if (neigen(iky) == 0) then
             nsegments(:,iky) = 1
             nsegments_poskx(:,iky) = 1
          else
             nsegments(:,iky) = (ntheta0-1)/neigen(iky)
             nsegments_poskx(:,iky) = (nakx-1)/neigen(iky)

             do ie = 1, mod(ntheta0-1,neigen(iky))+1
                nsegments(ie,iky) = nsegments(ie,iky) + 1
!                nsegments_poskx(ie,iky) = nsegments_poskx(ie,iky) + 1
             end do
             do ie = 1, mod(nakx-1,neigen(iky))+1
                nsegments_poskx(ie,iky) = nsegments_poskx(ie,iky) + 1
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
       
       neigen = ntheta0 ; neigen_max = ntheta0
       
       if (.not. allocated(ikx_shift_left_kypos)) then
          allocate (ikx_shift_left_kypos(neigen_max))
          allocate (ikx_shift_left_kyneg(neigen_max))
          allocate (ikx_shift(ntheta0,naky))
       end if
       ikx_shift = 0 ; ikx_shift_left_kypos = 0 ; ikx_shift_left_kyneg = 0
       
       if (.not. allocated(nsegments)) then
          allocate (nsegments(neigen_max,naky))
       end if
       
       ! this is the number of 2pi poloidal segments in the extended theta domain,
       ! which is needed in initializing the reponse matrix and doing the implicit sweep
       nsegments = 2*(nperiod-1) + 1
       nsegments_poskx = nsegments

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
       ! initialize ikxmod to ntheta0
       ! should not be necessary but just in case one tries to access
       ! a value beyond nsegments(ie,iky)
       ikxmod = ntheta0
    end if
    do iky = 1, naky
       ! only do the following once for each independent set of theta0s
       ! the assumption here is that all kx are on processor and sequential
       do ie = 1, neigen(iky)
          if (aky(iky) < 0.) then
             ikxshiftend = ikx_shift_left_kyneg(ie)
          else
             ikxshiftend = ikx_shift_left_kypos(ie)
          end if
          ! remap to start at theta0 = theta0_max
          ! (so that theta-theta0 is most negative)
          ! for this set of connected theta0s
          iseg = 1
          ikxmod(iseg,ie,iky) = ie + ikxshiftend
          if (nsegments(ie,iky) > 1) then
             do iseg = 2, nsegments(ie,iky)
                ikxmod(iseg,ie,iky) = ikxmod(iseg-1,ie,iky) + ikx_shift(ikxmod(iseg-1,ie,iky),iky)
             end do
          end if
       end do
    end do

    ! quick hack
    if (.not. allocated(nsegments_poskx)) then
       allocate (nsegments_poskx(neigen_max,naky))
    end if
    nsegments_poskx = 0
    do iky = 1, naky
       do ie = 1, neigen(iky)
          do iseg = 1, nsegments(ie,iky)
             if (ikxmod(iseg,ie,iky) <= nakx) nsegments_poskx(ie,iky) = nsegments_poskx(ie,iky) + 1
          end do
       end do
    end do
    
    if (allocated(ikx_shift_left_kypos)) deallocate (ikx_shift_left_kypos)
    if (allocated(ikx_shift_left_kyneg)) deallocate (ikx_shift_left_kyneg)
    if (allocated(ikx_shift)) deallocate (ikx_shift)

    ! this is the number of unique zed values in all segments but the first
    ! the first has one extra unique zed value (all others have one grid common
    ! with the previous segment due to periodicity)
    nzed_segment = iz_up(1)-iz_low(1)

  end subroutine init_extended_zgrid

  subroutine finish_extended_zgrid

    implicit none

    if (allocated(neigen)) deallocate (neigen)
    if (allocated(periodic)) deallocate (periodic)
    if (allocated(nsegments)) deallocate (nsegments)
    if (allocated(nsegments_poskx)) deallocate (nsegments_poskx)
    if (allocated(iz_low)) deallocate (iz_low, iz_mid, iz_up)
    if (allocated(ikxmod)) deallocate (ikxmod)

    extended_zgrid_initialized = .false.

  end subroutine finish_extended_zgrid

end module extended_zgrid
