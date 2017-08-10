module fields

  implicit none

  public :: init_fields, finish_fields

  private

  logical :: fields_initialized = .false.
  logical :: exist

  logical :: debug = .false.

contains

  subroutine init_fields

    use mp, only: proc0
    use fields_arrays, only: phi, apar
    use dist_fn_arrays, only: gvmu
    use stella_layouts, only: init_stella_layouts
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    use run_parameters, only: init_run_parameters
    use dist_fn, only: init_dist_fn
    use dist_fn, only: init_get_fields, get_fields
    use dist_fn, only: init_gxyz
    use init_g, only: ginit, init_init_g
!    use nonlinear_terms, only: nl_finish_init => finish_init

    implicit none

    logical :: restarted

    if (fields_initialized) return
    fields_initialized = .true.

    debug = debug .and. proc0
    
    if (debug) write(6,*) "fields::init_fields::init_theta_grid"
    call init_theta_grid
    if (debug) write(6,*) "fields::init_fields::init_init_g"
    call init_init_g
    if (debug) write(6,*) "fields::init_fields::init_run_parameters"
    call init_run_parameters
    if (debug) write(6,*) "fields::init_fields::init_dist_fn"
    call init_dist_fn
    if (debug) write(6,*) "fields::init_fields::allocate_arrays"
    call allocate_arrays
    if (debug) write(6,*) 'fields::init_fields::init_stella_layouts'
    call init_stella_layouts
    if (debug) write(6,*) 'fields::init_fields::init_kt_grids'
    call init_kt_grids

! Turn on nonlinear terms.
!    if (debug) write(6,*) "init_fields::nl_finish_init"
!    call nl_finish_init

    if (debug) write(*,*) "fields::init_fields::ginit"
    call ginit (restarted)
    if (debug) write(*,*) "fields::init_fields::init_gxyz"
    call init_gxyz

    if (restarted) return

    if (debug) write (*,*) 'fields::init_fields::init_get_fields'
    ! initialize get_fields subroutine
    call init_get_fields
    if (debug) write (*,*) 'fields::init_fields::get_fields'
    ! get initial field from initial distribution function
    call get_fields (gvmu, phi, apar)

  end subroutine init_fields

  subroutine allocate_arrays

    use fields_arrays, only: phi, apar
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0

    implicit none

    if (.not.allocated(phi)) then
       allocate (phi(naky,ntheta0,-ntgrid:ntgrid))
       phi = 0.
    end if
    if (.not. allocated(apar)) then
       allocate (apar(naky,ntheta0,-ntgrid:ntgrid))
       apar = 0.
    end if

  end subroutine allocate_arrays

!   subroutine phinorm (phifnc, aparfnc, bparfnc, phitot)
!     use theta_grid, only: delthet, ntgrid
!     use kt_grids, only: naky, ntheta0
!     use constants
!     implicit none
!     complex, dimension (-ntgrid:,:,:), intent (in) :: phifnc, aparfnc, bparfnc
!     real, dimension (:,:), intent (out) :: phitot
!     integer :: ik, it

!     do ik = 1, naky
!        do it = 1, ntheta0
!           phitot(it,ik) = 0.5/pi &
!            *(sum((abs(phifnc(:,it,ik))**2 + abs(aparfnc(:,it,ik))**2 &
!                   + abs(bparfnc(:,it,ik))**2) &
!                  *delthet))
!        end do
!     end do
!   end subroutine phinorm

!   subroutine kperp (ntgrid_output, akperp)
!     use theta_grid, only: delthet
!     use kt_grids, only: naky, aky, ntheta0
!     use run_parameters, only: fphi, fapar, fbpar
!     use dist_fn_arrays, only: kperp2
!     implicit none
!     integer, intent (in) :: ntgrid_output
!     real, dimension (:,:), intent (out) :: akperp
!     real :: anorm
!     integer :: ik, it

!     do ik = 1, naky
!        do it = 1, ntheta0
!           anorm = sum(abs(phinew(-ntgrid_output:ntgrid_output,it,ik)*fphi &
!                          + aparnew(-ntgrid_output:ntgrid_output,it,ik)*fapar &
!                          + bparnew(-ntgrid_output:ntgrid_output,it,ik)*fbpar)**2 &
!                       *delthet(-ntgrid_output:ntgrid_output))
!           if (anorm < 2.0*epsilon(0.0) .or. aky(ik) < epsilon(0.0)) then
!              akperp(it,ik) = 0.0
!           else
!              akperp(it,ik) &
!                   = sqrt(sum(kperp2(-ntgrid_output:ntgrid_output,it,ik) &
!                      *abs(phinew(-ntgrid_output:ntgrid_output,it,ik)*fphi &
!                           + aparnew(-ntgrid_output:ntgrid_output,it,ik)*fapar &
!                           + bparnew(-ntgrid_output:ntgrid_output,it,ik)*fbpar)**2 &
!                      *delthet(-ntgrid_output:ntgrid_output))/anorm)
!           end if
!        end do
!     end do
!   end subroutine kperp

  !!> This generates a flux surface average of phi. 

  !subroutine flux_surface_average_phi (phi_in, phi_average)
    !use theta_grid, only: ntgrid, drhodpsi, gradpar, bmag, delthet
    !use kt_grids, only: ntheta0, naky

    !implicit none
    !complex, intent (in) :: phi_in
    !complex, intent (out) :: phi_average
    !complex, dimension(-ntgrid:ntgrid,1:ntheta0,1:naky) :: phi_fieldline_avg
    !integer it, ig

    !call fieldline_average_phi(phi_in, phi_fieldline_avg)
    !do it = 1,ntheta0
      !do ig = -ntgrid,ntgrid
        !phi_average(ig, it, :) = sum(phi_fieldline_avg(ig, it, :))/real(naky)
      !end do
    !end do

  !end subroutine fieldline_average_phi

  ! subroutine timer
    
  !   character (len=10) :: zdate, ztime, zzone
  !   integer, dimension(8) :: ival
  !   real, save :: told=0., tnew=0.
    
  !   call date_and_time (zdate, ztime, zzone, ival)
  !   tnew = ival(5)*3600.+ival(6)*60.+ival(7)+ival(8)/1000.
  !   if (told > 0.) then
  !      print *, 'Fields: Time since last called: ',tnew-told,' seconds'
  !   end if
  !   told = tnew
  ! end subroutine timer

  subroutine finish_fields

    use fields_arrays, only: phi
    use fields_arrays, only: apar
!    use fields_arrays, only: bparold, bparnew
    use dist_fn, only: finish_get_fields

    implicit none

    call finish_get_fields
    if (allocated(phi)) deallocate (phi)
    if (allocated(apar)) deallocate (apar)
!    if (allocated(bparold)) deallocate (bparold)
!    if (allocated(bparnew)) deallocate (bparnew)

    fields_initialized = .false.

  end subroutine finish_fields

end module fields
