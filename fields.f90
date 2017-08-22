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
    use geometry, only: init_geometry
    use zgrid, only: init_zgrid
    use zgrid, only: nzed, nzgrid
    use zgrid, only: delthet, theta
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
    
    if (debug) write(6,*) "fields::init_fields::init_zgrid"
    call init_zgrid
    if (debug) write(6,*) "fields::init_fields::init_geometry"
    call init_geometry (nzed, nzgrid, theta, delthet)
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
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx

    implicit none

    if (.not.allocated(phi)) then
       allocate (phi(naky,nakx,-nzgrid:nzgrid))
       phi = 0.
    end if
    if (.not. allocated(apar)) then
       allocate (apar(naky,nakx,-nzgrid:nzgrid))
       apar = 0.
    end if

  end subroutine allocate_arrays

  subroutine finish_fields

    use fields_arrays, only: phi
    use fields_arrays, only: apar
    use geometry, only: finish_geometry
    use zgrid, only: finish_zgrid
    use dist_fn, only: finish_get_fields

    implicit none

    call finish_get_fields
    call finish_geometry
    call finish_zgrid
    if (allocated(phi)) deallocate (phi)
    if (allocated(apar)) deallocate (apar)

    fields_initialized = .false.

  end subroutine finish_fields

end module fields
