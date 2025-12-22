! ================================================================================================================================================================================= !
! --------------------------------------------------------------------- Routines for dealing with perioidic splines. -------------------------------------------------------------- !
! ================================================================================================================================================================================= !

module periodic_splines

    implicit none

    private

    public :: periodic_spline, new_periodic_spline, delete_periodic_spline

    interface handle_spline_error
        module procedure handle_spline_error_logical
    end interface handle_spline_error


    ! Holds data representing a periodic spline. Should be set up by calling [[new_periodic_spline]].
    type :: periodic_spline
        private
        
        ! Length of the data arrays represented by the spline.
        integer :: n = 0
        ! The actual size of the periodic domain.
        real :: period = 0
        ! Holds the independent and dependent values of the splined data in `x` and `y`. The second derivative is held in `y2` and calculated automatically.
        real, dimension (:), allocatable :: x, y, y2
        ! Indicates if the spline corresponding to this data is valid and can be used with the spline evaluation routines.
        logical, public :: valid = .false.
        ! The tension used in computing the splined data, note this must be the value used in the initialisation when passed to the spline evaluation routines.
        real :: tension = 1.0
  
        contains
            procedure :: interpolate => periodic_spline_interp
            procedure :: derivative => periodic_spline_deriv
    end type periodic_spline

    ! Constructor for periodic_spline.

    interface periodic_spline
        module procedure new_periodic_spline
    end interface periodic_spline


contains

    ! Populates a periodic_spline instance `spl` representing the periodic data y(x) of length n and periodic on `period`. 
    ! Note that the spline library expects `period > x(n) - x(1)`, which means the input data shouldn't include the duplicate periodic point.  
    ! As a convenience the user can pass data with the duplicate point and set `drop_last_point = .true.` to automatically exclude the duplicate point.
  
    type(periodic_spline) function new_periodic_spline (x, y, period, drop_last_point, tension) result(spl)
        use optionals, only: get_option_with_default
        use splines, only: fitp_curvp1
    
        implicit none
        
        real, dimension (:), intent (in) :: x, y
        real, intent (in) :: period
        logical, intent(in), optional :: drop_last_point
        real, intent(in), optional :: tension
        logical :: drop_point
        real, dimension (:), allocatable :: temp
        integer :: ierr
        spl%valid = .false.
        drop_point = get_option_with_default(drop_last_point, .false.)
        spl%tension = get_option_with_default(tension, 1.0)
        spl%n = size(x)
    
        if (drop_point) spl%n = spl%n - 1
        allocate(spl%x, source = x(:spl%n)) ; allocate(spl%y, source = y(:spl%n))
        allocate (spl%y2(spl%n), temp(2*spl%n))
        spl%period = period
        call fitp_curvp1 (spl%n,spl%x,spl%y,spl%period,spl%y2,temp,spl%tension,ierr)
        spl%valid = ierr == 0
    end function new_periodic_spline


    ! Reset and deallocate variables in passed periodic spline.

    subroutine delete_periodic_spline (spl)
        implicit none
    
        type (periodic_spline), intent (in out) :: spl
        spl%n = 0
        spl%period = 0.0
        if (allocated(spl%x)) deallocate (spl%x,spl%y)
        if (allocated(spl%y2)) deallocate (spl%y2)
        spl%valid = .false.
        spl%tension = 1.0
    end subroutine delete_periodic_spline


  !> Bound wrapper to splint.
    real function periodic_spline_interp(self, x)
        implicit none
    
        class (periodic_spline), intent(in) :: self
        real, intent(in) :: x
        
        periodic_spline_interp = periodic_splint(x, self)
    end function periodic_spline_interp


  !> Bound wrapper to dsplint.
    real function periodic_spline_deriv(self, x)
        implicit none
    
        class (periodic_spline), intent(in) :: self
        real, intent(in) :: x
        
        periodic_spline_deriv = periodic_dsplint(x, self)
    end function periodic_spline_deriv


    real function periodic_splint (x, spl)
        use splines, only: fitp_curvp2
        implicit none
    
        real, intent (in) :: x
        type (periodic_spline), intent (in) :: spl
    
        call handle_spline_error(spl%valid, 'periodic_splint')
        periodic_splint = fitp_curvp2(x, spl%n, spl%x, spl%y, spl%period, spl%y2, spl%tension)
    end function periodic_splint


    real function periodic_dsplint (x, spl)
        use splines, only: fitp_curvpd     

        implicit none
    
        real, intent (in) :: x
        type (periodic_spline), intent (in) :: spl
    
        call handle_spline_error(spl%valid, 'periodic_dsplint')
        periodic_dsplint = fitp_curvpd(x, spl%n, spl%x, spl%y, spl%period, spl%y2, spl%tension)
    end function periodic_dsplint


    !> If not valid abort with error noting invalid spline and which method was invoked.
    subroutine handle_spline_error_logical(valid, routine_name)
        use mp, only: mp_abort
    
        implicit none
    
        logical, intent(in) :: valid
        character(len = *), intent(in) :: routine_name
    
        if (.not. valid) call mp_abort('Attempt to use invalid spline in '//routine_name)
    end subroutine handle_spline_error_logical

! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------------- fitp_curvp1. --------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------------- fitp_curvp2. --------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------------- End Module. ---------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

end module periodic_splines
