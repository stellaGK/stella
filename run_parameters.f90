!> This module is basically a store for the input parameters that are specified in the namelists \a knobs and \a parameters. In general, the names of the public variables in this module are the same as the name of the input parameter they correspond to.
 

module run_parameters

  implicit none

  public :: init_run_parameters, finish_run_parameters
  public :: beta, zeff, tite, nine
  public :: fphi, fapar, fbpar
  public :: nonlinear
  public :: code_delt_max, wunits, woutunits, tunits
  public :: nstep, wstar_units, eqzip, margin
  public :: secondary, tertiary, harris
  public :: ieqzip
  public :: k0, cfl_cushion
  public :: avail_cpu_time
  
  private

  real :: cfl_cushion
  real :: beta, zeff, tite, nine
  real :: fphi, fapar, fbpar
  real :: delt, code_delt_max, margin
  logical :: nonlinear
  real, dimension (:), allocatable :: wunits, woutunits, tunits
  real :: avail_cpu_time
  integer :: nstep
  logical :: wstar_units, eqzip
  logical :: secondary, tertiary, harris
  real :: k0
  integer, public :: delt_option_switch
  integer, public, parameter :: delt_option_hand = 1, delt_option_auto = 2
  logical :: initialized = .false.
  logical :: rpexist, knexist

  integer, allocatable :: ieqzip(:,:)
  integer :: eqzip_option_switch
  integer, parameter :: &
       eqzip_option_none = 1, &
       eqzip_option_secondary = 2, &
       eqzip_option_tertiary = 3, &
       eqzip_option_equilibrium = 4

contains

  subroutine init_run_parameters

    use kt_grids, only: init_kt_grids, naky, nakx
    use stella_time, only: init_delt
    
    implicit none
!    logical, save :: initialized = .false.

    if (initialized) return
    initialized = .true.

    call read_parameters

    call init_kt_grids
    call init_delt (delt)

    allocate (wunits(naky))
    allocate (woutunits(naky))
    allocate (tunits(naky))

    ! omega_* normalization of time: 
    call adjust_time_norm

    if(.not.allocated(ieqzip)) allocate(ieqzip(nakx,naky))
    ieqzip(1:nakx,1:naky)=1
    select case (eqzip_option_switch)
    case (eqzip_option_secondary)
       ! suppress evolution of secondary mode
       ieqzip(1,2) = 0
    case (eqzip_option_tertiary)
       ! suppress evolution of tertiary mode
       ieqzip(2,1) = 0
       ieqzip(nakx,1) = 0
    case (eqzip_option_equilibrium)
       ! suppress evolution of 1D equilibrium (x dependent)
       ieqzip(1:nakx,1) = 0
    end select
  end subroutine init_run_parameters

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use mp, only: proc0, broadcast
    use stella_save, only: init_dt
    use text_options, only: text_option, get_option_value
    implicit none
    type (text_option), dimension (4), parameter :: eqzipopts = &
         (/ text_option('none', eqzip_option_none), &
            text_option('secondary', eqzip_option_secondary), &
            text_option('tertiary', eqzip_option_tertiary), &
            text_option('equilibrium', eqzip_option_equilibrium) /)
    character (len=20) :: eqzip_option
    type (text_option), dimension (3), parameter :: deltopts = &
         (/ text_option('default', delt_option_hand), &
            text_option('set_by_hand', delt_option_hand), &
            text_option('check_restart', delt_option_auto) /)
    character(20) :: delt_option

    integer :: ierr, istatus, in_file
    real :: delt_saved

    real :: teti  ! for back-compatibility
    namelist /parameters/ beta, zeff, tite, nine, teti, k0
    namelist /knobs/ fphi, fapar, fbpar, delt, nstep, wstar_units, eqzip, &
         delt_option, margin, secondary, tertiary, harris, &
         avail_cpu_time, eqzip_option, nonlinear, cfl_cushion

    if (proc0) then
       fphi = 1.0
       fapar = 1.0
       fbpar = -1.0
       nonlinear = .false.
       beta = 0.0
       zeff = 1.0
       tite = 1.0
       teti = -100.0
       nine = 1.0
       wstar_units = .false.
       eqzip_option = 'none'
       eqzip = .false.
       secondary = .true.
       tertiary = .false.
       harris = .false.
       k0 = 1.
       delt_option = 'default'
       margin = 0.05
       avail_cpu_time = 1.e10
       cfl_cushion = 0.5

       in_file = input_unit_exist("parameters", rpexist)
       if (rpexist) read (unit=in_file,nml=parameters)

       in_file = input_unit_exist("knobs", knexist)
       if (knexist) read (unit=in_file, nml=knobs)

       if (abs(teti+100.) > epsilon(0.)) tite = teti

       if (eqzip) then
          if (secondary .and. tertiary) then
             ierr = error_unit()
             write (ierr, *) 'Forcing secondary = FALSE'
             write (ierr, *) 'because you have chosen tertiary = TRUE'
             secondary = .false.
          end if
          if (secondary .and. harris) then
             ierr = error_unit()
             write (ierr, *) 'Forcing secondary = FALSE'
             write (ierr, *) 'because you have chosen harris = TRUE'
             secondary = .false.
          end if
          if (tertiary .and. harris) then
             ierr = error_unit()
             write (ierr, *) 'Forcing tertiary = FALSE'
             write (ierr, *) 'because you have chosen harris = TRUE'
             tertiary = .false.
          end if
       endif

       ierr = error_unit()
       call get_option_value &
            (delt_option, deltopts, delt_option_switch, ierr, &
            "delt_option in knobs")

       call get_option_value ( &
            eqzip_option, eqzipopts, eqzip_option_switch, error_unit(), &
            "eqzip_option in knobs")

    end if

    call broadcast (delt_option_switch)
    call broadcast (delt)
    call broadcast (cfl_cushion)
    call broadcast (beta)
    call broadcast (zeff)
    call broadcast (tite)
    call broadcast (nine)
    call broadcast (fphi)
    call broadcast (fapar)
    call broadcast (fbpar)
    call broadcast (nonlinear)
    call broadcast (nstep)
    call broadcast (wstar_units)
    call broadcast (eqzip)
    call broadcast (secondary)
    call broadcast (tertiary)
    call broadcast (harris)
    call broadcast (margin)
    call broadcast (k0)
    call broadcast (avail_cpu_time)
    call broadcast (eqzip_option_switch)
    
    code_delt_max = delt

    delt_saved = delt
    if (delt_option_switch == delt_option_auto) then
       call init_dt (delt_saved, istatus)
       if (istatus == 0) delt  = delt_saved
    endif

  end subroutine read_parameters

  subroutine adjust_time_norm
    use file_utils, only: error_unit
    use mp, only: proc0
    use kt_grids, only: aky
    implicit none
!CMR: Sep 2010
! Attempt to understand time normalisation variables, which are arrays(naky)
!    TUNITS: DT(KY)=TUNITS(KY).CODE_DT
!            This is a generally very useful variable to store ky dependent 
!            timestep in the code time normalisation.
!            Used to multiply ky independent source terms on RHS of GKE.
!    WUNITS: WUNITS(KY)=AKY(KY)*TUNITS(KY)/2
!            Auxiliary variable.  Used to save compute operations when 
!            evaluating source terms on RHS of GKE that are proportional to ky.
!            !! The Mysterious factor 1/2 Explained !!
!            The factor of 1/2 arises because those source terms were first
!            specified using the normalisation Tref=mref vtref^2 
! [R Numata et al, "AstroGK: Astrophysical gyrokinetics code", JCP, 2010].
!CMRend
    if (wstar_units) then
       wunits = 1.0
       where (aky > epsilon(0.0))
          tunits = 2.0/aky
       elsewhere
          tunits = 0.0
       end where
       if (any(abs(tunits) < epsilon(0.)) .and. proc0) then
          write (error_unit(), *) &
               "WARNING: wstar_units=.true. and aky=0.0: garbage results"
          print *, &
               "WARNING: wstar_units=.true. and aky=0.0: garbage results"
       end if
!CMR: Sep 2010
!  Changes to allow wstar_units to be used
!CMRend
    else
       tunits = 1.0
       wunits = aky/2.0
    end if
    woutunits = 1.0/tunits
  end subroutine adjust_time_norm

  subroutine finish_run_parameters

    implicit none

    if (allocated(wunits)) deallocate (wunits, woutunits, tunits)

    initialized = .false.

  end subroutine finish_run_parameters

end module run_parameters
