module sfincs_interface

  implicit none
  
  public :: get_neo_from_sfincs

  private

  integer :: nproc_sfincs
  logical :: includeXDotTerm
  logical :: includeElectricFieldTermInXiDot
  integer :: magneticDriftScheme
  logical :: includePhi1
  logical :: includePhi1InKineticEquation
  integer :: geometryScheme
  integer :: coordinateSystem
  integer :: inputRadialCoordinateForGradients

contains

  subroutine get_neo_from_sfincs

# ifdef USE_SFINCS    
    use mp, only: proc0, iproc
    use mp, only: comm_split, comm_free
    use sfincs_main, only: init_sfincs, prepare_sfincs, run_sfincs
# else
    use mp, only: mp_abort
# endif
    
    implicit none
    
# ifdef USE_SFINCS
    integer :: sfincs_comm
    integer :: color, ierr

    if (proc0) call read_sfincs_parameters
    call broadcast_sfincs_parameters
    if (iproc < nproc_sfincs) then
       color = 0
    else
       color = 1
    end if
    call comm_split (color, sfincs_comm, ierr)
    if (iproc < nproc_sfincs) then
       call init_sfincs (sfincs_comm)
       call pass_geometry_to_sfincs
       call prepare_sfincs
       call run_sfincs
    end if
    call comm_free (sfincs_comm, ierr)
    ! NB: NEED TO BROADCAST SFINCS RESULTS
# else
    call mp_abort ('to run with include_neoclassical_terms=.true., &
         & USE_SFINCS must be defined at compilation time.  Aborting.')
# endif

  end subroutine get_neo_from_sfincs

# ifdef USE_SFINCS  
  subroutine read_sfincs_parameters

    use mp, only: nproc
    use file_utils, only: input_unit_exist
    use species, only: nspec

    implicit none

    namelist /sfincs_input/ nproc_sfincs, &
         includeXDotTerm, &
         includeElectricFieldTermInXiDot, &
         magneticDriftScheme, &
         includePhi1, &
         includePhi1InKineticEquation, &
         geometryScheme, &
         coordinateSystem, &
         inputRadialCoordinateForGradients

    logical :: exist
    integer :: in_file

    nproc_sfincs = 1
    includeXDotTerm = .false.
    includeElectricFieldTermInXiDot = .false.
    magneticDriftScheme = 0
    includePhi1 = .true.
    includePhi1InKineticEquation = .false.
    geometryScheme = 1
    coordinateSystem = 3
    ! option 1 corresponds to using toroidal flux
    ! normalized by toroidal flux enclosed by the LCFS
    ! when calculating gradients of density, temperature, and potential
    inputRadialCoordinateForGradients = 1
    
    in_file = input_unit_exist("sfincs_input", exist)
    if (exist) read (unit=in_file, nml=sfincs_input)

    if (nproc_sfincs > nproc) then
       write (*,*) 'requested number of processors for sfincs is greater &
            & than total processor count.'
       write (*,*) 'allocating ', nproc, ' processors for sfincs.'
    end if

    if (nspec == 1 .and. includePhi1) then
       write (*,*) 'includePhi1 = .true. is incompatible with a single-species run.'
       write (*,*) 'forcing includePhi1 = .false.'
       includePhi1 = .false.
    end if

  end subroutine read_sfincs_parameters

  subroutine broadcast_sfincs_parameters

    use mp, only: broadcast

    implicit none

    call broadcast (nproc_sfincs)
    call broadcast (includeXDotTerm)
    call broadcast (includeElectricFieldTermInXiDot)
    call broadcast (magneticDriftScheme)
    call broadcast (includePhi1)
    call broadcast (includePhi1InKineticEquation)
    call broadcast (geometryScheme)
    call broadcast (coordinateSystem)
    call broadcast (inputRadialCoordinateForGradients)
    
  end subroutine broadcast_sfincs_parameters

  subroutine pass_geometry_to_sfincs

    use species, only: spec, nspec
    use theta_grid, only: ntgrid, drhotor2dr
    use globalVariables, only: includeXDotTerm_sfincs => includeXDotTerm
    use globalVariables, only: includeElectricFieldTermInXiDot_sfincs => includeElectricFieldTermInXiDot
    use globalVariables, only: magneticDriftScheme_sfincs => magneticDriftScheme
    use globalVariables, only: includePhi1_sfincs => includePhi1
    use globalVariables, only: includePhi1InKineticEquation_sfincs => includePhi1InKineticEquation
    use globalVariables, only: geometryScheme_sfincs => geometryScheme
    use globalVariables, only: coordinateSystem_sfincs => coordinateSystem
    use globalVariables, only: RadialCoordinateForGradients => inputRadialCoordinateForGradients
    use globalVariables, only: Nspecies, nHats, THats, MHats, Zs
    use globalVariables, only: Nzeta, Ntheta
    use globalVariables, only: dnHatdpsiHats, dTHatdpsiHats, dPhiHatdpsiHat
    
    implicit none

    includeXDotTerm_sfincs = includeXDotTerm
    includeElectricFieldTermInXiDot_sfincs = includeElectricFieldTermInXiDot
    magneticDriftScheme_sfincs = magneticDriftScheme
    includePhi1_sfincs = includePhi1
    includePhi1InKineticEquation_sfincs = includePhi1InKineticEquation
    geometryScheme_sfincs = geometryScheme
    coordinateSystem_sfincs = coordinateSystem
!     RadialCoordinateForGradients = inputRadialCoordinateForGradients
    Nspecies = nspec
    nHats(:nspec) = spec%dens
    THats(:nspec) = spec%temp
    mHats(:nspec) = spec%mass
    Zs(:nspec) = spec%z
!     ! FLAG -- need to modify for stellarator simulations
!     ! I think nzeta will be 2*nzgrid+1
!     ! and ntheta will be ny_ffs
    Nzeta = 1
    Ntheta = 2*ntgrid+1
!     ! FLAG -- this currently assumes we are using Miller specification
!    if (inputRadialCoordinateForGradients == 1) then
!        ! radial density gradient with respect to psitor/psitor_LCFS
!        ! normalized by reference density (not species density)
!        dnHatdpsiHats(:nspec) = -spec%fprim*spec%dens/drhotor2dr
!        ! radial temperature gradient with respect to psitor/psitor_LCFS = rhotor2
!        ! normalized by reference tmperatures (not species temperature)
!        dTHatdpsiHats(:nspec) = -spec%tprim*spec%temp/drhotor2dr
!        ! radial electric field
!        dPhiHatdpsiHat = 0.0
!     end if

  end subroutine pass_geometry_to_sfincs
# endif

end module sfincs_interface
