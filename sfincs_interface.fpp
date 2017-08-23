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
  integer :: inputRadialCoordinate
  integer :: inputRadialCoordinateForGradients
  real :: aHat, psiAHat
  real :: nu_n

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
       call pass_inputoptions_to_sfincs
       call prepare_sfincs
       call pass_geometry_to_sfincs
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
         inputRadialCoordinate, &
         inputRadialCoordinateForGradients, &
         aHat, psiAHat, nu_N

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
    ! option 3 corresponds to using sqrt of toroidal flux
    ! normalized by toroidal flux enclosed by the LCFS
    inputRadialCoordinate = 3
    ! option 3 corresponds to same choice
    ! when calculating gradients of density, temperature, and potential
    inputRadialCoordinateForGradients = 3
    ! corresponds to r_LCFS as reference length in sfincs
    aHat = 1.0
    ! corresponds to psitor_LCFS = B_ref * a_ref^2
    psiAHat = 1.0
    ! nu_n = nu_ref * aref/vt_ref
    ! nu_ref = 4*sqrt(2*pi)*nref*e**4*loglam/(3*sqrt(mref)*Tref**3/2)
    ! (with nref, Tref, and mref in Gaussian units)
    nu_N = 0.01 ! FLAG -- should replace with collisionality from stella

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
    call broadcast (inputRadialCoordinate)
    call broadcast (inputRadialCoordinateForGradients)
    call broadcast (aHat)
    call broadcast (psiAHat)
    call broadcast (nu_N)

  end subroutine broadcast_sfincs_parameters

  subroutine pass_inputoptions_to_sfincs

    use mp, only: mp_abort
    use geometry, only: geo_surf
    use species, only: spec, nspec
    use zgrid, only: nzgrid
    use globalVariables, only: includeXDotTerm_sfincs => includeXDotTerm
    use globalVariables, only: includeElectricFieldTermInXiDot_sfincs => includeElectricFieldTermInXiDot
    use globalVariables, only: magneticDriftScheme_sfincs => magneticDriftScheme
    use globalVariables, only: includePhi1_sfincs => includePhi1
    use globalVariables, only: includePhi1InKineticEquation_sfincs => includePhi1InKineticEquation
    use globalVariables, only: geometryScheme_sfincs => geometryScheme
    use globalVariables, only: coordinateSystem_sfincs => coordinateSystem
    use globalVariables, only: RadialCoordinate => inputRadialCoordinate
    use globalVariables, only: RadialCoordinateForGradients => inputRadialCoordinateForGradients
    use globalVariables, only: rN_wish
    use globalVariables, only: Nspecies, nHats, THats, MHats, Zs
    use globalVariables, only: Nzeta, Ntheta
    use globalVariables, only: dnHatdrNs, dTHatdrNs, dPhiHatdrN
    use globalVariables, only: aHat_sfincs => aHat
    use globalVariables, only: psiAHat_sfincs => psiAHat
    use globalVariables, only: nu_n_sfincs => nu_n

    implicit none

    includeXDotTerm_sfincs = includeXDotTerm
    includeElectricFieldTermInXiDot_sfincs = includeElectricFieldTermInXiDot
    magneticDriftScheme_sfincs = magneticDriftScheme
    includePhi1_sfincs = includePhi1
    includePhi1InKineticEquation_sfincs = includePhi1InKineticEquation
    geometryScheme_sfincs = geometryScheme
    coordinateSystem_sfincs = coordinateSystem
    RadialCoordinate = inputRadialCoordinate
    RadialCoordinateForGradients = inputRadialCoordinateForGradients
    Nspecies = nspec
    nHats(:nspec) = spec%dens
    THats(:nspec) = spec%temp
    mHats(:nspec) = spec%mass
    Zs(:nspec) = spec%z
!     ! FLAG -- need to modify for stellarator simulations
!     ! I think nzeta will be 2*nzgrid+1
!     ! and ntheta will be ny_ffs
    Nzeta = 1
    Ntheta = 2*nzgrid+1
    aHat_sfincs = aHat
    psiAHat_sfincs = psiAHat
    nu_n_sfincs = nu_n

    if (inputRadialCoordinate == 3) then
       rN_wish = geo_surf%rhotor
    else
       call mp_abort ('only inputRadialCoordinate=3 currently supported. aborting.')
    end if
    if (inputRadialCoordinateForGradients == 3) then
       ! radial density gradient with respect to rhotor = sqrt(psitor/psitor_LCFS)
       ! normalized by reference density (not species density)
       dnHatdrNs(:nspec) = -spec%fprim*spec%dens/geo_surf%drhotordrho
       ! radial temperature gradient with respect to rhotor = sqrt(psitor/psitor_LCFS)
       ! normalized by reference tmperatures (not species temperature)
       dTHatdrNs(:nspec) = -spec%tprim*spec%temp/geo_surf%drhotordrho
       ! radial electric field
       dPhiHatdrN = 0.0
    else
       call mp_abort ('only inputRadialCoordinateForGradients=3 currently supported. aborting.')
    end if

  end subroutine pass_inputoptions_to_sfincs

  subroutine pass_geometry_to_sfincs

    use geometry, only: bmag, dbdthet, gradpar
    use geometry, only: geo_surf
    use globalVariables, only: BHat
    use globalVariables, only: dBHatdtheta
    use globalVariables, only: iota
    use globalVariables, only: DHat
    use globalVariables, only: BHat_sup_theta
    use globalVariables, only: BHat_sub_zeta

    implicit none

    integer :: nzeta = 1

    call init_zero_arrays

    ! FLAG -- needs to be changed for stellarator runs
    BHat = spread(bmag,2,nzeta)
    dBHatdtheta = spread(dbdthet,2,nzeta)
    iota = 1./geo_surf%qinp
    ! this is grad psitor . (grad theta x grad zeta)
    ! note that + sign below relies on B = I grad zeta + grad zeta x grad psi
    DHat = geo_surf%qinp*spread(bmag*gradpar,2,nzeta)
    ! this is bhat . grad theta
    BHat_sup_theta = spread(bmag*gradpar,2,nzeta)
    ! this is I(psi) / (aref*Bref)
    BHat_sub_zeta = geo_surf%rgeo

  end subroutine pass_geometry_to_sfincs

  subroutine init_zero_arrays
    use globalVariables, only: dBHatdzeta
    use globalVariables, only: dBHatdpsiHat
    use globalVariables, only: BHat_sup_zeta
    use globalVariables, only: BHat_sub_psi
    use globalVariables, only: BHat_sub_theta
    use globalVariables, only: dBHat_sub_psi_dtheta
    use globalVariables, only: dBHat_sub_psi_dzeta
    use globalVariables, only: dBHat_sub_theta_dpsiHat
    use globalVariables, only: dBHat_sub_theta_dzeta
    use globalVariables, only: dBHat_sub_zeta_dpsiHat
    use globalVariables, only: dBHat_sub_zeta_dtheta
    use globalVariables, only: dBHat_sup_theta_dpsiHat
    use globalVariables, only: dBHat_sup_theta_dzeta
    use globalVariables, only: dBHat_sup_zeta_dpsiHat
    use globalVariables, only: dBHat_sup_zeta_dtheta
    implicit none
    dBHatdzeta = 0.
    dBHatdpsiHat = 0.
    BHat_sup_zeta = 0.
    BHat_sub_psi = 0.
    BHat_sub_theta = 0.
    dBHat_sub_psi_dtheta = 0.
    dBHat_sub_psi_dzeta = 0.
    dBHat_sub_theta_dpsiHat = 0.
    dBHat_sub_theta_dzeta = 0.
    dBHat_sub_zeta_dpsiHat = 0.
    dBHat_sub_zeta_dtheta = 0.
    dBHat_sup_theta_dpsiHat = 0.
    dBHat_sup_theta_dzeta = 0.
    dBHat_sup_zeta_dpsiHat = 0.
    dBHat_sup_zeta_dtheta = 0.
  end subroutine init_zero_arrays
# endif

end module sfincs_interface
