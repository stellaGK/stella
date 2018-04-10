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
  real :: aHat, psiAHat, Delta
  real :: nu_n
  integer :: nxi, nx, ntheta

contains

  subroutine get_neo_from_sfincs (nradii, drho, f_neoclassical, phi_neoclassical)

# ifdef USE_SFINCS
    use geometry, only: geo_surf
    use mp, only: proc0, iproc
    use mp, only: comm_split, comm_free
    use sfincs_main, only: init_sfincs, prepare_sfincs, run_sfincs, finish_sfincs
# else
    use mp, only: mp_abort
# endif
    use zgrid, only: nzgrid
    
    implicit none

    integer, intent (in) :: nradii
    real, intent (in) :: drho
    real, dimension (-nzgrid:,:,:,:,-nradii/2:), intent (out) :: f_neoclassical
    real, dimension (-nzgrid:,-nradii/2:), intent (out) :: phi_neoclassical

# ifdef USE_SFINCS
    integer :: sfincs_comm
    integer :: color, ierr
    integer :: irad
    real :: rhoc_neighbor

    if (proc0) call read_sfincs_parameters
    call broadcast_sfincs_parameters
    if (iproc < nproc_sfincs) then
       color = 0
    else
       color = 1
    end if
    call comm_split (color, sfincs_comm, ierr)
    if (iproc < nproc_sfincs) then
       do irad = -nradii/2, nradii/2
          rhoc_neighbor = geo_surf%rhoc + irad*drho
          call init_sfincs (sfincs_comm)
          call pass_inputoptions_to_sfincs (irad*drho)
          call pass_outputoptions_to_sfincs
          call prepare_sfincs
          call pass_geometry_to_sfincs (irad*drho)
          call run_sfincs
          if (proc0) call get_sfincs_output &
               (f_neoclassical(:,:,:,:,irad), phi_neoclassical(:,irad))
!          call broadcast_sfincs_output &
!               (f_neoclassical(:,:,:,:,irad), phi_neoclassical(:,irad))
          call finish_sfincs
       end do
    end if
    call comm_free (sfincs_comm, ierr)

    ! NB: NEED TO CHECK THIS BROADCAST OF SFINCS RESULTS 
    do irad = -nradii/2, nradii/2
       call broadcast_sfincs_output &
            (f_neoclassical(:,:,:,:,irad), phi_neoclassical(:,irad))
    end do
 
# else
    f_neoclassical = 0 ; phi_neoclassical = 0.
    call mp_abort ('to run with include_neoclassical_terms=.true., &
         & USE_SFINCS must be defined at compilation time.  Aborting.')
# endif

  end subroutine get_neo_from_sfincs

# ifdef USE_SFINCS
  subroutine read_sfincs_parameters

    use constants, only: pi
    use mp, only: nproc
    use file_utils, only: input_unit_exist
    use species, only: nspec
    use physics_parameters, only: rhostar, vnew_ref

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
         aHat, psiAHat, nu_N, nxi, nx, Delta, &
         ntheta

    logical :: exist
    integer :: in_file

    nproc_sfincs = 1
    ! do not include radial electric field term
    includeXDotTerm = .false.
    includeElectricFieldTermInXiDot = .false.
    ! no poloidal or toroidal magnetic drifts
    magneticDriftScheme = 0
    ! combo of next two variables means
    ! phi1 will be calculated via quasineutrality
    includePhi1 = .true.
    includePhi1InKineticEquation = .false.
    ! will be overridden by direct input of geometric quantities
    geometryScheme = 1
    ! seems to be a nonsensical option
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
    ! Delta is rho* = mref*vt_ref/(e*Bref*aref), with reference
    ! quantities given in SI units
    Delta = rhostar
    ! nu_n = nu_ref * aref/vt_ref
    ! nu_ref = 4*sqrt(2*pi)*nref*e**4*loglam/(3*sqrt(mref)*Tref**3/2)
    ! (with nref, Tref, and mref in Gaussian units)
    nu_N = vnew_ref*(4./(3.*pi))
    ! number of spectral coefficients in pitch angle
    nxi = 48
    ! number of speeds
    nx = 12
    ! number of poloidal angles
    Ntheta = 65

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
    
    ! ensure that ntheta is odd for SFINCS
    ntheta = 2*(ntheta/2)+1

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
    call broadcast (Delta)
    call broadcast (nu_N)
    call broadcast (nxi)
    call broadcast (nx)
    call broadcast (ntheta)

  end subroutine broadcast_sfincs_parameters

  subroutine pass_inputoptions_to_sfincs (delrho)

    use mp, only: mp_abort
    use geometry, only: geo_surf
    use species, only: spec, nspec
    use zgrid, only: nzed
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
    use globalVariables, only: Nzeta
    use globalVariables, only: nxi_sfincs => Nxi
    use globalVariables, only: nx_sfincs => Nx
    use globalVariables, only: ntheta_sfincs => Ntheta
    use globalVariables, only: dnHatdrNs, dTHatdrNs, dPhiHatdrN
    use globalVariables, only: aHat_sfincs => aHat
    use globalVariables, only: psiAHat_sfincs => psiAHat
    use globalVariables, only: Delta_sfincs => Delta
    use globalVariables, only: nu_n_sfincs => nu_n
    use globalVariables, only: Er_sfincs => Er

    implicit none

    real, intent (in) :: delrho

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
    nHats(:nspec) = spec%dens*(1.0-delrho*spec%fprim)
    THats(:nspec) = spec%temp*(1.0-delrho*spec%tprim)
    mHats(:nspec) = spec%mass
    Zs(:nspec) = spec%z
!     ! FLAG -- need to modify for stellarator simulations
!     ! I think nzeta will be 2*nzgrid+1
!     ! and ntheta will be ny_ffs
    Nzeta = 1
    ntheta_sfincs = ntheta
!    Ntheta = 2*nzgrid+1
    nx_sfincs = nx
    nxi_sfincs = nxi
    aHat_sfincs = aHat
    psiAHat_sfincs = psiAHat
    Delta_sfincs = Delta
    nu_n_sfincs = nu_n

    if (inputRadialCoordinate == 3) then
       rN_wish = geo_surf%rhotor + delrho*geo_surf%drhotordrho
    else
       call mp_abort ('only inputRadialCoordinate=3 currently supported. aborting.')
    end if
    if (inputRadialCoordinateForGradients == 3) then
       ! radial density gradient with respect to rhotor = sqrt(psitor/psitor_LCFS)
       ! normalized by reference density (not species density)
       dnHatdrNs(:nspec) = -spec%dens/geo_surf%drhotordrho*(spec%fprim - delrho*spec%d2ndr2)
       ! radial temperature gradient with respect to rhotor = sqrt(psitor/psitor_LCFS)
       ! normalized by reference tmperatures (not species temperature)
       dTHatdrNs(:nspec) = -spec%temp/geo_surf%drhotordrho*(spec%tprim - delrho*spec%d2Tdr2)
       ! radial electric field
       dPhiHatdrN = 0.0
    else
       call mp_abort ('only inputRadialCoordinateForGradients=3 currently supported. aborting.')
    end if

  end subroutine pass_inputoptions_to_sfincs

  subroutine pass_outputoptions_to_sfincs
    use export_f, only: export_f_theta_option
    use export_f, only: export_f_zeta_option
    use export_f, only: export_f_xi_option
    use export_f, only: export_f_x_option
    use export_f, only: export_delta_f
    implicit none
    export_f_theta_option = 0
    export_f_zeta_option = 0
    export_f_xi_option = 0
    export_f_x_option = 0
    export_delta_f = .true.
  end subroutine pass_outputoptions_to_sfincs

  subroutine pass_geometry_to_sfincs (delrho)

    use constants, only: pi
    use splines, only: linear_interp_periodic
    use zgrid, only: nz2pi, zed
    use geometry, only: bmag, dbdzed, gradpar
    use geometry, only: dBdrho, d2Bdrdth, dgradpardrho, dIdrho
    use geometry, only: geo_surf
    use globalVariables, only: BHat
    use globalVariables, only: dBHatdtheta
    use globalVariables, only: iota
    use globalVariables, only: DHat
    use globalVariables, only: BHat_sup_theta
    use globalVariables, only: BHat_sub_zeta
    use export_f, only: export_f_theta

    implicit none

    real, intent (in) :: delrho

    integer :: nzeta = 1
    integer :: nzpi, iz
    real :: q_local
    real, dimension (:), allocatable :: B_local, dBdz_local, gradpar_local
    real, dimension (:), allocatable :: zed_stella, zed_sfincs

    nzpi = nz2pi/2
    allocate (B_local(-nzpi:nzpi))
    allocate (dBdz_local(-nzpi:nzpi))
    allocate (gradpar_local(-nzpi:nzpi))
    allocate (zed_sfincs(ntheta))
    allocate (zed_stella(-nzpi:nzpi))

    call init_zero_arrays

    ! first get some geometric quantities at this radius
    ! for theta from -pi to pi
    q_local = geo_surf%qinp*(1.0+delrho*geo_surf%shat/geo_surf%rhoc)
    B_local = bmag(1,-nzpi:nzpi) + delrho*dBdrho(-nzpi:nzpi)
    dBdz_local = dbdzed(1,-nzpi:nzpi) + delrho*d2Bdrdth(-nzpi:nzpi)
    gradpar_local = gradpar(1,-nzpi:nzpi) + delrho*dgradpardrho(-nzpi:nzpi)

    zed_stella = zed(-nzpi:nzpi)+pi
    zed_sfincs = export_f_theta(:ntheta)

    iota = 1./q_local

    ! interpolate from stella zed-grid to sfincs theta grid
    ! point at -pi (stella) is same as point at 0 (sfincs)
    BHat(1,1) = B_local(-nzpi)
    call linear_interp_periodic (zed_stella, B_local, zed_sfincs(2:), BHat(2:,1))
    ! FLAG -- needs to be changed for stellarator runs
    BHat = spread(BHat(:,1),2,nzeta)
    
    dBHatdtheta(1,1) = dBdz_local(-nzpi)
    call linear_interp_periodic (zed_stella, dBdz_local, zed_sfincs(2:), dBHatdtheta(2:,1))
    dBHatdtheta = spread(dBHatdtheta(:,1),2,nzeta)

    ! this is bhat . grad theta
    BHat_sup_theta(1,1) = B_local(-nzpi)*gradpar_local(-nzpi)
    call linear_interp_periodic (zed_stella, B_local*gradpar_local, zed_sfincs(2:), BHat_sup_theta(2:,1))
    BHat_sup_theta = spread(BHat_sup_theta(:,1),2,nzeta)
    ! this is I(psi) / (aref*Bref)
    BHat_sub_zeta = geo_surf%rgeo + delrho*dIdrho
    ! this is grad psitor . (grad theta x grad zeta)
    ! note that + sign below relies on B = I grad zeta + grad zeta x grad psi
    DHat = q_local*BHat_sup_theta

    deallocate (B_local, dBdz_local, gradpar_local)
    deallocate (zed_sfincs, zed_stella)

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

  subroutine get_sfincs_output (f_neoclassical, phi_neoclassical)
    
    use constants, only: pi
    use splines, only: linear_interp_periodic
    use mp, only: mp_abort
    use species, only: nspec, spec
    use zgrid, only: nzgrid, nz2pi, nperiod, zed
    use vpamu_grids, only: nvgrid, nmu
    use vpamu_grids, only: vpa, mu, ztmax, vperp2, maxwell_mu
    use export_f, only: h_sfincs => delta_f
    use export_f, only: zed_sfincs => export_f_theta
    use globalVariables, only: nxi_sfincs => nxi
    use globalVariables, only: nx_sfincs => nx
    use globalVariables, only: x_sfincs => x
    use globalVariables, only: phi_sfincs => Phi1Hat
    use xGrid, only: xGrid_k
    use geometry, only: bmag

    implicit none

    real, dimension (-nzgrid:,:,:,:), intent (out) :: f_neoclassical
    real, dimension (-nzgrid:), intent (out) :: phi_neoclassical

    integer :: iz, iv, imu, is, ixi, ip, i, j

    integer :: nzpi, iz_low, iz_up
    integer :: nxi_stella
    real, dimension (1) :: x_stella
    integer, dimension (2) :: sgnvpa
    real, dimension (:), allocatable :: zed_stella, phi_stella
    real, dimension (:), allocatable :: xi_stella, hstella
    real, dimension (:), allocatable :: htmp, dhtmp_dx
    real, dimension (:), allocatable :: hdum
    real, dimension (:,:,:), allocatable :: h_stella
    real, dimension (:,:), allocatable :: xsfincs_to_xstella, legpoly
    real, dimension (:,:), allocatable :: hsfincs

    allocate (htmp(nxi_sfincs))
    allocate (dhtmp_dx(nxi_sfincs))
    allocate (hsfincs(nxi_sfincs,nx_sfincs))
    allocate (xsfincs_to_xstella(1,nx_sfincs))

    allocate (zed_stella(nz2pi))
    allocate (phi_stella(nz2pi))

    nzpi = nz2pi/2

    zed_stella = zed(-nzpi:nzpi) + pi

    phi_stella(1) = phi_sfincs(1,1)
    phi_stella(nz2pi) = phi_stella(1)
    call linear_interp_periodic (zed_sfincs(:ntheta), phi_sfincs(:ntheta,1), zed_stella(2:nz2pi-1), phi_stella(2:nz2pi-1))

    iz_low = -nzgrid
    iz_up = -nzgrid+nz2pi-1
    phi_neoclassical(iz_low:iz_up) = phi_stella
    ! if nperiod > 1 need to make copies of
    ! neoclassical potential for other 2pi segments
    if (nperiod > 1) then
       do ip = 2, 2*nperiod-1
          iz_low = iz_up + 1
          iz_up = iz_low + nz2pi -2
          phi_neoclassical(iz_low:iz_up) = phi_stella(2:)
       end do
    end if

    deallocate (phi_stella)
    allocate (h_stella(-nzgrid:nzgrid,size(h_sfincs,4),size(h_sfincs,5)))
    allocate (hdum(nz2pi))

    sgnvpa(1) = 1 ; sgnvpa(2) = -1
    do is = 1, nspec
       do i = 1, size(h_sfincs,5)
          do j = 1, size(h_sfincs,4)
             hdum(1) = h_sfincs(is,1,1,j,i)
             hdum(nz2pi) = hdum(1)
             call linear_interp_periodic (zed_sfincs(:ntheta), h_sfincs(is,:ntheta,1,j,i), &
                  zed_stella(2:nz2pi-1), hdum(2:nz2pi-1))
             iz_low = -nzgrid
             iz_up = -nzgrid+nz2pi-1
             h_stella(iz_low:iz_up,j,i) = hdum
             if (nperiod > 1) then
                do ip = 2, 2*nperiod-1
                   iz_low = iz_up + 1
                   iz_up = iz_low + nz2pi -2
                   h_stella(iz_low:iz_up,j,i) = hdum(2:)
                end do
             end if
          end do
       end do
       do iz = -nzgrid, nzgrid
          ! hsfincs is on the sfincs energy grid
          ! but is spectral in pitch-angle
          hsfincs = h_stella(iz,:,:)

          do imu = 1, nmu
             do iv = 1, nvgrid
                ! x_stella is the speed 
                ! corresponding to this (vpa,mu) grid point
                ! FLAG -- NEED TO EXTEND SFINCS TREATMENT TO INCLUDE MULTIPLE ALPHAS
                x_stella = sqrt(vpa(iv)**2+vperp2(1,iz,imu))
                ! note that with exception of vpa=0
                ! can use symmetry of vpa grid to see that
                ! each speed arc has two pitch angles on it
                ! correspondong to +/- vpa
                nxi_stella = 2
                allocate (xi_stella(nxi_stella))
                allocate (hstella(nxi_stella))
                allocate (legpoly(nxi_stella,0:nxi_sfincs-1))
                ! xi_stella is the pitch angle (vpa/v)
                ! corresponding to this (vpa,mu) grid point
                xi_stella = sgnvpa(:nxi_stella)*vpa(iv)/x_stella(1)

                ! set up matrix that interpolates from sfincs speed grid
                ! to the speed corresponding to this (vpa,mu) grid point
                call polynomialInterpolationMatrix (nx_sfincs, 1, &
                     x_sfincs, x_stella, exp(-x_sfincs*x_sfincs)*(x_sfincs**xGrid_k), &
                     exp(-x_stella*x_stella)*(x_stella**xGrid_k), xsfincs_to_xstella)
                
                ! do the interpolation
                do ixi = 1, nxi_sfincs
                   htmp(ixi) = sum(hsfincs(ixi,:)*xsfincs_to_xstella(1,:))
                end do

                ! next need to Legendre transform in pitch-angle
                ! first evaluate Legendre polynomials at requested pitch angles
                call legendre (xi_stella, legpoly)

                ! then do the transforms
                call legendre_transform (legpoly, htmp, hstella)

                f_neoclassical(iz,nvgrid-iv+1,imu,is) = hstella(2)
                f_neoclassical(iz,iv+nvgrid,imu,is) = hstella(1)

                deallocate (xi_stella, hstella, legpoly)
             end do
             ! h_sfincs is H_nc / (nref/vt_ref^3), with H_nc the non-Boltzmann part of F_nc
             ! NB: n_ref, etc. is fixed in stella to be the reference density
             ! at the central sfincs simulation; i.e., it does not vary with radius
             ! similarly, bmag below is the normalized B-field at the central radial location
             ! to be consistent with stella distribution functions,
             ! want H_nc / (n_s / vt_s^3 * pi^(3/2))
             f_neoclassical(iz,:,imu,is) = f_neoclassical(iz,:,imu,is) &
                  * pi**1.5 * spec(is)%stm**3/spec(is)%dens
             
             ! phi_sfincs is e phi / Tref as long as alpha=1 (default)
             ! need to multiply by Z_s * Tref/T_s * exp(-v^2)
             f_neoclassical(iz,:,imu,is) = f_neoclassical(iz,:,imu,is) &
                  - phi_neoclassical(iz)*ztmax(:,is)*maxwell_mu(1,iz,imu)
          end do
       end do
    end do

    deallocate (hdum, h_stella)
    deallocate (zed_stella)
    deallocate (htmp, dhtmp_dx)
    deallocate (hsfincs)
    deallocate (xsfincs_to_xstella)

  end subroutine get_sfincs_output

  ! returns the Legendre polynomials (legp)
  ! on requested grid (x)
  subroutine legendre (x, legp)
    
    implicit none
    
    real, dimension (:), intent (in) :: x
    real, dimension (:,0:), intent (out) :: legp
    
    integer :: n, idx
    
    n = size(legp,2)-1
    
    legp(:,0) = 1.0
    legp(:,1) = x
    
    do idx = 2, n
       legp(:,idx) = ((2.*idx-1.)*x*legp(:,idx-1) + (1.-idx)*legp(:,idx-2))/idx
    end do
    
  end subroutine legendre
  
  subroutine legendre_transform (legp, coefs, func)

    implicit none
    
    real, dimension (:,0:), intent (in) :: legp
    real, dimension (:), intent (in) :: coefs
    real, dimension (:), intent (out) :: func
    
    integer :: i

    func = 0.
    do i = 1, size(coefs)
       func = func + legp(:,i-1)*coefs(i)
    end do

  end subroutine legendre_transform

  subroutine broadcast_sfincs_output (fneo, phineo)
    use mp, only: broadcast
    use zgrid, only: nzgrid
    implicit none
    real, dimension (-nzgrid:,:,:,:), intent (in out) :: fneo
    real, dimension (-nzgrid:), intent (in out) :: phineo
    call broadcast (fneo)
    call broadcast (phineo)
  end subroutine broadcast_sfincs_output

# endif

end module sfincs_interface
