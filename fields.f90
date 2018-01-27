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
    use fields_arrays, only: phi_old, apar_old
    use dist_fn_arrays, only: gvmu
    use stella_layouts, only: init_stella_layouts
    use species, only: init_species
    use geometry, only: init_geometry
    use zgrid, only: init_zgrid
!    use kt_grids, only: init_kt_grids
    use run_parameters, only: init_run_parameters
    use physics_parameters, only: init_physics_parameters
    use dist_fn, only: init_dist_fn
    use dist_fn, only: init_get_fields, get_fields
    use dist_fn, only: stream_implicit
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
    call init_geometry
    if (debug) write(6,*) "fields::init_fields::init_physics_parameters"
    call init_physics_parameters
    if (debug) write (*,*) 'fields::init_fields::init_species'
    call init_species
    if (debug) write(6,*) "fields::init_fields::init_init_g"
    call init_init_g
    if (debug) write(6,*) "fields::init_fields::init_run_parameters"
    call init_run_parameters
    if (debug) write(6,*) "fields::init_fields::init_dist_fn"
    call init_dist_fn
    if (debug) write(6,*) "fields::init_fields::allocate_arrays"
    call allocate_arrays
    if (debug) write(*,*) "fields::init_fields::ginit"
    call ginit (restarted)
    if (debug) write(*,*) "fields::init_fields::init_gxyz"
    call init_gxyz
    ! initialize get_fields subroutine
    if (debug) write (*,*) 'fields::init_fields::init_get_fields'
    call init_get_fields
    if (debug) write(*,*) "fields::init_fields::init_response_matrix"
    if (stream_implicit) call init_response_matrix

    if (restarted) return

    if (debug) write (*,*) 'fields::init_fields::get_fields'
    ! get initial field from initial distribution function
    call get_fields (gvmu, phi, apar, dist='gbar')
    phi_old = phi ; apar_old = apar

  end subroutine init_fields

  subroutine init_response_matrix
    
    use linear_solve, only: lu_decomposition
    use fields_arrays, only: response_matrix
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, is_idx
    use species, only: nspec
    use kt_grids, only: naky, nakx
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: neigen, ikxmod
    use extended_zgrid, only: nsegments, nsegments_poskx
    use extended_zgrid, only: nzed_segment

    implicit none

    integer :: iky, ie, iseg, iz
    integer :: ikx
    integer :: nresponse, nz_ext
    integer :: llim, ulim
    integer :: idx
    integer :: izl_offset
    real :: dum
    complex, dimension (:), allocatable :: phiext
    complex, dimension (:,:), allocatable :: gext

    ! for a given ky and set of connected kx values
    ! give a unit impulse to phi at each zed location
    ! in the extended domain and solve for h(zed_extended,(vpa,mu,s))

    do iky = 1, naky

       ! the response matrix for each ky has neigen(ky)
       ! independent sets of connected kx values
       if (.not.associated(response_matrix(iky)%eigen)) &
            allocate (response_matrix(iky)%eigen(neigen(iky)))

       ! loop over the sets of connected kx values
       do ie = 1, neigen(iky)

          ! number of zeds x number of segments
          nz_ext = nsegments(ie,iky)*nzed_segment+1
          ! number of zeds x number of segments with kx >= 0
          ! note that we do not need to consider response to
          ! perturbations with kx < 0 as these are not explicitly evolved
          nresponse = nsegments_poskx(ie,iky)*nzed_segment+1

          ! llim and ulim indicate which chunk of the larger array
          ! of size nz_ext contains the information needed for the response matrix
          ! i.e., the non-negative kx values

          ! if the kx for the first segment is negative
          ! then the non-negative kx segments will be 
          ! at the end of the array
          if (ikxmod(1,ie,iky) > nakx) then
             llim = nz_ext-nresponse+1
             ulim = nz_ext
          ! otherwise, non-negative kx segments will be
          ! at the beginning
          else
             llim = 1
             ulim = nresponse
          end if

          ! for each ky and set of connected kx values,
          ! must have a response matrix that is N x N
          ! with N = number of zeds x number of 2pi segments with kx >= 0
          if (.not.associated(response_matrix(iky)%eigen(ie)%zloc)) &
               allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse,nresponse))

          ! response_matrix%idx is needed to keep track of permutations
          ! to the response matrix made during LU decomposition
          ! it will be input to LU back substitution during linear solve
          if (.not.associated(response_matrix(iky)%eigen(ie)%idx)) &
               allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))

          allocate (gext(nz_ext,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
          allocate (phiext(nresponse))
          ! idx is the index in the extended zed domain 
          ! that we are giving a unit impulse
          ! the offset with llim is to take into account
          ! the fact that we only want to perturb non-negative kx segments
          idx = llim-1

          ! loop over segments, starting with 1
          ! first segment is special because it has 
          ! one more unique zed value than all others
          ! since domain is [z0-pi:z0+pi], including both endpoints
          ! i.e., one endpoint is shared with the previous segment
          iseg = 1
          ! ikxmod gives the kx corresponding to iseg,ie,iky
          ikx = ikxmod(iseg,ie,iky)
          izl_offset = 0
          ! no need to obtain response to impulses at negative kx values
          if (ikx <= nakx) then
             do iz = iz_low(iseg), iz_up(iseg)
                idx = idx + 1
                call get_response_matrix_column (iky, ikx, iz, ie, idx, llim, ulim, phiext, gext)
             end do
             ! once we have used one segments, remaining segments
             ! have one fewer unique zed point
             izl_offset = 1
          end if
          if (nsegments(ie,iky) > 1) then
             do iseg = 2, nsegments(ie,iky)
                ikx = ikxmod(iseg,ie,iky)
                ! no need to treat negative kx values
                if (ikx > nakx) cycle
                do iz = iz_low(iseg)+izl_offset, iz_up(iseg)
                   idx = idx + 1
                   call get_response_matrix_column (iky, ikx, iz, ie, idx, llim, ulim, phiext, gext)
                end do
                if (izl_offset == 0) izl_offset = 1
             end do
          end if

          if (any(ikxmod(:nsegments(ie,iky),ie,iky) <= nakx)) then
             ! now that we have the reponse matrix for this ky and set of connected kx values
             ! get the LU decomposition so we are ready to solve the linear system
             call lu_decomposition (response_matrix(iky)%eigen(ie)%zloc,response_matrix(iky)%eigen(ie)%idx,dum)
          end if

          deallocate (gext, phiext)
       end do
    end do

  end subroutine init_response_matrix

  subroutine get_response_matrix_column (iky, ikx, iz, ie, idx, llim, ulim, phiext, gext)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, is_idx
    use stella_time, only: code_dt
    use zgrid, only: delzed
    use species, only: spec
    use geometry, only: gradpar
    use vpamu_grids, only: ztmax, vpa
    use fields_arrays, only: response_matrix
    use dist_fn_arrays, only: aj0x
    use dist_fn, only: stream_tridiagonal_solve
    use dist_fn, only: zed_upwind, time_upwind

    implicit none

    integer, intent (in) :: iky, ikx, iz, ie, idx, llim, ulim
    complex, dimension (:), intent (in out) :: phiext
    complex, dimension (:,vmu_lo%llim_proc:), intent (in out) :: gext

    integer :: ivmu, iv, is, idxp
    integer :: nz_ext
    integer :: sgn
    real :: fac

    nz_ext = size(gext,1)

    ! get -vpa*b.gradz*Ze/T*F0*d<phi>/dz corresponding to unit impulse in phi
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       ! initialize g to zero everywhere along extended zed domain
       gext(:,ivmu) = 0.0
       iv = iv_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)

       ! give unit impulse to phi at this zed location
       ! and compute -vpa*b.gradz*Ze/T*d<phi>/dz*F0 (RHS of streaming part of GKE)

       ! note -- assuming equal spacing in zed below
       fac = -0.5*(1.+time_upwind)*code_dt*vpa(iv)*gradpar(1,iz)*spec(is)%stm &
            *aj0x(iky,ikx,iz,ivmu)*ztmax(iv,is)/delzed(0)
       ! e.g., if centered differences, get RHS(i) = fac*(phi(i+1)-phi(i-1))*0.5
       ! so RHS(i-1) = fac*(phi(i)-phi(i-2))*0.5
       ! since only gave phi(i)=1 and all other phi=0,
       ! RHS(i-1) = fac*0.5
       ! similarly, RHS(i+1) = fac*(phi(i+2)-phi(i))*0.5
       ! = -0.5*fac*phi(i)
       ! if upwinded and vpa > 0, RHS(i) = fac*(phi(i)-phi(i-1))
       ! and so RHS(i) = fac*phi(i) and RHS(i+1) = -fac*phi(i)

       ! test for sign of vpa (vpa(iv<0) < 0, vpa(iv>0) > 0),
       ! as this affects upwinding of d<phi>/dz source term
       if (iv < 0) then
          sgn = -1
       else if (iv > 0) then
          sgn = 1
       else
          sgn = 0
       end if

       gext(idx,ivmu) = sgn*zed_upwind*fac
       if (idx > 1) gext(idx-1,ivmu) = 0.5*(1.-sgn*zed_upwind)*fac
       if (idx < nz_ext) gext(idx+1,ivmu) = -0.5*(1.+sgn*zed_upwind)*fac

       ! invert parallel streaming equation to get g^{n+1} on extended zed grid
       ! (I + (1+alph)/2*dt*vpa)*g_{inh}^{n+1} = RHS = gext
       call stream_tridiagonal_solve (iky, ie, iv, is, gext(:,ivmu))
    end do

    ! we now have g on the extended zed domain at this ky and set of connected kx values
    ! corresponding to a unit impulse in phi at this location
    ! obtain the fields associated with this g, but only need it for kx >= 0
    call get_fields_for_response_matrix (gext(llim:ulim,:), phiext, iky, ie)

    ! next need to create column in response matrix from phiext
    ! negative sign because matrix to be inverted in streaming equation
    ! is identity matrix - response matrix
    ! add in contribution from identity matrix
    idxp = idx-llim+1
    phiext(idxp) = phiext(idxp)-1.0
    response_matrix(iky)%eigen(ie)%zloc(:,idxp) = -phiext

  end subroutine get_response_matrix_column

  subroutine get_fields_for_response_matrix (g, phi, iky, ie)

    use stella_layouts, only: vmu_lo
    use species, only: nspec, spec
    use species, only: has_electron_species
    use geometry, only: dl_over_b
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: nsegments
    use kt_grids, only: aky
    use kt_grids, only: nakx
    use vpamu_grids, only: integrate_species
    use dist_fn_arrays, only: aj0x
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg
    use dist_fn, only: gamtot, gamtot3

    implicit none

    complex, dimension (:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:), intent (out) :: phi
    integer, intent (in) :: iky, ie
    
    integer :: idx, iseg, ikx, iz
    integer :: izl_offset
    real, dimension (nspec) :: wgt
    complex, dimension (:), allocatable :: g0
    complex :: tmp

    allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    wgt = spec%z*spec%dens
    phi = 0.

    ! FLAG -- IS THIS CORRECT FOR TWIST AND SHIFT
    ! BECAUSE ONLY KX>=0 PASSED IN WHICH MAY MEAN 
    ! ISEG = 1, ETC. ARE NOT INCLUDED

    idx = 0 ; izl_offset = 0
    iseg = 1
    ikx = ikxmod(iseg,ie,iky)
    if (ikx <= nakx) then
       do iz = iz_low(iseg), iz_up(iseg)
          idx = idx + 1
          g0 = aj0x(iky,ikx,iz,:)*g(idx,:)
          call integrate_species (g0, iz, wgt, phi(idx))
          phi(idx) = phi(idx)/gamtot(iky,ikx,iz)
       end do
       izl_offset = 1
    end if
    if (nsegments(ie,iky) > 1) then
       do iseg = 2, nsegments(ie,iky)
          ikx = ikxmod(iseg,ie,iky)
          if (ikx > nakx) cycle
          do iz = iz_low(iseg)+izl_offset, iz_up(iseg)
             idx = idx + 1
             g0 = aj0x(iky,ikx,iz,:)*g(idx,:)
             call integrate_species (g0, iz, wgt, phi(idx))
             phi(idx) = phi(idx)/gamtot(iky,ikx,iz)
          end do
          if (izl_offset == 0) izl_offset = 1
       end do
    end if

    if (.not.has_electron_species(spec) .and. &
         adiabatic_option_switch == adiabatic_option_fieldlineavg) then
       if (abs(aky(iky)) < epsilon(0.)) then
          ! no connections for ky = 0
          iseg = 1 
          tmp = sum(dl_over_b*phi)
          phi = phi + tmp*gamtot3(ikxmod(1,ie,iky),:)
       end if
    end if

    deallocate (g0)

  end subroutine get_fields_for_response_matrix

  subroutine allocate_arrays

    use fields_arrays, only: phi, apar
    use fields_arrays, only: phi_old, apar_old
    use fields_arrays, only: response_matrix
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx
    use dist_fn, only: stream_implicit

    implicit none

    if (.not.allocated(phi)) then
       allocate (phi(naky,nakx,-nzgrid:nzgrid))
       phi = 0.
    end if
    if (.not. allocated(apar)) then
       allocate (apar(naky,nakx,-nzgrid:nzgrid))
       apar = 0.
    end if
    if (.not.allocated(phi_old)) then
       allocate (phi_old(naky,nakx,-nzgrid:nzgrid))
       phi_old = 0.
    end if
    if (.not. allocated(apar_old)) then
       allocate (apar_old(naky,nakx,-nzgrid:nzgrid))
       apar_old = 0.
    end if
    if (.not.allocated(response_matrix)) then
       if (stream_implicit) then
          allocate (response_matrix(naky))
       else
          allocate (response_matrix(1))
       end if
!       response_matrix = 0.
    end if

  end subroutine allocate_arrays

  subroutine finish_response_matrix

    use fields_arrays, only: response_matrix

    implicit none

    if (allocated(response_matrix)) deallocate (response_matrix)

  end subroutine finish_response_matrix

  subroutine finish_fields

    use fields_arrays, only: phi, phi_old
    use fields_arrays, only: apar, apar_old
    use species, only: finish_species
    use geometry, only: finish_geometry
    use zgrid, only: finish_zgrid
    use dist_fn, only: finish_get_fields

    implicit none

    call finish_response_matrix
    call finish_get_fields
    call finish_geometry
    call finish_zgrid
    call finish_species
    if (allocated(phi)) deallocate (phi)
    if (allocated(phi_old)) deallocate (phi_old)
    if (allocated(apar)) deallocate (apar)
    if (allocated(apar_old)) deallocate (apar_old)

    fields_initialized = .false.

  end subroutine finish_fields

end module fields
