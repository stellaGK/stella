module response_matrix

  implicit none

  public :: init_response_matrix, finish_response_matrix
  public :: response_matrix_initialized

  private

  logical :: response_matrix_initialized = .false.

contains

  subroutine init_response_matrix
    
    use linear_solve, only: lu_decomposition
    use fields_arrays, only: response_matrix
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, is_idx
    use species, only: nspec
    use kt_grids, only: naky, zonal_mode
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: neigen, ikxmod
    use extended_zgrid, only: nsegments
    use extended_zgrid, only: nzed_segment

    implicit none

    integer :: iky, ie, iseg, iz
    integer :: ikx
    integer :: nz_ext, nresponse
    integer :: idx
    integer :: izl_offset, izup
    real :: dum
    complex, dimension (:), allocatable :: phiext
    complex, dimension (:,:), allocatable :: gext

    if (response_matrix_initialized) return
    response_matrix_initialized = .true.
    
    if (.not.allocated(response_matrix)) allocate (response_matrix(naky))
    
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
          
          ! treat zonal mode specially to avoid double counting
          ! as it is periodic
          if (zonal_mode(iky)) then
             nresponse = nz_ext-1
          else
             nresponse = nz_ext
          end if

          ! for each ky and set of connected kx values,
          ! must have a response matrix that is N x N
          ! with N = number of zeds per 2pi segment x number of 2pi segments
          if (.not.associated(response_matrix(iky)%eigen(ie)%zloc)) &
               allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse,nresponse))

          ! response_matrix%idx is needed to keep track of permutations
          ! to the response matrix made during LU decomposition
          ! it will be input to LU back substitution during linear solve
          if (.not.associated(response_matrix(iky)%eigen(ie)%idx)) &
               allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))

          allocate (gext(nz_ext,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
          allocate (phiext(nz_ext))
          ! idx is the index in the extended zed domain
          ! that we are giving a unit impulse
          idx = 0

          ! loop over segments, starting with 1
          ! first segment is special because it has 
          ! one more unique zed value than all others
          ! since domain is [z0-pi:z0+pi], including both endpoints
          ! i.e., one endpoint is shared with the previous segment
          iseg = 1
          ! ikxmod gives the kx corresponding to iseg,ie,iky
          ikx = ikxmod(iseg,ie,iky)
          izl_offset = 0
          ! avoid double-counting of periodic points for zonal mode
          if (zonal_mode(iky)) then
             izup = iz_up(iseg)-1
          else
             izup = iz_up(iseg)
          end if
          ! no need to obtain response to impulses at negative kx values
          do iz = iz_low(iseg), izup
             idx = idx + 1
             call get_response_matrix_column (iky, ikx, iz, ie, idx, nz_ext, nresponse, phiext, gext)
          end do
          ! once we have used one segments, remaining segments
          ! have one fewer unique zed point
          izl_offset = 1
          if (nsegments(ie,iky) > 1) then
             do iseg = 2, nsegments(ie,iky)
                ikx = ikxmod(iseg,ie,iky)
                do iz = iz_low(iseg)+izl_offset, iz_up(iseg)
                   idx = idx + 1
                   call get_response_matrix_column (iky, ikx, iz, ie, idx, nz_ext, nresponse, phiext, gext)
                end do
                if (izl_offset == 0) izl_offset = 1
             end do
          end if

          ! now that we have the reponse matrix for this ky and set of connected kx values
          ! get the LU decomposition so we are ready to solve the linear system
          call lu_decomposition (response_matrix(iky)%eigen(ie)%zloc,response_matrix(iky)%eigen(ie)%idx,dum)

          deallocate (gext, phiext)
       end do
    end do

  end subroutine init_response_matrix

  subroutine get_response_matrix_column (iky, ikx, iz, ie, idx, nz_ext, nresponse, phiext, gext)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use zgrid, only: delzed, nzgrid
    use kt_grids, only: zonal_mode
    use species, only: spec
    use stella_geometry, only: gradpar, dbdzed
    use vpamu_grids, only: vpa, mu
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use fields_arrays, only: response_matrix
    use gyro_averages, only: aj0x
    use run_parameters, only: driftkinetic_implicit
    use run_parameters, only: maxwellian_inside_zed_derivative
    use parallel_streaming, only: stream_tridiagonal_solve
    use parallel_streaming, only: stream_sign
    use run_parameters, only: zed_upwind, time_upwind

    implicit none

    integer, intent (in) :: iky, ikx, iz, ie, idx, nz_ext, nresponse
    complex, dimension (:), intent (in out) :: phiext
    complex, dimension (:,vmu_lo%llim_proc:), intent (in out) :: gext

    integer :: ivmu, iv, imu, is, ia
    integer :: izp, izm
    real :: mu_dbdzed_p, mu_dbdzed_m
    real :: fac, fac0, fac1, gyro_fac
    real, dimension (:), allocatable :: gradpar_fac

    ia = 1

    if (.not.allocated(gradpar_fac)) allocate (gradpar_fac(-nzgrid:nzgrid))
    gradpar_fac = gradpar

    if (.not.maxwellian_inside_zed_derivative) then
       ! get -vpa*b.gradz*Ze/T*F0*d<phi>/dz corresponding to unit impulse in phi
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          ! initialize g to zero everywhere along extended zed domain
          gext(:,ivmu) = 0.0
          iv = iv_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          is = is_idx(vmu_lo,ivmu)
          
          ! give unit impulse to phi at this zed location
          ! and compute -vpa*b.gradz*Ze/T*d<phi>/dz*F0 (RHS of streaming part of GKE)
          
          ! NB:  assuming equal spacing in zed below
          ! here, fac = -dt*(1+alph_t)/2*vpa*Ze/T*F0*J0/dz
          ! b.gradz left out because needs to be centred in zed
          if (driftkinetic_implicit) then
             gyro_fac = 1.0
          else
             gyro_fac = aj0x(iky,ikx,iz,ivmu)
          end if
          
          fac = -0.25*(1.+time_upwind)*code_dt*vpa(iv)*spec(is)%stm_psi0 &
               *gyro_fac*spec(is)%zt/delzed(0)*maxwell_vpa(iv,is)
          
          gradpar_fac = gradpar*maxwell_mu(ia,:,imu,is)*maxwell_fac(is)
          
          ! stream_sign < 0 corresponds to positive advection speed
          if (stream_sign(iv)<0) then
             if (iz > -nzgrid) then
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar_fac(iz) &
                     + (1.-zed_upwind)*gradpar_fac(iz-1))
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the right of
                ! this one
                if (iz < nzgrid) then
                   fac1 = fac*((1.+zed_upwind)*gradpar_fac(iz+1) &
                        + (1.-zed_upwind)*gradpar_fac(iz))
                else
                   fac1 = fac*((1.+zed_upwind)*gradpar_fac(-nzgrid+1) &
                        + (1.-zed_upwind)*gradpar_fac(nzgrid))
                end if
             else
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar_fac(iz) &
                     + (1.-zed_upwind)*gradpar_fac(nzgrid-1))
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the right of
                ! this one
                fac1 = fac*((1.+zed_upwind)*gradpar_fac(iz+1) &
                     + (1.-zed_upwind)*gradpar_fac(iz))
             end if
             gext(idx,ivmu) = fac0
             if (idx < nz_ext) gext(idx+1,ivmu) = -fac1
             ! zonal mode BC is periodic instead of zero, so must
             ! treat specially
             if (zonal_mode(iky)) then
                if (idx == 1) then
                   gext(nz_ext,ivmu) = fac0
                else if (idx == nz_ext-1) then
                   gext(1,ivmu) = -fac1
                end if
             end if
          else
             if (iz < nzgrid) then
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar_fac(iz) &
                     + (1.-zed_upwind)*gradpar_fac(iz+1))
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the left of
                ! this one
                if (iz > -nzgrid) then
                   fac1 = fac*((1.+zed_upwind)*gradpar_fac(iz-1) &
                        + (1.-zed_upwind)*gradpar_fac(iz))
                else
                   fac1 = fac*((1.+zed_upwind)*gradpar_fac(nzgrid-1) &
                        + (1.-zed_upwind)*gradpar_fac(iz))
                end if
             else
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar_fac(iz) &
                     + (1.-zed_upwind)*gradpar_fac(-nzgrid+1))
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the left of
                ! this one
                fac1 = fac*((1.+zed_upwind)*gradpar_fac(iz-1) &
                     + (1.-zed_upwind)*gradpar_fac(iz))
             end if
             gext(idx,ivmu) = -fac0
             if (idx > 1) gext(idx-1,ivmu) = fac1
             ! zonal mode BC is periodic instead of zero, so must
             ! treat specially
             if (zonal_mode(iky)) then
                if (idx == 1) then
                   gext(nz_ext,ivmu) = -fac0
                   gext(nz_ext-1,ivmu) = fac1
                else if (idx == 2) then
                   gext(nz_ext,ivmu) = fac1
                end if
             end if
          end if

          ! hack for now (duplicates much of the effort from sweep_zed_zonal)
          if (zonal_mode(iky)) then
             call sweep_zed_zonal_response (iv, is, stream_sign(iv), gext(:,ivmu))
          else
             ! invert parallel streaming equation to get g^{n+1} on extended zed grid
             ! (I + (1+alph)/2*dt*vpa)*g_{inh}^{n+1} = RHS = gext
             call stream_tridiagonal_solve (iky, ie, iv, is, gext(:,ivmu))
          end if

       end do
    else
       ! get -vpa*b.gradz*Ze/T*F0*d<phi>/dz corresponding to unit impulse in phi
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          ! initialize g to zero everywhere along extended zed domain
          gext(:,ivmu) = 0.0
          iv = iv_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          is = is_idx(vmu_lo,ivmu)
          
          ! give unit impulse to phi at this zed location
          ! and compute -vpa*b.gradz*Ze/T*d<phi>/dz*F0 (RHS of streaming part of GKE)
          
          ! NB:  assuming equal spacing in zed below
          ! here, fac = -dt*(1+alph_t)/2*vpa*Ze/T*F0*J0/dz
          ! b.gradz left out because needs to be centred in zed
          if (driftkinetic_implicit) then
             gyro_fac = 1.0
          else
             gyro_fac = aj0x(iky,ikx,iz,ivmu)
          end if
          
          fac = -0.25*(1.+time_upwind)*code_dt*vpa(iv)*spec(is)%stm_psi0 &
               *gyro_fac*spec(is)%zt*maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is)

          mu_dbdzed_p = 1./delzed(0)+mu(imu)*dbdzed(ia,iz)*(1.+zed_upwind)
          mu_dbdzed_m = 1./delzed(0)+mu(imu)*dbdzed(ia,iz)*(1.-zed_upwind)
          
          ! stream_sign < 0 corresponds to positive advection speed
          if (stream_sign(iv)<0) then
             if (iz > -nzgrid) then
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar(iz) &
                     + (1.-zed_upwind)*gradpar(iz-1))*mu_dbdzed_p
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the right of
                ! this one
                if (iz < nzgrid) then
                   izp = iz+1
                else
                   izp = -nzgrid+1
                end if
                fac1 = fac*((1.+zed_upwind)*gradpar(izp) &
                     + (1.-zed_upwind)*gradpar(iz))*mu_dbdzed_m
             else
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar(iz) &
                     + (1.-zed_upwind)*gradpar(nzgrid-1))*mu_dbdzed_p
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the right of
                ! this one
                fac1 = fac*((1.+zed_upwind)*gradpar(iz+1) &
                     + (1.-zed_upwind)*gradpar(iz))*mu_dbdzed_m
             end if
             gext(idx,ivmu) = fac0
             if (idx < nz_ext) gext(idx+1,ivmu) = -fac1
             ! zonal mode BC is periodic instead of zero, so must
             ! treat specially
             if (zonal_mode(iky)) then
                if (idx == 1) then
                   gext(nz_ext,ivmu) = fac0
                else if (idx == nz_ext-1) then
                   gext(1,ivmu) = -fac1
                end if
             end if
          else
             if (iz < nzgrid) then
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar(iz) &
                     + (1.-zed_upwind)*gradpar(iz+1))*mu_dbdzed_p
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the left of
                ! this one
                if (iz > -nzgrid) then
                   izm = iz-1
                else
                   izm = nzgrid-1
                end if
                fac1 = fac*((1.+zed_upwind)*gradpar(izm) &
                     + (1.-zed_upwind)*gradpar(iz))*mu_dbdzed_m
             else
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar(iz) &
                     + (1.-zed_upwind)*gradpar(-nzgrid+1))*mu_dbdzed_p
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the left of
                ! this one
                fac1 = fac*((1.+zed_upwind)*gradpar(iz-1) &
                     + (1.-zed_upwind)*gradpar(iz))*mu_dbdzed_m
             end if
             gext(idx,ivmu) = -fac0
             if (idx > 1) gext(idx-1,ivmu) = fac1
             ! zonal mode BC is periodic instead of zero, so must
             ! treat specially
             if (zonal_mode(iky)) then
                if (idx == 1) then
                   gext(nz_ext,ivmu) = -fac0
                   gext(nz_ext-1,ivmu) = fac1
                else if (idx == 2) then
                   gext(nz_ext,ivmu) = fac1
                end if
             end if
          end if

          ! hack for now (duplicates much of the effort from sweep_zed_zonal)
          if (zonal_mode(iky)) then
             call sweep_zed_zonal_response (iv, is, stream_sign(iv), gext(:,ivmu))
          else
             ! invert parallel streaming equation to get g^{n+1} on extended zed grid
             ! (I + (1+alph)/2*dt*vpa)*g_{inh}^{n+1} = RHS = gext
             call stream_tridiagonal_solve (iky, ie, iv, is, gext(:,ivmu))
          end if
          
       end do
    end if

    ! we now have g on the extended zed domain at this ky and set of connected kx values
    ! corresponding to a unit impulse in phi at this location
    call get_fields_for_response_matrix (gext, phiext, iky, ie)

    ! next need to create column in response matrix from phiext
    ! negative sign because matrix to be inverted in streaming equation
    ! is identity matrix - response matrix
    ! add in contribution from identity matrix
    phiext(idx) = phiext(idx)-1.0
    response_matrix(iky)%eigen(ie)%zloc(:,idx) = -phiext(:nresponse)

    if (allocated(gradpar_fac)) deallocate (gradpar_fac)

  end subroutine get_response_matrix_column

  subroutine sweep_zed_zonal_response (iv, is, sgn, g)

    use zgrid, only: nzgrid, delzed, nztot
    use run_parameters, only: zed_upwind, time_upwind
    use parallel_streaming, only: stream_c

    implicit none

    integer, intent (in) :: iv, is, sgn
    complex, dimension (:), intent (in out) :: g

    integer :: iz, iz1, iz2
    real :: fac1, fac2
    complex, dimension (:), allocatable :: gcf, gpi

    allocate (gpi(-nzgrid:nzgrid))
    allocate (gcf(-nzgrid:nzgrid))
    ! ky=0 is 2pi periodic (no extended zgrid)
    ! decompose into complementary function + particular integral
    ! zero BC for particular integral
    ! unit BC for complementary function (no source)
    if (sgn < 0) then
       iz1 = -nzgrid ; iz2 = nzgrid
    else
       iz1 = nzgrid ; iz2 = -nzgrid
    end if
    gpi(iz1) = 0. ; gcf(iz1) = 1.
    do iz = iz1-sgn, iz2, -sgn
       fac1 = 1.0+zed_upwind+sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
       fac2 = 1.0-zed_upwind-sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
       gpi(iz) = (-gpi(iz+sgn)*fac2 + 2.0*g(iz+nzgrid+1))/fac1
       gcf(iz) = -gcf(iz+sgn)*fac2/fac1
    end do
    ! g = g_PI + (g_PI(pi)/(1-g_CF(pi))) * g_CF
    g = gpi + (spread(gpi(iz2),1,nztot)/(1.-gcf(iz2)))*gcf
    deallocate (gpi, gcf)

  end subroutine sweep_zed_zonal_response

  subroutine get_fields_for_response_matrix (g, phi, iky, ie)

    use stella_layouts, only: vmu_lo
    use species, only: nspec, spec
    use species, only: has_electron_species
    use stella_geometry, only: dl_over_b
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: nsegments
    use kt_grids, only: zonal_mode, akx
    use vpamu_grids, only: integrate_species
    use gyro_averages, only: gyro_average
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg
    use fields, only: gamtot, gamtot3

    implicit none

    complex, dimension (:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:), intent (out) :: phi
    integer, intent (in) :: iky, ie
    
    integer :: idx, iseg, ikx, iz, ia
    integer :: izl_offset
    real, dimension (nspec) :: wgt
    complex, dimension (:), allocatable :: g0
    complex :: tmp

    ia = 1

    allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    wgt = spec%z*spec%dens
    phi = 0.

    idx = 0 ; izl_offset = 0
    iseg = 1
    ikx = ikxmod(iseg,ie,iky)
    if(zonal_mode(iky).and.abs(akx(ikx)) < epsilon(0.)) then
      phi(:) = 0.0
      return
    endif
    do iz = iz_low(iseg), iz_up(iseg)
       idx = idx + 1
       call gyro_average (g(idx,:), iky, ikx, iz, g0)
       call integrate_species (g0, iz, wgt, phi(idx))
       phi(idx) = phi(idx)/gamtot(iky,ikx,iz)
    end do
    izl_offset = 1
    if (nsegments(ie,iky) > 1) then
       do iseg = 2, nsegments(ie,iky)
          ikx = ikxmod(iseg,ie,iky)
          do iz = iz_low(iseg)+izl_offset, iz_up(iseg)
             idx = idx + 1
             call gyro_average (g(idx,:), iky, ikx, iz, g0)
             call integrate_species (g0, iz, wgt, phi(idx))
             phi(idx) = phi(idx)/gamtot(iky,ikx,iz)
          end do
          if (izl_offset == 0) izl_offset = 1
       end do
    end if

    if (.not.has_electron_species(spec) .and. &
         adiabatic_option_switch == adiabatic_option_fieldlineavg) then
       if (zonal_mode(iky)) then
          ! no connections for ky = 0
          iseg = 1 
          tmp = sum(dl_over_b(ia,:)*phi)
          phi = phi + tmp*gamtot3(ikxmod(1,ie,iky),:)
       end if
    end if

    deallocate (g0)

  end subroutine get_fields_for_response_matrix

  subroutine finish_response_matrix

    use fields_arrays, only: response_matrix

    implicit none
    
    if (allocated(response_matrix)) deallocate (response_matrix)

    response_matrix_initialized = .false.

  end subroutine finish_response_matrix

end module response_matrix
