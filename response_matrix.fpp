module response_matrix

  use netcdf

  implicit none

  public :: init_response_matrix, finish_response_matrix
  public :: read_response_matrix
  public :: response_matrix_initialized

  private

  logical :: response_matrix_initialized = .false.
  integer, parameter :: mat_unit = 70

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
    use job_manage, only: time_message
    use mp, only: proc0, iproc
    use run_parameters, only: mat_gen

    implicit none

    integer :: iky, ie, iseg, iz
    integer :: ikx
    integer :: nz_ext, nresponse
    integer :: idx
    integer :: izl_offset, izup
    real :: dum
    complex, dimension (:), allocatable :: phiext
    complex, dimension (:,:), allocatable :: gext
    logical :: debug = .false.
    character(100) :: message_lu
    real, dimension (2) :: time_response_matrix_lu

    ! Related to the saving of the the matrices in netcdf format
    character(len=15) :: fmt, proc_str 
    character(len=100) :: file_name
    integer :: istatus
    istatus = 0

!   All matrices handled by processor i_proc are stored
!   on a single file named: response_mat.iproc
    fmt = '(I5.5)'
    if (mat_gen) THEN
       call check_directories

       write (proc_str, fmt) iproc
       file_name = './mat/response_mat.'//trim(proc_str)

       open(unit=mat_unit, status='replace', file=file_name, &
            position='rewind', action='write', form='unformatted')
       write(unit=mat_unit) naky
    end if

    if (response_matrix_initialized) return
    response_matrix_initialized = .true.
    
    if (.not.allocated(response_matrix)) allocate (response_matrix(naky))
    
    ! for a given ky and set of connected kx values
    ! give a unit impulse to phi at each zed location
    ! in the extended domain and solve for h(zed_extended,(vpa,mu,s))

    
    do iky = 1, naky
       
       if (mat_gen) THEN
          write(unit=mat_unit) iky, neigen(iky)
       end if
       
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
          
          if (mat_gen) then
             write(unit=mat_unit) ie, nresponse
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
             call get_dgdphi_matrix_column (iky, ikx, iz, ie, idx, nz_ext, nresponse, phiext, gext)
          end do
          ! once we have used one segments, remaining segments
          ! have one fewer unique zed point
          izl_offset = 1
          if (nsegments(ie,iky) > 1) then
             do iseg = 2, nsegments(ie,iky)
                ikx = ikxmod(iseg,ie,iky)
                do iz = iz_low(iseg)+izl_offset, iz_up(iseg)
                   idx = idx + 1
                   call get_dgdphi_matrix_column (iky, ikx, iz, ie, idx, nz_ext, nresponse, phiext, gext)
                end do
                if (izl_offset == 0) izl_offset = 1
             end do
          end if
          deallocate (gext,phiext)
       enddo
       !DSO - This ends parallelization over velocity space.
       !      At this point every processor has int dv dgdphi for a given ky
       !      and so the quasineutrality solve and LU decomposition can be
       !      parallelized locally if need be.
       !      This is preferable to parallelization over ky as the LU 
       !      decomposition (and perhaps QN) will be dominated by the 
       !      ky with the most connections

       !      (perhaps further parallelization can be achieved by diviing up
       !       jtwist across nodes and communicating the results)
   


       ! solve quasineutrality 
       ! for local stella, this is a diagonal process, but global stella
       ! may require something more sophisticated

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

          allocate (phiext(nz_ext))

          do idx = 1, nresponse
             phiext(nz_ext) = 0.0
             phiext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(:,idx)
             call get_fields_for_response_matrix (phiext, iky, ie)

             ! next need to create column in response matrix from phiext
             ! negative sign because matrix to be inverted in streaming equation
             ! is identity matrix - response matrix
             ! add in contribution from identity matrix
             phiext(idx) = phiext(idx)-1.0
             response_matrix(iky)%eigen(ie)%zloc(:,idx) = -phiext(:nresponse)

          end do
       enddo

       !now we have the full response matrix. Finally, perform its LU decomposition
       do ie = 1, neigen(iky)
          ! now that we have the reponse matrix for this ky and set of connected kx values
          ! get the LU decomposition so we are ready to solve the linear system
          call lu_decomposition (response_matrix(iky)%eigen(ie)%zloc,response_matrix(iky)%eigen(ie)%idx,dum)
          
          if(proc0.and.debug) then
             call time_message(.false., time_response_matrix_lu, message_lu)
          end if
          
          if (mat_gen) then
             write(unit=mat_unit) response_matrix(iky)%eigen(ie)%idx
             write(unit=mat_unit) response_matrix(iky)%eigen(ie)%zloc
          end if
          
          deallocate (phiext)
       end do
    end do
    if (mat_gen) then
       close(unit=mat_unit)
    end if

  end subroutine init_response_matrix

  subroutine read_response_matrix
    
    use fields_arrays, only: response_matrix
    use common_types, only: response_matrix_type
    use mp, only: iproc

    implicit none

    integer :: iky, ie
    integer :: iky_dump, neigen_dump, naky_dump
    integer :: nresponse
    character(len=15) :: fmt, proc_str 
    character(len=100) :: file_name
    integer :: istatus, ie_dump, istat
    logical, parameter :: debug=.false.
    istatus = 0

!   All matrices handled by the processor i_proc are read
!   from a single file named: responst_mat.iproc
    fmt = '(I5.5)'
    write (proc_str, fmt) iproc
    file_name = './mat/response_mat.'//trim(proc_str)

    open(unit=mat_unit, status='old', file=file_name, &
         action='read', form='unformatted', iostat=istat)
    if (istat /= 0) then
       print *, 'Error opening response_matrix by processor', proc_str
    end if
!
    read(unit=mat_unit) naky_dump
!   
    if (response_matrix_initialized) return
    response_matrix_initialized = .true.
    
    if (.not.allocated(response_matrix)) allocate (response_matrix(naky_dump))
    
    do iky = 1, naky_dump
       read(unit=mat_unit) iky_dump, neigen_dump
       
       if (.not.associated(response_matrix(iky)%eigen)) &
            allocate (response_matrix(iky)%eigen(neigen_dump))
       
       do ie = 1, neigen_dump
          read(unit=mat_unit) ie_dump, nresponse
          
          if (.not.associated(response_matrix(iky)%eigen(ie)%zloc)) &
               allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse,nresponse))
          
          if (.not.associated(response_matrix(iky)%eigen(ie)%idx)) &
               allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))
          
          read(unit=mat_unit) response_matrix(iky)%eigen(ie)%idx
          read(unit=mat_unit) response_matrix(iky)%eigen(ie)%zloc
       end do
    end do
    close (mat_unit)
    
    if (debug) then
       print *, 'File', file_name, ' successfully read by proc: ', proc_str
    end if

  end subroutine read_response_matrix

  subroutine get_dgdphi_matrix_column (iky, ikx, iz, ie, idx, nz_ext, nresponse, phiext, gext)

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
    ! now integrate over velocities to get a square response matrix
    ! (this ends the parallelization over velocity space, so every core should have a 
    !  copy of phiext)
    call integrate_over_velocity (gext, phiext, iky, ie)

    response_matrix(iky)%eigen(ie)%zloc(:,idx) = phiext(:nresponse)

    if (allocated(gradpar_fac)) deallocate (gradpar_fac)

  end subroutine get_dgdphi_matrix_column

! subroutine get_phi_matrix
! end subroutine get_phi_matrix

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

  subroutine integrate_over_velocity(g,phi,iky,ie)

    use stella_layouts, only: vmu_lo
    use species, only: nspec, spec
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: nsegments
    use kt_grids, only: zonal_mode, akx
    use vpamu_grids, only: integrate_species
    use gyro_averages, only: gyro_average

    implicit none

    complex, dimension (:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:), intent (out) :: phi
    integer, intent (in) :: iky, ie
 
    integer :: idx, iseg, ikx, iz, ia
    integer :: izl_offset
    real, dimension (nspec) :: wgt
    complex, dimension (:), allocatable :: g0

    ia = 1

    allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    wgt = spec%z*spec%dens_psi0
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
    end do
    izl_offset = 1
    if (nsegments(ie,iky) > 1) then
       do iseg = 2, nsegments(ie,iky)
          ikx = ikxmod(iseg,ie,iky)
          do iz = iz_low(iseg)+izl_offset, iz_up(iseg)
             idx = idx + 1
             call gyro_average (g(idx,:), iky, ikx, iz, g0)
             call integrate_species (g0, iz, wgt, phi(idx))
          end do
          if (izl_offset == 0) izl_offset = 1
       end do
    end if

  end subroutine integrate_over_velocity

  subroutine get_fields_for_response_matrix (phi, iky, ie)

    use stella_layouts, only: vmu_lo
    use species, only: nspec, spec
    use species, only: has_electron_species
    use stella_geometry, only: dl_over_b
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: nsegments
    use kt_grids, only: zonal_mode, akx
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg
    use fields, only: gamtot, gamtot3

    implicit none

    complex, dimension (:), intent (inout) :: phi
    integer, intent (in) :: iky, ie
    
    integer :: idx, iseg, ikx, iz, ia
    integer :: izl_offset
    complex, dimension (:), allocatable :: g0
    complex :: tmp

    ia = 1

    allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    idx = 0 ; izl_offset = 0
    iseg = 1
    ikx = ikxmod(iseg,ie,iky)
    if(zonal_mode(iky).and.abs(akx(ikx)) < epsilon(0.)) then
      phi(:) = 0.0
      return
    endif
    do iz = iz_low(iseg), iz_up(iseg)
       idx = idx + 1
       phi(idx) = phi(idx)/gamtot(iky,ikx,iz)
    end do
    izl_offset = 1
    if (nsegments(ie,iky) > 1) then
       do iseg = 2, nsegments(ie,iky)
          ikx = ikxmod(iseg,ie,iky)
          do iz = iz_low(iseg)+izl_offset, iz_up(iseg)
             idx = idx + 1
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


!----------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------!


  subroutine parallel_LU_decomposition (iky)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use fields_arrays, only: response_matrix
    use mp, only: barrier, broadcast, scope, job, mp_comm
    use mp, only: sharedprocs, subprocs, iproc, nproc
    use job_manage, only: njobs
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: neigen, ikxmod
    use extended_zgrid, only: nsegments
    use extended_zgrid, only: nzed_segment
    use mpi
    use linear_solve, only: imaxloc

    implicit none

    integer, intent (in) :: iky
    integer, dimension (:), allocatable :: job_list
    integer, dimension (:), allocatable :: row_limits

    real, parameter :: zero = 1.0e-20
    real, dimension (:), allocatable :: vv
    complex, dimension (:), allocatable :: dum
    type(C_PTR) :: bptr

    complex, dimension (:,:), pointer :: lu

    integer :: idx, ijob, i,j,k,l,ie,n
    integer :: imax, jroot, neig, win_size, win, ierr
    integer :: rdiv, rmod
    integer :: disp_unit = 1 

    call scope (sharedprocs)

    allocate (job_list(nproc))
    allocate (row_limits(0:nproc))

    job_list(iproc) = job

    call barrier

    do ijob = 0, njobs-1
      jroot = -1

      do j = 1, nproc
        if(job_list(j).eq.job) then
          jroot = j !the first processor on this job will be the root process
          neig = neigen(iky)
          continue
        endif
      enddo
      
      if(jroot.eq.-1) continue !no processors on this core are on this job
      
      ! broadcast number of matrices
      call broadcast(neig,jroot)
      
      do ie = 1, neig

        win_size = 0
        if(iproc.eq.jroot) then
          n = size(response_matrix(iky)%eigen(ie)%idx)
          win_size = int(n*n,MPI_ADDRESS_KIND)*2*8_MPI_ADDRESS_KIND
        endif

        !broadcast size of matrix
        call broadcast(n, jroot)

        !allocate the window
        call mpi_win_allocate_share(win_size,disp_unit,MPI_INFO_NULL,mp_comm,bptr,win,ierr)

        if(iproc.ne.jroot) then
          !make sure all the procs have the right memory address
          call mpi_shared_query(win, jroot, win_size, disp_unit, bptr, ierr)
        endif

        ! bind this c_ptr to our fortran matrix
        call c_f_pointer(bptr,lu,(/n,n/))

        !load the matrix
        if(iproc.eq.jroot) lu = response_matrix(iky)%eigen(ie)%zloc

        !syncronize the processors
        call mpi_win_fence(0,win,ierr)

        ! All the processors have the matrix. 
        ! Now perform LU decomposition

        vv = maxval(cabs(lu),dim=2)
        if (any(vv==0.0)) &
           write (*,*) 'singular matrix in lu_decomposition'
        vv = 1.0/vv
        do j = 1, n
           !divide up the work using row_limits
           rdiv = (n-j)/nproc
           rmod = mod(n-j,nproc)
           row_limits(0) = j+1
           if(rdiv.eq.0) then
             row_limits(rmod+1:) = -1
             do k=1,rmod
               row_limits(k)  = row_limits(k-1) + 1
             enddo
           else
             do k=1,nproc
               row_limits(k) = row_limits(k-1) + rdiv
               if(k.le.rmod) row_limits(k) = row_limits(k) + 1
             enddo
           endif

           !pivot if needed
           imax = (j-1) + imaxloc(vv(j:n)*cabs(lu(j:n,j)))
           if (j /= imax) then
              dum = lu(imax,:)
              lu(imax,:) = lu(j,:)
              lu(j,:) = dum
              vv(imax) = vv(j)
           end if
           if(job.eq.ijob) response_matrix(iky)%eigen(ie)%idx(j) = imax

           !get the lead multiplier
           if (lu(j,j)==0.0) lu(j,j) = zero
           do i = row_limits(iproc-1), row_limits(iproc)
             lu(i,j) = lu(i,j)/lu(j,j)
           enddo

           do k=row_limits(iproc-1), row_limits(iproc)
             do i = j+1,n
                lu(i,k) = lu(i,k) - lu(i,j)*lu(j,k)
             enddo
           enddo
        enddo

        !copy the decomposed matrix over
        if(job.eq.ijob) response_matrix(iky)%eigen(ie)%zloc = lu

      enddo
    enddo

    deallocate (job_list, row_limits)



    call scope(subprocs)

  end subroutine parallel_LU_decomposition
  
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!! CHECK_DIRECTORIES !!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>if directories mat does not exist it is generated.
  subroutine check_directories
     
# if defined(__INTEL_COMPILER)
!   for system call with intel compiler
    use ifport, only: system
# endif
     
    integer :: ierr
    logical :: mat_exists

# if defined(__INTEL_COMPILER)    
    inquire(directory='mat', exist=mat_exists)
# else
    inquire(file="./mat/.", exist=mat_exists)
# endif

    if (.not. mat_exists) then
       ierr=system('mkdir -p mat')
    end if
  end subroutine check_directories

end module response_matrix
