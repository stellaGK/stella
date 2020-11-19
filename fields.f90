module fields

  use common_types, only: eigen_type

  implicit none

  public :: init_fields, finish_fields
  public :: advance_fields, get_fields
  public :: get_radial_correction
  public :: enforce_reality_field
  public :: get_fields_by_spec
  public :: gamtot, gamtot3
  public :: dgamtotdr, dgamtot3dr
  public :: time_field_solve
  public :: fields_updated
  public :: get_dchidy, get_dchidx

  private

  real, dimension (:,:,:), allocatable :: gamtot, apar_denom
  real, dimension (:,:), allocatable :: gamtot3
  real :: gamtot_h, gamtot3_h

  real, dimension (:,:,:), allocatable :: dgamtotdr
  real, dimension (:,:), allocatable :: dgamtot3dr

  complex, dimension (:,:), allocatable :: save1, save2, b_mat

  real, dimension (2,2) :: time_field_solve

  type (eigen_type), dimension (:,:), allocatable :: phi_solve
  type (eigen_type) :: phizf_solve
  logical :: fields_updated = .false.
  logical :: fields_initialized = .false.
  logical :: debug = .false.
  integer :: zm !this, along with zmi, keeps track of whether we want to 
                !include the zero mode in our radial solver or not

  interface get_dchidy
     module procedure get_dchidy_4d
     module procedure get_dchidy_2d
  end interface


contains

  subroutine init_fields

    use mp, only: sum_allreduce
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, onlY: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
    use dist_fn_arrays, only: kperp2, dkperp2dr
    use gyro_averages, only: aj0v, aj1v
    use run_parameters, only: fphi, fapar
    use run_parameters, only: ky_solve_radial
    use physics_parameters, only: tite, nine, beta
    use physics_flags, only: radial_variation
    use species, only: spec, has_electron_species, ion_species
    use stella_geometry, only: dl_over_b, d_dl_over_b_drho, dBdrho, bmag
    use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: nvpa, nmu, mu
    use vpamu_grids, only: vpa, vperp2
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use vpamu_grids, only: integrate_vmu
    use species, only: spec
    use kt_grids, only: naky, nakx, akx, nx
    use kt_grids, only: zonal_mode, rho_d
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg
    use linear_solve, only: lu_decomposition, lu_inverse
    

    implicit none

    integer :: ikxkyz, iz, it, ikx, iky, is, ia, jkx, zmi
    real :: efac, efacp, tmp, tmp2, wgt, dum
    real, dimension (:,:), allocatable :: g0
    real, dimension (:), allocatable :: g1

    complex, dimension (:,:), allocatable :: g0k, g0x, a_inv, a_fsa

    ia = 1

    ! do not see why this is before fields_initialized check below
    call allocate_arrays

    if (fields_initialized) return
    fields_initialized = .true.

    ! could move these array allocations to allocate_arrays to clean up code
    if (.not.allocated(gamtot)) allocate (gamtot(naky,nakx,-nzgrid:nzgrid)) ; gamtot = 0.
    if (.not.allocated(gamtot3)) then
       if (.not.has_electron_species(spec) &
            .and. adiabatic_option_switch==adiabatic_option_fieldlineavg) then
          allocate (gamtot3(nakx,-nzgrid:nzgrid)) ; gamtot3 = 0.
       else
          allocate (gamtot3(1,1)) ; gamtot3 = 0.
       end if
    end if
    if (.not.allocated(apar_denom)) then
       if (fapar > epsilon(0.0)) then
          allocate (apar_denom(naky,nakx,-nzgrid:nzgrid)) ; apar_denom = 0.
       else
          allocate (apar_denom(1,1,1)) ; apar_denom = 0.
       end if
    end if

    if (radial_variation) then
      if (.not.allocated(dgamtotdr)) allocate(dgamtotdr(naky,nakx,-nzgrid:nzgrid)) ; dgamtotdr=0.
      if (.not.allocated(dgamtot3dr)) then
        if (.not.has_electron_species(spec) &
            .and. adiabatic_option_switch==adiabatic_option_fieldlineavg) then
          allocate (dgamtot3dr(nakx,-nzgrid:nzgrid)) ; dgamtot3dr = 0.
          allocate (save1(nakx,ntubes)) ; save1 = 0.
          allocate (save2(nakx,ntubes)) ; save2 = 0.
        else
          allocate (dgamtot3dr(1,1)) ; dgamtot3dr = 0.
        endif
      endif
    else
      if (.not.allocated(dgamtotdr))  allocate(dgamtotdr(1,1,1)) ; dgamtotdr = 0.
      if (.not.allocated(dgamtot3dr)) allocate (dgamtot3dr(1,1)) ; dgamtot3dr = 0.
    endif

    if (fphi > epsilon(0.0)) then
       allocate (g0(nvpa,nmu))
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          it = it_idx(kxkyz_lo,ikxkyz)
          ! gamtot does not depend on flux tube index,
          ! so only compute for one flux tube index
          if (it /= 1) cycle
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          g0 = spread((1.0 - aj0v(:,ikxkyz)**2),1,nvpa) &
               * spread(maxwell_vpa(:,is),2,nmu)*spread(maxwell_mu(ia,iz,:,is),1,nvpa)*maxwell_fac(is)
          wgt = spec(is)%z*spec(is)%z*spec(is)%dens_psi0/spec(is)%temp
          call integrate_vmu (g0, iz, tmp)
          gamtot(iky,ikx,iz) = gamtot(iky,ikx,iz) + tmp*wgt
       end do
       call sum_allreduce (gamtot)
       ! avoid divide by zero when kx=ky=0
       ! do not evolve this mode, so value is irrelevant
       !if (zonal_mode(1).and.akx(1)<epsilon(0.).and.has_electron_species(spec)) gamtot(1,1,:) = 1.0
       if (zonal_mode(1).and.akx(1)<epsilon(0.)) gamtot(1,1,:) = 0.0

       gamtot_h = sum(spec%z*spec%z*spec%dens/spec%temp)

       if (radial_variation) then
         allocate (g1(nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
           it = it_idx(kxkyz_lo,ikxkyz)
           ! gamtot does not depend on flux tube index,
           ! so only compute for one flux tube index
           if (it /= 1) cycle
           iky = iky_idx(kxkyz_lo,ikxkyz)
           ikx = ikx_idx(kxkyz_lo,ikxkyz)
           iz = iz_idx(kxkyz_lo,ikxkyz)
           is = is_idx(kxkyz_lo,ikxkyz)
           g1 = aj0v(:,ikxkyz)*aj1v(:,ikxkyz)*(spec(is)%smz)**2 &
              * (kperp2(iky,ikx,ia,iz)*vperp2(ia,iz,:)/bmag(ia,iz)**2) &
              * (dkperp2dr(iky,ikx,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
              / (1.0 - aj0v(:,ikxkyz)**2 + 100.*epsilon(0.0))

           g0 = spread((1.0 - aj0v(:,ikxkyz)**2),1,nvpa) &
              * spread(maxwell_vpa(:,is),2,nmu)*spread(maxwell_mu(ia,iz,:,is),1,nvpa)*maxwell_fac(is) &
              * (-spec(is)%tprim*(spread(vpa**2,2,nmu)+spread(vperp2(ia,iz,:),1,nvpa)-2.5) &
                 -spec(is)%fprim + (dBdrho(iz)/bmag(ia,iz))*(1.0 - 2.0*spread(mu,1,nvpa)*bmag(ia,iz)) &
                 + spread(g1,1,nvpa))
           wgt = spec(is)%z*spec(is)%z*spec(is)%dens/spec(is)%temp
           call integrate_vmu (g0, iz, tmp)
           dgamtotdr(iky,ikx,iz) = dgamtotdr(iky,ikx,iz) + tmp*wgt
         end do
         call sum_allreduce (dgamtotdr)

         deallocate (g1)

       endif
       if (zonal_mode(1).and.akx(1)<epsilon(0.).and.has_electron_species(spec)) then
         gamtot(1,1,:)    = 0.0
         dgamtotdr(1,1,:) = 0.0
         zm = 1
       endif

       if (.not.has_electron_species(spec)) then
          efac = tite/nine * (spec(ion_species)%dens/spec(ion_species)%temp) 
          efacp = efac*(spec(ion_species)%tprim - spec(ion_species)%fprim)
          gamtot   = gamtot   + efac
          gamtot_h = gamtot_h + efac
          if(radial_variation) dgamtotdr = dgamtotdr + efacp
          if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
             if (zonal_mode(1)) then
                zm = 1
                gamtot3_h = tite/(nine*sum(spec%zt*spec%z*spec%dens))
                do ikx = 1, nakx
                   ! avoid divide by zero for kx=ky=0 mode,
                   ! which we do not need anyway
                   !if (abs(akx(ikx)) < epsilon(0.)) cycle
                   tmp = 1./efac - sum(dl_over_b(ia,:)/gamtot(1,ikx,:))
                   gamtot3(ikx,:) = 1./(gamtot(1,ikx,:)*tmp)
                   if (radial_variation) then
                     tmp2 = (spec(ion_species)%tprim - spec(ion_species)%fprim)/efac &
                            + sum(d_dl_over_b_drho(ia,:)/gamtot(1,ikx,:)) &
                            - sum(dl_over_b(ia,:)*dgamtotdr(1,ikx,:) &
                                / gamtot(1,ikx,:)**2)
                     dgamtot3dr(ikx,:)  = gamtot3(ikx,:) & 
                                        * (-dgamtotdr(1,ikx,:)/gamtot(1,ikx,:) + tmp2/tmp)
                   endif
                end do
                if(akx(1)<epsilon(0.)) dgamtot3dr(1,:) = 0.0
             end if
          end if
       end if

       if(radial_variation) then

         allocate (g0k(1,nakx))
         allocate (g0x(1,nakx))


         if(.not.allocated(phi_solve)) allocate(phi_solve(min(ky_solve_radial,naky),-nzgrid:nzgrid))

         do iky = 1, min(ky_solve_radial,naky)
           zmi = 0
           if(iky.eq.1) zmi=zm !zero mode may or may not be included in matrix
           do iz = -nzgrid, nzgrid
             if(.not.associated(phi_solve(iky,iz)%zloc)) &
               allocate(phi_solve(iky,iz)%zloc(nakx-zmi,nakx-zmi))
             if(.not.associated(phi_solve(iky,iz)%idx))  &
               allocate(phi_solve(iky,iz)%idx(nakx-zmi))

             phi_solve(iky,iz)%zloc = 0.0
             phi_solve(iky,iz)%idx = 0
             do ikx = 1+zmi, nakx
               g0k(1,:) = 0.0
               g0k(1,ikx) = dgamtotdr(1,ikx,iz)

               call transform_kx2x_unpadded (g0k,g0x)
               g0x(1,:) = rho_d*g0x(1,:)
               call transform_x2kx_unpadded(g0x,g0k)

               !column row
               phi_solve(iky,iz)%zloc(:,ikx-zmi) = g0k(1,(1+zmi):)
               phi_solve(iky,iz)%zloc(ikx-zmi,ikx-zmi) = phi_solve(iky,iz)%zloc(ikx-zmi,ikx-zmi) &
                                                     + gamtot(1,ikx,iz)
             enddo

             call lu_decomposition(phi_solve(iky,iz)%zloc, phi_solve(iky,iz)%idx, dum)
!          call zgetrf(nakx-zmi,nakx-zmi,phi_solve(iky,iz)%zloc,nakx-zmi,phi_solve(iky,iz)%idx,info)
           enddo
         enddo

         if (.not.has_electron_species(spec).and. & 
              adiabatic_option_switch == adiabatic_option_fieldlineavg) then

           if(.not.allocated(b_mat)) allocate(b_mat(nakx-zm,nakx-zm));
           
           allocate(a_inv(nakx-zm,nakx-zm))
           allocate(a_fsa(nakx-zm,nakx-zm)); a_fsa = 0.0

           if(.not.associated(phizf_solve%zloc)) &
               allocate(phizf_solve%zloc(nakx-zm,nakx-zm));
           phizf_solve%zloc = 0.0

           if(.not.associated(phizf_solve%idx)) allocate(phizf_solve%idx(nakx-zm));

           do ikx = 1+zm, nakx
             g0k(1,:) = 0.0
             g0k(1,ikx) = 1.0

             call transform_kx2x_unpadded (g0k,g0x)
             g0x(1,:) = (efac + efacp*rho_d)*g0x(1,:)
             call transform_x2kx_unpadded(g0x,g0k)

             !column row
             b_mat(:,ikx-zm) = g0k(1,(1+zm):) 
           enddo

           !get inverse of A
           do iz = -nzgrid, nzgrid

             call lu_inverse(phi_solve(iky,iz)%zloc, phi_solve(iky,iz)%idx, &
                             a_inv)

             !flux surface average it
             do ikx = 1, nakx-zm
               g0k(1,1) = 0
               g0k(1,(1+zm):) = a_inv(:,ikx)

               call transform_kx2x_unpadded (g0k,g0x)
               g0x(1,:) = (dl_over_b(ia,iz) + d_dl_over_b_drho(ia,iz)*rho_d)*g0x(1,:)
               call transform_x2kx_unpadded(g0x,g0k)

               a_fsa(:,ikx) = a_fsa(:,ikx) + g0k(1,(1+zm):) 
             enddo
           enddo

           ! calculate I - <A^-1>B
           do ikx = 1, nakx-zm
             do jkx = 1, nakx-zm
               phizf_solve%zloc(ikx,jkx) = -sum(a_fsa(ikx,:)*b_mat(:,jkx))
             enddo
           enddo
           do ikx = 1, nakx-zm
             phizf_solve%zloc(ikx,ikx) = 1.0 + phizf_solve%zloc(ikx,ikx)
           enddo

           call lu_decomposition(phizf_solve%zloc,phizf_solve%idx, dum)

           deallocate(a_inv,a_fsa)

         endif

         deallocate(g0k,g0x)
         deallocate (g0)
       endif
    end if

    if (fapar > epsilon(0.)) then
       allocate (g0(nvpa,nmu))
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          it = it_idx(kxkyz_lo,ikxkyz)
          ! apar_denom does not depend on flux tube index,
          ! so only compute for one flux tube index
          if (it /= 1) cycle
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          g0 = spread(maxwell_vpa(:,is)*vpa**2,2,nmu)*maxwell_fac(is) &
               * spread(maxwell_mu(ia,iz,:,is)*aj0v(:,ikxkyz)**2,1,nvpa)
          wgt = 2.0*beta*spec(is)%z*spec(is)%z*spec(is)%dens/spec(is)%mass
          call integrate_vmu (g0, iz, tmp)
          apar_denom(iky,ikx,iz) = apar_denom(iky,ikx,iz) + tmp*wgt
       end do
       call sum_allreduce (apar_denom)
       apar_denom = apar_denom + kperp2(:,:,ia,:)

       deallocate (g0)
    end if

!    filename=trim(run_name)//".gamtot"
!    open(3636,file=trim(filename),status='unknown')
!    do iky = 1, naky
!      do ikx = 1, nakx
!        do iz = -nzgrid,nzgrid
!          write(3636,'(4e25.8)') gamtot(iky,ikx,iz), dgamtotdr(iky,ikx,iz), &
!                                 gamtot3(ikx,iz), dgamtot3dr(ikx,iz)
!        enddo
!      enddo
!    enddo
!    close(3636)


!    if (wstar_implicit) call init_get_fields_wstar

  end subroutine init_fields


  subroutine allocate_arrays

    use fields_arrays, only: phi, apar, phi_old
    use fields_arrays, only: phi_corr_QN, phi_corr_GA
    use fields_arrays, only: apar_corr_QN, apar_corr_GA
    use zgrid, only: nzgrid, ntubes
    use stella_layouts, only: vmu_lo
    use physics_flags, only: radial_variation
    use kt_grids, only: naky, nakx

    implicit none

    if (.not.allocated(phi)) then
       allocate (phi(naky,nakx,-nzgrid:nzgrid,ntubes))
       phi = 0.
    end if
    if (.not. allocated(apar)) then
       allocate (apar(naky,nakx,-nzgrid:nzgrid,ntubes))
       apar = 0.
    end if
    if (.not.allocated(phi_old)) then
       allocate (phi_old(naky,nakx,-nzgrid:nzgrid,ntubes))
       phi_old = 0.
    end if
    if (.not.allocated(phi_corr_QN) .and. radial_variation) then
       allocate (phi_corr_QN(naky,nakx,-nzgrid:nzgrid,ntubes))
       phi_corr_QN = 0.
    end if
    if (.not.allocated(phi_corr_GA) .and. radial_variation) then
       allocate (phi_corr_GA(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       phi_corr_GA = 0.
    end if
    if (.not.allocated(apar_corr_QN) .and. radial_variation) then
       !allocate (apar_corr(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       allocate (apar_corr_QN(1,1,1,1))
       apar_corr_QN = 0.
    end if
    if (.not.allocated(apar_corr_GA) .and. radial_variation) then
       !allocate (apar_corr(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       allocate (apar_corr_GA(1,1,1,1,1))
       apar_corr_GA = 0.
    end if

  end subroutine allocate_arrays

  subroutine enforce_reality_field(fin)

!DSO> while most of the modes in the box have reality built in (as we 
!     throw out half the kx-ky plane, modes with ky=0 do not have
!     this enforcement built in. In theory this shouldn't be a problem
!     as these modes should be stable, but I made this function (and 
!     its relative in the dist file) just in case

    use kt_grids, only: nakx
    use zgrid, only: nzgrid
    
    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (inout) :: fin

    integer ikx

    fin(1,1,:,:) = real(fin(1,1,:,:))
    do ikx = 2, nakx/2+1
      fin(1,ikx,:,:) = 0.5*(fin(1,ikx,:,:) + conjg(fin(1,nakx-ikx+2,:,:)))
      fin(1,nakx-ikx+2,:,:) = conjg(fin(1,ikx,:,:))
    enddo

  end subroutine enforce_reality_field

  subroutine advance_fields (g, phi, apar, dist)

    use mp, only: proc0
    use stella_layouts, only: vmu_lo
    use job_manage, only: time_message
    use redistribute, only: scatter
    use dist_fn_arrays, only: gvmu
    use zgrid, only: nzgrid
    use dist_redistribute, only: kxkyz2vmu
    use run_parameters, only: fields_kxkyz

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: phi, apar
    character (*), intent (in) :: dist

    if (fields_updated) return

    ! time the communications + field solve
    if (proc0) call time_message(.false.,time_field_solve(:,1),' fields')
    if (fields_kxkyz) then
       ! first gather (vpa,mu) onto processor for v-space operations
       ! v-space operations are field solve, dg/dvpa, and collisions
       if (debug) write (*,*) 'dist_fn::advance_stella::scatter'
       if (proc0) call time_message(.false.,time_field_solve(:,2),' fields_redist')
       call scatter (kxkyz2vmu, g, gvmu)
       if (proc0) call time_message(.false.,time_field_solve(:,2),' fields_redist')
       ! given gvmu with vpa and mu local, calculate the corresponding fields
       if (debug) write (*,*) 'dist_fn::advance_stella::get_fields'
       call get_fields (gvmu, phi, apar, dist)
    else
       call get_fields_vmulo (g, phi, apar, dist)
    end if

    ! set a flag to indicate that the fields have been updated
    ! this helps avoid unnecessary field solves
    fields_updated = .true.
    ! time the communications + field solve
    if (proc0) call time_message(.false.,time_field_solve(:,1),' fields')

  end subroutine advance_fields

  subroutine get_fields (g, phi, apar, dist)

    use mp, only: proc0
    use mp, only: sum_allreduce, mp_abort
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
    use dist_fn_arrays, only: kperp2
    use gyro_averages, only: gyro_average
    use run_parameters, only: fphi, fapar
    use physics_parameters, only: beta
    use physics_flags, only: radial_variation
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: vpa
    use vpamu_grids, only: integrate_vmu
    use kt_grids, only: nakx
    use species, only: spec, has_electron_species

    implicit none
    
    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: phi, apar
    character (*), intent (in) :: dist
    complex :: tmp

    real :: wgt
    complex, dimension (:,:), allocatable :: g0
    integer :: ikxkyz, iz, it, ikx, iky, is, ia

    ia = 1

    phi = 0.
    if (fphi > epsilon(0.0)) then
       allocate (g0(nvpa,nmu))
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iky = iky_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          call gyro_average (g(:,:,ikxkyz), ikxkyz, g0)
          wgt = spec(is)%z*spec(is)%dens_psi0
          call integrate_vmu (g0, iz, tmp)
          phi(iky,ikx,iz,it) = phi(iky,ikx,iz,it) + wgt*tmp
       end do
       deallocate (g0)
       call sum_allreduce (phi)

       call get_phi (phi, dist)

    end if

    apar = 0.
    if (fapar > epsilon(0.0)) then
       allocate (g0(nvpa,nmu))
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iky = iky_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          call gyro_average (spread(vpa,2,nmu)*g(:,:,ikxkyz), ikxkyz, g0)
          wgt = 2.0*beta*spec(is)%z*spec(is)%dens*spec(is)%stm
          call integrate_vmu (g0, iz, tmp)
          apar(iky,ikx,iz,it) = apar(iky,ikx,iz,it) + tmp*wgt
       end do
       call sum_allreduce (apar)
       if (dist == 'h') then
          apar = apar/spread(kperp2(:,:,ia,:),4,ntubes)
       else if (dist == 'gbar') then
          apar = apar/spread(apar_denom,4,ntubes)
       else if (dist == 'gstar') then
          write (*,*) 'APAR NOT SETUP FOR GSTAR YET. aborting.'
          call mp_abort('APAR NOT SETUP FOR GSTAR YET. aborting.')
       else
          if (proc0) write (*,*) 'unknown dist option in get_fields. aborting'
          call mp_abort ('unknown dist option in get_fields. aborting')
       end if
       deallocate (g0)
    end if
    
  end subroutine get_fields

  subroutine get_fields_vmulo (g, phi, apar, dist)

    use mp, only: mp_abort
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: imu_idx, is_idx
    use gyro_averages, only: gyro_average, aj0x, aj1x
    use run_parameters, only: fphi, fapar
    use stella_geometry, only: dBdrho, bmag
    use physics_flags, only: radial_variation
    use dist_fn_arrays, only: kperp2, dkperp2dr
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: integrate_species, vperp2
    use kt_grids, only: nakx, naky, nx, multiply_by_rho
    use species, only: spec

    implicit none
    
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: phi, apar
    character (*), intent (in) :: dist

    integer :: ivmu, iz, it, ia, imu, is
    complex, dimension (:,:,:), allocatable :: gyro_g
    complex, dimension (:,:), allocatable :: g0k

    ia = 1

    phi = 0.
    if (fphi > epsilon(0.0)) then
       allocate (g0k(naky,nakx))
       allocate (gyro_g(naky,nakx,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       do it = 1, ntubes
         do iz = -nzgrid, nzgrid
           do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
             is = is_idx(vmu_lo,ivmu)
             imu = imu_idx(vmu_lo,ivmu)
             call gyro_average (g(:,:,iz,it,ivmu), iz, ivmu, gyro_g(:,:,ivmu))
             g0k = 0.0
             if(radial_variation) then
               g0k = gyro_g(:,:,ivmu) &
                   * (-0.5*aj1x(:,:,iz,ivmu)/aj0x(:,:,iz,ivmu)*(spec(is)%smz)**2 & 
                   * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                   * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                   + dBdrho(iz)/bmag(ia,iz))

               !g0k(1,1) = 0.
               call multiply_by_rho(g0k)

             endif

             gyro_g(:,:,ivmu) = gyro_g(:,:,ivmu) + g0k

           end do
           call integrate_species (gyro_g, iz, spec%z*spec%dens_psi0, phi(:,:,iz,it))
         end do
       end do
       deallocate (gyro_g)

       call get_phi(phi, dist)

    end if
    
    apar = 0.
    if (fapar > epsilon(0.0)) then
       ! FLAG -- NEW LAYOUT NOT YET SUPPORTED !!
       call mp_abort ('APAR NOT YET SUPPORTED FOR NEW FIELD SOLVE. ABORTING.')
!        allocate (g0(-nvgrid:nvgrid,nmu))
!        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
!           iz = iz_idx(kxkyz_lo,ikxkyz)
!           ikx = ikx_idx(kxkyz_lo,ikxkyz)
!           iky = iky_idx(kxkyz_lo,ikxkyz)
!           is = is_idx(kxkyz_lo,ikxkyz)
!           g0 = spread(aj0v(:,ikxkyz),1,nvpa)*spread(vpa,2,nmu)*g(:,:,ikxkyz)
!           wgt = 2.0*beta*spec(is)%z*spec(is)%dens*spec(is)%stm
!           call integrate_vmu (g0, iz, tmp)
!           apar(iky,ikx,iz) = apar(iky,ikx,iz) + tmp*wgt
!        end do
!        call sum_allreduce (apar)
!        if (dist == 'h') then
!           apar = apar/kperp2
!        else if (dist == 'gbar') then
!           apar = apar/apar_denom
!        else if (dist == 'gstar') then
!           write (*,*) 'APAR NOT SETUP FOR GSTAR YET. aborting.'
!           call mp_abort('APAR NOT SETUP FOR GSTAR YET. aborting.')
!        else
!           if (proc0) write (*,*) 'unknown dist option in get_fields. aborting'
!           call mp_abort ('unknown dist option in get_fields. aborting')
!        end if
!        deallocate (g0)
    end if
    
  end subroutine get_fields_vmulo

  subroutine get_fields_by_spec (g, fld)

    use mp, only: sum_allreduce
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
    use gyro_averages, only: gyro_average
    use run_parameters, only: fphi
    use stella_geometry, only: dl_over_b
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: integrate_vmu
    use kt_grids, only: nakx
    use kt_grids, only: zonal_mode
    use species, only: spec, nspec, has_electron_species
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg

    implicit none
    
    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,:), intent (out) :: fld

    real :: wgt
    complex, dimension (:,:), allocatable :: g0
    integer :: ikxkyz, iz, it, ikx, iky, is, ia
    complex, dimension (nspec) :: tmp

    ia = 1

    fld = 0.
    if (fphi > epsilon(0.0)) then
       allocate (g0(nvpa,nmu))
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iky = iky_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          wgt = spec(is)%z*spec(is)%dens_psi0
          call gyro_average (g(:,:,ikxkyz), ikxkyz, g0)
          g0 = g0*wgt
          call integrate_vmu (g0, iz, fld(iky,ikx,iz,it,is))
       end do
       call sum_allreduce (fld)

       fld = fld/gamtot_h

       if (.not.has_electron_species(spec) .and. &
            adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          if (zonal_mode(1)) then
             do ikx = 1, nakx
                do it = 1, ntubes
                   do is = 1, nspec
                      tmp(is) = sum(dl_over_b(ia,:)*fld(1,ikx,:,it,is))
                      fld(1,ikx,:,it,is) = fld(1,ikx,:,it,is) + tmp(is)*gamtot3_h
                   end do
                end do
             end do
          end if
       end if

       deallocate (g0)
    end if

  end subroutine get_fields_by_spec

  subroutine get_phi (phi, dist)

    use mp, only: proc0, mp_abort, job
    use physics_flags, only: full_flux_surface, radial_variation
    use run_parameters, only: ky_solve_radial
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: swap_kxky_ordered, nakx, naky, rho_d, zonal_mode
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
    use stella_geometry, only: dl_over_b, d_dl_over_b_drho
    use linear_solve, only: lu_back_substitution
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg
    use species, only: spec, has_electron_species

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi
    integer :: ia, it, iz, ikx, iky, zmi
    complex, dimension (:,:), allocatable :: g0k, g0x
    complex, dimension (:), allocatable :: g_fsa
    complex :: tmp

    character (*), intent (in) :: dist

    ia = 1

    if (dist == 'h') then
      phi = phi/gamtot_h
    else if (dist == 'gbar') then
      if (full_flux_surface) then
!        ! need to invert 1-gamma0 operator
!        allocate (phi_swap(naky_all,ikx_max))
!        do it = 1, ntubes
!           do iz = -nzgrid, nzgrid
!              call swap_kxky_ordered (phi, phi_swap)
!              do ikx = 1, ikx_max
!                 ! if ky values uncoupled, then obtain phi by simple divide
!                 if (maxval(ia_max_gam0a(:,ikx,iz)) == 1) then
!                    phi_swap(:,ikx) = phi_swap(:,ikx)/gam0a(1,:,ikx,iz)
!                 else
!                    ! for ky values that require info from other ky values
!                    ! must solve linear system
!                    call 
!                 end if
!              end do
!           end do
!        end do
!        deallocate (phi_swap)
      else if (radial_variation) then
        allocate (g0k(1,nakx))
        allocate (g0x(1,nakx))

         do it = 1, ntubes
           do iz = -nzgrid, nzgrid
             do iky = 1, naky
               if(iky > ky_solve_radial) then
                 phi(iky,:,iz,it) = phi(iky,:,iz,it)/gamtot(iky,:,iz)
               else
                 zmi=0
                 if(iky.eq.1) zmi=zm
                 call lu_back_substitution(phi_solve(iky,iz)%zloc, &
                                          phi_solve(iky,iz)%idx, phi(iky,(1+zmi):,iz,it))
!                call zgetrs('n',nakx-zm,1,phizf_solve(iz)%zloc,nakx-zm,phizf_solve(iz)%idx, &
!                            phi(1,(1+zm):,iz,it),nakx-zm,info)
               endif
             enddo

           enddo

           !if(zm.eq.1) phi(1,zm,:,it) = 0.0

         enddo
!         
!        g0k(1,:) = dgamtotdr(1,:,0)*phi(1,:,0,1)
!        call transform_kx2x_unpadded (g0k,g0x)
!        g0x(1,:) = rho_d*g0x(1,:)
!        call transform_x2kx_unpadded(g0x,g0k)
!        g0k(1,:) = g0k(1,:) + gamtot(1,:,0)*phi(1,:,0,1)
!        write (*,*) 'dsdft', g0k(1,:)

         deallocate (g0k,g0x)
      else
        phi = phi/spread(gamtot,4,ntubes)
        if(any(gamtot(1,1,:).lt.epsilon(0.))) phi(1,1,:,:) = 0.0
      end if

    else 
      if (proc0) write (*,*) 'unknown dist option in get_fields. aborting'
      call mp_abort ('unknown dist option in get_fields. aborting')
      return 
    end if


    if (.not.has_electron_species(spec) .and. &
      adiabatic_option_switch == adiabatic_option_fieldlineavg) then
      if (zonal_mode(1)) then
        if (dist == 'h') then
          do it = 1, ntubes
            do ikx = 1, nakx
              tmp = sum(dl_over_b(ia,:)*phi(1,ikx,:,it))
              phi(1,ikx,:,it) = phi(1,ikx,:,it) + tmp*gamtot3_h
            end do
          end do
        else if (dist == 'gbar'.and..not.(radial_variation.and.ky_solve_radial.gt.0)) then
          if(radial_variation) then
            do it = 1, ntubes
              do ikx = 1, nakx
                ! DSO - this is sort of hack in order to avoid extra communications
                !       However, get_radial_correction should be called immediately 
                !       after advance_fields, so it should be ok...
                save1(ikx,it) = sum(dl_over_b(ia,:)*phi(1,ikx,:,it))
                save2(ikx,it) = sum(d_dl_over_b_drho(ia,:)*phi(1,ikx,:,it))
              enddo
            enddo
          endif
          do ikx = 1, nakx
            do it = 1, ntubes
               tmp = sum(dl_over_b(ia,:)*phi(1,ikx,:,it))
               phi(1,ikx,:,it) = phi(1,ikx,:,it) + tmp*gamtot3(ikx,:)
            end do
          end do
        else if (dist == 'gbar'.and.(radial_variation.and.ky_solve_radial.gt.0)) then
          allocate (g0k(1,nakx))
          allocate (g0x(1,nakx))
          allocate (g_fsa(nakx-zm))

          do it = 1, ntubes
            g_fsa = 0.0
            do iz = -nzgrid, nzgrid
              g0k(1,:) = phi(1,:,iz,it)
              call transform_kx2x_unpadded (g0k,g0x)
              g0x(1,:) = (dl_over_b(ia,iz) + d_dl_over_b_drho(ia,iz)*rho_d)*g0x(1,:)
              call transform_x2kx_unpadded(g0x,g0k)

              g_fsa = g_fsa + g0k(1,(1+zm):)
            enddo
        
            call lu_back_substitution(phizf_solve%zloc,phizf_solve%idx, g_fsa)

            do ikx = 1,nakx-zm
              g0k(1,ikx+zm) = sum(b_mat(ikx,:)*g_fsa(:))
            enddo

            do iz = -nzgrid, nzgrid
              g_fsa = g0k(1,(1+zm):)
              call lu_back_substitution(phi_solve(1,iz)%zloc,phi_solve(1,iz)%idx, g_fsa)

              phi(1,(1+zm):,iz,it) = phi(1,(1+zm):,iz,it) + g_fsa
            enddo
          enddo

          deallocate(g0k,g0x,g_fsa)
        else 
          if (proc0) write (*,*) 'unknown dist option in get_fields. aborting'
          call mp_abort ('unknown dist option in get_fields. aborting')
        end if
      end if
    end if
    
    if(zm.eq.1) phi(1,zm,:,:) = 0.0

  end subroutine get_phi



  ! the following routine gets the correction in phi both from gyroaveraging and quasineutrality
  ! the output, phi, 
  subroutine get_radial_correction (g, phi_in, dist)

    use mp, only: proc0, mp_abort
    use stella_layouts, only: vmu_lo
    use gyro_averages, only: gyro_average, gyro_average_j1
    use gyro_averages, only: aj0x, aj1x
    use run_parameters, only: fphi, fapar, ky_solve_radial
    use stella_geometry, only: dl_over_b, bmag, dBdrho
    use stella_layouts, only: imu_idx, is_idx
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: integrate_species, vperp2
    use kt_grids, only: nakx, nx, naky
    use kt_grids, only: zonal_mode, multiply_by_rho
    use species, only: spec, has_electron_species
    use fields_arrays, only: phi_corr_QN, phi_corr_GA
    use dist_fn_arrays, only: kperp2, dkperp2dr
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg

    implicit none
    
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phi_in
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    character (*), intent (in) :: dist

    integer :: ikx, iky, ivmu, iz, it, ia, is, imu
    complex :: tmp
    complex, dimension (:,:,:,:), allocatable :: phi
    complex, dimension (:,:,:), allocatable :: gyro_g
    complex, dimension (:,:), allocatable :: g0k, g0x

    ia = 1

    if (fphi > epsilon(0.0)) then
       allocate (gyro_g(naky,nakx,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       allocate (g0k(naky,nakx))
       allocate (g0x(naky,nx))
       allocate (phi(naky,nakx,-nzgrid:nzgrid,ntubes))
       phi = 0.
       do it = 1, ntubes
         do iz = -nzgrid, nzgrid
           do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
             is = is_idx(vmu_lo,ivmu)
             imu = imu_idx(vmu_lo,ivmu)

             g0k = g(:,:,iz,it,ivmu) &
                 * (-0.5*aj1x(:,:,iz,ivmu)/aj0x(:,:,iz,ivmu) & 
                    * (spec(is)%smz)**2 & 
                    * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                    * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                 + dBdrho(iz)/bmag(ia,iz) - dgamtotdr(:,:,iz)/gamtot(:,:,iz))

             g0k(1,1) = 0.

             call gyro_average (g0k, iz, ivmu, gyro_g(:,:,ivmu))
           end do
           call integrate_species (gyro_g, iz, spec%z*spec%dens_psi0, phi(:,:,iz,it))
         end do
       end do


       if (dist == 'gbar') then
          !call get_phi (phi)
          phi = phi/spread(gamtot,4,ntubes)
          phi(1,1,:,:) = 0.0
       else if (dist == 'h') then
          if (proc0) write (*,*) 'dist option "h" not implemented in radial_correction. aborting'
          call mp_abort ('dist option "h" in radial_correction. aborting')
       else
          if (proc0) write (*,*) 'unknown dist option in radial_correction. aborting'
          call mp_abort ('unknown dist option in radial_correction. aborting')
       end if

       if (.not.has_electron_species(spec) .and. &
            adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          if (zonal_mode(1)) then
             if (dist == 'gbar') then
                do it = 1, ntubes
                   do ikx = 1, nakx
                      tmp = sum(dl_over_b(ia,:)*phi(1,ikx,:,it))
                      phi(1,ikx,:,it) = phi(1,ikx,:,it) &
                                      + tmp*gamtot3(ikx,:) &
                                      + dgamtot3dr(ikx,:)*save1(ikx,it) &
                                      + gamtot3(ikx,:)*save2(ikx,it)
                   end do
                end do
             else
                if (proc0) write (*,*) 'unknown dist option in radial_correction. aborting'
                call mp_abort ('unknown dist option in radial_correction. aborting')
             end if
          end if
       end if

       !collect quasineutrality corrections in wavenumber space
       do it = 1, ntubes
         do iz = -nzgrid, nzgrid
           g0k = phi(:,:,iz,it)
           call multiply_by_rho(g0k)
           phi_corr_QN(:,:,iz,it) = g0k
         enddo
       enddo
       do iky = 1, ky_solve_radial
         phi_corr_QN(iky,:,:,:) = 0.0
       enddo

       deallocate(phi)

       !collect gyroaveraging corrections in wavenumber space
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         is = is_idx(vmu_lo,ivmu)
         imu = imu_idx(vmu_lo,ivmu)
         do it = 1, ntubes
           do iz = -nzgrid, nzgrid
             call gyro_average_j1 (phi_in(:,:,iz,it), iz, ivmu, g0k)
             g0k = -g0k*(spec(is)%smz)**2 & 
                 * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                 * 0.5*(dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz))

             call multiply_by_rho(g0k)
             phi_corr_GA(:,:,iz,it,ivmu) = g0k
           enddo
         enddo
       enddo

       deallocate(g0x,g0k)
       deallocate (gyro_g)

    end if
    
  end subroutine get_radial_correction

  subroutine get_dchidy_4d (phi, apar, dchidy)

    use constants, only: zi
    use gyro_averages, only: gyro_average
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: is_idx, iv_idx
    use run_parameters, only: fphi, fapar
    use species, only: spec
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: vpa
    use kt_grids, only: nakx, aky, naky

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phi, apar
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (out) :: dchidy

    integer :: ivmu, iv, is
    complex, dimension (:,:,:,:), allocatable :: field

    allocate (field(naky,nakx,-nzgrid:nzgrid,ntubes))

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       iv = iv_idx(vmu_lo,ivmu)
       field = zi*spread(spread(spread(aky,2,nakx),3,2*nzgrid+1),4,ntubes) &
            * ( fphi*phi - fapar*vpa(iv)*spec(is)%stm*apar )
       call gyro_average (field, ivmu, dchidy(:,:,:,:,ivmu))
    end do

    deallocate (field)

  end subroutine get_dchidy_4d

  subroutine get_dchidy_2d (iz, ivmu, phi, apar, dchidy)

    use constants, only: zi
    use gyro_averages, only: gyro_average
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: is_idx, iv_idx
    use run_parameters, only: fphi, fapar
    use species, only: spec
    use vpamu_grids, only: vpa
    use kt_grids, only: nakx, aky, naky

    implicit none

    integer, intent (in) :: ivmu, iz
    complex, dimension (:,:), intent (in) :: phi, apar
    complex, dimension (:,:), intent (out) :: dchidy

    integer :: iv, is
    complex, dimension (:,:), allocatable :: field

    allocate (field(naky,nakx))

    is = is_idx(vmu_lo,ivmu)
    iv = iv_idx(vmu_lo,ivmu)
    field = zi*spread(aky,2,nakx) &
         * ( fphi*phi - fapar*vpa(iv)*spec(is)%stm*apar )
    call gyro_average (field, iz, ivmu, dchidy)
    
    deallocate (field)

  end subroutine get_dchidy_2d

  subroutine get_dchidx (iz, ivmu, phi, apar, dchidx)

    use constants, only: zi
    use gyro_averages, only: gyro_average
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: is_idx, iv_idx
    use run_parameters, only: fphi, fapar
    use species, only: spec
    use vpamu_grids, only: vpa
    use kt_grids, only: akx, naky, nakx

    implicit none

    integer, intent (in) :: ivmu, iz
    complex, dimension (:,:), intent (in) :: phi, apar
    complex, dimension (:,:), intent (out) :: dchidx

    integer :: iv, is
    complex, dimension (:,:), allocatable :: field

    allocate (field(naky,nakx))

    is = is_idx(vmu_lo,ivmu)
    iv = iv_idx(vmu_lo,ivmu)
    field = zi*spread(akx,1,naky) &
         * ( fphi*phi - fapar*vpa(iv)*spec(is)%stm*apar )
    call gyro_average (field, iz, ivmu, dchidx)
    
    deallocate (field)

  end subroutine get_dchidx

  subroutine finish_fields

    use fields_arrays, only: phi, phi_old
    use fields_arrays, only: phi_corr_QN, phi_corr_GA
    use fields_arrays, only: apar, apar_corr_QN, apar_corr_GA

    implicit none

    if (allocated(phi)) deallocate (phi)
    if (allocated(phi_old)) deallocate (phi_old)
    if (allocated(phi_corr_QN)) deallocate (phi_corr_QN)
    if (allocated(phi_corr_GA)) deallocate (phi_corr_GA)
    if (allocated(apar)) deallocate (apar)
    if (allocated(apar_corr_QN)) deallocate (apar_corr_QN)
    if (allocated(apar_corr_GA)) deallocate (apar_corr_GA)
    if (allocated(gamtot)) deallocate (gamtot)
    if (allocated(gamtot3)) deallocate (gamtot3)
    if (allocated(dgamtotdr)) deallocate (dgamtotdr)
    if (allocated(dgamtot3dr)) deallocate (dgamtot3dr)
    if (allocated(apar_denom)) deallocate (apar_denom)
    if (allocated(save1)) deallocate(save1)
    if (allocated(save2)) deallocate(save2)

    fields_initialized = .false.

  end subroutine finish_fields

end module fields
