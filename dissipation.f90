module dissipation

  implicit none

  public :: init_dissipation, finish_dissipation
  public :: include_collisions
  public :: include_krook_operator, update_delay_krook
  public :: remove_zero_projection, project_out_zero
  public :: advance_collisions_explicit, advance_collisions_implicit
  public :: time_collisions
  public :: hyper_dissipation
  public :: advance_hyper_dissipation
  public :: add_krook_operator
  public :: collisions_implicit
  public :: delay_krook, int_krook, int_proj

  private

  logical :: include_collisions, vpa_operator, mu_operator
  logical :: collisions_implicit, include_krook_operator
  logical :: momentum_conservation, energy_conservation
  logical :: hyper_dissipation, remove_zero_projection
  logical :: krook_odd
  real :: D_hyper, nu_krook, delay_krook, int_krook, int_proj
  integer:: ikxmax_source

  character(30) :: collision_model
  integer :: nresponse_vpa = 1
  integer :: nresponse_mu = 1
  real :: cfac

  real, dimension (:,:), allocatable :: aa_vpa, bb_vpa, cc_vpa
  real, dimension (:,:,:), allocatable :: aa_mu, cc_mu
  real, dimension (:,:), allocatable :: bb_mu
  complex, dimension (:,:,:), allocatable :: vpadiff_response
  integer, dimension (:,:), allocatable :: vpadiff_idx
  complex, dimension (:,:,:), allocatable :: mudiff_response
  integer, dimension (:,:), allocatable :: mudiff_idx

  complex, dimension (:,:,:,:,:), allocatable :: aa_blcs, cc_blcs
  complex, dimension (:,:,:,:), allocatable :: bb_blcs
  complex, dimension (:,:,:,:), allocatable :: dd_vecs
  complex, dimension (:,:,:), allocatable :: cdiffmat_band
  complex, dimension (:,:,:), allocatable :: blockmatrix
  integer, dimension (:,:), allocatable :: ipiv
  real, dimension (:,:,:), allocatable :: nus, nuD, nupa, nux, mw
  integer :: info

  logical :: collisions_initialized = .false.
  real, dimension (2,2) :: time_collisions = 0.

contains

  subroutine init_dissipation

    use mp, only: proc0
    use kt_grids, only: nakx
    use zgrid, only: ntubes
    use stella_layouts, only: vmu_lo
    use dist_fn_arrays, only: g_krook, g_proj

    implicit none

    call read_parameters
    if (include_collisions) then
        call init_collisions
    else
        if (proc0) then
           write (*,'(A)') "############################################################"
           write (*,'(A)') "                         COLLISIONS"
           write (*,'(A)') "############################################################"
           write (*,*) 'Coll. model:     None'
           write (*,*)
        end if
    end if

    if(include_krook_operator.and..not.allocated(g_krook)) then
      allocate (g_krook(nakx,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g_krook = 0.
    endif

    if(remove_zero_projection.and..not.allocated(g_proj)) then
      allocate (g_proj(nakx,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g_proj = 0.
    endif
    int_krook = 0.
    int_proj  = 0.


  end subroutine init_dissipation

  subroutine read_parameters

    use file_utils, only: input_unit_exist
    use mp, only: proc0, broadcast
    use kt_grids, only: ikx_max

    implicit none

    namelist /dissipation/ collision_model, hyper_dissipation, D_hyper, &
         include_collisions, collisions_implicit, &
         momentum_conservation, energy_conservation, &
         vpa_operator, mu_operator, include_krook_operator, &
         nu_krook, delay_krook, remove_zero_projection, &
         ikxmax_source, cfac, krook_odd

    integer :: in_file
    logical :: dexist

    if (proc0) then
       include_collisions = .false.
       include_krook_operator = .false.
       collisions_implicit = .true.
       collision_model = "dougherty"
       momentum_conservation = .true.
       energy_conservation = .true.
       vpa_operator = .true.
       mu_operator = .true.
       hyper_dissipation = .false.
       remove_zero_projection = .false.
       D_hyper = 0.05
       nu_krook = 0.05
       delay_krook =0.02
       ikxmax_source = 2 ! kx=0 and kx=1
       krook_odd = .true. ! damp only the odd mode that can affect profiles
       cfac = 1

       in_file = input_unit_exist("dissipation", dexist)
       if (dexist) read (unit=in_file, nml=dissipation)
    end if

    ikxmax_source = min(ikxmax_source,ikx_max)

    call broadcast (include_collisions)
    call broadcast (include_krook_operator)
    call broadcast (collisions_implicit)
    call broadcast (collision_model)
    call broadcast (momentum_conservation)
    call broadcast (energy_conservation)
    call broadcast (vpa_operator)
    call broadcast (mu_operator)
    call broadcast (hyper_dissipation)
    call broadcast (D_hyper)
    call broadcast (nu_krook)
    call broadcast (delay_krook)
    call broadcast (ikxmax_source)
    call broadcast (remove_zero_projection)
    call broadcast (cfac)
    call broadcast (krook_odd)

    if (.not.include_collisions) collisions_implicit = .false.

  end subroutine read_parameters

  subroutine init_collisions

    use species, only: spec, nspec
    use vpamu_grids, only: dvpa, dmu, mu, nmu
    use stella_geometry, only: bmag
    use stella_layouts

    implicit none

    integer :: is
    real :: cfl_dt_vpadiff, cfl_dt_mudiff
    real :: vnew_max


    if (collisions_initialized) return
    collisions_initialized = .true.

    if (collision_model == "dougherty") then
        write(*,*)
        write(*,*) 'Coll. model:     Dougherty'
        if (collisions_implicit) then
            write(*,*) 'Coll. algorithm: implicit'
           if (vpa_operator) then
              call init_vpadiff_matrix
              call init_vpadiff_conserve
           end if
           if (mu_operator) then
              call init_mudiff_matrix
              call init_mudiff_conserve
           end if
        else
            write(*,*) 'Coll. algorithm: explicit'
           vnew_max = 0.0
           do is = 1, nspec
              vnew_max = max(vnew_max,maxval(spec(is)%vnew))
           end do
           cfl_dt_vpadiff = 2.0*dvpa**2/vnew_max
           cfl_dt_mudiff = minval(bmag)/(vnew_max*maxval(mu(2:)/dmu(:nmu-1)**2))
        end if
    end if

    if (collision_model == "fokker-planck") then
        write(*,*) 'Coll. model:     linearized Fokker-Planck'
        write(*,*) 'Note:            tested for linear CBC w. adiabatic e-'
        call init_nusDpa
        if (collisions_implicit) then
            write(*,*) 'Coll. algorithm: implicit'
           call init_fp_diffmatrix
           call init_fp_conserve
        else
            write(*,*) 'Coll. algorithm: explicit'
           vnew_max = 0.0
           do is = 1, nspec
              vnew_max = max(vnew_max,maxval(spec(is)%vnew))
           end do
           cfl_dt_vpadiff = 2.0*dvpa**2/vnew_max
           cfl_dt_mudiff = minval(bmag)/(vnew_max*maxval(mu(2:)/dmu(:nmu-1)**2))
        end if
    end if
    write(*,*)
  end subroutine init_collisions

  subroutine init_nusDpa

          ! AVB: compute the collision frequencies nuD, nus and nupa

          use zgrid, only: nzgrid
          use vpamu_grids, only: nmu, mu, dmu, vpa, nvpa
          use stella_geometry, only: bmag
          use species, only: spec
          use spfunc, only: erf => erf_ext
          use finite_differences, only: fd3pt
          use vpamu_grids, only: maxwell_mu, maxwell_vpa
          use constants, only: pi

          implicit none

          integer :: ia, imu, iv, iz, is
          real :: x, Gf

          if (.not.allocated(nus)) allocate (nus(nvpa,nmu,-nzgrid:nzgrid))
          if (.not.allocated(nuD)) allocate (nuD(nvpa,nmu,-nzgrid:nzgrid))
          if (.not.allocated(nupa)) allocate (nupa(nvpa,nmu,-nzgrid:nzgrid))
          if (.not.allocated(nux)) allocate (nux(nvpa,nmu,-nzgrid:nzgrid))
          if (.not.allocated(mw)) allocate (mw(nvpa,nmu,-nzgrid:nzgrid))

          ia = 1
          is = 1

          do iz = -nzgrid, nzgrid
              do iv = 1, nvpa
                  do imu = 1, nmu
                      x = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mu(imu))
                      Gf = (erf(x)-x*(2/sqrt(pi))*exp(-x**2)) / (2*x**2)
                      nuD(iv,imu,iz) = spec(is)%vnew(is)*(erf(x)-Gf)/x**3
                      nus(iv,imu,iz) = spec(is)%vnew(is)*4*Gf/x
                      nupa(iv,imu,iz)= spec(is)%vnew(is)*2*Gf/x**3
                      mw(iv,imu,iz)  = maxwell_vpa(iv,is)*maxwell_mu(1,iz,imu,is)
                  end do
              end do
          end do

          nux = nupa - nuD

      end subroutine init_nusDpa

      subroutine finish_nusDpa

              implicit none

              if (allocated(nus)) deallocate (nus)
              if (allocated(nuD)) deallocate (nuD)
              if (allocated(nupa)) deallocate (nupa)
              if (allocated(nux)) deallocate (nux)
              if (allocated(mw)) deallocate (mw)

      end subroutine finish_nusDpa

      subroutine init_fp_diffmatrix

        use stella_time, only: code_dt
        use species, only: nspec, spec
        use zgrid, only: nzgrid
        use vpamu_grids, only: dvpa, vpa, nvpa, mu, nmu, maxwell_mu, maxwell_vpa, dmu, equally_spaced_mu_grid
        use stella_layouts, only: kxkyz_lo
        use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
        use stella_geometry, only: bmag
        use dist_fn_arrays, only: kperp2
        use constants, only: pi

        implicit none

        integer :: ikxkyz, iky, ikx, iz, is
        integer :: imu, ia, iv
        integer :: nc, nb, lldab
        real :: vpap, vpam, vfac, mum, mup
        real :: xpv, xmv, nupapv, nupamv, nuDpv, nuDmv, mwpv, mwmv, gam_mu, gam_mum, gam_mup
        real :: mwm, mwp, nuDm, nuDp, nupam, nupap, xm, xp

        if (.not.allocated(aa_blcs)) allocate (aa_blcs(nvpa,nmu,nmu,-nzgrid:nzgrid,nspec))
        if (.not.allocated(bb_blcs)) allocate (bb_blcs(nvpa,nmu,nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
        if (.not.allocated(cc_blcs)) allocate (cc_blcs(nvpa,nmu,nmu,-nzgrid:nzgrid,nspec))
        if (.not.allocated(dd_vecs)) allocate (dd_vecs(nvpa,nmu,-nzgrid:nzgrid,nspec))
        if (.not.allocated(cdiffmat_band)) allocate (cdiffmat_band(3*(nmu+1)+1,nmu*nvpa,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
        if (.not.allocated(blockmatrix)) allocate (blockmatrix(nvpa*nmu,nvpa*nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
        if (.not.allocated(ipiv)) allocate (ipiv(nvpa*nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

        ia = 1
        vfac = 1 ! zero vpa-operator, in beta
        aa_blcs = 0.
        bb_blcs = 0.
        cc_blcs = 0.
        dd_vecs = 0.

        ! AVB: construct difference matrix stored in block form
        ! AVB: this matrix is nvpa*nmu x nvpa*nmu
        ! AVB: aa_blcs stores subdiagonal blocks, bb_blcs diagonal blocks and cc_blcs the superdiagonal blocks
        ! AVB: aa_blcs(1,:,:) and cc_blcs(nmu,:,:) are never used

        do imu = 1, nmu
            do iv = 1, nvpa
                do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                   iky = iky_idx(kxkyz_lo,ikxkyz)
                   ikx = ikx_idx(kxkyz_lo,ikxkyz)
                   iz = iz_idx(kxkyz_lo,ikxkyz)
                   is = is_idx(kxkyz_lo,ikxkyz)

                   if (iv == 1) then
                        vpap   = 0.5*(vpa(iv)+vpa(iv+1))
                        mwpv   = exp(-vpap**2)*maxwell_mu(1,iz,imu,is)
                        xpv    = sqrt(vpap**2 + 2*bmag(ia,iz)*mu(imu))
                        nupapv = vfac*spec(is)%vnew(is)*2*(erf(xpv)-xpv*(2/sqrt(pi))*exp(-xpv**2)) / (2*xpv**2)/xpv**3
                        nuDpv  = vfac*spec(is)%vnew(is)*(erf(xpv)-(erf(xpv)-xpv*(2/sqrt(pi))*exp(-xpv**2)) / (2*xpv**2))/xpv**3

                        if (imu == 1) then
                            ! one-sided difference for mu-derivative at imu=1:
                            cc_blcs(iv,imu,imu+1,iz,is) = -vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*mw(iv+1,imu,iz)/mw(iv+1,imu+1,iz) / (2*dvpa) / dmu(imu)
                            cc_blcs(iv,imu,imu,iz,is)   = -code_dt*0.5*(nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz) &
                                                                +vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*mw(iv+1,imu,iz)/mw(iv+1,imu,iz) / (2*dvpa) / dmu(imu)
                            bb_blcs(iv,imu,imu,ikxkyz)  = 1 + code_dt*(0.5*(nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv) / dvpa**2 / mw(iv,imu,iz)
                            ! mu operator
                            if (mu_operator) then
                                ! use ghost cell at mu_{0} = 0, where term vanishes, so dmu(0) = mu(1).
                                mup = 0.5*(mu(imu)+mu(imu+1))
                                mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
                                xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
                                nuDp = spec(is)%vnew(is)*(erf(xp)-(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2)) / xp**3
                                nupap = spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
                                gam_mu = 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)
                                gam_mup = 2*(nupap*mup**2+nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp

                                bb_blcs(iv,imu,imu,ikxkyz)  = bb_blcs(iv,imu,imu,ikxkyz) - code_dt*(gam_mu*-1/dmu(imu) * dmu(imu)/2./mu(imu) &
                                                                +(gam_mup*-1/dmu(imu) - gam_mu*-1/dmu(imu)) * mu(imu)/(dmu(imu)/2.)) / mw(iv,imu,iz) / (dmu(imu)/2.+mu(imu))
                                bb_blcs(iv,imu,imu+1,ikxkyz)= bb_blcs(iv,imu,imu+1,ikxkyz) - code_dt*(gam_mu * 1/dmu(imu) * dmu(imu)/2./mu(imu) &
                                                                +(gam_mup*1/dmu(imu) - gam_mu*1/dmu(imu)) * mu(imu)/(dmu(imu)/2.)) / mw(iv,imu+1,iz) / (dmu(imu)/2.+mu(imu))
                                ! mixed derivative:
                                cc_blcs(iv,imu,imu,iz,is)   = cc_blcs(iv,imu,imu,iz,is) - code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv+1,imu,iz)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu))) / (mu(imu)+dmu(imu)) / dvpa
                                cc_blcs(iv,imu,imu+1,iz,is) = cc_blcs(iv,imu,imu+1,iz,is) - code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*1/mw(iv+1,imu+1,iz)* mu(imu)/dmu(imu)) / (mu(imu)+dmu(imu)) / dvpa
                                bb_blcs(iv,imu,imu,ikxkyz)  = bb_blcs(iv,imu,imu,ikxkyz) + code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv,imu,iz)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu))) / (mu(imu)+dmu(imu)) / dvpa
                                bb_blcs(iv,imu,imu+1,ikxkyz)= bb_blcs(iv,imu,imu+1,ikxkyz) + code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*1/mw(iv,imu+1,iz)* mu(imu)/dmu(imu)) / (mu(imu)+dmu(imu)) / dvpa
                            end if
                        else if (imu == nmu) then
                            ! AVB: one-sided difference for mu-derivative at imu=nmu:
                            cc_blcs(iv,imu,imu-1,iz,is) = +vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*mw(iv+1,imu,iz)/mw(iv+1,imu-1,iz) / (2*dvpa) / dmu(nmu-1)
                            cc_blcs(iv,imu,imu,iz,is)   = -code_dt*0.5*(nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz) &
                                                            -vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*mw(iv+1,imu,iz)/mw(iv+1,imu,iz) / (2*dvpa) / dmu(nmu-1)
                            bb_blcs(iv,imu,imu,ikxkyz)  = 1 + code_dt*(0.5*(nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv) / dvpa**2 / mw(iv,imu,iz)
                             ! mu operator
                            if (mu_operator) then
                                mum = 0.5*(mu(imu)+mu(imu-1))
                                mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
                                xm  = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
                                nuDm = spec(is)%vnew(is)*(erf(xm)-(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2)) / xm**3
                                nupam = spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
                                gam_mu = 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)
                                gam_mum = 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm

                                bb_blcs(iv,imu,imu,ikxkyz)  = bb_blcs(iv,imu,imu,ikxkyz) - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                                                + (-gam_mu/dmu(imu-1))*dmu(imu-1)/2./dmu(imu-1)) / mw(iv,imu,iz) / (dmu(imu-1)/2.+dmu(imu-1))
                                bb_blcs(iv,imu,imu-1,ikxkyz)= bb_blcs(iv,imu,imu-1,ikxkyz) - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * -1/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                                                -gam_mu * -1/dmu(imu-1)*dmu(imu-1)/2./dmu(imu-1)) / mw(iv,imu-1,iz) / (dmu(imu-1)/2.+dmu(imu-1))
                                ! mixed derivative:
                                cc_blcs(iv,imu,imu,iz,is)   = cc_blcs(iv,imu,imu,iz,is) - code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv+1,imu,iz)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                cc_blcs(iv,imu,imu-1,iz,is) = cc_blcs(iv,imu,imu-1,iz,is) + code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*1/mw(iv+1,imu-1,iz)* dmu(imu-1)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                bb_blcs(iv,imu,imu,ikxkyz)  = bb_blcs(iv,imu,imu,ikxkyz) + code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv,imu,iz)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                bb_blcs(iv,imu,imu-1,ikxkyz)= bb_blcs(iv,imu,imu-1,ikxkyz) - code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*1/mw(iv,imu-1,iz)* dmu(imu-1)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                            end if
                        else
                            cc_blcs(iv,imu,imu-1,iz,is) = +vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*mw(iv+1,imu,iz)/mw(iv+1,imu-1,iz)*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                            cc_blcs(iv,imu,imu+1,iz,is) = -vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*mw(iv+1,imu,iz)/mw(iv+1,imu+1,iz)*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                            cc_blcs(iv,imu,imu,iz,is)   = -code_dt*0.5*(nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz) &
                                                            -vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                            bb_blcs(iv,imu,imu,ikxkyz)  = 1 + code_dt*(0.5*(nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv) / dvpa**2 / mw(iv,imu,iz)
                            ! mu operator
                            if (mu_operator) then
                                ! quantities at mu_{i+1/2} and mu_{i-1/2}:
                                mup = 0.5*(mu(imu)+mu(imu+1))
                                mum = 0.5*(mu(imu)+mu(imu-1))
                                mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
                                mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
                                xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
                                xm = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
                                nuDp = spec(is)%vnew(is)*(erf(xp)-(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2)) / xp**3
                                nuDm = spec(is)%vnew(is)*(erf(xm)-(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2)) / xm**3
                                nupap = spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
                                nupam = spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
                                gam_mu = 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)
                                gam_mum = 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm
                                gam_mup = 2*(nupap*mup**2+nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp

                                bb_blcs(iv,imu,imu,ikxkyz)   = bb_blcs(iv,imu,imu,ikxkyz) - code_dt*((gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mum / dmu(imu-1))*dmu(imu)/dmu(imu-1) &
                                                                +(-gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mup / dmu(imu))* dmu(imu-1)/dmu(imu)) / mw(iv,imu,iz) * 2./(dmu(imu-1)+dmu(imu))
                                bb_blcs(iv,imu,imu-1,ikxkyz) = bb_blcs(iv,imu,imu-1,ikxkyz) - code_dt*((gam_mu * -1*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) - gam_mum* -1/dmu(imu-1)) * dmu(imu)/dmu(imu-1) &
                                                                -gam_mu * -1*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) * dmu(imu-1)/dmu(imu)) / mw(iv,imu-1,iz) * 2./(dmu(imu-1)+dmu(imu))
                                bb_blcs(iv,imu,imu+1,ikxkyz) = bb_blcs(iv,imu,imu+1,ikxkyz) - code_dt*( gam_mu*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) * dmu(imu)/dmu(imu-1) &
                                                                + (gam_mup/dmu(imu) - gam_mu*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu))) * dmu(imu-1)/dmu(imu)) / mw(iv,imu+1,iz) * 2/(dmu(imu-1)+dmu(imu))
                                ! mixed derivative:
                                cc_blcs(iv,imu,imu,iz,is)    = cc_blcs(iv,imu,imu,iz,is) - code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv+1,imu,iz)*(dmu(imu)/dmu(imu-1)-dmu(imu-1)/dmu(imu))) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                cc_blcs(iv,imu,imu-1,iz,is)  = cc_blcs(iv,imu,imu-1,iz,is) + code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*1/mw(iv+1,imu-1,iz)* dmu(imu)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                cc_blcs(iv,imu,imu+1,iz,is)  = cc_blcs(iv,imu,imu+1,iz,is) - code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*1/mw(iv+1,imu+1,iz)* dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                bb_blcs(iv,imu,imu,ikxkyz)   = bb_blcs(iv,imu,imu,ikxkyz) + code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv,imu,iz)*(dmu(imu)/dmu(imu-1)-dmu(imu-1)/dmu(imu))) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                bb_blcs(iv,imu,imu-1,ikxkyz) = bb_blcs(iv,imu,imu-1,ikxkyz) - code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*1/mw(iv,imu-1,iz)* dmu(imu)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                bb_blcs(iv,imu,imu+1,ikxkyz) = bb_blcs(iv,imu,imu+1,ikxkyz) + code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*1/mw(iv,imu+1,iz)* dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                            end if
                        end if

                    else if (iv == nvpa) then
                        vpam = 0.5*(vpa(iv)+vpa(iv-1))
                        mwmv = exp(-vpam**2)*maxwell_mu(1,iz,imu,is)
                        xmv = sqrt(vpam**2 + 2*bmag(ia,iz)*mu(imu))
                        nupamv = vfac*spec(is)%vnew(is)*2*(erf(xmv)-xmv*(2/sqrt(pi))*exp(-xmv**2)) / (2*xmv**2)/xmv**3
                        nuDmv = vfac*spec(is)%vnew(is)*(erf(xmv)-(erf(xmv)-xmv*(2/sqrt(pi))*exp(-xmv**2)) / (2*xmv**2))/xmv**3

                        if (imu == 1) then
                            ! one-sided difference for mu-derivative at imu=1:
                            aa_blcs(iv,imu,imu+1,iz,is) = +vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*mw(iv-1,imu,iz)/mw(iv-1,imu+1,iz) / (2*dvpa) / dmu(imu)
                            aa_blcs(iv,imu,imu,iz,is)   = -code_dt*0.5*(nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2  / mw(iv-1,imu,iz) &
                                                            -vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*mw(iv-1,imu,iz)/mw(iv-1,imu,iz) / (2*dvpa) / dmu(imu)
                            bb_blcs(iv,imu,imu,ikxkyz)  = 1 + code_dt*(0.5*(nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv) / dvpa**2  / mw(iv,imu,iz)
                            ! mu operator
                            if (mu_operator) then
                                ! use ghost cell at mu_{0} = 0, where term vanishes, so dmu(0) = mu(1).
                                mup = 0.5*(mu(imu)+mu(imu+1))
                                mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
                                xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
                                nuDp = spec(is)%vnew(is)*(erf(xp)-(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2)) / xp**3
                                nupap = spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
                                gam_mu = 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)
                                gam_mup = 2*(nupap*mup**2+nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp

                                bb_blcs(iv,imu,imu,ikxkyz)   = bb_blcs(iv,imu,imu,ikxkyz) - code_dt*(gam_mu *-1/dmu(imu) * dmu(imu)/2./mu(imu) &
                                                                +(gam_mup*-1/dmu(imu) - gam_mu*-1/dmu(imu)) * mu(imu)/(dmu(imu)/2.)) / mw(iv,imu,iz) / (dmu(imu)/2.+mu(imu))
                                bb_blcs(iv,imu,imu+1,ikxkyz) = bb_blcs(iv,imu,imu+1,ikxkyz) - code_dt*(gam_mu * 1/dmu(imu) * dmu(imu)/2./mu(imu) &
                                                                +(gam_mup* 1/dmu(imu) - gam_mu*1/dmu(imu)) * mu(imu)/(dmu(imu)/2.))  / mw(iv,imu+1,iz) / (dmu(imu)/2.+mu(imu))
                                ! mixed derivative:
                                aa_blcs(iv,imu,imu,iz,is)    = aa_blcs(iv,imu,imu,iz,is) + code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv-1,imu,iz)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu))) / (mu(imu)+dmu(imu)) / (dvpa)
                                aa_blcs(iv,imu,imu+1,iz,is)  = aa_blcs(iv,imu,imu+1,iz,is) + code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*1/mw(iv-1,imu+1,iz)* mu(imu)/dmu(imu)) / (mu(imu)+dmu(imu)) / (dvpa)
                                bb_blcs(iv,imu,imu,ikxkyz)   = bb_blcs(iv,imu,imu,ikxkyz) - code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv,imu,iz)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu))) / (mu(imu)+dmu(imu)) / (dvpa)
                                bb_blcs(iv,imu,imu+1,ikxkyz) = bb_blcs(iv,imu,imu+1,ikxkyz) - code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*1/mw(iv,imu+1,iz)* mu(imu)/dmu(imu)) / (mu(imu)+dmu(imu)) / (dvpa)
                            end if
                        else if (imu == nmu) then
                            ! one-sided difference for mu-derivative at imu=nmu:
                            aa_blcs(iv,imu,imu-1,iz,is) = -vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*mw(iv-1,imu,iz)/mw(iv-1,imu-1,iz) / (2*dvpa) / dmu(nmu-1)
                            aa_blcs(iv,imu,imu,iz,is)   = -code_dt*0.5*(nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2  / mw(iv-1,imu,iz) &
                                                            +vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*mw(iv-1,imu,iz)/mw(iv-1,imu,iz) / (2*dvpa) / dmu(nmu-1)
                            bb_blcs(iv,imu,imu,ikxkyz)  = 1 + code_dt*(0.5*(nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv) / dvpa**2 / mw(iv,imu,iz)
                            ! mu operator
                            if (mu_operator) then
                                mum = 0.5*(mu(imu)+mu(imu-1))
                                mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
                                xm = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
                                nuDm = spec(is)%vnew(is)*(erf(xm)-(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2)) / xm**3
                                nupam = spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
                                gam_mu  = 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)
                                gam_mum = 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm

                                bb_blcs(iv,imu,imu,ikxkyz)   = bb_blcs(iv,imu,imu,ikxkyz) - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                                                + (-gam_mu/dmu(imu-1))*dmu(imu-1)/2./dmu(imu-1)) / mw(iv,imu,iz) / (dmu(imu-1)/2.+dmu(imu-1))
                                bb_blcs(iv,imu,imu-1,ikxkyz) = bb_blcs(iv,imu,imu-1,ikxkyz) - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * -1/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                                                -gam_mu * -1/dmu(imu-1)*dmu(imu-1)/2./dmu(imu-1)) / mw(iv,imu-1,iz) / (dmu(imu-1)/2.+dmu(imu-1))
                                ! mixed derivative:
                                aa_blcs(iv,imu,imu,iz,is)    = aa_blcs(iv,imu,imu,iz,is) + code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv-1,imu,iz)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                aa_blcs(iv,imu,imu-1,iz,is)  = aa_blcs(iv,imu,imu-1,iz,is) - code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*1/mw(iv-1,imu-1,iz)* dmu(imu-1)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                bb_blcs(iv,imu,imu,ikxkyz)   = bb_blcs(iv,imu,imu,ikxkyz) - code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv,imu,iz)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                bb_blcs(iv,imu,imu-1,ikxkyz) = bb_blcs(iv,imu,imu-1,ikxkyz) + code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*1/mw(iv,imu-1,iz)* dmu(imu-1)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                            end if
                        else
                            aa_blcs(iv,imu,imu-1,iz,is)  = -vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*mw(iv-1,imu,iz)/mw(iv-1,imu-1,iz)* dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                            aa_blcs(iv,imu,imu+1,iz,is)  = +vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*mw(iv-1,imu,iz)/mw(iv-1,imu+1,iz)* dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                            aa_blcs(iv,imu,imu,iz,is)    = -code_dt*0.5*(nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz) &
                                                            +vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                            bb_blcs(iv,imu,imu,ikxkyz)   = 1 + code_dt*(0.5*(nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv) / dvpa**2 / mw(iv,imu,iz)
                            ! mu operator
                            if (mu_operator) then
                                ! quantities at mu_{i+1/2} and mu_{i-1/2}:
                                mup = 0.5*(mu(imu)+mu(imu+1))
                                mum = 0.5*(mu(imu)+mu(imu-1))
                                mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
                                mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
                                xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
                                xm = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
                                nuDp = spec(is)%vnew(is)*(erf(xp)-(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2)) / xp**3
                                nuDm = spec(is)%vnew(is)*(erf(xm)-(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2)) / xm**3
                                nupap = spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
                                nupam = spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
                                gam_mu = 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)
                                gam_mum = 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm
                                gam_mup = 2*(nupap*mup**2+nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp

                                bb_blcs(iv,imu,imu,ikxkyz)   = bb_blcs(iv,imu,imu,ikxkyz) - code_dt*((gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mum / dmu(imu-1))*dmu(imu)/dmu(imu-1) &
                                                                +(-gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mup / dmu(imu))* dmu(imu-1)/dmu(imu)) / mw(iv,imu,iz) * 2./(dmu(imu-1)+dmu(imu))
                                bb_blcs(iv,imu,imu-1,ikxkyz) = bb_blcs(iv,imu,imu-1,ikxkyz) - code_dt*((gam_mu * -1*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) - gam_mum* -1/dmu(imu-1)) * dmu(imu)/dmu(imu-1) &
                                                                -gam_mu * -1*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) * dmu(imu-1)/dmu(imu)) / mw(iv,imu-1,iz) * 2./(dmu(imu-1)+dmu(imu))
                                bb_blcs(iv,imu,imu+1,ikxkyz) = bb_blcs(iv,imu,imu+1,ikxkyz) - code_dt*( gam_mu*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) * dmu(imu)/dmu(imu-1) &
                                                                + (gam_mup/dmu(imu) - gam_mu*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu))) * dmu(imu-1)/dmu(imu)) / mw(iv,imu+1,iz) * 2/(dmu(imu-1)+dmu(imu))
                                ! mixed derivative, one-sided difference in vpa at iv = nvpa:
                                aa_blcs(iv,imu,imu,iz,is)    = aa_blcs(iv,imu,imu,iz,is) + code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv-1,imu,iz)*(dmu(imu)/dmu(imu-1)-dmu(imu-1)/dmu(imu))) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                aa_blcs(iv,imu,imu-1,iz,is)  = aa_blcs(iv,imu,imu-1,iz,is) - code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*1/mw(iv-1,imu-1,iz)* dmu(imu)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                aa_blcs(iv,imu,imu+1,iz,is)  = aa_blcs(iv,imu,imu+1,iz,is) + code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*1/mw(iv-1,imu+1,iz)* dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                bb_blcs(iv,imu,imu,ikxkyz)   = bb_blcs(iv,imu,imu,ikxkyz) - code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv,imu,iz)*(dmu(imu)/dmu(imu-1)-dmu(imu-1)/dmu(imu))) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                bb_blcs(iv,imu,imu-1,ikxkyz) = bb_blcs(iv,imu,imu-1,ikxkyz) + code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*1/mw(iv,imu-1,iz)* dmu(imu)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                bb_blcs(iv,imu,imu+1,ikxkyz) = bb_blcs(iv,imu,imu+1,ikxkyz) - code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*1/mw(iv,imu+1,iz)* dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                            end if
                        end if
                    else
                        vpam = 0.5*(vpa(iv)+vpa(iv-1))
                        vpap = 0.5*(vpa(iv)+vpa(iv+1))
                        mwmv = exp(-vpam**2)*maxwell_mu(1,iz,imu,is)
                        mwpv = exp(-vpap**2)*maxwell_mu(1,iz,imu,is)
                        xpv  = sqrt(vpap**2 + 2*bmag(ia,iz)*mu(imu))
                        xmv = sqrt(vpam**2 + 2*bmag(ia,iz)*mu(imu))
                        nupamv = vfac*spec(is)%vnew(is)*2*(erf(xmv)-xmv*(2/sqrt(pi))*exp(-xmv**2)) / (2*xmv**2)/xmv**3
                        nupapv = vfac*spec(is)%vnew(is)*2*(erf(xpv)-xpv*(2/sqrt(pi))*exp(-xpv**2)) / (2*xpv**2)/xpv**3
                        nuDmv = vfac*spec(is)%vnew(is)*(erf(xmv)-(erf(xmv)-xmv*(2/sqrt(pi))*exp(-xmv**2)) / (2*xmv**2))/xmv**3
                        nuDpv = vfac*spec(is)%vnew(is)*(erf(xpv)-(erf(xpv)-xpv*(2/sqrt(pi))*exp(-xpv**2)) / (2*xpv**2))/xpv**3

                        if (imu == 1) then
                            ! one-sided difference for mu-derivative at imu=1:
                             aa_blcs(iv,imu,imu+1,iz,is) = +vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*mw(iv-1,imu,iz)/mw(iv-1,imu+1,iz) / (2*dvpa) / dmu(imu)
                             aa_blcs(iv,imu,imu,iz,is)   = -code_dt*0.5*(nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz) -vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz) / (2*dvpa) / dmu(imu)
                             cc_blcs(iv,imu,imu+1,iz,is) = -vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*mw(iv+1,imu,iz)/mw(iv+1,imu+1,iz) / (2*dvpa) / dmu(imu)
                             cc_blcs(iv,imu,imu,iz,is)   = -code_dt*0.5*(nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz) +vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz) / (2*dvpa) / dmu(imu)
                             bb_blcs(iv,imu,imu,ikxkyz)  = 1 + code_dt*(0.5*(nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv + 0.5*(nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv) / dvpa**2 / mw(iv,imu,iz)
                             ! mu operator
                             if (mu_operator) then
                                 ! use ghost cell at mu_{0} = 0, where term vanishes, so dmu(0) = mu(1).
                                 mup = 0.5*(mu(imu)+mu(imu+1))
                                 mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
                                 xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
                                 nuDp = spec(is)%vnew(is)*(erf(xp)-(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2)) / xp**3
                                 nupap = spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
                                 gam_mu = 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)
                                 gam_mup = 2*(nupap*mup**2+nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp

                                 bb_blcs(iv,imu,imu,ikxkyz)   = bb_blcs(iv,imu,imu,ikxkyz) - code_dt*(gam_mu *-1/dmu(imu) * dmu(imu)/2./mu(imu) &
                                                                +(gam_mup*-1/dmu(imu) - gam_mu*-1/dmu(imu)) * mu(imu)/(dmu(imu)/2.)) / mw(iv,imu,iz) / (dmu(imu)/2.+mu(imu))
                                 bb_blcs(iv,imu,imu+1,ikxkyz) = bb_blcs(iv,imu,imu+1,ikxkyz)  - code_dt*(gam_mu * 1/dmu(imu) * dmu(imu)/2./mu(imu) &
                                                                +(gam_mup* 1/dmu(imu) - gam_mu*1/dmu(imu)) * mu(imu)/(dmu(imu)/2.))  / mw(iv,imu+1,iz) / (dmu(imu)/2.+mu(imu))
                                 ! mixed derivative:
                                 aa_blcs(iv,imu,imu,iz,is)    = aa_blcs(iv,imu,imu,iz,is) + code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv-1,imu,iz)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu))) / (mu(imu)+dmu(imu)) / (2*dvpa)
                                 aa_blcs(iv,imu,imu+1,iz,is)  = aa_blcs(iv,imu,imu+1,iz,is) + code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*1/mw(iv-1,imu+1,iz)* mu(imu)/dmu(imu)) / (mu(imu)+dmu(imu)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu,iz,is)    = cc_blcs(iv,imu,imu,iz,is) - code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv+1,imu,iz)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu))) / (mu(imu)+dmu(imu)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu+1,iz,is)  = cc_blcs(iv,imu,imu+1,iz,is) - code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*1/mw(iv+1,imu+1,iz)* mu(imu)/dmu(imu)) / (mu(imu)+dmu(imu)) / (2*dvpa)
                             end if
                        else if (imu == nmu) then
                             ! one-sided difference for mu-derivative at imu=nmu:
                             aa_blcs(iv,imu,imu-1,iz,is) = -vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*mw(iv-1,imu,iz)/mw(iv-1,imu-1,iz) / (2*dvpa) / dmu(nmu-1)
                             aa_blcs(iv,imu,imu,iz,is)   = -code_dt*0.5*(nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz) + vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz) / (2*dvpa) / dmu(nmu-1)
                             cc_blcs(iv,imu,imu-1,iz,is) = +vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*mw(iv+1,imu,iz)/mw(iv+1,imu-1,iz) / (2*dvpa) / dmu(nmu-1)
                             cc_blcs(iv,imu,imu,iz,is)   = -code_dt*0.5*(nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz) - vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz) / (2*dvpa) / dmu(nmu-1)
                             bb_blcs(iv,imu,imu,ikxkyz)  = 1 + code_dt*(0.5*(nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv + 0.5*(nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv) / dvpa**2 / mw(iv,imu,iz)
                             ! mu operator
                             if (mu_operator) then
                                 mum = 0.5*(mu(imu)+mu(imu-1))
                                 mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
                                 xm = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
                                 nuDm = spec(is)%vnew(is)*(erf(xm)-(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2)) / xm**3
                                 nupam = spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
                                 gam_mu = 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)
                                 gam_mum = 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm

                                 bb_blcs(iv,imu,imu,ikxkyz)   = bb_blcs(iv,imu,imu,ikxkyz) - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                                                + (-gam_mu/dmu(imu-1))*dmu(imu-1)/2./dmu(imu-1)) / mw(iv,imu,iz) / (dmu(imu-1)/2.+dmu(imu-1))
                                 bb_blcs(iv,imu,imu-1,ikxkyz) = bb_blcs(iv,imu,imu-1,ikxkyz) - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * -1/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                                                -gam_mu * -1/dmu(imu-1)*dmu(imu-1)/2./dmu(imu-1)) / mw(iv,imu-1,iz) / (dmu(imu-1)/2.+dmu(imu-1))
                                 ! mixed derivative:
                                 aa_blcs(iv,imu,imu,iz,is)    = aa_blcs(iv,imu,imu,iz,is) + code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv-1,imu,iz)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) / (dmu(imu-1)+dmu(imu-1)) / (2*dvpa)
                                 aa_blcs(iv,imu,imu-1,iz,is)  = aa_blcs(iv,imu,imu-1,iz,is) - code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*1/mw(iv-1,imu-1,iz)* dmu(imu-1)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu-1)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu,iz,is)    = cc_blcs(iv,imu,imu,iz,is) - code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv+1,imu,iz)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) / (dmu(imu-1)+dmu(imu-1)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu-1,iz,is)  = cc_blcs(iv,imu,imu-1,iz,is) + code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*1/mw(iv+1,imu-1,iz)* dmu(imu-1)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu-1)) / (2*dvpa)
                             end if
                        else
                            ! vpa operator (interior treatment):
                            aa_blcs(iv,imu,imu,iz,is)   = -code_dt*0.5*(nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz) &
                                                            +vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz) * (dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                            cc_blcs(iv,imu,imu,iz,is)   = -code_dt*0.5*(nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz) &
                                                            -vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz) * (dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                            bb_blcs(iv,imu,imu,ikxkyz)  = 1 + code_dt*(0.5*(nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv + 0.5*(nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv) / dvpa**2 / mw(iv,imu,iz)
                            ! vpa operator, mixed (interior treatment):
                            aa_blcs(iv,imu,imu-1,iz,is) = -vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*mw(iv-1,imu,iz)/mw(iv-1,imu-1,iz)*  dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                            aa_blcs(iv,imu,imu+1,iz,is) = +vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*mw(iv-1,imu,iz)/mw(iv-1,imu+1,iz)*  dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                            cc_blcs(iv,imu,imu-1,iz,is) = +vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*mw(iv+1,imu,iz)/mw(iv+1,imu-1,iz)*  dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                            cc_blcs(iv,imu,imu+1,iz,is) = -vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*mw(iv+1,imu,iz)/mw(iv+1,imu+1,iz)*  dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)

                            ! mu operator
                            if (mu_operator) then
                                ! quantities at mu_{i+1/2} and mu_{i-1/2}:
                                mup = 0.5*(mu(imu)+mu(imu+1))
                                mum = 0.5*(mu(imu)+mu(imu-1))
                                mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
                                mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
                                xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
                                xm = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
                                nuDp = spec(is)%vnew(is)*(erf(xp)-(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2)) / xp**3
                                nuDm = spec(is)%vnew(is)*(erf(xm)-(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2)) / xm**3
                                nupap = spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
                                nupam = spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
                                gam_mu = 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)
                                gam_mum = 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm
                                gam_mup = 2*(nupap*mup**2+nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp
                                ! mu_operator (interior treatment):
                                bb_blcs(iv,imu,imu,ikxkyz)   = bb_blcs(iv,imu,imu,ikxkyz) - code_dt*((gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mum / dmu(imu-1))*dmu(imu)/dmu(imu-1) &
                                                                +(-gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mup / dmu(imu))* dmu(imu-1)/dmu(imu)) / mw(iv,imu,iz) * 2./(dmu(imu-1)+dmu(imu))
                                bb_blcs(iv,imu,imu-1,ikxkyz) = bb_blcs(iv,imu,imu-1,ikxkyz) - code_dt*((gam_mu * -1*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) - gam_mum* -1/dmu(imu-1)) * dmu(imu)/dmu(imu-1) &
                                                                -gam_mu * -1*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) * dmu(imu-1)/dmu(imu)) / mw(iv,imu-1,iz) * 2./(dmu(imu-1)+dmu(imu))
                                bb_blcs(iv,imu,imu+1,ikxkyz) = bb_blcs(iv,imu,imu+1,ikxkyz) - code_dt*( gam_mu * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) * dmu(imu)/dmu(imu-1) &
                                                                + (gam_mup/dmu(imu) - gam_mu*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu))) * dmu(imu-1)/dmu(imu)) / mw(iv,imu+1,iz) * 2/(dmu(imu-1)+dmu(imu))
                                ! mu operator, mixed (interior treatment):
                                aa_blcs(iv,imu,imu,iz,is)    = aa_blcs(iv,imu,imu,iz,is)   + code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv-1,imu,iz)*(dmu(imu)/dmu(imu-1)-dmu(imu-1)/dmu(imu))) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                aa_blcs(iv,imu,imu-1,iz,is)  = aa_blcs(iv,imu,imu-1,iz,is) - code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*1/mw(iv-1,imu-1,iz)* dmu(imu)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                aa_blcs(iv,imu,imu+1,iz,is)  = aa_blcs(iv,imu,imu+1,iz,is) + code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*1/mw(iv-1,imu+1,iz)* dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                cc_blcs(iv,imu,imu,iz,is)    = cc_blcs(iv,imu,imu,iz,is)   - code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*1/mw(iv+1,imu,iz)*(dmu(imu)/dmu(imu-1)-dmu(imu-1)/dmu(imu))) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                cc_blcs(iv,imu,imu-1,iz,is)  = cc_blcs(iv,imu,imu-1,iz,is) + code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*1/mw(iv+1,imu-1,iz)* dmu(imu)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                cc_blcs(iv,imu,imu+1,iz,is)  = cc_blcs(iv,imu,imu+1,iz,is) - code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*1/mw(iv+1,imu+1,iz)* dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                           end if

                        end if
                    end if

                    ! add gyro-diffusive term:
                    bb_blcs(iv,imu,imu,ikxkyz) = bb_blcs(iv,imu,imu,ikxkyz) + code_dt*cfac*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*(nupa(iv,imu,iz)*bmag(ia,iz)*mu(imu) + nuD(iv,imu,iz)*(vpa(iv)**2 + bmag(ia,iz)*mu(imu)))

                end do
            end do
        end do

        ! switch from block storage to full matrix
        blockmatrix = 0.
        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz  = iz_idx(kxkyz_lo,ikxkyz)
            is  = is_idx(kxkyz_lo,ikxkyz)
            do iv = 1, nvpa
                ! diagonal blocks:
                blockmatrix(nmu*(iv-1)+1:nmu*iv, nmu*(iv-1)+1:nmu*iv, ikxkyz) = bb_blcs(iv,:,:,ikxkyz)
                if (iv < nvpa) then
                    ! subdiagonal blocks:
                    blockmatrix(nmu*iv+1:nmu*(iv+1), nmu*(iv-1)+1:nmu*iv, ikxkyz) = aa_blcs(iv+1,:,:,iz,is)
                    ! superdiagonal blocks:
                    blockmatrix(nmu*(iv-1)+1:nmu*iv, nmu*iv+1:nmu*(iv+1), ikxkyz) = cc_blcs(iv,:,:,iz,is)
                end if
            end do
        end do

        ! switch from full matrix to band-storage, for LAPACK banded solver routines:
        ! a_ij is stored in aband(ku+1+i-j,j) for $\max(1,j-ku) \leq i \leq \min(m,j+kl)$
        cdiffmat_band = 0.
        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz  = iz_idx(kxkyz_lo,ikxkyz)
            is  = is_idx(kxkyz_lo,ikxkyz)
            do iv = 1, nmu*nvpa
                do imu = 1, nmu*nvpa
                    if ((max(1, iv-(nmu+1)) .le. imu) .and. (imu .le. min(nvpa*nmu, iv+(nmu+1)))) then
                        cdiffmat_band(nmu+1 + nmu+1 + 1+imu-iv, iv, ikxkyz) = blockmatrix(imu, iv, ikxkyz)
                    end if
                end do
            end do
        end do

        ! AVB: LU factorise cdiffmat, using LAPACK's zgbtrf routine for factorising banded matrices:
        nc = nvpa*nmu
        nb = nmu+1
        lldab = 3*(nmu+1)+1
        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            call zgbtrf(nc, nc, nb, nb, cdiffmat_band(:,:,ikxkyz), lldab, ipiv(:,ikxkyz), info)
        end do

    end subroutine init_fp_diffmatrix

    subroutine init_fp_conserve

      use finite_differences, only: tridag
      use linear_solve, only: lu_decomposition
      use stella_time, only: code_dt
      use species, only: nspec, spec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: ztmax, maxwell_vpa, maxwell_mu
      use vpamu_grids, only: nmu, nvpa,vpa, vperp2
      use vpamu_grids, only: set_vpa_weights, wgts_vpa, wgts_mu
      use kt_grids, only: naky, nakx
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, it_idx
      use dist_fn_arrays, only: gvmu
      use gyro_averages, only: aj0v
      use fields, only: get_fields, get_fields_by_spec
      use stella_geometry, only: bmag
      use job_manage, only: time_message, timer_local

      implicit none

      integer :: ikxkyz, iky, ikx, iz, is, it
      integer :: imu, idx, ix, iv
      logical :: conservative_wgts
      real :: dum2, dum3
      complex, dimension (:,:,:,:), allocatable :: dum1
      complex, dimension (:,:,:,:,:), allocatable :: field
      complex, dimension (:,:), allocatable :: sumdelta
      complex, dimension (:,:), allocatable :: gvmutr
      complex, dimension (:), allocatable :: ghrs

      if (.not.allocated(vpadiff_response)) then
         nresponse_vpa = 1
         if (momentum_conservation) nresponse_vpa = nresponse_vpa + nspec
         if (energy_conservation) nresponse_vpa = nresponse_vpa + nspec
         allocate (vpadiff_response(nresponse_vpa,nresponse_vpa,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
         vpadiff_response = 0.
         allocate (vpadiff_idx(nresponse_vpa,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      end if

      allocate (dum1(naky,nakx,-nzgrid:nzgrid,ntubes))
      allocate (field(naky,nakx,-nzgrid:nzgrid,ntubes,nspec))
      allocate (sumdelta(nmu,nvpa)); sumdelta = 0.

      allocate (gvmutr(nvpa,nmu))
      allocate (ghrs(nmu*nvpa))

      ! set wgts to be equally spaced to ensure exact conservation properties
      conservative_wgts = .true.
      call set_vpa_weights (conservative_wgts)

      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo,ikxkyz)
         ikx = ikx_idx(kxkyz_lo,ikxkyz)
         iz = iz_idx(kxkyz_lo,ikxkyz)
         is = is_idx(kxkyz_lo,ikxkyz)
         ! AVB: calculate Green's function: supply unit impulse to rhs at location (iv,imu), solve for response
         sumdelta = 0.
         do iv = 1, nvpa
            do imu = 1, nmu
                ghrs = 0.
                ghrs(nmu*(iv-1)+imu) = 1.
                ! solve for response:
                call zgbtrs('No transpose', nvpa*nmu, nmu+1, nmu+1, 1, cdiffmat_band(:,:,ikxkyz), 3*(nmu+1)+1, ipiv(:,ikxkyz), ghrs, nvpa*nmu, info)
                do ix = 1, nvpa
                    gvmutr(ix,:) = ghrs(nmu*(ix-1)+1 : nmu*(ix-1)+nmu)
                end do
                sumdelta = sumdelta + ztmax(iv,is)*maxwell_mu(1,iz,imu,is)*aj0v(imu,ikxkyz)*transpose(gvmutr)
            end do
        end do
        gvmu(:,:,ikxkyz) = transpose(sumdelta) ! AVB: gvmu(:,:,ikxkyz) now contains: sum_iv sum_imu [ response(vpar,mu)_{iv,imu} ] / phi^(n+1) = h_{hom}^n+1(vpa,mu)_ikxkyz / phi^{n+1}_ikxkyz
      end do

      ! gvmu contains dhs/dphi
      ! for phi equation, need 1-P[dhs/dphi]
      call get_fields (gvmu, field(:,:,:,:,1), dum1, dist='h')

      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo,ikxkyz)
         ikx = ikx_idx(kxkyz_lo,ikxkyz)
         iz = iz_idx(kxkyz_lo,ikxkyz)
         it = it_idx(kxkyz_lo,ikxkyz)
         vpadiff_response(1,1,ikxkyz) = 1.0-field(iky,ikx,iz,it,1) ! AVB: ravelled response
      end do

      ! now get LU decomposition for vpadiff_response
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         call lu_decomposition (vpadiff_response(:,:,ikxkyz),vpadiff_idx(:,ikxkyz),dum2)
      end do

      ! reset wgts to default setting
      conservative_wgts = .false.
      call set_vpa_weights (conservative_wgts)

      deallocate (dum1, field)

  end subroutine init_fp_conserve

  subroutine init_vpadiff_matrix

    use stella_time, only: code_dt
    use species, only: nspec, spec
    use vpamu_grids, only: dvpa, vpa, nvpa
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
    use stella_geometry, only: bmag
    use dist_fn_arrays, only: kperp2

    implicit none

    integer :: ikxkyz, iky, ikx, iz, is
    integer :: ia

    if (.not.allocated(aa_vpa)) allocate (aa_vpa(nvpa,nspec))
!    if (.not.allocated(bb_vpa)) allocate (bb_vpa(nvpa,nspec))
    if (.not.allocated(bb_vpa)) allocate (bb_vpa(nvpa,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    if (.not.allocated(cc_vpa)) allocate (cc_vpa(nvpa,nspec))

    ! deal with boundary points (BC is f(vpa)=0 beyond +/- vpa_max)
    aa_vpa(1,:) = 0.0 ; cc_vpa(nvpa,:) = 0.0
    ! 2nd order centered differences for d/dvpa (1/2 dh/dvpa + vpa h)
    do is = 1, nspec
       aa_vpa(2:,is) = -code_dt*spec(is)%vnew(is)*0.5*(1.0/dvpa-vpa(:nvpa-1))/dvpa
!       bb_vpa(:,is) = 1.0+code_dt*spec(is)%vnew(is)/dvpa**2
       cc_vpa(:nvpa-1,is) = -code_dt*spec(is)%vnew(is)*0.5*(1.0/dvpa+vpa(2:))/dvpa
    end do

    ia = 1
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       bb_vpa(:,ikxkyz) = 1.0 + code_dt*spec(is)%vnew(is) &
          *  (0.25*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2 + 1./dvpa**2)
    end do

  end subroutine init_vpadiff_matrix

  subroutine init_mudiff_matrix

    use stella_time, only: code_dt
    use species, only: nspec, spec
    use zgrid, only: nzgrid
    use stella_geometry, only: bmag
    use vpamu_grids, only: dmu, mu, nmu
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
    use dist_fn_arrays, only: kperp2

    implicit none

    integer :: ikxkyz, iky, ikx, iz, is
    integer :: ia
    ! TMP FOR TESTING -- MAB
!    integer :: imu

    real, dimension (:), allocatable :: dmu_ghost, dmu_cell, mu_cell

    ! add ghost cell at mu=0 and beyond mu_max for purposes of differentiation
    ! note assuming here that grid spacing for ghost cell is equal to
    ! grid spacing for last non-ghost cell
    allocate (dmu_ghost(nmu))
    dmu_ghost(:nmu-1) = dmu ; dmu_ghost(nmu) = dmu(nmu-1)
    ! this is mu at cell centres (including to left and right of mu grid boundary points)
    allocate (mu_cell(nmu))
    mu_cell(:nmu-1) = 0.5*(mu(:nmu-1)+mu(2:))
    mu_cell(nmu) = mu(nmu)+0.5*dmu(nmu-1)
    ! this is mu_{j+1/2} - mu_{j-1/2}
    allocate (dmu_cell(nmu))
    dmu_cell(1) = mu_cell(1)
    dmu_cell(2:) = mu_cell(2:)-mu_cell(:nmu-1)

    if (.not.allocated(aa_mu)) allocate (aa_mu(-nzgrid:nzgrid,nmu,nspec))
    if (.not.allocated(bb_mu)) allocate (bb_mu(nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    if (.not.allocated(cc_mu)) allocate (cc_mu(-nzgrid:nzgrid,nmu,nspec))

    ia = 1

    ! deal with boundary points (BC is f(mu)=0 beyond mu_max and collision operator vanishes for mu -> 0)
    aa_mu(:,1,:) = 0.0 ; cc_mu(:,nmu,:) = 0.0
    ! 2nd order centered differences for dt * nu * d/dmu (mu/B*dh/dmu + 2*mu*h)
    do is = 1, nspec
       do iz = -nzgrid, nzgrid
          aa_mu(iz,2:,is) = -code_dt*spec(is)%vnew(is)*mu_cell(:nmu-1)*(1.0/(bmag(ia,iz)*dmu)-1.0)/dmu_cell(2:)
          cc_mu(iz,:nmu-1,is) = -code_dt*spec(is)%vnew(is)*mu_cell(:nmu-1)*(1.0/(bmag(ia,iz)*dmu)+1.0)/dmu_cell(:nmu-1)
       end do
    end do

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       bb_mu(1,ikxkyz) = 1.0 + code_dt*spec(is)%vnew(is) &
          * (0.25*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2 &
            + 1.0/(dmu(1)*bmag(ia,iz)) - 1.0)
       bb_mu(2:nmu-1,ikxkyz) = 1.0 + code_dt*spec(is)%vnew(is) &
          * (0.25*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2 &
            + (mu_cell(2:nmu-1)/dmu(2:)+mu_cell(:nmu-2)/dmu(:nmu-2)) &
            /(dmu_cell(2:nmu-1)*bmag(ia,iz)) - 1.0)
       bb_mu(nmu,ikxkyz) = 1.0 + code_dt*spec(is)%vnew(is) &
          * (0.25*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2 &
            + mu_cell(nmu-1)*(1.0/(dmu(nmu-1)*bmag(ia,iz)) + 1.0)/dmu_cell(nmu))
    end do

    deallocate (dmu_ghost, dmu_cell, mu_cell)

  end subroutine init_mudiff_matrix

  subroutine init_vpadiff_conserve

    use finite_differences, only: tridag
    use linear_solve, only: lu_decomposition
    use stella_time, only: code_dt
    use species, only: nspec, spec
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: ztmax, maxwell_vpa, maxwell_mu
    use vpamu_grids, only: nmu, vpa, vperp2
    use vpamu_grids, only: set_vpa_weights
    use kt_grids, only: naky, nakx
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use dist_fn_arrays, only: gvmu
    use gyro_averages, only: aj0v
    use fields, only: get_fields, get_fields_by_spec

    implicit none

    integer :: ikxkyz, iky, ikx, iz, it, is
    integer :: imu
    integer :: idx
    logical :: conservative_wgts
    real :: dum2
    complex, dimension (:,:,:,:), allocatable :: dum1
    complex, dimension (:,:,:,:,:), allocatable :: field

    if (.not.allocated(vpadiff_response)) then
       nresponse_vpa = 1
       if (momentum_conservation) nresponse_vpa = nresponse_vpa + nspec
       if (energy_conservation) nresponse_vpa = nresponse_vpa + nspec
       allocate (vpadiff_response(nresponse_vpa,nresponse_vpa,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       vpadiff_response = 0.
       allocate (vpadiff_idx(nresponse_vpa,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    end if
    allocate (dum1(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (field(naky,nakx,-nzgrid:nzgrid,ntubes,nspec))

    ! set wgts to be equally spaced to ensure exact conservation properties
    conservative_wgts = .true.
    call set_vpa_weights (conservative_wgts)

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          gvmu(:,imu,ikxkyz) = ztmax(:,is)*maxwell_mu(1,iz,imu,is)*aj0v(imu,ikxkyz)
          call tridag (1, aa_vpa(:,is), bb_vpa(:,ikxkyz), cc_vpa(:,is), gvmu(:,imu,ikxkyz))
       end do
    end do

    ! gvmu contains dhs/dphi
    ! for phi equation, need 1-P[dhs/dphi]
    ! for upar equations, need -Us[dhs/dphi]
    ! for energy conservation, need -Qs[dhs/dphi]
    call get_fields (gvmu, field(:,:,:,:,1), dum1, dist='h')

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       vpadiff_response(1,1,ikxkyz) = 1.0-field(iky,ikx,iz,it,1)
    end do
    idx = 2
    if (momentum_conservation) then
       call get_upar (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          vpadiff_response(idx:idx+nspec-1,1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do
       idx = idx + nspec
    end if
    if (energy_conservation) then
       call get_temp (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          vpadiff_response(idx:idx+nspec-1,1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do
    end if
    idx = 2

    if (momentum_conservation) then
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iz = iz_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          do imu = 1, nmu
             gvmu(:,imu,ikxkyz) = 2.*code_dt*spec(is)%vnew(is)*vpa*aj0v(imu,ikxkyz)*maxwell_vpa(:,is)*maxwell_mu(1,iz,imu,is)
             call tridag (1, aa_vpa(:,is), bb_vpa(:,ikxkyz), cc_vpa(:,is), gvmu(:,imu,ikxkyz))
          end do
       end do
       ! gvmu now contains dhs/dupars
       ! need to get -Ps[dhs/dupars] for phi equation
       ! need to get 1-Us[dhs/dupars] for momentum conservation
       ! need to get -Qs[dhs/dupars] for energy conservation
       call get_fields_by_spec (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          vpadiff_response(1,idx:idx+nspec-1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do

       call get_upar (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          do is = 1, nspec
             vpadiff_response(idx+is-1,idx+is-1,ikxkyz) = 1.0-field(iky,ikx,iz,it,is)
          end do
       end do

       if (energy_conservation) then
          call get_temp (gvmu, field)
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
             iky = iky_idx(kxkyz_lo,ikxkyz)
             ikx = ikx_idx(kxkyz_lo,ikxkyz)
             iz = iz_idx(kxkyz_lo,ikxkyz)
             it = it_idx(kxkyz_lo,ikxkyz)
             do is = 1, nspec
                vpadiff_response(idx+is+nspec-1,idx+is-1,ikxkyz) = -field(iky,ikx,iz,it,is)
             end do
          end do
       end if
       idx = idx + nspec
    end if

    if (energy_conservation) then
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iz = iz_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          do imu = 1, nmu
             gvmu(:,imu,ikxkyz) = 2.*code_dt*spec(is)%vnew(is)*(vpa**2+vperp2(1,iz,imu)-1.5) &
                  *aj0v(imu,ikxkyz)*maxwell_vpa(:,is)*maxwell_mu(1,iz,imu,is)
             call tridag (1, aa_vpa(:,is), bb_vpa(:,ikxkyz), cc_vpa(:,is), gvmu(:,imu,ikxkyz))
          end do
       end do
       ! gvmu now contains dhs/dQs
       ! need to get -Ps[dhs/dQs] for phi equation
       ! need to get 1-Us[dhs/dQs] for momentum conservation
       ! need to get -Qs[dhs/dQs] for energy conservation
       call get_fields_by_spec (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          vpadiff_response(1,idx:idx+nspec-1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do

       if (momentum_conservation) then
          call get_upar (gvmu, field)
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
             iky = iky_idx(kxkyz_lo,ikxkyz)
             ikx = ikx_idx(kxkyz_lo,ikxkyz)
             iz = iz_idx(kxkyz_lo,ikxkyz)
             it = it_idx(kxkyz_lo,ikxkyz)
             do is = 1, nspec
                vpadiff_response(idx+is-1-nspec,idx+is-1,ikxkyz) = -field(iky,ikx,iz,it,is)
             end do
          end do
       end if

       call get_temp (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          do is = 1, nspec
             vpadiff_response(idx+is-1,idx+is-1,ikxkyz) = 1.0-field(iky,ikx,iz,it,is)
          end do
       end do
    end if

    ! now get LU decomposition for vpadiff_response
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       call lu_decomposition (vpadiff_response(:,:,ikxkyz),vpadiff_idx(:,ikxkyz),dum2)
    end do

    ! reset wgts to default setting
    conservative_wgts = .false.
    call set_vpa_weights (conservative_wgts)

    deallocate (dum1, field)

  end subroutine init_vpadiff_conserve

  subroutine init_mudiff_conserve

    use finite_differences, only: tridag
    use linear_solve, only: lu_decomposition
    use stella_time, only: code_dt
    use species, only: nspec, spec
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: ztmax, maxwell_vpa, maxwell_mu
    use vpamu_grids, only: nvpa, vpa, vperp2
    use kt_grids, only: naky, nakx
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use dist_fn_arrays, only: gvmu, kperp2
    use gyro_averages, only: aj0v, aj1v
    use fields, only: get_fields, get_fields_by_spec
    use stella_geometry, only: bmag

    implicit none

    integer :: ikxkyz, iky, ikx, iz, it, is, ia
    integer :: iv
    integer :: idx
    real :: dum2
    complex, dimension (:,:,:,:), allocatable :: dum1
    complex, dimension (:,:,:,:,:), allocatable :: field

    if (.not.allocated(mudiff_response)) then
       nresponse_mu = 1
       if (momentum_conservation) nresponse_mu = nresponse_mu + nspec
       if (energy_conservation) nresponse_mu = nresponse_mu + nspec
       allocate (mudiff_response(nresponse_mu,nresponse_mu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       mudiff_response = 0.
       allocate (mudiff_idx(nresponse_mu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    end if
    allocate (dum1(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (field(naky,nakx,-nzgrid:nzgrid,ntubes,nspec))

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do iv = 1, nvpa
          gvmu(iv,:,ikxkyz) = ztmax(iv,is)*maxwell_mu(1,iz,:,is)*aj0v(:,ikxkyz)
          call tridag (1, aa_mu(iz,:,is), bb_mu(:,ikxkyz), cc_mu(iz,:,is), gvmu(iv,:,ikxkyz))
       end do
    end do

    ! gvmu contains dhs/dphi
    ! for phi equation, need 1-P[dhs/dphi]
    ! for uperp equations, need -Us[dhs/dphi]
    ! for energy conservation, need -Qs[dhs/dphi]
    call get_fields (gvmu, field(:,:,:,:,1), dum1, dist='h')

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       mudiff_response(1,1,ikxkyz) = 1.0-field(iky,ikx,iz,it,1)
    end do
    idx = 2
    if (momentum_conservation) then
       call get_uperp (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          mudiff_response(idx:idx+nspec-1,1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do
       idx = idx + nspec
    end if
    if (energy_conservation) then
       call get_temp_mu (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          mudiff_response(idx:idx+nspec-1,1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do
    end if
    idx = 2

    if (momentum_conservation) then
       ia = 1
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          do iv = 1, nvpa
             gvmu(iv,:,ikxkyz) = 2.*code_dt*spec(is)%vnew(is)*kperp2(iky,ikx,ia,iz)*vperp2(ia,iz,:) &
                  *(spec(is)%smz/bmag(ia,iz))**2*aj1v(:,ikxkyz)*maxwell_vpa(iv,is)*maxwell_mu(ia,iz,:,is)
             call tridag (1, aa_mu(iz,:,is), bb_mu(:,ikxkyz), cc_mu(iz,:,is), gvmu(iv,:,ikxkyz))
          end do
       end do
       ! gvmu now contains dhs/dupars
       ! need to get -Ps[dhs/dupars] for phi equation
       ! need to get 1-Us[dhs/dupars] for momentum conservation
       ! need to get -Qs[dhs/dupars] for energy conservation
       call get_fields_by_spec (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          mudiff_response(1,idx:idx+nspec-1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do

       call get_uperp (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          do is = 1, nspec
             mudiff_response(idx+is-1,idx+is-1,ikxkyz) = 1.0-field(iky,ikx,iz,it,is)
          end do
       end do

       if (energy_conservation) then
          call get_temp_mu (gvmu, field)
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
             iky = iky_idx(kxkyz_lo,ikxkyz)
             ikx = ikx_idx(kxkyz_lo,ikxkyz)
             iz = iz_idx(kxkyz_lo,ikxkyz)
             it = it_idx(kxkyz_lo,ikxkyz)
             do is = 1, nspec
                mudiff_response(idx+is+nspec-1,idx+is-1,ikxkyz) = -field(iky,ikx,iz,it,is)
             end do
          end do
       end if
       idx = idx + nspec
    end if

    if (energy_conservation) then
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iz = iz_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          do iv = 1, nvpa
             gvmu(iv,:,ikxkyz) = 2.0*code_dt*spec(is)%vnew(is)*(vpa(iv)**2+vperp2(1,iz,:)-1.5) &
                  *aj0v(:,ikxkyz)*maxwell_vpa(iv,is)*maxwell_mu(1,iz,:,is)
             call tridag (1, aa_mu(iz,:,is), bb_mu(:,ikxkyz), cc_mu(iz,:,is), gvmu(iv,:,ikxkyz))
          end do
       end do
       ! gvmu now contains dhs/dQs
       ! need to get -Ps[dhs/dQs] for phi equation
       ! need to get 1-Us[dhs/dQs] for momentum conservation
       ! need to get -Qs[dhs/dQs] for energy conservation
       call get_fields_by_spec (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          mudiff_response(1,idx:idx+nspec-1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do

       if (momentum_conservation) then
          call get_uperp (gvmu, field)
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
             iky = iky_idx(kxkyz_lo,ikxkyz)
             ikx = ikx_idx(kxkyz_lo,ikxkyz)
             iz = iz_idx(kxkyz_lo,ikxkyz)
             it = it_idx(kxkyz_lo,ikxkyz)
             do is = 1, nspec
                mudiff_response(idx+is-1-nspec,idx+is-1,ikxkyz) = -field(iky,ikx,iz,it,is)
             end do
          end do
       end if

       call get_temp_mu (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          do is = 1, nspec
             mudiff_response(idx+is-1,idx+is-1,ikxkyz) = 1.0-field(iky,ikx,iz,it,is)
          end do
       end do
    end if

    ! now get LU decomposition for mudiff_response
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       call lu_decomposition (mudiff_response(:,:,ikxkyz),mudiff_idx(:,ikxkyz),dum2)
    end do

    deallocate (dum1, field)

  end subroutine init_mudiff_conserve

  subroutine get_upar (g, fld)

    use mp, only: sum_allreduce
    use zgrid, only: nzgrid
    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: vpa
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use gyro_averages, only: aj0v

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,:), intent (out) :: fld

    integer :: ikxkyz, iky, ikx, iz, it, is
    complex, dimension (:,:), allocatable :: g0

    allocate (g0(nvpa,nmu))

    fld = 0.
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       g0 = g(:,:,ikxkyz)*spread(vpa,2,nmu)*spread(aj0v(:,ikxkyz),1,nvpa)
       call integrate_vmu (g0,iz,fld(iky,ikx,iz,it,is))
    end do
    deallocate (g0)

    call sum_allreduce (fld)

  end subroutine get_upar

  subroutine get_uperp (g, fld)

    use mp, only: sum_allreduce
    use zgrid, only: nzgrid
    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: vperp2
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use gyro_averages, only: aj1v

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,:), intent (out) :: fld

    integer :: ikxkyz, iky, ikx, iz, it, is
    complex, dimension (:,:), allocatable :: g0

    allocate (g0(nvpa,nmu))

    fld = 0.
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
!       g0 = 2.0*g(:,:,ikxkyz)*spread((vperp2(1,iz,:)-0.5)*aj1v(:,ikxkyz),1,nvpa)
       g0 = g(:,:,ikxkyz)*spread((vperp2(1,iz,:)-0.5)*aj1v(:,ikxkyz),1,nvpa)
       call integrate_vmu (g0,iz,fld(iky,ikx,iz,it,is))
    end do
    deallocate (g0)

    call sum_allreduce (fld)

  end subroutine get_uperp

  subroutine get_temp (g, fld)

    use mp, only: sum_allreduce
    use zgrid, only: nzgrid
    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: vpa
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use gyro_averages, only: aj0v

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,:), intent (out) :: fld

    integer :: ikxkyz, iky, ikx, iz, it, is
    complex, dimension (:,:), allocatable :: g0

    allocate (g0(nvpa,nmu))

    fld = 0.
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       g0 = g(:,:,ikxkyz)*(spread(vpa**2,2,nmu)-0.5) &
          * spread(aj0v(:,ikxkyz),1,nvpa)/1.5
       call integrate_vmu (g0,iz,fld(iky,ikx,iz,it,is))
    end do
    deallocate (g0)

    call sum_allreduce (fld)

  end subroutine get_temp

  subroutine get_temp_mu (g, fld)

    use mp, only: sum_allreduce
    use zgrid, only: nzgrid
    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: nvpa, nmu, vperp2
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use gyro_averages, only: aj0v

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,:), intent (out) :: fld

    integer :: ikxkyz, iky, ikx, iz, it, is
    complex, dimension (:,:), allocatable :: g0

    allocate (g0(nvpa,nmu))

    fld = 0.
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       g0 = g(:,:,ikxkyz)*(spread(vperp2(1,iz,:),1,nvpa)-1.0) &
          * spread(aj0v(:,ikxkyz),1,nvpa)/1.5
       call integrate_vmu (g0,iz,fld(iky,ikx,iz,it,is))
    end do
    deallocate (g0)

    call sum_allreduce (fld)

  end subroutine get_temp_mu

  subroutine finish_dissipation

    implicit none

    call finish_collisions

  end subroutine finish_dissipation

  subroutine finish_collisions

    use dist_fn_arrays, only: g_krook, g_proj

    implicit none

    if (collisions_implicit) then
       call finish_vpadiff_matrix
       call finish_mudiff_matrix
       call finish_vpadiff_response
       call finish_mudiff_response
    end if

    if(allocated(g_krook)) deallocate(g_krook)
    if(allocated(g_proj))  deallocate(g_proj)

    if (collision_model == "fokker-planck") then
        call finish_nusDpa
        call finish_fp_diffmatrix
    end if

    collisions_initialized = .false.

  end subroutine finish_collisions

  subroutine finish_fp_diffmatrix

    implicit none

    if (allocated(aa_vpa)) deallocate (aa_vpa)
    if (allocated(bb_vpa)) deallocate (bb_vpa)
    if (allocated(cc_vpa)) deallocate (cc_vpa)
    if (allocated(cdiffmat_band)) deallocate (cdiffmat_band)

  end subroutine finish_fp_diffmatrix

  subroutine finish_vpadiff_matrix

    implicit none

    if (allocated(aa_vpa)) deallocate (aa_vpa)
    if (allocated(bb_vpa)) deallocate (bb_vpa)
    if (allocated(cc_vpa)) deallocate (cc_vpa)

  end subroutine finish_vpadiff_matrix

  subroutine finish_mudiff_matrix

    implicit none

    if (allocated(aa_mu)) deallocate (aa_mu)
    if (allocated(bb_mu)) deallocate (bb_mu)
    if (allocated(cc_mu)) deallocate (cc_mu)

  end subroutine finish_mudiff_matrix

  subroutine finish_vpadiff_response

    implicit none

    if (allocated(vpadiff_response)) deallocate (vpadiff_response)
    if (allocated(vpadiff_idx)) deallocate (vpadiff_idx)

  end subroutine finish_vpadiff_response

  subroutine finish_mudiff_response

    implicit none

    if (allocated(mudiff_response)) deallocate (mudiff_response)
    if (allocated(mudiff_idx)) deallocate (mudiff_idx)

  end subroutine finish_mudiff_response

  subroutine add_krook_operator (g, gke_rhs)

    use zgrid, only: nzgrid, ntubes
    use constants, only: zi
    use kt_grids, only: akx, nakx, zonal_mode
    use stella_layouts, only: vmu_lo
    use stella_time, only: code_dt
    use dist_fn_arrays, only: g_krook
    use stella_geometry, only: dl_over_b

    implicit none

    real :: exp_fac
    complex :: tmp
    integer :: ikx, it, ia, ivmu

    !complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), optional, intent (in) :: f0
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gke_rhs

    ia = 1

    if(.not.zonal_mode(1)) return

    !TODO: add number and momentum conservation
    if(delay_krook.le.epsilon(0.)) then
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do ikx = 1, nakx
            if(abs(akx(ikx)).gt.akx(ikxmax_source)) cycle
            tmp = sum(dl_over_b(ia,:)*g(1,ikx,:,it,ivmu))
            if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) tmp = zi*aimag(tmp)
            gke_rhs(1,ikx,:,it,ivmu) = gke_rhs(1,ikx,:,it,ivmu) - code_dt*nu_krook*tmp
          enddo
        enddo
      enddo
    else
      exp_fac = exp(-code_dt/delay_krook)
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do ikx = 1, nakx
            if(abs(akx(ikx)).gt.akx(ikxmax_source)) cycle
            tmp = sum(dl_over_b(ia,:)*g(1,ikx,:,it,ivmu))
            if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) tmp = zi*aimag(tmp)
            gke_rhs(1,ikx,:,it,ivmu) = gke_rhs(1,ikx,:,it,ivmu) - code_dt*nu_krook &
                                     * (code_dt*tmp + exp_fac*int_krook*g_krook(ikx,it,ivmu)) &
                                     / (code_dt     + exp_fac*int_krook)
          enddo
        enddo
      enddo
    endif

  end subroutine add_krook_operator 

  subroutine update_delay_krook (g)

    use constants, only: zi
    use dist_fn_arrays, only: g_krook
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: akx, nakx, zonal_mode
    use stella_layouts, only: vmu_lo
    use stella_geometry, only: dl_over_b
    use stella_time, only: code_dt

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (in) :: g

    integer :: ivmu, it, ikx, ia
    real :: int_krook_old, exp_fac
    complex :: tmp

    if(.not.zonal_mode(1)) return

    exp_fac = exp(-code_dt/delay_krook)
    
    ia = 1

    int_krook_old = int_krook
    int_krook =  code_dt + exp_fac*int_krook_old

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      do it = 1, ntubes
        do ikx = 1, nakx
          tmp = sum(dl_over_b(ia,:)*g(1,ikx,:,it,ivmu))
          if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) tmp = zi*aimag(tmp)
          g_krook(ikx,it,ivmu) = (code_dt*tmp + exp_fac*int_krook_old*g_krook(ikx,it,ivmu))/int_krook
        enddo
      enddo
    enddo

    !g_krook   = (code_dt*g + exp_fac*int_krook_old*g_krook)/int_krook

  end subroutine update_delay_krook

  subroutine project_out_zero (g)

    use zgrid, only: nzgrid, ntubes
    use constants, only: zi
    use kt_grids, only: zonal_mode, akx, nakx
    use stella_layouts, only: vmu_lo
    use stella_time, only: code_dt
    use dist_fn_arrays, only: g_proj
    use stella_geometry, only: dl_over_b

    implicit none

    real :: exp_fac
    complex :: tmp
    integer :: ikx, it, ia, ivmu

    complex, dimension (:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (inout) :: g

    ia = 1
    if(.not.zonal_mode(1)) return

    exp_fac = exp(-code_dt/delay_krook)

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      do it = 1, ntubes
        do ikx = 1, nakx
          if(abs(akx(ikx)).gt.akx(ikxmax_source)) then
            g(ikx,:,it,ivmu) = 0.0
          else
            tmp = sum(dl_over_b(ia,:)*g(ikx,:,it,ivmu))
            if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) tmp = zi*aimag(tmp)
            if(delay_krook.le.epsilon(0.)) then
              g(ikx,:,it,ivmu) = tmp
            else
              g(ikx,:,it,ivmu) = (code_dt*tmp + exp_fac*int_proj*g_proj(ikx,it,ivmu)) &
                               / (code_dt     + exp_fac*int_proj)
            endif
          endif
          if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) then
            g_proj(ikx,it,ivmu) = zi*aimag(sum(dl_over_b(ia,:)*g(ikx,:,it,ivmu)))
          else
            g_proj(ikx,it,ivmu) = sum(dl_over_b(ia,:)*g(ikx,:,it,ivmu))
          endif
        enddo
      enddo
    enddo

    int_proj = code_dt + exp_fac*int_proj

  end subroutine project_out_zero

  subroutine advance_collisions_explicit (g, phi, gke_rhs)

    use mp, only: proc0
    use job_manage, only: time_message
    use redistribute, only: scatter, gather
    use stella_time, only: code_dt
    use zgrid, only: nzgrid, ntubes
    use species, only: spec
    use run_parameters, only: fphi
    use kt_grids, only: naky, nakx
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: set_vpa_weights
    use stella_geometry, only: bmag
    use stella_layouts, only: vmu_lo, kxkyz_lo
    use stella_layouts, only: is_idx, iky_idx, ikx_idx, iz_idx
    use dist_redistribute, only: kxkyz2vmu
    use dist_fn_arrays, only: gvmu, kperp2
    use g_tofrom_h, only: g_to_h

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phi
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gke_rhs

    integer :: is, ikxkyz, imu, iv, ivmu, ikx, iky, iz, ia
    logical :: conservative_wgts
    complex, dimension (:), allocatable :: mucoll
    complex, dimension (:,:,:), allocatable :: coll
    complex, dimension (:,:,:,:,:), allocatable :: tmp_vmulo

    complex, dimension (:,:,:), allocatable :: mucoll_fp
    complex, dimension (:,:,:), allocatable :: coll_fp

    if (proc0) call time_message(.false.,time_collisions(:,1),' collisions')

    allocate (tmp_vmulo(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    ! want exact conservation properties for collision operator
    conservative_wgts = .true.
    call set_vpa_weights (conservative_wgts)

    ! switch from g = <f> to h = f + Z*e*phi/T * F0
    tmp_vmulo = g
    call g_to_h (tmp_vmulo, phi, fphi)

    ! remap so that (vpa,mu) local
    if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')
    call scatter (kxkyz2vmu, tmp_vmulo, gvmu)
    if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')

    deallocate (tmp_vmulo)

    ia = 1

    ! take vpa derivatives
    if (collision_model=="dougherty") then
        allocate (coll(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
        allocate (mucoll(nmu))
        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
           iky = iky_idx(kxkyz_lo,ikxkyz)
           ikx = ikx_idx(kxkyz_lo,ikxkyz)
           iz = iz_idx(kxkyz_lo,ikxkyz)
           is = is_idx(kxkyz_lo,ikxkyz)
           if (vpa_operator) then
              do imu = 1, nmu
                 call vpa_differential_operator (gvmu(:,imu,ikxkyz), coll(:,imu,ikxkyz))
              end do
           end if
           if (mu_operator) then
              do iv = 1, nvpa
                 call mu_differential_operator (iz, ia, gvmu(iv,:,ikxkyz), mucoll)
                 coll(iv,:,ikxkyz) = coll(iv,:,ikxkyz) + mucoll
              end do
           end if
           if (momentum_conservation) call conserve_momentum (iky, ikx, iz, is, ikxkyz, gvmu(:,:,ikxkyz), coll(:,:,ikxkyz))
           if (energy_conservation) call conserve_energy (iz, is, ikxkyz, gvmu(:,:,ikxkyz), coll(:,:,ikxkyz))
           ! save memory by using gvmu and deallocating coll below
           ! before re-allocating tmp_vmulo
           gvmu(:,:,ikxkyz) = coll(:,:,ikxkyz) - 0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*gvmu(:,:,ikxkyz)
       end do
       deallocate (coll, mucoll)

       allocate (tmp_vmulo(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       ! remap so that (ky,kx,z,tube) local
       call gather (kxkyz2vmu, gvmu, tmp_vmulo)

       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          gke_rhs(:,:,:,:,ivmu) =  gke_rhs(:,:,:,:,ivmu) + code_dt*spec(is)%vnew(is)*tmp_vmulo(:,:,:,:,ivmu)
       end do
       deallocate (tmp_vmulo)
    end if

    if (collision_model=="fokker-planck") then
        allocate (coll_fp(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc)); coll_fp = 0.0
        allocate (mucoll_fp(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc)); mucoll_fp = 0.0

        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo,ikxkyz)
            ikx = ikx_idx(kxkyz_lo,ikxkyz)
            iz = iz_idx(kxkyz_lo,ikxkyz)
            is = is_idx(kxkyz_lo,ikxkyz)
            if (vpa_operator) then
              do imu = 1, nmu
                  call vpa_differential_operator_fp (gvmu(:,:,ikxkyz), coll_fp(:,:,ikxkyz), imu, iz, is, ia)
              end do
            end if
            if (mu_operator) then
              do iv = 1, nvpa
                 call mu_differential_operator_fp (gvmu(:,:,ikxkyz), mucoll_fp(:,:,ikxkyz), iv, iz, is, ia, iky, ikx, cfac)
              end do
            end if
            gvmu(:,:,ikxkyz) = coll_fp(:,:,ikxkyz) + mucoll_fp(:,:,ikxkyz)
        end do
        deallocate (coll_fp, mucoll_fp)

        allocate (tmp_vmulo(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        call gather (kxkyz2vmu, gvmu, tmp_vmulo)

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
           gke_rhs(:,:,:,:,ivmu) = gke_rhs(:,:,:,:,ivmu) + code_dt*tmp_vmulo(:,:,:,:,ivmu)
        end do
        deallocate (tmp_vmulo)
    end if

    ! reset to default integration wgts
    conservative_wgts = .false.
    call set_vpa_weights (conservative_wgts)

    if (proc0) call time_message(.false.,time_collisions(:,1),' collisions')

  end subroutine advance_collisions_explicit

  subroutine vpa_differential_operator_fp (h, Dh, imu, iz, is, ia)

      use vpamu_grids, only: nvpa, vpa, dvpa, mu, dmu, nmu, equally_spaced_mu_grid, maxwell_mu, maxwell_vpa
      use stella_geometry, only: bmag
      use constants, only: pi
      use species, only: spec

      implicit none

      complex, dimension (:,:), intent (out) :: Dh
      complex, dimension (:,:), intent (in) :: h
      integer, intent (in) :: imu, iz, ia, is
      integer :: iv
      complex :: Dhmu, Dhmu_u, Dhmu_l, dmuhp, dmuhm, dvpah, dvpahp, dvpahm
      real :: xp, xm, vpap, vpam, nupap, nupam, nuDp, nuDm, mwp, mwm, a, b, c

      iv = 1
      vpap = 0.5*(vpa(iv)+vpa(iv+1))
      xp   = sqrt(vpap**2 + 2*bmag(ia,iz)*mu(imu))
      nuDp = spec(is)%vnew(is)*(erf(xp) - (erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) )/xp**3
      nupap= spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
      mwp  = exp(-vpap**2)*maxwell_mu(1,iz,imu,is)
      dvpahp = (h(iv+1,imu)/mw(iv+1,imu,iz) - h(iv,imu)/mw(iv,imu,iz))/dvpa

      if (imu == 1) then
          ! second-order accurate
          !a = -(2.*dmu(1) + dmu(2))/(dmu(1)*(dmu(1)+dmu(2)))
          !b = (dmu(1)+dmu(2))/(dmu(1)*dmu(2))
          !c = -dmu(1)/(dmu(2)*(dmu(1)+dmu(2)))
          !dmuhp = a*h(iv+1,imu)/mw(iv+1,imu,iz)+ b*h(iv+1,imu+1)/mw(iv+1,imu+1,iz) + c*h(iv+1,imu+2)/mw(iv+1,imu+2,iz)
          ! first-order accurate, as in implicit routine:
          dmuhp = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz) - h(iv+1,imu)/mw(iv+1,imu,iz))/dmu(imu)
      else if (imu == nmu) then
          ! second-order accurate
          !a = dmu(nmu-1)/(dmu(nmu-2)*(dmu(nmu-2)+dmu(nmu-1)))
          !b = -(dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-2)*dmu(nmu-1))
          !c = (2.*dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-1)*(dmu(nmu-1)+dmu(nmu-2)))
          !dmuhp = a*h(iv+1,imu-2)/mw(iv+1,imu-2,iz) + b*h(iv+1,imu-1)/mw(iv+1,imu-1,iz) + c*h(iv+1,imu)/mw(iv+1,imu,iz)
          ! first-order accurate, as in implicit routine:
          dmuhp = (h(iv+1,imu)/mw(iv+1,imu,iz) - h(iv+1,imu-1)/mw(iv+1,imu-1,iz))/dmu(imu-1)
      else
          dmuhp = ((h(iv+1,imu)/mw(iv+1,imu,iz)-h(iv+1,imu-1)/mw(iv+1,imu-1,iz))*dmu(imu)/dmu(imu-1) &
                  + (h(iv+1,imu+1)/mw(iv+1,imu+1,iz)-h(iv+1,imu)/mw(iv+1,imu,iz))*dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu))
      end if
      Dh(iv,imu) = (2*0.5*(nupap*vpap**2 + 2*nuDp*bmag(ia,iz)*mu(imu))*mwp*dvpahp + vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*mw(iv+1,imu,iz)*dmuhp)/(2*dvpa)

      iv = nvpa
      vpam = 0.5*(vpa(iv)+vpa(iv-1))
      xm   = sqrt(vpam**2 + 2*bmag(ia,iz)*mu(imu))
      nuDm = spec(is)%vnew(is)*(erf(xm) - (erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) )/xm**3
      nupam= spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
      mwm  = exp(-vpam**2)*maxwell_mu(1,iz,imu,is)
      dvpahm = (h(iv,imu)/mw(iv,imu,iz) - h(iv-1,imu)/mw(iv-1,imu,iz))/dvpa

      if (imu == 1) then
          ! second-order accurate
          !a = -(2.*dmu(1) + dmu(2))/(dmu(1)*(dmu(1)+dmu(2)))
          !b = (dmu(1)+dmu(2))/(dmu(1)*dmu(2))
          !c = -dmu(1)/(dmu(2)*(dmu(1)+dmu(2)))
          !dmuhm = a*h(iv-1,imu)/mw(iv-1,imu,iz) + b*h(iv-1,imu+1)/mw(iv-1,imu+1,iz) + c*h(iv-1,imu+2)/mw(iv-1,imu+2,iz)
          ! first-order accurate, as in implicit routine:
          dmuhm = (h(iv-1,imu+1)/mw(iv-1,imu+1,iz) - h(iv-1,imu)/mw(iv-1,imu,iz))/dmu(imu)
      else if (imu == nmu) then
          ! second-order accurate
          !a = dmu(nmu-1)/(dmu(nmu-2)*(dmu(nmu-2)+dmu(nmu-1)))
          !b = -(dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-2)*dmu(nmu-1))
          !c = (2.*dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-1)*(dmu(nmu-1)+dmu(nmu-2)))
          !dmuhm = a*h(iv-1,imu-2)/mw(iv-1,imu-2,iz) + b*h(iv-1,imu-1)/mw(iv-1,imu-1,iz) + c*h(iv-1,imu)/mw(iv-1,imu,iz)
          ! first-order accurate, as in implicit routine:
          dmuhm = (h(iv-1,imu)/mw(iv-1,imu,iz) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz))/dmu(imu-1)
      else
          dmuhm = ((h(iv-1,imu)/mw(iv-1,imu,iz)-h(iv-1,imu-1)/mw(iv-1,imu-1,iz))*dmu(imu)/dmu(imu-1) + (h(iv-1,imu+1)/mw(iv-1,imu+1,iz)-h(iv-1,imu)/mw(iv-1,imu,iz))*dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu))
      end if
      Dh(iv,imu) = (-2*0.5*(nupam*vpam**2 + 2*nuDm*bmag(ia,iz)*mu(imu))*mwm*dvpahm - vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*mw(iv-1,imu,iz)*dmuhm) / (2*dvpa)

      do iv = 2, nvpa-1
          ! AVB: interior nodes:
          ! quantities at half-grid-points:
          vpap = 0.5*(vpa(iv)+vpa(iv+1))
          vpam = 0.5*(vpa(iv)+vpa(iv-1))
          xp   = sqrt(vpap**2 + 2*bmag(ia,iz)*mu(imu))
          xm   = sqrt(vpam**2 + 2*bmag(ia,iz)*mu(imu))
          nuDp = spec(is)%vnew(is)*(erf(xp) - (erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) )/xp**3
          nuDm = spec(is)%vnew(is)*(erf(xm) - (erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) )/xm**3
          nupap= spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
          nupam= spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
          mwp  = exp(-vpap**2)*maxwell_mu(1,iz,imu,is)
          mwm  = exp(-vpam**2)*maxwell_mu(1,iz,imu,is)
          dvpahp = (h(iv+1,imu)/mw(iv+1,imu,iz) - h(iv  ,imu)/mw(iv,imu,iz))/dvpa
          dvpahm = (h(iv  ,imu)/mw(iv,imu,iz) - h(iv-1,imu)/mw(iv-1,imu,iz))/dvpa

          if (imu == 1) then
              ! second-order accurate:
              !a = -(2.*dmu(1) + dmu(2))/(dmu(1)*(dmu(1)+dmu(2)))
              !b = (dmu(1)+dmu(2))/(dmu(1)*dmu(2))
              !c = -dmu(1)/(dmu(2)*(dmu(1)+dmu(2)))
              !dmuhp = a*h(iv+1,imu)/mw(iv+1,imu,iz) + b*h(iv+1,imu+1)/mw(iv+1,imu+1,iz) + c*h(iv+1,imu+2)/mw(iv+1,imu+2,iz)
              !dmuhm = a*h(iv-1,imu)/mw(iv-1,imu,iz) + b*h(iv-1,imu+1)/mw(iv-1,imu+1,iz) + c*h(iv-1,imu+2)/mw(iv-1,imu+2,iz)
              ! or first-order accurate, as in implicit routine:
              dmuhp = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz) - h(iv+1,imu)/mw(iv+1,imu,iz))/dmu(imu)
              dmuhm = (h(iv-1,imu+1)/mw(iv-1,imu+1,iz) - h(iv-1,imu)/mw(iv-1,imu,iz))/dmu(imu)
          else if (imu == nmu) then
              ! second-order accurate:
              !a = dmu(nmu-1)/(dmu(nmu-2)*(dmu(nmu-2)+dmu(nmu-1)))
              !b = -(dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-2)*dmu(nmu-1))
              !c = (2.*dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-1)*(dmu(nmu-1)+dmu(nmu-2)))
              !dmuhp = a*h(iv+1,imu-2)/mw(iv+1,imu-2,iz) + b*h(iv+1,imu-1)/mw(iv+1,imu-1,iz) + c*h(iv+1,imu)/mw(iv+1,imu,iz)
              !dmuhm = a*h(iv-1,imu-2)/mw(iv-1,imu-2,iz) + b*h(iv-1,imu-1)/mw(iv-1,imu-1,iz) + c*h(iv-1,imu)/mw(iv-1,imu,iz)
              ! or first-order accurate, as in implicit routine:
              dmuhp = (h(iv+1,imu)/mw(iv+1,imu,iz) - h(iv+1,imu-1)/mw(iv+1,imu-1,iz))/dmu(imu-1)
              dmuhm = (h(iv-1,imu)/mw(iv-1,imu,iz) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz))/dmu(imu-1)
          else
              dmuhp = ((h(iv+1,imu)/mw(iv+1,imu,iz)-h(iv+1,imu-1)/mw(iv+1,imu-1,iz))*dmu(imu)/dmu(imu-1) + (h(iv+1,imu+1)/mw(iv+1,imu+1,iz)-h(iv+1,imu)/mw(iv+1,imu,iz))*dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu))
              dmuhm = ((h(iv-1,imu)/mw(iv-1,imu,iz)-h(iv-1,imu-1)/mw(iv-1,imu-1,iz))*dmu(imu)/dmu(imu-1) + (h(iv-1,imu+1)/mw(iv-1,imu+1,iz)-h(iv-1,imu)/mw(iv-1,imu,iz))*dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu))
          end if
          Dh(iv,imu) = (2*0.5*(nupap*vpap**2 + 2*nuDp*bmag(ia,iz)*mu(imu))*mwp*dvpahp + vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz)*mw(iv+1,imu,iz)*dmuhp &
                      - 2*0.5*(nupam*vpam**2 + 2*nuDm*bmag(ia,iz)*mu(imu))*mwm*dvpahm - vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz)*mw(iv-1,imu,iz)*dmuhm) / (2*dvpa)
      end do

  end subroutine vpa_differential_operator_fp

  subroutine mu_differential_operator_fp (h, Dh, iv, iz, is, ia, iky, ikx, cfac)

      use vpamu_grids, only: nmu, mu, dmu, vpa, dvpa, nvpa, maxwell_mu, maxwell_vpa, equally_spaced_mu_grid
      use stella_geometry, only: bmag
      use species, only: spec
      use dist_fn_arrays, only: kperp2
      use constants, only: pi
      use job_manage, only: timer_local, time_message
      use mp, only: proc0

      implicit none

      complex, dimension (:,:), intent (in) :: h
      complex, dimension (:,:), intent (out) :: Dh
      integer, intent (in) :: iv, iz, is, ia, iky, ikx
      real, intent (in) :: cfac
      complex ::  Dvpah, Dvpah_p, Dvpah_m, Dmuh, Dmuh_m, Dmuh_p, Dmuh1, Dmuh2
      real :: nuDp, nuDm, nupap, nupam, mup, mum, nusp, nuxp, xp, xm, mwm, mwp, a, b, c
      integer :: imu

      imu = 1
      ! vpa-differential terms:
      if (iv == 1) then
          ! AVB: first order accurate:
          Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz) - h(iv,imu)/mw(iv,imu,iz))/dvpa
          Dvpah_p = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz) - h(iv,imu+1)/mw(iv,imu+1,iz))/dvpa
      else if (iv == nvpa) then
          ! AVB: first order accurate:
          Dvpah   = (h(iv,imu)/mw(iv,imu,iz) - h(iv-1,imu)/mw(iv-1,imu,iz))/dvpa
          Dvpah_p = (h(iv,imu+1)/mw(iv,imu+1,iz) - h(iv-1,imu+1)/mw(iv-1,imu+1,iz))/dvpa
      else
          Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz) - h(iv-1,imu)/mw(iv-1,imu,iz))/(2*dvpa)
          Dvpah_p = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz) - h(iv-1,imu+1)/mw(iv-1,imu+1,iz))/(2*dvpa)
      end if

      ! first mu-derivative term at mu_{i}:
      ! use ghost cell at mu_{0} = 0, where term vanishes, so dmu(0) = mu(1).
      Dmuh1 = ((vpa(iv)*mu(imu  )*nux(iv,imu  ,iz)*mw(iv,imu,iz)*Dvpah)*dmu(imu)/mu(imu) &
              +(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*Dvpah_p-vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz)*Dvpah)*mu(imu)/dmu(imu)) / (mu(imu)+dmu(imu))
      ! first derivative of h, at mu_{i+1/2}, and at mu_i:
      Dmuh  = (h(iv,imu+1)/mw(iv,imu+1,iz) - h(iv,imu)/mw(iv,imu,iz))/dmu(imu) ! first-order accurate, as used in implicit routine
      ! for second-order accuracy:
      !a = -(2.*dmu(1) + dmu(2))/(dmu(1)*(dmu(1)+dmu(2)))
      !b = (dmu(1)+dmu(2))/(dmu(1)*dmu(2))
      !c = -dmu(1)/(dmu(2)*(dmu(1)+dmu(2)))
      !Dmuh   = a*h(iv,imu)/mw(iv,imu,iz) + b*h(iv,imu+1)/mw(iv,imu+1,iz) + c*h(iv,imu+2)/mw(iv,imu+2,iz) ! second order accurate
      Dmuh_p = (h(iv,imu+1)/mw(iv,imu+1,iz) - h(iv,imu)/mw(iv,imu,iz))/dmu(imu) ! second-order accurate
      ! quantities at mu_{i+1/2}:
      mup = 0.5*(mu(imu)+mu(imu+1))
      mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
      xp  = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
      nuDp  = spec(is)%vnew(is)*( erf(xp)-(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) )/xp**3
      nupap = spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
      ! second mu-derivative term at mu_{i}:
      ! use d/dmu[...]_{1} = ([...]_{1+1/2} - [...]_{0})/(dmu_{1}/2+mu(1)), where [...]_{0} is a ghost cell at mu_{0} = 0, with [...]_{0} = 0.
      Dmuh2 = ( (2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)*Dmuh )*dmu(imu)/2./mu(imu) &
               +(2*(nupap*mup**2+nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp*Dmuh_p - 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)*Dmuh)*mu(imu)/(dmu(imu)/2.) )/(mu(imu)+dmu(imu)/2.)
      ! add differential terms and gyro-diffusive term:
      Dh(iv,imu) = Dmuh2 + Dmuh1 - cfac*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*(nupa(iv,imu,iz)*bmag(ia,iz)*mu(imu) + nuD(iv,imu,iz)*(vpa(iv)**2 + bmag(ia,iz)*mu(imu)))*h(iv,imu)

      imu = nmu
      ! vpa-differential terms:
      if (iv == 1) then
          ! AVB: first order accurate:
          Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz) - h(iv,imu)/mw(iv,imu,iz))/dvpa
          Dvpah_m = (h(iv+1,imu-1)/mw(iv+1,imu-1,iz) - h(iv,imu-1)/mw(iv,imu-1,iz))/dvpa
      else if (iv == nvpa) then
          ! AVB: first order accurate:
          Dvpah   = (h(iv,imu)/mw(iv,imu,iz) - h(iv-1,imu)/mw(iv-1,imu,iz))/dvpa
          Dvpah_m = (h(iv,imu-1)/mw(iv,imu-1,iz) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz))/dvpa
      else
          Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz) - h(iv-1,imu)/mw(iv-1,imu,iz))/(2*dvpa)
          Dvpah_m = (h(iv+1,imu-1)/mw(iv+1,imu-1,iz) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz))/(2*dvpa)
      end if

      ! first mu-derivative term at mu_{nmu}:
      Dmuh1 = (( vpa(iv)*mu(imu )*nux(iv,imu  ,iz)*mw(iv,imu,iz)*Dvpah - vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*Dvpah_m)*dmu(imu-1)/dmu(imu-1) &
              +( - vpa(iv)*mu(imu  )*nux(iv,imu  ,iz)*mw(iv,imu  ,iz)*Dvpah )*dmu(imu-1)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu-1))
      ! first derivative of h, at mu_{nmu} and mu_{nmu-1/2}:
      Dmuh_m = (h(iv,imu)/mw(iv,imu,iz) - h(iv,imu-1)/mw(iv,imu-1,iz))/dmu(imu-1)
      Dmuh   = (h(iv,imu)/mw(iv,imu,iz) - h(iv,imu-1)/mw(iv,imu-1,iz))/dmu(imu-1) ! first-order accurate
      ! for second-order accuracy:
      !a = dmu(nmu-1)/(dmu(nmu-2)*(dmu(nmu-2)+dmu(nmu-1)))
      !b = -(dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-2)*dmu(nmu-1))
      !c = (2.*dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-1)*(dmu(nmu-1)+dmu(nmu-2)))
      !Dmuh = a*h(iv,imu-2)/mw(iv,imu-2,iz) + b*h(iv,imu-1)/mw(iv,imu-1,iz) + c*h(iv,imu)/mw(iv,imu,iz)
      ! quantities at mu_{nmu-1/2}:
      mum = 0.5*(mu(imu)+mu(imu-1))
      mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
      xm  = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
      nuDm  = spec(is)%vnew(is)*( erf(xm)-(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) )/xm**3
      nupam = spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
      ! second mu-derivative term at mu_{nmu}:
      ! use d/dmu[...]_{nmu} = ([...]_{nmu+1} - [...]_{nmu-1/2})/(dmu_{nmu-1}/2+dmu(nmu-1)), where [...]_{nmu+1} is a ghost cell at mu = mu_{nmu} + dmu(nmu-1), with [...]_{nmu+1} = 0.
      Dmuh2 = ( ( 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)*Dmuh - 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm*Dmuh_m)*dmu(imu-1)/(dmu(imu-1)/2) &
               +(-2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)*Dmuh)*dmu(imu-1)/2./dmu(imu-1) )/(dmu(imu-1)/2.+dmu(imu-1))
      ! add differential terms and gyro-diffusive term:
      Dh(iv,imu) = Dmuh2 + Dmuh1 - cfac*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*(nupa(iv,imu,iz)*bmag(ia,iz)*mu(imu) + nuD(iv,imu,iz)*(vpa(iv)**2 + bmag(ia,iz)*mu(imu)))*h(iv,imu)

      do imu = 2, nmu-1
          ! vpa-differential terms:
          if (iv == 1) then
              ! AVB: first order accurate:
              Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz) - h(iv,imu)/mw(iv,imu,iz))/dvpa
              Dvpah_p = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz) - h(iv,imu+1)/mw(iv,imu+1,iz))/dvpa
              Dvpah_m = (h(iv+1,imu-1)/mw(iv+1,imu-1,iz) - h(iv,imu-1)/mw(iv,imu-1,iz))/dvpa
          else if (iv == nvpa) then
              ! AVB: first order accurate:
              Dvpah   = (h(iv,imu)/mw(iv,imu,iz) - h(iv-1,imu)/mw(iv-1,imu,iz))/dvpa
              Dvpah_p = (h(iv,imu+1)/mw(iv,imu+1,iz) - h(iv-1,imu+1)/mw(iv-1,imu+1,iz))/dvpa
              Dvpah_m = (h(iv,imu-1)/mw(iv,imu-1,iz) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz))/dvpa
          else
              Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz) - h(iv-1,imu)/mw(iv-1,imu,iz))/(2*dvpa)
              Dvpah_p = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz) - h(iv-1,imu+1)/mw(iv-1,imu+1,iz))/(2*dvpa)
              Dvpah_m = (h(iv+1,imu-1)/mw(iv+1,imu-1,iz) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz))/(2*dvpa)
          end if
          ! first mu-derivative of vpa-derivative term, at mu_{i}:
          Dmuh1 = ( (vpa(iv)*mu(imu  )*nux(iv,imu  ,iz)*mw(iv,imu,iz)*Dvpah     - vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz)*mw(iv,imu-1,iz)*Dvpah_m)*dmu(imu)/dmu(imu-1) &
                   +(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz)*Dvpah_p - vpa(iv)*mu(imu  )*nux(iv,imu  ,iz)*mw(iv,imu,iz)*Dvpah    )*dmu(imu-1)/dmu(imu) ) / (dmu(imu-1)+dmu(imu))
          ! first mu-derivatives of h, at mu_i, mu_{i+1/2} and mu_{i-1/2}:
          Dmuh   = ( (h(iv,imu)/mw(iv,imu,iz)-h(iv,imu-1)/mw(iv,imu-1,iz))*dmu(imu)/dmu(imu-1) + (h(iv,imu+1)/mw(iv,imu+1,iz)-h(iv,imu)/mw(iv,imu,iz))*dmu(imu-1)/dmu(imu) ) / (dmu(imu-1)+dmu(imu))
          Dmuh_m = (h(iv,imu)/mw(iv,imu,iz) - h(iv,imu-1)/mw(iv,imu-1,iz))/dmu(imu-1)
          Dmuh_p = (h(iv,imu+1)/mw(iv,imu+1,iz) - h(iv,imu)/mw(iv,imu,iz))/dmu(imu)
          ! quantities at mu_{i+1/2} and mu_{i-1/2}:
          mup = 0.5*(mu(imu)+mu(imu+1))
          mum = 0.5*(mu(imu)+mu(imu-1))
          mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
          mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
          xp  = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
          xm  = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
          nuDp  = spec(is)%vnew(is)*( erf(xp)-(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2) ) / (2*xp**2) )/xp**3
          nuDm  = spec(is)%vnew(is)*( erf(xm)-(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2) ) / (2*xm**2) )/xm**3
          nupap = spec(is)%vnew(is)*2*( erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2) ) / (2*xp**2) / xp**3
          nupam = spec(is)%vnew(is)*2*( erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2) ) / (2*xm**2) / xm**3
          ! second mu-derivative term at mu_{i}:
          Dmuh2 = ( ( 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)*Dmuh - 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm*Dmuh_m )*dmu(imu)/dmu(imu-1) &
                   +( 2*(nupap*mup**2+nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp*Dmuh_p - 2*(nupa(iv,imu,iz)*mu(imu)**2+nuD(iv,imu,iz)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz)*Dmuh )*dmu(imu-1)/dmu(imu) )*2/(dmu(imu-1)+dmu(imu))
          ! add differential terms and gyro-diffusive term:
          Dh(iv,imu) = Dmuh2 + Dmuh1 - cfac*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*(nupa(iv,imu,iz)*bmag(ia,iz)*mu(imu) + nuD(iv,imu,iz)*(vpa(iv)**2 + bmag(ia,iz)*mu(imu)))*h(iv,imu)
      end do

  end subroutine mu_differential_operator_fp

  subroutine vpa_differential_operator (h, Dh)

    use vpamu_grids, only: nvpa, vpa, dvpa

    implicit none

    complex, dimension (:), intent (in) :: h
    complex, dimension (:), intent (out) :: Dh

    integer :: iv

    ! use h = 0 at ghost cells beyond +/- vpa_max
    iv = 1
    Dh(iv) = (0.5*h(iv+1)*(1.0/dvpa+vpa(iv+1))-h(iv)/dvpa)/dvpa
    iv = nvpa
    Dh(iv) = (-h(iv)/dvpa+0.5*h(iv-1)*(1.0/dvpa-vpa(iv-1)))/dvpa
    do iv = 2, nvpa-1
       Dh(iv) = (0.5*h(iv+1)*(1.0/dvpa+vpa(iv+1))-h(iv)/dvpa+0.5*h(iv-1)*(1.0/dvpa-vpa(iv-1)))/dvpa
    end do

  end subroutine vpa_differential_operator

  subroutine mu_differential_operator (iz, ia, h, Dh)

    use vpamu_grids, only: nmu, mu, dmu
    use vpamu_grids, only: equally_spaced_mu_grid
    use finite_differences, only: d2_3pt, fd3pt
    use stella_geometry, only: bmag

    implicit none

    integer, intent (in) :: iz, ia
    complex, dimension (:), intent (in) :: h
    complex, dimension (:), intent (out) :: Dh

    integer :: imu
    real :: mup, mum
    complex, dimension (:), allocatable :: h_ghost, Dh_ghost
    real, dimension (:), allocatable :: dmu_ghost

    allocate (h_ghost(nmu+1))
    allocate (Dh_ghost(nmu+1))
    allocate (dmu_ghost(nmu))

    if (equally_spaced_mu_grid) then
       ! use mu_{i-1/2} = 0 for i = 1
       imu = 1
       mup = 0.5*(mu(imu+1)+mu(imu))/(bmag(ia,iz)*dmu(1))
       Dh(imu) = (h(imu+1)*(mup+mu(imu+1)) &
            -h(imu)*(mup-mu(imu)))/dmu(1)
       ! use h = 0 at ghost cells beyond mu_max
       imu = nmu
       mup = 0.5*(2.*mu(imu)+dmu(1))/(bmag(ia,iz)*dmu(1))
       mum = 0.5*(mu(imu)+mu(imu-1))/(bmag(ia,iz)*dmu(1))
       Dh(imu) = (-h(imu)*(mup+mum) + h(imu-1)*(mum-mu(imu-1)))/dmu(1)
       do imu = 2, nmu-1
          mup = 0.5*(mu(imu+1)+mu(imu))/(bmag(ia,iz)*dmu(1))
          mum = 0.5*(mu(imu)+mu(imu-1))/(bmag(ia,iz)*dmu(1))
          Dh(imu) = (h(imu+1)*(mup+mu(imu+1)) &
               -h(imu)*(mup+mum) + h(imu-1)*(mum-mu(imu-1)))/dmu(1)
       end do
    else
       ! pad h_ghost array with ghost cell beyond max(mu) with zero BC
       h_ghost(:nmu) = h ; h_ghost(nmu+1) = 0.
       ! assign extra dmu value at nmu (beyond mu grid)
       ! because it will be accessed (but not later used)
       ! by generic subroutine d2_3pt
       dmu_ghost(:nmu-1) = dmu(:nmu-1) ; dmu_ghost(nmu) = 1.0

       call d2_3pt (h_ghost, Dh_ghost, dmu_ghost)
       Dh = Dh_ghost(:nmu)*mu/bmag(ia,iz)

       ! next add (1/B + 2*mu)*dh/dmu + 2*h
       call fd3pt (h_ghost, Dh_ghost, dmu_ghost)
       Dh = Dh + (1./bmag(ia,iz) + 2.*mu)*Dh_ghost(:nmu) + 2.*h
    end if

    deallocate (h_ghost, Dh_ghost, dmu_ghost)

  end subroutine mu_differential_operator

  subroutine conserve_momentum (iky, ikx, iz, is, ikxkyz, h, Ch)

    use species, only: spec
    use stella_geometry, only: bmag
    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: vpa, nvpa, nmu, vperp2
    use vpamu_grids, only: maxwell_vpa, maxwell_mu
    use dist_fn_arrays, only: kperp2
    use gyro_averages, only: aj0v, aj1v

    implicit none

    integer, intent (in) :: iky, ikx, iz, is, ikxkyz
    complex, dimension (:,:), intent (in) :: h
    complex, dimension (:,:), intent (in out) :: Ch

    complex, dimension (:,:), allocatable :: u_fac
    complex :: integral
    integer :: ia

    allocate (u_fac(nvpa,nmu))

    ia = 1

    u_fac = spread(aj0v(:,ikxkyz),1,nvpa)*spread(vpa,2,nmu)
    call integrate_vmu (u_fac*h,iz,integral)

    Ch = Ch + 2.0*u_fac*integral*spread(maxwell_mu(1,iz,:,is),1,nvpa)*spread(maxwell_vpa(:,is),2,nmu)

    u_fac = spread(vperp2(ia,iz,:)*aj1v(:,ikxkyz),1,nvpa)*sqrt(kperp2(iky,ikx,ia,iz))*spec(is)%smz/bmag(ia,iz)
    call integrate_vmu (u_fac*h,iz,integral)
    
    Ch = Ch + 2.0*u_fac*integral*spread(maxwell_mu(1,iz,:,is),1,nvpa)*spread(maxwell_vpa(:,is),2,nmu)

    deallocate (u_fac)

  end subroutine conserve_momentum

  subroutine conserve_energy (iz, is, ikxkyz, h, Ch)

    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: vpa, nvpa, nmu, vperp2
    use vpamu_grids, only: maxwell_vpa, maxwell_mu
    use gyro_averages, only: aj0v

    implicit none

    integer, intent (in) :: iz, is, ikxkyz
    complex, dimension (:,:), intent (in) :: h
    complex, dimension (:,:), intent (in out) :: Ch

    complex, dimension (:,:), allocatable :: T_fac
    complex :: integral

    allocate (T_fac(nvpa,nmu))

    T_fac = spread(aj0v(:,ikxkyz),1,nvpa)*(spread(vpa**2,2,nmu)+spread(vperp2(1,iz,:),1,nvpa)-1.5)
    call integrate_vmu (T_fac*h,iz,integral)

    Ch = Ch + 4.0*T_fac*integral*spread(maxwell_mu(1,iz,:,is),1,nvpa)*spread(maxwell_vpa(:,is),2,nmu)/3.0

    deallocate (T_fac)

  end subroutine conserve_energy

  subroutine advance_collisions_implicit (mirror_implicit, phi, apar, g)

    use mp, only: proc0
    use redistribute, only: gather, scatter
    use dist_redistribute, only: kxkyz2vmu
    use job_manage, only: time_message
    use zgrid, only: nzgrid
    use vpamu_grids, only: set_vpa_weights
    use stella_layouts, only: vmu_lo
    use dist_fn_arrays, only: gvmu

    implicit none

    logical, intent (in) :: mirror_implicit
    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g

    logical :: conservative_wgts

    if (proc0) call time_message(.false.,time_collisions(:,1),' collisions')

    conservative_wgts = .true.
    call set_vpa_weights (conservative_wgts)

    if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')
    call scatter (kxkyz2vmu, g, gvmu)
    if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')

    if (collision_model == "dougherty") then
        if (vpa_operator) call advance_vpadiff_implicit (phi, apar, gvmu)
        if (mu_operator) call advance_mudiff_implicit (phi, apar, gvmu)
    end if
    if (collision_model == "fokker-planck") then
        call advance_implicit_fp (phi, apar, gvmu)
    end if

    if (.not.mirror_implicit) then
       ! then take the results and remap again so ky,kx,z local.
       if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')
       call gather (kxkyz2vmu, gvmu, g)
       if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')
    end if

    conservative_wgts = .false.
    call set_vpa_weights (conservative_wgts)

    if (proc0) call time_message(.false.,time_collisions(:,1),' collisions')

  end subroutine advance_collisions_implicit

  subroutine advance_implicit_fp (phi, apar, g)

    use mp, only: sum_allreduce
    use finite_differences, only: tridag
    use linear_solve, only: lu_back_substitution
    use stella_time, only: code_dt
    use run_parameters, only: fphi
    use species, only: nspec, spec
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: nmu, nvpa
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa, vperp2
    use vpamu_grids, only: set_vpa_weights
    use kt_grids, only: naky, nakx
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, it_idx
    use g_tofrom_h, only: g_to_h
    use gyro_averages, only: aj0v
    use fields, only: get_fields

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar
    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in out) :: g

    real, dimension (:,:), allocatable :: tmp
    complex, dimension (:,:,:,:,:), allocatable :: flds
    complex, dimension (:,:,:), allocatable :: g_in
    complex, dimension (:,:), allocatable :: gvmutr
    complex, dimension (:), allocatable :: ghrs

    integer :: ikxkyz, iky, ikx, iz, is, imu, iv, it
    integer :: idx

    ! store input g for use later, as g will be overwritten below
    allocate (g_in(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    g_in = g

    allocate (gvmutr(nvpa,nmu))
    allocate (ghrs(nmu*nvpa))

    !! since backwards difference in time, (I-dt*D)h_inh^{n+1} = g^{***}
    !! g = g^{***}. tridag below inverts above equation to get h_inh^{n+1}
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       is = is_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       do iv = 1, nvpa
             ghrs(nmu*(iv-1)+1 : nmu*iv) = g(iv,:,ikxkyz)
       end do
       ! using lapack:
       call zgbtrs('No transpose', nvpa*nmu, nmu+1, nmu+1, 1, cdiffmat_band(:,:,ikxkyz), 3*(nmu+1)+1, ipiv(:,ikxkyz), ghrs, nvpa*nmu, info)
       do iv = 1, nvpa
             gvmutr(iv,:) = ghrs(nmu*(iv-1)+1 : nmu*iv)
       end do
       g(:,:,ikxkyz) = gvmutr
    end do

    allocate (flds(naky,nakx,-nzgrid:nzgrid,ntubes,nresponse_vpa))

    ! need to obtain phi^{n+1} and conservation terms using response matrix approach
    ! first get phi_inh^{n+1}
    call get_fields (g, phi, apar, dist='h')
    flds(:,:,:,:,1) = phi

    ! AVB: obtain phi^{n+1} from response matrix
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       ! all is indices inside ikxkyz super-index have same info
       ! no need to compute multiple times
       is = is_idx(kxkyz_lo,ikxkyz) ; if (is /= 1) cycle
       call lu_back_substitution (vpadiff_response(:,:,ikxkyz), vpadiff_idx(:,ikxkyz), flds(iky,ikx,iz,it,:))
       phi(iky,ikx,iz,it) = flds(iky,ikx,iz,it,1)
    end do

    call sum_allreduce (phi)

    g = g_in

    ! RHS is g^{***} + Ze/T*<phi^{n+1}>*F0 + 2*dt*nu*J0*F0*(vpa*upar+(v^2-3/2)*temp)
    ! first two terms added via g_to_h subroutine
    call g_to_h (g, phi, fphi)

    !allocate (tmp(nvpa,nmu))
    !deallocate (tmp, flds)

    ! now invert system to get h^{n+1}
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       is = is_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)

       ! using lapack:
       do iv = 1, nvpa
           ghrs(nmu*(iv-1)+1 : nmu*iv) = g(iv,:,ikxkyz)
       end do
       call zgbtrs('No transpose', nvpa*nmu, nmu+1, nmu+1, 1, cdiffmat_band(:,:,ikxkyz), 3*(nmu+1)+1, ipiv(:,ikxkyz), ghrs, nvpa*nmu, info)
       do iv = 1, nvpa
           gvmutr(iv,:) = ghrs(nmu*(iv-1)+1 : nmu*iv)
       end do
       g(:,:,ikxkyz) = gvmutr
    end do

    ! now get g^{n+1} from h^{n+1} and phi^{n+1}
    call g_to_h (g, phi, -fphi)

    deallocate (g_in)

  end subroutine advance_implicit_fp

  subroutine advance_vpadiff_implicit (phi, apar, g)

    use mp, only: sum_allreduce
    use finite_differences, only: tridag
    use linear_solve, only: lu_back_substitution
    use stella_time, only: code_dt
    use run_parameters, only: fphi
    use species, only: nspec, spec
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: nmu, nvpa
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa, vperp2
    use vpamu_grids, only: set_vpa_weights
    use kt_grids, only: naky, nakx
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use g_tofrom_h, only: g_to_h
    use gyro_averages, only: aj0v
    use fields, only: get_fields

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar
    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in out) :: g

    integer :: ikxkyz, iky, ikx, iz, it, is
    integer :: imu
    integer :: idx
    real, dimension (:,:), allocatable :: tmp
    complex, dimension (:,:,:,:,:), allocatable :: flds
    complex, dimension (:,:,:), allocatable :: g_in

    ! store input g for use later, as g will be overwritten below
    allocate (g_in(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    g_in = g

    ! since backwards difference in time, (I-dt*D)h_inh^{n+1} = g^{***}
    ! g = g^{***}.  tridag below inverts above equation to get h_inh^{n+1}
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          call tridag (1, aa_vpa(:,is), bb_vpa(:,ikxkyz), cc_vpa(:,is), g(:,imu,ikxkyz))
       end do
    end do

    allocate (flds(naky,nakx,-nzgrid:nzgrid,ntubes,nresponse_vpa))

    ! need to obtain phi^{n+1} and conservation terms using response matrix approach
    ! first get phi_inh^{n+1}
    call get_fields (g, phi, apar, dist='h')
    flds(:,:,:,:,1) = phi

    idx = 2
    ! get upar_inh^{n+1}
    if (momentum_conservation) then
       call get_upar (g, flds(:,:,:,:,idx:idx+nspec-1))
       idx = idx + nspec
    end if

    ! get temp_inh^{n+1}
    if (energy_conservation) call get_temp (g, flds(:,:,:,:,idx:idx+nspec-1))

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       ! all is indices inside ikxkyz super-index have same info
       ! no need to compute multiple times
       is = is_idx(kxkyz_lo,ikxkyz) ; if (is /= 1) cycle
       call lu_back_substitution (vpadiff_response(:,:,ikxkyz), vpadiff_idx(:,ikxkyz), &
            flds(iky,ikx,iz,it,:))
       phi = flds(iky,ikx,iz,it,1)
    end do
    call sum_allreduce (phi)

    g = g_in

    ! RHS is g^{***} + Ze/T*<phi^{n+1}>*F0 + 2*dt*nu*J0*F0*(vpa*upar+(v^2-3/2)*temp)
    ! first two terms added via g_to_h subroutine
    call g_to_h (g, phi, fphi)

    allocate (tmp(nvpa,nmu))

    if (momentum_conservation .or. energy_conservation) then
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          tmp = 2.0*code_dt*spec(is)%vnew(is) &
               *spread(maxwell_vpa(:,is),2,nmu)*spread(aj0v(:,ikxkyz)*maxwell_mu(1,iz,:,is),1,nvpa)
          if (momentum_conservation) &
               g(:,:,ikxkyz) = g(:,:,ikxkyz) + tmp*spread(vpa*flds(iky,ikx,iz,it,is+1),2,nmu)
          if (energy_conservation) &
               g(:,:,ikxkyz) = g(:,:,ikxkyz) &
               + tmp*(spread(vpa**2,2,nmu)+spread(vperp2(1,iz,:),1,nvpa)-1.5)*flds(iky,ikx,iz,it,idx+is-1)
       end do
    end if

    deallocate (tmp, flds)

    ! now invert system to get h^{n+1}
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          call tridag (1, aa_vpa(:,is), bb_vpa(:,ikxkyz), cc_vpa(:,is), g(:,imu,ikxkyz))
       end do
    end do

    ! now get g^{n+1} from h^{n+1} and phi^{n+1}
    call g_to_h (g, phi, -fphi)

    deallocate (g_in)

  end subroutine advance_vpadiff_implicit

  subroutine advance_mudiff_implicit (phi, apar, g)

    use mp, only: sum_allreduce
    use finite_differences, only: tridag
    use linear_solve, only: lu_back_substitution
    use stella_time, only: code_dt
    use run_parameters, only: fphi
    use species, only: nspec, spec
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: nmu, nvpa
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa, vperp2
    use vpamu_grids, only: set_vpa_weights
    use kt_grids, only: naky, nakx
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use dist_fn_arrays, only: kperp2
    use gyro_averages, only: aj0v, aj1v
    use g_tofrom_h, only: g_to_h
    use fields, only: get_fields
    use stella_geometry, only: bmag

    ! TMP FOR TESTING
!    use vpamu_grids, only: mu

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar
    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in out) :: g

    integer :: ikxkyz, iky, ikx, iz, it, is, ia
    integer :: iv
    integer :: idx

    ! TMP FOR TESTING
!    integer :: imu

    real, dimension (:,:), allocatable :: tmp
    complex, dimension (:,:,:,:,:), allocatable :: flds
    complex, dimension (:,:,:), allocatable :: g_in

    ! store input g for use later, as g will be overwritten below
    allocate (g_in(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    g_in = g

    ! since backwards difference in time, (I-dt*D)h_inh^{n+1} = g^{***}
    ! g = g^{***}.  tridag below inverts above equation to get h_inh^{n+1}
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       ! TMP FOR TESTING
!       do imu = 1, nmu
!          g(:,imu,ikxkyz) = maxwell_vpa*maxwell_mu(1,iz,imu)
!       end do
       do iv = 1, nvpa
          call tridag (1, aa_mu(iz,:,is), bb_mu(:,ikxkyz), cc_mu(iz,:,is), g(iv,:,ikxkyz))
       end do
       ! TMP FOR TESTING
!       iv = nvpa/2
!       do imu = 1, nmu
!          write (*,*) 'ggg', mu(imu), real(g(iv,imu,ikxkyz)), aimag(g(iv,imu,ikxkyz)), maxwell_vpa(iv,is)*maxwell_mu(1,iz,imu)
!       end do
    end do

    allocate (flds(naky,nakx,-nzgrid:nzgrid,ntubes,nresponse_mu))

    ! need to obtain phi^{n+1} and conservation terms using response matrix approach
    ! first get phi_inh^{n+1}
    call get_fields (g, phi, apar, dist='h')
    flds(:,:,:,:,1) = phi

    idx = 2
    ! get upar_inh^{n+1}
    if (momentum_conservation) then
       call get_uperp (g, flds(:,:,:,:,idx:idx+nspec-1))
       idx = idx + nspec
    end if

    ! get temp_inh^{n+1}
    if (energy_conservation) call get_temp_mu (g, flds(:,:,:,:,idx:idx+nspec-1))

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       ! all is indices inside ikxkyz super-index have same info
       ! no need to compute multiple times
       is = is_idx(kxkyz_lo,ikxkyz) ; if (is /= 1) cycle
       call lu_back_substitution (mudiff_response(:,:,ikxkyz), mudiff_idx(:,ikxkyz), &
            flds(iky,ikx,iz,it,:))
       phi = flds(iky,ikx,iz,it,1)
    end do
    call sum_allreduce (phi)

    g = g_in

    ! RHS is g^{***} + Ze/T*<phi^{n+1}>*F0 + ...
    ! first two terms added via g_to_h subroutine
    call g_to_h (g, phi, fphi)

    allocate (tmp(nvpa,nmu))

    ia = 1

    if (momentum_conservation .or. energy_conservation) then
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          tmp = 2.0*code_dt*spec(is)%vnew(is) &
               *spread(maxwell_vpa(:,is),2,nmu)*spread(maxwell_mu(1,iz,:,is),1,nvpa)
          if (momentum_conservation) &
               g(:,:,ikxkyz) = g(:,:,ikxkyz) + tmp*kperp2(iky,ikx,ia,iz) &
             * spread(vperp2(ia,iz,:)*aj1v(:,ikxkyz),1,nvpa)*(spec(is)%smz/bmag(ia,iz))**2 &
             * flds(iky,ikx,iz,it,is+1)
          if (energy_conservation) &
               g(:,:,ikxkyz) = g(:,:,ikxkyz) &
               + flds(iky,ikx,iz,it,idx+is-1)*tmp*spread(aj0v(:,ikxkyz),1,nvpa) &
             * (spread(vpa**2,2,nmu)+spread(vperp2(1,iz,:),1,nvpa)-1.5)
       end do
    end if

    deallocate (tmp, flds)

    ! now invert system to get h^{n+1}
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do iv = 1, nvpa
          call tridag (1, aa_mu(iz,:,is), bb_mu(:,ikxkyz), cc_mu(iz,:,is), g(iv,:,ikxkyz))
       end do
    end do

    ! now get g^{n+1} from h^{n+1} and phi^{n+1}
    call g_to_h (g, phi, -fphi)

    deallocate (g_in)

  end subroutine advance_mudiff_implicit

  subroutine advance_hyper_dissipation (g)

    use stella_time, only: code_dt
    use physics_flags, only: full_flux_surface, radial_variation
    use zgrid, only: nzgrid, ntubes, nztot
    use stella_layouts, only: vmu_lo
    use dist_fn_arrays, only: kperp2
    use kt_grids, only: naky, nakx
    use kt_grids, only: aky, akx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g

    integer :: ia
    integer :: ivmu
    real :: k2max

    if (full_flux_surface.or.radial_variation) then
       ! avoid spatially dependent kperp
       k2max = aky(nakx)**2 + aky(naky)**2
       ! add in hyper-dissipation of form dg/dt = -D*(k/kmax)^4*g
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          g(:,:,:,:,ivmu) = g(:,:,:,:,ivmu)/(1.+code_dt &
             * (spread(spread(spread(akx**2,1,naky)+spread(aky**2,2,nakx),3,nztot),4,ntubes)/k2max)**2*D_hyper)
       end do
    else
       k2max = maxval(kperp2)
       ia = 1
       ! add in hyper-dissipation of form dg/dt = -D*(k/kmax)^4*g
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          g(:,:,:,:,ivmu) = g(:,:,:,:,ivmu)/(1.+code_dt*(spread(kperp2(:,:,ia,:),4,ntubes)/k2max)**2*D_hyper)
       end do
    end if

  end subroutine advance_hyper_dissipation

end module dissipation
