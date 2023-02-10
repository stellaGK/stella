module adjoint_run

   implicit none

   public :: init_adjoint

   public :: allocate_adjoint_variables
   public :: get_save
   public :: get_adjoint_save
   public :: finish_adjoint, finish_pert_terms
   ! public :: allocate_derivative_terms
!  public :: method_manufactured_solutions

   private

contains

   subroutine init_adjoint(np)

      use adjoint_convergence, only: init_convergence
      use adjoint_p_derivatives, only: allocate_unpert

      implicit none

      integer, intent(in) :: np

      call init_convergence
      call allocate_adjoint_variables
      call allocate_unpert
      call allocate_derivative_terms(np)

   end subroutine init_adjoint

!   subroutine method_manufactured_solutions

!     use dist_fn_arrays, only: gnew
!     use fields, only: advance_fields, fields_updated
!     use vpamu_grids, only:  maxwell_vpa, maxwell_mu, maxwell_fac

!     use kt_grids, only: nakx, naky, nalpha, aky
!     use zgrid, only: nzgrid, ntubes, zed
!     use vpamu_grids, only: nvpa, nmu
!     use species, only: spec, nspec

!     use stella_layouts, only: vmu_lo
!     use stella_layouts, only: iv_idx, imu_idx, is_idx
!     use stella_layouts, only: kxkyz_lo
!     use stella_layouts, onlY: iz_idx, it_idx, ikx_idx, iky_idx

!     use mp, only: sum_allreduce
!     use vpamu_grids, only: integrate_species, integrate_vmu

!     use gyro_averages, only: aj0v, aj0x
!     use adjoint_field_arrays, only: omega_g
!     use stella_geometry, only: dbdzed

!     use vpamu_grids, only: vpa, mu

!     use constants, only: zi
!     use dist_fn_arrays, only: wdrifty_g

!     use fields_arrays, only: phi, apar

!     use adjoint_distfn_arrays, only: source_adjoint

!     use stella_time, only: code_dt
!     use dist_fn_arrays, only: wstar

!     implicit none

!     complex, dimension (:,:,:,:,:), allocatable :: maxwell
!     integer :: ivmu,iv,imu,is
!     integer :: iz, ia

!     real, dimension(:, :), allocatable :: g0
!     real :: tmp,  wgt
!     real, dimension (nspec) :: fac
!     integer :: ikxkyz, it, ikx, iky
!     real, dimension (:,:,:), allocatable :: gamm

!     complex, dimension (:,:,:), allocatable :: gamm2
!     complex, dimension (:,:,:), allocatable :: g1

!     logical :: adjoint = .True.
!     logical :: restart_time_step

!     allocate(gamm(naky, nakx, -nzgrid:nzgrid)) ; gamm = 0.0
!     allocate(maxwell(nalpha,-nzgrid:nzgrid,nmu,nvpa,nspec))

!     allocate(gamm2(naky,nakx,-nzgrid:nzgrid)) ; gamm2 = 0.0

!     if(.not.allocated(source_adjoint)) allocate(source_adjoint (naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
!     source_adjoint = 0.0

!     ia = 1
!     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!        iv = iv_idx(vmu_lo,ivmu)
!        imu = imu_idx(vmu_lo,ivmu)
!        is = is_idx(vmu_lo,ivmu)
!        do iz = -nzgrid, nzgrid
!           maxwell(ia,iz,imu,iv,is) =  maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is)
!           gnew(:,:,iz,:,ivmu) = maxwell(ia,iz,imu,iv,is)* exp(-zed(iz)**2)
!        end do
!     end do

!     fields_updated = .false.
!     call advance_fields (gnew, phi, apar, dist='gbar', adjoint=adjoint)

!     allocate (g0(nvpa, nmu))
!     do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
!        it = it_idx(kxkyz_lo, ikxkyz)
!        iky = iky_idx(kxkyz_lo, ikxkyz)
!        ikx = ikx_idx(kxkyz_lo, ikxkyz)
!        iz = iz_idx(kxkyz_lo, ikxkyz)
!        is = is_idx(kxkyz_lo, ikxkyz)
!        g0 = spread(aj0v(:, ikxkyz), 1, nvpa) &
!             * (spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is))**2
!        wgt = spec(is)%zt
!        call integrate_vmu(g0, iz, tmp)
!        gamm(iky, ikx, iz) = gamm(iky, ikx, iz) + tmp * wgt
!     end do
!     deallocate(g0)
!     call sum_allreduce(gamm)

!     allocate(g1(naky,nakx,vmu_lo%llim_proc:vmu_lo%ulim_alloc)) ; g1 = 0.0

!     do iz = -nzgrid, nzgrid
!        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!           iv = iv_idx(vmu_lo,ivmu)
!           imu = imu_idx(vmu_lo,ivmu)
!           is = is_idx(vmu_lo,ivmu)

!           g1(:,:,ivmu) = -zi* spread(aky,2,nakx)* aj0x(:,:,iz,ivmu) * wstar(ia,iz,ivmu)/(code_dt*spec(is)%zt) * &
!                maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
!        end do
!        call integrate_species (g1, iz, spec%zt, gamm2(:,:,iz), reduce_in=.false.)
!     end do

!     call sum_allreduce(gamm2)
!     deallocate(g1)

!     do is = 1, nspec
!        fac(is) = spec(is)%z * spec(is)%z * spec(is)%dens_psi0 / spec(is)%temp_psi0
!     end do
!     call sum_allreduce(fac)

!     gamm = gamm / sum(fac)
!     gamm2 = gamm2 /sum(fac)

!     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!        iv = iv_idx(vmu_lo,ivmu)
!        imu = imu_idx(vmu_lo,ivmu)
!        is = is_idx(vmu_lo,ivmu)
!        do iz = -nzgrid, nzgrid
!           source_adjoint(1,1,iz,ia,ivmu) = (omega_g(1,1) * ( maxwell(ia,iz,imu,iv,is) - aj0x(1,1,iz,ivmu)*spec(is)%z*spec(is)%dens_psi0 * gamm(1,1,iz)) &
!                - 2*zed(iz) * maxwell(ia,iz,imu,iv,is) * b_dot_grad_z(ia, iz) * vpa(iv) * spec(is)%stm_psi0  &
!                - zi*aky(1)* wdrifty_g(ia,iz, ivmu)*maxwell(ia,iz,imu,iv,is)/code_dt &
!                + aj0x(1,1,iz,ivmu)*gamm2(1,1,iz) *spec(is)%z * spec(is)%dens_psi0) * exp(-zed(iz)**2 ) !! minus sign in ydrift !!
!        end do
!     end do

!     deallocate(maxwell, gamm)

!     ! open (4, file='text_file_outputs/source.txt', status="unknown",action="write",position="append")
!     ! ia = 1
!     ! do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!     !    iv = iv_idx(vmu_lo,ivmu)
!     !    imu = imu_idx(vmu_lo,ivmu)
!     !    is = is_idx(vmu_lo,ivmu)
!     !    do iz = -nzgrid, nzgrid
!     !       write(4,*) 'iz',iz,'iv',iv,'imu',imu, source_adjoint(1,1,iz,ia,ivmu), gnew(1,1,iz,ia,ivmu)
!     !    end do
!     ! end do
!     ! close(4)
! !    deallocate(adjoint_out, adjoint_out_check)

!   end subroutine method_manufactured_solutions

  !! Allocate Arrays
   subroutine allocate_adjoint_variables

      use adjoint_distfn_arrays, only: gsave, g_omega, g_omega2
      use adjoint_distfn_arrays, only: lam_save
      use adjoint_field_arrays, only: phi_save
      use adjoint_field_arrays, only: chi_save

      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, naky
      use stella_layouts, only: vmu_lo

!    use adjoint_distfn_arrays, only: gyrog_save, gyrophi_save
      implicit none

      if (.not. allocated(g_omega)) then
         allocate (g_omega(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         g_omega = 0.0
      end if
      if (.not. allocated(g_omega2)) then
         allocate (g_omega2(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         g_omega2 = 0.
      end if

    !! variables needed for adjoint method
      if (.not. allocated(gsave)) then
         allocate (gsave(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         gsave = 0.0
      end if
      if (.not. allocated(lam_save)) then
         allocate (lam_save(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         lam_save = 0.0
      end if
      if (.not. allocated(phi_save)) then
         allocate (phi_save(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi_save = 0.0
      end if
      if (.not. allocated(chi_save)) then
         allocate (chi_save(naky, nakx, -nzgrid:nzgrid, ntubes))
         chi_save = 0.0
      end if

   end subroutine allocate_adjoint_variables

   subroutine allocate_derivative_terms(np)

      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx

      use adjoint_distfn_arrays, only: g_store
      use adjoint_field_arrays, only: q_store

      implicit none

      integer, intent(in) :: np

      if (.not. allocated(g_store)) then
         allocate (g_store(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         g_store = 0.0
      end if

      if (.not. allocated(q_store)) then
         allocate (q_store(naky, nakx, -nzgrid:nzgrid, ntubes))
         q_store = 0.0
      end if

   end subroutine allocate_derivative_terms

  !! Get Adjoint Variables
   subroutine get_save

      use adjoint_distfn_arrays, only: gsave, g_omega, g_omega2
      use adjoint_field_arrays, only: phi_save
      use dist_fn_arrays, only: gnew
      use fields_arrays, only: apar, phi
      use adjoint_field_arrays, only: omega_g

      use fields, only: advance_fields, fields_updated
      use stella_time, only: code_dt, code_time
      use kt_grids, only: nakx, naky

      use stella_layouts, only: vmu_lo

      use stella_diagnostics, only: omega_vs_time, navg
      use constants, only: zi

      use zgrid, only: nzgrid

      use mp, only: broadcast, proc0
      implicit none

      integer :: iky, ikx
      integer :: ivmu, iz

    !! Adjoint - Store frequency: omega_g = gamma + i*omega
      omega_g = -zi * sum(omega_vs_time, dim=1) / real(navg)

      call broadcast(omega_g)
      if (proc0) write (*, *) 'omega_g', omega_g(1, 1)

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do iky = 1, naky
            do ikx = 1, nakx
               gsave(iky, ikx, :, :, ivmu) = (3 * gnew(iky, ikx, :, :, ivmu) - 4 * g_omega(iky, ikx, :, :, ivmu) + &
                                             g_omega2(iky, ikx, :, :, ivmu)) / (2 * code_dt * omega_g(iky, ikx)) * exp(-omega_g(iky, ikx) * code_time)
            end do
         end do
      end do

      fields_updated = .false.
      call advance_fields(gsave, phi_save, apar, dist='gbar')

   end subroutine get_save

   subroutine get_adjoint_save

      use dist_fn_arrays, only: gnew
      use adjoint_distfn_arrays, only: lam_save
      use adjoint_field_arrays, only: chi_save
      use fields_arrays, only: apar

      use fields, only: advance_fields, fields_updated

      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid
      implicit none

      integer :: ivmu, iz
      logical :: adjoint = .True.

      lam_save = gnew

      fields_updated = .false.
      call advance_fields(lam_save, chi_save, apar, dist='gbar', adjoint=adjoint)

      lam_save = lam_save
      chi_save = chi_save
      call get_chi_lam

   end subroutine get_adjoint_save

   subroutine get_chi_lam

      use kt_grids, only: naky, nakx
      use vpamu_grids, only: nvpa, nmu
      use stella_layouts, only: kxkyz_lo

      use redistribute, only: scatter, gather
      use dist_redistribute, only: kxkyz2vmu

      use adjoint_distfn_arrays, only: lam_save

      use vpamu_grids, only: maxwell_mu, maxwell_vpa, maxwell_fac
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use zgrid, only: nzgrid, ntubes

      use stella_geometry, only: bmag
      use adjoint_field_arrays, only: chi_save
      implicit none

      complex, dimension(:, :, :), allocatable :: lamv, lamvmu

      integer :: iv, ia, ivmu, imu, iz, is

      allocate (lamv(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc))
      allocate (lamvmu(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc))

      ia = 1.
      ! do iz = -nzgrid,nzgrid
      !    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
      !       iv = iv_idx(vmu_lo,ivmu)
      !       imu = imu_idx(vmu_lo,ivmu)
      !       is = is_idx(vmu_lo,ivmu)
      !       lam_save(:,:,iz,:,ivmu) = lam_save(:,:,iz,:,ivmu)/&
      !            (maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)&
      !            *maxwell_fac(is))
      !    end do
      ! end do

      call scatter(kxkyz2vmu, lam_save, lamv)

      do iv = 1, nvpa
         lamvmu(iv, :, :) = lamv(nvpa - iv + 1, :, :)
      end do

      call gather(kxkyz2vmu, lamvmu, lam_save)

      do iz = -nzgrid, nzgrid
         lam_save(:, :, iz, ia, :) = lam_save(:, :, iz, ia, :)
         chi_save(:, :, iz, ia) = chi_save(:, :, iz, ia)
      end do

      deallocate (lamv)
      deallocate (lamvmu)

   end subroutine get_chi_lam

  !! Deallocate Arrays

   subroutine finish_adjoint

      use adjoint_p_derivatives, only: deallocate_p_derivatives
      use adjoint_convergence, only: deallocate_convergence

      implicit none

      call finish_pert_terms
      call deallocate_p_derivatives
      call deallocate_adjoint_variables
      call deallocate_convergence

   end subroutine finish_adjoint

   subroutine deallocate_adjoint_variables

      use adjoint_distfn_arrays, only: gsave, lam_save
      use adjoint_distfn_arrays, only: g_omega, g_omega2
      use adjoint_field_arrays, only: phi_save, chi_save

      implicit none

      deallocate (g_omega)
      deallocate (g_omega2)

      deallocate (gsave)
      deallocate (lam_save)
      deallocate (phi_save)
      deallocate (chi_save)

   end subroutine deallocate_adjoint_variables

   subroutine finish_pert_terms

      use adjoint_distfn_arrays, only: g_store
      use adjoint_field_arrays, only: q_store
      implicit none

      if (allocated(g_store)) deallocate (g_store)
      if (allocated(q_store)) deallocate (q_store)
   end subroutine finish_pert_terms

end module adjoint_run
