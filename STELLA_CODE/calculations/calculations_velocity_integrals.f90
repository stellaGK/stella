!###############################################################################
!                            VELOCITY INTEGRATIONS                          
!###############################################################################
! 
! Perform integrations over v-space. Integrating over the perpendicular velocity
! is equivalent to taking weighted sums over mu = [1 --> nmu] of the quantity <g>
! with weights wgts_mu(ia,iz,imu). Integrating over the parallel velocity is
! equivalent to taking weighted sums over v = [-nvgrid --> nvgrid] of the quantity
! <g> with weights wgts_vpa(iv). Integrating over the species is equivalent to
! summing the contributions of each species with the weights ws(is).
! 
! The interface "integrate_mu" integrates over the perpendicular velocity.
! The interface "integrate_vmu" integrates over both velocities.
! The interface "integrate_species" integrates over both velocities and the species.
! 
! Each interface has a subroutine to deal with complex and real distribution
! functions <g>. Moreover the integration can be performed locally by summing
! over all (imu,ivmu,is) or non-locally by summing over the ivmu points and then
! summing the contributions of all processors. Finally there are "single" and
! "block" routines where <g> in the block routines has indices (kx,ky) as well.
! 
!###############################################################################
module calculations_velocity_integrals
   
   ! Import the velocity integration weights for all subroutines
   use grids_velocity, only: nmu, nvpa
   use grids_velocity, only: wgts_mu, wgts_vpa
   use grids_velocity, only: wgts_mu_bare

   implicit none
   
   ! Integrations
   public :: integrate_vmu, integrate_vpa
   public :: integrate_mu, integrate_species
   public :: integrate_species_ffs, integrate_vmu_ffs
   public :: integrate_species_ffs_rm

   private

   ! Integrations over both velocities and species
   interface integrate_species
      module procedure integrate_species_vmu
      module procedure integrate_species_vmu_single
      module procedure integrate_species_vmu_single_real
      module procedure integrate_species_vmu_block_complex
      module procedure integrate_species_vmu_block_real
   end interface

   ! Integrations over both velocities
   interface integrate_vmu
      module procedure integrate_vmu_local_real
      module procedure integrate_vmu_local_complex
      module procedure integrate_vmu_vmulo_complex
      module procedure integrate_vmu_vmulo_ivmu_only_real
   end interface

   ! Integrations over the magnetic moment (perpendicular velocity)
   interface integrate_mu
      module procedure integrate_mu_local
      module procedure integrate_mu_nonlocal
   end interface

   ! Integrations over the parallel velocity
   interface integrate_vpa
      module procedure integrate_vpa_nonlocal 
   end interface

contains
   
!###############################################################################
!############################## INTEGRALS OVER MU ##############################
!###############################################################################
  
   !-------------- Full velocity space is local to each processor --------------
   subroutine integrate_mu_local(iz, g, total)

      use grids_species, only: nspec

      implicit none

      ! Arguments
      integer, intent(in) :: iz
      real, dimension(:, :), intent(in) :: g
      real, dimension(:), intent(out) :: total
      
      ! Local variables
      integer :: is, imu, ia
      
      !-------------------------------------------------------------------------

      ! Initialise sum
      total = 0.

      ! Assume we only have a single flux tube
      ia = 1
      
      ! Integrate over mu
      do is = 1, nspec
         do imu = 1, nmu
            total(is) = total(is) + wgts_mu(ia, iz, imu) * g(imu, is)
         end do
      end do

   end subroutine integrate_mu_local

   !---------------- Each processor has a number of ivmu points ----------------
   subroutine integrate_mu_nonlocal(iz, g, total)

      use mp, only: nproc, sum_reduce
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, imu_idx, iv_idx

      implicit none

      ! Arguments
      integer, intent(in) :: iz
      real, dimension(vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(:, :), intent(out) :: total

      ! Local variables
      integer :: is, imu, iv, ivmu, ia
      
      !-------------------------------------------------------------------------

      ! Initialise sum
      total = 0.

      ! Assume we only have a single flux tube
      ia = 1
      
      ! Iterate over the i[vpa, mu, s] points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Obtain the [ivpa, imu, is] indices from [ivmus]
         is = is_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         iv = iv_idx(vmu_lo, ivmu)
         
         ! Integrate over mu
         total(iv, is) = total(iv, is) + wgts_mu(ia, iz, imu) * g(ivmu)
         
      end do

      ! Each processor has a few [ivmu] points, so sum all calculations
      if (nproc > 1) call sum_reduce(total, 0)

   end subroutine integrate_mu_nonlocal
 
!###############################################################################
!############################## INTEGRALS OVER VPA #############################
!###############################################################################

   !---------------- Each processor has a number of ivmu points ----------------
   subroutine integrate_vpa_nonlocal(g_vs_ivmus, g_vs_mus)

      use mp, only: nproc, sum_reduce
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, imu_idx, iv_idx

      implicit none

      ! Arguments
      real, dimension(vmu_lo%llim_proc:), intent(in) :: g_vs_ivmus
      real, dimension(:, :), intent(out) :: g_vs_mus
      
      ! Local variables
      integer :: is, imu, iv, ivmus
      
      !-------------------------------------------------------------------------

      ! Initialize sum
      g_vs_mus = 0.

      ! Iterate over the i[vpa, mu, s] points
      do ivmus = vmu_lo%llim_proc, vmu_lo%ulim_proc

         ! Obtain the [ivpa, imu, is] indices from [ivmus]
         is = is_idx(vmu_lo, ivmus)
         imu = imu_idx(vmu_lo, ivmus)
         iv = iv_idx(vmu_lo, ivmus)

         ! Integrate over the parallel velocity
         g_vs_mus(imu, is) = g_vs_mus(imu, is) + wgts_vpa(iv) * g_vs_ivmus(ivmus)

      end do

      ! Each processor has a few [ivmu] points, so sum all calculations
      if (nproc > 1) call sum_reduce(g_vs_mus, 0)

   end subroutine integrate_vpa_nonlocal

!###############################################################################
!########################### INTEGRALS OVER (MU, VPA) ##########################
!###############################################################################
   
   !-------------- Full velocity space is local to each processor --------------
   subroutine integrate_vmu_local_real(g, iz, total)

      implicit none

      ! Arguments
      real, dimension(:, :), intent(in) :: g
      integer, intent(in) :: iz
      real, intent(out) :: total

      ! Local variables
      integer :: iv, imu, ia
      
      !-------------------------------------------------------------------------

      ! Initialize sum
      total = 0.

      ! Assume we only have a single flux tube
      ia = 1
      
      ! Integrate over mu and vpa
      do imu = 1, nmu
         do iv = 1, nvpa
            total = total + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g(iv, imu)
         end do
      end do

   end subroutine integrate_vmu_local_real

   !-------------- Full velocity space is local to each processor --------------
   subroutine integrate_vmu_local_complex(g, iz, total)

      implicit none

      ! Arguments
      complex, dimension(:, :), intent(in) :: g
      integer, intent(in) :: iz
      complex, intent(out) :: total

      ! Local variables
      integer :: iv, imu, ia
      
      !-------------------------------------------------------------------------

      ! Initialize sum
      total = 0.

      ! Assume we only have a single flux tube
      ia = 1
      
      ! Integrate over mu and vpa
      do imu = 1, nmu
         do iv = 1, nvpa
            total = total + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g(iv, imu)
         end do
      end do

   end subroutine integrate_vmu_local_complex

   !---------------- Each processor has a number of ivmu points ----------------
   subroutine integrate_vmu_vmulo_complex(g, weights, total)

      use mp, only: sum_allreduce
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(:), intent(in) :: weights
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: total

      ! Local variables
      integer :: ivmu, iv, iz, is, imu, ia
      
      !-------------------------------------------------------------------------

      ! Initialize sum
      total = 0.

      ! Assume we only have a single flux tube
      ia = 1
      
      ! Iterate over the i[vpa, mu, s] points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Obtain the [ivpa, imu, is] indices from [ivmus]
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         
         ! Integrate over mu, vpa, species
         do iz = -nzgrid, nzgrid
            total(:, :, iz, :, is) = total(:, :, iz, :, is) + &
               wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g(:, :, iz, :, ivmu) * weights(is)
         end do
         
      end do

      ! Each processor has a few [ivmu] points, so sum all calculations
      call sum_allreduce(total)

   end subroutine integrate_vmu_vmulo_complex

   !---------------- Each processor has a number of ivmu points ----------------
   subroutine integrate_vmu_vmulo_ivmu_only_real(g, ia, iz, total)

      use mp, only: sum_allreduce
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

      implicit none

      ! Arguments
      real, dimension(vmu_lo%llim_proc:), intent(in) :: g
      integer, intent(in) :: ia, iz
      real, dimension(:), intent(out) :: total

      ! Local variables
      integer :: ivmu, iv, is, imu
      
      !-------------------------------------------------------------------------

      ! Initialize sum
      total = 0.

      ! Iterate over the i[vpa, mu, s] points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Obtain the [ivpa, imu, is] indices from [ivmus]
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         
         ! Integrate over mu and vpa
         total(is) = total(is) + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g(ivmu)
         
      end do

      ! Each processor has a few [ivmu] points, so sum all calculations
      call sum_allreduce(total)

   end subroutine integrate_vmu_vmulo_ivmu_only_real

   !---------------- Each processor has a number of ivmu points ----------------
   subroutine integrate_species_vmu(g, weights, total, ia_in)

      use mp, only: sum_allreduce
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use grids_z, only: nzgrid, ntubes

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(:), intent(in) :: weights
      integer, intent(in), optional :: ia_in
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: total

      ! Local variables
      integer :: ivmu, iv, it, iz, is, imu, ia
      
      !-------------------------------------------------------------------------
      
      ! Initialize sum
      total = 0.
      
      ! Number of field lines
      if (present(ia_in)) then
         ia = ia_in
      else
         ia = 1
      end if

      ! Iterate over the i[vpa, mu, s] points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Obtain the [ivpa, imu, is] indices from [ivmus]
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         
         ! Integrate over mu, vpa, species
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               total(:, :, iz, it) = total(:, :, iz, it) + &
                  wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g(:, :, iz, it, ivmu) * weights(is)
            end do
         end do
         
      end do

      ! Each processor has a few [ivmu] points, so sum all calculations
      call sum_allreduce(total)

   end subroutine integrate_species_vmu

   !---------------- Each processor has a number of ivmu points ----------------
   subroutine integrate_species_vmu_single(g, iz, weights, total, ia_in, reduce_in)

      use mp, only: sum_allreduce
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

      implicit none

      ! Arguments
      complex, dimension(vmu_lo%llim_proc:), intent(in) :: g
      integer, intent(in) :: iz
      real, dimension(:), intent(in) :: weights
      complex, intent(out) :: total
      integer, intent(in), optional :: ia_in
      logical, intent(in), optional :: reduce_in

      ! Local variables
      integer :: ivmu, iv, is, imu, ia
      logical :: reduce
      
      !-------------------------------------------------------------------------

      ! Initialize sum
      total = 0.

      ! Number of field lines
      if (present(ia_in)) then
         ia = ia_in
      else
         ia = 1
      end if
      
      ! Sum over processors or not
      if (present(reduce_in)) then
         reduce = reduce_in
      else
         reduce = .true.
      end if

      ! Iterate over the i[vpa, mu, s] points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Obtain the [ivpa, imu, is] indices from [ivmus]
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         
         ! Integrate over mu, vpa, species
         total = total + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g(ivmu) * weights(is)
         
      end do

      ! Each processor has a few [ivmu] points, so sum all calculations
      if (reduce) call sum_allreduce(total)

   end subroutine integrate_species_vmu_single

   !---------------- Each processor has a number of ivmu points ----------------
   subroutine integrate_species_vmu_single_real(g, iz, weights, total, ia_in, reduce_in)

      use mp, only: sum_allreduce
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

      implicit none

      ! Arguments
      real, dimension(vmu_lo%llim_proc:), intent(in) :: g
      integer, intent(in) :: iz
      real, dimension(:), intent(in) :: weights
      real, intent(out) :: total
      integer, intent(in), optional :: ia_in
      logical, intent(in), optional :: reduce_in

      ! Local variables
      integer :: ivmu, iv, is, imu, ia
      logical :: reduce
      
      !-------------------------------------------------------------------------

      ! Initialize sum
      total = 0.

      ! Number of field lines
      if (present(ia_in)) then
         ia = ia_in
      else
         ia = 1
      end if
      
      ! Sum over processors or not
      if (present(reduce_in)) then
         reduce = reduce_in
      else
         reduce = .true.
      end if

      ! Iterate over the i[vpa, mu, s] points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Obtain the [ivpa, imu, is] indices from [ivmus]
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         
         ! Integrate over mu, vpa, species
         total = total + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g(ivmu) * weights(is)
         
      end do

      ! Each processor has a few [ivmu] points, so sum all calculations
      if (reduce) call sum_allreduce(total)

   end subroutine integrate_species_vmu_single_real

   !---------------- Each processor has a number of ivmu points ----------------
   subroutine integrate_species_vmu_block_complex(g, iz, weights, pout, ia_in, reduce_in)

      use mp, only: sum_allreduce
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

      implicit none

      ! Arguments
      complex, dimension(:, :, vmu_lo%llim_proc:), intent(in) :: g
      integer, intent(in) :: iz
      integer, intent(in), optional :: ia_in
      logical, intent(in), optional :: reduce_in
      real, dimension(:), intent(in) :: weights
      complex, dimension(:, :), intent(out) :: pout

      ! Local variables
      integer :: ivmu, iv, is, imu, ia
      logical :: reduce
      
      !-------------------------------------------------------------------------

      ! Initialize sum
      pout = 0.

      ! Number of field lines
      if (present(ia_in)) then
         ia = ia_in
      else
         ia = 1
      end if
      
      ! Sum over processors or not
      if (present(reduce_in)) then
         reduce = reduce_in
      else
         reduce = .true.
      end if

      ! Iterate over the i[vpa, mu, s] points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Obtain the [ivpa, imu, is] indices from [ivmus]
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         
         ! Integrate over mu, vpa, species
         pout = pout + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g(:, :, ivmu) * weights(is)
         
      end do

      ! Each processor has a few [ivmu] points, so sum all calculations
      if (reduce) call sum_allreduce(pout)

   end subroutine integrate_species_vmu_block_complex

   !---------------- Each processor has a number of ivmu points ----------------
   subroutine integrate_species_vmu_block_real(g, iz, weights, pout, ia_in, reduce_in)

      use mp, only: sum_allreduce
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

      implicit none

      ! Arguments
      real, dimension(:, :, vmu_lo%llim_proc:), intent(in) :: g
      integer, intent(in) :: iz
      integer, intent(in), optional :: ia_in
      logical, intent(in), optional :: reduce_in
      real, dimension(:), intent(in) :: weights
      real, dimension(:, :), intent(out) :: pout

      ! Local variables
      integer :: ivmu, iv, is, imu, ia
      logical :: reduce
      
      !-------------------------------------------------------------------------

      ! Initialize sum
      pout = 0.

      ! Number of field lines
      if (present(ia_in)) then
         ia = ia_in
      else
         ia = 1
      end if
      
      ! Sum over processors or not
      if (present(reduce_in)) then
         reduce = reduce_in
      else
         reduce = .true.
      end if

      ! Iterate over the i[vpa, mu, s] points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Obtain the [ivpa, imu, is] indices from [ivmus]
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         
         ! Integrate over mu, vpa, species
         pout = pout + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g(:, :, ivmu) * weights(is)
         
      end do

      ! Each processor has a few [ivmu] points, so sum all calculations
      if (reduce) call sum_allreduce(pout)

   end subroutine integrate_species_vmu_block_real

   !---------------- Each processor has a number of ivmu points ----------------
   subroutine integrate_species_ffs(g, weights, pout, reduce_in)

      use mp, only: sum_allreduce
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

      implicit none

      ! Arguments
      complex, dimension(:, :, vmu_lo%llim_proc:), intent(in) :: g
      logical, intent(in), optional :: reduce_in
      real, dimension(:), intent(in) :: weights
      complex, dimension(:, :), intent(out) :: pout

      ! Local variables
      integer :: ivmu, iv, is, imu
      logical :: reduce
      
      !-------------------------------------------------------------------------

      ! Initialize sum
      pout = 0.

      ! Sum over processors or not
      if (present(reduce_in)) then
         reduce = reduce_in
      else
         reduce = .true.
      end if

      ! Iterate over the i[vpa, mu, s] points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Obtain the [ivpa, imu, is] indices from [ivmus]
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         
         ! Integrate over mu, vpa, species
         pout = pout + 2.0 * wgts_mu_bare(imu) * wgts_vpa(iv) * g(:, :, ivmu) * weights(is)
         
      end do

      ! Each processor has a few [ivmu] points, so sum all calculations
      if (reduce) call sum_allreduce(pout)

   end subroutine integrate_species_ffs

   !---------------- Each processor has a number of ivmu points ----------------
   subroutine integrate_species_ffs_rm(g, weights, pout, reduce_in)
      use mp, only: sum_allreduce
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

      implicit none

      ! Arguments
      complex, dimension(vmu_lo%llim_proc:), intent(in) :: g
      logical, intent(in), optional :: reduce_in
      real, dimension(:), intent(in) :: weights
      complex, intent(out) :: pout
      
      ! Local variables
      integer :: ivmu, iv, is, imu
      logical :: reduce
      
      !-------------------------------------------------------------------------

      ! Initialize sum
      pout = 0.
      
      ! Sum over processors or not
      if (present(reduce_in)) then
         reduce = reduce_in
      else
         reduce = .true.
      end if

      ! Iterate over the i[vpa, mu, s] points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Obtain the [ivpa, imu, is] indices from [ivmus]
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         
         ! Integrate over mu, vpa, species
         pout = pout + 2.0 * wgts_mu_bare(imu) * wgts_vpa(iv) * g(ivmu) * weights(is)
         
      end do

      ! Each processor has a few [ivmu] points, so sum all calculations
      if (reduce) call sum_allreduce(pout)

   end subroutine integrate_species_ffs_rm

   !-------------- Full velocity space is local to each processor --------------
   subroutine integrate_vmu_ffs(g, weights, ia, iz, pout, reduce_in)

      use mp, only: sum_allreduce
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

      implicit none

      ! Arguments
      complex, dimension(vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(:), intent(in) :: weights
      integer, intent(in) :: ia, iz
      complex, dimension(:), intent(out) :: pout
      logical, intent(in), optional :: reduce_in

      ! Local variables
      integer :: ivmu, iv, is, imu
      logical :: reduce
      
      !-------------------------------------------------------------------------

      ! Initialize sum
      pout = 0.

      ! Sum over processors or not
      if (present(reduce_in)) then
         reduce = reduce_in
      else
         reduce = .true.
      end if

      ! NB: for FFS, assume that there is only one flux annulus
      ! the inclusion of the Maxwellian term below is due to the fact that
      ! g/F is evolved for FFS
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         pout(is) = pout(is) + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g(ivmu) * weights(is)
      end do

      ! Each processor has a few [ivmu] points, so sum all calculations
      if (reduce) call sum_allreduce(pout)

   end subroutine integrate_vmu_ffs
   
end module calculations_velocity_integrals
