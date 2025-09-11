!###############################################################################
!                 CONVERT DISTRIBUTION FUNCTION BETWEEN (F,G,H)                 
!###############################################################################
! 
! The distribution function can be given as,
!     - The first-order turbulent distribution function delta f_s
!     - The first-order gyro-averaged turbulent distribution function g_s
!     - The non-adiabatic part h_s of the first-order turbulent distribution function
! 
!                                  MATHEMATICS                                  
! 
! One can convert between f, g and h using:
!     delta f_s = h_s - Z_s/T_s * phi * F_s
!     g_s = h_s - Z_s/T_s * <phi>_theta * F_s
!     h_s = delta f_s + Z_s/T_s * phi * F_s = g_s + Z_s/T_s * <phi>_theta * F_s
! 
! For electromagnetic simulations we redefine g as gbar
!     <chi>_theta = J_0*phi - 2*vpa*J_0*apar + 4*mu*(T_s/Z_s)*bpar*J_1/b_s
!     gbar = h - (Z_s/T_s) * <chi>_theta * F_s
!     gbar = h - (Z_s/T_s)*J_0*phi*F_s + 2*(Z_s/T_s)*J_0*vpa*apar*F_s - 4*mu*bpar*F_s*J_1/b_s
!     g = h - Z_s/T_s * <phi>_theta * F_s
!     gbar = g + 2*(Z_s/T_s)*J_0*vpa*apar*F_s - 4*mu*bpar*F_s*J_1/b_s
! 
! Note that a factor 2 appears in front of <apar> and <bpar> due to the 
! normalisation of the velocities with v_th = sqrt(2T/m).
! 
! Here F_s is a Maxwellian distribution. J_0 is the zeroth order Bessel function
! of the first kind, coming from the gyro-averages, i.e. <phi>_theta --> J_0 * phi.
! 
! Note that f, g and h are given in guiding-center coordinates (R,mu,vpa,t)
! whereas the electrostatic potential phi is given in (r,t).
! 
!###############################################################################
module calculations_tofrom_ghf
   
   ! Make routines available to other modules
   public :: gbar_to_g
   public :: g_to_h
   public :: g_to_f

   private

   interface gbar_to_g
      module procedure gbar_to_g_kxkyz
      module procedure gbar_to_g_1d_vpa
      module procedure gbar_to_g_vmu
      module procedure gbar_to_g_vmu_single
   end interface

   interface g_to_h
      module procedure g_to_h_kxkyz
      module procedure g_to_h_vmu
      module procedure g_to_h_vmu_single 
   end interface

   interface g_to_f
      module procedure g_to_f_kxkyz
      module procedure g_to_f_vmu
   end interface g_to_f

contains

!###############################################################################
!############################## CONVERT GBAR TO G ##############################
!###############################################################################
! 
! Recall the definition of <gbar>
!    <gbar> = <g> + 2 * (Z_s/T_s) * J_0 * vpa * <apar> * F_s
! 
! To obtain <g> we distract [2*(Z_s/T_s)*J_0*vpa*apar*F_s] from <gbar>
!    <gyro_averaged_field> = 2*(Z_s/T_s)*J_0*vpa*<apar>*F_s
!    g(kx,ky,z,mu,vpa,s) = gbar(kx,ky,z,mu,vpa,s) - gyro_averaged_field
! 
! To convert <g> to <gbar> simply supply <g> as an argument instead of <gbar>, 
! and set <facapar> to -1 which will result in:
!    gbar(kx,ky,z,mu,vpa,s) = g(kx,ky,z,mu,vpa,s) + gyro_averaged_field
! 
!###############################################################################

   !****************************************************************************
   !                 Convert <gbar> to <g> using (mu,vpa,kxkyzs)                
   !****************************************************************************
   subroutine gbar_to_g_kxkyz(g, apar, facapar)

      ! Parallelisation
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average
      
      ! Grids
      use grids_species, only: spec
      use grids_z, only: nzgrid
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: maxwell_vpa, maxwell_mu, vpa
      use parameters_numerical, only: maxwellian_normalization
      
      implicit none

      ! Arguments
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: apar
      real, intent(in) :: facapar

      ! Local variables
      integer :: ikxkyz, iz, it, iky, ikx, is, ia
      complex, dimension(:, :), allocatable :: field, gyro_averaged_field
      
      !-------------------------------------------------------------------------

      ! Allocate local arrays
      allocate (field(nvpa, nmu))
      allocate (gyro_averaged_field(nvpa, nmu))

      ! Assume we only have one field line
      ia = 1
      
      ! Iterate over the (kx,ky,z,mu,vpa,s) grid
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iky = iky_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         
         ! Calculate <gyro_averaged_field> = 2*(Z_s/T_s)*J_0*vpa*<apar>*F_s 
         ! First calculate [2*apar*(Z_s/T_s)*vpa], then add F_s, and then J_0
         field = 2.0 * facapar * apar(iky, ikx, iz, it) * spec(is)%zt * spec(is)%stm_psi0 * spread(vpa, 2, nmu)
         if (.not. maxwellian_normalization) then
            field = field * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa)
         end if
         call gyro_average(field, ikxkyz, gyro_averaged_field)
         
         ! Calculate <g>  = <gbar> - 2*(Z_s/T_s)*J_0*vpa*<apar>*F_s 
         g(:, :, ikxkyz) = g(:, :, ikxkyz) - gyro_averaged_field
         
      end do
      
      ! Deallocate local arrays
      deallocate (field, gyro_averaged_field)
      
   end subroutine gbar_to_g_kxkyz

   !****************************************************************************
   ! Convert <gbar> to <g> using (mu,vpa,kxkyzs) for a specific imu and ikxkyzs 
   !****************************************************************************
   subroutine gbar_to_g_1d_vpa(g, apar, imu, ikxkyz, facapar)

      ! Parallelisation
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: iz_idx, is_idx
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average
      
      ! Grids
      use grids_species, only: spec
      use grids_velocity, only: nvpa
      use grids_velocity, only: maxwell_vpa, maxwell_mu, vpa
      use parameters_numerical, only: maxwellian_normalization
      
      implicit none

      ! Arguments
      complex, dimension(:), intent(in out) :: g
      complex, intent(in) :: apar
      integer, intent(in) :: imu, ikxkyz
      real, intent(in) :: facapar

      ! Local variables
      integer :: iz, is, ia
      complex, dimension(:), allocatable :: field, gyro_averaged_field
      
      !-------------------------------------------------------------------------

      ! Allocate local arrays
      allocate (field(nvpa))
      allocate (gyro_averaged_field(nvpa))

      ! Assume we only have one field line
      ia = 1
      
      ! Get the specific (kx,ky,z,s,mu) point
      iz = iz_idx(kxkyz_lo, ikxkyz)
      is = is_idx(kxkyz_lo, ikxkyz)
      
      ! Calculate <gyro_averaged_field> = 2*(Z_s/T_s)*J_0*vpa*<apar>*F_s 
      ! First calculate [2*apar*(Z_s/T_s)*vpa], then add F_s, and then J_0
      field = 2.0 * facapar * apar * spec(is)%zt * spec(is)%stm_psi0 * vpa
      if (.not. maxwellian_normalization) then
         field = field * maxwell_vpa(:, is) * maxwell_mu(ia, iz, imu, is)
      end if
      call gyro_average(field, imu, ikxkyz, gyro_averaged_field)
      
      ! Calculate <g>  = <gbar> - 2*(Z_s/T_s)*J_0*vpa*<apar>*F_s 
      g = g - gyro_averaged_field
      
      ! Deallocate local arrays
      deallocate (field, gyro_averaged_field)
      
   end subroutine gbar_to_g_1d_vpa


   !****************************************************************************
   !               Convert <gbar> to <g> using (kx,ky,z,ivpamus)                
   !****************************************************************************
   subroutine gbar_to_g_vmu(g, apar, facapar)

      use grids_z, only: nzgrid
      use parallelisation_layouts, only: vmu_lo

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: apar
      real, intent(in) :: facapar

      ! Local variables
      integer :: ivmu
      
      !-------------------------------------------------------------------------
   
      ! Convert <gbar> to <g> for each ivpamus point
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         call gbar_to_g_vmu_single(ivmu, g(:, :, :, :, ivmu), apar, facapar)
      end do

   end subroutine gbar_to_g_vmu

   !****************************************************************************
   ! Convert <gbar> to <g> using (kx,ky,z,ivpamus) for a specific ivpamus point 
   !****************************************************************************
   subroutine gbar_to_g_vmu_single(ivmu, g0, apar, facapar)

      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      
      ! Calculations
      use calculations_kxky, only: multiply_by_rho
      use calculations_gyro_averages, only: gyro_average
      
      ! Grids
      use grids_species, only: spec
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, nakx
      use grids_velocity, only: vpa
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use parameters_numerical, only: maxwellian_normalization
      
      implicit none

      ! Arguments
      integer, intent(in) :: ivmu
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: g0
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: apar
      real, intent(in) :: facapar

      ! Local variables
      integer :: iv, imu, is
      integer :: it, iz, ia
      complex, dimension(:, :), allocatable :: field, gyro_averaged_field
      
      !-------------------------------------------------------------------------
      
      ! Allocate local arrays
      allocate (field(naky, nakx))
      allocate (gyro_averaged_field(naky, nakx))
      
      ! Get the indices of this ivpamus point
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      ! Assume we only have one field line
      ia = 1
      
      ! Iterate over the (it,iz) points
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
         
            ! Calculate <gyro_averaged_field> = 2*(Z_s/T_s)*J_0*vpa*<apar>*F_s 
            ! First calculate [2*apar*(Z_s/T_s)*vpa], then add F_s, and then J_0
            field = 2.0 * spec(is)%zt * spec(is)%stm_psi0 * vpa(iv) * facapar * apar(:, :, iz, it)
            if (.not. maxwellian_normalization) then
               field = field * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
            end if
            call gyro_average(field, iz, ivmu, gyro_averaged_field)
            
            ! Calculate <g>  = <gbar> - 2*(Z_s/T_s)*J_0*vpa*<apar>*F_s
            g0(:, :, iz, it) = g0(:, :, iz, it) - gyro_averaged_field
            
         end do
      end do
      
      ! Deallocate local arrays
      deallocate (field, gyro_averaged_field)

   end subroutine gbar_to_g_vmu_single
   
   
!###############################################################################
!################################ CONVERT G TO H ###############################
!###############################################################################
! 
! Recall the definition of <g> and <h>
!     <g> = <h> - Z_s/T_s * <phi>_theta * F_s
!     <h> = <g> + Z_s/T_s * <phi>_theta * F_s
! 
! To obtain <h> we add [Z_s/T_s * <phi>_theta * F_s] to <g>
!    <gyro_averaged_field> = Z_s/T_s * <phi>_theta * F_s
!    h(kx,ky,z,mu,vpa,s) = g(kx,ky,z,mu,vpa,s) + <gyro_averaged_field>
! 
! To convert <h> to <g> simply supply <h> as an argument instead of <g>,
! and set <facphi> to -1 which will result in:
!    g(kx,ky,z,mu,vpa,s) = h(kx,ky,z,mu,vpa,s) - <gyro_averaged_field>
! 
!###############################################################################



   !****************************************************************************
   !                 Convert <gbar> to <g> using (mu,vpa,kxkyzs)                
   !****************************************************************************
   subroutine g_to_h_kxkyz(g, phi, bpar, facphi)

      ! Parallelisation
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: iky_idx, ikx_idx
      use parallelisation_layouts, only: iz_idx, it_idx, is_idx
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average
      use calculations_gyro_averages, only: gyro_average_j1
      
      ! Flags
      use parameters_numerical, only: maxwellian_normalization
      use parameters_physics, only: include_bpar
      
      ! Grids
      use grids_species, only: spec
      use grids_z, only: nzgrid
      use grids_velocity, only: nvpa, nmu, mu
      use grids_velocity, only: maxwell_vpa, maxwell_mu
      
      implicit none

      ! Arguments
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
      real, intent(in) :: facphi

      ! Local variables
      integer :: ikxkyz, iz, it, iky, ikx, is, ia
      complex, dimension(:, :), allocatable :: field, gyro_averaged_field
      real :: facbpar
      
      !-------------------------------------------------------------------------
      
      allocate (field(nvpa, nmu))
      allocate (gyro_averaged_field(nvpa, nmu))

      ! Assume we only have one field line
      ia = 1
      
      ! Iterate over the (kx,ky,z,mu,vpa,s) grid
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iky = iky_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         
         ! Calculate <gyro_averaged_field> = Z_s/T_s * <phi>_theta * F_s
         ! First calculate [(Z_s/T_s)*<phi>], then add F_s and J_0
         field = facphi * phi(iky, ikx, iz, it) * spec(is)%zt
         if (.not. maxwellian_normalization) then
            field = field * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa)
         end if
         call gyro_average(field, ikxkyz, gyro_averaged_field)
         
         ! Calculate <h> = <g> + (Z_s/T_s)*J_0*phi*F_s
         g(:, :, ikxkyz) = g(:, :, ikxkyz) + gyro_averaged_field
         
      end do
      
      
      ! Add electromagnetic terms: 4*mu*J_1*<bpar>/b_s
      if (include_bpar) then
      
         ! This factor determines whether we add or substract the electromagnetic term
         facbpar = facphi
         
         ! Iterate over the (it,iz) points
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            
            ! Calculate <gyro_averaged_field> = 4*mu*<bpar>*Fs*J_1/b_s
            ! First calculate [4*mu*<bpar>], then add F_s, and then J_1/b_s
            field = 4.0 * facbpar * spread(mu, 1, nvpa) * bpar(iky, ikx, iz, it)
            if (.not. maxwellian_normalization) then
               field = field * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa)
            end if
            call gyro_average_j1(field, ikxkyz, gyro_averaged_field)
            
            ! Calculate <h> = <g> + Z_s/T_s*J_0*phi*Fs + 4*mu*<bpar>*Fs*J_1/b_s
            g(:, :, ikxkyz) = g(:, :, ikxkyz) + gyro_averaged_field
            
         end do
         
      end if
      
      ! Deallocate local arrays
      deallocate (field, gyro_averaged_field)

   end subroutine g_to_h_kxkyz
   
   !****************************************************************************
   !                 Convert <g> to <h> using (kx,ky,z,ivpamus)                 
   !****************************************************************************
   subroutine g_to_h_vmu(g, phi, bpar, facphi, phi_corr)

      use grids_z, only: nzgrid
      use parallelisation_layouts, only: vmu_lo

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
      complex, dimension(:, :, -nzgrid:, :), optional, intent(in) :: phi_corr
      real, intent(in) :: facphi

      ! Local variables
      integer :: ivmu
      
      !-------------------------------------------------------------------------

      ! Convert <g> to <h> for each ivpamus point
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         call g_to_h_vmu_single(ivmu, g(:, :, :, :, ivmu), phi, bpar, facphi, phi_corr)
      end do

   end subroutine g_to_h_vmu

   !****************************************************************************
   !   Convert <g> to <h> using (kx,ky,z,ivpamus) for a specific ivpamus point  
   !****************************************************************************
   subroutine g_to_h_vmu_single(ivmu, g0, phi, bpar, facphi, phi_corr)

      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      use arrays_gyro_averages, only: aj0x, aj1x
      use calculations_kxky, only: multiply_by_rho
      
      ! Flags
      use parameters_physics, only: radial_variation
      use parameters_physics, only: include_bpar
      use parameters_numerical, only: maxwellian_normalization
      
      ! Geometry
      use geometry, only: bmag, dBdrho
      use arrays, only: kperp2, dkperp2dr
      
      ! Grids
      use grids_species, only: spec
      use grids_kxky, only: naky, nakx
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: vpa, vperp2, mu
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac

      implicit none

      ! Arguments
      integer, intent(in) :: ivmu
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: g0
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
      real, intent(in) :: facphi
      complex, dimension(:, :, -nzgrid:, :), optional, intent(in) :: phi_corr

      ! Local variables
      real :: facbpar
      integer :: iv, imu, is
      integer :: it, iz, ia
      complex, dimension(:, :), allocatable :: field, gyro_averaged_field, g0k
      
      !-------------------------------------------------------------------------

      ! Allocate local arrays
      allocate (field(naky, nakx))
      allocate (gyro_averaged_field(naky, nakx))
      if (radial_variation) allocate (g0k(naky, nakx))
      
      ! Get the indices of this ivpamus point
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      ! Assume we only have one field line
      ia = 1
      
      ! Iterate over the (it,iz) points
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
         
            ! Calculate <gyro_averaged_field> = Z_s/T_s * <phi>_theta * F_s
            ! First calculate [(Z_s/T_s)*<phi>] and add F_s
            field = spec(is)%zt * facphi * phi(:, :, iz, it)
            if (.not. maxwellian_normalization) then
               field = field * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
            end if
            
            ! Add radial variation corrections
            if (radial_variation .and. present(phi_corr)) then
               g0k = field * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                  - spec(is)%fprim - 2.0 * dBdrho(iz) * mu(imu) &
                  - 0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                  * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                  * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)))

               call multiply_by_rho(g0k)

               field = field + g0k + phi_corr(:, :, iz, it) * spec(is)%zt * facphi &
                  * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
            end if
            
            ! Finally add <...>_theta, which is equivalent to adding J_0
            call gyro_average(field, iz, ivmu, gyro_averaged_field)
            
            ! Calculate <h> = <g> + Z_s/T_s * <phi>_theta * F_s
            g0(:, :, iz, it) = g0(:, :, iz, it) + gyro_averaged_field
            
         end do
      end do
      
      ! Add electromagnetic terms: 4*mu*<bpar>*Fs*J_1/b_s
      if (include_bpar) then
      
         ! This factor determines whether we add or substract the electromagnetic term
         facbpar = facphi
         
         ! Iterate over the (it,iz) points
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
            
               ! Calculate <gyro_averaged_field> = 4*mu*<bpar>*Fs*J_1/b_s
               ! First calculate [4*mu*<bpar>], then add F_s, and then J_1/b_s
               field = 4.0 * mu(imu) * facbpar * bpar(:,:,iz,it) 
               if (.not. maxwellian_normalization) then
                  field = field * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
               end if
               call gyro_average_j1(field, iz, ivmu, gyro_averaged_field)
               
               ! Calculate <h> = <g> + Z_s/T_s*J_0*phi*Fs + 4*mu*<bpar>*Fs*J_1/b_s
               g0(:, :, iz, it) = g0(:, :, iz, it) + gyro_averaged_field
               
            end do
         end do
         
      end if
      
      ! Deallocate local arrays
      deallocate (field, gyro_averaged_field)
      if (allocated(g0k)) deallocate (g0k)

   end subroutine g_to_h_vmu_single

!###############################################################################
!################################ CONVERT G TO F ###############################
!###############################################################################
! 
! Recall the definitions of the distribution functions
!     delta f_s = h_s - Z_s/T_s * phi * F_s
!     g_s = h_s - Z_s/T_s * <phi>_theta * F_s
!     h_s = delta f_s + Z_s/T_s * phi * F_s = g_s + Z_s/T_s * <phi>_theta * F_s
! 
! To convert <g> to <f> we use,
!     <f> = <g> + (Z_s/T_s)*<phi>_theta*F_s - (Z_s/T_s)*phi*F_s
! 
! To obtain <f> we add [(Z_s/T_s)*<phi>_theta*F_s - (Z_s/T_s)*phi*F_s] to <g>
!    <field> = (Z_s/T_s)*phi*F_s
!    <gyro_averaged_field> = (Z_s/T_s)*<phi>_theta*F_s
!    f(kx,ky,z,mu,vpa,s) = g(kx,ky,z,mu,vpa,s) + <gyro_averaged_field> - <field>
! 
! To convert <g> to <f> simply supply <f> as an argument instead of <g>,
! and set <facphi> to -1 which will result in:
!    g(kx,ky,z,mu,vpa,s) = f(kx,ky,z,mu,vpa,s) - <gyro_averaged_field> + <field>
! 
!###############################################################################


   !****************************************************************************
   !                   Convert <g> to <f> using (mu,vpa,kxkyzs)                 
   !****************************************************************************
   subroutine g_to_f_kxkyz(g, phi, facphi)

      ! Parallelisation
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: iky_idx, ikx_idx
      use parallelisation_layouts, only: iz_idx, it_idx, is_idx
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average
      
      ! Flags
      use parameters_numerical, only: maxwellian_normalization
      
      ! Grids
      use grids_species, only: spec
      use grids_z, only: nzgrid
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac

      implicit none

      ! Arguments
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      real, intent(in) :: facphi

      ! Local variables
      integer :: ikxkyz, iz, it, iky, ikx, is, ia
      complex, dimension(:, :), allocatable :: field, gyro_averaged_field
      
      !-------------------------------------------------------------------------

      ! Allocate arrays
      allocate (field(nvpa, nmu))
      allocate (gyro_averaged_field(nvpa, nmu))

      ! Assume we only have one field line
      ia = 1

      ! Iterate over the (kx,ky,z,mu,vpa,s) grid
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iky = iky_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         
         ! Calculate <gyro_averaged_field> = Z_s/T_s * <phi>_theta * F_s
         ! First calculate [(Z_s/T_s)*<phi>], add F_s and J_0
         field = facphi * phi(iky, ikx, iz, it) * spec(is)%zt
         if (.not. maxwellian_normalization) then
            field = field * spread(maxwell_vpa(:, is), 2, nmu) * &
                 spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
         end if
         call gyro_average(field, ikxkyz, gyro_averaged_field)
         
         ! Calculate <f> = <g> + (Z_s/T_s)*<phi>_theta*F_s - (Z_s/T_s)*phi*F_s
         g(:, :, ikxkyz) = g(:, :, ikxkyz) + gyro_averaged_field - field
         
      end do

      ! Deallocate local arrays
      deallocate (field, gyro_averaged_field)

   end subroutine g_to_f_kxkyz
   
   !****************************************************************************
   !                 Convert <g> to <f> using (kx,ky,z,ivpamus)                 
   !****************************************************************************
   subroutine g_to_f_vmu(g, phi, facphi, phi_corr)

      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      
      ! Calculations
      use calculations_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
      use calculations_gyro_averages, only: gyro_average
      use arrays_gyro_averages, only: aj0x, aj1x, j0_ffs
      use calculations_kxky, only: multiply_by_rho
      
      ! Flags
      use parameters_physics, only: radial_variation
      use parameters_physics, only: full_flux_surface
      use parameters_numerical, only: maxwellian_normalization
      
      ! Geometry
      use geometry, only: bmag, dBdrho
      use arrays, only: kperp2, dkperp2dr
      
      ! Grids
      use grids_species, only: spec
      use grids_kxky, only: naky, nakx
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac, vperp2, mu, vpa

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :), optional, intent(in) :: phi_corr
      real, intent(in) :: facphi

      ! Local variables
      integer :: ivmu, iz, it, is, imu, iv, ia
      complex, dimension(:, :), allocatable :: field, gyro_averaged_field, g0k
      
      !-------------------------------------------------------------------------

      ! Allocate local arrays
      allocate (field(naky, nakx))
      allocate (gyro_averaged_field(naky, nakx))
      if (radial_variation) allocate (g0k(naky, nakx))

      ! Assume we only have one field line
      ia = 1
      
      ! Iterate over the (kx,ky,z,mu,vpa,s) points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
            
               ! Full-flux-surface
               if (full_flux_surface) then
               
                  !!FLAG!! Need to gyro_averaged_field ia = 1 for ffs
                  field = spec(is)%zt * facphi * phi(:, :, iz, it) * &
                          maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                  call gyro_average(field, gyro_averaged_field, j0_ffs(:, :, iz, ivmu))

               ! Flux-tube
               else
               
                  ! Calculate <gyro_averaged_field> = Z_s/T_s * <phi>_theta * F_s
                  ! First calculate [(Z_s/T_s)*<phi>] and add F_s
                  field = spec(is)%zt * facphi * phi(:, :, iz, it)
                  if (.not. maxwellian_normalization) then
                     field = field * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                  end if
                          
                  ! Add radial variation corrections
                  if (radial_variation .and. present(phi_corr)) then
                     g0k = field * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                        - spec(is)%fprim - 2.0 * dBdrho(iz) * mu(imu) &
                        - 0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                        * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                        * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)))
                     call multiply_by_rho(g0k)
                     field = field + g0k + phi_corr(:, :, iz, it) * spec(is)%zt * facphi &
                       * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                  end if
                  
                  ! Finally add <...>_theta, which is equivalent to adding J_0
                  call gyro_average(field, iz, ivmu, gyro_averaged_field)
                  
               end if
               
               ! Calculate <f> = <g> + (Z_s/T_s)*<phi>_theta*F_s - (Z_s/T_s)*phi*F_s
               g(:, :, iz, it, ivmu) = g(:, :, iz, it, ivmu) + gyro_averaged_field - field
               
            end do
         end do
      end do

      ! Deallocate local arrays
      deallocate (field, gyro_averaged_field)
      if (allocated(g0k)) deallocate (g0k)

   end subroutine g_to_f_vmu

end module calculations_tofrom_ghf
