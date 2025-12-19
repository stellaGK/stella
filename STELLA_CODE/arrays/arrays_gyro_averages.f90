!###############################################################################
!                            ARRAYS GYRO AVERAGES
!###############################################################################
! This module computes the gyro-average of the fields.
! 
! Gyro-averaging in real space is an integral over theta:
!    < phi >_theta = 1/(2*pi) * int_0^(2*pi) phi dtheta
! In reciprocal space this is equivalent with multiplying with the Bessel function:
!    Fourier < phi >_theta = J_0(k_perp*rho_s) * phi
! 
! There are gyro-average routines each layout of the grid on the processors:
!       kxkyz_layout    (naky*nakx*nzed*ntubes*nspec - 1) points
!       kxyz_layout     (ny*nakx*nzed*ntubes*nspec - 1) points
!       xyz_layout      (ny*nx*nzed*ntubes*nspec - 1) points
!       vmu_layout      (2*nvgrid*nmu*nspec - 1) points
! 
! The Bessel functions J0(x) and J1(x) are stored only for the points <ivmu> and
! <ikxkyz> that are located on the current processor <iproc>.
! 
! For the (v,mu,s)-grid we have aj0x and aj1x(iky, ikx, iz, ivmu)
! For the (kx,ky,z)-grid we have aj0v and aj1v(imu, ikxkyz)
! 
! The argument of the bessel function is a_k = k_perp*rho_s     see Eq.(19)
! if we use v_th = sqrt(2*T/m) and w = q*B/m we find:
!    a_k = k_perpN/rho_r*rho_s = k_perpN/(v_thr/w_r)*(v_perps/w_s)
!        = k_perpN*v_perpsN*(v_ths/v_thr)*(w_r/w_s)
!        = k_perpN*v_perpsN*sqrt(T_s/T_r)*sqrt(m_r/m_s)*(Z_r/Z_s)*(m_s/m_r)*B_r/B_s
!        = k_perpN*v_perpsN*sqrt(T_N*m_N)/Z_N/bmag(ia,iz)
!        = sqrt(kperp2*vperp2)*spec(is)%smz_psi0/bmag
!
!###############################################################################
module arrays_gyro_averages

   ! Load debug flags
   use common_types, only: coupled_alpha_type
   use debug_flags, only: debug => gyro_averages_debug
   use debug_flags, only: debug_test_gyro_average

   implicit none
   
   ! Make routines available to other modules
   public :: init_arrays_bessel_functions
   public :: finish_arrays_bessel_functions
   public :: find_max_required_kalpha_index

   ! Make arrays available
   public :: aj0x, aj0v, aj1x, aj1v
   public :: j0_B_ffs, j0_ffs 
   public :: j1_ffs
   public :: j0_const, j0_B_const
   public :: j0max_const

   private

   ! Dimension (naky, nakx, nalpha, -nzgrid:nzgrid, imuvpaspecies) on vmu-layout[llim_proc:ulim_alloc]
   real, dimension(:, :, :, :), allocatable :: aj0x, aj1x
   
   ! Dimension (nmu, ikxkyzspecies) on kxkyz-layout[llim_proc:ulim_alloc])
   real, dimension(:, :), allocatable :: aj0v, aj1v

   ! Variables for the full flux surface simulations
   type(coupled_alpha_type), dimension(:, :, :, :), allocatable :: j0_ffs, j0_B_ffs
   type(coupled_alpha_type), dimension(:, :, :, :), allocatable :: j1_ffs
   real, dimension(:, :, :, :), allocatable :: j0_const, j0_B_const, j0max_const

   ! Only initialize the Bessel functions once
   logical :: initialised_bessels = .false.

contains

!###############################################################################
!                              BESSEL FUNCTIONS
!###############################################################################
! Store the Bessel functions J0(x) and J1(x) only for the points <ivmu> and
! <ikxkyz> that are located on the current processor <iproc>.
! 
! For the (v,mu,s)-grid we have aj0x and aj1x(iky, ikx, iz, ivmu)
! For the (kx,ky,z)-grod we have aj0v and aj1v(imu, ikxkyz)
! 
! The argument of the bessel function is a_k = k_perp*rho_s     Eq.(19)
! if we use v_th = sqrt(2*T/m) and w = q*B/m we find:
!    a_k = k_perpN/rho_r*rho_s = k_perpN/(v_thr/w_r)*(v_perps/w_s)
!        = k_perpN*v_perpsN*(v_ths/v_thr)*(w_r/w_s)
!        = k_perpN*v_perpsN*sqrt(T_s/T_r)*sqrt(m_r/m_s)*(Z_r/Z_s)*(m_s/m_r)*B_r/B_s
!        = k_perpN*v_perpsN*sqrt(T_N*m_N)/Z_N/bmag(ia,iz)
!        = sqrt(kperp2*vperp2)*spec(is)%smz_psi0/bmag
!###############################################################################

   subroutine init_arrays_bessel_functions

      implicit none
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_bessels) return
      initialised_bessels = .true.
      
      ! Initialise Bessel functions J_0(x) and J_1(x) on the (nmu) grid
      call init_bessel_versus_mu_aj0v_aj1v
      
      ! Initialise Bessel functions J_0(x) and J_1(x) on the (kx,ky) grid
      call init_bessel_versus_kxky_aj0x_aj1x
      
      ! Flag -- remove this
      !if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::test_gyro_average'
      !if(debug_test_gyro_average) call test_gyro_average
      
   end subroutine init_arrays_bessel_functions
   
   !****************************************************************************
   !                    Initialise Bessel functions versus mu                   
   !****************************************************************************
   subroutine init_bessel_versus_mu_aj0v_aj1v

      ! Import Bessel functions
      use spfunc, only: j0, j1
      
      ! Parallelisation
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
      
      ! Grids and geometry
      use geometry, only: bmag
      use grids_velocity, only: nmu
      use grids_velocity, only: vperp2
      use arrays, only: kperp2
      
      ! Species parameters
      use grids_species, only: spec
      
      implicit none
      
      ! Local variables
      integer :: iz, iky, ikx, imu, is, ia
      integer :: ikxkyz
      real :: arg
      
      !-------------------------------------------------------------------------

      ! Allocate the arrays
      if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::allocate_aj0v_aj1v'
      if (.not. allocated(aj0v)) then
         allocate (aj0v(nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
         aj0v = 0.
      end if
      if (.not. allocated(aj1v)) then
         allocate (aj1v(nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
         aj1v = 0.
      end if
      
      ! Assume we have a single field line
      ia = 1

      ! Calculate the Bessel functions J_0(x) and J_1(x) on the (nmu) grid
      ! Note that j1 returns J_1(x)/x, not J_1(x)
      if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::calculate_aj0v_aj1v'
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
         
            ! The argument of the bessel function is a_k = k_perp*rho_s
            arg = spec(is)%bess_fac * spec(is)%smz_psi0 * sqrt(vperp2(ia, iz, imu) * kperp2(iky, ikx, ia, iz)) / bmag(ia, iz)
            
            ! Calculate the Bessel functions with arg = a_k = k_perp*rho_s
            aj0v(imu, ikxkyz) = j0(arg)
            aj1v(imu, ikxkyz) = j1(arg)
            
         end do
      end do
      
   end subroutine init_bessel_versus_mu_aj0v_aj1v

   !****************************************************************************
   !                 Initialise Bessel functions versus (kx,ky)                 
   !****************************************************************************
   subroutine init_bessel_versus_kxky_aj0x_aj1x

      ! Import Bessel functions
      use spfunc, only: j0, j1
      
      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: is_idx, imu_idx
      
      ! Grids and geometry
      use geometry, only: bmag
      use grids_z, only: nzgrid
      use grids_kxky, only: naky, nakx
      use grids_velocity, only: vperp2
      use arrays, only: kperp2
      
      ! Species parameters
      use grids_species, only: spec
      
      ! Flags
      use parameters_physics, only: full_flux_surface
      
      implicit none
      
      ! Local variables
      integer :: iz, iky, ikx, imu, is, ia
      integer :: ivmu
      real :: arg
      
      !-------------------------------------------------------------------------
 
      ! Calculate the Bessel functions for full-flux-surface simulations
      if (full_flux_surface) then
         call init_bessel_ffs 
         
      ! Calculate the Bessel functions for flux-tube simulations
      else
      
         ! Allocate the arrays
         if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::allocate_aj0x_aj1x'
         if (.not. allocated(aj0x)) then
               allocate (aj0x(naky, nakx, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
               aj0x = 0.
         end if 
         if (.not. allocated(aj1x)) then
               allocate (aj1x(naky, nakx, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
               aj1x = 0.
         end if

         ! Assume we have a single field line
         ia = 1
         
         ! Calculate the Bessel functions J_0(x) and J_1(x) on the (kx,ky) grid
         ! Note that j1 returns J_1(x)/x, not J_1(x)
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  do iky = 1, naky
                   
                     ! The argument of the bessel function is a_k = k_perp*rho_s
                     arg = spec(is)%bess_fac * spec(is)%smz_psi0 * sqrt(vperp2(ia, iz, imu) * kperp2(iky, ikx, ia, iz)) / bmag(ia, iz)
                     
                     ! Calculate the Bessel functions with arg = a_k = k_perp*rho_s
                     aj0x(iky, ikx, iz, ivmu) = j0(arg)
                     aj1x(iky, ikx, iz, ivmu) = j1(arg)
                     
                  end do
               end do
            end do
         end do
         
      end if
      
   end subroutine init_bessel_versus_kxky_aj0x_aj1x
   
   !****************************************************************************
   !              Initialise Bessel function for full-flux-surface              
   !****************************************************************************
   subroutine init_bessel_ffs

      use mp, only: sum_allreduce, proc0
      use spfunc, only: j0
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      use calculations_transforms, only: transform_alpha2kalpha
      use grids_species, only: nspec, spec
      use geometry, only: bmag
      use grids_z, only: nzgrid, nztot
      use grids_velocity, only: nmu, nvpa
      use grids_velocity, only: vperp2, maxwell_vpa, maxwell_mu
      use grids_kxky, only: nalpha, naky, naky_all, ikx_max
      use calculations_kxky, only: swap_kxky_ordered
      use arrays, only: kperp2
      use grids_kxky, only: nakx
      use calculations_kxky, only: swap_kxky_back_ordered
      use spfunc, only: j1

      implicit none

      ! integer :: j0_ffs_unit, j0_B_ffs_unit
      integer :: iky, ikx, ia, iz
      integer :: ivmu, iv, imu, is
      integer :: ia_max_j0_count, ia_max_j0_B_count
      real :: arg, rtmp
      real :: ia_max_j0_reduction_factor, ia_max_j0_B_reduction_factor
      real, dimension(:), allocatable :: wgts
      real, dimension(:), allocatable :: aj0_alpha, j0_B
      real, dimension(:, :, :), allocatable :: kperp2_swap
      complex, dimension(:), allocatable :: aj0_kalpha, j0_B_kalpha

      real, dimension(:), allocatable :: aj1_alpha
      complex, dimension(:), allocatable :: aj1_kalpha
      integer :: ia_max_j1_count
      real :: ia_max_j1_reduction_factor

      complex, dimension(:, :), allocatable :: j0_const_in_kalpha, j0_B_const_in_kalpha
      complex, dimension(:, :), allocatable :: j0_const_c, j0_B_const_c

      real, dimension(:), allocatable :: j0max
      complex, dimension(:, :), allocatable :: j0max_const_in_kalpha, j0max_const_c
      
      !-------------------------------------------------------------------------
      
      
      if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::full_flux_surface'

      ! if (debug_test_gyro_average) call open_output_file (j0_ffs_unit, '.j0_ffs')
      ! if (debug_test_gyro_average) call open_output_file (j0_B_ffs_unit, '.j0_over_B_ffs')

      ! wgts are species-dependent factors appearing in Gamma0 factor
      allocate (wgts(nspec))
      wgts = spec%dens * spec%z**2 / spec%temp

      if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::full_flux_surface::allocate_arrays'
      ! aj0_alpha will contain J_0 as a function of k_alpha and alpha
      allocate (aj0_alpha(nalpha))
      allocate (aj0_kalpha(naky))
      ! j0_B will contain J_0*B*exp(-v^2) as a function of k_alpha and alpha
      allocate (j0_B(nalpha))
      allocate (j0_B_kalpha(naky))
      allocate (kperp2_swap(naky_all, ikx_max, nalpha))
      if (.not. allocated(j0_ffs)) then
         allocate (j0_ffs(naky_all, ikx_max, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      end if
      if (.not. allocated(j0_B_ffs)) then
         allocate (j0_B_ffs(naky_all, ikx_max, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      end if

      allocate (aj1_alpha(nalpha))
      allocate (aj1_kalpha(naky))
      if (.not. allocated(j1_ffs)) then
         allocate (j1_ffs(naky_all, ikx_max, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      end if

      allocate (j0_const_in_kalpha(naky_all, ikx_max)); j0_const_in_kalpha = 0.0
      allocate (j0_B_const_in_kalpha(naky_all, ikx_max)); j0_B_const_in_kalpha = 0.0
      allocate (j0_const_c(naky, nakx)); j0_const_c = 0.0
      allocate (j0_B_const_c(naky, nakx)); j0_B_const_c = 0.0
      allocate (j0_const(naky, nakx, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); j0_const = 0.0
      allocate (j0_B_const(naky, nakx, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); j0_B_const = 0.0

      allocate (j0max(nalpha))
      allocate (j0max_const_in_kalpha(naky_all, ikx_max)); j0max_const_in_kalpha = 0.0
      allocate (j0max_const_c(naky, nakx)); j0max_const_c = 0.0
      allocate (j0max_const(naky, nakx, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); j0max_const = 0.0

      ia_max_j0_count = 0; ia_max_j0_B_count = 0
      do iz = -nzgrid, nzgrid
   !         if (proc0) write (*, *) 'calculating Fourier coefficients needed for gyro-averaging with alpha variation; zed index: ', iz
         ! for each value of alpha, take kperp^2 calculated on domain kx = [-kx_max, kx_max] and ky = [0, ky_max]
         ! and use symmetry to obtain kperp^2 on domain kx = [0, kx_max] and ky = [-ky_max, ky_max]
         ! this makes later convolutions involving sums over all ky more straightforward
         if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::full_flux_surface::swap_kxky'
         do ia = 1, nalpha
               call swap_kxky_ordered(kperp2(:, :, ia, iz), kperp2_swap(:, :, ia))
         end do
         if (debug) write (*, *) 'arrays_gyro_averages:init_bessel::full_flux_surface::j0_loop'
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               is = is_idx(vmu_lo, ivmu)
               iv = iv_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               do ikx = 1, ikx_max
               do iky = 1, naky_all
                  do ia = 1, nalpha
                     ! calculate the argument of the Bessel function, which depends on both alpha and k_alpha
                     arg = spec(is)%bess_fac * spec(is)%smz_psi0 * sqrt(vperp2(ia, iz, imu) * kperp2_swap(iky, ikx, ia)) / bmag(ia, iz)
                     ! compute the value of the Bessel function J0 corresponding to argument arg
                     aj0_alpha(ia) = j0(arg)
                     ! compute J_0*B*exp(-v^2), needed when integrating g over v-space in Maxwell's equations,
                     ! due to B in v-space Jacobian
                     j0_B(ia) = aj0_alpha(ia) * bmag(ia, iz)
                     j0max = aj0_alpha(ia) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is)
                     aj1_alpha(ia) = j1(arg)
                  end do

                  ! fourier transform aj0_alpha and j0_B.
                  ! note that fourier coefficients aj0_kalpha and j0_B_kalpha have
                  ! been filtered to avoid aliasing
                  call transform_alpha2kalpha(aj0_alpha, aj0_kalpha)
                  call transform_alpha2kalpha(j0_B, j0_B_kalpha)
                  call transform_alpha2kalpha(aj1_alpha, aj1_kalpha)

                  j0_const_in_kalpha(iky, ikx) = aj0_kalpha(1)
                  j0_B_const_in_kalpha(iky, ikx) = j0_B_kalpha(1)
                  j0max_const_in_kalpha(iky, ikx) = j0max(1)
                     
                  ! given the Fourier coefficients aj0_kalpha, calculate the minimum number of coefficients needed,
                  ! called j0_ffs%max_idx, to ensure that the relative error in the total spectral energy is below a specified tolerance
                  !if (debug) write (*,*) 'arrays_gyro_averages:init_bessel::full_flux_surface::find_max_required_kalpha_index'
                  !                ! TMP FOR TESTING
                  !                j0_ffs(iky,ikx,iz,ivmu)%max_idx = naky
                  call find_max_required_kalpha_index(aj0_kalpha, j0_ffs(iky, ikx, iz, ivmu)%max_idx, imu, iz, is)
                  ! given the Fourier coefficients j0_B_kalpha, calculate the minimum number of coefficients needed,
                  ! called j0_B_ffs%max_idx, to ensure that the relative error in the total spectral energy is below a specified tolerance
                  call find_max_required_kalpha_index(j0_B_kalpha, j0_B_ffs(iky, ikx, iz, ivmu)%max_idx, imu, iz, is)
                  call find_max_required_kalpha_index(aj1_kalpha, j1_ffs(iky, ikx, iz, ivmu)%max_idx, imu, iz, is)
                  ! keep track of the total number of coefficients that must be retained across different phase space points
                  ia_max_j0_count = ia_max_j0_count + j0_ffs(iky, ikx, iz, ivmu)%max_idx
                  ! keep track of the total number of coefficients that must be retained across different phase space points
                  ia_max_j0_B_count = ia_max_j0_B_count + j0_B_ffs(iky, ikx, iz, ivmu)%max_idx
                  ia_max_j1_count = ia_max_j1_count + j1_ffs(iky, ikx, iz, ivmu)%max_idx

                  ! allocate array to hold the reduced number of Fourier coefficients
                  if (.not. associated(j0_ffs(iky, ikx, iz, ivmu)%fourier)) &
                     allocate (j0_ffs(iky, ikx, iz, ivmu)%fourier(j0_ffs(iky, ikx, iz, ivmu)%max_idx))
                  ! fill the array with the requisite coefficients
                  j0_ffs(iky, ikx, iz, ivmu)%fourier = aj0_kalpha(:j0_ffs(iky, ikx, iz, ivmu)%max_idx)
                  ! if (debug_test_gyro_average) call test_ffs_bessel_coefs (j0_ffs(iky,ikx,iz,ivmu)%fourier, aj0_alpha, iky, ikx, iz, j0_ffs_unit, ivmu)
                  if (.not. associated(j0_B_ffs(iky, ikx, iz, ivmu)%fourier)) &
                     allocate (j0_B_ffs(iky, ikx, iz, ivmu)%fourier(j0_B_ffs(iky, ikx, iz, ivmu)%max_idx))
                  ! fill the array with the requisite coefficients
                  j0_B_ffs(iky, ikx, iz, ivmu)%fourier = j0_B_kalpha(:j0_B_ffs(iky, ikx, iz, ivmu)%max_idx)
                  ! if (debug_test_gyro_average) call test_ffs_bessel_coefs (j0_B_ffs(iky,ikx,iz,ivmu)%fourier, j0_B, iky, ikx, iz, j0_B_ffs_unit, ivmu)
                  if (.not. associated(j1_ffs(iky, ikx, iz, ivmu)%fourier)) &
                     allocate (j1_ffs(iky, ikx, iz, ivmu)%fourier(j1_ffs(iky, ikx, iz, ivmu)%max_idx))
                  j1_ffs(iky, ikx, iz, ivmu)%fourier = aj1_kalpha(:j1_ffs(iky, ikx, iz, ivmu)%max_idx)
               end do
               end do
               call swap_kxky_back_ordered(j0_const_in_kalpha, j0_const_c)
               j0_const(:, :, iz, ivmu) = real(j0_const_c)
               call swap_kxky_back_ordered(j0_B_const_in_kalpha, j0_B_const_c)
               j0_B_const(:, :, iz, ivmu) = real(j0_B_const_c)
               call swap_kxky_back_ordered(j0max_const_in_kalpha, j0max_const_c)
               j0max_const(:, :, iz, ivmu) = real(j0max_const_c)
         end do
      end do

      deallocate (j0_B, j0_B_kalpha)
      deallocate (aj0_alpha)
      deallocate (j0_const_in_kalpha, j0_const_c)
      deallocate (j0_B_const_in_kalpha, j0_B_const_c)
      deallocate (j0max_const_in_kalpha, j0max_const_c, j0max)

      ! calculate the reduction factor of Fourier modes
      ! used to represent J0
      ! avoid overflow by converting integers to reals before multiplying
      rtmp = real(naky) * real(naky_all) * real(ikx_max) * real(nztot) * real(nmu) * real(nvpa) * real(nspec)
      call sum_allreduce(ia_max_j0_count)
      ia_max_j0_reduction_factor = real(ia_max_j0_count) / rtmp
      call sum_allreduce(ia_max_j0_B_count)
      ia_max_j0_B_reduction_factor = real(ia_max_j0_B_count) / rtmp

      call sum_allreduce(ia_max_j1_count)
      ia_max_j1_reduction_factor = real(ia_max_j1_count) / rtmp

      if (proc0) then
         write (*, *) 'average number of k-alphas needed to represent J0(kperp(alpha))=', ia_max_j0_reduction_factor * naky, 'out of ', naky
         write (*, *) 'average number of k-alphas needed to represent J0(kperp(alpha))*B(alpha)*exp(-v^2)=', &
               ia_max_j0_B_reduction_factor * naky, 'out of ', naky
         write (*, *)
      end if

      deallocate (wgts)
      deallocate (aj0_kalpha)
      deallocate (kperp2_swap)
      deallocate (aj1_alpha, aj1_kalpha)

      ! if (debug_test_gyro_average) call close_output_file (j0_ffs_unit)
      ! if (debug_test_gyro_average) call close_output_file (j0_B_ffs_unit)
      
   end subroutine init_bessel_ffs

!###############################################################################
!################################## FUNCTIONS ##################################
!###############################################################################

   ! subroutine takes a set of Fourier coefficients (ft)
   ! and returns the minimum number of coeffients that must be retained (idx)
   ! to ensure that the relative error in the total spectral energy is
   ! below a specified tolerance (tol_floor)
   subroutine find_max_required_kalpha_index(ft, idx, imu, iz, is, tol_in)

      use grids_velocity, only: maxwell_mu

      implicit none

      complex, dimension(:), intent(in) :: ft
      integer, intent(out) :: idx
      integer, intent(in), optional :: imu, iz, is
      real, intent(in), optional :: tol_in

      real, parameter :: tol_floor = 1.0e-8
      integer :: i, n
      real :: subtotal, total
      real :: tol
      real, dimension(:), allocatable :: ftmod2

      n = size(ft)

      ! use conservative estimate
      ! when deciding number of modes to retain
      if (present(tol_in)) then
         tol = tol_in
      elseif (present(imu) .and. present(iz) .and. present(is)) then
         !       tol = min(0.1,tol_floor/maxval(maxwell_mu(:,iz,imu,is)))
         tol = min(1.0e-6, tol_floor / maxval(maxwell_mu(:, iz, imu, is)))
      else
         tol = tol_floor
      end if

      allocate (ftmod2(n))
      ! get spectral energy associated with each mode
      ftmod2 = sqrt(real(ft * conjg(ft)))
      ! get total spectral energy
      total = sqrt(sum(ftmod2))
      subtotal = 0.

      ! find minimum spectral index for which
      ! desired percentage of spectral energy contained
      ! in modes with indices at or below it
      if (total > 0.) then
         i = 1
         do while (subtotal < total * (1.0 - tol))
               idx = i
               subtotal = sqrt(sum(ftmod2(:i)))
               i = i + 1
         end do
      else
         idx = 1
      end if

      deallocate (ftmod2)

   end subroutine find_max_required_kalpha_index

!###############################################################################
!################### FINALISE BESSEL FUNCTIONS FOR FLUX TUBE ###################
!###############################################################################

   subroutine finish_arrays_bessel_functions

      implicit none

      if (allocated(aj0v)) deallocate (aj0v)
      if (allocated(aj1v)) deallocate (aj1v)
      if (allocated(aj0x)) deallocate (aj0x)
      if (allocated(aj1x)) deallocate (aj1x)
      if (allocated(j0_ffs)) deallocate (j0_ffs)
      if (allocated(j0_B_ffs)) deallocate (j0_B_ffs)
      if (allocated(j0_B_const)) deallocate (j0_B_const)
      if (allocated(j0_const)) deallocate (j0_const)
      if (allocated(j0max_const)) deallocate (j0max_const)

      initialised_bessels = .false.

   end subroutine finish_arrays_bessel_functions

end module arrays_gyro_averages
