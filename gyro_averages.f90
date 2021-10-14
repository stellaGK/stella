module gyro_averages

  use common_types, only: coupled_alpha_type

  public :: aj0x, aj0v, aj1x, aj1v
  public :: init_bessel, finish_bessel
  public :: gyro_average
  public :: gyro_average_j1
  public :: j0_over_B_ffs, j0_ffs, gam0_ffs
  
  private

  interface gyro_average
     module procedure gyro_average_kxky_local
     module procedure gyro_average_kxkyz_local
     module procedure gyro_average_vmu_local
     module procedure gyro_average_vmus_nonlocal
     module procedure gyro_average_ffs_kxky_local
     module procedure gyro_average_ffs_kxkyz_local
  end interface gyro_average

  interface gyro_average_j1
     module procedure gyro_average_j1_kxky_local
     module procedure gyro_average_j1_kxkyz_local
     module procedure gyro_average_j1_vmu_local
  end interface

  real, dimension (:,:,:,:), allocatable :: aj0x, aj1x
  ! (naky, nakx, nalpha, -nzgrid:nzgrid, -vmu-layout-)

  real, dimension (:,:), allocatable :: aj0v, aj1v
  ! (nmu, -kxkyz-layout-)

!  integer, dimension (:,:,:,:), allocatable :: ia_max_aj0a
!  complex, dimension (:,:,:,:,:), allocatable :: aj0a

  type (coupled_alpha_type), dimension (:,:,:,:), allocatable :: j0_ffs, j0_over_B_ffs
  type (coupled_alpha_type), dimension (:,:,:), allocatable :: gam0_ffs
  
!  complex, dimension (:,:,:,:), allocatable :: lu_gam0a
!  integer, dimension (:), allocatable :: lu_gam0a_idx

  logical :: bessinit = .false.

  logical :: debug = .false.
  
contains

  subroutine init_bessel

    use mp, only: sum_allreduce, proc0
    use dist_fn_arrays, only: kperp2
    use physics_flags, only: full_flux_surface
    use species, only: spec, nspec
    use stella_geometry, only: bmag
    use zgrid, only: nzgrid, nztot
    use vpamu_grids, only: vperp2, nmu, nvpa
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
!    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: integrate_species
    use kt_grids, only: naky, nakx, nalpha
    use kt_grids, only: naky_all, ikx_max
    use kt_grids, only: swap_kxky_ordered
    use stella_layouts, only: kxkyz_lo, vmu_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, imu_idx, iv_idx
    use spfunc, only: j0, j1
    use stella_transforms, only: transform_alpha2kalpha
    use file_utils, only: open_output_file, close_output_file
    
    implicit none

    integer :: iz, iky, ikx, imu, is, ia, iv
    integer :: ikxkyz, ivmu
    real :: arg!, dum
    integer :: ia_max_j0_count, ia_max_j0_over_B_count, ia_max_gam0_count
    real :: ia_max_j0_reduction_factor, ia_max_j0_over_B_reduction_factor, ia_max_gam0_reduction_factor

    real, dimension (:), allocatable :: wgts
    real, dimension (:), allocatable :: aj0_alpha, j0_over_B
    complex, dimension (:), allocatable :: aj0_kalpha, j0_over_B_kalpha, gam0_kalpha
    real, dimension (:), allocatable :: gam0_alpha
    real, dimension (:,:,:), allocatable :: kperp2_swap

!    integer :: j0_ffs_unit, j0_over_B_ffs_unit, gam0_ffs_unit

    if (bessinit) return
    bessinit = .true.

    if (debug) write (*,*) 'gyro_averages::init_bessel::allocate_aj0v_aj1v'
    if (.not.allocated(aj0v)) then
       allocate (aj0v(nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       aj0v = 0.
    end if
    if (.not.allocated(aj1v)) then
       allocate (aj1v(nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       aj1v = 0.
    end if

    if (debug) write (*,*) 'gyro_averages::init_bessel::calculate_aj0v_aj1v'
    ia = 1
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          arg = spec(is)%bess_fac*spec(is)%smz_psi0*sqrt(vperp2(ia,iz,imu)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
          aj0v(imu,ikxkyz) = j0(arg)
          ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
          aj1v(imu,ikxkyz) = j1(arg)
       end do
    end do

    if (debug) write (*,*) 'gyro_averages::init_bessel::full_flux_surface'
    if (full_flux_surface) then
!       call open_output_file (j0_ffs_unit, '.j0_ffs')
!       call open_output_file (j0_over_B_ffs_unit, '.j0_over_B_ffs')
!       call open_output_file (gam0_ffs_unit, '.gam0_ffs')
       
       ! wgts are species-dependent factors appearing in Gamma0 factor
       allocate (wgts(nspec))
       wgts = spec%dens*spec%z**2/spec%temp

       if (debug) write (*,*) 'gyro_averages::init_bessel::full_flux_surface::allocate_arrays'
       allocate (aj0_kalpha(naky))
       ! j0_over_B will contain J_0/B as a function of k_alpha and alpha
       allocate (j0_over_B(nalpha))
       allocate (j0_over_B_kalpha(naky))
       allocate (gam0_kalpha(naky))
       allocate (kperp2_swap(naky_all,ikx_max,nalpha))
       if (.not.allocated(j0_ffs)) then
          allocate(j0_ffs(naky_all,ikx_max,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       end if
       if (.not.allocated(j0_over_B_ffs)) then
          allocate(j0_over_B_ffs(naky_all,ikx_max,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       end if
       if (.not.allocated(gam0_ffs)) then
          allocate(gam0_ffs(naky_all,ikx_max,-nzgrid:nzgrid))
       end if

       ia_max_j0_count = 0 ; ia_max_j0_over_B_count = 0 ; ia_max_gam0_count = 0
       do iz = -nzgrid, nzgrid
          write (*,*) 'calculating Fourier coefficients needed for gyro-averaging with alpha variation; zed index: ', iz
          ! for each value of alpha, take kperp^2 calculated on domain kx = [-kx_max, kx_max] and ky = [0, ky_max]
          ! and use symmetry to obtain kperp^2 on domain kx = [0, kx_max] and ky = [-ky_max, ky_max]
          ! this makes later convolutions involving sums over all ky more straightforward
          if (debug) write (*,*) 'gyro_averages::init_bessel::full_flux_surface::swap_kxky'
          do ia = 1, nalpha
             call swap_kxky_ordered (kperp2(:,:,ia,iz), kperp2_swap(:,:,ia))
          end do
          ! aj0_alpha will contain J_0 as a function of k_alpha and alpha
          allocate (aj0_alpha(nalpha))
          if (debug) write (*,*) 'gyro_averages::init_bessel::full_flux_surface::j0_loop'
          do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
             is = is_idx(vmu_lo,ivmu)
             imu = imu_idx(vmu_lo,ivmu)
             do ikx = 1, ikx_max
                do iky = 1, naky_all
                   do ia = 1, nalpha
                      ! calculate the argument of the Bessel function, which depends on both alpha and k_alpha
                      arg = spec(is)%bess_fac*spec(is)%smz_psi0*sqrt(vperp2(ia,iz,imu)*kperp2_swap(iky,ikx,ia))/bmag(ia,iz)
                      ! compute the value of the Bessel function J0 corresponding to argument arg
                      aj0_alpha(ia) = j0(arg)
                      ! compute J_0/B, needed when integrating g over v-space in Maxwell's equations, due to 1/B in v-space Jacobian
                      j0_over_B(ia) = aj0_alpha(ia)/bmag(ia,iz)
                   end do
                   ! fourier transform aj0_alpha and j0_over_B.
                   ! note that fourier coefficients aj0_kalpha and j0_over_B_kalpha have
                   ! been filtered to avoid aliasing
                   !if (debug) write (*,*) 'gyro_averages::init_bessel::full_flux_surface::transform_alpha2kalpha'
                   call transform_alpha2kalpha (aj0_alpha, aj0_kalpha)
                   call transform_alpha2kalpha (j0_over_B, j0_over_B_kalpha)
                   ! given the Fourier coefficients aj0_kalpha, calculate the minimum number of coefficients needed,
                   ! called j0_ffs%max_idx, to ensure that the relative error in the total spectral energy is below a specified tolerance
                   !if (debug) write (*,*) 'gyro_averages::init_bessel::full_flux_surface::find_max_required_kalpha_index'
!                   ! TMP FOR TESTING
!                   j0_ffs(iky,ikx,iz,ivmu)%max_idx = naky
                   call find_max_required_kalpha_index (aj0_kalpha, j0_ffs(iky,ikx,iz,ivmu)%max_idx, imu, iz, is)
                   ! given the Fourier coefficients j0_over_B_kalpha, calculate the minimum number of coefficients needed,
                   ! called j0_over_B_ffs%max_idx, to ensure that the relative error in the total spectral energy is below a specified tolerance
                   call find_max_required_kalpha_index (j0_over_B_kalpha, j0_over_B_ffs(iky,ikx,iz,ivmu)%max_idx, imu, iz, is)
                   ! keep track of the total number of coefficients that must be retained across different phase space points
                   ia_max_j0_count = ia_max_j0_count + j0_ffs(iky,ikx,iz,ivmu)%max_idx
                   ! keep track of the total number of coefficients that must be retained across different phase space points
                   ia_max_j0_over_B_count = ia_max_j0_over_B_count + j0_over_B_ffs(iky,ikx,iz,ivmu)%max_idx
                   ! allocate array to hold the reduced number of Fourier coefficients
                   if (.not.associated(j0_ffs(iky,ikx,iz,ivmu)%fourier)) &
                        allocate (j0_ffs(iky,ikx,iz,ivmu)%fourier(j0_ffs(iky,ikx,iz,ivmu)%max_idx))
                   ! fill the array with the requisite coefficients
                   j0_ffs(iky,ikx,iz,ivmu)%fourier = aj0_kalpha(:j0_ffs(iky,ikx,iz,ivmu)%max_idx)
!                   call test_ffs_bessel_coefs (j0_ffs(iky,ikx,iz,ivmu)%fourier, aj0_alpha, iky, ikx, iz, j0_ffs_unit, ivmu)
                   ! allocate array to hold the reduced number of Fourier coefficients
                   if (.not.associated(j0_over_B_ffs(iky,ikx,iz,ivmu)%fourier)) &
                        allocate (j0_over_B_ffs(iky,ikx,iz,ivmu)%fourier(j0_over_B_ffs(iky,ikx,iz,ivmu)%max_idx))
                   ! fill the array with the requisite coefficients
                   j0_over_B_ffs(iky,ikx,iz,ivmu)%fourier = j0_over_B_kalpha(:j0_over_B_ffs(iky,ikx,iz,ivmu)%max_idx)
!                   call test_ffs_bessel_coefs (j0_over_B_ffs(iky,ikx,iz,ivmu)%fourier, j0_over_B, iky, ikx, iz, j0_over_B_ffs_unit, ivmu)
                end do
             end do
          end do
          ! aj0_alpha will be re-used below with a different size so deallocate and re-allocate below
          deallocate (aj0_alpha)

          ! calculate the reduced set of Fourier coefficients needed to accurately represent \Gamma_0,
          ! the term associated with the polarization density that appears in quasineutrality

          ! re-use aj0_alpha array, and allocate new gam0_alpha array
          allocate (aj0_alpha(vmu_lo%llim_proc:vmu_lo%ulim_alloc))
          allocate (gam0_alpha(nalpha))
          ! in calculating the Fourier coefficients for Gamma_0, change loop orders
          ! so that inner loop is over ivmu super-index;
          ! this is done because we must integrate over v-space and sum over species,
          ! and we want to minimise memory usage where possible (so, e.g., aj0_alpha need
          ! only be a function of ivmu and can be over-written for each (ia,iky,ikx)).
          if (debug) write (*,*) 'gyro_averages::init_bessel::full_flux_surface::gam0_loop'
          do ikx = 1, ikx_max
             do iky = 1, naky_all
                do ia = 1, nalpha
                   !if (debug) write (*,*) 'gyro_averages::init_bessel::full_flux_surface::gam0_loop::aj0_alpha'
                   ! get J0 for all vpar, mu, spec values
                   do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                      is = is_idx(vmu_lo,ivmu)
                      imu = imu_idx(vmu_lo,ivmu)
                      iv = iv_idx(vmu_lo,ivmu)
                      ! calculate the argument of the Bessel function J0
                      arg = spec(is)%bess_fac*spec(is)%smz_psi0*sqrt(vperp2(ia,iz,imu)*kperp2_swap(iky,ikx,ia))/bmag(ia,iz)
                      ! compute J0 corresponding to the given argument arg
                      aj0_alpha(ivmu) = j0(arg)
                      ! form coefficient needed to calculate 1-Gamma_0
                      aj0_alpha(ivmu) = (1.0-aj0_alpha(ivmu)**2) &
                           * maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is)
                   end do

                   !if (debug) write (*,*) 'gyro_averages::init_bessel::full_flux_surface::gam0_loop::integrate_species'
                   ! calculate gamma0(kalpha,alpha,...) = sum_s Zs^2 * ns / Ts int d3v (1-J0^2)*F_{Maxwellian}
                   ! note that v-space Jacobian contains alpha-dependent factor, B(z,alpha),
                   ! but this is not a problem as we have yet to transform from alpha to k_alpha
                   call integrate_species (aj0_alpha, iz, wgts, gam0_alpha(ia), ia)
                end do
                !if (debug) write (*,*) 'gyro_averages::init_bessel::full_flux_surface::gam0_loop::transform_alpha2kalpha'
                ! fourier transform Gamma_0(alpha) from alpha to k_alpha space
                call transform_alpha2kalpha (gam0_alpha, gam0_kalpha)
                ! find the minimum number of Fourier coefficients contained in gam0_kalpha
                ! that must be retained to ensure that the relative error in the total spectral energy
                ! is below a specified tolerance (~1 percent)
                call find_max_required_kalpha_index (gam0_kalpha, gam0_ffs(iky,ikx,iz)%max_idx)
                ! keep track of the total number of coefficients that must be retained across different phase space points
                ia_max_gam0_count = ia_max_gam0_count + gam0_ffs(iky,ikx,iz)%max_idx
                ! allocate array to hold the reduced number of Fourier coefficients
                if (.not.associated(gam0_ffs(iky,ikx,iz)%fourier)) &
                     allocate (gam0_ffs(iky,ikx,iz)%fourier(gam0_ffs(iky,ikx,iz)%max_idx))
                ! fill the array with the requisite coefficients
                gam0_ffs(iky,ikx,iz)%fourier = gam0_kalpha(:gam0_ffs(iky,ikx,iz)%max_idx)
!                call test_ffs_bessel_coefs (gam0_ffs(iky,ikx,iz)%fourier, gam0_alpha, iky, ikx, iz, gam0_ffs_unit)
             end do
          end do
          ! no longer neeed aj0_alpha and gam0_alpha, so deallocate to free up memory
          deallocate (aj0_alpha, gam0_alpha)
       end do
       deallocate (j0_over_B, j0_over_B_kalpha)

!       lu_gam0a = gam0a
!       call lu_decomposition (lu_gam0a, lu_gam0a_idx, dum)

       ! calculate the reduction factor of Fourier modes
       ! used to represent J0
       call sum_allreduce (ia_max_j0_count)
       !       ia_max_j0_reduction_factor = real(ia_max_j0_count)/real(naky*nakx*nztot*nmu*nvpa*nspec*naky)
       ia_max_j0_reduction_factor = real(ia_max_j0_count)/real(naky*ikx_max*nztot*nmu*nvpa*nspec*naky_all)
       call sum_allreduce (ia_max_j0_over_B_count)
       ia_max_j0_over_B_reduction_factor = real(ia_max_j0_over_B_count)/real(naky*ikx_max*nztot*nmu*nvpa*nspec*naky_all)
       call sum_allreduce (ia_max_gam0_count)
       ia_max_gam0_reduction_factor = real(ia_max_gam0_count)/real(naky*ikx_max*nztot*naky_all)

       if (proc0) then
          write (*,*) 'average number of k-alphas needed to represent J0(kperp(alpha))=', ia_max_j0_reduction_factor*naky, 'out of ', naky
          write (*,*) 'average number of k-alphas needed to represent J0(kperp(alpha))/B(alpha)=', ia_max_j0_over_B_reduction_factor*naky, 'out of ', naky
          write (*,*) 'average number of k-alphas needed to represent Gamma0(kperp(alpha))=', ia_max_gam0_reduction_factor*naky, 'out of ', naky
          write (*,*)
       end if

       deallocate (wgts)
       deallocate (aj0_kalpha, gam0_kalpha)
       deallocate (kperp2_swap)

       !       call close_output_file (j0_ffs_unit)
       !       call close_output_file (gam0_ffs_unit)
!       call close_output_file (j0_over_B_ffs_unit)
    else
       if (.not.allocated(aj0x)) then
          allocate (aj0x(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
          aj0x = 0.
       end if

       if (.not.allocated(aj1x)) then
          allocate (aj1x(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
          aj1x = 0.
       end if

       ia = 1
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          do iz = -nzgrid, nzgrid
             do ikx = 1, nakx
                do iky = 1, naky
                   arg = spec(is)%bess_fac*spec(is)%smz_psi0*sqrt(vperp2(ia,iz,imu)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
                   aj0x(iky,ikx,iz,ivmu) = j0(arg)
                   ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
                   aj1x(iky,ikx,iz,ivmu) = j1(arg)
                end do
             end do
          end do
       end do
    end if
    if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_average'
!    call test_gyro_average

  contains

    ! inverse fourier transform coefs%fourier for several phase space points and compare with
    ! unfiltered version in alpha-space
    subroutine test_ffs_bessel_coefs (coefs, f_alpha, iky, ikx, iz, unit, ivmu)

      use stella_layouts, only: vmu_lo, iv_idx, is_idx, imu_idx
      
      implicit none

      complex, dimension (:), intent (in) :: coefs
      real, dimension (:), intent (in) :: f_alpha
      integer, intent (in) :: iky, ikx, iz
      integer, intent (in) :: unit
      integer, intent (in), optional :: ivmu
      
      integer :: iv, imu, is

      if (present(ivmu)) then
         ! coefficients should all be independent of vpa, so only do comparison for one vpa point
         iv = iv_idx(vmu_lo,ivmu)
         if (iv == 1) then
            ! only sample subset of mu locations
            imu = imu_idx(vmu_lo,ivmu)
            if (mod(imu-1,nmu/2-1)==0) then
               is = is_idx(vmu_lo,ivmu)
               call test_ffs_bessel_coefs_subset (coefs, f_alpha, iky, ikx, iz, unit, iv, imu, is)
            end if
         end if
      else
         call test_ffs_bessel_coefs_subset (coefs, f_alpha, iky, ikx, iz, unit)
      end if
         
    end subroutine test_ffs_bessel_coefs

    subroutine test_ffs_bessel_coefs_subset (coefs, f_alpha, iky, ikx, iz, unit, iv, imu, is)

      use constants, only: pi
      use zgrid, only: nzgrid, zed
      use kt_grids, only: naky, nalpha, aky_all_ordered
      use vpamu_grids, only: mu
      use stella_transforms, only: transform_kalpha2alpha
      use stella_geometry, only: alpha
      
      implicit none

      complex, dimension (:), intent (in) :: coefs
      real, dimension (:), intent (in) :: f_alpha
      integer, intent (in) :: iky, ikx, iz
      integer, intent (in) :: unit
      integer, intent (in), optional :: iv, imu, is
      
      complex, dimension (:), allocatable :: coefs_padded
      real, dimension (:), allocatable :: f_alpha_approx

      integer :: ia
      integer :: max_idx
      real :: relative_error
      real, parameter :: minval = 1.0e-3
      
      ! only sample a subset of z locations
      if (mod(iz,nzgrid/2)==0) then
         ! consider only kx = 0
         if (ikx == 1) then
            allocate (coefs_padded(naky))
            allocate (f_alpha_approx(nalpha))
            ! initialize the padded coefficient array to zero
            coefs_padded = 0.0
            ! fill in non-zero entries with truncated Fourier coefficients
            max_idx = size(coefs)
            coefs_padded(:max_idx) = coefs
            ! inverse Fourier transform to get alpha-dependent function
            call transform_kalpha2alpha (coefs_padded, f_alpha_approx)
            if (present(iv)) then
               do ia = 1, nalpha
                  relative_error = 2.0*abs(f_alpha(ia)-f_alpha_approx(ia))/(abs(f_alpha(ia)) + abs(f_alpha_approx(ia)))
                  write (unit,*) alpha(ia), f_alpha(ia), f_alpha_approx(ia), &
                       relative_error, aky_all_ordered(iky), ikx, iz, zed(iz), iv, imu, mu(imu), is
               end do
               ! user 2*pi periodicity in alpha to fill in final point (for visualization purposes)
               ia = 1
               relative_error = 2.0*abs(f_alpha(ia)-f_alpha_approx(ia))/(abs(f_alpha(ia)) + abs(f_alpha_approx(ia)))
               !if (abs(f_alpha(ia)) > minval) then
               !   relative_error = abs((f_alpha(ia)-f_alpha_approx(ia))/f_alpha(ia))
               !else
               !   relative_error = abs(f_alpha(ia)-f_alpha_approx(ia))
               !end if
               write (unit,*) 2.0*pi, f_alpha(1), f_alpha_approx(1), &
                    relative_error, aky_all_ordered(iky), ikx, iz, zed(iz), iv, imu, mu(imu), is
            else
               do ia = 1, nalpha
                  relative_error = 2.0*abs(f_alpha(ia)-f_alpha_approx(ia))/(abs(f_alpha(ia)) + abs(f_alpha_approx(ia)))
                  !if (abs(f_alpha(ia)) > minval) then
                  !   relative_error = abs((f_alpha(ia)-f_alpha_approx(ia))/f_alpha(ia))
                  !else
                  !   relative_error = abs(f_alpha(ia)-f_alpha_approx(ia))
                  !end if
                  write (unit,*) alpha(ia), f_alpha(ia), f_alpha_approx(ia), &
                       relative_error, aky_all_ordered(iky), ikx, iz, zed(iz)
               end do
               ! user 2*pi periodicity in alpha to fill in final point (for visualization purposes)
               ia = 1
               relative_error = 2.0*abs(f_alpha(ia)-f_alpha_approx(ia))/(abs(f_alpha(ia)) + abs(f_alpha_approx(ia)))
               !if (abs(f_alpha(ia)) > minval) then
               !   relative_error = abs((f_alpha(ia)-f_alpha_approx(ia))/f_alpha(ia))
               !else
               !   relative_error = abs(f_alpha(ia)-f_alpha_approx(ia))
               !end if
               write (unit,*) 2.0*pi, f_alpha(1), f_alpha_approx(1), &
                    relative_error, aky_all_ordered(iky), ikx, iz, zed(iz)
            end if
            write (unit,*)
            deallocate (coefs_padded, f_alpha_approx)
         end if
      end if
      
    end subroutine test_ffs_bessel_coefs_subset

    ! set up field that varies as x^2 = rho^2 * cos(angle)^2,
    ! with rho the distance from the origin, and 'angle' is the angle made with the horizontal
    ! if considering a particle at x=0, then rho is thee gyro-radius and angle is the gyro-angle
    ! the gyro-average should theen be 1/(2pi) * int_0^2pi dangle rho^2 * cos(angle)^2 = rho^2/2
    subroutine test_gyro_average

      use constants, only: pi
      use kt_grids, only: ny, nx, x, x0, y, y0
      use kt_grids, only: nakx, ikx_max, naky, naky_all
      use kt_grids, only: swap_kxky, swap_kxky_back
      use stella_transforms, only: transform_x2kx, transform_y2ky
      use stella_transforms, only: transform_kx2x, transform_ky2y
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use vpamu_grids, only: nmu
      use species, only: nspec, spec
      use stella_geometry, only: alpha, bmag, x_displacement_fac
      use spfunc, only: bessi0
      ! TMP FOR TESTING
      use kt_grids, only: akx
      
      implicit none

      real, dimension (:,:), allocatable :: fld_yx
      complex, dimension (:,:), allocatable :: fld_ykx
      complex, dimension (:,:), allocatable :: fld_kykx_swapped
      complex, dimension (:,:), allocatable :: fld_kykx, gyro_fld
      real, dimension (:,:,:,:), allocatable :: gyro_fld_yx
      real :: gyroradius
      integer :: iy, ix, ivmu, iv, imu, is
      ! TMP FOR TESTING
      integer :: iky, ikx
      
      integer, parameter :: iz = 0
      
      allocate (fld_yx(ny,nx))
      allocate (fld_ykx(ny,ikx_max))
      allocate (fld_kykx_swapped(naky_all,ikx_max))
      allocate (fld_kykx(naky,nakx))
      allocate (gyro_fld(naky,nakx))
      allocate (gyro_fld_yx(ny,nx,nmu,nspec))

      if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_averages::fld_yx'
      ! set up field that varies as x^2 = rho^2 * cos(angle)^2 and is constant in y
!      fld_yx = spread(0.1*(x-pi*x0),1,ny)**2
      fld_yx = spread(cos(50.0*(x/x0-pi)),1,ny)
!      fld_yx = spread(exp(-0.1*(x-pi*x0)**2),1,ny)
      if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_averages::transform_x2kx'
      ! transform from (y,x) to (y,kx), with kx going from 0 to kxmax
      call transform_x2kx (fld_yx, fld_ykx)
      if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_averages::transform_y2ky'
      ! transform from (y,kx) to (ky,kx), with ky going from (0,...,kymax,-kymax,...,-dky)
      ! and kx going from 0 to kxmax
      call transform_y2ky (fld_ykx, fld_kykx_swapped)
      if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_averages::swap_kxky_back'
      ! use reality condition to re-arrange array so that ky goes from 0 to kymax
      ! and kx goes from (0,...,kxmax,-kxmax,...,-dkx)
      call swap_kxky_back (fld_kykx_swapped, fld_kykx)

!      do iky = 1, naky
!         do ikx = 1, nakx
!            write (*,*) 'fld_kykx: ', iky, ikx, akx(ikx), real(fld_kykx(iky,ikx)), aimag(fld_kykx(iky,ikx))
!         end do
!         write (*,*)
!      end do
!      stop
      
      ! gyro-average the field at z=0 for different values of mu
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! get the vpa index
         iv = iv_idx(vmu_lo,ivmu)
         ! as J0 independent of vpa, pick only one vpa to test
         if (iv /= 1) cycle
         ! get the species index
         is = is_idx(vmu_lo,ivmu)
         ! get the mu index
         imu = imu_idx(vmu_lo,ivmu)
!         if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_averages::gyro_average'
         if (full_flux_surface) then
            call gyro_average (fld_kykx, gyro_fld, j0_ffs(:,:,iz,ivmu))
         else
            call gyro_average (fld_kykx, iz, ivmu, gyro_fld)
         end if
         ! use reality to re-arrange array entries so that ky goes from (0,...,kymax,-kymax,...,-dky)
         ! and kx goes from 0 to kxmax
         call swap_kxky (gyro_fld, fld_kykx_swapped)
         ! transform from (ky,kx) to (y,kx)
         call transform_ky2y (fld_kykx_swapped, fld_ykx)
         ! transform from (y,kx) to (y,x)
         call transform_kx2x (fld_ykx, gyro_fld_yx(:,:,imu,is))
      end do
      if (debug) write (*,*) 'gyro_averages::init_bessel::test_gyro_averages::write_to_screen'
      ! NB: this is only set up to work on a single processor at the moment
      ! NB: to extend, must move information about gyro_fld onto proc0
      do is = 1, nspec
         do imu = 1, nmu
            do ix = 1, nx
               do iy = 1, ny
                  ! gyro-radius/reference gyro-radius is v_perp/Omega/rho_ref = (v_perp/vths)*(rho_s/rho_ref)
                  ! = vperp * sqrt(T_s/T_ref*m_s/m_ref)(B_ref/Z*B) = vperp / (spec%zstm*bmag)
                  gyroradius = sqrt(vperp2(iy,iz,imu))/(spec(is)%zstm*bmag(iy,iz))
!                  gyroradius = sqrt(vperp2(1,iz,imu))/(spec(is)%zstm*bmag(1,iz))
                  write (42,*) 'y: ', y(iy), 'x: ', x(ix)-x0*pi, 'gyro_fld: ', gyro_fld_yx(iy,ix,imu,is), 'gyroradius: ', gyroradius, 'spec: ', is, &
                       'alpha: ', alpha(iy), 'x_displacement_fac: ', x_displacement_fac(iy,iz)
               end do
               write (42,*)
            end do
            write (42,*)
         end do
         write (42,*)
      end do
      do iy = 1, 1000
         gyroradius = (iy-1)*15.0/999.
         !         write (43,*) 'gyroradius: ', gyroradius, 'bes: ', bessi0(0.1*0.5*(gyroradius/x_displacement_fac(1,iz))**2)*exp(-0.1*0.5*(gyroradius/x_displacement_fac(1,iz))**2)
         !         write (43,*) 'gyroradius: ', gyroradius, 'analytical: ', 0.5*(0.1*gyroradius/x_displacement_fac(1,iz))**2
         write (43,*) 'gyroradius: ', gyroradius, 'analytical: ', j0(50.*gyroradius/(x0*x_displacement_fac(1,iz)))
      end do
      ! TMP FOR TESTING
      stop

      deallocate (fld_yx, fld_ykx)
      deallocate (fld_kykx_swapped, fld_kykx)
      deallocate (gyro_fld, gyro_fld_yx)
      
    end subroutine test_gyro_average
      
  end subroutine init_bessel

  ! subroutine takes a set of Fourier coefficients (ft)
  ! and returns the minimum number of coeffients that must be retained (idx)
  ! to ensure that the relative error in the total spectral energy is
  ! below a specified tolerance (tol_floor)
  subroutine find_max_required_kalpha_index (ft, idx, imu, iz, is)

    use vpamu_grids, only: maxwell_mu

    implicit none

    complex, dimension (:), intent (in) :: ft
    integer, intent (out) :: idx
    integer, intent (in), optional :: imu, iz, is

    real, parameter :: tol_floor = 1.0e-8
    integer :: i, n
    real :: subtotal, total
    real :: tol
    real, dimension (:), allocatable :: ftmod2

    n = size(ft)

    ! use conservative estimate
    ! when deciding number of modes to retain
    if (present(imu) .and. present(iz).and.present(is)) then
       !       tol = min(0.1,tol_floor/maxval(maxwell_mu(:,iz,imu,is)))
       tol = min(1.0e-6,tol_floor/maxval(maxwell_mu(:,iz,imu,is)))
    else
       tol = tol_floor
    end if

    allocate (ftmod2(n))
    ! get spectral energy associated with each mode
    ftmod2 = sqrt(real(ft*conjg(ft)))
    ! get total spectral energy
    total = sqrt(sum(ftmod2))
    subtotal = 0.

    ! find minimum spectral index for which
    ! desired percentage of spectral energy contained
    ! in modes with indices at or below it
    if (total > 0.) then
       i = 1
       do while (subtotal < total*(1.0-tol))
          idx = i
          subtotal = sqrt(sum(ftmod2(:i)))
          i = i + 1
       end do
    else
       idx = 1
    end if

    deallocate (ftmod2)

  end subroutine find_max_required_kalpha_index

  subroutine finish_bessel

    implicit none

    if (allocated(aj0v)) deallocate (aj0v)
    if (allocated(aj1v)) deallocate (aj1v)
    if (allocated(aj0x)) deallocate (aj0x)
    if (allocated(aj1x)) deallocate (aj1x)
    if (allocated(j0_ffs)) deallocate (j0_ffs)
    if (allocated(gam0_ffs)) deallocate (gam0_ffs)
    if (allocated(j0_over_B_ffs)) deallocate (j0_over_B_ffs)
!    if (allocated(lu_gam0a)) deallocate (lu_gam0a)

    bessinit = .false.

  end subroutine finish_bessel

  subroutine gyro_average_kxky_local (field, iz, ivmu, gyro_field)

    implicit none

    complex, dimension (:,:), intent (in) :: field
    integer, intent (in) :: iz, ivmu
    complex, dimension (:,:), intent (out) :: gyro_field

    gyro_field = aj0x(:,:,iz,ivmu)*field

  end subroutine gyro_average_kxky_local

  subroutine gyro_average_kxkyz_local (field, ivmu, gyro_field)

    use zgrid, only: nzgrid, ntubes

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: field
    integer, intent (in) :: ivmu
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: gyro_field

    gyro_field = spread(aj0x(:,:,:,ivmu),4,ntubes)*field

  end subroutine gyro_average_kxkyz_local

  subroutine gyro_average_ffs_kxky_local (field, gyro_field, coefs)

    use kt_grids, only: naky, nakx
    use kt_grids, only: naky_all, ikx_max
    use kt_grids, only: swap_kxky_ordered, swap_kxky_back_ordered

    implicit none

    complex, dimension (:,:), intent (in) :: field
    complex, dimension (:,:), intent (out) :: gyro_field
    type (coupled_alpha_type), dimension (:,:), intent (in) :: coefs
    
    integer :: iky, ikx, ikyp
    integer :: idx
    complex, dimension (:,:), allocatable :: field_kyall, gyro_field_kyall

    ! need to switch from ky>=0 and all kx
    ! to kx>=0 and all ky (using reality condition)
    allocate (field_kyall(naky_all,ikx_max))
    allocate (gyro_field_kyall(naky_all,ikx_max)) ; gyro_field_kyall = 0.
    call swap_kxky_ordered (field, field_kyall)
    ! NB: J0(kx,ky) = J0(-kx,-ky) and Gamma0(kx,ky) = Gamma0(-kx,-ky)
    do ikx = 1, ikx_max
       do iky = 1, naky_all
          ! account for contributions from less positive ky values (and this ky itself)
          do ikyp = 1, min(naky, iky)
             ! idx is the index corresponding to k_alpha - k_alpha'
             ! runs from iky down to 1
             idx = iky - ikyp + 1
             ! if the Fourier coefficient corresponding to this value of (k_alpha-k_alpha',k_alpha')
             ! is sufficiently small, then it will not have been included in the truncated version
             ! of the coefficients; in this case, it makes no contribution to the gyro-average sum
             if (coefs(idx,ikx)%max_idx >= ikyp) then
                gyro_field_kyall(iky,ikx) = gyro_field_kyall(iky,ikx) &
                     + coefs(idx,ikx)%fourier(ikyp)*field_kyall(idx,ikx)
             end if
          end do
          ! if iky = naky_all, then already at max positive ky, so no contributions
          ! from more positive ky value possible
          if (iky == naky_all) cycle
          ! account for contributions from more positive ky values (but not this ky itself,
          ! as already accounted for above
          do ikyp = 2, min(naky,naky_all-iky+1)
             ! idx is the index corresponding to k_alpha - k_alpha'
             ! runs from iky + 1 up to iky + naky (or until naky_all, if it is reached first)
             idx = iky + ikyp - 1
             ! if the Fourier coefficient corresponding to this value of (k_alpha-k_alpha',k_alpha')
             ! is sufficiently small, then it will not have been included in the truncated version
             ! of the coefficients; in this case, it makes no contribution to the gyro-average sum
             if (coefs(idx,ikx)%max_idx >= ikyp) then
                ! the k_alpha' values considered in this loop are negative, but only have
                ! Fourier coefficients for positive ky values;
                ! must use the reality condition to convert this to the equivalent coefficients for negative ky
                gyro_field_kyall(iky,ikx) = gyro_field_kyall(iky,ikx) &
                     + conjg(coefs(idx,ikx)%fourier(ikyp))*field_kyall(idx,ikx)
             end if
          end do
       end do
    end do
             
    !       ! account for contributions from less positive ky values (and this ky itself)
    !       do ia = 1, min(coefs(iky,ikx)%max_idx,iky)
    !          idx = iky-ia+1
    !          gyro_field_kyall(iky,ikx) = gyro_field_kyall(iky,ikx) &
    !               + coefs(idx,ikx)%fourier(ia)*field_kyall(idx,ikx)
    !       end do
    !       ! account for contributions from more positive ky values
    !       if (coefs(iky,ikx)%max_idx > 1 .and. iky /= naky_all) then
    !          do ia = 2, min(coefs(iky,ikx)%max_idx,naky_all-iky+1)
    !             idx = iky+ia-1
    !             gyro_field_kyall(iky,ikx) = gyro_field_kyall(iky,ikx) &
    !                  + coefs(idx,ikx)%fourier(ia)*field_kyall(idx,ikx)
    !          end do
    !       end if
    !    end do
    ! end do

    call swap_kxky_back_ordered (gyro_field_kyall, gyro_field)
    deallocate (field_kyall, gyro_field_kyall)

  end subroutine gyro_average_ffs_kxky_local

  subroutine gyro_average_ffs_kxkyz_local (field, gyro_field, coefs)
    
    use zgrid, only: nzgrid, ntubes
    
    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: field
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: gyro_field
    type (coupled_alpha_type), dimension (:,:,-nzgrid:), intent (in) :: coefs
    
    integer :: iz, it

    do it = 1, ntubes
       do iz = -nzgrid, nzgrid
          call gyro_average (field(:,:,iz,it), gyro_field(:,:,iz,it), coefs(:,:,iz))
       end do
    end do

  end subroutine gyro_average_ffs_kxkyz_local
  
  subroutine gyro_average_vmu_local (distfn, ikxkyz, gyro_distfn)

    use vpamu_grids, only: nvpa

    implicit none

    complex, dimension (:,:), intent (in) :: distfn
    integer, intent (in) :: ikxkyz
    complex, dimension (:,:), intent (out) :: gyro_distfn

    gyro_distfn = spread(aj0v(:,ikxkyz),1,nvpa)*distfn

  end subroutine gyro_average_vmu_local

  subroutine gyro_average_vmus_nonlocal (field, iky, ikx, iz, gyro_field)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (vmu_lo%llim_proc:), intent (in) :: field
    integer, intent (in) :: iky, ikx, iz
    complex, dimension (vmu_lo%llim_proc:), intent (out) :: gyro_field

    gyro_field = aj0x(iky,ikx,iz,:)*field

  end subroutine gyro_average_vmus_nonlocal

  subroutine gyro_average_j1_kxky_local (field, iz, ivmu, gyro_field)

    implicit none

    complex, dimension (:,:), intent (in) :: field
    integer, intent (in) :: iz, ivmu
    complex, dimension (:,:), intent (out) :: gyro_field

    gyro_field = aj1x(:,:,iz,ivmu)*field

  end subroutine gyro_average_j1_kxky_local

  subroutine gyro_average_j1_kxkyz_local (field, ivmu, gyro_field)

    use zgrid, only: nzgrid, ntubes

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: field
    integer, intent (in) :: ivmu
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: gyro_field

    integer :: iz, it

    do it = 1, ntubes
       do iz = -nzgrid, nzgrid
          call gyro_average_j1 (field(:,:,iz,it), iz, ivmu, gyro_field(:,:,iz,it))
       end do
    end do

  end subroutine gyro_average_j1_kxkyz_local

  subroutine gyro_average_j1_vmu_local (distfn, ikxkyz, gyro_distfn)

    use vpamu_grids, only: nvpa

    implicit none

    complex, dimension (:,:), intent (in) :: distfn
    integer, intent (in) :: ikxkyz
    complex, dimension (:,:), intent (out) :: gyro_distfn

    gyro_distfn = spread(aj1v(:,ikxkyz),1,nvpa)*distfn

  end subroutine gyro_average_j1_vmu_local

end module gyro_averages
