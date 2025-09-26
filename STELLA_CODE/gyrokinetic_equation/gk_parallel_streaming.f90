!###############################################################################
!############# PARALLEL STREAMING TERM OF THE GYROKINETIC EQUATION #############
!###############################################################################
! 
! This module evolves the parallel streaming term:
!     v_{th,s} v_parallel b . ∇z ( ∂g_{k,s} / ∂z + Z_s/T_s ∂(J_0 ϕ_k) / ∂z exp(-v*2_s))
! 
!###############################################################################
module gk_parallel_streaming

   ! Load the debug flags
   use debug_flags, only: debug => parallel_streaming_debug 
   
   implicit none

   ! Make routines available to other modules
   public :: init_parallel_streaming, finish_parallel_streaming
   public :: advance_parallel_streaming_explicit
   public :: add_parallel_streaming_radial_variation
   public :: stream_tridiagonal_solve
   public :: initialised_parallel_streaming
   public :: stream, stream_c, stream_sign, gradpar_c
   public :: time_parallel_streaming
   public :: stream_rad_var1
   public :: stream_rad_var2
   public :: center_zed, get_dzed
   public :: get_zed_derivative_extended_domain
   public :: stream_correction
   public :: stream_correction_sign
   public :: stream_store_full, stream_full_sign
   public :: get_dgdz_centered

   private

   interface center_zed
      module procedure center_zed_segment_real
      module procedure center_zed_segment_complex
      module procedure center_zed_extended
   end interface center_zed

   integer, dimension(:), allocatable :: stream_sign
   real, dimension(:, :, :, :), allocatable :: stream
   real, dimension(:, :, :), allocatable :: stream_c
   real, dimension(:, :, :), allocatable :: stream_rad_var1
   real, dimension(:, :, :), allocatable :: stream_rad_var2
   real, dimension(:, :), allocatable :: stream_tri_a1, stream_tri_a2
   real, dimension(:, :), allocatable :: stream_tri_b1, stream_tri_b2
   real, dimension(:, :), allocatable :: stream_tri_c1, stream_tri_c2
   real, dimension(:, :), allocatable :: gradpar_c

   real, dimension(2, 3) :: time_parallel_streaming = 0.
   integer, dimension(:,:), allocatable :: stream_correction_sign, stream_full_sign
   real, dimension(:, :, :, :), allocatable :: stream_correction, stream_store_full

   logical :: initialised_parallel_streaming = .false.

contains

   !****************************************************************************
   !                          Initialise Parallel streaming
   !****************************************************************************
   subroutine init_parallel_streaming

      use calculations_finite_differences, only: fd3pt
      use grids_time, only: code_dt
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      use grids_species, only: spec, nspec, pfac
      use grids_velocity, only: nvpa, nvpa
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use grids_velocity, only: vperp2, vpa, mu
      use grids_kxky, only: nalpha
      use grids_z, only: nzgrid, nztot
      use geometry, only: gradpar, dgradpardrho, dBdrho, gfac, b_dot_grad_z
      use parameters_numerical, only: stream_implicit, driftkinetic_implicit
      use parameters_physics, only: include_parallel_streaming, radial_variation
      use parameters_physics, only: full_flux_surface

      implicit none

      integer :: iv, imu, is, ivmu
      integer :: ia, iz

      real, dimension(:), allocatable :: energy
      real, dimension(:, :, :), allocatable :: stream_store

      !-------------------------------------------------------------------------
      
      if (initialised_parallel_streaming) return
      initialised_parallel_streaming = .true.
      if (debug) write (6, *) 'parallel_streaming::init_parallel_streaming'

      ! ========================================================================
      if (.not. allocated(stream)) allocate (stream(nalpha, -nzgrid:nzgrid, nvpa, nspec)); stream = 0.
      if (.not. allocated(stream_sign)) allocate (stream_sign(nvpa)); stream_sign = 0
      !-------------------------------------------------------------------------
      !                        Full Flux Surface simulation
      !-------------------------------------------------------------------------
      if(full_flux_surface) then 
         if (.not. allocated(stream_store_full)) allocate (stream_store_full(nalpha,-nzgrid:nzgrid, nvpa, nspec)); stream_store_full = 0.
         if (.not. allocated(stream_full_sign)) allocate(stream_full_sign(nalpha, nvpa)); stream_full_sign = 0.0
         if (driftkinetic_implicit) then
            if (.not. allocated(stream_correction)) allocate (stream_correction(nalpha, -nzgrid:nzgrid, nvpa, nspec)); stream_correction = 0.
            if (.not. allocated(stream_store)) allocate (stream_store(-nzgrid:nzgrid, nvpa, nspec)); stream_store = 0.
            if (.not. allocated(stream_correction_sign)) allocate (stream_correction_sign(nalpha, nvpa)); stream_correction_sign = 0.0
         end if
      end if

      ! ========================================================================
      ! Get the streaming coefficient. Note that the sign of 'stream; corresponds
      ! to the sign of the coefficient appearing on RHS of GK equation -> i.e. 
      ! the sign of ' - b . grad z vpar '.
      ! This is the factor multiplying dg/dz on RHS of equation
      if (include_parallel_streaming) then
         do iv = 1, nvpa
            do iz = -nzgrid, nzgrid
               do ia = 1, nalpha
                  stream(ia, iz, iv, :) = -code_dt * b_dot_grad_z(ia, iz) * vpa(iv) * spec%stm_psi0
               end do
               !----------------------------------------------------------------
               !               Full Flux Surface simulation
               !----------------------------------------------------------------
               if (driftkinetic_implicit) then
                  stream_store(iz, iv, :) = stream(1, iz, iv, :)
               end if
               !----------------------------------------------------------------
            end do
         end do
      else
         stream = 0.0
      end if

      ! ========================================================================
      ! Get additional arrays for FFS or Radial Variation 
      !-------------------------------------------------------------------------
      !                        Full Flux Surface simulation
      !-------------------------------------------------------------------------
      if(full_flux_surface) then 
         stream_store_full = stream
         if (driftkinetic_implicit) then
            stream_correction = stream - spread(stream_store, 1, nalpha)
            stream = spread(stream_store, 1, nalpha)
            deallocate (stream_store)
         end if
      end if
      !-------------------------------------------------------------------------
      !                           Radial Variation
      !-------------------------------------------------------------------------
      if (radial_variation) then
         allocate (energy(-nzgrid:nzgrid))

         if (.not. allocated(stream_rad_var1)) then
            allocate (stream_rad_var1(-nzgrid:nzgrid, nvpa, nspec))
         end if
         if (.not. allocated(stream_rad_var2)) then
            allocate (stream_rad_var2(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            stream_rad_var2 = 0.0
         end if
         ia = 1
         stream_rad_var1 = -code_dt * spread(spread(spec%stm_psi0, 1, nztot), 2, nvpa) &
                           * gfac * spread(spread(vpa, 1, nztot) * spread(dgradpardrho, 2, nvpa), 3, nspec)
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            energy = (vpa(iv)**2 + vperp2(ia, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
            stream_rad_var2(ia, :, ivmu) = &
               +code_dt * spec(is)%stm_psi0 * vpa(iv) * gradpar &
               * spec(is)%zt * maxwell_vpa(iv, is) * maxwell_mu(ia, :, imu, is) * maxwell_fac(is) &
               * (pfac * (spec(is)%fprim + spec(is)%tprim * (energy - 2.5)) &
                  + gfac * 2 * mu(imu) * dBdrho)
         end do
         deallocate (energy)
      end if

      ! ========================================================================
      ! Create an array that gives the sign of the streaming coeffient. 
      ! This is needed becuase the upwinding factor is dependent on the direction of
      ! advection. 
      ! Here, stream_sign set to +/- 1 depending on the sign of the parallel streaming term.
      ! NB: stream_sign = -1 corresponds to positive advection velocity
      ! only need to consider ia=1, iz=0 and is=1 because alpha, z and species dependences
      ! do not lead to change in sign of the streaming pre-factor
      do iv = 1, nvpa
         stream_sign(iv) = int(sign(1.0, stream(1, 0, iv, 1)))
         !----------------------------------------------------------------------
         !                      Full Flux Surface simulation
         !----------------------------------------------------------------------
         if(full_flux_surface) then 
            do ia = 1, nalpha
               stream_full_sign(ia,:) = int(sign(1.0, stream_store_full(ia, 0, iv, 1)))
               if (driftkinetic_implicit) then         
                   stream_correction_sign(ia,:) =  int(sign(1.0, stream_correction(ia, 0, iv, 1)))
               end if
            end do 
         end if
         !----------------------------------------------------------------------
      end do

      ! ========================================================================
      ! Create centered arrays, which are needed for implicit advance
      if (stream_implicit .or. driftkinetic_implicit) then
         call init_invert_stream_operator
         if (.not. allocated(stream_c)) allocate (stream_c(-nzgrid:nzgrid, nvpa, nspec))
         stream_c = stream(1, :, :, :)
         do is = 1, nspec
            do iv = 1, nvpa
               call center_zed(iv, stream_c(:, iv, is), -nzgrid)
            end do
         end do
         if (.not. allocated(gradpar_c)) allocate (gradpar_c(-nzgrid:nzgrid, -1:1))
         gradpar_c = spread(gradpar, 2, 3)
         ! get gradpar centred in zed for negative vpa (affects upwinding) 
         call center_zed(1, gradpar_c(:, -stream_sign(1)), -nzgrid)
         ! get gradpar centred in zed for positive vpa (affects upwinding)
         call center_zed(nvpa, gradpar_c(:, -stream_sign(nvpa)), -nzgrid)
         stream = spread(stream_c, 1, nalpha)
      end if
      ! ========================================================================

   end subroutine init_parallel_streaming

   !****************************************************************************
   !      Initialise the tridiagonal matrix needed for the stream inversion 
   !****************************************************************************
   ! Only needed for implicit advance of parallel streaming 
   subroutine init_invert_stream_operator

      use grids_z, only: delzed
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: nsegments, neigen_max
      use parameters_numerical, only: zed_upwind_plus, zed_upwind_minus, time_upwind_plus

      implicit none

      integer :: nz, nseg_max

      !-------------------------------------------------------------------------
      
      ! <iz_up> and <iz_low> are not defined when we dont have any connections  
      if (neigen_max>0) then 
          nz = maxval(iz_up - iz_low)
          nseg_max = maxval(nsegments)
      else 
          nz = 0
          nseg_max = 0
      end if

      if (.not. allocated(stream_tri_a1)) then
         allocate (stream_tri_a1(nz * nseg_max + 1, -1:1)); stream_tri_a1 = 0.
         allocate (stream_tri_a2(nz * nseg_max + 1, -1:1)); stream_tri_a2 = 0.
         allocate (stream_tri_b1(nz * nseg_max + 1, -1:1)); stream_tri_b1 = 1.
         allocate (stream_tri_b2(nz * nseg_max + 1, -1:1)); stream_tri_b2 = 0.
         allocate (stream_tri_c1(nz * nseg_max + 1, -1:1)); stream_tri_c1 = 0.
         allocate (stream_tri_c2(nz * nseg_max + 1, -1:1)); stream_tri_c2 = 0.
      end if

      ! Corresponds to sign of stream term positive on RHS of equation
      ! i.e., negative parallel advection speed
      ! NB: assumes equal spacing in zed
      stream_tri_b1(:, 1) = zed_upwind_plus
      stream_tri_b2(:, 1) = -1.0 / delzed(0)
      stream_tri_c1(:nz * nseg_max, 1) = zed_upwind_minus
      stream_tri_c2(:nz * nseg_max, 1) = 1.0 / delzed(0)
      ! Corresponds to sign of stream term negative on RHS of equation
      ! NB: assumes equal spacing in zed
      stream_tri_b1(:, -1) = zed_upwind_plus
      stream_tri_b2(:, -1) = 1.0 / delzed(0)
      stream_tri_a1(2:, -1) = zed_upwind_minus
      stream_tri_a2(2:, -1) = -1.0 / delzed(0)

      stream_tri_a2 = time_upwind_plus * stream_tri_a2
      stream_tri_b2 = time_upwind_plus * stream_tri_b2
      stream_tri_c2 = time_upwind_plus * stream_tri_c2

   end subroutine init_invert_stream_operator

   !****************************************************************************
   !                      Advance Parallel Streaming Explicitly 
   !****************************************************************************
   subroutine advance_parallel_streaming_explicit(g, phi, bpar, gout)

      use mp, only: proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      
      use parameters_numerical, only: maxwellian_normalization
      use parameters_physics, only: full_flux_surface, include_bpar
	
      use grids_velocity, only: mu
      use grids_species, only: spec
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, naky_all, nakx, ikx_max, ny
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac

      use calculations_kxky, only: swap_kxky
      use calculations_transforms, only: transform_ky2y
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1

	   ! For FFS 
      use arrays_gyro_averages, only: j0_ffs

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      integer :: ivmu, iv, imu, is, ia, iz, it
      complex, dimension(:, :, :, :), allocatable :: g0, dgphi_dz, dgbpar_dz
      complex, dimension(:, :, :, :), allocatable :: g0y, g1y
      complex, dimension(:, :), allocatable :: g0_swap

      ! ------------------------------------------------------------------------
      ! Start the timer for the parallel streaming part of the time advance
      if (proc0) call time_message(.false., time_parallel_streaming(:, 1), ' Stream advance')

      ! ========================================================================
      ! Allocate arrays needed for intermmediate calculations
      allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dgphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dgbpar_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
      ! ------------------------------------------------------------------------
      !                            Full Flux Surface
      ! ------------------------------------------------------------------------
      ! If simulating a full flux surface, will also need version of the above arrays
      ! that is Fourier transformed to y-space
      if (full_flux_surface) then
         allocate (g0_swap(naky_all, ikx_max))
         allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes))
         allocate (g1y(ny, ikx_max, -nzgrid:nzgrid, ntubes))
      end if

      ! ========================================================================
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! Get (iv,imu,is) indices corresponding to ivmu super-index
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         ! =====================================================================
         ! Obtain <phi> 
         if (full_flux_surface) then
            ! ------------------------------------------------------------------
            !                           Full Flux Surface
            ! ------------------------------------------------------------------
            call gyro_average(phi, g0(:, :, :, :), j0_ffs(:, :, :, ivmu))
         else
            ! ------------------------------------------------------------------
            !                             Flux Tube
            ! ------------------------------------------------------------------
            call gyro_average(phi, ivmu, g0(:, :, :, :))
         end if
         
         ! =====================================================================
         ! Get d<phi>/dz, with z the parallel coordinate and store in dgphi_dz .
         ! Mote that this should be a centered difference to avoid numerical
         ! unpleasantness to do with inexact cancellations in later velocity integration.
         ! See Appendix of the stella JCP 2019 for details
         call get_dgdz_centered(g0, ivmu, dgphi_dz)
         ! ------------------------------------------------------------------
         !                           Electromagnetic
         ! ------------------------------------------------------------------
         if (include_bpar) then
            call gyro_average_j1(bpar, ivmu, g0(:, :, :, :))
            call get_dgdz_centered(g0, ivmu, dgbpar_dz)
         else
            dgbpar_dz = 0.
         end if 

         ! =====================================================================
         ! Compute dg/dz in k-space and store in g0
         call get_dgdz(g(:, :, :, :, ivmu), ivmu, g0)

         ! =====================================================================
         ! Add the stream term to advance parallel streaming in time
         ! ---------------------------------------------------------------------
         !                              Full Flux Surface
         ! ---------------------------------------------------------------------
         ! If simulating a full flux surface, need to obtain the contribution 
         ! from parallel streaming in y-space, so FFT d(g/F)/dz from ky to y.
         if (full_flux_surface) then
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  ! Get dg/dz in real space
                  call swap_kxky(g0(:, :, iz, it), g0_swap)
                  call transform_ky2y(g0_swap, g0y(:, :, iz, it))

                  ! Get d<phi>/dz in real space
                  call swap_kxky(dgphi_dz(:, :, iz, it), g0_swap)
                  call transform_ky2y(g0_swap, g1y(:, :, iz, it))
               end do
            end do

            g0y(:, :, :, :) = g0y(:, :, :, :) + g1y(:, :, :, :) * spec(is)%zt * maxwell_fac(is) &
                 * maxwell_vpa(iv, is) * spread(spread(maxwell_mu(:, :, imu, is), 2, ikx_max), 4, ntubes) * maxwell_fac(is)
               
            ! Multiply d(g/F)/dz and d<phi>/dz terms with vpa*(b . grad z) and add to source (RHS of GK equation)
            call add_stream_term_full_ffs(g0y, ivmu, gout(:, :, :, :, ivmu))
         else
            ! ------------------------------------------------------------------
            !                           Flux Tube
            ! ------------------------------------------------------------------
            ! Bpar term is zero unless include_bpar = T, see if statement above.
            ia = 1
            if (maxwellian_normalization) then
               g0(:, :, :, :) = g0(:, :, :, :) + dgphi_dz(:, :, :, :) * spec(is)%zt & 
                                               + 4.*mu(imu)*dgbpar_dz(:, :, :, :) 
            else
               g0(:, :, :, :) = g0(:, :, :, :) + (dgphi_dz(:, :, :, :) * spec(is)%zt &
                                               + 4.*mu(imu)*dgbpar_dz(:, :, :, :)) &
                                               * maxwell_fac(is) * maxwell_vpa(iv, is) & 
               * spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx), 4, ntubes)
            end if

            ! Multiply dg/dz with vpa*(b . grad z) and add to source (RHS of GK equation)
            call add_stream_term(g0, ivmu, gout(:, :, :, :, ivmu))
         end if

      end do

      ! ========================================================================
      ! Deallocate intermediate arrays used in this subroutine
      deallocate (g0, dgphi_dz, dgbpar_dz)
      if (full_flux_surface) deallocate (g0y, g1y, g0_swap)
      ! ========================================================================
      ! Finish timing the subroutine
      if (proc0) call time_message(.false., time_parallel_streaming(:, 1), ' Stream advance')

   end subroutine advance_parallel_streaming_explicit

   !****************************************************************************
   !           Radial Variation Contribution for Parallel Streaming 
   !****************************************************************************
   subroutine add_parallel_streaming_radial_variation(g, gout, rhs)

      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      use job_manage, only: time_message
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, nakx
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use grids_species, only: spec
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      use arrays_fields, only: phi, phi_corr_QN, phi_corr_GA

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      ! The next input/output is for quasineutrality and gyroaveraging corrections
      ! that go directly in the RHS (since they don't require further FFTs)
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: rhs

      integer :: ivmu, iv, imu, is, it, ia, iz

      complex, dimension(:, :, :, :), allocatable :: g0, g1, g2, g3
      complex, dimension(:, :), allocatable :: g0k

      ! ========================================================================

      allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g1(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g2(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g3(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g0k(naky, nakx)); g0k = 0

      ! Parallel streaming stays in ky,kx,z space with ky,kx,z local
      ia = 1
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! Obtain <phi>
         call gyro_average(phi, ivmu, g0)
         ! Get d<phi>/dz, with z the parallel coordinate and store in g1
         call get_dgdz_centered(g0, ivmu, g1)
         ! Get variation in gyroaveraging and store in g2
         call get_dgdz_centered(phi_corr_GA(:, :, :, :, ivmu), ivmu, g2)

         ! Get variation in quasineutrality and store in g3
         call gyro_average(phi_corr_QN, ivmu, g0)
         call get_dgdz_centered(g0, ivmu, g3)

         call get_dgdz(g(:, :, :, :, ivmu), ivmu, g0)

         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid

               ! Add variation in gradpar
               g0k = g0(:, :, iz, it) &
                     + g1(:, :, iz, it) * spec(is)%zt * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)

               g0k = g0k * stream_rad_var1(iz, iv, is)

               ! Add variation in F_s/T_s
               g0k = g0k + g1(:, :, iz, it) * stream_rad_var2(ia, iz, ivmu)

               gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + g0k

               ! Add variation in the gyroaveraging and quasineutrality of phi
               ! These variations already have the linear part calculated, so
               ! add it into the rhs directly
               g0k = spec(is)%zt * stream(1, iz, iv, is) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is) &
                     * (g2(:, :, iz, it) + g3(:, :, iz, it))

               rhs(:, :, iz, it, ivmu) = rhs(:, :, iz, it, ivmu) + g0k

            end do
         end do
      end do
      deallocate (g0, g1, g2, g3, g0k)

   end subroutine add_parallel_streaming_radial_variation

   !****************************************************************************
   !                                Add Streaming term
   !****************************************************************************
   subroutine add_stream_term(g, ivmu, src)

      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, is_idx
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: src
      integer, intent(in) :: ivmu

      integer :: iz, iv, is

      !-------------------------------------------------------------------------

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      do iz = -nzgrid, nzgrid
         src(:, :, iz, :) = src(:, :, iz, :) + stream(1, iz, iv, is) * g(:, :, iz, :)
      end do

   end subroutine add_stream_term

   !****************************************************************************
   !                       Add Streaming term - Full Flux Surface
   !****************************************************************************
   subroutine add_stream_term_full_ffs(g, ivmu, src)

      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, is_idx
      use grids_z, only: nzgrid
      use grids_kxky, only: ny

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: src
      integer, intent(in) :: ivmu

      integer :: iz, iy, iv, is

      !-------------------------------------------------------------------------

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      do iz = -nzgrid, nzgrid
         do iy = 1, ny
            src(iy, :, iz, :) = src(iy, :, iz, :) + stream_store_full(iy, iz, iv, is) * g(iy, :, iz, :)
         end do
      end do

   end subroutine add_stream_term_full_ffs

   !****************************************************************************
   !                        Tridiagonal Solve for Streaming 
   !****************************************************************************
   subroutine stream_tridiagonal_solve(iky, ie, iv, is, g)

      use calculations_finite_differences, only: tridag
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: nsegments
      use grids_extended_zgrid, only: nzed_segment

      implicit none

      integer, intent(in) :: iky, ie, iv, is
      complex, dimension(:), intent(in out) :: g

      integer :: iseg, llim, ulim, n
      integer :: nz, nseg_max, sgn, n_ext
      integer :: ia
      real, dimension(:), allocatable :: a, b, c

      !-------------------------------------------------------------------------

      ia = 1

      ! Avoid double-counting at boundaries between 2pi segments
      nz = nzed_segment
      nseg_max = nsegments(ie, iky)
      sgn = stream_sign(iv)

      n_ext = nseg_max * nz + 1
      allocate (a(n_ext))
      allocate (b(n_ext))
      allocate (c(n_ext))

      iseg = 1
      llim = 1; ulim = nz + 1
      a(llim:ulim) = stream_tri_a1(llim:ulim, sgn) &
                     - stream(ia, iz_low(iseg):iz_up(iseg), iv, is) * stream_tri_a2(llim:ulim, sgn)
      b(llim:ulim) = stream_tri_b1(llim:ulim, sgn) &
                     - stream(ia, iz_low(iseg):iz_up(iseg), iv, is) * stream_tri_b2(llim:ulim, sgn)
      c(llim:ulim) = stream_tri_c1(llim:ulim, sgn) &
                     - stream(ia, iz_low(iseg):iz_up(iseg), iv, is) * stream_tri_c2(llim:ulim, sgn)

      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            llim = ulim + 1
            ulim = llim + nz - 1
            a(llim:ulim) = stream_tri_a1(llim:ulim, sgn) &
                           - stream(ia, iz_low(iseg) + 1:iz_up(iseg), iv, is) * stream_tri_a2(llim:ulim, sgn)
            b(llim:ulim) = stream_tri_b1(llim:ulim, sgn) &
                           - stream(ia, iz_low(iseg) + 1:iz_up(iseg), iv, is) * stream_tri_b2(llim:ulim, sgn)
            c(llim:ulim) = stream_tri_c1(llim:ulim, sgn) &
                           - stream(ia, iz_low(iseg) + 1:iz_up(iseg), iv, is) * stream_tri_c2(llim:ulim, sgn)
         end do
      end if
      n = size(stream_tri_a1, 1)
      a(ulim) = stream_tri_a1(n, sgn) - stream(ia, iz_up(nsegments(ie, iky)), iv, is) * stream_tri_a2(n, sgn)
      b(ulim) = stream_tri_b1(n, sgn) - stream(ia, iz_up(nsegments(ie, iky)), iv, is) * stream_tri_b2(n, sgn)
      c(ulim) = 0. ! This line should not be necessary, as c(ulim) should not be accessed by tridag
      call tridag(1, a(:ulim), b(:ulim), c(:ulim), g)

      deallocate (a, b, c)

   end subroutine stream_tridiagonal_solve

   !****************************************************************************
   !                               Get dg/dz
   !****************************************************************************
   ! Get third order accurate upwinded dg/dz, assuming delta zed is equally spaced
   subroutine get_dgdz(g, ivmu, dgdz)

      use calculations_finite_differences, only: third_order_upwind_zed
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx
      use grids_z, only: nzgrid, delzed, ntubes
      use grids_extended_zgrid, only: neigen, nsegments
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: fill_zed_ghost_zones
      use grids_extended_zgrid, only: periodic
      use grids_kxky, only: naky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdz
      integer, intent(in) :: ivmu

      integer :: iseg, ie, it, iky, iv
      complex, dimension(2) :: gleft, gright

      !-------------------------------------------------------------------------

      iv = iv_idx(vmu_lo, ivmu)
      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! First fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g(:, :, :, :), gleft, gright)
                  ! Now get dg/dz
                  call third_order_upwind_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
                     g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                     delzed(0), stream_sign(iv), gleft, gright, periodic(iky), &
                     dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
      end do

   end subroutine get_dgdz

   !****************************************************************************
   !                             Get centered dg/dz
   !****************************************************************************
   ! Get second order accurate centered dg/dz, assuming delta zed is equally spaced
   subroutine get_dgdz_centered(g, ivmu, dgdz)

      use calculations_finite_differences, only: second_order_centered_zed
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx
      use grids_z, only: nzgrid, delzed, ntubes
      use grids_extended_zgrid, only: neigen, nsegments
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: fill_zed_ghost_zones
      use grids_extended_zgrid, only: periodic
      use grids_kxky, only: naky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdz
      integer, intent(in) :: ivmu

      integer :: iseg, ie, iky, iv, it
      complex, dimension(2) :: gleft, gright

      !-------------------------------------------------------------------------

      iv = iv_idx(vmu_lo, ivmu)
      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! First fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g(:, :, :, :), gleft, gright)
                  ! Now get dg/dz
                  call second_order_centered_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
                     g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                     delzed(0), stream_sign(iv), gleft, gright, periodic(iky), &
                     dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
      end do

   end subroutine get_dgdz_centered

   !****************************************************************************
   !                            Get d./dz at cell centres
   !****************************************************************************
   subroutine get_dzed(iv, g, dgdz)

      use calculations_finite_differences, only: fd_cell_centres_zed
      use grids_kxky, only: naky
      use grids_z, only: nzgrid, delzed, ntubes
      use grids_extended_zgrid, only: neigen, nsegments
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: fill_zed_ghost_zones

      implicit none

      integer, intent(in) :: iv
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdz

      integer :: iky, ie, iseg, it
      complex, dimension(2) :: gleft, gright

      !-------------------------------------------------------------------------

      do it = 1, ntubes
         do iky = 1, naky
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! First fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g, gleft, gright)
                  ! Get finite difference approximation for dg/dz at cell centres
                  ! iv > nvgrid corresponds to positive vpa, iv <= nvgrid to negative vpa
                  call fd_cell_centres_zed(iz_low(iseg), &
                     g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                     delzed(0), stream_sign(iv), gleft(2), gright(1), &
                     dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
      end do

   end subroutine get_dzed

   !****************************************************************************
   !                      Get d./dz on the extended zed-domain
   !****************************************************************************
   subroutine get_zed_derivative_extended_domain(iv, f, f_left, f_right, df_dz)

      use grids_z, only: delzed
      use calculations_finite_differences, only: fd_cell_centres_zed

      implicit none

      integer, intent(in) :: iv
      complex, dimension(:), intent(in) :: f
      complex, intent(in) :: f_left, f_right
      complex, dimension(:), intent(out) :: df_dz

      !-------------------------------------------------------------------------

      call fd_cell_centres_zed(1, f, delzed(0), stream_sign(iv), f_left, f_right, df_dz)

   end subroutine get_zed_derivative_extended_domain

   !****************************************************************************
   !                  Get centered d./dz on the extended zed-domain
   !****************************************************************************
   subroutine center_zed_extended(iv, g)

      use calculations_finite_differences, only: cell_centres_zed
      use grids_kxky, only: naky, nakx
      use grids_z, only: nzgrid, ntubes
      use grids_extended_zgrid, only: neigen, nsegments
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: fill_zed_ghost_zones
      use parameters_numerical, only: zed_upwind

      implicit none

      integer, intent(in) :: iv
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: g

      integer :: iky, ie, iseg, it
      complex, dimension(2) :: gleft, gright
      complex, dimension(:, :, :), allocatable :: gc

      !-------------------------------------------------------------------------

      allocate (gc(nakx, -nzgrid:nzgrid, ntubes))

      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! first fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g, gleft, gright)
                  ! get cell centres values
                  call cell_centres_zed(iz_low(iseg), &
                     g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                     zed_upwind, stream_sign(iv), gleft(2), gright(1), &
                     gc(ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
         g(iky, :, :, :) = gc
      end do

      deallocate (gc)

   end subroutine center_zed_extended

   !****************************************************************************
   !                              Center in zed - Real
   !****************************************************************************
   ! center_zed_segment_real takes as arguments the vpa index (iv)
   ! the z-depenendent real function f, and the starting iz index for the array f (llim),
   ! and overwrites f with the cell-centered version
   ! it is assumed that any function passed to this subroutine is periodic
   !****************************************************************************
   subroutine center_zed_segment_real(iv, f, llim)

      use parameters_numerical, only: zed_upwind_plus, zed_upwind_minus

      integer, intent(in) :: iv, llim
      real, dimension(llim:), intent(in out) :: f

      integer :: ulim

      !-------------------------------------------------------------------------

      ulim = llim + size(f) - 1

      if (stream_sign(iv) > 0) then
         f(:ulim - 1) = zed_upwind_plus * f(:ulim - 1) + zed_upwind_minus * f(llim + 1:)
         f(ulim) = f(llim)
      else
         f(llim + 1:) = zed_upwind_minus * f(:ulim - 1) + zed_upwind_plus * f(llim + 1:)
         f(llim) = f(ulim)
      end if

   end subroutine center_zed_segment_real

   !****************************************************************************
   !                          Center in zed - Complex
   !****************************************************************************
   ! center_zed_segment_complex takes as arguments the vpa index (iv)
   ! the z-depenendent conplex function f, and the starting iz index for the array f (llim),
   ! and overwrites f with the cell-centered version;
   !****************************************************************************
   subroutine center_zed_segment_complex(iv, f, llim, periodic)

      use parameters_numerical, only: zupwnd_p => zed_upwind_plus
      use parameters_numerical, only: zupwnd_m => zed_upwind_minus

      integer, intent(in) :: iv, llim
      complex, dimension(llim:), intent(in out) :: f
      logical, intent(in) :: periodic

      integer :: ulim

      !-------------------------------------------------------------------------

      ulim = llim + size(f) - 1

      ! stream_sign > 0 means negative advection velocity
      if (stream_sign(iv) > 0) then
         f(:ulim - 1) = zupwnd_p * f(:ulim - 1) + zupwnd_m * f(llim + 1:)
         if (periodic) then
            f(ulim) = f(llim)
         else
            f(ulim) = zupwnd_p * f(ulim)
         end if
      else
         f(llim + 1:) = zupwnd_m * f(:ulim - 1) + zupwnd_p * f(llim + 1:)
         if (periodic) then
            f(llim) = f(ulim)
         else
            f(llim) = zupwnd_p * f(llim)
         end if
      end if

   end subroutine center_zed_segment_complex

   !****************************************************************************
   !                          Finish Parallel Streaming 
   !****************************************************************************
   subroutine finish_parallel_streaming

      use parameters_numerical, only: stream_implicit, driftkinetic_implicit
      use parameters_physics, only: full_flux_surface
      
      implicit none

      !-------------------------------------------------------------------------

      if (allocated(stream)) deallocate (stream)
      if (allocated(stream_c)) deallocate (stream_c)
      if (allocated(stream_sign)) deallocate (stream_sign)
      if (allocated(gradpar_c)) deallocate (gradpar_c)
      if (allocated(stream_rad_var1)) deallocate (stream_rad_var1)
      if (allocated(stream_rad_var2)) deallocate (stream_rad_var2)

      if (stream_implicit .or. driftkinetic_implicit) call finish_invert_stream_operator
      if(full_flux_surface) then 
         if(allocated(stream_store_full)) deallocate(stream_store_full)
         if(allocated(stream_full_sign)) deallocate(stream_full_sign)
         if (driftkinetic_implicit) then 
            if(allocated(stream_correction)) deallocate(stream_correction)
            if(allocated(stream_correction_sign)) deallocate(stream_correction_sign)
         end if
      end if

      initialised_parallel_streaming = .false.

   end subroutine finish_parallel_streaming

   !****************************************************************************
   !                            Finish Tridiagonal solve
   !****************************************************************************
   ! Only for implicit parallel streaming
   subroutine finish_invert_stream_operator

      implicit none

      if (allocated(stream_tri_a1)) then
         deallocate (stream_tri_a1)
         deallocate (stream_tri_a2)
         deallocate (stream_tri_b1)
         deallocate (stream_tri_b2)
         deallocate (stream_tri_c1)
         deallocate (stream_tri_c2)
      end if

   end subroutine finish_invert_stream_operator

end module gk_parallel_streaming
