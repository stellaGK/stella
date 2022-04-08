module mirror_radial_implicit

   implicit none

   public :: mirror_radial_initialized
   public :: init_mirror_radial, finish_mirror_radial
   public :: advance_mirror_radially_implicit

   private

   logical :: mirror_radial_initialized = .false.

   integer, dimension(:, :), allocatable :: mirror_sign
   real, dimension(:, :, :, :), allocatable :: mirror
   real, dimension(:, :, :, :), allocatable :: mirror_interp_loc
   integer, dimension(:, :, :, :), allocatable :: mirror_interp_idx_shift

contains

   subroutine init_mirror_radial

      use stella_time, only: code_dt
      use species, only: spec, nspec
      use vpamu_grids, only: nmu, mu
      use zgrid, only: nzgrid
      use kt_grids, only: nakx, rho_d_clamped
      use stella_geometry, only: dbdzed, b_dot_grad_z
      use stella_geometry, only: d2Bdrdth, dgradpardrho
      use physics_flags, only: include_mirror

      implicit none

      integer :: ia, ix, iz, imu

      if (mirror_radial_initialized) return
      mirror_radial_initialized = .true.

      ia = 1

      if (.not. allocated(mirror)) allocate (mirror(nakx, -nzgrid:nzgrid, nmu, nspec)); mirror = 0.
      if (.not. allocated(mirror_sign)) allocate (mirror_sign(nakx, -nzgrid:nzgrid)); mirror_sign = 0

      !> mirror has sign consistent with being on RHS of GKE;
      !> it is the factor multiplying dg/dvpa in the mirror term
      if (include_mirror) then
         do imu = 1, nmu
            do iz = -nzgrid, nzgrid
               do ix = 1, nakx
                  mirror(ix, iz, imu, :) = code_dt * spec%stm_psi0 * mu(imu) &
                                           * (b_dot_grad_z(ia, iz) * dbdzed(ia, iz) + &
                                              rho_d_clamped(ix) * (dgradpardrho(iz) * dbdzed(ia, iz) &
                                                                   + b_dot_grad_z(ia, iz) * d2Bdrdth(iz)))
               end do
            end do
         end do
      else
         mirror = 0.
      end if

      do ix = 1, nakx
         !> mirror_sign set to +/- 1 depending on the sign of the mirror term.
         !> NB: mirror_sign = -1 corresponds to positive advection velocity
         do iz = -nzgrid, nzgrid
            mirror_sign(ix, iz) = int(sign(1.0, mirror(ix, iz, 1, 1)))
         end do
      end do

      call init_mirror_semi_lagrange

   end subroutine init_mirror_radial

   subroutine init_mirror_semi_lagrange

      use zgrid, only: nzgrid
      use vpamu_grids, only: nmu, dvpa
      use species, only: nspec
      use kt_grids, only: nakx

      implicit none

      if (.not. allocated(mirror_interp_idx_shift)) &
         allocate (mirror_interp_idx_shift(nakx, -nzgrid:nzgrid, nmu, nspec))
      if (.not. allocated(mirror_interp_loc)) &
         allocate (mirror_interp_loc(nakx, -nzgrid:nzgrid, nmu, nspec))

      mirror_interp_idx_shift = int(mirror / dvpa)
      mirror_interp_loc = abs(mod(mirror, dvpa)) / dvpa

      ! f at shifted vpa
      ! is f(iv+idx_shift)*(1-mirror_interp_loc)
      ! + f(iv+idx_shift + mirror_sign)*mirror_interp_loc

   end subroutine init_mirror_semi_lagrange

   ! advance mirror implicit solve dg/dt = mu/m * bhat . grad B (dg/dvpa + m*vpa/T * g)
   subroutine advance_mirror_radially_implicit(time_mirror, g)

      use mp, only: proc0
      use job_manage, only: time_message
      use redistribute, only: gather, scatter
      use stella_layouts, only: vmu_lo, kxkyz_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use vpamu_grids, only: nvpa, nmu
      use dist_fn_arrays, only: gvmu
      use dist_redistribute, only: kxkyz2vmu
      use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      real, dimension(:, :), intent (inout) :: time_mirror 

      integer :: ivmu, iz, it
      complex, dimension(:, :, :), allocatable :: g0v
      complex, dimension(:, :), allocatable :: g0k, g0x


      if (proc0) call time_message(.false., time_mirror(:, 1), ' Mirror advance')

      allocate (g0k(naky, nakx))
      allocate (g0x(naky, nakx))

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0k = g(:, :, iz, it, ivmu)
               call transform_kx2x_unpadded(g0k, g0x)
               g(:, :, iz, it, ivmu) = g0x
            enddo
         enddo
      enddo

      if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
      call scatter(kxkyz2vmu, g, gvmu)
      if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')

      allocate (g0v(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

      call vpa_interpolation(gvmu, g0v)

      ! then take the results and remap again so ky,kx,z local.
      if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
      call gather(kxkyz2vmu, g0v, g)
      if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')

      deallocate (g0v)

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0x = g(:, :, iz, it, ivmu)
               call transform_x2kx_unpadded(g0x, g0k)
               g(:, :, iz, it, ivmu) = g0k
            enddo
         enddo
      enddo

      deallocate (g0k, g0x)

      if (proc0) call time_message(.false., time_mirror, ' Mirror advance')

   end subroutine advance_mirror_radially_implicit

   subroutine vpa_interpolation(grid, interp)

      use vpamu_grids, only: nvpa, nmu
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: ikx_idx, iz_idx, is_idx

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: grid
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(out) :: interp

      integer :: ikxkyz, ix, iz, is, iv, imu
      integer :: shift, sgn, llim, ulim
      real :: fac0, fac1, fac2, fac3
      real :: tmp0, tmp1, tmp2, tmp3

      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         ix = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            tmp0 = mirror_interp_loc(ix, iz, imu, is)
            tmp1 = tmp0 - 2.0
            tmp2 = tmp0 - 1.0
            tmp3 = tmp0 + 1.0

            shift = mirror_interp_idx_shift(ix, iz, imu, is)
            sgn = mirror_sign(ix, iz)

            if (shift == 0) then
            ! deal with boundary point near outgoing BC
               ! using 2-point (linear) interpolation
               ! could do 3-point to improve accuracy
               fac1 = 1.0 - tmp0
               fac2 = tmp0
               if (sgn > 0) then
                  llim = 2
                  ulim = nvpa - 2 - shift
               else
                  llim = nvpa - 1
                  ulim = 3 - shift
               end if
               iv = llim - sgn
               interp(iv, imu, ikxkyz) = grid(iv + shift, imu, ikxkyz) * fac1 &
                                         + grid(iv + shift + sgn, imu, ikxkyz) * fac2
            else
               if (sgn > 0) then
                  llim = 1
                  ulim = nvpa - 2 - shift
               else
                  llim = nvpa
                  ulim = 3 - shift
               end if
            end if

            ! if llim > ulim (for sgn > 0) or llim < ulim (for sgn < 0)
            ! then there are no elements to be interpolated
            if (sgn * ulim < sgn * llim) then
               interp(:, imu, ikxkyz) = 0.
            else
               ! coefficient multiplying g_{iv-1}
               fac0 = -tmp0 * tmp1 * tmp2
               ! coefficient multiplying g_{iv}
               fac1 = 3.*tmp1 * tmp2 * tmp3
               ! coefficient multiplying g_{iv+1}
               fac2 = -3.*tmp0 * tmp1 * tmp3
               ! coefficient multiplying g_{iv+2}
               fac3 = tmp0 * tmp2 * tmp3
               do iv = llim, ulim, sgn
                  interp(iv, imu, ikxkyz) = (grid(iv + shift - sgn, imu, ikxkyz) * fac0 &
                                             + grid(iv + shift, imu, ikxkyz) * fac1 &
                                             + grid(iv + shift + sgn, imu, ikxkyz) * fac2 &
                                             + grid(iv + shift + 2 * sgn, imu, ikxkyz) * fac3) / 6.
               end do

               ! at boundary cell, use zero incoming BC for point just beyond boundary
               interp(ulim + sgn, imu, ikxkyz) = (grid(ulim + shift, imu, ikxkyz) * fac0 &
                                                  + grid(ulim + shift + sgn, imu, ikxkyz) * fac1 &
                                                  + grid(ulim + shift + 2 * sgn, imu, ikxkyz) * fac2) / 6.
               interp(ulim + 2 * sgn, imu, ikxkyz) = (grid(ulim + shift + sgn, imu, ikxkyz) * fac0 &
                                                      + grid(ulim + shift + 2 * sgn, imu, ikxkyz) * fac1) / 6.
               ! use zero incoming BC for cells beyond +/- nvgrid
               if (shift > 0) then
                  interp(nvpa - shift + 1:, imu, ikxkyz) = 0.
               else if (shift < 0) then
                  interp(:-shift, imu, ikxkyz) = 0.
               end if
            end if
         end do
      end do

   end subroutine vpa_interpolation

   subroutine finish_mirror_radial

      use run_parameters, only: mirror_implicit

      implicit none

      if (allocated(mirror)) deallocate (mirror)
      if (allocated(mirror_sign)) deallocate (mirror_sign)

      call finish_mirror_semi_lagrange

      mirror_radial_initialized = .false.

   end subroutine finish_mirror_radial

   subroutine finish_mirror_semi_lagrange

      implicit none

      if (allocated(mirror_interp_loc)) deallocate (mirror_interp_loc)
      if (allocated(mirror_interp_idx_shift)) deallocate (mirror_interp_idx_shift)

   end subroutine finish_mirror_semi_lagrange

end module mirror_radial_implicit
