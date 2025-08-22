module multibox

   use fft_work, only: fft_type

   implicit none

   public :: init_multibox
   public :: finish_multibox
   public :: multibox_communicate
   public :: apply_radial_boundary_conditions
   public :: init_mb_get_phi
   public :: mb_get_phi
   public :: communicate_multibox_parameters
   public :: add_multibox_krook
   public :: bs_fullgrid
   public :: xL, xR
   public :: rhoL, rhoR
   public :: kx0_L, kx0_R
   public :: time_multibox
   public :: phi_buffer0, phi_buffer1

   private

   complex, dimension(:), allocatable :: g_buffer0, g_buffer1, phi_buffer0, phi_buffer1
   complex, dimension(:), allocatable :: fsa_x
   real, dimension(:), allocatable :: copy_mask_left, copy_mask_right
   real, dimension(:), allocatable :: krook_mask_left, krook_mask_right
   real, dimension(:), allocatable :: krook_fac
   real, dimension(:), allocatable :: b_mat
   real, dimension(:), allocatable :: x_mb, rho_mb, rho_mb_clamped

   real, dimension(2, 2) :: time_multibox = 0.

   real :: dx_mb

   complex, dimension(:, :), allocatable :: fft_kxky, fft_xky
   real, dimension(:, :), allocatable :: fft_xy

   ! for the unpadded FFTs
   type(fft_type) :: yf_fft, yb_fft
   type(fft_type) :: xf_fft, xb_fft

   complex, dimension(:), allocatable :: fft_x_k, fft_x_x
   complex, dimension(:), allocatable :: fft_y_k
   real, dimension(:), allocatable :: fft_y_y

   logical :: mb_transforms_initialized = .false.
   logical :: get_phi_initialized = .false.
   logical :: use_multibox
   integer :: temp_ind = 0
   integer :: bs_fullgrid
   ! integer :: mb_debug_step
   integer :: x_fft_size
   integer :: phi_bound, phi_pow
   integer :: ikymin

   real :: xL = 0., xR = 0.
   real :: rhoL = 0., rhoR = 0.
   real :: kx0_L, kx0_R
   !real :: efac_l, efacp_l

contains


   subroutine init_multibox

      use constants, only: pi
      use stella_layouts, only: vmu_lo
      use geometry, only: geo_surf, q_as_x, get_x_to_rho
      use geometry, only: drhodpsi, dxdpsi
      use z_grid, only: nzgrid, ntubes
      use parameters_kxky_grid, only: nakx, naky, nx, x0
      use grids_kxky, only: akx, aky
      use grids_kxky, only: x, x_d
      use grids_kxky, only: rho_clamped, rho_d, rho_d_clamped
      use parameters_kxky_grid, only: centered_in_rho
      use parameters_kxky_grid, only: periodic_variation
      use parameters_multibox, only: boundary_size, krook_size
      use file_utils, only: runtype_option_switch, runtype_multibox
      use job_manage, only: njobs
      use parameters_physics, only: rhostar
      use mp, only: scope, crossdomprocs, subprocs, &
                    send, receive, job

      use parameters_multibox, only: lr_debug_switch, lr_debug_option_default, & 
                                     lr_debug_option_L, lr_debug_option_R

      use parameters_multibox, only: krook_option_switch, &
                                     krook_option_flat, krook_option_linear, &
                                     krook_option_exp, krook_option_exp_rev

      use parameters_multibox, only: krook_efold, krook_exponent
      use parameters_multibox, only: use_dirichlet_bc
      use parameters_multibox, only: include_multibox_krook
      use parameters_multibox, only: mb_zf_option_switch, mb_zf_option_skip_ky0, &
                                     mb_zf_option_zero_ky0, mb_zf_option_zero_fsa
      use parameters_multibox, only: nu_krook_mb

      implicit none

      integer :: g_buff_size
      integer :: phi_buff_size
      integer :: i, pfac

      real, dimension(:), allocatable :: x_clamped, x_d_clamped, x_mb_clamped

      integer :: nakxl, nxl, nakxr, nxr, fac
      integer :: copy_size
      real :: db, x_shift, dqdrho

      if (runtype_option_switch == runtype_multibox) then
         call scope(crossdomprocs)

         if (job == 1) then
               call receive(nakxl, 0)
               call receive(nxl, 0)
               call receive(nakxr, njobs - 1)
               call receive(nxr, njobs - 1)

               ! the following assumes nx in the center domain is some
               ! integer multiple of nx in the left or right domain.
               ! Also assumes dx is the same in every domain, which should
               ! be the case
               fac = nx / nxl
               x_fft_size = nakxl * fac
         else
               call send(nakx, 1)
               call send(nx, 1)
               x_fft_size = nakx
         end if

         call scope(subprocs)
      else
         x_fft_size = nakx
      end if

      if ((runtype_option_switch == runtype_multibox) .or. use_dirichlet_bc) then
         use_multibox = .true.
      else
         use_multibox = .false.
      end if

      if (abs(nu_krook_mb) > epsilon(0.0) .and. use_multibox) then
         include_multibox_krook = .true.
      end if

      copy_size = boundary_size - krook_size

      if (.not. use_multibox) return

      bs_fullgrid = nint((3.0 * boundary_size) / 2.0)

      ikymin = 1
      if (mb_zf_option_switch == mb_zf_option_skip_ky0) ikymin = 2

      pfac = 1
      if (periodic_variation) pfac = 2
      phi_buff_size = pfac * boundary_size * naky * ntubes * (2 * nzgrid + 1)
      g_buff_size = phi_buff_size * (vmu_lo%ulim_alloc - vmu_lo%llim_proc + 1)

      if (.not. allocated(g_buffer0)) allocate (g_buffer0(g_buff_size))
      if (.not. allocated(g_buffer1)) allocate (g_buffer1(g_buff_size))
      if (.not. allocated(phi_buffer0)) allocate (phi_buffer0(phi_buff_size))
      if (.not. allocated(phi_buffer1)) allocate (phi_buffer1(phi_buff_size))
      if (.not. allocated(fsa_x) .and. (mb_zf_option_switch == mb_zf_option_zero_fsa)) then
         allocate (fsa_x(nakx)); fsa_x = 0.0
      end if
      if (.not. allocated(copy_mask_left)) allocate (copy_mask_left(pfac * boundary_size)); copy_mask_left = 1.0
      if (.not. allocated(copy_mask_right)) allocate (copy_mask_right(pfac * boundary_size)); copy_mask_right = 1.0
      if (.not. allocated(krook_mask_left)) allocate (krook_mask_left(pfac * boundary_size)); krook_mask_left = 0.0
      if (.not. allocated(krook_mask_right)) allocate (krook_mask_right(pfac * boundary_size)); krook_mask_right = 0.0

      g_buffer0 = 0.
      g_buffer1 = 0.
      phi_buffer0 = 0.
      phi_buffer1 = 0.

      if (krook_size > 0) then
         select case (krook_option_switch)
         case (krook_option_flat)
            do i = 1, krook_size
               krook_mask_right(i) = 1.0
               copy_mask_right(i) = 0.0
            end do
         case (krook_option_linear)
            db = 1.0 / krook_size
            do i = 1, krook_size
               krook_mask_right(i) = i * db
               copy_mask_right(i) = 0.0
            end do
         case (krook_option_exp)
            db = krook_efold / krook_size
            do i = 1, krook_size
               krook_mask_right(i) = 1.0 - (1.0 - exp(-(krook_size - i) * db)) / (1.0 - exp(-krook_efold))
               copy_mask_right(i) = 0.0
            end do
         case (krook_option_exp_rev)
            db = krook_efold / krook_size
            do i = 1, krook_size
               krook_mask_right(i) = (1.0 - exp(-i * db)) / (1.0 - exp(-krook_efold))
               copy_mask_right(i) = 0.0
            end do
         end select
      end if

      if (periodic_variation) then
         do i = 1, boundary_size
            copy_mask_right(i + boundary_size) = copy_mask_right(boundary_size - i + 1)
            krook_mask_right(i + boundary_size) = krook_mask_right(boundary_size - i + 1)
         end do
         copy_mask_left = copy_mask_right
         krook_mask_left = krook_mask_right
      else
         do i = 1, boundary_size
            copy_mask_left(i) = copy_mask_right(boundary_size - i + 1)
            krook_mask_left(i) = krook_mask_right(boundary_size - i + 1)
         end do
      end if

      if (.not. allocated(krook_fac)) allocate (krook_fac(naky))

      krook_fac = 1.0
      do i = 2, naky
         krook_fac(i) = (aky(i) / aky(2))**krook_exponent
      end do

      call init_mb_transforms

      if (.not. allocated(x_mb)) allocate (x_mb(x_fft_size))
      if (.not. allocated(rho_mb)) allocate (rho_mb(x_fft_size))
      if (.not. allocated(rho_mb_clamped)) allocate (rho_mb_clamped(x_fft_size))

      if (runtype_option_switch /= runtype_multibox) then
         x_mb = x_d
         rho_mb = rho_d
         rho_mb_clamped = rho_d_clamped
         return
      end if

      call scope(crossdomprocs)

      dx_mb = (2 * pi * x0) / x_fft_size

      if (job == 1) then
         x_shift = pi * x0
         if (centered_in_rho) then
            if (q_as_x) then
               dqdrho = geo_surf%shat * geo_surf%qinp / geo_surf%rhoc
               x_shift = pi * x0 * (1.0 &
                                    - 0.5 * rhostar * pi * x0 * geo_surf%d2qdr2 / (pfac * dqdrho**2 * dxdpsi))
            else
               x_shift = pi * x0 * (1.0 &
                                    - 0.5 * rhostar * pi * x0 * geo_surf%d2psidr2 * drhodpsi**2 / (pfac * dxdpsi))
            end if
         end if
         do i = 1, x_fft_size
            if (periodic_variation) then
               if (i <= (x_fft_size / 2)) then
                  x_mb(i) = (i - 1) * dx_mb - 0.5 * x_shift
               else
                  x_mb(i) = x_mb(x_fft_size - i + 1)
               end if
            else
               x_mb(i) = (i - 0.5) * dx_mb - x_shift
            end if
         end do
         call get_x_to_rho(1, x_mb, rho_mb)

         xL = x_mb(boundary_size)
         xR = x_mb(x_fft_size / pfac - boundary_size + 1)
         rhoL = rho_mb(boundary_size)
         rhoR = rho_mb(x_fft_size / pfac - boundary_size + 1)

         allocate (x_clamped(nx))
         allocate (x_d_clamped(nakx))
         allocate (x_mb_clamped(x_fft_size))

         if (lr_debug_switch == lr_debug_option_L) then
            x_clamped = xL
            x_d_clamped = xL
            x_mb_clamped = xL
         else if (lr_debug_switch == lr_debug_option_R) then
            x_clamped = xR
            x_d_clamped = xR
            x_mb_clamped = xR
         else
            x_clamped = x
            x_d_clamped = x_d
            x_mb_clamped = x_mb
            do i = 1, nx
               if (x_clamped(i) < xL) x_clamped(i) = xL
               if (x_clamped(i) > xR) x_clamped(i) = xR
            end do
            do i = 1, x_fft_size
               if (x_mb_clamped(i) < xL) x_mb_clamped(i) = xL
               if (x_mb_clamped(i) > xR) x_mb_clamped(i) = xR
            end do
            do i = 1, nakx
               if (x_d_clamped(i) < xL) x_d_clamped(i) = xL
               if (x_d_clamped(i) > xR) x_d_clamped(i) = xR
            end do
         end if

         call get_x_to_rho(1, x_clamped, rho_clamped)
         call get_x_to_rho(1, x_d_clamped, rho_d_clamped)
         call get_x_to_rho(1, x_mb_clamped, rho_mb_clamped)

         deallocate (x_clamped, x_d_clamped)

      elseif (job == 0) then
         do i = 1, x_fft_size
            x_mb(i) = (i - 1) * dx_mb
         end do

         call receive(xL, 1)
         call send(akx(2), 1)
      elseif (job == njobs - 1) then
         do i = 1, x_fft_size
            x_mb(i) = (i - 1) * dx_mb
         end do

         call receive(xR, 1)
         call send(akx(2), 1)
      end if

      call scope(subprocs)

   end subroutine init_multibox

   subroutine communicate_multibox_parameters
      use job_manage, only: njobs
      use mp, only: scope, crossdomprocs, subprocs, &
                    send, receive, job

      implicit none

      if (job == 1) then
         call scope(crossdomprocs)

         call send(xL, 0)
         call send(xR, njobs - 1)

         call receive(kx0_L, 0)
         call receive(kx0_R, njobs - 1)

         call scope(subprocs)

      end if

   end subroutine communicate_multibox_parameters

   subroutine finish_multibox

      implicit none

      if (allocated(g_buffer0)) deallocate (g_buffer0)
      if (allocated(g_buffer1)) deallocate (g_buffer1)
      if (allocated(phi_buffer0)) deallocate (phi_buffer0)
      if (allocated(phi_buffer1)) deallocate (phi_buffer1)
      if (allocated(fsa_x)) deallocate (fsa_x)

      if (allocated(copy_mask_left)) deallocate (copy_mask_left)
      if (allocated(copy_mask_right)) deallocate (copy_mask_right)
      if (allocated(krook_mask_left)) deallocate (krook_mask_left)
      if (allocated(krook_mask_right)) deallocate (krook_mask_right)
      if (allocated(krook_fac)) deallocate (krook_fac)

      if (allocated(fft_kxky)) deallocate (fft_kxky)
      if (allocated(fft_xky)) deallocate (fft_xky)
      if (allocated(fft_xy)) deallocate (fft_xy)

      if (allocated(x_mb)) deallocate (x_mb)
      if (allocated(rho_mb)) deallocate (rho_mb)
      if (allocated(rho_mb_clamped)) deallocate (rho_mb_clamped)

      call finish_mb_transforms

   end subroutine finish_multibox

   subroutine multibox_communicate(gin)

      use constants, only: zi
      use parameters_kxky_grid, only: nakx, naky, naky_all, nx, ny
      use grids_kxky, only: akx, aky, dx, dy, zonal_mode
      use parameters_kxky_grid, only: periodic_variation
      use file_utils, only: runtype_option_switch, runtype_multibox
      use file_utils, only: get_unused_unit
      use store_arrays_fields, only: phi, phi_corr_QN, shift_state
      use job_manage, only: njobs, time_message
      use parameters_physics, only: radial_variation, prp_shear_enabled, hammett_flow_shear
      use parameters_physics, only: g_exb, g_exbfac
      use stella_layouts, only: vmu_lo
      use geometry, only: dl_over_b
      use z_grid, only: nzgrid
      use mp, only: job, scope, mp_abort, &
                    crossdomprocs, subprocs, allprocs, &
                    ssend, receive, proc0
      use parameters_multibox, only: lr_debug_switch, lr_debug_option_default, & 
                                     lr_debug_option_L, lr_debug_option_R
      use parameters_multibox, only: mb_debug_step
      use parameters_multibox, only: phi_bound, phi_pow
      use parameters_multibox, only: boundary_size
      use parameters_multibox, only: mb_zf_option_switch, mb_zf_option_skip_ky0, &
                                     mb_zf_option_zero_ky0, mb_zf_option_zero_fsa

      implicit none

      integer :: num, ia, ix, iix, iky, iz, it, iv
      integer :: offset
      integer :: ii, jj, temp_unit, pfac
      real :: afacx, afacy
      character(len=512) :: filename

      complex, dimension(:, :), allocatable :: prefac
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(inout) :: gin

      if (runtype_option_switch /= runtype_multibox) return
      if (lr_debug_switch /= lr_debug_option_default) return
      if (njobs /= 3) call mp_abort("Multibox only supports 3 domains at the moment.")

      if (proc0) call time_message(.false., time_multibox(:, 1), ' mb_comm')

      allocate (prefac(naky, x_fft_size)); prefac = 1.0

      pfac = 1
      if (periodic_variation) pfac = 2

      if (prp_shear_enabled .and. hammett_flow_shear) then
         prefac = exp(-zi * g_exb * g_exbfac * spread(x_mb, 1, naky) * spread(aky * shift_state, 2, x_fft_size))
      end if

      if (mb_debug_step > 0 .and. mod(temp_ind, mb_debug_step) == 0 .and. proc0) then
         ! call get_unused_unit(temp_unit)
         temp_unit = 3023 + job
         afacx = real(nx) / real(x_fft_size)
         afacy = real(ny) / real(2 * naky - 1)
         fft_kxky = phi(:, :, 0, 1)
         if (radial_variation) then
            fft_kxky = fft_kxky + phi_corr_QN(:, :, 0, 1)
         end if
         call transform_kx2x(fft_kxky, fft_xky)
         fft_xky = fft_xky * prefac
         call transform_ky2y(fft_xky, fft_xy)
         write (filename, "(A,I1,A,I0.6)") "phiout", job, "_", temp_ind
         open (unit=temp_unit, file=filename, status="replace", &
               action="write", form="unformatted", access="stream")
         write (temp_unit) real(x_fft_size, 4)
         do ii = 1, x_fft_size
            write (temp_unit) real(afacx * dx * (ii - 1), 4)
         end do
         do ii = 1, naky_all
            write (temp_unit) real(afacy * dy * (ii - 1), 4)
            do jj = 1, x_fft_size
               write (temp_unit) real(fft_xy(ii, jj), 4)
            end do
         end do
         close (unit=temp_unit)
      end if

! DSO - change communicator
      call scope(crossdomprocs)

      ia = 1
      if (job == 0 .or. job == (njobs - 1)) then
         if (periodic_variation) then
            offset = -boundary_size
         else
            offset = 0; 
            ! DSO the next line might seem backwards, but this makes it easier to stitch together imaages
            ! FLAG DSO - might do something weird with magnetic shear
            if (job == njobs - 1) offset = nakx - boundary_size
         end if
         !first g
         num = 1
         do iv = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, vmu_lo%ntubes

               !this is where the FSA goes
               if (zonal_mode(1) .and. mb_zf_option_switch == mb_zf_option_zero_fsa) then
                  do ix = 1, nakx
                     fsa_x(ix) = sum(dl_over_b(ia, :) * gin(1, ix, :, it, iv))
                  end do
               end if

               do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
                  fft_kxky = gin(:, :, iz, it, iv)

                  if (zonal_mode(1)) then
                     if (mb_zf_option_switch == mb_zf_option_zero_ky0) then
                        fft_kxky(1, :) = 0.0
                     elseif (mb_zf_option_switch == mb_zf_option_zero_fsa) then
                        fft_kxky(1, :) = fft_kxky(1, :) - fsa_x
                     end if
                  end if

                  call transform_kx2x(fft_kxky, fft_xky)
                  fft_xky = fft_xky * prefac
                  do ix = 1, pfac * boundary_size
                     iix = ix + offset
                     if (iix <= 0) iix = iix + nakx
                     do iky = ikymin, naky
                        !DSO if in the future the grids can have different naky, one will
                        !have to divide by naky here, and multiply on the receiving end
                        g_buffer0(num) = fft_xky(iky, iix)
                        num = num + 1
                     end do
                  end do
               end do
            end do
         end do
         call ssend(g_buffer0, 1, 43 + job)
         !now phi
         num = 1
         do it = 1, vmu_lo%ntubes
            !this is where the FSA goes
            if (zonal_mode(1) .and. mb_zf_option_switch == mb_zf_option_zero_fsa) then
               do ix = 1, nakx
                  fsa_x(ix) = sum(dl_over_b(ia, :) * phi(1, ix, :, it))
               end do
            end if

            do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
               fft_kxky = spread((zi * akx)**phi_pow, 1, naky) * phi(:, :, iz, it)

               if (zonal_mode(1)) then
                  if (mb_zf_option_switch == mb_zf_option_zero_ky0) then
                     fft_kxky(1, :) = 0.0
                  elseif (mb_zf_option_switch == mb_zf_option_zero_fsa) then
                     fft_kxky(1, :) = fft_kxky(1, :) - fsa_x
                  end if
               end if

               call transform_kx2x(fft_kxky, fft_xky)
               fft_xky = fft_xky * prefac
               do iky = ikymin, naky
                  do ix = 1, pfac * boundary_size
                     iix = ix + offset
                     if (iix <= 0) iix = iix + nakx
                     !DSO if in the future the grids can have different naky, one will
                     !have to divide by naky here, and multiply on the receiving end
                     phi_buffer0(num) = fft_xky(iky, iix)
                     num = num + 1
                  end do
               end do
            end do
         end do
! DSO - send data
         call ssend(phi_buffer0, 1, 143 + job)
      else
! DSO - receive the data
! left
         call receive(g_buffer0, 0, 43)
         call receive(phi_buffer0, 0, 143)
! right
         call receive(g_buffer1, njobs - 1, 43 + njobs - 1)
         call receive(phi_buffer1, njobs - 1, 143 + njobs - 1)

         !apply the BCs
         call apply_radial_boundary_conditions(gin)
      end if

! DSO - change communicator
      call scope(subprocs)

      temp_ind = temp_ind + 1

      deallocate (prefac)

      if (proc0) call time_message(.false., time_multibox(:, 1), ' mb_comm')

   end subroutine multibox_communicate

   subroutine apply_radial_boundary_conditions(gin)
 
      use grids_kxky, only: zonal_mode
      use parameters_multibox, only: boundary_size
      use parameters_kxky_grid, only: periodic_variation, naky
      use stella_layouts, only: vmu_lo
      use z_grid, only: nzgrid

      use parameters_multibox, only: smooth_zf

      implicit none

      integer :: num, ix, iky, iz, it, iv, iL, iR
      integer :: offsetL, offsetR, pfac
      complex :: dzm, dzp

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(inout) :: gin

      pfac = 1
      if (periodic_variation) pfac = 2

      offsetL = 0
      offsetR = x_fft_size - boundary_size
      if (periodic_variation) then
         offsetL = -boundary_size
         offsetR = x_fft_size / 2 - boundary_size + 1
      end if

      num = 1
      do iv = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, vmu_lo%ntubes
            do iz = -vmu_lo%nzgrid, vmu_lo%nzgrid
               call transform_kx2x(gin(:, :, iz, it, iv), fft_xky)
               do ix = 1, pfac * boundary_size
                  iL = ix + offsetL
                  iR = ix + offsetR
                  if (iL <= 0) iL = iL + x_fft_size
                  if (iR <= 0) iR = iR + x_fft_size
                  do iky = ikymin, naky
                     fft_xky(iky, iL) = fft_xky(iky, iL) * (1 - copy_mask_left(ix)) &
                                        + g_buffer0(num) * copy_mask_left(ix)

                     fft_xky(iky, iR) = fft_xky(iky, iR) * (1 - copy_mask_right(ix)) &
                                        + g_buffer1(num) * copy_mask_right(ix)
                     num = num + 1
                  end do
               end do
               if (smooth_zf) then
                  dzm = fft_xky(1, boundary_size + 1) - fft_xky(1, boundary_size)
                  dzp = fft_xky(1, x_fft_size - boundary_size + 1) - fft_xky(1, x_fft_size - boundary_size)
                  do ix = 1, pfac * boundary_size
                     iL = ix + offsetL
                     iR = ix + offsetR
                     if (iL <= 0) iL = iL + x_fft_size
                     if (iR <= 0) iR = iR + x_fft_size
                     fft_xky(1, iL) = fft_xky(1, iL) + dzm
                     fft_xky(1, iR) = fft_xky(1, iR) - dzp
                  end do
               end if
               if (zonal_mode(1)) fft_xky(1, :) = real(fft_xky(1, :))
               call transform_x2kx(fft_xky, gin(:, :, iz, it, iv))
            end do
         end do
      end do

   end subroutine apply_radial_boundary_conditions

   subroutine add_multibox_krook(g, rhs)

      use stella_time, only: code_dt
      use stella_layouts, only: vmu_lo
      use parameters_kxky_grid, only: nakx, naky, periodic_variation
      use parameters_multibox, only: boundary_size
      use z_grid, only: nzgrid, ntubes
      use mp, only: job, proc0
      use job_manage, only: time_message
      use parameters_multibox, only: nu_krook_mb
      use parameters_multibox, only: use_dirichlet_BC

      implicit none

      integer :: iky, ix, iL, iR, iz, it, ivmu, num, offsetL, offsetR
      integer :: pfac

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: rhs

      complex, allocatable, dimension(:, :) :: g0x, g0k
      if (job /= 1 .and. .not. use_dirichlet_BC) return

      if (proc0) call time_message(.false., time_multibox(:, 2), ' mb_krook')

      allocate (g0k(naky, nakx))
      allocate (g0x(naky, x_fft_size))

      if (periodic_variation) then
         offsetL = -boundary_size
         offsetR = x_fft_size / 2 - boundary_size + 1
      else
         offsetL = 0
         offsetR = x_fft_size - boundary_size
      end if

      pfac = 1
      if (periodic_variation) pfac = 2

      num = 1
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0x = 0.0
               call transform_kx2x(g(:, :, iz, it, ivmu), fft_xky)
               do ix = 1, pfac * boundary_size
                  iL = ix + offsetL
                  iR = ix + offsetR
                  if (iL <= 0) iL = iL + x_fft_size
                  if (iR <= 0) iR = iR + x_fft_size
                  do iky = ikymin, naky
                     g0x(iky, iL) = (fft_xky(iky, iL) - g_buffer0(num)) * krook_mask_left(ix)
                     g0x(iky, iR) = (fft_xky(iky, iR) - g_buffer1(num)) * krook_mask_right(ix)
                     num = num + 1
                  end do
               end do
               call transform_x2kx(g0x, g0k)
               rhs(:, :, iz, it, ivmu) = rhs(:, :, iz, it, ivmu) &
                                         - code_dt * nu_krook_mb * spread(krook_fac, 2, nakx) * g0k
            end do
         end do
      end do

      deallocate (g0k, g0x)

      if (proc0) call time_message(.false., time_multibox(:, 2), ' mb_krook')

   end subroutine add_multibox_krook

!!
!!>DSO - The following subroutines solve for phi in the _physical_ region of space
!!       It is done here because the radial grid may include an extra point

   subroutine init_mb_get_phi(has_elec, adiabatic_elec, efac, efacp)
      use parameters_kxky_grid, only: nakx, naky
      use parameters_multibox, only: boundary_size
      use z_grid, only: nzgrid
      use parameters_physics, only: radial_variation
      use geometry, only: dl_over_b, d_dl_over_b_drho
      use parameters_multibox, only: ky_solve_radial
      use store_arrays_fields, only: phi_solve, phizf_solve, gamtot, dgamtotdr
      use linear_solve, only: lu_decomposition, lu_inverse

      use parameters_multibox, only: phi_bound
      implicit none

      logical, intent(in) :: has_elec, adiabatic_elec
      real, intent(in) :: efac, efacp
      integer :: ia, iz, iky, ikx, b_solve
      real :: dum
      complex, dimension(:, :), allocatable :: g0k, g1k, g0x, a_inv, a_fsa

      if (.not. radial_variation) return

      !this does not depend on the timestep, so only do once
      if (get_phi_initialized) return

      get_phi_initialized = .true.

      !efac_l  = efac
      !efacp_l = efacp_l

      ia = 1
      b_solve = boundary_size - phi_bound

      allocate (g0k(1, nakx))
      allocate (g1k(1, nakx))
      allocate (g0x(1, x_fft_size))

      if (.not. allocated(phi_solve)) allocate (phi_solve(min(ky_solve_radial, naky), -nzgrid:nzgrid))

      do iky = 1, min(ky_solve_radial, naky)
         do iz = -nzgrid, nzgrid
            if (.not. associated(phi_solve(iky, iz)%zloc)) &
               allocate (phi_solve(iky, iz)%zloc(x_fft_size - 2 * b_solve, x_fft_size - 2 * b_solve))
            if (.not. associated(phi_solve(iky, iz)%idx)) &
               allocate (phi_solve(iky, iz)%idx(x_fft_size - 2 * b_solve))

            phi_solve(iky, iz)%zloc = 0.0
            phi_solve(iky, iz)%idx = 0
            do ikx = 1 + b_solve, x_fft_size - b_solve
               g0x(1, :) = 0.0
               g0x(1, ikx) = 1.0
               call transform_x2kx(g0x, g0k)

               g1k(1, :) = g0k(1, :) * gamtot(iky, :, iz)
               call transform_kx2x(g1k, g0x)

               !row column
               phi_solve(iky, iz)%zloc(:, ikx - b_solve) = g0x(1, (1 + b_solve):(x_fft_size - b_solve))

               g1k(1, :) = g0k(1, :) * dgamtotdr(iky, :, iz)
               call transform_kx2x(g1k, g0x)
               g0x(1, :) = rho_mb_clamped * g0x(1, :)

               phi_solve(iky, iz)%zloc(:, ikx - b_solve) = phi_solve(iky, iz)%zloc(:, ikx - b_solve) &
                                                           + g0x(1, (1 + b_solve):(x_fft_size - b_solve))
            end do

            call lu_decomposition(phi_solve(iky, iz)%zloc, phi_solve(iky, iz)%idx, dum)
!       call zgetrf(nakx,nakx,phi_solve(iky,iz)%zloc,nakx,phi_solve(iky,iz)%idx,info)
         end do
      end do

      if (.not. has_elec) then
         if (.not. allocated(b_mat)) allocate (b_mat(x_fft_size - 2 * b_solve)); b_mat = 0.0
         do ikx = 1 + b_solve, x_fft_size - b_solve
            !row column
            b_mat(ikx - b_solve) = efac + efacp * rho_mb_clamped(ikx)
         end do
      end if

      if (adiabatic_elec) then

         allocate (a_inv(x_fft_size - 2 * b_solve, x_fft_size - 2 * b_solve))
         allocate (a_fsa(x_fft_size - 2 * b_solve, x_fft_size - 2 * b_solve)); a_fsa = 0.0

         if (.not. associated(phizf_solve%zloc)) &
            allocate (phizf_solve%zloc(x_fft_size - 2 * b_solve, x_fft_size - 2 * b_solve)); 
         phizf_solve%zloc = 0.0

         if (.not. associated(phizf_solve%idx)) allocate (phizf_solve%idx(x_fft_size - 2 * b_solve)); 
         !get inverse of A
         do iz = -nzgrid, nzgrid

            call lu_inverse(phi_solve(1, iz)%zloc, phi_solve(1, iz)%idx, a_inv)

            !flux surface average it
            do ikx = 1, x_fft_size - 2 * b_solve
               a_fsa(:, ikx) = a_fsa(:, ikx) + (dl_over_b(ia, iz) + d_dl_over_b_drho(ia, iz) &
                                                * rho_mb_clamped((1 + b_solve):(x_fft_size - b_solve))) * a_inv(:, ikx)
            end do
         end do

         ! calculate I - <A^-1>B
         do ikx = 1, x_fft_size - 2 * b_solve
            phizf_solve%zloc(ikx, :) = -a_fsa(ikx, :) * b_mat
            phizf_solve%zloc(ikx, ikx) = 1.0 + phizf_solve%zloc(ikx, ikx)
         end do

         call lu_decomposition(phizf_solve%zloc, phizf_solve%idx, dum)

         deallocate (a_inv, a_fsa)

      end if

      deallocate (g0k, g1k, g0x)
   end subroutine init_mb_get_phi

   subroutine mb_get_phi(phi, has_elec, adiabatic_elec)
      use constants, only: zi
      use parameters_kxky_grid, only: nakx, naky
      use parameters_multibox, only: boundary_size
      use grids_kxky, only: zonal_mode, akx
      use z_grid, only: nzgrid, ntubes
      use geometry, only: dl_over_b, d_dl_over_b_drho
      use parameters_multibox, only: ky_solve_radial
      use store_arrays_fields, only: gamtot, dgamtotdr, phi_solve, phizf_solve
      use linear_solve, only: lu_back_substitution

      use parameters_multibox, only: phi_pow, phi_bound

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
      logical, intent(in) :: has_elec, adiabatic_elec

      integer :: ia, it, iz, ix, iky, ind, b_solve
      complex, dimension(:, :), allocatable :: g0k, g1k, g0x, g1x, g0z
      complex, dimension(:), allocatable :: g_fsa, pb_fsa
      real :: tmp

      ia = 1
      b_solve = boundary_size - phi_bound

      allocate (g0k(1, nakx))
      allocate (g1k(1, nakx))
      allocate (g0x(1, x_fft_size))
      allocate (g1x(1, x_fft_size))
      if (adiabatic_elec .and. zonal_mode(1)) then
         allocate (g0z(x_fft_size, -nzgrid:nzgrid))
         allocate (g_fsa(x_fft_size - 2 * b_solve))
         allocate (pb_fsa(x_fft_size - 2 * b_solve))
      end if

      do it = 1, ntubes
         if (adiabatic_elec .and. zonal_mode(1)) pb_fsa = 0.0
         do iz = -nzgrid, nzgrid
            do iky = 1, naky
               if (iky > ky_solve_radial) then
                  phi(iky, :, iz, it) = phi(iky, :, iz, it) / gamtot(iky, :, iz)
               else
                  g0x = 0.0
                  tmp = 0
                  ind = boundary_size * (iky - 1 + naky * (iz + nzgrid + (2 * nzgrid + 1) * (it - 1)))
                  do ix = 1, b_solve
                     g0x(1, ix) = phi_buffer0(ind + ix)
                     g0x(1, x_fft_size + 1 - ix) = phi_buffer1(ind + boundary_size + 1 - ix)
                  end do

                  if (iky == 1) tmp = sum(real(g0x(1, :)))
                  call transform_x2kx(g0x, g0k)
                  if (phi_pow /= 0) then
                     g0k(1, :) = g0k(1, :) / ((zi * akx)**phi_pow)
                     if (iky == 1) g0k(1, 1) = 0.0
                  end if

                  g0x = 0.0
                  if (.not. has_elec .and. phi_pow /= 0) then
                     g1k(1, :) = g0k(1, :)
                     call transform_kx2x(g1k, g0x)

                     g0x(1, (b_solve + 1):(x_fft_size - b_solve)) = &
                        -g0x(1, (b_solve + 1):(x_fft_size - b_solve)) * b_mat

                     if (adiabatic_elec .and. iky == 1) then
                        pb_fsa = pb_fsa + (dl_over_b(ia, iz) + d_dl_over_b_drho(ia, iz) &
                                           * rho_mb_clamped((1 + b_solve):(x_fft_size - b_solve))) &
                                 * g0x(1, (1 + b_solve):(x_fft_size - b_solve))
                     end if
                  end if

                  g1k(1, :) = g0k(1, :) * gamtot(iky, :, iz)
                  call transform_kx2x(g1k, g1x)
                  g0x = g0x + g1x
                  g1k(1, :) = g0k(1, :) * dgamtotdr(iky, :, iz)
                  call transform_kx2x(g1k, g1x)
                  g1x(1, :) = rho_mb_clamped * g1x(1, :) + g0x(1, :)

                  g0k(1, :) = phi(iky, :, iz, it)
                  call transform_kx2x(g0k, g0x)
                  g0x = g0x - g1x

                  call lu_back_substitution(phi_solve(iky, iz)%zloc, phi_solve(iky, iz)%idx, &
                                            g0x(1, (b_solve + 1):(x_fft_size - b_solve)))

                  if (iky == 1) then
                     tmp = (tmp + &
                            sum(real(g0x(1, (b_solve + 1):(x_fft_size - b_solve))))) &
                           / (2 * b_solve)
                  end if

                  if (phi_pow /= 0) then
                     do ix = 1, b_solve
                        g0x(1, ix) = 0.0
                        g0x(1, x_fft_size + 1 - ix) = 0.0
                     end do
                     call transform_x2kx(g0x, g0k)
                     g0k(1, :) = g0k(1, :) * (zi * akx)**phi_pow
                     call transform_kx2x(g0k, g0x)
                  end if

                  do ix = 1, b_solve
                     g0x(1, ix) = phi_buffer0(ind + ix) !- tmp
                     g0x(1, x_fft_size + 1 - ix) = phi_buffer1(ind + boundary_size + 1 - ix) !- tmp
                  end do

                  call transform_x2kx(g0x, g0k)
                  if (phi_pow /= 0) then
                     g0k(1, :) = g0k(1, :) / (zi * akx)**phi_pow
                     if (iky == 1) g0k(1, 1) = 0.
                  end if
                  phi(iky, :, iz, it) = g0k(1, :)

               end if
            end do
         end do

         if (ky_solve_radial == 0 .and. any(gamtot(1, 1, :) < epsilon(0.))) phi(1, 1, :, it) = 0.0

         if (adiabatic_elec .and. zonal_mode(1)) then
            !get A_p^-1.(g - A_b.phi_b) in real space
            do iz = -nzgrid, nzgrid
               g0k(1, :) = phi(1, :, iz, it)
               call transform_kx2x(g0k, g0x)
               g0z(:, iz) = g0x(1, :)

               if (phi_pow /= 0) then
                  call lu_back_substitution(phi_solve(1, iz)%zloc, phi_solve(1, iz)%idx, pb_fsa)
                  g0z((1 + b_solve):(x_fft_size - b_solve), iz) = &
                     g0z((1 + b_solve):(x_fft_size - b_solve), iz) + pb_fsa
               end if
            end do

            ! get <A_p^-1.(g- - A_b.phi_b)>_psi
            g_fsa = 0.0
            do iz = -nzgrid, nzgrid
               g_fsa = g_fsa + (dl_over_b(ia, iz) + d_dl_over_b_drho(ia, iz) &
                                * rho_mb_clamped((1 + b_solve):(x_fft_size - b_solve))) &
                       * g0z((1 + b_solve):(x_fft_size - b_solve), iz)
            end do

            call lu_back_substitution(phizf_solve%zloc, phizf_solve%idx, g_fsa)
            g1x(1, (1 + b_solve):(x_fft_size - b_solve)) = b_mat * g_fsa

            do iz = -nzgrid, nzgrid
               g_fsa = g1x(1, (1 + b_solve):(x_fft_size - b_solve))
               call lu_back_substitution(phi_solve(1, iz)%zloc, phi_solve(1, iz)%idx, g_fsa)

               g0z((1 + b_solve):(x_fft_size - b_solve), iz) = &
                  g0z((1 + b_solve):(x_fft_size - b_solve), iz) + g_fsa

               g0x(1, :) = g0z(:, iz)
               call transform_x2kx(g0x, g0k)
               phi(1, :, iz, it) = g0k(1, :)
            end do
         end if
      end do

      deallocate (g0k, g1k, g0x, g1x)

      if (allocated(g0z)) deallocate (g0z)
      if (allocated(g_fsa)) deallocate (g_fsa)
      if (allocated(pb_fsa)) deallocate (pb_fsa)

   end subroutine mb_get_phi

!!>DSO - The following subroutines are the _ analogues of the ones found in
! stella_transforms.f90.
! The ones uses here ensure that the grid spacing in real space is consistent between
! domains (since we do not keep the checkboard mode)

   subroutine init_mb_transforms

      use stella_layouts, only: init_stella_layouts
      use parameters_kxky_grid, only: nakx, naky, naky_all

      implicit none

      if (mb_transforms_initialized) return
      mb_transforms_initialized = .true.

      if (.not. allocated(fft_kxky)) allocate (fft_kxky(naky, nakx))
      if (.not. allocated(fft_xky)) allocate (fft_xky(naky, x_fft_size))
      if (.not. allocated(fft_xy)) allocate (fft_xy(naky_all, x_fft_size))

      call init_x_fft
      call init_y_fft

   end subroutine init_mb_transforms

   subroutine init_x_fft

      use fft_work, only: init_ccfftw, FFT_BACKWARD, FFT_FORWARD

      implicit none

      if (.not. allocated(fft_x_k)) allocate (fft_x_k(x_fft_size))
      if (.not. allocated(fft_x_x)) allocate (fft_x_x(x_fft_size))

      call init_ccfftw(xf_fft, FFT_BACKWARD, x_fft_size, fft_x_k, fft_x_x)
      call init_ccfftw(xb_fft, FFT_FORWARD, x_fft_size, fft_x_x, fft_x_k)

   end subroutine init_x_fft

   subroutine init_y_fft

      use parameters_kxky_grid, only: naky, naky_all
      use fft_work, only: init_crfftw, init_rcfftw, FFT_BACKWARD, FFT_FORWARD

      implicit none

      if (.not. allocated(fft_y_k)) allocate (fft_y_k(naky))
      if (.not. allocated(fft_y_y)) allocate (fft_y_y(naky_all))

      call init_crfftw(yf_fft, FFT_BACKWARD, naky_all, fft_y_k, fft_y_y)
      call init_rcfftw(yb_fft, FFT_FORWARD, naky_all, fft_y_y, fft_y_k)

   end subroutine init_y_fft

!
!> transform routines start here
!

   subroutine transform_kx2x(gkx, gx)

      use parameters_kxky_grid, only: ikx_max

      implicit none

      complex, dimension(:, :), intent(in)  :: gkx
      complex, dimension(:, :), intent(out) :: gx

      integer :: iy

      do iy = 1, size(gkx, 1)
         fft_x_k = 0.
         fft_x_k(:ikx_max) = gkx(iy, :ikx_max)
         fft_x_k((x_fft_size - ikx_max + 2):) = gkx(iy, ikx_max + 1:)
         call dfftw_execute_dft(xf_fft%plan, fft_x_k, fft_x_x)
         gx(iy, :) = fft_x_x * xf_fft%scale
      end do

   end subroutine transform_kx2x

   subroutine transform_x2kx(gx, gkx)

      use parameters_kxky_grid, only: ikx_max

      implicit none

      complex, dimension(:, :), intent(in)  :: gx
      complex, dimension(:, :), intent(out) :: gkx

      integer :: iy

      do iy = 1, size(gx, 1)
         fft_x_x = gx(iy, :)
         call dfftw_execute_dft(xb_fft%plan, fft_x_x, fft_x_k)
         gkx(iy, :ikx_max) = fft_x_k(:ikx_max) * xb_fft%scale
         gkx(iy, ikx_max + 1:) = fft_x_k((x_fft_size - ikx_max + 2):) * xb_fft%scale
      end do

   end subroutine transform_x2kx

   subroutine transform_ky2y(gky, gy)

      implicit none

      complex, dimension(:, :), intent(in) :: gky
      real, dimension(:, :), intent(out) :: gy

      integer :: ikx

      do ikx = 1, size(gky, 2)
         fft_y_k = gky(:, ikx)
         call dfftw_execute_dft_c2r(yf_fft%plan, fft_y_k, fft_y_y)
         gy(:, ikx) = fft_y_y * yf_fft%scale
      end do

   end subroutine transform_ky2y

!   subroutine transform_y2ky (gy, gky)
!
!    implicit none
!
!    real, dimension (:,:), intent (in out) :: gy
!    complex, dimension (:,:), intent (out) :: gky
!
!    integer :: ikx
!
!    do ikx = 1, size(gy,2)
!       fft_y_k = gy(:,ikx)
!       call dfftw_execute_dft_r2c(yb_fft%plan, fft_y_y, fft_y_k)
!       gky(:,ikx) = fft_y_y*yb_fft%scale
!    end do
!
!  end subroutine transform_y2ky

   subroutine finish_mb_transforms

      implicit none
      if (.not. use_multibox) return
      call dfftw_destroy_plan(yf_fft%plan)
      call dfftw_destroy_plan(yb_fft%plan)
      call dfftw_destroy_plan(xf_fft%plan)
      call dfftw_destroy_plan(xb_fft%plan)
      if (allocated(fft_y_k)) deallocate (fft_y_k)
      if (allocated(fft_y_y)) deallocate (fft_y_y)
      if (allocated(fft_x_k)) deallocate (fft_x_k)
      if (allocated(fft_x_x)) deallocate (fft_x_x)
      if (allocated(fft_xky)) deallocate (fft_xky)
      if (allocated(fft_xy)) deallocate (fft_xy)

      mb_transforms_initialized = .false.

   end subroutine finish_mb_transforms

end module multibox
