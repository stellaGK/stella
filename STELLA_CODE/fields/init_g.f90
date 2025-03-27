!> This module contains the subroutines which set the initial value of the
!! fields and the distribution function.

module init_g
   implicit none

   public :: ginit
   public :: init_init_g, finish_init_g
   public :: width0
   public :: scale_to_phiinit, phiinit
   public :: tstart
   public :: reset_init

   private

   ! Choose the initalization option for the potential
   integer :: init_distribution_switch

   real :: width0, phiinit, imfac, refac, zf_init
   real :: den0, upar0, tpar0, tperp0
   real :: den1, upar1, tpar1, tperp1
   real :: den2, upar2, tpar2, tperp2
   real :: tstart, scale, kxmax, kxmin
   logical :: chop_side, left, scale_to_phiinit, oddparity
   character(300), public :: restart_file
   character(len=150) :: restart_dir

   logical :: initialized = .false.
   logical :: exist

contains

   subroutine init_init_g

      use stella_save, only: init_save, read_many
      use stella_layouts, only: init_stella_layouts
      use system_fortran, only: systemf
      use mp, only: proc0, broadcast
      
      ! Read namelist from input file
      use input_file, only: read_namelist_initialize_distribution
      use input_file, only: read_namelist_initialize_distribution_maxwellian
      
      ! Load the <init_distribution_switch> parameters
      use input_file, only: init_distribution_option_maxwellian
      use input_file, only: init_distribution_option_noise
      use input_file, only: init_distribution_option_restart_many
      use input_file, only: init_distribution_option_kpar
      use input_file, only: init_distribution_option_rh
      use input_file, only: init_distribution_option_remap

      implicit none

      integer :: ind_slash

      if (initialized) return
      initialized = .true.

      call init_stella_layouts
      
      ! Read <initialize_distribution> namelist
      if (proc0) call read_namelist_initialize_distribution(init_distribution_switch, &
         phiinit, left, chop_side, scale_to_phiinit)
         
      ! Broadcast to all processors
      call broadcast(init_distribution_switch)
      call broadcast(phiinit)
      call broadcast(left)
      call broadcast(chop_side)
      call broadcast(scale_to_phiinit)
      
      ! Read <initialize_distribution_maxwellian> namelist
      if (init_distribution_switch==init_distribution_option_maxwellian) then
         if (proc0) call read_namelist_initialize_distribution_maxwellian(width0, den0, upar0, oddparity)
         call broadcast(width0)
         call broadcast(den0)
         call broadcast(upar0)
         call broadcast(oddparity)
      end if

      if (proc0) call read_parameters

      ! prepend restart_dir to restart_file
      ! append trailing slash if not exists
      if (restart_dir(len_trim(restart_dir):) /= "/") &
         restart_dir = trim(restart_dir)//"/"

      if (proc0) call systemf('mkdir -p '//trim(restart_dir))

      !Determine if restart file contains "/" if so split on this point to give DIR//FILE
      !so restart files are created in DIR//restart_dir//FILE
      ind_slash = index(restart_file, "/", .true.)
      if (ind_slash == 0) then !No slash present
         restart_file = trim(restart_dir)//trim(restart_file)
      else !Slash present
         restart_file = trim(restart_file(1:ind_slash))//trim(restart_dir)//trim(restart_file(ind_slash + 1:))
      end if

      call broadcast(refac)
      call broadcast(imfac)
      call broadcast(tpar0)
      call broadcast(tperp0)
      call broadcast(den1)
      call broadcast(upar1)
      call broadcast(tpar1)
      call broadcast(tperp1)
      call broadcast(den2)
      call broadcast(upar2)
      call broadcast(tpar2)
      call broadcast(tperp2)
      call broadcast(zf_init)
      call broadcast(kxmax)
      call broadcast(kxmin)
      call broadcast(tstart)
      call broadcast(restart_file)
      call broadcast(read_many)
      call broadcast(scale)

      call init_save(restart_file)

   end subroutine init_init_g

   subroutine ginit(restarted, istep0)

      use stella_save, only: init_tstart
      use parameters_numerical, only: maxwellian_normalization
      
      ! Load the <init_distribution_switch> parameters
      use input_file, only: init_distribution_option_maxwellian
      use input_file, only: init_distribution_option_noise
      use input_file, only: init_distribution_option_restart_many
      use input_file, only: init_distribution_option_kpar
      use input_file, only: init_distribution_option_rh
      use input_file, only: init_distribution_option_remap

      logical, intent(out) :: restarted
      integer, intent(out) :: istep0
      integer :: istatus

      restarted = .false.
      istep0 = 0
      select case (init_distribution_switch)
      case (init_distribution_option_maxwellian)
         call ginit_maxwellian
      case (init_distribution_option_noise)
         call ginit_noise
      case (init_distribution_option_kpar)
         call ginit_kpar
      case (init_distribution_option_rh)
         call ginit_rh
      case (init_distribution_option_remap)
         call ginit_remap
      case (init_distribution_option_restart_many)
         call ginit_restart_many
         call init_tstart(tstart, istep0, istatus)
         restarted = .true.
         scale = 1.
      end select

      !> if maxwwellian_normalization = .true., the pdf is normalized by F0 (which is not the case otherwise)
      !> unless reading in g from a restart file, normalise g by F0 for a full flux surface simulation
      if (maxwellian_normalization .and. init_distribution_switch /= init_distribution_option_restart_many) then
         call normalize_by_maxwellian
      end if

   end subroutine ginit

   subroutine read_parameters
      use file_utils, only: input_unit, error_unit, run_name, input_unit_exist
      use stella_save, only: read_many

      implicit none

      namelist /init_g_knobs/ width0, &
         restart_file, restart_dir, read_many, scale, tstart, zf_init, &
         den0, upar0, tpar0, tperp0, imfac, refac, &
         den1, upar1, tpar1, tperp1, &
         den2, upar2, tpar2, tperp2, &
         kxmax, kxmin
      integer :: in_file

      width0 = -3.5 ! Used for <ginit_options> = {default, kpar}
      den0 = 1. ! Used for <ginit_options> = {default, kpar}
      upar0 = 0. ! Used for <ginit_options> = {default, kpar}
      tstart = 0. ! Used for restarted simulations
      scale = 1.0 ! Used for restarted simulations
      refac = 1. ! Used for <ginit_options> = {kpar}
      imfac = 0. ! Used for <ginit_options> = {kpar}
      tpar0 = 0. ! Used for <ginit_options> = {kpar}
      tperp0 = 0. ! Used for <ginit_options> = {kpar}
      den1 = 0. ! Used for <ginit_options> = {kpar}
      upar1 = 0. ! Used for <ginit_options> = {kpar}	
      tpar1 = 0. ! Used for <ginit_options> = {kpar}
      tperp1 = 0. ! Used for <ginit_options> = {kpar}
      den2 = 0. ! Used for <ginit_options> = {kpar}
      upar2 = 0. ! Used for <ginit_options> = {kpar}
      tpar2 = 0. ! Used for <ginit_options> = {kpar}
      tperp2 = 0. ! Used for <ginit_options> = {kpar}
      zf_init = 1.0 ! Used for <ginit_options> = {noise}
      kxmax = 1.e100 ! Used for <ginit_options> = {rh}
      kxmin = 0. ! Used for <ginit_options> = {rh}

      restart_file = trim(run_name)//".nc"
      restart_dir = "./"
      in_file = input_unit_exist("init_g_knobs", exist)
!    if (exist) read (unit=input_unit("init_g_knobs"), nml=init_g_knobs)
      if (exist) read (unit=in_file, nml=init_g_knobs)

   end subroutine read_parameters

   subroutine ginit_maxwellian

      use constants, only: zi
      use species, only: spec
      use zgrid, only: nzgrid, zed
      use parameters_kxky_grids, only: naky, nakx, ikx_max
      use grids_kxky, only: theta0, akx, zonal_mode
      use parameters_kxky_grids, only: reality
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: vpa
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use arrays_dist_fn, only: gvmu
      use stella_layouts, only: kxkyz_lo, iz_idx, ikx_idx, iky_idx, is_idx
      use ran, only: ranf

      implicit none

      complex, dimension(naky, nakx, -nzgrid:nzgrid) :: phi
      logical :: right
      integer :: ikxkyz
      integer :: iz, iky, ikx, is, ia

      right = .not. left

      do iz = -nzgrid, nzgrid
         phi(:, :, iz) = exp(-((zed(iz) - theta0) / width0)**2) * cmplx(1.0, 1.0)
      end do

      ! this is a messy way of doing things
      ! could tidy it up a bit
      if (sum(cabs(phi)) < epsilon(0.)) then
         do iz = -nzgrid, nzgrid
            phi(:, :, iz) = exp(-(zed(iz) / width0)**2) * cmplx(1.0, 1.0)
         end do
      end if

      if (oddparity) then
      ! make phi an odd function of zed
        do iz = -nzgrid,nzgrid
            phi(:, :, iz) = zed(iz)*phi(:, :, iz)
        end do
      end if

      if (chop_side) then
         if (left) phi(:, :, :-1) = 0.0
         if (right) phi(:, :, 1:) = 0.0
      end if

      if (zonal_mode(1)) then
         ! zero out kx = ky = 0 mode
         if (abs(akx(1)) < epsilon(0.0)) then
            phi(1, 1, :) = 0.0
         end if

         ! force the reality condition on the zonal mode; i.e. phi_(-kx,ky=0) = conjugate(phi_(kx,ky=0))
         if (reality) then
            ! ikx_max is the index corresponding to max k_x value
            do ikx = 1, nakx - ikx_max
               phi(1, nakx - ikx + 1, :) = conjg(phi(1, ikx + 1, :))
            end do
         end if
      end if

      ! need better way to initialise for full flux surface cases
      ia = 1

      gvmu = 0.
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iz = iz_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iky = iky_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         gvmu(:, :, ikxkyz) = phiinit * phi(iky, ikx, iz) / abs(spec(is)%z) &
                              * (den0 + 2.0 * zi * spread(vpa, 2, nmu) * upar0) &
                              * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * spread(maxwell_vpa(:, is), 2, nmu) * maxwell_fac(is)
      end do

   end subroutine ginit_maxwellian

   subroutine ginit_noise

      use mp, only: proc0, broadcast
      use arrays_dist_fn, only: kperp2
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use extended_zgrid, only: ikxmod, nsegments, neigen
      use extended_zgrid, only: it_right
      use extended_zgrid, only: periodic, phase_shift
      use parameters_kxky_grids, only: naky, nakx, reality
      use grids_kxky, only: zonal_mode
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use arrays_dist_fn, only: gvmu
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use mp, only: proc0, broadcast, max_allreduce
      use mp, only: scope, crossdomprocs, subprocs
      use file_utils, only: runtype_option_switch, runtype_multibox
      use parameters_physics, only: nonlinear 
      use ran

      implicit none

      complex, dimension(naky, nakx, -nzgrid:nzgrid, ntubes) :: phi
      real :: a, b, kmin
      integer :: ikxkyz, iz, it, iky, ikx, is, ie, iseg, ia
      integer :: itmod

      if ((naky == 1 .and. nakx == 1) .or. (.not. nonlinear)) then
         if (proc0) then
            write (*, *) 'Noise initialization option is not suited for single mode simulations,'
            write (*, *) 'or linear simulations, using default initialization option instead.'
            write (*, *)
         end if
         call ginit_maxwellian
         return
      else
         ! zero out ky=kx=0 mode
         phi(1, 1, :, :) = 0.0
      end if

      ia = 1
      if (proc0) then
         phi(1, 1, :, :) = 0.0
         kmin = 1.e6
         if (naky > 1) kmin = minval(kperp2(2, 1, ia, :))
         if (nakx > 1) kmin = min(kmin, minval(kperp2(1, 2, ia, :)))

         if (runtype_option_switch == runtype_multibox) then
            call scope(crossdomprocs)
            call max_allreduce(kmin)
            call scope(subprocs)
         end if

         ! keep old (ikx, iky) loop order to get old results exactly:
         !Fill phi with random (complex) numbers between -0.5 and 0.5
         do ikx = 1, nakx
            do iky = 1, naky
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid

                     a = ranf() - 0.5
                     b = ranf() - 0.5
                     ! do not populate high k modes with large amplitudes
                     if ((ikx > 1 .or. iky > 1) .and. (kperp2(iky, ikx, ia, iz) >= kmin)) then
                        !the following as an extra factor of kmin to offset the Gamma-1 in quasineutrality
                        phi(iky, ikx, iz, it) = cmplx(a, b) * kmin * kmin / kperp2(iky, ikx, ia, iz)
                     else
                        phi(iky, ikx, iz, it) = 0.0
                     end if
                  end do
                  if (chop_side) then
                     if (left) then
                        phi(iky, ikx, :-1, it) = 0.0
                     else
                        phi(iky, ikx, 1:, it) = 0.0
                     end if
                  end if
               end do
            end do
         end do

         ! enforce periodicity where required
         do iky = 1, naky
            if (periodic(iky)) then
               phi(iky, :, nzgrid, :) = phi(iky, :, -nzgrid, :) / phase_shift(iky)
            end if
         end do

         ! zero out the kx=ky=0 mode and apply optional
         ! scaliing factor to all zonal modes
         if (zonal_mode(1)) then
            !Apply scaling factor
            phi(1, :, :, :) = phi(1, :, :, :) * zf_init

            !Set ky=kx=0.0 mode to zero in amplitude
            phi(1, 1, :, :) = 0.0
         end if

         !Apply reality condition (i.e. -kx mode is conjugate of +kx mode)
         if (reality) then
            do ikx = nakx / 2 + 2, nakx
               phi(1, ikx, :, :) = conjg(phi(1, nakx - ikx + 2, :, :))
            end do
         end if

      end if

      do iky = 1, naky
         do ie = 1, neigen(iky)
            ! enforce zero BC at ends of domain, unless periodic
            if (.not. periodic(iky)) then
               phi(iky, ikxmod(1, ie, iky), -nzgrid, :) = 0.0
               phi(iky, ikxmod(nsegments(ie, iky), ie, iky), nzgrid, :) = 0.0
            end if
            ! enforce equality of g values at duplicate zed points
            if (nsegments(ie, iky) > 1) then
               do it = 1, ntubes
                  itmod = it
                  do iseg = 2, nsegments(ie, iky)
                     phi(iky, ikxmod(iseg, ie, iky), -nzgrid, it_right(itmod)) = phi(iky, ikxmod(iseg - 1, ie, iky), nzgrid, itmod)
                     itmod = it_right(itmod)
                  end do
               end do
            end if
         end do
      end do

      call broadcast(phi)

      !Now set g using data in phi
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iky = iky_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         gvmu(:, :, ikxkyz) = spec(is)%z * phiinit * phi(iky, ikx, iz, it) &
                              * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
      end do

   end subroutine ginit_noise

   subroutine ginit_kpar

!    use species, only: spec, has_electron_species
      use zgrid, only: nzgrid, zed
      use parameters_kxky_grids, only: naky, nakx
      use grids_kxky, only: theta0
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: vpa, vperp2
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use arrays_dist_fn, only: gvmu
      use stella_layouts, only: kxkyz_lo, iky_idx, ikx_idx, iz_idx, is_idx
      use constants, only: zi

      implicit none

      complex, dimension(naky, nakx, -nzgrid:nzgrid) :: phi, odd
      real, dimension(-nzgrid:nzgrid) :: dfac, ufac, tparfac, tperpfac
      integer :: ikxkyz
      integer :: iz, iky, ikx, imu, iv, ia, is

      phi = 0.
      odd = 0.
      if (width0 > 0.) then
         do iz = -nzgrid, nzgrid
            phi(:, :, iz) = exp(-((zed(iz) - theta0) / width0)**2) * cmplx(refac, imfac)
         end do
      else
         do iz = -nzgrid, nzgrid
            phi(:, :, iz) = cmplx(refac, imfac)
         end do
      end if
      if (chop_side) then
         if (left) then
            phi(:, :, :-1) = 0.0
         else
            phi(:, :, 1:) = 0.0
         end if
      end if

      odd = zi * phi

      dfac = den0 + den1 * cos(zed) + den2 * cos(2.*zed)
      ufac = upar0 + upar1 * sin(zed) + upar2 * sin(2.*zed)
      tparfac = tpar0 + tpar1 * cos(zed) + tpar2 * cos(2.*zed)
      tperpfac = tperp0 + tperp1 * cos(zed) + tperp2 * cos(2.*zed)

      ia = 1
      ! charge dependence keeps initial Phi from being too small
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            do iv = 1, nvpa
               gvmu(iv, imu, ikxkyz) = phiinit &
                                       * (dfac(iz) * phi(iky, ikx, iz) &
                                          + 2.0 * vpa(iv) * ufac(iz) * odd(iky, ikx, iz) &
                                          + (vpa(iv)**2 - 0.5) * tparfac(iz) * phi(iky, ikx, iz) &
                                          + tperpfac(iz) * (vperp2(ia, iz, imu) - 1.) * phi(iky, ikx, iz))
            end do
         end do
         gvmu(:, :, ikxkyz) = gvmu(:, :, ikxkyz) &
                              * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
      end do

! FLAG -- should be uncommented, which means I need to fix flae
!    if (has_electron_species(spec)) then
!       call flae (gold, gnew)
!       gold = gold - gnew
!    end if
!    gnew = gold

   end subroutine ginit_kpar

   subroutine ginit_rh

      use species, only: spec
      use arrays_dist_fn, only: gvmu, kperp2
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: nvpa, nmu
      use grids_kxky, only: akx

      implicit none

      integer :: ikxkyz, iky, ikx, iz, is, ia

      ! initialize g to be a Maxwellian with a constant density perturbation

      gvmu = 0.

      ia = 1
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         ! only set the first ky mode to be non-zero
         ! this is because this is meant to test the damping of zonal flow (ky=0)
         iky = iky_idx(kxkyz_lo, ikxkyz); if (iky /= 1) cycle
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)

         if (abs(akx(ikx)) < kxmax .and. abs(akx(ikx)) > kxmin) then
            gvmu(:, :, ikxkyz) = spec(is)%z * 0.5 * phiinit * kperp2(iky, ikx, ia, iz) &
                                 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
         end if
      end do

   end subroutine ginit_rh

   subroutine ginit_remap

      use species, only: spec
      use arrays_dist_fn, only: gvmu
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: nvpa, nmu

      implicit none

      integer :: ikxkyz, iky, ikx, iz, is, ia

      ! initialize g to be a Maxwellian with a constant density perturbation

      gvmu = 0.

      ia = 1
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)

         !if((ikx.eq.15.and.iky.eq.5).or.((ikx-nakx).eq.-12.and.iky.eq.3)) then
         if ((ikx == 1 .and. iky == 2)) then
            gvmu(:, :, ikxkyz) = spec(is)%z * phiinit &
                                 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
         end if
      end do

   end subroutine ginit_remap

   subroutine ginit_restart_many

      use arrays_dist_fn, only: gvmu
      use stella_save, only: stella_restore
      use mp, only: proc0
      use file_utils, only: error_unit

      implicit none

      integer :: istatus, ierr

      ! should really check if profile_variation=T here but need
      ! to move profile_variation to module that is accessible here
      call stella_restore(gvmu, scale, istatus)

      if (istatus /= 0) then
         ierr = error_unit()
         if (proc0) write (ierr, *) "Error reading file: ", trim(restart_file)
         gvmu = 0.
      end if

   end subroutine ginit_restart_many

   subroutine normalize_by_maxwellian

      use stella_layouts, only: kxkyz_lo, is_idx, iz_idx
      use arrays_dist_fn, only: gvmu
      use vpamu_grids, only: nmu
      use vpamu_grids, only: maxwell_mu, maxwell_vpa, maxwell_fac

      implicit none

      integer :: ia, imu
      integer :: ikxkyz, iz, is

      !> gvmu is initialised with a Maxwellian weighting for flux tube simulations,
      !> with the Maxwellian evaluated at ia = 1
      !> we are undoing that weighting here, so also need to use ia = 1
      ia = 1
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            gvmu(:, imu, ikxkyz) = gvmu(:, imu, ikxkyz) / (maxwell_mu(ia, iz, imu, is) * maxwell_vpa(:, is) * maxwell_fac(is))
         end do
      end do

   end subroutine normalize_by_maxwellian

   subroutine reset_init
   
      use input_file, only: init_distribution_option_restart_many
      
      implicit none

      init_distribution_switch = init_distribution_option_restart_many

   end subroutine reset_init

!   subroutine flae (g, gavg)

!     use species, only: spec, electron_species
!     use zgrid, only: nzgrid, delthet, jacob
!     use kt_grids, only: aky, ntheta0
!     use vpamu_grids, only: nvgrid
!     use stella_layouts, only: gxyz_lo, is_idx
!     complex, dimension (-nzgrid:,:,:,gxyz_lo%llim_proc:), intent (in) :: g
!     complex, dimension (-nzgrid:,:,:,gxyz_lo%llim_proc:), intent (out) :: gavg

!     real :: wgt
!     integer :: iglo, it, ik

!     gavg = 0.
!     wgt = 1./sum(delthet*jacob)

!     do iglo = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
!        if (spec(is_idx(gxyz_lo, iglo))%type /= electron_species) cycle
!        ik = 1
!        if (aky(ik) > epsilon(0.)) cycle
!        do it = 1, ntheta0
!           gavg(:,it,ik,iglo) = sum(g(:,it,ik,iglo)*delthet*jacob)*wgt
!        end do
!     end do

!   end subroutine flae

   subroutine finish_init_g

      use stella_save, only: finish_save

      implicit none

      initialized = .false.

      call finish_save

   end subroutine finish_init_g

end module init_g
