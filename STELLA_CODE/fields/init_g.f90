!> This module contains the subroutines which set the initial value of the
!! fields and the distribution function.

module init_g

   implicit none
   
   ! Public routines
   public :: ginit, reset_init
   public :: init_init_g, finish_init_g
   
   ! The stella.f90 script will check if we want to call rescale_fields()
   public :: phiinit, scale_to_phiinit
   
   ! When we restart a simulation, we need to access <tstart>
   public :: tstart

   private
   
   !----------------------------- public variables -----------------------------

   ! Initialization parameters used in all options
   ! Moreover, these are used in rescale_fields() in the fields.fpp module
   logical :: scale_to_phiinit
   real :: phiinit
   
   ! When we restart a simulation, we need to access <tstart>
   real :: tstart
   
   !----------------------------- module variables -----------------------------

   ! Remember whether the module has already been initialised
   logical :: initialised = .false.

   ! Choose the initalization option for the potential
   integer :: init_distribution_switch
   
   ! This variable is read when the module is initialised and used later in stella_restore()
   real :: scale
   
   ! During the initialization of this module we set the restart path
   ! and we will parse it to stella_save.fpp through save_init()
   character(len=300) :: restart_file
   character(len=150) :: restart_dir

contains

   !****************************************************************************
   !                          INITIALIZE THIS MODULE                           !
   !****************************************************************************
   subroutine init_init_g

      use stella_save, only: init_save, read_many
      use stella_layouts, only: init_stella_layouts
      use system_fortran, only: systemf
      use mp, only: proc0, broadcast
      use stella_save, only: read_many
      
      ! Read namelist from input file
      use input_file_fields, only: read_namelist_initialise_distribution
      use input_file_fields, only: read_namelist_restart_options
      
      ! Load the <init_distribution_switch> parameters
      use input_file_fields, only: init_distribution_option_maxwellian
      use input_file_fields, only: init_distribution_option_noise
      use input_file_fields, only: init_distribution_option_restart_many
      use input_file_fields, only: init_distribution_option_kpar
      use input_file_fields, only: init_distribution_option_rh
      use input_file_fields, only: init_distribution_option_remap

      implicit none

      integer :: ind_slash
      
      !-------------------------------------------------------------------------

      if (initialised) return
      initialised = .true.

      call init_stella_layouts
      
      ! Read <initialise_distribution> namelist 
      if (proc0) call read_namelist_initialise_distribution(init_distribution_switch, phiinit, scale_to_phiinit)
         
      ! Broadcast to all processors
      call broadcast(init_distribution_switch)
      call broadcast(scale_to_phiinit)
      call broadcast(phiinit)

      ! Read <restart_options> namelist
      ! Most of these options will be parsed to other stella modules
      ! Except the <scale> variable which is used in init_distribution_switch = 'many'
      if (proc0) call read_namelist_restart_options(tstart, scale, restart_file, restart_dir, read_many)
         
      ! Broadcast to all processors
      call broadcast(tstart)
      call broadcast(scale)
      call broadcast(restart_file)
      call broadcast(restart_dir)
      call broadcast(read_many)

      ! Prepend restart_dir to restart_file, and append trailing slash if not exists
      if (restart_dir(len_trim(restart_dir):) /= "/") &
         restart_dir = trim(restart_dir)//"/" 
      if (proc0) call systemf('mkdir -p '//trim(restart_dir))

      ! Determine if restart file contains "/" if so split on this point to give DIR//FILE
      ! so restart files are created in DIR//restart_dir//FILE
      ind_slash = index(restart_file, "/", .true.)
      if (ind_slash == 0) then !No slash present
         restart_file = trim(restart_dir)//trim(restart_file)
      else !Slash present
         restart_file = trim(restart_file(1:ind_slash))//trim(restart_dir)//trim(restart_file(ind_slash + 1:))
      end if 
      
      ! Initialize the netcdf saving
      call init_save(restart_file)

   end subroutine init_init_g

   !****************************************************************************
   !                   INITIALIZE THE DISTRIBUTION FUNCTION                    !
   !****************************************************************************
   subroutine ginit(restarted, istep0)

      use stella_save, only: init_tstart
      use parameters_numerical, only: maxwellian_normalization
      
      ! Load the <init_distribution_switch> parameters
      use input_file_fields, only: init_distribution_option_maxwellian
      use input_file_fields, only: init_distribution_option_noise
      use input_file_fields, only: init_distribution_option_restart_many
      use input_file_fields, only: init_distribution_option_kpar
      use input_file_fields, only: init_distribution_option_rh
      use input_file_fields, only: init_distribution_option_remap

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

   !****************************************************************************
   !                     INITIALIZE POTENTIAL: MAXWELLIAN                      !
   !****************************************************************************
   subroutine ginit_maxwellian

      use mp, only: proc0, broadcast
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
      use input_file_fields, only: read_namelist_initialise_distribution_maxwellian

      implicit none

      complex, dimension(naky, nakx, -nzgrid:nzgrid) :: phi
      logical :: right
      integer :: ikxkyz
      integer :: iz, iky, ikx, is, ia
      
      ! Read the following variables from the input file
      real :: width0, den0, upar0
      logical :: oddparity, left, chop_side
      
      !-------------------------------------------------------------------------
      
      ! Read <initialise_distribution_maxwellian> namelist
      if (proc0) call read_namelist_initialise_distribution_maxwellian(width0, den0, upar0, oddparity, left, chop_side)
         
      ! Broadcast to all processors
      call broadcast(width0)
      call broadcast(den0)
      call broadcast(upar0)
      call broadcast(oddparity)
      call broadcast(left)
      call broadcast(chop_side)

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

   !****************************************************************************
   !                       INITIALIZE POTENTIAL: NOISE                         !
   !****************************************************************************
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
      use physics_parameters, only: include_nonlinear
      use ran
      use input_file_fields, only: read_namelist_initialise_distribution_noise

      implicit none

      complex, dimension(naky, nakx, -nzgrid:nzgrid, ntubes) :: phi
      real :: a, b, kmin
      integer :: ikxkyz, iz, it, iky, ikx, is, ie, iseg, ia
      integer :: itmod 
      
      ! Read the following variables from the input file
      real :: zf_init
      logical :: left, chop_side
      
      !-------------------------------------------------------------------------
      
      ! Read <initialise_distribution_noise> namelist
      if (proc0) call read_namelist_initialise_distribution_noise(zf_init, left, chop_side)
         
      ! Broadcast to all processors
      call broadcast(zf_init)
      call broadcast(left)
      call broadcast(chop_side)

      if ((naky == 1 .and. nakx == 1) .or. (.not. include_nonlinear)) then
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

   !****************************************************************************
   !                       INITIALIZE POTENTIAL: KPAR                          !
   !****************************************************************************
   subroutine ginit_kpar 
   
      use mp, only: proc0, broadcast
      use zgrid, only: nzgrid, zed
      use parameters_kxky_grids, only: naky, nakx
      use grids_kxky, only: theta0
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: vpa, vperp2
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use arrays_dist_fn, only: gvmu
      use stella_layouts, only: kxkyz_lo, iky_idx, ikx_idx, iz_idx, is_idx
      use input_file_fields, only: read_namelist_initialise_distribution_kpar
      use constants, only: zi

      implicit none

      complex, dimension(naky, nakx, -nzgrid:nzgrid) :: phi, odd
      real, dimension(-nzgrid:nzgrid) :: dfac, ufac, tparfac, tperpfac
      integer :: ikxkyz
      integer :: iz, iky, ikx, imu, iv, ia, is
      
      ! Read the following variables from the input file
      real :: width0, imfac, refac
      real :: den0, upar0, tpar0, tperp0
      real :: den1, upar1, tpar1, tperp1
      real :: den2, upar2, tpar2, tperp2
      logical :: left, chop_side
      
      !-------------------------------------------------------------------------
      
      ! Read <initialise_distribution_noise> namelist
      if (proc0) call read_namelist_initialise_distribution_kpar(&
         width0, refac, imfac, den0, upar0, tpar0, tperp0, &
         den1, upar1, tpar1, tperp1, den2, upar2, tpar2, tperp2, left, chop_side)
         
      ! Broadcast to all processors
      call broadcast(width0)
      call broadcast(refac)
      call broadcast(imfac)
      call broadcast(den0)
      call broadcast(upar0)
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
      call broadcast(left)
      call broadcast(chop_side)

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

   !****************************************************************************
   !                        INITIALIZE POTENTIAL: RH                           !
   !****************************************************************************
   subroutine ginit_rh

      use mp, only: proc0, broadcast
      use species, only: spec
      use arrays_dist_fn, only: gvmu, kperp2
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: nvpa, nmu
      use grids_kxky, only: akx
      use input_file_fields, only: read_namelist_initialise_distribution_rh

      implicit none

      integer :: ikxkyz, iky, ikx, iz, is, ia
      
      ! Read the following variables from the input file
      real :: imfac, refac
      real :: kxmax, kxmin
      
      !-------------------------------------------------------------------------
      
      ! Read <initialise_distribution_rh> namelist
      if (proc0) call read_namelist_initialise_distribution_rh(kxmin, kxmax, imfac, refac)
      
      ! Broadcast to all processors
      call broadcast(refac)
      call broadcast(imfac)
      call broadcast(kxmax)
      call broadcast(kxmin)

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

   !****************************************************************************
   !                      INITIALIZE POTENTIAL: REMAP                          !
   !****************************************************************************
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

   !****************************************************************************
   !                       INITIALIZE POTENTIAL: MANY                          !
   !****************************************************************************
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

   !****************************************************************************
   !                        NORMALIZE BY MAXWELLIAN                            !
   !****************************************************************************
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

   !****************************************************************************
   !                                                                           !
   !****************************************************************************
   subroutine reset_init
   
      use input_file_fields, only: init_distribution_option_restart_many
      
      implicit none

      init_distribution_switch = init_distribution_option_restart_many

   end subroutine reset_init

   !****************************************************************************
   !                                                                           !
   !****************************************************************************
   subroutine finish_init_g

      use stella_save, only: finish_save

      implicit none

      initialised = .false.

      call finish_save

   end subroutine finish_init_g

end module init_g
